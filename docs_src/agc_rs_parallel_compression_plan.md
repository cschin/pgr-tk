# Implementation Plan: Multi-threaded Compression for agc-rs

**Date:** 2026-04-04
**Branch:** `aln_map_improvment`
**File:** `agc-rs/src/compressor.rs` (primary), `agc-rs/src/lz_diff.rs` (minor)

---

## Current State

`rayon` is already a dependency and `par_iter()` is already imported.
Two places already use parallel iteration:

| Location | What is parallel | What is not |
|---|---|---|
| `add_as_reference` line 641 | `split_contig + compress_reference` per record | k-mer collection, splitter scan, DB writes |
| `add_as_delta` Phase 2 line 767 | per-record split + group lookup + `lz_from_ref_blob` | Phase 1 DB reads, Phase 3 ZSTD+DB writes |

The remaining single-threaded work accounts for most wall-clock time on
human-genome scale inputs.

---

## Bottleneck Analysis

### A — `determine_splitters` (reference path, one-time cost)

Called once for the first sample. Three sequential passes over all records:

1. `collect_singleton_kmers` — counts every k-mer across all records  
2. `find_splitters_in_contigs` — scans each record for singleton k-mers

Both are embarrassingly parallel per record. Currently single-threaded.

### B — `build_fallback_index` + `collect_ref_kmers` (delta path, Phase 1)

Both functions:
1. Issue a single `SELECT … FROM segment_group` query (must be serial)
2. For each row: `decompress_reference` (ZSTD) + scan k-mers (CPU-bound)

Step 2 is independent per row and can be parallelized. For a reference
with ~50 k segment groups, this is the dominant Phase 1 cost.

### C — `lz_from_ref_blob` called per segment in Phase 2 (biggest win)

Current code (line 838–840):
```rust
let mut lz = segment::lz_from_ref_blob(ref_blob, params_ref)?;
let raw_delta = segment::compress_delta(&mut lz, &seg.seq)?;
```

`lz_from_ref_blob` does:
1. ZSTD-decompress the ref blob (~60 KB after decompression)
2. `LzDiff::prepare()` — build hash table over the reference O(ref_len)

This is called **once per segment**, even when many segments map to the
same group.  For a 3 GB genome with 50 k segments mapping to ~10 k
groups, the LZ index is rebuilt ~5× per group on average, entirely
redundant.

### D — `batch_compress_deltas` inside the serial transaction (Phase 3)

Current code (line 998–1024): for each group with new deltas, inside the
transaction:
1. Load existing `delta_blob`
2. `extract_delta_from_batch` for each existing entry (ZSTD decompress + walk)
3. Append new raw deltas
4. `batch_compress_deltas` (ZSTD compress at level 19)

ZSTD at level 19 is expensive.  Doing this inside the transaction blocks
all writes until every group's blob is recomputed.

---

## Proposed Changes

### Change 1 — Parallelize `collect_singleton_kmers`

**File:** `compressor.rs` lines 36–52

Replace the sequential HashMap accumulation with rayon fold+reduce:

```rust
fn collect_singleton_kmers(records: &[FastaRecord], k: usize) -> HashSet<u64> {
    records
        .par_iter()
        .fold(
            || HashMap::<u64, u32>::new(),
            |mut counts, rec| {
                let mut km = Kmer::new(k as u8);
                for &b in &rec.seq {
                    if km.push_bits(b) && km.full() {
                        *counts.entry(km.canonical()).or_insert(0) =
                            counts.get(&km.canonical()).copied().unwrap_or(0)
                                .saturating_add(1);
                    }
                }
                counts
            },
        )
        .reduce(
            || HashMap::new(),
            |mut a, b| {
                for (k, v) in b {
                    *a.entry(k).or_insert(0) =
                        a.get(&k).copied().unwrap_or(0).saturating_add(v);
                }
                a
            },
        )
        .into_iter()
        .filter(|(_, cnt)| *cnt == 1)
        .map(|(kmer, _)| kmer)
        .collect()
}
```

**Note on correctness:** fold+reduce sums counts across threads.
A k-mer that appears exactly once across all records is still a singleton.
The semantics are identical to the sequential version.

### Change 2 — Parallelize `find_splitters_in_contigs`

**File:** `compressor.rs` lines 107–118

```rust
fn find_splitters_in_contigs(
    records: &[FastaRecord],
    singletons: &HashSet<u64>,
    k: usize,
    segment_size: usize,
) -> HashSet<u64> {
    records
        .par_iter()
        .map(|rec| find_splitters_in_seq(&rec.seq, singletons, k, segment_size))
        .reduce(HashSet::new, |mut a, b| { a.extend(b); a })
}
```

`find_splitters_in_seq` is pure (read-only `singletons`).  The union of
per-record splitter sets equals the result of the sequential version.

### Change 3 — Parallelize `build_fallback_index`

**File:** `compressor.rs` lines 464–497

Load rows from DB serially (unchanged), then process in parallel:

```rust
fn build_fallback_index(
    conn: &Connection,
    kmer_len: usize,
    sample_rate: u64,
    max_bucket: usize,
) -> Result<(HashMap<u64, Vec<i64>>, HashMap<i64, Vec<u8>>, HashMap<i64, usize>)> {
    // Serial DB load (unchanged)
    let rows: Vec<(i64, Vec<u8>)> = { ... };

    // Parallel processing: decompress + sample k-mers per group
    type PerGroup = (i64, Vec<u64>, Vec<u8>, usize);
    let processed: Vec<PerGroup> = rows
        .par_iter()
        .map(|(group_id, ref_blob)| {
            let seq_2bit = segment::decompress_reference(ref_blob)
                .expect("decompress_reference");
            let kmers = sample_kmers_2bit(&seq_2bit, kmer_len, sample_rate);
            (*group_id, kmers, ref_blob.clone(), seq_2bit.len())
        })
        .collect();

    // Serial merge into output maps
    let mut kmer_index: HashMap<u64, Vec<i64>> = HashMap::new();
    let mut ref_data_map: HashMap<i64, Vec<u8>> = HashMap::with_capacity(processed.len());
    let mut ref_len_map: HashMap<i64, usize> = HashMap::with_capacity(processed.len());
    for (group_id, kmers, blob, len) in processed {
        for kmer in kmers {
            let bucket = kmer_index.entry(kmer).or_default();
            if bucket.len() < max_bucket {
                bucket.push(group_id);
            }
        }
        ref_data_map.insert(group_id, blob);
        ref_len_map.insert(group_id, len);
    }

    Ok((kmer_index, ref_data_map, ref_len_map))
}
```

Apply the same parallel pattern to `collect_ref_kmers` (lines 180–205):
load blobs from DB serially, then `par_iter()` to decompress + scan k-mers,
finally reduce into a single `HashSet<u64>`.

### Change 4 — Pre-build LzDiff cache before Phase 2 (biggest win)

**File:** `compressor.rs`, `add_as_delta`, before the `par_iter` at line 767

Instead of calling `lz_from_ref_blob` once per segment (redundantly for
the same group), build every needed LzDiff instance **once** before Phase 2,
in parallel, and share via `Mutex`:

```rust
// After build_fallback_index returns ref_data_map:
use std::sync::Mutex;

let lz_cache: HashMap<i64, Mutex<LzDiff>> = ref_data_map
    .par_iter()
    .map(|(&group_id, blob)| {
        let lz = segment::lz_from_ref_blob(blob, params_ref)
            .expect("lz_from_ref_blob");
        (group_id, Mutex::new(lz))
    })
    .collect();
```

Then in Phase 2, replace lines 838–840:

```rust
// Before (per-segment rebuild):
let mut lz = segment::lz_from_ref_blob(ref_blob, params_ref)?;
let raw_delta = segment::compress_delta(&mut lz, &seg.seq)?;

// After (cached, lock per group):
let raw_delta = {
    let mut lz = lz_cache[&group_id].lock().expect("lz_cache lock");
    segment::compress_delta(&mut lz, &seg.seq)?
};
```

**Contention analysis:** The lock is held only during `encode()`, which
runs O(seg_len) ≈ 60 KB.  Two threads encoding against the same group
will briefly serialize.  In practice, with ~50 k segments mapping to
~10 k groups and 16 threads, contention is low.

**Alternative (no locking):** If `LzDiff::encode_v1/v2` does not actually
mutate struct fields (the hash tables are read-only after `prepare()`),
change `encode(&mut self, ...)` to `encode(&self, ...)` and store LzDiff
in `Arc<LzDiff>` directly — zero contention.  This requires verifying
that `encode_v1`/`encode_v2` have no side effects on `self`.  If they
don't, add `unsafe impl Sync for LzDiff {}` or restructure the
interior-mutability (move `index_ready` into a `Cell<bool>`).

### Change 5 — Move ZSTD compression out of the transaction (Phase 3)

**File:** `compressor.rs` lines 998–1024

Currently the expensive ZSTD level-19 compression runs inside the SQLite
transaction, blocking all writes.  Move it before the transaction:

```rust
// ── Step A: load existing delta_blobs BEFORE the transaction ─────────────
let existing_blobs: HashMap<i64, Option<Vec<u8>>> = new_deltas_per_group
    .keys()
    .map(|gid| {
        let blob: Option<Vec<u8>> = self.db.conn().query_row(
            "SELECT delta_blob FROM segment_group WHERE id = ?1",
            params![gid],
            |r| r.get(0),
        ).unwrap_or(None);
        (*gid, blob)
    })
    .collect();

// ── Step B: parallel ZSTD recompression ──────────────────────────────────
let combined_blobs: HashMap<i64, Vec<u8>> = new_deltas_per_group
    .par_iter()
    .map(|(group_id, new_raws)| {
        let mut all_raws: Vec<Vec<u8>> = Vec::new();
        if let Some(Some(blob)) = existing_blobs.get(group_id) {
            let existing_count =
                existing_counts.get(group_id).copied().unwrap_or(0) as usize;
            for idx in 0..existing_count {
                let raw = segment::extract_delta_from_batch(blob, idx)
                    .expect("extract_delta_from_batch");
                all_raws.push(raw);
            }
        }
        all_raws.extend_from_slice(new_raws);
        let blob = segment::batch_compress_deltas(&all_raws)
            .expect("batch_compress_deltas");
        (*group_id, blob)
    })
    .collect();

// ── Step C: inside transaction — fast UPDATE only ─────────────────────────
// (replace the old lines 998–1024 with:)
for (group_id, combined_blob) in &combined_blobs {
    tx.execute(
        "UPDATE segment_group SET delta_blob = ?1 WHERE id = ?2",
        params![combined_blob, group_id],
    )?;
}
```

The transaction now contains only fast `INSERT` + `UPDATE` statements,
not ZSTD decompression+recompression.

---

## Implementation Order

| Step | Change | Expected speedup |
|---|---|---|
| 1 | `collect_singleton_kmers` parallel | ~Nx on reference path |
| 2 | `find_splitters_in_contigs` parallel | ~Nx on reference path |
| 3 | `build_fallback_index` + `collect_ref_kmers` parallel | ~Nx Phase 1 of append |
| 4 | Pre-build LzDiff cache | Largest win: eliminates redundant LZ index builds |
| 5 | Move ZSTD out of transaction | Shortens critical section; ~Nx on groups |

N = number of rayon threads (default = hardware threads).

---

## Additional Notes

### `LzDiff::encode` mutability

`encode(&mut self, ...)` takes `&mut self` only to allow lazy
`prepare_index()` on first call.  After `lz_from_ref_blob()` the index
is already built, so the `&mut` is never actually needed.
A follow-up refactor: move `index_ready` to `Cell<bool>`, change
`encode` to `&self`, and implement `Sync` — this enables zero-contention
sharing with `Arc<LzDiff>` (no Mutex needed at all).

### Re-compression cost on repeated appends

Each `add_as_delta` call re-decompresses + re-compresses the full
`delta_blob` for every touched group.  For N appends this is O(N²) in
compressed blob I/O.  A future optimization: switch to an
append-only framing (multiple independent ZSTD frames concatenated) so
old deltas are never re-read.  This is orthogonal to the changes above.

### No new dependencies needed

All changes use `rayon` (already in `Cargo.toml`) and `std::sync::Mutex`
(stdlib).  No `dashmap` or other crate is required.
