# pgr-mdb Memory Reduction: Partition-and-Merge Plan

## Problem

`pgr-mdb` consumes >300 GB of RAM when building a minimizer database from an
AGC archive that contains 100+ haplotypes of a human-scale genome.  Two
independent memory peaks drive this:

### Peak 1 — `par_fetch_seqs()` decompresses everything at once

`AGCFile::par_fetch_seqs` (`agc_io.rs`) opens one SQLite read-only connection
per rayon task and decompresses **all** contigs simultaneously:

```
100 haplotypes × 3 Gbp × 1 byte/bp ≈ 300 GB
```

The resulting `Vec<SeqRec>` is passed directly to `load_index_from_seq_vec`,
which also holds the data while building shimmer pairs.

### Peak 2 — `frag_map` accumulates before any flush

`ShmmrToFrags = FxHashMap<(u64,u64), Vec<FragmentSignature>>` grows without
bound until every sequence has been indexed:

```
~940 M shimmer pairs × 17 bytes/FragmentSignature
+ FxHashMap bucket overhead
≈ 16–20 GB (plus Vec heap allocations add substantially)
```

Both peaks overlap in time, giving the observed >300 GB total.

---

## Solution: Partition-and-Merge

Process `B` haplotypes at a time.  Each batch writes a **shard** (sorted
`.shard_N.mdbi` + `.shard_N.mdbv`).  After all batches, a k-way external merge
combines the sorted shards into the final `.mdbi` / `.mdbv`.  The `.midx`
SQLite seq-index is written once at the end from the accumulated (tiny)
`CompactSeq` metadata.

### Phase 1: Batched Indexing

```
for batch in sample_ctg_list.chunks(B):
    recs  = par_fetch_seqs_batch(batch)   # decompress only B haplotypes
    seqs  = assign sids continuing from previous batch
    self.load_index_from_seq_vec(&seqs)   # build frag_map for this batch
    write_shmmr_map_shard(prefix, shard_id)
    self.frag_map.clear(); self.frag_map.shrink_to_fit()
```

Each shard's `.mdbi` has keys sorted in ascending `(k1, k2)` order (same
invariant as `write_shmmr_map_split`).

### Phase 2: K-Way External Merge

```
open streaming key readers for each shard .mdbi
open seekable file handles for each shard .mdbv
initialize min-heap with first key from each shard

while heap is non-empty:
    pop minimum (k1, k2) from shard S
    collect fragments from S (and drain any other shards with same key)
    write concatenated fragment bytes to final .mdbv
    emit idx entry → final .mdbi
    advance each contributing shard in heap

delete shard files
```

Memory during merge ≈ O(num_shards) key entries + I/O buffers — negligible.

### Seq Metadata

`CompactSeq` (name, source, len, seq_frag_range) is tiny (~100 bytes).
It is accumulated in `self.seqs` across all batches and written to `.midx`
(SQLite) once after the merge.

---

## Memory Profile

| Mode | Decompress peak | frag_map peak | Total peak |
|---|---|---|---|
| Current (100 haplotypes) | ~300 GB | ~16–20 GB | **>300 GB** |
| B=10 (10 shards) | ~30 GB | ~2 GB | **~32 GB** |
| B=5  (20 shards) | ~15 GB | ~1 GB | **~16 GB** |

---

## Implementation

### New API surface

**`agc_io.rs`**

```rust
impl AGCFile {
    /// Expose the full (sample, contig) list for batch slicing.
    pub fn sample_ctg_list(&self) -> &[(String, String)]

    /// Decompress only the given (sample, contig) pairs in parallel.
    pub fn par_fetch_seqs_batch(&self, batch: &[(String, String)]) -> Vec<SeqRec>
}
```

**`seq_db.rs`**

```rust
/// Write one shard: same binary layout as write_shmmr_map_split.
pub fn write_shmmr_map_shard(
    shmmr_spec: &ShmmrSpec,
    shmmr_map:  &ShmmrToFrags,
    prefix:     &str,
    shard_id:   usize,
) -> Result<(), io::Error>

/// K-way external merge of sorted shards → final .mdbi / .mdbv.
pub fn merge_shmmr_map_shards(
    shmmr_spec: &ShmmrSpec,
    prefix:     &str,
    num_shards: usize,
) -> Result<(), io::Error>

impl CompactSeqDB {
    /// Load index from AGC file in memory-bounded batches.
    pub fn load_index_from_agcfile_batched(
        &mut self,
        agcfile:    &AGCFile,
        batch_size: usize,
        prefix:     &str,
        agc_path:   Option<&str>,
    ) -> Result<(), io::Error>
}
```

**`pgr-mdb.rs`**

```
--batch-size <N>    Haplotypes per processing batch [default: 0 = all at once]
```

When `--batch-size > 0`, call `load_index_from_agcfile_batched` and skip the
separate `write_shmmr_map_index` call (batched function handles everything).

---

## Implementation Order

1. `agc_io.rs` — `sample_ctg_list` + `par_fetch_seqs_batch` (trivial)
2. `seq_db.rs` — `write_shmmr_map_shard` (copy of split writer, different paths)
3. `seq_db.rs` — `merge_shmmr_map_shards` (streaming k-way merge)
4. `seq_db.rs` — `load_index_from_agcfile_batched` (orchestrator)
5. `pgr-mdb.rs` — `--batch-size` flag + call site

---

## Shard File Format

Shards use the identical binary layout as the final `.mdbi`/`.mdbv`:

- `{prefix}.shard_{N}.mdbi` — sorted key-index (28 bytes/key)
- `{prefix}.shard_{N}.mdbv` — fragment values (17 bytes/record)

Headers embed the ShmmrSpec so each shard is self-describing.  Shard files are
deleted after a successful merge.

---

## Correctness Notes

- `sid` values are assigned globally across batches (counter continues across
  batches), so `seq_id` in `FragmentSignature` records remain unique.
- `frag_id` within each `FragmentSignature` is sequence-local (0-based index
  within that sequence's shimmer pairs).  The AGC-backed lookup path only uses
  `(seq_id, bgn, end, orientation)`, so `frag_id` pass-through is safe.
- Each shard's keys are already sorted by `write_shmmr_map_shard`, so the
  k-way merge is a standard external merge with no re-sorting step.
- The merge preserves the `(k1 ≤ k2)` invariant established by `seq_to_index`.
