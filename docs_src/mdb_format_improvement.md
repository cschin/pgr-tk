# pgr-tk `.mdb` / `.midx` Format Improvement Proposal

## Background

The shimmer database used by `pgr-mdb` and `pgr-query` consists of four files written
by `CompactSeqDB::write_shmmr_map_index` and read back by `SeqIndexDB::load_from_agc_index`
/ `load_from_frg_index`.

| File | Encoding | Compression | Purpose |
|------|----------|-------------|---------|
| `.mdb`  | Custom `byteorder` LE binary | None | Shimmer pair → fragment locations |
| `.midx` | Tab-delimited text → **SQLite** (implemented) | None | Contig metadata (sid, len, name, source) |
| `.sdx`  | `bincode` v2 | None | Chunk offsets + `Vec<CompactSeq>` |
| `.frg`  | `bincode` v2 | Deflate per 256-chunk | Raw fragment sequences |

This document covers improvements to the `.mdb` (proposed) and `.midx` (already
implemented) files.

---

## `.midx` — Implemented: SQLite Database

### Previous format (tab-delimited text)

```
{sid}\t{len}\t{ctg_name}\t{source}\n
```

No types, no schema version, no indices, not queryable from external tools.

### New format (SQLite, `PRAGMA user_version = 1`)

Written by `write_seq_index_sqlite` / read by `read_seq_index_sqlite` in
`seq_db.rs`.  The file keeps the `.midx` extension so no call-site renames
are needed.

```sql
PRAGMA user_version = 1;

CREATE TABLE shmmr_spec (
    w        INTEGER NOT NULL,
    k        INTEGER NOT NULL,
    r        INTEGER NOT NULL,
    min_span INTEGER NOT NULL,
    sketch   INTEGER NOT NULL
);

CREATE TABLE seq_index (
    sid      INTEGER PRIMARY KEY,
    len      INTEGER NOT NULL,
    ctg_name TEXT    NOT NULL,
    source   TEXT
);

CREATE INDEX idx_ctg_source ON seq_index (ctg_name, source);
```

Benefits over the old text file:

| Property | Old tab-delimited | New SQLite |
|----------|-------------------|------------|
| Types | All strings, parsed at runtime | INTEGER / TEXT, enforced by schema |
| Schema version | None | `PRAGMA user_version = 1` |
| External tooling | `awk`, `cut` only | `sqlite3`, DuckDB, Python `sqlite3`, R |
| Name lookup | Full file scan | O(log n) via `idx_ctg_source` index |
| `shmmr_spec` stored | No (only in `.mdb`) | Yes (deduplicated source of truth) |
| Corruption detection | Silent | SQLite integrity checks |

The `shmmr_spec` table provides a second copy of the compression parameters
alongside the `.mdb` header, which will become the primary copy once the
`.mdb` → `.mdbi`+`.mdbv` split below is implemented.

---

## Current `.mdb` Binary Layout

Written by `write_shmmr_map_file` (`seq_db.rs:1320`):

```
[3 bytes]  magic          ASCII "mdb"  — no version field
[4 bytes]  w              u32 LE
[4 bytes]  k              u32 LE
[4 bytes]  r              u32 LE
[4 bytes]  min_span       u32 LE
[4 bytes]  sketch         u32 LE  (boolean packed as full u32)
[8 bytes]  n_pairs        u64 LE  — number of unique shimmer pairs

for each shimmer pair (order = FxHashMap iteration, non-deterministic):
    [8 bytes]  k1          u64 LE
    [8 bytes]  k2          u64 LE
    [8 bytes]  vec_len     u64 LE
    for each fragment (17 bytes, hard-coded stride):
        [4 bytes]  frag_id    u32 LE
        [4 bytes]  seq_id     u32 LE
        [4 bytes]  bgn        u32 LE
        [4 bytes]  end        u32 LE
        [1 byte]   orient     u8
```

---

## How pgr-query Uses the Format

Understanding the actual access pattern is essential for evaluating any
replacement.

### Startup (`load_from_agc_index`, `ext.rs:87`)

1. **Linear scan** — `read_mdb_file_to_frag_locations` (`seq_db.rs:1438`) reads
   the entire `.mdb` sequentially, key by key, recording each pair's byte offset.
   Returns `ShmmrIndexFileLocation: Vec<((u64,u64),(usize,usize))>`.

2. **Build FxHashMap** — the offset vector is converted to
   `frag_location_map: FxHashMap<(u64,u64),(usize,usize)>`.
   For a human genome (~80 M unique shimmer pairs): **~3 GB of RAM** just for
   the key-to-offset map.

3. **Memory-map the file** — `Mmap::map(&fmap_file)` maps the full `.mdb` so
   fragment records can be read by byte offset without an extra `read()` call.

### Hot query path (`raw_query_fragment_from_mmap_midx`, `seq_db.rs:1259`)

For every shimmer pair `(s0, s1)` extracted from the query sequence (parallelised
with `par_iter`):

```
frag_location_map.get(&(s0, s1))         →  (byte_offset, vec_len)
frag_map_mmap_file[byte_offset ..]       →  decode 17-byte fragment records
```

Access pattern: **O(1) FxHashMap lookup + random mmap read** per shimmer pair.
A typical query produces hundreds to thousands of shimmer pairs.

---

## Problems with the Current Format

### 1. Startup requires scanning the full file
`read_mdb_file_to_frag_locations` must read all ~6 GB (human genome) sequentially
to learn where each key's data starts.  The file interleaves keys and values, so
there is no way to jump directly to the key table.  This scan dominates load time.

### 2. Key-offset map consumes ~3 GB of RAM
`frag_location_map` holds one `(u64, u64) → (usize, usize)` entry per unique
shimmer pair.  For 80 M pairs at ~40 bytes per FxHashMap slot: ~3.2 GB of heap
allocation that lives for the entire lifetime of `pgr-query`.

### 3. No schema version
The 3-byte magic `"mdb"` carries no version number.  An incompatible format
change silently corrupts reads — there is no way to detect old vs new files.

### 4. Non-deterministic output
`shmmr_map` is iterated as an `FxHashMap`, so two identical genomes produce
bit-for-bit different `.mdb` files.  This prevents content-based checksums and
reproducible builds.

### 5. Hard-coded 17-byte stride
`read_mdb_file_to_frag_locations` line 1494:
```rust
let advance = 17 * vec_len;
```
The record size is a magic constant in the seek logic.  Any field addition
silently breaks existing readers.

### 6. No compression
Fragment coordinate arrays are highly compressible (delta-coded coordinates,
mostly small values), but the file is stored as raw binary.  For a human genome
the data section alone is ~4 GB uncompressed.

### 7. u32 overflow risk
`frag_id`, `seq_id`, `bgn`, `end` are all `u32`.  A large pan-genome or a
single chromosome at single-base resolution approaches the 4 billion limit.

---

## Why Some Alternatives Are Unsuitable

**SQLite** — ruled out for the hot query path.  A B-tree lookup is 50–100×
slower than `FxHashMap::get` (~10–50 µs vs ~100 ns per key).  A query with
1,000 shimmer pairs would add 10–50 ms overhead per query vs 0.1 ms today.
SQLite is already used for agc-rs sequence storage, where query counts are low;
it is not appropriate here.

**Apache Arrow IPC (columnar)** — good for analytical access (scan all pairs
with a given `seq_id`, range queries, etc.) but the hot path is random
point-lookup by `(k1, k2)`, not a columnar scan.  Arrow's row-group structure
does not help random mmap reads.  Worth considering for an offline analysis
layer, but not as the primary query format.

---

## Proposed Format

The core idea is to **decouple the sorted key index from the fragment data**.
Interleaving them forces a full linear scan at startup; separating them allows
the key index to be read as a single sequential block, and enables binary search
on the mmapped key table as a future zero-RAM option.

### File 1: `{prefix}.mdbi` — sorted key index

Replaces the startup linear scan.

```
Offset  Size    Field
------  ------  -----
0       4       magic+version   b"mdbi"  (last byte = 0x01 for version 1)
4       4       w               u32 LE
8       4       k               u32 LE
12      4       r               u32 LE
16      4       min_span        u32 LE
20      4       sketch          u32 LE
24      8       n_keys          u64 LE

28 + i*24  (for i in 0..n_keys, SORTED by (k1, k2)):
    k1:           u64 LE    shimmer hash 1
    k2:           u64 LE    shimmer hash 2
    data_offset:  u64 LE    byte offset into .mdbv file
    vec_len:      u32 LE    number of fragment records at that offset
```

For 80 M shimmer pairs: **80 M × 24 = 1.92 GB** (vs scanning 6 GB today).

The sorted layout enables two modes:

- **HashMap mode (default)**: load `.mdbi` sequentially into `FxHashMap` in one
  pass — startup I/O drops from ~6 GB → ~1.9 GB, roughly 3× faster.
- **Binary-search mode (low-RAM option)**: `mmap` the `.mdbi` file and binary-
  search the key table directly (O(log 80 M) = ~27 comparisons, ~3–5 µs per
  lookup).  Eliminates the 3 GB FxHashMap from RAM entirely at the cost of
  ~30–50× slower per-key lookup — acceptable when memory is the bottleneck.

### File 2: `{prefix}.mdbv` — fragment values

Contains only fragment records, in the same order as the key index.

```
Offset  Size    Field
------  ------  -----
0       5       magic+version   b"mdbv\x01"
5       8       n_records       u64 LE   (total fragment records across all keys)

9 + i*17  (for i in 0..n_records):
    frag_id:   u32 LE
    seq_id:    u32 LE
    bgn:       u32 LE
    end:       u32 LE
    orient:    u8
```

The file is memory-mapped exactly as `.mdb` is today.  The hot-path read code
(`get_fragment_signatures_from_mmap_file`, `seq_db.rs:1502`) requires no change
beyond pointing at the `.mdbv` mmap instead of `.mdb`.

### File 3: `{prefix}.midx` — keep as-is

Tab-delimited text (~5 KB for a human genome).  Not a bottleneck.

---

## Optional Future Step: Block Compression

Once the key/value split is in place, Zstd block compression can be added to
`.mdbv` without changing the query path:

- Divide fragment records into fixed-size blocks (e.g. 4 096 records each).
- Zstd-compress each block independently.
- Store block byte offsets in `.mdbi` by changing `data_offset` to encode
  `(block_id: u32, within_block_offset: u16)`.
- On lookup: decompress the relevant block (cached in an LRU), read the records.

Expected compression ratio for fragment coordinate arrays: 3–4× (delta-coded
`bgn`/`end` with small increments), reducing `.mdbv` from ~4 GB to ~1–1.5 GB.

---

## Implementation Plan

All changes are confined to `pgr-db/src/seq_db.rs` and
`pgr-db/src/frag_file_io.rs`.  No changes to `pgr-query` or `pgr-mdb` binaries
are needed beyond updating the prefix constants.

### Step 1 — Write side (`seq_db.rs`)

Replace `write_shmmr_map_file`:

```rust
pub fn write_shmmr_map_index_v2(
    shmmr_spec: &ShmmrSpec,
    shmmr_map: &ShmmrToFrags,
    prefix: &str,
) -> Result<(), io::Error> {
    // 1. Collect and sort keys
    let mut keys: Vec<(u64, u64)> = shmmr_map.keys().copied().collect();
    keys.sort_unstable();

    // 2. Write .mdbv (values only, in key-sorted order)
    let mut val_file = BufWriter::new(File::create(format!("{prefix}.mdbv"))?);
    val_file.write_all(b"mdbv\x01")?;
    let total_recs: u64 = shmmr_map.values().map(|v| v.len() as u64).sum();
    val_file.write_u64::<LittleEndian>(total_recs)?;
    let mut offset: u64 = 13; // header size

    // 3. Write .mdbi (sorted key index)
    let mut idx_file = BufWriter::new(File::create(format!("{prefix}.mdbi"))?);
    idx_file.write_all(b"mdbi\x01")?;
    // ... write shmmr_spec fields ...
    idx_file.write_u64::<LittleEndian>(keys.len() as u64)?;

    for &(k1, k2) in &keys {
        let frags = &shmmr_map[&(k1, k2)];
        idx_file.write_u64::<LittleEndian>(k1)?;
        idx_file.write_u64::<LittleEndian>(k2)?;
        idx_file.write_u64::<LittleEndian>(offset)?;
        idx_file.write_u32::<LittleEndian>(frags.len() as u32)?;
        for f in frags {
            val_file.write_u32::<LittleEndian>(f.0)?;
            val_file.write_u32::<LittleEndian>(f.1)?;
            val_file.write_u32::<LittleEndian>(f.2)?;
            val_file.write_u32::<LittleEndian>(f.3)?;
            val_file.write_u8(f.4)?;
        }
        offset += frags.len() as u64 * 17;
    }
    Ok(())
}
```

### Step 2 — Read side (`seq_db.rs` / `frag_file_io.rs`)

Replace `read_mdb_file_to_frag_locations`:

```rust
pub fn read_mdbi_file(
    filepath: &str,
) -> Result<(ShmmrSpec, ShmmrToIndexFileLocation), io::Error> {
    let mut f = BufReader::new(File::open(filepath)?);
    // validate magic + version, read shmmr_spec, n_keys
    // sequential read of key table into FxHashMap — one pass, no seek
    let mut map = ShmmrToIndexFileLocation::default();
    for _ in 0..n_keys {
        let k1 = f.read_u64::<LittleEndian>()?;
        let k2 = f.read_u64::<LittleEndian>()?;
        let offset = f.read_u64::<LittleEndian>()? as usize;
        let vec_len = f.read_u32::<LittleEndian>()? as usize;
        map.insert((k1, k2), (offset, vec_len));
    }
    Ok((shmmr_spec, map))
}
```

### Step 3 — Update callers

- `load_from_agc_index` (`ext.rs:87`): call `read_mdbi_file(prefix + ".mdbi")`
  and `Mmap` the `.mdbv` file.
- `CompactSeqFragFileStorage::new` (`frag_file_io.rs:29`): same.
- `write_shmmr_map_index` (`seq_db.rs:821`): call `write_shmmr_map_index_v2`.

### Step 4 — Backward compatibility shim (optional)

Keep `read_mdb_file_to_frag_locations` for reading legacy `.mdb` files.
Auto-detect by checking magic bytes: `b"mdb"` → old path, `b"mdbi"` → new path.

---

## Summary

| Property | Current `.mdb` | Proposed `.mdbi`+`.mdbv` |
|----------|---------------|--------------------------|
| Startup I/O to build offset map | ~6 GB linear scan | ~1.9 GB sequential read |
| RAM for key-offset map | ~3 GB FxHashMap | ~3 GB (HashMap mode) or 0 (binary-search mode) |
| Hot query path | O(1) HashMap + mmap | O(1) HashMap + mmap (unchanged) |
| Disk usage (data) | ~6 GB | ~6 GB (same; Zstd optional later) |
| Schema version | 3-char magic, no version | 4-byte magic with version byte |
| Output determinism | Non-deterministic | Deterministic (sorted keys) |
| Hard-coded stride | Yes (`17 * vec_len`) | Yes (same; documented) |
| Implementation scope | — | `seq_db.rs` + `frag_file_io.rs` only |
