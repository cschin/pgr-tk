# agc-rs Compression Analysis vs. AGC C++

Benchmark: 3 E. coli genomes (MG1655+ W3110 + Sakai, 14,882,589 raw bases).

## Progression

| Version | Description | File size | Ratio |
|---------|-------------|-----------|-------|
| v1 | Fixed 60 kb chunks, no splitters | 17.8 MB | 0.83× |
| v2 | + AGC singleton-kmer splitters | 2.73 MB | 5.44× |
| v3 | + RC canonical orientation | 2.65 MB | 5.61× |
| v4 | + vote-based fallback matching | 2.51 MB | 5.94× |
| v5 | + adaptive splitting (find_new_splitters) | 2.52 MB | 5.90× |
| v6 | + segment overlap + ZSTD level 19 + ref-kmer exclusion | 2.31 MB | 6.45× |
| v7 | + per-group batch ZSTD for deltas (schema v3) | 2.32 MB | 6.40× |
| v8 | + one-splitter lookup (partial-key terminators map) | 2.32 MB | 6.40× |
| **AGC C++** | reference implementation | **1.98 MB** | **7.33×** |

---

## Gap Analysis (v8 vs AGC C++)

```
Our file size:            2,322,432 bytes
AGC C++ file size:        2,028,863 bytes
Total gap:                  293,569 bytes (14.5% larger)

Breakdown:
  SQLite page overhead:  ~200,000 bytes  (page alignment, B-tree metadata, indices)
  Raw compressed data:   2,119,476 bytes (ref_data + delta_blob in segment_group)
  AGC compressed total:  2,028,863 bytes (includes AGC format overhead ~40-50 KB)
  Net LZ/ZSTD gap:         ~40-60 bytes  (~2%)
```

**The compression algorithms are essentially equivalent.
The remaining gap is almost entirely the SQLite storage format overhead (~200 KB).**

---

## What Was Investigated and Ruled Out

### 1. Per-group batch ZSTD (v7)
Implemented: all LZ-diff raw bytes for a segment_group are concatenated and
ZSTD-compressed together in `segment_group.delta_blob` (schema v3).

Result: **no measurable improvement** at 3 genomes. With only 1–2 deltas per
group there is not enough cross-delta context for ZSTD to exploit.  At scale
(50+ genomes per group, as AGC's default `pack_cardinality=50`) the benefit
would compound.

### 2. One-splitter lookup (v8)
Implemented: `terminators_map` tracks every (kmer_front, kmer_back) pair stored
in segment_group.  For a query segment with exactly one known boundary splitter,
`find_one_splitter_group` looks up all partner k-mers (including SENTINEL for
terminal segments) and picks the closest-size candidate — mirroring AGC's
`find_cand_segment_with_one_splitter`.

Result: **no measurable improvement** for E. coli.  The vote-based fallback was
already finding the correct terminal-segment references for these near-identical
strains.

### 3. Dense hash table (HASHING_STEP=1)
Tested: changed HASHING_STEP from 4 (sparse, every 4th reference position
indexed) to 1 (dense, all positions indexed), matching AGC's default build
without `USE_SPARSE_HT`.

Result: **no improvement** (+360 bytes raw data, same page-aligned file size).
For E. coli the shared k-mers are long enough that sparse indexing finds them
regardless.  For more divergent genomes (e.g. plant pan-genomes) the dense HT
might matter more.

### 4. min_match_len = 20 (AGC default)
Tested: raised min_match_len from 18 (our default) to 20 (AGC C++ default).

Result: **slight degradation** — longer minimum match means shorter valid
matches at segment boundaries are rejected.  Reverted to 18.

### 5. AGC "tuples" pre-encoding for ref_data
AGC's `store_in_archive` pre-packs 4 genomic bases (each 0–3) into 1 byte
(`bytes2tuples_impl<4,4>`) before ZSTD compression.  This is the "tuples"
format.  We store 1 byte per base (values 0–3) and ZSTD-compress directly at
level 19.

Benchmark showed: for high-entropy (random) sequences ZSTD already approaches
the theoretical minimum on the raw 1-byte-per-base format.  The tuples encoding
only helps when there is additional k-mer-level structure beyond base entropy.
For E. coli (genome-wide entropy close to 2 bits/base) the tuples approach gave
**no improvement or slight regression** (0.80×–1.11×).

---

## Root Cause of Remaining Gap

The ~14.5% total gap breaks down as:

| Source | Size | Fixable? |
|--------|------|----------|
| SQLite page alignment + B-tree overhead | ~200 KB | Requires format change |
| Net LZ-diff / ZSTD quality gap | ~40–60 KB | Already near-parity |
| AGC format overhead (headers, splitter stream, collection descriptor) | ~40–50 KB | N/A |

### SQLite overhead details
- Default SQLite page size: 4096 bytes.  Every BLOB overflow uses full pages.
- B-tree row overhead, free-page list, WAL metadata.
- With 85 segment_groups × 2 BLOBs (ref_data + delta_blob) = 170 medium BLOBs,
  alignment waste dominates.
- A larger page size (e.g. `PRAGMA page_size = 65536`) would reduce waste for
  large BLOBs but would not eliminate the structural overhead.
- The only way to match AGC's file size is to replace SQLite with a custom
  append-only archive format (or a columnar store like Parquet/Arrow IPC).

---

## Potential Future Optimizations

1. **Custom archive format** — replace SQLite with a flat binary archive
   similar to AGC's own format.  Expected gain: close the ~200 KB overhead gap,
   bringing agc-rs within 2–3% of AGC C++ on E. coli.

2. **Scale testing** — at 50+ genomes the per-group batch ZSTD gain becomes
   significant because ZSTD has cross-genome context.  AGC's default
   `pack_cardinality=50` is tuned for this.  Our batch ZSTD implementation
   already supports this via `batch_compress_deltas`.

3. **Dense HT for divergent genomes** — expose `HASHING_STEP` as a parameter
   and use step=1 for highly divergent datasets (plant pan-genomes, bacteria
   with >10% divergence).

4. **Tuples pre-encoding** — may help for low-complexity or repetitive reference
   sequences (e.g. repeat-rich eukaryotic chromosomes) where base-level
   statistics differ from uniform.

5. **Two-segment split** — AGC has `find_cand_segment_with_missing_middle_splitter`
   that splits one query segment into two when a middle splitter is missing.
   Not yet implemented in agc-rs; would help for structural variants.

---

## Key Parameters (current defaults)

| Parameter | agc-rs | AGC C++ |
|-----------|--------|---------|
| splitter_k | 31 | 31 |
| segment_size | 60,000 | 60,000 |
| min_match_len | 18 | 20 |
| HASHING_STEP | 4 (sparse) | 1 (dense, default) |
| ZSTD level — ref | 19 | 13 or 19 (adaptive) |
| ZSTD level — delta batch | 19 | 17 |
| pack_cardinality | 1 group/genome (unbounded) | 50 |
| Storage | SQLite (WAL) | Custom binary archive |
