# pgr-alnmap Output Files

The command-line tool `pgr-alnmap` aligns a set of assembled contigs (the *query*) against a
reference genome (the *target*) using SHIMMER (sparse hierarchical minimizer matching) anchoring
followed by Wavefront Alignment (WFA) or Smith-Waterman (SW) for base-level variant calling.

For each run it writes the following output files:

| File | Description |
|------|-------------|
| `<prefix>.alnmap` | Alignment chains with match, variant, and SV candidate records |
| `<prefix>.ctgmap.bed` | Contig-to-reference mapping in BED format |
| `<prefix>.ctgmap.json` | Contig-to-reference mapping in JSON (input to `pgr-generate-diploid-vcf`) |
| `<prefix>.target_len.json` | Reference sequence lengths in JSON |
| `<prefix>.query_len.json` | Query contig lengths in JSON |
| `<prefix>.svcnd.bed` | Structural variant candidate regions in BED format |
| `<prefix>.ctgsv.bed` | Contig-level SV summary in BED format |
| `<prefix>.svcnd.seqs` | FASTA sequences spanning SV candidates (omitted with `--skip-uncalled-sv-seq-file`) |

---

## The `.alnmap` File

The `.alnmap` file encodes the full alignment of each contig as a sequence of tab-separated
records grouped into *alignment chains*. Each chain covers one contiguous mapping of a query
contig to a target sequence.

All coordinates are **0-based, half-open** `[start, end)`.

### Common fields (all record types)

Every line starts with these nine tab-separated fields:

| Col | Field | Description |
|-----|-------|-------------|
| 1 | `aln_idx` | Alignment chain index (zero-padded integer) |
| 2 | `rec_type` | Record type (see below) |
| 3 | `target_name` | Reference sequence name |
| 4 | `target_start` | Start position in the reference (0-based) |
| 5 | `target_end` | End position in the reference (exclusive) |
| 6 | `query_name` | Query contig name |
| 7 | `query_start` | Start position in the contig (0-based) |
| 8 | `query_end` | End position in the contig (exclusive) |
| 9 | `orientation` | Alignment strand: `0` = forward, `1` = reverse complement |

### Record types

The second field identifies the record type. Types come in three groups based on
whether the target region is uniquely mapped, duplicated, or part of an overlapping chain.

| Type | Suffix | Meaning |
|------|--------|---------|
| `B` | — | Begin of an alignment chain |
| `E` | — | End of an alignment chain |
| `M` | — | Match block: query and target are identical over this region |
| `M_D` | `_D` | Match block, but another query contig also maps to the same target region (duplicated) |
| `M_O` | `_O` | Match block from an overlapping alignment chain |
| `V` | — | Variant block: base-level variants were called between query and target |
| `V_D` | `_D` | Variant block in a duplicated target region |
| `V_O` | `_O` | Variant block from an overlapping alignment chain |
| `S` | — | SV candidate: the aligner could not produce a clean base-level alignment |
| `S_D` | `_D` | SV candidate in a duplicated target region |
| `S_O` | `_O` | SV candidate from an overlapping alignment chain |

Records with the `_D` suffix indicate that two or more query contigs map to the same
target region — a signal of potential segmental duplication. Records with the `_O` suffix
arise when alignment chains overlap each other on the target.

---

### B (Begin) record — extra fields

```
aln_idx  B  target_name  ts  te  query_name  qs  qe  orientation  query_length  ctg_orientation  target_dup  target_ovlp  query_dup  query_ovlp
```

| Field | Description |
|-------|-------------|
| `query_length` | Total length of the query contig (bases) |
| `ctg_orientation` | Dominant strand of the contig across the full chain: `0` = forward, `1` = reverse |
| `target_dup` | `1` if the target region is duplicated, `0` otherwise |
| `target_ovlp` | `1` if the target region is part of an overlapping chain, `0` otherwise |
| `query_dup` | `1` if the query region is duplicated, `0` otherwise |
| `query_ovlp` | `1` if the query region is part of an overlapping chain, `0` otherwise |

---

### E (End) record — extra fields

```
aln_idx  E  target_name  ts  te  query_name  qs  qe  orientation  query_length  ctg_orientation
```

The last two fields are the same as in the B record. The `ts`/`te` and `qs`/`qe` coordinates
here correspond to the **last** alignment block of the chain.

---

### M / M_D / M_O (Match) record

No extra fields beyond the common nine. The query and reference sequences are identical
over `[target_start, target_end)`.

---

### S / S_D / S_O (SV Candidate) record — extra fields

```
aln_idx  S  target_name  ts  te  query_name  qs  qe  orientation  ctg_orientation  diff_type
```

| Field | Description |
|-------|-------------|
| `ctg_orientation` | Dominant strand of the contig: `0` = forward, `1` = reverse |
| `diff_type` | Reason the aligner could not call variants (see table below) |

**`diff_type` values:**

| Code | Name | Meaning |
|------|------|---------|
| `A` | `FailAln` | WFA/SW alignment returned no result |
| `E` | `FailEndMatch` | The first or last 16 bases of the two sequences did not match (anchor check failed) |
| `S` | `FailShortSeq` | One or both sequences are shorter than 16 bases |
| `L` | `FailLengthDiff` | The length difference between query and target exceeds the `max_sw_aln_size` threshold |
| `U` | `Unknown` | Unclassified failure |

SV candidate regions are written to `<prefix>.svcnd.bed` and, unless
`--skip-uncalled-sv-seq-file` is set, to `<prefix>.svcnd.seqs` as FASTA sequences for
downstream structural variant analysis.

---

### V / V_D / V_O (Variant) record — extra fields

```
aln_idx  V  target_name  ts  te  query_name  qs  qe  orientation  target_diff  query_diff  target_coord  variant_type  ref_seq  alt_seq
```

| Field | Description |
|-------|-------------|
| `target_diff` | Offset of the variant from `target_start` within this alignment block |
| `query_diff` | Offset of the variant from `query_start` within this alignment block |
| `target_coord` | Absolute position of the variant in the target sequence (`target_start + target_diff`) |
| `variant_type` | Single character: `X` = SNP/mismatch, `I` = insertion, `D` = deletion |
| `ref_seq` | Reference (target) sequence at the variant site |
| `alt_seq` | Query sequence at the variant site |

**`variant_type` characters:**

| Code | Meaning | `ref_seq` vs `alt_seq` |
|------|---------|------------------------|
| `X` | SNP or multi-nucleotide mismatch | equal length, different sequence |
| `I` | Insertion in the query | `ref_seq` is shorter than `alt_seq` |
| `D` | Deletion from the query | `ref_seq` is longer than `alt_seq` |

Each V record corresponds to a single variant event. A single alignment block may produce
multiple consecutive V records. These are merged into a single VCF record by
`pgr-generate-diploid-vcf` when the variants overlap on the reference.

---

## Practical Recipes

### Extract uniquely mapped blocks only

```bash
awk '$2 == "M" || $2 == "V" || $2 == "S"' sample.alnmap \
  | cut -f1-9 \
  | sort -k3,3 -k4,4n -u \
  > sample_unique_blocks.alnmap
```

### Extract duplicated-region blocks

```bash
awk '$2 == "M_D" || $2 == "V_D" || $2 == "S_D"' sample.alnmap \
  | cut -f1-9 \
  | sort -k3,3 -k4,4n -u \
  > sample_dup_blocks.alnmap
```

### Extract all called variants (for inspection or custom VCF building)

```bash
awk '$2 == "V"' sample.alnmap \
  | awk 'BEGIN{OFS="\t"} {print $3, $12, $12+length($14), $14, $15, $13}' \
  | sort -k1,1 -k2,2n \
  > sample_variants.tsv
# Columns: target_name, target_coord, target_coord_end, ref_seq, alt_seq, variant_type
```

### Count SV candidates per chromosome

```bash
awk '$2 == "S"' sample.alnmap | cut -f3 | sort | uniq -c | sort -rn
```

### Summarise alignment chain lengths

```bash
awk '$2 == "B" {print $3, $5-$4}' sample.alnmap \
  | awk '{sum[$1]+=$2; n[$1]++} END{for(c in sum) print c, n[c], sum[c]}' \
  | sort -k1,1V
# Columns: target_name, number_of_chains, total_bases_covered
```

---

## Relationship to Other Output Files

- **`.ctgmap.json`** — used as the `target_len_json_path` argument to
  `pgr-generate-diploid-vcf`; contains contig-to-reference mappings and target sequence
  lengths needed to write the VCF header.
- **`.svcnd.bed`** — BED file of target regions that are SV candidates; can be loaded
  directly into a genome browser alongside the alignment.
- **`.svcnd.seqs`** — FASTA sequences spanning SV candidate regions; input to
  `pgr-generate-sv-analysis` for principal-bundle-based structural variant decomposition.

---

## Design Limitations and Potential Future Formats

### Current design is clunky

The `.alnmap` format has several practical pain points that accumulate at whole-genome scale:

- **Variable column count per record type.** B records have 15 columns, E records have 11,
  M records have 9, S records have 11, and V records have 15. Every downstream parser must
  dispatch on column 2 before it knows how many fields to expect. This makes the format
  fragile — a single extra tab anywhere silently shifts all subsequent field indices.

- **Mixed concerns in one file.** Chain boundaries (B/E), match blocks (M), variant calls
  (V), and SV candidates (S) are interleaved in output order rather than grouped by type.
  Extracting just variants requires scanning the entire file even if the caller only needs
  V records.

- **No random access.** A human genome alignment against GRCh38 can produce hundreds of
  millions of lines. There is no index, so any query — "give me all variants on chr7
  between 100 Mb and 110 Mb" — requires a full linear scan and `grep`/`awk` filtering.

- **String-encoded numbers.** Every coordinate, length, and flag is serialised as ASCII
  text. At scale this wastes I/O bandwidth and CPU time compared with binary-encoded
  integers, and makes column statistics (mean block length, total aligned bases, etc.)
  expensive to compute.

- **Sequence data mixed with positional data.** The `ref_seq` and `alt_seq` strings in V
  records can be arbitrarily long (large indels, MNPs). Storing them inline in the same
  flat file as the coordinate data prevents efficient columnar access to either.

---

### Alternative 1 — Apache Parquet with separate tables

[Apache Parquet](https://parquet.apache.org) is a columnar binary format with built-in
compression and predicate pushdown. The `.alnmap` content maps naturally onto three
separate Parquet tables that can be stored alongside each other or inside a single
[Delta Lake](https://delta.io) or [Apache Iceberg](https://iceberg.apache.org) table set:

**`chains` table** — one row per alignment chain (B/E pair)

| Column | Type | Notes |
|--------|------|-------|
| `aln_idx` | `uint32` | chain ID |
| `target_name` | `string` (dictionary-encoded) | chromosome/contig name |
| `target_start` | `uint32` | |
| `target_end` | `uint32` | |
| `query_name` | `string` (dictionary-encoded) | |
| `query_start` | `uint32` | |
| `query_end` | `uint32` | |
| `orientation` | `uint8` | 0/1 |
| `ctg_orientation` | `uint8` | 0/1 |
| `query_length` | `uint32` | |
| `target_dup` | `bool` | |
| `target_ovlp` | `bool` | |
| `query_dup` | `bool` | |
| `query_ovlp` | `bool` | |

**`blocks` table** — one row per M/V/S record

| Column | Type | Notes |
|--------|------|-------|
| `aln_idx` | `uint32` | foreign key into `chains` |
| `block_type` | `uint8` | 0=M, 1=V, 2=S |
| `dup_flag` | `bool` | |
| `ovlp_flag` | `bool` | |
| `target_name` | `string` (dictionary-encoded) | |
| `target_start` | `uint32` | |
| `target_end` | `uint32` | |
| `query_name` | `string` (dictionary-encoded) | |
| `query_start` | `uint32` | |
| `query_end` | `uint32` | |
| `orientation` | `uint8` | |
| `sv_diff_type` | `uint8` (nullable) | only for S records |

**`variants` table** — one row per V record variant event

| Column | Type | Notes |
|--------|------|-------|
| `aln_idx` | `uint32` | foreign key into `chains` |
| `target_name` | `string` (dictionary-encoded) | |
| `target_coord` | `uint32` | absolute position |
| `variant_type` | `uint8` | 0=X, 1=I, 2=D |
| `ref_seq` | `large_binary` | reference bases |
| `alt_seq` | `large_binary` | query bases |

**Advantages:**

- Any Parquet-aware tool (DuckDB, Polars, Pandas, Spark) can query variants on a specific
  chromosome range without reading the full file, using Parquet row-group statistics and
  predicate pushdown.
- Dictionary encoding on `target_name` and `query_name` compresses chromosome strings to
  a single byte per record.
- The `variants` table can be sorted by `(target_name, target_coord)` independently of
  the `blocks` table, enabling fast positional lookups.
- Parquet files are self-describing (schema embedded in the footer) — no separate
  documentation is needed to parse them.

**Quick conversion sketch using Python/Polars:**

```python
import polars as pl

# Read the current alnmap format
df = pl.read_csv("sample.alnmap", separator="\t", has_header=False,
                 infer_schema_length=0)  # all strings initially

variants = (
    df.filter(pl.col("column_2").str.starts_with("V"))
    .select([
        pl.col("column_1").cast(pl.UInt32).alias("aln_idx"),
        pl.col("column_3").alias("target_name"),
        pl.col("column_12").cast(pl.UInt32).alias("target_coord"),
        pl.col("column_13").alias("variant_type"),
        pl.col("column_14").alias("ref_seq"),
        pl.col("column_15").alias("alt_seq"),
    ])
)

variants.write_parquet("sample_variants.parquet",
                       compression="zstd",
                       statistics=True)
```

---

### Alternative 2 — Apache Arrow IPC with an index

[Apache Arrow IPC](https://arrow.apache.org/docs/format/IPC.html) (the "Feather v2"
on-disk format) stores data in the same columnar layout that Arrow uses in memory. The
critical addition for random-access use is a lightweight companion index.

**Layout:**

```
sample.arrow          # Arrow IPC file, record batches sorted by (target_name, target_coord)
sample.arrow.idx      # Small index: one entry per record batch → (target_name, min_coord, max_coord, byte_offset)
```

Each record batch covers one chromosome (or one fixed-size coordinate window). The index
fits entirely in memory (a few kilobytes for a human genome), so any positional query
resolves to a single `pread()` call into the Arrow file.

**Schema (single unified table):**

```
record_type:    uint8         # 0=chain_bgn, 1=chain_end, 2=match, 3=variant, 4=sv_cnd
aln_idx:        uint32
dup_flag:       bool
ovlp_flag:      bool
target_name:    dictionary<uint16, utf8>
target_start:   uint32
target_end:     uint32
query_name:     dictionary<uint16, utf8>
query_start:    uint32
query_end:      uint32
orientation:    uint8
# nullable fields — set only for relevant record types
query_length:   uint32        (chain_bgn / chain_end only)
ctg_orientation:uint8         (chain_bgn / chain_end only)
target_coord:   uint32        (variant only)
variant_type:   uint8         (variant only — 0=X, 1=I, 2=D)
ref_seq:        large_utf8    (variant only)
alt_seq:        large_utf8    (variant only)
sv_diff_type:   uint8         (sv_cnd only — 0=FailAln … 4=Unknown)
```

**Advantages over Parquet for this use case:**

- Arrow IPC is a streaming format — a writer can append record batches incrementally as
  `pgr-alnmap` processes each contig, without buffering the entire output in memory first.
  Parquet requires knowing the full row group before it can flush.
- The in-process representation is identical to the on-disk layout, so there is zero
  deserialisation cost when loading into Rust (`arrow2`/`arrow-rs`), Python (`pyarrow`),
  or R (`arrow`).
- The companion index makes random-access region queries `O(log n)` in the number of
  record batches rather than `O(n)` in the number of lines.

**Reading a region in Python:**

```python
import pyarrow as pa
import pyarrow.ipc as ipc
import struct, mmap

# Load the lightweight index
with open("sample.arrow.idx", "rb") as f:
    index = [struct.unpack("QII", f.read(16)) for _ in range(f.read())]

# Find record batches overlapping chr7:100_000_000-110_000_000
target_offsets = [off for off, mn, mx in index if mn < 110_000_000 and mx > 100_000_000]

with open("sample.arrow", "rb") as f:
    mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    for offset in target_offsets:
        reader = ipc.open_stream(pa.BufferReader(mm[offset:]))
        batch = reader.read_next_batch()
        # filter locally to exact range
        mask = (pa.compute.greater(batch["target_end"],   100_000_000) &
                pa.compute.less   (batch["target_start"], 110_000_000))
        print(batch.filter(mask).to_pandas())
```

---

### Summary comparison

| Criterion | Current `.alnmap` | Parquet (3 tables) | Arrow IPC + index |
|-----------|:-----------------:|:------------------:|:-----------------:|
| Human-readable | Yes | No | No |
| Random access by position | No (full scan) | Partial (row-group stats) | Yes (byte-level index) |
| Streaming write | Yes | No (row-group buffering) | Yes |
| Columnar compression | No | Yes (Zstd/Snappy) | Yes (LZ4/Zstd) |
| Zero-copy in-process read | No | No | Yes |
| Schema self-describing | No | Yes | Yes |
| Tooling (grep/awk) | Yes | No | No |
| Inter-op (Python/R/SQL) | Limited | Excellent | Excellent |

The Parquet approach is the better fit if the primary use case is analytical (joining
variants against annotation databases, computing statistics across chromosomes). The Arrow
IPC approach is better if the primary use case is streaming processing or if `pgr-alnmap`
itself is to be refactored to write binary output incrementally without buffering.
