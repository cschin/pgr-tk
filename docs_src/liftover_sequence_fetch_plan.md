# Implementation Plan: Sequence Fetching for `pgr align liftover-gtf`

**Date:** 2026-04-04
**Author:** Jason Chin
**Branch:** `aln_map_improvment`

---

## Motivation

`pgr align liftover-gtf` currently maps GTF transcript annotations from GRCh38 to
haplotype contigs and writes coordinates (GTF, SQLite, anomaly reports).  It
does not retrieve the underlying nucleotide sequences.  Adding sequence output
enables downstream analyses such as:

- Comparing reference vs. haplotype transcript sequences to detect SNPs/indels
  within exons
- Validating liftover quality by aligning the two sequences
- Providing input to transcript-level expression tools that require
  haplotype-aware FASTA references

Both the **spliced reference transcript** sequence and the **spliced liftover
contig** sequence are needed.  Sequences are fetched from AGC (Assembled
Genomes Compressor) archives using the existing `pgr_db::agc_io::AGCFile` API
already used by other `pgr-bin` binaries.

---

## Scope

| In scope | Out of scope |
|---|---|
| Fetch spliced ref sequences from a GRCh38 AGC archive | Alignment of ref vs. contig sequences |
| Fetch spliced contig sequences from a haplotype AGC archive | CDS / protein translation |
| Write FASTA files for all hits and for HQ hits only | Variant calling |
| Store sequences in the output SQLite database | Multi-threaded fetch |
| Auto-detect sample name from single-sample archives | |

---

## AGC API

The project provides `pgr_db::agc_io::AGCFile` (gated by the `with_agc`
feature, which is enabled by default in `pgr-bin`).  The relevant calls are:

```rust
let agc = pgr_db::agc_io::AGCFile::new(path: String) -> Result<AGCFile, io::Error>;
let seq: Vec<u8> = agc.get_sub_seq(
    sample_name: String,   // e.g. "hg38" or "HG002"
    contig_name: String,   // e.g. "chr1" or "HG002#2#chr1"
    bgn: usize,            // 0-based inclusive
    end: usize,            // 0-based exclusive
);
```

- Returns uppercase ASCII bytes (A/C/G/T/N).  No built-in strand awareness.
- Reverse-complement is handled by `pgr_db::fasta_io::reverse_complement`.
- Sample name is auto-detected from `agc.samples[0].name` when the archive
  has exactly one sample; a `--*-agc-sample` override flag handles
  multi-sample archives.

No `Cargo.toml` changes are required — `pgr-bin` already declares
`default = ["with_agc"]`.

---

## New CLI Flags

Added to `CmdOptions` in `pgr align liftover-gtf.rs`:

| Flag | Type | Description |
|---|---|---|
| `--ref-agc <path>` | `Option<String>` | AGC archive containing GRCh38 reference sequences |
| `--ref-agc-sample <name>` | `Option<String>` | Sample name in ref archive (auto-detected if omitted and archive has 1 sample) |
| `--query-agc <path>` | `Option<String>` | AGC archive containing haplotype contig sequences |
| `--query-agc-sample <name>` | `Option<String>` | Sample name in query archive (auto-detected if omitted) |

Both `--ref-agc` and `--query-agc` are independently optional.  Omitting both
leaves current behavior completely unchanged.

---

## Output Files

| File | Written when | Content |
|---|---|---|
| `<prefix>.ref_tx.fa` | `--ref-agc` provided | Spliced reference transcript sequences, one record per transcript, strand-corrected |
| `<prefix>.contig_tx.fa` | `--query-agc` provided | Spliced contig sequences for every liftover hit ≥ `--min-coverage`, strand-corrected |
| `<prefix>.hq_contig_tx.fa` | `--query-agc` provided | Same, filtered to hits ≥ `--full-coverage` (default 90%) |

FASTA headers carry enough metadata to be self-describing:

```
# ref_tx.fa
>NM_000059.4 gene_id=BRCA2 chrom=chr13:32315507-32400268 strand=+

# contig_tx.fa / hq_contig_tx.fa
>NM_000059.4 contig=HG002#2#chr13:31261907-31345680 strand=+ coverage=71.3
```

---

## SQLite Schema Additions

Two new tables are added by `init_output_db`.  They are populated only when
the corresponding AGC flag is provided; empty tables add no overhead otherwise.

```sql
CREATE TABLE IF NOT EXISTS ref_sequences (
    transcript_pk  INTEGER PRIMARY KEY
                       REFERENCES transcripts(transcript_pk)
                       ON DELETE CASCADE,
    ref_seq        TEXT NOT NULL   -- spliced, strand-corrected ASCII
);

CREATE TABLE IF NOT EXISTS contig_sequences (
    liftover_pk    INTEGER PRIMARY KEY
                       REFERENCES liftover(liftover_pk)
                       ON DELETE CASCADE,
    contig_seq     TEXT NOT NULL   -- block-concatenated, strand-corrected ASCII
);
```

Sequences are stored in separate tables (not as columns on `liftover`) to
avoid bloating coordinate-only queries with multi-kilobase TEXT fields.

---

## New Helper Functions

All sequence-fetching functions are gated with `#[cfg(feature = "with_agc")]`.

### `open_agc_and_sample`

```rust
fn open_agc_and_sample(
    path: &str,
    sample_override: Option<&str>,
) -> Result<(AGCFile, String), io::Error>
```

Opens the archive and resolves the sample name.  Returns an error if
`sample_override` is `None` and the archive contains more than one sample.

### `fetch_ref_tx_seq`

```rust
fn fetch_ref_tx_seq(
    agc: &AGCFile,
    sample: &str,
    chrom: &str,
    exons: &[(u32, u32)],   // 0-based half-open, sorted ascending
    strand: char,
) -> Vec<u8>
```

Iterates exon intervals, fetches each chunk, concatenates, then
reverse-complements the whole splice if `strand == '-'`.

### `fetch_contig_tx_seq`

```rust
fn fetch_contig_tx_seq(
    agc: &AGCFile,
    sample: &str,
    hit: &LiftoverResult,
) -> Vec<u8>
```

Iterates `hit.block_sizes` / `hit.block_starts` to reconstruct the exon
blocks on the contig, concatenates, then reverse-complements if
`hit.strand == '-'`.

### `write_fasta_record`

```rust
fn write_fasta_record(
    out: &mut BufWriter<File>,
    header: &str,
    seq: &[u8],
) -> io::Result<()>
```

Writes a standard 80-column wrapped FASTA record.  No feature gate — uses
`std` only.

---

## Integration in `main()`

### Placement

The sequence pass runs as a **separate second sweep** after the liftover
transaction commits.  It iterates `&transcripts` again and re-calls
`liftover_transcript` (which is pure; `target_intervals` is still in memory).
This keeps the existing per-transcript loop clean.

### Control flow sketch

```
parse args
open ref_agc_opt, query_agc_opt  (if flags provided)
...existing liftover logic...
liftover tx.commit()

#[cfg(with_agc)]
if ref_agc_opt.is_some() || query_agc_opt.is_some() {
    open output FASTA files
    open DB sequence insert statements
    for t in &transcripts {
        if ref_agc_opt → fetch_ref_tx_seq → write FASTA + insert DB
        if query_agc_opt {
            re-run liftover_transcript to get hits
            for each hit → fetch_contig_tx_seq → write contig_tx.fa
                if coverage >= full_coverage → write hq_contig_tx.fa
                insert contig_sequences table (liftover_pk via DB lookup)
        }
    }
    seq_tx.commit()
    eprintln output paths
}
```

The `liftover_pk` for `contig_sequences` is looked up with:

```sql
SELECT liftover_pk FROM liftover
WHERE transcript_pk=? AND contig=? AND contig_start=? AND contig_end=?
```

covered by the existing `idx_liftover_contig` index.

---

## Implementation Steps

1. **`CmdOptions`** — add 4 new `Option<String>` fields with `#[clap(long)]`
2. **`init_output_db`** — add `CREATE TABLE IF NOT EXISTS` for `ref_sequences`
   and `contig_sequences`
3. **`write_fasta_record`** — add unconditionally (std only)
4. **`#[cfg(feature = "with_agc")]` helpers** — add `open_agc_and_sample`,
   `fetch_ref_tx_seq`, `fetch_contig_tx_seq`
5. **`use` declarations** — add `pgr_db::agc_io::AGCFile` and
   `pgr_db::fasta_io::reverse_complement` inside `#[cfg(feature = "with_agc")]`
6. **`main()`** — insert sequence pass after liftover commit (inside
   `#[cfg(feature = "with_agc")]` block)
7. **`eprintln!` summary** — extend to report new output file paths

---

## Example Usage

```bash
# Both reference and haplotype sequences
../target/release/pgr align liftover-gtf \
    hg002_hap0.alndb \
    hg38.ncbiRefSeq.gtf.gz \
    hg002_hap0_liftover.db \
    --target-chr-prefix "GRCh38#0#" \
    --min-coverage 0.5 \
    --ref-agc   hg38.agcrs \
    --query-agc hg002_hap0.agcrs

# Expected additional output:
#   hg002_hap0_liftover.ref_tx.fa        — 207 k spliced ref transcripts
#   hg002_hap0_liftover.contig_tx.fa     — 174 k spliced contig transcripts
#   hg002_hap0_liftover.hq_contig_tx.fa  — 79 k high-quality (≥90%) transcripts
```

---

## Open Questions

1. **Memory vs. re-compute**: Re-running `liftover_transcript` in the second
   pass is clean but does O(n_transcripts) extra interval queries.  For the
   current dataset (~207 k transcripts) this is fast.  If datasets grow
   significantly, caching `Vec<LiftoverResult>` in the first pass would be
   preferable.

2. **Large sequences in SQLite**: Storing full sequences in SQLite TEXT columns
   is convenient for queries but may produce large DB files for whole-genome
   datasets.  A `--no-store-seq` flag could suppress DB insertion while still
   writing FASTA files.

3. **Multi-sample AGC archives**: Auto-detection only works for single-sample
   archives.  The `--*-agc-sample` flags handle the multi-sample case but
   require the user to know the sample name.  A future `--list-agc-samples`
   subcommand could help discovery.
