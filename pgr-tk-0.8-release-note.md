# PGR-TK 0.8.0 Release Notes

**Branch:** `aln_map_improvment` → `main`  
**Version:** 0.6.0 → 0.8.0  
**Date:** 2026-04-11  
**Commits since 0.6.0:** 104

---

## Overview

PGR-TK 0.8.0 is a major release that rearchitects the toolkit from the ground
up. The C++ AGC submodule is gone, replaced by a pure-Rust `agc-rs` crate with
a new SQLite-backed storage format. Twenty-two standalone binaries are
consolidated into a single `pgr` executable with named subcommand groups. Index
files move to a transparent, queryable format. Alignment output gains a
structured SQLite database. GTF liftover and SV annotation are now first-class
pipeline steps. End-to-end example pipelines—from a 5-minute E. coli demo to a
full 100+ haplotype HPRC R2 build—ship with the repository. The web server and
frontend are removed; a redesigned interface will follow in a future release.

---

## 1. agc-rs: Pure-Rust Genome Compression (replaces C++ AGC)

The C++ AGC submodule is fully removed. `agc-rs`, a faithful Rust port, is now
the exclusive compression backend and has graduated from prototype to production.

### What changed

- **Comparable compression performance:** `agc-rs` implements the same
  delta-encoding algorithm as AGC C++ and has been verified to produce
  equivalent compression on E. coli benchmarks. Formal benchmarks against
  AGC C++ on human haplotype assemblies are in progress.
- **SQLite storage:** Instead of an opaque binary blob, sequences are stored
  in a SQLite database (`.agcrs` file). This means the archive is inspectable
  with standard tools (`sqlite3`), importable from Python/R, and
  straightforward to back up, query, or diff.
- **New subcommands:**
  - `agc-rs create` / `append` — single-sample archive operations (unchanged API)
  - `agc-rs batch-append` — compress multiple samples in one invocation, sharing
    the delta-index build across all contigs in the batch (significant I/O
    savings for large batches)
  - `agc-rs merge` — combine pre-built sub-archives that share the same reference
    and splitter set into a single archive, enabling a parallelisable
    batch-then-merge strategy for very large collections
- **Parallel reads:** Each rayon task holds its own read-only SQLite connection
  (WAL mode), eliminating the mutex serialisation bottleneck that affected
  multi-threaded decompression in earlier versions.
- **`with_agc` feature flag removed:** AGC is no longer optional; all pgr-db
  and pgr-bin code unconditionally links against `agc-rs`.

### Batch-and-merge workflow for large haplotype collections

```bash
# Build one sub-archive per batch of 24 haplotypes
bash examples/hprc_r2/01_build_archive.sh

# Merge all sub-archives into the final pangenome archive
bash examples/hprc_r2/02_merge_batches.sh
```

The HPRC R2 example (`examples/hprc_r2/`) demonstrates this strategy for
300+ haplotypes with a `--test` mode (first 10 haplotypes) for quick
validation.

---

## 2. Index Format v3: Transparent, Queryable `.mdbi`/`.mdbv`/`.midx`

The monolithic `.mdb` binary shimmer index is replaced by three files that
are individually inspectable and composable.

| File | Content | Format |
|------|---------|--------|
| `.mdbi` | Sorted shimmer-pair keys | Flat binary (little-endian) |
| `.mdbv` | Fragment location values | Flat binary (little-endian) |
| `.midx` | Sequence metadata + shimmer parameters | **SQLite** |

The `.midx` SQLite schema stores contig name, source, length, and the shimmer
specification (`w`, `k`, `r`, `min_span`, `sketch`) in typed tables with
indices, replacing the old tab-delimited text file that had no schema and
was not queryable.

### Batch-and-merge indexing for memory efficiency

`pgr index mdb` now supports a `--batch-size N` flag (default: 16 whole
haplotypes per batch). Each batch is decompressed, indexed, and written to a
partial shard; shards are then merged in a second pass. This caps peak memory
regardless of archive size and makes it practical to build an index for 100+
haplotype archives on a workstation.

```bash
pgr index mdb --agcrs-input archive.agcrs --prefix out --batch-size 16
```

### Breaking change

Old `.mdb` + tab-delimited `.midx` files are not forward-compatible. Existing
indices must be rebuilt with `pgr index mdb`.

---

## 3. CLI Consolidation: 22 Binaries → Single `pgr` Executable

All standalone `pgr-*` binaries are replaced by a single `pgr` binary with
six subcommand groups. The tool surface is now self-documenting (`pgr --help`
lists all groups; `pgr <group> --help` lists subcommands).

```
pgr
├── index        pgr index mdb          (was pgr-mdb)
│                pgr index shmmr-count  (was pgr-shmmr-count)
│
├── align        pgr align alnmap       (was pgr-alnmap)
│                pgr align map-coord    (was pgr-map-coordinate)
│                pgr align liftover-gtf (new — see §4)
│
├── query        pgr query seqs         (was pgr-query)
│                pgr query fetch        (was pgr-fetch-seqs)
│                pgr query cov          (was pgr-compare-cov)
│                pgr query cov2         (was pgr-compare-cov2)
│
├── bundle       pgr bundle decomp      (was pgr-pbundle-decomp)
│                pgr bundle aln/dist/offset/sort/svg/shmmr-dist
│
├── variant      pgr variant diploid-vcf   (was pgr-generate-diploid-vcf)
│                pgr variant sv-analysis   (was pgr-generate-sv-analysis)
│                pgr variant merge-sv      (was pgr-merge-svcnd-bed)
│                pgr variant annotate-vcf  (was pgr-annotate-vcf-file)
│                pgr variant annotate-bed  (was pgr-annotate-bed-file)
│
└── plot         pgr plot chr-aln       (was pgr-generate-chr-aln-plot)
```

All positional arguments have been replaced with explicit `--long-flags` (with
short aliases) for unambiguous scripting and better shell completion.

### Memory mode for query

`pgr query seqs` gains a `--memory-mode` flag:

```bash
pgr query seqs --memory-mode low|moderate|high
```

`low` reduces thread count to bound peak RSS; `high` maximises parallelism on
machines with ample RAM. Default is `moderate`.

---

## 4. Alignment, Liftover, and End-to-End Annotation Pipeline

### SQLite alignment database (`.alndb`)

`pgr align alnmap` now writes a structured SQLite database alongside the
legacy text outputs (`.alnmap`, `.ctgmap.*`, `.svcnd.*`). The `.alndb` file
contains:

| Table | Content |
|-------|---------|
| `run_params` | Command-line parameters used for this run |
| `sequences` | Reference and query sequence metadata |
| `chains` | Sparse shimmer-pair alignment chains |
| `blocks` | Alignment block details (duplication / overlap flags) |
| `variants` | SNP and indel calls |
| `ctgmap` | Contig-to-reference mappings |
| `sv_candidates` | Structural variant candidate records |
| `alignment_flags` | Miscellaneous flags and metadata |

The legacy text outputs are preserved for backward compatibility.

### GTF liftover (`pgr align liftover-gtf`)

A new subcommand lifts GTF transcript annotations from reference coordinates
to haplotype-contig coordinates using the alignment output. Liftover results
are stored in a SQLite database that records per-transcript and per-exon
coordinates, strand, coverage, and mapping quality.

The hg002 example runs liftover as step 04:

```bash
bash examples/hg002/04_liftover_gtf.sh
```

### SV candidate analysis and gene impact annotation

`pgr variant sv-analysis` is substantially expanded. The end-to-end report
(`generate_e2e_report.py`) now includes a gene-impact subtab for each SV
candidate:

- SVs overlapping a gene body are classified as **exon disrupted**, **exon
  partial**, or **intronic**
- Intergenic SVs report the nearest upstream (5′) and downstream (3′) genes
  with distances
- A genome-wide ideogram plots SV positions by type across all chromosomes

### Performance improvements to alnmap

- Restored `into_par_iter()` on alignment block processing, recovering a
  multi-thread regression present in 0.6.x
- `get_sub_seq_by_id` regression fixed
- Correct native toolchain selection on Apple Silicon (aarch64)
- Technical algorithm description: `AlnMapTechReport/main.tex` (compile with
  `pdflatex` for the full document covering shimmer-based sparse alignment, SV
  candidate detection, and diploid assembly annotation)

---

## 5. End-to-End Example Pipelines

Three self-contained example pipelines cover the full range from quick
evaluation to production pangenome builds.

### E. coli (~5 min, no downloads)

```bash
cd examples/ecoli
bash 01_agcrs_basics.sh   # archive create/append/verify
bash 02_build_index.sh    # shimmer index build
bash 03_query_seqs.sh     # pangenome query with --memory-mode demo
```

### HG002 diploid human genome (~2.5 GB download)

A complete diploid analysis pipeline with six numbered steps:

```
00_download.sh         — fetch reference + assemblies
01_align_alnmap.sh     — sparse alignment for both haplotypes
02_variant_vcf.sh      — diploid VCF generation
03_annotate_vcf.sh     — ClinVar and gene annotation
04_liftover_gtf.sh     — GTF liftover → liftover_report.html
05_query_mhc.sh        — pangenome query of the MHC region
06_generate_report.sh  — composite e2e HTML report
run_all.sh             — orchestrate all steps with timing
```

`run_all.sh` is resumable: each step checks for output sentinels and skips if
already complete. On macOS systems with less than 32 GB RAM it automatically
switches to a chr6-only mode. A timing log (`run_all_timings.tsv`) records
wall time, CPU time, peak RSS, and CPU % per step; entries are upserted on
re-runs so partial reruns produce accurate totals.

The HTML report (`e2e_report.html`) includes:
- Alignment summary and contig mapping overview
- Variant statistics and ClinVar overlap
- SV candidates with genome-wide ideogram and gene-impact breakdown
- GTF liftover status with NM gene funnel, exon-length scatter plots, and
  per-chromosome coverage

### HPRC R2 (~300 GB, large compute)

Batch-append + merge strategy for 100+ haplotype pangenome:

```bash
bash examples/hprc_r2/01_build_archive.sh   # batch-append, 24 haplotypes/batch
bash examples/hprc_r2/02_merge_batches.sh   # merge sub-archives into final
```

Batches are skipped if already built (resume support). A `--test` flag
processes the first 10 haplotypes for quick pipeline validation.

---

## 6. pgr-tk Python Bindings: Maturin + uv + Python 3.13

The Python binding build system is updated:

- **maturin** replaces the old setuptools approach
- **uv** manages the virtual environment and dev dependencies
- **Python 3.13** is the minimum (and recommended) version
- **pyo3 0.22** enables Python 3.13 support with the updated `Bound`-based
  module API

### Build

```bash
cd pgr-tk
uv venv --python 3.13
uv pip install --python .venv/bin/python maturin numpy
env -u CONDA_PREFIX \
    VIRTUAL_ENV=$(pwd)/.venv \
    PYO3_PYTHON=$(pwd)/.venv/bin/python \
    PYTHON_SYS_EXECUTABLE=$(pwd)/.venv/bin/python \
    .venv/bin/maturin develop
```

### Known issue

`load_from_agc_index()` and `load_from_frg_index()` in the Python bindings
are currently broken following the index format migration (`.mdb` →
`.mdbi`/`.mdbv`/`.midx`). `load_from_fastx()` and `load_from_seq_list()`
remain functional. Updated bindings for the new index format are planned for
the next release.

---

## 7. Removed: pgr-web Frontend and Server

`pgr-web/frontend` (Dioxus/WASM) and `pgr-web/pgr-server` are removed from
the repository and workspace. The prototype served its purpose in exploring
what a browser-based interface could look like. A redesigned web interface
will be developed as a separate project once the core toolkit stabilises.

---

## Breaking Changes Summary

| Area | 0.6.x | 0.8.0 |
|------|-------|-------|
| CLI | 22 separate `pgr-*` binaries | Single `pgr <group> <subcommand>` |
| All CLI args | Mix of positional and flags | All explicit `--long-flags` |
| Index files | `.mdb` + tab-delimited `.midx` | `.mdbi` + `.mdbv` + SQLite `.midx` |
| AGC backend | C++ submodule (optional) | `agc-rs` Rust crate (required) |
| Python: min version | 3.6+ | 3.13+ |
| Python: `load_from_agc_index()` | Working | Broken (fix pending) |

Old index files must be rebuilt. Old shell scripts using `pgr-alnmap`,
`pgr-query`, etc. must be updated to the new subcommand syntax.

---

## Dependency Updates

| Crate | 0.6.x | 0.8.0 |
|-------|-------|-------|
| `pyo3` | 0.18.3 | 0.22 |
| `rusqlite` | 0.31 | 0.39 |
| `rustc-hash` | 1 | 2 |
| `thiserror` | 1 | 2 |
| `flate2` | 1.x | 1.1.9 |
| `voracious_radix_sort` | — | 1.2.0 (new) |

Minimum Rust version: 1.80 or later (edition 2021).
