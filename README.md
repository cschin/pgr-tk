# PGR-TK: A PanGenomic Research Tool Kit

[![test_and_build](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml/badge.svg)](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml)

PGR-TK provides pangenome assembly management, query, and Minimizer Anchored
Pangenome (MAP) Graph generation with principal bundle decomposition. It is
designed for the analysis of large diploid and pangenomic assemblies, from
sequence archiving through structural variant calling and gene annotation
liftover.

Research preprint:
[Multiscale Analysis of Pangenome Enables Improved Representation of Genomic Diversity For Repetitive And Clinically Relevant Genes](https://www.biorxiv.org/content/10.1101/2022.08.05.502980v2)

![Pangenome Data Management and Minimizer Anchored Pangenome Graph Generation](/images/PGR_TK_Sketch_MAPG_construction.png)

With the MAP graph, principal bundle decomposition can reveal complicated
structural variants and genome rearrangements across a population.

![AMY1A Example](/images/AMY1A_example.png)

---

## What's New in 0.8.0

### agc-rs: Pure-Rust genome compression replaces C++ AGC

The C++ AGC submodule is removed. `agc-rs`, a faithful Rust port with
comparable compression performance, is now the exclusive backend. Sequences
are stored in a SQLite-backed `.agcrs` archive that is inspectable with
standard tools and accessible from any SQLite client.

New subcommands:
- `agc-rs batch-append` — compress multiple samples in one pass, sharing the
  delta-index build across all contigs in the batch
- `agc-rs merge` — combine pre-built sub-archives into a single archive,
  enabling a parallelisable batch-then-merge strategy for very large
  haplotype collections (see `examples/hprc_r2/`)

### Index format v3: transparent `.mdbi` / `.mdbv` / SQLite `.midx`

The monolithic `.mdb` binary shimmer index is replaced by three files:
`.mdbi` (sorted keys), `.mdbv` (fragment location values), and `.midx`
(SQLite metadata with typed schema). The new format is externally queryable,
composable, and supports a batch-and-merge memory-efficient build strategy
(`--batch-size` flag on `pgr index mdb`).

> **Breaking change:** old `.mdb` and tab-delimited `.midx` files are not
> compatible. Rebuild existing indices with `pgr index mdb`.

### Unified `pgr` CLI: 22 binaries consolidated into one

All `pgr-*` standalone binaries are replaced by a single `pgr` executable
with six subcommand groups (`index`, `align`, `query`, `bundle`, `variant`,
`plot`). All arguments are explicit `--long-flags`. Run `pgr --help` for
the full command map.

### GTF liftover and end-to-end annotation pipeline

New `pgr align liftover-gtf` subcommand lifts reference transcript
annotations to haplotype-contig coordinates. Combined with the alignment
(`.alndb`), variant calling, ClinVar annotation, and SV analysis steps, this
provides a complete diploid annotation pipeline with an interactive HTML
report (see `examples/hg002/`).

`pgr align alnmap` now writes a queryable SQLite `.alndb` database alongside
legacy text outputs.

### SV annotation with gene impact

The end-to-end report classifies SV candidates by gene impact (exon
disrupted / exon partial / intronic / intergenic with nearest-neighbour
genes) and includes a genome-wide ideogram.

### Python bindings: maturin + uv + Python 3.13

The Python binding build system migrates to maturin and uv; Python 3.13 is
now required. See the [Python bindings section](#python-bindings) below.

### pgr-web removed

The prototype web frontend and server are removed. A redesigned interface
will follow in a future release.

---

## Quick Start

### Prerequisites

```bash
cargo build --release -p agc-rs -p pgr-bin
# add target/release to PATH, or set:
export PGR=$(pwd)/target/release/pgr
export AGC_RS=$(pwd)/target/release/agc-rs
```

### Try it in 5 minutes (E. coli, no downloads)

```bash
cd examples/ecoli
bash 01_agcrs_basics.sh   # create archive, append, round-trip verify
bash 02_build_index.sh    # build shimmer index
bash 03_query_seqs.sh     # pangenome query with --memory-mode demo
```

### Typical workflow

```bash
# 1. Build an AGC archive
agc-rs create pangenome.agcrs --sample GRCh38 GRCh38.fa
agc-rs append pangenome.agcrs --sample HG002_mat HG002_mat.fa
agc-rs append pangenome.agcrs --sample HG002_pat HG002_pat.fa

# 2. Build the shimmer index (batch-size 16 haplotypes per pass)
pgr index mdb --agcrs-input pangenome.agcrs --batch-size 16

# 3. Query a region of interest
pgr query seqs --pgr-db-prefix pangenome --query-fastx-path query.fa \
    --output-prefix output --max-count 128 --min-anchor-count 10
```

### Large pangenome (100+ haplotypes, batch-and-merge)

```bash
# Build sub-archives in batches of 24, then merge
bash examples/hprc_r2/01_build_archive.sh
bash examples/hprc_r2/02_merge_batches.sh
```

---

## Command Line Reference

### Sequence archive management (`agc-rs`)

```
agc-rs create  <archive.agcrs> --sample <name> <sequences.fa>
agc-rs append  <archive.agcrs> --sample <name> <sequences.fa>
agc-rs batch-append <archive.agcrs> name1:file1.fa name2:file2.fa ...
agc-rs merge   --output merged.agcrs input1.agcrs input2.agcrs ...
agc-rs info    <archive.agcrs>
agc-rs list    <archive.agcrs>
agc-rs get     <archive.agcrs> --sample <name> [--contig <ctg>]
```

### `pgr index` — build shimmer indices

| Subcommand | Description |
|-----------|-------------|
| `pgr index mdb` | Build `.mdbi`/`.mdbv`/`.midx` shimmer index from an `.agcrs` archive |
| `pgr index shmmr-count` | Count shimmer occurrences across three sequence sets |

```bash
pgr index mdb --agcrs-input archive.agcrs [--prefix out] [--batch-size 16]
```

### `pgr align` — genome alignment and liftover

| Subcommand | Description |
|-----------|-------------|
| `pgr align alnmap` | Align long contigs to a reference; produces `.alndb`, `.alnmap`, `.ctgmap.*`, `.svcnd.*` |
| `pgr align map-coord` | Map query coordinates to target via an alnmap/alndb file |
| `pgr align liftover-gtf` | Lift GTF transcript annotations from reference to haplotype contigs |

```bash
pgr align alnmap -R ref.fa -Q assembly.fa -O output_prefix
pgr align liftover-gtf --alndb hap0.alndb --gtf annotation.gtf --output hap0_liftover.db
```

### `pgr query` — search and fetch

| Subcommand | Description |
|-----------|-------------|
| `pgr query seqs` | Query a pangenome index; outputs hit summary and FASTA |
| `pgr query fetch` | Fetch sequences from an archive by region |
| `pgr query cov` / `cov2` | Compare shimmer-pair coverage between two sequence sets |

```bash
pgr query seqs --pgr-db-prefix pangenome --query-fastx-path query.fa \
    --output-prefix out --memory-mode moderate
```

`--memory-mode low|moderate|high` controls the parallelism/RSS trade-off.

### `pgr bundle` — principal bundle decomposition

| Subcommand | Description |
|-----------|-------------|
| `pgr bundle decomp` | Generate principal bundle decomposition via MAP Graph |
| `pgr bundle svg` | Generate interactive SVG from a bundle BED file |
| `pgr bundle sort` | Generate contig sort order from bundle decomposition |
| `pgr bundle dist` | Compute pairwise alignment scores from a bundle BED file |
| `pgr bundle offset` | Compute bundle offsets |
| `pgr bundle shmmr-dist` | Shimmer-pair distance between bundle entries |
| `pgr bundle aln` | Bundle-guided alignment |

### `pgr variant` — variant calling and annotation

| Subcommand | Description |
|-----------|-------------|
| `pgr variant diploid-vcf` | Merge two haplotype `.alndb` files into a phased diploid VCF |
| `pgr variant annotate-vcf` | Annotate a VCF with gene names from a GTF file |
| `pgr variant annotate-bed` | Annotate BED regions with gene annotation features |
| `pgr variant sv-analysis` | Analyse SV candidates with principal bundle decomposition |
| `pgr variant merge-sv` | Merge SV candidate BED records from multiple haplotypes |

### `pgr plot` — visualisation

| Subcommand | Description |
|-----------|-------------|
| `pgr plot chr-aln` | Generate chromosome alignment SVG plots from an alnmap JSON file |

Run `pgr <group> <subcommand> --help` for detailed usage of any subcommand.

---

## Examples

| Directory | Data | What it covers |
|-----------|------|----------------|
| `examples/ecoli/` | ~5 MB (included) | agc-rs archive CRUD · shimmer index · `pgr query seqs` |
| `examples/hg002/` | ~2.5 GB download | Full diploid pipeline: align → VCF → annotate → liftover → query → HTML report |
| `examples/hprc_r2/` | ~300 GB download | Batch-append + merge for 100+ haplotype pangenome |

### HG002 end-to-end pipeline steps

```
00_download.sh         fetch reference + assemblies
01_align_alnmap.sh     sparse alignment for both haplotypes
02_variant_vcf.sh      diploid VCF generation
03_annotate_vcf.sh     ClinVar and gene annotation
04_liftover_gtf.sh     GTF liftover → liftover_report.html
05_query_mhc.sh        pangenome query of the MHC region
06_generate_report.sh  composite e2e HTML report
run_all.sh             orchestrate all steps with timing log
```

`run_all.sh` is resumable (sentinel-based step skipping), auto-switches to
chr6-only mode on machines with less than 32 GB RAM, and records wall time,
CPU time, peak RSS, and CPU % per step in `run_all_timings.tsv`.

---

## Python Bindings

The Python bindings are built with [maturin](https://github.com/PyO3/maturin)
and managed with [uv](https://github.com/astral-sh/uv). Python 3.13 or later
is required.

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

> **Note:** `load_from_agc_index()` and `load_from_frg_index()` are currently
> broken following the index format migration. `load_from_fastx()` and
> `load_from_seq_list()` remain functional. Updated bindings are planned for
> the next release.

API documentation: https://genedx.github.io/pgr-tk/  
Jupyter Notebooks: https://github.com/genedx/pgr-tk-notebooks/

---

## Build

### Native build

```bash
cargo build --release -p agc-rs -p pgr-bin
```

Binaries land in `target/release/`. On Apple Silicon, the build system
automatically selects the native aarch64 toolchain.

### Docker (Linux / Ubuntu 22.04)

```bash
git clone git@github.com:cschin/pgr-tk.git
cd pgr-tk/docker
ln -s Dockerfile.build_env-20.04 Dockerfile
docker build -t pgr-tk-build .
docker run -it --rm -v $PWD:/wd/pgr-tk pgr-tk-build /bin/bash
# inside container:
cd /wd/pgr-tk && bash build.sh
```

### Singularity

```bash
# 1. Commit and push the Docker container
docker commit <container_id> <image_name>:<version>
docker push <image_name>:<version>

# 2. Build Singularity image
singularity build ./pgr-tk.<version>.sif docker://<docker_repo>/<image_name>:<version>

# 3. Run
singularity exec --fakeroot -B <host_path>:/<container_path> \
    ./pgr-tk.<version>.sif pgr index mdb --agcrs-input pangenome.agcrs
```

---

## Release Notes

See [pgr-tk-0.8-release-note.md](pgr-tk-0.8-release-note.md) for the full
0.8.0 changelog, breaking changes, and migration guide.
