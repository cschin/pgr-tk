# pgr-tk Examples

Choose the starting point that fits your situation:

```
Have ~10 minutes and just want to try pgr-tk?
  → examples/ecoli/        (data already in repo, no downloads)

Want a complete diploid human-genome analysis?
  → examples/hg002/        (download ~2.5 GB, then run end-to-end)

Building a full 100+ haplotype pangenome from HPRC Release 2?
  → examples/hprc_r2/      (download ~300 GB, large compute required)
```

## Directory map

| Directory | Data size | What it covers |
|---|---|---|
| `ecoli/` | ~5 MB (included) | agc-rs archive CRUD · shimmer index · `pgr query seqs` |
| `hg002/` | ~2.5 GB download | full diploid pipeline: align → VCF → annotate → liftover → query |
| `hprc_r2/` | ~300 GB download | build a pangenome archive from HPRC R2 via batch-append + merge |

## Prerequisites

All examples assume you have built the release binaries:

```bash
cargo build --release -p agc-rs -p pgr-bin
```

Override the binary paths at any time by setting environment variables:

```bash
export PGR=/path/to/pgr
export AGC_RS=/path/to/agc-rs
```
