# HPRC Release 2 — full pangenome archive

Build a compressed pangenome archive from all HPRC Release 2 assemblies
(~100+ haplotypes, ~3 GB each) using batched compression and archive merging.

## Requirements

| Resource | Estimate |
|---|---|
| Disk (peak) | `BATCH_SIZE × 3 GB` + final archive (~60 GB compressed) |
| Disk (default, BATCH_SIZE=24) | ~130 GB peak, ~60 GB final |
| RAM | ~16 GB during batch-append |
| Network | ~300 GB download total |
| Time | 12–48 h depending on bandwidth and core count |

## Quick start

```bash
cd examples/hprc_r2
bash 01_build_archive.sh hprc_r2.agcrs
```

**Test mode** (first 10 haplotypes only, ~30 GB download, ~10 min):

```bash
bash 01_build_archive.sh --test hprc_r2_test.agcrs
```

## Tuning

| Variable | Default | Effect |
|---|---|---|
| `BATCH_SIZE` | 24 | haplotypes per sub-archive; lower = less peak disk |
| `KEEP_BATCHES=1` | 0 | keep sub-archives after merge (useful for debugging) |

Example — smaller batches for tighter disk budget:

```bash
BATCH_SIZE=8 bash 01_build_archive.sh hprc_r2.agcrs
```

## After building the archive

Index the archive for sequence queries:

```bash
../../target/release/pgr index mdb \
    --agcrs-input hprc_r2.agcrs \
    --prefix hprc_r2 \
    --batch-size 16
```

Then query any region — see `examples/hg002/05_query_mhc.sh` for an example.
