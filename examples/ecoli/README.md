# E. coli examples

Three self-contained scripts that run in order and take about 5 minutes total.
Test data (~5 MB) is already in `test_data/ecoli/` — no downloads needed.

## Steps

| Script | Tool | What it shows |
|---|---|---|
| `01_agcrs_basics.sh` | `agc-rs` | create archive · append · info · list · get · round-trip verify |
| `02_build_index.sh` | `pgr index mdb` | build shimmer index (.mdbi/.mdbv/.midx) from the archive |
| `03_query_seqs.sh` | `pgr query seqs` | query a 50 kb region across all three strains; demonstrates `--memory-mode` |

## Quick start

```bash
cd examples/ecoli
bash 01_agcrs_basics.sh
bash 02_build_index.sh
bash 03_query_seqs.sh
```

Each script is standalone — you can run them individually once the archive
(`ecoli_demo.agcrs`) and index (`ecoli_demo.mdbi` / `.mdbv` / `.midx`) exist.
