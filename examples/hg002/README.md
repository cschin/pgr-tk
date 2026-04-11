# HG002 end-to-end example

A complete diploid human-genome analysis pipeline using HG002 (GIAB/HPRC)
and GRCh38 as the reference.

## Prerequisites

- Binaries built: `cargo build --release -p agc-rs -p pgr-bin`
- ~2.5 GB free disk for downloads; ~20 GB for intermediate files
- `curl`, `python3`, `sqlite3`
- Steps 03 and 05: `bgzip`, `tabix`, `bcftools` (htslib suite)

## Steps

| Script | Input | What it does |
|---|---|---|
| `00_download.sh` | internet | download GRCh38, HG002 mat/pat, RefSeq GTF, ClinVar from HPRC S3 / NCBI |
| `01_align_alnmap.sh` | FASTAs | `pgr align alnmap` — align each haplotype to GRCh38, produce `.alndb` |
| `02_variant_vcf.sh` | `.alndb` × 2 | `pgr variant diploid-vcf` + `pgr plot chr-aln` |
| `03_annotate_vcf.sh` | VCF + GTF | gene annotation, ClinVar cross-reference, bgzip/tabix |
| `04_liftover_gtf.sh` | `.alndb` + GTF | `pgr align liftover-gtf` — map transcripts to HG002 contigs |
| `05_query_mhc.sh` | downloaded FASTAs | build a 2-sample pangenome index, query the MHC region |

## Quick start

```bash
cd examples/hg002
bash 00_download.sh          # ~30 min depending on bandwidth
bash 01_align_alnmap.sh      # ~1 h per haplotype
bash 02_variant_vcf.sh
bash 03_annotate_vcf.sh
bash 04_liftover_gtf.sh
bash 05_query_mhc.sh
```

Or run all steps in sequence with timing instrumentation:

```bash
bash run_all.sh
```

## Disk layout

Downloaded files land in the same directory (gitignored).  Intermediate and
output files (`.alndb`, `.vcf.gz`, `.db`, etc.) also accumulate here — you
can delete them after reviewing the results.

## Reference

See `lift_over_examples.md` for the liftover database schema and example SQL
queries.
