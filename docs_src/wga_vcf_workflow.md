# Whole Genome Alignment and VCF Generation with PGR-TK

This guide explains how to use PGR-TK to align a new genome assembly against a reference,
generate a diploid VCF file, and prepare it for annotation with ClinVar and other databases
(e.g., OMIM, gnomAD, dbSNP).

---

## Overview

The workflow has four steps:

```
New Assembly (FASTA, two haplotypes)
         |
         v
[1] pgr-alnmap          -- align each haplotype against the reference
         |
         v
[2] pgr-generate-diploid-vcf   -- merge both alnmap files into a diploid VCF
         |
         v
[3] (optional) pgr-annotate-vcf-file  -- annotate with gene names from GTF
         |
         v
[4] External annotation  -- ClinVar, gnomAD, dbSNP via bcftools annotate
```

---

## Prerequisites

Build the PGR-TK binaries:

```bash
cargo build --release
```

The tools used in this workflow are:

- `pgr-alnmap` — whole genome alignment, produces `.alnmap` files
- `pgr-generate-diploid-vcf` — merges two haplotype `.alnmap` files into a VCF
- `pgr-annotate-vcf-file` — optional gene-name annotation using a GTF file

All binaries are in `target/release/`.

---

## Step 1: Align Each Haplotype Against the Reference

Run `pgr-alnmap` once per haplotype. It uses SHIMMER (sparse hierarchical minimizer matching)
anchoring followed by Wavefront Alignment (WFA) or Smith-Waterman to call base-level variants.

```bash
# Haplotype 0
pgr-alnmap \
    reference.fasta \
    assembly_hap0.fasta \
    sample_hap0 \
    --preset default

# Haplotype 1
pgr-alnmap \
    reference.fasta \
    assembly_hap1.fasta \
    sample_hap1 \
    --preset default
```

### Preset options

| Preset    | w  | k  | r | Use case                          |
|-----------|----|----|---|-----------------------------------|
| `fast`    | 80 | 55 | 4 | Quick pass, lower sensitivity     |
| `default` | 48 | 55 | 2 | Balanced (recommended)            |
| `detail`  | 48 | 55 | 2 | Higher sensitivity, slower        |

Parameters can also be set manually with `--w`, `--k`, `--r`, `--min-span`.

### Output files

| File                     | Contents                                              |
|--------------------------|-------------------------------------------------------|
| `sample_hap0.alnmap`     | Alignment blocks and called variants (M, V, S records)|
| `sample_hap0.ctgmap.json`| Contig-to-reference mapping with sequence lengths     |
| `sample_hap0.svcnd.bed`  | Structural variant candidate regions                  |
| `sample_hap0_svcnd/`     | FASTA sequences spanning SV candidates (if enabled)   |

The `.ctgmap.json` from either haplotype (they reference the same target) is used in Step 2.

---

## Step 2: Generate the Diploid VCF

`pgr-generate-diploid-vcf` reads the two `.alnmap` files, groups overlapping variant records,
resolves alleles per haplotype, deduplicates, and emits VCFv4.2.

```bash
pgr-generate-diploid-vcf \
    sample_hap0.alnmap \
    sample_hap1.alnmap \
    sample_hap0.ctgmap.json \
    sample \
    --sample-name SAMPLE_ID
```

### Output files

| File           | Contents                                                    |
|----------------|-------------------------------------------------------------|
| `sample.vcf`   | Diploid VCF with phased genotypes                           |
| `sample.bed`   | BED file of regions covered by both haplotypes              |

### VCF format details

The VCF is VCFv4.2. Each record has:

- `CHROM`, `POS` (1-based), `REF`, `ALT` — standard fields
- `QUAL` — 40 for PASS variants, 30 for flagged records
- `FILTER` — one of:
  - `PASS` — clean, uniquely aligned variant
  - `DUP` — variant from a duplicated alignment block
  - `OVLP` — variant from an overlapping alignment chain
  - `NC` — no confident diploid call (one haplotype not covered)
- `FORMAT/GT` — phased genotype, e.g. `0|1`, `1|1`

Example records:

```
#CHROM  POS     ID  REF   ALT   QUAL  FILTER  INFO  FORMAT  SAMPLE_ID
chr1    925952  .   A     G     40    PASS    .     GT      0|1
chr1    931279  .   ACGT  A     40    PASS    .     GT      1|1
chr1    982736  .   T     TGCA  30    DUP     .     GT      .|1
```

---

## Step 3: (Optional) Annotate with Gene Names

`pgr-annotate-vcf-file` adds a `GN` INFO field with the gene name(s) overlapping each variant.
It accepts a gzip-compressed GTF file (NCBI RefSeq or GENCODE format).

Download a GTF annotation for GRCh38:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
```

Run annotation:

```bash
pgr-annotate-vcf-file \
    sample.vcf \
    hg38.ncbiRefSeq.gtf.gz \
    sample.annotated.vcf
```

The output adds `INFO/GN=<gene_name>` to any variant overlapping a transcript. Variants outside
annotated transcripts are dropped from the output. To keep all variants, use this output only
for gene-based filtering and keep the original VCF for database annotation.

---

## Step 4: Annotate Against ClinVar, gnomAD, and dbSNP

The VCF from Step 2 (or Step 3) is standard VCFv4.2 and can be processed with any standard
annotation pipeline.

### Sort and compress the VCF

```bash
bcftools sort sample.vcf -O z -o sample.sorted.vcf.gz
bcftools index -t sample.sorted.vcf.gz
```

### ClinVar annotation

Download the ClinVar VCF:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

Annotate:

```bash
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
    sample.sorted.vcf.gz \
    -O z -o sample.clinvar.vcf.gz
```

Key ClinVar INFO fields added:

| Field         | Description                                              |
|---------------|----------------------------------------------------------|
| `CLNSIG`      | Clinical significance (Pathogenic, Benign, VUS, etc.)    |
| `CLNDN`       | Disease name                                             |
| `CLNREVSTAT`  | Review status (number of submitters, expert panel, etc.) |

### dbSNP annotation (rsID)

```bash
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi

bcftools annotate \
    -a 00-All.vcf.gz \
    -c ID \
    sample.sorted.vcf.gz \
    -O z -o sample.dbsnp.vcf.gz
```

### gnomAD annotation

```bash
# gnomAD v4 genomes (adjust chromosome files as needed)
bcftools annotate \
    -a gnomad.genomes.v4.vcf.bgz \
    -c INFO/AF,INFO/AF_popmax \
    sample.sorted.vcf.gz \
    -O z -o sample.gnomad.vcf.gz
```

### Combining all annotations

```bash
bcftools annotate -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN \
    sample.sorted.vcf.gz | \
bcftools annotate -a 00-All.vcf.gz \
    -c ID | \
bcftools annotate -a gnomad.genomes.v4.vcf.bgz \
    -c INFO/AF,INFO/AF_popmax \
    -O z -o sample.fully_annotated.vcf.gz

bcftools index -t sample.fully_annotated.vcf.gz
```

---

## Filtering Recommendations

After annotation, apply filters to focus on clinically or functionally relevant variants:

```bash
# Keep only PASS variants with pathogenic or likely-pathogenic ClinVar classification
bcftools filter \
    -i 'FILTER="PASS" && INFO/CLNSIG~"Pathogenic"' \
    sample.fully_annotated.vcf.gz

# Keep rare variants (gnomAD AF < 0.01) that PASS alignment QC
bcftools filter \
    -i 'FILTER="PASS" && INFO/AF < 0.01' \
    sample.fully_annotated.vcf.gz

# Exclude duplicated or overlapped alignment artifacts
bcftools filter \
    -e 'FILTER="DUP" || FILTER="OVLP"' \
    sample.fully_annotated.vcf.gz
```

---

## Notes on Haploid Assemblies

If only a single (haploid) assembly is available, run `pgr-alnmap` once and pass the same
`.alnmap` file for both haplotype inputs to `pgr-generate-diploid-vcf`:

```bash
pgr-alnmap reference.fasta assembly.fasta sample sample --preset default

pgr-generate-diploid-vcf \
    sample.alnmap \
    sample.alnmap \
    sample.ctgmap.json \
    sample \
    --sample-name SAMPLE_ID
```

Genotypes will be homozygous (`1|1`) for all called variants. Use the `sample.bed` output
to restrict downstream analysis to covered regions.

---

## Summary of Commands

```bash
# 1. Align both haplotypes
pgr-alnmap reference.fasta hap0.fasta out_hap0 --preset default
pgr-alnmap reference.fasta hap1.fasta out_hap1 --preset default

# 2. Generate diploid VCF
pgr-generate-diploid-vcf \
    out_hap0.alnmap out_hap1.alnmap \
    out_hap0.ctgmap.json \
    sample --sample-name MY_SAMPLE

# 3. (Optional) Annotate with gene names
pgr-annotate-vcf-file sample.vcf hg38.ncbiRefSeq.gtf.gz sample.gene_annotated.vcf

# 4. Sort, compress, index
bcftools sort sample.vcf -O z -o sample.vcf.gz && bcftools index -t sample.vcf.gz

# 5. Annotate with ClinVar
bcftools annotate -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
    sample.vcf.gz -O z -o sample.clinvar.vcf.gz
```
