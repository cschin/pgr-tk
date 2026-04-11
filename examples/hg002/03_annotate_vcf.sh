#!/usr/bin/env bash
# 03_annotate_vcf.sh — annotate variants with gene names and ClinVar significance.
#
# Pipeline:
#   1. Strip PanSN prefix (GRCh38#0#) from CHROM so chromosome names are
#      compatible with standard annotation databases.
#   2. pgr variant annotate-vcf — add gene names from RefSeq GTF.
#   3. bgzip + tabix — index the annotated VCF.
#   4. bcftools sort + rename chromosomes (add "chr" prefix for ClinVar).
#   5. bcftools annotate — add CLNSIG / CLNDN from ClinVar.
#
# Requires:
#   hg002.vcf            (from 02_variant_vcf.sh)
#   hg38.ncbiRefSeq.gtf.gz + clinvar.vcf.gz + clinvar.vcf.gz.tbi
#                        (from 00_download.sh)
#   bgzip, tabix, bcftools (htslib suite)
#
# Usage:
#   bash examples/hg002/03_annotate_vcf.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"

if [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr binary not found at $PGR" >&2; exit 1
fi
for tool in bgzip tabix bcftools; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not found"; exit 1; }
done
for f in hg002.vcf hg38.ncbiRefSeq.gtf.gz clinvar.vcf.gz; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found"; exit 1; }
done

# ---------------------------------------------------------------------------
# 1. Strip PanSN prefix from CHROM column (GRCh38#0#chrN → chrN)
# ---------------------------------------------------------------------------
echo "=== [1] Strip PanSN prefix from hg002.vcf ==="
grep  "^#" hg002.vcf  > hg002.stripped.vcf
grep -v "^#" hg002.vcf | sed 's/^[^#]*#[^#]*#//' >> hg002.stripped.vcf

# ---------------------------------------------------------------------------
# 2. Gene annotation with RefSeq GTF
# ---------------------------------------------------------------------------
echo "=== [2] pgr variant annotate-vcf ==="
"$PGR" variant annotate-vcf \
    --vcf-path hg002.stripped.vcf \
    --annotation-path hg38.ncbiRefSeq.gtf.gz \
    --output-path hg002.annotated.vcf

# ---------------------------------------------------------------------------
# 3. bgzip + tabix
# ---------------------------------------------------------------------------
echo "=== [3] bgzip + tabix ==="
bgzip -f hg002.annotated.vcf
tabix -p vcf hg002.annotated.vcf.gz

# ---------------------------------------------------------------------------
# 4. Sort and rename chromosomes for ClinVar compatibility
# ---------------------------------------------------------------------------
echo "=== [4] bcftools sort + rename chromosomes ==="
bcftools sort hg002.annotated.vcf.gz -O z -o hg002.annotated.sorted.vcf.gz
bcftools index -t hg002.annotated.sorted.vcf.gz

# Build chrN → N rename table (strip "chr" prefix to match ClinVar)
bcftools query -f '%CHROM\n' hg002.annotated.sorted.vcf.gz \
    | sort -u \
    | grep "^chr" \
    | awk '{print $0"\t"substr($0,4)}' > chr_rename.txt

bcftools annotate --rename-chrs chr_rename.txt \
    hg002.annotated.sorted.vcf.gz \
    -O z -o hg002.nochr.vcf.gz
bcftools index --tbi hg002.nochr.vcf.gz

# ---------------------------------------------------------------------------
# 5. ClinVar annotation
# ---------------------------------------------------------------------------
echo "=== [5] bcftools annotate with ClinVar ==="
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
    hg002.nochr.vcf.gz \
    -O z -o hg002.annotated.sorted.clinvar.vcf.gz
bcftools index --tbi hg002.annotated.sorted.clinvar.vcf.gz

# ---------------------------------------------------------------------------
# Show annotated ClinVar variants
# ---------------------------------------------------------------------------
echo
echo "=== ClinVar-annotated variants (first 20) ==="
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
    -i 'INFO/CLNSIG!="."' \
    hg002.annotated.sorted.clinvar.vcf.gz \
    | head -20

echo
echo "Output: hg002.annotated.sorted.clinvar.vcf.gz"
echo "Next: bash 04_liftover_gtf.sh"
