#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

REPORT="e2e_report.html"
TIMELOG="e2e_timings.tsv"
GTIME_FMT="%e\t%U\t%S\t%M\t%P"

> "$TIMELOG"

# Always generate the HTML report on exit (even on failure)
trap 'generate_report' EXIT

generate_report() {
    python3 "$(dirname "$0")/generate_e2e_report.py" \
        --timelog "$TIMELOG" \
        --out "$REPORT"
}

run_step() {
    local step_name="$1"; shift
    local t_file; t_file=$(mktemp)
    echo "=== $step_name ==="
    gtime -f "$GTIME_FMT" -o "$t_file" "$@"
    printf "%s\t%s\n" "$step_name" "$(cat "$t_file")" >> "$TIMELOG"
    rm "$t_file"
}

PGR=../target/release/pgr

# ---------------------------------------------------------------------------
# Step 1 — alignment map
# ---------------------------------------------------------------------------
run_step "pgr align alnmap hap0" \
    "$PGR" align alnmap \
        --reference-fasta-path GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.b.fa.gz \
        --assembly-contig-path hg002v1.1.mat_MT.PanSN.b.fa.gz \
        --output-prefix hg002_hap0 \
        --preset default

run_step "pgr align alnmap hap1" \
    "$PGR" align alnmap \
        --reference-fasta-path GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.b.fa.gz \
        --assembly-contig-path hg002v1.1.pat.PanSN.b.fa.gz \
        --output-prefix hg002_hap1 \
        --preset default

# ---------------------------------------------------------------------------
# Step 2 — generate diploid VCF
# ---------------------------------------------------------------------------
run_step "pgr variant diploid-vcf" \
    "$PGR" variant diploid-vcf \
        --hap0-path hg002_hap0.alndb \
        --hap1-path hg002_hap1.alndb \
        --output-prefix hg002 \
        --sample-name hg002

# ---------------------------------------------------------------------------
# Step 2b — whole genome alignment plots
# ---------------------------------------------------------------------------
run_step "pgr plot chr-aln hap0" \
    "$PGR" plot chr-aln \
        --ctgmap-json-path hg002_hap0.ctgmap.json \
        --output-prefix hg002_hap0_aln_plot

run_step "pgr plot chr-aln hap1" \
    "$PGR" plot chr-aln \
        --ctgmap-json-path hg002_hap1.ctgmap.json \
        --output-prefix hg002_hap1_aln_plot

# ---------------------------------------------------------------------------
# Step 3a — strip PanSN prefix (GRCh38#0#) from CHROM column
# ---------------------------------------------------------------------------
echo "=== strip PanSN prefix from hg002.vcf ==="
grep "^#" hg002.vcf > hg002.stripped.vcf
grep -v "^#" hg002.vcf | sed 's/^[^#]*#[^#]*#//' >> hg002.stripped.vcf

# ---------------------------------------------------------------------------
# Step 3b — annotate VCF with GTF
# ---------------------------------------------------------------------------
run_step "pgr variant annotate-vcf" \
    "$PGR" variant annotate-vcf \
        --vcf-path hg002.stripped.vcf \
        --annotation-path hg38.ncbiRefSeq.gtf.gz \
        --output-path hg002.annotated.vcf

# ---------------------------------------------------------------------------
# Step 3c — bgzip and tabix index the annotated VCF
# ---------------------------------------------------------------------------
run_step "bgzip annotated VCF" \
    bgzip -f hg002.annotated.vcf

run_step "tabix annotated VCF" \
    tabix -p vcf hg002.annotated.vcf.gz

# ---------------------------------------------------------------------------
# Step 4 — sort, strip chr, annotate with ClinVar
# ---------------------------------------------------------------------------
run_step "bcftools sort" \
    bcftools sort hg002.annotated.vcf.gz -O z -o hg002.annotated.sorted.vcf.gz

run_step "bcftools index sorted" \
    bcftools index -t hg002.annotated.sorted.vcf.gz

echo "=== building chr_rename.txt ==="
bcftools query -f '%CHROM\n' hg002.annotated.sorted.vcf.gz \
    | sort -u \
    | grep "^chr" \
    | awk '{print $0"\t"substr($0,4)}' > chr_rename.txt

run_step "bcftools rename-chrs" \
    bcftools annotate --rename-chrs chr_rename.txt hg002.annotated.sorted.vcf.gz \
        -O z -o hg002.nochr.vcf.gz

run_step "bcftools index nochr" \
    bcftools index --tbi hg002.nochr.vcf.gz

run_step "bcftools annotate ClinVar" \
    bcftools annotate \
        -a clinvar.vcf.gz \
        -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
        hg002.nochr.vcf.gz \
        -O z -o hg002.annotated.sorted.clinvar.vcf.gz

run_step "bcftools index clinvar" \
    bcftools index --tbi hg002.annotated.sorted.clinvar.vcf.gz

# ---------------------------------------------------------------------------
# Show annotated variants
# ---------------------------------------------------------------------------
echo ""
echo "=== Annotated ClinVar variants (first 20) ==="
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
    -i 'INFO/CLNSIG!="."' hg002.annotated.sorted.clinvar.vcf.gz | head -20

# ---------------------------------------------------------------------------
# Step 5 — transcript liftover + liftover report
# ---------------------------------------------------------------------------
echo "=== transcript liftover ==="
bash "$(dirname "$0")/get_tx_seqs.sh"

echo "Done."
