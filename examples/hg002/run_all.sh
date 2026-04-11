#!/usr/bin/env bash
# run_all.sh — run the complete HG002 end-to-end pipeline with timing.
#
# Each step is skipped automatically if its key output file already exists.
# Delete example_output/ (or individual files) to re-run specific steps.
#
# Steps:
#   01  pgr align alnmap        → example_output/hg002_hap{0,1}.alndb
#   02  pgr variant diploid-vcf → example_output/hg002.vcf + aln plots
#   03  VCF annotation          → example_output/hg002.annotated.sorted.clinvar.vcf.gz
#   04  GTF liftover + report   → example_output/liftover_report.html
#   05  MHC pangenome query     → example_output/mhc_hits.000.hit
#
# Run 00_download.sh first if you haven't already downloaded the input data.
#
# Usage:
#   bash examples/hg002/run_all.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT="example_output"
mkdir -p "$OUT"

TIMELOG="$OUT/run_all_timings.tsv"

# Use GNU time if available; fall back to plain execution (no timing data).
if [[ "$(uname)" == "Darwin" ]] && command -v gtime &>/dev/null; then
    _time() { gtime -f "%e\t%U\t%S\t%M\t%P" -a -o "$TIMELOG" "$@"; }
elif command -v /usr/bin/time &>/dev/null && /usr/bin/time --version &>/dev/null 2>&1; then
    _time() { /usr/bin/time -f "%e\t%U\t%S\t%M\t%P" --append -o "$TIMELOG" "$@"; }
else
    _time() { "$@"; }
fi

# run_step NAME SCRIPT SENTINEL...
#   Runs SCRIPT unless all SENTINEL files already exist.
run_step() {
    local name="$1" script="$2"; shift 2
    local sentinels=("$@")

    # Check if all sentinel files exist
    local skip=true
    for f in "${sentinels[@]}"; do
        [[ -f "$f" ]] || { skip=false; break; }
    done

    if $skip; then
        echo
        echo "╌╌ $name — skipped (outputs already exist)"
        return
    fi

    echo
    echo "╔══ $name ══"
    local t0=$SECONDS
    local rc=0
    _time bash "$script" 2>&1 | sed 's/^/│ /' || rc=${PIPESTATUS[0]}
    local elapsed=$(( SECONDS - t0 ))
    if [[ $rc -ne 0 ]]; then
        echo "╚══ $name — FAILED (exit $rc) after ${elapsed}s"
    else
        echo "╚══ $name — ${elapsed}s elapsed"
    fi
    printf "%s\t%s\n" "$name" "${rc:-0}" >> "$TIMELOG"
}

run_step "01 align alnmap" \
    01_align_alnmap.sh \
    "$OUT/hg002_hap0.alndb" \
    "$OUT/hg002_hap1.alndb"

run_step "02 variant vcf + aln plots" \
    02_variant_vcf.sh \
    "$OUT/hg002.vcf" \
    "$OUT/hg002_hap0_aln_plot.html" \
    "$OUT/hg002_hap1_aln_plot.html"

run_step "03 annotate vcf" \
    03_annotate_vcf.sh \
    "$OUT/hg002.annotated.sorted.clinvar.vcf.gz"

run_step "04 liftover gtf + report" \
    04_liftover_gtf.sh \
    "$OUT/hg002_hap0_liftover.db" \
    "$OUT/hg002_hap1_liftover.db" \
    "$OUT/liftover_report.html"

run_step "05 MHC pangenome query" \
    05_query_mhc.sh \
    "$OUT/mhc_hits.000.hit"

echo
echo "=== All steps complete ==="
if [[ -s "$TIMELOG" ]]; then
    echo "Timing log: $TIMELOG"
    column -t "$TIMELOG" 2>/dev/null || cat "$TIMELOG"
fi
