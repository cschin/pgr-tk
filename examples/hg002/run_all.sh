#!/usr/bin/env bash
# run_all.sh — run the complete HG002 end-to-end pipeline with timing.
#
# Each step is skipped automatically if its key output file already exists.
# Delete example_output/ (or individual files) to re-run specific steps.
# The timing log is updated per-step (upserted) so re-runs keep accurate
# timings for steps that actually ran and preserve old times for skipped ones.
#
# Steps:
#   01  pgr align alnmap        → example_output/hg002_hap{0,1}.alndb
#   02  pgr variant diploid-vcf → example_output/hg002.vcf + aln plots
#   03  VCF annotation          → example_output/hg002.annotated.sorted.clinvar.vcf.gz
#   04  GTF liftover + report   → example_output/hg002_hap{0,1}_liftover.db + liftover_report.html
#   05  MHC pangenome query     → example_output/mhc_hits.000.hit + mhc_bundle.html
#   06  Final e2e report        → example_output/e2e_report.html (always runs)
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

# Temp file to capture gtime/GNU time output; cleaned up on exit.
_GTIME=$(mktemp)
trap 'rm -f "$_GTIME"' EXIT

# Use GNU time if available; write to temp file (not directly to TIMELOG).
_HAS_TIME=0
if [[ "$(uname)" == "Darwin" ]] && command -v gtime &>/dev/null; then
    _time() { gtime -f "%e\t%U\t%S\t%M\t%P" -o "$_GTIME" "$@"; }
    _HAS_TIME=1
elif command -v /usr/bin/time &>/dev/null && /usr/bin/time --version &>/dev/null 2>&1; then
    _time() { /usr/bin/time -f "%e\t%U\t%S\t%M\t%P" -o "$_GTIME" "$@"; }
    _HAS_TIME=1
else
    _time() { "$@"; }
fi

# _upsert_timelog NAME WALL USER SYS RSS CPU
#   Write (or replace) a single 6-column row keyed by NAME in TIMELOG.
#   Skipped steps are not touched, so their previous timing is preserved.
_upsert_timelog() {
    local name="$1" wall="$2" user="$3" sys_="$4" rss="$5" cpu="$6"
    local tmp
    tmp=$(mktemp)
    # Remove any existing row for this step name (fixed-string tab anchor).
    { [[ -f "$TIMELOG" ]] && grep -Fv "$(printf '%s\t' "$name")" "$TIMELOG"; } \
        > "$tmp" 2>/dev/null || true
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$wall" "$user" "$sys_" "$rss" "$cpu" >> "$tmp"
    mv "$tmp" "$TIMELOG"
}

# run_step NAME SCRIPT SENTINEL...
#   Runs SCRIPT unless all SENTINEL files already exist.
#   Pass no sentinels to always run (e.g. the final report step).
run_step() {
    local name="$1" script="$2"; shift 2
    local sentinels=("$@")

    # No sentinels → always run.  Otherwise skip if all outputs exist.
    local skip=false
    if (( ${#sentinels[@]} > 0 )); then
        skip=true
        for f in "${sentinels[@]}"; do
            [[ -f "$f" ]] || { skip=false; break; }
        done
    fi

    if $skip; then
        echo
        echo "╌╌ $name — skipped (outputs already exist)"
        return
    fi

    echo
    echo "╔══ $name ══"
    local t0=$SECONDS
    local rc=0
    > "$_GTIME"   # clear previous gtime output
    _time bash "$script" 2>&1 | sed 's/^/│ /' || rc=${PIPESTATUS[0]}
    local elapsed=$(( SECONDS - t0 ))
    if [[ $rc -ne 0 ]]; then
        echo "╚══ $name — FAILED (exit $rc) after ${elapsed}s"
    else
        echo "╚══ $name — ${elapsed}s elapsed"
    fi

    # Upsert timing row (one row per step, updated on each actual run).
    if (( _HAS_TIME )) && [[ -s "$_GTIME" ]]; then
        IFS=$'\t' read -r _wall _user _sys _rss _cpu < "$_GTIME"
        _upsert_timelog "$name" "$_wall" "$_user" "$_sys" "$_rss" "$_cpu"
    else
        # gtime unavailable: record wall-clock seconds, zeroes for the rest.
        _upsert_timelog "$name" "$elapsed" "0" "0" "0" "N/A"
    fi
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

run_step "06 final e2e report" \
    06_generate_report.sh

echo
echo "=== All steps complete ==="
if [[ -s "$TIMELOG" ]]; then
    echo "Timing log: $TIMELOG"
    column -t "$TIMELOG" 2>/dev/null || cat "$TIMELOG"
fi
