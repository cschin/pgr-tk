#!/usr/bin/env bash
# run_all.sh — run the complete HG002 end-to-end pipeline with timing.
#
# Calls steps 01–05 in order.  Run 00_download.sh first if you haven't
# already downloaded the input data.
#
# All output lands in example_output/
#
# Usage:
#   bash examples/hg002/run_all.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT="example_output"
mkdir -p "$OUT"

TIMELOG="$OUT/run_all_timings.tsv"
> "$TIMELOG"

# Use GNU time (gtime on macOS via brew coreutils, /usr/bin/time on Linux)
if [[ "$(uname)" == "Darwin" ]]; then
    _time() { gtime -f "%e\t%U\t%S\t%M\t%P" -a -o "$TIMELOG" "$@"; }
else
    _time() { /usr/bin/time -f "%e\t%U\t%S\t%M\t%P" --append -o "$TIMELOG" "$@"; }
fi

run_step() {
    local name="$1"; shift
    echo
    echo "╔══ $name ══"
    local t0=$SECONDS
    _time bash "$@" 2>&1 | sed 's/^/│ /'
    echo "╚══ $name — $((SECONDS - t0))s elapsed"
    printf "%s\t" "$name" >> "$TIMELOG"
}

run_step "01 align alnmap"    01_align_alnmap.sh
run_step "02 variant vcf"     02_variant_vcf.sh
run_step "03 annotate vcf"    03_annotate_vcf.sh
run_step "04 liftover gtf"    04_liftover_gtf.sh
run_step "05 query mhc"       05_query_mhc.sh

echo
echo "=== All steps complete ==="
echo "Timing log: $TIMELOG"
column -t "$TIMELOG" 2>/dev/null || cat "$TIMELOG"
