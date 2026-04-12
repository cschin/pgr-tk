#!/usr/bin/env bash
# 02_build_index.sh — build a PGR-TK shimmer index from the E. coli archive.
#
# Requires: example_output/ecoli_demo.agcrs (created by 01_agcrs_basics.sh)
#
# Produces in example_output/:
#   ecoli_demo.mdbi   — sorted shimmer-pair → value-file offset index
#   ecoli_demo.mdbv   — packed fragment-signature values
#   ecoli_demo.midx   — SQLite sequence metadata (names, lengths, AGC path)
#
# Usage:
#   bash examples/ecoli/02_build_index.sh

set -euo pipefail
cd "$(dirname "$0")"

# NOTE: This script requires pgr to be installed in your PATH.
# Install from the repository root with:
#   cargo install --path ../../pgr-bin
PGR="${PGR:-pgr}"
OUT="example_output"
ARCHIVE="$OUT/ecoli_demo.agcrs"
PREFIX="$OUT/ecoli_demo"

if ! command -v "$PGR" &>/dev/null && [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr not found in PATH" >&2
    echo "       Install with: cargo install --path <repo-root>/pgr-bin" >&2
    exit 1
fi

if [[ ! -f "$ARCHIVE" ]]; then
    echo "ERROR: $ARCHIVE not found — run 01_agcrs_basics.sh first" >&2
    exit 1
fi

mkdir -p "$OUT"

echo "=== Building shimmer index from $ARCHIVE ==="
echo "    Output prefix: $PREFIX"
echo "    Default parameters: w=80 k=56 r=4 min_span=64"
echo

# For a 3-strain E. coli archive the entire index fits in one batch.
# (--batch-size 0 disables sharding; fine for small archives.)
"$PGR" index mdb \
    --agcrs-input "$ARCHIVE" \
    --prefix "$PREFIX" \
    --batch-size 0

echo
echo "Index files written:"
ls -lh "${PREFIX}.mdbi" "${PREFIX}.mdbv" "${PREFIX}.midx"
echo
echo "Done."
