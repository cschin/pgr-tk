#!/usr/bin/env bash
# 02_merge_batches.sh — merge batch sub-archives produced by 01_build_archive.sh
# into a single output archive.
#
# Use this script when the batch phase completed but the merge step still
# needs to run (e.g. after fixing a merge issue, or when KEEP_BATCHES=1 was
# set during the original build).
#
# Resume support: if the output archive already exists the script exits
# immediately after printing info.
#
# Usage:
#   bash examples/hprc_r2/02_merge_batches.sh <archive.agcrs> [agc-rs-bin] [work-dir]
#
# Arguments:
#   archive.agcrs  — final merged output archive path
#   agc-rs-bin     — path to the agc-rs binary (default: ../../target/release/agc-rs)
#   work-dir       — directory containing batch_NNNN.agcrs files
#                    (default: ./hprc_r2_work)

set -euo pipefail

ARCHIVE="${1:?usage: $0 <archive.agcrs> [agc-rs-binary] [work-dir]}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
AGC_RS="${2:-${SCRIPT_DIR}/../../target/release/agc-rs}"
WORK_DIR="${3:-$(pwd)/hprc_r2_work}"

if [[ ! -x "$AGC_RS" ]]; then
    echo "ERROR: agc-rs binary not found at $AGC_RS" >&2
    echo "       Build with: cargo build -p agc-rs --release" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Short-circuit: if the final archive already exists, just print info and exit
# ---------------------------------------------------------------------------
if [[ -f "$ARCHIVE" ]]; then
    echo "[INFO] Archive already exists: $ARCHIVE"
    "$AGC_RS" info "$ARCHIVE"
    exit 0
fi

# ---------------------------------------------------------------------------
# Discover batch sub-archives in order
# ---------------------------------------------------------------------------
mapfile -t BATCH_ARCHIVES < <(ls -1 "${WORK_DIR}"/batch_[0-9]*.agcrs 2>/dev/null | sort)

if (( ${#BATCH_ARCHIVES[@]} == 0 )); then
    echo "ERROR: no batch_NNNN.agcrs files found in $WORK_DIR" >&2
    exit 1
fi

echo "[MERGE] Found ${#BATCH_ARCHIVES[@]} sub-archive(s) in $WORK_DIR:"
for f in "${BATCH_ARCHIVES[@]}"; do
    echo "         $f"
done

# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------
echo ""
echo "[MERGE] Merging ${#BATCH_ARCHIVES[@]} sub-archive(s) → $ARCHIVE ..."
"$AGC_RS" merge --output "$ARCHIVE" "${BATCH_ARCHIVES[@]}"

echo ""
echo "[DONE] Archive: $ARCHIVE"
"$AGC_RS" info "$ARCHIVE"
