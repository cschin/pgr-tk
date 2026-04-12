#!/usr/bin/env bash
# 01_build_archive.sh — build an agc-rs pangenome archive from the HPRC
# release-2 assembly table by processing haplotypes in batches, then merging.
#
# Strategy:
#   - The first haplotype in the work list is the common reference for every
#     batch sub-archive.  It is downloaded once and kept until the final merge.
#   - For each batch of BATCH_SIZE haplotypes a sub-archive is created:
#       1. agc-rs create  (seeds it with the common reference)
#       2. agc-rs batch-append  (compresses the batch with a shared index,
#          avoiding N redundant index-build passes)
#       3. Downloaded FASTAs are deleted immediately after the batch is done.
#   - All sub-archives are merged into the final output archive.
#   - Sub-archives and the cached reference FA are removed after a successful
#     merge (set KEEP_BATCHES=1 to disable cleanup).
#
# Disk headroom needed during a batch: BATCH_SIZE × genome_size (e.g. 24 ×
# 3 GB ≈ 72 GB for human).  Adjust BATCH_SIZE to fit available disk space
# (larger batches amortize the shared-index build across more samples).
#
# Resume support: if a batch sub-archive already exists it is skipped; if the
# final archive already exists the script exits immediately after printing info.
#
# Usage:
#   bash examples/hprc_r2/01_build_archive.sh [--test] <archive.agcrs> [agc-rs-bin] [work-dir]
#
# Options:
#   --test         Process only the first 10 haplotypes (quick smoke test)
#
# Arguments:
#   archive.agcrs  — final merged output archive path
#   agc-rs-bin     — path to the agc-rs binary (default: ../../target/release/agc-rs)
#   work-dir       — scratch directory for downloads and sub-archives
#                    (default: ./hprc_r2_work)

set -euo pipefail

BATCH_SIZE=24
KEEP_BATCHES="${KEEP_BATCHES:-0}"   # set to 1 to skip cleanup of sub-archives

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
TEST_MODE=0
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=1
    shift
fi

# NOTE: This script requires agc-rs to be installed in your PATH.
# Install from the repository root with:
#   cargo install --path <repo-root>/agc-rs
ARCHIVE="${1:?usage: $0 [--test] <archive.agcrs> [agc-rs-binary] [work-dir]}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
AGC_RS="${2:-agc-rs}"
WORK_DIR="${3:-$(pwd)/hprc_r2_work}"

if ! command -v "$AGC_RS" &>/dev/null && [[ ! -x "$AGC_RS" ]]; then
    echo "ERROR: agc-rs not found in PATH" >&2
    echo "       Install with: cargo install --path <repo-root>/agc-rs" >&2
    exit 1
fi

if (( TEST_MODE )); then
    echo "[INFO] TEST MODE — limiting to first 10 haplotypes"
fi

S3_HTTP_BASE="https://human-pangenomics.s3.us-west-2.amazonaws.com"
CSV_URL="https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv"
CSV_FILE="${WORK_DIR}/assemblies_release2_v1.0.index.csv"

mkdir -p "$WORK_DIR"

# ---------------------------------------------------------------------------
# Short-circuit: if the final archive already exists, just print info and exit
# ---------------------------------------------------------------------------
if [[ -f "$ARCHIVE" ]]; then
    echo "[INFO] Archive already exists: $ARCHIVE"
    "$AGC_RS" info "$ARCHIVE"
    exit 0
fi

# ---------------------------------------------------------------------------
# Download assembly index if not already cached
# ---------------------------------------------------------------------------
if [[ ! -f "$CSV_FILE" ]]; then
    echo "[INDEX] Downloading assembly index ..."
    curl -fsSL "$CSV_URL" -o "$CSV_FILE"
    echo "[INDEX] Saved: $CSV_FILE"
fi

# ---------------------------------------------------------------------------
# Build work list: "sample_name\thttps_url" pairs, sorted by sample_id/hap
# ---------------------------------------------------------------------------
WORK_LIST="${WORK_DIR}/work_list.tsv"
HAP_LIMIT=$(( TEST_MODE ? 10 : 0 ))
python3 - "$CSV_FILE" "$WORK_LIST" "$S3_HTTP_BASE" "$HAP_LIMIT" <<'PYEOF'
import csv, sys

csv_path, out_path, http_base = sys.argv[1], sys.argv[2], sys.argv[3]

rows = []
with open(csv_path, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        name  = row["assembly_name"].strip()
        s3url = row["assembly"].strip()
        sid   = row["sample_id"].strip()

        if "_mat_" in name:
            hap = "mat"
        elif "_pat_" in name:
            hap = "pat"
        elif "_hap1_" in name:
            hap = "hap1"
        elif "_hap2_" in name:
            hap = "hap2"
        else:
            phasing   = row.get("phasing", "").strip()
            haplotype = row["haplotype"].strip()
            if phasing in ("trio", "hic"):
                hap = "mat" if haplotype == "2" else "pat"
            else:
                hap = f"hap{haplotype}"

        s3_prefix = "s3://human-pangenomics/"
        key = s3url[len(s3_prefix):]
        http_url = f"{http_base}/{key}"

        rows.append((sid, hap, f"{sid}_{hap}", http_url))

rows.sort(key=lambda r: (r[0], r[1]))

limit = int(sys.argv[4]) if len(sys.argv) > 4 else 0
if limit:
    rows = rows[:limit]

with open(out_path, "w") as out:
    for _, _, sample_name, http_url in rows:
        out.write(f"{sample_name}\t{http_url}\n")

print(f"[INDEX] {len(rows)} assemblies written to work list", flush=True)
PYEOF

# ---------------------------------------------------------------------------
# Load work list into arrays
# ---------------------------------------------------------------------------
mapfile -t SAMPLE_NAMES < <(cut -f1 "$WORK_LIST")
mapfile -t HTTP_URLS    < <(cut -f2 "$WORK_LIST")
TOTAL=${#SAMPLE_NAMES[@]}

if (( TOTAL == 0 )); then
    echo "ERROR: work list is empty" >&2
    exit 1
fi

echo "[BUILD] $TOTAL haplotypes to process (batch size: $BATCH_SIZE)"

# ---------------------------------------------------------------------------
# Phase 1 — Download the common reference (first haplotype), keep it
# ---------------------------------------------------------------------------
REF_SAMPLE="${SAMPLE_NAMES[0]}"
REF_URL="${HTTP_URLS[0]}"
REF_FA="${WORK_DIR}/$(basename "$REF_URL")"

if [[ ! -f "$REF_FA" ]]; then
    echo "[REF] Downloading reference: $REF_SAMPLE"
    curl -fL --retry 3 --retry-delay 5 -o "$REF_FA" "$REF_URL"
    echo "[REF] Saved: $REF_FA"
else
    echo "[REF] Using cached reference: $REF_FA"
fi

# ---------------------------------------------------------------------------
# Phase 2 — Process remaining haplotypes in batches of BATCH_SIZE
# ---------------------------------------------------------------------------
# Batch 0 covers indices 1 .. BATCH_SIZE   (index 0 = reference, not repeated)
# Batch 1 covers indices BATCH_SIZE+1 .. 2*BATCH_SIZE
# ... and so on.
#
# Each batch sub-archive is seeded with the common reference so that all
# sub-archives share the same splitter set and can be merged cleanly.

BATCH_ARCHIVES=()
BATCH_NUM=0
BATCH_START=1   # index 0 is the reference handled above

while (( BATCH_START < TOTAL )); do
    BATCH_END=$(( BATCH_START + BATCH_SIZE ))
    (( BATCH_END > TOTAL )) && BATCH_END=$TOTAL
    BATCH_COUNT=$(( BATCH_END - BATCH_START ))

    BATCH_LABEL=$(printf "%04d" "$BATCH_NUM")
    BATCH_ARCHIVE="${WORK_DIR}/batch_${BATCH_LABEL}.agcrs"
    BATCH_ARCHIVES+=("$BATCH_ARCHIVE")

    if [[ -f "$BATCH_ARCHIVE" ]]; then
        echo "[BATCH ${BATCH_LABEL}] sub-archive exists — skipping (hap ${BATCH_START}–$((BATCH_END-1)))"
        BATCH_NUM=$(( BATCH_NUM + 1 ))
        BATCH_START=$BATCH_END
        continue
    fi

    echo ""
    echo "[BATCH ${BATCH_LABEL}] === haplotypes $BATCH_START–$((BATCH_END-1)) ($BATCH_COUNT samples) ==="

    # Scratch directory for this batch's FASTAs
    BATCH_DIR="${WORK_DIR}/batch_${BATCH_LABEL}_fa"
    mkdir -p "$BATCH_DIR"

    # Seed sub-archive with the common reference
    echo "[BATCH ${BATCH_LABEL}] Creating sub-archive (reference: $REF_SAMPLE) ..."
    "$AGC_RS" create --output "$BATCH_ARCHIVE" --sample "$REF_SAMPLE" "$REF_FA"

    # Download all FASTAs for this batch
    echo "[BATCH ${BATCH_LABEL}] Downloading $BATCH_COUNT FASTAs ..."
    BATCH_INPUT_ARGS=()
    BATCH_LOCAL_FILES=()
    for (( i = BATCH_START; i < BATCH_END; i++ )); do
        SNAME="${SAMPLE_NAMES[$i]}"
        URL="${HTTP_URLS[$i]}"
        LOCAL_FA="${BATCH_DIR}/$(basename "$URL")"
        echo "[BATCH ${BATCH_LABEL}]   DL $(( i - BATCH_START + 1 ))/${BATCH_COUNT}  $SNAME"
        curl -fL --retry 3 --retry-delay 5 -o "$LOCAL_FA" "$URL"
        BATCH_INPUT_ARGS+=("${SNAME}:${LOCAL_FA}")
        BATCH_LOCAL_FILES+=("$LOCAL_FA")
    done

    # Compress all batch samples with a single shared index
    echo "[BATCH ${BATCH_LABEL}] batch-append $BATCH_COUNT samples (shared index) ..."
    "$AGC_RS" batch-append "$BATCH_ARCHIVE" "${BATCH_INPUT_ARGS[@]}"

    # Remove downloaded FASTAs and scratch dir
    rm -f "${BATCH_LOCAL_FILES[@]}"
    rmdir "$BATCH_DIR" 2>/dev/null || true
    echo "[BATCH ${BATCH_LABEL}] DONE  (downloaded FASTAs removed)"

    BATCH_NUM=$(( BATCH_NUM + 1 ))
    BATCH_START=$BATCH_END
done

echo ""
echo "[BUILD] All ${#BATCH_ARCHIVES[@]} sub-archive(s) complete"

# ---------------------------------------------------------------------------
# Phase 3 — Merge all sub-archives into the final archive
# ---------------------------------------------------------------------------
echo ""
echo "[MERGE] Merging ${#BATCH_ARCHIVES[@]} sub-archive(s) → $ARCHIVE ..."
"$AGC_RS" merge --output "$ARCHIVE" "${BATCH_ARCHIVES[@]}"

echo ""
echo "[DONE] Archive: $ARCHIVE"
"$AGC_RS" info "$ARCHIVE"

# ---------------------------------------------------------------------------
# Phase 4 — Cleanup (sub-archives and cached reference FA)
# ---------------------------------------------------------------------------
if (( ! KEEP_BATCHES )); then
    echo "[CLEAN] Removing sub-archives and cached reference FA ..."
    rm -f "${BATCH_ARCHIVES[@]}"
    rm -f "$REF_FA"
    echo "[CLEAN] Done"
else
    echo "[CLEAN] KEEP_BATCHES=1 — sub-archives and reference FA kept in $WORK_DIR"
fi
