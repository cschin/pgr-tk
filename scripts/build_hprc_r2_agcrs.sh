#!/usr/bin/env bash
# build_hprc_r2_agcrs.sh — incrementally build an agc-rs pangenome archive
# from the HPRC release-2 assembly table.
#
# For each row the script:
#   1. Derives a sample name:  <sample_id>_mat  or  <sample_id>_pat
#   2. Skips samples already present in the archive (resume support)
#   3. Downloads the .fa.gz via HTTPS (public bucket, no credentials needed)
#   4. Creates or appends to the archive
#   5. Deletes the local .fa.gz immediately after a successful add
#
# Usage:
#   bash scripts/build_hprc_r2_agcrs.sh [--test] <archive.agcrs> [agc-rs-bin] [work-dir]
#
# Options:
#   --test         Process only the first 10 haplotypes (quick smoke test)
#
# Arguments:
#   archive.agcrs  — output archive path (created on first run, appended thereafter)
#   agc-rs-bin     — path to the agc-rs binary (default: target/release/agc-rs)
#   work-dir       — scratch directory for downloads (default: ./hprc_r2_downloads)

set -euo pipefail

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
TEST_MODE=0
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=1
    shift
fi

ARCHIVE="${1:?usage: $0 [--test] <archive.agcrs> [agc-rs-binary] [work-dir]}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
AGC_RS="${2:-${SCRIPT_DIR}/../target/release/agc-rs}"
WORK_DIR="${3:-$(pwd)/hprc_r2_downloads}"

if [[ ! -x "$AGC_RS" ]]; then
    echo "ERROR: agc-rs binary not found at $AGC_RS" >&2
    echo "       Build with: cargo build -p agc-rs --release" >&2
    exit 1
fi

if (( TEST_MODE )); then
    echo "[TEST MODE] limiting to first 10 haplotypes"
fi

# S3 → HTTPS base for the public human-pangenomics bucket
S3_HTTP_BASE="https://human-pangenomics.s3.us-west-2.amazonaws.com"

CSV_URL="https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv"
CSV_FILE="${WORK_DIR}/assemblies_release2_v1.0.index.csv"

mkdir -p "$WORK_DIR"

# ---------------------------------------------------------------------------
# Download the assembly index if not already cached
# ---------------------------------------------------------------------------
if [[ ! -f "$CSV_FILE" ]]; then
    echo "[INDEX] Downloading assembly index ..."
    curl -fsSL "$CSV_URL" -o "$CSV_FILE"
    echo "[INDEX] Saved: $CSV_FILE"
fi

# ---------------------------------------------------------------------------
# Build a work list: "sample_name\thttps_url" pairs.
# Sort by sample_id so pat/mat pairs are adjacent (better delta compression).
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
        s3url = row["assembly"].strip()          # s3://human-pangenomics/...
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
            # No mat/pat/hap1/hap2 in assembly name — fall back to haplotype
            # column (HPRC: 1=pat, 2=mat) but keep the generic label if phasing
            # is unknown (neither trio nor hic).
            phasing = row.get("phasing", "").strip()
            haplotype = row["haplotype"].strip()
            if phasing in ("trio", "hic"):
                hap = "mat" if haplotype == "2" else "pat"
            else:
                hap = f"hap{haplotype}"

        # Convert s3://bucket/key  →  https://bucket.s3.region.amazonaws.com/key
        s3_prefix = "s3://human-pangenomics/"
        key = s3url[len(s3_prefix):]
        http_url = f"{http_base}/{key}"

        rows.append((sid, hap, f"{sid}_{hap}", http_url))

# Sort by sample_id then hap so pat/mat pairs are adjacent
rows.sort(key=lambda r: (r[0], r[1]))

import sys as _sys
limit = int(_sys.argv[4]) if len(_sys.argv) > 4 else 0
if limit:
    rows = rows[:limit]

with open(out_path, "w") as out:
    for _, _, sample_name, http_url in rows:
        out.write(f"{sample_name}\t{http_url}\n")

print(f"[INDEX] {len(rows)} assemblies written to work list", flush=True)
PYEOF

# ---------------------------------------------------------------------------
# Collect already-present samples (for resume support)
# ---------------------------------------------------------------------------
EXISTING_SAMPLES_FILE="${WORK_DIR}/existing_samples.txt"
if [[ -f "$ARCHIVE" ]]; then
    "$AGC_RS" list "$ARCHIVE" > "$EXISTING_SAMPLES_FILE" 2>/dev/null || true
    echo "[RESUME] $(wc -l < "$EXISTING_SAMPLES_FILE") samples already in archive"
else
    > "$EXISTING_SAMPLES_FILE"
    echo "[RESUME] Archive does not exist yet; will create on first sample"
fi

# ---------------------------------------------------------------------------
# Process each assembly
# ---------------------------------------------------------------------------
TOTAL=$(wc -l < "$WORK_LIST")
IDX=0

while IFS=$'\t' read -r SAMPLE_NAME HTTP_URL; do
    IDX=$((IDX + 1))
    BASENAME="$(basename "$HTTP_URL")"
    LOCAL_FA="${WORK_DIR}/${BASENAME}"

    # Skip if already in archive
    if grep -qx "$SAMPLE_NAME" "$EXISTING_SAMPLES_FILE" 2>/dev/null; then
        echo "[${IDX}/${TOTAL}] SKIP  $SAMPLE_NAME — already in archive"
        continue
    fi

    echo "[${IDX}/${TOTAL}] ===== $SAMPLE_NAME ====="

    # ---- Download -----------------------------------------------------------
    echo "[${IDX}/${TOTAL}] DL    $HTTP_URL"
    curl -fL --retry 3 --retry-delay 5 -o "$LOCAL_FA" "$HTTP_URL"

    # ---- Create or append ---------------------------------------------------
    if [[ ! -f "$ARCHIVE" ]]; then
        echo "[${IDX}/${TOTAL}] CREATE  $ARCHIVE  (first sample: $SAMPLE_NAME)"
        "$AGC_RS" create "$LOCAL_FA" --output "$ARCHIVE" --sample "$SAMPLE_NAME"
    else
        echo "[${IDX}/${TOTAL}] APPEND  $SAMPLE_NAME → $ARCHIVE"
        "$AGC_RS" append "$ARCHIVE" --sample "$SAMPLE_NAME" "$LOCAL_FA"
    fi

    # ---- Cleanup ------------------------------------------------------------
    rm -f "$LOCAL_FA"
    echo "[${IDX}/${TOTAL}] DONE  $SAMPLE_NAME  (local copy removed)"

    # Update cache so resume within the same run works correctly
    echo "$SAMPLE_NAME" >> "$EXISTING_SAMPLES_FILE"

done < "$WORK_LIST"

echo ""
echo "All done.  Archive: $ARCHIVE"
"$AGC_RS" info "$ARCHIVE"
