#!/usr/bin/env bash
# 01_agcrs_basics.sh — build an .agcrs archive from the E. coli test set and
# exercise the core agc-rs commands: create, append, info, list, get, subrange.
#
# Run from the examples/ecoli/ directory or anywhere — the script resolves
# all paths relative to itself.
#
# Output is written to examples/ecoli/example_output/
#
# Usage:
#   bash examples/ecoli/01_agcrs_basics.sh

set -euo pipefail
cd "$(dirname "$0")"

# NOTE: This script requires agc-rs to be installed in your PATH.
# Install from the repository root with:
#   cargo install --path ../../agc-rs
AGC_RS="${AGC_RS:-agc-rs}"
TESTDATA="../../test_data/ecoli"
OUT="example_output"
ARCHIVE="$OUT/ecoli_demo.agcrs"

if ! command -v "$AGC_RS" &>/dev/null && [[ ! -x "$AGC_RS" ]]; then
    echo "ERROR: agc-rs not found in PATH" >&2
    echo "       Install with: cargo install --path <repo-root>/agc-rs" >&2
    exit 1
fi

mkdir -p "$OUT"

# ---- 1. Create archive with MG1655 as the reference sample -----------------
rm -f "$ARCHIVE" "$ARCHIVE-wal" "$ARCHIVE-shm"
echo
echo "=== [1] Creating archive: $ARCHIVE (reference: MG1655) ==="
"$AGC_RS" create \
    --output "$ARCHIVE" \
    --sample MG1655 \
    "$TESTDATA/ecoli_k12_mg1655.fna.gz"

# ---- 2. Append W3110 -------------------------------------------------------
echo
echo "=== [2] Appending W3110 ==="
"$AGC_RS" append "$ARCHIVE" \
    --sample W3110 \
    "$TESTDATA/ecoli_k12_w3110.fna.gz"

# ---- 3. Append O157:H7 Sakai (chromosome + 2 plasmids) --------------------
echo
echo "=== [3] Appending Sakai ==="
"$AGC_RS" append "$ARCHIVE" \
    --sample Sakai \
    "$TESTDATA/ecoli_o157h7_sakai.fna.gz"

# ---- 4. Archive statistics -------------------------------------------------
echo
echo "=== [4] Archive info ==="
"$AGC_RS" info "$ARCHIVE"

# ---- 5. List samples -------------------------------------------------------
echo
echo "=== [5] Samples ==="
"$AGC_RS" list "$ARCHIVE"

# ---- 6. List contigs for each sample ---------------------------------------
echo
echo "=== [6] Contigs: MG1655 ==="
"$AGC_RS" list --contigs MG1655 "$ARCHIVE"

echo
echo "=== Contigs: W3110 ==="
"$AGC_RS" list --contigs W3110 "$ARCHIVE"

echo
echo "=== Contigs: Sakai (chromosome + 2 plasmids) ==="
"$AGC_RS" list --contigs Sakai "$ARCHIVE"

# ---- 7. Get a full contig (first 200 bp shown) -----------------------------
echo
echo "=== [7] Get W3110/NC_007779.1 (first 200 bp) ==="
"$AGC_RS" get "$ARCHIVE" "W3110/NC_007779.1" \
    | awk 'NR==1{print} NR==2{print substr($0,1,200)"..."; exit}'

# ---- 8. Subrange queries with round-trip verification ----------------------
check_range() {
    local sample="$1" contig="$2" start="$3" end="$4" fasta="$5"
    local len=$(( end - start ))
    echo
    echo "=== Subrange: $sample/$contig:$start-$end ==="
    local got ref
    got=$("$AGC_RS" get "$ARCHIVE" "$sample/$contig:$start-$end" | tail -n +2 | tr -d '\n')
    ref=$(gunzip -c "$fasta" \
        | awk -v id="$contig" '/^>/{found=($0 ~ id); next} found{printf "%s",$0}' \
        | cut -c$(( start + 1 ))-$end)
    if [ "$got" = "$ref" ] && [ "${#got}" -eq "$len" ]; then
        echo "  PASS: $len bases match"
    else
        echo "  FAIL: mismatch (got ${#got} bases, expected $len)"
        exit 1
    fi
}

echo
echo "=== [8] Subrange queries ==="
check_range MG1655  NC_000913.3  1000000  1000060  "$TESTDATA/ecoli_k12_mg1655.fna.gz"
check_range MG1655  NC_000913.3  2300000  2300080  "$TESTDATA/ecoli_k12_mg1655.fna.gz"
check_range W3110   NC_007779.1  4500000  4500100  "$TESTDATA/ecoli_k12_w3110.fna.gz"
check_range Sakai   NC_002695.2   500000   500060  "$TESTDATA/ecoli_o157h7_sakai.fna.gz"
check_range Sakai   NC_002128.1    10000    10060  "$TESTDATA/ecoli_o157h7_sakai.fna.gz"

# ---- 9. Full round-trip check ----------------------------------------------
echo
echo "=== [9] Full round-trip: MG1655/NC_000913.3 ==="
ORIG=$(gunzip -c "$TESTDATA/ecoli_k12_mg1655.fna.gz" | awk 'NR>1{printf "%s",$0} END{print ""}')
GOT=$("$AGC_RS" get "$ARCHIVE" "MG1655/NC_000913.3" | tail -n +2 | tr -d '\n')
if [ "$ORIG" = "$GOT" ]; then
    echo "PASS: $(echo -n "$GOT" | wc -c | tr -d ' ') bases match"
else
    echo "FAIL: sequence mismatch"; exit 1
fi

# ---- 10. Direct SQLite query (no binary needed) ----------------------------
echo
echo "=== [10] Direct SQLite query ==="
sqlite3 "$ARCHIVE" \
    "SELECT sm.name AS sample, c.name AS contig, c.length
     FROM sample sm JOIN contig c ON c.sample_id = sm.id
     ORDER BY sm.name, c.name;"

echo
echo "Archive size: $(du -sh "$ARCHIVE" | cut -f1)"
echo "Done. Archive: $ARCHIVE"
