#!/usr/bin/env bash
# Demo: build an .agcrs archive from the E. coli test set and run queries.
#
# Usage: bash scripts/ecoli_agcrs_demo.sh
#
# The script assumes it is run from the pgr-tk workspace root.

set -euo pipefail

ARCHIVE=ecoli_demo.agcrs
BIN=target/release/agc-rs
TESTDATA=test_data/ecoli

# ---- 0. Build the binary ---------------------------------------------------
echo "=== Building agc-rs (release) ==="
cargo build --release -p agc-rs 2>&1 | tail -3

# ---- 1. Create archive with MG1655 as the reference sample -----------------
rm -f "$ARCHIVE" "$ARCHIVE-wal" "$ARCHIVE-shm"
echo
echo "=== Creating archive: $ARCHIVE ==="
"$BIN" create \
    --output "$ARCHIVE" \
    --sample MG1655 \
    "$TESTDATA/ecoli_k12_mg1655.fna.gz"

# ---- 2. Append W3110 -------------------------------------------------------
echo
echo "=== Appending W3110 ==="
"$BIN" append "$ARCHIVE" \
    --sample W3110 \
    "$TESTDATA/ecoli_k12_w3110.fna.gz"

# ---- 3. Append O157:H7 Sakai (chromosome + 2 plasmids) --------------------
echo
echo "=== Appending Sakai ==="
"$BIN" append "$ARCHIVE" \
    --sample Sakai \
    "$TESTDATA/ecoli_o157h7_sakai.fna.gz"

# ---- 4. Archive statistics -------------------------------------------------
echo
echo "=== Info ==="
"$BIN" info "$ARCHIVE"

# ---- 5. List samples -------------------------------------------------------
echo
echo "=== Samples ==="
"$BIN" list "$ARCHIVE"

# ---- 6. List contigs for each sample ---------------------------------------
echo
echo "=== Contigs: MG1655 ==="
"$BIN" list --contigs MG1655 "$ARCHIVE"

echo
echo "=== Contigs: Sakai (chromosome + 2 plasmids) ==="
"$BIN" list --contigs Sakai "$ARCHIVE"

# ---- 7. Get a full contig --------------------------------------------------
echo
echo "=== Full contig: W3110/NC_007779.1 (first 200 bp shown) ==="
{ "$BIN" get "$ARCHIVE" "W3110/NC_007779.1" 2>/dev/null || true; } \
    | awk 'NR==1{print} NR==2{print substr($0,1,200)"..."; exit}'

# ---- 8. Subrange query -----------------------------------------------------
echo
echo "=== Subrange: MG1655/NC_000913.3:1000000-1000060 ==="
"$BIN" get "$ARCHIVE" "MG1655/NC_000913.3:1000000-1000060"

echo
echo "=== Subrange: Sakai/NC_002695.2:500000-500060 ==="
"$BIN" get "$ARCHIVE" "Sakai/NC_002695.2:500000-500060"

# ---- 9. Verify round-trip with gunzip + diff ------------------------------
echo
echo "=== Round-trip check: MG1655/NC_000913.3 ==="
ORIG=$(gunzip -c "$TESTDATA/ecoli_k12_mg1655.fna.gz" | awk 'NR>1{printf "%s", $0} END{print ""}')
GOT=$("$BIN" get "$ARCHIVE" "MG1655/NC_000913.3" | tail -n +2 | tr -d '\n')
if [ "$ORIG" = "$GOT" ]; then
    echo "PASS: sequences match ($(echo -n "$GOT" | wc -c) bases)"
else
    echo "FAIL: sequence mismatch"
    exit 1
fi

# ---- 10. Verify SQLite is readable without agc-rs -------------------------
echo
echo "=== Raw SQLite query (no agc-rs needed) ==="
sqlite3 "$ARCHIVE" \
    "SELECT sm.name AS sample, c.name AS contig, c.length
     FROM sample sm JOIN contig c ON c.sample_id = sm.id
     ORDER BY sm.name, c.name;"

echo
echo "Archive size: $(stat -f%z "$ARCHIVE" | awk '{printf "%.2f MB", $1/1048576}')  (uncompressed input: ~$(gunzip -c "$TESTDATA"/*.fna.gz | wc -c | awk '{printf "%.0f MB", $1/1048576}'))"
echo
echo "Done. Archive written to: $ARCHIVE"
