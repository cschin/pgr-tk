#!/usr/bin/env bash
# pgr-server/tests/run_integration_tests.sh
#
# Integration tests for pgr-server.  Starts a server using the small hg002
# chr6 example database, exercises every endpoint, and reports PASS / FAIL.
#
# Usage:
#   bash run_integration_tests.sh [OPTIONS]
#
# OPTIONS
#   --server-url  URL    test against an already-running server (skip start)
#   --server-bin  PATH   pgr-server binary (default: ../../target/release/pgr-server)
#   --pgr-bin     PATH   pgr binary for /bundle tests (default: pgr)
#   --port        N      port to bind when starting server (default: 13737)
#   -h, --help           show this help
#
# The test databases are taken from:
#   ../../examples/hg002/example_output/hg002_chr6_pan.*
#   ../../examples/hg002/hg38.ncbiRefSeq.db

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ── defaults ──────────────────────────────────────────────────────────────────
SERVER_BIN="${REPO_ROOT}/target/release/pgr-server"
PGR_BIN="pgr"
PORT=13737
SERVER_URL=""
SKIP_START=0

# ── option parsing ─────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --server-url) SERVER_URL="$2"; SKIP_START=1; shift 2 ;;
        --server-bin) SERVER_BIN="$2"; shift 2 ;;
        --pgr-bin)    PGR_BIN="$2";    shift 2 ;;
        --port)       PORT="$2";       shift 2 ;;
        -h|--help)
            sed -n '2,/^set /p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$SERVER_URL" ]] && SERVER_URL="http://127.0.0.1:${PORT}"

# The .midx file stores relative paths to .agcrs / .mdbv; the server must be
# started from the same working directory the database was built in.
HG002_DIR="${REPO_ROOT}/examples/hg002"
DB_PREFIX="example_output/hg002_chr6_pan"   # relative to HG002_DIR
GENE_DB="${HG002_DIR}/hg38.ncbiRefSeq.db"

# ── prerequisites check ───────────────────────────────────────────────────────
for f in "${HG002_DIR}/${DB_PREFIX}.agcrs" "${HG002_DIR}/${DB_PREFIX}.mdbi" \
          "${HG002_DIR}/${DB_PREFIX}.mdbv"  "${HG002_DIR}/${DB_PREFIX}.midx"; do
    if [[ ! -f "$f" ]]; then
        echo "SKIP: test database not found: $f"
        echo "      Run examples/hg002/run_all.sh first to build it."
        exit 0
    fi
done

if [[ ! -f "$GENE_DB" ]]; then
    echo "SKIP: gene annotation DB not found: $GENE_DB"
    exit 0
fi

# ── helpers ───────────────────────────────────────────────────────────────────
PASS=0; FAIL=0; SKIP_COUNT=0
SERVER_PID=""

cleanup() {
    if [[ -n "$SERVER_PID" ]] && kill -0 "$SERVER_PID" 2>/dev/null; then
        kill "$SERVER_PID" 2>/dev/null || true
        wait "$SERVER_PID" 2>/dev/null || true   # reap so bash doesn't print "Terminated"
    fi
    rm -f "/tmp/pgr_test_$$.yaml" "/tmp/pgr_test_$$.log"
}
trap cleanup EXIT

check() {
    local label="$1" result="$2"
    if [[ "$result" == "ok" ]]; then
        printf "  PASS  %s\n" "$label"
        (( PASS++ )) || true
    else
        printf "  FAIL  %s — %s\n" "$label" "$result"
        (( FAIL++ )) || true
    fi
}

http_get() {   # http_get PATH → stdout; exit 1 on HTTP error
    python3 -c "
import urllib.request, sys
try:
    with urllib.request.urlopen('${SERVER_URL}' + sys.argv[1], timeout=10) as r:
        sys.stdout.buffer.write(r.read())
except urllib.error.HTTPError as e:
    print(e.read().decode(errors='replace'), file=sys.stderr)
    sys.exit(1)
" "$1"
}

http_post() {  # http_post PATH JSON → stdout; exit 1 on HTTP error
    python3 -c "
import urllib.request, sys, json
path, body = sys.argv[1], sys.argv[2].encode()
req = urllib.request.Request('${SERVER_URL}' + path, data=body,
      headers={'Content-Type': 'application/json'}, method='POST')
try:
    with urllib.request.urlopen(req, timeout=30) as r:
        sys.stdout.buffer.write(r.read())
except urllib.error.HTTPError as e:
    print(e.read().decode(errors='replace'), file=sys.stderr)
    sys.exit(1)
" "$1" "$2"
}

jq_field() {   # jq_field JSON KEY → value string
    python3 -c "import json,sys; d=json.loads(sys.argv[1]); print(d$(echo "$2" | sed "s/\./']['/g;s/^/['/;s/$/']/" | sed "s/\['\([0-9]*\)'\]/[\1]/g"))" "$1" 2>/dev/null || echo ""
}

jq_len() {     # jq_len JSON KEY → length of list at KEY
    python3 -c "import json,sys; d=json.loads(sys.argv[1]); print(len(d['$2']))" "$1" 2>/dev/null || echo "0"
}

# ── start server (unless --server-url given) ───────────────────────────────────
if [[ "$SKIP_START" -eq 0 ]]; then
    if [[ ! -x "$SERVER_BIN" ]]; then
        echo "ERROR: pgr-server binary not found: $SERVER_BIN"
        echo "       Build with: cargo build -p pgr-server --release"
        exit 1
    fi

    CFG_FILE="/tmp/pgr_test_$$.yaml"
    LOG_FILE="/tmp/pgr_test_$$.log"
    : > "$CFG_FILE"   # ensure it exists before writing

    cat > "$CFG_FILE" <<YAML
server:
  host: "127.0.0.1"
  port: ${PORT}
  log_level: "info"
  pgr_binary: "${PGR_BIN}"

databases:
  - name: "default"
    db_prefix: "${DB_PREFIX}"
    gene_db: "${GENE_DB}"
    memory_mode: "moderate"
    ref_sample_hint: "GRCh38"

query_defaults:
  flank: 100000
  max_count: 128
  max_query_count: 128
  max_target_count: 128
  merge_range_tol: 100000
  min_anchor_count: 5
YAML

    # Kill any stale process already listening on the port.
    stale=$(lsof -ti ":${PORT}" 2>/dev/null || true)
    if [[ -n "$stale" ]]; then
        echo "     Killing stale process on port ${PORT}: PID $stale"
        kill "$stale" 2>/dev/null || true
        sleep 1
    fi

    echo "Starting pgr-server (port ${PORT}, cwd=${HG002_DIR})..."
    # Start from HG002_DIR so the relative paths stored in .midx resolve correctly.
    # exec replaces the subshell, so $! == actual pgr-server PID.
    (cd "$HG002_DIR" && exec "$SERVER_BIN" --config "$CFG_FILE") \
        >> "$LOG_FILE" 2>&1 &
    SERVER_PID=$!

    for i in $(seq 1 30); do
        sleep 1
        if curl -sf "${SERVER_URL}/api/v1/health" > /dev/null 2>&1; then
            echo "Server ready after ${i}s  (PID ${SERVER_PID})"
            break
        fi
        if ! kill -0 "$SERVER_PID" 2>/dev/null; then
            echo "ERROR: pgr-server exited unexpectedly. Log:"
            cat "$LOG_FILE"
            exit 1
        fi
        if [[ $i -eq 30 ]]; then
            echo "ERROR: server did not become healthy within 30s. Log:"
            cat "$LOG_FILE"
            exit 1
        fi
    done
    echo
fi

# ══════════════════════════════════════════════════════════════════════════════
echo "=== Running integration tests against ${SERVER_URL} ==="
echo

# ── T01: GET /health ──────────────────────────────────────────────────────────
R=$(http_get /api/v1/health 2>/dev/null) || R='{"status":"error"}'
STATUS=$(python3 -c "import json,sys; print(json.loads(sys.argv[1]).get('status',''))" "$R" 2>/dev/null)
[[ "$STATUS" == "ok" ]] && check "T01 GET /health → {status:ok}" "ok" \
                         || check "T01 GET /health → {status:ok}" "got: $R"

# ── T02: GET /info ────────────────────────────────────────────────────────────
R=$(http_get /api/v1/info 2>/dev/null) || R='{}'
N=$(jq_len "$R" databases)
[[ "$N" -ge 1 ]] && check "T02 GET /info → ≥1 database" "ok" \
                  || check "T02 GET /info → ≥1 database" "databases=$N"

# ── T03: GET /sequences ───────────────────────────────────────────────────────
SEQS=$(http_get /api/v1/sequences 2>/dev/null) || SEQS='[]'
N=$(python3 -c "import json,sys; print(len(json.loads(sys.argv[1])))" "$SEQS" 2>/dev/null || echo 0)
[[ "$N" -ge 2 ]] && check "T03 GET /sequences → ≥2 entries" "ok" \
                  || check "T03 GET /sequences → ≥2 entries" "got $N"

# ── T04: resolve reference ctg (PanSN-aware, done here not in server) ────────
REF_CTG=$(python3 - "$SEQS" "GRCh38" "chr6" <<'PYEOF'
import json, sys
seqs, hint, chrom = json.loads(sys.argv[1]), sys.argv[2], sys.argv[3]
for s in seqs:
    ctg, src = s["ctg"], s.get("src", "")
    if hint.lower() in src.lower() and (ctg.rsplit("#",1)[-1] == chrom or ctg == chrom):
        print(ctg); raise SystemExit(0)
print("NOT_FOUND"); raise SystemExit(1)
PYEOF
) || REF_CTG="NOT_FOUND"
[[ "$REF_CTG" != "NOT_FOUND" ]] && check "T04 resolve chr6 → PanSN ctg ($REF_CTG)" "ok" \
                                 || check "T04 resolve chr6 → PanSN ctg" "not found"

# ── T05: POST /query/region (hits only) ───────────────────────────────────────
if [[ "$REF_CTG" != "NOT_FOUND" ]]; then
    REQ=$(python3 -c "
import json
print(json.dumps({'chrom':'${REF_CTG}','start':31218748,'end':31322092,
  'include_sequences':False,'min_anchor_count':5,'merge_range_tol':100000,
  'max_count':64,'max_query_count':64,'max_target_count':64}))
")
    R=$(http_post /api/v1/query/region "$REQ" 2>/dev/null) || R='{"hits":[]}'
    N=$(jq_len "$R" hits)
    [[ "$N" -ge 2 ]] && check "T05 POST /query/region HLA-C region → ≥2 hits" "ok" \
                      || check "T05 POST /query/region HLA-C region → ≥2 hits" "hits=$N resp=$R"
else
    printf "  SKIP  T05 POST /query/region (ctg not resolved)\n"; (( SKIP_COUNT++ )) || true
fi

# ── T06: POST /query/region with sequences ────────────────────────────────────
if [[ "$REF_CTG" != "NOT_FOUND" ]]; then
    REQ=$(python3 -c "
import json
print(json.dumps({'chrom':'${REF_CTG}','start':31268748,'end':31272092,
  'include_sequences':True,'min_anchor_count':3,'merge_range_tol':50000,
  'max_count':64,'max_query_count':64,'max_target_count':64}))
")
    R=$(http_post /api/v1/query/region "$REQ" 2>/dev/null) || R='{}'
    FASTA=$(python3 -c "import json,sys; print(json.loads(sys.argv[1]).get('sequences_fasta') or '')" "$R" 2>/dev/null)
    [[ "${#FASTA}" -gt 10 ]] && check "T06 POST /query/region include_sequences → non-empty FASTA" "ok" \
                              || check "T06 POST /query/region include_sequences → non-empty FASTA" "fasta len=${#FASTA}"
else
    printf "  SKIP  T06 (ctg not resolved)\n"; (( SKIP_COUNT++ )) || true
fi

# ── T07: POST /bundle ─────────────────────────────────────────────────────────
if [[ -n "${FASTA:-}" && "${#FASTA}" -gt 10 ]]; then
    BUNDLE_REQ=$(python3 -c "
import json, sys
fasta = sys.argv[1]
print(json.dumps({'sequences_fasta': fasta, 'generate_html': True}))
" "$FASTA")
    R=$(http_post /api/v1/bundle "$BUNDLE_REQ" 2>/dev/null) || R='{}'
    BED=$(python3 -c "import json,sys; print(json.loads(sys.argv[1]).get('bed',''))" "$R" 2>/dev/null)
    HTML=$(python3 -c "import json,sys; print(json.loads(sys.argv[1]).get('html') or '')" "$R" 2>/dev/null)
    if [[ "${#BED}" -gt 5 && "${#HTML}" -gt 100 ]]; then
        check "T07 POST /bundle → BED + HTML" "ok"
    else
        check "T07 POST /bundle → BED + HTML" "bed=${#BED}b html=${#HTML}b resp=$R"
    fi
else
    printf "  SKIP  T07 POST /bundle (no FASTA from T06)\n"; (( SKIP_COUNT++ )) || true
fi

# ── T08: GET /api/v1/openapi.json ────────────────────────────────────────────
R=$(http_get /api/v1/openapi.json 2>/dev/null) || R='{}'
TITLE=$(python3 -c "import json,sys; print(json.loads(sys.argv[1]).get('info',{}).get('title',''))" "$R" 2>/dev/null)
[[ "$TITLE" == "pgr-server" ]] && check "T08 GET /openapi.json → title=pgr-server" "ok" \
                                || check "T08 GET /openapi.json → title=pgr-server" "title='$TITLE'"

# ── summary ───────────────────────────────────────────────────────────────────
echo
TOTAL=$(( PASS + FAIL + SKIP_COUNT ))
echo "Results: ${PASS} passed, ${FAIL} failed, ${SKIP_COUNT} skipped  (${TOTAL} total)"
[[ "$FAIL" -eq 0 ]]
