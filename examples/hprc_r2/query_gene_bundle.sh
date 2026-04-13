#!/usr/bin/env bash
# query_gene_bundle.sh — query any gene locus from a pangenome index,
# run principal bundle decomposition on all haplotype hits, and generate
# an interactive HTML visualisation.
#
# Gene coordinates are derived from either:
#   • a SQLite gene database (.db) built by fetch_refseq_gtf_db.sh  [recommended]
#   • a raw GTF annotation file (plain or .gz)
#
# The reference sequence for the query is fetched directly from the agcrs
# archive via pgr query fetch — no UCSC API call needed.
#
# ── Modes ────────────────────────────────────────────────────────────────────
#
# local  (default, -m local)
#   All work is done by spawning pgr sub-commands inline.  The shimmer index
#   is loaded fresh for every query, which is slow on large databases.
#   Steps: [1] gene lookup → [2] seq-list cache → [3] fetch ref seq →
#          [4] pgr query seqs → [5] pgr bundle decomp → [6] pgr bundle svg
#
# server (-m server)
#   Queries are routed to a pgr-server HTTP process.  On first invocation the
#   script writes a server config and starts the server in the background; on
#   subsequent calls the existing process is reused (index already in RAM).
#   Steps: [S0] check / start server → [S1] POST /query/gene (or /query/region)
#          → [S2] POST /bundle
#   Requires: pgr-server binary in PATH (or SERVER_BIN env var).
#
# ── Requirements ─────────────────────────────────────────────────────────────
#   local mode : pgr in PATH
#   server mode: pgr-server in PATH (or SERVER_BIN=<path>), python3
#   both modes : python3, sqlite3 (for .db annotation)
#
# ── Usage ────────────────────────────────────────────────────────────────────
#   bash query_gene_bundle.sh [OPTIONS] <gene_name> <annotation> [db_prefix] [flank_bp] [out_dir]
#
#   OPTIONS
#     -m, --mode  local|server   execution mode (default: local)
#     -h, --help                 show this help and exit
#
#   POSITIONAL
#     gene_name  — gene symbol or gene_id  (e.g. C9orf72, BRCA1, HLA-A)
#     annotation — SQLite .db from fetch_refseq_gtf_db.sh  OR  .gtf / .gtf.gz
#     db_prefix  — prefix for .agcrs/.mdbi/.mdbv/.midx files   (default: hprc_r2)
#     flank_bp   — flanking bases added on each side            (default: 100000)
#     out_dir    — output directory                             (default: <gene>_bundle)
#
# ── Environment overrides ────────────────────────────────────────────────────
#   PGR=<path>             path to pgr binary         (local mode)
#   REF_SAMPLE_HINT=<s>    substring identifying the reference sample (GRCh38)
#   SERVER_URL=<url>       base URL of the server     (default: http://127.0.0.1:3000)
#   SERVER_PORT=<n>        port to bind when starting (default: 3000)
#   SERVER_BIN=<path>      pgr-server binary          (default: pgr-server)
#   SERVER_MEMORY_MODE=<m> moderate | high            (default: moderate)

set -euo pipefail

# ── Option parsing ─────────────────────────────────────────────────────────────
MODE="local"   # default; overridden by -m / --mode

usage() {
    cat >&2 <<EOF
Usage: $0 [OPTIONS] <gene_name> <annotation> [db_prefix] [flank_bp] [out_dir]

OPTIONS
  -m, --mode  local|server   execution mode (default: local)
  -h, --help                 show this help and exit

POSITIONAL
  gene_name   gene symbol or gene_id (e.g. HLA-A, BRCA1, C9orf72)
  annotation  SQLite .db from fetch_refseq_gtf_db.sh  OR  .gtf / .gtf.gz
  db_prefix   prefix for .agcrs/.mdbi/.mdbv/.midx     (default: hprc_r2)
  flank_bp    flanking bases on each side              (default: 100000)
  out_dir     output directory                         (default: <gene>_bundle)

ENV (fine-tuning)
  PGR=<path>              pgr binary              (local mode)
  REF_SAMPLE_HINT=<s>     reference sample hint   (default: GRCh38)
  SERVER_URL=<url>        server base URL         (default: http://127.0.0.1:3000)
  SERVER_PORT=<n>         port to bind on start   (default: 3000)
  SERVER_BIN=<path>       pgr-server binary       (default: pgr-server)
  SERVER_MEMORY_MODE=<m>  moderate|high           (default: moderate)
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mode)
            MODE="${2:?--mode requires an argument: local|server}"
            shift 2
            ;;
        --mode=*)
            MODE="${1#--mode=}"
            shift
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "ERROR: unknown option: $1" >&2
            usage
            ;;
        *)
            break   # first positional argument
            ;;
    esac
done

# ── Positional arguments ───────────────────────────────────────────────────────
GENE_NAME="${1:?$(usage)}"
ANNOTATION="${2:?$(usage)}"
DB_PREFIX="${3:-hprc_r2}"
FLANK="${4:-100000}"
OUT_DIR="${5:-${GENE_NAME}_bundle}"

# ── Environment overrides (fine-tuning) ───────────────────────────────────────
PGR="${PGR:-pgr}"
REF_SAMPLE_HINT="${REF_SAMPLE_HINT:-GRCh38}"
SERVER_URL="${SERVER_URL:-http://127.0.0.1:3000}"
SERVER_PORT="${SERVER_PORT:-3000}"
SERVER_BIN="${SERVER_BIN:-pgr-server}"
SERVER_MEMORY_MODE="${SERVER_MEMORY_MODE:-moderate}"

# ── Validate inputs ───────────────────────────────────────────────────────────
[[ -f "$ANNOTATION" ]] || { echo "ERROR: annotation file not found: $ANNOTATION" >&2; exit 1; }

if [[ "$MODE" == "local" ]]; then
    if ! command -v "$PGR" &>/dev/null && [[ ! -x "$PGR" ]]; then
        echo "ERROR: pgr not found in PATH" >&2
        echo "       Install with: cargo install --path <repo-root>/pgr-bin" >&2
        exit 1
    fi
elif [[ "$MODE" == "server" ]]; then
    : # SERVER_BIN checked lazily in check_or_start_server
else
    echo "ERROR: MODE must be 'local' or 'server' (got: $MODE)" >&2
    exit 1
fi

# ── Setup ─────────────────────────────────────────────────────────────────────
mkdir -p "$OUT_DIR"

PERF_LOG="$OUT_DIR/perf.log"
: > "$PERF_LOG"

GENE_COORDS="$OUT_DIR/gene_coords.tsv"
QUERY_FA="$OUT_DIR/${GENE_NAME}_query.fa"
HIT_PREFIX="$OUT_DIR/${GENE_NAME}_hits"
HIT_FA="${HIT_PREFIX}.000.fa"
BUNDLE_PREFIX="$OUT_DIR/${GENE_NAME}_bundle"

# ── Performance logging ───────────────────────────────────────────────────────
# run_timed LABEL COMMAND [ARGS...]  — wraps a shell command with GNU time.
run_timed() {
    local label="$1"; shift
    local t0=$SECONDS
    printf '\n### %s\n' "$label" >> "$PERF_LOG"
    /usr/bin/time -v -a -o "$PERF_LOG" -- "$@"
    local rc=$?
    printf 'Wall time: %ds\n' $(( SECONDS - t0 )) >> "$PERF_LOG"
    return $rc
}

# run_timed_simple LABEL — records only wall time (for server HTTP calls).
run_timed_simple() {
    local label="$1"
    printf '\n### %s\n' "$label" >> "$PERF_LOG"
    printf 'Wall time: %ds\n' 0 >> "$PERF_LOG"   # placeholder; caller overwrites
}

# log_wall LABEL SECONDS — append a wall-time-only perf entry (server mode).
log_wall() {
    local label="$1" elapsed="$2"
    printf '\n### %s\n' "$label" >> "$PERF_LOG"
    printf 'Wall time: %ds\n' "$elapsed" >> "$PERF_LOG"
}

# ── Performance summary (shared by both modes) ────────────────────────────────
print_perf_summary() {
python3 - "$PERF_LOG" <<'PYEOF'
import sys, re

log = open(sys.argv[1]).read()
sections = re.split(r'^### ', log, flags=re.MULTILINE)

hdr = f"{'Step':<36} {'Wall':>6} {'User':>7} {'Sys':>7} {'MaxRSS(MB)':>12} {'MajFlt':>8} {'MinFlt':>10}"
print(hdr)
print('-' * len(hdr))

for sec in sections[1:]:
    lines = sec.splitlines()
    label = lines[0].strip()

    def find(pat):
        for l in lines:
            m = re.search(pat, l)
            if m: return m.group(1)
        return '?'

    wall_s = find(r'Wall time: (\d+)s')
    user_s = find(r'User time \(seconds\): ([\d.]+)')
    sys_s  = find(r'System time \(seconds\): ([\d.]+)')
    rss_kb = find(r'Maximum resident set size \(kbytes\): (\d+)')
    maj    = find(r'Major \(requiring I/O\) page faults: (\d+)')
    minor  = find(r'Minor \(reclaiming a frame\) page faults: (\d+)')

    wall_str = f"{wall_s}s"               if wall_s != '?' else '?'
    user_str = f"{float(user_s):.1f}s"   if user_s != '?' else '?'
    sys_str  = f"{float(sys_s):.1f}s"    if sys_s  != '?' else '?'
    rss_str  = str(int(rss_kb) // 1024)  if rss_kb != '?' else '?'

    print(f"{label[:36]:<36} {wall_str:>6} {user_str:>7} {sys_str:>7} {rss_str:>12} {maj:>8} {minor:>10}")
PYEOF
}

# ─────────────────────────────────────────────────────────────────────────────
# SERVER MODE helpers
# ─────────────────────────────────────────────────────────────────────────────

# generate_server_config OUT_FILE ANNOTATION_DB_OR_EMPTY
generate_server_config() {
    local cfg="$1"
    local gene_db_line=""
    if [[ "$ANNOTATION" == *.db ]]; then
        gene_db_line="    gene_db: \"${ANNOTATION}\""
    fi

    cat > "$cfg" <<YAML
# Auto-generated by query_gene_bundle.sh — do not edit while server is running.
server:
  host: "127.0.0.1"
  port: ${SERVER_PORT}
  log_level: "info"
  pgr_binary: "${PGR}"

databases:
  - name: "default"
    db_prefix: "${DB_PREFIX}"
${gene_db_line}
    memory_mode: "${SERVER_MEMORY_MODE}"
    ref_sample_hint: "${REF_SAMPLE_HINT}"

query_defaults:
  flank: ${FLANK}
  max_count: 128
  max_query_count: 128
  max_target_count: 128
  merge_range_tol: 100000
  min_anchor_count: 10
YAML
}

# check_server_health — returns 0 if the server responds to /api/v1/health
check_server_health() {
    python3 -c "
import urllib.request, sys
try:
    urllib.request.urlopen('${SERVER_URL}/api/v1/health', timeout=3)
    sys.exit(0)
except Exception:
    sys.exit(1)
" 2>/dev/null
}

# check_or_start_server — starts pgr-server in background if not already running.
# Config, PID, and log are written beside the DB prefix (persistent across gene calls).
check_or_start_server() {
    local cfg="${DB_PREFIX}.pgr-server.yaml"
    local pid_file="${DB_PREFIX}.pgr-server.pid"
    local log_file="${DB_PREFIX}.pgr-server.log"

    if check_server_health; then
        echo "[S0] pgr-server already running at ${SERVER_URL}"
        return 0
    fi

    # Ensure the binary is available
    if ! command -v "$SERVER_BIN" &>/dev/null && [[ ! -x "$SERVER_BIN" ]]; then
        echo "ERROR: pgr-server binary not found: $SERVER_BIN" >&2
        echo "       Build with: cargo build -p pgr-server --release" >&2
        echo "       Or set SERVER_BIN=<path>" >&2
        exit 1
    fi

    echo "[S0] pgr-server not running — generating config and starting..."
    generate_server_config "$cfg"
    echo "     Config : $cfg"
    echo "     Log    : $log_file"

    # Start server in background; suppress stdin so it doesn't block
    "$SERVER_BIN" --config "$cfg" >> "$log_file" 2>&1 < /dev/null &
    local pid=$!
    echo "$pid" > "$pid_file"
    echo "     PID    : $pid  (saved to $pid_file)"
    echo "     Waiting for index to load (memory_mode=${SERVER_MEMORY_MODE})..."

    local i
    for i in $(seq 1 180); do   # up to 15 min (5 s × 180)
        sleep 5
        if check_server_health; then
            echo "     Ready after $((i * 5))s"
            return 0
        fi
        # Check if process is still alive
        if ! kill -0 "$pid" 2>/dev/null; then
            echo "ERROR: pgr-server (pid $pid) exited unexpectedly." >&2
            echo "       Check $log_file for details." >&2
            exit 1
        fi
        # Progress indicator every 30 s
        if (( i % 6 == 0 )); then
            echo "     ... still loading (${i}*5s elapsed)"
        fi
    done

    echo "ERROR: pgr-server did not become healthy within 15 minutes." >&2
    echo "       Check $log_file for details." >&2
    exit 1
}

# server_post ENDPOINT REQUEST_JSON_FILE — HTTP POST, print response body to stdout.
server_post() {
    local endpoint="$1"
    local req_file="$2"
    python3 -c "
import urllib.request, sys

req_file = sys.argv[1]
url = '${SERVER_URL}${endpoint}'

with open(req_file, 'rb') as f:
    data = f.read()

req = urllib.request.Request(url, data=data,
      headers={'Content-Type': 'application/json'}, method='POST')
try:
    with urllib.request.urlopen(req) as r:
        sys.stdout.buffer.write(r.read())
except urllib.error.HTTPError as e:
    body = e.read().decode(errors='replace')
    print(f'HTTP {e.code} from {url}: {body}', file=sys.stderr)
    sys.exit(1)
" "$req_file"
}

# ─────────────────────────────────────────────────────────────────────────────
# ── GENE LOOKUP helper (shared between local .db path and server GTF path) ───
# ─────────────────────────────────────────────────────────────────────────────

lookup_gene_sqlite() {
    echo "=== [1] SQLite lookup for gene: ${GENE_NAME} ==="
    local table_check
    table_check=$(sqlite3 "$ANNOTATION" \
        "SELECT name FROM sqlite_master WHERE type='table' AND name='genes';" 2>/dev/null || true)
    if [[ -z "$table_check" ]]; then
        echo "ERROR: $ANNOTATION has no 'genes' table — run fetch_refseq_gtf_db.sh first" >&2
        exit 1
    fi
    local row
    row=$(sqlite3 "$ANNOTATION" \
        "SELECT chrom, start, end, strand FROM genes WHERE gene_name='${GENE_NAME}'
         ORDER BY CASE WHEN chrom NOT GLOB '*_*' THEN 0 ELSE 1 END, end-start DESC LIMIT 1;")
    if [[ -z "$row" ]]; then
        row=$(sqlite3 "$ANNOTATION" \
            "SELECT chrom, start, end, strand FROM genes WHERE gene_id='${GENE_NAME}'
             ORDER BY CASE WHEN chrom NOT GLOB '*_*' THEN 0 ELSE 1 END, end-start DESC LIMIT 1;")
    fi
    if [[ -z "$row" ]]; then
        echo "ERROR: gene '${GENE_NAME}' not found in $ANNOTATION" >&2
        local similar
        similar=$(sqlite3 "$ANNOTATION" \
            "SELECT gene_name FROM genes WHERE gene_name LIKE '${GENE_NAME}%' LIMIT 10;" \
            2>/dev/null || true)
        if [[ -n "$similar" ]]; then
            echo "       Similar names in db:" >&2
            echo "$similar" | sed 's/^/         /' >&2
        fi
        exit 1
    fi
    CHROM=$(echo "$row"      | cut -d'|' -f1)
    GENE_START=$(echo "$row" | cut -d'|' -f2)
    GENE_END=$(echo "$row"   | cut -d'|' -f3)
    STRAND=$(echo "$row"     | cut -d'|' -f4)
    printf '%s\t%s\t%s\t%s\n' "$CHROM" "$GENE_START" "$GENE_END" "$STRAND" > "$GENE_COORDS"
    echo "  ${CHROM}:${GENE_START}-${GENE_END}  strand=${STRAND}  ($(( GENE_END - GENE_START )) bp)"
}

lookup_gene_gtf() {
    echo "=== [1] Parsing GTF for gene: ${GENE_NAME} ==="
    python3 - "$ANNOTATION" "$GENE_NAME" "$GENE_COORDS" <<'PYEOF'
import sys, gzip, re

gtf_path, gene_name, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
opener = gzip.open if gtf_path.endswith(".gz") else open
pattern = re.compile(r'gene_name\s+"([^"]+)"|gene_id\s+"([^"]+)"')

found = []
with opener(gtf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        for m in pattern.finditer(fields[8]):
            name = m.group(1) or m.group(2)
            if name == gene_name:
                found.append((fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]))
                break

if not found:
    print(f"ERROR: gene '{gene_name}' not found in GTF", file=sys.stderr)
    sys.exit(1)

chrom  = found[0][0]
start  = min(r[1] for r in found)
end    = max(r[2] for r in found)
strand = found[0][3]

with open(out_path, "w") as out:
    out.write(f"{chrom}\t{start}\t{end}\t{strand}\n")
print(f"  {chrom}:{start}-{end}  strand={strand}  ({end-start:,} bp)")
PYEOF
    read -r CHROM GENE_START GENE_END STRAND < "$GENE_COORDS"
}

# =============================================================================
# ── SERVER MODE ───────────────────────────────────────────────────────────────
# =============================================================================

if [[ "$MODE" == "server" ]]; then

# ---------------------------------------------------------------------------
# S0. Check / start pgr-server
# ---------------------------------------------------------------------------
t0=$SECONDS
check_or_start_server
log_wall "[S0] check/start pgr-server" $(( SECONDS - t0 ))

# ---------------------------------------------------------------------------
# S1. Query via the HTTP API
#
#   .db annotation → POST /query/gene  (server resolves gene + query in one shot)
#   GTF annotation → local gene lookup + POST /query/region
# ---------------------------------------------------------------------------
RESP_FILE="$OUT_DIR/${GENE_NAME}_server_response.json"
REQ_FILE="$OUT_DIR/${GENE_NAME}_server_request.json"

if [[ "$ANNOTATION" == *.db ]]; then
    echo
    echo "=== [S1] POST /api/v1/query/gene — ${GENE_NAME} ==="
    python3 -c "
import json, sys
print(json.dumps({
    'gene_name': '${GENE_NAME}',
    'flank': ${FLANK},
    'include_sequences': True,
    'min_anchor_count': 10,
    'merge_range_tol': 100000,
    'max_count': 128,
    'max_query_count': 128,
    'max_target_count': 128,
}))
" > "$REQ_FILE"

    t0=$SECONDS
    server_post /api/v1/query/gene "$REQ_FILE" > "$RESP_FILE"
    log_wall "[S1] POST /query/gene" $(( SECONDS - t0 ))

else
    # GTF: resolve gene locally, then POST /query/region
    if [[ "$ANNOTATION" == *.gtf || "$ANNOTATION" == *.gtf.gz ]]; then
        lookup_gene_gtf
    else
        echo "ERROR: in server mode, annotation must be a .db, .gtf, or .gtf.gz file" >&2
        exit 1
    fi

    QUERY_START=$(( GENE_START - FLANK ))
    QUERY_END=$(( GENE_END   + FLANK ))
    (( QUERY_START < 0 )) && QUERY_START=0
    echo "    Query region with ${FLANK} bp flanks: ${CHROM}:${QUERY_START}-${QUERY_END}"

    echo
    echo "=== [S1] POST /api/v1/query/region — ${CHROM}:${QUERY_START}-${QUERY_END} ==="
    python3 -c "
import json
print(json.dumps({
    'chrom': '${CHROM}',
    'start': ${QUERY_START},
    'end':   ${QUERY_END},
    'ref_sample_hint': '${REF_SAMPLE_HINT}',
    'include_sequences': True,
    'min_anchor_count': 10,
    'merge_range_tol': 100000,
    'max_count': 128,
    'max_query_count': 128,
    'max_target_count': 128,
}))
" > "$REQ_FILE"

    t0=$SECONDS
    server_post /api/v1/query/region "$REQ_FILE" > "$RESP_FILE"
    log_wall "[S1] POST /query/region" $(( SECONDS - t0 ))
fi

# Parse the query response → gene_coords.tsv, hit summary, hit FASTA
python3 - "$RESP_FILE" "$GENE_COORDS" "$HIT_PREFIX" "$HIT_FA" <<'PYEOF'
import json, sys

resp_file, coords_file, hit_prefix, hit_fa = sys.argv[1:]

with open(resp_file) as f:
    resp = json.load(f)

qr   = resp['query_region']
hits = resp['hits']
fasta = resp.get('sequences_fasta') or ''

# Gene metadata (present only in /query/gene responses)
gene = resp.get('gene')
if gene:
    chrom, gstart, gend, strand = gene['chrom'], gene['start'], gene['end'], gene['strand']
    print(f"  Gene   : {gene['gene_name']}  {chrom}:{gstart}-{gend}  strand={strand}")
    with open(coords_file, 'w') as f:
        f.write(f"{chrom}\t{gstart}\t{gend}\t{strand}\n")
else:
    chrom = qr['ctg'].split('#')[-1] if '#' in qr['ctg'] else qr['ctg']

print(f"  Region : {qr['ctg']}:{qr['bgn']}-{qr['end']}")
print(f"  Hits   : {len(hits)}")

# Hit summary (mimics local mode .000.hit format)
with open(f"{hit_prefix}.000.hit", 'w') as f:
    f.write("# idx\tq_ctg\tq_bgn\tq_end\tanchors\tsrc\tctg\tctg_bgn\tctg_end\torientation\tname\n")
    for i, h in enumerate(hits):
        f.write(f"{i:03d}\t{qr['ctg']}\t{qr['bgn']}\t{qr['end']}\t"
                f"{h['anchor_count']}\t{h['src']}\t{h['ctg']}\t"
                f"{h['bgn']}\t{h['end']}\t{h['orientation']}\t{h['name']}\n")

# Hit FASTA
with open(hit_fa, 'w') as f:
    f.write(fasta)

print(f"  Written: {hit_prefix}.000.hit  ({len(hits)} hits)")
print(f"  Written: {hit_fa}  ({len(fasta)} bytes)")
PYEOF

echo
echo "=== Hit summary (${HIT_PREFIX}.000.hit) ==="
echo "# idx  q_ctg  q_bgn  q_end  anchors  src  ctg  ctg_bgn  ctg_end  orientation  name"
cat "${HIT_PREFIX}.000.hit"

if [[ ! -s "$HIT_FA" ]]; then
    echo
    echo "NOTE: no hit sequences — skipping bundle decomposition."
    echo
    echo "=== Performance summary ($PERF_LOG) ==="
    print_perf_summary
    echo
    echo "=== Done ==="
    exit 0
fi

# ---------------------------------------------------------------------------
# S2. Principal bundle decomposition + HTML via /bundle
# ---------------------------------------------------------------------------
BUNDLE_REQ="$OUT_DIR/${GENE_NAME}_bundle_request.json"
BUNDLE_RESP="$OUT_DIR/${GENE_NAME}_bundle_response.json"

echo
echo "=== [S2] POST /api/v1/bundle ==="
python3 -c "
import json, sys

hit_fa = sys.argv[1]
req_out = sys.argv[2]

with open(hit_fa) as f:
    fasta = f.read()

if not fasta.strip():
    print('NOTE: empty hit FASTA — nothing to bundle', file=sys.stderr)
    sys.exit(0)

with open(req_out, 'w') as f:
    json.dump({'sequences_fasta': fasta, 'generate_html': True}, f)
" "$HIT_FA" "$BUNDLE_REQ"

t0=$SECONDS
server_post /api/v1/bundle "$BUNDLE_REQ" > "$BUNDLE_RESP"
log_wall "[S2] POST /bundle" $(( SECONDS - t0 ))

python3 - "$BUNDLE_RESP" "$BUNDLE_PREFIX" <<'PYEOF'
import json, sys

resp_file, prefix = sys.argv[1], sys.argv[2]

with open(resp_file) as f:
    resp = json.load(f)

bed  = resp.get('bed', '')
html = resp.get('html', '')

with open(f'{prefix}.bed', 'w') as f:
    f.write(bed)
print(f"  Written: {prefix}.bed  ({len(bed.splitlines())} lines)")

if html:
    with open(f'{prefix}.html', 'w') as f:
        f.write(html)
    print(f"  Written: {prefix}.html")
else:
    print("  No HTML generated (generate_html was false or bundle empty)")
PYEOF

# ---------------------------------------------------------------------------
# Performance summary and done
# ---------------------------------------------------------------------------
echo
echo "=== Performance summary (wall time; server processes not measured) ==="
print_perf_summary

echo
echo "=== Done (server mode) ==="
echo "  Hit FASTA   : $HIT_FA"
echo "  Bundle BED  : ${BUNDLE_PREFIX}.bed"
echo "  Bundle HTML : ${BUNDLE_PREFIX}.html"
echo "  Perf log    : $PERF_LOG"
echo
echo "  Server      : $SERVER_URL  (still running — reuse with MODE=server)"
echo "  Server PID  : $(cat "${DB_PREFIX}.pgr-server.pid" 2>/dev/null || echo n/a)"
echo "  Server log  : ${DB_PREFIX}.pgr-server.log"

exit 0   # end of server mode
fi  # [[ "$MODE" == "server" ]]

# =============================================================================
# ── LOCAL MODE (original pipeline) ───────────────────────────────────────────
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Look up gene coordinates → chromosome, start (0-based), end
# ---------------------------------------------------------------------------
if [[ "$ANNOTATION" == *.db ]]; then
    lookup_gene_sqlite
else
    lookup_gene_gtf
    read -r CHROM GENE_START GENE_END STRAND < "$GENE_COORDS"
fi

QUERY_START=$(( GENE_START - FLANK ))
QUERY_END=$(( GENE_END   + FLANK ))
(( QUERY_START < 0 )) && QUERY_START=0
echo "    Query region with ${FLANK} bp flanks: ${CHROM}:${QUERY_START}-${QUERY_END}"

# ---------------------------------------------------------------------------
# 2. List all sequences in the DB and find the matching reference contig.
#    Cache is stored beside the DB prefix (shared across all gene queries).
# ---------------------------------------------------------------------------
SEQ_LIST="${DB_PREFIX}.seq_list.tsv"
if [[ ! -s "$SEQ_LIST" ]]; then
    echo
    echo "=== [2] Listing sequences in ${DB_PREFIX}.midx ==="
    run_timed "[2] pgr query fetch --list" \
        "$PGR" query fetch \
            --pgr-db-prefix "$DB_PREFIX" \
            --list \
            --output-file "$SEQ_LIST"
    echo "    $(wc -l < "$SEQ_LIST") sequences indexed"
else
    echo
    echo "[SKIP] ${DB_PREFIX}.seq_list.tsv already cached"
fi

REF_REGION=$(python3 - "$SEQ_LIST" "$REF_SAMPLE_HINT" "$CHROM" "$QUERY_START" "$QUERY_END" <<'PYEOF'
import sys

seq_list, hint, chrom, q_start, q_end = \
    sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

with open(seq_list) as fh:
    for line in fh:
        sid, src, ctg, length = line.rstrip().split("\t")
        length = int(length)
        if hint.lower() in src.lower() and ctg.split("#")[-1] == chrom:
            bgn = max(0, q_start)
            end = min(length, q_end)
            print(f"{chrom}:{bgn}-{end}_{src}\t{src}\t{ctg}\t{bgn}\t{end}\t0")
            break
    else:
        print(f"ERROR: no reference contig matching '{hint}' / '{chrom}' found",
              file=sys.stderr)
        sys.exit(1)
PYEOF
)
echo "    Reference region: $REF_REGION"

# ---------------------------------------------------------------------------
# 3. Fetch the reference sequence for the locus from the agcrs
# ---------------------------------------------------------------------------
if [[ ! -s "$QUERY_FA" ]]; then
    echo
    echo "=== [3] Fetching reference sequence from agcrs ==="
    REF_REGION_FILE="$OUT_DIR/ref_region.tsv"
    echo "$REF_REGION" > "$REF_REGION_FILE"
    run_timed "[3] pgr query fetch (ref seq)" \
        "$PGR" query fetch \
            --pgr-db-prefix "$DB_PREFIX" \
            --region-file   "$REF_REGION_FILE" \
            --output-file   "$QUERY_FA"
    echo "    Written: $QUERY_FA  ($(wc -c < "$QUERY_FA" | tr -d ' ') bytes)"
else
    echo
    echo "[SKIP] $QUERY_FA already exists"
fi

# ---------------------------------------------------------------------------
# 4. Query all haplotype contigs covering the locus
# ---------------------------------------------------------------------------
echo
echo "=== [4] pgr query seqs — ${GENE_NAME} vs pangenome ==="
run_timed "[4] pgr query seqs" \
    "$PGR" query seqs \
        --pgr-db-prefix    "$DB_PREFIX" \
        --query-fastx-path "$QUERY_FA" \
        --output-prefix    "$HIT_PREFIX" \
        --memory-mode      low \
        --merge-range-tol  100000 \
        --max-count        128 \
        --max-query-count  128 \
        --max-target-count 128 \
        --min-anchor-count 10

echo
echo "=== Hit summary (${HIT_PREFIX}.000.hit) ==="
echo "# columns: idx  query  q_bgn  q_end  q_len  anchors  src  contig  bgn  end  orient  name"
cat "${HIT_PREFIX}.000.hit"

if [[ ! -s "$HIT_FA" ]]; then
    echo
    echo "NOTE: no hit sequences found — skipping bundle decomposition."
    exit 0
fi

# ---------------------------------------------------------------------------
# 5. Principal bundle decomposition
# ---------------------------------------------------------------------------
echo
echo "=== [5] pgr bundle decomp ==="
run_timed "[5] pgr bundle decomp" \
    "$PGR" bundle decomp \
        --fastx-path    "$HIT_FA" \
        --output-prefix "$BUNDLE_PREFIX"

# ---------------------------------------------------------------------------
# 6. Interactive HTML bundle visualisation
# ---------------------------------------------------------------------------
echo
echo "=== [6] pgr bundle svg --html ==="
run_timed "[6] pgr bundle svg --html" \
    "$PGR" bundle svg \
        --bed-file-path "${BUNDLE_PREFIX}.bed" \
        --output-prefix "$BUNDLE_PREFIX" \
        --html

# ---------------------------------------------------------------------------
# Performance summary
# ---------------------------------------------------------------------------
echo
echo "=== Performance summary ($PERF_LOG) ==="
print_perf_summary

echo
echo "=== Done (local mode) ==="
echo "  Query FASTA : $QUERY_FA"
echo "  Hit FASTA   : $HIT_FA"
echo "  Bundle BED  : ${BUNDLE_PREFIX}.bed"
echo "  Bundle HTML : ${BUNDLE_PREFIX}.html"
echo "  Perf log    : $PERF_LOG"
