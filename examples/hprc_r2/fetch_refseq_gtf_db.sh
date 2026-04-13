#!/usr/bin/env bash
# fetch_refseq_gtf_db.sh — download hg38.ncbiRefSeq.gtf.gz from UCSC and
# load it into a SQLite gene-annotation database for fast coordinate lookups.
#
# The resulting database is consumed by query_gene_bundle.sh in place of the
# raw (gzipped) GTF, avoiding repeated decompression and parsing.
#
# Schema
# ------
#   genes        — one row per gene (gene_id, gene_name, chrom, start, end, strand)
#   transcripts  — one row per transcript linked to its gene
#   exons        — one row per exon linked to its transcript
#
# Indexes on gene_name and gene_id support O(log n) coordinate lookup.
#
# Usage:
#   bash fetch_refseq_gtf_db.sh [out_db] [gtf_gz]
#
#   out_db  — output SQLite file           (default: hg38.ncbiRefSeq.db)
#   gtf_gz  — local GTF.gz to use instead of downloading  (optional)

set -euo pipefail

OUT_DB="${1:-hg38.ncbiRefSeq.db}"
GTF_GZ="${2:-}"

UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"

# ---------------------------------------------------------------------------
# 1. Obtain the GTF
# ---------------------------------------------------------------------------
if [[ -n "$GTF_GZ" ]]; then
    [[ -f "$GTF_GZ" ]] || { echo "ERROR: $GTF_GZ not found" >&2; exit 1; }
    echo "[INFO] Using local GTF: $GTF_GZ"
else
    GTF_GZ="hg38.ncbiRefSeq.gtf.gz"
    if [[ ! -f "$GTF_GZ" ]]; then
        echo "=== Downloading $GTF_GZ from UCSC ==="
        curl -fL --retry 3 --retry-delay 5 -o "$GTF_GZ" "$UCSC_URL"
        echo "    Downloaded: $GTF_GZ  ($(du -sh "$GTF_GZ" | cut -f1))"
    else
        echo "[SKIP] $GTF_GZ already exists"
    fi
fi

# ---------------------------------------------------------------------------
# 2. Parse GTF and populate SQLite database
# ---------------------------------------------------------------------------
echo
echo "=== Building SQLite database: $OUT_DB ==="

python3 - "$GTF_GZ" "$OUT_DB" <<'PYEOF'
import sys, gzip, re, sqlite3

gtf_gz, db_path = sys.argv[1], sys.argv[2]

# ── schema ──────────────────────────────────────────────────────────────────
SCHEMA = """
PRAGMA journal_mode = WAL;

CREATE TABLE IF NOT EXISTS genes (
    gene_pk   INTEGER PRIMARY KEY,
    gene_id   TEXT    NOT NULL,
    gene_name TEXT    NOT NULL,
    chrom     TEXT    NOT NULL,
    start     INTEGER NOT NULL,   -- 0-based
    end       INTEGER NOT NULL,   -- exclusive
    strand    TEXT    NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS idx_genes_id_chrom ON genes(gene_id, chrom);
CREATE        INDEX IF NOT EXISTS idx_genes_name     ON genes(gene_name);
CREATE        INDEX IF NOT EXISTS idx_genes_loc      ON genes(chrom, start, end);

CREATE TABLE IF NOT EXISTS transcripts (
    tx_pk       INTEGER PRIMARY KEY,
    gene_pk     INTEGER NOT NULL REFERENCES genes(gene_pk),
    tx_id       TEXT    NOT NULL,
    chrom       TEXT    NOT NULL,
    start       INTEGER NOT NULL,
    end         INTEGER NOT NULL,
    strand      TEXT    NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS idx_tx_id  ON transcripts(tx_id);
CREATE        INDEX IF NOT EXISTS idx_tx_gene ON transcripts(gene_pk);

CREATE TABLE IF NOT EXISTS exons (
    exon_pk  INTEGER PRIMARY KEY,
    tx_pk    INTEGER NOT NULL REFERENCES transcripts(tx_pk),
    start    INTEGER NOT NULL,
    end      INTEGER NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_exon_tx ON exons(tx_pk);
"""

attr_re = re.compile(r'(\w+)\s+"([^"]+)"')

def parse_attrs(s):
    return dict(attr_re.findall(s))

con = sqlite3.connect(db_path)
con.executescript(SCHEMA)

opener = gzip.open if gtf_gz.endswith(".gz") else open

gene_rows = {}   # (gene_id, chrom) -> (gene_name, start, end, strand)
tx_rows   = {}   # tx_id -> (gene_id, chrom, start, end, strand)
exon_rows = []   # (tx_id, start, end)

n_lines = 0
with opener(gtf_gz, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip().split("\t")
        if len(f) < 9:
            continue
        chrom, feature, start, end, strand, attrs = \
            f[0], f[2], int(f[3])-1, int(f[4]), f[6], f[8]
        a = parse_attrs(attrs)
        gene_id   = a.get("gene_id", "")
        gene_name = a.get("gene_name", gene_id)
        tx_id     = a.get("transcript_id", "")

        key = (gene_id, chrom)
        if feature == "gene":
            if key not in gene_rows:
                gene_rows[key] = (gene_name, start, end, strand)
            else:
                g = gene_rows[key]
                gene_rows[key] = (gene_name, min(g[1], start), max(g[2], end), strand)
        elif feature == "transcript" and tx_id:
            tx_rows[tx_id] = (gene_id, chrom, start, end, strand)
            # ncbiRefSeq GTF has no "gene" lines — derive gene span from transcripts,
            # keyed by (gene_id, chrom) so alt contigs stay separate
            if key not in gene_rows:
                gene_rows[key] = (gene_name, start, end, strand)
            else:
                g = gene_rows[key]
                gene_rows[key] = (gene_name, min(g[1], start), max(g[2], end), strand)
        elif feature == "exon" and tx_id:
            exon_rows.append((tx_id, start, end))

        n_lines += 1
        if n_lines % 500_000 == 0:
            print(f"  … {n_lines:,} lines parsed", flush=True)

print(f"  {n_lines:,} GTF lines  |  "
      f"{len(gene_rows):,} gene-chrom entries  |  "
      f"{len(tx_rows):,} transcripts  |  "
      f"{len(exon_rows):,} exons", flush=True)

# ── insert ───────────────────────────────────────────────────────────────────
print("  Inserting genes …", flush=True)
gene_key_to_pk = {}   # (gene_id, chrom) -> pk
with con:
    for (gene_id, chrom), (gene_name, start, end, strand) in gene_rows.items():
        cur = con.execute(
            "INSERT OR IGNORE INTO genes(gene_id,gene_name,chrom,start,end,strand) "
            "VALUES(?,?,?,?,?,?)",
            (gene_id, gene_name, chrom, start, end, strand)
        )
        pk = cur.lastrowid or \
            con.execute("SELECT gene_pk FROM genes WHERE gene_id=? AND chrom=?",
                        (gene_id, chrom)).fetchone()[0]
        gene_key_to_pk[(gene_id, chrom)] = pk

print("  Inserting transcripts …", flush=True)
tx_id_to_pk = {}
with con:
    for tx_id, (gene_id, chrom, start, end, strand) in tx_rows.items():
        gene_pk = gene_key_to_pk.get((gene_id, chrom))
        if gene_pk is None:
            continue
        cur = con.execute(
            "INSERT OR IGNORE INTO transcripts(gene_pk,tx_id,chrom,start,end,strand) "
            "VALUES(?,?,?,?,?,?)",
            (gene_pk, tx_id, chrom, start, end, strand)
        )
        tx_id_to_pk[tx_id] = cur.lastrowid or \
            con.execute("SELECT tx_pk FROM transcripts WHERE tx_id=?",
                        (tx_id,)).fetchone()[0]

print("  Inserting exons …", flush=True)
with con:
    con.executemany(
        "INSERT INTO exons(tx_pk,start,end) VALUES(?,?,?)",
        [(tx_id_to_pk[tid], s, e)
         for tid, s, e in exon_rows if tid in tx_id_to_pk]
    )

con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
con.close()
print("  Done.", flush=True)
PYEOF

echo
echo "=== Database ready: $OUT_DB ==="
sqlite3 "$OUT_DB" "
SELECT 'genes       : ' || COUNT(*) FROM genes;
SELECT 'transcripts : ' || COUNT(*) FROM transcripts;
SELECT 'exons       : ' || COUNT(*) FROM exons;
"
echo
echo "Example query:"
echo "  sqlite3 $OUT_DB \"SELECT chrom,start,end,strand FROM genes WHERE gene_name='C9orf72'\""
