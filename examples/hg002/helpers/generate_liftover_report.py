#!/usr/bin/env python3
"""
generate_liftover_report.py — produce a standalone GTF liftover HTML report
from hg002_hap0_liftover.db and hg002_hap1_liftover.db.

Tabs:
  Summary        — headline KPI cards
  Liftover Status — transcript status breakdown (hap0 / hap1)
  NM Genes       — NM gene funnel + coverage distribution chart
  Scatter Plots  — exon-length and genomic-span scatter (linear + log)
  Per-Chromosome — per-chrom NM gene coverage table

Usage:
    python3 generate_liftover_report.py [--base-dir DIR] [--out FILE]

Arguments:
    --base-dir DIR   Directory containing the liftover .db files
                     (default: directory of this script)
    --out FILE       Output HTML path (default: <base-dir>/liftover_report.html)
"""

import sys, sqlite3, json, html, os, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--base-dir", default=os.path.dirname(os.path.abspath(__file__)),
                help="Directory containing hg002_hap{0,1}_liftover.db")
ap.add_argument("--out", default=None,
                help="Output HTML path (default: <base-dir>/liftover_report.html)")
args = ap.parse_args()

base_dir = os.path.abspath(args.base_dir)
out_path = args.out if args.out else os.path.join(base_dir, "liftover_report.html")

db0_path = os.path.join(base_dir, "hg002_hap0_liftover.db")
db1_path = os.path.join(base_dir, "hg002_hap1_liftover.db")

for p in (db0_path, db1_path):
    if not os.path.exists(p):
        print(f"ERROR: {p} not found — run 04_liftover_gtf.sh first", file=sys.stderr)
        sys.exit(1)


def query(db_path, sql, params=()):
    conn = sqlite3.connect(db_path)
    cur = conn.execute(sql, params)
    rows = cur.fetchall()
    conn.close()
    return rows


# ── 1. Transcript-level status counts ───────────────────────────────────────
def tx_status(db):
    rows = query(db, """
        SELECT status, COUNT(*) AS n
        FROM transcript_summary
        GROUP BY status
        ORDER BY n DESC
    """)
    total = query(db, "SELECT COUNT(*) FROM transcript_summary")[0][0]
    return total, {r[0]: r[1] for r in rows}


total0, st0 = tx_status(db0_path)
total1, st1 = tx_status(db1_path)

STATUS_ORDER = ["single_full", "single_partial", "multi_contig", "multi_location", "no_hit"]
STATUS_COLOR = {
    "single_full":    "#2ecc71",
    "single_partial": "#f39c12",
    "multi_contig":   "#3498db",
    "multi_location": "#9b59b6",
    "no_hit":         "#e74c3c",
}

# ── 2. NM curated mRNA gene funnel ─────────────────────────────────────────
NM_FUNNEL_SQL = """
WITH nm_primary AS (
  SELECT DISTINCT gene_name FROM transcripts
  WHERE SUBSTR(transcript_id,1,2)='NM'
    AND transcript_id NOT GLOB '*_[0-9]'
    AND transcript_id NOT GLOB '*_[0-9][0-9]'
    AND gene_name != ''
),
best AS (
  SELECT t.gene_name, MAX(l.coverage_pct) AS bc
  FROM transcripts t JOIN liftover l ON t.transcript_pk=l.transcript_pk
  WHERE SUBSTR(t.transcript_id,1,2)='NM'
    AND t.transcript_id NOT GLOB '*_[0-9]'
    AND t.transcript_id NOT GLOB '*_[0-9][0-9]'
    AND t.gene_name!=''
  GROUP BY t.gene_name
)
SELECT
  COUNT(*),
  SUM(CASE WHEN b.bc IS NOT NULL THEN 1 ELSE 0 END),
  SUM(CASE WHEN b.bc IS NULL     THEN 1 ELSE 0 END),
  SUM(CASE WHEN b.bc>=90         THEN 1 ELSE 0 END),
  SUM(CASE WHEN b.bc>=50 AND b.bc<90 THEN 1 ELSE 0 END)
FROM nm_primary p LEFT JOIN best b ON p.gene_name=b.gene_name
"""


def nm_funnel(db):
    r = query(db, NM_FUNNEL_SQL)[0]
    return dict(total=r[0], any_hit=r[1], no_hit=r[2], hq=r[3], partial=r[4])


f0 = nm_funnel(db0_path)
f1 = nm_funnel(db1_path)

# ── 3. HQ gene counts (all types) ───────────────────────────────────────────
HQ_SQL = """
SELECT COUNT(DISTINCT t.gene_name)
FROM transcripts t JOIN liftover l ON t.transcript_pk=l.transcript_pk
WHERE l.coverage_pct>=90 AND t.gene_name!=''
"""
hq0 = query(db0_path, HQ_SQL)[0][0]
hq1 = query(db1_path, HQ_SQL)[0][0]
total_genes = query(db0_path,
    "SELECT COUNT(DISTINCT gene_name) FROM transcripts WHERE gene_name!=''"
)[0][0]

# ── 4. HQ by transcript type ────────────────────────────────────────────────
TX_TYPE_SQL = """
SELECT
  CASE SUBSTR(t.transcript_id,1,2)
    WHEN 'NM' THEN 'NM curated mRNA'
    WHEN 'XM' THEN 'XM predicted mRNA'
    WHEN 'NR' THEN 'NR curated ncRNA'
    WHEN 'XR' THEN 'XR predicted ncRNA'
    ELSE            'other'
  END AS tx_type,
  COUNT(DISTINCT t.gene_name) AS hq_genes
FROM transcripts t JOIN liftover l ON t.transcript_pk=l.transcript_pk
WHERE l.coverage_pct>=90 AND t.gene_name!=''
  AND t.transcript_id NOT GLOB '*_[0-9]'
  AND t.transcript_id NOT GLOB '*_[0-9][0-9]'
GROUP BY tx_type ORDER BY hq_genes DESC
"""
tx_types0 = query(db0_path, TX_TYPE_SQL)
tx_types1 = query(db1_path, TX_TYPE_SQL)

# ── 5. Per-chromosome HQ NM gene counts ─────────────────────────────────────
CHR_SQL = """
WITH best AS (
  SELECT t.gene_name, t.ref_chrom, MAX(l.coverage_pct) AS bc
  FROM transcripts t JOIN liftover l ON t.transcript_pk=l.transcript_pk
  WHERE SUBSTR(t.transcript_id,1,2)='NM'
    AND t.transcript_id NOT GLOB '*_[0-9]'
    AND t.transcript_id NOT GLOB '*_[0-9][0-9]'
    AND t.gene_name!=''
    AND t.ref_chrom NOT LIKE '%_alt%'
    AND t.ref_chrom NOT LIKE '%_fix%'
    AND t.ref_chrom NOT LIKE '%_random%'
    AND t.ref_chrom NOT LIKE 'chrUn%'
  GROUP BY t.gene_name, t.ref_chrom
),
totals AS (
  SELECT ref_chrom,
         COUNT(DISTINCT gene_name)                           AS total,
         SUM(CASE WHEN bc>=90 THEN 1 ELSE 0 END)            AS hq,
         SUM(CASE WHEN bc>=50 AND bc<90 THEN 1 ELSE 0 END)  AS partial,
         SUM(CASE WHEN bc IS NULL        THEN 1 ELSE 0 END)  AS no_hit
  FROM (
    SELECT DISTINCT p.gene_name, p.ref_chrom, b.bc
    FROM (SELECT DISTINCT gene_name, ref_chrom FROM transcripts
          WHERE SUBSTR(transcript_id,1,2)='NM'
            AND transcript_id NOT GLOB '*_[0-9]'
            AND transcript_id NOT GLOB '*_[0-9][0-9]'
            AND gene_name!=''
            AND ref_chrom NOT LIKE '%_alt%'
            AND ref_chrom NOT LIKE '%_fix%'
            AND ref_chrom NOT LIKE '%_random%'
            AND ref_chrom NOT LIKE 'chrUn%') p
    LEFT JOIN best b ON p.gene_name=b.gene_name AND p.ref_chrom=b.ref_chrom
  )
  GROUP BY ref_chrom
)
SELECT ref_chrom, total, hq, partial, no_hit,
       ROUND(100.0*hq/total,1) AS hq_pct
FROM totals
ORDER BY
  CASE ref_chrom
    WHEN 'chr1'  THEN 1  WHEN 'chr2'  THEN 2  WHEN 'chr3'  THEN 3
    WHEN 'chr4'  THEN 4  WHEN 'chr5'  THEN 5  WHEN 'chr6'  THEN 6
    WHEN 'chr7'  THEN 7  WHEN 'chr8'  THEN 8  WHEN 'chr9'  THEN 9
    WHEN 'chr10' THEN 10 WHEN 'chr11' THEN 11 WHEN 'chr12' THEN 12
    WHEN 'chr13' THEN 13 WHEN 'chr14' THEN 14 WHEN 'chr15' THEN 15
    WHEN 'chr16' THEN 16 WHEN 'chr17' THEN 17 WHEN 'chr18' THEN 18
    WHEN 'chr19' THEN 19 WHEN 'chr20' THEN 20 WHEN 'chr21' THEN 21
    WHEN 'chr22' THEN 22 WHEN 'chrX'  THEN 23 WHEN 'chrY'  THEN 24
    WHEN 'chrM'  THEN 25 ELSE 99
  END
"""
chr0 = query(db0_path, CHR_SQL)
chr1 = query(db1_path, CHR_SQL)

# ── 6. Coverage distribution of NM genes (hap0) ─────────────────────────────
COV_DIST_SQL = """
WITH best AS (
  SELECT t.gene_name, MAX(l.coverage_pct) AS bc
  FROM transcripts t JOIN liftover l ON t.transcript_pk=l.transcript_pk
  WHERE SUBSTR(t.transcript_id,1,2)='NM'
    AND t.transcript_id NOT GLOB '*_[0-9]'
    AND t.transcript_id NOT GLOB '*_[0-9][0-9]'
    AND t.gene_name!=''
  GROUP BY t.gene_name
)
SELECT ROUND(bc/5)*5 AS bin, COUNT(*) AS n
FROM best WHERE bc IS NOT NULL
GROUP BY bin ORDER BY bin
"""
cov_dist = query(db0_path, COV_DIST_SQL)

# ── 7. Scatter plot data ─────────────────────────────────────────────────────
SCATTER_SQL = """
WITH ref_exon AS (
  SELECT e.transcript_pk, SUM(e.exon_end - e.exon_start) AS ref_exon_len
  FROM exons e GROUP BY e.transcript_pk
),
best_hit AS (
  SELECT l.transcript_pk,
         (l.contig_end - l.contig_start) AS contig_genomic_span,
         l.block_sizes
  FROM liftover l
  JOIN transcript_summary ts ON l.transcript_pk = ts.transcript_pk
  WHERE l.contig = ts.best_contig
    AND ts.status IN ('single_full', 'single_partial')
)
SELECT t.ref_end - t.ref_start  AS ref_genomic_span,
       r.ref_exon_len,
       b.contig_genomic_span,
       b.block_sizes,
       t.transcript_id,
       t.gene_name
FROM transcripts t
JOIN ref_exon r ON t.transcript_pk = r.transcript_pk
JOIN best_hit b ON t.transcript_pk = b.transcript_pk
ORDER BY t.transcript_pk
"""


def parse_block_sizes(bs):
    return sum(int(x) for x in bs.split(',') if x.strip())


print("Querying liftover databases ...", flush=True)
scatter_raw0 = query(db0_path, SCATTER_SQL)
scatter_raw1 = query(db1_path, SCATTER_SQL)

sc_exon0 = [{"x": round(r[1]/1000, 2), "y": round(parse_block_sizes(r[3])/1000, 2),
             "n": f"{r[5]} ({r[4]})"}
            for r in scatter_raw0]
sc_exon1 = [{"x": round(r[1]/1000, 2), "y": round(parse_block_sizes(r[3])/1000, 2),
             "n": f"{r[5]} ({r[4]})"}
            for r in scatter_raw1]
sc_span0 = [{"x": round(r[0]/1000, 2), "y": round(r[2]/1000, 2),
             "n": f"{r[5]} ({r[4]})"}
            for r in scatter_raw0]
sc_span1 = [{"x": round(r[0]/1000, 2), "y": round(r[2]/1000, 2),
             "n": f"{r[5]} ({r[4]})"}
            for r in scatter_raw1]

sc_exon0_js = json.dumps(sc_exon0)
sc_exon1_js = json.dumps(sc_exon1)
sc_span0_js = json.dumps(sc_span0)
sc_span1_js = json.dumps(sc_span1)


# ── Helpers ──────────────────────────────────────────────────────────────────
def bar(pct_val, color, width=120):
    w = int(pct_val / 100 * width)
    return (f'<div style="display:inline-block;width:{w}px;height:14px;'
            f'background:{color};border-radius:2px;vertical-align:middle"></div> '
            f'<span style="font-variant-numeric:tabular-nums">{pct_val:.1f}%</span>')


def tx_status_rows(total, st):
    rows = ""
    for s in STATUS_ORDER:
        n = st.get(s, 0)
        rows += (f"<tr><td>{html.escape(s)}</td><td class='num'>{n:,}</td>"
                 f"<td>{bar(n/total*100, STATUS_COLOR[s])}</td></tr>")
    return rows


def nm_funnel_html(f, label):
    items = [
        ("Total NM genes",  f['total'],   "#7f8c8d"),
        ("Has any hit",     f['any_hit'], "#3498db"),
        ("No hit at all",   f['no_hit'],  "#e74c3c"),
        ("HQ ≥ 90%",       f['hq'],      "#2ecc71"),
        ("Partial 50–89%", f['partial'], "#f39c12"),
    ]
    s = (f"<h3>{html.escape(label)}</h3>"
         "<table class='data'><thead><tr>"
         "<th>Category</th><th>Count</th><th>% of total NM</th>"
         "</tr></thead><tbody>")
    for name, n, color in items:
        s += (f"<tr><td>{name}</td><td class='num'>{n:,}</td>"
              f"<td>{bar(n/f['total']*100, color)}</td></tr>")
    s += "</tbody></table>"
    return s


def chr_table(rows, label):
    s = f"<h3>{html.escape(label)}</h3>"
    s += ("<table class='data'><thead><tr>"
          "<th>Chrom</th><th>Total NM</th><th>HQ ≥90%</th>"
          "<th>Partial</th><th>No hit</th><th>HQ%</th>"
          "</tr></thead><tbody>")
    for chrom, total, hq, partial, nohit, hq_pct in rows:
        hq_pct_v = hq_pct if hq_pct is not None else 0.0
        color = ("#2ecc71" if hq_pct_v >= 80 else
                 "#f39c12" if hq_pct_v >= 50 else "#e74c3c")
        s += (f"<tr><td>{html.escape(chrom)}</td>"
              f"<td class='num'>{total:,}</td>"
              f"<td class='num'>{hq:,}</td>"
              f"<td class='num'>{partial:,}</td>"
              f"<td class='num'>{nohit:,}</td>"
              f"<td>{bar(hq_pct_v, color, 80)}</td></tr>")
    s += "</tbody></table>"
    return s


def tx_type_table(rows, label):
    s = f"<h3>{html.escape(label)}</h3>"
    s += "<table class='data'><thead><tr><th>Transcript type</th><th>HQ genes</th></tr></thead><tbody>"
    for name, n in rows:
        s += f"<tr><td>{html.escape(name)}</td><td class='num'>{n:,}</td></tr>"
    s += "</tbody></table>"
    return s


cov_labels_js = json.dumps([str(int(b)) for b, _ in cov_dist])
cov_values_js = json.dumps([n for _, n in cov_dist])

# ── SV candidates ─────────────────────────────────────────────────────────────
SV_TYPE_LABEL = {
    "SV": "Structural variant (size-discordant gap)",
    "TG": "Ref gap (unaligned region between contigs)",
    "TD": "Contig duplicate (contained alignment)",
    "TO": "Contig overlap (partially overlapping)",
    "QG": "Contig gap (unaligned between ref alignments)",
    "QD": "Ref duplicate (ref block contained)",
    "QO": "Ref overlap (ref block partially overlapping)",
}
SV_TYPE_COLOR = {
    "SV": "#e74c3c", "TG": "#e67e22", "TD": "#f39c12", "TO": "#f1c40f",
    "QG": "#3498db", "QD": "#9b59b6", "QO": "#1abc9c",
}

def parse_svcnd_bed(path):
    records = []
    if not os.path.exists(path):
        return records
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            sv_type = parts[3][:2]
            size = end - start
            records.append((chrom, start, end, sv_type, size))
    return records

def svcnd_type_rows(records, codes):
    from collections import Counter
    total = len(records)
    counts = Counter(r[3] for r in records)
    rows = ""
    for code in codes:
        n = counts.get(code, 0)
        color = SV_TYPE_COLOR.get(code, "#95a5a6")
        label = SV_TYPE_LABEL.get(code, code)
        pct_v = n / total * 100 if total else 0
        w = int(pct_v / 100 * 120)
        rows += (f"<tr><td><b>{html.escape(code)}</b></td>"
                 f"<td>{html.escape(label)}</td>"
                 f"<td class='num'>{n:,}</td>"
                 f"<td><div style='display:inline-block;width:{w}px;height:14px;"
                 f"background:{color};border-radius:2px;vertical-align:middle'></div>"
                 f" {pct_v:.1f}%</td></tr>")
    return rows, total

svcnd0 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap0.svcnd.bed"))
svcnd1 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap1.svcnd.bed"))
ctgsv0 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap0.ctgsv.bed"))
ctgsv1 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap1.ctgsv.bed"))

if svcnd0 or svcnd1:
    ref_rows0, ref_total0 = svcnd_type_rows(svcnd0, ["SV", "TG", "TD", "TO"])
    ref_rows1, ref_total1 = svcnd_type_rows(svcnd1, ["SV", "TG", "TD", "TO"])
    ctg_rows0, ctg_total0 = svcnd_type_rows(ctgsv0, ["QG", "QD", "QO"])
    ctg_rows1, ctg_total1 = svcnd_type_rows(ctgsv1, ["QG", "QD", "QO"])
    sv_tab_html = f"""
<h2>SV Candidate Summary</h2>
<p class="note">
  <b>Ref view</b> (.svcnd.bed) — anomalies on the reference axis &nbsp;&middot;&nbsp;
  <b>Contig view</b> (.ctgsv.bed) — anomalies on the contig axis
</p>
<div class="col2">
  <div>
    <h3>Hap0 &mdash; Ref view ({ref_total0:,} records)</h3>
    <table class="data"><thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th></tr></thead>
    <tbody>{ref_rows0}</tbody></table>
    <h3>Hap0 &mdash; Contig view ({ctg_total0:,} records)</h3>
    <table class="data"><thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th></tr></thead>
    <tbody>{ctg_rows0}</tbody></table>
  </div>
  <div>
    <h3>Hap1 &mdash; Ref view ({ref_total1:,} records)</h3>
    <table class="data"><thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th></tr></thead>
    <tbody>{ref_rows1}</tbody></table>
    <h3>Hap1 &mdash; Contig view ({ctg_total1:,} records)</h3>
    <table class="data"><thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th></tr></thead>
    <tbody>{ctg_rows1}</tbody></table>
  </div>
</div>
"""
    sv_tab_button = '<button onclick="showTab(\'svcnd\')" id="btn-svcnd">SV Candidates</button>'
    sv_tab_panel  = f'<div class="tabpanel" id="tab-svcnd">{sv_tab_html}</div>'
else:
    sv_tab_button = ""
    sv_tab_panel  = ""

# ── Compose page ─────────────────────────────────────────────────────────────
doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PGR-TK GTF Liftover Report — HG002</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
  body  {{ font-family: sans-serif; margin: 0; background: #f9f9f9; color: #222; }}
  h2    {{ color: #2c3e50; }}
  h3    {{ color: #34495e; margin-bottom: .4em; }}
  h4    {{ color: #34495e; margin-bottom: .4em; }}
  .summary-grid {{
    display: grid; grid-template-columns: repeat(auto-fit,minmax(180px,1fr));
    gap: 1em; margin-bottom: 2em;
  }}
  .card {{ background:#fff; border:1px solid #ddd; border-radius:8px;
           padding:1em 1.2em; text-align:center; }}
  .card .val {{ font-size:2em; font-weight:700; color:#2980b9; }}
  .card .lbl {{ font-size:.85em; color:#7f8c8d; margin-top:.2em; }}
  .col2 {{ display:grid; grid-template-columns:1fr 1fr; gap:2em; }}
  table.data {{ border-collapse:collapse; width:100%; background:#fff;
                border:1px solid #ddd; border-radius:6px; overflow:hidden;
                margin-bottom:1.5em; }}
  table.data th {{ background:#2c3e50; color:#fff; padding:.6em 1em; text-align:left; }}
  table.data td {{ padding:.5em 1em; border-bottom:1px solid #eee; }}
  table.data tr:last-child td {{ border-bottom:none; }}
  table.data tr:hover td {{ background:#f0f4f8; }}
  .num  {{ text-align:right; font-variant-numeric:tabular-nums; }}
  .chart-wrap {{ background:#fff; border:1px solid #ddd; border-radius:8px;
                 padding:1em; margin-bottom:1.5em; }}
  .note {{ font-size:.82em; color:#7f8c8d; margin:.3em 0 1em; }}
  .tabbar {{ display: flex; gap: 0; background: #ecf0f1;
             border-bottom: 2px solid #bdc3c7; margin-bottom: 0; }}
  .tabbar button {{
    background: transparent; border: none; color: #7f8c8d;
    padding: .55em 1.3em; font-size: .92em; cursor: pointer;
    border-bottom: 3px solid transparent; margin-bottom: -2px;
    transition: color .15s;
  }}
  .tabbar button:hover  {{ color: #2c3e50; }}
  .tabbar button.active {{ color: #2980b9; border-bottom-color: #2980b9; font-weight: 600; }}
  .tabpanel {{ display: none; padding: 1.5em 2em; }}
  .tabpanel.active {{ display: block; }}
</style>
</head>
<body>
<div class="tabbar">
  <button class="active" onclick="showTab('summary')"  id="btn-summary">Summary</button>
  <button onclick="showTab('status')"   id="btn-status">Liftover Status</button>
  <button onclick="showTab('nmgenes')"  id="btn-nmgenes">NM Genes</button>
  <button onclick="showTab('scatter')"  id="btn-scatter">Scatter Plots</button>
  <button onclick="showTab('perchrom')" id="btn-perchrom">Per-Chromosome</button>
  {sv_tab_button}
</div>

<!-- Tab: Summary -->
<div class="tabpanel active" id="tab-summary">
  <p class="note">
    Reference: GRCh38 &middot; Annotation: hg38.ncbiRefSeq.gtf.gz &middot;
    Hap0 = maternal &middot; Hap1 = paternal &middot;
    HQ threshold: &ge;90% exon coverage on a single contig
  </p>
  <div class="summary-grid">
    <div class="card"><div class="val">{total_genes:,}</div><div class="lbl">total annotated gene names</div></div>
    <div class="card"><div class="val">{hq0:,}</div><div class="lbl">HQ genes &mdash; hap0</div></div>
    <div class="card"><div class="val">{hq1:,}</div><div class="lbl">HQ genes &mdash; hap1</div></div>
    <div class="card"><div class="val">{f0['hq']:,}</div><div class="lbl">HQ NM curated genes &mdash; hap0</div></div>
    <div class="card"><div class="val">{f1['hq']:,}</div><div class="lbl">HQ NM curated genes &mdash; hap1</div></div>
    <div class="card"><div class="val">{f0['no_hit']:,}</div><div class="lbl">NM genes &mdash; no hit &mdash; hap0</div></div>
  </div>
</div>

<!-- Tab: Liftover Status -->
<div class="tabpanel" id="tab-status">
  <h2>Transcript-level Liftover Status</h2>
  <div class="col2">
    <div>
      <h3>Hap0 &mdash; {total0:,} transcripts</h3>
      <table class="data">
        <thead><tr><th>Status</th><th>Count</th><th>Fraction</th></tr></thead>
        <tbody>{tx_status_rows(total0, st0)}</tbody>
      </table>
    </div>
    <div>
      <h3>Hap1 &mdash; {total1:,} transcripts</h3>
      <table class="data">
        <thead><tr><th>Status</th><th>Count</th><th>Fraction</th></tr></thead>
        <tbody>{tx_status_rows(total1, st1)}</tbody>
      </table>
    </div>
  </div>
</div>

<!-- Tab: NM Genes -->
<div class="tabpanel" id="tab-nmgenes">
  <h2>NM Curated mRNA Gene Funnel</h2>
  <div class="col2">
    {nm_funnel_html(f0, "Hap0")}
    {nm_funnel_html(f1, "Hap1")}
  </div>
  <h2>HQ Genes by Transcript Type</h2>
  <div class="col2">
    {tx_type_table(tx_types0, "Hap0 &mdash; HQ gene count by transcript type")}
    {tx_type_table(tx_types1, "Hap1 &mdash; HQ gene count by transcript type")}
  </div>
  <h2>NM Gene Coverage Distribution &mdash; Hap0</h2>
  <div class="chart-wrap">
    <canvas id="covChart" height="80"></canvas>
  </div>
  <p class="note">
    Coverage = best single-contig exon coverage across all NM transcripts for each gene.
    Genes with no hit are excluded (shown in funnel above).
  </p>
</div>

<!-- Tab: Scatter Plots -->
<div class="tabpanel" id="tab-scatter">
  <h2>Transcript Length: Contig vs Reference</h2>
  <p class="note">
    Spliced exon bases aligned to the best-matching contig vs total exon bases in the
    reference annotation. Points near the diagonal indicate complete, undistorted liftover.
    Single-full and single-partial transcripts only. Hap0 (blue) and Hap1 (red) overlaid.
  </p>
  <div class="col2">
    <div class="chart-wrap"><canvas id="scExon"    height="300"></canvas></div>
    <div class="chart-wrap"><canvas id="scExonLog" height="300"></canvas></div>
  </div>
  <h2>Genomic Span: Contig vs Reference</h2>
  <p class="note">
    Genomic span on the contig vs the reference gene body span (including all introns).
    Deviations from the diagonal indicate structural differences.
    Hap0 (blue) and Hap1 (red) overlaid.
  </p>
  <div class="col2">
    <div class="chart-wrap"><canvas id="scSpan"    height="300"></canvas></div>
    <div class="chart-wrap"><canvas id="scSpanLog" height="300"></canvas></div>
  </div>
</div>

<!-- Tab: Per-Chromosome -->
<div class="tabpanel" id="tab-perchrom">
  <h2>Per-Chromosome NM Gene Coverage</h2>
  <div class="col2">
    {chr_table(chr0, "Hap0")}
    {chr_table(chr1, "Hap1")}
  </div>
</div>

<!-- Tab: SV Candidates -->
{sv_tab_panel}

<script>
const _chartsDone = {{}};
function initChartsForTab(id) {{
  if (_chartsDone[id]) return;
  _chartsDone[id] = true;
  if (id === 'nmgenes') initCovChart();
  if (id === 'scatter') initScatterCharts();
}}
function showTab(id) {{
  document.querySelectorAll('.tabpanel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tabbar button').forEach(b => b.classList.remove('active'));
  document.getElementById('tab-' + id).classList.add('active');
  document.getElementById('btn-' + id).classList.add('active');
  initChartsForTab(id);
}}
initChartsForTab('summary');

function scatterCfg(canvasId, data0, data1, title) {{
  let maxVal = 0;
  for (const p of data0.concat(data1)) {{
    if (p.x > maxVal) maxVal = p.x;
    if (p.y > maxVal) maxVal = p.y;
  }}
  new Chart(document.getElementById(canvasId), {{
    type: 'scatter',
    data: {{
      datasets: [
        {{ label: 'Hap0', data: data0, backgroundColor: 'rgba(41,128,185,0.25)',
           pointRadius: 2, pointHoverRadius: 4 }},
        {{ label: 'Hap1', data: data1, backgroundColor: 'rgba(231,76,60,0.20)',
           pointRadius: 2, pointHoverRadius: 4 }},
        {{ label: 'y = x', data: [{{x:0,y:0}},{{x:maxVal,y:maxVal}}],
           type: 'line', borderColor: '#7f8c8d', borderDash: [4,4],
           borderWidth: 1, pointRadius: 0, fill: false }}
      ]
    }},
    options: {{
      animation: false,
      plugins: {{
        title: {{ display: true, text: title }},
        tooltip: {{ callbacks: {{ label: ctx => {{
          const p = ctx.raw;
          return p.n ? [p.n, `Ref: ${{p.x}} kb  Contig: ${{p.y}} kb`]
                     : `(${{p.x}}, ${{p.y}})`;
        }} }} }}
      }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Reference (kb)' }}, beginAtZero: true }},
        y: {{ title: {{ display: true, text: 'Contig (kb)' }},   beginAtZero: true }}
      }}
    }}
  }});
}}

function scatterLogCfg(canvasId, data0, data1, title) {{
  const pos0 = data0.filter(p => p.x > 0 && p.y > 0);
  const pos1 = data1.filter(p => p.x > 0 && p.y > 0);
  let minV = Infinity, maxV = 0;
  for (const p of pos0.concat(pos1)) {{
    if (p.x < minV) minV = p.x;
    if (p.y < minV) minV = p.y;
    if (p.x > maxV) maxV = p.x;
    if (p.y > maxV) maxV = p.y;
  }}
  new Chart(document.getElementById(canvasId), {{
    type: 'scatter',
    data: {{
      datasets: [
        {{ label: 'Hap0', data: pos0, backgroundColor: 'rgba(41,128,185,0.15)',
           pointRadius: 1.5, pointHoverRadius: 4 }},
        {{ label: 'Hap1', data: pos1, backgroundColor: 'rgba(231,76,60,0.12)',
           pointRadius: 1.5, pointHoverRadius: 4 }},
        {{ label: 'y = x', data: [{{x:minV,y:minV}},{{x:maxV,y:maxV}}],
           type: 'line', borderColor: '#7f8c8d', borderDash: [4,4],
           borderWidth: 1, pointRadius: 0, fill: false }}
      ]
    }},
    options: {{
      animation: false,
      plugins: {{
        title: {{ display: true, text: title }},
        tooltip: {{ callbacks: {{ label: ctx => {{
          const p = ctx.raw;
          return p.n ? [p.n, `Ref: ${{p.x}} kb  Contig: ${{p.y}} kb`]
                     : `(${{p.x}}, ${{p.y}})`;
        }} }} }}
      }},
      scales: {{
        x: {{ type: 'logarithmic', title: {{ display: true, text: 'Reference (kb)' }} }},
        y: {{ type: 'logarithmic', title: {{ display: true, text: 'Contig (kb)' }} }}
      }}
    }}
  }});
}}

function initScatterCharts() {{
  scatterCfg(   'scExon',    {sc_exon0_js}, {sc_exon1_js}, 'Spliced exon length (linear)');
  scatterLogCfg('scExonLog', {sc_exon0_js}, {sc_exon1_js}, 'Spliced exon length (log scale)');
  scatterCfg(   'scSpan',    {sc_span0_js}, {sc_span1_js}, 'Genomic span (linear)');
  scatterLogCfg('scSpanLog', {sc_span0_js}, {sc_span1_js}, 'Genomic span (log scale)');
}}

function initCovChart() {{
  new Chart(document.getElementById('covChart'), {{
    type: 'bar',
    data: {{
      labels: {cov_labels_js},
      datasets: [{{
        label: 'NM genes',
        data: {cov_values_js},
        backgroundColor: '#2980b9',
        borderRadius: 2,
      }}]
    }},
    options: {{
      plugins: {{ legend: {{ display: false }},
                 tooltip: {{ callbacks: {{ title: ctx => 'Coverage ~' + ctx[0].label + '%' }} }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Best exon coverage (%)' }} }},
        y: {{ title: {{ display: true, text: 'NM genes' }}, beginAtZero: true }}
      }}
    }}
  }});
}}
</script>
</body>
</html>
"""

with open(out_path, "w") as f:
    f.write(doc)
print(f"Liftover report written to {out_path}", flush=True)
