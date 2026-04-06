#!/usr/bin/env python3
"""
generate_e2e_report.py — produce e2e_report.html from existing output files.

Usage:
    python3 generate_e2e_report.py [--timelog e2e_timings.tsv] [--out e2e_report.html]

All paths are resolved relative to the directory containing this script,
so the script can be run from anywhere.
"""

import sys, html, os, subprocess, shutil, argparse

script_dir = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument("--timelog", default=os.path.join(script_dir, "e2e_timings.tsv"))
ap.add_argument("--out",     default=os.path.join(script_dir, "e2e_report.html"))
args = ap.parse_args()

tsv_path  = args.timelog
html_path = args.out
base_dir  = script_dir

# ── Step timings ──────────────────────────────────────────────────────────────
rows = []
if os.path.exists(tsv_path):
    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) == 6:
                rows.append(parts)

max_wall   = max((float(r[1]) for r in rows), default=1)
total_wall = sum(float(r[1]) for r in rows)
peak_rss   = max((float(r[4]) for r in rows), default=0)

timing_rows = ""
for r in rows:
    name, wall, user, sys_, rss, cpu = r
    pct = min(100, float(wall) / max_wall * 100)
    try:
        rss_mb = f"{float(rss)/1024:.1f}"
    except ValueError:
        rss_mb = rss
    timing_rows += f"""
      <tr>
        <td class="name">{html.escape(name)}</td>
        <td class="bar-cell"><div class="bar" style="width:{pct:.1f}%"></div></td>
        <td class="num">{float(wall):.1f}s</td>
        <td class="num">{float(user):.1f}s</td>
        <td class="num">{float(sys_):.1f}s</td>
        <td class="num">{rss_mb} MB</td>
        <td class="num">{html.escape(cpu)}</td>
      </tr>"""

# ── Alignment plots (embedded iframes) ───────────────────────────────────────
PLOT_DEFS = [
    ("hg002_hap0_aln_plot", "Haplotype 0", "plot-hap0"),
    ("hg002_hap1_aln_plot", "Haplotype 1", "plot-hap1"),
]

def embed_plot(prefix):
    fname = os.path.join(base_dir, prefix + ".html")
    if not os.path.exists(fname):
        return ""
    with open(fname) as f:
        content = f.read()
    escaped = content.replace("&", "&amp;").replace('"', "&quot;")
    return (f'<iframe srcdoc="{escaped}"'
            f' style="width:100%;height:600px;border:1px solid #ddd;'
            f'border-radius:6px;background:#fff;" sandbox="allow-scripts"></iframe>')

# Build subtab bar and panels
_plot_buttons = ""
_plot_panels  = ""
_first_plot   = True
for prefix, label, tid in PLOT_DEFS:
    iframe = embed_plot(prefix)
    if not iframe:
        continue
    active_btn   = ' class="active"'   if _first_plot else ""
    active_panel = ' class="subpanel active"' if _first_plot else ' class="subpanel"'
    _plot_buttons += (f'<button{active_btn} '
                      f'onclick="showSubTab(\'plots\',\'{tid}\')" '
                      f'id="subbtn-plots-{tid}">{html.escape(label)}</button>\n')
    _plot_panels  += f'<div{active_panel} id="subpanel-plots-{tid}">{iframe}</div>\n'
    _first_plot    = False

if _plot_buttons:
    plots_html = (f'<div class="subtabbar">{_plot_buttons}</div>'
                  f'{_plot_panels}')
else:
    plots_html = ""

# ── ClinVar annotation summary ────────────────────────────────────────────────
CLNSIG_CATEGORIES = [
    ("Pathogenic/Likely_pathogenic",
     lambda s: ("Pathogenic" in s or "Likely_pathogenic" in s)
               and "Benign" not in s and "Likely_benign" not in s),
    ("Uncertain_significance",
     lambda s: "Uncertain_significance" in s),
    ("Benign/Likely_benign",
     lambda s: ("Benign" in s or "Likely_benign" in s)),
    ("Drug response / association / risk factor",
     lambda s: any(k in s for k in ("drug_response", "association", "risk_factor",
                                    "protective", "Likely_risk_allele"))),
    ("Other / no classification",
     lambda s: True),  # catch-all
]
CAT_COLOR = {
    "Pathogenic/Likely_pathogenic":              "#e74c3c",
    "Uncertain_significance":                    "#e67e22",
    "Benign/Likely_benign":                      "#27ae60",
    "Drug response / association / risk factor": "#2980b9",
    "Other / no classification":                 "#95a5a6",
}

clinvar_vcf = os.path.join(base_dir, "hg002.annotated.sorted.clinvar.vcf.gz")
bcftools    = shutil.which("bcftools")

clinvar_section = ""
if bcftools and os.path.exists(clinvar_vcf):
    print("Querying ClinVar VCF ...", flush=True)

    total_vars = int(subprocess.check_output(
        [bcftools, "view", "--no-header", clinvar_vcf],
        stderr=subprocess.DEVNULL).decode().count('\n'))

    raw = subprocess.check_output(
        [bcftools, "query",
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/GN\n",
         "-i", 'INFO/CLNSIG!="."',
         clinvar_vcf],
        stderr=subprocess.DEVNULL).decode()

    records = []
    for line in raw.splitlines():
        parts = line.split('\t')
        if len(parts) == 9:
            records.append(parts)
    total_annotated = len(records)

    from collections import defaultdict
    cat_counts   = defaultdict(int)
    cat_by_chrom = defaultdict(lambda: defaultdict(int))
    pathogenic_rows = []

    # Detect whether VCF uses "chr1" or bare "1" style
    _sample_chrom = records[0][0] if records else "chr1"
    _chr_prefix   = "chr" if _sample_chrom.startswith("chr") else ""
    CHROMS = [f"{_chr_prefix}{i}" for i in range(1, 23)] + \
             [f"{_chr_prefix}X", f"{_chr_prefix}Y", f"{_chr_prefix}M"]
    # Display label always uses "chr" prefix
    def chrom_label(c):
        return c if c.startswith("chr") else "chr" + c

    for chrom, pos, ref, alt, vtype, clnsig, clndn, revstat, gene in records:
        for cat_name, predicate in CLNSIG_CATEGORIES:
            if predicate(clnsig):
                cat_counts[cat_name] += 1
                cat_by_chrom[chrom][cat_name] += 1
                if cat_name == "Pathogenic/Likely_pathogenic":
                    diseases = clndn.replace('_', ' ').replace('|', '; ')
                    coord    = f"{chrom_label(chrom)}:{pos}"
                    def _trunc(s, n=8):
                        return html.escape(s[:n]) + "…" if len(s) > n else html.escape(s)
                    variant  = f"{_trunc(ref)} &gt; {_trunc(alt)}"
                    pathogenic_rows.append(
                        (coord, variant, gene, clnsig, diseases[:80],
                         revstat.replace('_', ' ')))
                break

    def card(val, lbl, color="#2980b9"):
        return (f'<div class="cv-card">'
                f'<div class="cv-val" style="color:{color}">{val:,}</div>'
                f'<div class="cv-lbl">{html.escape(lbl)}</div></div>')

    n_path = cat_counts["Pathogenic/Likely_pathogenic"]
    n_vus  = cat_counts["Uncertain_significance"]
    n_ben  = cat_counts["Benign/Likely_benign"]

    cards = (card(total_vars,      "total variants in VCF") +
             card(total_annotated, "ClinVar annotated",           "#8e44ad") +
             card(n_ben,           "Benign / Likely benign",      "#27ae60") +
             card(n_vus,           "Uncertain significance",      "#e67e22") +
             card(n_path,          "Pathogenic / Likely pathogenic", "#e74c3c"))

    max_cat = max(cat_counts.values(), default=1)
    breakdown_rows = ""
    for cat_name, _ in CLNSIG_CATEGORIES:
        n = cat_counts.get(cat_name, 0)
        color = CAT_COLOR[cat_name]
        w = int(n / max_cat * 180)
        pct_s = f"{n/total_annotated*100:.1f}%" if total_annotated else "0%"
        breakdown_rows += (
            f"<tr><td>{html.escape(cat_name)}</td>"
            f"<td class='num'>{n:,}</td>"
            f"<td class='num'>{pct_s}</td>"
            f"<td><div style='display:inline-block;width:{w}px;height:14px;"
            f"background:{color};border-radius:2px;vertical-align:middle'></div></td></tr>")

    chrom_rows = ""
    for chrom in CHROMS:
        d = cat_by_chrom.get(chrom, {})
        if not d:
            continue
        p = d.get("Pathogenic/Likely_pathogenic", 0)
        v = d.get("Uncertain_significance", 0)
        b = d.get("Benign/Likely_benign", 0)
        total_chr = sum(d.values())
        chrom_rows += (f"<tr><td>{html.escape(chrom_label(chrom))}</td>"
                       f"<td class='num'>{total_chr:,}</td>"
                       f"<td class='num' style='color:#e74c3c'>{p}</td>"
                       f"<td class='num' style='color:#e67e22'>{v}</td>"
                       f"<td class='num' style='color:#27ae60'>{b}</td></tr>")

    path_rows = ""
    for coord, variant, gene, clnsig, diseases, rev in pathogenic_rows:
        path_rows += (f"<tr><td style='white-space:nowrap'>{html.escape(coord)}</td>"
                      f"<td style='font-family:monospace'>{variant}</td>"
                      f"<td><strong>{html.escape(gene)}</strong></td>"
                      f"<td>{html.escape(clnsig)}</td>"
                      f"<td>{html.escape(diseases)}</td>"
                      f"<td style='font-size:.8em;color:#7f8c8d'>{html.escape(rev)}</td></tr>")

    clinvar_section = f"""
<h2>ClinVar Annotation Summary</h2>
<div class="cv-cards">{cards}</div>

<div class="subtabbar">
  <button class="active" onclick="showSubTab('clinvar','cv-breakdown')" id="subbtn-clinvar-cv-breakdown">Classification</button>
  <button onclick="showSubTab('clinvar','cv-chrom')"     id="subbtn-clinvar-cv-chrom">Per-chromosome</button>
  <button onclick="showSubTab('clinvar','cv-pathogenic')" id="subbtn-clinvar-cv-pathogenic">Pathogenic / LP ({n_path})</button>
</div>

<div class="subpanel active" id="subpanel-clinvar-cv-breakdown">
  <table>
    <thead><tr><th>ClinVar classification</th><th>Variants</th><th>%</th><th></th></tr></thead>
    <tbody>{breakdown_rows}</tbody>
  </table>
</div>

<div class="subpanel" id="subpanel-clinvar-cv-chrom">
  <table style="width:auto;min-width:340px">
    <thead><tr><th>Chrom</th><th>Total</th>
      <th style="color:#e74c3c">Path/LP</th>
      <th style="color:#e67e22">VUS</th>
      <th style="color:#27ae60">Benign</th></tr></thead>
    <tbody>{chrom_rows}</tbody>
  </table>
</div>

<div class="subpanel" id="subpanel-clinvar-cv-pathogenic">
  <table style="width:100%;table-layout:fixed">
    <colgroup>
      <col style="width:14%">
      <col style="width:12%">
      <col style="width:8%">
      <col style="width:18%">
      <col style="width:30%">
      <col style="width:18%">
    </colgroup>
    <thead><tr><th>Coordinate</th><th>Variant</th><th>Gene</th>
      <th>Classification</th><th>Disease</th><th>Review status</th></tr></thead>
    <tbody style="word-break:break-word">{path_rows if path_rows else "<tr><td colspan='6'>None found</td></tr>"}</tbody>
  </table>
</div>
"""
elif not bcftools:
    clinvar_section = ("<h2>ClinVar Annotation Summary</h2>"
                       "<p><em>bcftools not found in PATH — skipping.</em></p>")
else:
    clinvar_section = ("<h2>ClinVar Annotation Summary</h2>"
                       "<p><em>hg002.annotated.sorted.clinvar.vcf.gz not found — "
                       "run all pipeline steps first.</em></p>")

# ── SV candidates ─────────────────────────────────────────────────────────────
SV_TYPE_LABEL = {
    "SV": "Structural variant",
    "TG": "Ref gap",
    "TD": "Contig duplicate (ref view)",
    "TO": "Contig overlap (ref view)",
    "QG": "Contig gap",
    "QD": "Ref duplicate (contig view)",
    "QO": "Ref overlap (contig view)",
}
SV_TYPE_COLOR = {
    "SV": "#e74c3c",
    "TG": "#e67e22",
    "TD": "#3498db",
    "TO": "#9b59b6",
    "QG": "#16a085",
    "QD": "#2980b9",
    "QO": "#8e44ad",
}
SIZE_BINS = [
    ("<200 bp",    lambda s: s < 200),
    ("200 bp–1 kb", lambda s: 200 <= s < 1_000),
    ("1 kb–10 kb",  lambda s: 1_000 <= s < 10_000),
    ("≥10 kb",     lambda s: s >= 10_000),
]
CHROMS_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

def parse_svcnd_bed(path):
    """Return list of (chrom, start, end, sv_type) stripping PanSN prefix."""
    records = []
    if not os.path.exists(path):
        return records
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0].split('#')[-1]   # strip GRCh38#0# prefix
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            sv_type = parts[3][:2]             # first two chars: SV/TG/TD/TO
            records.append((chrom, start, end, sv_type))
    return records

def svcnd_summary(records):
    from collections import defaultdict
    by_type    = defaultdict(int)
    by_size    = {label: 0 for label, _ in SIZE_BINS}
    by_chrom   = defaultdict(lambda: defaultdict(int))
    for chrom, start, end, sv_type in records:
        size = end - start
        by_type[sv_type] += 1
        for label, pred in SIZE_BINS:
            if pred(size):
                by_size[label] += 1
                break
        by_chrom[chrom][sv_type] += 1
    return by_type, by_size, by_chrom

svcnd0 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap0.svcnd.bed"))
svcnd1 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap1.svcnd.bed"))
ctgsv0 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap0.ctgsv.bed"))
ctgsv1 = parse_svcnd_bed(os.path.join(base_dir, "hg002_hap1.ctgsv.bed"))
sv_type0, sv_size0, sv_chrom0 = svcnd_summary(svcnd0)
sv_type1, sv_size1, sv_chrom1 = svcnd_summary(svcnd1)
q_type0, q_size0, _ = svcnd_summary(ctgsv0)
q_type1, q_size1, _ = svcnd_summary(ctgsv1)

def sv_type_rows(by_type, total, codes):
    rows = ""
    for code in codes:
        n = by_type.get(code, 0)
        color = SV_TYPE_COLOR[code]
        label = SV_TYPE_LABEL[code]
        pct = n / total * 100 if total else 0
        w = int(pct / 100 * 160)
        rows += (f"<tr><td><span style='font-family:monospace;font-weight:bold;"
                 f"color:{color}'>{code}</span></td>"
                 f"<td>{html.escape(label)}</td>"
                 f"<td class='num'>{n:,}</td>"
                 f"<td class='num'>{pct:.1f}%</td>"
                 f"<td><div style='display:inline-block;width:{w}px;height:14px;"
                 f"background:{color};border-radius:2px;vertical-align:middle'></div></td></tr>")
    return rows

def sv_size_rows(by_size0, by_size1):
    rows = ""
    for label, _ in SIZE_BINS:
        n0 = by_size0.get(label, 0)
        n1 = by_size1.get(label, 0)
        rows += (f"<tr><td>{html.escape(label)}</td>"
                 f"<td class='num'>{n0:,}</td>"
                 f"<td class='num'>{n1:,}</td></tr>")
    return rows

def sv_chrom_rows(by_chrom0, by_chrom1):
    rows = ""
    seen = set()
    for chrom in CHROMS_ORDER:
        d0 = by_chrom0.get(chrom, {})
        d1 = by_chrom1.get(chrom, {})
        if not d0 and not d1:
            continue
        seen.add(chrom)
        t0 = sum(d0.values()); t1 = sum(d1.values())
        sv0 = d0.get("SV", 0); sv1 = d1.get("SV", 0)
        rows += (f"<tr><td>{html.escape(chrom)}</td>"
                 f"<td class='num'>{t0:,}</td>"
                 f"<td class='num' style='color:#e74c3c'>{sv0}</td>"
                 f"<td class='num'>{t1:,}</td>"
                 f"<td class='num' style='color:#e74c3c'>{sv1}</td></tr>")
    return rows

total_sv0 = len(svcnd0)
total_sv1 = len(svcnd1)
total_q0  = len(ctgsv0)
total_q1  = len(ctgsv1)

sv_type_rows0 = sv_type_rows(sv_type0, total_sv0, ["SV", "TG", "TD", "TO"])
sv_type_rows1 = sv_type_rows(sv_type1, total_sv1, ["SV", "TG", "TD", "TO"])
q_type_rows0  = sv_type_rows(q_type0,  total_q0,  ["QG", "QD", "QO"])
q_type_rows1  = sv_type_rows(q_type1,  total_q1,  ["QG", "QD", "QO"])
sv_size_table = sv_size_rows(sv_size0, sv_size1)
sv_chrom_table = sv_chrom_rows(sv_chrom0, sv_chrom1)

if svcnd0 or svcnd1:
    sv_section = f"""
<h2>SV Candidate Summary</h2>
<p style="font-size:.85em;color:#7f8c8d">
  <b>Ref view</b> (.svcnd.bed) — gaps and anomalies on the reference coordinate axis:<br>
  <b>SV</b> = structural variant (size-discordant gap) &nbsp;·&nbsp;
  <b>TG</b> = ref gap (unaligned region between consecutive contigs) &nbsp;·&nbsp;
  <b>TD</b> = contig duplicate (contig alignment contained within a previous one) &nbsp;·&nbsp;
  <b>TO</b> = contig overlap (contig alignment partially overlapping a previous one)<br><br>
  <b>Contig view</b> (.ctgsv.bed) — gaps and anomalies on the contig coordinate axis:<br>
  <b>QG</b> = contig gap (unaligned region between consecutive ref alignments) &nbsp;·&nbsp;
  <b>QD</b> = ref duplicate (ref block contained within a previous one) &nbsp;·&nbsp;
  <b>QO</b> = ref overlap (ref block partially overlapping a previous one)
</p>

<h3>Ref view (svcnd.bed)</h3>
<div style="display:grid;grid-template-columns:1fr 1fr;gap:2em;margin-bottom:1.5em">
  <div>
    <h4>Hap0 — {total_sv0:,} records</h4>
    <table>
      <thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th><th></th></tr></thead>
      <tbody>{sv_type_rows0}</tbody>
    </table>
  </div>
  <div>
    <h4>Hap1 — {total_sv1:,} records</h4>
    <table>
      <thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th><th></th></tr></thead>
      <tbody>{sv_type_rows1}</tbody>
    </table>
  </div>
</div>

<h3>Contig view (ctgsv.bed)</h3>
<div style="display:grid;grid-template-columns:1fr 1fr;gap:2em;margin-bottom:1.5em">
  <div>
    <h4>Hap0 — {total_q0:,} records</h4>
    <table>
      <thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th><th></th></tr></thead>
      <tbody>{q_type_rows0}</tbody>
    </table>
  </div>
  <div>
    <h4>Hap1 — {total_q1:,} records</h4>
    <table>
      <thead><tr><th>Code</th><th>Type</th><th>Count</th><th>%</th><th></th></tr></thead>
      <tbody>{q_type_rows1}</tbody>
    </table>
  </div>
</div>

<div style="display:grid;grid-template-columns:auto 1fr;gap:2em;align-items:start">
  <div>
    <h3>Size distribution (ref view)</h3>
    <table style="width:auto">
      <thead><tr><th>Size bin</th><th>Hap0</th><th>Hap1</th></tr></thead>
      <tbody>{sv_size_table}</tbody>
    </table>
  </div>
  <div>
    <h3>Per-chromosome (primary, ref view)</h3>
    <table>
      <thead><tr>
        <th>Chrom</th>
        <th>Hap0 total</th><th style="color:#e74c3c">Hap0 SV</th>
        <th>Hap1 total</th><th style="color:#e74c3c">Hap1 SV</th>
      </tr></thead>
      <tbody>{sv_chrom_table}</tbody>
    </table>
  </div>
</div>
"""
else:
    sv_section = ("<h2>SV Candidate Summary</h2>"
                  "<p><em>hg002_hap*.svcnd.bed not found — run pgr-alnmap first.</em></p>")

# ── Liftover report (embedded iframe) ────────────────────────────────────────
liftover_iframe = ""
_liftover_path = os.path.join(base_dir, "liftover_report.html")
if os.path.exists(_liftover_path):
    with open(_liftover_path) as _f:
        _liftover_content = _f.read()
    _liftover_escaped = _liftover_content.replace("&", "&amp;").replace('"', "&quot;")
    liftover_iframe = (f'<iframe srcdoc="{_liftover_escaped}"'
                       f' style="width:100%;height:calc(100vh - 120px);border:none;" '
                       f'sandbox="allow-scripts"></iframe>')

# ── Determine which tabs to show ──────────────────────────────────────────────
tabs = []
tabs.append(("timings", "Step Timings"))
if plots_html:
    tabs.append(("plots", "Alignment Plots"))
if clinvar_section:
    tabs.append(("clinvar", "ClinVar"))
if liftover_iframe:
    tabs.append(("liftover", "Liftover"))
if svcnd0 or svcnd1 or ctgsv0 or ctgsv1:
    tabs.append(("svcnd", "SV Candidates"))

tab_buttons = ""
for i, (tid, tlabel) in enumerate(tabs):
    active = ' class="active"' if i == 0 else ""
    tab_buttons += f'<button{active} onclick="showTab(\'{tid}\')" id="btn-{tid}">{html.escape(tlabel)}</button>\n'

# ── Compose page ──────────────────────────────────────────────────────────────
doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>DGI Post Assembly Report (mockup, developed with PGR-TK + Claude Code)</title>
<style>
  body  {{ font-family: sans-serif; margin: 0; background: #f9f9f9; color: #222; }}
  h1    {{ color: #2c3e50; margin: 0; padding: .6em 2em; background: #fff;
           border-bottom: 1px solid #ddd; font-size: 1.4em; }}
  h2    {{ color: #2c3e50; }}
  h3    {{ color: #34495e; margin-bottom: .4em; }}
  /* ── tab bar ── */
  .tabbar {{ display: flex; gap: 0; background: #2c3e50; padding: 0 2em; }}
  .tabbar button {{
    background: transparent; border: none; color: #bdc3c7;
    padding: .7em 1.4em; font-size: .95em; cursor: pointer;
    border-bottom: 3px solid transparent; transition: color .15s;
  }}
  .tabbar button:hover  {{ color: #ecf0f1; }}
  .tabbar button.active {{ color: #fff; border-bottom-color: #2980b9; font-weight: 600; }}
  /* ── tab panels ── */
  .tabpanel {{ display: none; padding: 1.5em 2em; }}
  .tabpanel.active {{ display: block; }}
  /* ── summary strip ── */
  .summary {{ background: #fff; border: 1px solid #ddd; border-radius: 6px;
              padding: .8em 1.4em; display: inline-block; margin-bottom: 1.5em; }}
  .summary span {{ font-weight: bold; color: #2980b9; }}
  /* ── tables ── */
  table {{ border-collapse: collapse; width: 100%; background: #fff;
           border: 1px solid #ddd; border-radius: 6px; overflow: hidden;
           margin-bottom: 1.5em; }}
  th {{ background: #2c3e50; color: #fff; padding: .6em 1em; text-align: left; }}
  td {{ padding: .5em 1em; border-bottom: 1px solid #eee; }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: #f0f4f8; }}
  .name {{ white-space: nowrap; font-weight: 500; }}
  .num  {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .bar-cell {{ width: 35%; padding: 0 1em; }}
  .bar  {{ height: 18px; background: #2980b9; border-radius: 3px; min-width: 2px; }}
  /* ── subtab bar (used inside Alignment Plots tab) ── */
  .subtabbar {{ display: flex; gap: 0; background: #ecf0f1;
                border-bottom: 2px solid #bdc3c7; margin-bottom: 1em; }}
  .subtabbar button {{
    background: transparent; border: none; color: #7f8c8d;
    padding: .5em 1.2em; font-size: .9em; cursor: pointer;
    border-bottom: 3px solid transparent; margin-bottom: -2px;
    transition: color .15s;
  }}
  .subtabbar button:hover  {{ color: #2c3e50; }}
  .subtabbar button.active {{ color: #2980b9; border-bottom-color: #2980b9; font-weight: 600; }}
  .subpanel {{ display: none; }}
  .subpanel.active {{ display: block; }}
  /* ── clinvar cards ── */
  .cv-cards {{ display: flex; flex-wrap: wrap; gap: 1em; margin-bottom: 1.5em; }}
  .cv-card  {{ background: #fff; border: 1px solid #ddd; border-radius: 8px;
               padding: .8em 1.2em; min-width: 160px; text-align: center; }}
  .cv-val   {{ font-size: 1.8em; font-weight: 700; }}
  .cv-lbl   {{ font-size: .8em; color: #7f8c8d; margin-top: .2em; }}
</style>
</head>
<body>
<h1>DGI Post Assembly Report (mockup, developed with PGR-TK + Claude Code)</h1>
<div class="tabbar">
{tab_buttons}</div>

<!-- Tab: Step Timings -->
<div class="tabpanel active" id="tab-timings">
  <div class="summary">
    Total wall time: <span>{total_wall:.1f}s</span> &nbsp;|&nbsp;
    Peak RSS: <span>{peak_rss/1024:.1f} MB</span> &nbsp;|&nbsp;
    Steps: <span>{len(rows)}</span>
  </div>
  <h2>Step timings</h2>
  <table>
    <thead>
      <tr>
        <th>Step</th><th>Wall time</th><th>Wall (s)</th>
        <th>User (s)</th><th>Sys (s)</th><th>Max RSS</th><th>CPU%</th>
      </tr>
    </thead>
    <tbody>{timing_rows}
    </tbody>
  </table>
</div>

<!-- Tab: Alignment Plots -->
<div class="tabpanel" id="tab-plots">
  <h2>Alignment plots</h2>
  {plots_html}
</div>

<!-- Tab: ClinVar -->
<div class="tabpanel" id="tab-clinvar">
  {clinvar_section}
</div>

<!-- Tab: Liftover -->
<div class="tabpanel" id="tab-liftover" style="padding:0">
  {liftover_iframe}
</div>

<!-- Tab: SV Candidates -->
<div class="tabpanel" id="tab-svcnd">
  {sv_section}
</div>

<script>
function showTab(id) {{
  document.querySelectorAll('.tabpanel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tabbar button').forEach(b => b.classList.remove('active'));
  document.getElementById('tab-' + id).classList.add('active');
  document.getElementById('btn-' + id).classList.add('active');
}}
function showSubTab(parentId, subId) {{
  const panel = document.getElementById('tab-' + parentId);
  panel.querySelectorAll('.subpanel').forEach(p => p.classList.remove('active'));
  panel.querySelectorAll('.subtabbar button').forEach(b => b.classList.remove('active'));
  document.getElementById('subpanel-' + parentId + '-' + subId).classList.add('active');
  document.getElementById('subbtn-' + parentId + '-' + subId).classList.add('active');
}}
</script>
</body>
</html>
"""

with open(html_path, "w") as f:
    f.write(doc)
print(f"Report written to {html_path}")
