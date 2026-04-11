#!/usr/bin/env python3
"""
generate_e2e_report.py — produce e2e_report.html from existing output files.

Usage:
    python3 generate_e2e_report.py [--timelog e2e_timings.tsv] [--out e2e_report.html]

All paths are resolved relative to the directory containing this script,
so the script can be run from anywhere.
"""

import sys, html, os, subprocess, shutil, argparse, urllib.request, re

script_dir = os.path.dirname(os.path.abspath(__file__))

# ── CDN script cache (fetch once, reuse) ──────────────────────────────────────
_cdn_cache = {}
def inline_ext_scripts(html_str):
    """Replace all <script src="..."> with fully inlined content for offline use."""
    def _fetch(m):
        url = m.group(1)
        if url not in _cdn_cache:
            try:
                print(f"Fetching {url} ...", flush=True)
                with urllib.request.urlopen(url, timeout=30) as r:
                    _cdn_cache[url] = r.read().decode('utf-8')
            except Exception as e:
                print(f"Warning: could not fetch {url}: {e}", flush=True)
                return m.group(0)   # keep original tag on failure
        return f'<script>{_cdn_cache[url]}</script>'
    return re.sub(r'<script\s+src="([^"]+)"[^>]*></script>', _fetch, html_str)

import json as _json, random as _random

def subsample_scatter_html(html_str, max_pts, seed=42):
    """Replace large scatter-plot JSON arrays with random subsamples of max_pts points.

    Targets arrays of the form [{"x":...,"y":...,"n":"..."},...] embedded in JS.
    Uses a cache so identical arrays (shared between linear/log charts) get the
    same subsample, and a fixed seed for reproducibility.
    """
    rng   = _random.Random(seed)
    cache = {}   # original JSON str -> replacement str
    # Each scatter array element contains only x/y/n keys with no nested [ or ]
    _array_re = re.compile(r'\[[^\[\]]+\]')
    orig_total = [0]   # track pre-subsample size for the note (first large array found)

    def _replace(m):
        orig = m.group(0)
        # Only process arrays that look like scatter data
        if not orig.startswith('[{"x"'):
            return orig
        if orig in cache:
            return cache[orig]
        try:
            arr = _json.loads(orig)
        except Exception:
            cache[orig] = orig
            return orig
        if len(arr) <= max_pts:
            cache[orig] = orig
            return orig
        if orig_total[0] == 0:
            orig_total[0] = len(arr)
        sampled = rng.sample(arr, max_pts)
        result = _json.dumps(sampled, separators=(', ', ': '))
        cache[orig] = result
        return result

    subsampled = _array_re.sub(_replace, html_str)

    # Inject a visible note after each <canvas> tag inside the scatter section
    if orig_total[0] > 0:
        note = (f'<p style="font-size:.8em;color:#e67e22;margin:.2em 0 .8em">'
                f'&#9888; Lite mode: showing {max_pts:,} randomly sampled points '
                f'out of {orig_total[0]:,} total. Run without --lite for full data.</p>')
        subsampled = re.sub(
            r'(<canvas\s[^>]*id="sc[A-Za-z]+[^"]*"[^>]*>)',
            r'\1' + note,
            subsampled)
    return subsampled

ap = argparse.ArgumentParser()
ap.add_argument("--base-dir", default=script_dir,
                help="Directory containing input files (default: script directory)")
ap.add_argument("--timelog", default=None)
ap.add_argument("--out",     default=None)
ap.add_argument("--lite-pts", metavar="N", type=int, default=5000,
                help="Max scatter points in the lite report (default: 5000)")
args = ap.parse_args()

base_dir  = os.path.abspath(args.base_dir)
tsv_path  = args.timelog if args.timelog else os.path.join(base_dir, "e2e_timings.tsv")
html_path = args.out     if args.out     else os.path.join(base_dir, "e2e_report.html")

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
         "-f", "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/GN\t[%GT]\n",
         "-i", 'INFO/CLNSIG!="."',
         clinvar_vcf],
        stderr=subprocess.DEVNULL).decode()

    def gt_zygosity(gt):
        alleles = re.split(r'[|/]', gt)
        non_ref = [a for a in alleles if a not in ('.', '0')]
        if not non_ref:
            return "ref"
        return "hom" if len(set(non_ref)) == 1 and len(non_ref) > 1 else "het"

    records = []
    for line in raw.splitlines():
        parts = line.split('\t')
        if len(parts) == 10:
            records.append(parts)
    total_annotated = len(records)

    from collections import defaultdict
    cat_counts   = defaultdict(int)
    cat_by_chrom = defaultdict(lambda: defaultdict(int))
    zyg_counts   = {"het": 0, "hom": 0}   # overall het/hom tally
    pathogenic_rows = []

    # Detect whether VCF uses "chr1" or bare "1" style
    _sample_chrom = records[0][0] if records else "chr1"
    _chr_prefix   = "chr" if _sample_chrom.startswith("chr") else ""
    CHROMS = [f"{_chr_prefix}{i}" for i in range(1, 23)] + \
             [f"{_chr_prefix}X", f"{_chr_prefix}Y", f"{_chr_prefix}M"]
    # Display label always uses "chr" prefix
    def chrom_label(c):
        return c if c.startswith("chr") else "chr" + c

    def _trunc(s, n=8):
        return html.escape(s[:n]) + "…" if len(s) > n else html.escape(s)

    ZYG_LABEL = {"het": "Heterozygous", "hom": "Homozygous", "ref": "Ref"}
    ZYG_COLOR = {"het": "#2980b9", "hom": "#e74c3c", "ref": "#95a5a6"}

    for chrom, pos, ref, alt, vtype, clnsig, clndn, revstat, gene, gt in records:
        zyg = gt_zygosity(gt)
        if zyg in zyg_counts:
            zyg_counts[zyg] += 1
        for cat_name, predicate in CLNSIG_CATEGORIES:
            if predicate(clnsig):
                cat_counts[cat_name] += 1
                cat_by_chrom[chrom][cat_name] += 1
                if cat_name == "Pathogenic/Likely_pathogenic":
                    diseases = clndn.replace('_', ' ').replace('|', '; ')
                    coord    = f"{chrom_label(chrom)}:{pos}"
                    variant  = f"{_trunc(ref)} &gt; {_trunc(alt)}"
                    pathogenic_rows.append(
                        (coord, variant, gene, zyg, clnsig, diseases[:80],
                         revstat.replace('_', ' ')))
                break

    def card(val, lbl, color="#2980b9"):
        return (f'<div class="cv-card">'
                f'<div class="cv-val" style="color:{color}">{val:,}</div>'
                f'<div class="cv-lbl">{html.escape(lbl)}</div></div>')

    n_path = cat_counts["Pathogenic/Likely_pathogenic"]
    n_vus  = cat_counts["Uncertain_significance"]
    n_ben  = cat_counts["Benign/Likely_benign"]
    n_het  = zyg_counts["het"]
    n_hom  = zyg_counts["hom"]

    cards = (card(total_vars,      "total variants in VCF") +
             card(total_annotated, "ClinVar annotated",           "#8e44ad") +
             card(n_het,           "Heterozygous",                "#2980b9") +
             card(n_hom,           "Homozygous",                  "#e74c3c") +
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
    for coord, variant, gene, zyg, clnsig, diseases, rev in pathogenic_rows:
        zyg_color = ZYG_COLOR.get(zyg, "#95a5a6")
        zyg_label = ZYG_LABEL.get(zyg, zyg)
        path_rows += (f"<tr><td style='white-space:nowrap'>{html.escape(coord)}</td>"
                      f"<td style='font-family:monospace'>{variant}</td>"
                      f"<td><strong>{html.escape(gene)}</strong></td>"
                      f"<td><span style='color:{zyg_color};font-weight:600'>{zyg_label}</span></td>"
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
      <col style="width:12%">
      <col style="width:10%">
      <col style="width:7%">
      <col style="width:10%">
      <col style="width:16%">
      <col style="width:27%">
      <col style="width:18%">
    </colgroup>
    <thead><tr><th>Coordinate</th><th>Variant</th><th>Gene</th>
      <th>Zygosity</th><th>Classification</th><th>Disease</th><th>Review status</th></tr></thead>
    <tbody style="word-break:break-word">{path_rows if path_rows else "<tr><td colspan='7'>None found</td></tr>"}</tbody>
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

# GRCh38 primary chromosome lengths (bp)
CHROM_LENGTHS = {
    "chr1":  248956422, "chr2":  242193529, "chr3":  198295559,
    "chr4":  190214555, "chr5":  181538259, "chr6":  170805979,
    "chr7":  159345973, "chr8":  145138636, "chr9":  138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16":  90338345, "chr17":  83257441, "chr18":  80373285,
    "chr19":  58617616, "chr20":  64444167, "chr21":  46709983,
    "chr22":  50818468, "chrX":  156040895, "chrY":   57227415,
    "chrM":      16569,
}

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

# ── SV gene annotation (Steps 1–4) ────────────────────────────────────────────
import sqlite3, bisect
from collections import defaultdict, Counter

IMPACT_COLOR = {
    "exon disrupted": "#e74c3c",
    "exon partial":   "#e67e22",
    "intronic":       "#3498db",
    "intergenic":     "#95a5a6",
}

def load_gene_exon_index(db_path):
    """Return (gene_index, exon_index) built from liftover.db."""
    gene_index = defaultdict(list)   # chrom → [(g_start, g_end, gene_name, strand)]
    exon_index = defaultdict(list)   # (chrom, gene_name) → [(es, ee)]

    if not os.path.exists(db_path):
        return gene_index, exon_index

    conn = sqlite3.connect(db_path)
    # Gene bodies — one span per gene_name+chrom
    for row in conn.execute("""
        SELECT gene_name, ref_chrom,
               MIN(ref_start) AS g_start, MAX(ref_end) AS g_end, ref_strand
        FROM transcripts
        WHERE gene_name != ''
          AND ref_chrom NOT LIKE '%_alt%'
          AND ref_chrom NOT LIKE '%_fix%'
          AND ref_chrom NOT LIKE '%_random%'
          AND ref_chrom NOT LIKE 'chrUn%'
        GROUP BY gene_name, ref_chrom, ref_strand
    """):
        gene_name, chrom, g_start, g_end, strand = row
        gene_index[chrom].append((g_start, g_end, gene_name, strand))

    # Exons — collapse to unique ranges per gene
    for row in conn.execute("""
        SELECT DISTINCT t.gene_name, t.ref_chrom, e.exon_start, e.exon_end
        FROM exons e JOIN transcripts t ON e.transcript_pk = t.transcript_pk
        WHERE t.gene_name != ''
          AND t.ref_chrom NOT LIKE '%_alt%'
          AND t.ref_chrom NOT LIKE '%_fix%'
          AND t.ref_chrom NOT LIKE '%_random%'
          AND t.ref_chrom NOT LIKE 'chrUn%'
    """):
        gene_name, chrom, es, ee = row
        exon_index[(chrom, gene_name)].append((es, ee))

    conn.close()

    # Sort for bisect queries
    for chrom in gene_index:
        gene_index[chrom].sort()
    for key in exon_index:
        exon_index[key].sort()

    return gene_index, exon_index

def overlapping_genes(chrom, sv_start, sv_end, gene_index):
    entries = gene_index.get(chrom, [])
    # entries sorted by g_start; find insertion point for sv_end
    idx = bisect.bisect_left(entries, (sv_end,))
    hits = []
    for i in range(idx - 1, -1, -1):
        g_start, g_end, gene_name, strand = entries[i]
        if g_end <= sv_start:
            break
        hits.append((gene_name, strand, g_start, g_end))
    return hits

def exon_impact(chrom, gene_name, sv_start, sv_end, exon_index):
    fully, partially = [], []
    for (es, ee) in exon_index.get((chrom, gene_name), []):
        if sv_start <= es and ee <= sv_end:
            fully.append((es, ee))
        elif sv_start < ee and sv_end > es:
            partially.append((es, ee))
    return fully, partially

def nearest_genes(chrom, sv_start, sv_end, gene_index):
    entries = gene_index.get(chrom, [])
    if not entries:
        return None, None
    starts = [e[0] for e in entries]
    idx = bisect.bisect_left(starts, sv_start)

    upstream = None
    for i in range(idx - 1, -1, -1):
        g_start, g_end, gene_name, strand = entries[i]
        if g_end <= sv_start:
            upstream = (gene_name, sv_start - g_end, strand)
            break

    downstream = None
    for i in range(idx, len(entries)):
        g_start, g_end, gene_name, strand = entries[i]
        if g_start >= sv_end:
            downstream = (gene_name, g_start - sv_end, strand)
            break

    return upstream, downstream

def sv_ideogram_html(svcnd0, svcnd1):
    """Return HTML+JS for a genome-wide canvas ideogram of SV positions."""
    import json
    type_to_idx = {"SV": 0, "TG": 1, "TD": 2, "TO": 3}
    sv_data = {"hap0": {}, "hap1": {}}
    for chrom, start, end, sv_type in svcnd0:
        sv_data["hap0"].setdefault(chrom, []).append([(start + end) // 2,
                                                       type_to_idx.get(sv_type, 0)])
    for chrom, start, end, sv_type in svcnd1:
        sv_data["hap1"].setdefault(chrom, []).append([(start + end) // 2,
                                                       type_to_idx.get(sv_type, 0)])
    sv_json        = json.dumps(sv_data)
    chrom_len_json = json.dumps({c: CHROM_LENGTHS[c] for c in CHROMS_ORDER if c in CHROM_LENGTHS})
    chroms_json    = json.dumps(CHROMS_ORDER)
    # colours match SV_TYPE_COLOR: SV, TG, TD, TO
    colors_json    = json.dumps(["#e74c3c", "#3498db", "#2ecc71", "#e67e22"])
    labels_json    = json.dumps([
        "SV — structural variant",
        "TG — ref gap",
        "TD — contig duplicate",
        "TO — contig overlap",
    ])
    n_rows  = len(CHROMS_ORDER)
    row_h   = 28          # px per chromosome row
    legend  = 55          # px for legend area
    top_pad = 10
    canvas_h = top_pad + n_rows * row_h + legend

    return f"""
<p style="font-size:.85em;color:#666">
  Each horizontal bar represents one chromosome (GRCh38 primary assembly, to scale).
  Ticks <b>above</b> a bar = Hap0 SV candidates; ticks <b>below</b> = Hap1.
  Colors indicate SV type.
</p>
<canvas id="sv-ideogram" width="1100" height="{canvas_h}"
        style="display:block;max-width:100%;cursor:crosshair"></canvas>
<div id="sv-ideogram-tip"
     style="position:fixed;background:rgba(0,0,0,.75);color:#fff;padding:4px 8px;
            border-radius:4px;font-size:12px;pointer-events:none;display:none"></div>
<script>
(function(){{
  var svData    = {sv_json};
  var chromLens = {chrom_len_json};
  var chroms    = {chroms_json};
  var colors    = {colors_json};
  var typeNames = ["SV","TG","TD","TO"];
  var typeLabels= {labels_json};

  var LABEL_W = 52, RIGHT_PAD = 12;
  var ROW_H   = {row_h};
  var TOP     = {top_pad};
  var LEGEND_H= {legend};
  var BAR_THICK = 6, MARK_H = 7;
  var MAX_LEN = 248956422;  // chr1

  function getScale(canvas) {{
    return (canvas.width - LABEL_W - RIGHT_PAD) / MAX_LEN;
  }}

  function drawIdeogram() {{
    var canvas = document.getElementById('sv-ideogram');
    if (!canvas) return;
    var ctx = canvas.getContext('2d');
    var W = canvas.width;
    ctx.clearRect(0, 0, W, canvas.height);
    var scale = getScale(canvas);
    ctx.font = '11px sans-serif';

    for (var i = 0; i < chroms.length; i++) {{
      var chrom = chroms[i];
      var clen  = chromLens[chrom] || 0;
      if (!clen) continue;
      var bw = clen * scale;
      var cy = TOP + i * ROW_H + ROW_H / 2;

      // chromosome bar
      ctx.fillStyle = '#c8c8c8';
      ctx.fillRect(LABEL_W, cy - BAR_THICK/2, bw, BAR_THICK);

      // label
      ctx.fillStyle = '#333';
      ctx.textAlign  = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText(chrom, LABEL_W - 4, cy);

      // Hap0 ticks above
      var h0 = (svData.hap0 || {{}})[chrom] || [];
      ctx.lineWidth = 1;
      for (var j = 0; j < h0.length; j++) {{
        ctx.strokeStyle  = colors[h0[j][1]];
        ctx.globalAlpha  = 0.55;
        var xp = LABEL_W + h0[j][0] * scale;
        ctx.beginPath();
        ctx.moveTo(xp, cy - BAR_THICK/2 - MARK_H);
        ctx.lineTo(xp, cy - BAR_THICK/2);
        ctx.stroke();
      }}

      // Hap1 ticks below
      var h1 = (svData.hap1 || {{}})[chrom] || [];
      for (var j = 0; j < h1.length; j++) {{
        ctx.strokeStyle  = colors[h1[j][1]];
        ctx.globalAlpha  = 0.55;
        var xp = LABEL_W + h1[j][0] * scale;
        ctx.beginPath();
        ctx.moveTo(xp, cy + BAR_THICK/2);
        ctx.lineTo(xp, cy + BAR_THICK/2 + MARK_H);
        ctx.stroke();
      }}
      ctx.globalAlpha = 1.0;
    }}

    // legend
    ctx.textAlign    = 'left';
    ctx.textBaseline = 'middle';
    var ly = TOP + chroms.length * ROW_H + 12;
    ctx.fillStyle = '#555';
    ctx.fillText('\\u25b2 Hap0  \\u25bc Hap1', LABEL_W, ly);
    ly += 18;
    var lx = LABEL_W;
    for (var t = 0; t < typeNames.length; t++) {{
      ctx.fillStyle = colors[t];
      ctx.fillRect(lx, ly - 6, 12, 12);
      ctx.fillStyle = '#333';
      ctx.fillText(typeLabels[t], lx + 16, ly);
      lx += 240;
    }}
  }}

  // Tooltip on mousemove
  var canvas = document.getElementById('sv-ideogram');
  var tip    = document.getElementById('sv-ideogram-tip');
  if (canvas && tip) {{
    canvas.addEventListener('mousemove', function(e) {{
      var rect  = canvas.getBoundingClientRect();
      var mx    = (e.clientX - rect.left) * (canvas.width / rect.width);
      var my    = (e.clientY - rect.top)  * (canvas.height / rect.height);
      var scale = getScale(canvas);
      var found = null;
      for (var i = 0; i < chroms.length; i++) {{
        var cy = TOP + i * ROW_H + ROW_H / 2;
        if (Math.abs(my - cy) > ROW_H) continue;
        var chrom = chroms[i];
        var pos   = (mx - LABEL_W) / scale;
        var hap   = my < cy ? 'hap0' : 'hap1';
        var arr   = (svData[hap] || {{}})[chrom] || [];
        var best  = null, bestD = 8 / scale;  // 8px tolerance
        for (var j = 0; j < arr.length; j++) {{
          var d = Math.abs(arr[j][0] - pos);
          if (d < bestD) {{ bestD = d; best = arr[j]; }}
        }}
        if (best) {{
          found = chrom + ':' + best[0].toLocaleString() +
                  '  ' + typeNames[best[1]] + '  (' + hap + ')';
        }}
        break;
      }}
      if (found) {{
        tip.textContent = found;
        tip.style.display = 'block';
        tip.style.left = (e.clientX + 12) + 'px';
        tip.style.top  = (e.clientY - 20) + 'px';
      }} else {{
        tip.style.display = 'none';
      }}
    }});
    canvas.addEventListener('mouseleave', function() {{
      tip.style.display = 'none';
    }});
  }}

  // Lazy draw: called by showSubTab when this panel becomes visible
  window._svIdeogramDraw = drawIdeogram;
}})();
</script>
"""

def annotate_svcnd(records, gene_index, exon_index):
    ann = []
    for chrom, start, end, sv_type in records:
        hits = overlapping_genes(chrom, start, end, gene_index)
        if hits:
            for gene_name, strand, g_start, g_end in hits:
                fully, partially = exon_impact(chrom, gene_name, start, end, exon_index)
                impact = ("exon disrupted" if fully else
                          "exon partial"   if partially else
                          "intronic")
                ann.append(dict(
                    chrom=chrom, start=start, end=end, sv_type=sv_type,
                    context="genic", gene=gene_name, strand=strand,
                    impact=impact,
                    n_exon_full=len(fully), n_exon_partial=len(partially),
                    upstream=None, downstream=None,
                ))
        else:
            up, dn = nearest_genes(chrom, start, end, gene_index)
            ann.append(dict(
                chrom=chrom, start=start, end=end, sv_type=sv_type,
                context="intergenic", gene=None, strand=None,
                impact="intergenic",
                n_exon_full=0, n_exon_partial=0,
                upstream=up, downstream=dn,
            ))
    return ann

def sv_gene_impact_html(ann, hap_label):
    if not ann:
        return "<p><em>No data.</em></p>"

    # ── summary counts ──
    total     = len(ann)
    n_genic   = sum(1 for r in ann if r['context'] == 'genic')
    n_interg  = sum(1 for r in ann if r['context'] == 'intergenic')
    n_exon_d  = sum(1 for r in ann if r['impact'] == 'exon disrupted')
    n_exon_p  = sum(1 for r in ann if r['impact'] == 'exon partial')
    n_intron  = sum(1 for r in ann if r['impact'] == 'intronic')

    # unique genes with any exon impact (potential protein truncation / splicing change)
    genes_exon_d   = {r['gene'] for r in ann if r['impact'] == 'exon disrupted' and r['gene']}
    genes_exon_p   = {r['gene'] for r in ann if r['impact'] == 'exon partial'   and r['gene']}
    genes_exon_any = genes_exon_d | genes_exon_p

    def scard(val, lbl, color="#2980b9", subtitle=""):
        sub = f'<div style="font-size:.72em;color:#aaa;margin-top:.1em">{subtitle}</div>' if subtitle else ""
        return (f'<div style="background:#fff;border:1px solid #ddd;border-radius:8px;'
                f'padding:.6em 1em;text-align:center;min-width:120px">'
                f'<div style="font-size:1.6em;font-weight:700;color:{color}">{val:,}</div>'
                f'<div style="font-size:.8em;color:#7f8c8d">{lbl}</div>{sub}</div>')

    cards = (scard(total,    "total records") +
             scard(n_genic,  "genic",            "#e67e22") +
             scard(n_interg, "intergenic",        "#95a5a6") +
             scard(n_exon_d, "exon disrupted",    "#e74c3c") +
             scard(n_exon_p, "exon partial",      "#e67e22") +
             scard(n_intron, "intronic",          "#3498db") +
             scard(len(genes_exon_any), "genes with exon impact", "#8e44ad",
                   subtitle=f"{len(genes_exon_d)} fully disrupted · {len(genes_exon_p)} partial"))

    # ── impact × SV type cross-tab ──
    sv_types_present = sorted({r['sv_type'] for r in ann})
    impacts = ["exon disrupted", "exon partial", "intronic", "intergenic"]
    cross = defaultdict(lambda: defaultdict(int))
    for r in ann:
        cross[r['impact']][r['sv_type']] += 1

    cross_head = "".join(f"<th>{t}</th>" for t in sv_types_present)
    cross_rows = ""
    for imp in impacts:
        color = IMPACT_COLOR.get(imp, "#222")
        row = f"<tr><td style='font-weight:600;color:{color}'>{imp}</td>"
        for t in sv_types_present:
            row += f"<td class='num'>{cross[imp].get(t, 0):,}</td>"
        cross_rows += row + "</tr>"

    # ── top 20 affected genes ──
    gene_hits = Counter(r['gene'] for r in ann if r['gene'])
    gene_impacts = defaultdict(set)
    for r in ann:
        if r['gene']:
            gene_impacts[r['gene']].add(r['impact'])

    top_gene_rows = ""
    for gene_name, count in gene_hits.most_common(20):
        imps = ", ".join(sorted(gene_impacts[gene_name]))
        top_gene_rows += (f"<tr><td><strong>{html.escape(gene_name)}</strong></td>"
                          f"<td class='num'>{count}</td>"
                          f"<td style='font-size:.85em'>{html.escape(imps)}</td></tr>")

    # ── genes with exon impact (potential protein truncation / splicing change) ──
    # Collect per-gene SV count and worst impact for sorting
    exon_gene_info = {}   # gene -> {sv_count, has_full, has_partial, chroms}
    for r in ann:
        if r['impact'] not in ('exon disrupted', 'exon partial'):
            continue
        g = r['gene']
        if g not in exon_gene_info:
            exon_gene_info[g] = {'sv_count': 0, 'has_full': False, 'has_partial': False,
                                 'chroms': set(), 'strand': r['strand']}
        exon_gene_info[g]['sv_count'] += 1
        if r['impact'] == 'exon disrupted':
            exon_gene_info[g]['has_full'] = True
        else:
            exon_gene_info[g]['has_partial'] = True
        exon_gene_info[g]['chroms'].add(r['chrom'])

    # Sort: fully-disrupted first, then by SV count desc
    exon_gene_sorted = sorted(
        exon_gene_info.items(),
        key=lambda kv: (not kv[1]['has_full'], -kv[1]['sv_count']))

    exon_gene_rows = ""
    for gene_name, info in exon_gene_sorted:
        worst = "exon disrupted" if info['has_full'] else "exon partial"
        wcolor = "#e74c3c" if info['has_full'] else "#e67e22"
        note = ("protein truncation / loss likely" if info['has_full']
                else "altered splicing possible")
        chroms_str = ", ".join(sorted(info['chroms'],
                                      key=lambda c: CHROMS_ORDER.index(c)
                                      if c in CHROMS_ORDER else 99))
        exon_gene_rows += (
            f"<tr>"
            f"<td><strong>{html.escape(gene_name)}</strong> ({info['strand']})</td>"
            f"<td>{html.escape(chroms_str)}</td>"
            f"<td class='num'>{info['sv_count']}</td>"
            f"<td style='color:{wcolor};font-weight:600'>{worst}</td>"
            f"<td style='font-size:.85em;color:#555'>{note}</td>"
            f"</tr>")

    # ── full annotation table (first 200 rows) ──
    def fmt_neighbor(n):
        if n is None:
            return "—"
        gene, dist, strand = n
        return f"{html.escape(gene)} ({strand}, {dist:,} bp)"

    ann_rows = ""
    for r in sorted(ann, key=lambda x: (x['chrom'], x['start']))[:200]:
        size = r['end'] - r['start']
        impact_color = IMPACT_COLOR.get(r['impact'], "#222")
        if r['context'] == 'genic':
            gene_cell  = f"<strong>{html.escape(r['gene'])}</strong> ({r['strand']})"
            n_exon_cell = f"{r['n_exon_full']} full / {r['n_exon_partial']} partial"
            up_cell    = "—"
            dn_cell    = "—"
        else:
            gene_cell  = "<span style='color:#95a5a6'>—</span>"
            n_exon_cell = "—"
            up_cell    = fmt_neighbor(r['upstream'])
            dn_cell    = fmt_neighbor(r['downstream'])

        sv_color = SV_TYPE_COLOR.get(r['sv_type'], "#222")
        ann_rows += (
            f"<tr>"
            f"<td style='white-space:nowrap'>{html.escape(r['chrom'])}:{r['start']:,}</td>"
            f"<td class='num'>{size:,}</td>"
            f"<td><span style='color:{sv_color};font-weight:bold'>{r['sv_type']}</span></td>"
            f"<td style='color:{impact_color};font-weight:600'>{r['impact']}</td>"
            f"<td>{gene_cell}</td>"
            f"<td class='num' style='font-size:.85em'>{n_exon_cell}</td>"
            f"<td style='font-size:.8em'>{up_cell}</td>"
            f"<td style='font-size:.8em'>{dn_cell}</td>"
            f"</tr>"
        )
    note = "" if len(ann) <= 200 else f"showing first 200 of {len(ann):,} records"

    return f"""
<h4>{html.escape(hap_label)}</h4>

<h4 style="color:#8e44ad;margin:.4em 0 .4em">Genes with exon impact — potential protein truncation or splicing change
  <span style="font-weight:400;font-size:.8em;color:#7f8c8d">({len(genes_exon_any)} unique genes:
  {len(genes_exon_d)} with fully disrupted exon · {len(genes_exon_p)} with partial exon overlap)</span>
</h4>
<div style="background:#f8f8f8;border:1px solid #ddd;border-radius:6px;
            padding:.7em 1.2em;margin:.3em 0 .8em;font-size:.85em;color:#444">
  <div style="margin-bottom:.5em"><strong>How impact is classified:</strong></div>
  <table style="border:none;background:none;margin:0;width:auto">
    <tr style="border:none">
      <td style="border:none;padding:.2em 1.2em .2em 0;white-space:nowrap;vertical-align:top">
        <span style="color:#e74c3c;font-weight:600">exon disrupted</span>
      </td>
      <td style="border:none;padding:.2em 0;font-family:monospace;white-space:pre;line-height:1.4">SV:    |&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;|
Exon:      |&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;|          &#x2190; fully inside SV
&#x21d2; entire exon deleted/replaced; reading frame likely broken</td>
    </tr>
    <tr style="border:none">
      <td style="border:none;padding:.4em 1.2em .2em 0;white-space:nowrap;vertical-align:top">
        <span style="color:#e67e22;font-weight:600">exon partial</span>
      </td>
      <td style="border:none;padding:.4em 0;font-family:monospace;white-space:pre;line-height:1.4">SV:    |&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;&#x2550;|
Exon:              |&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;&#x2500;|  &#x2190; only partly overlapped
&#x21d2; splice site may be disrupted; partial coding-sequence change</td>
    </tr>
  </table>
</div>
<h4 style="color:#34495e;margin:1.2em 0 .4em">Summary</h4>
<div style="display:flex;flex-wrap:wrap;gap:.8em;margin-bottom:1.2em">{cards}</div>

<details>
  <summary style="cursor:pointer;font-weight:600;color:#8e44ad;margin:.4em 0 .6em;
                  list-style:none;display:flex;align-items:center;gap:.4em">
    <span id="exon-gene-toggle-{hap_label}">&#9654;</span>
    Show / hide gene table ({len(exon_gene_info)} genes)
  </summary>
  <table>
    <thead><tr>
      <th>Gene (strand)</th><th>Chromosome(s)</th><th>SV records</th>
      <th>Worst impact</th><th>Consequence</th>
    </tr></thead>
    <tbody>{exon_gene_rows}</tbody>
  </table>
</details>

<details>
  <summary style="cursor:pointer;font-weight:600;color:#34495e;margin:.6em 0 .4em">
    Impact × SV type
  </summary>
  <table>
    <thead><tr><th>Impact</th>{cross_head}</tr></thead>
    <tbody>{cross_rows}</tbody>
  </table>
</details>

<details>
  <summary style="cursor:pointer;font-weight:600;color:#34495e;margin:.6em 0 .4em">
    Top 20 affected genes
  </summary>
  <table>
    <thead><tr><th>Gene</th><th>SV records</th><th>Impact types</th></tr></thead>
    <tbody>{top_gene_rows}</tbody>
  </table>
</details>

<details>
  <summary style="cursor:pointer;font-weight:600;color:#34495e;margin:.6em 0 .4em">
    Per-SV annotation (first 200){" — " + note.strip() if note.strip() else ""}
  </summary>
  <table style="table-layout:fixed;width:100%">
    <colgroup>
      <col style="width:13%"><col style="width:7%"><col style="width:5%">
      <col style="width:11%"><col style="width:13%"><col style="width:11%">
      <col style="width:20%"><col style="width:20%">
    </colgroup>
    <thead><tr>
      <th>Position</th><th>Size</th><th>Type</th><th>Impact</th>
      <th>Gene</th><th>Exons</th><th>5′ neighbor</th><th>3′ neighbor</th>
    </tr></thead>
    <tbody style="word-break:break-word">{ann_rows}</tbody>
  </table>
</details>
"""

# Load gene/exon index (from hap0 DB; both haps share the same reference annotation)
_db0 = os.path.join(base_dir, "hg002_hap0_liftover.db")
_db1 = os.path.join(base_dir, "hg002_hap1_liftover.db")
_have_db = os.path.exists(_db0) or os.path.exists(_db1)

sv_ideogram_section = sv_ideogram_html(svcnd0, svcnd1) if (svcnd0 or svcnd1) else ""

sv_gene_section0 = ""
sv_gene_section1 = ""
if (svcnd0 or svcnd1) and _have_db:
    print("Building gene/exon index ...", flush=True)
    _ref_db = _db0 if os.path.exists(_db0) else _db1
    gene_index, exon_index = load_gene_exon_index(_ref_db)
    print(f"  {sum(len(v) for v in gene_index.values()):,} gene spans, "
          f"{sum(len(v) for v in exon_index.values()):,} exon ranges loaded", flush=True)

    print("Annotating SV candidates ...", flush=True)
    ann0 = annotate_svcnd(svcnd0, gene_index, exon_index)
    ann1 = annotate_svcnd(svcnd1, gene_index, exon_index)
    print(f"  hap0: {len(ann0):,} records  hap1: {len(ann1):,} records", flush=True)

    sv_gene_section0 = sv_gene_impact_html(ann0, "Hap0")
    sv_gene_section1 = sv_gene_impact_html(ann1, "Hap1")
elif svcnd0 or svcnd1:
    _no_db_msg = "<p><em>Liftover DB not found — run get_tx_seqs.sh first to enable gene impact annotation.</em></p>"
    sv_gene_section0 = _no_db_msg
    sv_gene_section1 = _no_db_msg

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

<div class="subtabbar">
  <button class="active" onclick="showSubTab('svcnd','sv-types')"   id="subbtn-svcnd-sv-types">SV Types</button>
  <button onclick="showSubTab('svcnd','sv-sizechr')"  id="subbtn-svcnd-sv-sizechr">Size &amp; Chromosome</button>
  <button onclick="showSubTab('svcnd','sv-genome')"   id="subbtn-svcnd-sv-genome">Genome View</button>
  <button onclick="showSubTab('svcnd','sv-gene0')"    id="subbtn-svcnd-sv-gene0">Gene Impact — Hap0</button>
  <button onclick="showSubTab('svcnd','sv-gene1')"    id="subbtn-svcnd-sv-gene1">Gene Impact — Hap1</button>
</div>

<!-- SV Types (ref + contig combined) -->
<div class="subpanel active" id="subpanel-svcnd-sv-types">
  <h3>Ref view — SV candidate type breakdown (svcnd.bed)</h3>
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
  <h3>Contig view — SV candidate type breakdown (ctgsv.bed)</h3>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:2em">
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
</div>

<!-- Size & Chromosome -->
<div class="subpanel" id="subpanel-svcnd-sv-sizechr">
  <h3>Size distribution &amp; Per-chromosome breakdown (ref view)</h3>
  <div style="display:grid;grid-template-columns:auto 1fr;gap:2em;align-items:start">
    <div>
      <h4>Size distribution</h4>
      <table style="width:auto">
        <thead><tr><th>Size bin</th><th>Hap0</th><th>Hap1</th></tr></thead>
        <tbody>{sv_size_table}</tbody>
      </table>
    </div>
    <div>
      <h4>Per-chromosome (primary chromosomes)</h4>
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
</div>

<!-- Genome View -->
<div class="subpanel" id="subpanel-svcnd-sv-genome">
  <h3>Genome-wide SV candidate distribution</h3>
  {sv_ideogram_section}
</div>

<!-- Gene Impact Hap0 -->
<div class="subpanel" id="subpanel-svcnd-sv-gene0">
  <h3>Gene Impact Annotation — Hap0</h3>
  {sv_gene_section0}
</div>

<!-- Gene Impact Hap1 -->
<div class="subpanel" id="subpanel-svcnd-sv-gene1">
  <h3>Gene Impact Annotation — Hap1</h3>
  {sv_gene_section1}
</div>
"""
else:
    sv_section = ("<h2>SV Candidate Summary</h2>"
                  "<p><em>hg002_hap*.svcnd.bed not found — run pgr-alnmap first.</em></p>")

# ── Liftover summary (queried directly from SQLite databases) ─────────────────
import sqlite3 as _sqlite3

def _liftover_summary_html(db_path, hap_label):
    """Generate an HTML summary section from a liftover SQLite database."""
    if not os.path.exists(db_path):
        return f'<p><em>{hap_label}: {db_path} not found.</em></p>'
    try:
        con = _sqlite3.connect(db_path)
        cur = con.cursor()

        # Overall gene stats
        cur.execute("""
            SELECT SUM(total), SUM(single_full), SUM(single_partial),
                   SUM(multi_contig), SUM(multi_location), SUM(no_hit)
            FROM gene_summary
        """)
        total_g, sf, sp, mc, ml, nh = cur.fetchone()
        total_g = total_g or 0

        # Transcript status breakdown
        cur.execute("""
            SELECT status, COUNT(*) FROM transcript_summary GROUP BY status ORDER BY COUNT(*) DESC
        """)
        tx_rows = cur.fetchall()
        total_tx = sum(c for _, c in tx_rows)

        # Per-chromosome gene summary (top chroms by total transcripts)
        cur.execute("""
            SELECT ref_chrom,
                   SUM(total)          AS total,
                   SUM(single_full)    AS sf,
                   SUM(single_partial) AS sp,
                   SUM(multi_contig)   AS mc,
                   SUM(multi_location) AS ml,
                   SUM(no_hit)         AS nh
            FROM gene_summary
            GROUP BY ref_chrom
            ORDER BY total DESC
        """)
        chrom_rows = cur.fetchall()
        con.close()
    except Exception as e:
        return f'<p><em>Error reading {db_path}: {html.escape(str(e))}</em></p>'

    def _pct(n, d):
        return f"{100*n/d:.1f}%" if d else "—"

    def _bar(n, d, color="#2980b9"):
        pct = 100*n/d if d else 0
        return (f'<div style="height:14px;background:{color};border-radius:3px;'
                f'min-width:2px;width:{pct:.1f}%"></div>')

    # Gene-level summary cards
    cards = ""
    for label, val, color in [
        ("Total transcripts", total_g, "#2c3e50"),
        ("Single full",       sf,      "#27ae60"),
        ("Single partial",    sp,      "#f39c12"),
        ("Multi-contig",      mc,      "#8e44ad"),
        ("Multi-location",    ml,      "#e67e22"),
        ("No hit",            nh,      "#e74c3c"),
    ]:
        cards += (f'<div class="cv-card">'
                  f'<div class="cv-val" style="color:{color}">{val:,}</div>'
                  f'<div class="cv-lbl">{label}<br>({_pct(val, total_g)})</div>'
                  f'</div>\n')

    # Transcript status table
    tx_table_rows = ""
    for status, cnt in tx_rows:
        tx_table_rows += (f'<tr><td class="name">{html.escape(status)}</td>'
                          f'<td class="num">{cnt:,}</td>'
                          f'<td class="num">{_pct(cnt, total_tx)}</td>'
                          f'<td class="bar-cell">{_bar(cnt, total_tx)}</td></tr>\n')

    # Per-chromosome table
    chrom_table_rows = ""
    for chrom, tot, sf_c, sp_c, mc_c, ml_c, nh_c in chrom_rows:
        chrom_table_rows += (
            f'<tr><td class="name">{html.escape(chrom)}</td>'
            f'<td class="num">{tot:,}</td>'
            f'<td class="num">{sf_c:,} ({_pct(sf_c, tot)})</td>'
            f'<td class="num">{sp_c:,}</td>'
            f'<td class="num">{mc_c:,}</td>'
            f'<td class="num">{ml_c:,}</td>'
            f'<td class="num">{nh_c:,} ({_pct(nh_c, tot)})</td>'
            f'</tr>\n')

    return f"""
<h3>{html.escape(hap_label)}</h3>
<p style="color:#555;font-size:.9em">Database: <code>{html.escape(db_path)}</code></p>
<div class="cv-cards">{cards}</div>

<h4>Transcript status breakdown ({total_tx:,} transcripts)</h4>
<table>
  <thead><tr><th>Status</th><th>Count</th><th>Percent</th><th>Bar</th></tr></thead>
  <tbody>{tx_table_rows}</tbody>
</table>

<h4>Per-chromosome gene summary</h4>
<table>
  <thead><tr>
    <th>Chrom</th><th>Total</th><th>Single full</th>
    <th>Single partial</th><th>Multi-contig</th><th>Multi-location</th><th>No hit</th>
  </tr></thead>
  <tbody>{chrom_table_rows}</tbody>
</table>
"""

_db0 = os.path.join(base_dir, "hg002_hap0_liftover.db")
_db1 = os.path.join(base_dir, "hg002_hap1_liftover.db")
_lo0 = _liftover_summary_html(_db0, "Haplotype 0 (maternal)")
_lo1 = _liftover_summary_html(_db1, "Haplotype 1 (paternal)")

liftover_inline = f"""
<h2>GTF Liftover Summary</h2>
<p>Transcripts lifted from GRCh38 RefSeq annotation onto HG002 haplotype contigs.</p>
{_lo0}
<hr>
{_lo1}
"""

# ── Determine which tabs to show ──────────────────────────────────────────────
tabs = []
tabs.append(("timings", "Step Timings"))
if plots_html:
    tabs.append(("plots", "Alignment Plots"))
if clinvar_section:
    tabs.append(("clinvar", "ClinVar"))
if os.path.exists(_db0) or os.path.exists(_db1):
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
<div class="tabpanel" id="tab-liftover">
  {liftover_inline}
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
  if (parentId === 'svcnd' && subId === 'sv-genome' && window._svIdeogramDraw) {{
    window._svIdeogramDraw();
    window._svIdeogramDraw = null;
  }}
}}
</script>
</body>
</html>
"""

# ── Derive lite filename: e2e_report.html → e2e_report_lite.html ─────────────
_base, _ext = os.path.splitext(html_path)
lite_path = _base + "_lite" + _ext
lite_pts  = args.lite_pts

# Build lite doc: subsample scatter arrays, add a banner at the top
print(f"Building lite report (≤{lite_pts:,} scatter points) ...", flush=True)
_lite_doc = subsample_scatter_html(doc, lite_pts)

# Cross-link banners (inserted just after <body> opening)
_full_banner = (
    f'<div style="background:#eaf4ff;border-bottom:1px solid #b3d4f5;padding:.5em 1.5em;'
    f'font-size:.85em">Full report &mdash; scatter plots contain all data points. '
    f'<a href="{os.path.basename(lite_path)}">Switch to lite version</a> '
    f'(faster to load, {lite_pts:,} points per chart).</div>')

_lite_banner = (
    f'<div style="background:#fff8e1;border-bottom:1px solid #ffe082;padding:.5em 1.5em;'
    f'font-size:.85em">&#9888; Lite report &mdash; scatter plots comparing transcripts only show {lite_pts:,} randomly '
    f'sampled points. '
    f'<a href="{os.path.basename(html_path)}">Switch to full version</a> '
    f'(all data, larger file).</div>')

doc       = doc.replace('<body>', '<body>' + _full_banner, 1)
_lite_doc = _lite_doc.replace('<body>', '<body>' + _lite_banner, 1)

with open(html_path, "w") as f:
    f.write(doc)
print(f"Report written to {html_path}")

with open(lite_path, "w") as f:
    f.write(_lite_doc)
print(f"Lite report written to {lite_path}")
