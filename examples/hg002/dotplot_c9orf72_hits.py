#!/usr/bin/env python3
"""
dotplot_c9orf72_hits.py — pairwise sequence dotplot of C9orf72 hit contigs
using the pgr-tk Python interface (shimmer-pair matches).

Usage:
    python dotplot_c9orf72_hits.py [hit_fasta] [output_png]

Defaults:
    hit_fasta  = example_output/c9orf72_hits.000.fa
    output_png = example_output/c9orf72_dotplot.png
"""

import sys
import gzip
import os

# pgr-tk lives one level up in the repo; add its site-packages to sys.path
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_VENV_SP = os.path.join(_REPO, "pgr-tk", ".venv", "lib")
for _d in os.listdir(_VENV_SP):
    sp = os.path.join(_VENV_SP, _d, "site-packages")
    if os.path.isdir(sp) and sp not in sys.path:
        sys.path.insert(0, sp)

import pgrtk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HIT_FA  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SCRIPT_DIR, "example_output", "c9orf72_hits.000.fa")
OUT_PNG = sys.argv[2] if len(sys.argv) > 2 else os.path.join(SCRIPT_DIR, "example_output", "c9orf72_dotplot.png")

# ---------------------------------------------------------------------------
# Read FASTA
# ---------------------------------------------------------------------------
def read_fasta(path):
    seqs = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        name, buf = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs.append((name, "".join(buf)))
    return seqs

seqs = read_fasta(HIT_FA)
print(f"Loaded {len(seqs)} sequences from {HIT_FA}")
for name, seq in seqs:
    print(f"  {name}  ({len(seq):,} bp)")

# ---------------------------------------------------------------------------
# Build SeqIndexDB from all hit sequences
# ---------------------------------------------------------------------------
W, K, R, MIN_SPAN = 48, 56, 4, 12

sdb = pgrtk.SeqIndexDB()
seq_list = [(name, list(seq.upper().encode())) for name, seq in seqs]
sdb.load_from_seq_list(seq_list, source=None, w=W, k=K, r=R, min_span=MIN_SPAN)

seq_info = sdb.seq_info          # id -> (ctg_name, src, length)
id_to_name = {sid: info[0] for sid, info in seq_info.items()}
id_to_len  = {sid: info[2] for sid, info in seq_info.items()}
name_to_id = {v: k for k, v in id_to_name.items()}

n = len(seqs)
short_labels = []
for name, _ in seqs:
    # keep last component of the pangenome-style name for readability
    parts = name.split("::")
    short_labels.append(parts[-1] if len(parts) > 1 else name)

# ---------------------------------------------------------------------------
# For each sequence, query the database and collect (qpos, tpos, orient)
#   get_match_positions_with_fragment returns: {target_sid: [(qpos, tpos, orient), ...]}
# ---------------------------------------------------------------------------
dots = {}   # (qi, ti) -> list of (qpos, tpos, orient)
for qi, (qname, qseq) in enumerate(seqs):
    qbytes = list(qseq.upper().encode())
    hits = sdb.get_match_positions_with_fragment(qbytes)
    for target_sid, positions in hits.items():
        ti = list(id_to_name.keys()).index(target_sid)
        key = (qi, ti)
        dots.setdefault(key, []).extend(positions)

# ---------------------------------------------------------------------------
# Plot: n×n grid of dot plots
# ---------------------------------------------------------------------------
FIG_SIZE = max(4 * n, 8)
fig = plt.figure(figsize=(FIG_SIZE, FIG_SIZE))
gs  = gridspec.GridSpec(n, n, hspace=0.05, wspace=0.05)

for qi in range(n):
    qname, qseq = seqs[qi]
    qsid = name_to_id[qname]
    qlen = len(qseq)

    for ti in range(n):
        tname, tseq = seqs[ti]
        tsid = name_to_id[tname]
        tlen = len(tseq)

        ax = fig.add_subplot(gs[ti, qi])  # row=target, col=query

        key = (qi, ti)
        if key in dots and dots[key]:
            pts = dots[key]
            fwd = [(p[0], p[1]) for p in pts if p[2] == 0]
            rev = [(p[0], p[1]) for p in pts if p[2] == 1]
            if fwd:
                xs, ys = zip(*fwd)
                ax.scatter(xs, ys, s=0.5, c="steelblue", linewidths=0, rasterized=True)
            if rev:
                xs, ys = zip(*rev)
                ax.scatter(xs, ys, s=0.5, c="tomato",    linewidths=0, rasterized=True)

        ax.set_xlim(0, qlen)
        ax.set_ylim(0, tlen)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        # Labels on outer edges only
        if ti == n - 1:
            ax.set_xlabel(short_labels[qi], fontsize=7, rotation=15, ha="right")
        if qi == 0:
            ax.set_ylabel(short_labels[ti], fontsize=7, rotation=0, ha="right", va="center")

# Legend
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue", markersize=6, label="forward"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="tomato",    markersize=6, label="reverse"),
]
fig.legend(handles=legend_handles, loc="lower right", fontsize=8, framealpha=0.8)
fig.suptitle("C9orf72 hit contigs — pairwise shimmer dotplot", fontsize=11, y=1.01)

plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"\nDotplot written: {OUT_PNG}")
