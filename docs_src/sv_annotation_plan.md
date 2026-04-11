# SV Candidate Annotation Plan

**Date:** 2026-04-07  
**Branch:** `aln_map_improvment`  
**Target file:** `test_genomes/generate_e2e_report.py`

---

## Goals

1. **Gene/exon impact** — for each SV candidate that overlaps a gene body, report which genes are hit and whether any exons are fully or partially disrupted.
2. **Nearest neighbors** — for intergenic SVs, report the closest upstream (5') and downstream (3') gene and the distance to each.
3. *(Future)* Repeat element overlap.

---

## Input Data

| File | Coordinate system | Notes |
|------|-------------------|-------|
| `hg002_hap{0,1}.svcnd.bed` | GRCh38, CHROM = `GRCh38#0#chr1` | ref-view SVs |
| `hg002_hap{0,1}.ctgsv.bed` | Contig coordinates | contig-view; skip for gene annotation |
| `hg002_hap{0,1}_liftover.db` | `transcripts.ref_chrom` uses `chr` prefix | SQLite |

**Coordinate normalisation:** strip the `GRCh38#0#` PanSN prefix from `svcnd.bed` CHROM → `chr1`, `chr2`, … to match the liftover DB.

---

## Implementation Plan

### Step 1 — Build in-memory gene/exon interval index

Load once from the liftover DB (hap0; genes are reference-based so both haps share the same annotation):

```sql
-- Gene bodies (one row per unique gene span)
SELECT DISTINCT gene_name, ref_chrom, MIN(ref_start) AS g_start, MAX(ref_end) AS g_end, ref_strand
FROM transcripts
WHERE gene_name != ''
  AND ref_chrom NOT LIKE '%_alt%'
  AND ref_chrom NOT LIKE '%_fix%'
  AND ref_chrom NOT LIKE '%_random%'
  AND ref_chrom NOT LIKE 'chrUn%'
GROUP BY gene_name, ref_chrom, ref_strand;

-- Exons (for impact classification)
SELECT t.gene_name, t.ref_chrom, e.exon_start, e.exon_end
FROM exons e JOIN transcripts t ON e.transcript_pk = t.transcript_pk
WHERE t.gene_name != ''
  AND t.ref_chrom NOT LIKE '%_alt%'
  ...
```

**Data structure — simple sorted list + bisect (no external deps):**

```python
from collections import defaultdict
import bisect

# gene_index[chrom] = sorted list of (start, end, gene_name, strand)
gene_index = defaultdict(list)   # sorted by start after loading

# exon_index[(chrom, gene_name)] = sorted list of (exon_start, exon_end)
exon_index = defaultdict(list)
```

A full interval-overlap query over ~60 k gene bodies with `bisect` is fast enough:

```python
def overlapping_genes(chrom, sv_start, sv_end, gene_index):
    entries = gene_index[chrom]
    # entries sorted by g_start; find first where g_start < sv_end
    idx = bisect.bisect_left(entries, (sv_end,))
    hits = []
    for i in range(idx - 1, -1, -1):
        g_start, g_end, gene_name, strand = entries[i]
        if g_end <= sv_start:
            break
        hits.append((gene_name, strand, g_start, g_end))
    return hits
```

### Step 2 — Classify exon impact

For each overlapping gene, check its exons:

```python
def exon_impact(chrom, gene_name, sv_start, sv_end, exon_index):
    fully, partially = [], []
    for (es, ee) in exon_index[(chrom, gene_name)]:
        if sv_start <= es and ee <= sv_end:
            fully.append((es, ee))          # SV fully contains exon
        elif sv_start < ee and sv_end > es:
            partially.append((es, ee))      # partial overlap
    return fully, partially
```

**Impact label:**

| Condition | Label |
|-----------|-------|
| ≥1 exon fully contained | `exon disrupted` |
| ≥1 exon partially overlapped | `exon partial` |
| Gene body only (no exon overlap) | `intronic` |

### Step 3 — Nearest neighbor for intergenic SVs

For SVs with no gene overlap, binary-search the sorted gene list per chromosome:

```python
def nearest_genes(chrom, sv_start, sv_end, gene_index):
    entries = gene_index[chrom]   # sorted by g_start
    starts = [e[0] for e in entries]

    # upstream (5'): last gene whose end <= sv_start
    idx = bisect.bisect_left(starts, sv_start)
    upstream = None
    for i in range(idx - 1, -1, -1):
        g_start, g_end, gene_name, strand = entries[i]
        if g_end <= sv_start:
            upstream = (gene_name, sv_start - g_end, strand)
            break

    # downstream (3'): first gene whose start >= sv_end
    downstream = None
    for i in range(idx, len(entries)):
        g_start, g_end, gene_name, strand = entries[i]
        if g_start >= sv_end:
            downstream = (gene_name, g_start - sv_end, strand)
            break

    return upstream, downstream
```

### Step 4 — Annotate each SV record

```python
SV_ANN = []   # list of dicts, one per SV candidate

for chrom, start, end, sv_type in svcnd0:
    hits = overlapping_genes(chrom, start, end, gene_index)
    if hits:
        for gene_name, strand, g_start, g_end in hits:
            fully, partially = exon_impact(chrom, gene_name, start, end, exon_index)
            impact = ("exon disrupted" if fully else
                      "exon partial"   if partially else
                      "intronic")
            SV_ANN.append(dict(
                chrom=chrom, start=start, end=end, sv_type=sv_type,
                context="genic", gene=gene_name, strand=strand,
                impact=impact,
                n_exon_full=len(fully), n_exon_partial=len(partially),
                upstream=None, downstream=None
            ))
    else:
        up, dn = nearest_genes(chrom, start, end, gene_index)
        SV_ANN.append(dict(
            chrom=chrom, start=start, end=end, sv_type=sv_type,
            context="intergenic", gene=None, strand=None, impact=None,
            n_exon_full=0, n_exon_partial=0,
            upstream=up, downstream=dn
        ))
```

### Step 5 — Summary statistics

From `SV_ANN` compute:

| Metric | How |
|--------|-----|
| Total genic | `sum(r['context']=='genic')` |
| Total intergenic | `sum(r['context']=='intergenic')` |
| Exon-disrupting | `sum(r['impact']=='exon disrupted')` |
| Exon-partial | `sum(r['impact']=='exon partial')` |
| Intronic | `sum(r['impact']=='intronic')` |
| Top genes by SV count | `Counter(r['gene'] for r in SV_ANN if r['gene'])` |

### Step 6 — HTML output (new subtab inside SV Candidates tab)

Add a third subtab **"Gene Impact"** alongside the existing Ref view / Contig view subtabs:

```
[Ref view]  [Contig view]  [Gene Impact]
```

**Gene Impact subtab layout:**

1. **Summary cards** — Total SVs / Genic / Intergenic / Exon-disrupting / Exon-partial / Intronic
2. **Impact breakdown table** — rows: SV type × impact category (cross-tabulation)
3. **Top affected genes table** — gene name | strand | # SVs | impact types
4. **Full annotation table** (collapsible / paginated) — chrom | pos | size | type | context | gene | impact | upstream neighbor | downstream neighbor

---

## Performance Estimates

| Operation | Scale | Expected time |
|-----------|-------|---------------|
| Load gene index from SQLite | ~60 k genes | <1 s |
| Load exon index from SQLite | ~500 k exons | ~2 s |
| Annotate 14 k SVs (bisect) | 14 k × O(log 60k) | <1 s |
| HTML generation | — | negligible |

Total added runtime: **~3 s** per run.

---

## Future: Repeat Element Overlap (Phase 3)

- Source: UCSC RepeatMasker BED (`hg38.repeatmasker.bed.gz`) — downloadable, ~170 MB
- Same interval-index approach; add `repeat_family` field to `SV_ANN`
- Report: fraction of SV sequence covered by repeats; breakdown by repeat class (SINE, LINE, LTR, Simple_repeat, …)
- Implementation identical to gene overlap — load once into `repeat_index[chrom]`, query with `bisect`

---

## Files to Modify

| File | Change |
|------|--------|
| `generate_e2e_report.py` | Add Steps 1–6 after existing SV parsing; add Gene Impact subtab HTML |
| `get_tx_seqs.sh` | No change needed (liftover DB already built) |

No new dependencies — pure Python stdlib + SQLite (already used elsewhere in the script).
