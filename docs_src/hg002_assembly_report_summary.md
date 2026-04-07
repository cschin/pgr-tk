# HG002 Post-Assembly Analysis — Presentation Summary

**Date:** 2026-04-07  
**Sample:** HG002 (GIAB Ashkenazi trio — son)  
**Reference:** GRCh38 (no-alt analysis set)  
**Assembly:** hg002v1.1 (maternal hap0 + paternal hap1, T2T-quality)  
**Pipeline:** PGR-TK + Claude Code (mockup report)

---

## Slide 1 — What We Did

> **Goal:** Characterise how well a diploid human assembly matches the GRCh38 reference — at the level of alignments, structural variants, transcript liftover, and clinical variant annotation.

**Pipeline steps:**

```
Assembly (hap0 + hap1)
       │
       ▼
pgr-alnmap            →  alignment chains, M-blocks, SV candidates  (.alndb)
       │
       ├──► pgr-generate-diploid-vcf  →  diploid VCF
       │         │
       │         └──► bcftools + ClinVar  →  clinical annotation
       │
       └──► pgr-liftover-gtf          →  transcript liftover  (_liftover.db)
                   │
                   └──► generate_e2e_report.py  →  HTML report
```

---

## Slide 2 — Transcript Liftover (Hap0)

| Metric | Count | % |
|--------|------:|--:|
| Total transcripts | 207,289 | 100% |
| Fully mapped (≥ 90% exon coverage) | 188,095 | **90.7%** |
| Partial (50–89%) | 370 | 0.2% |
| Multi-contig | 1,062 | 0.5% |
| No hit | 17,599 | **8.5%** |

**Gene-level HQ (≥ 90% coverage):**

| Gene set | HQ genes |
|----------|--------:|
| All gene symbols | 57,268 |
| NM curated mRNA genes | 19,187 / 19,421 (**98.8%**) |

> 98.8% of curated protein-coding gene loci are fully covered in the maternal haplotype —  
> a hallmark of a near-complete T2T-quality assembly.

**No-hit root causes (hap0, NM genes):**

| Cause | Genes |
|-------|------:|
| Annotated only on alt/fix/unplaced locus | ~1,598 |
| Centromere / heterochromatin / repeat-dense | ~1,113 |
| chrY (absent in maternal haplotype) | 57 |

---

## Slide 3 — Structural Variant Candidates

SV candidates are alignment anomalies flagged by `pgr-alnmap` on two coordinate axes.

### Ref view (`.svcnd.bed`) — anomalies on the reference axis

| Type | Hap0 | Hap1 | Meaning |
|------|-----:|-----:|---------|
| **SV** | 7,298 | 7,120 | Size-discordant gap — true structural variant |
| **TG** | 1,542 | 1,725 | Reference gap — unaligned region between contigs |
| **TD** | 5,747 | 5,659 | Contig duplicate — alignment contained in a prior one |
| **TO** |   164 |   167 | Contig overlap — alignment partially overlapping a prior one |
| **Total** | **14,751** | **14,671** | |

### Contig view (`.ctgsv.bed`) — anomalies on the assembly axis

| Type | Meaning |
|------|---------|
| **QG** | Contig gap — unaligned region between consecutive ref alignments |
| **QD** | Ref duplicate — reference block contained in a prior one |
| **QO** | Ref overlap — reference block partially overlapping a prior one |

---

## Slide 4 — SV Gene Impact

SV candidates annotated against GRCh38 gene/exon coordinates (from `_liftover.db`).

### Impact classification

```
exon disrupted   SV:  |══════════════════════|
                 Exon:    |────────|            ← fully inside SV
                 → entire exon deleted; reading frame likely broken

exon partial     SV:  |══════════════|
                 Exon:          |────────────|  ← only partly overlapped
                 → splice site disrupted; altered splicing possible

intronic         SV inside gene body but no exon overlap
intergenic       SV between genes (nearest 5′/3′ neighbors reported)
```

### Summary (hap0 / hap1)

| Impact | Hap0 records | Hap1 records |
|--------|-------------:|-------------:|
| exon disrupted | — | — |
| exon partial | — | — |
| intronic | — | — |
| intergenic | — | — |
| **Unique genes with any exon impact** | **3,139** | **5,108** |

> Full per-gene tables (sorted by impact severity) are available in the HTML report.

---

## Slide 5 — ClinVar Annotation

**Pipeline:** Diploid VCF → GTF annotation → ClinVar 2024 → bcftools annotate

| Category | Variants |
|----------|--------:|
| Pathogenic / Likely pathogenic | see report |
| Uncertain significance | see report |
| Benign / Likely benign | see report |
| Drug response / risk factor | see report |

**Zygosity** (from diploid genotype field):
- **Heterozygous** — one ref allele, one alt allele (carrier or single-copy variant)
- **Homozygous** — both alleles are the same alternate (full penetrance expected)

> Pathogenic variants table in the report includes coordinate, REF→ALT, gene, zygosity, disease name, and ClinVar review status.

---

## Slide 6 — Interactive HTML Report

Two versions generated automatically per run:

| Version | File | Size | Use |
|---------|------|------|-----|
| Full | `e2e_report.html` | ~162 MB | Archive / detailed analysis |
| Lite | `e2e_report_lite.html` | ~7 MB | Sharing / OneDrive / presentations |

**Tabs in the report:**

| Tab | Content |
|-----|---------|
| Step Timings | Wall/CPU time for each pipeline step |
| Alignment Plots | Whole-genome dot-plot per chromosome, hap0 and hap1 |
| ClinVar | Classification breakdown, per-chromosome table, pathogenic variant list |
| Liftover | Transcript status, NM gene funnel, scatter plots (exon length, genomic span) |
| SV Candidates | SV type counts, size distribution, genome-wide ideogram, gene impact |

---

## Slide 7 — Key Takeaways

1. **Near-complete assembly** — 98.8% of curated NM protein-coding genes fully covered (≥ 90% exon bases) in hap0.

2. **~7,200 true SVs per haplotype** detected from alignment anomalies on the reference axis.

3. **3,139 / 5,108 genes** have at least one exon overlapped by an SV candidate in hap0/hap1 — a subset of these may cause protein truncation or altered splicing.

4. **ClinVar pathogenic variants** annotated with zygosity — allows prioritisation of homozygous pathogenic findings.

5. **Self-contained HTML report** — fully offline-capable (CDN scripts inlined), sharable via OneDrive without a server.

---

## Appendix — Tool Versions & Reference Data

| Item | Value |
|------|-------|
| Reference | GCA_000001405.15 GRCh38 no-alt analysis set |
| Assembly (hap0) | hg002v1.1 maternal + MT |
| Assembly (hap1) | hg002v1.1 paternal |
| GTF annotation | hg38.ncbiRefSeq.gtf.gz |
| ClinVar VCF | clinvar.vcf.gz (NCBI) |
| HQ liftover threshold | ≥ 90% exon bases covered |
| Partial liftover threshold | ≥ 50% |
| SV size bins | < 200 bp / 200 bp–1 kb / 1–10 kb / ≥ 10 kb |
| Lite report scatter sample | 5,000 points per chart (random, seed = 42) |
