# Liftover Database (`_liftover.db`) — Reference Guide

**Date:** 2026-04-07  
**Produced by:** `pgr align liftover-gtf`  
**Example files:** `hg002_hap0_liftover.db`, `hg002_hap1_liftover.db`  
**Annotation source:** `hg38.ncbiRefSeq.gtf.gz`

---

## Overview

One SQLite database is created per haplotype.  It records:

1. Every transcript from the input GTF (reference coordinates)
2. How well each transcript's exons map onto the assembled contigs
3. Summary classifications per transcript and per gene
4. Optionally: the spliced sequences (reference and contig)

**Row counts (hap0, HG002):**

| Table | Rows |
|-------|-----:|
| `transcripts` | 207,289 |
| `exons` | 2,214,932 |
| `liftover` | 192,988 |
| `transcript_summary` | 207,289 |
| `gene_summary` | 59,251 |
| `ref_sequences` | 207,289 |
| `contig_sequences` | 192,988 |

---

## Schema

### `transcripts`

One row per transcript from the GTF.  The primary key (`transcript_pk`) is the join key used by all other tables.

```sql
CREATE TABLE transcripts (
    transcript_pk  INTEGER PRIMARY KEY,
    transcript_id  TEXT    NOT NULL UNIQUE,  -- e.g. "NM_000001.3"
    gene_id        TEXT    NOT NULL,          -- NCBI gene ID
    gene_name      TEXT    NOT NULL,          -- gene symbol, e.g. "BRCA1"
    ref_chrom      TEXT    NOT NULL,          -- "chr1" .. "chrM"
    ref_start      INTEGER NOT NULL,          -- gene body start, 0-based, inclusive
    ref_end        INTEGER NOT NULL,          -- gene body end, exclusive
    ref_strand     TEXT    NOT NULL,          -- "+" or "-"
    exon_count     INTEGER NOT NULL           -- number of exons in this transcript
);
CREATE INDEX idx_transcripts_chrom ON transcripts(ref_chrom, ref_start, ref_end);
```

**Transcript ID prefixes present (hap0):**

| Prefix | Type | Count |
|--------|------|------:|
| `NM_` | NCBI curated mRNA | 71,867 |
| `XM_` | NCBI predicted mRNA | 66,380 |
| `XR_` | NCBI predicted ncRNA | 28,636 |
| `NR_` | NCBI curated ncRNA | 21,515 |
| other | pseudogenes, snRNA, rRNA, LOC*, OR*, … | ~18,891 |

**Coverage by chromosome:** chr1 (18,362 transcripts) through chrM (37 transcripts).

---

### `exons`

One row per exon–transcript pair.  A transcript with `exon_count = N` has exactly N rows here.

```sql
CREATE TABLE exons (
    exon_pk        INTEGER PRIMARY KEY,
    transcript_pk  INTEGER NOT NULL REFERENCES transcripts ON DELETE CASCADE,
    exon_start     INTEGER NOT NULL,   -- 0-based, inclusive (GRCh38 coordinates)
    exon_end       INTEGER NOT NULL    -- exclusive
);
CREATE INDEX idx_exons_tx ON exons(transcript_pk);
```

**Statistics (hap0):**  2,214,932 exon records.  Average 10.7 exons per transcript; range 1–363.

---

### `liftover`

One row per transcript–contig alignment hit.  A transcript can have zero (no hit) or more rows (multiple contig hits above the minimum coverage threshold of 50%).

```sql
CREATE TABLE liftover (
    liftover_pk    INTEGER PRIMARY KEY,
    transcript_pk  INTEGER NOT NULL REFERENCES transcripts ON DELETE CASCADE,
    contig         TEXT    NOT NULL,    -- query contig name, e.g. "HG002#1#chr1"
    contig_start   INTEGER NOT NULL,   -- mapped span start on contig, 0-based
    contig_end     INTEGER NOT NULL,   -- mapped span end on contig
    strand         TEXT    NOT NULL,   -- "+" or "-" relative to contig
    coverage_pct   REAL    NOT NULL,   -- exon bases covered / total exon bases × 100
    block_count    INTEGER NOT NULL,   -- number of mapped exon blocks
    block_sizes    TEXT    NOT NULL,   -- comma-separated block lengths (bp)
    block_starts   TEXT    NOT NULL    -- comma-separated offsets from contig_start
);
CREATE INDEX idx_liftover_tx     ON liftover(transcript_pk);
CREATE INDEX idx_liftover_contig ON liftover(contig, contig_start, contig_end);
```

**Key fields:**

| Field | Notes |
|-------|-------|
| `coverage_pct` | Range 50.0–200.0 in practice (>100 possible when assembly duplicates a region). Average ~99.8 for mapped transcripts. |
| `block_sizes` | Mirrors the exon structure on the contig; parse with `SUM(int(x) for x in block_sizes.split(',') if x)` to get total mapped exon bases. |
| `block_starts` | 0-based offsets from `contig_start`; add `contig_start` to recover absolute contig coordinates. |
| `contig` | 24 distinct contigs in hap0 (one per primary chromosome + MT). |
| span range | 8 bp (tiny ncRNA) to 237,913,878 bp (very large gene on chr1). |

**Minimum coverage threshold:** 50% (`--min-coverage 0.5` in `get_tx_seqs.sh`).  
**HQ threshold:** 90% (`--full-coverage 0.9`).

---

### `transcript_summary`

One row per transcript (same set as `transcripts`).  Summarises the best alignment hit.

```sql
CREATE TABLE transcript_summary (
    transcript_pk  INTEGER PRIMARY KEY REFERENCES transcripts ON DELETE CASCADE,
    status         TEXT    NOT NULL,
    hit_count      INTEGER NOT NULL,   -- alignments above min-coverage threshold
    contig_count   INTEGER NOT NULL,   -- distinct contigs with hits
    best_coverage  REAL    NOT NULL,   -- coverage_pct of the best hit
    best_contig    TEXT    NOT NULL    -- contig name of the best hit
);
CREATE INDEX idx_summary_status ON transcript_summary(status);
```

**Status codes and counts (hap0):**

| Status | Count | % | Meaning |
|--------|------:|--:|---------|
| `single_full` | 188,095 | 90.7% | One hit, ≥ 90% exon coverage — **high quality** |
| `no_hit` | 17,599 | 8.5% | No alignment above 50% threshold |
| `multi_contig` | 1,062 | 0.5% | Best hit spans multiple contigs (fragmented assembly) |
| `single_partial` | 370 | 0.2% | One hit, 50–89% exon coverage |
| `multi_location` | 163 | 0.1% | Multiple hits on different contigs (segmental duplication) |

**Note on `no_hit` root causes:**
- Transcripts annotated only on alt/fix/unplaced loci (~1,598 NM genes) — those chromosomes are absent from the alndb
- Genes in centromere / heterochromatin / repeat-dense regions (~1,113 NM genes)
- chrY transcripts in the maternal haplotype (57 NM genes)

---

### `gene_summary`

One row per unique `gene_id` in the GTF.  Aggregates transcript status counts per gene.

```sql
CREATE TABLE gene_summary (
    gene_id        TEXT    PRIMARY KEY,
    gene_name      TEXT    NOT NULL,
    ref_chrom      TEXT    NOT NULL,
    total          INTEGER NOT NULL,
    single_full    INTEGER NOT NULL,
    single_partial INTEGER NOT NULL,
    multi_contig   INTEGER NOT NULL,
    multi_location INTEGER NOT NULL,
    no_hit         INTEGER NOT NULL
);
CREATE INDEX idx_gene_summary_chrom ON gene_summary(ref_chrom);
```

59,251 unique gene symbols in hap0.  A gene is considered **HQ** if any of its transcripts achieves `coverage_pct ≥ 90`.

---

### `ref_sequences` *(written with `--ref-fa`)*

Spliced exon sequence for each transcript as extracted from the GRCh38 reference FASTA.

```sql
CREATE TABLE ref_sequences (
    transcript_pk  INTEGER PRIMARY KEY REFERENCES transcripts ON DELETE CASCADE,
    ref_seq        TEXT NOT NULL   -- concatenated exon sequence (5'→3', coding strand)
);
```

---

### `contig_sequences` *(written with `--query-fa`)*

Spliced exon sequence as extracted from the best-hit assembly contig.

```sql
CREATE TABLE contig_sequences (
    liftover_pk    INTEGER PRIMARY KEY REFERENCES liftover ON DELETE CASCADE,
    contig_seq     TEXT NOT NULL   -- concatenated mapped exon sequence
);
```

Join with `liftover` to get the contig name and coverage, then with `transcripts` for the gene/transcript ID.

---

## Common Queries

### All HQ transcripts for a gene

```sql
SELECT t.transcript_id, t.ref_chrom, t.ref_start, t.ref_end,
       l.contig, l.contig_start, l.contig_end, l.coverage_pct
FROM transcripts t
JOIN liftover l ON t.transcript_pk = l.transcript_pk
JOIN transcript_summary ts ON t.transcript_pk = ts.transcript_pk
WHERE t.gene_name = 'BRCA1'
  AND ts.status = 'single_full'
ORDER BY l.coverage_pct DESC;
```

---

### NM gene liftover funnel (HQ threshold = 90%)

```sql
WITH nm_genes AS (
    SELECT DISTINCT gene_name
    FROM transcripts
    WHERE SUBSTR(transcript_id, 1, 2) = 'NM'
      AND gene_name != ''
      AND ref_chrom NOT LIKE '%_alt%'
      AND ref_chrom NOT LIKE '%_fix%'
      AND ref_chrom NOT LIKE '%_random%'
      AND ref_chrom NOT LIKE 'chrUn%'
),
best AS (
    SELECT t.gene_name, MAX(l.coverage_pct) AS best_cov
    FROM transcripts t
    JOIN liftover l ON t.transcript_pk = l.transcript_pk
    WHERE SUBSTR(t.transcript_id, 1, 2) = 'NM'
    GROUP BY t.gene_name
)
SELECT
    COUNT(*)                                                        AS total_nm_genes,
    SUM(CASE WHEN b.best_cov >= 90 THEN 1 ELSE 0 END)             AS hq_ge90,
    SUM(CASE WHEN b.best_cov >= 50 AND b.best_cov < 90 THEN 1 END) AS partial_50_90,
    COUNT(*) - COUNT(b.gene_name)                                  AS no_hit
FROM nm_genes n
LEFT JOIN best b ON n.gene_name = b.gene_name;
```

---

### Compare spliced length: reference vs contig (scatter plot data)

```sql
WITH ref_exon AS (
    SELECT t.transcript_pk,
           SUM(e.exon_end - e.exon_start) AS ref_exon_len
    FROM transcripts t
    JOIN exons e ON t.transcript_pk = e.transcript_pk
    GROUP BY t.transcript_pk
),
best_hit AS (
    SELECT l.transcript_pk,
           l.contig_end - l.contig_start AS contig_genomic_span,
           l.block_sizes
    FROM liftover l
    JOIN transcript_summary ts ON l.transcript_pk = ts.transcript_pk
    WHERE l.contig = ts.best_contig
      AND ts.status IN ('single_full', 'single_partial')
)
SELECT t.ref_end - t.ref_start  AS ref_genomic_span,   -- full gene body incl. introns
       r.ref_exon_len,                                  -- spliced exon bases
       b.contig_genomic_span,                           -- span on contig incl. introns
       b.block_sizes,                                   -- parse to get contig exon bases
       t.transcript_id,
       t.gene_name
FROM transcripts t
JOIN ref_exon r  ON t.transcript_pk = r.transcript_pk
JOIN best_hit b  ON t.transcript_pk = b.transcript_pk
ORDER BY t.transcript_pk;
```

**Interpreting the scatter plots:**

| Point position | Interpretation |
|----------------|----------------|
| On the diagonal | Perfect liftover — assembly matches reference |
| Above diagonal  | Contig is longer (assembly has an insertion or duplication) |
| Below diagonal  | Contig is shorter (assembly has a deletion) |
| Far below diagonal | Partial liftover — exons missed or contig break |

---

### SV candidates overlapping partial-liftover genes

```sql
ATTACH 'hg002_hap0.alndb' AS aln;

SELECT COUNT(DISTINCT t.gene_name)
FROM transcripts t
JOIN liftover l      ON t.transcript_pk = l.transcript_pk
JOIN aln.sv_candidates sv
  ON sv.target_name = 'GRCh38#0#' || t.ref_chrom
  AND sv.target_end   > t.ref_start
  AND sv.target_start < t.ref_end
WHERE SUBSTR(t.transcript_id, 1, 2) = 'NM'
GROUP BY t.gene_name
HAVING MAX(l.coverage_pct) >= 50 AND MAX(l.coverage_pct) < 90;
-- Returns ~1,030 partially-mapped NM genes that overlap at least one SV candidate
```

---

### Per-chromosome HQ gene summary

```sql
SELECT t.ref_chrom,
       COUNT(DISTINCT t.gene_name)                                       AS total_genes,
       COUNT(DISTINCT CASE WHEN l.coverage_pct >= 90 THEN t.gene_name END) AS hq_genes,
       ROUND(100.0 * COUNT(DISTINCT CASE WHEN l.coverage_pct >= 90
             THEN t.gene_name END) / COUNT(DISTINCT t.gene_name), 1)    AS pct_hq
FROM transcripts t
LEFT JOIN liftover l ON t.transcript_pk = l.transcript_pk
WHERE t.gene_name != ''
  AND t.ref_chrom NOT LIKE '%_alt%'
  AND t.ref_chrom NOT LIKE '%_fix%'
  AND t.ref_chrom NOT LIKE '%_random%'
  AND t.ref_chrom NOT LIKE 'chrUn%'
GROUP BY t.ref_chrom
ORDER BY t.ref_chrom;
```

---

## Notes on Coordinate Systems

- All coordinates in the liftover DB are **0-based, half-open** intervals `[start, end)`,
  consistent with BED format and Python slice notation.
- `ref_chrom` uses the bare `chr` prefix (`chr1`, `chrX`, `chrM`), **not** the PanSN format
  (`GRCh38#0#chr1`) used inside `.alndb`.  When joining across the two databases,
  prepend `'GRCh38#0#'` to `ref_chrom`.
- `contig` names use the PanSN format as written by the assembly
  (e.g., `HG002#1#chr1` for maternal chr1, `HG002#2#chr1` for paternal).

---

## Files Written by `pgr align liftover-gtf`

| File | Contents |
|------|----------|
| `*_liftover.db` | This SQLite database |
| `*.gtf` | All transcripts with `coverage_pct ≥ min_coverage` in GTF format |
| `*.hq.gtf` | Transcripts with `coverage_pct ≥ full_coverage` (≥ 90%) |
| `*.anomaly.tsv` | Transcript-level anomalies (multi-contig, low coverage, etc.) |
| `*.gene.anomaly.tsv` | Gene-level anomaly summary |
| `*_ref_tx.fa` | Spliced reference sequences (with `--ref-fa`) |
| `*_contig_tx.fa` | All mapped contig sequences (with `--query-fa`) |
| `*_hq_contig_tx.fa` | HQ-only contig sequences |
