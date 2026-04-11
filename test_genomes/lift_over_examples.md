# pgr align liftover-gtf — usage and query examples

`pgr align liftover-gtf` maps GTF transcript annotations from a reference genome (the
*target* side of a `pgr align alnmap` alignment) to haplotype contigs (the *query* side)
using the alignment data stored in an `.alndb` file.  Results are written to a
SQLite database with three linked tables.

---

## Running the tool

```bash
cd test_genomes

../target/release/pgr align liftover-gtf \
    --alndb-path hg002_hap0.alndb \
    --gtf-path hg38.ncbiRefSeq.gtf.gz \
    --output-db hg002_hap0_liftover.db \
    --target-chr-prefix "GRCh38#0#" \
    --min-coverage 0.5
```

| Argument / flag | Description |
|---|---|
| `--alndb-path` | Alignment database produced by `pgr align alnmap` |
| `--gtf-path` | Reference annotation (plain or `.gz`) |
| `--output-db` | Output SQLite database (overwritten if it exists) |
| `--target-chr-prefix` | Prefix to strip from alndb target names so they match GTF chromosome names.  Use `"GRCh38#0#"` for PanSN-formatted references. |
| `--min-coverage` | Minimum fraction of exon bases that must land on one contig (default `0.5`). |

---

## Output database schema

```
transcripts (transcript_pk PK, transcript_id UNIQUE, gene_id, gene_name,
             ref_chrom, ref_start, ref_end, ref_strand, exon_count)
    │
    ├─── exons    (exon_pk PK, transcript_pk FK → transcripts, exon_start, exon_end)
    │
    └─── liftover (liftover_pk PK, transcript_pk FK → transcripts,
                   contig, contig_start, contig_end, strand,
                   coverage_pct, block_count, block_sizes, block_starts)
```

Both `exons` and `liftover` reference `transcripts.transcript_pk` with
`ON DELETE CASCADE`.  Rows in `transcripts` preserve the order they appear in the
input GTF.

### Row counts (hg002 haplotype-0 example)

| Table | Rows |
|---|---|
| `transcripts` | 207,289 |
| `exons` | 2,214,932 |
| `liftover` | 174,267 |
| Transcripts with no hit | 33,149 |

---

## Example queries

### 1 — Summary statistics

```sql
-- transcripts lifted vs. not lifted
SELECT
    count(*) FILTER (WHERE transcript_pk     IN (SELECT transcript_pk FROM liftover)) AS lifted,
    count(*) FILTER (WHERE transcript_pk NOT IN (SELECT transcript_pk FROM liftover)) AS no_hit
FROM transcripts;
```

```sql
-- lifted transcripts per contig, most populated first
SELECT contig, count(*) AS n
FROM liftover
GROUP BY contig
ORDER BY n DESC
LIMIT 10;
-- HG002#2#chr1   16900
-- HG002#2#chr2   13254
-- HG002#2#chr3   11173
-- HG002#2#chr11   9722
-- HG002#2#chr12   9317
```

---

### 2 — Look up a gene (BRCA2)

```sql
SELECT t.transcript_id,
       t.ref_chrom, t.ref_start, t.ref_end, t.ref_strand,
       l.contig, l.contig_start, l.contig_end, l.strand,
       round(l.coverage_pct, 1) AS cov,
       l.block_count
FROM liftover l
JOIN transcripts t USING(transcript_pk)
WHERE t.gene_name = 'BRCA2';
-- NM_000059.4 | chr13 | 32315507 | 32400268 | + | HG002#2#chr13 | 31261907 | 31345680 | + | 71.3 | 26
```

```sql
-- exon coordinates for BRCA2 on the reference
SELECT t.transcript_id, e.exon_start, e.exon_end
FROM transcripts t
JOIN exons e USING(transcript_pk)
WHERE t.gene_name = 'BRCA2'
ORDER BY e.exon_start;
-- NM_000059.4 | 32315507 | 32315667
-- NM_000059.4 | 32316421 | 32316527
-- ...
```

---

### 3 — Best liftover hit per transcript (highest coverage)

```sql
SELECT t.gene_name, t.transcript_id,
       l.contig, l.contig_start, l.contig_end, l.strand,
       round(l.coverage_pct, 1) AS cov
FROM liftover l
JOIN transcripts t USING(transcript_pk)
WHERE l.coverage_pct = (
    SELECT max(coverage_pct)
    FROM liftover
    WHERE transcript_pk = t.transcript_pk
)
ORDER BY t.ref_chrom, t.ref_start
LIMIT 10;
```

---

### 4 — Clinically relevant genes (BRCA1, TP53, MYC)

```sql
SELECT t.gene_name, l.contig,
       l.contig_start, l.contig_end, l.strand,
       round(l.coverage_pct, 1) AS cov
FROM liftover l
JOIN transcripts t USING(transcript_pk)
WHERE t.gene_name IN ('BRCA1', 'TP53', 'MYC')
ORDER BY t.gene_name, l.coverage_pct DESC;
-- BRCA1 | HG002#2#chr17 | 43405160 | 43484981 | + | 100.0
-- MYC   | HG002#2#chr8  | 129342625| 129349346 | + | 100.0
-- TP53  | HG002#2#chr17 | 7574699  | 7580570   | + |  78.7
```

---

### 5 — High-confidence liftovers (≥ 90 % coverage) with exon structure

```sql
SELECT t.gene_name, t.transcript_id,
       l.contig, l.contig_start, l.contig_end, l.strand,
       round(l.coverage_pct, 1) AS cov,
       l.block_count,
       l.block_sizes,
       l.block_starts
FROM liftover l
JOIN transcripts t USING(transcript_pk)
WHERE l.coverage_pct >= 90
ORDER BY t.ref_chrom, t.ref_start
LIMIT 5;
```

---

### 6 — Transcripts that map to multiple contigs

```sql
SELECT t.transcript_id, t.gene_name, count(DISTINCT l.contig) AS n_contigs
FROM liftover l
JOIN transcripts t USING(transcript_pk)
GROUP BY l.transcript_pk
HAVING n_contigs > 1
ORDER BY n_contigs DESC
LIMIT 10;
```

---

### 7 — Export lifted transcripts as BED-like TSV

```bash
sqlite3 -separator $'\t' hg002_hap0_liftover.db \
  "SELECT l.contig, l.contig_start, l.contig_end,
          t.transcript_id, round(l.coverage_pct,1), l.strand,
          l.contig_start, l.contig_end, '0',
          l.block_count, l.block_sizes, l.block_starts
   FROM liftover l
   JOIN transcripts t USING(transcript_pk)
   WHERE l.coverage_pct >= 90
   ORDER BY l.contig, l.contig_start;" \
  > hg002_hap0_liftover_90pct.bed12
```

---

### 8 — Running on haplotype 1

```bash
../target/release/pgr align liftover-gtf \
    --alndb-path hg002_hap1.alndb \
    --gtf-path hg38.ncbiRefSeq.gtf.gz \
    --output-db hg002_hap1_liftover.db \
    --target-chr-prefix "GRCh38#0#" \
    --min-coverage 0.5
```

---

### 9 — Generate a GTF file for haplotype contigs from the database

The script below queries the database and reconstructs a GTF annotation on the
haplotype contigs.  For each transcript it emits one `transcript` record and one
`exon` record per alignment block (derived from `block_sizes` / `block_starts`).
By default only the **best hit** (highest `coverage_pct`) per transcript is kept;
pass `--all-hits` to write every hit above the coverage threshold.

Save as `liftover_to_gtf.py` and run:

```bash
python3 liftover_to_gtf.py hg002_hap0_liftover.db hg002_hap0_liftover.gtf
# or with a custom coverage threshold:
python3 liftover_to_gtf.py hg002_hap0_liftover.db hg002_hap0_liftover.gtf --min-cov 90
# include all hits above the threshold (not just the best per transcript):
python3 liftover_to_gtf.py hg002_hap0_liftover.db hg002_hap0_liftover.gtf --all-hits
```

```python
#!/usr/bin/env python3
"""Generate a GTF file for haplotype contigs from a pgr align liftover-gtf SQLite database.

Usage:
    python3 liftover_to_gtf.py <db> <output.gtf> [--min-cov FLOAT] [--all-hits]
"""
import sqlite3, sys, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db",         help="liftover SQLite database")
    ap.add_argument("output_gtf", help="output GTF file")
    ap.add_argument("--min-cov",  type=float, default=50.0,
                    help="minimum coverage_pct to include (default 50)")
    ap.add_argument("--all-hits", action="store_true",
                    help="write every hit above --min-cov instead of only the best per transcript")
    args = ap.parse_args()

    conn = sqlite3.connect(args.db)

    if args.all_hits:
        # All liftover rows above the coverage threshold
        query = """
            SELECT l.contig, l.contig_start, l.contig_end, l.strand,
                   l.block_sizes, l.block_starts, l.coverage_pct,
                   t.transcript_id, t.gene_id, t.gene_name,
                   t.ref_chrom, t.ref_start, t.ref_end, t.ref_strand
            FROM liftover l
            JOIN transcripts t USING(transcript_pk)
            WHERE l.coverage_pct >= ?
            ORDER BY l.contig, l.contig_start
        """
        params = (args.min_cov,)
    else:
        # Best hit (highest coverage_pct) per transcript, preserving GTF order
        # via the transcript_pk (which reflects GTF insertion order).
        query = """
            SELECT l.contig, l.contig_start, l.contig_end, l.strand,
                   l.block_sizes, l.block_starts, l.coverage_pct,
                   t.transcript_id, t.gene_id, t.gene_name,
                   t.ref_chrom, t.ref_start, t.ref_end, t.ref_strand
            FROM liftover l
            JOIN transcripts t USING(transcript_pk)
            WHERE l.coverage_pct >= ?
              AND l.liftover_pk = (
                  SELECT liftover_pk FROM liftover l2
                  WHERE l2.transcript_pk = l.transcript_pk
                  ORDER BY l2.coverage_pct DESC LIMIT 1
              )
            ORDER BY l.contig, l.contig_start
        """
        params = (args.min_cov,)

    source = "pgr-liftover"
    written = 0

    with open(args.output_gtf, "w") as out:
        # GTF uses no version header — comment lines only
        out.write(f"# source db: {args.db}\n")
        out.write(f"# min_coverage_pct: {args.min_cov}\n")

        for row in conn.execute(query, params):
            (contig, cstart, cend, strand,
             block_sizes_str, block_starts_str, cov,
             tx_id, gene_id, gene_name,
             ref_chrom, ref_start, ref_end, ref_strand) = row

            attrs = (
                f'gene_id "{gene_id}"; transcript_id "{tx_id}"; '
                f'gene_name "{gene_name}"; '
                f'ref_chrom "{ref_chrom}"; ref_start {ref_start}; '
                f'ref_end {ref_end}; ref_strand "{ref_strand}"; '
                f'coverage_pct {cov:.1f};'
            )

            # transcript record — GTF is 1-based inclusive
            out.write(
                f"{contig}\t{source}\ttranscript\t"
                f"{cstart + 1}\t{cend}\t.\t{strand}\t.\t{attrs}\n"
            )

            # exon records — reconstruct from block_sizes / block_starts
            sizes  = [int(x) for x in block_sizes_str.split(",")]
            starts = [int(x) for x in block_starts_str.split(",")]
            for size, rel in zip(sizes, starts):
                exon_s = cstart + rel
                exon_e = cstart + rel + size
                out.write(
                    f"{contig}\t{source}\texon\t"
                    f"{exon_s + 1}\t{exon_e}\t.\t{strand}\t.\t{attrs}\n"
                )

            written += 1

    conn.close()
    print(f"Wrote {written} transcript records to {args.output_gtf}")

if __name__ == "__main__":
    main()
```

**Expected output (first few lines of `hg002_hap0_liftover.gtf`):**

```
# source db: hg002_hap0_liftover.db
# min_coverage_pct: 50.0
HG002#2#chrM  pgr-liftover  transcript  15379  15446  .  -  .  gene_id "TRNP"; transcript_id "rna-TRNP"; gene_name "TRNP"; ...
HG002#2#chrM  pgr-liftover  exon        15379  15446  .  -  .  gene_id "TRNP"; transcript_id "rna-TRNP"; gene_name "TRNP"; ...
HG002#2#chrM  pgr-liftover  transcript  15311  15376  .  +  .  gene_id "TRNT"; transcript_id "rna-TRNT"; ...
```
