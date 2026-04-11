# alnmap SQLite improvement plan

## Motivation

`pgr align alnmap` currently writes 8–9 separate files per haplotype run. Downstream tools
consume overlapping subsets of those files, requiring fragile multi-argument command lines
and full linear scans of large text files. This document describes a plan to replace all
per-run outputs with a single SQLite database file (`.alndb`) and update the downstream
tools accordingly.

---

## Current outputs and downstream consumers

| File | Downstream consumer in pipeline |
|------|----------------------------------|
| `<prefix>.alnmap` | `pgr variant diploid-vcf` (V and M records) |
| `<prefix>.target_len.json` | `pgr variant diploid-vcf` (VCF header contig lines) |
| `<prefix>.ctgmap.json` | `pgr plot chr-aln` |
| `<prefix>.vcf` | not used downstream (redundant with diploid VCF) |
| `<prefix>.ctgmap.bed` | not used in current pipeline |
| `<prefix>.query_len.json` | not used in current pipeline |
| `<prefix>.svcnd.bed` | not used in current pipeline (future SV analysis) |
| `<prefix>.ctgsv.bed` | not used in current pipeline (future SV analysis) |
| `<prefix>.svcnd.seqs` | `pgr variant sv-analysis` (not yet in pipeline) |

---

## Proposed SQLite schema — `<prefix>.alndb`

All original names and coordinates are stored verbatim. No normalisation at the DB layer.

```sql
-- Run provenance: shimmer spec and command-line parameters
CREATE TABLE run_params (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

-- Reference (target) and query sequence lengths
-- Replaces target_len.json and query_len.json
CREATE TABLE sequences (
    seq_id   INTEGER PRIMARY KEY,
    seq_name TEXT NOT NULL,
    seq_type TEXT NOT NULL,    -- 'target' or 'query'
    length   INTEGER NOT NULL
);

-- One row per alignment chain (B/E record pair)
-- Replaces B/E records in .alnmap
CREATE TABLE chains (
    aln_idx         INTEGER PRIMARY KEY,
    target_name     TEXT NOT NULL,
    target_start    INTEGER NOT NULL,
    target_end      INTEGER NOT NULL,
    query_name      TEXT NOT NULL,
    query_start     INTEGER NOT NULL,
    query_end       INTEGER NOT NULL,
    orientation     INTEGER NOT NULL,    -- 0=fwd 1=rc
    ctg_orientation INTEGER NOT NULL,
    query_length    INTEGER NOT NULL,
    target_dup      INTEGER NOT NULL,    -- 0/1
    target_ovlp     INTEGER NOT NULL,
    query_dup       INTEGER NOT NULL,
    query_ovlp      INTEGER NOT NULL
);

-- Match (M) and SV candidate (S) records
-- Replaces M/S records in .alnmap
CREATE TABLE blocks (
    block_id        INTEGER PRIMARY KEY,
    aln_idx         INTEGER NOT NULL REFERENCES chains(aln_idx),
    block_type      INTEGER NOT NULL,    -- 0=M 1=S
    dup_flag        INTEGER NOT NULL,
    ovlp_flag       INTEGER NOT NULL,
    target_name     TEXT NOT NULL,
    target_start    INTEGER NOT NULL,
    target_end      INTEGER NOT NULL,
    query_name      TEXT NOT NULL,
    query_start     INTEGER NOT NULL,
    query_end       INTEGER NOT NULL,
    orientation     INTEGER NOT NULL,
    ctg_orientation INTEGER,             -- S records only
    sv_diff_type    INTEGER              -- S records only: 0=FailAln 1=FailEndMatch
                                         --   2=FailShortSeq 3=FailLengthDiff 4=Unknown
);
CREATE INDEX blocks_pos ON blocks(target_name, target_start, target_end);

-- Variant (V) records
-- Replaces V records in .alnmap
CREATE TABLE variants (
    variant_id   INTEGER PRIMARY KEY,
    aln_idx      INTEGER NOT NULL REFERENCES chains(aln_idx),
    dup_flag     INTEGER NOT NULL,
    ovlp_flag    INTEGER NOT NULL,
    target_name  TEXT NOT NULL,
    target_start INTEGER NOT NULL,    -- alignment block start
    target_end   INTEGER NOT NULL,    -- alignment block end
    query_name   TEXT NOT NULL,
    query_start  INTEGER NOT NULL,
    query_end    INTEGER NOT NULL,
    orientation  INTEGER NOT NULL,
    target_diff  INTEGER NOT NULL,
    query_diff   INTEGER NOT NULL,
    target_coord INTEGER NOT NULL,    -- absolute position = target_start + target_diff
    variant_type INTEGER NOT NULL,    -- 0=X(SNP/MNP) 1=I(insertion) 2=D(deletion)
    ref_seq      TEXT NOT NULL,
    alt_seq      TEXT NOT NULL
);
CREATE INDEX variants_pos ON variants(target_name, target_coord);

-- Contig-to-reference mapping
-- Replaces ctgmap.bed and ctgmap.json
CREATE TABLE ctgmap (
    ctgmap_id       INTEGER PRIMARY KEY,
    target_name     TEXT NOT NULL,
    target_start    INTEGER NOT NULL,
    target_end      INTEGER NOT NULL,
    query_name      TEXT NOT NULL,
    query_start     INTEGER NOT NULL,
    query_end       INTEGER NOT NULL,
    ctg_len         INTEGER NOT NULL,
    orientation     INTEGER NOT NULL,
    ctg_orientation INTEGER NOT NULL,
    target_dup      INTEGER NOT NULL,
    target_ovlp     INTEGER NOT NULL,
    query_dup       INTEGER NOT NULL,
    query_ovlp      INTEGER NOT NULL
);
CREATE INDEX ctgmap_tgt ON ctgmap(target_name, target_start, target_end);
CREATE INDEX ctgmap_qry ON ctgmap(query_name, query_start, query_end);

-- SV candidate regions
-- Replaces svcnd.bed
CREATE TABLE sv_candidates (
    svcnd_id     INTEGER PRIMARY KEY,
    target_name  TEXT NOT NULL,
    target_start INTEGER NOT NULL,
    target_end   INTEGER NOT NULL,
    sv_type      TEXT NOT NULL
);
CREATE INDEX svcnd_pos ON sv_candidates(target_name, target_start, target_end);

-- Contig-level SV summary
-- Replaces ctgsv.bed
CREATE TABLE ctgsv (
    ctgsv_id    INTEGER PRIMARY KEY,
    query_name  TEXT NOT NULL,
    query_start INTEGER NOT NULL,
    query_end   INTEGER NOT NULL,
    sv_type     TEXT NOT NULL
);

-- Sequences spanning SV candidates (optional — omitted with --skip-uncalled-sv-seq-file)
-- Replaces svcnd.seqs
CREATE TABLE sv_sequences (
    svcnd_id    INTEGER NOT NULL REFERENCES sv_candidates(svcnd_id),
    target_seq  TEXT NOT NULL,
    query_seq   TEXT NOT NULL
);
```

### Design rules

- **All original names stored verbatim.** PanSN prefixes (`assembly#hap#chr`) are preserved
  in `target_name` and `query_name` columns throughout. No normalisation at the DB layer.
- **Normalisation is handled by third-party tools in the pipeline.** PanSN stripping and
  chromosome renaming are done with standard bioinformatics tools (`sed`, `bcftools
  annotate --rename-chrs`) in the pipeline script, not inside any pgr-tk binary. This
  keeps the Rust code free of format-specific post-processing logic and makes the
  normalisation steps explicit and auditable in the shell script.
- **Indices created after all inserts.** `pgr align alnmap` writes records in contig order.
  Creating indices after the final insert is faster than maintaining them incrementally.
- **WAL mode enabled.** Allows `pgr variant diploid-vcf` to open both `hap0.alndb` and
  `hap1.alndb` concurrently without write contention.

---

## Changes required per tool

### `pgr align alnmap`

- Open a `rusqlite::Connection` to `<prefix>.alndb` at startup (WAL mode).
- Create the schema (DDL above) before processing begins.
- Write records into the appropriate table as they are produced — streaming inserts
  wrapped in a single transaction per contig for throughput.
- Create indices after all contigs are processed.
- Store shimmer spec and preset in `run_params`.
- **Retain** `.vcf` as a separate file output (standard format, consumed by bcftools
  and other VCF-aware tools). All other current file outputs are replaced by the DB.

### `pgr variant diploid-vcf`

- Accept `.alndb` paths instead of `.alnmap` + `.target_len.json`.
- Replace the full text-scan + in-memory sort with two indexed queries:

  ```sql
  -- variants, pre-sorted by position
  SELECT aln_idx, target_name, target_coord, dup_flag, ovlp_flag,
         target_start, target_end, ref_seq, alt_seq, variant_type
    FROM variants
   ORDER BY target_name, target_coord;

  -- blocks for interval building (M and V only)
  SELECT aln_idx, target_name, target_start, target_end,
         query_name, query_start, query_end, orientation, dup_flag, ovlp_flag
    FROM blocks
   WHERE block_type IN (0, 1);

  -- contig lengths for VCF header
  SELECT seq_name, length FROM sequences WHERE seq_type = 'target';
  ```

- `target_name` values written to the VCF `CHROM` column are emitted verbatim (PanSN
  names intact). PanSN stripping is done in the pipeline script after the fact using
  `sed` or `bcftools`, not inside this tool.
- The `variant_records.sort()` call in the current code is eliminated — the ORDER BY
  in the query delivers records already sorted.
- Remove the `target_len_json_path` positional argument; contig lengths come from the
  `sequences` table of either `.alndb` (both haplotypes share the same reference).

### `pgr plot chr-aln`

- Accept an `.alndb` path instead of `.ctgmap.json`.
- Query `ctgmap` and `sequences` tables directly:

  ```sql
  SELECT t_name, ts, te, q_name, qs, qe, ctg_len,
         orientation, ctg_orientation, target_dup, target_ovlp, query_dup, query_ovlp
    FROM ctgmap;

  SELECT seq_name, length FROM sequences WHERE seq_type = 'target';
  SELECT seq_name, length FROM sequences WHERE seq_type = 'query';
  ```

- No JSON parsing, no separate file argument for lengths.

### `run_e2e.sh`

- Step 1 outputs change from 8 files to 1 `.alndb` + 1 `.vcf` per haplotype.
- Step 2 (`pgr variant diploid-vcf`) argument list simplifies:
  ```bash
  pgr variant diploid-vcf hg002_hap0.alndb hg002_hap1.alndb hg002 --sample-name hg002
  ```
- Step 2b (`pgr plot chr-aln`) argument list simplifies:
  ```bash
  pgr plot chr-aln hg002_hap0.alndb hg002_hap0_aln_plot
  ```
- Steps 2 and 2b are independent and can run in parallel (both depend only on Step 1).
- Step 3a (`sed` PanSN strip) is **retained** — PanSN names pass through verbatim from
  the DB into the VCF, and the `sed` step in the pipeline script is the designated place
  to strip them before downstream VCF tools.
- Step 4 (`bcftools annotate --rename-chrs`) is retained for the `chr` → bare rename
  required to match ClinVar chromosome naming; this is a separate concern from PanSN.

---

## File count comparison

| | Current | After |
|-|---------|-------|
| Files produced per haplotype by `pgr align alnmap` | 8–9 | 2 (`.alndb` + `.vcf`) |
| Arguments to `pgr variant diploid-vcf` | 4 positional | 3 positional |
| Arguments to `pgr plot chr-aln` | 2 positional | 2 positional (type changes) |
| Intermediate files in pipeline | 10+ | 4 |
| sed/awk post-processing steps | 1 | 1 (retained — PanSN strip stays in script) |

---

---

## Phase 2 — combined diploid sample database (`.sampledb`)

### Concept

Rather than keeping two separate `.alndb` files (one per haplotype), merge them into a
single **sample database** (`<sample>.sampledb`). This file becomes the single source of
truth for a diploid sample: both haplotype alignments, all derived variant calls, and all
external annotations (ClinVar, gene models, regulatory regions, custom tracks) stored in
one place.

```
hg002_hap0.alndb  ─┐
                    ├──► hg002.sampledb
hg002_hap1.alndb  ─┘
```

A new tool `pgr-build-sampledb` merges the two `.alndb` files and initialises the
annotation tables. Subsequent annotation tools append rows to the annotation tables
without touching the alignment data.

---

### Schema additions over the per-haplotype `.alndb`

All alignment tables from the single-haplotype schema are retained with one extra column:

```sql
-- Sample-level metadata
CREATE TABLE sample (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
    -- e.g. sample_name, reference_name, reference_build, creation_date
);
```

Every alignment table gains a `haplotype` column:

```sql
-- haplotype: 0 = hap0 (maternal), 1 = hap1 (paternal)
-- Added to: chains, blocks, variants, ctgmap, sv_candidates, ctgsv, sv_sequences
ALTER TABLE chains       ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
ALTER TABLE blocks       ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
ALTER TABLE variants     ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
ALTER TABLE ctgmap       ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
ALTER TABLE sv_candidates ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
ALTER TABLE ctgsv        ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
ALTER TABLE sv_sequences ADD COLUMN haplotype INTEGER NOT NULL DEFAULT 0;
```

`aln_idx` is re-keyed per haplotype to avoid collision:

```sql
-- Composite primary key on alignment tables
-- aln_idx is unique within a haplotype, (haplotype, aln_idx) is globally unique
CREATE INDEX chains_hap   ON chains(haplotype, aln_idx);
CREATE INDEX variants_hap ON variants(haplotype, target_name, target_coord);
CREATE INDEX blocks_hap   ON blocks(haplotype, target_name, target_start, target_end);
CREATE INDEX ctgmap_hap   ON ctgmap(haplotype, target_name, target_start, target_end);
```

---

### Diploid variant table

The merged diploid calls (output of `pgr variant diploid-vcf`) are stored back into the
sample DB rather than only as a VCF file:

```sql
CREATE TABLE diploid_variants (
    dv_id        INTEGER PRIMARY KEY,
    target_name  TEXT NOT NULL,    -- PanSN name, verbatim
    target_coord INTEGER NOT NULL, -- 0-based
    ref_seq      TEXT NOT NULL,
    alt_seq      TEXT NOT NULL,    -- comma-separated if multi-allelic
    genotype     TEXT NOT NULL,    -- e.g. "0|1", "1|2"
    filter       TEXT NOT NULL,    -- PASS / DUP / OVLP / NC
    quality      INTEGER NOT NULL
);
CREATE INDEX dv_pos ON diploid_variants(target_name, target_coord);
```

---

### Annotation layer

Annotations are stored in typed tables. Each annotation table links to
`diploid_variants` via a match on `(target_name, target_coord, ref_seq, alt_seq)`, or to
genomic regions independently of any called variant.

#### Variant-level annotations

```sql
-- ClinVar annotations (loaded from clinvar.vcf.gz via bcftools query or direct VCF parse)
CREATE TABLE annot_clinvar (
    dv_id       INTEGER REFERENCES diploid_variants(dv_id),
    target_name TEXT NOT NULL,
    target_coord INTEGER NOT NULL,
    ref_seq     TEXT NOT NULL,
    alt_seq     TEXT NOT NULL,
    clnsig      TEXT,    -- CLNSIG field: Pathogenic, Benign, etc.
    clndn       TEXT,    -- CLNDN field: disease name
    clnrevstat  TEXT     -- CLNREVSTAT field: review status
);
CREATE INDEX annot_clinvar_pos ON annot_clinvar(target_name, target_coord);

-- Gene-level variant effect (loaded from pgr variant annotate-vcf output or VEP/SnpEff)
CREATE TABLE annot_gene_effect (
    dv_id         INTEGER REFERENCES diploid_variants(dv_id),
    target_name   TEXT NOT NULL,
    target_coord  INTEGER NOT NULL,
    gene_name     TEXT,
    transcript_id TEXT,
    consequence   TEXT,    -- e.g. missense_variant, synonymous_variant, intron_variant
    strand        INTEGER  -- 0=fwd 1=rev
);
CREATE INDEX annot_gene_pos ON annot_gene_effect(target_name, target_coord);
```

#### Region-level annotations

```sql
-- Gene models from GTF (loaded once, shared across samples if using a common reference)
CREATE TABLE annot_genes (
    gene_id     INTEGER PRIMARY KEY,
    target_name TEXT NOT NULL,
    gene_start  INTEGER NOT NULL,
    gene_end    INTEGER NOT NULL,
    gene_name   TEXT,
    gene_type   TEXT,    -- protein_coding, lncRNA, etc.
    strand      INTEGER
);
CREATE INDEX annot_genes_pos ON annot_genes(target_name, gene_start, gene_end);

-- Exon / CDS features
CREATE TABLE annot_features (
    feature_id    INTEGER PRIMARY KEY,
    gene_id       INTEGER REFERENCES annot_genes(gene_id),
    target_name   TEXT NOT NULL,
    feature_start INTEGER NOT NULL,
    feature_end   INTEGER NOT NULL,
    feature_type  TEXT NOT NULL,   -- exon / CDS / UTR
    transcript_id TEXT
);
CREATE INDEX annot_feat_pos ON annot_features(target_name, feature_start, feature_end);

-- Generic extensible track (BED-style, for regulatory regions, repeat masker, etc.)
CREATE TABLE annot_regions (
    region_id    INTEGER PRIMARY KEY,
    track_name   TEXT NOT NULL,    -- e.g. "ENCODE_enhancers", "RepeatMasker"
    target_name  TEXT NOT NULL,
    region_start INTEGER NOT NULL,
    region_end   INTEGER NOT NULL,
    name         TEXT,
    score        REAL,
    strand       INTEGER
);
CREATE INDEX annot_regions_pos ON annot_regions(track_name, target_name, region_start, region_end);
```

---

### Annotation loading tools

Each annotation source gets a small dedicated loader that appends rows to the appropriate
table. No pgr-tk binary needs to understand annotation formats — loaders are thin scripts
or standalone tools:

| Source | Loader | Populates |
|--------|--------|-----------|
| `clinvar.vcf.gz` | `pgr-load-clinvar` (or Python script via `bcftools query`) | `annot_clinvar` |
| GTF file | `pgr-load-gtf` | `annot_genes`, `annot_features` |
| `pgr variant annotate-vcf` output | `pgr-load-vcf-annot` | `annot_gene_effect` |
| Any BED file | `pgr-load-bed-track` | `annot_regions` |
| `diploid_variants` VCF | `pgr-load-diploid-vcf` | `diploid_variants` |

All loaders use the same approach: open the `.sampledb`, begin a transaction, insert
rows, commit. They are idempotent with a `DELETE FROM <table> WHERE track_name = ?`
before reload, or use `INSERT OR REPLACE`.

---

### Representative queries enabled by the combined DB

```sql
-- All pathogenic ClinVar variants in the sample, with genotype
SELECT dv.target_name, dv.target_coord, dv.ref_seq, dv.alt_seq,
       dv.genotype, c.clnsig, c.clndn
  FROM diploid_variants dv
  JOIN annot_clinvar c USING (target_name, target_coord, ref_seq, alt_seq)
 WHERE c.clnsig LIKE '%Pathogenic%';

-- Heterozygous variants overlapping protein-coding exons
SELECT dv.target_name, dv.target_coord, dv.ref_seq, dv.alt_seq,
       dv.genotype, g.gene_name, f.feature_type
  FROM diploid_variants dv
  JOIN annot_features f ON dv.target_name = f.target_name
                       AND dv.target_coord BETWEEN f.feature_start AND f.feature_end
  JOIN annot_genes g USING (gene_id)
 WHERE dv.genotype LIKE '0|%' OR dv.genotype LIKE '%|0'
   AND f.feature_type = 'CDS';

-- SV candidates in both haplotypes overlapping a gene
SELECT sc.target_name, sc.target_start, sc.target_end,
       sc.haplotype, g.gene_name
  FROM sv_candidates sc
  JOIN annot_genes g ON sc.target_name = g.target_name
                    AND sc.target_start < g.gene_end
                    AND sc.target_end   > g.gene_start;
```

---

### Updated pipeline sketch

```
pgr align alnmap hap0  ──►  hg002_hap0.alndb
pgr align alnmap hap1  ──►  hg002_hap1.alndb
                         │
                   pgr-build-sampledb
                         │
                   hg002.sampledb   ◄── pgr-load-diploid-vcf  (hg002.vcf → diploid_variants)
                         │           ◄── pgr-load-gtf          (hg38.ncbiRefSeq.gtf.gz → annot_genes/features)
                         │           ◄── pgr-load-clinvar      (clinvar.vcf.gz → annot_clinvar)
                         │           ◄── pgr-load-vcf-annot    (pgr variant annotate-vcf output → annot_gene_effect)
                         │           ◄── pgr-load-bed-track    (any BED → annot_regions)
                         ▼
               single queryable file: all alignments +
               diploid calls + all annotations
```

---

### Design principles

- **Append-only annotation.** Alignment tables are written once by `pgr align alnmap` and never
  modified. Annotation tables are appended to by separate loaders and can be reloaded
  independently.
- **No annotation logic in Rust alignment code.** Loaders can be Python scripts, R
  scripts, or small standalone Rust binaries — the interface is just SQLite inserts.
- **PanSN names preserved throughout.** All `target_name` values in both alignment and
  annotation tables use the original PanSN form. Name normalisation (stripping PanSN,
  mapping `chr1` ↔ `1`) is done at query time or in the loader's pre-processing step,
  not in the schema. A `name_map` table can optionally record the mapping:

  ```sql
  CREATE TABLE name_map (
      stored_name TEXT PRIMARY KEY,   -- e.g. "GRCh38#0#chr1"
      canonical   TEXT NOT NULL       -- e.g. "chr1"
  );
  ```

- **Reference annotations are shareable.** `annot_genes`, `annot_features`, and
  `annot_regions` tracks that depend only on the reference (not the sample) can be loaded
  into a separate read-only `reference.annotdb` and ATTACHed at query time, avoiding
  duplication across many sample databases:

  ```sql
  ATTACH 'GRCh38.annotdb' AS ref;
  SELECT * FROM diploid_variants dv
  JOIN ref.annot_genes g ON dv.target_name = g.target_name ...;
  ```

---

## Implementation order

### Phase 1 — per-haplotype `.alndb`

1. Add SQLite write to `pgr align alnmap` alongside existing file outputs (no tool breakage).
2. Add SQLite read to `pgr variant diploid-vcf`; keep `.alnmap` read path under a
   `--legacy` flag during transition. No PanSN logic added to Rust code — the pipeline
   script handles it with `sed` as before.
3. Add SQLite read to `pgr plot chr-aln`; keep `.ctgmap.json` read path under
   `--legacy` flag.
4. Update `run_e2e.sh`: new `.alndb` arguments for steps 2 and 2b; `sed` PanSN strip
   and `bcftools rename-chrs` steps unchanged.
5. Remove legacy file outputs from `pgr align alnmap` and legacy read paths from consumers
   once pipeline is validated end-to-end.

### Phase 2 — combined `.sampledb`

6. Implement `pgr-build-sampledb`: merges `hap0.alndb` + `hap1.alndb` into
   `<sample>.sampledb`, adding the `haplotype` column and composite indices.
7. Implement `pgr-load-diploid-vcf`: reads `hg002.vcf`, populates `diploid_variants`.
8. Implement `pgr-load-gtf`: reads GTF, populates `annot_genes` + `annot_features`.
   Designed to write to a shared `reference.annotdb` as well as per-sample.
9. Implement `pgr-load-clinvar`: reads ClinVar VCF, populates `annot_clinvar`.
   Handles `chr` ↔ bare name mapping in the loader, not in the schema.
10. Implement `pgr-load-vcf-annot`: reads `pgr variant annotate-vcf` output, populates
    `annot_gene_effect`.
11. Implement `pgr-load-bed-track`: generic BED loader for `annot_regions`.
12. Update `run_e2e.sh` to drive the full Phase 2 pipeline with a single
    `hg002.sampledb` as the terminal artefact.
