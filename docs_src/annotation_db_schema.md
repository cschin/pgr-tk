# PGR-TK Annotation Database Schema & Summary Report Guide

**Date:** 2026-04-07  
**Branch:** `aln_map_improvment`  
**Tools covered:** `pgr align alnmap`, `pgr align liftover-gtf`, `generate_e2e_report.py`, `get_tx_seqs.sh`

---

## 1. Alignment Database (`.alndb`)

Produced by **`pgr align alnmap`**.  One `.alndb` per haplotype (e.g., `hg002_hap0.alndb`).

### `run_params`
| Column | Type | Description |
|--------|------|-------------|
| `key`   | TEXT PK | Parameter name |
| `value` | TEXT    | Parameter value |

Stores run metadata: reference name, query name, preset, version, etc.

---

### `sequences`
| Column | Type | Description |
|--------|------|-------------|
| `seq_id`   | INTEGER PK | Internal ID |
| `seq_name` | TEXT       | Sequence name |
| `seq_type` | TEXT       | `"target"` (reference) or `"query"` (assembly contig) |
| `length`   | INTEGER    | Sequence length (bp) |

---

### `chains`
Top-level alignment records between one reference chromosome and one query contig.

| Column | Type | Description |
|--------|------|-------------|
| `aln_idx`        | INTEGER PK | Alignment index |
| `target_name`    | TEXT       | Reference chromosome |
| `target_start`   | INTEGER    | Start on reference (0-based) |
| `target_end`     | INTEGER    | End on reference |
| `query_name`     | TEXT       | Query contig name |
| `query_start`    | INTEGER    | Start on contig (0-based) |
| `query_end`      | INTEGER    | End on contig |
| `orientation`    | INTEGER    | 0 = forward, 1 = reverse |
| `ctg_orientation`| INTEGER    | Contig orientation flag |
| `query_length`   | INTEGER    | Full contig length |
| `target_dup`     | INTEGER    | 1 if this alignment is a duplicate on the reference |
| `target_ovlp`    | INTEGER    | 1 if this alignment overlaps a previous one on the reference |
| `query_dup`      | INTEGER    | 1 if this alignment is a duplicate on the contig |
| `query_ovlp`     | INTEGER    | 1 if this alignment overlaps a previous one on the contig |

---

### `blocks`
Individual aligned segments within a chain (M-blocks = matching regions).

| Column | Type | Description |
|--------|------|-------------|
| `block_id`       | INTEGER PK | Block ID |
| `aln_idx`        | INTEGER FK | → `chains.aln_idx` |
| `block_type`     | INTEGER    | Block classification |
| `dup_flag`       | INTEGER    | Duplication flag |
| `ovlp_flag`      | INTEGER    | Overlap flag |
| `target_name`    | TEXT       | Reference chromosome |
| `target_start`   | INTEGER    | Start on reference |
| `target_end`     | INTEGER    | End on reference |
| `query_name`     | TEXT       | Contig name |
| `query_start`    | INTEGER    | Start on contig |
| `query_end`      | INTEGER    | End on contig |
| `orientation`    | INTEGER    | Strand |
| `ctg_orientation`| INTEGER    | Contig orientation |
| `sv_diff_type`   | INTEGER    | SV classification (0=match, 1=insertion, 2=deletion, …) |

---

### `variants`
Small variants (SNPs and short indels) discovered from the alignment.

| Column | Type | Description |
|--------|------|-------------|
| `variant_id`   | INTEGER PK | |
| `aln_idx`      | INTEGER FK | → `chains.aln_idx` |
| `dup_flag`     | INTEGER    | |
| `ovlp_flag`    | INTEGER    | |
| `target_name`  | TEXT       | Reference chromosome |
| `target_start` | INTEGER    | Position on reference |
| `target_end`   | INTEGER    | |
| `query_name`   | TEXT       | Contig |
| `query_start`  | INTEGER    | Position on contig |
| `query_end`    | INTEGER    | |
| `orientation`  | INTEGER    | Strand |
| `target_diff`  | INTEGER    | Size change on reference |
| `query_diff`   | INTEGER    | Size change on contig |
| `target_coord` | INTEGER    | 0-based reference coordinate |
| `variant_type` | INTEGER    | SNP / indel code |
| `ref_seq`      | TEXT       | Reference allele |
| `alt_seq`      | TEXT       | Alternate allele |

---

### `ctgmap`
Contig-level alignment summary (one row per contig placed on the reference).

| Column | Type | Description |
|--------|------|-------------|
| `ctgmap_id`      | INTEGER PK | |
| `target_name`    | TEXT       | Reference chromosome |
| `target_start`   | INTEGER    | Aligned span start |
| `target_end`     | INTEGER    | Aligned span end |
| `query_name`     | TEXT       | Contig name |
| `query_start`    | INTEGER    | |
| `query_end`      | INTEGER    | |
| `ctg_len`        | INTEGER    | Full contig length |
| `orientation`    | INTEGER    | Strand |
| `ctg_orientation`| INTEGER    | |
| `target_dup`     | INTEGER    | Duplication flags |
| `target_ovlp`    | INTEGER    | Overlap flags |
| `query_dup`      | INTEGER    | |
| `query_ovlp`     | INTEGER    | |

---

### `sv_candidates`
Structural variant candidates on the **reference** coordinate axis (exported to `.svcnd.bed`).

| Column | Type | Description |
|--------|------|-------------|
| `svcnd_id`     | INTEGER PK | |
| `target_name`  | TEXT       | Reference chromosome (PanSN format, e.g. `GRCh38#0#chr1`) |
| `target_start` | INTEGER    | SV start on reference |
| `target_end`   | INTEGER    | SV end on reference |
| `sv_type`      | TEXT       | Type code — see table below |

**SV type codes (ref view):**

| Code | Meaning |
|------|---------|
| `SV` | Size-discordant gap — true structural variant |
| `TG` | Reference gap — unaligned region between consecutive contig alignments |
| `TD` | Contig duplicate — a contig alignment fully contained within a previous one |
| `TO` | Contig overlap — a contig alignment partially overlapping a previous one |

---

### `ctgsv`
Structural variant candidates on the **contig** coordinate axis (exported to `.ctgsv.bed`).

| Column | Type | Description |
|--------|------|-------------|
| `ctgsv_id`    | INTEGER PK | |
| `query_name`  | TEXT       | Contig name |
| `query_start` | INTEGER    | SV start on contig |
| `query_end`   | INTEGER    | SV end on contig |
| `sv_type`     | TEXT       | Type code — see table below |

**SV type codes (contig view):**

| Code | Meaning |
|------|---------|
| `QG` | Contig gap — unaligned region between consecutive reference alignments |
| `QD` | Reference duplicate — a reference block fully contained within a previous one |
| `QO` | Reference overlap — a reference block partially overlapping a previous one |

---

### `sv_sequences`
Reference and query sequences flanking each SV candidate.

| Column | Type | Description |
|--------|------|-------------|
| `svsq_id`      | INTEGER PK | |
| `target_name`  | TEXT       | Reference chromosome |
| `target_start` | INTEGER    | |
| `target_end`   | INTEGER    | |
| `query_name`   | TEXT       | Contig |
| `query_start`  | INTEGER    | |
| `query_end`    | INTEGER    | |
| `orientation`  | INTEGER    | Strand |
| `target_seq`   | TEXT       | Reference sequence around SV |
| `query_seq`    | TEXT       | Contig sequence around SV |

---

## 2. Liftover Database (`_liftover.db`)

Produced by **`pgr align liftover-gtf`**.  One database per haplotype.

---

### `transcripts`
One row per transcript from the input GTF (GRCh38 annotation).

| Column | Type | Description |
|--------|------|-------------|
| `transcript_pk` | INTEGER PK | Internal ID |
| `transcript_id` | TEXT UNIQUE | NCBI transcript accession (e.g. `NM_000001.3`) |
| `gene_id`       | TEXT       | NCBI gene ID |
| `gene_name`     | TEXT       | Human-readable gene symbol |
| `ref_chrom`     | TEXT       | Reference chromosome (`chr1`, `chr2`, …) |
| `ref_start`     | INTEGER    | Gene body start, 0-based |
| `ref_end`       | INTEGER    | Gene body end |
| `ref_strand`    | TEXT       | `"+"` or `"-"` |
| `exon_count`    | INTEGER    | Number of exons |

---

### `exons`
One row per exon for every transcript.

| Column | Type | Description |
|--------|------|-------------|
| `exon_pk`       | INTEGER PK | |
| `transcript_pk` | INTEGER FK | → `transcripts.transcript_pk` |
| `exon_start`    | INTEGER    | Exon start on reference, 0-based |
| `exon_end`      | INTEGER    | Exon end on reference |

---

### `liftover`
One row per transcript–contig alignment hit (a transcript may hit multiple contigs).

| Column | Type | Description |
|--------|------|-------------|
| `liftover_pk`   | INTEGER PK | |
| `transcript_pk` | INTEGER FK | → `transcripts.transcript_pk` |
| `contig`        | TEXT       | Query contig name |
| `contig_start`  | INTEGER    | Mapped span start on contig, 0-based |
| `contig_end`    | INTEGER    | Mapped span end on contig |
| `strand`        | TEXT       | Strand on contig |
| `coverage_pct`  | REAL       | Fraction of exon bases covered (0–100) |
| `block_count`   | INTEGER    | Number of mapped exon blocks |
| `block_sizes`   | TEXT       | Comma-separated exon block sizes (bp) |
| `block_starts`  | TEXT       | Comma-separated block offsets from `contig_start` |

---

### `transcript_summary`
Best-hit classification for each transcript (one row per transcript).

| Column | Type | Description |
|--------|------|-------------|
| `transcript_pk` | INTEGER PK FK | → `transcripts.transcript_pk` |
| `status`        | TEXT          | See status codes below |
| `hit_count`     | INTEGER       | Alignments above min-coverage threshold |
| `contig_count`  | INTEGER       | Number of distinct contigs with hits |
| `best_coverage` | REAL          | Coverage % of best hit |
| `best_contig`   | TEXT          | Name of best-matching contig |

**Status codes:**

| Status | Meaning |
|--------|---------|
| `single_full`    | Single hit ≥ 90% exon coverage (high-quality / HQ) |
| `single_partial` | Single hit ≥ 50%, < 90% coverage |
| `multi_contig`   | Best hit spans multiple contigs (fragmented assembly) |
| `multi_location` | Multiple hits on different contigs (segmental duplication) |
| `no_hit`         | No alignment above min-coverage threshold |

---

### `gene_summary`
Aggregated per-gene liftover status.

| Column | Type | Description |
|--------|------|-------------|
| `gene_id`       | TEXT PK | NCBI gene ID |
| `gene_name`     | TEXT    | Gene symbol |
| `ref_chrom`     | TEXT    | Reference chromosome |
| `total`         | INTEGER | Total transcripts for this gene |
| `single_full`   | INTEGER | |
| `single_partial`| INTEGER | |
| `multi_contig`  | INTEGER | |
| `multi_location`| INTEGER | |
| `no_hit`        | INTEGER | |

---

### `ref_sequences` *(optional — written with `--ref-fa`)*
| Column | Type | Description |
|--------|------|-------------|
| `transcript_pk` | INTEGER PK FK | |
| `ref_seq`       | TEXT          | Spliced exon sequence from GRCh38 |

---

### `contig_sequences` *(optional — written with `--query-fa`)*
| Column | Type | Description |
|--------|------|-------------|
| `liftover_pk` | INTEGER PK FK | |
| `contig_seq`  | TEXT          | Spliced exon sequence mapped from assembly contig |

---

## 3. SV Gene Annotation (`generate_e2e_report.py`)

The report script annotates each SV candidate against the gene/exon index loaded from the liftover database.

### Algorithm

1. **Load index** from `_liftover.db`:
   - `gene_index[chrom]` — sorted list of `(g_start, g_end, gene_name, strand)` per chromosome
   - `exon_index[(chrom, gene_name)]` — sorted list of `(exon_start, exon_end)` per gene
   - Primary chromosomes only (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`); alt/fix/unplaced excluded

2. **Overlap query** — binary search (`bisect`) on sorted gene starts:
   - If ≥1 gene overlaps → classify as **genic**
   - Otherwise → classify as **intergenic**

3. **Exon impact** (genic SVs):

   ```
   exon disrupted  SV:  |══════════════════════|
                   Exon:    |────────|            ← exon fully inside SV
                   → entire exon deleted/replaced; reading frame likely broken

   exon partial    SV:  |══════════════|
                   Exon:          |────────────|  ← exon only partly overlapped
                   → splice site may be disrupted; partial coding-sequence change

   intronic        SV overlaps gene body but no exon
                   → minimal coding impact; possible regulatory effect
   ```

4. **Nearest-neighbor** (intergenic SVs):
   - `upstream` — closest gene whose end ≤ SV start (5′ neighbor)
   - `downstream` — closest gene whose start ≥ SV end (3′ neighbor)
   - Reports gene name, strand, and distance (bp)

### Per-SV annotation record

| Field | Description |
|-------|-------------|
| `chrom`          | Chromosome (prefix-stripped, e.g. `chr1`) |
| `start`, `end`   | SV coordinates |
| `sv_type`        | `SV`, `TG`, `TD`, `TO` |
| `context`        | `"genic"` or `"intergenic"` |
| `impact`         | `"exon disrupted"`, `"exon partial"`, `"intronic"`, or `None` |
| `gene`           | Overlapping gene symbol (genic only) |
| `strand`         | Gene strand `"+"` or `"-"` (genic only) |
| `n_exon_full`    | Count of exons fully contained in SV |
| `n_exon_partial` | Count of exons partially overlapped |
| `upstream`       | `(gene, distance_bp, strand)` or `None` |
| `downstream`     | `(gene, distance_bp, strand)` or `None` |

---

## 4. Generated Report Files

| File | Produced by | Description |
|------|-------------|-------------|
| `hg002_hap{0,1}.alndb`              | `pgr align alnmap`         | Alignment database |
| `hg002_hap{0,1}.svcnd.bed`          | `pgr align alnmap`         | SV candidates, ref coordinates |
| `hg002_hap{0,1}.ctgsv.bed`          | `pgr align alnmap`         | SV candidates, contig coordinates |
| `hg002_hap{0,1}_liftover.db`        | `pgr align liftover-gtf`   | Transcript liftover database |
| `hg002_hap{0,1}.hq.gtf`             | `pgr align liftover-gtf`   | High-quality transcript annotations |
| `hg002_hap{0,1}_ref_tx.fa`          | `pgr align liftover-gtf`   | Reference transcript FASTA |
| `hg002_hap{0,1}_hq_contig_tx.fa`    | `pgr align liftover-gtf`   | HQ contig-mapped transcript FASTA |
| `hg002.annotated.sorted.clinvar.vcf.gz` | bcftools pipeline | VCF with ClinVar annotations |
| `liftover_report.html`              | `get_tx_seqs.sh`     | Transcript liftover HTML report |
| `e2e_report.html`                   | `generate_e2e_report.py` | Full end-to-end HTML report |
| `e2e_report_lite.html`              | `generate_e2e_report.py` | Lite report (≤5,000 scatter points) |

---

## 5. ClinVar Annotation Pipeline

### Steps

1. **`pgr align alnmap`** — produces `hg002_hap{0,1}.alndb`
2. **`pgr variant diploid-vcf`** — merges hap0/hap1 alndb into `hg002.vcf`
3. Strip PanSN prefix from CHROM column (`GRCh38#0#chr1` → `chr1`)
4. **`pgr variant annotate-vcf`** — annotates with `hg38.ncbiRefSeq.gtf.gz` → `hg002.annotated.vcf`
5. **bgzip / tabix** — compress and index
6. **bcftools sort + rename** — sort, strip `chr` prefix for ClinVar matching → `hg002.nochr.vcf.gz`
7. **bcftools annotate** — add `CLNSIG`, `CLNDN`, `CLNREVSTAT` from `clinvar.vcf.gz` → `hg002.annotated.sorted.clinvar.vcf.gz`

### Report fields extracted per variant

| Field | Source | Description |
|-------|--------|-------------|
| `CHROM`, `POS` | VCF | Chromosomal position |
| `REF`, `ALT`   | VCF | Alleles (indels truncated to 20 chars + `…`) |
| `TYPE`         | VCF INFO | `snp`, `indel`, etc. |
| `CLNSIG`       | ClinVar | Significance (e.g. `Pathogenic`, `Likely_benign`) |
| `CLNDN`        | ClinVar | Associated disease/phenotype |
| `CLNREVSTAT`   | ClinVar | Review status (number of stars) |
| `GN`           | VCF INFO | Gene name from GTF annotation |
| `GT`           | VCF FORMAT | Genotype → zygosity (het / hom) |

### Significance categories

| Category | Criteria |
|----------|----------|
| Pathogenic / Likely pathogenic | Contains `Pathogenic` or `Likely_pathogenic` and does NOT contain `Benign` |
| Uncertain significance | Contains `Uncertain_significance` |
| Benign / Likely benign | Contains `Benign` or `Likely_benign` |
| Drug response / risk factor | Contains `drug_response`, `association`, `risk_factor`, `protective`, or `Likely_risk_allele` |
| Other | Everything else |
