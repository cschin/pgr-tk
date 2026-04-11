use clap::{self, Parser};
use flate2::bufread::MultiGzDecoder;
use iset::IntervalMap;
use rusqlite::params;
use rusqlite::Connection;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use pgr_db::fasta_io::reverse_complement;

/// Lift over GTF transcript annotations from a reference (e.g. GRCh38) to haplotype
/// contigs using an alndb alignment file produced by pgr-alnmap.
///
/// The reference is the *target* side of the alignment (what pgr-alnmap was run against),
/// and the haplotype contigs are the *query* side.
///
/// Output: one TSV line per successfully lifted transcript (BED12-like with extra columns).
#[derive(Parser, Debug)]
#[clap(about, long_about = None)]
pub struct Args {
    /// path to the alndb file (reference → haplotype contig alignment)
    #[clap(long, short)]
    pub alndb_path: String,

    /// path to the GTF annotation file (plain or gzip-compressed)
    #[clap(long, short)]
    pub gtf_path: String,

    /// path for the output SQLite database
    #[clap(long, short)]
    pub output_db: String,

    /// minimum fraction [0.0–1.0] of exon bases that must be covered by M-blocks on
    /// a single contig for the transcript to be reported
    #[clap(long, default_value_t = 0.5)]
    pub min_coverage: f64,

    /// strip this prefix string from target chromosome names in the alndb before
    /// matching against GTF chromosome names (e.g. "GRCh38#0#" to convert PanSN
    /// names to plain chromosome names).  Leave unset to use names as-is.
    #[clap(long)]
    pub target_chr_prefix: Option<String>,

    /// coverage threshold [0.0–1.0] above which a single hit is classified as
    /// "full" rather than "partial" in the anomaly summary (default 0.9)
    #[clap(long, default_value_t = 0.9)]
    pub full_coverage: f64,

    /// path to the reference FASTA file (plain or .gz) containing the sequences
    /// used as the alignment target (e.g. GRCh38).  When provided together with
    /// --query-fa, spliced transcript sequences are written to
    /// <output>.ref_tx.fa / <output>.contig_tx.fa / <output>.hq_contig_tx.fa
    /// and stored in the ref_sequences / contig_sequences tables.
    #[clap(long)]
    pub ref_fa: Option<String>,

    /// path to the haplotype FASTA file (plain or .gz) containing the query
    /// contig sequences (e.g. hg002_hap0).
    #[clap(long)]
    pub query_fa: Option<String>,
}

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Per-block data stored in the inner interval map: (qs, qe, orientation)
type QBlock = (u32, u32, u32);

/// target_chrom → query_name → IntervalMap<QBlock>
///
/// Keeping a separate IntervalMap per query contig prevents collisions when two
/// different contigs produce M-blocks with identical (target_start, target_end)
/// coordinates (common in repetitive / paralogous regions).
type TargetIntervals = FxHashMap<String, FxHashMap<String, IntervalMap<u32, QBlock>>>;

struct TranscriptInfo {
    transcript_id: String,
    gene_id: String,
    gene_name: String,
    chrom: String,
    start: u32, // 0-based
    end: u32,   // 0-based exclusive
    strand: char,
    exons: Vec<(u32, u32)>, // (start, end) 0-based exclusive, sorted
}

struct LiftoverResult {
    query_name: String,
    query_start: u32,
    query_end: u32,
    strand: char,
    block_count: usize,
    block_sizes: Vec<u32>,
    block_starts: Vec<u32>, // relative to query_start
    coverage_pct: f64,
}

// ---------------------------------------------------------------------------
// Load M-blocks from alndb
// ---------------------------------------------------------------------------

fn load_target_intervals(db_path: &str, strip_prefix: Option<&str>) -> TargetIntervals {
    let conn = Connection::open_with_flags(db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .expect("can't open alndb file");

    // ── Pass 1: collect all M-blocks grouped by chain ──────────────────────
    // Key: (target_name, query_name, aln_idx, orientation)
    // Value: Vec<(ts, te, qs, qe)>
    type ChainKey = (String, String, u32, u32);
    let mut chain_blocks: FxHashMap<ChainKey, Vec<(u32, u32, u32, u32)>> = FxHashMap::default();

    let mut stmt = conn
        .prepare(
            "SELECT target_name, target_start, target_end,
                    query_name, query_start, query_end, orientation, aln_idx
             FROM blocks WHERE block_type = 0",
        )
        .expect("prepare blocks query");

    let mut rows = stmt.query([]).expect("query blocks");
    while let Some(row) = rows.next().expect("blocks row") {
        let raw_t_name: String = row.get(0).unwrap();
        let t_name = match strip_prefix {
            Some(pfx) => raw_t_name
                .strip_prefix(pfx)
                .unwrap_or(&raw_t_name)
                .to_string(),
            None => raw_t_name,
        };
        let ts: u32 = row.get(1).unwrap();
        let te: u32 = row.get(2).unwrap();
        let q_name: String = row.get(3).unwrap();
        let qs: u32 = row.get(4).unwrap();
        let qe: u32 = row.get(5).unwrap();
        let orientation: u32 = row.get::<_, i64>(6).unwrap() as u32;
        let aln_idx: u32 = row.get::<_, i64>(7).unwrap() as u32;
        if te > ts {
            chain_blocks
                .entry((t_name, q_name, aln_idx, orientation))
                .or_default()
                .push((ts, te, qs, qe));
        }
    }

    // ── Pass 2: for each chain, merge overlapping M-blocks then fill gaps ──
    let mut target_intervals = TargetIntervals::default();

    for ((t_name, q_name, _aln_idx, orientation), mut blocks) in chain_blocks {
        // Sort by target_start
        blocks.sort_unstable_by_key(|&(ts, _, _, _)| ts);

        // Merge overlapping / adjacent M-blocks.
        // For orientation=0: qs increases with ts, so take min(qs) / max(qe).
        // For orientation=1: qs decreases with ts, so take max(qs) / min(qe).
        let mut merged: Vec<(u32, u32, u32, u32)> = Vec::new(); // (ts, te, qs, qe)
        for (ts, te, qs, qe) in blocks {
            if let Some(last) = merged.last_mut() {
                if ts < last.1 {
                    // Overlapping: extend
                    last.1 = last.1.max(te);
                    if orientation == 0 {
                        last.2 = last.2.min(qs);
                        last.3 = last.3.max(qe);
                    } else {
                        last.2 = last.2.max(qs);
                        last.3 = last.3.min(qe);
                    }
                    continue;
                }
            }
            merged.push((ts, te, qs, qe));
        }

        let qmap = target_intervals
            .entry(t_name)
            .or_default()
            .entry(q_name)
            .or_default();

        // Insert merged M-blocks
        for &(ts, te, qs, qe) in &merged {
            if te > ts {
                qmap.insert(ts..te, (qs, qe, orientation));
            }
        }

        // Fill intra-chain gaps between consecutive merged blocks.
        // The gap block carries the bounding query coordinates so that
        // map_interval() can project exon positions linearly within the gap.
        // This lets exons that fall in locally divergent (but chain-aligned)
        // regions receive a coordinate mapping even without an M-block anchor.
        for w in merged.windows(2) {
            let (_, te_prev, _, qe_prev) = w[0];
            let (ts_next, _, qs_next, _) = w[1];
            if ts_next > te_prev {
                // Gap exists on the target side
                let (gap_qs, gap_qe) = if orientation == 0 {
                    // forward: query gap is [qe_prev, qs_next)
                    (qe_prev, qs_next)
                } else {
                    // reverse: query gap is [qe_next_block, qs_prev_block)
                    // stored as w[1].qe .. w[0].qs  (lower .. higher query)
                    (w[1].3, w[0].2)
                };
                if gap_qe > gap_qs && ts_next > te_prev {
                    qmap.insert(te_prev..ts_next, (gap_qs, gap_qe, orientation));
                }
            }
        }
    }

    target_intervals
}

// ---------------------------------------------------------------------------
// GTF parsing
// ---------------------------------------------------------------------------

fn parse_attribute(attributes: &str, key: &str) -> Option<String> {
    // Matches: key "value"; or key "value"
    for field in attributes.split(';') {
        let field = field.trim();
        if let Some(rest) = field.strip_prefix(key) {
            let rest = rest.trim();
            // rest should be  "value"
            let value = rest.trim_matches('"');
            if !value.is_empty() {
                return Some(value.to_string());
            }
            // handle: key "value"  (space separated)
            if let Some(idx) = rest.find('"') {
                let inner = &rest[idx + 1..];
                if let Some(end) = inner.find('"') {
                    return Some(inner[..end].to_string());
                }
            }
        }
    }
    None
}

fn parse_gtf<R: BufRead>(reader: R) -> Vec<TranscriptInfo> {
    // Use a Vec to preserve GTF order and a separate index map for O(1) lookup.
    let mut transcripts: Vec<TranscriptInfo> = Vec::new();
    let mut tx_index = FxHashMap::<String, usize>::default(); // transcript_id → index in transcripts

    for line in reader.lines() {
        let Ok(line) = line else { continue };
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        let chrom = fields[0];
        let feature = fields[2];
        let start_1based: u32 = match fields[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end_1based: u32 = match fields[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        // Convert from GTF 1-based inclusive to 0-based half-open
        let start = start_1based - 1;
        let end = end_1based;
        let strand = fields[6].chars().next().unwrap_or('.');
        let attributes = fields[8];

        match feature {
            "transcript" => {
                let transcript_id = parse_attribute(attributes, "transcript_id")
                    .unwrap_or_else(|| "unknown_tx".to_string());
                if tx_index.contains_key(&transcript_id) {
                    continue; // deduplicate
                }
                let gene_id = parse_attribute(attributes, "gene_id")
                    .unwrap_or_else(|| "unknown_gene".to_string());
                let gene_name = parse_attribute(attributes, "gene_name")
                    .or_else(|| parse_attribute(attributes, "gene"))
                    .unwrap_or_else(|| gene_id.clone());
                let idx = transcripts.len();
                tx_index.insert(transcript_id.clone(), idx);
                transcripts.push(TranscriptInfo {
                    transcript_id,
                    gene_id,
                    gene_name,
                    chrom: chrom.to_string(),
                    start,
                    end,
                    strand,
                    exons: Vec::new(),
                });
            }
            "exon" => {
                let transcript_id = match parse_attribute(attributes, "transcript_id") {
                    Some(id) => id,
                    None => continue,
                };
                let idx = if let Some(&i) = tx_index.get(&transcript_id) {
                    i
                } else {
                    // exon appeared before its transcript record — create a stub
                    let gene_id = parse_attribute(attributes, "gene_id")
                        .unwrap_or_else(|| "unknown_gene".to_string());
                    let gene_name = parse_attribute(attributes, "gene_name")
                        .or_else(|| parse_attribute(attributes, "gene"))
                        .unwrap_or_else(|| gene_id.clone());
                    let i = transcripts.len();
                    tx_index.insert(transcript_id.clone(), i);
                    transcripts.push(TranscriptInfo {
                        transcript_id,
                        gene_id,
                        gene_name,
                        chrom: chrom.to_string(),
                        start,
                        end,
                        strand,
                        exons: Vec::new(),
                    });
                    i
                };
                transcripts[idx].exons.push((start, end));
            }
            _ => {}
        }
    }

    // For transcripts with no explicit exon records use the transcript span.
    for tx in &mut transcripts {
        if tx.exons.is_empty() {
            tx.exons.push((tx.start, tx.end));
        } else {
            tx.exons.sort_unstable();
            tx.exons.dedup();
        }
    }

    transcripts
}

// ---------------------------------------------------------------------------
// Coordinate mapping helpers
// ---------------------------------------------------------------------------

/// Map a half-open target interval [exon_s, exon_e) that overlaps block [ts, te)
/// to query coordinates.  Returns (qs_mapped, qe_mapped) in 0-based half-open.
fn map_interval(
    exon_s: u32,
    exon_e: u32,
    ts: u32,
    te: u32,
    qs: u32,
    _qe: u32,
    orientation: u32,
) -> (u32, u32) {
    let s = exon_s.max(ts);
    let e = exon_e.min(te);
    if s >= e {
        return (0, 0);
    }
    if orientation == 0 {
        (qs + (s - ts), qs + (e - ts))
    } else {
        // reverse: ts↔qe, te↔qs
        (qs + (te - e), qs + (te - s))
    }
}

fn flip_strand(s: char) -> char {
    match s {
        '+' => '-',
        '-' => '+',
        c => c,
    }
}

// ---------------------------------------------------------------------------
// Liftover logic for one transcript
// ---------------------------------------------------------------------------

fn liftover_transcript(
    tx: &TranscriptInfo,
    target_intervals: &TargetIntervals,
    min_coverage: f64,
) -> Vec<LiftoverResult> {
    let query_maps = match target_intervals.get(&tx.chrom) {
        Some(m) => m,
        None => return Vec::new(),
    };

    let total_exon_bases: u32 = tx.exons.iter().map(|(s, e)| e - s).sum();
    if total_exon_bases == 0 {
        return Vec::new();
    }

    // Collect all mapped segments grouped by (query_name, orientation)
    // key: (query_name, orientation)  value: Vec<(query_start, query_end)>
    let mut groups: FxHashMap<(String, u32), Vec<(u32, u32)>> = FxHashMap::default();

    for &(exon_s, exon_e) in &tx.exons {
        for (q_name, qmap) in query_maps {
            for (range, (qs, qe, orientation)) in qmap.iter(exon_s..exon_e) {
                let ts = range.start;
                let te = range.end;
                let (ms, me) = map_interval(exon_s, exon_e, ts, te, *qs, *qe, *orientation);
                if me > ms {
                    groups
                        .entry((q_name.clone(), *orientation))
                        .or_default()
                        .push((ms, me));
                }
            }
        }
    }

    if groups.is_empty() {
        return Vec::new();
    }

    // For each group, merge overlapping/adjacent segments, compute covered bases
    let mut results = Vec::new();

    for ((q_name, orientation), mut segs) in groups {
        segs.sort_unstable();

        // Merge overlapping / adjacent intervals
        let mut merged: Vec<(u32, u32)> = Vec::new();
        for (s, e) in segs {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }

        let covered_bases: u32 = merged.iter().map(|(s, e)| e - s).sum();
        let coverage = covered_bases as f64 / total_exon_bases as f64;
        if coverage < min_coverage {
            continue;
        }

        let query_start = merged.first().unwrap().0;
        let query_end = merged.last().unwrap().1;

        let strand = if orientation == 0 {
            tx.strand
        } else {
            flip_strand(tx.strand)
        };

        // For orientation=1 the segments are already in ascending order on the contig,
        // but they represent the reverse-complement of the transcript.  In BED12 the
        // blocks are always listed left-to-right on the chromosome regardless of strand.
        let block_count = merged.len();
        let block_sizes: Vec<u32> = merged.iter().map(|(s, e)| e - s).collect();
        let block_starts: Vec<u32> = merged.iter().map(|(s, _)| s - query_start).collect();

        results.push(LiftoverResult {
            query_name: q_name,
            query_start,
            query_end,
            strand,
            block_count,
            block_sizes,
            block_starts,
            coverage_pct: coverage * 100.0,
        });
    }

    // Sort results by coverage descending so the best hit comes first
    results.sort_by(|a, b| b.coverage_pct.partial_cmp(&a.coverage_pct).unwrap());
    results
}

// ---------------------------------------------------------------------------
// Output database schema
// ---------------------------------------------------------------------------

fn init_output_db(conn: &Connection) {
    conn.execute_batch(
        "PRAGMA journal_mode = WAL;
         PRAGMA foreign_keys = ON;

         CREATE TABLE IF NOT EXISTS transcripts (
             transcript_pk  INTEGER PRIMARY KEY,
             transcript_id  TEXT    NOT NULL UNIQUE,
             gene_id        TEXT    NOT NULL,
             gene_name      TEXT    NOT NULL,
             ref_chrom      TEXT    NOT NULL,
             ref_start      INTEGER NOT NULL,
             ref_end        INTEGER NOT NULL,
             ref_strand     TEXT    NOT NULL,
             exon_count     INTEGER NOT NULL
         );

         CREATE TABLE IF NOT EXISTS exons (
             exon_pk        INTEGER PRIMARY KEY,
             transcript_pk  INTEGER NOT NULL
                                REFERENCES transcripts(transcript_pk)
                                ON DELETE CASCADE,
             exon_start     INTEGER NOT NULL,
             exon_end       INTEGER NOT NULL
         );

         CREATE TABLE IF NOT EXISTS liftover (
             liftover_pk    INTEGER PRIMARY KEY,
             transcript_pk  INTEGER NOT NULL
                                REFERENCES transcripts(transcript_pk)
                                ON DELETE CASCADE,
             contig         TEXT    NOT NULL,
             contig_start   INTEGER NOT NULL,
             contig_end     INTEGER NOT NULL,
             strand         TEXT    NOT NULL,
             coverage_pct   REAL    NOT NULL,
             block_count    INTEGER NOT NULL,
             block_sizes    TEXT    NOT NULL,
             block_starts   TEXT    NOT NULL
         );

         CREATE INDEX IF NOT EXISTS idx_transcripts_chrom
             ON transcripts(ref_chrom, ref_start, ref_end);
         CREATE INDEX IF NOT EXISTS idx_exons_tx
             ON exons(transcript_pk);
         CREATE INDEX IF NOT EXISTS idx_liftover_tx
             ON liftover(transcript_pk);
         CREATE INDEX IF NOT EXISTS idx_liftover_contig
             ON liftover(contig, contig_start, contig_end);

         CREATE TABLE IF NOT EXISTS transcript_summary (
             transcript_pk  INTEGER PRIMARY KEY
                                REFERENCES transcripts(transcript_pk)
                                ON DELETE CASCADE,
             status         TEXT    NOT NULL,
             hit_count      INTEGER NOT NULL,
             contig_count   INTEGER NOT NULL,
             best_coverage  REAL    NOT NULL,
             best_contig    TEXT    NOT NULL
         );

         CREATE INDEX IF NOT EXISTS idx_summary_status
             ON transcript_summary(status);

         CREATE TABLE IF NOT EXISTS gene_summary (
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

         CREATE INDEX IF NOT EXISTS idx_gene_summary_chrom
             ON gene_summary(ref_chrom);

         CREATE TABLE IF NOT EXISTS ref_sequences (
             transcript_pk  INTEGER PRIMARY KEY
                                REFERENCES transcripts(transcript_pk)
                                ON DELETE CASCADE,
             ref_seq        TEXT NOT NULL
         );

         CREATE TABLE IF NOT EXISTS contig_sequences (
             liftover_pk    INTEGER PRIMARY KEY
                                REFERENCES liftover(liftover_pk)
                                ON DELETE CASCADE,
             contig_seq     TEXT NOT NULL
         );",
    )
    .expect("failed to initialise output database schema");
}

// ---------------------------------------------------------------------------
// Sequence loading and fetch helpers (direct FASTA)
// ---------------------------------------------------------------------------

/// Load all sequences from a FASTA file (plain or .gz) into a HashMap keyed
/// by the first whitespace-delimited token of the header line.
fn load_fasta_seqs(path: &str) -> std::io::Result<FxHashMap<String, Vec<u8>>> {
    use pgr_db::fasta_io::FastaReader;
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut fa = FastaReader::new(reader, &path.to_string(), 1 << 24, false, false)?;
    let mut map: FxHashMap<String, Vec<u8>> = FxHashMap::default();
    while let Some(Ok(rec)) = fa.next_rec() {
        let name = String::from_utf8_lossy(&rec.id).to_string();
        map.insert(name, rec.seq);
    }
    Ok(map)
}

/// Write a FASTA record with 80-column line wrapping.
fn write_fasta_record(out: &mut BufWriter<File>, header: &str, seq: &[u8]) -> std::io::Result<()> {
    writeln!(out, ">{}", header)?;
    for chunk in seq.chunks(80) {
        out.write_all(chunk)?;
        writeln!(out)?;
    }
    Ok(())
}

/// Splice reference exon sequences from a loaded FASTA map and apply strand.
/// `chrom` may use PanSN format ("GRCh38#0#chr1"); the map is tried with the
/// full name first, then with just the last '#'-delimited field.
fn fetch_ref_tx_seq(
    seqs: &FxHashMap<String, Vec<u8>>,
    chrom: &str,
    exons: &[(u32, u32)], // 0-based half-open, sorted ascending
    strand: char,
) -> Vec<u8> {
    let seq = seqs.get(chrom).or_else(|| {
        // Try bare chromosome name (last field after '#')
        let bare = chrom.split('#').next_back().unwrap_or(chrom);
        seqs.get(bare)
    });
    let seq = match seq {
        Some(s) => s,
        None => return Vec::new(),
    };
    let mut spliced: Vec<u8> = Vec::new();
    for &(es, ee) in exons {
        let es = es as usize;
        let ee = (ee as usize).min(seq.len());
        if es < ee {
            spliced.extend_from_slice(&seq[es..ee]);
        }
    }
    if strand == '-' {
        reverse_complement(&spliced)
    } else {
        spliced
    }
}

/// Splice contig block sequences from a loaded FASTA map and apply strand.
fn fetch_contig_tx_seq(seqs: &FxHashMap<String, Vec<u8>>, hit: &LiftoverResult) -> Vec<u8> {
    let seq = match seqs.get(&hit.query_name) {
        Some(s) => s,
        None => return Vec::new(),
    };
    let mut spliced: Vec<u8> = Vec::new();
    for (&size, &rel) in hit.block_sizes.iter().zip(hit.block_starts.iter()) {
        let bs = (hit.query_start + rel) as usize;
        let be = (bs + size as usize).min(seq.len());
        if bs < be {
            spliced.extend_from_slice(&seq[bs..be]);
        }
    }
    if hit.strand == '-' {
        reverse_complement(&spliced)
    } else {
        spliced
    }
}

// ---------------------------------------------------------------------------
// Anomaly classification
// ---------------------------------------------------------------------------

/// Classification of a transcript's liftover result.
#[derive(Default)]
struct LiftoverSummary {
    /// Anomaly category string stored in the DB.
    status: &'static str,
    hit_count: usize,
    contig_count: usize,
    best_coverage: f64,  // 0.0 when no_hit
    best_contig: String, // empty when no_hit
}

/// Classify the set of hits for one transcript.
fn classify_hits(hits: &[LiftoverResult], full_cov_threshold: f64) -> LiftoverSummary {
    if hits.is_empty() {
        return LiftoverSummary {
            status: "no_hit",
            ..Default::default()
        };
    }

    // hits are already sorted coverage-desc by liftover_transcript()
    let best = &hits[0];
    let contig_count = hits
        .iter()
        .map(|h| h.query_name.as_str())
        .collect::<std::collections::HashSet<_>>()
        .len();

    let status = if hits.len() == 1 {
        if best.coverage_pct >= full_cov_threshold * 100.0 {
            "single_full"
        } else {
            "single_partial"
        }
    } else if contig_count > 1 {
        "multi_contig" // best hit spans or maps to different contigs
    } else {
        "multi_location" // multiple hits on the same contig (e.g. segmental dup)
    };

    LiftoverSummary {
        status,
        hit_count: hits.len(),
        contig_count,
        best_coverage: best.coverage_pct,
        best_contig: best.query_name.clone(),
    }
}

// ---------------------------------------------------------------------------
// GTF writing helper
// ---------------------------------------------------------------------------

fn write_gtf_records(
    out: &mut BufWriter<File>,
    tx: &TranscriptInfo,
    hit: &LiftoverResult,
) -> std::io::Result<()> {
    let attrs = format!(
        "gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; \
         ref_chrom \"{}\"; ref_start {}; ref_end {}; ref_strand \"{}\"; \
         coverage_pct {:.1};",
        tx.gene_id,
        tx.transcript_id,
        tx.gene_name,
        tx.chrom,
        tx.start,
        tx.end,
        tx.strand,
        hit.coverage_pct,
    );

    // transcript record — GTF is 1-based inclusive
    writeln!(
        out,
        "{}\tpgr-liftover\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
        hit.query_name,
        hit.query_start + 1,
        hit.query_end,
        hit.strand,
        attrs
    )?;

    // exon records — one per alignment block
    for (size, rel) in hit.block_sizes.iter().zip(hit.block_starts.iter()) {
        let exon_s = hit.query_start + rel;
        let exon_e = exon_s + size;
        writeln!(
            out,
            "{}\tpgr-liftover\texon\t{}\t{}\t.\t{}\t.\t{}",
            hit.query_name,
            exon_s + 1,
            exon_e,
            hit.strand,
            attrs
        )?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

pub fn run(args: Args) -> Result<(), std::io::Error> {
    // Load alignment blocks from alndb
    let target_intervals =
        load_target_intervals(&args.alndb_path, args.target_chr_prefix.as_deref());

    // Parse GTF (supports .gz)
    let transcripts: Vec<TranscriptInfo> = {
        let file = File::open(Path::new(&args.gtf_path))?;
        if args.gtf_path.ends_with(".gz") {
            parse_gtf(BufReader::new(MultiGzDecoder::new(BufReader::new(file))))
        } else {
            parse_gtf(BufReader::new(file))
        }
    };

    // Open / create output database
    let db_path = Path::new(&args.output_db);
    if db_path.exists() {
        std::fs::remove_file(db_path)?;
    }
    let conn = Connection::open(db_path).expect("failed to open output database");
    init_output_db(&conn);

    // Open output GTF (derive path by replacing / appending .gtf to the db path)
    let gtf_path = db_path.with_extension("gtf");
    let mut gtf_out = BufWriter::new(File::create(&gtf_path)?);
    // GTF format uses no version header — comment lines only
    writeln!(gtf_out, "# source-db: {}", args.output_db)?;
    writeln!(gtf_out, "# alndb: {}", args.alndb_path)?;
    writeln!(
        gtf_out,
        "# min_coverage_pct: {:.0}",
        args.min_coverage * 100.0
    )?;

    // High-quality GTF: only hits >= full_coverage threshold
    let hq_gtf_path = db_path.with_extension("hq.gtf");
    let mut hq_gtf_out = BufWriter::new(File::create(&hq_gtf_path)?);
    writeln!(hq_gtf_out, "# source-db: {}", args.output_db)?;
    writeln!(hq_gtf_out, "# alndb: {}", args.alndb_path)?;
    writeln!(
        hq_gtf_out,
        "# min_coverage_pct: {:.0}",
        args.full_coverage * 100.0
    )?;

    // --- Insert GTF data (transcripts + exons) in one transaction ---
    let tx = conn.unchecked_transaction().expect("begin transaction");
    let mut tx_stmt = tx
        .prepare(
            "INSERT INTO transcripts
                 (transcript_id, gene_id, gene_name, ref_chrom,
                  ref_start, ref_end, ref_strand, exon_count)
             VALUES (?1,?2,?3,?4,?5,?6,?7,?8)",
        )
        .expect("prepare transcripts insert");
    let mut ex_stmt = tx
        .prepare(
            "INSERT INTO exons (transcript_pk, exon_start, exon_end)
             VALUES (?1,?2,?3)",
        )
        .expect("prepare exons insert");

    // transcript_id → transcript_pk for the liftover phase
    let mut tx_pk: FxHashMap<String, i64> = FxHashMap::default();

    for t in &transcripts {
        tx_stmt
            .execute(params![
                t.transcript_id,
                t.gene_id,
                t.gene_name,
                t.chrom,
                t.start,
                t.end,
                t.strand.to_string(),
                t.exons.len() as i64,
            ])
            .expect("insert transcript");
        let pk = tx.last_insert_rowid();
        tx_pk.insert(t.transcript_id.clone(), pk);
        for &(es, ee) in &t.exons {
            ex_stmt.execute(params![pk, es, ee]).expect("insert exon");
        }
    }
    drop(tx_stmt);
    drop(ex_stmt);
    tx.commit().expect("commit GTF data");

    // --- Liftover + summary in one transaction ---
    let tx = conn
        .unchecked_transaction()
        .expect("begin liftover transaction");
    let mut lo_stmt = tx
        .prepare(
            "INSERT INTO liftover
                 (transcript_pk, contig, contig_start, contig_end,
                  strand, coverage_pct, block_count, block_sizes, block_starts)
             VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9)",
        )
        .expect("prepare liftover insert");
    let mut sum_stmt = tx
        .prepare(
            "INSERT INTO transcript_summary
                 (transcript_pk, status, hit_count, contig_count,
                  best_coverage, best_contig)
             VALUES (?1,?2,?3,?4,?5,?6)",
        )
        .expect("prepare summary insert");

    // Unique gene IDs with at least one high-quality liftover
    let mut hq_genes: FxHashSet<String> = FxHashSet::default();

    // Aggregate counts per status
    let mut cnt_no_hit: u64 = 0;
    let mut cnt_single_full: u64 = 0;
    let mut cnt_single_partial: u64 = 0;
    let mut cnt_multi_contig: u64 = 0;
    let mut cnt_multi_location: u64 = 0;

    // Per-gene accumulators: gene_id → (gene_name, ref_chrom, [sf, sp, mc, ml, nh])
    let mut gene_counts: FxHashMap<String, (String, String, [u64; 5])> = FxHashMap::default();

    // Collect anomaly detail rows for the summary file
    // (transcript_id, gene_id, gene_name, ref_chrom, ref_start, ref_end, status,
    //  hit_count, contig_count, best_coverage, best_contig, contigs)
    let mut anomaly_rows: Vec<(
        String,
        String,
        String,
        String,
        u32,
        u32,
        &'static str,
        usize,
        usize,
        f64,
        String,
        String,
    )> = Vec::new();

    for t in &transcripts {
        let pk = tx_pk[&t.transcript_id];
        let hits = liftover_transcript(t, &target_intervals, args.min_coverage);

        let s = classify_hits(&hits, args.full_coverage);
        let status_idx = match s.status {
            "single_full" => {
                cnt_single_full += 1;
                0
            }
            "single_partial" => {
                cnt_single_partial += 1;
                1
            }
            "multi_contig" => {
                cnt_multi_contig += 1;
                2
            }
            "multi_location" => {
                cnt_multi_location += 1;
                3
            }
            _ => {
                cnt_no_hit += 1;
                4
            }
        };
        // Accumulate into gene-level counts
        let ge = gene_counts
            .entry(t.gene_id.clone())
            .or_insert_with(|| (t.gene_name.clone(), t.chrom.clone(), [0u64; 5]));
        ge.2[status_idx] += 1;

        // Collect anomaly rows (everything except single_full)
        if s.status != "single_full" {
            let contigs = hits
                .iter()
                .map(|h| h.query_name.as_str())
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
                .join(";");
            anomaly_rows.push((
                t.transcript_id.clone(),
                t.gene_id.clone(),
                t.gene_name.clone(),
                t.chrom.clone(),
                t.start,
                t.end,
                s.status,
                s.hit_count,
                s.contig_count,
                s.best_coverage,
                s.best_contig.clone(),
                contigs,
            ));
        }

        sum_stmt
            .execute(params![
                pk,
                s.status,
                s.hit_count as i64,
                s.contig_count as i64,
                s.best_coverage,
                s.best_contig,
            ])
            .expect("insert summary");

        // GTF: best hit only (above min_coverage)
        if let Some(best) = hits.first() {
            write_gtf_records(&mut gtf_out, t, best)?;
            // High-quality GTF: best hit must also meet full_coverage threshold
            if best.coverage_pct >= args.full_coverage * 100.0 {
                write_gtf_records(&mut hq_gtf_out, t, best)?;
                hq_genes.insert(t.gene_id.clone());
            }
        }

        // SQLite liftover: all hits above min_coverage
        for hit in &hits {
            let block_sizes_str = hit
                .block_sizes
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",");
            let block_starts_str = hit
                .block_starts
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",");
            lo_stmt
                .execute(params![
                    pk,
                    hit.query_name,
                    hit.query_start,
                    hit.query_end,
                    hit.strand.to_string(),
                    hit.coverage_pct,
                    hit.block_count as i64,
                    block_sizes_str,
                    block_starts_str,
                ])
                .expect("insert liftover");
        }
    }
    drop(lo_stmt);
    drop(sum_stmt);
    tx.commit().expect("commit liftover data");

    // --- Insert gene summary into DB ---
    {
        let tx = conn
            .unchecked_transaction()
            .expect("begin gene_summary transaction");
        let mut gs_stmt = tx
            .prepare(
                "INSERT INTO gene_summary
                     (gene_id, gene_name, ref_chrom,
                      total, single_full, single_partial,
                      multi_contig, multi_location, no_hit)
                 VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9)",
            )
            .expect("prepare gene_summary insert");
        for (gene_id, (gene_name, chrom, c)) in &gene_counts {
            let total = c.iter().sum::<u64>();
            gs_stmt
                .execute(params![
                    gene_id,
                    gene_name,
                    chrom,
                    total as i64,
                    c[0] as i64,
                    c[1] as i64,
                    c[2] as i64,
                    c[3] as i64,
                    c[4] as i64,
                ])
                .expect("insert gene_summary");
        }
        drop(gs_stmt);
        tx.commit().expect("commit gene_summary");
    }

    // --- Write anomaly summary file ---
    let summary_path = db_path.with_extension("anomaly.tsv");
    let mut sum_out = BufWriter::new(File::create(&summary_path)?);
    let total = transcripts.len() as u64;

    // Section 1: aggregate counts
    writeln!(sum_out, "## liftover anomaly summary")?;
    writeln!(sum_out, "## alndb      : {}", args.alndb_path)?;
    writeln!(sum_out, "## gtf        : {}", args.gtf_path)?;
    writeln!(sum_out, "## min_cov    : {:.0}%", args.min_coverage * 100.0)?;
    writeln!(
        sum_out,
        "## full_cov   : {:.0}%",
        args.full_coverage * 100.0
    )?;
    writeln!(sum_out, "#")?;
    writeln!(sum_out, "## status\tcount\tpct")?;
    for (label, n) in [
        ("total", total),
        ("single_full", cnt_single_full),
        ("single_partial", cnt_single_partial),
        ("multi_contig", cnt_multi_contig),
        ("multi_location", cnt_multi_location),
        ("no_hit", cnt_no_hit),
    ] {
        writeln!(sum_out, "##   {label:<20}\t{n}\t{:.2}%", pct(n, total))?;
    }
    writeln!(sum_out, "#")?;

    // Section 2: per-transcript anomaly list
    writeln!(
        sum_out,
        "transcript_id\tgene_id\tgene_name\tref_chrom\tref_start\tref_end\t\
         status\thit_count\tcontig_count\tbest_coverage\tbest_contig\tall_contigs"
    )?;
    for (
        tx_id,
        gene_id,
        gene_name,
        chrom,
        start,
        end,
        status,
        hit_count,
        contig_count,
        best_cov,
        best_contig,
        contigs,
    ) in &anomaly_rows
    {
        writeln!(
            sum_out,
            "{tx_id}\t{gene_id}\t{gene_name}\t{chrom}\t{start}\t{end}\t\
             {status}\t{hit_count}\t{contig_count}\t{best_cov:.1}\t{best_contig}\t{contigs}"
        )?;
    }

    // --- Gene-level summary — separate file ---
    let gene_summary_path = db_path.with_extension("gene.anomaly.tsv");
    let mut gene_out = BufWriter::new(File::create(&gene_summary_path)?);

    writeln!(gene_out, "## gene-level liftover summary")?;
    writeln!(gene_out, "## alndb      : {}", args.alndb_path)?;
    writeln!(gene_out, "## gtf        : {}", args.gtf_path)?;
    writeln!(
        gene_out,
        "## min_cov    : {:.0}%",
        args.min_coverage * 100.0
    )?;
    writeln!(
        gene_out,
        "## full_cov   : {:.0}%",
        args.full_coverage * 100.0
    )?;
    writeln!(gene_out, "#")?;
    writeln!(
        gene_out,
        "gene_id\tgene_name\tref_chrom\ttotal\tsingle_full\tsingle_partial\
         \tmulti_contig\tmulti_location\tno_hit\tpct_full"
    )?;

    // Sort by ref_chrom then by total anomalies descending
    let mut gene_vec: Vec<(&String, &(String, String, [u64; 5]))> = gene_counts.iter().collect();
    gene_vec.sort_by(|a, b| {
        let chrom_cmp = a.1 .1.cmp(&b.1 .1);
        if chrom_cmp != std::cmp::Ordering::Equal {
            return chrom_cmp;
        }
        // within same chrom: most anomalies (least single_full) first
        let a_anom = a.1 .2[1] + a.1 .2[2] + a.1 .2[3] + a.1 .2[4];
        let b_anom = b.1 .2[1] + b.1 .2[2] + b.1 .2[3] + b.1 .2[4];
        b_anom.cmp(&a_anom)
    });

    for (gene_id, (gene_name, chrom, c)) in &gene_vec {
        let total_g = c.iter().sum::<u64>();
        let pct_full = pct(c[0], total_g);
        writeln!(
            gene_out,
            "{gene_id}\t{gene_name}\t{chrom}\t{total_g}\t\
             {}\t{}\t{}\t{}\t{}\t{pct_full:.1}",
            c[0], c[1], c[2], c[3], c[4]
        )?;
    }

    eprintln!("Liftover anomaly summary  (total transcripts: {total})");
    eprintln!(
        "  single_full     {:>8}  ({:.1}%)",
        cnt_single_full,
        pct(cnt_single_full, total)
    );
    eprintln!(
        "  single_partial  {:>8}  ({:.1}%)",
        cnt_single_partial,
        pct(cnt_single_partial, total)
    );
    eprintln!(
        "  multi_contig    {:>8}  ({:.1}%)",
        cnt_multi_contig,
        pct(cnt_multi_contig, total)
    );
    eprintln!(
        "  multi_location  {:>8}  ({:.1}%)",
        cnt_multi_location,
        pct(cnt_multi_location, total)
    );
    eprintln!(
        "  no_hit          {:>8}  ({:.1}%)",
        cnt_no_hit,
        pct(cnt_no_hit, total)
    );
    eprintln!("Outputs:");
    eprintln!("  SQLite       : {}", args.output_db);
    eprintln!("  GTF          : {}", gtf_path.display());
    eprintln!(
        "  HQ GTF       : {} ({} genes, {:.0}% coverage threshold)",
        hq_gtf_path.display(),
        hq_genes.len(),
        args.full_coverage * 100.0
    );
    eprintln!("  Anomalies    : {}", summary_path.display());
    eprintln!("  Gene summary : {}", gene_summary_path.display());

    // --- Optional sequence pass (requires --ref-fa and --query-fa) ---
    if args.ref_fa.is_some() || args.query_fa.is_some() {
        let ref_seqs = if let Some(ref path) = args.ref_fa {
            eprintln!("Loading ref FASTA: {path}");
            load_fasta_seqs(path).unwrap_or_else(|e| panic!("cannot load ref FASTA '{path}': {e}"))
        } else {
            FxHashMap::default()
        };
        let query_seqs = if let Some(ref path) = args.query_fa {
            eprintln!("Loading query FASTA: {path}");
            load_fasta_seqs(path)
                .unwrap_or_else(|e| panic!("cannot load query FASTA '{path}': {e}"))
        } else {
            FxHashMap::default()
        };

        // Output FASTA paths
        let ref_fa_path = db_path.with_extension("ref_tx.fa");
        let ctg_fa_path = db_path.with_extension("contig_tx.fa");
        let hq_ctg_fa_path = db_path.with_extension("hq_contig_tx.fa");

        let mut ref_fa_out = BufWriter::new(File::create(&ref_fa_path)?);
        let mut ctg_fa_out = BufWriter::new(File::create(&ctg_fa_path)?);
        let mut hq_fa_out = BufWriter::new(File::create(&hq_ctg_fa_path)?);

        // Sequence insert statements
        let seq_tx = conn.unchecked_transaction().expect("begin seq transaction");
        let mut ref_seq_stmt = seq_tx
            .prepare("INSERT INTO ref_sequences (transcript_pk, ref_seq) VALUES (?1,?2)")
            .expect("prepare ref_seq insert");
        let mut ctg_seq_stmt = seq_tx
            .prepare("INSERT INTO contig_sequences (liftover_pk, contig_seq) VALUES (?1,?2)")
            .expect("prepare contig_seq insert");

        for t in &transcripts {
            let pk = tx_pk[&t.transcript_id];

            // Reference spliced sequence
            if !ref_seqs.is_empty() {
                let ref_seq = fetch_ref_tx_seq(&ref_seqs, &t.chrom, &t.exons, t.strand);
                let ref_header = format!(
                    "{} gene_id={} chrom={}:{}-{} strand={}",
                    t.transcript_id, t.gene_id, t.chrom, t.start, t.end, t.strand
                );
                write_fasta_record(&mut ref_fa_out, &ref_header, &ref_seq)?;
                ref_seq_stmt
                    .execute(params![pk, std::str::from_utf8(&ref_seq).unwrap_or("")])
                    .expect("insert ref_seq");
            }

            // Liftover contig sequences
            if !query_seqs.is_empty() {
                let hits = liftover_transcript(t, &target_intervals, args.min_coverage);
                for hit in &hits {
                    let ctg_seq = fetch_contig_tx_seq(&query_seqs, hit);
                    let ctg_header = format!(
                        "{} gene_id={} contig={}:{}-{} strand={} coverage={:.1}",
                        t.transcript_id,
                        t.gene_id,
                        hit.query_name,
                        hit.query_start,
                        hit.query_end,
                        hit.strand,
                        hit.coverage_pct
                    );
                    write_fasta_record(&mut ctg_fa_out, &ctg_header, &ctg_seq)?;
                    if hit.coverage_pct >= args.full_coverage * 100.0 {
                        write_fasta_record(&mut hq_fa_out, &ctg_header, &ctg_seq)?;
                    }

                    // Look up liftover_pk for the contig_sequences FK
                    let lo_pk: i64 = seq_tx
                        .query_row(
                            "SELECT liftover_pk FROM liftover
                             WHERE transcript_pk=?1 AND contig=?2
                               AND contig_start=?3 AND contig_end=?4",
                            params![pk, hit.query_name, hit.query_start, hit.query_end],
                            |r| r.get(0),
                        )
                        .expect("liftover_pk lookup");
                    ctg_seq_stmt
                        .execute(params![lo_pk, std::str::from_utf8(&ctg_seq).unwrap_or("")])
                        .expect("insert contig_seq");
                }
            }
        }
        drop(ref_seq_stmt);
        drop(ctg_seq_stmt);
        seq_tx.commit().expect("commit sequence data");

        eprintln!("  Ref seqs      : {}", ref_fa_path.display());
        eprintln!("  Contig seqs   : {}", ctg_fa_path.display());
        eprintln!("  HQ contig seqs: {}", hq_ctg_fa_path.display());
    }

    Ok(())
}

#[inline]
fn pct(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        n as f64 / total as f64 * 100.0
    }
}
