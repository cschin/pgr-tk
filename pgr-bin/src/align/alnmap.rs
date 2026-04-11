use clap::{self, Parser};
use iset::set::IntervalSet;
use pgr_db::aln;
use pgr_db::ext::{get_fastx_reader, GZFastaReader, SeqIndexDB};
use pgr_db::fasta_io::{reverse_complement, SeqRec};
use rayon::prelude::*;
use rusqlite::{params, Connection};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

#[derive(Clone, Copy, clap::ValueEnum, Default, Debug)]
enum OptPreset {
    Fast,
    #[default]
    Default,
    Detail,
    Overwrite,
}

/// Align long contigs and identify potential SV regions with respect to the reference fasta file
#[derive(Parser, Debug)]
#[clap(about, long_about = None)]
pub struct Args {
    /// path to the reference fasta file
    #[clap(long, short = 'R')]
    pub reference_fasta_path: String,

    /// the path to the query assembly contig file
    #[clap(long, short)]
    pub assembly_contig_path: String,

    /// the prefix of the output files
    #[clap(long, short)]
    pub output_prefix: String,

    /// use preset parameters ( (w,k,r,min_span,max_sw_aln_size) = (80, 55, 4, 64, 1024) for fast, (48, 55, 2, 16, 32864) for detail)
    #[clap(long, default_value_t, value_enum)]
    pub preset: OptPreset,

    ///number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    pub number_of_thread: usize,

    /// overwrite the preset, minimizer window size: w 
    #[clap(long, short, default_value_t = 48)]
    pub w: u32,

    /// overwrite the preset, minimizer k-mer size
    #[clap(long, short, default_value_t = 55)]
    pub k: u32,

    /// overwrite the preset, sparse minimizer (shimmer) reduction factor
    #[clap(long, short, default_value_t = 2)]
    pub r: u32,

    /// overwrite the preset, min span for neighboring minimizers
    #[clap(long, short, default_value_t = 16)]
    pub min_span: u32,

    /// overwrite the preset, max size to do SW for calling structure variants
    #[clap(long, short, default_value_t = 1024)]
    pub max_sw_aln_size: u32,

    /// the gap penalty factor for sparse alignments in the SHIMMER space
    #[clap(long, default_value_t = 0.025)]
    pub gap_penalty_factor: f32,

    /// the max gap length allowed in the alignment blocks
    #[clap(long, default_value_t = 100000)]
    pub max_gap: u32,

    /// the span of the chain for building the sparse alignment directed acyclic graph
    #[clap(long, default_value_t = 8)]
    pub max_aln_chain_span: u32,

    /// if specified, generate fasta files for the sequence covering the SV candidates
    #[clap(long, short, default_value_t = false)]
    pub skip_uncalled_sv_seq_file: bool,
}

struct Parameters {
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    max_sw_aln_size: u32,
}

type ShimmerMatchBlock = (u32, u32, u32, u32, u32, u32, u32); //t_idx, ts, ted, q_idx, ts, te, orientation

#[derive(Clone)]
enum Record {
    Bgn(ShimmerMatchBlock, u32, u32), // MatchBlock, q_len, ctg_aln_orientation
    End(ShimmerMatchBlock, u32, u32), // MatchBlock, q_len, ctg_aln_orientation
    Match(ShimmerMatchBlock),
    SvCnd((ShimmerMatchBlock, AlnDiff, u32)), // MatchBlock, diff_type, ctg_aln_orientation
    Variant(ShimmerMatchBlock, u32, u32, u32, char, String, String),
}

// ((q_smp_start, q_smp_end, q_smp_orientation), (t_smp_start, t_smp_end, t_smp_orientation))
type AlignSegment = ((u32, u32, u8), (u32, u32, u8));

type AlignSegments = Vec<AlignSegment>;

type AlignmentResult = Vec<(u32, u32, char, String, String)>;

#[derive(Clone)]
enum AlnDiff {
    Aligned(AlignmentResult),
    FailAln,
    FailEndMatch,
    FailLengthDiff,
    FailShortSeq,
}

#[derive(Serialize)]
struct CtgMapRec {
    t_name: String,
    ts: u32,
    te: u32,
    q_name: String,
    qs: u32,
    qe: u32,
    ctg_len: u32,
    orientation: u32,
    ctg_orientation: u32,
    t_dup: bool,
    t_ovlp: bool,
    q_dup: bool,
    q_ovlp: bool,
}

#[derive(Serialize)]
struct CtgMapSet {
    records: Vec<CtgMapRec>,
    target_length: Vec<(u32, String, u32)>,
    query_length: Vec<(u32, String, u32)>,
}

fn init_alndb(conn: &Connection) {
    conn.execute_batch(
        "PRAGMA journal_mode=WAL;
         PRAGMA synchronous=NORMAL;
         CREATE TABLE IF NOT EXISTS run_params (
             key   TEXT PRIMARY KEY,
             value TEXT NOT NULL
         );
         CREATE TABLE IF NOT EXISTS sequences (
             seq_id   INTEGER PRIMARY KEY,
             seq_name TEXT NOT NULL,
             seq_type TEXT NOT NULL,
             length   INTEGER NOT NULL
         );
         CREATE TABLE IF NOT EXISTS chains (
             aln_idx         INTEGER PRIMARY KEY,
             target_name     TEXT NOT NULL,
             target_start    INTEGER NOT NULL,
             target_end      INTEGER NOT NULL,
             query_name      TEXT NOT NULL,
             query_start     INTEGER NOT NULL,
             query_end       INTEGER NOT NULL,
             orientation     INTEGER NOT NULL,
             ctg_orientation INTEGER NOT NULL,
             query_length    INTEGER NOT NULL,
             target_dup      INTEGER NOT NULL,
             target_ovlp     INTEGER NOT NULL,
             query_dup       INTEGER NOT NULL,
             query_ovlp      INTEGER NOT NULL
         );
         CREATE TABLE IF NOT EXISTS blocks (
             block_id        INTEGER PRIMARY KEY,
             aln_idx         INTEGER NOT NULL,
             block_type      INTEGER NOT NULL,
             dup_flag        INTEGER NOT NULL,
             ovlp_flag       INTEGER NOT NULL,
             target_name     TEXT NOT NULL,
             target_start    INTEGER NOT NULL,
             target_end      INTEGER NOT NULL,
             query_name      TEXT NOT NULL,
             query_start     INTEGER NOT NULL,
             query_end       INTEGER NOT NULL,
             orientation     INTEGER NOT NULL,
             ctg_orientation INTEGER,
             sv_diff_type    INTEGER
         );
         CREATE TABLE IF NOT EXISTS variants (
             variant_id   INTEGER PRIMARY KEY,
             aln_idx      INTEGER NOT NULL,
             dup_flag     INTEGER NOT NULL,
             ovlp_flag    INTEGER NOT NULL,
             target_name  TEXT NOT NULL,
             target_start INTEGER NOT NULL,
             target_end   INTEGER NOT NULL,
             query_name   TEXT NOT NULL,
             query_start  INTEGER NOT NULL,
             query_end    INTEGER NOT NULL,
             orientation  INTEGER NOT NULL,
             target_diff  INTEGER NOT NULL,
             query_diff   INTEGER NOT NULL,
             target_coord INTEGER NOT NULL,
             variant_type INTEGER NOT NULL,
             ref_seq      TEXT NOT NULL,
             alt_seq      TEXT NOT NULL
         );
         CREATE TABLE IF NOT EXISTS ctgmap (
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
         CREATE TABLE IF NOT EXISTS sv_candidates (
             svcnd_id     INTEGER PRIMARY KEY,
             target_name  TEXT NOT NULL,
             target_start INTEGER NOT NULL,
             target_end   INTEGER NOT NULL,
             sv_type      TEXT NOT NULL
         );
         CREATE TABLE IF NOT EXISTS ctgsv (
             ctgsv_id    INTEGER PRIMARY KEY,
             query_name  TEXT NOT NULL,
             query_start INTEGER NOT NULL,
             query_end   INTEGER NOT NULL,
             sv_type     TEXT NOT NULL
         );
         CREATE TABLE IF NOT EXISTS sv_sequences (
             svsq_id      INTEGER PRIMARY KEY,
             target_name  TEXT NOT NULL,
             target_start INTEGER NOT NULL,
             target_end   INTEGER NOT NULL,
             query_name   TEXT NOT NULL,
             query_start  INTEGER NOT NULL,
             query_end    INTEGER NOT NULL,
             orientation  INTEGER NOT NULL,
             target_seq   TEXT NOT NULL,
             query_seq    TEXT NOT NULL
         );",
    )
    .expect("failed to initialise alndb schema");
}

fn filter_aln(aln_segs: &AlignSegments) -> Vec<((u32, u32), (u32, u32))> {
    // the aln_segs should be sorted already
    let aln_segs = aln_segs.clone();

    let mut last_ts = aln_segs[0].1 .0;
    let mut last_te = aln_segs[0].1 .1;

    let mut last_qs = aln_segs[0].0 .0;
    let mut last_qe = aln_segs[0].0 .1;

    let mut rtn = Vec::<((u32, u32), (u32, u32))>::new();
    rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    for ((_qs, qe, qo), (ts, te, to)) in aln_segs {
        if te < ts {
            continue;
        };
        if qo != to {
            continue;
        };
        if ts > last_te {
            last_ts = last_te;
            last_te = te;

            last_qs = last_qe;
            last_qe = qe;
            if last_ts == last_te {
                continue;
            }
            rtn.push(((last_ts, last_te), (last_qs, last_qe)));
        }
    }
    rtn
}

fn filter_aln_rev(aln_segs: &AlignSegments) -> Vec<((u32, u32), (u32, u32))> {
    // the aln_segs should be sorted already
    let mut aln_segs = aln_segs.clone();
    aln_segs.reverse();

    let mut last_ts = aln_segs[0].1 .0;
    let mut last_te = aln_segs[0].1 .1;

    let mut last_qs = aln_segs[0].0 .0;
    let mut last_qe = aln_segs[0].0 .1;

    let mut rtn = Vec::<((u32, u32), (u32, u32))>::new();
    rtn.push(((last_ts, last_te), (last_qs, last_qe)));
    for ((qs, _qe, qo), (ts, te, to)) in aln_segs {
        if te < ts {
            continue;
        };
        if qo == to {
            continue;
        };
        if ts >= last_te {
            last_ts = last_te;
            last_te = te;

            last_qe = last_qs;
            last_qs = qs;
            if last_ts == last_te {
                continue;
            }
            rtn.push(((last_ts, last_te), (last_qs, last_qe)));
        }
    }
    rtn
}

pub fn run(args: Args) -> Result<(), std::io::Error> {

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let mut ref_seq_index_db = SeqIndexDB::new();

    let parameters = match args.preset {
        OptPreset::Fast => Parameters {
            w: 80,
            k: 55,
            r: 4,
            min_span: 64,
            max_sw_aln_size: 1 << 10,
        },
        OptPreset::Default => Parameters {
            w: 48,
            k: 55,
            r: 2,
            min_span: 16,
            max_sw_aln_size: 1 << 10,
        },
        OptPreset::Detail => Parameters {
            w: 48,
            k: 55,
            r: 2,
            min_span: 16,
            max_sw_aln_size: 1 << 15,
        },
        OptPreset::Overwrite => Parameters {
            w: args.w,
            k: args.k,
            r: args.r,
            min_span: args.min_span,
            max_sw_aln_size: args.max_sw_aln_size,
        },
    };

    ref_seq_index_db.load_from_fastx(
        args.reference_fasta_path,
        parameters.w,
        parameters.k,
        parameters.r,
        parameters.min_span,
        true,
    )?;

    let mut out_alnmap = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("alnmap")).unwrap(),
    );

    let mut out_vcf =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("vcf")).unwrap());

    let mut out_ctgmap = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("ctgmap.bed")).unwrap(),
    );

    let mut out_ctgmap_json = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("ctgmap.json")).unwrap(),
    );

    let mut out_target_len = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("target_len.json")).unwrap(),
    );

    let mut out_query_len = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("query_len.json")).unwrap(),
    );

    let mut out_svcnd = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("svcnd.bed")).unwrap(),
    );

    let mut out_ctgsv = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("ctgsv.bed")).unwrap(),
    );
    let mut out_sv_seq_file = if !args.skip_uncalled_sv_seq_file {
        Some(BufWriter::new(
            File::create(Path::new(&args.output_prefix).with_extension("svcnd.seqs")).unwrap(),
        ))
    } else {
        None
    };

    // --- SQLite output (alongside legacy file outputs) ---
    let db_path = Path::new(&args.output_prefix).with_extension("alndb");
    if db_path.exists() {
        std::fs::remove_file(&db_path).expect("failed to remove existing alndb");
    }
    let conn = Connection::open(&db_path).expect("failed to open alndb");
    init_alndb(&conn);

    let mut query_seqs: Vec<SeqRec> = vec![];
    let mut add_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                query_seqs.push(r);
            };
        });
    };
   

    match get_fastx_reader(args.assembly_contig_path, true)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => add_seqs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => add_seqs(&mut reader.into_iter()),
    };
    
    let kmer_size = parameters.k;

    let query_name = query_seqs
        .iter()
        .enumerate()
        .map(|(idx, seq_rec)| {
            (
                idx as u32,
                String::from_utf8_lossy(&seq_rec.id[..]).to_string(),
            )
        })
        .collect::<FxHashMap<_, _>>();

    let query_len = query_seqs
        .iter()
        .enumerate()
        .map(|(idx, seq_rec)| (idx as u32, seq_rec.seq.len()))
        .collect::<FxHashMap<_, _>>();

    let target_name = ref_seq_index_db
        .seq_info
        .as_ref()
        .unwrap()
        .iter()
        .map(|(k, v)| (*k, v.0.clone()))
        .collect::<FxHashMap<_, _>>();

    let target_len = ref_seq_index_db
        .seq_info
        .as_ref()
        .unwrap()
        .iter()
        .map(|(k, v)| (*k, v.2))
        .collect::<FxHashMap<_, _>>();

    // --- run_params and sequences ---
    {
        let tx = conn.unchecked_transaction().expect("begin transaction");
        let preset_str = format!("{:?}", args.preset);
        for (k, v) in &[
            ("w", parameters.w.to_string()),
            ("k", parameters.k.to_string()),
            ("r", parameters.r.to_string()),
            ("min_span", parameters.min_span.to_string()),
            ("max_sw_aln_size", parameters.max_sw_aln_size.to_string()),
            ("preset", preset_str),
        ] {
            tx.execute(
                "INSERT INTO run_params(key, value) VALUES(?1, ?2)",
                params![k, v],
            )
            .expect("insert run_params");
        }
        for (id, name) in &target_name {
            let len = *target_len.get(id).unwrap();
            tx.execute(
                "INSERT INTO sequences(seq_id, seq_name, seq_type, length) VALUES(?1,?2,'target',?3)",
                params![id, name, len],
            )
            .expect("insert target sequence");
        }
        for (id, name) in &query_name {
            let len = *query_len.get(id).unwrap() as u32;
            tx.execute(
                "INSERT INTO sequences(seq_id, seq_name, seq_type, length) VALUES(?1,?2,'query',?3)",
                params![*id + 1_000_000u32, name, len],
            )
            .expect("insert query sequence");
        }
        tx.commit().expect("commit sequences");
    }

    let all_records = query_seqs
        .par_iter()
        .enumerate()
        .map(|(q_idx, seq_rec)| {
            // let q_name = String::from_utf8_lossy(&seq_rec.id);
            let query_seq = seq_rec.seq.clone();
            //let q_len = query_seq.len();
            let max_gap = args.max_gap;
            let query_results = ref_seq_index_db.query_fragment_to_hps(
                &query_seq,
                args.gap_penalty_factor,
                Some(1),
                Some(1),
                Some(1),
                Some(args.max_aln_chain_span),
                Some(max_gap),
                true,
            );
            (q_idx, seq_rec, query_results)
        })
        .flat_map(|(q_idx, seq_rec, query_results)| {
            if let Some(qr) = query_results {
                let query_seq = &seq_rec.seq;
                let q_len: usize = query_seq.len();
                let mut target_id_to_mapped_regions = FxHashMap::default();
                let mut target_id_to_orientation_len_count = FxHashMap::default();
                qr.into_iter().for_each(|(t_idx, mapped_segments)| {
                    let mut aln_lens = vec![];
                    let mut ctg_orientation_count = (0_usize, 0_usize); // ctg level orientation count: (fwd_count, rev_count)
                    mapped_segments.into_iter().for_each(|(_score, aln)| {
                        let mut segment_orientation_count = (0_usize, 0_usize); // ctg level orientation count: (fwd_count, rev_count)
                        if aln.len() > 2 {
                            aln_lens.push(aln.len());
                            for hp in &aln {
                                let seg_len = (hp.0 .1 - hp.0 .0) as usize;
                                if hp.0 .2 == hp.1 .2 {
                                    ctg_orientation_count.0 += seg_len;
                                    segment_orientation_count.0 += seg_len;
                                } else {
                                    ctg_orientation_count.1 += seg_len;
                                    segment_orientation_count.1 += seg_len;
                                }
                            }
                            let seg_orientation =
                                if segment_orientation_count.0 > segment_orientation_count.1 {
                                    0_u32
                                } else {
                                    1_u32
                                };

                            let e = target_id_to_mapped_regions
                                .entry(t_idx)
                                .or_insert_with(Vec::new);
                            e.push((aln, seg_orientation))
                        }
                        let ctg_orientation = if ctg_orientation_count.0 > ctg_orientation_count.1 {
                            0_u32
                        } else {
                            1_u32
                        };
                        let e = target_id_to_orientation_len_count
                            .entry(t_idx)
                            .or_insert(((0, 0), 0));
                        *e = (ctg_orientation_count, ctg_orientation);
                    })
                });

                let rtn = target_id_to_mapped_regions
                    .into_iter()
                    .flat_map(|(t_idx, mapped_regions)| {
                        let mapped_region_aln = mapped_regions
                            .into_iter()
                            .map(|(aln_segs, orientation)| {
                                let aln_segs = if orientation == 0 {
                                    filter_aln(&aln_segs)
                                } else {
                                    filter_aln_rev(&aln_segs)
                                };

                                aln_segs
                                    .into_iter()
                                    .map(|((ts, te), (qs, qe))| {
                                        let ts = ts - kmer_size; // add one to ensure a match base if the first call is deletion
                                                                 //let te = te;
                                        let qs = if orientation == 0 { qs - kmer_size } else { qs };
                                        let qe = if orientation == 0 { qe } else { qe + kmer_size };
                                        // Fetch only the needed subregion instead of the full chromosome.
                                        // This avoids O(n_query_contigs) full chromosome reconstructions
                                        // and the associated memory allocation pressure.
                                        let s0str = ref_seq_index_db
                                            .get_sub_seq_by_id(t_idx, ts as usize, te as usize)
                                            .unwrap();
                                        let s1str = if orientation == 0 {
                                            query_seq[qs as usize..qe as usize].to_vec()
                                        } else {
                                            reverse_complement(
                                                &query_seq[(qs - kmer_size) as usize
                                                    ..(qe - kmer_size) as usize],
                                            )
                                        };

                                        let wf_aln_diff: AlnDiff = if s0str.len() <= 16
                                            || s1str.len() <= 16
                                        {
                                            AlnDiff::FailShortSeq
                                        } else if s0str[..16] != s1str[..16]
                                            || s0str[s0str.len() - 16..]
                                                != s1str[s1str.len() - 16..]
                                        {
                                            AlnDiff::FailEndMatch
                                        } else if (s0str.len() as isize - s1str.len() as isize)
                                            .abs()
                                            >= 128
                                        {
                                            // AlnDiff::FailLengthDiff
                                            if s0str.len() < parameters.max_sw_aln_size as usize
                                                && s1str.len() < parameters.max_sw_aln_size as usize
                                            {
                                                if let Some(aln_res) = aln::get_sw_variant_segments(
                                                    &s0str, &s1str, 1, 4, 4, 1,
                                                ) {
                                                    AlnDiff::Aligned(aln_res)
                                                } else {
                                                    AlnDiff::FailAln
                                                }
                                            } else {
                                                AlnDiff::FailLengthDiff
                                            }
                                        } else if let Some(aln_res) = aln::get_wfa_variant_segments(
                                            &s0str,
                                            &s1str,
                                            1,
                                            Some(384),
                                            4,
                                            4,
                                            1,
                                        ) {
                                            AlnDiff::Aligned(aln_res)
                                        } else {
                                            AlnDiff::FailAln
                                        };
                                        ((ts, te), (qs, qe), orientation, wf_aln_diff)
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .filter(|v| !v.is_empty())
                            .collect::<Vec<_>>();

                        let (_, ctg_orientation) =
                            target_id_to_orientation_len_count.get(&t_idx).unwrap();

                        mapped_region_aln
                            .into_iter()
                            .map(|v| {
                                let mut output_records = Vec::<Record>::new();
                                let ((ts, te), (qs, qe), orientation, _diff) = v[0].clone();
                                let qs = if orientation == 0 { qs } else { qs - kmer_size };
                                let qe = if orientation == 0 { qe } else { qe - kmer_size };
                                output_records.push(Record::Bgn(
                                    (t_idx, ts, te, q_idx as u32, qs, qe, orientation),
                                    q_len as u32,
                                    *ctg_orientation,
                                ));
                                let v_last = v.last().unwrap().clone();
                                v.into_iter().for_each(
                                    |((ts, te), (qs, qe), orientation, diff)| {
                                        let qs = if orientation == 0 { qs } else { qs - kmer_size };
                                        let qe = if orientation == 0 { qe } else { qe - kmer_size };
                                        if let AlnDiff::Aligned(diff) = diff {
                                            if diff.is_empty() {
                                                output_records.push(Record::Match((
                                                    t_idx,
                                                    ts,
                                                    te,
                                                    q_idx as u32,
                                                    qs,
                                                    qe,
                                                    orientation,
                                                )))
                                            } else {
                                                diff.into_iter().for_each(
                                                    |(td, qd, vt, t_str, q_str)| {
                                                        output_records.push(Record::Variant(
                                                            (
                                                                t_idx,
                                                                ts,
                                                                te,
                                                                q_idx as u32,
                                                                qs,
                                                                qe,
                                                                orientation,
                                                            ),
                                                            td,
                                                            qd,
                                                            ts + td,
                                                            vt,
                                                            t_str,
                                                            q_str,
                                                        ));
                                                    },
                                                )
                                            }
                                        } else {
                                            output_records.push(Record::SvCnd((
                                                (t_idx, ts, te, q_idx as u32, qs, qe, orientation),
                                                diff,
                                                *ctg_orientation,
                                            )));
                                        }
                                    },
                                );

                                let ((ts, te), (qs, qe), orientation, _diff) = v_last;
                                let qs = if orientation == 0 { qs } else { qs - kmer_size };
                                let qe = if orientation == 0 { qe } else { qe - kmer_size };
                                output_records.push(Record::End(
                                    (t_idx, ts, te, q_idx as u32, qs, qe, orientation),
                                    q_len as u32,
                                    *ctg_orientation,
                                ));
                                output_records
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                Some(rtn)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    let mut in_aln_sv_cnd_records = Vec::<(ShimmerMatchBlock, char, u32)>::new();
    let mut target_aln_blocks =
        FxHashMap::<u32, Vec<(usize, ShimmerMatchBlock, u32, u32)>>::default();
    let mut query_aln_blocks =
        FxHashMap::<u32, Vec<(usize, ShimmerMatchBlock, u32, u32)>>::default();

    // the first round loop through all_records for computing duplicated / overlapped match blocks
    all_records
        .iter()
        .flatten()
        .enumerate()
        .for_each(|(aln_idx, vr)| {
            let mut bgn_rec: Option<(ShimmerMatchBlock, u32, u32)> = None;
            let mut end_rec: Option<(ShimmerMatchBlock, u32, u32)> = None;
            vr.iter().for_each(|r| {
                match r.clone() {
                    Record::Bgn(match_block, q_len, ctg_orientation) => {
                        bgn_rec = Some((match_block, q_len, ctg_orientation));
                    }
                    Record::SvCnd((
                        (t_idx, ts, te, q_idx, qs, qe, orientation),
                        diff,
                        ctg_orientation,
                    )) => {
                        let diff_type = match diff {
                            AlnDiff::FailAln => 'A',
                            AlnDiff::FailEndMatch => 'E',
                            AlnDiff::FailShortSeq => 'S',
                            AlnDiff::FailLengthDiff => 'L',
                            _ => 'U',
                        };
                        in_aln_sv_cnd_records.push((
                            (t_idx, ts + 1, te + 1, q_idx, qs + 1, qe + 1, orientation),
                            diff_type,
                            ctg_orientation,
                        ));
                    }
                    Record::End(match_block, q_len, ctg_orientation) => {
                        end_rec = Some((match_block, q_len, ctg_orientation));
                    }
                    _ => {}
                };
                //writeln!(out_alnmap, "{}", rec_out).expect("fail to write the output file");
            });
            //aln_block.push( (aln_idx, bgn_rec.unwrap(), end_rec.unwrap()) );
            let (
                (b_t_idx, b_ts, _b_te, b_q_idx, b_qs, b_qe, b_orientation),
                _ctg_len,
                _ctg_orientation,
            ) = bgn_rec.unwrap();
            let (
                (e_t_idx, _e_ts, e_te, e_q_idx, e_qs, e_qe, e_orientation),
                ctg_len,
                ctg_orientation,
            ) = end_rec.unwrap();
            assert_eq!(b_orientation, e_orientation);
            assert_eq!(b_t_idx, e_t_idx);
            assert_eq!(b_q_idx, e_q_idx);
            let t_entry = target_aln_blocks.entry(b_t_idx).or_insert_with(Vec::new);
            let q_entry = query_aln_blocks.entry(b_q_idx).or_insert_with(Vec::new);
            if b_orientation == 0 {
                t_entry.push((
                    aln_idx,
                    (b_t_idx, b_ts, e_te, b_q_idx, b_qs, e_qe, b_orientation),
                    ctg_len,
                    ctg_orientation,
                ));
                q_entry.push((
                    aln_idx,
                    (b_t_idx, b_ts, e_te, b_q_idx, b_qs, e_qe, b_orientation),
                    ctg_len,
                    ctg_orientation,
                ));
            } else {
                t_entry.push((
                    aln_idx,
                    (b_t_idx, b_ts, e_te, b_q_idx, e_qs, b_qe, b_orientation),
                    ctg_len,
                    ctg_orientation,
                ));
                q_entry.push((
                    aln_idx,
                    (b_t_idx, b_ts, e_te, b_q_idx, e_qs, b_qe, b_orientation),
                    ctg_len,
                    ctg_orientation,
                ));
            }
        });

    let mut target_aln_blocks = target_aln_blocks.into_iter().collect::<Vec<_>>();
    target_aln_blocks.sort();

    let mut target_aln_bed_records = Vec::<(String, u32, u32, String)>::new();
    let mut target_duplicate_blocks = FxHashSet::<ShimmerMatchBlock>::default();
    let mut target_overlap_blocks = FxHashSet::<ShimmerMatchBlock>::default();
    target_aln_blocks
        .iter_mut()
        .for_each(|(t_idx, match_blocks)| {
            match_blocks.sort_by_key(|v| v.1 .1);
            let mut cts = 0_u32;
            let mut cte = 0_u32;
            let mut c_ctg = &String::from("BGN");
            let t_name = target_name.get(t_idx).unwrap();
            match_blocks
                .iter()
                .for_each(|&(_aln_idx, match_block, ctg_len, ctg_orientation)| {
                    let (t_idx, ts, te, q_idx, qs, qe, orientation) = match_block;
                    //println!("T {} {} {} {} {} {} {}", t_name, ts, te, q_idx, qs, qe, orientation);
                    let next_ctg = query_name.get(&q_idx).unwrap();
                    if ts > cte {
                        let bed_annotation = format!(
                            "TG:{}>{}:{}:{}:{}:{}:{}",
                            c_ctg, next_ctg, qs, qe, ctg_len, orientation, ctg_orientation
                        );
                        target_aln_bed_records.push((t_name.clone(), cte, ts, bed_annotation));
                        //println!("G {} {} {} {} {}", t_name, cte, ts, c_ctg, next_ctg);
                        c_ctg = next_ctg;
                        cts = ts;
                        cte = te;
                    } else if te <= cte {
                        let bed_annotation = format!(
                            "TD:{}>{}:{}:{}:{}:{}:{}",
                            c_ctg, next_ctg, qs, qe, ctg_len, orientation, ctg_orientation
                        );
                        target_duplicate_blocks.insert(match_block);
                        target_aln_bed_records.push((t_name.clone(), ts, te, bed_annotation));
                        //println!("D {} {} {} {} {}", t_name, cts, te, c_ctg, next_ctg);
                    } else {
                        let bed_annotation = format!(
                            "TO:{}>{}:{}:{}:{}:{}:{}",
                            c_ctg, next_ctg, qs, qe, ctg_len, orientation, ctg_orientation
                        );
                        target_overlap_blocks.insert((t_idx, ts, cte, q_idx, qs, qe, orientation));
                        target_aln_bed_records.push((t_name.clone(), ts, cte, bed_annotation));
                        //println!("O {} {} {} {} {}", t_name, ts, cte, c_ctg, next_ctg);
                        c_ctg = next_ctg;
                        cte = te;
                    };
                });
            let next_ctg = &String::from("END");
            let t_len = *target_len.get(t_idx).unwrap();
            let bed_annotation = format!("TG:{}>{}", c_ctg, next_ctg);
            target_aln_bed_records.push((t_name.clone(), cte, t_len, bed_annotation));
        });

    let mut query_aln_bed_records = Vec::<(String, u32, u32, String)>::new();
    let mut query_duplicate_blocks = FxHashSet::<ShimmerMatchBlock>::default();
    let mut query_overlap_blocks = FxHashSet::<ShimmerMatchBlock>::default();
    query_aln_blocks
        .iter_mut()
        .for_each(|(q_idx, match_blocks)| {
            match_blocks.sort_by_key(|v| v.1 .4);
            let mut cqs = 0_u32;
            let mut cqe = 0_u32;
            let mut c_target = &String::from("BGN");
            let q_name = query_name.get(q_idx).unwrap();
            match_blocks
                .iter()
                .for_each(|&(_aln_idx, match_block, ctg_len, ctg_orientation)| {
                    //println!("Q {} {} {} {} {} {} {}", t_name, ts, te, q_idx, qs, qe, orientation);
                    let (t_idx, ts, te, _q_idx, qs, qe, orientation) = match_block;
                    let next_target = target_name.get(&t_idx).unwrap();
                    if qs > cqe {
                        let bed_annotation = format!(
                            "QG:{}>{}:{}:{}:{}:{}:{}",
                            c_target, next_target, ts, te, ctg_len, orientation, ctg_orientation
                        );
                        query_aln_bed_records.push((q_name.clone(), cqe, qs, bed_annotation));
                        //println!("G {} {} {} {}", t_name, ts, te, p_target);
                        c_target = next_target;
                        cqs = qs;
                        cqe = qe;
                    } else if qe <= cqe {
                        let bed_annotation = format!(
                            "QD:{}>{}:{}:{}:{}:{}:{}",
                            c_target, next_target, ts, te, ctg_len, orientation, ctg_orientation
                        );
                        query_duplicate_blocks.insert(match_block);
                        query_aln_bed_records.push((q_name.clone(), qs, qe, bed_annotation));
                        //println!("D {} {} {} {}", t_name, ts, te, p_target);
                    } else {
                        let bed_annotation = format!(
                            "QO:{}>{}:{}:{}:{}:{}:{}",
                            c_target, next_target, ts, te, ctg_len, orientation, ctg_orientation
                        );
                        query_overlap_blocks.insert(match_block);
                        query_aln_bed_records.push((q_name.clone(), qs, cqe, bed_annotation));
                        //println!("O {} {} {} {}", t_name, ts, te, p_target);
                        c_target = next_target;
                        cqe = qe;
                    }
                });
            let next_target = &String::from("END");
            let bed_annotation = format!("QG:{}>{}", c_target, next_target);
            let q_len = *query_len.get(q_idx).unwrap() as u32;
            query_aln_bed_records.push((q_name.clone(), cqe, q_len, bed_annotation));
        });

    let mut target_duplicate_intervals = FxHashMap::<u32, IntervalSet<u32>>::default();
    target_duplicate_blocks
        .iter()
        .for_each(|block: &ShimmerMatchBlock| {
            let e = target_duplicate_intervals.entry(block.0).or_default();
            if block.2 > block.1 {
                e.insert(block.1..block.2);
            }
        });

    let mut target_overlap_intervals = FxHashMap::<u32, IntervalSet<u32>>::default();
    target_overlap_blocks
        .iter()
        .for_each(|block: &ShimmerMatchBlock| {
            let e = target_overlap_intervals.entry(block.0).or_default();
            if block.2 > block.1 {
                e.insert(block.1..block.2);
            }
        });

    let mut in_aln_sv_and_bed_records = Vec::<(String, u32, u32, String)>::new();
    in_aln_sv_cnd_records.sort();
    in_aln_sv_cnd_records.iter().for_each(
        |((t_idx, ts, te, q_idx, qs, qe, orientation), diff_type, ctg_orientation)| {
            let q_name = query_name.get(q_idx).unwrap();
            let dup =
                if let Some(target_duplicate_intervals) = target_duplicate_intervals.get(t_idx) {
                    if te > ts {
                        target_duplicate_intervals.has_overlap(ts..te)
                    } else {
                        false
                    }
                } else {
                    false
                };

            let ovlp = if let Some(target_overlap_intervals) = target_overlap_intervals.get(t_idx) {
                if te > ts {
                    target_overlap_intervals.has_overlap(ts..te)
                } else {
                    false
                }
            } else {
                false
            };
            let svc_type = if dup {
                "SVC_D"
            } else if ovlp {
                "SVC_O"
            } else {
                "SVC"
            };

            let bed_annotation = format!(
                "{}:{}:{}-{}:{}:{}:{}",
                svc_type, q_name, qs, qe, orientation, ctg_orientation, diff_type
            );
            let t_name = target_name.get(t_idx).unwrap();
            in_aln_sv_and_bed_records.push((t_name.clone(), ts + 1, te + 1, bed_annotation));
        },
    );

    let mut all_bed_records = Vec::<_>::new();
    all_bed_records.extend(in_aln_sv_and_bed_records);
    all_bed_records.extend(target_aln_bed_records);
    //all_bed_record.extend(query_aln_bed_records);
    all_bed_records.sort();

    {
        let tx = conn.unchecked_transaction().expect("begin svcnd tx");
        let mut stmt = tx
            .prepare(
                "INSERT INTO sv_candidates(target_name,target_start,target_end,sv_type)
                 VALUES(?1,?2,?3,?4)",
            )
            .expect("prepare sv_candidates insert");
        all_bed_records.iter().for_each(|r| {
            writeln!(out_svcnd, "{}\t{}\t{}\t{}", r.0, r.1, r.2, r.3)
                .expect("fail to write the 'in-alignment' sv candidate bed file");
            stmt.execute(params![r.0, r.1, r.2, r.3])
                .expect("insert sv_candidates");
        });
        drop(stmt);
        tx.commit().expect("commit sv_candidates");
    }

    // output ctgmap file

    let mut ctgmap_records = Vec::<CtgMapRec>::new();
    target_aln_blocks
        .into_iter()
        .for_each(|(t_idx, match_blocks)| {
            let t_name = target_name.get(&t_idx).unwrap();
            match_blocks
                .iter()
                .for_each(|&(_aln_idx, match_block, ctg_len, ctg_orientation)| {
                    let (_t_idx, ts, te, q_idx, qs, qe, orientation) = match_block;
                    let q_name = query_name.get(&q_idx).unwrap();
                    let t_dup = if target_duplicate_blocks.contains(&match_block) {
                        1
                    } else {
                        0
                    };
                    let t_ovlp = if target_overlap_blocks.contains(&match_block) {
                        1
                    } else {
                        0
                    };
                    let q_dup = if query_duplicate_blocks.contains(&match_block) {
                        1
                    } else {
                        0
                    };
                    let q_ovlp = if query_overlap_blocks.contains(&match_block) {
                        1
                    } else {
                        0
                    };
                    ctgmap_records.push(CtgMapRec {
                        t_name: t_name.clone(),
                        ts,
                        te,
                        q_name: q_name.clone(),
                        qs,
                        qe,
                        ctg_len,
                        orientation,
                        ctg_orientation,
                        t_dup: t_dup == 1,
                        t_ovlp: t_ovlp == 1,
                        q_dup: q_dup == 1,
                        q_ovlp: q_ovlp == 1,
                    });
                    writeln!(
                        out_ctgmap,
                        "{}\t{}\t{}\t{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
                        t_name,
                        ts,
                        te,
                        q_name,
                        qs,
                        qe,
                        ctg_len,
                        orientation,
                        ctg_orientation,
                        t_dup,
                        t_ovlp,
                        q_dup,
                        q_ovlp
                    )
                    .expect("can't write ctgmap file");
                });
        });

    let query_length = query_len
        .into_iter()
        .map(|(id, length)| (id, query_name.get(&id).unwrap().clone(), length as u32))
        .collect::<Vec<_>>();
    let target_length = target_len
        .into_iter()
        .map(|(id, length)| (id, target_name.get(&id).unwrap().clone(), length))
        .collect::<Vec<_>>();
    let ctg_map_set = CtgMapSet {
        records: ctgmap_records,
        query_length,
        target_length,
    };

    let ctgmap_json =
        serde_json::to_string(&ctg_map_set).expect("fail to construct json for ctg map");
    writeln!(out_ctgmap_json, "{}", ctgmap_json).expect("fail to write ctg map json file");

    // --- ctgmap inserts ---
    {
        let tx = conn.unchecked_transaction().expect("begin ctgmap tx");
        let mut stmt = tx
            .prepare(
                "INSERT INTO ctgmap(target_name,target_start,target_end,
                  query_name,query_start,query_end,ctg_len,orientation,ctg_orientation,
                  target_dup,target_ovlp,query_dup,query_ovlp)
                 VALUES(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13)",
            )
            .expect("prepare ctgmap insert");
        for r in &ctg_map_set.records {
            stmt.execute(params![
                r.t_name, r.ts, r.te, r.q_name, r.qs, r.qe, r.ctg_len,
                r.orientation, r.ctg_orientation,
                r.t_dup as i32, r.t_ovlp as i32, r.q_dup as i32, r.q_ovlp as i32
            ])
            .expect("insert ctgmap");
        }
        drop(stmt);
        tx.commit().expect("commit ctgmap");
    }

    let target_length_json = serde_json::to_string(&ctg_map_set.target_length)
        .expect("fail to construct json for ctg map");
    writeln!(out_target_len, "{}", target_length_json).expect("fail to write ctg map json file");

    let query_length_json = serde_json::to_string(&ctg_map_set.query_length)
        .expect("fail to construct json for ctg map");
    writeln!(out_query_len, "{}", query_length_json).expect("fail to write ctg map json file");

    query_aln_bed_records.sort();
    {
        let tx = conn.unchecked_transaction().expect("begin ctgsv tx");
        let mut stmt = tx
            .prepare(
                "INSERT INTO ctgsv(query_name,query_start,query_end,sv_type)
                 VALUES(?1,?2,?3,?4)",
            )
            .expect("prepare ctgsv insert");
        query_aln_bed_records.iter().for_each(|r| {
            writeln!(out_ctgsv, "{}\t{}\t{}\t{}", r.0, r.1, r.2, r.3)
                .expect("fail to write the 'in-alignment' sv candidate bed file");
            stmt.execute(params![r.0, r.1, r.2, r.3])
                .expect("insert ctgsv");
        });
        drop(stmt);
        tx.commit().expect("commit ctgsv");
    }

    let mut vcf_records = Vec::<(u32, u32, String, String, ShimmerMatchBlock)>::new();

    // the second round loop through all_records to output and tagged variant from duplicate / overlapped blocks
    let tx = conn.unchecked_transaction().expect("begin main tx");
    let mut chain_stmt = tx
        .prepare(
            "INSERT INTO chains(aln_idx,target_name,target_start,target_end,
              query_name,query_start,query_end,orientation,ctg_orientation,query_length,
              target_dup,target_ovlp,query_dup,query_ovlp)
             VALUES(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14)",
        )
        .expect("prepare chains insert");
    let mut block_stmt = tx
        .prepare(
            "INSERT INTO blocks(aln_idx,block_type,dup_flag,ovlp_flag,
              target_name,target_start,target_end,query_name,query_start,query_end,
              orientation,ctg_orientation,sv_diff_type)
             VALUES(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13)",
        )
        .expect("prepare blocks insert");
    let mut variant_stmt = tx
        .prepare(
            "INSERT INTO variants(aln_idx,dup_flag,ovlp_flag,
              target_name,target_start,target_end,query_name,query_start,query_end,
              orientation,target_diff,query_diff,target_coord,variant_type,ref_seq,alt_seq)
             VALUES(?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16)",
        )
        .expect("prepare variants insert");
    let mut svsq_stmt = tx
        .prepare(
            "INSERT INTO sv_sequences(target_name,target_start,target_end,
              query_name,query_start,query_end,orientation,target_seq,query_seq)
             VALUES(?1,?2,?3,?4,?5,?6,?7,?8,?9)",
        )
        .expect("prepare sv_sequences insert");

    for (aln_idx, vr) in all_records.into_iter().flatten().enumerate() {
        for r in vr.into_iter() {
            let rec_out = match r.clone() {
                    Record::Bgn(match_block, q_len, ctg_orientation) => {
                        let (t_idx, ts, te, q_idx, qs, qe, orientation) = match_block;
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();
                        let t_dup = if target_duplicate_blocks.contains(&match_block) { 1 } else { 0 };
                        let t_ovlp = if target_overlap_blocks.contains(&match_block) { 1 } else { 0 };
                        let q_dup = if query_duplicate_blocks.contains(&match_block) { 1 } else { 0 };
                        let q_ovlp = if query_overlap_blocks.contains(&match_block) { 1 } else { 0 };
                        chain_stmt.execute(params![
                            aln_idx as u32, tn, ts, te, qn, qs, qe,
                            orientation, ctg_orientation, q_len,
                            t_dup, t_ovlp, q_dup, q_ovlp
                        ]).expect("insert chain");
                        format!(
                            "{:06}\tB\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            aln_idx, tn, ts, te, qn, qs, qe, orientation,
                            q_len, ctg_orientation, t_dup, t_ovlp, q_dup, q_ovlp
                        )
                    }
                    Record::End(match_block, q_len, ctg_orientation) => {
                        let (t_idx, ts, te, q_idx, qs, qe, orientation) = match_block;
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();
                        format!(
                            "{:06}\tE\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            aln_idx, tn, ts, te, qn, qs, qe, orientation, q_len, ctg_orientation
                        )
                    }
                    Record::Match((t_idx, ts, te, q_idx, qs, qe, orientation)) => {
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();
                        let dup = if let Some(target_duplicate_intervals) =
                            target_duplicate_intervals.get(&t_idx)
                        {
                            if te > ts {
                                target_duplicate_intervals.has_overlap(ts..te)
                            } else {
                                false
                            }
                        } else {
                            false
                        };

                        let ovlp = if let Some(target_overlap_intervals) =
                            target_overlap_intervals.get(&t_idx)
                        {
                            if te > ts {
                                target_overlap_intervals.has_overlap(ts..te)
                            } else {
                                false
                            }
                        } else {
                            false
                        };
                        let match_type = if dup {
                            "M_D"
                        } else if ovlp {
                            "M_O"
                        } else {
                            "M"
                        };

                        block_stmt.execute(params![
                            aln_idx as u32, 0i32, dup as i32, ovlp as i32,
                            tn, ts, te, qn, qs, qe, orientation,
                            rusqlite::types::Null, rusqlite::types::Null
                        ]).expect("insert M block");
                        format!(
                            "{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            aln_idx, match_type, tn, ts, te, qn, qs, qe, orientation
                        )
                    }
                    Record::SvCnd((
                        (t_idx, ts, te, q_idx, qs, qe, orientation),
                        diff,
                        ctg_orientation,
                    )) => {
                        let diff_type = match diff {
                            AlnDiff::FailAln => 'A',
                            AlnDiff::FailEndMatch => 'E',
                            AlnDiff::FailShortSeq => 'S',
                            AlnDiff::FailLengthDiff => 'L',
                            _ => 'U',
                        };

                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();
                        let dup = if let Some(target_duplicate_intervals) =
                            target_duplicate_intervals.get(&t_idx)
                        {
                            target_duplicate_intervals.has_overlap(ts..te)
                        } else {
                            false
                        };

                        let ovlp = if let Some(target_overlap_intervals) =
                            target_overlap_intervals.get(&t_idx)
                        {
                            if te > ts {
                                target_overlap_intervals.has_overlap(ts..te)
                            } else {
                                false
                            }
                        } else {
                            false
                        };

                        let svc_type = if dup {
                            "S_D"
                        } else if ovlp {
                            "S_O"
                        } else {
                            "S"
                        };

                        let out = format!(
                            "{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            aln_idx,
                            svc_type,
                            tn,
                            ts,
                            te,
                            qn,
                            qs,
                            qe,
                            orientation,
                            ctg_orientation,
                            diff_type
                        );

                        let sv_diff_code = match diff_type {
                            'A' => 0i32, 'E' => 1, 'S' => 2, 'L' => 3, _ => 4,
                        };
                        block_stmt.execute(params![
                            aln_idx as u32, 1i32, dup as i32, ovlp as i32,
                            tn, ts, te, qn, qs, qe, orientation,
                            ctg_orientation, sv_diff_code
                        ]).expect("insert S block");

                        if let Some(out_sv_seq_file) = out_sv_seq_file.as_mut() {
                            let t_seq_slice = &ref_seq_index_db
                                .get_sub_seq_by_id(t_idx, ts as usize, te as usize)
                                .unwrap()[..];
                            let t_seq = String::from_utf8_lossy(t_seq_slice);
                            let q_seq = if orientation == 0 {
                                query_seqs[q_idx as usize].seq[(qs as usize)..(qe as usize)]
                                    .to_vec()
                            } else {
                                reverse_complement(
                                    &query_seqs[q_idx as usize].seq[(qs as usize)..(qe as usize)],
                                )
                            };
                            let q_seq = String::from_utf8_lossy(&q_seq[..]);
                            writeln!(out_sv_seq_file, "{}\t{}\t{}", out, t_seq, q_seq)
                                .expect("writing fasta for SV candidate fail");
                            svsq_stmt.execute(params![
                                tn, ts, te, qn, qs, qe, orientation,
                                t_seq.as_ref(), q_seq.as_ref()
                            ]).expect("insert sv_sequences");
                        };

                        out
                    }
                    Record::Variant(match_block, td, qd, tc, vt, tvs, qvs) => {
                        let (t_idx, ts, te, q_idx, qs, qe, orientation) = match_block;
                        vcf_records.push((t_idx, tc + 1, tvs.clone(), qvs.clone(), match_block));
                        let tn = target_name.get(&t_idx).unwrap();
                        let qn = query_name.get(&q_idx).unwrap();

                        let dup = if let Some(target_duplicate_intervals) =
                            target_duplicate_intervals.get(&t_idx)
                        {
                            target_duplicate_intervals.has_overlap(ts..te)
                        } else {
                            false
                        };

                        let ovlp = if let Some(target_overlap_intervals) =
                            target_overlap_intervals.get(&t_idx)
                        {
                            if te > ts {
                                target_overlap_intervals.has_overlap(ts..te)
                            } else {
                                false
                            }
                        } else {
                            false
                        };

                        let variant_type = if dup { "V_D" } else if ovlp { "V_O" } else { "V" };
                        let vt_code: i32 = match vt { 'X' => 0, 'I' => 1, _ => 2 };
                        variant_stmt.execute(params![
                            aln_idx as u32, dup as i32, ovlp as i32,
                            tn, ts, te, qn, qs, qe, orientation,
                            td, qd, tc, vt_code, tvs, qvs
                        ]).expect("insert variant");
                        format!(
                            "{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            aln_idx, variant_type, tn, ts, te, qn, qs, qe,
                            orientation, td, qd, tc, vt, tvs, qvs
                        )
                    }
                };
                writeln!(out_alnmap, "{}", rec_out).expect("fail to write the output file");
        }
    }
    drop(chain_stmt);
    drop(block_stmt);
    drop(variant_stmt);
    drop(svsq_stmt);
    tx.commit().expect("commit main alignment records");

    writeln!(out_vcf, "##fileformat=VCFv4.2").expect("fail to write the vcf file");
    ctg_map_set
        .target_length
        .into_iter()
        .for_each(|(_, t_name, t_len)| {
            writeln!(out_vcf, r#"##contig=<ID={},length={}>"#, t_name, t_len)
                .expect("fail to write the vcf file");
        });
    writeln!(
        out_vcf,
        r#"##FILTER=<ID=td,Description="variant from duplicated contig alignment on target">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FILTER=<ID=to,Description="variant from overlapped contig alignment on query">"#
    )
    .expect("fail to write the vcf file");
    writeln!(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        .expect("fail to write the vcf file");

    vcf_records.sort();
    vcf_records
        .into_iter()
        .for_each(|(t_idx, tc, tvs, qvs, match_block)| {
            let tn = target_name.get(&t_idx).unwrap();

            let dup =
                if let Some(target_duplicate_intervals) = target_duplicate_intervals.get(&t_idx) {
                    if match_block.2 > match_block.1 {
                        target_duplicate_intervals.has_overlap(match_block.1..match_block.2)
                    } else {
                        false
                    }
                } else {
                    false
                };

            let ovlp = if let Some(target_overlap_intervals) = target_overlap_intervals.get(&t_idx)
            {
                if match_block.2 > match_block.1 {
                    target_overlap_intervals.has_overlap(match_block.1..match_block.2)
                } else {
                    false
                }
            } else {
                false
            };
            let filter = if dup {
                "DUP"
            } else if ovlp {
                "OVLP"
            } else {
                "PASS"
            };
            let qv: u32 = if filter != "PASS" { 10 } else { 60 };
            writeln!(
                out_vcf,
                "{}\t{}\t.\t{}\t{}\t{}\t{}\t.",
                tn,
                tc,
                tvs.trim_end_matches('-'),
                qvs.trim_end_matches('-'),
                qv,
                filter
            )
            .expect("fail to write the vcf file");
        });

    conn.execute_batch(
        "CREATE INDEX IF NOT EXISTS blocks_pos
             ON blocks(target_name, target_start, target_end);
         CREATE INDEX IF NOT EXISTS variants_pos
             ON variants(target_name, target_coord);
         CREATE INDEX IF NOT EXISTS ctgmap_tgt
             ON ctgmap(target_name, target_start, target_end);
         CREATE INDEX IF NOT EXISTS ctgmap_qry
             ON ctgmap(query_name, query_start, query_end);
         CREATE INDEX IF NOT EXISTS svcnd_pos
             ON sv_candidates(target_name, target_start, target_end);",
    )
    .expect("failed to create alndb indices");

    Ok(())
}
