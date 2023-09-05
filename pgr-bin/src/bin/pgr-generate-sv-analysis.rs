const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
// use rayon::prelude::*;
use pgr_db::aln;
use pgr_db::ext::{get_principal_bundle_decomposition, SeqIndexDB};
use pgr_db::fasta_io::reverse_complement;
use rustc_hash::FxHashMap;
use serde::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// perform structural variation principle bundle decomposition
#[derive(Parser, Debug)]
#[clap(name = "pgr-sv-decompse")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the data file containing the reference and query sequence of SV candidates
    sv_candidate_seq_path: String,
    /// the prefix of the output files
    output_prefix: String,
    /// the prefix of the output files
    #[clap(long, default_value = "Sample")]
    sample_name: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

#[derive(Debug)]
struct CandidateRecord {
    aln_id: u32,
    svc_type: String,
    target_name: String,
    ts: u32,
    te: u32,
    query_name: String,
    qs: u32,
    qe: u32,
    orientation: u8,
    ctg_orientation: u8,
    aln_type: String,
    target_sequence: Vec<u8>,
    query_sequence: Vec<u8>,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug, Serialize, Deserialize)]
struct BundleSegment {
    bgn: u32,
    end: u32,
    bundle_id: u32,
    bundle_v_count: u32,
    bundle_dir: u32,
    bundle_v_bgn: u32,
    bundle_v_end: u32,
    is_repeat: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
enum AlnType {
    Match,
    Insertion,
    Deletion,
    Begin,
}
type AlnPathElement = (usize, usize, AlnType, u32, u32, usize, usize);
type AlnPath = Vec<AlnPathElement>;

struct ShimmerConfig {
    w: u32,
    k: u32,
    r: u32,
    min_span: u32,
    min_cov: u32,
    min_branch_size: u32,
}

type AlignmentResult = Vec<(u32, u32, char, String, String)>;

#[derive(Clone, Debug)]
enum AlnDiff {
    Aligned(AlignmentResult),
    FailAln,
    _FailEndMatch,
    FailLengthDiff,
    FailShortSeq,
}

type ShimmerMatchBlock = (String, u32, u32, String, u32, u32, u32); //t_name, ts, ted, q_name, ts, te, orientation
#[derive(Clone, Debug)]
enum Record {
    _Bgn(ShimmerMatchBlock, u32, u32), // MatchBlock, q_len, ctg_aln_orientation
    _End(ShimmerMatchBlock, u32, u32), // MatchBlock, q_len, ctg_aln_orientation
    Match(ShimmerMatchBlock, String, String), // MatchBlock, target_bundle_path, query_bundle_path
    SvCnd((ShimmerMatchBlock, AlnDiff, u32, String, String)), // MatchBlock, diff_type, ctg_aln_orientation, target_bundle_path, query_bundle_path
    Variant(
        ShimmerMatchBlock,
        u32,
        u32,
        u32,
        char,
        String,
        String,
        String,
        String,
    ), // MatchBlock, Variant Info..., target_bundle_path, query_bundle_path
}

fn align_bundles(
    q_bundles: &Vec<BundleSegment>,
    t_bundles: &Vec<BundleSegment>,
) -> (f32, usize, usize, AlnPath) {
    let q_count = q_bundles.len();
    let t_count = t_bundles.len();
    let mut s_map = FxHashMap::<(usize, usize), i64>::default();
    let mut t_map = FxHashMap::<(usize, usize), AlnType>::default();

    let mut get_aln_direction_with_best_score =
        |q_idx: usize, t_idx: usize, s_map: &FxHashMap<(usize, usize), i64>| -> (AlnType, i64) {
            let mut best = (AlnType::Match, i64::MIN);
            let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
            let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
            let min_len = if q_len > t_len { t_len } else { q_len };
            let q_b_seg = q_bundles[q_idx];
            let t_b_seg = t_bundles[t_idx];
            if q_idx == 0
                && t_idx == 0
                && (q_b_seg.bundle_id == t_b_seg.bundle_id)
                && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
            {
                best = (AlnType::Match, 2 * min_len)
            } else if q_idx == 0 && t_idx == 0 {
                best = (AlnType::Begin, 0)
            };

            if q_idx > 0
                && t_idx > 0
                && q_b_seg.bundle_id == t_b_seg.bundle_id
                && (q_b_seg.bundle_dir == t_b_seg.bundle_dir)
            {
                best = (
                    AlnType::Match,
                    2 * min_len + s_map.get(&(q_idx - 1, t_idx - 1)).unwrap(),
                )
            };
            if t_idx > 0 {
                let score = -2 * q_len + s_map.get(&(q_idx, t_idx - 1)).unwrap();
                if score > best.1 {
                    best = (AlnType::Deletion, score)
                };
            };
            if q_idx > 0 {
                let score = -2 * t_len + s_map.get(&(q_idx - 1, t_idx)).unwrap();
                if score > best.1 {
                    best = (AlnType::Insertion, score)
                }
            }
            t_map.insert((q_idx, t_idx), best.0);
            best
        };

    let mut aln_path = AlnPath::new();

    (0..t_count)
        .flat_map(|t_idx| (0..q_count).map(move |q_idx| (q_idx, t_idx)))
        .for_each(|(q_idx, t_idx)| {
            //println!("{} {}", q_idx, t_idx);
            let (_, score) = get_aln_direction_with_best_score(q_idx, t_idx, &s_map);
            s_map.insert((q_idx, t_idx), score);
        });
    let mut q_idx = q_count - 1;
    let mut t_idx = t_count - 1;
    let mut diff_len = 0_usize;
    let mut max_len = 1_usize;
    while let Some(aln_type) = t_map.get(&(q_idx, t_idx)) {
        let qq_idx = q_idx;
        let tt_idx = t_idx;
        let (diff_len_delta, max_len_delta) = match aln_type {
            AlnType::Match => {
                let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                let diff_len_delta = (q_len - t_len).unsigned_abs() as usize;
                let max_len_delata = if q_len > t_len {
                    q_len as usize
                } else {
                    t_len as usize
                };
                q_idx -= 1;
                t_idx -= 1;
                (diff_len_delta, max_len_delata)
            }
            AlnType::Insertion => {
                let q_len = (q_bundles[q_idx].end as i64 - q_bundles[q_idx].bgn as i64).abs();
                q_idx -= 1;
                (q_len as usize, q_len as usize)
            }
            AlnType::Deletion => {
                let t_len = (t_bundles[t_idx].end as i64 - t_bundles[t_idx].bgn as i64).abs();
                t_idx -= 1;
                (t_len as usize, t_len as usize)
            }
            AlnType::Begin => break,
        };
        diff_len += diff_len_delta;
        max_len += max_len_delta;
        aln_path.push((
            qq_idx,
            tt_idx,
            *aln_type,
            q_bundles[qq_idx].bundle_id,
            t_bundles[tt_idx].bundle_id,
            diff_len_delta,
            max_len_delta,
        ));
    }
    aln_path.reverse();
    (
        diff_len as f32 / max_len as f32,
        diff_len,
        max_len,
        aln_path,
    )
}

#[allow(clippy::type_complexity)]
fn group_smps_by_principle_bundle_id(
    smps: &[((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>)],
    bundle_length_cutoff: usize,
    bundle_merge_distance: usize,
) -> Vec<Vec<((u64, u64, u32, u32, u8), usize, u32, usize)>> {
    let mut pre_bundle_id: Option<usize> = None;
    let mut pre_direction: Option<u32> = None;
    let mut all_partitions = vec![];
    let mut new_partition = vec![];
    smps.iter().for_each(|&(smp, bundle_info)| {
        if bundle_info.is_none() {
            return;
        };
        let bundle_info = bundle_info.unwrap();
        let d = if smp.4 == bundle_info.1 { 0_u32 } else { 1_u32 };
        let bid = bundle_info.0;
        let bpos = bundle_info.2;
        if pre_bundle_id.is_none() {
            new_partition.clear();
            new_partition.push((smp, bid, d, bpos));
            pre_bundle_id = Some(bid);
            pre_direction = Some(d);
            return;
        };
        if bid != pre_bundle_id.unwrap() || d != pre_direction.unwrap() {
            let l = new_partition.len();
            if new_partition[l - 1].0 .3 as usize - new_partition[0].0 .2 as usize
                > bundle_length_cutoff
            {
                all_partitions.push(new_partition.clone());
                new_partition.clear();
            } else {
                new_partition.clear();
            };
            pre_bundle_id = Some(bid);
            pre_direction = Some(d);
        };
        new_partition.push((smp, bid, d, bpos));
    });
    let l = new_partition.len();
    if l > 0
        && new_partition[l - 1].0 .3 as usize - new_partition[0].0 .2 as usize
            > bundle_length_cutoff
    {
        all_partitions.push(new_partition);
    };

    let mut rtn_partitions = vec![];

    if all_partitions.is_empty() {
        return rtn_partitions;
    }
    let mut partition = all_partitions[0].clone();
    (1..all_partitions.len()).for_each(|idx| {
        let p = all_partitions[idx].clone();
        let p_len = partition.len();
        let p_end = partition[p_len - 1].0 .3;
        let p_bid = partition[p_len - 1].1;
        let p_d = partition[p_len - 1].2;
        let np_bgn = p[0].0 .2;
        let np_bid = p[0].1;
        let np_d = p[0].2;
        if p_bid == np_bid
            && p_d == np_d
            && (np_bgn as i64 - p_end as i64).abs() < bundle_merge_distance as i64
        {
            partition.extend(p);
        } else {
            rtn_partitions.push(partition.clone());
            partition = p;
        }
    });
    if !partition.is_empty() {
        rtn_partitions.push(partition);
    }
    rtn_partitions
}

fn get_aln_diff(s0str: &[u8], s1str: &[u8]) -> AlnDiff {
    let wf_aln_diff: AlnDiff = if s0str.is_empty() || s1str.is_empty() {
        AlnDiff::FailShortSeq
    } else if (s0str.len() as isize - s1str.len() as isize).abs() >= 128 {
        AlnDiff::FailLengthDiff
    } else if let Some(aln_res) = aln::get_variant_segments(s0str, s1str, 1, Some(384), 3, 3, 1) {
        AlnDiff::Aligned(aln_res)
    } else {
        AlnDiff::FailAln
    };
    wf_aln_diff
}

fn aln_segments(
    ts: usize,
    te: usize,
    qs: usize,
    qe: usize,
    rec: &CandidateRecord,
    target_bundle_path: &str,
    query_bundle_path: &str,
) -> Vec<Record> {
    let target_name = &rec.target_name;
    let query_name = &rec.query_name;
    let target_seg_sequence = &rec.target_sequence[ts..te];
    let query_seg_sequence = &rec.query_sequence[qs..qe];
    let diff = get_aln_diff(target_seg_sequence, query_seg_sequence);

    let ts = ts as u32 + rec.ts;
    let te = te as u32 + rec.ts;

    let (qs, qe) = if rec.orientation == 0 {
        (qs as u32 + rec.qs, qe as u32 + rec.qs)
    } else {
        (rec.qe - qe as u32, rec.qe - qs as u32)
    };
    let mut aln_block_records = Vec::<Record>::new();
    if let AlnDiff::Aligned(diff) = diff {
        if diff.is_empty() {
            aln_block_records.push(Record::Match(
                (
                    target_name.clone(),
                    ts,
                    te,
                    query_name.clone(),
                    qs,
                    qe,
                    rec.orientation as u32,
                ),
                target_bundle_path.to_string(),
                query_bundle_path.to_string(),
            ))
        } else {
            diff.into_iter().for_each(|(td, qd, vt, t_str, q_str)| {
                aln_block_records.push(Record::Variant(
                    (
                        target_name.clone(),
                        ts,
                        te,
                        query_name.clone(),
                        qs,
                        qe,
                        rec.orientation as u32,
                    ),
                    td,
                    qd,
                    ts + td,
                    vt,
                    t_str,
                    q_str,
                    target_bundle_path.to_string(),
                    query_bundle_path.to_string(),
                ));
            })
        }
    } else {
        aln_block_records.push(Record::SvCnd((
            (
                target_name.clone(),
                ts,
                te,
                query_name.clone(),
                qs,
                qe,
                rec.orientation as u32,
            ),
            diff,
            rec.ctg_orientation as u32,
            target_bundle_path.to_string(),
            query_bundle_path.to_string(),
        )));
    }
    aln_block_records
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let seq_pair_file = BufReader::new(File::open(Path::new(&args.sv_candidate_seq_path)).unwrap());
    let mut paired_seq_records = Vec::new();

    seq_pair_file.lines().for_each(|line| {
        if let Ok(line) = line {
            let fields = line.split('\t').collect::<Vec<&str>>();
            let paser_err_msg = "can't parse the input file";
            assert!(fields.len() == 13);
            let aln_id = fields[0].parse::<u32>().expect(paser_err_msg);
            let svc_type = fields[1].to_string();
            let target_name = fields[2].to_string();
            let ts = fields[3].parse::<u32>().expect(paser_err_msg);
            let te = fields[4].parse::<u32>().expect(paser_err_msg);
            let query_name = fields[5].to_string();
            let qs = fields[6].parse::<u32>().expect(paser_err_msg);
            let qe = fields[7].parse::<u32>().expect(paser_err_msg);
            let orientation = fields[8].parse::<u8>().expect(paser_err_msg);
            let ctg_orientation = fields[9].parse::<u8>().expect(paser_err_msg);
            let aln_type = fields[10].to_string();
            let target_sequence = fields[11].into();
            let query_sequence = if orientation == 0 {
                fields[12].into()
            } else {
                reverse_complement(fields[12].as_bytes())
            };

            let rec = CandidateRecord {
                aln_id,
                svc_type,
                target_name,
                ts,
                te,
                query_name,
                qs,
                qe,
                orientation,
                ctg_orientation,
                aln_type,
                target_sequence,
                query_sequence,
            };
            paired_seq_records.push(rec);
        }
    });

    let mut outpu_alnmap_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("svcnd.alnmap"))
            .expect("can't create the output file"),
    );

    paired_seq_records.into_iter().for_each(|rec| {
        let mut sdb = SeqIndexDB::new();
        let seq_list = vec![
            (rec.target_name.clone(), rec.target_sequence.clone()),
            (rec.query_name.clone(), rec.query_sequence.clone()),
        ];

        let shmmr_cfg = ShimmerConfig {
            w: 31,
            k: 23,
            r: 1,
            min_span: 13,
            min_cov: 0,
            min_branch_size: 0,
        };
        sdb.load_from_seq_list(
            seq_list,
            None,
            shmmr_cfg.w,
            shmmr_cfg.k,
            shmmr_cfg.r,
            shmmr_cfg.min_span,
        )
        .expect("can't load the sequences");
        let (principal_bundles_with_id, vertex_to_bundle_id_direction_pos) = sdb
            .get_principal_bundles_with_id(
                shmmr_cfg.min_cov as usize,
                shmmr_cfg.min_branch_size as usize,
                Some(vec![0, 1]),
            );
        let sid_smps = get_principal_bundle_decomposition(&vertex_to_bundle_id_direction_pos, &sdb);
        let sid_smps: FxHashMap<u32, Vec<_>> = sid_smps.into_iter().collect();
        let bid_to_size = principal_bundles_with_id
            .iter()
            .map(|v| (v.0, v.2.len()))
            .collect::<FxHashMap<usize, usize>>();
        let seq_info = sdb
            .seq_info
            .unwrap()
            .into_iter()
            .map(|(k, v)| (k, v))
            .collect::<Vec<_>>();

        let bundle_length_cutoff = 0;
        let bundle_merge_distance = 0;

        let mut repeat_count = FxHashMap::<u32, Vec<u32>>::default();
        let mut non_repeat_count = FxHashMap::<u32, Vec<u32>>::default();
        let mut sid_to_bundle_segs = FxHashMap::<u32, Vec<BundleSegment>>::default();

        seq_info.iter().for_each(|(sid, sdata)| {
            let (_ctg, _src, _len) = sdata;
            let smps = sid_smps.get(sid).unwrap();
            let smp_partitions = group_smps_by_principle_bundle_id(
                smps,
                bundle_length_cutoff,
                bundle_merge_distance,
            );
            let mut ctg_bundle_count = FxHashMap::<usize, usize>::default();
            smp_partitions.iter().for_each(|p| {
                let bid = p[0].1;
                *ctg_bundle_count.entry(bid).or_insert_with(|| 0) += 1;
            });
            let mut bundle_segs = Vec::<BundleSegment>::new();
            smp_partitions.into_iter().for_each(|p| {
                let b = p[0].0 .2 - shmmr_cfg.k;
                let e = p[p.len() - 1].0 .3;
                let bid = p[0].1;
                let direction = p[0].2;
                let is_repeat = if *ctg_bundle_count.get(&bid).unwrap_or(&0) > 1 {
                    repeat_count
                        .entry(*sid)
                        .or_insert_with(Vec::new)
                        .push(e - b - shmmr_cfg.k);
                    true
                } else {
                    non_repeat_count
                        .entry(*sid)
                        .or_insert_with(Vec::new)
                        .push(e - b - shmmr_cfg.k);
                    false
                };
                bundle_segs.push(BundleSegment {
                    bgn: b,
                    end: e,
                    bundle_id: bid as u32,
                    bundle_v_count: bid_to_size[&bid] as u32,
                    bundle_dir: direction,
                    bundle_v_bgn: p[0].3 as u32,
                    bundle_v_end: p[p.len() - 1].3 as u32,
                    is_repeat,
                });
            });

            sid_to_bundle_segs.insert(*sid, bundle_segs);
        });

        let target_bundles = sid_to_bundle_segs.get(&0).unwrap();
        let query_bundles = sid_to_bundle_segs.get(&1).unwrap();
        let (_dist0, _diff_len0, _max_len0, aln_path) =
            align_bundles(query_bundles, target_bundles);

        writeln!(
            outpu_alnmap_file,
            "## {:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            rec.aln_id,
            rec.svc_type,
            rec.target_name,
            rec.ts,
            rec.te,
            rec.query_name,
            rec.qs,
            rec.qe,
            rec.orientation,
            rec.ctg_orientation,
            rec.aln_type
        )
        .expect("can't write the alnmap output file");

        let mut cur_target_seg_bgn = 0u32;
        let mut cur_query_seg_bgn = 0u32;
        let mut aln_block_records = Vec::<Record>::new();
        let mut pre_aln_type = None;
        let mut pre_target_bundles = Vec::<(u32, u32, u8)>::new();
        let mut pre_query_bundles = Vec::<(u32, u32, u8)>::new();
        aln_path.into_iter().for_each(|elm| {
            let (qb_idx, tb_idx, aln_type, _qb_bid, _tb_bid, _diff_delta, _max_diff) = elm;
            let t_seg = target_bundles[tb_idx];
            let q_seg = query_bundles[qb_idx];

            // writeln!(outpu_alnmap_file, "{:?}", elm);

            match aln_type {
                AlnType::Match => {
                    match pre_aln_type {
                        Some(AlnType::Match) => (), // if the previous is a match, there is no sequence between two match blocks
                        _ => {
                            // process any thing after the previous match block before the current match block
                            let ts = cur_target_seg_bgn as usize;
                            let te = (t_seg.bgn + shmmr_cfg.k) as usize;
                            let qs = cur_query_seg_bgn as usize;
                            let qe = (q_seg.bgn + shmmr_cfg.k) as usize;
                            let target_path = pre_target_bundles
                                .iter()
                                .map(|(id, dir, rep)| format!("{}:{}:{}", id, dir, rep))
                                .collect::<Vec<_>>()
                                .join("-");
                            let query_path = pre_query_bundles
                                .iter()
                                .map(|(id, dir, rep)| format!("{}:{}:{}", id, dir, rep))
                                .collect::<Vec<_>>()
                                .join("-");
                            let target_path = if target_path.is_empty() { "*" } else { &target_path[..] }; 
                            let query_path = if query_path.is_empty() { "*" } else { &query_path[..] }; 

                            if ts != te || qs != qe {
                                // println!("target_segment: {} {} {}", ts, te, te - ts);
                                // println!("query_segment: {} {} {}", qs, qe, qe - qs);
                                aln_block_records.extend(aln_segments(
                                    ts,
                                    te,
                                    qs,
                                    qe,
                                    &rec,
                                    target_path,
                                    query_path,
                                ))
                            };
                        }
                    }

                    pre_target_bundles.clear();
                    pre_query_bundles.clear();

                    // process the current match block
                    let ts = t_seg.bgn as usize;
                    let te = t_seg.end as usize;
                    let qs = q_seg.bgn as usize;
                    let qe = q_seg.end as usize;
                    let target_bundle_info = format!(
                        "{}:{}:{}",
                        t_seg.bundle_id,
                        t_seg.bundle_dir,
                        if t_seg.is_repeat { 1 } else { 0 }
                    );
                    let query_bundle_info = format!(
                        "{}:{}:{}",
                        q_seg.bundle_id,
                        q_seg.bundle_dir,
                        if q_seg.is_repeat { 1 } else { 0 }
                    );

                    aln_block_records.extend(aln_segments(
                        ts,
                        te,
                        qs,
                        qe,
                        &rec,
                        &target_bundle_info[..],
                        &query_bundle_info[..],
                    ));

                    // println!("target_m_segment: {} {} {}", ts, te, te - ts);
                    // println!("query_m_segment: {} {} {}", qs, qe, qe - qs);

                    cur_target_seg_bgn = t_seg.end - shmmr_cfg.k;
                    cur_query_seg_bgn = q_seg.end - shmmr_cfg.k;
                }
                AlnType::Deletion => {
                    pre_target_bundles.push((
                        t_seg.bundle_id,
                        t_seg.bundle_dir,
                        if t_seg.is_repeat { 1 } else { 0 },
                    ));
                }
                AlnType::Insertion => {
                    pre_query_bundles.push((
                        q_seg.bundle_id,
                        q_seg.bundle_dir,
                        if q_seg.is_repeat { 1 } else { 0 },
                    ));
                }
                _ => (),
            };

            pre_aln_type = Some(aln_type);
        });
        // the last segment
        let ts = cur_target_seg_bgn as usize;
        let te = rec.target_sequence.len();
        let qs = cur_query_seg_bgn as usize;
        let qe = rec.query_sequence.len();
        if ts != te && qs != qe {
            //println!("target_e_segment: {} {} {}", ts, te, te - ts);
            //println!("query_e_segment: {} {} {}", qs, qe, qe - qs);
            aln_block_records.extend(aln_segments(ts, te, qs, qe, &rec, "*", "*"))
        };

        // generate the alnmap records
        aln_block_records.iter().for_each(|record| {
            let rec_out = match record.clone() {
                Record::Match((tn, ts, te, qn, qs, qe, orientation), target_path, query_path) => {
                    let match_type = if rec.svc_type.ends_with('D') {
                        "M_D"
                    } else if rec.svc_type.ends_with('O') {
                        "M_O"
                    } else {
                        "M"
                    };

                    Some(format!(
                        "{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        rec.aln_id,
                        match_type,
                        tn,
                        ts,
                        te,
                        qn,
                        qs,
                        qe,
                        orientation,
                        target_path,
                        query_path
                    ))
                }
                Record::SvCnd((
                    (tn, ts, te, qn, qs, qe, orientation),
                    diff,
                    ctg_orientation,
                    target_path,
                    query_path,
                )) => {
                    let diff_type = match diff {
                        AlnDiff::FailAln => 'A',
                        AlnDiff::_FailEndMatch => 'E',
                        AlnDiff::FailShortSeq => 'S',
                        AlnDiff::FailLengthDiff => 'L',
                        _ => 'U',
                    };

                    let svc_type = if rec.svc_type.ends_with('D') {
                        "S_D"
                    } else if rec.svc_type.ends_with('O') {
                        "S_O"
                    } else {
                        "S"
                    };

                    Some(format!(
                        "{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        rec.aln_id,
                        svc_type,
                        tn,
                        ts,
                        te,
                        qn,
                        qs,
                        qe,
                        orientation,
                        ctg_orientation,
                        diff_type,
                        target_path,
                        query_path
                    ))
                }
                Record::Variant(match_block, td, qd, tc, vt, tvs, qvs, target_path, query_path) => {
                    let (tn, ts, te, qn, qs, qe, orientation) = match_block;
                    let variant_type = if rec.svc_type.ends_with('D') {
                        "V_D"
                    } else if rec.svc_type.ends_with('O') {
                        "V_O"
                    } else {
                        "V"
                    };
                    Some(format!(
                        "{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        rec.aln_id,
                        variant_type,
                        tn,
                        ts,
                        te,
                        qn,
                        qs,
                        qe,
                        orientation,
                        td,
                        qd,
                        tc,
                        vt,
                        tvs,
                        qvs,
                        target_path,
                        query_path
                    ))
                }
                _ => None,
            };
            if let Some(rec_out) = rec_out {
                writeln!(outpu_alnmap_file, "{}", rec_out)
                    .expect("can't write the alnmap output file");
            }
        });
    });

    Ok(())
}
