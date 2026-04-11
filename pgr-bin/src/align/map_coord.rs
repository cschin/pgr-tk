use clap::{self, Parser};
use iset::IntervalMap;
use pgr_db::aln::{aln_pair_map, wfa_align_bases};
// use rayon::prelude::*;
use pgr_db::ext::{get_fastx_reader, GZFastaReader};
use pgr_db::fasta_io::{reverse_complement, SeqRec};
use rusqlite::Connection;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

fn get_aln_blocks_from_db(db_path: &str) -> FxHashMap<u32, Vec<ShimmerMatchBlock>> {
    let conn = Connection::open_with_flags(db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .expect("can't open alndb file");
    let mut aln_blocks = FxHashMap::<u32, Vec<ShimmerMatchBlock>>::default();

    // M blocks (block_type = 0)
    {
        let mut stmt = conn
            .prepare(
                "SELECT aln_idx, target_name, target_start, target_end,
                        query_name, query_start, query_end, orientation
                 FROM blocks WHERE block_type = 0",
            )
            .expect("prepare blocks query");
        let mut rows = stmt.query([]).expect("query blocks");
        while let Some(row) = rows.next().expect("blocks row") {
            let aln_idx: u32 = row.get::<_, i64>(0).unwrap() as u32;
            let t_name: String = row.get(1).unwrap();
            let ts: u32 = row.get(2).unwrap();
            let te: u32 = row.get(3).unwrap();
            let q_name: String = row.get(4).unwrap();
            let qs: u32 = row.get(5).unwrap();
            let qe: u32 = row.get(6).unwrap();
            let orientation: u32 = row.get::<_, i64>(7).unwrap() as u32;
            let e = aln_blocks.entry(aln_idx).or_default();
            e.push((t_name, ts, te, q_name, qs, qe, orientation, "M".to_string()));
        }
    }

    // V records from variants table
    {
        let mut stmt = conn
            .prepare(
                "SELECT aln_idx, dup_flag, ovlp_flag, target_name, target_start, target_end,
                        query_name, query_start, query_end, orientation
                 FROM variants",
            )
            .expect("prepare variants query");
        let mut rows = stmt.query([]).expect("query variants");
        while let Some(row) = rows.next().expect("variants row") {
            let aln_idx: u32 = row.get::<_, i64>(0).unwrap() as u32;
            let dup_flag: i32 = row.get(1).unwrap();
            let ovlp_flag: i32 = row.get(2).unwrap();
            let t_name: String = row.get(3).unwrap();
            let ts: u32 = row.get(4).unwrap();
            let te: u32 = row.get(5).unwrap();
            let q_name: String = row.get(6).unwrap();
            let qs: u32 = row.get(7).unwrap();
            let qe: u32 = row.get(8).unwrap();
            let orientation: u32 = row.get::<_, i64>(9).unwrap() as u32;
            let rec_type = if dup_flag != 0 {
                "V_D".to_string()
            } else if ovlp_flag != 0 {
                "V_O".to_string()
            } else {
                "V".to_string()
            };
            let e = aln_blocks.entry(aln_idx).or_default();
            e.push((t_name, ts, te, q_name, qs, qe, orientation, rec_type));
        }
    }

    aln_blocks
}

/// given an alnmap or alndb file, two sequence files and a list of coordinates in the query sequence,
/// map those coordinates in the target sequence according to the alnmap
#[derive(Parser, Debug)]
#[clap(about, long_about = None)]
pub struct Args {
    /// path to the alnmap or alndb file
    #[clap(long, short)]
    pub alnmap_path: String,
    /// path to the target fasta file
    #[clap(long, short)]
    pub target_fasta_path: String,
    /// the path to the query fasta file
    #[clap(long, short)]
    pub query_fasta_path: String,
    /// path to query coordinate file
    #[clap(long, short)]
    pub coorindate_file_path: String,
    /// the prefix of the output files
    #[clap(long, short)]
    pub output_path: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    pub number_of_thread: usize,
}

type ShimmerMatchBlock = (String, u32, u32, String, u32, u32, u32, String);

pub fn run(args: Args) -> Result<(), std::io::Error> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let aln_blocks: FxHashMap<u32, Vec<ShimmerMatchBlock>> = if args.alnmap_path.ends_with(".alndb")
    {
        get_aln_blocks_from_db(&args.alnmap_path)
    } else {
        let f = BufReader::new(File::open(Path::new(&args.alnmap_path)).unwrap());
        let mut aln_blocks = FxHashMap::<u32, Vec<ShimmerMatchBlock>>::default();
        f.lines().for_each(|line| {
            if let Ok(line) = line {
                if line.trim().starts_with('#') {
                    return;
                };
                let fields = line.split('\t').collect::<Vec<&str>>();
                assert!(fields.len() > 3);
                let rec_type = fields[1];
                let err_msg = format!("fail to parse on {}", line);
                let aln_block_id = fields[0].parse::<u32>().expect(&err_msg);
                let t_name = fields[2];
                let ts = fields[3].parse::<u32>().expect(&err_msg);
                let te = fields[4].parse::<u32>().expect(&err_msg);
                let q_name = fields[5];
                let qs = fields[6].parse::<u32>().expect(&err_msg);
                let qe = fields[7].parse::<u32>().expect(&err_msg);
                let orientation = fields[8].parse::<u32>().expect(&err_msg);
                let e = aln_blocks.entry(aln_block_id).or_default();
                e.push((
                    t_name.to_string(),
                    ts,
                    te,
                    q_name.to_string(),
                    qs,
                    qe,
                    orientation,
                    rec_type.to_string(),
                ));
            }
        });
        aln_blocks
    };

    let blocks_to_intervals =
        |blocks: FxHashMap<u32, Vec<ShimmerMatchBlock>>| -> FxHashMap<String, IntervalMap<u32, ShimmerMatchBlock>> {
            let mut aln_intervals = FxHashMap::<String, IntervalMap<u32, ShimmerMatchBlock>>::default();
            blocks.into_iter().for_each(|(_block_id, records)| {
                records.into_iter().for_each(|rec| {
                    let t_name = rec.3.clone();
                    let bgn = rec.4;
                    let end = rec.5;
                    let interval_set = aln_intervals.entry(t_name).or_default();
                    interval_set.insert(bgn..end, rec);
                })
            });
            aln_intervals
        };

    let aln_intervals = blocks_to_intervals(aln_blocks);

    let mut target_seqs: Vec<SeqRec> = vec![];
    let mut add_target_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                target_seqs.push(r);
            };
        });
    };

    match get_fastx_reader(args.target_fasta_path, true)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => add_target_seqs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => add_target_seqs(&mut reader.into_iter()),
    };

    let target_seqs = target_seqs
        .into_iter()
        .map(|srec| (String::from_utf8_lossy(&srec.id[..]).to_string(), srec))
        .collect::<FxHashMap<String, SeqRec>>();

    let mut query_seqs: Vec<SeqRec> = vec![];
    let mut add_query_seqs = |seq_iter: &mut dyn Iterator<Item = io::Result<SeqRec>>| {
        seq_iter.into_iter().for_each(|r| {
            if let Ok(r) = r {
                query_seqs.push(r);
            };
        });
    };

    match get_fastx_reader(args.query_fasta_path, true)? {
        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::GZFile(reader) => add_query_seqs(&mut reader.into_iter()),

        #[allow(clippy::useless_conversion)] // the into_iter() is necessary for dyn patching
        GZFastaReader::RegularFile(reader) => add_query_seqs(&mut reader.into_iter()),
    };

    let query_seqs = query_seqs
        .into_iter()
        .map(|srec| (String::from_utf8_lossy(&srec.id[..]).to_string(), srec))
        .collect::<FxHashMap<String, SeqRec>>();

    // TODO: this is not an ideal cache mechanism, but it is easy to implement.
    // Find a better way in the future
    let mut cached_map = FxHashMap::<
        (String, u32, u32, String, u32, u32, u32),
        Option<Box<FxHashMap<u32, u32>>>,
    >::default();

    let mut get_target_position_map = |t_name: &String,
                                       ts: &u32,
                                       te: &u32,
                                       q_name: &String,
                                       qs: &u32,
                                       qe: &u32,
                                       orientation: &u32|
     -> Option<Box<FxHashMap<u32, u32>>> {
        let e = cached_map
            .entry((
                t_name.clone(),
                *ts,
                *ts,
                q_name.clone(),
                *qs,
                *qs,
                *orientation,
            ))
            .or_insert_with(|| {
                let t_seq = &target_seqs.get(t_name).unwrap().seq;
                let t_sub_seq = t_seq[(*ts as usize)..(*te as usize)].to_vec();

                let q_seq = &query_seqs.get(q_name).unwrap().seq;

                let q_sub_seq = if *orientation == 0 {
                    q_seq[(*qs as usize)..(*qe as usize)].to_vec()
                } else {
                    reverse_complement(&q_seq[(*qs as usize)..(*qe as usize)])
                };
                let t_str = String::from_utf8_lossy(&t_sub_seq[..]);
                let q_str = String::from_utf8_lossy(&q_sub_seq[..]);
                if let Some((aln_target_str, aln_query_str)) =
                    wfa_align_bases(&t_str, &q_str, 384, 4, 4, 1)
                {
                    let mut q_pos_to_t_pos_map = FxHashMap::<u32, u32>::default();
                    aln_pair_map(&aln_target_str, &aln_query_str)
                        .into_iter()
                        .for_each(|(tp, qp, _)| {
                            q_pos_to_t_pos_map.entry(qp).or_insert(tp);
                        });
                    Some(Box::new(q_pos_to_t_pos_map))
                } else {
                    None
                }
            });
        e.clone()
    };

    let mut position_of_interests = FxHashMap::<String, Vec<u32>>::default();

    let coorindate_file =
        BufReader::new(File::open(Path::new(&args.coorindate_file_path)).unwrap());

    coorindate_file.lines().for_each(|line| {
        if let Ok(line) = line {
            if line.trim().starts_with('#') {
                return;
            };
            let fields = line.split('\t').collect::<Vec<&str>>();

            let q_name = fields[0].to_string();
            let q_coordinate = fields[1].parse::<u32>().expect("parsing error");
            let e = position_of_interests.entry(q_name).or_default();
            e.push(q_coordinate);
        }
    });

    let mut out_file = BufWriter::new(File::create(Path::new(&args.output_path)).unwrap());

    position_of_interests
        .iter_mut()
        .for_each(|(q_name, q_coordiates)| {
            if let Some(interval_map) = aln_intervals.get(q_name) {
                q_coordiates.sort();
                q_coordiates.iter().for_each(|coordinate| {
                    let mut overlap_records = Vec::<(&String, &u32, &ShimmerMatchBlock)>::new();
                    interval_map.values_overlap(*coordinate).for_each(|block| {
                        overlap_records.push((q_name, coordinate, block));
                    });
                    if overlap_records.is_empty() {
                        writeln!(out_file, "{}\t{}\t*\t*\t*\t*\t0", q_name, coordinate)
                            .expect("can't write the output file");
                    } else {
                        let mut target_collection = FxHashSet::<(
                            String,
                            u32,
                            Option<String>,
                            Option<u32>,
                            u32,
                            String,
                        )>::default();
                        let mut unique_targets = FxHashSet::<(String, u32)>::default();
                        overlap_records
                            .into_iter()
                            .for_each(|(q_name, coordinate, block)| {
                                let (t_name, ts, te, _, qs, qe, orientation, btype) = block;
                                if btype.starts_with('M') {
                                    if *orientation == 0 {
                                        let t_name = t_name.clone();
                                        let t_coordinate = coordinate - qs + ts;
                                        target_collection.insert((
                                            q_name.clone(),
                                            *coordinate,
                                            Some(t_name.clone()),
                                            Some(t_coordinate),
                                            *orientation,
                                            btype.clone(),
                                        ));
                                        unique_targets.insert((t_name, t_coordinate));
                                    } else {
                                        let t_name = t_name.clone();
                                        let t_coordinate = (qe - coordinate) + ts; //TODO: check
                                        target_collection.insert((
                                            q_name.clone(),
                                            *coordinate,
                                            Some(t_name.clone()),
                                            Some(t_coordinate),
                                            *orientation,
                                            btype.clone(),
                                        ));
                                        unique_targets.insert((t_name, t_coordinate));
                                    }
                                } else if btype.starts_with('V') {
                                    let q_pos_to_t_pos_map = get_target_position_map(
                                        t_name,
                                        ts,
                                        te,
                                        q_name,
                                        qs,
                                        qe,
                                        orientation,
                                    );

                                    if let Some(q_pos_to_t_pos_map) = q_pos_to_t_pos_map {
                                        let q_pos = if *orientation == 0 {
                                            coordinate - qs
                                        } else {
                                            qe - coordinate
                                        };
                                        if let Some(t_pos) = q_pos_to_t_pos_map.get(&q_pos) {
                                            target_collection.insert((
                                                q_name.clone(),
                                                *coordinate,
                                                Some(t_name.clone()),
                                                Some(t_pos + ts),
                                                *orientation,
                                                btype.clone(),
                                            ));
                                            unique_targets.insert((t_name.clone(), t_pos + ts));
                                        } else {
                                            target_collection.insert((
                                                q_name.clone(),
                                                *coordinate,
                                                Some(t_name.clone()),
                                                None,
                                                *orientation,
                                                btype.clone(),
                                            ));
                                        }
                                    } else {
                                        target_collection.insert((
                                            q_name.clone(),
                                            *coordinate,
                                            Some(t_name.clone()),
                                            None,
                                            *orientation,
                                            btype.clone(),
                                        ));
                                    }
                                } else {
                                    target_collection.insert((
                                        q_name.clone(),
                                        *coordinate,
                                        Some(t_name.clone()),
                                        None,
                                        *orientation,
                                        btype.clone(),
                                    ));
                                };
                            });
                        let hit_count = unique_targets.len();
                        target_collection.into_iter().for_each(
                            |(q_name, q_pos, t_name, t_pos, orientation, btype)| {
                                let t_name = if let Some(t_name) = t_name {
                                    t_name
                                } else {
                                    "*".to_string()
                                };
                                let t_pos = if let Some(t_pos) = t_pos {
                                    format!("{}", t_pos)
                                } else {
                                    "*".to_string()
                                };
                                writeln!(
                                    out_file,
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    q_name, q_pos, t_name, t_pos, orientation, btype, hit_count
                                )
                                .expect("can't write the output file")
                            },
                        );
                    }
                });
            }
        });

    Ok(())
}
