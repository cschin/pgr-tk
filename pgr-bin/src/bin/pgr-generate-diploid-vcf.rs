const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use iset::set::IntervalSet;
// use rayon::prelude::*;
use rusqlite::Connection;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Generate diploid VCF field from paired alnmap file from two haplotype assembly
#[derive(Parser, Debug)]
#[clap(name = "pgr-generate-diploid-vcf")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the first haplotype alnmap or alndb file
    hap0_path: String,
    /// path to the second haplotype alnmap or alndb file
    hap1_path: String,
    /// the prefix of the output files
    output_prefix: String,
    /// path to target_len.json (required when using .alnmap input; inferred from .alndb)
    #[clap(long)]
    target_len_json: Option<String>,
    /// the sample name in the VCF
    #[clap(long, default_value = "Sample")]
    sample_name: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

type TargetSeqLength = Vec<(u32, String, u32)>;

type ShimmerMatchBlock = (String, u32, u32, String, u32, u32, u32);
type VariantRecord = (String, u32, u32, u64, u8, String, String, String); //t_name, tc, tl, aln_block_id, hap_type, tvs, qvs, rec_type

fn get_variant_recs_from_db(
    db_path: &str,
    hap_type: u8,
) -> (
    Vec<VariantRecord>,
    FxHashMap<u64, Vec<ShimmerMatchBlock>>,
    FxHashMap<u64, Vec<ShimmerMatchBlock>>,
) {
    let conn = Connection::open_with_flags(db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .expect("can't open alndb file");

    let mut variant_records = Vec::<VariantRecord>::new();
    let mut aln_blocks = FxHashMap::<u64, Vec<ShimmerMatchBlock>>::default();
    let mut unique_aln_blocks = FxHashMap::<u64, Vec<ShimmerMatchBlock>>::default();

    // Variant records from the variants table
    {
        let mut stmt = conn
            .prepare(
                "SELECT aln_idx, dup_flag, ovlp_flag, target_name, target_coord,
                        ref_seq, alt_seq
                 FROM variants ORDER BY target_name, target_coord",
            )
            .expect("prepare variants query");
        let mut rows = stmt.query([]).expect("query variants");
        while let Some(row) = rows.next().expect("variants row") {
            let aln_idx: u64 = row.get::<_, i64>(0).unwrap() as u64;
            let dup_flag: i32 = row.get(1).unwrap();
            let ovlp_flag: i32 = row.get(2).unwrap();
            let t_name: String = row.get(3).unwrap();
            let tc: u32 = row.get(4).unwrap();
            let tvs: String = row.get(5).unwrap();
            let qvs: String = row.get(6).unwrap();
            let rec_type = if dup_flag != 0 {
                "V_D".to_string()
            } else if ovlp_flag != 0 {
                "V_O".to_string()
            } else {
                "V".to_string()
            };
            let tl = tvs.len() as u32;
            variant_records.push((t_name, tc, tl, aln_idx, hap_type, tvs, qvs, rec_type));
        }
    }

    // Block-level coords from both M blocks and V records contribute to aln_blocks
    let add_blocks = |stmt_sql: &str,
                      aln_blocks: &mut FxHashMap<u64, Vec<ShimmerMatchBlock>>,
                      unique_aln_blocks: &mut FxHashMap<u64, Vec<ShimmerMatchBlock>>| {
        let mut stmt = conn.prepare(stmt_sql).expect("prepare blocks query");
        let mut rows = stmt.query([]).expect("query blocks");
        while let Some(row) = rows.next().expect("blocks row") {
            let aln_idx: u64 = row.get::<_, i64>(0).unwrap() as u64;
            let dup_flag: i32 = row.get(1).unwrap();
            let ovlp_flag: i32 = row.get(2).unwrap();
            let t_name: String = row.get(3).unwrap();
            let ts: u32 = row.get(4).unwrap();
            let te: u32 = row.get(5).unwrap();
            let q_name: String = row.get(6).unwrap();
            let qs: u32 = row.get(7).unwrap();
            let qe: u32 = row.get(8).unwrap();
            let orientation: u32 = row.get(9).unwrap();
            let smb: ShimmerMatchBlock = (t_name, ts, te, q_name, qs, qe, orientation);
            aln_blocks.entry(aln_idx).or_default().push(smb.clone());
            if dup_flag == 0 && ovlp_flag == 0 {
                unique_aln_blocks.entry(aln_idx).or_default().push(smb);
            }
        }
    };

    add_blocks(
        "SELECT aln_idx, dup_flag, ovlp_flag, target_name, target_start, target_end,
                query_name, query_start, query_end, orientation
         FROM blocks WHERE block_type = 0",
        &mut aln_blocks,
        &mut unique_aln_blocks,
    );
    add_blocks(
        "SELECT aln_idx, dup_flag, ovlp_flag, target_name, target_start, target_end,
                query_name, query_start, query_end, orientation
         FROM variants",
        &mut aln_blocks,
        &mut unique_aln_blocks,
    );

    (variant_records, aln_blocks, unique_aln_blocks)
}

#[allow(clippy::type_complexity)]
fn get_variant_recs_from_alnmap(
    f: BufReader<File>,
    hap_type: u8,
) -> (
    Vec<VariantRecord>,
    FxHashMap<u64, Vec<ShimmerMatchBlock>>,
    FxHashMap<u64, Vec<ShimmerMatchBlock>>,
) {
    let mut variant_records = Vec::<VariantRecord>::new();
    let mut aln_blocks = FxHashMap::<u64, Vec<ShimmerMatchBlock>>::default();
    let mut unique_aln_blocks = FxHashMap::<u64, Vec<ShimmerMatchBlock>>::default();

    f.lines().for_each(|line| {
        if let Ok(line) = line {
            if line.trim().starts_with('#') {
                return;
            };
            let fields = line.split('\t').collect::<Vec<&str>>();
            assert!(fields.len() > 3);
            let rec_type = fields[1];
            if rec_type.starts_with('V') {
                assert!(fields.len() == 15 || fields.len() == 17);
                let err_msg = format!("fail to parse on {}", line);
                let aln_block_id = fields[0].parse::<u64>().expect(&err_msg);
                let t_name = fields[2];
                let tc = fields[11].parse::<u32>().expect(&err_msg);
                let tvs = fields[13];
                let qvs = fields[14];
                variant_records.push((
                    t_name.to_string(),
                    tc,
                    tvs.len() as u32,
                    aln_block_id,
                    hap_type,
                    tvs.to_string(),
                    qvs.to_string(),
                    rec_type.to_string(),
                ));
            };

            if rec_type.starts_with('M') || rec_type.starts_with('V') {
                let err_msg = format!("fail to parse on {}", line);
                let aln_block_id = fields[0].parse::<u64>().expect(&err_msg);
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
                ));
                if rec_type == "M" || rec_type == "V" {
                    let e = unique_aln_blocks.entry(aln_block_id).or_default();
                    e.push((
                        t_name.to_string(),
                        ts,
                        te,
                        q_name.to_string(),
                        qs,
                        qe,
                        orientation,
                    ));
                }
            }
        }
    });
    (variant_records, aln_blocks, unique_aln_blocks)
}

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let use_db = args.hap0_path.ends_with(".alndb");

    let mut target_length: TargetSeqLength = if use_db {
        let conn = Connection::open_with_flags(
            &args.hap0_path,
            rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
        )
        .expect("can't open hap0 alndb");
        let mut stmt = conn
            .prepare(
                "SELECT seq_id, seq_name, length FROM sequences WHERE seq_type='target'",
            )
            .expect("prepare sequences query");
        stmt.query_map([], |row| {
            Ok((
                row.get::<_, u32>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, u32>(2)?,
            ))
        })
        .expect("query sequences")
        .map(|r| r.expect("sequence row"))
        .collect()
    } else {
        let json_path = args
            .target_len_json
            .as_ref()
            .expect("--target-len-json is required when not using .alndb input");
        let mut json_file =
            BufReader::new(File::open(Path::new(json_path)).expect("can't open the input file"));
        let mut buffer = Vec::new();
        json_file.read_to_end(&mut buffer)?;
        serde_json::from_str(&String::from_utf8_lossy(&buffer[..]))
            .expect("can't parse the target_len.json file")
    };

    target_length.sort();

    let (hap0_recs, hap0_aln_blocks, hap0_unique_aln_blocks) = if use_db {
        get_variant_recs_from_db(&args.hap0_path, 0)
    } else {
        let f = BufReader::new(File::open(Path::new(&args.hap0_path)).unwrap());
        get_variant_recs_from_alnmap(f, 0)
    };
    let (hap1_recs, hap1_aln_blocks, hap1_unique_aln_blocks) = if use_db {
        get_variant_recs_from_db(&args.hap1_path, 1)
    } else {
        let f = BufReader::new(File::open(Path::new(&args.hap1_path)).unwrap());
        get_variant_recs_from_alnmap(f, 1)
    };

    let blocks_to_intervals =
        |blocks: FxHashMap<u64, Vec<ShimmerMatchBlock>>| -> FxHashMap<String, IntervalSet<u32>> {
            let mut aln_intervals = FxHashMap::<String, IntervalSet<u32>>::default();
            blocks.into_iter().for_each(|(_block_id, records)| {
                records.into_iter().for_each(|rec| {
                    let t_name = rec.0;
                    let bgn = rec.1;
                    let end = rec.2;
                    let interval_set = aln_intervals.entry(t_name).or_default();
                    interval_set.insert(bgn..end);
                })
            });
            aln_intervals
        };

    let hap0_aln_intervals = blocks_to_intervals(hap0_aln_blocks);
    let hap1_aln_intervals = blocks_to_intervals(hap1_aln_blocks);
    let hap0_unique_aln_intervals = blocks_to_intervals(hap0_unique_aln_blocks);
    let hap1_unique_aln_intervals = blocks_to_intervals(hap1_unique_aln_blocks);

    let mut out_vcf =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("vcf")).unwrap());
    let mut out_bed =
        BufWriter::new(File::create(Path::new(&args.output_prefix).with_extension("bed")).unwrap());
    writeln!(out_vcf, "##fileformat=VCFv4.2").expect("fail to write the vcf file");
    target_length.into_iter().for_each(|(_, t_name, t_len)| {
        writeln!(out_vcf, r#"##contig=<ID={},length={}>"#, t_name, t_len)
            .expect("fail to write the vcf file");
    });
    writeln!(
        out_vcf,
        r#"##FILTER=<ID=DUP,Description="duplicated alignment block">"#
    )
    .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FILTER=<ID=OVLP,Description="overlapped alignment block">"#
    )
    .expect("fail to write the vcf file");
    writeln!(out_vcf, r#"##FILTER=<ID=NC,Description="no diploid call">"#)
        .expect("fail to write the vcf file");
    writeln!(
        out_vcf,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )
    .expect("fail to write the vcf file");

    writeln!(
        out_vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        args.sample_name
    )
    .expect("fail to write the vcf file");

    let convert_to_vcf_record = |records: &mut Vec<VariantRecord>| {
        records.sort_by_key(|v| (v.4, v.1, v.3)); // sorted by haplotype index, start reference start coordinate, aln_block
        let mut ref_bases = FxHashSet::<(u32, char)>::default();
        let mut h0alleles = FxHashMap::<u32, Vec<VariantRecord>>::default();
        let mut h1alleles = FxHashMap::<u32, Vec<VariantRecord>>::default();

        let mut al_idx_map = FxHashMap::<(u8, u64), u32>::default();
        let mut al_idx = 0_u32;

        let ref_name = records.first().unwrap().0.clone();
        let mut rec_type = Option::<String>::None;
        records.iter().for_each(|rec| {
            let (_t_name, ts, tl, aln_block_id, ht, vts, _vqs, rt) = rec;

            if rec_type.is_none() && (rt == "V_D" || rt == "V_O") {
                rec_type = Some(rt.clone());
            }

            (0..*tl).for_each(|t_pos| {
                let vts = vts.chars().collect::<Vec<_>>();
                ref_bases.insert((*ts + t_pos, vts[t_pos as usize]));
            });

            let key = (*ht, *aln_block_id);

            let al_idx2 = al_idx_map.entry(key).or_insert_with(|| {
                al_idx += 1;
                al_idx
            });

            if *ht == 0 {
                h0alleles.entry(*al_idx2).or_default().push(rec.clone());
            };
            if *ht == 1 {
                h1alleles.entry(*al_idx2).or_default().push(rec.clone());
            };
        });

        let mut ref_bases = ref_bases.into_iter().collect::<Vec<_>>();
        ref_bases.sort();
        let ref_str = String::from_iter(ref_bases.iter().map(|(_, c)| *c).collect::<Vec<_>>());
        assert!(ref_str.len() == ref_bases.len()); // make sure all bases at the same t_pos are the same, if not the vectors will have different lengths
        let ts0 = ref_bases.first().unwrap().0;
        let tl0 = ref_str.len() as u32;

        let mut query_alleles = al_idx_map
            .iter()
            .map(|(&(ht, _block_id), &al_idx)| {
                let alleles = if ht == 0 {
                    h0alleles.get(&al_idx).unwrap().clone()
                } else {
                    h1alleles.get(&al_idx).unwrap().clone()
                };
                let mut allele_str = Vec::<String>::new();
                let mut offset = 0usize;
                alleles
                    .iter()
                    .for_each(|(_t_name, ts, tl, _aln_block_id, _ht, _vts, vqs, _rt)| {
                        let end = (*ts - ts0) as usize;
                        allele_str.push(ref_str[offset..end].to_string());
                        allele_str.push(vqs.clone());
                        offset = end + *tl as usize;
                    });
                allele_str.push(ref_str[offset..].to_string());

                (al_idx, allele_str.join(""))
            })
            .collect::<Vec<_>>();

        // deduplicate query_alleles
        let mut al_idx_map = FxHashMap::<u32, u32>::default();
        let mut unique_query_alleles = FxHashMap::<String, u32>::default();
        al_idx_map.insert(0, 0);
        unique_query_alleles.entry(ref_str.clone()).or_insert(0);
        query_alleles.sort_by_key(|v| v.1.len());
        let mut new_idx = 1u32;
        query_alleles.iter().for_each(|(idx, allele)| {
            if !unique_query_alleles.contains_key(allele) {
                unique_query_alleles.insert(allele.clone(), new_idx);
                al_idx_map.insert(*idx, new_idx);
                new_idx += 1;
            } else {
                al_idx_map.insert(*idx, *unique_query_alleles.get(allele).unwrap());
            }
        });

        let query_alleles = unique_query_alleles
            .into_iter()
            .filter(|(_, v)| *v != 0)
            .map(|(qs, _)| qs.clone())
            .collect::<Vec<_>>()
            .join(",");

        let h0_al_idx = if let Some(i_set) = hap0_aln_intervals.get(&ref_name) {
            if i_set.has_overlap(ts0..ts0 + tl0) {
                let h0_al_idx = if h0alleles.is_empty() {
                    vec!["0".to_string()]
                } else {
                    let mut allele_count = FxHashMap::<u32, u32>::default();
                    h0alleles.keys().for_each(|&idx| {
                        *allele_count
                            .entry(*al_idx_map.get(&idx).unwrap())
                            .or_default() += 1
                    });
                    //println!("allele_count: {:?}", allele_count);
                    allele_count
                        .keys()
                        .map(|idx| format!("{}", idx))
                        .collect::<Vec<_>>()
                };
                if h0_al_idx.len() == 1 {
                    h0_al_idx[0].clone()
                } else {
                    ".".to_string()
                }
            } else {
                ".".to_string()
            }
        } else {
            ".".to_string()
        };
        let h1_al_idx = if let Some(i_set) = hap1_aln_intervals.get(&ref_name) {
            if i_set.has_overlap(ts0..ts0 + tl0) {
                let h1_al_idx = if h1alleles.is_empty() {
                    vec!["0".to_string()]
                } else {
                    let mut allele_count = FxHashMap::<u32, u32>::default();
                    h1alleles.keys().for_each(|&idx| {
                        *allele_count
                            .entry(*al_idx_map.get(&idx).unwrap())
                            .or_default() += 1
                    });
                    allele_count
                        .keys()
                        .map(|idx| format!("{}", idx))
                        .collect::<Vec<_>>()
                };
                if h1_al_idx.len() == 1 {
                    h1_al_idx[0].clone()
                } else {
                    ".".to_string()
                }
            } else {
                ".".to_string()
            }
        } else {
            ".".to_string()
        };
        let gt = [h0_al_idx, h1_al_idx].join("|");
        (ref_name, ts0, ref_str, query_alleles, gt, rec_type)
    };

    let mut variant_records = Vec::<VariantRecord>::new();
    variant_records.extend(hap0_recs);
    variant_records.extend(hap1_recs);

    let mut variant_records = variant_records.into_iter().collect::<Vec<_>>();
    // variant_group: represent a group of overlapped variants, ref_id, ref_start, len, REF, ALT
    let mut variant_group = Vec::<VariantRecord>::new();
    // currrent_vg_end: represent the end coordinate of the current variant group
    let mut current_vg_end = Option::<(String, u32)>::None;
    variant_records.sort();
    variant_records.into_iter().for_each(
        |(ref_name, ts, tl, block_id, ht, vts, vqs, rec_type)| {
            if let Some(current_vg_end) = current_vg_end.clone() {
                //println!("{} {} {} {} {:?} {}", ref_name, ts, tl, ts + tl ,  current_vg_end, variant_group.len()  );
                if ref_name == current_vg_end.0 && ts < current_vg_end.1 {
                    variant_group.push((
                        ref_name.clone(),
                        ts,
                        tl,
                        block_id,
                        ht,
                        vts,
                        vqs,
                        rec_type,
                    ));
                } else if !variant_group.is_empty() {
                    //println!("X {} {} {} {} {:?}", ref_name, ts, tl, variant_group.len(), variant_group);
                    let (vcf_rec_ref_name, ts0, ref_str, query_alleles, gt, g_rec_type) =
                        convert_to_vcf_record(&mut variant_group);
                    let rt = if let Some(g_rec_type) = g_rec_type {
                        if g_rec_type == "V_D" {
                            "DUP"
                        } else if g_rec_type == "V_O" {
                            "OVLP"
                        } else {
                            "PASS"
                        }
                    } else {
                        "PASS"
                    };
                    let rt = if rt == "PASS" && gt.contains('.') {
                        "NC"
                    } else {
                        rt
                    };
                    let qv: u32 = if rt != "PASS" { 30 } else { 40 };
                    //writeln!(
                    //    out_vcf, "{:?}", variant_group
                    //);
                    writeln!(
                        out_vcf,
                        "{}\t{}\t.\t{}\t{}\t{}\t{}\t.\tGT\t{}",
                        vcf_rec_ref_name,
                        ts0 + 1,
                        ref_str,
                        query_alleles,
                        qv,
                        rt,
                        gt,
                    )
                    .expect("fail to write the vcf file");
                    variant_group.clear();
                    variant_group.push((
                        ref_name.clone(),
                        ts,
                        tl,
                        block_id,
                        ht,
                        vts,
                        vqs,
                        rec_type,
                    ));
                }
            } else {
                variant_group.push((
                    ref_name.clone(),
                    ts,
                    tl,
                    block_id,
                    ht,
                    vts,
                    vqs,
                    rec_type,
                ));
                current_vg_end = Some((ref_name, ts + tl));
                return;
            }
            current_vg_end = Some((ref_name.clone(), ts + tl));
        },
    );
    if !variant_group.is_empty() {
        // println!("X {} {} {} {}", ref_name, ts, tl, variant_group.len());
        let (vcf_rec_ref_name, ts0, ref_str, query_alleles, gt, g_rec_type) =
            convert_to_vcf_record(&mut variant_group);
        let rt = if let Some(g_rec_type) = g_rec_type {
            if g_rec_type == "V_D" {
                "DUP"
            } else if g_rec_type == "V_O" {
                "OVLP"
            } else {
                "PASS"
            }
        } else {
            "PASS"
        };
        let qv: u32 = if rt != "PASS" { 30 } else { 40 };
        writeln!(
            out_vcf,
            "{}\t{}\t.\t{}\t{}\t{}\t{}\t.\tGT\t{}",
            vcf_rec_ref_name,
            ts0 + 1,
            ref_str,
            query_alleles,
            qv,
            rt,
            gt,
        )
        .expect("fail to write the vcf file");
    };
    variant_group.clear();

    let merge_intervals =
        |intervals: FxHashMap<String, IntervalSet<u32>>| -> FxHashMap<String, IntervalSet<u32>> {
            intervals
                .into_iter()
                .flat_map(|(t_name, i_set)| {
                    let mut intervals = i_set
                        .unsorted_iter()
                        .map(|range| (range.start, range.end))
                        .collect::<Vec<_>>();
                    if intervals.is_empty() {
                        return None;
                    }
                    intervals.sort();
                    let mut merged_intervals = IntervalSet::<u32>::new();
                    let mut current_range = intervals.first().unwrap().to_owned();
                    intervals.into_iter().for_each(|(bgn, end)| {
                        if bgn <= current_range.1 && end > current_range.1 {
                            current_range.1 = end;
                        } else if bgn > current_range.1 {
                            merged_intervals.insert(current_range.0..current_range.1);
                            current_range = (bgn, end);
                        }
                    });
                    merged_intervals.insert(current_range.0..current_range.1);
                    Some((t_name, merged_intervals))
                })
                .collect()
        };

    let hap0_aln_merged_intervals: FxHashMap<String, IntervalSet<u32>> =
        merge_intervals(hap0_unique_aln_intervals);

    let hap1_aln_merged_intervals: FxHashMap<String, IntervalSet<u32>> =
        merge_intervals(hap1_unique_aln_intervals);

    let mut t_names = hap0_aln_merged_intervals
        .keys()
        .cloned()
        .collect::<Vec<_>>();
    t_names.sort();

    t_names.into_iter().for_each(|t_name| {
        let hap0_aln_merged_intervals = hap0_aln_merged_intervals.get(&t_name).unwrap();
        if let Some(hap1_aln_merged_intervals) = hap1_aln_merged_intervals.get(&t_name) {
            let mut intervals = hap0_aln_merged_intervals
                .unsorted_iter()
                .map(|range| (range.start, range.end))
                .collect::<Vec<_>>();
            if intervals.is_empty() {
                return;
            }
            intervals.sort();
            let msg = "can't write the output bed file";
            intervals.into_iter().for_each(|(bgn, end)| {
                hap1_aln_merged_intervals.iter(bgn..end).for_each(|range| {
                    let (bgn1, end1) = (range.start, range.end);
                    if bgn1 < bgn && end1 < end {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn, end1).expect(msg);
                    } else if bgn1 < bgn && end <= end1 {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn, end).expect(msg);
                    } else if bgn <= bgn1 && end1 < end {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn1, end1).expect(msg);
                    } else if bgn <= bgn1 && end <= end1 {
                        writeln!(out_bed, "{}\t{}\t{}", t_name, bgn1, end).expect(msg);
                    }
                });
            });
        }
    });

    Ok(())
}
