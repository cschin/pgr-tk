const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
// use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Merge svcnd from multiple *.svcnd.bed files into one and compute the merged regions
/// It is useful to identify unique bed regions to one specific haplotype
#[derive(Parser, Debug)]
#[clap(name = "pgr-merge-svcnd-bed")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// path to the file contain the input bed files, each line should be "label<tab>input file path"
    #[clap(long, short)]
    input_files: String,
    /// the path of the output files
    #[clap(long, short)]
    output_path: String,
    /// number of threads used in parallel (more memory usage), default to "0" using all CPUs available or the number set by RAYON_NUM_THREADS
    #[clap(long, default_value_t = 0)]
    number_of_thread: usize,
}

type Interval = ((u32, u32), (String, String));
fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.number_of_thread)
        .build_global()
        .unwrap();

    let input_files = BufReader::new(File::open(Path::new(&args.input_files)).unwrap());

    let input_files = input_files
        .lines()
        .flat_map(|line| {
            if let Ok(line) = line {
                let rec = line.trim().split('\t').collect::<Vec<&str>>();
                assert!(rec.len() >= 2);
                Some((rec[0].to_string(), rec[1].to_string()))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    let mut interval_collection =
        FxHashMap::<String, Vec<((u32, u32), (String, String))>>::default();
    input_files.iter().for_each(|(label, path)| {
        let bed_reader = BufReader::new(File::open(Path::new(path)).unwrap());
        bed_reader.lines().for_each(|line| {
            if let Ok(line) = line {
                if line.starts_with('#') {
                    return;
                };
                let err_msg = format!("fail to parse on {}", line);
                let fields = line.split('\t').collect::<Vec<&str>>();
                let chr = fields[0].to_string();
                let bgn = fields[1].parse::<u32>().expect(&err_msg);
                let end = fields[2].parse::<u32>().expect(&err_msg);
                let annotation = fields[3].to_string();
                let e = interval_collection.entry(chr).or_insert_with(Vec::new);
                e.push(((bgn, end), (label.clone(), annotation)));
            }
        });
    });

    let group_intervals = |intervals: &mut Vec<Interval>| -> Vec<(u32, u32, Vec<Interval>)> {
        let mut interval_groups = Vec::<(u32, u32, Vec<Interval>)>::new();
        if intervals.is_empty() {
            return interval_groups;
        }

        intervals.sort();
        let (mut current_bgn, mut current_end) = intervals.first().unwrap().0;

        let mut current_groups = Vec::<Interval>::new();
        intervals.iter().for_each(|(interval, payload)| {
            if current_end < interval.0 {
                interval_groups.push((current_bgn, current_end, current_groups.clone()));
                current_groups.clear();
                current_groups.push((*interval, payload.clone()));
                (current_bgn, current_end) = *interval;
            } else {
                current_groups.push((*interval, payload.clone()));
                if current_end < interval.1 {
                    current_end = interval.1;
                } 
            }
        });
        if !current_groups.is_empty() {
            interval_groups.push((current_bgn, current_end, current_groups.clone()));
        }
        interval_groups
    };

    let mut out_bed = BufWriter::new(File::create(Path::new(&args.output_path)).unwrap());
    let mut keys = interval_collection.keys().cloned().collect::<Vec<_>>();
    keys.sort();
    keys.into_iter().for_each(|key| {
        let intervals = interval_collection.get_mut(&key).unwrap();
        let interval_groups = group_intervals(intervals);
        interval_groups.into_iter().for_each(|intervals| {
            if intervals.2.is_empty() {
                return;
            }
            let itvl_group_bgn = intervals.0;
            let itvl_group_end = intervals.1;
            if itvl_group_bgn > itvl_group_end {
                return;
            };

            let mut label_count = FxHashMap::<String, u32>::default();
            let mut total_interval_counts = 0u32;
            intervals.2.iter().for_each(|(_interval, payload)| {
                let e = label_count.entry(payload.0.clone()).or_default();
                *e += 1;
                total_interval_counts += 1;
            });

            writeln!(
                out_bed,
                "{}\t{}\t{}\tmerged:{}:{}",
                key,
                itvl_group_bgn,
                itvl_group_end,
                label_count.len(),
                total_interval_counts
            )
            .expect("unable to write the output file");

            intervals.2.iter().for_each(|(interval, payload)| {
                let number_haplotype = label_count.len();
                let e = label_count.entry(payload.0.clone()).or_default();
                writeln!(
                    out_bed,
                    "{}\t{}\t{}\t{}:{}:{}-{}:{}:{}",
                    key,
                    interval.0,
                    interval.1,
                    payload.0,
                    payload.1,
                    itvl_group_bgn,
                    itvl_group_end,
                    number_haplotype,
                    *e,
                )
                .expect("unable to write the output file");
            });
        });
    });
}
