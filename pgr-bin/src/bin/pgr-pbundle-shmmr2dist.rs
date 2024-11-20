const VERSION_STRING: &str = env!("VERSION_STRING");
use clap::{self, CommandFactory, Parser};
use kodama::{linkage, Method};
use rustc_hash::{FxHashMap, FxHashSet};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::{fs::File, path};

/// Generate alignment scores between sequences using bundle decomposition from a principal bundle bed file
#[derive(Parser, Debug)]
#[clap(name = "pgr-pbundle-bed2dist")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// the path to the principal bundle bed file
    idx_file_path: String,
    /// the prefix of the output file
    output_prefix: String,
}

type Smps = Vec<(String, u32, u32, u8)>; // shmmr_string, bgn, end, orientation

fn align_smps(smps0: &Smps, smps1: &Smps) -> (f32, usize, usize, i64, isize) {
    // return: dist, diff_len, max_len, best_score, best_offset
    let mut smp_to_frags0 = FxHashMap::<(String, u8), Vec<(u32, u32)>>::default();
    let mut smp_to_frags1 = FxHashMap::<(String, u8), Vec<(u32, u32)>>::default();
    let mut all_smps = FxHashSet::<(String, u8)>::default();
    let mut length0 = 0_u32;
    let mut length1 = 0_u32;
    smps0.iter().for_each(|(frag_id, bgn, end, orientation)| {
        let e = smp_to_frags0
            .entry((frag_id.clone(), *orientation))
            .or_default();
        e.push((*bgn, *end));
        all_smps.insert((frag_id.clone(), *orientation));
        length0 += *end - *bgn;
    });

    smps1.iter().for_each(|(frag_id, bgn, end, orientation)| {
        let e = smp_to_frags1
            .entry((frag_id.clone(), *orientation))
            .or_default();
        e.push((*bgn, *end));
        all_smps.insert((frag_id.clone(), *orientation));
        length1 += *end - *bgn;
    });

    let mut match_score = 0_i32;
    let mut diff_len = 0_u32;
    let mut offsets = Vec::<(i32, u32)>::new();
    for smp in all_smps {
        if smp_to_frags0.contains_key(&smp) && smp_to_frags1.contains_key(&smp) {
            let frags0 = &smp_to_frags0[&smp];
            let frags1 = &smp_to_frags1[&smp];
            let l0 = frags0.iter().map(|v| v.1 - v.0).sum::<u32>();
            let l1 = frags1.iter().map(|v| v.1 - v.0).sum::<u32>();

            if frags0.len() == frags1.len() {
                match_score += (l0 + l1) as i32;
                for i in 0..frags0.len() {
                    let (bgn0, _end0) = frags0[i];
                    let (bgn1, _end1) = frags1[i];
                    if frags0.len() == 1 {
                        // only use unique ones
                        offsets.push((bgn1 as i32 - bgn0 as i32, l0 + l1));
                    }
                }
            } else {
                match_score += u32::min(l0, l1) as i32 - l0.abs_diff(l1) as i32;
                diff_len += l0.abs_diff(l1);
            };
        } else if smp_to_frags0.contains_key(&smp) {
            let frags0 = &smp_to_frags0[&smp];
            let l0 = frags0.iter().map(|v| v.1 - v.0).sum::<u32>();
            match_score -= l0 as i32;
            diff_len += l0;
        } else if smp_to_frags1.contains_key(&smp) {
            let frags1 = &smp_to_frags1[&smp];
            let l1 = frags1.iter().map(|v| v.1 - v.0).sum::<u32>();
            match_score -= l1 as i32;
            diff_len += l1;
        }
    }

    offsets.sort();
    let mut offset_clusters = Vec::<Vec<(i32, u32)>>::new();
    const MERGE_LENGTH: i32 = 16;

    let mut current_cluster = Vec::<(i32, u32)>::new();
    let mut current_offset: Option<i32> = None;
    for (offset, length) in offsets {
        if let Some(last_offset) = current_offset {
            if offset - last_offset < MERGE_LENGTH {
                current_offset = Some(offset);
                current_cluster.push((offset, length));
            } else {
                offset_clusters.push(current_cluster.clone());
                current_cluster.clear();
                current_offset = Some(offset);
                current_cluster.push((offset, length));
            }
        } else {
            current_offset = Some(offset);
            current_cluster.push((offset, length));
        };
    }
    if !current_cluster.is_empty() {
        offset_clusters.push(current_cluster);
    };
    if offset_clusters.is_empty() {
        return (1.0, diff_len as usize, usize::MAX, match_score as i64, 0);
    };
    offset_clusters.sort_by_key(|v| -(v.len() as isize));
    let best_cluster = &offset_clusters[0];
    let ave_offset = best_cluster.iter().map(|v| v.0).sum::<i32>() / best_cluster.len() as i32;
    let max_len = best_cluster.iter().map(|v| v.1).sum::<u32>();
    let dist = 1.0 - 0.5 * (match_score as f32 / (length0 + length1) as f32 + 1.0_f32);
    // return: dist, diff_len, max_len, best_score, best_offset
    (
        dist,
        diff_len as usize,
        max_len as usize,
        match_score as i64,
        ave_offset as isize,
    )
}

type Contigs = FxHashMap<u32, (String, String, u32)>; // contig_id -> contig_name, source, length
type FragMap = FxHashMap<String, (u32, u32, u32, u8)>; // shmmr string -> seq_id, bgn, end, orientation
type CtgToFrags = FxHashMap<String, Smps>; // contig_id -> shmmr_string, bgn, end, orientation

fn main() -> Result<(), std::io::Error> {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();
    let shmmr_idx_filename = path::Path::new(&args.idx_file_path);
    let shmmr_idx_file =
        BufReader::new(File::open(shmmr_idx_filename).expect("can't open the bed file"));
    //let mut ctg_data = FxHashMap::<String, Vec<_>>::default();
    let shmmr_idx_file_parse_err_msg = "shmmr idx file parsing error";
    let mut contigs = Contigs::default();
    let mut frag_map = FragMap::default();
    let mut ctg_to_frags = CtgToFrags::default();
    shmmr_idx_file.lines().for_each(|line| {
        let line = line.unwrap().trim().to_string();
        if line.is_empty() {
            return;
        }
        let t = &line[0..1];
        match t {
            "#" | "K" => {}
            "C" => {
                let idx_fields = line.split('\t').collect::<Vec<&str>>();
                let ctg_id: u32 = idx_fields[1].parse().expect(shmmr_idx_file_parse_err_msg);
                let ctg_name: String = idx_fields[2].to_string();
                let source = idx_fields[3].to_string();
                let length: u32 = idx_fields[4].parse().expect(shmmr_idx_file_parse_err_msg);
                contigs.entry(ctg_id).or_insert((ctg_name, source, length));
            }
            "F" => {
                let idx_fields = line.split('\t').collect::<Vec<&str>>();
                let frag_id = idx_fields[1].to_string();
                let seq_id: u32 = idx_fields[3].parse().expect(shmmr_idx_file_parse_err_msg);
                let bgn: u32 = idx_fields[4].parse().expect(shmmr_idx_file_parse_err_msg);
                let end: u32 = idx_fields[5].parse().expect(shmmr_idx_file_parse_err_msg);
                let orientation: u8 = idx_fields[6].parse().expect(shmmr_idx_file_parse_err_msg);
                frag_map
                    .entry(frag_id.clone())
                    .or_insert((seq_id, bgn, end, orientation));
                let ctg = contigs
                    .get(&seq_id)
                    .expect(shmmr_idx_file_parse_err_msg)
                    .0
                    .clone();
                let e = ctg_to_frags.entry(ctg.clone()).or_default();
                e.push((frag_id, bgn, end, orientation));
            }
            _ => panic!("{}", shmmr_idx_file_parse_err_msg),
        }
    });

    let ctg_to_frags = ctg_to_frags
        .into_iter()
        .map(|(k, mut v)| {
            v.sort_by_key(|k| k.1);
            (k, v)
        })
        .collect::<Vec<_>>();

    let n_ctg = ctg_to_frags.len();

    let out_path = Path::new(&args.output_prefix).with_extension("dist");
    let mut out_file = BufWriter::new(File::create(out_path).expect("can't create the dist file"));

    let mut dist_map = FxHashMap::<(usize, usize), f32>::default();
    let mut offset_map = FxHashMap::<(usize, usize), isize>::default();
    let mut min_dist = 0.0_f32;
    let mut max_dist = 1.0_f32;
    (0..n_ctg)
        .flat_map(|ctg_idx0| (0..n_ctg).map(move |ctg_idx1| (ctg_idx0, ctg_idx1)))
        .for_each(|(ctg_idx0, ctg_idx1)| {
            if ctg_idx0 > ctg_idx1 {
                return;
            };
            let (ctg0, ctg0_smps) = &ctg_to_frags[ctg_idx0];
            let (ctg1, ctg1_smps) = &ctg_to_frags[ctg_idx1];
            let (dist, diff_len, max_len, best_score, best_offset) =
                align_smps(ctg0_smps, ctg1_smps);
            writeln!(
                out_file,
                "{} {} {} {} {} {} {}",
                ctg0, ctg1, dist, diff_len, max_len, best_score, best_offset
            )
            .expect("writing error");

            if ctg_idx1 != ctg_idx0 {
                writeln!(
                    out_file,
                    "{} {} {} {} {} {} {}",
                    ctg1, ctg0, dist, diff_len, max_len, best_score, -best_offset
                )
                .expect("writing error");
                min_dist = if dist < min_dist { dist } else { min_dist };
                max_dist = if dist > max_dist { dist } else { max_dist };
                dist_map.insert((ctg_idx0, ctg_idx1), dist);
                offset_map.insert((ctg_idx0, ctg_idx1), best_offset);
                offset_map.insert((ctg_idx1, ctg_idx0), -best_offset);
            }
        });

    let w = max_dist - min_dist + 0.01;
    dist_map.iter_mut().for_each(|(_k, v)| {
        *v = (*v - min_dist + 0.01) / w;
    });
    let mut dist_mat = vec![];
    (0..n_ctg - 1).for_each(|i| {
        (i + 1..n_ctg).for_each(|j| {
            dist_mat.push(*dist_map.get(&(i, j)).unwrap());
        })
    });
    let dend = linkage(&mut dist_mat, n_ctg, Method::Average);

    let steps = dend.steps().to_vec();
    let mut node_data = FxHashMap::<usize, (String, Vec<usize>, f32)>::default();
    (0..n_ctg).for_each(|ctg_idx| {
        node_data.insert(ctg_idx, (format!("{}", ctg_idx), vec![ctg_idx], 0.0_f32));
    });

    let mut last_node_id = 0_usize;
    steps.iter().enumerate().for_each(|(c, s)| {
        let (node_string1, nodes1, height1) = node_data.remove(&s.cluster1).unwrap();
        let (node_string2, nodes2, height2) = node_data.remove(&s.cluster2).unwrap();
        let new_node_id = c + n_ctg;
        let mut nodes = Vec::<usize>::new();
        let new_node_string = if nodes1.len() > nodes2.len() {
            nodes.extend(nodes1);
            nodes.extend(nodes2);
            format!(
                "({}:{}, {}:{})",
                node_string1,
                s.dissimilarity - height1,
                node_string2,
                s.dissimilarity - height2
            )
        } else {
            nodes.extend(nodes2);
            nodes.extend(nodes1);
            format!(
                "({}:{}, {}:{})",
                node_string2,
                s.dissimilarity - height2,
                node_string1,
                s.dissimilarity - height1
            )
        };
        node_data.insert(new_node_id, (new_node_string, nodes, s.dissimilarity));
        last_node_id = new_node_id;
    });

    let mut tree_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("nwk"))
            .expect("can't create the nwk file"),
    );

    let emptyp_string = ("".to_string(), vec![], 0.0);
    let (tree_string, nodes, _) = node_data.get(&last_node_id).unwrap_or(&emptyp_string);
    writeln!(tree_file, "{};", tree_string).expect("can't write the nwk file");

    let mut dendrogram_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("ddg"))
            .expect("can't create the dendrogram file"),
    );

    let mut offset_file = BufWriter::new(
        File::create(Path::new(&args.output_prefix).with_extension("offset"))
            .expect("can't create the alignment offset file"),
    );

    let mut node_position_size = FxHashMap::<usize, ((f32, f32), usize)>::default();
    let mut position = 0.0_f32;
    let mut offset = 0_isize;
    let mut p_idx: Option<usize> = None;
    let mut offset_group = Vec::<_>::new();
    let mut group_min_offset = 100000_isize;
    nodes.iter().for_each(|&ctg_idx| {
        node_position_size.insert(ctg_idx, ((position, 0.0), 1));
        writeln!(
            dendrogram_file,
            "L\t{}\t{}",
            ctg_idx, ctg_to_frags[ctg_idx].0
        )
        .expect("can't write the dendrogram file");
        position += 1.0;
        if let Some(p_idx) = p_idx {
            let (idx0, idx1) = if p_idx < ctg_idx {
                (p_idx, ctg_idx)
            } else {
                (ctg_idx, p_idx)
            };
            if *dist_map.get(&(idx0, idx1)).unwrap_or(&1.0) < 0.25 {
                offset += *offset_map.get(&(p_idx, ctg_idx)).unwrap_or(&0);
                offset_group.push((ctg_idx, offset));
                if offset < group_min_offset {
                    group_min_offset = offset;
                };
            } else {
                offset_group.iter().for_each(|&(ctg_idx, offset)| {
                    writeln!(
                        offset_file,
                        "{}\t{}",
                        ctg_to_frags[ctg_idx].0,
                        offset - group_min_offset
                    )
                    .expect("can't write the offset file");
                });
                group_min_offset = 100000_isize;
                offset_group.clear();
                offset = 0;
            }
        } else {
            offset_group.push((ctg_idx, offset));
        };
        p_idx = Some(ctg_idx)
    });

    offset_group.iter().for_each(|&(ctg_idx, offset)| {
        writeln!(
            offset_file,
            "{}\t{}",
            ctg_to_frags[ctg_idx].0,
            offset - group_min_offset
        )
        .expect("can't write the offset file");
    });

    steps.into_iter().enumerate().for_each(|(c, s)| {
        let ((pos0, _), size0) = *node_position_size.get(&s.cluster1).unwrap();
        let ((pos1, _), size1) = *node_position_size.get(&s.cluster2).unwrap();

        let pos = ((size0 as f32) * pos0 + (size1 as f32) * pos1) / ((size0 + size1) as f32);
        writeln!(
            dendrogram_file,
            "I\t{}\t{}\t{}\t{}\t{}",
            c + n_ctg,
            s.cluster1,
            s.cluster2,
            s.size,
            s.dissimilarity,
        )
        .expect("can't write the dendrogram file");
        node_position_size.insert(c + n_ctg, ((pos, s.dissimilarity), s.size));
    });
    let mut node_positions = node_position_size
        .into_iter()
        .collect::<Vec<(usize, ((f32, f32), usize))>>();
    node_positions.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    node_positions
        .into_iter()
        .for_each(|(vid, ((pos, h), size))| {
            writeln!(dendrogram_file, "P\t{}\t{}\t{}\t{}", vid, pos, h, size)
                .expect("can't write the dendrogram file");
        });
    Ok(())
}
