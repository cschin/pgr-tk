pub mod seq_index_db {
    use pgr_db::aln::{self, HitPair};
    use pgr_db::seq_db;
    use pgr_db::shmmrutils::{sequence_to_shmmrs, DeltaPoint, ShmmrSpec};
    use pgr_db::{agc_io, fasta_io};
    use rayon::prelude::*;
    use rayon::vec;
    use rustc_hash::{FxHashMap, FxHashSet};
    use serde::{Deserialize, Serialize};
    use serde_json::{json, Value};
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader, BufWriter, Stderr, Write};
    use std::num;

    pub type SmpBundleTuple = ((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>);
    pub type SmpsWithBundleLabel = Vec<SmpBundleTuple>;

    pub struct SeqIndexDB {
        /// Rust internal: store the specification of the shmmr specifcation
        pub shmmr_spec: Option<ShmmrSpec>,
        /// Rust internal: store the sequences
        pub seq_db: Option<seq_db::CompactSeqDB>,
        /// Rust internal: store the agc file and the index
        pub agc_db: Option<(agc_io::AGCFile, seq_db::ShmmrToFrags)>,
        /// a dictionary maps (ctg_name, source) -> (id, len)
        pub seq_index: Option<HashMap<(String, Option<String>), (u32, u32)>>,
        /// a dictionary maps id -> (ctg_name, source, len)
        pub seq_info: Option<HashMap<u32, (String, Option<String>, u32)>>,
    }

    impl SeqIndexDB {
        /// constructor, take no argument
        pub fn new() -> Self {
            SeqIndexDB {
                seq_db: None,
                agc_db: None,
                shmmr_spec: None,
                seq_index: None,
                seq_info: None,
            }
        }

        pub fn load_from_agc_index(&mut self, prefix: String) -> Result<(), std::io::Error> {
            // let (shmmr_spec, new_map) = seq_db::read_mdb_file(prefix.to_string() + ".mdb").unwrap();
            let (shmmr_spec, new_map) =
                seq_db::read_mdb_file_parallel(prefix.to_string() + ".mdb").unwrap();
            let agc_file = agc_io::AGCFile::new(prefix.to_string() + ".agc")?;
            self.agc_db = Some((agc_file, new_map));
            self.seq_db = None;
            self.shmmr_spec = Some(shmmr_spec);

            let mut seq_index = HashMap::<(String, Option<String>), (u32, u32)>::new();
            let mut seq_info = HashMap::<u32, (String, Option<String>, u32)>::new();
            let midx_file = BufReader::new(File::open(prefix.to_string() + ".midx")?);
            midx_file
                .lines()
                .into_iter()
                .try_for_each(|line| -> Result<(), std::io::Error> {
                    let line = line.unwrap();
                    let mut line = line.as_str().split("\t");
                    let sid = line.next().unwrap().parse::<u32>().unwrap();
                    let len = line.next().unwrap().parse::<u32>().unwrap();
                    let ctg_name = line.next().unwrap().to_string();
                    let source = line.next().unwrap().to_string();
                    seq_index.insert((ctg_name.clone(), Some(source.clone())), (sid, len));
                    seq_info.insert(sid, (ctg_name, Some(source), len));
                    Ok(())
                })?;
            self.seq_index = Some(seq_index);
            self.seq_info = Some(seq_info);
            Ok(())
        }

        pub fn load_from_seq_list(
            &mut self,
            seq_list: Vec<(String, Vec<u8>)>,
            source: Option<&str>,
            w: u32,
            k: u32,
            r: u32,
            min_span: u32,
        ) -> () {
            let spec = ShmmrSpec {
                w,
                k,
                r,
                min_span,
                sketch: false,
            };
            let source = Some(source.unwrap().to_string());
            let mut sdb = seq_db::CompactSeqDB::new(spec.clone());
            let seq_vec = seq_list
                .into_iter()
                .enumerate()
                .map(|(sid, v)| (sid as u32, source.clone(), v.0, v.1))
                .collect::<Vec<(u32, Option<String>, String, Vec<u8>)>>();
            sdb.load_seqs_from_seq_vec(&seq_vec);

            self.shmmr_spec = Some(spec);
            let mut seq_index = HashMap::<(String, Option<String>), (u32, u32)>::new();
            let mut seq_info = HashMap::<u32, (String, Option<String>, u32)>::new();
            sdb.seqs.iter().for_each(|v| {
                seq_index.insert((v.name.clone(), v.source.clone()), (v.id, v.len as u32));
                seq_info.insert(v.id, (v.name.clone(), v.source.clone(), v.len as u32));
            });
            self.seq_index = Some(seq_index);
            self.seq_info = Some(seq_info);
            self.seq_db = Some(sdb);
            self.agc_db = None;
            ()
        }

        fn get_vertex_map_from_priciple_bundles(
            &self,
            pb: Vec<Vec<(u64, u64, u8)>>,
        ) -> FxHashMap<(u64, u64), (usize, u8, usize)> {
            // conut segment for filtering, some undirectional seg may have both forward and reverse in the principle bundles
            let mut seg_count = FxHashMap::<(u64, u64), usize>::default();
            pb.iter().for_each(|bundle| {
                bundle.iter().for_each(|v| {
                    *seg_count.entry((v.0, v.1)).or_insert(0) += 1;
                })
            });

            pb.iter()
                .enumerate()
                .flat_map(|(bundle_id, path)| {
                    path.iter()
                        .enumerate()
                        .filter(|(_, &v)| *seg_count.get(&(v.0, v.1)).unwrap_or(&0) == 1)
                        .map(|(p, v)| ((v.0, v.1), (bundle_id, v.2, p)))
                        .collect::<Vec<((u64, u64), (usize, u8, usize))>>()
                })
                .collect()
        }

        pub fn get_seq(&self, sample_name: String, ctg_name: String) -> Vec<u8> {
            if self.agc_db.is_some() {
                self.agc_db
                    .as_ref()
                    .unwrap()
                    .0
                    .get_seq(sample_name, ctg_name)
            } else {
                let &(sid, _) = self
                    .seq_index
                    .as_ref()
                    .unwrap()
                    .get(&(ctg_name, Some(sample_name)))
                    .unwrap();
                self.seq_db.as_ref().unwrap().get_seq_by_id(sid)
            }
        }

        pub fn get_principal_bundles(
            &self,
            min_count: usize,
            path_len_cutoff: usize,
        ) -> Vec<Vec<(u64, u64, u8)>> {
            let frag_map = if self.agc_db.is_some() {
                &self.agc_db.as_ref().unwrap().1
            } else {
                &self.seq_db.as_ref().unwrap().frag_map
            };

            let adj_list = seq_db::frag_map_to_adj_list(frag_map, min_count as usize);

            seq_db::get_principal_bundles_from_adj_list(frag_map, &adj_list, path_len_cutoff)
                .0
                .into_iter()
                .map(|p| p.into_iter().map(|v| (v.0, v.1, v.2)).collect())
                .collect::<Vec<Vec<(u64, u64, u8)>>>()
        }

        pub fn get_principal_bundle_decomposition(
            &self,
            min_count: usize,
            path_len_cutoff: usize,
        ) -> (
            Vec<(usize, usize, Vec<(u64, u64, u8)>)>,
            Vec<(
                u32,
                SmpsWithBundleLabel,
            )>,
        ) {
            fn get_smps(seq: Vec<u8>, shmmr_spec: &ShmmrSpec) -> Vec<(u64, u64, u32, u32, u8)> {
                let shmmrs = sequence_to_shmmrs(0, &seq, &shmmr_spec, false);
                seq_db::pair_shmmrs(&shmmrs)
                    .par_iter()
                    .map(|(s0, s1)| {
                        let p0 = s0.pos() + 1;
                        let p1 = s1.pos() + 1;
                        let s0 = s0.x >> 8;
                        let s1 = s1.x >> 8;
                        if s0 < s1 {
                            (s0, s1, p0, p1, 0_u8)
                        } else {
                            (s1, s0, p0, p1, 1_u8)
                        }
                    })
                    .collect::<Vec<(u64, u64, u32, u32, u8)>>()
            }

            let pb = self.get_principal_bundles(min_count, path_len_cutoff);
            //println!("DBG: # bundles {}", pb.len());

            let mut vertex_to_bundle_id_direction_pos =
                self.get_vertex_map_from_priciple_bundles(pb.clone()); //not efficient but it is PyO3 limit now

            let seqid_smps: Vec<(u32, Vec<(u64, u64, u32, u32, u8)>)> = self
                .seq_info
                .clone()
                .unwrap_or_default()
                .iter()
                .map(|(sid, data)| {
                    let (ctg_name, source, _) = data;
                    let source = source.clone().unwrap();
                    let seq = self.get_seq(source.clone(), ctg_name.clone());
                    (*sid, get_smps(seq, &self.shmmr_spec.clone().unwrap()))
                })
                .collect();

            // data for reordering the bundles and for re-ordering them along the sequnces
            let mut bundle_id_to_directions = FxHashMap::<usize, Vec<u32>>::default();
            let mut bundle_id_to_orders = FxHashMap::<usize, Vec<f32>>::default();
            seqid_smps.iter().for_each(|(_sid, smps)| {
                let mut bundle_visited = FxHashSet::<usize>::default();
                smps.iter().enumerate().for_each(|(order, v)| {
                    if let Some(bid) = vertex_to_bundle_id_direction_pos.get(&(v.0, v.1)) {
                        if !bundle_visited.contains(&bid.0) {
                            bundle_id_to_orders
                                .entry(bid.0)
                                .or_insert(vec![])
                                .push(order as f32);
                            bundle_visited.insert(bid.0);
                        }
                        let direction = match bid.1 == v.4 {
                            true => 0,
                            false => 1,
                        };
                        bundle_id_to_directions
                            .entry(bid.0)
                            .or_insert(vec![])
                            .push(direction);
                    }
                })
            });

            // determine the bundles' overall orders and directions by consensus voting
            let mut bundle_mean_order_direction = (0..pb.len())
                .into_iter()
                .map(|bid| {
                    if let Some(orders) = bundle_id_to_orders.get(&bid) {
                        let sum: f32 = orders.into_iter().sum();
                        let mean_ord = sum / (orders.len() as f32);
                        let mean_ord = mean_ord as usize;
                        let directions = bundle_id_to_directions.get(&bid).unwrap();
                        let dir_sum = directions.iter().sum::<u32>() as usize;
                        let direction = if dir_sum < (directions.len() >> 1) {
                            0_u8
                        } else {
                            1_u8
                        };
                        (mean_ord, bid, direction)
                    } else {
                        let mean_ord = usize::MAX;
                        (mean_ord, bid, 0)
                    }
                })
                .collect::<Vec<(usize, usize, u8)>>();

            //println!("DBG: length of bundle_mean_order_direction: {}", bundle_mean_order_direction.len());

            bundle_mean_order_direction.sort();
            // re-order the principal bundles
            let principal_bundles = bundle_mean_order_direction
                .iter()
                .map(|(ord, bid, direction)| {
                    let bundle = if *direction == 1 {
                        let rpb = pb[*bid]
                            .iter()
                            .rev()
                            .map(|v| (v.0, v.1, 1 - v.2))
                            .collect::<Vec<(u64, u64, u8)>>();
                        rpb.iter().enumerate().for_each(|(p, v)| {
                            vertex_to_bundle_id_direction_pos.insert((v.0, v.1), (*bid, v.2, p));
                            // overide what in the hashpmap
                        });
                        rpb
                    } else {
                        pb[*bid].clone()
                    };

                    (*bid, *ord, bundle)
                })
                .collect::<Vec<(usize, usize, Vec<(u64, u64, u8)>)>>();

            // loop through each sequnece and generate the decomposition for the sequence
            let seqid_smps_with_bundle_id_seg_direction = seqid_smps
                .iter()
                .map(|(sid, smps)| {
                    let smps = smps
                        .into_iter()
                        .map(|v| {
                            let seg_match = if let Some(m) =
                                vertex_to_bundle_id_direction_pos.get(&(v.0, v.1))
                            {
                                Some(*m)
                            } else {
                                None
                            };
                            (*v, seg_match)
                        })
                        .collect::<SmpsWithBundleLabel>();
                    (*sid, smps)
                })
                .collect::<Vec<(u32, SmpsWithBundleLabel)>>();

            (principal_bundles, seqid_smps_with_bundle_id_seg_direction)
        }
    }

    pub fn query_fragment_to_hps(
        seq_db: &SeqIndexDB,
        seq: Vec<u8>,
        penality: f32,
        max_count: Option<u32>,
        max_count_query: Option<u32>,
        max_count_target: Option<u32>,
        max_aln_span: Option<u32>,
    ) -> Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)> {
        let shmmr_spec = seq_db.shmmr_spec.as_ref().unwrap();
        let shmmr_to_frags = &seq_db.agc_db.as_ref().unwrap().1;
        let res = aln::query_fragment_to_hps(
            shmmr_to_frags,
            &seq,
            shmmr_spec,
            penality,
            max_count,
            max_count_query,
            max_count_target,
            max_aln_span,
        );
        res
    }

    pub fn group_smps_by_principle_bundle_id(
        smps: &SmpsWithBundleLabel,
        length_cutoff: Option<u32>,
        merge_length: Option<u32>,
    ) -> Vec<SmpsWithBundleLabel> {
        let length_cutoff = if length_cutoff.is_some() {
            length_cutoff.unwrap()
        } else {
            2500
        };
        let merge_length = if merge_length.is_some() {
            merge_length.unwrap()
        } else {
            5000
        };

        let mut pre_bundle_id: Option<usize> = None;
        let mut pre_direction: Option<u8> = None;
        let mut all_partitions: Vec<SmpsWithBundleLabel> = vec![];
        let mut new_partition: Vec<SmpBundleTuple> = vec![];

        smps.iter().for_each(|(smp, bundle_info)| {
            if let Some(bundle_info) = bundle_info {
                let direction = if smp.4 == bundle_info.1 { 0_u8 } else { 1_u8 };
                let bundle_id = bundle_info.0;
                if pre_bundle_id.is_none() || pre_direction.is_none() {
                    new_partition.push((*smp, Some(*bundle_info)));
                    pre_bundle_id = Some(bundle_id);
                    pre_direction = Some(direction);

                    return;
                };
                if pre_bundle_id.unwrap() != bundle_id || pre_direction.unwrap() != direction {
                    if new_partition[new_partition.len() - 1].0 .3 - new_partition[0].0 .2
                        > length_cutoff
                    {
                        all_partitions.push(new_partition.clone());
                        new_partition.clear();
                    } else {
                        new_partition.clear();
                    }
                    pre_bundle_id = Some(bundle_id);
                    pre_direction = Some(direction);
                };
                new_partition.push((*smp, Some(*bundle_info)));
            }
        });

        if new_partition[new_partition.len() - 1].0 .3 - new_partition[0].0 .2 > length_cutoff {
            all_partitions.push(new_partition.clone());
        }

        let partition = &mut all_partitions[0].clone();
        let mut partitions: Vec<SmpsWithBundleLabel> = vec![];
        all_partitions[1..].iter().for_each(|p| {
            let p_end = partition[partition.len() - 1].0 .3;

            //can unwrap as unmatch segement is filtered out
            let (p_bundle_id, p_direction, _) = partition[partition.len() - 1].1.unwrap();
            let p_direction = if partition[partition.len() - 1].0 .4 == p_direction {
                0_u8
            } else {
                1_u8
            };

            let n_bgn = p[0].0 .2;
            let (n_budle_id, n_direction, _) = p[0].1.unwrap();
            let n_direction = if p[0].0 .4 == n_direction { 0_u8 } else { 1_u8 };

            if p_bundle_id == n_budle_id
                && p_direction == n_direction
                && i64::abs(n_bgn as i64 - p_end as i64) < merge_length as i64
            {
                partition.extend(p.clone());
            } else {
                partitions.push(partition.clone());
                partition.clear();
                partition.extend(p.clone());
            };
        });
        if partition.len() > 0 {
            partitions.push(partition.clone());
        }
        partitions
    }
}
