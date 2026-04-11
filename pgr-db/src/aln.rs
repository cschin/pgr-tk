// use rayon::prelude::*;
use crate::seq_db::{self, FragmentHit};
use crate::shmmrutils::{self, ShmmrSpec};
use log::debug;
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::collections::HashSet;
use wavefront_aln::*;

pub type HitPair = ((u32, u32, u8), (u32, u32, u8)); //(bgn1, end1, orientation1),  (bgn2, end2, orientation2)

pub fn sparse_aln(
    sp_hits: &mut Vec<HitPair>,
    max_span: u32,
    penalty: f32,
    max_gap: Option<u32>,
    orientated: bool,
) -> Vec<(f32, Vec<HitPair>)> {
    // given a set of hits in the form of (bgn1, end1, orientation1),  (bgn2, end2, orientation2)
    // perform (banded) dynamic programming to group them into list of hit chains
    sp_hits.sort_by(|a, b| a.0 .0.partial_cmp(&b.0 .0).unwrap());
    let mut v_s = FxHashMap::<HitPair, f32>::default(); // score for each vertex
    let mut best_pre_v = FxHashMap::<HitPair, Option<HitPair>>::default(); // look up for the best pre-vertex
    assert!(sp_hits.len() > 1);
    let first_hp = sp_hits[0];
    v_s.insert(first_hp, first_hp.0 .1 as f32 - first_hp.0 .0 as f32); // the score of the first node is just its length
    best_pre_v.insert(first_hp, None);

    (1..sp_hits.len()).for_each(|i| {
        let hp = sp_hits[i];
        let mut best_v = Option::<HitPair>::None;
        let mut best_s = 0_f32;
        let mut j = i;
        let mut span_set = HashSet::<(u32, u32, u8)>::new();
        loop {
            if j == 0 {
                break;
            };
            j -= 1;

            let pre_hp = sp_hits[j];

            if orientated {
                // don't connect if orientations are not agreed if orientated == true
                let p_orientation = pre_hp.0 .2 ^ pre_hp.1 .2; // pre_hp.0.2  = 0 or 1 and pre_hp.1.2 = 0 or 1
                let orientation = hp.0 .2 ^ hp.1 .2; // hp.0.2  = 0 or 1 and hp.1.2 = 0 or 1
                if p_orientation != orientation {
                    continue;
                }
            }

            if let Some(max_gap) = max_gap {
                let max_gap = max_gap as f32;
                if hp.0 .2 == hp.1 .2 {
                    if (hp.0 .0 as f32 - pre_hp.0 .1 as f32).abs() > max_gap
                        || (hp.1 .0 as f32 - pre_hp.1 .1 as f32).abs() > max_gap
                    {
                        continue;
                    }
                } else if (hp.0 .0 as f32 - pre_hp.0 .1 as f32).abs() > max_gap
                    || (hp.1 .1 as f32 - pre_hp.1 .0 as f32).abs() > max_gap
                {
                    continue;
                }
            }

            if pre_hp.0 == hp.0 {
                continue;
            }; // don't connect node with the same left coordinate
            span_set.insert(pre_hp.0);
            let p_s = v_s.get(&pre_hp).unwrap_or(&0_f32);
            let mut s: f32 = *p_s + (hp.0 .1 as f32 - hp.0 .0 as f32);

            if hp.0 .2 == hp.1 .2 {
                // same orientation
                s -= penalty
                    * ((hp.0 .0 as f32 - pre_hp.0 .1 as f32).abs()
                        + (hp.1 .0 as f32 - pre_hp.1 .1 as f32).abs());
            } else {
                // opposite orientation
                s -= penalty
                    * ((hp.0 .0 as f32 - pre_hp.0 .1 as f32).abs()
                        + (hp.1 .1 as f32 - pre_hp.1 .0 as f32).abs());
            }

            if s > best_s {
                best_s = s;
                best_v = Some(pre_hp);
            }

            if span_set.len() >= max_span as usize {
                break;
            };
        }

        if best_s > 0_f32 {
            v_s.insert(hp, best_s);
            best_pre_v.insert(hp, best_v);
        } else {
            v_s.insert(hp, hp.0 .1 as f32 - hp.0 .0 as f32);
            best_pre_v.insert(hp, None);
        }
    });

    let mut unvisited_v = FxHashSet::<HitPair>::default();
    unvisited_v.extend(sp_hits.iter());
    let mut out = Vec::<(f32, Vec<HitPair>)>::new();
    while !unvisited_v.is_empty() {
        let mut best_s = 0_f32; // global best score
        let mut best_v: Option<HitPair> = None; // global best vertex
                                                // println!("DBG un-visit len; {}", unvisited_v.len());
        unvisited_v.iter().for_each(|hp| {
            let s = v_s.get(hp).unwrap_or(&0_f32);
            if *s > best_s {
                best_s = *s;
                best_v = Some(*hp);
            }
        });
        let mut track = Vec::<HitPair>::new();
        let mut v = best_v;
        while v.is_some() {
            let hp = v.unwrap();
            if !unvisited_v.contains(&hp) {
                break;
            };
            track.push(hp);
            v = *best_pre_v.get(&hp).unwrap_or(&None);
        }
        if track.is_empty() {
            continue;
        };
        track.reverse();
        track.iter().for_each(|hp| {
            // let s = v_s.get(hp).unwrap_or(&0_f32);
            // println!("H {} {} {} {} {} {} {}", hp.0.0, hp.0.1, hp.0.2, hp.1.0, hp.1.1, hp.1.2, s );
            unvisited_v.remove(hp);
        });
        let bgn_s = v_s.get(&track[0]).unwrap_or(&0_f32);
        out.push((best_s - bgn_s, track));
    }
    out
}

pub type TargetHitPairLists = Vec<(u32, Vec<(f32, Vec<HitPair>)>)>; // target_id, Vec<(score, HitPairs)>

#[allow(clippy::too_many_arguments)]
pub fn query_fragment_to_hps(
    raw_query_hits: Vec<FragmentHit>,
    frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
    penalty: f32,
    max_count: Option<u32>,
    query_max_count: Option<u32>,
    target_max_count: Option<u32>,
    max_aln_span: Option<u32>,
    max_gap: Option<u32>,
    oriented: bool,
) -> TargetHitPairLists {
    let mut shmmr_pair_hash_count = FxHashMap::<(u64, u64), u32>::default();
    let mut query_shmmr_pair_hash_count = FxHashMap::<(u64, u64), u32>::default();
    let mut target_shmer_pair_count = FxHashMap::<(u64, u64, u32), u32>::default();

    seq_db::pair_shmmrs(&shmmrutils::sequence_to_shmmrs(0, frag, shmmr_spec, false))
        .iter()
        .for_each(|shmmr_pair| {
            let entry = query_shmmr_pair_hash_count
                .entry((shmmr_pair.0.hash(), shmmr_pair.1.hash()))
                .or_insert(0);
            *entry += 1;
        });

    raw_query_hits.iter().for_each(
        |(shmmr_pair_hash, _query_position, frag_signature): &(
            (u64, u64),
            _,
            Vec<seq_db::FragmentSignature>,
        )| {
            //let sp = d.0;
            // count shimmer pair hits
            let entry = shmmr_pair_hash_count.entry(*shmmr_pair_hash).or_insert(0);
            *entry += 1;

            frag_signature
                .iter()
                .for_each(|(_frg_id, seq_id, _bgn, _end, _orientation)| {
                    // count shimmer pair on target hits
                    // v = frg_id, seq_id, bgn, end, orientation(to shimmer pair)
                    let key = (shmmr_pair_hash.0, shmmr_pair_hash.1, *seq_id);
                    let entry = target_shmer_pair_count.entry(key).or_insert(0);
                    *entry += 1;
                })
        },
    );

    let mut target_squence_id_to_hits =
        FxHashMap::<u32, Vec<((u32, u32, u8), (u32, u32, u8))>>::default();
    raw_query_hits.into_iter().for_each(
        |(shmmr_pair, query_position, frag_signature): (
            (u64, u64),
            _,
            Vec<seq_db::FragmentSignature>,
        )| {
            let count = *shmmr_pair_hash_count.get(&shmmr_pair).unwrap_or(&0);
            let max_count = max_count.unwrap_or(128);
            if count > max_count {
                return;
            };
            let max_count_query = query_max_count.unwrap_or(128);
            if count > max_count_query {
                return;
            };
            let left_frag_coordinate = query_position;
            frag_signature
                .iter()
                .for_each(|&(_frg_id, sid, pos0, pos1, orientation)| {
                    let count = *target_shmer_pair_count
                        .get(&(shmmr_pair.0, shmmr_pair.1, sid))
                        .unwrap_or(&0);
                    let max_count_target = target_max_count.unwrap_or(128);
                    if count > max_count_target {
                        return;
                    };
                    let e = target_squence_id_to_hits.entry(sid).or_default();
                    let right_frag_coordinate = (pos0, pos1, orientation);
                    e.push((left_frag_coordinate, right_frag_coordinate));
                });
        },
    );

    let max_aln_span = max_aln_span.unwrap_or(8);

    target_squence_id_to_hits
        .into_iter()
        .filter(|(_sid, hps)| hps.len() > 1)
        .map(|(sid, mut hps)| {
            (
                sid,
                sparse_aln(&mut hps, max_aln_span, penalty, max_gap, oriented),
            )
        })
        .collect::<Vec<_>>()
}

pub fn wfa_align_bases(
    target_str: &str,
    query_str: &str,
    max_wf_length: u32,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
) -> Option<(String, String)> {
    let capacity = std::cmp::max(1024, std::cmp::max(target_str.len(), query_str.len()) >> 5);
    let mut wfs = WaveFronts::new_with_capacity(
        target_str,
        query_str,
        max_wf_length,
        mismatch_penalty,
        open_penalty,
        extension_penalty,
        capacity,
    );
    if wfs.step_all(Some(1024)) == WaveFrontStepResult::ReachEnd {
        Some(wfs.backtrace())
    } else {
        None
    }
}

pub fn aln_pair_map(aln_target_str: &str, aln_query_str: &str) -> Vec<(u32, u32, char)> {
    let paired = std::iter::zip(aln_target_str.as_bytes(), aln_query_str.as_bytes());
    let mut t_pos = 0_u32;
    let mut q_pos = 0_u32;
    paired
        .into_iter()
        .map(|(&tb, &qb)| {
            let mut t = '-';
            let new_t_pos = if tb == b'-' {
                t = 'I';
                t_pos
            } else {
                t_pos + 1
            };
            let new_q_pos = if qb == b'-' {
                t = 'D';
                q_pos
            } else {
                q_pos + 1
            };
            if tb == qb {
                t = 'M';
            };
            if tb != qb && tb != b'-' && qb != b'-' {
                t = 'X';
            };
            let out = (t_pos, q_pos, t);
            t_pos = new_t_pos;
            q_pos = new_q_pos;
            out
        })
        .collect::<Vec<_>>()
}

pub fn get_variants_from_aln_pair_map(
    aln_pairs: &[(u32, u32, char)],
    target_str: &str,
    query_str: &str,
) -> Vec<(u32, u32, char, String, String)> {
    let mut current_variant = Vec::<(char, char, char)>::new();

    let aggregate_variants = |previous_match: &(u32, u32, char, char, char),
                              current_variant: &Vec<(char, char, char)>|
     -> Option<(u32, u32, char, String, String)> {
        let t_variant_segment = String::from_iter(current_variant.iter().map(|v| v.0));
        let q_variant_segment = String::from_iter(current_variant.iter().map(|v| v.1));
        let t_variant_segment = t_variant_segment.replace('-', "").trim().to_string();
        let q_variant_segment = q_variant_segment.replace('-', "").trim().to_string();
        let t_len = t_variant_segment.len();
        let q_len = q_variant_segment.len();
        let v_type = match t_len.cmp(&q_len) {
            Ordering::Greater => 'D',
            Ordering::Less => 'I',
            Ordering::Equal => 'X',
        };

        match v_type {
            'X' => Some((
                previous_match.0 + 1,
                previous_match.1 + 1,
                'X',
                t_variant_segment,
                q_variant_segment,
            )),
            'D' => Some((
                previous_match.0,
                previous_match.1,
                'D',
                [previous_match.3.to_string(), t_variant_segment].join(""),
                [previous_match.4.to_string(), q_variant_segment].join(""),
            )),
            'I' => Some((
                previous_match.0,
                previous_match.1,
                'I',
                [previous_match.3.to_string(), t_variant_segment].join(""),
                [previous_match.4.to_string(), q_variant_segment].join(""),
            )),
            _ => None,
        }
    };
    let mut previous_match = (0_u32, 0_u32, 'U', '-', '-');
    let mut variants = Vec::<Option<(u32, u32, char, String, String)>>::new();

    aln_pairs.iter().for_each(|&(t_pos, q_pos, t)| match t {
        'M' => {
            let t_char = target_str.as_bytes()[t_pos as usize] as char;
            let q_char = query_str.as_bytes()[q_pos as usize] as char;
            if !current_variant.is_empty() {
                variants.push(aggregate_variants(&previous_match, &current_variant));
            };

            current_variant.clear();
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, t_char, q_char, t);
            previous_match = (t_pos, q_pos, 'M', t_char, q_char);
        }
        'X' => {
            let t_char = target_str.as_bytes()[t_pos as usize] as char;
            let q_char = query_str.as_bytes()[q_pos as usize] as char;
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, t_char, q_char, t);
            current_variant.push((t_char, q_char, t));
        }
        'I' => {
            let q_char = query_str.as_bytes()[q_pos as usize] as char;
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, '-', q_char, t);
            current_variant.push(('-', q_char, t));
        }
        'D' => {
            let t_char = target_str.as_bytes()[t_pos as usize] as char;
            debug!("{} {} {:1} {:1} {}", t_pos, q_pos, t_char, '-', t);
            current_variant.push((t_char, '-', t));
        }
        _ => {}
    });
    if !current_variant.is_empty() {
        variants.push(aggregate_variants(&previous_match, &current_variant));
    };
    variants.into_iter().flatten().collect::<Vec<_>>()
}

type AlignmentResult = Vec<(u32, u32, char, String, String)>;
pub fn get_wfa_variant_segments(
    target_str: &[u8],
    query_str: &[u8],
    left_padding: usize,
    max_wf_length: Option<u32>,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
) -> Option<AlignmentResult> {
    let set_len_diff = (query_str.len() as i64 - target_str.len() as i64).unsigned_abs() as u32;
    let max_wf_length = if let Some(max_wf_length) = max_wf_length {
        max_wf_length
    } else {
        std::cmp::max(2 * set_len_diff, 128_u32)
    };

    // we need to reverse the string for alignment such the the gaps are on the left
    // maybe we can do this in the stack for performance in the future
    // we assume the left_padding base on the left side are identical
    let mut r_t_str = target_str[left_padding..].to_vec();
    let mut r_q_str = query_str[left_padding..].to_vec();
    r_t_str.reverse();
    r_q_str.reverse();
    let r_t_str = String::from_utf8_lossy(&r_t_str[..]);
    let r_q_str = String::from_utf8_lossy(&r_q_str[..]);
    let t_len_minus_one = left_padding as u32 + r_t_str.len() as u32 - 1;
    let q_len_minus_one = left_padding as u32 + r_q_str.len() as u32 - 1;

    if let Some((aln_target_str, aln_query_str)) = wfa_align_bases(
        &r_t_str,
        &r_q_str,
        max_wf_length,
        mismatch_penalty,
        open_penalty,
        extension_penalty,
    ) {
        /*
        // print out the alignment string for debugging

        // let mut r_aln_target_str = aln_target_str.clone().as_bytes().to_owned();
        // let mut r_aln_query_str = aln_query_str.clone().as_bytes().to_owned();
        // r_aln_target_str.reverse();
        // r_aln_query_str.reverse();
        // let r_aln_target_str = String::from_utf8_lossy(&r_aln_target_str[..]);
        // let r_aln_query_str = String::from_utf8_lossy(&r_aln_query_str[..]);
        // println!("XX: {}", r_aln_target_str);
        // println!("XX: {}", r_aln_query_str);
        */

        let mut aln_pairs = aln_pair_map(&aln_target_str, &aln_query_str);
        // assume the base on the left are identical  ( # of base = left_padding)
        (0..left_padding).for_each(|delta| {
            aln_pairs.push((
                (r_t_str.len() + delta) as u32,
                (r_q_str.len() + delta) as u32,
                'M',
            ));
        });
        // convert the coordinate from the reverse to the forward sequence
        aln_pairs.iter_mut().for_each(|(t_pos, q_pos, _c)| {
            *t_pos = t_len_minus_one - *t_pos;
            *q_pos = q_len_minus_one - *q_pos;
        });
        aln_pairs.reverse();

        // compute the VCF like variant representation
        let target_str = String::from_utf8_lossy(target_str);
        let query_str = String::from_utf8_lossy(query_str);
        Some(get_variants_from_aln_pair_map(
            &aln_pairs,
            &target_str,
            &query_str,
        ))
    } else {
        None
    }
}

pub fn sw_align_bases(
    target_str: &str,
    query_str: &str,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
) -> Option<(String, String)> {
    let mut target_str = (*target_str).as_bytes().to_vec();
    let mut query_str = (*query_str).as_bytes().to_vec();
    target_str.reverse();
    query_str.reverse();
    let t_len = target_str.len();
    let q_len = query_str.len();

    // initial condition for j = 0
    let mut match_scores = (0..t_len + 1)
        .map(|i| {
            if i == 0 {
                0
            } else {
                -open_penalty - (i as i32) * extension_penalty
            }
        })
        .collect::<Vec<i32>>();

    let mut e_scores = (0..t_len + 1)
        .map(|i| {
            if i == 0 {
                i32::MIN
            } else {
                -open_penalty - (i as i32) * extension_penalty
            }
        })
        .collect::<Vec<i32>>();
    let mut f_scores = vec![i32::MIN; t_len + 1];

    let mut trace_back = (0..t_len + 1)
        .map(|_| vec![(0_i8, 0_i8); q_len + 1])
        .collect::<Vec<_>>();

    for row in trace_back.iter_mut().skip(1).take(t_len) {
        row[0] = (-1, 0);
    }

    for j in 1..q_len + 1 {
        // for i = 0
        let p_match_score = match_scores.clone();
        match_scores[0] = -open_penalty - (j as i32) * extension_penalty;
        e_scores[0] = i32::MIN;
        f_scores[0] = -open_penalty - (j as i32) * extension_penalty;
        trace_back[0][j] = (0, -1);

        for i in 1..t_len + 1 {
            let s = p_match_score[i - 1]
                - (if target_str[i - 1] == query_str[j - 1] {
                    0
                } else {
                    mismatch_penalty
                });

            let e = if e_scores[i - 1] == i32::MIN {
                i32::MIN
            } else {
                e_scores[i - 1] - extension_penalty
            };

            let f = if f_scores[i] == i32::MIN {
                i32::MIN
            } else {
                f_scores[i] - extension_penalty
            };

            (trace_back[i][j], match_scores[i]) = if s > e && s > f {
                ((-1, -1), s)
            } else if e > f {
                ((-1, 0), e)
            } else {
                ((0, -1), f)
            };

            let o = match_scores[i] - open_penalty;

            e_scores[i] = if o > e { o } else { e };

            f_scores[i] = if o > f { o } else { f }
        }
    }
    let mut t_pos = t_len;
    let mut q_pos = q_len;
    let mut aln_t = Vec::<u8>::new();
    let mut aln_q = Vec::<u8>::new();

    while t_pos != 0 || q_pos != 0 {
        let d = trace_back[t_pos][q_pos];
        if d.0 != 0 {
            t_pos -= 1;
            aln_t.push(target_str[t_pos]);
        } else {
            aln_t.push(b'-');
        };
        if d.1 != 0 {
            q_pos -= 1;
            aln_q.push(query_str[q_pos]);
        } else {
            aln_q.push(b'-');
        }
    }
    //aln_t.reverse();
    //aln_q.reverse();

    Some((
        String::from_utf8_lossy(&aln_t[..]).to_string(),
        String::from_utf8_lossy(&aln_q[..]).to_string(),
    ))
}

pub fn get_sw_variant_segments(
    target_str: &[u8],
    query_str: &[u8],
    left_padding: usize,
    mismatch_penalty: i32,
    open_penalty: i32,
    extension_penalty: i32,
) -> Option<AlignmentResult> {
    let t_str = target_str[left_padding..].to_vec();
    let q_str = query_str[left_padding..].to_vec();
    let t_str = String::from_utf8_lossy(&t_str[..]);
    let q_str = String::from_utf8_lossy(&q_str[..]);

    if let Some((aln_target_str, aln_query_str)) = sw_align_bases(
        &t_str,
        &q_str,
        mismatch_penalty,
        open_penalty,
        extension_penalty,
    ) {
        /*
        // print out the alignment string for debugging

        // let mut r_aln_target_str = aln_target_str.clone().as_bytes().to_owned();
        // let mut r_aln_query_str = aln_query_str.clone().as_bytes().to_owned();
        // r_aln_target_str.reverse();
        // r_aln_query_str.reverse();
        // let r_aln_target_str = String::from_utf8_lossy(&r_aln_target_str[..]);
        // let r_aln_query_str = String::from_utf8_lossy(&r_aln_query_str[..]);
        // println!("XX: {}", r_aln_target_str);
        // println!("XX: {}", r_aln_query_str);
        */

        let mut aln_pairs = Vec::<_>::new();
        // assume the base on the left are identical  ( # of base = left_padding)
        (0..left_padding).for_each(|delta| {
            aln_pairs.push((delta as u32, delta as u32, 'M'));
        });
        aln_pairs.extend(
            aln_pair_map(&aln_target_str, &aln_query_str)
                .into_iter()
                .map(|v| (v.0 + left_padding as u32, v.1 + left_padding as u32, v.2)),
        );

        // compute the VCF like variant representation
        let target_str = String::from_utf8_lossy(target_str);
        let query_str = String::from_utf8_lossy(query_str);
        Some(get_variants_from_aln_pair_map(
            &aln_pairs,
            &target_str,
            &query_str,
        ))
    } else {
        None
    }
}

#[cfg(test)]
mod test {

    #[test]

    fn sparse_aln_test() {
        use crate::aln::{sparse_aln, HitPair};
        use std::fs::File;
        use std::io::{BufRead, BufReader};
        let f = BufReader::new(File::open("./test/test_data/test_hits").unwrap());
        let mut hp = Vec::<HitPair>::new();
        f.lines().for_each(|s| {
            if let Ok(s) = s {
                let s = s.split_ascii_whitespace();
                let out = s
                    .into_iter()
                    .map(|s| s.parse::<u32>().unwrap())
                    .collect::<Vec<u32>>();
                assert_eq!(out.len(), 6);
                hp.push((
                    (out[0], out[1], out[2] as u8),
                    (out[3], out[4], out[5] as u8),
                ));
            }
        });
        let oriented = false;
        let max_gap = None;
        let out = sparse_aln(&mut hp, 8, 0.5_f32, max_gap, oriented);
        out.iter().for_each(|(s, v)| println!("{} {}", s, v.len()));
        // TODO: Test the output properly
    }

    #[test]
    fn test_wfa_align_bases() {
        use crate::aln::{aln_pair_map, get_variants_from_aln_pair_map, wfa_align_bases};
        use log::debug;
        //use simple_logger::SimpleLogger;
        //SimpleLogger::new().init().unwrap();
        let t_str = "ACATACATGTGTGTGAAAAATATATAAGTAAAAAAAATGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCATGAAACCCCATACATGAAAGTTGCATGAA";
        let q_str = "ACATACATGTGAAATATAATAAAAGTTGCATGAAAAAACATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCAAAAGTTGCATGAAACATACATGAAAATGAAAAAACATACATGAAAGTTGCATGAA";

        // Test that alignment works and produces a result
        let result = wfa_align_bases(t_str, q_str, 20, 2, 2, 1);
        assert!(result.is_some(), "WFA alignment should produce a result");

        if let Some((t_aln_str, q_aln_str)) = result {
            // Verify alignment strings have equal length
            assert_eq!(
                t_aln_str.len(),
                q_aln_str.len(),
                "Alignment strings should have equal length"
            );

            // Alignment strings should only contain valid characters
            assert!(
                t_aln_str
                    .chars()
                    .all(|c| c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == '-'),
                "Target alignment string should only contain ACGT or gaps"
            );
            assert!(
                q_aln_str
                    .chars()
                    .all(|c| c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == '-'),
                "Query alignment string should only contain ACGT or gaps"
            );

            // Generate alignment pairs
            let aln_pairs = aln_pair_map(&t_aln_str, &q_aln_str);
            assert!(!aln_pairs.is_empty(), "Alignment pairs should not be empty");

            // Get variants
            let variants = get_variants_from_aln_pair_map(&aln_pairs, t_str, q_str);

            // Verify variants make sense
            for (t_pos, q_pos, variant_type, t_segment, q_segment) in &variants {
                // Positions should be valid
                assert!(
                    *t_pos < t_str.len() as u32,
                    "Target position should be valid"
                );
                assert!(
                    *q_pos < q_str.len() as u32,
                    "Query position should be valid"
                );

                // Variant type should be one of the expected values
                assert!(
                    *variant_type == 'X' || *variant_type == 'I' || *variant_type == 'D',
                    "Variant type should be X, I, or D"
                );

                // For debug output
                println!(
                    "{} {} {} {} {}",
                    t_pos, q_pos, variant_type, t_segment, q_segment
                );
            }
        }
    }

    #[test]
    fn test_wfa_aggreate_variant() {
        use crate::aln::{
            aln_pair_map, get_variants_from_aln_pair_map, get_wfa_variant_segments, wfa_align_bases,
        };
        use log::debug;
        //use simple_logger::SimpleLogger;
        //SimpleLogger::new().init().unwrap();
        let t_str =
            "ACGGAGGTGAGCCTGGGAGCATAGAGGTGGGCCTGGGAGCATGGCGGCGGGGGGGGGGCCTGGGAGCACAGGGCGGGCC";
        let q_str =
            "ACGGAGGTGAGCCTGGGAGCATAGAGGTGGGCCTGGGAGCATGGCGGTGGGGGGGGGCCTGGGAGCACAGGGCGGGCC";

        if let Some(aln_res) =
            get_wfa_variant_segments(t_str.as_bytes(), q_str.as_bytes(), 1, Some(128), 3, 3, 1)
        {
            aln_res
                .into_iter()
                .for_each(|v| println!("{} {} {} {} {}", v.0, v.1, v.2, v.3, v.4));
        };
        // TODO: Test the output properly
    }

    #[test]
    fn test_sw_align_bases() {
        use crate::aln::{aln_pair_map, get_variants_from_aln_pair_map, sw_align_bases};
        use log::debug;
        //use simple_logger::SimpleLogger;
        //SimpleLogger::new().init().unwrap();
        let t_str = "ACATACATGTGTGTGAAAAATATATAAGTAAAAAAAATGCATGAAACCCCAAAAGTTGCATGAAACATACATGAAAATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCATGAAACCCCATACATGAAAGTTGCATGAA";
        let q_str = "ACATACATGTGAAATATAATAAAAGTTGCATGAAAAAACATACATGAAAGTTGCATGAAACATACATGAAAAAAGTTGCAAAAGTTGCATGAAACATACATGAAAATGAAAAAACATACATGAAAGTTGCATGAA";
        if let Some((t_aln_str, q_aln_str)) = sw_align_bases(t_str, q_str, 2, 2, 1) {
            println!("{}", t_aln_str);
            println!("{}", q_aln_str);
            let aln_pairs = aln_pair_map(&t_aln_str, &q_aln_str);
            let variants = get_variants_from_aln_pair_map(&aln_pairs, t_str, q_str);
            variants.into_iter().for_each(|(t_pos, q_pos, t, s1, s2)| {
                println!("{} {} {} {} {}", t_pos, q_pos, t, s1, s2);
            });
        }
        // TODO: Test the output properly
    }

    #[test]
    fn test_sw_aggreate_variant() {
        use crate::aln::{
            aln_pair_map, get_sw_variant_segments, get_variants_from_aln_pair_map, sw_align_bases,
        };
        use log::debug;
        //use simple_logger::SimpleLogger;
        //SimpleLogger::new().init().unwrap();
        let t_str =
            "ACGGAGGTGAGCCTGGGAGCATAGAGGTGGGCCTGGGAGCATGGCGGCGGGGGGGGGGCCTGGGAGCACAGGGCGGGCC";
        let q_str =
            "ACGGAGGTGAGCCTGGGAGCATAGAGGTGGGCCTGGGAGCATGGCGGTGGGGGGGGGCCTGGGAGCACAGGGCGGGCC";

        if let Some(aln_res) =
            get_sw_variant_segments(t_str.as_bytes(), q_str.as_bytes(), 1, 3, 3, 1)
        {
            aln_res
                .into_iter()
                .for_each(|v| println!("{} {} {} {} {}", v.0, v.1, v.2, v.3, v.4));
        };
        // TODO: Test the output properly
    }
}
