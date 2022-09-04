use axum::{
    extract::{Extension, Path, Query},
    http::{
        header::{self, HeaderMap, HeaderName},
        HeaderValue, Method,
    },
    response::{Html, IntoResponse},
    routing::{get, post},
    Json, Router,
};
use pgr_db::{
    aln::{self, HitPair},
    fasta_io::reverse_complement,
};
use pgr_server::seq_index_db::*;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::net::SocketAddr;
use std::sync::Arc;
use tower::ServiceBuilder;
use tower_http::cors::Any;
use tower_http::cors::CorsLayer;
use tower_http::trace::TraceLayer;
use tracing::Span;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

#[derive(Deserialize)]

struct SequenceQuerySpec {
    source: String,
    ctg: String,
    bgn: usize,
    end: usize,
    padding: usize,
    merge_range_tol: usize,
    full_match: bool,
}

#[derive(Serialize)]
struct TargetRanges {
    query_src_ctg: (String, String),
    matches: Vec<(u32, Vec<(f32, Vec<aln::HitPair>)>)>,
    sid_ctg_src: Vec<(u32, String, String)>,
    principal_bundle_decomposition: Vec<(u32, Vec<SmpsWithBundleLabel>)>,
}

#[derive(Serialize)]
struct TargetRangesSimplified {
    query_src_ctg: (String, String),
    match_summary: Vec<(u32, Vec<(u32, u32, u32, u32, usize, bool)>)>, // (q_bgn, q_end, t_bgn, t_end, num_hits)
    sid_ctg_src: Vec<(u32, String, String)>,
    principal_bundle_decomposition: Vec<(u32, String, Vec<(u32, u32, u32, u8)>)>, //bgn, end, bundle_id, bundle_direction
}

#[tokio::main]
async fn main() {
    tracing_subscriber::registry()
        .with(tracing_subscriber::EnvFilter::new(
            std::env::var("RUST_LOG")
                .unwrap_or_else(|_| "example_tracing_aka_logging=debug,tower_http=debug".into()),
        ))
        .with(tracing_subscriber::fmt::layer())
        .init();

    let mut seq_db = SeqIndexDB::new();
    let _ = seq_db.load_from_agc_index(
        //"/wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-small_panel".to_string(),
        "/wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-v0".to_string(),
    );
    let seq_db = Arc::new(seq_db);
    // build our application with a route
    let app = Router::new()
        .route(
            "/",
            get({
                let seq_db = seq_db.clone();
                move || handler(seq_db)
            }),
        )
        .route(
            "/query_sdb",
            post({
                let seq_db = seq_db.clone();
                move |params| query_sdb_with(params, seq_db)
            }),
        )
        .layer(
            CorsLayer::new()
                .allow_origin(Any)
                //.allow_origin("http://127.0.0.1:8080".parse::<HeaderValue>().unwrap())
                .allow_methods(Any)
                .allow_headers(Any),
        )
        .layer(ServiceBuilder::new().layer(TraceLayer::new_for_http()));

    // run it
    let addr = SocketAddr::from(([127, 0, 0, 1], 3000));
    println!("listening on {}", addr);
    axum::Server::bind(&addr)
        .serve(app.into_make_service())
        .await
        .unwrap();
}

/*
async fn handler(seq_db: Arc<SeqIndexDB>) -> impl IntoResponse {
    let n_ctg = 0;
    let mut headers = HeaderMap::new();
    headers.insert(header::CONTENT_TYPE, "text/plain".parse().unwrap());
    headers.insert(header::URI, "http://127.0.0.1:3000".parse().unwrap());
    let rtn = format!("Hello, World! {}", n_ctg);
    (headers, rtn)
}
*/

async fn handler(seq_db: Arc<SeqIndexDB>) -> Json<usize> {
    let n_ctg = seq_db.seq_index.as_ref().unwrap().len();
    Json(n_ctg)
}

async fn query_sdb_with(
    Json(payload): Json<SequenceQuerySpec>,
    seq_db: Arc<SeqIndexDB>,
) -> Json<TargetRangesSimplified> {
    let agc_db = seq_db.agc_db.as_ref().unwrap();
    let sample_name = payload.source;
    let ctg_name = payload.ctg;
    let padding = payload.padding;
    let merge_range_tol = payload.merge_range_tol;
    let seq_len = match seq_db
        .seq_index
        .as_ref()
        .unwrap()
        .get(&(ctg_name.clone(), Some(sample_name.clone())))
    {
        None => 0,
        Some(value) => value.1,
    };

    let q_seq_len = payload.end - payload.bgn;
    let q_seq_bgn = if padding > payload.bgn {
        0
    } else {
        payload.bgn - padding
    };
    let q_seq_end = if payload.end + padding > seq_len as usize {
        seq_len as usize
    } else {
        payload.end + padding
    };

    let sub_seq =
        (&agc_db.0).get_sub_seq(sample_name.clone(), ctg_name.clone(), q_seq_bgn, q_seq_end);

    /*
    println!(
        "DBG: sub_seq_len {:?} {} {}",
        sub_seq.len(),
        q_seq_bgn,
        q_seq_end
    );
     */

    let matches = query_fragment_to_hps(
        &seq_db,
        sub_seq.clone(),
        0.25,
        Some(128),
        Some(128),
        Some(128),
        Some(0),
    );

    let mut sid_target_regions: Vec<_> = matches
        .iter()
        .map(|(sid, ms)| {
            let mut targegt_regions = ms
                .iter()
                .filter(|(_, m)| m.len() >= 4)
                .map(|(_, m)| {
                    let mut f_count = 0_u32;
                    let mut r_count = 0_u32;
                    let mut rgns: Vec<(u32, u32, u32, u32)> = vec![];
                    m.iter().for_each(|v| {
                        if v.0 .2 == v.1 .2 {
                            f_count += 1;
                        } else {
                            r_count += 1;
                        };
                        rgns.push((v.1 .0, v.1 .1, v.0 .0, v.0 .1));
                    });
                    rgns.sort();

                    let t_bgn = rgns[0].0;
                    let q_bgn = rgns[0].2;
                    let t_end = rgns[rgns.len() - 1].1;
                    let q_end = rgns[rgns.len() - 1].3;

                    if f_count > r_count {
                        (t_bgn, t_end, q_bgn, q_end, 0_u8, m)
                    } else {
                        (t_bgn, t_end, q_bgn, q_end, 1_u8, m)
                    }
                })
                .collect::<Vec<_>>();

            targegt_regions.sort();

            type Matches = Vec<((u32, u32, u8), (u32, u32, u8))>;

            let mut merged_regions: Vec<Vec<(u32, u32, u32, u32, u8, &Matches)>> = vec![];

            if targegt_regions.len() > 0 {
                //println!("DBG: targegt_regions count: {}", targegt_regions.len());

                let fwd_regions = targegt_regions
                    .iter()
                    .filter(|&r| r.4 == 0)
                    .collect::<Vec<_>>();
                //println!("DBG: fwd_regions count: {}:{}", sid, fwd_regions.len());
                let rev_regions = targegt_regions
                    .iter()
                    .filter(|&r| r.4 == 1)
                    .collect::<Vec<_>>();
                //println!("DBG: rev_regions count: {}:{}", sid, rev_regions.len());
                fwd_regions.into_iter().for_each(|v| {
                    if merged_regions.len() == 0 {
                        merged_regions.push(vec![v.clone()]);
                        return;
                    } else {
                        let last_idx = merged_regions.len() - 1;
                        let last_m_rgn = &mut merged_regions[last_idx];
                        let last_idx = last_m_rgn.len() - 1;
                        let last_rgn = last_m_rgn[last_idx];
                        //println!("mfDBG {} {} : {} {}", last_rgn.0, last_rgn.1, v.0, v.1);
                        if i64::abs((v.0 as i64) - (last_rgn.1 as i64)) < (merge_range_tol as i64) {
                            last_m_rgn.push(v.clone());
                        } else {
                            merged_regions.push(vec![v.clone()]);
                        }
                    }
                });
                rev_regions.into_iter().for_each(|v| {
                    if merged_regions.len() == 0 {
                        merged_regions.push(vec![v.clone()]);
                        return;
                    } else {
                        let last_idx = merged_regions.len() - 1;
                        let last_m_rgn = &mut merged_regions[last_idx];
                        let last_idx = last_m_rgn.len() - 1;
                        let last_rgn = last_m_rgn[last_idx];
                        //println!("mrDBG {} {} : {} {}", last_rgn.0, last_rgn.1, v.0, v.1);
                        if i64::abs((v.0 as i64) - (last_rgn.1 as i64)) < (merge_range_tol as i64) {
                            last_m_rgn.push(v.clone());
                        } else {
                            merged_regions.push(vec![v.clone()]);
                        }
                    }
                });
            }
            /*
            println!(
                "DBG: merged_regions count: {}:{}",
                sid,
                merged_regions.len()
            );
            */
            merged_regions.sort();
            (*sid, merged_regions)
        })
        .collect();

    let mut sid_ctg_src = sid_target_regions
        .iter()
        .map(|&(sid, _)| {
            let r = seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
            match &r.1 {
                Some(src) => (sid, r.0.clone(), src.clone()),
                None => (sid, r.0.clone(), "none".to_string()),
            }
        })
        .collect::<Vec<(u32, String, String)>>();
    sid_ctg_src.sort();
    sid_target_regions.sort_by_key(|v| v.0);

    let match_summary: Vec<(u32, Vec<(u32, u32, u32, u32, usize, bool)>)> = sid_target_regions
        .iter()
        .map(|(sid, h)| {
            let summary = h
                .iter()
                .map(|m| {
                    let n_hits = m.iter().map(|v| v.5.len()).sum();

                    let mut q_list = m.iter().map(|v| (v.2, v.3)).collect::<Vec<(u32, u32)>>();
                    q_list.sort();

                    let t_min_bgn = m[0].0;
                    let t_max_end = m[m.len() - 1].1;
                    let reversed = if m[0].2 > m[m.len() - 1].3 {
                        true
                    } else {
                        false
                    };
                    let q_min_bgn = q_list[0].0;
                    let q_max_end = q_list[q_list.len() - 1].1;

                    (q_min_bgn, q_max_end, t_min_bgn, t_max_end, n_hits, reversed)
                })
                .filter(|v| {
                    let (q_bgn, q_end) = if (v.0 < v.1) { (v.0, v.1) } else { (v.1, v.0) };
                    (q_bgn as usize) <= (padding as usize)
                        && (q_end as usize) >= q_seq_len + (padding as usize)
                        && ((v.3 - v.2) as f32) > ((q_seq_len + 2 * padding) as f32) * 0.5
                })
                .collect::<Vec<(u32, u32, u32, u32, usize, bool)>>();

            (*sid, summary)
        })
        .filter(|v| v.1.len() > 0)
        .collect();

    let seq_list = match_summary
        .iter()
        .flat_map(|v| {
            let sid = v.0;
            //let (ctg_name, sample_name, _) = seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
            //let sample_name = sample_name.as_ref().unwrap();
            // println!("DBG0: {}", ctg_name);
            v.1.iter()
                .map(|h| {
                    let t_bgn = h.2;
                    let t_end = h.3;
                    let reversed = h.5;
                    //if t_bgn > t_end {
                    //    (t_bgn, t_end) = (t_end, t_bgn);
                    //    reversed = true;
                    //}
                    let (ctg_name, sample_name, _) =
                        seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
                    let sample_name = sample_name.as_ref().unwrap();
                    //println!("DBG: {}", ctg_name);
                    let mut seq = (&agc_db.0).get_sub_seq(
                        sample_name.clone(),
                        ctg_name.clone(),
                        t_bgn as usize,
                        t_end as usize,
                    );
                    if reversed {
                        seq = reverse_complement(&seq);
                    }
                    (format!("{}_{}_{}", ctg_name, t_bgn, t_end), seq)
                })
                .collect::<Vec<(String, Vec<u8>)>>()
        })
        .collect::<Vec<(String, Vec<u8>)>>();

    let mut new_sdb = SeqIndexDB::new();
    new_sdb.load_from_seq_list(seq_list.clone(), Some(&"Memory".to_string()), 56, 56, 4, 28);

    let (_principal_bundles, seqid_smps_with_bundle_id_seg_direction) =
        new_sdb.get_principal_bundle_decomposition(0, 8);

    let principal_bundle_decomposition = seqid_smps_with_bundle_id_seg_direction
        .iter()
        .map(|(sid, smps_with_bundle_info)| {
            (
                *sid,
                group_smps_by_principle_bundle_id(smps_with_bundle_info, None, None),
            )
        })
        .collect::<Vec<(u32, Vec<SmpsWithBundleLabel>)>>();

    let mut principal_bundle_decomposition: Vec<(u32, String, Vec<(u32, u32, u32, u8)>)> =
        principal_bundle_decomposition
            .into_iter()
            .map(|(sid, bundles)| {
                let summary = bundles
                    .into_iter()
                    .map(|b| {
                        let bgn = b[0].0 .2;
                        let end = b[b.len() - 1].0 .3;
                        let bundle_id = b[0].1.unwrap().0;
                        let direction = if b[0].0 .4 == b[0].1.unwrap().1 {
                            0_u8
                        } else {
                            1_u8
                        };
                        (bgn, end, bundle_id as u32, direction)
                    })
                    .collect::<Vec<(u32, u32, u32, u8)>>();
                let ctg_name = new_sdb
                    .seq_info
                    .as_ref()
                    .unwrap()
                    .get(&sid)
                    .unwrap()
                    .0
                    .clone();
                (sid, ctg_name, summary)
            })
            .collect();

    principal_bundle_decomposition.sort();

    Json(TargetRangesSimplified {
        query_src_ctg: (sample_name, ctg_name),
        match_summary,
        sid_ctg_src,
        principal_bundle_decomposition,
    })
}
