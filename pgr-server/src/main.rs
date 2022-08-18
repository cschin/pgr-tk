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
use pgr_db::aln::{self, HitPair};
use pgr_server::seq_index_db::*;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::sync::Arc;
use std::{collections::HashMap, net::SocketAddr};
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
    match_summary: Vec<(u32, Vec<(u32, u32, u32, u32, usize)>)>, // (q_bgn, q_end, t_bgn, t_end, num_hits)
    sid_ctg_src: Vec<(u32, String, String)>,
    principal_bundle_decomposition: Vec<(String, Vec<(u32, u32, u32, u8)>)>, //bgn, end, bundle_id, bundle_direction
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
        "/wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-small_panel".to_string(),
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
    let seq_len = match seq_db
        .seq_index
        .as_ref()
        .unwrap()
        .get(&(ctg_name.clone(), Some(sample_name.clone())))
    {
        None => 0,
        Some(value) => value.1,
    };
    let bgn = if padding > payload.bgn {
        0
    } else {
        payload.bgn - padding
    };
    let end = if payload.end + padding > seq_len as usize {
        seq_len as usize
    } else {
        payload.end + padding
    };

    let sub_seq = (&agc_db.0).get_sub_seq(sample_name.clone(), ctg_name.clone(), bgn, end);
    let mut matches = query_fragment_to_hps(
        &seq_db,
        sub_seq.clone(),
        0.25,
        Some(128),
        Some(128),
        Some(128),
        Some(0),
    );

    matches = matches
        .into_iter()
        .map(|v| {
            let hits =
                v.1.into_iter()
                    .filter(|h| {
                        let q_bgn = h.1[0].0 .0;
                        let q_end = h.1[h.1.len() - 1].0 .1;
                        //println!("DBG: {} {} {} {}", q_bgn, q_end, padding, padding + (payload.end - payload.bgn) as usize);
                        (q_bgn as usize) < padding
                            && q_end as usize > padding + (payload.end - payload.bgn) as usize
                    })
                    .collect::<Vec<(f32, Vec<((u32, u32, u8), (u32, u32, u8))>)>>();
            (v.0, hits)
        })
        .collect::<Vec<(u32, Vec<(f32, Vec<((u32, u32, u8), (u32, u32, u8))>)>)>>();

    let mut sid_ctg_src = matches
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
    //println!("{:?}", src_ctg_sid);
    matches.sort_by_key(|v| v.0);

    let match_summary: Vec<(u32, Vec<(u32, u32, u32, u32, usize)>)> = matches
        .iter()
        .map(|(sid, h)| {
            let summary = h
                .iter()
                .map(|m| {
                    let n_hits = m.1.len();
                    let q_bgn = m.1[0].0 .0;
                    let q_end = m.1[n_hits - 1].0 .1;
                    let t_bgn = m.1[0].1 .0;
                    let t_end = m.1[n_hits - 1].1 .1;
                    (q_bgn, q_end, t_bgn, t_end, n_hits)
                })
                .collect::<Vec<(u32, u32, u32, u32, usize)>>();
            (*sid, summary)
        })
        .collect();

    let seq_list = matches
        .iter()
        .flat_map(|v| {
            let sid = v.0;
            v.1.iter()
                .map(|h| {
                    let mut t_bgn = h.1[0].1 .0;
                    let mut t_end = h.1[h.1.len() - 1].1 .1;
                    if t_bgn > t_end {
                        (t_bgn, t_end) = (t_end, t_bgn);
                    }
                    let (ctg_name, sample_name, _) =
                        seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
                    let sample_name = sample_name.as_ref().unwrap();
                    let seq = (&agc_db.0).get_sub_seq(
                        sample_name.clone(),
                        ctg_name.clone(),
                        t_bgn as usize,
                        t_end as usize,
                    );
                    (format!("{}_{}_{}", ctg_name, t_bgn, t_end), seq)
                })
                .collect::<Vec<(String, Vec<u8>)>>()
        })
        .collect::<Vec<(String, Vec<u8>)>>();

    let mut new_sdb = SeqIndexDB::new();
    new_sdb.load_from_seq_list(
        seq_list.clone(),
        Some(&"Memory".to_string()),
        56,
        56,
        1,
        26,
    );

    let (principal_bundles, seqid_smps_with_bundle_id_seg_direction) =
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

    let principal_bundle_decomposition: Vec<(String, Vec<(u32, u32, u32, u8)>)> =
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
                let ctg_name = new_sdb.seq_info.as_ref().unwrap().get(&sid).unwrap().0.clone();
                (ctg_name, summary)
            })
            .collect();

    Json(TargetRangesSimplified {
        query_src_ctg: (sample_name, ctg_name),
        match_summary,
        sid_ctg_src,
        principal_bundle_decomposition,
    })
}
