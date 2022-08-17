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
    principal_bundle_decomposition: Vec<(u32, Vec<( (u64,u64,u32,u32,u8), Option<(usize, u8, usize)> )> )> 
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
) -> Json<TargetRanges> {
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
    let mut target_ranges = query_fragment_to_hps(
        &seq_db,
        sub_seq.clone(),
        0.25,
        Some(128),
        Some(128),
        Some(128),
        Some(0),
    );

    target_ranges = target_ranges
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

    let mut src_ctg_sid = target_ranges
        .iter()
        .map(|&(sid, _)| {
            let r = seq_db.seq_info.as_ref().unwrap().get(&sid).unwrap();
            match &r.1 {
                Some(src) => (sid, r.0.clone(), src.clone()),
                None => (sid, r.0.clone(), "none".to_string()),
            }
        })
        .collect::<Vec<(u32, String, String)>>();
    src_ctg_sid.sort();
    //println!("{:?}", src_ctg_sid);
    target_ranges.sort_by_key(|v| v.0);

    let seq_list = target_ranges
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
        128,
        56,
        4,
        26,
    );

    let (principal_bundles, seqid_smps_with_bundle_id_seg_direction) =
        new_sdb.get_principal_bundle_decomposition(0, 8);

    Json(TargetRanges {
        query_src_ctg: (sample_name, ctg_name),
        matches: target_ranges,
        sid_ctg_src: src_ctg_sid,
        principal_bundle_decomposition: seqid_smps_with_bundle_id_seg_direction
    })
}
