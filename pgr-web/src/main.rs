// main.rs

use dioxus::prelude::*;
use reqwest;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use web_sys::console;
use web_sys::Document;
//use pgr_db::aln::{self, HitPair};
type HitPair = ((u32, u32, u8), (u32, u32, u8)); //(bgn1, end1, orientation1),  (bgn2, end2, orientation2)


type SmpBundleTuple = ((u64, u64, u32, u32, u8), Option<(usize, u8, usize)>);
type SmpsWithBundleLabel = Vec<SmpBundleTuple>;


#[derive(Deserialize)]
struct TargetRanges {
    query_src_ctg: (String, String),
    matches: Vec<(u32, Vec<(f32, Vec<HitPair>)>)>,
    sid_ctg_src: Vec<(u32, String, String)>,
    principal_bundle_decomposition: Vec<(u32, u32, Vec<SmpsWithBundleLabel>)>,

}

#[derive(Deserialize)]
struct TargetRangesSimplified {
    query_src_ctg: (String, String),
    match_summary: Vec<(u32, Vec<(u32, u32, u32, u32, usize)>)>, // (q_bgn, q_end, t_bgn, t_end, num_hits)
    sid_ctg_src: Vec<(u32, String, String)>,
    principal_bundle_decomposition: Vec<(u32, String, Vec<(u32, u32, u32, u8)>)>, //bgn, end, bundle_id, bundle_direction
}

#[derive(Serialize)]
struct SequenceQuerySpec {
    source: String,
    ctg: String,
    bgn: usize,
    end: usize,
    padding: usize,
    merge_range_tol: usize,
    full_match: bool,
}

fn main() {
    dioxus::web::launch(app);
}

fn app(cx: Scope) -> Element {

    let rgn = use_future(&cx, (), |_| async move {
        let client = reqwest::Client::new();
        client
            .get("http://127.0.0.1:3000/")
            .send()
            .await
            .unwrap()
            .text()
            .await
    });

    cx.render( match  rgn.value() {
        Some(Ok(val)) =>  rsx!( 
            div { link { href:"https://fonts.googleapis.com/icon?family=Material+Icons", rel:"stylesheet", }
            }
            //div {"{val}"}
            crate::query_results {}    ),
        Some(Err(err)) => rsx!("errored!"),
        None => rsx!("loading!"),

    } )

}

pub fn query_results(cx: Scope) -> Element {
    let query = SequenceQuerySpec {
        source: "hg38_tagged".to_string(),
        ctg: "chr6_hg38".to_string(),
        //bgn: 160952514,
        //end: 161087407,
        bgn: 32163513,
        end: 32992088,
        padding: 200000,
        merge_range_tol: 2000000,
        full_match: true,
    };

    let targets = use_future(&cx, (), |_| async move {
        let client = reqwest::Client::new();
        client
            .post("http://127.0.0.1:3000/query_sdb")
            .json(&query)
            .send()
            .await
            .unwrap()
            .json::<TargetRangesSimplified>()
            .await
    });

    cx.render( match targets.value() {
        Some(Ok(val)) => {

            let sid_to_ctg_src = val.sid_ctg_src.iter().map(|v| {
                let (sid, ctg_name, src) = v;
                (*sid, (ctg_name, src))
            }).collect::<HashMap<u32,(&String, &String)>>();

            rsx!{
                div { class: "content-center",
                    br {}
                    br {}
                    div {class: "p-8 content-center",
                        table { class: "p-8 table-auto border-collapse border-spacing-4 border border-solid",
                            thead {
                                tr{
                                    th {class: "p-1", "sid"} 
                                    th {class: "p-1", "contig"}
                                    th {class: "p-1", "source"}
                                    th {class: "p-1", "hit count"}
                                    th {class: "p-1", "query span"}
                                    th {class: "p-1", "query len"}
                                    th {class: "p-1", "target span"}
                                    th {class: "p-1", "target len"}
                                    th {class: "p-1", "ADD"}
                                
                                }
                            }
                            tbody {
                                rsx!(val.match_summary.iter().map(|v| {
                                    let sid = v.0;
                                    let (ctg, src) = *sid_to_ctg_src.get(&sid).unwrap();
                                    let style_classes = "p-1 border border-slate-900 text-center";
                                    let hit_summary = v.1.iter().map(move |(q_bgn, q_end, t_bgn, t_end, n_hits)| {

                                        let q_span = format!("{}-{}", q_bgn, q_end);
                                        let t_span = format!("{}-{}", t_bgn, t_end);
                                        let q_len = q_end - q_bgn;
                                        let t_len = if t_end > t_bgn {t_end - t_bgn} else { t_bgn - t_end};
                                        //let t_span: i32 = w.1[l-1].1.1 as i32 - w.1[0].1.0 as i32;
                                        rsx!( tr { class: "border-solid text-center" ,
                                            td { class: "{style_classes}", "{sid}"}  
                                            td { class: "{style_classes}", "{ctg}"} 
                                            td { class: "{style_classes}", "{src}"}
                                            td { class: "{style_classes}", "{n_hits}"} 
                                            td { class: "{style_classes}", "{q_span}"} 
                                            td { class: "{style_classes}", "{q_len}"} 
                                            td { class: "{style_classes}", "{t_span}"}
                                            td { class: "{style_classes}", "{t_len}"}
                                            td { class: "{style_classes}", button { class: "border border-slate-900", "ADD"}}
                                            } )
                                    });
                                        
                                rsx!( hit_summary)
                                }))
                            }
                        }
                    }
                    rsx!( 
                        br {} 
                        br {}   
                        val.principal_bundle_decomposition.iter().flat_map(|(sid, ctg_name, r)| {track(cx, (*sid, r.clone()))}
                    )) 
                }
            }
        },
        Some(Err(err)) => rsx!( "Err" ),
        None => rsx!("loading")
    })
}

pub fn track(cx: Scope, range:  (u32, Vec<(u32, u32, u32, u8)>) ) -> Element {
    console::log_1(&"Rendering the track".into());
    cx.render(
        rsx! {
            div { 
                class: "p-1",
                svg {
                    width: "1250",
                    height: "25",
                    view_box: "-10000 -20 1750000 40",
                    defs {
                        marker {
                            id: "arrow",
                            markerWidth: "10",
                            markerHeight: "10",
                            refX: "0",
                            refY: "3",
                            orient: "auto",
                            markerUnits: "strokWidth",
                            view_box: "0 0 20 20",
                            path {
                                d: "M0,0 L0,6 L4,3 z",
                                fill: "#888"
                            }
                        }
                    } 
                    
                    range.1.iter().map(|(bgn, end, bundle_id, direction)| {
                        let sid = range.0;
                        let mut bgn = bgn;
                        let mut end = end;
                        let mut y = "0";
                        let mut style = "stroke:rgb(255,0,0);stroke-width:5000";
                        if *direction == 1 {
                            (bgn ,end) = (end, bgn);
                            y = "-2400";
                            style = "stroke:rgb(0,0,255);stroke-width:5000";
                        }
                        let line_id = format!("s_{}_{}_{}_{}", sid, bundle_id, bgn, end);
                        let line_id2 = line_id.clone();
                        rsx! {
                            line {
                                id: "{line_id}",
                                onclick: move |evt| {
                                    let window = web_sys::window().expect("global window does not exists");    
	                                let document = window.document().expect("expecting a document on window");
                                    let el = document.get_element_by_id(&line_id2).unwrap();
                                    let style = el.attributes().get_named_item("style").unwrap();
                                    console::log_1(&line_id2.clone().into());
                                    console::log_1(&style.value().into());
                                    style.set_value("stroke:rgb(0,0,255);stroke-width:15000");
                                },
                                x1: "{bgn}",
                                y1: "{y}",
                                x2: "{end}",
                                y2: "{y}", 
                                style: "{style}",
                                marker_end: "url(#arrow)"
                            }
                        }
                    })
                } 
            }
            
        }
    )
}