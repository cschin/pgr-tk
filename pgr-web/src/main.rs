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

static cmap: [&str;97] = ["#870098","#00aaa5","#3bff00","#ec0000","#00a2c3","#00f400","#ff1500","#0092dd",
                          "#00dc00","#ff8100","#007ddd","#00c700","#ffb100","#0038dd","#00af00","#fcd200",
                          "#0000d5","#009a00","#f1e700","#0000b1","#00a55d","#d4f700","#4300a2","#00aa93",
                          "#a1ff00","#dc0000","#00aaab","#1dff00","#f40000","#009fcb","#00ef00","#ff2d00",
                          "#008ddd","#00d700","#ff9900","#0078dd","#00c200","#ffb900","#0025dd","#00aa00",
                          "#f9d700","#0000c9","#009b13","#efed00","#0300aa","#00a773","#ccf900","#63009e",
                          "#00aa98","#84ff00","#e10000","#00a7b3","#00ff00","#f90000","#009bd7","#00ea00",
                          "#ff4500","#0088dd","#00d200","#ffa100","#005ddd","#00bc00","#ffc100","#0013dd",
                          "#00a400","#f7dd00","#0000c1","#009f33","#e8f000","#1800a7","#00aa88","#c4fc00",
                          "#78009b","#00aaa0","#67ff00","#e60000","#00a4bb","#00fa00","#fe0000","#0098dd",
                          "#00e200","#ff5d00","#0082dd","#00cc00","#ffa900","#004bdd","#00b400","#ffc900",
                          "#0000dd","#009f00","#f4e200","#0000b9","#00a248","#dcf400","#2d00a4","#00aa8d",
                          "#bcff00"];

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
                        let line_color = cmap[(bundle_id % 97) as usize];
                        let line_with = 5000;
                        let style = format!("stroke:{};stroke-width:{}", line_color, line_with);
                        
                        if *direction == 1 {
                            (bgn ,end) = (end, bgn);
                            y = "-2400";
                        }
                        let line_id = format!("s_{}_{}_{}_{}", sid, bundle_id, bgn, end);
                        let line_class = format!("bdl_{}", bundle_id);
                        let line_class2 = line_class.clone();
                        rsx! {
                            line {
                                id: "{line_id}",
                                class: "{line_class} normal",
                                onclick: move |evt| {
                                    let window = web_sys::window().expect("global window does not exists");    
	                                let document = window.document().expect("expecting a document on window");
                                    let line_elements = document.get_elements_by_class_name(&line_class2);
                                    console::log_1(&line_class2.clone().into());
                                    (0..line_elements.length()).into_iter().for_each(|idx| {
                                        let el = line_elements.item(idx).unwrap(); 
                                        let classes = el.class_list();
                                        let mut line_with = line_with; 
                                        if classes.contains(&"normal") {
                                            line_with *= 2; 
                                            classes.remove_1(&"normal");
                                            classes.add_1(&"highlited");
                                        } else {
                                            classes.add_1(&"normal");
                                            classes.remove_1(&"highlited");
                                        };

                                        let style = el.attributes().get_named_item("style").unwrap();
                                        console::log_1(&style.value().into());
                                        let style_str = format!("stroke:{};stroke-width:{}", line_color, line_with);
                                        style.set_value(&style_str[..]);
                                    });
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