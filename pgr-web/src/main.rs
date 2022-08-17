// main.rs

use dioxus::prelude::*;
use reqwest;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
//use pgr_db::aln::{self, HitPair};
type HitPair = ((u32, u32, u8), (u32, u32, u8)); //(bgn1, end1, orientation1),  (bgn2, end2, orientation2)

#[derive(Deserialize, Serialize)]
struct Test {
    number: u32,
}

#[derive(Deserialize)]
struct TargetRanges {
    query_src_ctg: (String, String),
    matches: Vec<(u32, Vec<(f32, Vec<HitPair>)>)>,
    sid_ctg_src: Vec<(u32, String, String)>
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
        bgn: 160952514,
        end: 160952514+300000,
        padding: 10000,
        merge_range_tol: 100000,
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
            .json::<TargetRanges>()
            .await
    });

    cx.render( match targets.value() {
        Some(Ok(val)) => {
            let sid_to_ctg_src = val.sid_ctg_src.iter().map(|v| {
                let (sid, ctg_name, src) = v;
                (*sid, (ctg_name, src))
             
            }).collect::<HashMap<u32,(&String, &String)>>();
            rsx!{
                br {}
                br {}
                div {class: "p-8",
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
                        rsx!(val.matches.iter().map(|v| {
                            let sid = v.0;
                            let (ctg, src) = *sid_to_ctg_src.get(&sid).unwrap();
                            let style_classes = "p-1 border border-slate-900 text-center";
                            let hit_summary = v.1.iter().map(move |w| {
                                let l = w.1.len();
                                let qbgn = w.1[0].0.0;
                                let qend = w.1[l-1].0.1;
                                let tbgn = w.1[0].1.0;
                                let tend =  w.1[l-1].1.1;

                                let q_span = format!("{}-{}", qbgn, qend);
                                let t_span = format!("{}-{}", tbgn, tend);
                                let q_len = qend - qbgn;
                                let t_len = if tend > tbgn {tend - tbgn} else { tbgn - tend};
                                //let t_span: i32 = w.1[l-1].1.1 as i32 - w.1[0].1.0 as i32;
                                rsx!( tr { class: "border-solid text-center" ,
                                    td { class: "{style_classes}", "{sid}"}  
                                    td { class: "{style_classes}", "{ctg}"} 
                                    td { class: "{style_classes}", "{src}"}
                                    td { class: "{style_classes}", "{l}"} 
                                    td { class: "{style_classes}", "{q_span}"} 
                                    td { class: "{style_classes}", "{q_len}"} 
                                    td { class: "{style_classes}", "{t_span}"}
                                    td { class: "{style_classes}", "{t_len}"}
                                    td { class: "{style_classes}", button { class: "border border-slate-900", "ADD"}}
                                    } )
                            });         
                        rsx!(
                            hit_summary
                            )
                        }))
                        
                    }}
                }

            }
        },
        Some(Err(err)) => rsx!( "Err" ),
        None => rsx!("loading")
    })
}
