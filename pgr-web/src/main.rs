// main.rs

use dioxus::prelude::*;
use reqwest;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use web_sys::console;
use web_sys::Document;
use rustc_hash::FxHashMap;
use wasm_bindgen::JsCast;
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

#[derive(Deserialize, Clone)]
pub struct TargetRangesSimplified {
    query_src_ctg: (String, String),
    match_summary: Vec<(u32, Vec<(u32, u32, u32, u32, usize, bool)>)>, // (q_bgn, q_end, t_bgn, t_end, num_hits)
    sid_ctg_src: Vec<(u32, String, String)>,
    principal_bundle_decomposition: Vec<(u32, String, Vec<(u32, u32, u32, u8)>)>, //bgn, end, bundle_id, bundle_direction
}


#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ShmmrSpec {
    pub w: u32,
    pub k: u32,
    pub r: u32,
    pub min_span: u32,
    pub sketch: bool,
}

#[derive(Serialize, Deserialize, Clone, PartialEq)]
pub struct SequenceQuerySpec {
    source: String,
    ctg: String,
    bgn: usize,
    end: usize,
    padding: usize,
    merge_range_tol: usize,
    full_match: bool,
    pb_shmmr_spec: ShmmrSpec
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
    let ROI_json = include_str!("data/ROIs.json");
    let rois: FxHashMap<String, SequenceQuerySpec> = serde_json::from_str(ROI_json).unwrap();
    let rois2 = rois.clone();

    let query = use_state(&cx, || <Option<SequenceQuerySpec>>::None);
    let query_name = use_state(&cx, || <Option<String>>::None);
   
    //let q = query.current().as_ref().clone();

    let t = use_future(&cx, (query_name,), |(query_name)| async move {
        console::log_1(&"query".into());    
        let client = reqwest::Client::new();
        let qn = query_name.0.current().as_ref().clone();
        let q = if qn.is_none() {None} else {
            Some(rois2.get(&qn.unwrap()).unwrap())
        };
        client
            .post("http://127.0.0.1:3000/query_sdb")
            .json(&q)
            .send()
            .await
            .unwrap()
            .json::<TargetRangesSimplified>()
            .await
    });



    cx.render(
        rsx! {
            div { class: "p-8", 
                div {"PanGenome Research Tool Kit Principal Bundle Demo"}
           
                select { 
                    name: "ROI_selector",
                    id: "ROI_selector",
                    class: "px-6 py-2.5",
                    rois.iter().map(|(k, v)| {
                        rsx! { 
                            option {
                                class: "px-6 py-2.5",
                                value: "{k}",
                                selected: "true" ,
                                "{k}"
                            }
                        }
                    })
                }
                button { 
                    class: "inline-block px-6 py-2.5 bg-blue-600 text-white rounded",
                    onclick: move |_| {
                        console::log_1(&"clicked".into()); 
                        let window = web_sys::window().expect("global window does not exists");    
                        let document = window.document().expect("expecting a document on window");
                        
                        let roi_selector: web_sys::HtmlSelectElement = document.get_element_by_id(&"ROI_selector").unwrap().dyn_into().unwrap();
                        let options =  roi_selector.options();
                        console::log_1(&options.selected_index().unwrap().into());
                        let selected_value = options.get_with_index(options.selected_index().unwrap() as u32).unwrap().get_attribute("value").unwrap();
                        console::log_1(&selected_value.clone().into());
                        let new_query =rois.get(&selected_value).unwrap().clone(); 
                        query.modify(move |_| Some(new_query));
                        query_name.modify(move |_| Some(selected_value.clone()));
                    },
                    "Show" 
                }
            }
            div {[
                
                match t.value() {
                    Some(Ok(val)) => {
                        //target.modify(|_| Some(val.clone())); 
                        console::log_1(&"target modified".into());       
                        //cx.needs_update();
                        rsx! { div {[query_results2(cx, query.current().as_ref().clone(), 
                                                       Some(val.clone()))]} }
                    },
                    Some(Err(err)) => rsx! {div {"Err"}},
                    None => rsx! {div {"Loading"}},
                }
            ]}
            
                //query_results2(cx, query.current().as_ref().clone(), 
                //                     target.current().as_ref().clone())]}
        }
    )
}


pub fn query_results2(cx: Scope, query: Option<SequenceQuerySpec>, target: Option<TargetRangesSimplified>) -> Element {
    console::log_1(&"rendering query_results2".into()); 

    if target.is_none() {
        let r = rsx! { div {} };
        return cx.render(r)
    } 
    if query.is_none() {
        let r = rsx! { div {} };
        return cx.render(r)
    } 
    
    console::log_1(&"rendering query_results2, 2".into()); 
    let val = target.unwrap();

    console::log_1(&query.clone().unwrap().ctg.into()); 
    console::log_1(&val.match_summary.len().into()); 
    let sid_to_ctg_src = val.sid_ctg_src.iter().map(|v| {
        let (sid, ctg_name, src) = v;
        (*sid, (ctg_name, src))
    }).collect::<HashMap<u32,(&String, &String)>>();
    let q = query.unwrap().clone();
    let query = q.clone();
    let ctg = query.ctg;            
    let bgn = query.bgn;            
    let end = query.end;            
    let mut track_size = (query.end - query.bgn + 2 * query.padding);  
    track_size = track_size + (track_size >> 1);
    console::log_1(&"rendering query_results2, 3".into()); 
    cx.render (
    rsx!{
        div { class: "grid p-8  grid-cols-1 justify-center space-y-4",
            h2 {"Query: {ctg}:{bgn}-{end}"}
            div { class: "overflow-x-auto sm:-mx-6 lg:-mx-8",
                p {class: "px-8 py-2", "Query results"}
                div {class: "flex flex-col max-h-[250px]",
                div {class: "flex-grow overflow-auto",
                    table { class: "relative w-full",
                        thead {
                            tr{
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "sid"} 
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "contig"}
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "source"}
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "hit count"}
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "query span"}
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "query len"}
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "target span"}
                                th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "target len"}
                            }
                        }
                        tbody {
                            class: "divide-y",
                            rsx!(val.match_summary.iter().map(|v| {
                                let sid = v.0;
                                let (ctg, src) = *sid_to_ctg_src.get(&sid).unwrap();
                                let style_classes = "px-1 py-2 text-center";
                                let hit_summary = v.1.iter().map(move |(q_bgn, q_end, t_bgn, t_end, n_hits, reversed)| {

                                    let q_span = format!("{}-{}", q_bgn, q_end);
                                    let t_span = format!("{}-{}", t_bgn, t_end);
                                    let q_len = if q_end > q_bgn { q_end - q_bgn } else { q_bgn - q_end };
                                    let t_len = if t_end > t_bgn {t_end - t_bgn} else { t_bgn - t_end};
                                    rsx!( tr {
                                        td { class: "{style_classes}", "{sid}"}  
                                        td { class: "{style_classes}", "{ctg}"} 
                                        td { class: "{style_classes}", "{src}"}
                                        td { class: "{style_classes}", "{n_hits}"} 
                                        td { class: "{style_classes}", "{q_span}"} 
                                        td { class: "{style_classes}", "{q_len}"} 
                                        td { class: "{style_classes}", "{t_span}"}
                                        td { class: "{style_classes}", "{t_len}"}
                                        } )
                                });
                                    
                            rsx!( hit_summary)
                            }))
                        }
                    }
                }
                }
                rsx!( 
                    hr {class: "my-2 h-px bg-gray-700 border-0 dark:bg-gray-700"}
                    h2 {class: "px-8 py-2", "Principal Bundle Decomposition"}
                    div {
                        class: "px-8 content-center overflow-auto min-w-[1280px] max-h-[550px]",
                        val.principal_bundle_decomposition.iter().flat_map(|(sid, ctg_name, r)| {
                            track(cx, ctg_name.clone(), track_size, (*sid, r.clone()))
                        })
                    }) 
                }
            }
        }
    )
}



pub fn query_results(cx: Scope, query: Option<SequenceQuerySpec>) -> Element {
    if query.is_none() {
        let r = rsx! { div {} };
        return cx.render(r)
    } 

    let q = query.unwrap().clone();
    let query = q.clone();

    let targets = use_future(&cx, (), |_| async move {
        let client = reqwest::Client::new();
        client
            .post("http://127.0.0.1:3000/query_sdb")
            .json(&q)
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
            let ctg = query.ctg;            
            let bgn = query.bgn;            
            let end = query.end;            
            let mut track_size = (query.end - query.bgn + 2 * query.padding);  
            track_size = track_size + (track_size >> 1);
            rsx!{
                div { class: "grid p-8  grid-cols-1 justify-center space-y-4",
                    h2 {"Query: {ctg}:{bgn}-{end}"}
                    div { class: "overflow-x-auto sm:-mx-6 lg:-mx-8",
                        p {class: "px-8 py-2", "Query results"}
                        div {class: "flex flex-col max-h-[250px]",
                        div {class: "flex-grow overflow-auto",
                            table { class: "relative w-full",
                                thead {
                                    tr{
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "sid"} 
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "contig"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "source"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "hit count"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "query span"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "query len"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "target span"}
                                        th {class: "px-1 py-2 sticky top-0 text-blue-900 bg-blue-300", "target len"}
                                    }
                                }
                                tbody {
                                    class: "divide-y",
                                    rsx!(val.match_summary.iter().map(|v| {
                                        let sid = v.0;
                                        let (ctg, src) = *sid_to_ctg_src.get(&sid).unwrap();
                                        let style_classes = "px-1 py-2 text-center";
                                        let hit_summary = v.1.iter().map(move |(q_bgn, q_end, t_bgn, t_end, n_hits, reversed)| {

                                            let q_span = format!("{}-{}", q_bgn, q_end);
                                            let t_span = format!("{}-{}", t_bgn, t_end);
                                            let q_len = if q_end > q_bgn { q_end - q_bgn } else { q_bgn - q_end };
                                            let t_len = if t_end > t_bgn {t_end - t_bgn} else { t_bgn - t_end};
                                            rsx!( tr {
                                                td { class: "{style_classes}", "{sid}"}  
                                                td { class: "{style_classes}", "{ctg}"} 
                                                td { class: "{style_classes}", "{src}"}
                                                td { class: "{style_classes}", "{n_hits}"} 
                                                td { class: "{style_classes}", "{q_span}"} 
                                                td { class: "{style_classes}", "{q_len}"} 
                                                td { class: "{style_classes}", "{t_span}"}
                                                td { class: "{style_classes}", "{t_len}"}
                                                } )
                                        });
                                            
                                    rsx!( hit_summary)
                                    }))
                                }
                            }
                        }
                        }
                        rsx!( 
                            hr {class: "my-2 h-px bg-gray-700 border-0 dark:bg-gray-700"}
                            h2 {class: "px-8 py-2", "Principal Bundle Decomposition"}
                            div {
                                class: "px-8 content-center overflow-auto min-w-[1280px] max-h-[550px]",
                                val.principal_bundle_decomposition.iter().flat_map(|(sid, ctg_name, r)| {
                                    track(cx, ctg_name.clone(), track_size, (*sid, r.clone()))
                                })
                            }) 
                        }
                    }
                }
            },
        Some(Err(err)) => rsx!( "Err" ),
        None => rsx!("loading")
    })
}

pub fn track(cx: Scope, ctg_name: String, track_size: usize, range:  (u32, Vec<(u32, u32, u32, u8)>) ) -> Element {
    console::log_1(&"Rendering the track".into());
    let left_padding = track_size >> 4;
    cx.render(
        rsx! {
            div { 
                class: "p-1",
                p { "{ctg_name}"}
                svg {
                    width: "1250",
                    height: "50",
                    view_box: "-{left_padding} -80 {track_size} 120",
                    preserveAspectRatio: "none",
                    
                    range.1.iter().map(|(bgn, end, bundle_id, direction)| {
                        let sid = range.0;
                        let mut bgn = bgn;
                        let mut end = end;
                        if *direction == 1 {
                            (bgn, end) = (end, bgn);
                        }

                        let bundle_color = cmap[(bundle_id % 97) as usize];
                        let arror_end = *end as f32;
                        let end = *end as f32 - (*end as f32 - *bgn as f32) * 0.20;

                        let line_id = format!("s_{}_{}_{}_{}", sid, bundle_id, bgn, end);
                        let line_class = format!("bdl_{}", bundle_id);
                        let line_class2 = line_class.clone();
                        let path_str;
                        if *direction == 1 {
                            path_str = format!("M {bgn} -32 L {bgn} -16 L {end} -16 L {end} -12 L {arror_end} -24 L {end} -36 L {end} -32 Z");         
                        } else {
                            path_str = format!("M {bgn} -8 L {bgn} 8 L {end} 8 L {end} 12 L {arror_end} 0 L {end} -12 L {end} -8 Z");  
                        }
                        rsx! {
                            g {
                                id: "{line_id}",
                                class: "{line_class} normal",
                                onclick: move |_evt| {
                                    let window = web_sys::window().expect("global window does not exists");    
	                                let document = window.document().expect("expecting a document on window");
                                    let line_elements = document.get_elements_by_class_name(&line_class2);
                                    console::log_1(&line_class2.clone().into());
                                    (0..line_elements.length()).into_iter().for_each(|idx| {
                                        let el = line_elements.item(idx).unwrap(); 
                                        let classes = el.class_list();
                                        let transform_str;
                                        if classes.contains(&"normal") {
                                            let _ = classes.remove_1(&"normal");
                                            let _ = classes.add_1(&"highlited");
                                            transform_str = format!("scale(1,2)");
                                        } else {
                                            let _ = classes.add_1(&"normal");
                                            let _ = classes.remove_1(&"highlited");
                                            transform_str = format!("scale(1,1)");
                                        };

                                        let stroke = el.attributes().get_named_item("transform").unwrap();
                                        console::log_1(&stroke.value().into());
                                        stroke.set_value(&transform_str[..]);
                                    });
                                },
                                transform: "scale(1,1)",
                          
                                path {
                                    d: "{path_str}", 
                                    fill: "{bundle_color}",
                                    //stroke: "black",
                                    //stroke_width: "2",
                                    fill_opacity: "0.8",
                                }
                            }
                        }
                    })
                } 
            }
            
        }
    )
}