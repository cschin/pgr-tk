// main.rs

use dioxus::prelude::*;
use reqwest;
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashMap;
use web_sys::console;
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


#[derive(Clone)]
struct QueryState(String);


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
    let query_state = use_state(&cx, || "done".to_string()); 
   
    //let q = query.current().as_ref().clone();

    let targets = use_future(&cx, (query_name,), |(query_name)| async move {
        console::log_1(&"query".into());        
        let window = web_sys::window().expect("global window does not exists");    
        let document = window.document().expect("expecting a document on window");
        let query_result_div = document.get_element_by_id(&"query_results").unwrap();
        let _ = query_result_div.set_attribute("hidden", "true");

   
        let query_status_div = document.get_element_by_id(&"query_status").unwrap();
        let _ = query_status_div.remove_attribute("hidden");

        let query_button_div = document.get_element_by_id(&"query_button").unwrap();
        let _ = query_button_div.set_attribute("disabled", "true");


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
            .json::<Option<TargetRangesSimplified>>()
            .await
    });

    let mut kvs = rois.iter().map(|(k,v)| { (k.clone(), v.clone())} ).collect::<Vec<_>>();
    kvs.sort_by_key(|v| v.0.clone());
    cx.render(
        rsx! {
            div { class: "flex flex-row p-4", 
                div { class: "basis-2/4",
                    h2 {"PanGenome Research Tool Kit: Principal Bundle Decomposition Demo"}
                }
                div { class: "basis-1/4 mb-3 xl:w-96",
                 
                    select { 
                        name: "ROI_selector",
                        id: "ROI_selector",
                        class: "form-select appearance-none  w-full px-3 py-1.5 focus:text-gray-700 focus:bg-white focus:border-blue-600 focus:outline-none",
                        
                        kvs.iter().map(|(k, v)| {
                            rsx! { 
                                option {
                                    value: "{k}",
                                    "{k}"
                                }
                            }
                        })
                    }
                }

                div { class: "basis-1/4",
                    button { 
                        id: "query_button",
                        disabled: "false",
                        class: "inline-block px-6 py-1.5 bg-blue-600 text-white rounded",
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
            }
            div { id: "query_status",
                  class: "p-4",
                  "waiting"
            }

            div { id: "query_results",
                  hidden: "true",
                [
                match targets.value() {
                    Some(Ok(val)) => {
                        let window = web_sys::window().expect("global window does not exists");    
                        let document = window.document().expect("expecting a document on window");
                        let query_result_div = document.get_element_by_id(&"query_results").unwrap();
                        let _ = query_result_div.remove_attribute("hidden");
                        
                        let query_status_div = document.get_element_by_id(&"query_status").unwrap();
                        let _ = query_status_div.set_attribute("hidden", "true");

                        let query_button_div = document.get_element_by_id(&"query_button").unwrap();
                        let _ = query_button_div.remove_attribute("disabled");
                        
                        console::log_1(&"target modified xxx".into());  
                        
                        // reset all bundle class
                        let bundle_elements = document.get_elements_by_class_name(&"bundle");
                        (0..bundle_elements.length()).into_iter().for_each(|idx| {
                            let el = bundle_elements.item(idx).unwrap(); 
                            let classes = el.class_list();
                            let _ = classes.add_1(&"normal");
                            let _ = classes.remove_1(&"highlited");
                            let path = el.children().item(0).unwrap();
                            let stroke = path.attributes().get_named_item("stroke-width").unwrap();
                            console::log_1(&stroke.value().into());
                            stroke.set_value(&format!("0.5"));
                        });

                        rsx! { div {[query_results(cx, query.clone(), query_state.clone(), 
                                                       val.clone())]} }
                    },
                    Some(Err(err)) => rsx! {div {class: "p-4", "Err"}},
                    None => rsx! {div { class: "p-4", "No target yet"}},
                }
            ]}
            
        }
    )
}


pub fn query_results(cx: Scope, 
                     query: UseState<Option<SequenceQuerySpec>>, 
                     query_state: UseState<String>, 
                     target: Option<TargetRangesSimplified>) -> Element {

    console::log_1(&"rendering query_results2".into()); 
    
    //let query_state = query_state.current().as_ref().clone();
    //console::log_1(&query_state.clone().into()); 
    //if query_state == "requesting".to_string() {
    //    let r = rsx! { div { class: "p-4", "Requesting data" } };
    //    return cx.render(r)
    //} 


    let query = query.current().as_ref().clone();

    if query.is_none() {
        let r = rsx! { div { class: "p-4", "No query yet" } };
        return cx.render(r)
    } 

    if target.is_none() {
        let r = rsx! { div { class: "p-4", "Query sent, waiting for targets" } };
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
    
    let query = query.unwrap().clone();
    let ctg = query.ctg;            
    let bgn = query.bgn;            
    let end = query.end;            
    let mut track_size = (query.end - query.bgn + 2 * query.padding);  
    track_size = track_size + (track_size >> 1);
    console::log_1(&"rendering query_results2, 3".into()); 
    cx.render (
    rsx!{
        div { class: "grid p-2  grid-cols-1 justify-center space-y-2",
            div { class: "overflow-x-auto sm:-mx-6 lg:-mx-8",
                div {class: "flex flex-col min-w-[1280px]  max-h-screen",
                    rsx!( 
                        h2 {class: "px-8 py-2", "Principal Bundle Decomposition, Query: {ctg}:{bgn}-{end}"}
                        div {
                            class: "px-8 content-center overflow-auto min-w-[1280px] max-h-[450px]",
                            val.principal_bundle_decomposition.iter().flat_map(|(sid, ctg_name, r)| {
                                track(cx, ctg_name.clone(), track_size, (*sid, r.clone()))
                            })
                        }) 
                    }
                    hr {class: "my-2 h-px bg-gray-700 border-0 dark:bg-gray-700"}
                    div {class: "px-8 py-1",
                        div {class: "flex-grow overflow-auto max-h-[250px]",
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
                }
            }
        }
    )
}

pub fn track(cx: Scope, ctg_name: String, track_range: usize, range:  (u32, Vec<(u32, u32, u32, u8)>) ) -> Element {
    console::log_1(&"Rendering the track".into());
    let track_length = 1600;
    let left_padding = track_range >> 8;
    let scaling_factor = track_length as f32 / (track_range + 2*left_padding) as f32; 
    let left_padding = left_padding as f32 * scaling_factor as f32;
    let stroke_width = 0.5;
    let ctg_id = format!("ctg_{}", ctg_name);
    cx.render(
        rsx! {
            div { 
                class: "p-1",
                p { "{ctg_name}"}
                svg {
                    id: "{ctg_id}",
                    width: "{track_length}",
                    height: "40",
                    view_box: "0 -16 {track_length} 24",
                    preserveAspectRatio: "none",
                    
                    range.1.iter().map(|(bgn, end, bundle_id, direction)| {
                        let sid = range.0;
                        let mut bgn = *bgn as f32 * scaling_factor + left_padding;
                        let mut end = *end as f32 * scaling_factor + left_padding;
                        if *direction == 1 {
                            (bgn, end) = (end, bgn);
                        }

                        let bundle_color = cmap[(bundle_id % 97) as usize];
                        let stroke_color = cmap[((bundle_id * 47) % 43) as usize];
                        let arror_end = end as f32;
                        let end = if *direction == 0 {
                            if end as f32 - 5.0 < bgn { bgn } else { end as f32 - 5.0 }
                        } else {
                            if end as f32 + 5.0 > bgn { bgn } else { end as f32 + 5.0 } 
                        };  

                        let line_id = format!("s_{}_{}_{}_{}", sid, bundle_id, bgn, end);
                        let line_class = format!("bundle-{}", bundle_id);
                        let line_class2 = line_class.clone();
                        let path_str = format!("M {bgn} -3 L {bgn} 3 L {end} 3 L {end} 4 L {arror_end} 0 L {end} -4 L {end} -3 Z");  
                        rsx! {
                            g {
                                id: "{line_id}",
                                class: "{line_class} bundle normal",
                                onclick: move |_evt| {
                                    let window = web_sys::window().expect("global window does not exists");    
	                                let document = window.document().expect("expecting a document on window");
                                    let line_elements = document.get_elements_by_class_name(&line_class2);
                                    console::log_1(&line_class2.clone().into());
                                    (0..line_elements.length()).into_iter().for_each(|idx| {
                                        let el = line_elements.item(idx).unwrap(); 
                                        let classes = el.class_list();
                                        let stroke_width_str;
                                        if classes.contains(&"normal") {
                                            let _ = classes.remove_1(&"normal");
                                            let _ = classes.add_1(&"highlited");
                                            stroke_width_str =  format!("1.5");
                                        } else {
                                            let _ = classes.add_1(&"normal");
                                            let _ = classes.remove_1(&"highlited");
                                            stroke_width_str =  format!("0.5");
                                        };

                                        let path = el.children().item(0).unwrap();
                                        let stroke = path.attributes().get_named_item("stroke-width").unwrap();
                                        console::log_1(&stroke.value().into());
                                        stroke.set_value(&stroke_width_str[..]);
                                    });
                                },                                
                          
                                path {
                                    d: "{path_str}", 
                                    fill: "{bundle_color}",
                                    stroke: "{stroke_color}",
                                    stroke_width: "{stroke_width}",
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