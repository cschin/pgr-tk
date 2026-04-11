use clap::{self, Parser};
use rusqlite::Connection;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Deserialize;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{self, Path};
use svg::node::{element, Node};
use svg::Document;

#[allow(dead_code)] // need the standard names for deserialization if they are not use
#[derive(Deserialize, Clone, Debug)]
struct CtgMapRec {
    t_name: String,
    ts: u32,
    te: u32,
    q_name: String,
    qs: u32,
    qe: u32,
    ctg_len: u32,
    orientation: u32,
    ctg_orientation: u32,
    t_dup: bool,
    t_ovlp: bool,
    q_dup: bool,
    q_ovlp: bool,
}

#[derive(Deserialize)]
struct CtgMapSet {
    records: Vec<CtgMapRec>,
    target_length: Vec<(u32, String, u32)>,
    query_length: Vec<(u32, String, u32)>,
}

type CytoRecord = (u32, u32, String, String);
#[derive(Deserialize)]
struct CytoBands {
    cytobands: FxHashMap<String, Vec<CytoRecord>>,
}

/// generate align block plot from ctgmap.json file
#[derive(Parser, Debug)]
#[clap(about, long_about = None)]

pub struct Args {
    /// path to a ctgmap.json file
    #[clap(long, short)]
    pub ctgmap_json_path: String,

    /// the prefix of the output files
    #[clap(long, short)]
    pub output_prefix: String,

    /// if given, we will use this to determine the plot scale, this is useful for generate many plot in the same scale
    #[clap(long)]
    pub total_target_bases: Option<f64>,

    /// set the panel width
    #[clap(long, default_value_t = 1400.0)]
    pub panel_width: f64,

    /// draw the reference track with cytoband
    #[clap(long)]
    pub cytoband_json: Option<String>,

    /// if given, we will only generate plot for the specified contig in the reference
    #[clap(long)]
    pub ctg: Option<String>,

    /// if given, it will highlight regions specified by the bed file in the reference(target) track
    #[clap(long)]
    pub ref_annotation_bed: Option<String>,

    /// generate SVG instead of HTML
    #[clap(long)]
    pub svg: bool,
}

static CMAP: [&str; 97] = [
    "#870098", "#00aaa5", "#3bff00", "#ec0000", "#00a2c3", "#00f400", "#ff1500", "#0092dd",
    "#00dc00", "#ff8100", "#007ddd", "#00c700", "#ffb100", "#0038dd", "#00af00", "#fcd200",
    "#0000d5", "#009a00", "#f1e700", "#0000b1", "#00a55d", "#d4f700", "#4300a2", "#00aa93",
    "#a1ff00", "#dc0000", "#00aaab", "#1dff00", "#f40000", "#009fcb", "#00ef00", "#ff2d00",
    "#008ddd", "#00d700", "#ff9900", "#0078dd", "#00c200", "#ffb900", "#0025dd", "#00aa00",
    "#f9d700", "#0000c9", "#009b13", "#efed00", "#0300aa", "#00a773", "#ccf900", "#63009e",
    "#00aa98", "#84ff00", "#e10000", "#00a7b3", "#00ff00", "#f90000", "#009bd7", "#00ea00",
    "#ff4500", "#0088dd", "#00d200", "#ffa100", "#005ddd", "#00bc00", "#ffc100", "#0013dd",
    "#00a400", "#f7dd00", "#0000c1", "#009f33", "#e8f000", "#1800a7", "#00aa88", "#c4fc00",
    "#78009b", "#00aaa0", "#67ff00", "#e60000", "#00a4bb", "#00fa00", "#fe0000", "#0098dd",
    "#00e200", "#ff5d00", "#0082dd", "#00cc00", "#ffa900", "#004bdd", "#00b400", "#ffc900",
    "#0000dd", "#009f00", "#f4e200", "#0000b9", "#00a248", "#dcf400", "#2d00a4", "#00aa8d",
    "#bcff00",
];

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

fn load_ctgmap_from_db(db_path: &str) -> CtgMapSet {
    let conn = Connection::open_with_flags(db_path, rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY)
        .expect("can't open alndb file");

    let mut records = Vec::<CtgMapRec>::new();
    {
        let mut stmt = conn
            .prepare(
                "SELECT target_name, target_start, target_end,
                        query_name, query_start, query_end,
                        ctg_len, orientation, ctg_orientation,
                        target_dup, target_ovlp, query_dup, query_ovlp
                 FROM ctgmap",
            )
            .expect("prepare ctgmap query");
        let mut rows = stmt.query([]).expect("query ctgmap");
        while let Some(row) = rows.next().expect("ctgmap row") {
            records.push(CtgMapRec {
                t_name: row.get(0).unwrap(),
                ts: row.get(1).unwrap(),
                te: row.get(2).unwrap(),
                q_name: row.get(3).unwrap(),
                qs: row.get(4).unwrap(),
                qe: row.get(5).unwrap(),
                ctg_len: row.get(6).unwrap(),
                orientation: row.get(7).unwrap(),
                ctg_orientation: row.get(8).unwrap(),
                t_dup: row.get::<_, i32>(9).unwrap() != 0,
                t_ovlp: row.get::<_, i32>(10).unwrap() != 0,
                q_dup: row.get::<_, i32>(11).unwrap() != 0,
                q_ovlp: row.get::<_, i32>(12).unwrap() != 0,
            });
        }
    }

    let mut target_length = Vec::<(u32, String, u32)>::new();
    let mut query_length = Vec::<(u32, String, u32)>::new();
    {
        let mut stmt = conn
            .prepare("SELECT seq_id, seq_name, length FROM sequences")
            .expect("prepare sequences query");
        let mut rows = stmt.query([]).expect("query sequences");
        while let Some(row) = rows.next().expect("sequences row") {
            let seq_id: u32 = row.get(0).unwrap();
            let seq_name: String = row.get(1).unwrap();
            let length: u32 = row.get(2).unwrap();
            // query seqs have seq_id >= 1_000_000 (offset applied at write time)
            if seq_id >= 1_000_000 {
                query_length.push((seq_id - 1_000_000, seq_name, length));
            } else {
                target_length.push((seq_id, seq_name, length));
            }
        }
    }

    CtgMapSet {
        records,
        target_length,
        query_length,
    }
}

pub fn run(args: Args) -> Result<(), std::io::Error> {
    let mut ctgmap_set: CtgMapSet = if args.ctgmap_json_path.ends_with(".alndb") {
        load_ctgmap_from_db(&args.ctgmap_json_path)
    } else {
        let mut ctgmap_json_file = BufReader::new(
            File::open(Path::new(&args.ctgmap_json_path)).expect("can't open the input file"),
        );
        let mut buffer = Vec::new();
        ctgmap_json_file.read_to_end(&mut buffer)?;
        serde_json::from_str(&String::from_utf8_lossy(&buffer[..]))
            .expect("can't parse the ctgmap.json file")
    };

    let cytobands = if let Some(cytoband_path) = args.cytoband_json.clone() {
        let mut cytoband_file = BufReader::new(
            File::open(Path::new(&cytoband_path)).expect("can't open the cytoband json file"),
        );
        let mut buffer = Vec::new();
        cytoband_file.read_to_end(&mut buffer)?;
        let cytobands: CytoBands = serde_json::from_str(&String::from_utf8_lossy(&buffer[..]))
            .expect("can't parse the cytoband json file");
        Some(cytobands)
    } else {
        None
    };

    let ref_highlight = if let Some(ref_annotation_bed) = args.ref_annotation_bed.clone() {
        let bed_file_path = path::Path::new(&ref_annotation_bed);
        let bed_file = BufReader::new(File::open(bed_file_path).expect("can't open the bed file"));
        let mut ref_highlight = FxHashMap::<String, Vec<(u32, u32)>>::default();
        let bed_file_parse_err_msg = "bed file parsing error";
        bed_file.lines().for_each(|line| {
            let line = line.unwrap().trim().to_string();
            if line.is_empty() {
                return;
            }
            if &line[0..1] == "#" {
                return;
            }
            let bed_fields = line.split('\t').collect::<Vec<&str>>();
            let ctg: String = bed_fields[0].to_string();
            let bgn: u32 = bed_fields[1].parse().expect(bed_file_parse_err_msg);
            let end: u32 = bed_fields[2].parse().expect(bed_file_parse_err_msg);
            let e = ref_highlight.entry(ctg).or_insert_with(Vec::new);
            e.push((bgn, end));
        });
        Some(ref_highlight)
    } else {
        None
    };

    ctgmap_set.query_length.sort();
    ctgmap_set.target_length.sort();
    let mut ctg_target_hit_len = FxHashMap::<String, FxHashMap<String, u32>>::default();

    let query_length = ctgmap_set
        .query_length
        .iter()
        .map(|v| (v.1.clone(), v.2))
        .collect::<FxHashMap<_, _>>();

    ctgmap_set.records.iter().for_each(|r| {
        if r.q_dup {
            return;
        };
        let e = ctg_target_hit_len.entry(r.q_name.clone()).or_default();
        let e2 = e.entry(r.t_name.clone()).or_default();
        *e2 += (r.qe as i32 - r.qs as i32).unsigned_abs();
    });

    let mut ctg2tgt = FxHashMap::<String, String>::default();
    let mut tgt2ctg = FxHashMap::<String, String>::default();

    ctg_target_hit_len.into_iter().for_each(|(ctg, tgt_len)| {
        let mut tgt_len = tgt_len.into_iter().collect::<Vec<_>>();
        if !tgt_len.is_empty() {
            tgt_len.sort_by(|a, b| b.1.cmp(&a.1));
            let tgt = tgt_len[0].0.clone();
            //println!("DBG: {} {:?}", ctg, tgt_len);
            ctg2tgt.insert(ctg.clone(), tgt.clone());
            tgt2ctg.insert(tgt.clone(), ctg.clone());
        };
    });

    let mut tgt_to_records = FxHashMap::<String, Vec<CtgMapRec>>::default();
    let mut qry_to_alt_tgt_records = FxHashMap::<String, Vec<CtgMapRec>>::default();
    let mut tgt_to_alt_qry_records = FxHashMap::<String, Vec<CtgMapRec>>::default();
    ctgmap_set.records.iter().for_each(|r| {
        if r.q_dup {
            return;
        };
        if *ctg2tgt.get(&r.q_name).unwrap() != r.t_name {
            let e = qry_to_alt_tgt_records.entry(r.q_name.clone()).or_default();
            e.push((*r).clone());
            let e = tgt_to_alt_qry_records.entry(r.t_name.clone()).or_default();
            e.push((*r).clone());
            return;
        }
        let e = tgt_to_records.entry(r.t_name.clone()).or_default();
        e.push((*r).clone());
    });

    let target_padding = 1.5e6;
    let mut offset = 0_f64;
    let target_aln_blocks = ctgmap_set
        .target_length
        .iter()
        .flat_map(|(id, t_name, t_len)| {
            if let Some(target_ctg) = args.ctg.as_ref() {
                if target_ctg != "summary" && target_ctg != t_name {
                    return None;
                }
            };
            let mut q_len_sum = 0.0;
            let mut q_set = FxHashSet::<String>::default();
            tgt_to_records
                .get(t_name)
                .unwrap_or(&vec![])
                .iter()
                .for_each(|record| {
                    let q_len = query_length.get(&record.q_name).unwrap();
                    if !q_set.contains(&record.q_name) {
                        q_set.insert(record.q_name.clone());
                        q_len_sum += *q_len as f64;
                    };
                });

            if let Some(records) = tgt_to_records.get(t_name) {
                let out = (*id, t_name.clone(), *t_len, offset, records);
                offset += (*t_len as f64).max(q_len_sum) + target_padding;
                Some(out)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    let svg_box_height = if args.ctg.is_none() { 3500 } else { 180 };

    // start to construct the SVG element
    let mut document = Document::new()
        .set(
            "viewBox",
            (
                -args.panel_width * 0.05,
                -50,
                args.panel_width * 0.95 * 2.0,
                svg_box_height,
            ),
        )
        .set("width", args.panel_width * 2.0)
        .set("height", svg_box_height)
        .set("preserveAspectRatio", "none")
        .set("id", "WholeGenomeViwer")
        .set("overflow", "visible");

    let scaling_factor = if let Some(total_target_bases) = args.total_target_bases {
        args.panel_width * 0.8 / total_target_bases
    } else {
        args.panel_width * 0.8 / offset
    };

    let mut plot_overview = || {
        target_aln_blocks
            .iter()
            .for_each(|target_aln_block_records| {
                let t_offset = target_aln_block_records.3;
                let t_name = target_aln_block_records.1.clone();
                let mut group = element::Group::new().set("id", format!("overview_{}", t_name));
                let b = t_offset * scaling_factor;
                let e = (t_offset + target_aln_block_records.2 as f64) * scaling_factor;
                let w = 4.0 + ((target_aln_block_records.0 + 1) % 2) as f64 * 1.5;
                let path_str = format!("M {b:0.4} 6 L {e:0.4} 6");
                let path = element::Path::new()
                    .set("stroke", "#000")
                    .set("stroke-width", format!("{w}"))
                    .set("opacity", 0.7)
                    .set("stroke-opacity", 0.7)
                    .set("d", path_str);
                group.append(path);

                let text = element::Text::new(target_aln_block_records.1.clone())
                    .set("x", b)
                    .set("y", 0)
                    .set("font-size", "6px")
                    .set("font-family", "monospace");
                group.append(text);

                if let Some(ref_highlight) = &ref_highlight {
                    if let Some(regions) = ref_highlight.get(&t_name) {
                        regions.iter().for_each(|(bgn, end)| {
                            let b = (t_offset + *bgn as f64) * scaling_factor;
                            let e = (t_offset + *end as f64) * scaling_factor;
                            let path_str = format!("M {b:0.4} 3 L {e:0.4} 3");
                            let path = element::Path::new()
                                .set("stroke", "#F00")
                                .set("stroke-width", 6)
                                .set("opacity", 0.7)
                                .set("stroke-opacity", 0.7)
                                .set("d", path_str);
                            group.append(path);
                        });
                    }
                };

                let mut best_query_block = FxHashMap::<String, CtgMapRec>::default();
                target_aln_block_records.4.iter().for_each(|record| {
                    let e = best_query_block
                        .entry(record.q_name.clone())
                        .or_insert(record.clone());
                    if (e.qs as i32 - e.qe as i32).abs()
                        < (record.qs as i32 - record.qe as i32).abs()
                    {
                        *e = record.clone();
                    }
                });

                let mut best_query_block = best_query_block.values().collect::<Vec<_>>();
                best_query_block.sort_by_key(|&v| v.ts);
                let mut q_offset = 0.0;
                let mut q_offset_map = FxHashMap::<String, f64>::default();
                best_query_block.into_iter().for_each(|record| {
                    let q_len = query_length.get(&record.q_name).unwrap();
                    if !q_offset_map.contains_key(&record.q_name) {
                        q_offset_map.insert(record.q_name.clone(), q_offset);

                        let b = (t_offset + q_offset) * scaling_factor;
                        let e = (t_offset + q_offset + *q_len as f64) * scaling_factor;
                        let y = 95.0;
                        let path_str = format!("M {b:0.4} {y:0.4} L {e:0.4} {y:0.4}");
                        let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
                        let path = element::Path::new()
                            .set("stroke", color)
                            .set("stroke-width", "5")
                            .set("opacity", 0.7)
                            .set("stroke-opacity", 0.7)
                            .set("d", path_str);
                        group.append(path);

                        q_offset += *q_len as f64;
                    };
                });

                target_aln_block_records.4.iter().for_each(|record| {
                    if record.t_dup && record.q_dup {
                        return;
                    };

                    let q_len = query_length.get(&record.q_name).unwrap();

                    let ts = record.ts as f64 + t_offset;
                    let te = record.te as f64 + t_offset;

                    let (qs, qe) = if record.ctg_orientation == 1 {
                        (q_len - record.qe, q_len - record.qs)
                    } else {
                        (record.qs, record.qe)
                    };

                    // let qs = record.qs;
                    // let qe = record.qe;
                    let (qs, qe) = if record.orientation != record.ctg_orientation {
                        (qe, qs)
                    } else {
                        (qs, qe)
                    };
                    let offset = q_offset_map.get(&record.q_name).unwrap();
                    let qs = qs as f64 + t_offset + offset;
                    let qe = qe as f64 + t_offset + offset;
                    let ts = ts * scaling_factor;
                    let te = te * scaling_factor;
                    let qs = qs * scaling_factor;
                    let qe = qe * scaling_factor;
                    // println!("{:?}", record);
                    // println!("{} {} {} {}", ts, te, qs, qe);

                    let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];

                    let path_str =
                        format!("M {ts:0.4} 10 L {te:0.4} 10 L {qe:0.4} 90 L {qs:0.4} 90 Z");
                    let path = element::Path::new()
                        .set("fill", color)
                        .set("stroke", "#000")
                        .set("stroke-width", 0.25)
                        .set("opacity", 0.7)
                        .set("stroke-opacity", 0.4)
                        .set("d", path_str);
                    group.append(path);
                });
                document.append(group);
            });
    };

    if let Some(target_ctg) = args.ctg.as_ref() {
        if target_ctg.eq("summary") {
            plot_overview();
        };
    } else {
        plot_overview();
    };

    // per chromosome plot

    let mut y_offset = if args.ctg.is_none() { 200.0 } else { 0.0 };
    let scaling_factor = if args.ctg.is_some() {
        scaling_factor
    } else {
        scaling_factor * 12.0
    };

    target_aln_blocks
        .iter()
        .for_each(|target_aln_block_record| {
            let t_name = target_aln_block_record.1.clone();
            if let Some(target_ctg) = args.ctg.as_ref() {
                if t_name != target_ctg.clone() {
                    return;
                };
            };
            let group = match get_chr_svg_group(
                target_aln_block_record,
                scaling_factor,
                &cytobands,
                &ref_highlight,
                &tgt_to_alt_qry_records,
                &ctg2tgt,
                &query_length,
                &qry_to_alt_tgt_records,
            ) {
                Some(value) => value,
                None => return,
            };
            let mut sub_svg = Document::new()
                .set("viewBox", (0, -25, args.panel_width, 130))
                .set("width", args.panel_width)
                .set("height", 130)
                .set("preserveAspectRatio", "none")
                .set("y", y_offset)
                .set("id", t_name.clone())
                .set("class", "chr_view")
                .set("overflow", "visible");

            sub_svg.append(group);
            let text = element::Text::new(target_aln_block_record.1.clone())
                .set("x", 0.0)
                .set("y", y_offset + 20.0)
                .set("font-size", "20px")
                .set("font-family", "monospace");
            document.append(text);
            document.append(sub_svg);
            y_offset += 130.0;
        });

    let mut out_file = if args.svg {
        BufWriter::new(
            File::create(path::Path::new(&args.output_prefix).with_extension("svg"))
                .expect("can't create the HTML output file"),
        )
    } else {
        BufWriter::new(
            File::create(path::Path::new(&args.output_prefix).with_extension("html"))
                .expect("can't create the SVG output file"),
        )
    };
    let mut svg_elment = BufWriter::new(Vec::new());
    svg::write(&mut svg_elment, &document).unwrap();
    if !args.svg {
        let jscript = r#"
        <script>
        document.addEventListener('readystatechange', event => {
            if (event.target.readyState === "complete") {
                var views = document.getElementsByClassName("chr_view");
                for (let i = 0; i < views.length; i++) {
                    views[i].addEventListener('mousedown', function(event) {
                        event.preventDefault();
                        const viewBoxValues = views[i].getAttribute('viewBox').split(' ').map(val => parseFloat(val));
                        let viewBox = { x: viewBoxValues[0], y: viewBoxValues[1], width: viewBoxValues[2], height: viewBoxValues[3] };
                        console.log(event);
                        if (event.button != 0) {
                            return;
                        }
                        if (event.altKey) {
                            scalingFactor = 1.25;
                        } else {
                            scalingFactor = 0.8; 
                        };
                        viewBox.width *= scalingFactor;
                        views[i].setAttribute('viewBox', `${viewBox.x} ${viewBox.y} ${viewBox.width} ${viewBox.height}`);
                    });
                };
            }
        });
        </script>
        "#;
        writeln!(out_file, r#"<html><body>"#).expect("can't write the output html file");
        writeln!(out_file, "{}", jscript).expect("can't write the output html file");
        writeln!(out_file, r#"<div style="overflow:scroll;">"#)
            .expect("can't write the output html file");
    };

    writeln!(
        out_file,
        "{}",
        String::from_utf8_lossy(&svg_elment.into_inner().unwrap())
    )
    .expect("can't write the output HTML or SVG file");

    if !args.svg {
        writeln!(out_file, "</div></body></html>").expect("can't write the output html file");
    };

    Ok(())
}

fn get_chr_svg_group(
    target_aln_block_record: &(u32, String, u32, f64, &Vec<CtgMapRec>),
    scaling_factor: f64,
    cytobands: &Option<CytoBands>,
    ref_highlight: &Option<FxHashMap<String, Vec<(u32, u32)>>>,
    tgt_to_alt_qry_records: &FxHashMap<String, Vec<CtgMapRec>>,
    ctg2tgt: &FxHashMap<String, String>,
    query_length: &FxHashMap<String, u32>,
    qry_to_alt_tgt_records: &FxHashMap<String, Vec<CtgMapRec>>,
) -> Option<element::Group> {
    let t_name = target_aln_block_record.1.clone();
    let mut group = element::Group::new();
    let t_offset = 0.0;
    let t_len = target_aln_block_record.2;
    let y = 6.0;
    let mut draw_plain_ref_track = || {
        let b = t_offset * scaling_factor;
        let e = (t_offset + t_len as f64) * scaling_factor;
        // let w = 4.0 + ((target_aln_block_records.0 + 1) % 2) as f64 * 1.5;
        let path_str = format!("M {b:0.4} {y:0.4} L {e:0.4} {y:0.4}");
        let path = element::Path::new()
            .set("stroke", "#000")
            .set("stroke-width", 8)
            .set("opacity", 0.7)
            .set("stroke-opacity", 0.7)
            .set("d", path_str);
        group.append(path);
    };
    if let Some(cytobands) = cytobands.as_ref() {
        if let Some(cyto_records) = cytobands.cytobands.get(&t_name) {
            cyto_records.iter().for_each(|(cs, ce, c_name, band)| {
                let b = (t_offset + *cs as f64) * scaling_factor;
                let e = (t_offset + *ce as f64) * scaling_factor;
                let mut color = if band.starts_with("gpos") {
                    "#000"
                } else {
                    "#AAA"
                };
                if band == "acen" {
                    color = "#FF0";
                };
                let path_str = format!("M {b:0.4} {y:0.4} L {e:0.4} {y:0.4}");
                let mut path = element::Path::new()
                    .set("stroke", color)
                    .set("stroke-width", 8)
                    .set("opacity", 0.7)
                    .set("stroke-opacity", 0.7)
                    .set("d", path_str);
                path.append(element::Title::new(c_name.clone()));
                group.append(path);
            })
        } else {
            draw_plain_ref_track()
        };
    } else {
        draw_plain_ref_track()
    }
    if let Some(ref_highlight) = ref_highlight.as_ref() {
        if let Some(regions) = ref_highlight.get(&t_name) {
            let y2 = y - 8.0;
            regions.iter().for_each(|(bgn, end)| {
                let b = (t_offset + *bgn as f64) * scaling_factor;
                let e = (t_offset + *end as f64) * scaling_factor;
                let path_str = format!("M {b:0.4} {y2:0.4} L {e:0.4} {y2:0.4}");
                let mut path = element::Path::new()
                    .set("stroke", "#F00")
                    .set("stroke-width", 6)
                    .set("opacity", 0.7)
                    .set("stroke-opacity", 0.7)
                    .set("d", path_str);
                path.append(element::Title::new(format!("{}-{}", bgn, end)));
                group.append(path);
            });
        }
    };

    if let Some(tgt_to_alt_qry_records) = tgt_to_alt_qry_records.get(&target_aln_block_record.1) {
        let t_offset = 0.0;
        tgt_to_alt_qry_records.iter().for_each(|record| {
            let b = (t_offset + record.ts as f64) * scaling_factor;
            let e = (t_offset + record.te as f64) * scaling_factor;
            let y = 14.0;
            let path_str = format!("M {b:0.4} {y:0.4} L {e:0.4} {y:0.4}");
            let mut path = element::Path::new()
                .set("stroke", "#000")
                .set("stroke-width", 8)
                .set("opacity", 0.7)
                .set("stroke-opacity", 0.7)
                .set("d", path_str);
            let na = "N/A".to_string();
            let q_tgt = ctg2tgt.get(&record.q_name).unwrap_or(&na);
            path.append(element::Title::new(format!(
                "{} to {} with {}:{}-{}",
                record.t_name, q_tgt, record.q_name, record.qs, record.qe
            )));
            group.append(path);
        })
    };
    let mut best_query_block = FxHashMap::<String, CtgMapRec>::default();
    target_aln_block_record.4.iter().for_each(|record| {
        let e = best_query_block
            .entry(record.q_name.clone())
            .or_insert(record.clone());
        if (e.qs as i32 - e.qe as i32).abs() < (record.qs as i32 - record.qe as i32).abs() {
            *e = record.clone();
        }
    });
    let mut best_query_block = best_query_block.values().collect::<Vec<_>>();
    best_query_block.sort_by_key(|&v| v.ts);
    let mut q_offset = 0.0;
    let mut q_offset_map = FxHashMap::<String, f64>::default();
    best_query_block.into_iter().for_each(|record| {
        let q_len = query_length.get(&record.q_name).unwrap();
        if !q_offset_map.contains_key(&record.q_name) {
            let ctg_aln_orientation = record.ctg_orientation;
            q_offset_map.insert(record.q_name.clone(), q_offset);

            let b = (t_offset + q_offset) * scaling_factor;
            let e = (t_offset + q_offset + *q_len as f64) * scaling_factor;
            let y = 95.0;
            let path_str = format!("M {b:0.4} {y:0.4} L {e:0.4} {y:0.4}");
            let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
            let mut path = element::Path::new()
                .set("stroke", color)
                .set("stroke-width", 8)
                .set("opacity", 0.7)
                .set("stroke-opacity", 0.7)
                .set("d", path_str);
            path.append(element::Title::new(record.q_name.clone()));
            group.append(path);

            if let Some(qry_to_alt_tgt_records) = qry_to_alt_tgt_records.get(&record.q_name) {
                qry_to_alt_tgt_records.iter().for_each(|record| {
                    let qe = if ctg_aln_orientation == 0 {
                        record.qs
                    } else {
                        q_len - record.qs
                    };
                    let qs = if ctg_aln_orientation == 0 {
                        record.qe
                    } else {
                        q_len - record.qe
                    };
                    let b = (t_offset + q_offset + qs as f64) * scaling_factor;
                    let e = (t_offset + q_offset + qe as f64) * scaling_factor;
                    let y = 105.0;
                    let path_str = format!("M {b:0.4} {y:0.4} L {e:0.4} {y:0.4}");
                    let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
                    let mut path = element::Path::new()
                        .set("stroke", color)
                        .set("stroke-width", 8)
                        .set("opacity", 0.7)
                        .set("stroke-opacity", 0.7)
                        .set("d", path_str);
                    path.append(element::Title::new(format!(
                        "{}@{}:{}-{}",
                        record.q_name, record.t_name, record.ts, record.te
                    )));
                    group.append(path);
                });
            };

            q_offset += *q_len as f64;
        };
    });
    target_aln_block_record.4.iter().for_each(|record| {
        if record.t_dup && record.q_dup {
            return;
        };

        let q_len = query_length.get(&record.q_name).unwrap();

        let ts = record.ts as f64 + t_offset;
        let te = record.te as f64 + t_offset;

        let (qs, qe) = if record.ctg_orientation == 1 {
            (q_len - record.qe, q_len - record.qs)
        } else {
            (record.qs, record.qe)
        };

        // let qs = record.qs;
        // let qe = record.qe;
        let (qs, qe) = if record.orientation != record.ctg_orientation {
            (qe, qs)
        } else {
            (qs, qe)
        };
        let offset = q_offset_map.get(&record.q_name).unwrap();
        let qs = qs as f64 + t_offset + offset;
        let qe = qe as f64 + t_offset + offset;
        let ts = ts * scaling_factor;
        let te = te * scaling_factor;
        let qs = qs * scaling_factor;
        let qe = qe * scaling_factor;
        // println!("{:?}", record);
        // println!("{} {} {} {}", ts, te, qs, qe);

        let color = CMAP[(calculate_hash(&record.q_name) % 97) as usize];
        let y = 14.0;
        let y2 = 88.0;
        let path_str = format!(
            "M {ts:0.4} {y:0.4} L {te:0.4} {y:0.4} L {qe:0.4} {y2:0.4} L {qs:0.4} {y2:0.4} Z"
        );
        let mut path = element::Path::new()
            .set("fill", color)
            .set("stroke", "#000")
            .set("stroke-width", "0.25")
            .set("opacity", "0.7")
            .set("stroke-opacity", "0.4")
            .set("d", path_str);
        let orientation = if record.orientation == 0 { '+' } else { '-' };
        let t_dup_mark = if record.t_dup { 1 } else { 0 };
        let q_dup_mark = if record.q_dup { 1 } else { 0 };
        path.append(element::Title::new(format!(
            "{}:{}-{} @ {}:{}-{} {}:{}:{}",
            record.t_name,
            record.ts,
            record.te,
            record.q_name,
            record.qs,
            record.qe,
            orientation,
            t_dup_mark,
            q_dup_mark
        )));

        group.append(path);
    });
    Some(group)
}
