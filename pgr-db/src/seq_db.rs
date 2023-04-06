use crate::agc_io::AGCFile;
use crate::fasta_io::{reverse_complement, FastaReader, SeqRec};
use crate::graph_utils::{AdjList, AdjPair, ShmmrGraphNode};
use crate::shmmrutils::{match_reads, sequence_to_shmmrs, DeltaPoint, ShmmrSpec, MM128};
use anyhow::Result;
use bincode::{config, Decode, Encode};
use byteorder::{ByteOrder, LittleEndian, WriteBytesExt};
use flate2::bufread::MultiGzDecoder;
use flate2::write::DeflateEncoder;
use flate2::Compression;
//use libflate::deflate::{Decoder, EncodeOptions, Encoder};
//use libflate::lz77::DefaultLz77Encoder;
use petgraph::graphmap::DiGraphMap;
use petgraph::visit::Dfs;
use petgraph::EdgeDirection::{Incoming, Outgoing};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use zstd::stream::{decode_all, encode_all};

use std::fmt;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

pub const KMERSIZE: u32 = 56;
pub const SHMMRSPEC: ShmmrSpec = ShmmrSpec {
    w: 80,
    k: KMERSIZE,
    r: 4,
    min_span: 64,
    sketch: true,
};

pub type Bases = Vec<u8>;
pub type AlnSegments = (u32, bool, u32, Vec<AlnSegment>); //(refFragID, orientation, SeqLength, AlnSegments)

#[derive(Debug, Clone, Decode, Encode)]
pub enum AlnSegment {
    // this still use a lot of space, we will find way to reduce the memory footprint later
    FullMatch,
    // u16 should be enough, the max span should be less than 128 * 144 = 18423 * 2 < 2**16
    Match(u32, u32),
    Insertion(u8),
}
#[allow(clippy::large_enum_variant)]
enum GZFastaReader {
    GZFile(FastaReader<BufReader<MultiGzDecoder<BufReader<File>>>>),
    RegularFile(FastaReader<BufReader<BufReader<File>>>),
}

#[derive(Debug, Clone, Decode, Encode)]
pub enum Fragment {
    AlnSegments(AlnSegments),
    Prefix(Bases),
    Internal(Bases),
    Suffix(Bases),
}

pub const FRAG_SHIFT: usize = 4;
pub const FRAG_GROUP_MAX: usize = 1 << FRAG_SHIFT;
#[derive(Debug, Clone, Decode, Encode)]
pub struct FragmentGroup {
    pub seqs: Vec<Vec<u8>>,
    seq_len: Vec<usize>,
    total_len: usize,
    pub compressed_data: Vec<u8>,
    pub compressed: bool,
}

impl FragmentGroup {
    pub fn new() -> Self {
        let seqs = Vec::new();
        let compressed_data = Vec::new();
        let seq_len = Vec::new();
        FragmentGroup {
            seqs,
            seq_len,
            total_len: 0,            
            compressed_data,
            compressed: false,
        }
    }

    pub fn compress(&mut self) {
        if self.compressed == true {
            return;
        }
        let data = self.seqs.iter().flat_map(|v| v.clone()).collect::<Vec<u8>>();
        self.compressed_data = encode_all(&data[..], 1).unwrap();
        self.compressed = true;
        self.seqs.clear();
        /*
        println!(
            "compress ratio {}/{}={} {}/{}={}",
            self.total_len,
            self.compressed_data.len(),
            self.total_len as f32 / self.compressed_data.len() as f32,
            data.len(),
            self.compressed_data.len(),
            data.len() as f32 / self.compressed_data.len() as f32,
        );
        */
    }

    pub fn add_frag(&mut self, v: &[u8]) -> Option<usize> {
        if self.seqs.len() >= FRAG_GROUP_MAX {
            if !self.compressed {
                self.compress()
            };
            None
        } else if self.compressed {
            None
        } else {
            let length = self.seqs.len();
            self.total_len += v.len();
            let single_compressed_seq = v.to_vec();
            self.seq_len.push(single_compressed_seq.len());
            self.seqs.push(single_compressed_seq);
            Some(length)
        }

    }

    pub fn get_frag(&self, sub_idx: u32) -> Vec<u8> {
        if !self.compressed {
            self.seqs[sub_idx as usize].clone()
        } else {
            let decoded_data = decode_all(&self.compressed_data[..]).unwrap();
            let mut offset = 0;
            for sidx in 0..sub_idx as usize {
                offset += self.seq_len[sidx];
            };

            decoded_data
                [offset..offset + self.seq_len[sub_idx as usize]].to_vec()
            
        }
    }
}

impl fmt::Display for Fragment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Fragment::AlnSegments(d) => write!(
                f,
                "RefFragId:{} Orientation:{:?} len:{:?} AlignSegs:{:?}",
                d.0, d.1, d.2, d.3
            ),
            Fragment::Prefix(b) => write!(f, "Seq:{} AlignSegs:None", String::from_utf8_lossy(b)),
            Fragment::Internal(b) => write!(f, "Seq:{} AlignSegs:None", String::from_utf8_lossy(b)),
            Fragment::Suffix(b) => write!(f, "Seq:{} AlignSegs:None", String::from_utf8_lossy(b)),
        }
    }
}

pub type ShmmrPair = (u64, u64);

pub type Fragments = Vec<Fragment>;
pub type FragmentSignature = (u32, u32, u32, u32, u8); //frg_id, seq_id, bgn, end, orientation(to shimmer pair)
pub type ShmmrToFrags = FxHashMap<ShmmrPair, Vec<FragmentSignature>>;

pub trait GetSeq {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8>;
    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8>;
}

#[derive(Debug, Clone, Decode, Encode)]
pub struct CompactSeq {
    pub source: Option<String>,
    pub name: String,
    pub id: u32,
    pub seq_frags: Vec<u32>, // (start, len)
    pub len: usize,
}

#[derive(Debug, Clone)]
pub struct CompactSeqDB {
    pub shmmr_spec: ShmmrSpec,
    pub seqs: Vec<CompactSeq>,
    pub frag_map: ShmmrToFrags,
    pub frag_groups: Option<Vec<FragmentGroup>>,
}

pub fn pair_shmmrs(shmmrs: &Vec<MM128>) -> Vec<(&MM128, &MM128)> {
    if shmmrs.len() < 2 {
        return vec![];
    }
    let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
        .iter()
        .zip(shmmrs[1..shmmrs.len()].iter())
        .collect::<Vec<_>>();
    shmmr_pairs
}

pub fn deltas_to_aln_segs(
    deltas: &Vec<DeltaPoint>,
    endx: usize,
    endy: usize,
    base_frg: &Vec<u8>,
    frg: &Vec<u8>,
) -> Vec<AlnSegment> {
    let mut aln_segs = Vec::<AlnSegment>::new();
    if deltas.is_empty() && base_frg.len() == frg.len() {
        aln_segs.push(AlnSegment::FullMatch);
        //println!("aln_segs: {:?}", aln_segs);
        return aln_segs;
    }

    let mut x = endx;
    let mut y = endy;

    for yy in (y..frg.len()).rev() {
        aln_segs.push(AlnSegment::Insertion(frg[yy]));
    }

    for d in deltas.iter() {
        let x1 = d.x as usize;
        let y1 = d.y as usize;
        if x1 < x {
            aln_segs.push(AlnSegment::Match(x1 as u32, x as u32));
        }
        x = x1;
        y = y1;
        if d.dk > 0 {
            x -= d.dk as usize; // deletion from the base_frg
        } else {
            for yy in 0..(-d.dk) as usize {
                aln_segs.push(AlnSegment::Insertion(frg[y - yy - 1]));
            }
        }
    }
    if x != 0 {
        aln_segs.push(AlnSegment::Match(0, x as u32));
    };
    aln_segs.reverse();
    //println!("aln_segs: {:?}", aln_segs);
    aln_segs
}

pub fn reconstruct_seq_from_aln_segs(base_seq: &[u8], aln_segs: &[AlnSegment]) -> Vec<u8> {
    let mut seq = Vec::<u8>::new();
    for s in aln_segs.iter() {
        match s {
            AlnSegment::FullMatch => {
                seq.extend_from_slice(base_seq);
            }
            AlnSegment::Match(x1, x2) => {
                seq.extend_from_slice(&base_seq[*x1 as usize..*x2 as usize]);
            }
            AlnSegment::Insertion(c) => {
                seq.push(*c);
            }
        }
    }
    seq
}

impl CompactSeqDB {
    pub fn new(shmmr_spec: ShmmrSpec) -> Self {
        let seqs = Vec::<CompactSeq>::new();
        let frag_map = ShmmrToFrags::default();
        let frags = None;
        CompactSeqDB {
            shmmr_spec,
            seqs,
            frag_map,
            frag_groups: frags,
        }
    }

    pub fn seq_to_compressed(
        &mut self,
        source: Option<String>,
        name: String,
        id: u32,
        seq: &Vec<u8>,
        shmmrs: Vec<MM128>,
        try_compress: bool,
    ) -> CompactSeq {
        let mut seq_frags = Vec::<u32>::new();

        assert!(self.frag_groups.is_some());
        let frag_groups: &mut Vec<FragmentGroup> = self.frag_groups.as_mut().unwrap();

        let mut frag_group_id = frag_groups.len() as u32;

        //assert!(shmmrs.len() > 0);
        if shmmrs.is_empty() {
            let mut frag_group = FragmentGroup::new();
            let sub_idx = frag_group.add_frag(&seq[..]).unwrap(); // unwrap, first element
            assert!(sub_idx == 0);
            frag_groups.push(frag_group);
            seq_frags.push(((frag_group_id << FRAG_SHIFT) | (sub_idx as u32)) << 2 | 0b00);

            return CompactSeq {
                source,
                name,
                id,
                seq_frags,
                len: seq.len(),
            };
        }

        let mut seq_len = 0_usize;
        // prefix
        let end = (shmmrs[0].pos() + 1) as usize;

        let mut frag_group = FragmentGroup::new();
        let sub_idx = frag_group.add_frag(&seq[..end]).unwrap(); // unwrap, the 0th element
        assert!(sub_idx == 0);
        frag_groups.push(frag_group);
        seq_frags.push((frag_group_id << FRAG_SHIFT | sub_idx as u32) << 2 | 0b00);
        seq_len += end;
        frag_group_id += 1;

        pair_shmmrs(&shmmrs).iter().for_each(|(shmmr0, shmmr1)| {
            let s0 = shmmr0.hash();
            let s1 = shmmr1.hash();
            let (shmmr_pair, orientation) = if s0 <= s1 {
                ((s0, s1), 0_u8)
            } else {
                ((s1, s0), 1_u8)
            };
            let bgn = shmmr0.pos() + 1;
            let end = shmmr1.pos() + 1;
            let frag_len = end - bgn;
            let frag = &seq[(bgn - self.shmmr_spec.k) as usize..end as usize];
            let mut added = false;

            if self.frag_map.contains_key(&shmmr_pair) {
                let e = self.frag_map.get_mut(&shmmr_pair).unwrap();

                for t_frag in e.iter() {
                    if orientation != t_frag.4 {
                        continue;
                    };
                    let t_frag_id = t_frag.0;
                    let t_frag_group_id = t_frag_id >> FRAG_SHIFT >> 2;
                    if let Some(frag_group) = frag_groups.get_mut(t_frag_group_id as usize) {
                        if let Some(sub_idx) = frag_group.add_frag(frag) {
                            let frag_id =
                                (t_frag_group_id << FRAG_SHIFT | sub_idx as u32) << 2 | 0b01;
                            seq_frags.push(frag_id);
                            seq_len += frag_len as usize;
                            e.push((frag_id, id, bgn, end, orientation));
                            added = true;
                            break;
                        } else {
                            let mut frag_group = FragmentGroup::new();
                            let sub_idx = frag_group.add_frag(frag).unwrap(); // unwrap, first element
                            frag_groups.push(frag_group);
                            let frag_id =
                                (frag_group_id << FRAG_SHIFT | sub_idx as u32) << 2 | 0b01;
                            seq_frags.push(frag_id);
                            frag_group_id += 1;
                            seq_len += frag_len as usize;
                            e.push((frag_id, id, bgn, end, orientation));
                            added = true;
                            break;
                        }
                    }
                }
            };
            if !added {
                let mut frag_group = FragmentGroup::new();
                let sub_idx = frag_group.add_frag(frag).unwrap(); // unwrap, first element
                assert!(sub_idx == 0);
                frag_groups.push(frag_group);
                let frag_id = (frag_group_id << FRAG_SHIFT | sub_idx as u32) << 2 | 0b01;
                seq_frags.push(frag_id);
                self.frag_map
                    .insert(shmmr_pair, vec![(frag_id, id, bgn, end, orientation)]);
                frag_group_id += 1;
                seq_len += frag_len as usize;
            }
        });

        // suffix
        let bgn = (shmmrs[shmmrs.len() - 1].pos() + 1) as usize;
        let frag = &seq[bgn..];
        let mut frag_group = FragmentGroup::new();
        let sub_idx = frag_group.add_frag(frag).unwrap(); // unwrap, first element
        frag_groups.push(frag_group);
        seq_frags.push((frag_group_id << FRAG_SHIFT | sub_idx as u32) << 2 | 0b10);
        //frag_group_id += 1;
        seq_len += frag.len();

        assert_eq!(seq_len, seq.len());
        CompactSeq {
            source,
            name,
            id,
            seq_frags,
            len: seq.len(),
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn seq_to_index(
        source: Option<String>,
        name: String,
        id: u32,
        seqlen: usize,
        shmmrs: Vec<MM128>,
    ) -> (CompactSeq, Vec<((u64, u64), u32, u32, u8)>) {
        //assert!(shmmrs.len() > 0);
        if shmmrs.is_empty() {
            return (
                CompactSeq {
                    source,
                    name,
                    id,
                    seq_frags: Vec::<_>::new(),
                    len: seqlen,
                },
                vec![],
            );
        }

        let shmmr_pairs = shmmrs[0..shmmrs.len() - 1]
            .iter()
            .zip(shmmrs[1..shmmrs.len()].iter())
            .collect::<Vec<_>>();

        let internal_frags: Vec<((u64, u64), u32, u32, u8)> = shmmr_pairs
            .par_iter()
            .map(|(shmmr0, shmmr1)| {
                let s0 = shmmr0.hash();
                let s1 = shmmr1.hash();
                let (shmmr_pair, orientation) = if s0 <= s1 {
                    ((s0, s1), 0_u8)
                } else {
                    ((s1, s0), 1_u8)
                };
                let bgn = shmmr0.pos() + 1;
                let end = shmmr1.pos() + 1;
                (shmmr_pair, bgn, end, orientation)
            })
            .collect::<Vec<_>>();

        let seq_frags: Vec<u32> = (0..(shmmr_pairs.len() as u32)).collect();
        (
            CompactSeq {
                source,
                name,
                id,
                seq_frags,
                len: seqlen,
            },
            internal_frags,
        )
    }

    fn get_fastx_reader(&mut self, filepath: String) -> Result<GZFastaReader, std::io::Error> {
        let file = File::open(&filepath)?;
        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!("input file: {} detected as gz-compressed file", filepath);
                is_gzfile = true;
            }
        }
        drop(reader);

        let file = File::open(&filepath)?;
        let reader = BufReader::new(file);
        let gz_buf = BufReader::new(MultiGzDecoder::new(reader));

        let file = File::open(&filepath)?;
        let reader = BufReader::new(file);
        let std_buf = BufReader::new(reader);

        if is_gzfile {
            drop(std_buf);
            Ok(GZFastaReader::GZFile(
                FastaReader::new(gz_buf, &filepath, 1 << 14, true).unwrap(),
            ))
        } else {
            drop(gz_buf);
            Ok(GZFastaReader::RegularFile(
                FastaReader::new(std_buf, &filepath, 1 << 14, true).unwrap(),
            ))
        }
    }

    fn get_shmmrs_from_seqs(
        &mut self,
        seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>,
    ) -> Vec<(u32, Vec<MM128>)> {
        let all_shmmers = seqs
            .par_iter()
            .map(|(sid, _, _, seq)| {
                let shmmrs = sequence_to_shmmrs(*sid, seq, &self.shmmr_spec, false);
                //let shmmrs = sequence_to_shmmrs2(*sid, &seq, 80, KMERSIZE, 4);
                (*sid, shmmrs)
            })
            .collect::<Vec<(u32, Vec<MM128>)>>();
        all_shmmers
    }

    fn load_seq_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = self.seqs.len() as u32;
        if self.frag_groups.is_none() {
            self.frag_groups = Some(Vec::<FragmentGroup>::new());
        };

        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let source = rec.source.clone();
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, source, seqname, rec.seq));
                    sid += 1;
                } else {
                    end_ext_loop = true;
                    break;
                }
                count += 1;
                if count > 128 {
                    break;
                }
            }

            self.load_seqs_from_seq_vec(&seqs);
            if end_ext_loop {
                break;
            }
        }
    }

    pub fn load_seqs_from_seq_vec(&mut self, seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>) {
        if self.frag_groups.is_none() {
            self.frag_groups = Some(Vec::<FragmentGroup>::new());
        }
        let all_shmmers = self.get_shmmrs_from_seqs(seqs);
        seqs.iter()
            .zip(all_shmmers)
            .for_each(|((sid, source, seqname, seq), (_sid, shmmrs))| {
                let compress_seq = self.seq_to_compressed(
                    source.clone(),
                    seqname.clone(),
                    *sid,
                    seq,
                    shmmrs,
                    true,
                );
                self.seqs.push(compress_seq);
            });
    }

    pub fn load_seqs_from_fastx(&mut self, filepath: String) -> Result<(), std::io::Error> {
        match self.get_fastx_reader(filepath)? {
            #[allow(clippy::useless_conversion)] // the into_iter() is neceesay for dyn patching
            GZFastaReader::GZFile(reader) => self.load_seq_from_reader(&mut reader.into_iter()),

            #[allow(clippy::useless_conversion)] // the into_iter() is neceesay for dyn patching
            GZFastaReader::RegularFile(reader) => {
                self.load_seq_from_reader(&mut reader.into_iter())
            }
        };

        Ok(())
    }

    fn load_index_from_reader(&mut self, reader: &mut dyn Iterator<Item = io::Result<SeqRec>>) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let source = rec.source;
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, source, seqname, rec.seq));
                    sid += 1;
                } else {
                    end_ext_loop = true;
                    break;
                }
                count += 1;
                if count > 128 {
                    break;
                }
            }

            self.load_index_from_seq_vec(&seqs);
            if end_ext_loop {
                break;
            }
        }
    }

    pub fn load_index_from_seq_vec(&mut self, seqs: &Vec<(u32, Option<String>, String, Vec<u8>)>) {
        let all_shmmers = self.get_shmmrs_from_seqs(seqs);
        let seq_names = seqs
            .iter()
            .map(|(_sid, src, n, s)| (src.clone(), n.clone(), s.len()))
            .collect::<Vec<(Option<String>, String, usize)>>();

        /*
        seq_names.iter().zip(all_shmmers).for_each(
            |((source, seq_name, seqlen), (sid, shmmrs))| {
                let compress_seq =
                    self._seq_to_index(source.clone(), seq_name.clone(), sid, *seqlen, shmmrs);
                self.seqs.push(compress_seq);
            },
        );
        */

        seq_names
            .par_iter()
            .zip(all_shmmers)
            .map(|((source, seq_name, seqlen), (sid, shmmrs))| {
                let tmp = self::CompactSeqDB::seq_to_index(
                    source.clone(),
                    seq_name.clone(),
                    sid,
                    *seqlen,
                    shmmrs,
                );
                (sid, tmp.0, tmp.1)
            })
            .collect::<Vec<(u32, CompactSeq, Vec<_>)>>()
            .into_iter()
            .for_each(|(sid, cs, internal_frags)| {
                internal_frags.iter().zip(cs.seq_frags.clone()).for_each(
                    |((shmmr, bgn, end, orientation), frg_id)| {
                        let e = self.frag_map.entry(*shmmr).or_default();
                        e.push((frg_id, sid, *bgn, *end, *orientation));
                    },
                );
                self.seqs.push(cs);
            });
    }

    fn _write_shmmr_vec_from_reader(
        &mut self,
        reader: &mut dyn Iterator<Item = io::Result<SeqRec>>,
        writer: &mut Vec<u8>,
    ) {
        let mut seqs = <Vec<(u32, Option<String>, String, Vec<u8>)>>::new();
        let mut sid = 0;
        loop {
            let mut count = 0;
            let mut end_ext_loop = false;
            seqs.clear();

            loop {
                if let Some(rec) = reader.next() {
                    let rec = rec.unwrap();
                    let source = rec.source;
                    let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                    seqs.push((sid, source, seqname, rec.seq));
                    sid += 1;
                } else {
                    end_ext_loop = true;
                    break;
                }
                count += 1;
                if count > 128 {
                    break;
                }
            }

            self.get_shmmrs_from_seqs(&seqs)
                .iter()
                .map(|(_, v)| v)
                .for_each(|v| {
                    v.iter().for_each(|m| {
                        let _ = writer.write_u64::<LittleEndian>(m.x);
                        let _ = writer.write_u64::<LittleEndian>(m.y);
                    });
                });

            if end_ext_loop {
                break;
            }
        }
    }

    pub fn load_index_from_fastx(&mut self, filepath: String) -> Result<(), std::io::Error> {
        match self.get_fastx_reader(filepath)? {
            #[allow(clippy::useless_conversion)] // the into_iter() is neceesay for dyn patching
            GZFastaReader::GZFile(reader) => self.load_index_from_reader(&mut reader.into_iter()),

            #[allow(clippy::useless_conversion)] // the into_iter() is neceesay for dyn patching
            GZFastaReader::RegularFile(reader) => {
                self.load_index_from_reader(&mut reader.into_iter())
            }
        };

        Ok(())
    }

    pub fn load_index_from_agcfile(&mut self, agcfile: AGCFile) -> Result<(), std::io::Error> {
        //let agcfile = AGCFile::new(filepath);

        self.load_index_from_reader(&mut agcfile.into_iter());
        Ok(())
    }
}

impl CompactSeqDB {
    fn reconstruct_seq_from_frags<I: Iterator<Item = u32>>(&self, frag_ids: I) -> Vec<u8> {
        let mut reconstructed_seq = <Vec<u8>>::new();
        let frag_groups = self.frag_groups.as_ref().unwrap();
        // let mut _p = 0;
        frag_ids.for_each(|frag_id| {
            let t = frag_id & 0b11;
            let sub_idx = (frag_id >> 2) & ((0x01 << FRAG_SHIFT) - 1);
            let frag_group_id = frag_id >> 2 >> FRAG_SHIFT;
            let b = frag_groups
                .get(frag_group_id as usize)
                .unwrap()
                .get_frag(sub_idx);
            //println!("{}:{}", frg_id, sdb.frags[*frg_id as usize]);
            match t {
                0b00 => {
                    //prefix
                    reconstructed_seq.extend_from_slice(&b[..]);
                }
                0b10 => {
                    //suffix
                    reconstructed_seq.extend_from_slice(&b[..]);
                    //_p += b.len();
                }
                0b01 => {
                    // internal
                    reconstructed_seq.extend_from_slice(&b[self.shmmr_spec.k as usize..]);
                    //_p += b.len()-self.shmmr_spec.k as usize;
                }
                _ => (),
            }
        });

        reconstructed_seq
    }

    pub fn get_seq(&self, seq: &CompactSeq) -> Vec<u8> {
        self.reconstruct_seq_from_frags(seq.seq_frags.clone().into_iter())
    }

    /* TODO */
    /*
    pub fn get_sub_seq(&self, seq: &CompactSeq, b: usize, e:usize) -> Vec<u8> {
        vec![]
    }
    */
}

impl GetSeq for CompactSeqDB {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8> {
        let seq = self.seqs.get(sid as usize).unwrap();
        self.reconstruct_seq_from_frags(seq.seq_frags.clone().into_iter())
    }

    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let seq = self.get_seq_by_id(sid);
        seq[bgn as usize..end as usize].into()
    }
}

impl CompactSeqDB {
    pub fn write_shmr_map_index(&self, fp_prefix: String) -> Result<(), std::io::Error> {
        let seq_idx_fp = fp_prefix.clone() + ".midx";
        let data_fp = fp_prefix + ".mdb";
        write_shmr_map_file(&self.shmmr_spec, &self.frag_map, data_fp)?;
        let mut idx_file = BufWriter::new(File::create(seq_idx_fp).expect("file create error"));
        self.seqs
            .iter()
            .try_for_each(|s| -> Result<(), std::io::Error> {
                writeln!(
                    idx_file,
                    "{}\t{}\t{}\t{}",
                    s.id,
                    s.len,
                    s.name,
                    s.source.clone().unwrap_or_else(|| "-".to_string())
                )?;
                Ok(())
            })?;

        Ok(())
    }
}

impl CompactSeqDB {
    pub fn write_to_frag_files(&mut self, file_prefix: String) {
        let mut sdx_file = BufWriter::new(
            File::create(file_prefix.clone() + ".sdx").expect("sdx file creating fail\n"),
        );
        let mut frg_file =
            BufWriter::new(File::create(file_prefix + ".frg").expect("frg file creating fail\n"));

        let config = config::standard();

        //self.seqs.iter().for_each(|s| {
        //    println!("{:?} {:?} {} {}", s.id, s.seq_frags.len(), s.seq_frags[0], s.seq_frags[s.seq_frags.len()-1]);
        //});

        self
        .frag_groups.as_mut().unwrap().par_iter_mut().for_each(|f| f.compress());

        let compressed_frag_groups = self
            .frag_groups
            .as_ref()
            .unwrap()
            .par_iter()
            .map(|f| {
                
                let w = bincode::encode_to_vec(f, config).unwrap();
                let mut compressor = DeflateEncoder::new(Vec::new(), Compression::default());
                compressor.write_all(&w).unwrap();
                let compress_frag = compressor.finish().unwrap();
                compress_frag
            })
            .collect::<Vec<Vec<u8>>>();

        let mut frag_grpup_addr_offeset = vec![];
        let mut offset = 0_usize;
        compressed_frag_groups.iter().for_each(|v| {
            let l = v.len();
            frag_grpup_addr_offeset.push((offset, v.len()));
            offset += l;
            frg_file.write_all(v).expect("frag file writing error\n");
        });

        bincode::encode_into_std_write(
            (frag_grpup_addr_offeset, &self.seqs),
            &mut sdx_file,
            config,
        )
        .expect("sdx file writing error\n");
        //bincode::encode_into_std_write(compressed_frags, &mut frg_file, config)
        //    .expect(" frag file writing error");
    }
}

pub fn frag_map_to_adj_list(
    frag_map: &ShmmrToFrags,
    min_count: usize,
    keeps: Option<Vec<u32>>, // a list of sequence id that we like to keep the sequence in the adj list regardless the coverage
) -> AdjList {
    let mut out = frag_map
        .par_iter()
        .flat_map(|v| {
            v.1.iter()
                .map(|vv| (vv.1, vv.2, vv.3, ShmmrGraphNode(v.0 .0, v.0 .1, vv.4)))
                .collect::<Vec<(u32, u32, u32, ShmmrGraphNode)>>() //(seq_id, bgn, end, (hash0, hash1, orientation))
        })
        .collect::<Vec<(u32, u32, u32, ShmmrGraphNode)>>();
    if out.len() < 2 {
        return vec![];
    }
    out.par_sort();

    let out = if let Some(keeps) = keeps {
        let keeps = FxHashSet::<u32>::from_iter(keeps.into_iter());

        // more or less duplicate code, but this takes the hashset check out of the loop if keeps is None.
        out.into_par_iter()
            .map(|v| {
                if frag_map.get(&(v.3 .0, v.3 .1)).unwrap().len() >= min_count
                    || keeps.contains(&v.0)
                {
                    Some(v)
                } else {
                    None
                }
            })
            .collect::<Vec<Option<(u32, u32, u32, ShmmrGraphNode)>>>()
    } else {
        out.into_par_iter()
            .map(|v| {
                if frag_map.get(&(v.3 .0, v.3 .1)).unwrap().len() >= min_count {
                    Some(v)
                } else {
                    None
                }
            })
            .collect::<Vec<Option<(u32, u32, u32, ShmmrGraphNode)>>>()
    };

    (0..out.len() - 1)
        //.into_par_iter()
        .into_iter()
        .flat_map(|i| {
            if let (Some(v), Some(w)) = (out[i], out[i + 1]) {
                // println!("DBG v: {} {} {} {:?} w: {} {} {} {:?}", v.0, v.1, v.2, v.3, w.0, w.1, w.2, w.3); // XXX
                if v.0 != w.0 || v.2 != w.1 {
                    vec![None]
                } else {
                    vec![
                        Some((v.0, v.3, w.3)),
                        Some((
                            v.0,
                            ShmmrGraphNode(w.3 .0, w.3 .1, 1 - w.3 .2),
                            ShmmrGraphNode(v.3 .0, v.3 .1, 1 - v.3 .2),
                        )),
                    ]
                }
            } else {
                vec![None]
            }
        })
        .filter(|v| v.is_some())
        .map(|v| v.unwrap())
        .collect::<AdjList>() // seq_id, node0, node1
}

pub fn generate_smp_adj_list_for_seq(
    seq: &Vec<u8>,
    sid: u32,
    frag_map: &ShmmrToFrags,
    shmmr_spec: &ShmmrSpec,
    min_count: usize,
) -> AdjList {
    let shmmrs = sequence_to_shmmrs(0, seq, shmmr_spec, false);
    let res = pair_shmmrs(&shmmrs)
        .iter()
        .map(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.x >> 8;
            let s1 = s1.x >> 8;
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            }
        })
        .collect::<Vec<(u64, u64, u32, u32, u8)>>();

    if res.len() < 2 {
        vec![]
    } else {
        (0..res.len() - 1)
            .into_iter()
            .flat_map(|i| {
                let v = res[i];
                let w = res[i + 1];
                if (frag_map.get(&(v.0, v.1)).is_none() || frag_map.get(&(w.0, w.1)).is_none())
                    || (frag_map.get(&(v.0, v.1)).unwrap().len() < min_count
                        || frag_map.get(&(w.0, w.1)).unwrap().len() < min_count)
                    || v.3 != w.2
                {
                    vec![None]
                } else {
                    vec![
                        Some((
                            sid,
                            ShmmrGraphNode(v.0, v.1, v.4),
                            ShmmrGraphNode(w.0, w.1, w.4),
                        )),
                        Some((
                            sid,
                            ShmmrGraphNode(w.0, w.1, 1 - w.4),
                            ShmmrGraphNode(v.0, v.1, 1 - v.4),
                        )),
                    ]
                }
            })
            .flatten()
            .collect::<AdjList>()
    }
}

type PBundleNode = (
    // node, Option<previous_node>, node_weight, is_leaf, global_rank, branch, branch_rank
    ShmmrGraphNode,
    Option<ShmmrGraphNode>,
    u32,
    bool,
    u32,
    u32,
    u32,
);

pub fn sort_adj_list_by_weighted_dfs(
    frag_map: &ShmmrToFrags,
    adj_list: &[AdjPair],
    start: ShmmrGraphNode,
) -> Vec<PBundleNode> {
    use crate::graph_utils::BiDiGraphWeightedDfs;

    let mut g = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut score = FxHashMap::<ShmmrGraphNode, u32>::default();
    adj_list.iter().for_each(|&(_sid, v, w)| {
        let vv = (v.0, v.1);
        let ww = (w.0, w.1);
        let v = ShmmrGraphNode(v.0, v.1, v.2);
        let w = ShmmrGraphNode(w.0, w.1, w.2);
        g.add_edge(v, w, ());

        // println!("DBG: add_edge {:?} {:?}", v, w);
        score
            .entry(v)
            .or_insert_with(|| frag_map.get(&vv).unwrap().len() as u32);
        score
            .entry(w)
            .or_insert_with(|| frag_map.get(&ww).unwrap().len() as u32);
    });

    // println!("DBG: # node: {}, # edgg: {}", g.node_count(), g.edge_count());

    let start = ShmmrGraphNode(start.0, start.1, start.2);

    let mut wdfs_walker = BiDiGraphWeightedDfs::new(&g, start, &score);
    let mut out = vec![];
    while let Some((node, p_node, is_leaf, rank, branch_id, branch_rank)) = wdfs_walker.next(&g) {
        let node_count = *score.get(&node).unwrap();
        let p_node = p_node.map(|pnode| ShmmrGraphNode(pnode.0, pnode.1, pnode.2));
        out.push((
            ShmmrGraphNode(node.0, node.1, node.2),
            p_node,
            node_count,
            is_leaf,
            rank,
            branch_id,
            branch_rank,
        ));
        //println!("DBG, next node: {:?}", node);
    }
    out
}

pub fn get_principal_bundles_from_adj_list(
    frag_map: &ShmmrToFrags,
    adj_list: &[AdjPair],
    path_len_cutoff: usize,
) -> (Vec<Vec<ShmmrGraphNode>>, AdjList) {
    assert!(!adj_list.is_empty());
    // println!("DBG: adj_list[0]: {:?}", adj_list[0]);
    let s = adj_list[0].1;
    let sorted_adj_list = sort_adj_list_by_weighted_dfs(frag_map, adj_list, s);

    // println!("DGB: sorted_adj_list len: {}", sorted_adj_list.len());

    let mut paths: Vec<Vec<ShmmrGraphNode>> = vec![];
    let mut path: Vec<ShmmrGraphNode> = vec![];
    for v in sorted_adj_list.into_iter() {
        path.push(v.0);
        if v.3 {
            // it is a leaf node
            paths.push(path.clone());
            path.clear()
        }
    }

    let long_paths = paths
        .into_iter()
        .filter(|p| p.len() > path_len_cutoff as usize);

    let mut main_bundle_path_vertices = FxHashSet::<(u64, u64)>::default();

    long_paths.for_each(|p| {
        p.into_iter().for_each(|v| {
            main_bundle_path_vertices.insert((v.0, v.1));
        })
    });

    let mut g0 = DiGraphMap::<ShmmrGraphNode, ()>::new();
    let mut filtered_adj_list = AdjList::new();
    adj_list.iter().for_each(|&(sid, v, w)| {
        if main_bundle_path_vertices.contains(&(v.0, v.1))
            && main_bundle_path_vertices.contains(&(w.0, w.1))
        {
            g0.add_edge(
                ShmmrGraphNode(v.0, v.1, v.2),
                ShmmrGraphNode(w.0, w.1, w.2),
                (),
            );
            filtered_adj_list.push((sid, v, w));
        }
    });

    let mut g1 = g0.clone();
    let mut terminal_vertices = FxHashSet::<ShmmrGraphNode>::default();

    for (v, w, _) in g0.all_edges() {
        if g0.neighbors_directed(v, Outgoing).count() > 1 {
            terminal_vertices.insert(v);
        };
        if g0.neighbors_directed(w, Incoming).count() > 1 {
            terminal_vertices.insert(v);
        };
    }

    let mut starts = Vec::<ShmmrGraphNode>::default();
    for v in g1.nodes() {
        if g1.neighbors_directed(v, Incoming).count() == 0 {
            starts.push(v);
        }
    }
    // if the whole graph is a loop
    if starts.is_empty() {
        if let Some(v) = g1.nodes().next() {
            starts.push(v);
        }
    };

    let mut principal_bundles = Vec::<Vec<ShmmrGraphNode>>::new();

    while !starts.is_empty() {
        let s = starts.pop().unwrap();
        let mut dfs = Dfs::new(&g1, s);
        let mut path = Vec::<ShmmrGraphNode>::new();
        while let Some(v) = dfs.next(&g1) {
            if terminal_vertices.contains(&v) {
                path.push(v);
                break;
            } else {
                path.push(v);
            }
        }
        if !path.is_empty() {
            path.iter().for_each(|&v| {
                g1.remove_node(v);
                g1.remove_node(ShmmrGraphNode(v.0, v.1, 1 - v.2));
            });

            /*
            let v = path[path.len()-1];

            for w in g1.neighbors_directed(v, Outgoing) {
                if g1.neighbors_directed(w, Incoming).count() == 0 {
                    starts.push(w);
                }
            }
            */
            starts.clear();
            for v in g1.nodes() {
                if g1.neighbors_directed(v, Incoming).count() == 0 {
                    starts.push(v);
                }
            }

            principal_bundles.push(path);
        }

        // if the whole graph is a loop
        if starts.is_empty() {
            if let Some(v) = g1.nodes().next() {
                starts.push(v);
            }
        };
    }
    principal_bundles.sort_by(|a, b| b.len().partial_cmp(&(a.len())).unwrap());
    (principal_bundles, filtered_adj_list)
}

impl CompactSeqDB {
    pub fn generate_smp_adj_list_from_frag_map(
        &self,
        min_count: usize,
        keeps: Option<Vec<u32>>,
    ) -> AdjList {
        frag_map_to_adj_list(&self.frag_map, min_count, keeps)
    }
}

type FragmentHit = ((u64, u64), (u32, u32, u8), Vec<FragmentSignature>); // ((hash0, hash1), (pos0, pos1, orientation), fragments)

pub fn query_fragment(
    shmmr_map: &ShmmrToFrags,
    frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
) -> Vec<FragmentHit> {
    let shmmrs = sequence_to_shmmrs(0, frag, shmmr_spec, false);
    let query_results = pair_shmmrs(&shmmrs)
        .par_iter()
        .map(|(s0, s1)| {
            let p0 = s0.pos() + 1;
            let p1 = s1.pos() + 1;
            let s0 = s0.hash();
            let s1 = s1.hash();
            if s0 < s1 {
                (s0, s1, p0, p1, 0_u8)
            } else {
                (s1, s0, p0, p1, 1_u8)
            }
        })
        .map(|(s0, s1, p0, p1, orientation)| {
            if let Some(m) = shmmr_map.get(&(s0, s1)) {
                ((s0, s1), (p0, p1, orientation), m.clone())
            } else {
                ((s0, s1), (p0, p1, orientation), vec![])
            }
        })
        .collect::<Vec<_>>();
    query_results
}

pub fn get_match_positions_with_fragment(
    shmmr_map: &ShmmrToFrags,
    frag: &Vec<u8>,
    shmmr_spec: &ShmmrSpec,
) -> FxHashMap<u32, Vec<(u32, u32, u8)>> {
    let mut res = FxHashMap::<u32, Vec<(u32, u32, u8)>>::default();
    query_fragment(shmmr_map, frag, shmmr_spec)
        .into_iter()
        .for_each(|v| {
            let q_direction = v.1 .2;
            v.2.into_iter().for_each(|w| {
                let (_, sid, p0, p1, direction) = w;
                let direction = if direction == q_direction { 0 } else { 1 };
                res.entry(sid).or_default().push((p0, p1, direction));
            });
        });
    res.iter_mut().for_each(|(_k, v)| v.sort());
    res
}

pub fn write_shmr_map_file(
    shmmr_spec: &ShmmrSpec,
    shmmr_map: &ShmmrToFrags,
    filepath: String,
) -> Result<(), std::io::Error> {
    let mut out_file =
        File::create(filepath).expect("open fail while writing the SHIMMER map (.mdb) file\n");
    let mut buf = Vec::<u8>::new();

    buf.extend("mdb".to_string().into_bytes());

    buf.write_u32::<LittleEndian>(shmmr_spec.w as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.k as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.r as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.min_span as u32)?;
    buf.write_u32::<LittleEndian>(shmmr_spec.sketch as u32)?;

    buf.write_u64::<LittleEndian>(shmmr_map.len() as u64)?;
    shmmr_map
        .iter()
        .try_for_each(|(k, v)| -> Result<(), std::io::Error> {
            buf.write_u64::<LittleEndian>(k.0)?;
            buf.write_u64::<LittleEndian>(k.1)?;
            buf.write_u64::<LittleEndian>(v.len() as u64)?;
            v.iter().try_for_each(|r| -> Result<(), std::io::Error> {
                buf.write_u32::<LittleEndian>(r.0)?;
                buf.write_u32::<LittleEndian>(r.1)?;
                buf.write_u32::<LittleEndian>(r.2)?;
                buf.write_u32::<LittleEndian>(r.3)?;
                buf.write_u8(r.4)?;
                Ok(())
            })
        })?;
    let _ = out_file.write_all(&buf);
    Ok(())
}

pub fn read_mdb_file(filepath: String) -> Result<(ShmmrSpec, ShmmrToFrags), io::Error> {
    let mut in_file =
        File::open(filepath).expect("Error while opening the SHIMMER map file (.mdb) file");
    let mut buf = Vec::<u8>::new();

    let mut u64bytes = [0_u8; 8];
    let mut u32bytes = [0_u8; 4];
    in_file.read_to_end(&mut buf)?;
    let mut cursor = 0_usize;
    assert!(buf[0..3] == "mdb".to_string().into_bytes());
    cursor += 3; // skip "mdb"

    let w = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let k = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let r = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let min_span = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let flag = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let sketch = (flag & 0b01) == 0b01;

    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch,
    };
    u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
    let shmmr_key_len = usize::from_le_bytes(u64bytes);
    cursor += 8;
    let mut shmmr_map = ShmmrToFrags::default();
    (0..shmmr_key_len).into_iter().for_each(|_| {
        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k1 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k2 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let vec_len = usize::from_le_bytes(u64bytes);
        cursor += 8;

        let value = (0..vec_len)
            .into_iter()
            .map(|_| {
                let mut v = (0_u32, 0_u32, 0_u32, 0_u32, 0_u8);

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.0 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.1 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.2 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                v.3 = u32::from_le_bytes(u32bytes);
                cursor += 4;

                v.4 = buf[cursor..cursor + 1][0];
                cursor += 1;

                v
            })
            .collect::<Vec<(u32, u32, u32, u32, u8)>>();

        shmmr_map.insert((k1, k2), value);
    });

    Ok((shmmr_spec, shmmr_map))
}

pub fn read_mdb_file_parallel(filepath: String) -> Result<(ShmmrSpec, ShmmrToFrags), io::Error> {
    let mut in_file =
        File::open(filepath).expect("open fail while reading the SHIMMER map (.mdb) file");
    let mut buf = Vec::<u8>::new();

    let mut u64bytes = [0_u8; 8];

    in_file.read_to_end(&mut buf)?;
    let mut cursor = 0_usize;
    assert!(buf[0..3] == "mdb".to_string().into_bytes());
    cursor += 3; // skip "mdb"

    let w = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let k = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let r = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let min_span = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let flag = LittleEndian::read_u32(&buf[cursor..cursor + 4]);
    cursor += 4;
    let sketch = (flag & 0b01) == 0b01;

    let shmmr_spec = ShmmrSpec {
        w,
        k,
        r,
        min_span,
        sketch,
    };
    u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
    let shmmr_key_len = usize::from_le_bytes(u64bytes);
    cursor += 8;
    ShmmrToFrags::default();
    let mut rec_loc = Vec::<(u64, u64, usize, usize)>::new();
    for _ in 0..shmmr_key_len {
        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k1 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let k2 = u64::from_le_bytes(u64bytes);
        cursor += 8;

        u64bytes.clone_from_slice(&buf[cursor..cursor + 8]);
        let vec_len = usize::from_le_bytes(u64bytes);
        cursor += 8;

        let start = cursor;
        cursor += vec_len * 17;
        rec_loc.push((k1, k2, start, vec_len))
    }

    let shmmr_map = rec_loc
        .par_iter()
        .map(|&(k1, k2, start, vec_len)| {
            let mut cursor = start;
            let value = (0..vec_len)
                .into_iter()
                .map(|_| {
                    let mut u32bytes = [0_u8; 4];
                    let mut v = (0_u32, 0_u32, 0_u32, 0_u32, 0_u8);
                    u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                    v.0 = u32::from_le_bytes(u32bytes);
                    cursor += 4;

                    u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                    v.1 = u32::from_le_bytes(u32bytes);
                    cursor += 4;

                    u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                    v.2 = u32::from_le_bytes(u32bytes);
                    cursor += 4;

                    u32bytes.clone_from_slice(&buf[cursor..cursor + 4]);
                    v.3 = u32::from_le_bytes(u32bytes);
                    cursor += 4;

                    v.4 = buf[cursor..cursor + 1][0];
                    cursor += 1;
                    v
                })
                .collect::<Vec<(u32, u32, u32, u32, u8)>>();
            ((k1, k2), value)
        })
        .collect::<FxHashMap<(u64, u64), Vec<(u32, u32, u32, u32, u8)>>>();
    Ok((shmmr_spec, shmmr_map))
}
