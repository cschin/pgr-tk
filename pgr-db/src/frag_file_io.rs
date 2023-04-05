use crate::seq_db::{
    self, read_mdb_file_parallel, CompactSeq, Fragment, FragmentGroup, GetSeq, ShmmrToFrags,
    FRAG_SHIFT,
};
use crate::shmmrutils::ShmmrSpec;
use bincode::config;
use flate2::read::DeflateDecoder;
use memmap::Mmap;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

pub struct CompactSeqDBStorage {
    pub shmmr_spec: ShmmrSpec,
    pub seqs: Vec<CompactSeq>,
    pub frag_map: ShmmrToFrags,
    pub frag_file_prefix: String,
    pub frag_file: Mmap,
    pub frag_group_addr_offsets: Vec<(usize, usize)>, //offset, compress_chunk_size
    pub seq_index: FxHashMap<(String, Option<String>), (u32, u32)>,
    /// a dictionary maps id -> (ctg_name, source, len)
    pub seq_info: FxHashMap<u32, (String, Option<String>, u32)>,
}

impl CompactSeqDBStorage {
    pub fn new(prefix: String) -> Self {
        let frag_file_prefix = prefix;
        let (shmmr_spec, frag_map) =
            read_mdb_file_parallel(frag_file_prefix.clone() + ".mdb").unwrap();
        let mut sdx_file = BufReader::new(
            File::open(frag_file_prefix.clone() + ".sdx").expect("sdx file open error"),
        );
        let config = config::standard();
        let (frag_addr_offsets, seqs): (Vec<(usize, usize)>, Vec<CompactSeq>) =
            bincode::decode_from_std_read(&mut sdx_file, config).expect("read sdx file error");
        let f_file = File::open(frag_file_prefix.clone() + ".frg").expect("frag file open fail");
        let frag_file = unsafe { Mmap::map(&f_file).expect("frag mmap fail") };
        let mut seq_index = FxHashMap::<(String, Option<String>), (u32, u32)>::default();
        let mut seq_info = FxHashMap::<u32, (String, Option<String>, u32)>::default();

        let midx_file = BufReader::new(
            File::open(frag_file_prefix.clone() + ".midx").expect("open midx file fail"),
        );
        midx_file
            .lines()
            .into_iter()
            .try_for_each(|line| -> Result<(), std::io::Error> {
                let line = line.unwrap();
                let mut line = line.as_str().split('\t');
                let sid = line.next().unwrap().parse::<u32>().unwrap();
                let len = line.next().unwrap().parse::<u32>().unwrap();
                let ctg_name = line.next().unwrap().to_string();
                let source = line.next().unwrap().to_string();
                seq_index.insert((ctg_name.clone(), Some(source.clone())), (sid, len));
                seq_info.insert(sid, (ctg_name, Some(source), len));
                Ok(())
            })
            .expect("read midx file fail");

        Self {
            shmmr_spec,
            seqs,
            frag_map,
            frag_file_prefix,
            frag_file,
            frag_group_addr_offsets: frag_addr_offsets,
            seq_index,
            seq_info,
        }
    }

    fn get_seq_from_frag_ids<I: Iterator<Item = u32>>(&self, frag_ids: I) -> Vec<u8> {
        let mut reconstructed_seq = <Vec<u8>>::new();

        let mut _p = 0;
        frag_ids.for_each(|frag_id| {
            let t = frag_id & 0b11;
            let sub_idx = (frag_id >> 2) & 0b1111;
            let frag_group_id = frag_id >> 2 >> FRAG_SHIFT;
            let frag_group = fetch_frag(
                frag_group_id,
                &self.frag_group_addr_offsets,
                &self.frag_file,
            );
            let b = frag_group.get_uncompressed_frag(sub_idx);
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
}

impl GetSeq for CompactSeqDBStorage {
    fn get_seq_by_id(&self, sid: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let seq_frags = self.seqs[sid as usize].seq_frags.clone();
        self.get_seq_from_frag_ids(seq_frags.into_iter())
    }

    fn get_sub_seq_by_id(&self, sid: u32, bgn: u32, end: u32) -> Vec<u8> {
        assert!((sid as usize) < self.seqs.len());
        let seq = self.get_seq_by_id(sid);
        seq[bgn as usize..end as usize].into()
    }
}

fn fetch_frag(
    frag_group_id: u32,
    frag_group_addr_offsets: &[(usize, usize)],
    frag_file: &Mmap,
) -> FragmentGroup {
    let config = config::standard();
    let (offset, size) = frag_group_addr_offsets[frag_group_id as usize];
    let compress_chunk = frag_file[offset..(offset + size as usize)].to_vec();
    let mut deflater = DeflateDecoder::new(&compress_chunk[..]);
    let mut s: Vec<u8> = vec![];
    deflater.read_to_end(&mut s).expect("decompression error");
    let (frag_group, _size): (FragmentGroup, usize) =
        bincode::decode_from_slice::<FragmentGroup, bincode::config::Configuration>(&s[..], config)
            .unwrap();
    frag_group
}
