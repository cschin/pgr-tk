pub const VERSION_STRING: &str = env!("VERSION_STRING");

//pub mod agc_io;
pub mod aln;
//pub mod bindings;
pub mod ec;
pub mod fasta_io;
pub mod frag_file_io;
pub mod gff_db;
pub mod graph_utils;
pub mod kmer_filter;
pub mod seq_db;
//pub mod seqs2variants;
pub mod shmmrutils;

#[cfg(test)]
mod tests {
    use crate::fasta_io::FastaReader;
    use crate::shmmrutils::{self, match_reads, DeltaPoint};
    use flate2::bufread::MultiGzDecoder;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Read};

    use crate::seq_db::{self, deltas_to_aln_segs, reconstruct_seq_from_aln_segs};

    pub fn load_seqs() -> HashMap<String, Vec<u8>> {
        let mut seqs = HashMap::<String, Vec<u8>>::new();
        let filepath = "test/test_data/test_seqs.fa";
        let file = File::open(filepath.to_string()).unwrap();
        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!("input file detected as gz-compressed file",);
                is_gzfile = true;
            }
        }
        drop(reader);

        let file = File::open(&filepath).unwrap();
        let mut reader = BufReader::new(file);
        let gz_buf = &mut BufReader::new(MultiGzDecoder::new(&mut reader));

        let file = File::open(&filepath).unwrap();
        let reader = BufReader::new(file);
        let std_buf = &mut BufReader::new(reader);

        let fastx_buf: &mut dyn BufRead = if is_gzfile {
            drop(std_buf);
            gz_buf
        } else {
            drop(gz_buf);
            std_buf
        };

        let mut fastx_reader =
            FastaReader::new(fastx_buf, &filepath.to_string(), 1 << 12, true).unwrap();
        while let Some(rec) = fastx_reader.next_rec() {
            let rec = rec.unwrap();
            let seqname = String::from_utf8_lossy(&rec.id).into_owned();
            seqs.insert(seqname, rec.seq.clone());
        }
        seqs
    }

    #[test]
    pub fn gz_file_read_test() {
        let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
        let _ = sdb.load_seqs_from_fastx("test/test_data/test_seqs2.fa.gz".to_string());
        println!("{:?}", sdb.seqs[0].seq_frag_range);
    }

    #[test]
    fn load_seq_test() {
        let seqs = load_seqs();
        let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
        let _ = sdb.load_seqs_from_fastx("test/test_data/test_seqs2.fa.gz".to_string());
        //println!("test");
        for seq in sdb.seqs.iter() {
            //println!();
            //println!("{}", seq.name);
            let reconstruct_seq = sdb.get_seq(seq);
            let orig_seq = seqs.get(&seq.name).unwrap();
            if reconstruct_seq != *orig_seq {
                //println!("{}", seq.name);
                //println!("{:?}", reconstruct_seq);
                //println!("{:?}", orig_seq);
                for i in 0..reconstruct_seq.len() {
                    if orig_seq[i] != reconstruct_seq[i] {
                        println!("{} {} {} X", i, orig_seq[i], reconstruct_seq[i]);
                    } else {
                        println!("{} {} {}  ", i, orig_seq[i], reconstruct_seq[i]);
                    }
                }
            };
            assert_eq!(reconstruct_seq, *orig_seq);
        }
    }

    #[test]
    fn reconstruct_test1() {
        let base_frg = "TATTTATATTTATTTATATATATTTATATATTTATATATATATTTATATATAAATAT"
            .as_bytes()
            .to_vec();
        let frg = "TTTTTATTTTTTTAATTAATTAATTATTTATTTATTTATTTATTTATTTATTTATTT"
            .as_bytes()
            .to_vec();
        //let frg = "TTATATTTATTTATATATATTTATATAGTTTATATATATATTTATATATAAATATATA".as_bytes().to_vec();
        let m = match_reads(&base_frg, &frg, true, 0.1, 0, 0, 32);
        if let Some(m) = m {
            let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
            let aln_segs = deltas_to_aln_segs(&deltas, m.end0 as usize, m.end1 as usize, &frg);
            let re_seq = reconstruct_seq_from_aln_segs(&base_frg, &aln_segs);
            if frg != re_seq || true {
                println!("{} {}", String::from_utf8_lossy(&base_frg), base_frg.len());
                println!("{} {}", String::from_utf8_lossy(&frg), frg.len());
                println!("{} {} {} {}", m.bgn0, m.end0, m.bgn1, m.end1);
                println!("{:?}", deltas);
                println!(
                    "{}",
                    String::from_utf8_lossy(&reconstruct_seq_from_aln_segs(&base_frg, &aln_segs))
                );
                println!("{:?}", aln_segs);
            }
            assert_eq!(frg, reconstruct_seq_from_aln_segs(&base_frg, &aln_segs));
        }
    }

    #[test]
    fn reconstruct_test2() {
        let base_frg = "TATTTATATTTATTTATATATATTTATATATTTATATATATATTTATATATAAATAT"
            .as_bytes()
            .to_vec();
        let frg = "TTTTTTATTTTTTTAATTAATTAATTATTTATTTATTTATTTATTTATTTATTTATT"
            .as_bytes()
            .to_vec();
        //let frg = "TTATATTTATTTATATATATTTATATAGTTTATATATATATTTATATATAAATATATA".as_bytes().to_vec();
        let m = match_reads(&base_frg, &frg, true, 0.1, 0, 0, 32);
        if let Some(m) = m {
            let deltas: Vec<DeltaPoint> = m.deltas.unwrap();
            let aln_segs = deltas_to_aln_segs(&deltas, m.end0 as usize, m.end1 as usize, &frg);
            let re_seq = reconstruct_seq_from_aln_segs(&base_frg, &aln_segs);
            if frg != re_seq || true {
                println!("{} {}", String::from_utf8_lossy(&base_frg), base_frg.len());
                println!("{} {}", String::from_utf8_lossy(&frg), frg.len());
                println!("{} {} {} {}", m.bgn0, m.end0, m.bgn1, m.end1);
                println!("{:?}", deltas);
                println!(
                    "{}",
                    String::from_utf8_lossy(&reconstruct_seq_from_aln_segs(&base_frg, &aln_segs))
                );
                println!("{:?}", aln_segs);
            }
            assert_eq!(frg, reconstruct_seq_from_aln_segs(&base_frg, &aln_segs));
        }
    }

    #[test]
    fn rc_match() {
        use seq_db::GetSeq;
        let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
        let _ = sdb.load_seqs_from_fastx("test/test_data/test_rev.fa".to_string());
        let cs0 = sdb.get_seq_by_id(0);
        let cs1 = sdb.get_seq_by_id(1);
        let shmmr_spec = seq_db::SHMMRSPEC;
        let shmmr0 = shmmrutils::sequence_to_shmmrs(0, &cs0, &shmmr_spec, false);
        let shmmr1 = shmmrutils::sequence_to_shmmrs(0, &cs1, &shmmr_spec, false);
        let shmmr0 = shmmr0.iter().map(|m| m.hash()).collect::<Vec<u64>>();
        let shmmr1 = shmmr1.iter().rev().map(|m| m.hash()).collect::<Vec<u64>>();
        assert!(shmmr0.len() > 0);
        assert_eq!(shmmr0, shmmr1);
    }






    #[test]
    fn test_shmmr_reduction_boundary_condition() {
        // In some cases, the numnber of shimmers might be less than the number of the reducetion factor.
        // In such cases, we want to output the first and the last shimmers as if there are other bigger one surrounding them.
        let seq = b"CCAGTTGTATCCATGACAAAGATGAGGCCGCGAGGAGGGCGAGTGGGTTTGGGGGCAGGCAGAGTGCCTTGGAGAACTTACAGGTCCTGCCACAATCCTAATGCAAGGATGGAGCTGCAAGTTCAGTTTGGGAATCATCAGCCTGGATTGGTTTGGTGGAAGCCAGGGAGTGGTTGAGGACCCCCACAGGGGAGCTCTGAGGAAGGAAGTTCCGAAGGAGGGAACGTAAGAAATGACCAGGTCAGAACCAAGGGTGGTCCAGAAGCTAACCCTTAGCTTAGGGACAGTTTCACAGAGAACACGTCCATGATGCAAGACTCTGCTGAGGGCCTGGAGCAGTGAAGACTGGGGCAAGGTCACCCTCTGGGAAGTGAAGTCACCAGAGACCTTGCGGAGCAGCTTTGAGAGTTCTCTGAGTAGGAAGGTAACAGAATGTGAAGGACACTGGAGAGAAGGCCAATAGGAAGCAAACAAAAACAGGCCAAGGAAACCCAGTACAGGGGGCTGCAGGGCCCAGGGAGTGGGTCCCTCATCTCTCCTCCCCACGCTTGGCCAGGTCCCCACCTCCCCCGGGAGTGCGTGGGCTTTGAGGCTGTGCAGGAAGTGCCGGTGGGGCTGGTGCAGCCGGCCAGCGCAACCCTGTACGACTACTACAACCCCGGTGAGCACTGCAGGACACCCTGAAATTCAGGAGAACTTTGGCATAGGTGCCCTCCTATGGGACAATGGACACCGGGGTAGTGAGGGGGCAGAGAGCCCTGGGGCTCCCTGGGACTGAGGAGGCAGAATGGAGGGGCCTGTGCCCTAACTCCTCTCTGTTCTCCAGAGCGCAGATGTTCTGTGTTTTACGGGGCACCAAGTAAGAGCAGACTCTTGGCCACCTTGTGTTCTGCTGAAGTCTGCCAGTGTGCTGAGGGTGAGACTGAGGGCCTGGGGCGGGGCAGT";
        let seq2 = b"CCAGTTGTATCCATGACAAAGATGAGGCCGCGAGGAGGGCGAGTGGGTTTGGGGGCAGGCAGAGTGCCTTGGAGAACTTACAGGTCCTGCCACAATCCTAATGCAAGGATGGAGCTGCAAGTTCAGTTTGGGAATCATCAGCCTGGATTGGTTTGGTGGAAGCCAGGGAGTGGTTGAGACCCCCACAGGGGAGCTCTGAGGAAGGAAGTTCCGAAGGAGGGAACGTAAGAAATGACCAGGTCAGAACCAAGGGTGGTCCAGAAGCTAACCCTTAGCTTAGGGACAGTTTCACAGAGAACACGTCCATGATGCAAGACTCTGCTGAGGGCCTGGAGCAGTGAAGACTGGGGCAAGGTCACCCTCTGGGAAGTGAAGTCACCAGAGACCTTGCGGAGCAGCTTTGAGAGTTCTCTGAGTAGGAAGGTAACAGAATGTGAAGGACACTGGAGAGAAGGCCAATAGGAAGCAAACAAAAACAGGCCAAGGAAACCCAGTACAGGGGGCTGCAGGGCCCAGGGAGTGGGTCCCTCATCTCTCCTCCCCACGCTTGGCCAGGTCCCCACCTCCCGGGAGTGCGTGGGCTTTGAGGCTGTGCAGGAAGTGCCGGTGGGGCTGGTGCAGCCGGCCAGCGCAACCCTGTACGACTACTACAACCCCGGTGAGCACTGCAGGACACCCTGAAATTCAGGAGAACTTTGGCATAGGTGCCCTCCTATGGGACAATGGACACCGGGGTAGTGAGGGGGCAGAGAGCCCTGGGGCTCCCTGGGACTGAGGAGGCAGAATGGAGGGGCCTGTGCCCTAACTCCTCTCTGTTCTCCAGAGCGCAGATGTTCTGTGTTTTACGGGGCACCAAGTAAGAGCAGACTCTTGGCCACCTTGTGTTCTGCTGAAGTCTGCCAGTGTGCTGAGGGTGAGACTGAGGGCCTGGGGCGGGGCAGT";
        use shmmrutils::sequence_to_shmmrs;
        use shmmrutils::ShmmrSpec;
        let spec = ShmmrSpec {
            w: 24,
            k: 24,
            r: 12,
            min_span: 24,
            sketch: false,
        };
        let out1 = sequence_to_shmmrs(0, &seq.to_vec(), &spec, true);
        println!("out1: {} {:?}", out1.len(), out1);
        let out2 = sequence_to_shmmrs(0, &seq2.to_vec(), &spec, true);
        println!("out2: {} {:?}", out2.len(), out2);
        assert!(out1.len() == 2);
        assert!(out2.len() == 2);
    }

    #[test]
    fn test_open_compact_seq_db_storage() {
        use seq_db::GetSeq;
        use crate::frag_file_io::CompactSeqDBStorage;
        let seq_storage = CompactSeqDBStorage::new("test/test_data/test_seqs_frag".to_string());
        let seq = seq_storage.get_seq_by_id(0);
        println!("{}", String::from_utf8_lossy(&seq[..]));
        let seq = seq_storage.get_sub_seq_by_id(0, 100, 200);
        println!("{}", String::from_utf8_lossy(&seq[..]));
    }

    #[test]
    fn test_seq_db_storage_get_sub_read() {
        use crate::frag_file_io::CompactSeqDBStorage;
        use seq_db::GetSeq;
        let seq_storage = CompactSeqDBStorage::new("test/test_data/test_seqs_frag".to_string());
        let sid = 0;

        let seq = seq_storage.get_seq_by_id(sid);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 0, 105);
        assert_eq!(seq[0..105], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 105, 286);
        assert_eq!(seq[105..286], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 104, 286);
        assert_eq!(seq[104..286], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 105, 287);
        assert_eq!(seq[105..287], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 250, 1423);
        assert_eq!(seq[250..1423], sub_seq[..]);

        let sid = 5;

        let seq = seq_storage.get_seq_by_id(sid);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 0, 105);
        assert_eq!(seq[0..105], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 105, 286);
        assert_eq!(seq[105..286], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 104, 286);
        assert_eq!(seq[104..286], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 105, 287);
        assert_eq!(seq[105..287], sub_seq[..]);

        let sub_seq = seq_storage.get_sub_seq_by_id(sid, 250, 1423);
        assert_eq!(seq[250..1423], sub_seq[..]);
    }

    #[test]
    fn test_seq_db_get_sub_read() {
        use seq_db::GetSeq;
        let mut sdb = seq_db::CompactSeqDB::new(seq_db::SHMMRSPEC);
        let _ = sdb.load_seqs_from_fastx("test/test_data/test_seqs.fa".to_string());
        let sid = 0;
        let seq = sdb.get_seq_by_id(sid);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 0, 105);
        assert_eq!(seq[0..105], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 105, 286);
        assert_eq!(seq[105..286], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 104, 286);
        assert_eq!(seq[104..286], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 105, 287);
        assert_eq!(seq[105..287], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 250, 1423);
        assert_eq!(seq[250..1423], sub_seq[..]);

        let sid = 5;
        let seq = sdb.get_seq_by_id(sid);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 0, 105);
        assert_eq!(seq[0..105], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 105, 286);
        assert_eq!(seq[105..286], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 104, 286);
        assert_eq!(seq[104..286], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 105, 287);
        assert_eq!(seq[105..287], sub_seq[..]);

        let sub_seq = sdb.get_sub_seq_by_id(sid, 250, 1423);
        assert_eq!(seq[250..1423], sub_seq[..]);
    }
}
