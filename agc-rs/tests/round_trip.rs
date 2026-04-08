use agc_rs::compressor::Compressor;
use agc_rs::decompressor::{bits_to_ascii, AgcFile};
use agc_rs::fasta_io::read_fasta_gz;
use agc_rs::segment::Params;
use std::collections::HashMap;
use std::path::Path;
use tempfile::tempdir;

/// Read a FASTA (plain or .gz) file and return a name -> ASCII-sequence map.
/// Returns an empty map (rather than panicking) if the file does not exist,
/// so that tests skip gracefully when test data are absent.
fn load_fasta_as_map(path: &str) -> HashMap<String, Vec<u8>> {
    if !Path::new(path).exists() {
        return HashMap::new();
    }
    let records = read_fasta_gz(Path::new(path)).expect("read_fasta_gz");
    records
        .into_iter()
        .map(|r| (r.name, bits_to_ascii(&r.seq)))
        .collect()
}

#[test]
fn ecoli_three_genome_round_trip() {
    let mg1655 = "test_data/ecoli/ecoli_k12_mg1655.fna.gz";
    let w3110 = "test_data/ecoli/ecoli_k12_w3110.fna.gz";
    let sakai = "test_data/ecoli/ecoli_o157h7_sakai.fna.gz";

    if !Path::new(mg1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("ecoli.agcrs");

    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new(mg1655), "MG1655").unwrap();
    c.add_fasta(Path::new(w3110), "W3110").unwrap();
    c.add_fasta(Path::new(sakai), "Sakai").unwrap();
    c.finish().unwrap();

    let agc = AgcFile::open(&archive).unwrap();
    assert_eq!(agc.n_samples().unwrap(), 3);

    for (sample, path) in [("MG1655", mg1655), ("W3110", w3110), ("Sakai", sakai)] {
        let expected = load_fasta_as_map(path);
        for (contig, seq) in &expected {
            let got = agc.full_contig(sample, contig).unwrap();
            assert_eq!(
                got,
                *seq,
                "mismatch: {sample}/{contig} (first diff at {:?})",
                got.iter().zip(seq.iter()).position(|(a, b)| a != b)
            );
        }
    }
}

#[test]
fn ecoli_subrange_queries() {
    let mg1655 = "test_data/ecoli/ecoli_k12_mg1655.fna.gz";
    if !Path::new(mg1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("ecoli_sub.agcrs");
    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new(mg1655), "MG1655").unwrap();
    c.finish().unwrap();

    let agc = AgcFile::open(&archive).unwrap();
    let contigs = agc.list_contigs("MG1655").unwrap();
    let contig = &contigs[0];
    let full = agc.full_contig("MG1655", contig).unwrap();
    let len = full.len() as u64;

    for (s, e) in [
        (0u64, 1000u64),
        (500, 2000),
        (100_000, 200_000),
        (len - 1000, len),
        (0, len),
    ] {
        let sub = agc.contig_seq("MG1655", contig, s, e).unwrap();
        assert_eq!(
            sub,
            full[s as usize..e as usize],
            "mismatch at range {s}..{e}"
        );
    }
}

#[test]
fn ecoli_append_sample() {
    let mg1655 = "test_data/ecoli/ecoli_k12_mg1655.fna.gz";
    let w3110 = "test_data/ecoli/ecoli_k12_w3110.fna.gz";
    if !Path::new(mg1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("ecoli_append.agcrs");

    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new(mg1655), "MG1655").unwrap();
    c.finish().unwrap();

    let mut c = Compressor::append(&archive).unwrap();
    c.add_fasta(Path::new(w3110), "W3110").unwrap();
    c.finish().unwrap();

    let agc = AgcFile::open(&archive).unwrap();
    assert_eq!(agc.n_samples().unwrap(), 2);
    assert!(agc.list_samples().unwrap().contains(&"MG1655".to_string()));
    assert!(agc.list_samples().unwrap().contains(&"W3110".to_string()));
}

#[test]
fn ecoli_sakai_three_sequences() {
    let sakai = "test_data/ecoli/ecoli_o157h7_sakai.fna.gz";
    if !Path::new(sakai).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("sakai.agcrs");
    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new(sakai), "Sakai").unwrap();
    c.finish().unwrap();

    let agc = AgcFile::open(&archive).unwrap();
    // Sakai has chromosome + 2 plasmids = 3 sequences.
    assert_eq!(agc.n_contigs("Sakai").unwrap(), 3);
}
