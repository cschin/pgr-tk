use agc_rs::compressor::Compressor;
use agc_rs::decompressor::{bits_to_ascii, AgcFile};
use agc_rs::fasta_io::read_fasta_gz;
use agc_rs::merge::merge_archives;
use agc_rs::segment::Params;
use std::collections::HashMap;
use std::path::Path;
use tempfile::tempdir;

const MG1655: &str = "test_data/ecoli/ecoli_k12_mg1655.fna.gz";
const W3110: &str = "test_data/ecoli/ecoli_k12_w3110.fna.gz";
const SAKAI: &str = "test_data/ecoli/ecoli_o157h7_sakai.fna.gz";

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

/// Assert every contig in `expected_fasta` round-trips correctly from `agc`.
fn assert_sequences_match(agc: &AgcFile, sample: &str, expected_fasta: &str) {
    let expected = load_fasta_as_map(expected_fasta);
    for (contig, seq) in &expected {
        let got = agc
            .full_contig(sample, contig)
            .unwrap_or_else(|e| panic!("full_contig({sample}/{contig}) failed: {e}"));
        assert_eq!(
            got,
            *seq,
            "sequence mismatch: {sample}/{contig} (first diff at {:?})",
            got.iter().zip(seq.iter()).position(|(a, b)| a != b)
        );
    }
}

#[test]
fn ecoli_three_genome_round_trip() {
    let mg1655 = MG1655;
    let w3110 = W3110;
    let sakai = SAKAI;

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
    let mg1655 = MG1655;
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
    let mg1655 = MG1655;
    let w3110 = W3110;
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
    let sakai = SAKAI;
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

// ---------------------------------------------------------------------------
// batch-append tests
// ---------------------------------------------------------------------------

/// batch-append W3110 + Sakai in one call and verify every sequence round-trips
/// identically to the result of sequential add_fasta calls.
#[test]
fn batch_append_matches_sequential_append() {
    if !Path::new(MG1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let sequential = dir.path().join("sequential.agcrs");
    let batched = dir.path().join("batched.agcrs");

    // --- sequential reference ------------------------------------------------
    {
        let mut c = Compressor::create(&sequential, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
        c.add_fasta(Path::new(W3110), "W3110").unwrap();
        c.add_fasta(Path::new(SAKAI), "Sakai").unwrap();
        c.finish().unwrap();
    }

    // --- batch append --------------------------------------------------------
    {
        let mut c = Compressor::create(&batched, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap(); // reference
        c.finish().unwrap();
    }
    {
        let mut c = Compressor::append(&batched).unwrap();
        c.add_fasta_batch(&[(Path::new(W3110), "W3110"), (Path::new(SAKAI), "Sakai")])
            .unwrap();
        c.finish().unwrap();
    }

    // --- compare -------------------------------------------------------------
    let ref_agc = AgcFile::open(&sequential).unwrap();
    let test_agc = AgcFile::open(&batched).unwrap();

    assert_eq!(test_agc.n_samples().unwrap(), ref_agc.n_samples().unwrap());

    for (sample, path) in [("MG1655", MG1655), ("W3110", W3110), ("Sakai", SAKAI)] {
        assert_sequences_match(&test_agc, sample, path);
        // Also verify against the sequential reference to catch subtle divergence.
        let expected = load_fasta_as_map(path);
        for (contig, seq) in &expected {
            let from_ref = ref_agc.full_contig(sample, contig).unwrap();
            let from_test = test_agc.full_contig(sample, contig).unwrap();
            assert_eq!(
                from_ref, from_test,
                "batch vs sequential mismatch: {sample}/{contig}"
            );
            assert_eq!(
                from_test, *seq,
                "batch vs FASTA mismatch: {sample}/{contig}"
            );
        }
    }
}

/// batch-append to an empty archive must fail with a clear error.
#[test]
fn batch_append_empty_archive_returns_error() {
    if !Path::new(W3110).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("empty.agcrs");

    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.finish().unwrap();

    let mut c = Compressor::append(&archive).unwrap();
    let result = c.add_fasta_batch(&[(Path::new(W3110), "W3110")]);
    assert!(
        result.is_err(),
        "expected error when batch-appending to archive with no reference"
    );
}

/// batch-append must reject a duplicate sample name before doing any work.
#[test]
fn batch_append_duplicate_sample_returns_error() {
    if !Path::new(MG1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("dup.agcrs");

    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
    c.add_fasta(Path::new(W3110), "W3110").unwrap();
    c.finish().unwrap();

    // Try to batch-append "W3110" again.
    let mut c = Compressor::append(&archive).unwrap();
    let result = c.add_fasta_batch(&[(Path::new(SAKAI), "W3110")]);
    assert!(result.is_err(), "expected error for duplicate sample name");
}

// ---------------------------------------------------------------------------
// merge tests
// ---------------------------------------------------------------------------

/// Merging two archives that each hold ref + one haplotype produces an archive
/// with all three samples, and all sequences round-trip correctly.
#[test]
fn merge_two_archives_round_trip() {
    if !Path::new(MG1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive_a = dir.path().join("a.agcrs");
    let archive_b = dir.path().join("b.agcrs");
    let merged = dir.path().join("merged.agcrs");

    // archive_a: MG1655 (ref) + W3110 (delta)
    {
        let mut c = Compressor::create(&archive_a, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
        c.add_fasta(Path::new(W3110), "W3110").unwrap();
        c.finish().unwrap();
    }

    // archive_b: MG1655 (ref) + Sakai (delta)
    {
        let mut c = Compressor::create(&archive_b, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
        c.add_fasta(Path::new(SAKAI), "Sakai").unwrap();
        c.finish().unwrap();
    }

    merge_archives(&merged, &[&archive_a, &archive_b]).unwrap();

    let agc = AgcFile::open(&merged).unwrap();
    // MG1655 (from archive_a) + W3110 + Sakai
    assert_eq!(agc.n_samples().unwrap(), 3);

    let names = agc.list_samples().unwrap();
    assert!(names.contains(&"MG1655".to_string()));
    assert!(names.contains(&"W3110".to_string()));
    assert!(names.contains(&"Sakai".to_string()));

    assert_sequences_match(&agc, "MG1655", MG1655);
    assert_sequences_match(&agc, "W3110", W3110);
    assert_sequences_match(&agc, "Sakai", SAKAI);
}

/// Merging archives built from different references must fail.
#[test]
fn merge_different_references_returns_error() {
    if !Path::new(MG1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive_mg = dir.path().join("mg.agcrs");
    let archive_w3 = dir.path().join("w3.agcrs");
    let merged = dir.path().join("bad_merge.agcrs");

    // Two archives with different reference genomes → different splitter sets.
    {
        let mut c = Compressor::create(&archive_mg, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
        c.finish().unwrap();
    }
    {
        let mut c = Compressor::create(&archive_w3, Params::default()).unwrap();
        c.add_fasta(Path::new(W3110), "W3110").unwrap();
        c.finish().unwrap();
    }

    let result = merge_archives(&merged, &[&archive_mg, &archive_w3]);
    assert!(
        result.is_err(),
        "expected error when merging archives with different reference splitters"
    );
}

/// Merging archives with overlapping non-reference sample names must fail.
#[test]
fn merge_duplicate_sample_names_returns_error() {
    if !Path::new(MG1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive_a = dir.path().join("a_dup.agcrs");
    let archive_b = dir.path().join("b_dup.agcrs");
    let merged = dir.path().join("dup_merge.agcrs");

    // Both archives have the same non-reference sample name "W3110".
    {
        let mut c = Compressor::create(&archive_a, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
        c.add_fasta(Path::new(W3110), "W3110").unwrap();
        c.finish().unwrap();
    }
    {
        let mut c = Compressor::create(&archive_b, Params::default()).unwrap();
        c.add_fasta(Path::new(MG1655), "MG1655").unwrap();
        c.add_fasta(Path::new(W3110), "W3110").unwrap(); // same name
        c.finish().unwrap();
    }

    let result = merge_archives(&merged, &[&archive_a, &archive_b]);
    assert!(
        result.is_err(),
        "expected error when merging archives with duplicate non-reference sample names"
    );
}
