use agc_rs::compressor::Compressor;
use agc_rs::segment::Params;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn metadata_readable_without_agc_rs_api() {
    let mg1655 = "test_data/ecoli/ecoli_k12_mg1655.fna.gz";
    if !Path::new(mg1655).exists() {
        return;
    }

    let dir = tempdir().unwrap();
    let archive = dir.path().join("interop.agcrs");
    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new(mg1655), "MG1655").unwrap();
    c.finish().unwrap();

    // Open with raw rusqlite — no agc-rs API — to verify schema is plain SQLite.
    let conn = rusqlite::Connection::open(&archive).unwrap();

    let n: i64 = conn
        .query_row("SELECT COUNT(*) FROM sample", [], |r| r.get(0))
        .unwrap();
    assert_eq!(n, 1);

    let name: String = conn
        .query_row("SELECT name FROM sample", [], |r| r.get(0))
        .unwrap();
    assert_eq!(name, "MG1655");

    let nc: i64 = conn
        .query_row("SELECT COUNT(*) FROM contig", [], |r| r.get(0))
        .unwrap();
    assert!(nc >= 1, "expected at least one contig, got {nc}");
}
