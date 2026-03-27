# Plan: Porting AGC to Rust

## Why

AGC is integrated into pgr-tk via a C FFI wrapper calling a pre-built static library.
Problems: platform-specific `.a` files, panic-prone `CString::new(...).unwrap()` at the
FFI boundary, a C++ build toolchain requirement, and a binary format no other tool can
read without the AGC binary.

**The port does not reproduce the AGC file format.** SQLite is used as the container —
readable from Python, R, DuckDB, and the shell with no custom code. The LZ-diff + ZSTD
compression algorithm is preserved because it delivers the compression ratio. Everything
else (custom archive format, thread pool, hash set, pre-built static libraries) is
replaced by idiomatic Rust.

---

## Container: SQLite

```sql
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
-- keys: 'schema_version', 'kmer_len', 'min_match_len', 'segment_size', 'created_at'

CREATE TABLE sample (id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL);

CREATE TABLE contig (id INTEGER PRIMARY KEY,
                     sample_id INTEGER NOT NULL REFERENCES sample(id),
                     name TEXT NOT NULL, length INTEGER NOT NULL,
                     UNIQUE(sample_id, name));

CREATE TABLE segment_group (id INTEGER PRIMARY KEY,
                             ref_data BLOB NOT NULL,  -- ZSTD(reference sequence)
                             params   TEXT NOT NULL); -- JSON: key_len, lz_version, …

CREATE TABLE segment (id INTEGER PRIMARY KEY,
                      contig_id   INTEGER NOT NULL REFERENCES contig(id),
                      seg_order   INTEGER NOT NULL,
                      group_id    INTEGER NOT NULL REFERENCES segment_group(id),
                      in_group_id INTEGER NOT NULL,  -- 0 = reference of this group
                      is_rev_comp BOOLEAN NOT NULL,
                      raw_length  INTEGER NOT NULL,
                      delta_data  BLOB,              -- ZSTD(LZ-diff output); NULL for ref
                      UNIQUE(contig_id, seg_order));
CREATE INDEX idx_seg ON segment(contig_id, seg_order);
```

Readable immediately without agc-rs:
```python
import sqlite3
con = sqlite3.connect("archive.agcrs")
print(con.execute("SELECT name, length FROM contig").fetchall())
```
```bash
sqlite3 archive.agcrs "SELECT sm.name, c.name, c.length FROM sample sm JOIN contig c ON c.sample_id=sm.id"
```

---

## Crate Layout

```
agc-rs/
  src/
    lib.rs           # public AgcFile API
    error.rs         # AgcError, Result<T>
    db.rs            # SQLite schema, open/create/readonly
    kmer.rs          # Kmer (2-bit u64), push, canonical, rev_comp
    lz_diff.rs       # LZ-diff encode + decode
    segment.rs       # ZSTD(LZ-diff) compress/decompress per segment group
    compressor.rs    # full pipeline: FASTA → segment groups → SQLite
    decompressor.rs  # read-only AgcFile query API
    fasta_io.rs      # reuse pgr-db/src/fasta_io.rs
    bin/agc.rs       # CLI (clap): create, append, get, list, info
  tests/
    round_trip.rs    # compress then decompress, compare to input FASTA
    subrange.rs      # contig_seq(s, e) == full_seq[s..e]
    interop.rs       # schema readable via raw SQLite (no agc-rs API)
```

**Cargo.toml** — no C++ toolchain, no pre-built `.a`:
```toml
[dependencies]
rusqlite   = { version = "0.31", features = ["bundled"] }
zstd       = "0.13"
flate2     = "1"        # gzip FASTA
rayon      = "1"        # parallel segment compression
rustc-hash = "1"        # FxHashMap for k-mer tables
thiserror  = "1"
clap       = { version = "4", features = ["derive"] }
serde_json = "1"        # params column in segment_group
```

---

## Test Data

Three *E. coli* genomes are stored in `test_data/ecoli/`:

| File | Strain | Accession | Role |
|------|--------|-----------|------|
| `ecoli_k12_mg1655.fna.gz` | K-12 MG1655 | NC_000913.3 | **reference** |
| `ecoli_k12_w3110.fna.gz` | K-12 W3110 | NC_007779.1 | query 1 |
| `ecoli_o157h7_sakai.fna.gz` | O157:H7 Sakai | NC_002695.2 + plasmids | query 2 |

These are real complete genome assemblies (~4.6 MB each uncompressed). K-12 MG1655 and
W3110 share ~99% identity; O157:H7 Sakai is more divergent and also carries two plasmids,
testing multi-sequence FASTA handling.

---

## Modules and Tests

### `error.rs`

`AgcError` enum with variants for DB, I/O, ZSTD, sample-not-found, contig-not-found,
range-out-of-bounds, LZ-diff failure. `Result<T>` alias. `From` impls for
`rusqlite::Error` and `std::io::Error`.

---

### `db.rs`

Schema creation, `open` / `create` / `open_readonly`. WAL mode. `schema_version` check
on open (reject incompatible versions with a clear error).

**Tests:**
```rust
#[test] fn create_and_reopen() {
    // create in tempfile, insert one meta row, close, reopen, verify row
}
#[test] fn wrong_schema_version_is_rejected() {
    // write wrong version into meta, expect Err(AgcError::UnsupportedVersion)
}
#[test] fn readonly_rejects_writes() {
    // open_readonly, attempt INSERT, expect Err
}
```

---

### `kmer.rs`

`Kmer` struct: packed 2-bit u64 pair (forward + reverse complement), current length.
`push(base)`, `canonical() -> u64`, `is_forward_canonical() -> bool`.
Reverse complement via bitwise inversion + shift.

**Tests:**
```rust
#[test] fn canonical_is_min_of_fwd_and_rc() { /* ACGT k=4 */ }
#[test] fn push_wraps_correctly_at_max_k()  { /* k=32 boundary */ }
#[test] fn reverse_complement_known_seqs()  { /* AATTCG → CGAATT */ }
```

---

### `lz_diff.rs`

Port of `CLZDiff_V1` from the C++ AGC source. Reference-based LZ encoder/decoder.

**`prepare(reference)`** — index every 4th k-mer (length `key_len = min_match_len - 4`)
into a `FxHashMap<u64, Vec<u32>>` mapping kmer_hash → positions in reference.

**`encode(query) -> Vec<u8>`** — greedy match finding: scan query, look up k-mer in
table, extend match backward and forward, emit `Match(pos, len)` or `Literal(base)`.
Special token for N-runs. Serialise as compact byte stream (LEB128 for positions/lengths).

**`decode(encoded) -> Result<Vec<u8>>`** — reconstruct query from literal/match token
stream and the reference.

**Tests:**
```rust
#[test] fn round_trip_identical_seqs()     { /* query == ref */ }
#[test] fn round_trip_single_snp()         { /* one base differs */ }
#[test] fn round_trip_insertion()          { /* query longer */ }
#[test] fn round_trip_deletion()           { /* query shorter */ }
#[test] fn round_trip_no_shared_kmers()    { /* completely different sequences */ }
#[test] fn round_trip_n_run()              { /* 500-base N stretch */ }
#[test] fn round_trip_ecoli_mg1655_vs_w3110() {
    // load first 100 KB of each genome, encode W3110 against MG1655, decode, compare
    let ref_seq  = load_first_n("test_data/ecoli/ecoli_k12_mg1655.fna.gz", 100_000);
    let query    = load_first_n("test_data/ecoli/ecoli_k12_w3110.fna.gz",  100_000);
    let mut lz   = LzDiff::new(16);
    lz.prepare(&ref_seq);
    let encoded  = lz.encode(&query);
    let decoded  = lz.decode(&encoded).unwrap();
    assert_eq!(decoded, query);
}
```

---

### `segment.rs`

Per-segment-group compress/decompress. The reference sequence is stored as `ZSTD(ref)`.
Each other sequence is stored as `ZSTD(lz_diff.encode(query))`. Decompression reverses
both steps. Subrange reconstruction accumulates `raw_length` offsets across segments
and decompresses only the segments that overlap `[start, end)`.

**Tests:**
```rust
#[test] fn compress_decompress_group_of_three()  { /* ref + 2 queries, recover all */ }
#[test] fn subrange_matches_full_slice()          { /* [500..1500] == full[500..1500] */ }
#[test] fn single_sequence_group()               { /* only a reference, no deltas */ }
```

---

### `decompressor.rs`

`AgcFile` public API: `open()`, `open_readonly()`, `list_samples()`, `list_contigs()`,
`contig_len()`, `contig_seq(sample, contig, start, end)`, `full_contig()`.
All queries go through prepared SQLite statements cached on the connection.

**Tests (use in-memory SQLite with synthetic rows):**
```rust
#[test] fn list_samples_empty_db()     { }
#[test] fn contig_not_found_is_err()   { }
#[test] fn contig_seq_full_sequence()  { }
#[test] fn contig_seq_subrange()       { }
#[test] fn contig_len_correct()        { }
```

---

### `compressor.rs`

FASTA reading → k-mer-based segmentation at `~segment_size` boundaries → reference
selection (first genome) → parallel `rayon::par_iter()` over segment groups → insert
into SQLite in one transaction per sample.

Append mode: open existing archive, add new samples without disturbing existing data.

---

### `bin/agc.rs`

```
agc-rs create  -k 31 -l 20 -s 60000 -t 8 -o out.agcrs ref.fa [query.fa...]
agc-rs append  out.agcrs [query.fa...]
agc-rs get     out.agcrs sample/contig[:start-end]
agc-rs getset  out.agcrs sample1 sample2 ... -o output/
agc-rs list    out.agcrs [--samples | --contigs sample]
agc-rs info    out.agcrs
```

---

## Integration Tests (`tests/round_trip.rs`)

```rust
// helper: load all sequences from a gzipped FASTA into a HashMap<name, Vec<u8>>
fn load_fasta_gz(path: &str) -> HashMap<String, Vec<u8>> { ... }

#[test]
fn ecoli_three_genome_round_trip() {
    let dir = tempdir().unwrap();
    let archive = dir.path().join("ecoli.agcrs");

    // compress: MG1655 as reference, W3110 and O157:H7 Sakai as queries
    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new("test_data/ecoli/ecoli_k12_mg1655.fna.gz"), "MG1655").unwrap();
    c.add_fasta(Path::new("test_data/ecoli/ecoli_k12_w3110.fna.gz"),  "W3110").unwrap();
    c.add_fasta(Path::new("test_data/ecoli/ecoli_o157h7_sakai.fna.gz"), "Sakai").unwrap();
    c.finish().unwrap();

    // decompress and compare every contig byte-for-byte
    let agc = AgcFile::open(&archive).unwrap();
    for (sample, path) in [
        ("MG1655", "test_data/ecoli/ecoli_k12_mg1655.fna.gz"),
        ("W3110",  "test_data/ecoli/ecoli_k12_w3110.fna.gz"),
        ("Sakai",  "test_data/ecoli/ecoli_o157h7_sakai.fna.gz"),
    ] {
        let expected = load_fasta_gz(path);
        for (contig, seq) in &expected {
            let got = agc.full_contig(sample, contig).unwrap();
            assert_eq!(got, *seq, "mismatch: {sample}/{contig}");
        }
    }
}

#[test]
fn ecoli_subrange_queries() {
    // after compression, verify contig_seq(s, e) == full[s..e] for
    // 20 random (start, end) pairs across all three genomes
    let dir = tempdir().unwrap();
    let archive = dir.path().join("ecoli.agcrs");
    // ... compress ...
    let agc = AgcFile::open(&archive).unwrap();
    let full = agc.full_contig("MG1655", "NC_000913.3").unwrap();
    for (s, e) in [(0,1000), (500,1500), (100_000, 200_000),
                   (4_000_000, 4_100_000), (0, full.len() as u64)] {
        let sub = agc.contig_seq("MG1655", "NC_000913.3", s, e).unwrap();
        assert_eq!(sub, full[s as usize..e as usize]);
    }
}

#[test]
fn ecoli_append_sample() {
    let dir = tempdir().unwrap();
    let archive = dir.path().join("ecoli.agcrs");
    // create with only MG1655
    let mut c = Compressor::create(&archive, Params::default()).unwrap();
    c.add_fasta(Path::new("test_data/ecoli/ecoli_k12_mg1655.fna.gz"), "MG1655").unwrap();
    c.finish().unwrap();
    // append W3110
    let mut c = Compressor::append(&archive).unwrap();
    c.add_fasta(Path::new("test_data/ecoli/ecoli_k12_w3110.fna.gz"), "W3110").unwrap();
    c.finish().unwrap();
    // both samples present and correct
    let agc = AgcFile::open(&archive).unwrap();
    assert_eq!(agc.n_samples().unwrap(), 2);
}

#[test]
fn ecoli_sakai_plasmids_preserved() {
    // O157:H7 Sakai has 3 sequences (chromosome + 2 plasmids)
    // verify all three are stored and recovered correctly
    let agc = /* open compressed archive */;
    let contigs = agc.list_contigs("Sakai").unwrap();
    assert_eq!(contigs.len(), 3);
}
```

## Interop Test (`tests/interop.rs`)

```rust
#[test]
fn metadata_readable_without_agc_rs() {
    // open the .agcrs file with a raw rusqlite::Connection (no agc-rs API at all)
    // verify sample and contig tables are populated correctly
    let conn = rusqlite::Connection::open(&archive).unwrap();
    let n: i64 = conn.query_row(
        "SELECT COUNT(*) FROM sample", [], |r| r.get(0)).unwrap();
    assert_eq!(n, 3);
    let names: Vec<String> = conn.prepare("SELECT name FROM sample ORDER BY name")
        .unwrap().query_map([], |r| r.get(0)).unwrap()
        .map(|r| r.unwrap()).collect();
    assert_eq!(names, ["MG1655", "Sakai", "W3110"]);
}
```

---

## Updating pgr-db

Once `agc-rs` is complete, replace the C FFI in `pgr-db/src/agc_io.rs`:

```rust
// before: extern "C" { fn agc_open(...) -> *mut c_void; ... }

// after:
use agc_rs::AgcFile;

pub struct AgcDB { inner: AgcFile }

impl AgcDB {
    pub fn open(path: &str) -> crate::Result<Self> {
        AgcFile::open_readonly(Path::new(path))
            .map(|f| AgcDB { inner: f })
            .map_err(crate::PgrError::from)
    }
    pub fn get_seq(&self, sample: &str, contig: &str,
                   start: u64, end: u64) -> crate::Result<Vec<u8>> {
        self.inner.contig_seq(sample, contig, start, end)
            .map_err(crate::PgrError::from)
    }
}
```

Remove the `agc` git submodule and any `build.rs` C++ compilation steps.
The `py_agc_test.py` tests are adapted to the new format and run against the Rust backend.

---

## What Is Dropped

| C++ component | Replacement |
|---------------|-------------|
| Custom `CArchive` binary format | SQLite |
| `CCollection` metadata serialisation | SQL schema |
| `hash_set_lp` (linear-probe hash set) | `FxHashMap` |
| `CBoundedQueue` + manual thread pool | `rayon::par_iter` |
| `mimalloc`, `libraduls`, pre-built `.a` | not needed |
| AGC V1/V2 format compatibility | new format only |
| pybind11 | PyO3 (already in pgr-tk) |
