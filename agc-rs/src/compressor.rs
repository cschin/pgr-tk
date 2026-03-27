use std::path::Path;

use rayon::prelude::*;
use rusqlite::params;

use crate::db::AgcDb;
use crate::error::{AgcError, Result};
use crate::fasta_io;
use crate::segment::{self, Params};

/// Full compression pipeline: FASTA → segments → SQLite.
///
/// Phase 1: each contig is split into fixed-size chunks of `params.segment_size`
/// bases. Each chunk is stored as its own `segment_group` (the chunk itself is
/// the reference, `in_group_id = 0`, `delta_data = NULL`).  This is lossless
/// and correct but does not exploit cross-contig similarity via LZ-diff.
///
/// TODO: Phase 2 – group corresponding chunks from different contigs / samples
/// into the same `segment_group` and store subsequent contigs as LZ-diff
/// deltas against the group reference.
pub struct Compressor {
    db: AgcDb,
    params: Params,
}

/// Serialise `params` to the JSON string stored in `segment_group.params`.
fn params_json(p: &Params) -> String {
    format!(
        r#"{{"min_match_len":{},"segment_size":{}}}"#,
        p.min_match_len, p.segment_size
    )
}

impl Compressor {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new archive at `path`.
    ///
    /// Fails if the file already exists and already contains an agc-rs schema.
    pub fn create(path: &Path, params: Params) -> Result<Self> {
        Ok(Self {
            db: AgcDb::create(path)?,
            params,
        })
    }

    /// Open an existing archive for appending.
    ///
    /// The `params` stored in the existing segment groups are left untouched;
    /// new segments written during this session will use `Params::default()`.
    pub fn append(path: &Path) -> Result<Self> {
        Ok(Self {
            db: AgcDb::open(path)?,
            params: Params::default(),
        })
    }

    // -----------------------------------------------------------------------
    // Compression
    // -----------------------------------------------------------------------

    /// Compress all sequences in a FASTA (plain or `.gz`) file and store them
    /// under `sample_name`.
    ///
    /// Returns [`AgcError::SampleNotFound`] (repurposed as "already exists")
    /// if `sample_name` is already in the archive.  Uses a single SQLite
    /// transaction for the whole sample for durability and speed.
    pub fn add_fasta(&mut self, path: &Path, sample_name: &str) -> Result<()> {
        // Reject duplicate sample names early.
        let already: i64 = self.db.conn().query_row(
            "SELECT COUNT(*) FROM sample WHERE name = ?1",
            params![sample_name],
            |r| r.get(0),
        )?;
        if already > 0 {
            return Err(AgcError::SampleNotFound(format!(
                "sample '{}' already exists in the archive",
                sample_name
            )));
        }

        // Read all FASTA records.
        let records = fasta_io::read_fasta_gz(path)?;
        if records.is_empty() {
            return Ok(());
        }

        let segment_size = self.params.segment_size as usize;
        let params_json = params_json(&self.params);

        // Pre-compute all (chunk_bytes, raw_length) pairs in parallel using
        // rayon so ZSTD compression of ref_data blobs runs concurrently.
        //
        // Structure: for each record, a Vec of (ref_blob, raw_length) per chunk.
        let compressed_per_record: Vec<Vec<(Vec<u8>, usize)>> = records
            .par_iter()
            .map(|rec| {
                rec.seq
                    .chunks(segment_size)
                    .map(|chunk| {
                        let blob = segment::compress_reference(chunk)
                            .expect("zstd compress_reference failed");
                        (blob, chunk.len())
                    })
                    .collect()
            })
            .collect();

        // Insert everything inside a single transaction.
        let conn = self.db.conn();
        let tx = conn.unchecked_transaction()?;

        // Insert sample row.
        tx.execute(
            "INSERT INTO sample (name) VALUES (?1)",
            params![sample_name],
        )?;
        let sample_id: i64 = tx.last_insert_rowid();

        for (record, chunks) in records.iter().zip(compressed_per_record.iter()) {
            let total_len: usize = record.seq.len();

            // Insert contig row.
            tx.execute(
                "INSERT INTO contig (sample_id, name, length) VALUES (?1, ?2, ?3)",
                params![sample_id, &record.name, total_len as i64],
            )?;
            let contig_id: i64 = tx.last_insert_rowid();

            for (seg_order, (ref_blob, raw_length)) in chunks.iter().enumerate() {
                // One segment_group per chunk (Phase 1: no cross-contig delta).
                tx.execute(
                    "INSERT INTO segment_group (ref_data, params) VALUES (?1, ?2)",
                    params![ref_blob, &params_json],
                )?;
                let group_id: i64 = tx.last_insert_rowid();

                // The single segment in this group is the reference itself.
                tx.execute(
                    "INSERT INTO segment \
                     (contig_id, seg_order, group_id, in_group_id, is_rev_comp, raw_length, delta_data) \
                     VALUES (?1, ?2, ?3, 0, 0, ?4, NULL)",
                    params![contig_id, seg_order as i64, group_id, *raw_length as i64],
                )?;
            }
        }

        tx.commit()?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Finalisation
    // -----------------------------------------------------------------------

    /// Flush any pending changes and close the database.
    pub fn finish(self) -> Result<()> {
        // Connection is closed on drop; nothing else needed for SQLite WAL.
        drop(self.db);
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decompressor::AgcFile;
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Write an in-memory FASTA string to a temporary file and return the handle.
    fn write_tmp_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(content.as_bytes()).expect("write fasta");
        f
    }

    /// A temporary directory + archive path inside it.
    fn tmp_archive() -> (TempDir, std::path::PathBuf) {
        let dir = TempDir::new().expect("tempdir");
        let path = dir.path().join("test.agcrs");
        (dir, path)
    }

    // -----------------------------------------------------------------------
    // Tests
    // -----------------------------------------------------------------------

    #[test]
    fn compress_and_list_sample() {
        let fasta = ">chr1\nACGTACGTACGT\n>chr2\nTTTTGGGG\n";
        let fasta_file = write_tmp_fasta(fasta);
        let (_dir, archive) = tmp_archive();

        let mut c = Compressor::create(&archive, Params::default()).unwrap();
        c.add_fasta(fasta_file.path(), "sample1").unwrap();
        c.finish().unwrap();

        let agc = AgcFile::open(&archive).unwrap();
        assert_eq!(agc.n_samples().unwrap(), 1);

        let contigs = agc.list_contigs("sample1").unwrap();
        assert_eq!(contigs.len(), 2);
        assert!(contigs.contains(&"chr1".to_string()));
        assert!(contigs.contains(&"chr2".to_string()));

        assert_eq!(agc.contig_len("sample1", "chr1").unwrap(), 12);
        assert_eq!(agc.contig_len("sample1", "chr2").unwrap(), 8);
    }

    #[test]
    fn compress_and_decompress_sequence() {
        let seq = "ACGTACGTNNACGTTTTTGGGGCCCCAAAA";
        let fasta = format!(">contig1\n{}\n", seq);
        let fasta_file = write_tmp_fasta(&fasta);
        let (_dir, archive) = tmp_archive();

        let mut c = Compressor::create(&archive, Params::default()).unwrap();
        c.add_fasta(fasta_file.path(), "mysample").unwrap();
        c.finish().unwrap();

        let agc = AgcFile::open(&archive).unwrap();
        let recovered = agc.full_contig("mysample", "contig1").unwrap();
        // N stays as N in ASCII round-trip (ascii_to_2bit gives 4, bits_to_ascii gives N).
        let expected: Vec<u8> = seq.bytes().map(|b| {
            match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
                _ => b'N',
            }
        }).collect();
        assert_eq!(recovered, expected);
    }

    #[test]
    fn append_second_sample() {
        let fasta_a = ">chr1\nACGTACGT\n";
        let fasta_b = ">chr1\nTTTTCCCC\n";
        let fa = write_tmp_fasta(fasta_a);
        let fb = write_tmp_fasta(fasta_b);
        let (_dir, archive) = tmp_archive();

        // Create with sample A.
        let mut c = Compressor::create(&archive, Params::default()).unwrap();
        c.add_fasta(fa.path(), "sampleA").unwrap();
        c.finish().unwrap();

        // Append sample B.
        let mut c = Compressor::append(&archive).unwrap();
        c.add_fasta(fb.path(), "sampleB").unwrap();
        c.finish().unwrap();

        let agc = AgcFile::open(&archive).unwrap();
        assert_eq!(agc.n_samples().unwrap(), 2);

        let samples = agc.list_samples().unwrap();
        assert!(samples.contains(&"sampleA".to_string()));
        assert!(samples.contains(&"sampleB".to_string()));
    }

    #[test]
    fn duplicate_sample_returns_error() {
        let fasta = ">chr1\nACGT\n";
        let fa = write_tmp_fasta(fasta);
        let (_dir, archive) = tmp_archive();

        let mut c = Compressor::create(&archive, Params::default()).unwrap();
        c.add_fasta(fa.path(), "dup").unwrap();
        // Second add with same name must fail.
        let result = c.add_fasta(fa.path(), "dup");
        assert!(
            result.is_err(),
            "expected error for duplicate sample, got Ok"
        );
    }

    #[test]
    fn multi_segment_contig_round_trip() {
        // Use a very small segment_size to force multi-segment storage.
        let params = Params {
            segment_size: 10,
            ..Params::default()
        };
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bases → 4 segments of 10 + 1 of 2
        let fasta = format!(">longseq\n{}\n", seq);
        let fa = write_tmp_fasta(&fasta);
        let (_dir, archive) = tmp_archive();

        let mut c = Compressor::create(&archive, params).unwrap();
        c.add_fasta(fa.path(), "s").unwrap();
        c.finish().unwrap();

        let agc = AgcFile::open(&archive).unwrap();
        let recovered = agc.full_contig("s", "longseq").unwrap();
        assert_eq!(recovered, seq.as_bytes());
    }
}
