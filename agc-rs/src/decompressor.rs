use std::path::Path;

use rusqlite::params;

use crate::db::AgcDb;
use crate::error::{AgcError, Result};
use crate::segment;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Read-only (or read-write) view of an agc-rs archive.
///
/// Sequences are decompressed on demand; no data is cached between calls.
pub struct AgcFile {
    db: AgcDb,
}

impl AgcFile {
    /// Open the database at `path` for read-write access.
    pub fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            db: AgcDb::open(path)?,
        })
    }

    /// Open the database at `path` in read-only mode.
    pub fn open_readonly(path: &Path) -> Result<Self> {
        Ok(Self {
            db: AgcDb::open_readonly(path)?,
        })
    }

    // -----------------------------------------------------------------------
    // Metadata queries
    // -----------------------------------------------------------------------

    /// Return the number of samples in the archive.
    pub fn n_samples(&self) -> Result<usize> {
        let count: i64 = self
            .db
            .conn()
            .query_row("SELECT COUNT(*) FROM sample", [], |r| r.get(0))?;
        Ok(count as usize)
    }

    /// Return the number of contigs belonging to `sample`.
    pub fn n_contigs(&self, sample: &str) -> Result<usize> {
        let sample_id = self.sample_id(sample)?;
        let count: i64 = self.db.conn().query_row(
            "SELECT COUNT(*) FROM contig WHERE sample_id = ?1",
            params![sample_id],
            |r| r.get(0),
        )?;
        Ok(count as usize)
    }

    /// List all sample names in the archive.
    pub fn list_samples(&self) -> Result<Vec<String>> {
        let mut stmt = self
            .db
            .conn()
            .prepare("SELECT name FROM sample ORDER BY id")?;
        let names = stmt
            .query_map([], |r| r.get(0))?
            .collect::<std::result::Result<Vec<String>, _>>()?;
        Ok(names)
    }

    /// List the contig names for `sample`.
    pub fn list_contigs(&self, sample: &str) -> Result<Vec<String>> {
        let sample_id = self.sample_id(sample)?;
        let mut stmt = self
            .db
            .conn()
            .prepare("SELECT name FROM contig WHERE sample_id = ?1 ORDER BY id")?;
        let names = stmt
            .query_map(params![sample_id], |r| r.get(0))?
            .collect::<std::result::Result<Vec<String>, _>>()?;
        Ok(names)
    }

    /// Return the total length (in bases) of `contig` within `sample`.
    pub fn contig_len(&self, sample: &str, contig: &str) -> Result<u64> {
        let sample_id = self.sample_id(sample)?;
        let len: i64 = self
            .db
            .conn()
            .query_row(
                "SELECT length FROM contig WHERE sample_id = ?1 AND name = ?2",
                params![sample_id, contig],
                |r| r.get(0),
            )
            .map_err(|_| AgcError::ContigNotFound {
                sample: sample.to_string(),
                contig: contig.to_string(),
            })?;
        Ok(len as u64)
    }

    // -----------------------------------------------------------------------
    // Sequence extraction
    // -----------------------------------------------------------------------

    /// Extract the subsequence `[start, end)` (0-based, half-open, in bases)
    /// for `contig` within `sample`.
    ///
    /// Returns the sequence as uppercase ASCII bytes.
    pub fn contig_seq(&self, sample: &str, contig: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        // 1. Contig metadata.
        let total_len = self.contig_len(sample, contig)?;
        let contig_id = self.contig_id(sample, contig)?;

        // 2. Validate range.
        if start > end || end > total_len {
            return Err(AgcError::RangeOutOfBounds {
                start: start as usize,
                end: end as usize,
                len: total_len as usize,
            });
        }
        if start == end {
            return Ok(Vec::new());
        }

        // 3. Load all segments for this contig, ordered by seg_order.
        struct SegRow {
            group_id: i64,
            in_group_id: i64,
            is_rev_comp: bool,
            raw_length: u64,
            delta_data: Option<Vec<u8>>,
        }

        let mut stmt = self.db.conn().prepare(
            "SELECT group_id, in_group_id, is_rev_comp, raw_length, delta_data \
             FROM segment \
             WHERE contig_id = ?1 \
             ORDER BY seg_order",
        )?;

        let segs: Vec<SegRow> = stmt
            .query_map(params![contig_id], |r| {
                Ok(SegRow {
                    group_id: r.get(0)?,
                    in_group_id: r.get(1)?,
                    is_rev_comp: r.get::<_, bool>(2)?,
                    raw_length: r.get::<_, i64>(3)? as u64,
                    delta_data: r.get(4)?,
                })
            })?
            .collect::<std::result::Result<Vec<_>, _>>()?;

        // 4 & 5. Walk segments, decompress those that overlap [start, end).
        let mut result_bits: Vec<u8> = Vec::with_capacity((end - start) as usize);
        let mut offset: u64 = 0;

        for seg in &segs {
            let seg_start = offset;
            let seg_end = offset + seg.raw_length;

            if seg_end > start && seg_start < end {
                // This segment contributes to the output.
                let seq_bits =
                    self.decompress_segment(seg.group_id, seg.in_group_id, &seg.delta_data)?;

                // 6. Handle reverse complement.
                let seq_bits = if seg.is_rev_comp {
                    rev_comp_bits(&seq_bits)
                } else {
                    seq_bits
                };

                // 7. Strip the overlap prefix.
                //
                // Stored segments include a k-base overlap prefix from the
                // previous splitter k-mer (AGC's `split_pos = pos+1-k`
                // scheme).  `raw_length` is the actual contig contribution
                // (stored_len - k for non-first segments).  Deriving the
                // overlap from the stored/decoded length avoids needing to
                // know k explicitly, and is automatically zero for first
                // segments (where raw_length == stored_length).
                let overlap = seq_bits.len().saturating_sub(seg.raw_length as usize);
                let seq_bits = &seq_bits[overlap..];

                // Trim to the portion that falls within [start, end).
                let local_start = if start > seg_start {
                    (start - seg_start) as usize
                } else {
                    0
                };
                let local_end = if end < seg_end {
                    (end - seg_start) as usize
                } else {
                    seq_bits.len()
                };

                result_bits.extend_from_slice(&seq_bits[local_start..local_end]);
            }

            offset = seg_end;
            if offset >= end {
                break;
            }
        }

        // 7. Convert 2-bit to uppercase ASCII.
        Ok(bits_to_ascii(&result_bits))
    }

    /// Retrieve the full sequence of `contig` within `sample` as uppercase ASCII.
    pub fn full_contig(&self, sample: &str, contig: &str) -> Result<Vec<u8>> {
        let len = self.contig_len(sample, contig)?;
        self.contig_seq(sample, contig, 0, len)
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    fn sample_id(&self, sample: &str) -> Result<i64> {
        self.db
            .conn()
            .query_row(
                "SELECT id FROM sample WHERE name = ?1",
                params![sample],
                |r| r.get(0),
            )
            .map_err(|_| AgcError::SampleNotFound(sample.to_string()))
    }

    fn contig_id(&self, sample: &str, contig: &str) -> Result<i64> {
        let sample_id = self.sample_id(sample)?;
        self.db
            .conn()
            .query_row(
                "SELECT id FROM contig WHERE sample_id = ?1 AND name = ?2",
                params![sample_id, contig],
                |r| r.get(0),
            )
            .map_err(|_| AgcError::ContigNotFound {
                sample: sample.to_string(),
                contig: contig.to_string(),
            })
    }

    /// Decompress one segment.
    ///
    /// - `in_group_id == 0`: reference segment — decompress `ref_data`.
    /// - `in_group_id > 0`: delta segment — decompress `delta_data` (ZSTD),
    ///   then LZ-diff decode against the group's reference.
    fn decompress_segment(
        &self,
        group_id: i64,
        in_group_id: i64,
        delta_data: &Option<Vec<u8>>,
    ) -> Result<Vec<u8>> {
        if in_group_id == 0 {
            let ref_blob: Vec<u8> = self.db.conn().query_row(
                "SELECT ref_data FROM segment_group WHERE id = ?1",
                params![group_id],
                |r| r.get(0),
            )?;
            segment::decompress_reference(&ref_blob)
        } else {
            let blob = delta_data.as_ref().ok_or_else(|| {
                AgcError::LzDiff(format!("segment in group {} has no delta_data", group_id))
            })?;
            let ref_blob: Vec<u8> = self.db.conn().query_row(
                "SELECT ref_data FROM segment_group WHERE id = ?1",
                params![group_id],
                |r| r.get(0),
            )?;
            let raw = segment::decompress_delta_data(blob)?;
            let params = segment::Params::default();
            let lz = segment::lz_from_ref_blob(&ref_blob, &params)?;
            segment::decompress_delta(&lz, &raw)
        }
    }
}

// ---------------------------------------------------------------------------
// Helper functions (public so callers can reuse them)
// ---------------------------------------------------------------------------

/// Convert a 2-bit encoded sequence to uppercase ASCII.
///
/// | Bits | ASCII |
/// |------|-------|
/// |  0   | A     |
/// |  1   | C     |
/// |  2   | G     |
/// |  3   | T     |
/// |  4+  | N     |
pub fn bits_to_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect()
}

/// Reverse complement of a 2-bit encoded sequence.
///
/// Complement: A(0)<->T(3), C(1)<->G(2).  N(4) stays N.
pub fn rev_comp_bits(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            0 => 3, // A -> T
            1 => 2, // C -> G
            2 => 1, // G -> C
            3 => 0, // T -> A
            x => x, // N or invalid
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::db::AgcDb;
    use crate::fasta_io::FastaRecord;
    use crate::segment;
    use rusqlite::params;
    use tempfile::NamedTempFile;

    // -----------------------------------------------------------------------
    // Setup helpers
    // -----------------------------------------------------------------------

    fn temp_db() -> (NamedTempFile, AgcDb) {
        let f = NamedTempFile::new().expect("tempfile");
        let path = f.path().to_path_buf();
        std::fs::remove_file(&path).ok();
        let db = AgcDb::create(&path).expect("create db");
        (f, db)
    }

    /// Insert a single contig backed by a single segment group and segment.
    ///
    /// `sequence_ascii` must be uppercase ASCII (A/C/G/T/N).
    fn insert_synthetic_contig(db: &AgcDb, sample: &str, contig: &str, sequence_ascii: &[u8]) {
        let conn = db.conn();

        // sample
        conn.execute(
            "INSERT OR IGNORE INTO sample (name) VALUES (?1)",
            params![sample],
        )
        .expect("insert sample");
        let sample_id: i64 = conn
            .query_row(
                "SELECT id FROM sample WHERE name = ?1",
                params![sample],
                |r| r.get(0),
            )
            .expect("sample id");

        // contig
        conn.execute(
            "INSERT OR IGNORE INTO contig (sample_id, name, length) VALUES (?1, ?2, ?3)",
            params![sample_id, contig, sequence_ascii.len() as i64],
        )
        .expect("insert contig");
        let contig_id: i64 = conn
            .query_row(
                "SELECT id FROM contig WHERE sample_id = ?1 AND name = ?2",
                params![sample_id, contig],
                |r| r.get(0),
            )
            .expect("contig id");

        // 2-bit encode
        let seq_bits: Vec<u8> = sequence_ascii
            .iter()
            .map(|&b| match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            })
            .collect();

        let ref_blob = segment::compress_reference(&seq_bits).expect("compress ref");

        // segment_group: ref_data = compressed reference; params JSON
        conn.execute(
            "INSERT INTO segment_group (ref_data, params) VALUES (?1, ?2)",
            params![ref_blob, r#"{"min_match_len":18,"segment_size":60000}"#],
        )
        .expect("insert segment_group");
        let group_id: i64 = conn.last_insert_rowid();

        // segment: in_group_id = 0 means this segment IS the reference;
        // delta_data is NULL.
        conn.execute(
            "INSERT INTO segment (contig_id, seg_order, group_id, in_group_id, is_rev_comp, raw_length, delta_data) \
             VALUES (?1, 0, ?2, 0, 0, ?3, NULL)",
            params![contig_id, group_id, sequence_ascii.len() as i64],
        )
        .expect("insert segment");
    }

    fn open_agc_file(f: &NamedTempFile) -> AgcFile {
        AgcFile {
            db: AgcDb::open(f.path()).expect("open db"),
        }
    }

    // -----------------------------------------------------------------------
    // Tests
    // -----------------------------------------------------------------------

    #[test]
    fn list_samples_empty_db() {
        let (f, _db) = temp_db();
        let agc = open_agc_file(&f);
        let samples = agc.list_samples().expect("list_samples");
        assert!(samples.is_empty());
    }

    #[test]
    fn contig_not_found_returns_err() {
        let (f, db) = temp_db();
        insert_synthetic_contig(&db, "sampleA", "chr1", b"ACGTACGT");
        drop(db);
        let agc = open_agc_file(&f);
        match agc.contig_len("sampleA", "chrNONEXISTENT") {
            Err(AgcError::ContigNotFound { .. }) => {}
            other => panic!("expected ContigNotFound, got {:?}", other),
        }
        match agc.contig_len("NO_SUCH_SAMPLE", "chr1") {
            Err(AgcError::SampleNotFound(_)) => {}
            other => panic!("expected SampleNotFound, got {:?}", other),
        }
    }

    #[test]
    fn full_contig_round_trip() {
        let seq = b"ACGTACGTNNACGT".to_vec();
        let (f, db) = temp_db();
        insert_synthetic_contig(&db, "s1", "c1", &seq);
        drop(db);
        let agc = open_agc_file(&f);
        let result = agc.full_contig("s1", "c1").expect("full_contig");
        assert_eq!(result, seq);
    }

    #[test]
    fn subrange_query() {
        // Build a 2000-base synthetic sequence.
        let seq: Vec<u8> = (0u8..200u8)
            .cycle()
            .take(2000)
            .map(|i| b"ACGT"[(i as usize) % 4])
            .collect();

        let (f, db) = temp_db();
        insert_synthetic_contig(&db, "s2", "c2", &seq);
        drop(db);
        let agc = open_agc_file(&f);

        let full = agc.full_contig("s2", "c2").expect("full");
        let sub = agc.contig_seq("s2", "c2", 500, 1500).expect("subrange");

        assert_eq!(sub, &full[500..1500]);
    }
}
