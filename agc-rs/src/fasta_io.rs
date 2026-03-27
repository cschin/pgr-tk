use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

use crate::error::{AgcError, Result};

/// A single FASTA record.
pub struct FastaRecord {
    /// Sequence name (everything after `>` up to the first whitespace).
    pub name: String,
    /// 2-bit encoded sequence: A=0, C=1, G=2, T=3, N (and anything else)=4.
    pub seq: Vec<u8>,
    /// Original uppercase ASCII bases.
    pub seq_ascii: Vec<u8>,
}

/// Convert an uppercase ASCII nucleotide to its 2-bit code.
#[inline]
fn ascii_to_2bit(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // N or anything ambiguous
    }
}

/// Read a FASTA file at `path`.
///
/// Accepts both plain text and gzip-compressed files (detected by the `.gz`
/// extension).  All bases are uppercased.  Records with empty sequences are
/// skipped.  Returns an error on I/O failure or if the file contains no `>`
/// header at all.
pub fn read_fasta_gz(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path)?;
    let buf = BufReader::new(file);

    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gz {
        let decoder = MultiGzDecoder::new(buf);
        parse_fasta(BufReader::new(decoder))
    } else {
        parse_fasta(buf)
    }
}

/// Parse FASTA from any `BufRead` source.
fn parse_fasta<R: BufRead>(reader: R) -> Result<Vec<FastaRecord>> {
    let mut records: Vec<FastaRecord> = Vec::new();

    let mut current_name: Option<String> = None;
    let mut current_ascii: Vec<u8> = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        let line = line.trim_end().as_bytes().to_vec();

        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            // Flush previous record.
            if let Some(name) = current_name.take() {
                if !current_ascii.is_empty() {
                    let seq: Vec<u8> = current_ascii.iter().map(|&b| ascii_to_2bit(b)).collect();
                    records.push(FastaRecord {
                        name,
                        seq,
                        seq_ascii: current_ascii,
                    });
                }
                current_ascii = Vec::new();
            }

            // Parse the new header: name is everything after '>' up to the
            // first whitespace character.
            let header = std::str::from_utf8(&line[1..])
                .map_err(|e| AgcError::FastaParse(e.to_string()))?;
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
        } else if current_name.is_some() {
            // Sequence line: uppercase and append.
            let uppercased: Vec<u8> = line.iter().map(|b| b.to_ascii_uppercase()).collect();
            current_ascii.extend_from_slice(&uppercased);
        }
    }

    // Flush the last record.
    if let Some(name) = current_name {
        if !current_ascii.is_empty() {
            let seq: Vec<u8> = current_ascii.iter().map(|&b| ascii_to_2bit(b)).collect();
            records.push(FastaRecord {
                name,
                seq,
                seq_ascii: current_ascii,
            });
        }
    }

    Ok(records)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Write a string into a temporary plain-text file and return the handle.
    fn write_tmp_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(content.as_bytes()).expect("write");
        f
    }

    #[test]
    fn parse_in_memory_fasta() {
        let fasta = ">seq1 description\nACGTACGT\nACGT\n>seq2\nNNACGT\n";
        let tmp = write_tmp_fasta(fasta);
        let records = read_fasta_gz(tmp.path()).expect("read_fasta_gz");
        assert_eq!(records.len(), 2);

        let r1 = &records[0];
        assert_eq!(r1.name, "seq1");
        assert_eq!(r1.seq_ascii, b"ACGTACGTACGT");
        // 2-bit: A=0, C=1, G=2, T=3
        assert_eq!(r1.seq, vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]);

        let r2 = &records[1];
        assert_eq!(r2.name, "seq2");
        // N => 4
        assert_eq!(r2.seq[0], 4);
        assert_eq!(r2.seq[1], 4);
    }

    #[test]
    fn uppercase_lowercased_input() {
        let fasta = ">lower\nacgt\n";
        let tmp = write_tmp_fasta(fasta);
        let records = read_fasta_gz(tmp.path()).expect("read");
        assert_eq!(records[0].seq_ascii, b"ACGT");
        assert_eq!(records[0].seq, vec![0, 1, 2, 3]);
    }

    #[test]
    fn skips_empty_sequences() {
        let fasta = ">empty\n\n>real\nACGT\n";
        let tmp = write_tmp_fasta(fasta);
        let records = read_fasta_gz(tmp.path()).expect("read");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "real");
    }

    #[test]
    fn read_real_ecoli_if_present() {
        let p = std::path::Path::new(
            "/Users/cschin/Sandbox/repos/pgr-tk/test_data/ecoli/ecoli_k12_mg1655.fna.gz",
        );
        if !p.exists() {
            return; // skip if not present
        }
        let records = read_fasta_gz(p).expect("read ecoli");
        assert!(!records.is_empty(), "expected at least one record");
        assert!(!records[0].seq_ascii.is_empty(), "sequence must not be empty");
    }
}
