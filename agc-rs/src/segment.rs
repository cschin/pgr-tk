use crate::error::{AgcError, Result};
use crate::lz_diff::{LzDiff, LzVersion};

// ---------------------------------------------------------------------------
// Compression parameters
// ---------------------------------------------------------------------------

/// Parameters that control segment compression.
pub struct Params {
    /// Minimum match length for the LZ-diff index (default 18).
    pub min_match_len: u32,
    /// Target segment size in bases (default 60 000).
    pub segment_size: u32,
    /// K-mer length used to find splitter boundaries (AGC strategy, default 31).
    pub splitter_k: u32,
    /// LZ-diff encoding version (default V1).
    pub lz_version: LzVersion,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            min_match_len: 18,
            segment_size: 60_000,
            splitter_k: 31,
            lz_version: LzVersion::V1,
        }
    }
}

// ---------------------------------------------------------------------------
// Reference compression
// ---------------------------------------------------------------------------

/// ZSTD-compress a 2-bit encoded reference sequence for storage in the
/// `segment_group.ref_data` column.
///
/// Compression level 13 is used for a good size/speed tradeoff.
pub fn compress_reference(reference: &[u8]) -> Result<Vec<u8>> {
    zstd::encode_all(reference, 19).map_err(|e| AgcError::Zstd(e.to_string()))
}

/// Decompress a `ref_data` BLOB back to its 2-bit encoded sequence.
pub fn decompress_reference(blob: &[u8]) -> Result<Vec<u8>> {
    zstd::decode_all(blob).map_err(|e| AgcError::Zstd(e.to_string()))
}

// ---------------------------------------------------------------------------
// Delta encoding (raw — no ZSTD at this layer)
// ---------------------------------------------------------------------------

/// Encode `query` as an LZ-diff delta against the reference already loaded
/// into `lz` and return the raw byte stream.
///
/// The caller is responsible for compression: either store the bytes
/// temporarily in `segment.delta_data` and batch-compress later via
/// [`batch_compress_deltas`], or compress individually for testing.
pub fn compress_delta(lz: &mut LzDiff, query: &[u8]) -> Result<Vec<u8>> {
    Ok(lz.encode(query))
}

/// Decode a raw LZ-diff byte stream against the reference stored in `lz`.
pub fn decompress_delta(lz: &LzDiff, raw: &[u8]) -> Result<Vec<u8>> {
    lz.decode(raw)
}

// ---------------------------------------------------------------------------
// Per-group batch delta compression
// ---------------------------------------------------------------------------

/// Concatenate multiple raw LZ-diff byte streams and ZSTD-compress them as a
/// single blob for storage in `segment_group.delta_blob`.
///
/// Frame format (little-endian):
/// ```text
/// [n : u32]
/// [len_0 : u32][bytes_0 ...]
/// [len_1 : u32][bytes_1 ...]
/// ...
/// ```
///
/// `in_group_id` values are 1-based; element at index `i` (0-based) in
/// `raw_deltas` corresponds to `in_group_id = i + 1`.
pub fn batch_compress_deltas(raw_deltas: &[Vec<u8>]) -> Result<Vec<u8>> {
    let n = raw_deltas.len() as u32;
    let payload = 4 + raw_deltas.iter().map(|d| 4 + d.len()).sum::<usize>();
    let mut buf = Vec::with_capacity(payload);
    buf.extend_from_slice(&n.to_le_bytes());
    for d in raw_deltas {
        buf.extend_from_slice(&(d.len() as u32).to_le_bytes());
        buf.extend_from_slice(d);
    }
    zstd::encode_all(buf.as_slice(), 19).map_err(|e| AgcError::Zstd(e.to_string()))
}

/// ZSTD-decompress a `delta_blob` and return the raw LZ-diff bytes at
/// position `idx` (0-based, so `in_group_id - 1`).
pub fn extract_delta_from_batch(blob: &[u8], idx: usize) -> Result<Vec<u8>> {
    let frame = zstd::decode_all(blob).map_err(|e| AgcError::Zstd(e.to_string()))?;
    if frame.len() < 4 {
        return Err(AgcError::LzDiff("delta_blob frame too short".to_string()));
    }
    let n = u32::from_le_bytes(frame[..4].try_into().unwrap()) as usize;
    if idx >= n {
        return Err(AgcError::LzDiff(format!(
            "delta index {idx} out of range (batch has {n} entries)"
        )));
    }
    let mut pos = 4usize;
    for i in 0..=idx {
        if pos + 4 > frame.len() {
            return Err(AgcError::LzDiff("truncated delta_blob frame".to_string()));
        }
        let len = u32::from_le_bytes(frame[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        if pos + len > frame.len() {
            return Err(AgcError::LzDiff("truncated delta entry in batch".to_string()));
        }
        if i == idx {
            return Ok(frame[pos..pos + len].to_vec());
        }
        pos += len;
    }
    unreachable!()
}

// ---------------------------------------------------------------------------
// Convenience constructor
// ---------------------------------------------------------------------------

/// Build a fully prepared [`LzDiff`] from a stored `ref_data` BLOB.
///
/// This decompresses the blob, creates an `LzDiff` with the given `params`,
/// and calls `prepare()` so the instance is ready for encoding or decoding.
pub fn lz_from_ref_blob(blob: &[u8], params: &Params) -> Result<LzDiff> {
    let reference = decompress_reference(blob)?;
    let mut lz = LzDiff::new(params.lz_version, params.min_match_len);
    lz.prepare(&reference);
    Ok(lz)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic sequence of `len` bases with a deterministic pattern.
    fn synthetic_seq(len: usize, seed: u8) -> Vec<u8> {
        (0..len)
            .map(|i| (i as u8 ^ seed) & 0x3) // values in 0..=3
            .collect()
    }

    #[test]
    fn compress_decompress_reference() {
        let reference = synthetic_seq(10_000, 0x5A);
        let blob = compress_reference(&reference).expect("compress_reference");
        let recovered = decompress_reference(&blob).expect("decompress_reference");
        assert_eq!(reference, recovered);
    }

    #[test]
    fn compress_decompress_delta_three_seqs() {
        let params = Params::default();
        let reference = synthetic_seq(1_000, 0x00);

        // Three queries: first identical to ref, second with a few SNPs,
        // third with more differences.
        let mut query1 = reference.clone();
        let mut query2 = reference.clone();
        let mut query3 = reference.clone();
        // Introduce some single-base differences (SNPs in 2-bit space).
        for i in (0..1_000).step_by(50) {
            query2[i] = (query2[i] + 1) & 0x3;
        }
        for i in (0..1_000).step_by(20) {
            query3[i] = (query3[i] + 2) & 0x3;
        }

        let ref_blob = compress_reference(&reference).expect("compress ref");
        let mut lz = lz_from_ref_blob(&ref_blob, &params).expect("lz_from_ref_blob");

        let d1 = compress_delta(&mut lz, &query1).expect("compress delta1");
        let d2 = compress_delta(&mut lz, &query2).expect("compress delta2");
        let d3 = compress_delta(&mut lz, &query3).expect("compress delta3");

        // Rebuild lz for decoding (same reference).
        let lz_dec = lz_from_ref_blob(&ref_blob, &params).expect("lz_dec");
        let r1 = decompress_delta(&lz_dec, &d1).expect("decompress delta1");
        let r2 = decompress_delta(&lz_dec, &d2).expect("decompress delta2");
        let r3 = decompress_delta(&lz_dec, &d3).expect("decompress delta3");

        assert_eq!(r1, query1, "query1 round-trip failed");
        assert_eq!(r2, query2, "query2 round-trip failed");
        assert_eq!(r3, query3, "query3 round-trip failed");
    }

    #[test]
    fn subrange_raw_length_accounting() {
        // Scenario: a single contig of length 4_600_000 stored as one segment
        // (raw_lengths = [4_600_000, 0]).
        // Verify that boundary math correctly identifies the segment covering
        // the full range [0, 4_600_000).
        let raw_lengths: &[u64] = &[4_600_000, 0];

        let query_start: u64 = 0;
        let query_end: u64 = 4_600_000;

        // Walk segments to find those that overlap [query_start, query_end).
        let mut offset: u64 = 0;
        let mut overlapping_indices = Vec::new();
        for (idx, &rl) in raw_lengths.iter().enumerate() {
            if rl == 0 {
                break;
            }
            let seg_start = offset;
            let seg_end = offset + rl;
            if seg_end > query_start && seg_start < query_end {
                overlapping_indices.push(idx);
            }
            offset += rl;
        }

        assert_eq!(
            overlapping_indices,
            vec![0],
            "single segment should cover the full contig"
        );
        assert_eq!(
            offset, 4_600_000,
            "total accumulated length should equal contig length"
        );
    }
}
