use std::cell::RefCell;

use crate::error::{AgcError, Result};
use crate::lz_diff::{LzDiff, LzVersion};

// ---------------------------------------------------------------------------
// Thread-local ZSTD compressor — one context per rayon worker thread.
//
// Allocating a fresh ZSTD_CCtx (~1–8 MB) per `encode_all` call at level 19
// was triggering hundreds of millions of minor page faults when 50 000+
// segments were compressed concurrently (reference path).  Re-using a single
// context per thread cuts allocator pressure to O(n_threads) instead of
// O(n_segments).
//
// Compression level 9 gives a good size/speed tradeoff and matches the level
// commonly used in C++ AGC.  Level 19 has diminishing size returns but is
// orders of magnitude slower.
// ---------------------------------------------------------------------------

thread_local! {
    static ZSTD_CTX: RefCell<zstd::bulk::Compressor<'static>> = RefCell::new(
        zstd::bulk::Compressor::new(9).expect("zstd compressor init"),
    );
    static ZSTD_DCTX: RefCell<zstd::bulk::Decompressor<'static>> = RefCell::new(
        zstd::bulk::Decompressor::new().expect("zstd decompressor init"),
    );
}

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
/// Uses a thread-local compressor context (level 9) so that each rayon
/// worker thread reuses one ZSTD_CCtx across all segments it processes,
/// avoiding the per-call context allocation overhead that causes page-fault
/// pressure when tens of thousands of segments are compressed in parallel.
pub fn compress_reference(reference: &[u8]) -> Result<Vec<u8>> {
    ZSTD_CTX.with(|ctx| {
        ctx.borrow_mut()
            .compress(reference)
            .map_err(|e| AgcError::Zstd(e.to_string()))
    })
}

/// Decompress a `ref_data` BLOB back to its 2-bit encoded sequence.
///
/// Uses a thread-local decompressor context so that parallel calls in
/// `build_fallback_index` (one per segment group) do not each allocate a
/// new ZSTD_DCtx.  `bulk::Compressor` stores the content size in the frame
/// header, so we read the exact capacity from the frame rather than guessing.
pub fn decompress_reference(blob: &[u8]) -> Result<Vec<u8>> {
    let capacity = zstd::zstd_safe::get_frame_content_size(blob)
        .ok()
        .flatten()
        .unwrap_or(blob.len() as u64 * 20) as usize;
    ZSTD_DCTX.with(|ctx| {
        ctx.borrow_mut()
            .decompress(blob, capacity)
            .map_err(|e| AgcError::Zstd(e.to_string()))
    })
}

// ---------------------------------------------------------------------------
// Delta encoding (raw — no ZSTD at this layer)
// ---------------------------------------------------------------------------

/// Encode `query` as an LZ-diff delta against the reference already loaded
/// into `lz` and return the raw byte stream.
///
/// Compress the result with [`compress_delta_data`] before storing in
/// `segment.delta_data`.
pub fn compress_delta(lz: &LzDiff, query: &[u8]) -> Result<Vec<u8>> {
    Ok(lz.encode(query))
}

/// Decode a raw LZ-diff byte stream against the reference stored in `lz`.
pub fn decompress_delta(lz: &LzDiff, raw: &[u8]) -> Result<Vec<u8>> {
    lz.decode(raw)
}

/// ZSTD-compress a single raw LZ-diff delta for storage in `segment.delta_data`.
///
/// Uses the thread-local compressor so that parallel calls within the append
/// path do not each allocate a new ZSTD_CCtx.
pub fn compress_delta_data(raw: &[u8]) -> Result<Vec<u8>> {
    ZSTD_CTX.with(|ctx| {
        ctx.borrow_mut()
            .compress(raw)
            .map_err(|e| AgcError::Zstd(e.to_string()))
    })
}

/// Decompress a `segment.delta_data` BLOB back to the raw LZ-diff bytes.
pub fn decompress_delta_data(blob: &[u8]) -> Result<Vec<u8>> {
    let capacity = zstd::zstd_safe::get_frame_content_size(blob)
        .ok()
        .flatten()
        .unwrap_or(blob.len() as u64 * 20) as usize;
    ZSTD_DCTX.with(|ctx| {
        ctx.borrow_mut()
            .decompress(blob, capacity)
            .map_err(|e| AgcError::Zstd(e.to_string()))
    })
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
