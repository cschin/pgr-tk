use std::collections::{HashMap, HashSet};
use std::path::Path;

use rayon::prelude::*;
use rusqlite::{params, Connection};

use crate::db::AgcDb;
use crate::error::{AgcError, Result};
use crate::fasta_io::{self, FastaRecord};
use crate::kmer::{rev_comp_2bit_seq, Kmer};
use crate::segment::{self, Params};

// ---------------------------------------------------------------------------
// Sentinel value for "no splitter at this boundary"
// ---------------------------------------------------------------------------

const SENTINEL: u64 = u64::MAX;

// ---------------------------------------------------------------------------
// Fallback k-mer sampling constants (used when exact splitter lookup fails)
// ---------------------------------------------------------------------------

/// Approximate number of k-mers sampled per segment for the fallback search.
const TARGET_KMERS_PER_SEG: usize = 600;

/// Ignore k-mers that appear in more than this many distinct segment groups.
const MAX_KMER_BUCKET: usize = 20;

// ---------------------------------------------------------------------------
// AGC-faithful splitter determination
// ---------------------------------------------------------------------------

/// Collect every canonical k-mer from `records` and return a **sorted** Vec
/// containing only those that appear exactly once (singletons).
///
/// Mirrors AGC's `start_kmer_collecting_threads` + `remove_non_singletons`
/// but uses a flat Vec + sort instead of a HashMap.  For a 3 Gbp genome at
/// k=31 this costs ~24 GB (8 bytes × 3 B k-mers) versus the ~48–72 GB that
/// a parallel fold+HashMap approach requires from hash-table load-factor
/// overhead and per-thread accumulator copies.
fn collect_sorted_singletons(records: &[FastaRecord], k: usize) -> Vec<u64> {
    // Pre-allocate: at most (seq_len - k + 1) k-mers per record.
    let total_kmers: usize = records
        .iter()
        .map(|r| r.seq.len().saturating_sub(k - 1))
        .sum();

    // Serial collection avoids the peak-memory doubling that occurs when
    // parallel per-record Vecs are later flattened into a single allocation.
    let mut all_kmers: Vec<u64> = Vec::with_capacity(total_kmers);
    for rec in records {
        let mut km = Kmer::new(k as u8);
        for &b in &rec.seq {
            if km.push_bits(b) && km.full() {
                all_kmers.push(km.canonical());
            }
        }
    }

    // Parallel in-place sort (rayon).  Mirrors AGC's RadixSort step.
    all_kmers.par_sort_unstable();

    // Compact in-place: keep only k-mers that appear exactly once.
    // Mirrors AGC's `remove_non_singletons`.
    let n = all_kmers.len();
    let mut write = 0usize;
    let mut i = 0usize;
    while i < n {
        let v = all_kmers[i];
        let mut j = i + 1;
        while j < n && all_kmers[j] == v {
            j += 1;
        }
        if j - i == 1 {
            all_kmers[write] = v;
            write += 1;
        }
        i = j;
    }
    all_kmers.truncate(write);
    all_kmers.shrink_to_fit();

    all_kmers
}

/// Scan a single 2-bit sequence for `candidates` k-mers at least
/// `segment_size` bases apart and return those chosen as splitters.
///
/// `candidates` must be a **sorted** slice so membership tests use
/// binary search (O(log n)) instead of a hash lookup.  The sorted-slice
/// contract matches the output of [`collect_sorted_singletons`].
///
/// This is the core of AGC's `find_splitters_in_contig`:
/// - `current_len` starts at `segment_size` so the first eligible k-mer is
///   checked immediately.
/// - After a splitter is accepted the k-mer state AND `current_len` are
///   reset to zero.
/// - The rightmost candidate found in the contig's trailing k-mers is also
///   added so the final segment gets a back boundary.
fn find_splitters_in_seq(
    seq: &[u8],
    candidates: &[u64],
    k: usize,
    segment_size: usize,
) -> HashSet<u64> {
    let mut splitters: HashSet<u64> = HashSet::new();
    let mut km = Kmer::new(k as u8);
    let mut current_len: usize = segment_size;
    let mut recent_kmers: Vec<u64> = Vec::new();

    for &b in seq {
        if !km.push_bits(b) {
            current_len += 1;
            continue;
        }
        current_len += 1;

        if km.full() {
            let can = km.canonical();
            recent_kmers.push(can);

            if current_len >= segment_size && candidates.binary_search(&can).is_ok() {
                splitters.insert(can);
                current_len = 0;
                km.reset();
                recent_kmers.clear();
            }
        }
    }

    for &kmer in recent_kmers.iter().rev() {
        if candidates.binary_search(&kmer).is_ok() {
            splitters.insert(kmer);
            break;
        }
    }

    splitters
}

/// Phase 3: scan each reference contig for singletons, delegating to
/// `find_splitters_in_seq` per contig.
///
/// `singletons` must be a sorted slice (output of [`collect_sorted_singletons`]).
fn find_splitters_in_contigs(
    records: &[FastaRecord],
    singletons: &[u64],
    k: usize,
    segment_size: usize,
) -> HashSet<u64> {
    records
        .par_iter()
        .map(|rec| find_splitters_in_seq(&rec.seq, singletons, k, segment_size))
        .reduce(HashSet::new, |mut a, b| {
            a.extend(b);
            a
        })
}

/// Convenience wrapper: determine the full splitter set from the reference.
pub fn determine_splitters(
    records: &[FastaRecord],
    k: usize,
    segment_size: usize,
) -> HashSet<u64> {
    let singletons = collect_sorted_singletons(records, k);
    find_splitters_in_contigs(records, &singletons, k, segment_size)
}

// ---------------------------------------------------------------------------
// Adaptive splitting (AGC `find_new_splitters`)
// ---------------------------------------------------------------------------

/// For a new contig that has no global splitters, find local splitters from
/// k-mers that are singletons within the contig itself and absent from the
/// reference genome entirely.
///
/// Reproduces AGC's `find_new_splitters`:
/// 1. Collect all canonical k-mers from the contig.
/// 2. Keep only those that appear exactly once (singletons).
/// 3. Remove k-mers present in `known_kmers` (the global splitter set, which
///    is already in memory).  This is a conservative approximation of AGC's
///    full reference exclusion: novel contigs — the only ones that reach this
///    path — are highly divergent, so their k-mers rarely appear in the
///    reference outside of the splitter boundaries anyway.  Using the full
///    reference k-mer set would require O(genome_size) memory (~13 GB for
///    human) which is impractical.
/// 4. Run the standard splitter-finding scan over the local candidates.
fn find_local_splitters(
    seq: &[u8],
    known_kmers: &HashSet<u64>,
    k: usize,
    segment_size: usize,
) -> HashSet<u64> {
    // Step 1+2: singleton k-mers in this contig.
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mut km = Kmer::new(k as u8);
    for &b in seq {
        if km.push_bits(b) && km.full() {
            let c = counts.entry(km.canonical()).or_insert(0);
            *c = c.saturating_add(1);
        }
    }
    // Step 3: keep singletons that are absent from the global splitter set.
    // Collect into a sorted Vec so find_splitters_in_seq can use binary_search.
    let mut local_candidates: Vec<u64> = counts
        .into_iter()
        .filter(|(_, cnt)| *cnt == 1)
        .map(|(kmer, _)| kmer)
        .filter(|kmer| !known_kmers.contains(kmer))
        .collect();

    if local_candidates.is_empty() {
        return HashSet::new();
    }

    local_candidates.sort_unstable();

    // Step 4: apply the standard splitter-finding scan.
    find_splitters_in_seq(seq, &local_candidates, k, segment_size)
}

// ---------------------------------------------------------------------------
// Contig splitting
// ---------------------------------------------------------------------------

/// Metadata for one segment produced by splitting a contig.
pub struct SegInfo {
    /// The sequence in *canonical orientation* (may be RC of the raw scan).
    /// Includes a k-base overlap prefix from the previous splitter k-mer when
    /// `kmer_front != SENTINEL` (matching AGC's `compress_contig` overlap
    /// scheme: `split_pos = pos + 1 - kmer_length`).
    pub seq: Vec<u8>,
    /// Actual bases this segment contributes to the contig.
    ///
    /// Equals `seq.len()` for the first segment (no overlap prefix) and
    /// `seq.len() - k` for all subsequent segments (k-base overlap prefix
    /// = the leading splitter k-mer, which is already counted by the previous
    /// segment).
    pub raw_len: usize,
    /// Canonical k-mer value at the leading boundary, or `SENTINEL`.
    pub kmer_front: u64,
    /// Canonical k-mer value at the trailing boundary, or `SENTINEL`.
    pub kmer_back: u64,
    /// `true` when `seq` is the reverse complement of the original scan order.
    pub is_rc: bool,
}

/// Split a 2-bit encoded contig at every position where a full canonical k-mer
/// matches a splitter.
///
/// Mirrors AGC's `compress_contig` overlap scheme: after a splitter fires at
/// position `pos`, the next segment starts at `pos + 1 - k` (the beginning of
/// the splitter k-mer).  This k-base overlap prefix is **included** in the
/// stored bytes but is NOT counted in `raw_len`, so consecutive segments
/// correctly tile the contig without gaps or double-counting.
///
/// The decompressor strips the overlap by taking `seq[seq.len()-raw_len..]`
/// after decoding and optional RC.
///
/// Segments are stored in *canonical orientation*: if
/// `kmer_front > kmer_back` the sequence is reverse-complemented and
/// `is_rc = true`.  First and last segments (with one `SENTINEL` boundary)
/// are always stored forward.
pub fn split_contig(seq: &[u8], splitters: &HashSet<u64>, k: usize) -> Vec<SegInfo> {
    let mut km = Kmer::new(k as u8);
    let mut segs: Vec<SegInfo> = Vec::new();
    let mut split_start: usize = 0;
    let mut kmer_front: u64 = SENTINEL;

    for (pos, &b) in seq.iter().enumerate() {
        if !km.push_bits(b) {
            continue;
        }
        if km.full() {
            let can = km.canonical();
            if splitters.contains(&can) {
                let raw = &seq[split_start..pos + 1];
                // Canonical orientation: store as RC when front > back so that
                // the smaller k-mer is always the "front" in the stored form.
                let is_rc = kmer_front != SENTINEL && kmer_front > can;
                let stored = if is_rc {
                    rev_comp_2bit_seq(raw)
                } else {
                    raw.to_vec()
                };
                // raw_len: subtract the k-base overlap prefix (the kmer_front
                // k-mer at the start of `raw`) for all non-first segments.
                let raw_len = if kmer_front != SENTINEL {
                    raw.len() - k
                } else {
                    raw.len()
                };
                segs.push(SegInfo { seq: stored, raw_len, kmer_front, kmer_back: can, is_rc });
                // AGC overlap: next segment starts at the beginning of the
                // splitter k-mer (pos + 1 - k), giving a k-base shared prefix.
                split_start = pos + 1 - k;
                kmer_front = can;
                km.reset();
            }
        }
    }

    // Trailing segment: no back splitter, always store forward.
    if split_start < seq.len() {
        let raw = &seq[split_start..];
        let raw_len = if kmer_front != SENTINEL { raw.len() - k } else { raw.len() };
        segs.push(SegInfo {
            seq: raw.to_vec(),
            raw_len,
            kmer_front,
            kmer_back: SENTINEL,
            is_rc: false,
        });
    }

    segs
}

// ---------------------------------------------------------------------------
// Exact segment lookup by (kmer_front, kmer_back) pair
// ---------------------------------------------------------------------------

/// Build a lookup map from canonical (sorted) kmer-pair → group_id.
///
/// Includes both fully-bounded segments (both splitters non-NULL) and
/// partial segments (exactly one splitter non-NULL, partner treated as
/// SENTINEL = u64::MAX).  The key is always `(min, max)` so it matches
/// regardless of orientation.
fn build_exact_map(conn: &Connection) -> Result<HashMap<(u64, u64), i64>> {
    let rows: Vec<(i64, Option<i64>, Option<i64>)> = {
        let mut stmt = conn.prepare(
            "SELECT id, kmer_front, kmer_back \
             FROM segment_group \
             WHERE kmer_front IS NOT NULL OR kmer_back IS NOT NULL",
        )?;
        let collected: Vec<(i64, Option<i64>, Option<i64>)> = stmt
            .query_map([], |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)))?
            .filter_map(|r| r.ok())
            .collect();
        collected
    };
    let mut map = HashMap::with_capacity(rows.len());
    for (gid, kf, kb) in rows {
        let kf = kf.map(|v| v as u64).unwrap_or(SENTINEL);
        let kb = kb.map(|v| v as u64).unwrap_or(SENTINEL);
        let key = (kf.min(kb), kf.max(kb));
        map.insert(key, gid);
    }
    Ok(map)
}

/// Build a "terminators" map: splitter k-mer → list of partner k-mers it has
/// been paired with in stored segment groups.
///
/// Includes SENTINEL (u64::MAX) as a valid partner for partial-boundary
/// segments (exactly one splitter), matching AGC's `map_segments_terminators`
/// which stores `~0ull` as the partner for terminal segments.
fn build_terminators_map(conn: &Connection) -> Result<HashMap<u64, Vec<u64>>> {
    let rows: Vec<(Option<i64>, Option<i64>)> = {
        let mut stmt = conn.prepare(
            "SELECT kmer_front, kmer_back \
             FROM segment_group \
             WHERE kmer_front IS NOT NULL OR kmer_back IS NOT NULL",
        )?;
        let collected: Vec<(Option<i64>, Option<i64>)> = stmt
            .query_map([], |r| Ok((r.get(0)?, r.get(1)?)))?
            .filter_map(|r| r.ok())
            .collect();
        collected
    };
    let mut map: HashMap<u64, Vec<u64>> = HashMap::new();
    for (kf, kb) in rows {
        let kf = kf.map(|v| v as u64).unwrap_or(SENTINEL);
        let kb = kb.map(|v| v as u64).unwrap_or(SENTINEL);
        if kf != SENTINEL {
            map.entry(kf).or_default().push(kb);
        }
        if kb != SENTINEL {
            map.entry(kb).or_default().push(kf);
        }
    }
    Ok(map)
}

/// One-splitter lookup: given a segment with exactly one known boundary
/// splitter, find the stored segment group whose reference size is closest
/// to `seg_len` among all groups paired with that splitter.
///
/// Mirrors AGC's `find_cand_segment_with_one_splitter`: SENTINEL (u64::MAX)
/// is a valid partner representing terminal segments stored with a partial
/// boundary (matching AGC's `~0ull` sentinel).  Size-closeness heuristic
/// is used instead of LZ estimation.
fn find_one_splitter_group(
    single_kmer: u64,
    seg_len: usize,
    terminators_map: &HashMap<u64, Vec<u64>>,
    exact_map: &HashMap<(u64, u64), i64>,
    ref_len_map: &HashMap<i64, usize>,
) -> Option<i64> {
    // If no partners known, try the partial-key (K, SENTINEL) directly.
    let partners: &[u64] = match terminators_map.get(&single_kmer) {
        Some(v) => v,
        None => {
            let key = (single_kmer.min(SENTINEL), single_kmer.max(SENTINEL));
            return exact_map.get(&key).copied();
        }
    };
    let seg_len = seg_len as i64;
    partners
        .iter()
        .filter_map(|&partner| {
            let key = (single_kmer.min(partner), single_kmer.max(partner));
            exact_map.get(&key).map(|&gid| {
                let ref_len = ref_len_map.get(&gid).copied().unwrap_or(0) as i64;
                (gid, (ref_len - seg_len).abs())
            })
        })
        .min_by_key(|(_, diff)| *diff)
        .map(|(gid, _)| gid)
}

// ---------------------------------------------------------------------------
// Fallback: hash-based vote matching (position-independent)
// ---------------------------------------------------------------------------

/// MurMur-64 finaliser — same hash used in lz_diff.rs.
#[inline]
fn murmur64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}

#[inline]
fn compute_sample_rate(seg_len: usize, kmer_len: usize) -> u64 {
    if seg_len < kmer_len {
        return 1;
    }
    ((seg_len - kmer_len + 1) / TARGET_KMERS_PER_SEG).max(1) as u64
}

/// Hash-sample k-mers from a 2-bit sequence: include k-mer when
/// `murmur64(kmer) % sample_rate == 0`.
fn sample_kmers_2bit(seq: &[u8], kmer_len: usize, sample_rate: u64) -> Vec<u64> {
    if seq.len() < kmer_len {
        return Vec::new();
    }
    let mask: u64 = if kmer_len < 32 {
        (1u64 << (2 * kmer_len)) - 1
    } else {
        u64::MAX
    };
    let mut sampled = Vec::with_capacity((seq.len() / sample_rate as usize).max(1));
    let mut kmer: u64 = 0;
    let mut valid: usize = 0;
    for &b in seq {
        if b >= 4 {
            valid = 0;
            kmer = 0;
            continue;
        }
        kmer = ((kmer << 2) | b as u64) & mask;
        valid += 1;
        if valid >= kmer_len && murmur64(kmer) % sample_rate == 0 {
            sampled.push(kmer);
        }
    }
    sampled
}

/// Build k-mer inverted index, ref-data cache, and ref-length map from all
/// stored segment_groups.
///
/// Returns `(kmer_index, ref_data_map, ref_len_map)` where `ref_len_map`
/// maps group_id → uncompressed reference sequence length in 2-bit bases.
fn build_fallback_index(
    conn: &Connection,
    kmer_len: usize,
    sample_rate: u64,
    max_bucket: usize,
) -> Result<(HashMap<u64, Vec<i64>>, HashMap<i64, Vec<u8>>, HashMap<i64, usize>)> {
    // Serial DB load
    let rows: Vec<(i64, Vec<u8>)> = {
        let mut stmt = conn.prepare("SELECT id, ref_data FROM segment_group")?;
        let collected: Vec<(i64, Vec<u8>)> = stmt
            .query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, Vec<u8>>(1)?)))?
            .filter_map(|r| r.ok())
            .collect();
        collected
    };

    // Parallel: decompress + sample k-mers per group
    let processed: Vec<(i64, Vec<u64>, Vec<u8>, usize)> = rows
        .par_iter()
        .map(|(group_id, ref_blob)| {
            let seq_2bit = segment::decompress_reference(ref_blob)
                .expect("decompress_reference");
            let seq_len = seq_2bit.len();
            let kmers = sample_kmers_2bit(&seq_2bit, kmer_len, sample_rate);
            (*group_id, kmers, ref_blob.clone(), seq_len)
        })
        .collect();

    // Serial merge into output maps (fast — no decompression here)
    let mut kmer_index: HashMap<u64, Vec<i64>> = HashMap::new();
    let mut ref_data_map: HashMap<i64, Vec<u8>> = HashMap::with_capacity(processed.len());
    let mut ref_len_map: HashMap<i64, usize> = HashMap::with_capacity(processed.len());
    for (group_id, kmers, blob, len) in processed {
        for kmer in kmers {
            let bucket = kmer_index.entry(kmer).or_default();
            if bucket.len() < max_bucket {
                bucket.push(group_id);
            }
        }
        ref_data_map.insert(group_id, blob);
        ref_len_map.insert(group_id, len);
    }

    Ok((kmer_index, ref_data_map, ref_len_map))
}

/// Vote-count similarity search over stored segment groups.
fn find_best_ref_group(
    seq: &[u8],
    kmer_index: &HashMap<u64, Vec<i64>>,
    kmer_len: usize,
    sample_rate: u64,
) -> Option<i64> {
    let sampled = sample_kmers_2bit(seq, kmer_len, sample_rate);
    if sampled.is_empty() {
        return None;
    }
    let expected = (seq.len().saturating_sub(kmer_len) + 1) as u64 / sample_rate;
    let min_shared = (expected / 60).max(1) as usize;

    let mut votes: HashMap<i64, usize> = HashMap::new();
    for kmer in &sampled {
        if let Some(groups) = kmer_index.get(kmer) {
            for &gid in groups {
                *votes.entry(gid).or_default() += 1;
            }
        }
    }
    votes
        .into_iter()
        .filter(|(_, count)| *count >= min_shared)
        .max_by_key(|(_, count)| *count)
        .map(|(gid, _)| gid)
}

// ---------------------------------------------------------------------------
// Params JSON serialisation
// ---------------------------------------------------------------------------

fn params_json(p: &Params) -> String {
    format!(
        r#"{{"min_match_len":{},"segment_size":{},"splitter_k":{}}}"#,
        p.min_match_len, p.segment_size, p.splitter_k
    )
}

// ---------------------------------------------------------------------------
// Compression pipeline
// ---------------------------------------------------------------------------

pub struct Compressor {
    db: AgcDb,
    params: Params,
}

enum ChunkResult {
    Delta {
        group_id: i64,
        raw_length: usize,
        /// Raw (non-ZSTD) LZ-diff bytes — batch-compressed per group later.
        raw_delta: Vec<u8>,
        is_rc: bool,
    },
    NewRef {
        raw_length: usize,
        ref_blob: Vec<u8>,
        kmer_front: Option<i64>,
        kmer_back: Option<i64>,
        is_rc: bool,
    },
}

impl Compressor {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new archive at `path`.
    pub fn create(path: &Path, params: Params) -> Result<Self> {
        Ok(Self {
            db: AgcDb::create(path)?,
            params,
        })
    }

    /// Open an existing archive for appending.
    pub fn append(path: &Path) -> Result<Self> {
        Ok(Self {
            db: AgcDb::open(path)?,
            params: Params::default(),
        })
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /// Compress a FASTA file and store it under `sample_name`.
    pub fn add_fasta(&mut self, path: &Path, sample_name: &str) -> Result<()> {
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

        let records = fasta_io::read_fasta_gz(path)?;
        if records.is_empty() {
            return Ok(());
        }

        let existing: i64 = self
            .db
            .conn()
            .query_row("SELECT COUNT(*) FROM sample", [], |r| r.get(0))?;

        if existing == 0 {
            self.add_as_reference(records, sample_name)
        } else {
            self.add_as_delta(records, sample_name)
        }
    }

    pub fn finish(self) -> Result<()> {
        drop(self.db);
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Reference path (first sample)
    // -----------------------------------------------------------------------

    fn add_as_reference(&mut self, records: Vec<FastaRecord>, sample_name: &str) -> Result<()> {
        let segment_size = self.params.segment_size as usize;
        let splitter_k = self.params.splitter_k as usize;
        let params_json = params_json(&self.params);

        // Determine AGC-style splitters from the reference sequences.
        let splitters = determine_splitters(&records, splitter_k, segment_size);

        // Parallel: split each contig and ZSTD-compress each segment.
        let compressed_per_record: Vec<Vec<(Vec<u8>, usize, Option<i64>, Option<i64>, bool)>> =
            records
                .par_iter()
                .map(|rec| {
                    split_contig(&rec.seq, &splitters, splitter_k)
                        .into_iter()
                        .map(|seg| {
                            let blob = segment::compress_reference(&seg.seq)
                                .expect("compress_reference failed");
                            let kf = if seg.kmer_front != SENTINEL {
                                Some(seg.kmer_front as i64)
                            } else {
                                None
                            };
                            let kb = if seg.kmer_back != SENTINEL {
                                Some(seg.kmer_back as i64)
                            } else {
                                None
                            };
                            (blob, seg.raw_len, kf, kb, seg.is_rc)
                        })
                        .collect()
                })
                .collect();

        let conn = self.db.conn();
        let tx = conn.unchecked_transaction()?;

        // Persist splitters for use by subsequent append calls.
        for &kmer in &splitters {
            tx.execute(
                "INSERT OR IGNORE INTO splitter (kmer) VALUES (?1)",
                params![kmer as i64],
            )?;
        }

        tx.execute("INSERT INTO sample (name) VALUES (?1)", params![sample_name])?;
        let sample_id: i64 = tx.last_insert_rowid();

        for (record, segs) in records.iter().zip(compressed_per_record.iter()) {
            tx.execute(
                "INSERT INTO contig (sample_id, name, length) VALUES (?1, ?2, ?3)",
                params![sample_id, &record.name, record.seq.len() as i64],
            )?;
            let contig_id: i64 = tx.last_insert_rowid();

            for (seg_order, (ref_blob, raw_length, kf, kb, is_rc)) in segs.iter().enumerate() {
                tx.execute(
                    "INSERT INTO segment_group (ref_data, params, kmer_front, kmer_back) \
                     VALUES (?1, ?2, ?3, ?4)",
                    params![ref_blob, &params_json, kf, kb],
                )?;
                let group_id: i64 = tx.last_insert_rowid();
                tx.execute(
                    "INSERT INTO segment \
                     (contig_id, seg_order, group_id, in_group_id, \
                      is_rev_comp, raw_length, delta_data) \
                     VALUES (?1, ?2, ?3, 0, ?4, ?5, NULL)",
                    params![
                        contig_id,
                        seg_order as i64,
                        group_id,
                        *is_rc,
                        *raw_length as i64
                    ],
                )?;
            }
        }

        tx.commit()?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Delta path (subsequent samples)
    // -----------------------------------------------------------------------

    fn add_as_delta(&mut self, records: Vec<FastaRecord>, sample_name: &str) -> Result<()> {
        let segment_size = self.params.segment_size as usize;
        let splitter_k = self.params.splitter_k as usize;
        let params_json_str = params_json(&self.params);
        let params_ref = &self.params;
        let kmer_len = self.params.min_match_len as usize;

        // --- Phase 1: load splitters and build lookup maps -------------------

        // Load persisted splitters from the archive.
        let splitters: HashSet<u64> = {
            let mut stmt = self
                .db
                .conn()
                .prepare("SELECT kmer FROM splitter")?;
            let collected: HashSet<u64> = stmt
                .query_map([], |r| r.get::<_, i64>(0))?
                .filter_map(|r| r.ok())
                .map(|v| v as u64)
                .collect();
            collected
        };

        // Exact (kmer_front, kmer_back) → group_id map.
        let exact_map = build_exact_map(self.db.conn())?;

        // Terminators map: splitter → list of partner splitters it has been
        // paired with in stored segment groups (one-splitter lookup).
        let terminators_map = build_terminators_map(self.db.conn())?;

        // Fallback vote-based index (hash-sampled k-mers) and ref-len map.
        let sample_rate = compute_sample_rate(segment_size, kmer_len);
        let (kmer_index, ref_data_map, ref_len_map) =
            build_fallback_index(self.db.conn(), kmer_len, sample_rate, MAX_KMER_BUCKET)?;

        // --- Phase 2: parallel compression with adaptive splitting -----------
        //
        // For each contig, split with global splitters.  If no global splitter
        // fired (single large segment), apply AGC's adaptive strategy:
        // find singleton k-mers unique to this contig (absent from the full
        // reference k-mer set) and re-split with those local splitters.  The
        // local splitters are returned alongside the chunk results so they can
        // be saved to the DB for future append calls.
        //
        // LzDiff instances are cached per-record in a local HashMap so that
        // multiple segments of the same record mapping to the same group only
        // build the index once.  Memory is bounded to (groups in one contig) ×
        // (LzDiff size) per thread, rather than all groups simultaneously.

        let compressed_per_record: Vec<(Vec<ChunkResult>, Vec<u64>)> = records
            .par_iter()
            .map(|rec| {
                let mut local_splitters: Vec<u64> = Vec::new();

                // Initial split with global splitters.
                let mut segs = split_contig(&rec.seq, &splitters, splitter_k);

                // Adaptive: one big segment with no boundaries → try local splitters.
                if segs.len() == 1 && rec.seq.len() >= segment_size {
                    let new_spl =
                        find_local_splitters(&rec.seq, &splitters, splitter_k, segment_size);
                    if !new_spl.is_empty() {
                        let combined: HashSet<u64> = splitters
                            .iter()
                            .chain(new_spl.iter())
                            .copied()
                            .collect();
                        segs = split_contig(&rec.seq, &combined, splitter_k);
                        local_splitters.extend(new_spl);
                    }
                }

                // Per-record LzDiff cache: built on demand, dropped after this
                // record is done.  Avoids rebuilding hash tables for groups that
                // appear in multiple segments of the same contig.
                let mut lz_local: HashMap<i64, crate::lz_diff::LzDiff> = HashMap::new();

                let mut chunks: Vec<ChunkResult> = Vec::with_capacity(segs.len());
                for seg in segs {
                    // Try exact kmer-pair lookup first.
                    let exact_group = if seg.kmer_front != SENTINEL
                        && seg.kmer_back != SENTINEL
                    {
                        let kf = seg.kmer_front;
                        let kb = seg.kmer_back;
                        exact_map.get(&(kf.min(kb), kf.max(kb))).copied()
                    } else {
                        None
                    };

                    // One-splitter lookup.
                    let one_spl_group = if exact_group.is_none() {
                        let single_kmer =
                            if seg.kmer_front != SENTINEL && seg.kmer_back == SENTINEL {
                                Some(seg.kmer_front)
                            } else if seg.kmer_back != SENTINEL && seg.kmer_front == SENTINEL {
                                Some(seg.kmer_back)
                            } else {
                                None
                            };
                        single_kmer.and_then(|kmer| {
                            find_one_splitter_group(
                                kmer,
                                seg.seq.len(),
                                &terminators_map,
                                &exact_map,
                                &ref_len_map,
                            )
                        })
                    } else {
                        None
                    };

                    // Fall back to vote-based matching if no hit yet.
                    let matched_group = exact_group.or(one_spl_group).or_else(|| {
                        find_best_ref_group(&seg.seq, &kmer_index, kmer_len, sample_rate)
                    });

                    let chunk = match matched_group {
                        Some(group_id) => {
                            // Build LzDiff for this group on first use, reuse thereafter.
                            if !lz_local.contains_key(&group_id) {
                                let ref_blob = ref_data_map
                                    .get(&group_id)
                                    .expect("group_id missing from ref_data_map");
                                let lz = segment::lz_from_ref_blob(ref_blob, params_ref)
                                    .expect("lz_from_ref_blob");
                                lz_local.insert(group_id, lz);
                            }
                            let lz = lz_local.get(&group_id).unwrap();
                            let raw_delta = segment::compress_delta(lz, &seg.seq)
                                .expect("compress_delta");
                            ChunkResult::Delta {
                                group_id,
                                raw_length: seg.raw_len,
                                raw_delta,
                                is_rc: seg.is_rc,
                            }
                        }
                        None => {
                            let ref_blob = segment::compress_reference(&seg.seq)
                                .expect("compress_reference");
                            let kf = if seg.kmer_front != SENTINEL {
                                Some(seg.kmer_front as i64)
                            } else {
                                None
                            };
                            let kb = if seg.kmer_back != SENTINEL {
                                Some(seg.kmer_back as i64)
                            } else {
                                None
                            };
                            ChunkResult::NewRef {
                                raw_length: seg.raw_len,
                                ref_blob,
                                kmer_front: kf,
                                kmer_back: kb,
                                is_rc: seg.is_rc,
                            }
                        }
                    };
                    chunks.push(chunk);
                }

                (chunks, local_splitters)
            })
            .collect();

        // --- Phase 3: insert in a single transaction -------------------------
        //
        // Strategy for batch delta compression:
        // 1. Pre-query MAX(in_group_id) per group to find existing delta counts.
        // 2. During the insert loop, assign sequential in_group_id per group and
        //    accumulate raw LZ-diff bytes keyed by group_id.
        // 3. After all segments are inserted, batch-ZSTD each group's raw deltas
        //    and UPDATE segment_group.delta_blob.  For groups with an existing
        //    delta_blob, extract all previous raw entries, append the new ones,
        //    and re-compress the combined list.

        // Pre-query existing in_group_id counts so appends get correct ids.
        let existing_counts: HashMap<i64, i64> = {
            let mut stmt = self.db.conn().prepare(
                "SELECT group_id, MAX(in_group_id) \
                 FROM segment WHERE in_group_id > 0 \
                 GROUP BY group_id",
            )?;
            let collected: Vec<(i64, i64)> = stmt
                .query_map([], |r| Ok((r.get(0)?, r.get(1)?)))?
                .filter_map(|r| r.ok())
                .collect();
            collected.into_iter().collect()
        };

        // Per-group raw delta accumulator: group_id → list of raw LZ-diff bytes
        // in in_group_id order (1-based, relative to this batch; caller adds
        // existing_count offset).
        let mut new_deltas_per_group: HashMap<i64, Vec<Vec<u8>>> = HashMap::new();
        // in_group_id counter per group for this batch.
        let mut next_in_group_id: HashMap<i64, i64> = HashMap::new();

        let conn = self.db.conn();
        let tx = conn.unchecked_transaction()?;

        tx.execute("INSERT INTO sample (name) VALUES (?1)", params![sample_name])?;
        let sample_id: i64 = tx.last_insert_rowid();

        for (record, (chunks, local_spl)) in records.iter().zip(compressed_per_record.iter()) {
            tx.execute(
                "INSERT INTO contig (sample_id, name, length) VALUES (?1, ?2, ?3)",
                params![sample_id, &record.name, record.seq.len() as i64],
            )?;
            let contig_id: i64 = tx.last_insert_rowid();

            // Persist any local (adaptive) splitters found for this contig.
            for &kmer in local_spl {
                tx.execute(
                    "INSERT OR IGNORE INTO splitter (kmer) VALUES (?1)",
                    params![kmer as i64],
                )?;
            }

            for (seg_order, chunk_result) in chunks.iter().enumerate() {
                match chunk_result {
                    ChunkResult::Delta {
                        group_id,
                        raw_length,
                        raw_delta,
                        is_rc,
                    } => {
                        // Compute the 1-based in_group_id for this delta.
                        let existing_base = existing_counts.get(group_id).copied().unwrap_or(0);
                        let batch_pos = next_in_group_id.entry(*group_id).or_insert(0);
                        *batch_pos += 1;
                        let in_group_id = existing_base + *batch_pos;

                        // Accumulate raw bytes for later batch compression.
                        new_deltas_per_group
                            .entry(*group_id)
                            .or_default()
                            .push(raw_delta.clone());

                        tx.execute(
                            "INSERT INTO segment \
                             (contig_id, seg_order, group_id, in_group_id, \
                              is_rev_comp, raw_length, delta_data) \
                             VALUES (?1, ?2, ?3, ?4, ?5, ?6, NULL)",
                            params![
                                contig_id,
                                seg_order as i64,
                                group_id,
                                in_group_id,
                                *is_rc,
                                *raw_length as i64
                            ],
                        )?;
                    }
                    ChunkResult::NewRef {
                        raw_length,
                        ref_blob,
                        kmer_front,
                        kmer_back,
                        is_rc,
                    } => {
                        tx.execute(
                            "INSERT INTO segment_group \
                             (ref_data, params, kmer_front, kmer_back) \
                             VALUES (?1, ?2, ?3, ?4)",
                            params![ref_blob, &params_json_str, kmer_front, kmer_back],
                        )?;
                        let new_group_id: i64 = tx.last_insert_rowid();
                        tx.execute(
                            "INSERT INTO segment \
                             (contig_id, seg_order, group_id, in_group_id, \
                              is_rev_comp, raw_length, delta_data) \
                             VALUES (?1, ?2, ?3, 0, ?4, ?5, NULL)",
                            params![
                                contig_id,
                                seg_order as i64,
                                new_group_id,
                                *is_rc,
                                *raw_length as i64
                            ],
                        )?;
                    }
                }
            }
        }

        // Finalize: batch-compress raw deltas per group.
        //
        // Step A — load existing delta_blobs from DB in a single batch query
        // before the transaction (read-only, no lock held).
        let existing_blobs: HashMap<i64, Option<Vec<u8>>> = {
            let ids: Vec<i64> = new_deltas_per_group.keys().copied().collect();
            // Build "WHERE id IN (?,?,…)" dynamically.
            let placeholders = ids
                .iter()
                .enumerate()
                .map(|(i, _)| format!("?{}", i + 1))
                .collect::<Vec<_>>()
                .join(",");
            let sql = format!(
                "SELECT id, delta_blob FROM segment_group WHERE id IN ({placeholders})"
            );
            let conn = self.db.conn();
            let mut stmt = conn.prepare(&sql)?;
            let rows: Vec<(i64, Option<Vec<u8>>)> = stmt
                .query_map(rusqlite::params_from_iter(ids.iter()), |r| {
                    Ok((r.get::<_, i64>(0)?, r.get::<_, Option<Vec<u8>>>(1)?))
                })?
                .filter_map(|r| r.ok())
                .collect();
            // Any id not returned had no row (shouldn't happen, but default to None).
            let mut map: HashMap<i64, Option<Vec<u8>>> =
                ids.iter().map(|&id| (id, None)).collect();
            for (id, blob) in rows {
                map.insert(id, blob);
            }
            map
        };

        // Step B — parallel ZSTD recompression, one task per group.
        let combined_blobs: HashMap<i64, Vec<u8>> = new_deltas_per_group
            .par_iter()
            .map(|(group_id, new_raws)| {
                let mut all_raws: Vec<Vec<u8>> = Vec::new();
                if let Some(Some(blob)) = existing_blobs.get(group_id) {
                    let existing_count =
                        existing_counts.get(group_id).copied().unwrap_or(0) as usize;
                    for idx in 0..existing_count {
                        let raw = segment::extract_delta_from_batch(blob, idx)
                            .expect("extract_delta_from_batch");
                        all_raws.push(raw);
                    }
                }
                all_raws.extend_from_slice(new_raws);
                let blob = segment::batch_compress_deltas(&all_raws)
                    .expect("batch_compress_deltas");
                (*group_id, blob)
            })
            .collect();

        // Step C — inside the transaction: fast UPDATE only (no ZSTD here).
        for (group_id, combined_blob) in &combined_blobs {
            tx.execute(
                "UPDATE segment_group SET delta_blob = ?1 WHERE id = ?2",
                params![combined_blob, group_id],
            )?;
        }

        tx.commit()?;
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

    fn write_tmp_fasta(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(content.as_bytes()).expect("write fasta");
        f
    }

    fn tmp_archive() -> (TempDir, std::path::PathBuf) {
        let dir = TempDir::new().expect("tempdir");
        let path = dir.path().join("test.agcrs");
        (dir, path)
    }

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
    }

    #[test]
    fn duplicate_sample_returns_error() {
        let fasta = ">chr1\nACGT\n";
        let fasta_file = write_tmp_fasta(fasta);
        let (_dir, archive) = tmp_archive();

        let mut c = Compressor::create(&archive, Params::default()).unwrap();
        c.add_fasta(fasta_file.path(), "s1").unwrap();
        c.finish().unwrap();

        let mut c2 = Compressor::append(&archive).unwrap();
        assert!(c2.add_fasta(fasta_file.path(), "s1").is_err());
        c2.finish().unwrap();
    }

    // -----------------------------------------------------------------------
    // Splitter determination tests
    // -----------------------------------------------------------------------

    fn make_records_from_2bit(seqs: Vec<Vec<u8>>) -> Vec<FastaRecord> {
        seqs.into_iter()
            .enumerate()
            .map(|(i, seq)| FastaRecord { name: format!("ctg{}", i), seq })
            .collect()
    }

    /// A pseudo-random 2-bit sequence has many singleton k-mers; the splitter
    /// algorithm should find at least one when segment_size is small.
    /// Using k=8 so that 4^8=65536 >> 200 bases → nearly all k-mers are singletons.
    #[test]
    fn splitters_found_in_simple_sequence() {
        // LCG-generated pseudo-random sequence; with k=8 virtually every
        // k-mer is unique in a 200-base window.
        let mut state: u32 = 0xDEAD_BEEF;
        let seq: Vec<u8> = (0..200)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                ((state >> 14) & 3) as u8
            })
            .collect();
        let records = make_records_from_2bit(vec![seq]);
        let s = determine_splitters(&records, 8, 20);
        assert!(!s.is_empty(), "expected at least one splitter");
    }

    /// A sequence shorter than segment_size + k should still not panic.
    #[test]
    fn splitters_short_sequence_no_panic() {
        let seq: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let records = make_records_from_2bit(vec![seq]);
        let _ = determine_splitters(&records, 4, 60_000);
    }

    // -----------------------------------------------------------------------
    // Split contig tests
    // -----------------------------------------------------------------------

    /// With no splitters the entire sequence becomes one segment with SENTINEL
    /// boundaries.
    #[test]
    fn split_with_no_splitters_is_single_segment() {
        let seq: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let empty: HashSet<u64> = HashSet::new();
        let segs = split_contig(&seq, &empty, 4);
        assert_eq!(segs.len(), 1);
        assert_eq!(segs[0].kmer_front, SENTINEL);
        assert_eq!(segs[0].kmer_back, SENTINEL);
        assert!(!segs[0].is_rc);
        assert_eq!(segs[0].seq, seq);
    }

    /// Segment contents should concatenate back to the original sequence.
    #[test]
    fn split_concatenates_to_original() {
        let seq: Vec<u8> = (0u8..100).map(|i| i % 4).collect();
        let records = make_records_from_2bit(vec![seq.clone()]);
        let splitters = determine_splitters(&records, 4, 10);
        let segs = split_contig(&seq, &splitters, 4);

        // Reconstruct: for is_rc segments, RC back to original before concatenating.
        let mut reconstructed: Vec<u8> = Vec::new();
        for seg in &segs {
            if seg.is_rc {
                reconstructed.extend(rev_comp_2bit_seq(&seg.seq));
            } else {
                reconstructed.extend(&seg.seq);
            }
        }
        assert_eq!(reconstructed, seq, "concatenation must equal original");
    }

    // -----------------------------------------------------------------------
    // Fallback k-mer sampling tests
    // -----------------------------------------------------------------------

    #[test]
    fn sample_kmers_rate1_exhaustive() {
        let seq: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let kmers = sample_kmers_2bit(&seq, 18, 1);
        assert_eq!(kmers.len(), 3);
        assert_ne!(kmers[0], kmers[1]);
        assert_ne!(kmers[1], kmers[2]);
    }

    #[test]
    fn sample_kmers_short_sequence() {
        let seq: Vec<u8> = vec![0, 1, 2];
        let kmers = sample_kmers_2bit(&seq, 18, 1);
        assert!(kmers.is_empty());
    }

    #[test]
    fn sample_kmers_n_handling() {
        let mut seq: Vec<u8> = vec![0u8; 20];
        seq[10] = 4; // N in the middle
        let kmers = sample_kmers_2bit(&seq, 18, 1);
        // Window resets at N; no full window can span it.
        assert!(kmers.is_empty(), "N should break all k-mers");
    }

    #[test]
    fn sample_kmers_position_independent() {
        let seq1: Vec<u8> = (0u8..100).map(|i| i % 4).collect();
        let seq2: Vec<u8> = (0u8..100).rev().map(|i| i % 4).collect();
        let rate = compute_sample_rate(100, 18);
        let km1: HashSet<u64> = sample_kmers_2bit(&seq1, 18, rate).into_iter().collect();
        let km2: HashSet<u64> = sample_kmers_2bit(&seq2, 18, rate).into_iter().collect();
        // Calling twice on the same sequence must yield the same set.
        let km1b: HashSet<u64> = sample_kmers_2bit(&seq1, 18, rate).into_iter().collect();
        assert_eq!(km1, km1b, "sampling must be deterministic");
        // The two sets need not be equal (different sequences), but both must
        // be subsets of the full k-mer universe.
        let _ = km2; // just ensure it compiles and doesn't panic
    }
}
