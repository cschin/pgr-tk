use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

use rusqlite::params;

use crate::db::AgcDb;
use crate::error::{AgcError, Result};

/// Merge multiple agc-rs archives into a single output archive.
///
/// All source archives must have been built against the **same reference**
/// (identical splitter sets).  The reference sample and all non-reference
/// samples from `sources[0]` are taken as the canonical base.  Each
/// subsequent archive contributes only its non-reference samples (the first /
/// lowest-id sample in each archive is treated as the reference and skipped).
///
/// Segment groups are deduplicated by `(kmer_front, kmer_back)`: if the same
/// genomic region boundary pair is already present in the output, the existing
/// group is reused rather than duplicated.  Groups with NULL boundary kmers
/// (terminal / fallback-matched segments) have no reliable identity key and
/// are always copied as new groups.
///
/// No recompression is performed: `delta_data` blobs are copied byte-for-byte.
pub fn merge_archives(output_path: &Path, sources: &[&Path]) -> Result<()> {
    if sources.is_empty() {
        return Ok(());
    }

    // Trivial case: one source → just copy the file.
    if sources.len() == 1 {
        std::fs::copy(sources[0], output_path)?;
        return Ok(());
    }

    let t0 = Instant::now();
    eprintln!(
        "[{:7.2}s] merge: {} archives → {}",
        t0.elapsed().as_secs_f64(),
        sources.len(),
        output_path.display(),
    );

    // Open all sources read-only.
    let src_dbs: Vec<AgcDb> = sources
        .iter()
        .map(|p| AgcDb::open_readonly(p))
        .collect::<Result<Vec<_>>>()?;

    // --- Validate splitter sets are compatible and collect their union --------
    //
    // Two archives built from the same reference genome share the same base
    // splitter set.  However, appending novel (highly-divergent) sequences may
    // add local splitters that are unique to each archive.  We therefore allow
    // archives that share a large common core (Jaccard ≥ 0.90) and warn when
    // their sets differ.  The output archive receives the *union* of all
    // splitter sets so that subsequent appends can use any boundary k-mer that
    // was valid in any source.
    //
    // A Jaccard < 0.90 almost certainly means the archives were built from
    // different reference genomes, which would corrupt the merged output.
    eprintln!(
        "[{:7.2}s]   validating splitter sets ...",
        t0.elapsed().as_secs_f64()
    );

    fn load_splitters(db: &AgcDb) -> Result<std::collections::HashSet<i64>> {
        let mut stmt = db
            .conn()
            .prepare("SELECT kmer FROM splitter")?;
        let collected: std::collections::HashSet<i64> = stmt
            .query_map([], |r| r.get::<_, i64>(0))?
            .filter_map(|r| r.ok())
            .collect();
        Ok(collected)
    }

    let base_splitters = load_splitters(&src_dbs[0])?;
    let mut union_splitters: std::collections::HashSet<i64> = base_splitters.clone();

    for (i, src) in src_dbs.iter().enumerate().skip(1) {
        let splitters = load_splitters(src)?;
        let common = base_splitters.intersection(&splitters).count();
        let total  = base_splitters.union(&splitters).count();
        let jaccard = if total == 0 { 1.0_f64 } else { common as f64 / total as f64 };

        if jaccard < 0.90 {
            return Err(AgcError::SampleNotFound(format!(
                "source {} has an incompatible splitter set \
                 (Jaccard similarity with source 1 = {:.3}; expected ≥ 0.90) — \
                 archives do not appear to share the same reference genome",
                sources[i].display(),
                jaccard
            )));
        }

        let novel_in_src = splitters.difference(&base_splitters).count();
        let novel_in_base = base_splitters.difference(&splitters).count();
        if novel_in_src > 0 || novel_in_base > 0 {
            eprintln!(
                "[{:7.2}s]   note: source {} has {} novel splitters not in source 1 \
                 and source 1 has {} not in source {} (Jaccard={:.4}); \
                 these are likely from locally-split novel contigs — merging union",
                t0.elapsed().as_secs_f64(),
                i + 1,
                novel_in_src,
                novel_in_base,
                i + 1,
                jaccard
            );
        }

        union_splitters.extend(&splitters);
    }

    let ref_splitters: Vec<i64> = {
        let mut v: Vec<i64> = union_splitters.into_iter().collect();
        v.sort_unstable();
        v
    };

    eprintln!(
        "[{:7.2}s]   splitters compatible ({} kmers in union)",
        t0.elapsed().as_secs_f64(),
        ref_splitters.len(),
    );

    // --- Validate sample names are unique across all sources -----------------
    // Collect the names that will actually be written to the output:
    //   source 0 → all samples
    //   source 1..N → all samples except the first (reference)
    // Report every conflict before bailing so the user sees the full picture.
    eprintln!(
        "[{:7.2}s]   checking sample name uniqueness ...",
        t0.elapsed().as_secs_f64()
    );
    {
        let mut seen: HashMap<String, usize> = HashMap::new(); // name → first source index
        let mut conflicts: Vec<String> = Vec::new();
        for (src_idx, src_db) in src_dbs.iter().enumerate() {
            let names: Vec<String> = {
                let query = if src_idx == 0 {
                    "SELECT name FROM sample ORDER BY id"
                } else {
                    "SELECT name FROM sample ORDER BY id LIMIT -1 OFFSET 1"
                };
                let mut stmt = src_db.conn().prepare(query)?;
                let collected: Vec<String> = stmt
                    .query_map([], |r| r.get::<_, String>(0))?
                    .filter_map(|r| r.ok())
                    .collect();
                collected
            };
            for name in names {
                if let Some(&prev_src) = seen.get(&name) {
                    conflicts.push(format!(
                        "  sample '{}' appears in both source {} ({}) and source {} ({})",
                        name,
                        prev_src + 1,
                        sources[prev_src].display(),
                        src_idx + 1,
                        sources[src_idx].display(),
                    ));
                } else {
                    seen.insert(name, src_idx);
                }
            }
        }
        if !conflicts.is_empty() {
            return Err(AgcError::SampleNotFound(format!(
                "duplicate sample names across archives — cannot merge:\n{}",
                conflicts.join("\n")
            )));
        }
    }
    eprintln!(
        "[{:7.2}s]   sample names are unique",
        t0.elapsed().as_secs_f64()
    );

    // --- Create output archive -----------------------------------------------
    let out_db = AgcDb::create(output_path)?;

    // Write splitters.
    {
        let tx = out_db.conn().unchecked_transaction()?;
        for kmer in &ref_splitters {
            tx.execute(
                "INSERT OR IGNORE INTO splitter (kmer) VALUES (?1)",
                params![kmer],
            )?;
        }
        tx.commit()?;
    }

    // `kmer_key_to_out_group`: (kmer_front, kmer_back) → group_id in output.
    // Populated as we copy segment_groups; used to deduplicate groups from
    // subsequent sources that map to the same boundary pair.
    let mut kmer_key_to_out_group: HashMap<(i64, i64), i64> = HashMap::new();

    // --- Process each source -------------------------------------------------
    for (src_idx, src_db) in src_dbs.iter().enumerate() {
        eprintln!(
            "[{:7.2}s]   source {}/{}: {}",
            t0.elapsed().as_secs_f64(),
            src_idx + 1,
            sources.len(),
            sources[src_idx].display(),
        );

        // Build a group-id remapping for this source:
        //   src_group_id → out_group_id
        //
        // We preload ALL segment_groups from the source so that segment rows
        // can be rewritten with output-relative group IDs in O(1) per segment.
        let src_to_out_group =
            build_group_mapping(src_db, &out_db, &mut kmer_key_to_out_group, t0)?;

        // Determine which samples to copy from this source.
        // For the first source: ALL samples (including the reference).
        // For subsequent sources: skip the reference (the sample with the
        // smallest id, which was the first one inserted when the archive was
        // created).
        let all_samples: Vec<(i64, String)> = {
            let mut stmt = src_db
                .conn()
                .prepare("SELECT id, name FROM sample ORDER BY id")?;
            let collected: Vec<(i64, String)> = stmt
                .query_map([], |r| Ok((r.get::<_, i64>(0)?, r.get::<_, String>(1)?)))?
                .filter_map(|r| r.ok())
                .collect();
            collected
        };

        let samples_to_copy: Vec<(i64, String)> = if src_idx == 0 {
            all_samples
        } else {
            // Skip the first (reference) sample.
            all_samples.into_iter().skip(1).collect()
        };

        eprintln!(
            "[{:7.2}s]     copying {} samples ...",
            t0.elapsed().as_secs_f64(),
            samples_to_copy.len(),
        );

        let tx = out_db.conn().unchecked_transaction()?;
        for (src_sample_id, sample_name) in &samples_to_copy {
            tx.execute(
                "INSERT INTO sample (name) VALUES (?1)",
                params![sample_name],
            )?;
            let out_sample_id: i64 = tx.last_insert_rowid();

            // Copy contigs.
            let contigs: Vec<(i64, String, i64)> = {
                let mut stmt = src_db.conn().prepare(
                    "SELECT id, name, length FROM contig WHERE sample_id = ?1 ORDER BY id",
                )?;
                let collected: Vec<(i64, String, i64)> = stmt
                    .query_map(params![src_sample_id], |r| {
                        Ok((
                            r.get::<_, i64>(0)?,
                            r.get::<_, String>(1)?,
                            r.get::<_, i64>(2)?,
                        ))
                    })?
                    .filter_map(|r| r.ok())
                    .collect();
                collected
            };

            for (src_contig_id, contig_name, contig_len) in &contigs {
                tx.execute(
                    "INSERT INTO contig (sample_id, name, length) VALUES (?1, ?2, ?3)",
                    params![out_sample_id, contig_name, contig_len],
                )?;
                let out_contig_id: i64 = tx.last_insert_rowid();

                // Copy segments with remapped group_id and contig_id.
                let segments: Vec<(i64, i64, i64, bool, i64, Option<Vec<u8>>)> = {
                    let mut stmt = src_db.conn().prepare(
                        "SELECT seg_order, group_id, in_group_id, is_rev_comp, \
                                raw_length, delta_data \
                         FROM segment WHERE contig_id = ?1 ORDER BY seg_order",
                    )?;
                    let collected: Vec<(i64, i64, i64, bool, i64, Option<Vec<u8>>)> = stmt
                        .query_map(params![src_contig_id], |r| {
                            Ok((
                                r.get::<_, i64>(0)?,
                                r.get::<_, i64>(1)?,
                                r.get::<_, i64>(2)?,
                                r.get::<_, bool>(3)?,
                                r.get::<_, i64>(4)?,
                                r.get::<_, Option<Vec<u8>>>(5)?,
                            ))
                        })?
                        .filter_map(|r| r.ok())
                        .collect();
                    collected
                };

                for (seg_order, src_group_id, in_group_id, is_rc, raw_len, delta_data) in &segments
                {
                    let out_group_id = src_to_out_group
                        .get(src_group_id)
                        .copied()
                        .unwrap_or(*src_group_id);

                    tx.execute(
                        "INSERT INTO segment \
                         (contig_id, seg_order, group_id, in_group_id, \
                          is_rev_comp, raw_length, delta_data) \
                         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
                        params![
                            out_contig_id,
                            seg_order,
                            out_group_id,
                            in_group_id,
                            is_rc,
                            raw_len,
                            delta_data,
                        ],
                    )?;
                }
            }
        }
        tx.commit()?;

        eprintln!(
            "[{:7.2}s]     done ({} segment_group mappings)",
            t0.elapsed().as_secs_f64(),
            src_to_out_group.len(),
        );
    }

    eprintln!("[{:7.2}s] merge complete", t0.elapsed().as_secs_f64());
    Ok(())
}

// ---------------------------------------------------------------------------
// Helper: build segment_group id mapping for one source
// ---------------------------------------------------------------------------

/// Load all segment_groups from `src_db`, copy any that are not yet in
/// `out_db` (as determined by `kmer_key_to_out_group`), and return a
/// `HashMap<src_group_id, out_group_id>` covering every group in the source.
fn build_group_mapping(
    src_db: &AgcDb,
    out_db: &AgcDb,
    kmer_key_to_out_group: &mut HashMap<(i64, i64), i64>,
    t0: Instant,
) -> Result<HashMap<i64, i64>> {
    // Preload all groups from source.
    let src_groups: Vec<(i64, Vec<u8>, String, Option<i64>, Option<i64>)> = {
        let mut stmt = src_db.conn().prepare(
            "SELECT id, ref_data, params, kmer_front, kmer_back \
             FROM segment_group ORDER BY id",
        )?;
        let collected: Vec<(i64, Vec<u8>, String, Option<i64>, Option<i64>)> = stmt
            .query_map([], |r| {
                Ok((
                    r.get::<_, i64>(0)?,
                    r.get::<_, Vec<u8>>(1)?,
                    r.get::<_, String>(2)?,
                    r.get::<_, Option<i64>>(3)?,
                    r.get::<_, Option<i64>>(4)?,
                ))
            })?
            .filter_map(|r| r.ok())
            .collect();
        collected
    };

    eprintln!(
        "[{:7.2}s]     mapping {} segment_groups ...",
        t0.elapsed().as_secs_f64(),
        src_groups.len(),
    );

    let mut src_to_out: HashMap<i64, i64> = HashMap::with_capacity(src_groups.len());

    let tx = out_db.conn().unchecked_transaction()?;

    for (src_gid, ref_data, grp_params, kf, kb) in src_groups {
        // Can we deduplicate by (kmer_front, kmer_back)?
        if let (Some(kf_val), Some(kb_val)) = (kf, kb) {
            let key = (kf_val.min(kb_val), kf_val.max(kb_val));
            if let Some(&out_gid) = kmer_key_to_out_group.get(&key) {
                // Already in output — reuse without copying.
                src_to_out.insert(src_gid, out_gid);
                continue;
            }
            // Not yet in output — copy and register.
            tx.execute(
                "INSERT INTO segment_group (ref_data, params, kmer_front, kmer_back) \
                 VALUES (?1, ?2, ?3, ?4)",
                params![ref_data, grp_params, kf, kb],
            )?;
            let out_gid: i64 = tx.last_insert_rowid();
            kmer_key_to_out_group.insert(key, out_gid);
            src_to_out.insert(src_gid, out_gid);
        } else {
            // NULL-keyed group: no reliable dedup key — always copy as new.
            tx.execute(
                "INSERT INTO segment_group (ref_data, params, kmer_front, kmer_back) \
                 VALUES (?1, ?2, ?3, ?4)",
                params![ref_data, grp_params, kf, kb],
            )?;
            let out_gid: i64 = tx.last_insert_rowid();
            src_to_out.insert(src_gid, out_gid);
        }
    }

    tx.commit()?;
    Ok(src_to_out)
}
