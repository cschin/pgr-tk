use axum::{extract::State, Json};
use rusqlite::{Connection, OpenFlags};
use std::sync::Arc;

use crate::error::AppError;
use crate::models::{
    GeneInfo, GeneQueryRequest, GeneQueryResponse, HitRecord, QueryRegionInfo, RegionQueryRequest,
    RegionQueryResponse,
};
use crate::state::{AppState, SyncSeqIndexDB};

// ── helpers ──────────────────────────────────────────────────────────────────

/// Find the first reference sequence ID whose source contains `hint` and
/// whose contig name equals `chrom`.
fn find_ref_sid(
    db: &SyncSeqIndexDB,
    chrom: &str,
    hint: &str,
) -> Option<(u32, u32)> {
    let seq_info = db.0.seq_info.as_ref()?;
    for (&sid, (ctg, src, length)) in seq_info {
        if ctg == chrom {
            let src_str = src.as_deref().unwrap_or("");
            if src_str.contains(hint) {
                return Some((sid, *length));
            }
        }
    }
    None
}

/// Run a shimmer query against the database and convert raw hit-pair chains
/// into `HitRecord`s ready for the API response.
fn run_query(
    db: Arc<SyncSeqIndexDB>,
    query_seq: Vec<u8>,
    max_count: u32,
    max_query_count: u32,
    max_target_count: u32,
    merge_range_tol: usize,
    min_anchor_count: usize,
) -> anyhow::Result<Vec<HitRecord>> {
    let query_results = db
        .0
        .query_fragment_to_hps_from_mmap_file(
            &query_seq,
            0.025,
            Some(max_count),
            Some(max_query_count),
            Some(max_target_count),
            Some(8),
            None,
            false,
        )
        .ok_or_else(|| anyhow::anyhow!("shimmer query returned None (index not loaded?)"))?;

    let seq_info = db
        .0
        .seq_info
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("seq_info not loaded"))?;

    let mut hits: Vec<HitRecord> = Vec::new();

    for (sid, alns) in query_results {
        // Collect per-(sid, orientation) ranges
        let mut fwd_ranges: Vec<(u32, u32, usize)> = Vec::new(); // (bgn, end, anchor_count)
        let mut rev_ranges: Vec<(u32, u32, usize)> = Vec::new();

        for (_score, aln) in alns {
            if aln.len() < min_anchor_count {
                continue;
            }
            let mut f_count = 0usize;
            let mut r_count = 0usize;
            for hp in &aln {
                if hp.0 .2 == hp.1 .2 {
                    f_count += 1;
                } else {
                    r_count += 1;
                }
            }
            let orientation = if f_count >= r_count { 0u32 } else { 1u32 };

            let mut target_coords: Vec<(u32, u32)> = aln
                .iter()
                .map(|hp| (hp.1 .0, hp.1 .1))
                .collect();
            target_coords.sort();
            let bgn = target_coords[0].0;
            let end = target_coords[target_coords.len() - 1].1;
            let anchor_count = aln.len();

            if orientation == 0 {
                fwd_ranges.push((bgn, end, anchor_count));
            } else {
                rev_ranges.push((bgn, end, anchor_count));
            }
        }

        let (ctg, src, _ctg_len) = match seq_info.get(&sid) {
            Some(v) => v,
            None => continue,
        };
        let src_str = src.as_deref().unwrap_or("").to_string();

        for (ranges, orientation) in [(&mut fwd_ranges, 0u32), (&mut rev_ranges, 1u32)] {
            if ranges.is_empty() {
                continue;
            }
            ranges.sort_by_key(|r| r.0);
            let mut merged: Vec<(u32, u32, usize)> = Vec::new();
            for (bgn, end, ac) in ranges.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if (bgn as i64) - (last.1 as i64) <= merge_range_tol as i64 {
                        last.1 = last.1.max(end);
                        last.2 = last.2.max(ac);
                        continue;
                    }
                }
                merged.push((bgn, end, ac));
            }
            for (bgn, end, anchor_count) in merged {
                let base = std::path::Path::new(&src_str)
                    .file_stem()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_else(|| src_str.clone());
                let name = format!("{}::{}_{}_{}_{}", base, ctg, bgn, end, orientation);
                hits.push(HitRecord {
                    src: src_str.clone(),
                    ctg: ctg.clone(),
                    bgn,
                    end,
                    orientation,
                    anchor_count,
                    name,
                });
            }
        }
    }

    hits.sort_by(|a, b| {
        a.ctg
            .cmp(&b.ctg)
            .then(a.src.cmp(&b.src))
            .then(a.bgn.cmp(&b.bgn))
    });
    Ok(hits)
}

/// Fetch FASTA sequences for the listed hits.
fn hits_to_fasta(db: &SyncSeqIndexDB, hits: &[HitRecord]) -> anyhow::Result<String> {
    let seq_info = db
        .0
        .seq_info
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("seq_info not loaded"))?;

    // Build (ctg, src) → sid lookup
    let ctg_src_to_sid: std::collections::HashMap<(String, String), u32> = seq_info
        .iter()
        .map(|(&sid, (ctg, src, _))| {
            let src_str = src.as_deref().unwrap_or("").to_string();
            ((ctg.clone(), src_str), sid)
        })
        .collect();

    let mut fasta = String::new();
    for hit in hits {
        let sid = ctg_src_to_sid
            .get(&(hit.ctg.clone(), hit.src.clone()))
            .copied()
            .ok_or_else(|| {
                anyhow::anyhow!("sequence not found: ctg={} src={}", hit.ctg, hit.src)
            })?;
        let mut seq = db
            .0
            .get_sub_seq_by_id(sid, hit.bgn as usize, hit.end as usize)
            .map_err(|e| anyhow::anyhow!("get_sub_seq_by_id: {e}"))?;
        if hit.orientation == 1 {
            seq = pgr_db::fasta_io::reverse_complement(&seq);
        }
        fasta.push('>');
        fasta.push_str(&hit.name);
        fasta.push('\n');
        fasta.push_str(&String::from_utf8_lossy(&seq));
        fasta.push('\n');
    }
    Ok(fasta)
}

// ── /query/region ─────────────────────────────────────────────────────────────

/// Query all haplotype sequences that align to a reference genomic region.
///
/// Fetches the reference subsequence `[start, end)` from the sample identified
/// by `ref_sample_hint`, runs a shimmer-based pangenome alignment, merges
/// overlapping hits, and returns coordinate ranges for every matching haplotype.
/// Set `include_sequences` to `true` to also receive the hit sequences as
/// inline FASTA (may be large for wide regions).
#[utoipa::path(
    post,
    path = "/api/v1/query/region",
    request_body = RegionQueryRequest,
    responses(
        (status = 200, description = "Alignment hits and optional FASTA sequences",
         body = RegionQueryResponse),
        (status = 404, description = "Database or reference chromosome not found",
         body = ErrorBody),
        (status = 500, description = "Shimmer query failed",
         body = ErrorBody),
    ),
    tag = "query"
)]
pub async fn query_region(
    State(state): State<AppState>,
    Json(req): Json<RegionQueryRequest>,
) -> Result<Json<RegionQueryResponse>, AppError> {
    let db_name = req
        .db
        .as_deref()
        .unwrap_or(&state.default_db)
        .to_string();

    let db = state
        .databases
        .get(&db_name)
        .ok_or_else(|| AppError::NotFound(format!("database '{db_name}' not found")))?
        .clone();

    let ref_hint = req
        .ref_sample_hint
        .as_deref()
        .unwrap_or(&db.ref_sample_hint)
        .to_string();

    let max_count = req.max_count.unwrap_or(state.query_defaults.max_count);
    let max_query_count = req
        .max_query_count
        .unwrap_or(state.query_defaults.max_query_count);
    let max_target_count = req
        .max_target_count
        .unwrap_or(state.query_defaults.max_target_count);
    let merge_range_tol = req
        .merge_range_tol
        .unwrap_or(state.query_defaults.merge_range_tol);
    let min_anchor_count = req
        .min_anchor_count
        .unwrap_or(state.query_defaults.min_anchor_count);
    let include_sequences = req.include_sequences;

    let chrom = req.chrom.clone();
    let start = req.start;
    let end = req.end;

    let seq_db = Arc::clone(&db.seq_db);

    let result = tokio::task::spawn_blocking(move || {
        // 1. Find reference SID
        let (ref_sid, ref_len) = find_ref_sid(&seq_db, &chrom, &ref_hint).ok_or_else(|| {
            AppError::NotFound(format!(
                "reference sequence '{chrom}' not found (hint='{ref_hint}')"
            ))
        })?;

        let clamped_start = start.min(ref_len as u64) as usize;
        let clamped_end = end.min(ref_len as u64) as usize;

        // 2. Fetch the query subsequence
        let query_seq = seq_db
            .0
            .get_sub_seq_by_id(ref_sid, clamped_start, clamped_end)
            .map_err(|e| AppError::Internal(format!("get_sub_seq_by_id: {e}")))?;

        // 3. Run shimmer query
        let hits = run_query(
            Arc::clone(&seq_db),
            query_seq,
            max_count,
            max_query_count,
            max_target_count,
            merge_range_tol,
            min_anchor_count,
        )
        .map_err(|e| AppError::Internal(e.to_string()))?;

        // 4. Optionally fetch sequences
        let sequences_fasta = if include_sequences {
            Some(
                hits_to_fasta(&seq_db, &hits)
                    .map_err(|e| AppError::Internal(e.to_string()))?,
            )
        } else {
            None
        };

        Ok::<_, AppError>(RegionQueryResponse {
            query_region: QueryRegionInfo {
                src: ref_hint.clone(),
                ctg: chrom,
                bgn: clamped_start as u64,
                end: clamped_end as u64,
            },
            hits,
            sequences_fasta,
        })
    })
    .await
    .map_err(AppError::from)?;

    Ok(Json(result?))
}

// ── /query/gene ───────────────────────────────────────────────────────────────

fn lookup_gene(gene_db_path: &str, gene_name: &str) -> anyhow::Result<Option<GeneInfo>> {
    let conn = Connection::open_with_flags(
        gene_db_path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )?;

    let mut stmt = conn.prepare(
        "SELECT gene_name, chrom, start, end, strand FROM genes \
         WHERE gene_name = ?1 \
         ORDER BY CASE WHEN chrom NOT GLOB '*_*' THEN 0 ELSE 1 END, end-start DESC \
         LIMIT 1",
    )?;

    let row = stmt.query_row([gene_name], |r| {
        Ok(GeneInfo {
            gene_name: r.get::<_, String>(0)?,
            chrom: r.get::<_, String>(1)?,
            start: r.get::<_, i64>(2)?,
            end: r.get::<_, i64>(3)?,
            strand: r.get::<_, String>(4)?,
        })
    });

    match row {
        Ok(g) => Ok(Some(g)),
        Err(rusqlite::Error::QueryReturnedNoRows) => Ok(None),
        Err(e) => Err(e.into()),
    }
}

/// Look up a gene by symbol, then query all haplotype sequences covering
/// the gene span ± `flank` bp.
///
/// Requires `gene_db` to be configured for the selected database.
/// The gene annotation DB is queried with SQL; primary-assembly chromosomes
/// are preferred over alt/patch contigs; the longest span wins among ties.
#[utoipa::path(
    post,
    path = "/api/v1/query/gene",
    request_body = GeneQueryRequest,
    responses(
        (status = 200, description = "Gene metadata, alignment hits and optional FASTA",
         body = GeneQueryResponse),
        (status = 400, description = "Database has no gene_db configured",
         body = ErrorBody),
        (status = 404, description = "Database or gene symbol not found",
         body = ErrorBody),
        (status = 500, description = "Shimmer query or gene DB error",
         body = ErrorBody),
    ),
    tag = "query"
)]
pub async fn query_gene(
    State(state): State<AppState>,
    Json(req): Json<GeneQueryRequest>,
) -> Result<Json<GeneQueryResponse>, AppError> {
    let db_name = req
        .db
        .as_deref()
        .unwrap_or(&state.default_db)
        .to_string();

    let db = state
        .databases
        .get(&db_name)
        .ok_or_else(|| AppError::NotFound(format!("database '{db_name}' not found")))?
        .clone();

    let gene_db_path = db.gene_db_path.clone().ok_or_else(|| {
        AppError::BadRequest(format!(
            "database '{db_name}' has no gene_db configured"
        ))
    })?;

    let ref_hint = req
        .ref_sample_hint
        .as_deref()
        .unwrap_or(&db.ref_sample_hint)
        .to_string();

    let flank = req.flank.unwrap_or(state.query_defaults.flank);
    let max_count = req.max_count.unwrap_or(state.query_defaults.max_count);
    let max_query_count = req
        .max_query_count
        .unwrap_or(state.query_defaults.max_query_count);
    let max_target_count = req
        .max_target_count
        .unwrap_or(state.query_defaults.max_target_count);
    let merge_range_tol = req
        .merge_range_tol
        .unwrap_or(state.query_defaults.merge_range_tol);
    let min_anchor_count = req
        .min_anchor_count
        .unwrap_or(state.query_defaults.min_anchor_count);
    let include_sequences = req.include_sequences;
    let gene_name = req.gene_name.clone();

    let seq_db = Arc::clone(&db.seq_db);

    let result = tokio::task::spawn_blocking(move || {
        // 1. Look up gene
        let gene = lookup_gene(&gene_db_path, &gene_name)
            .map_err(|e| AppError::Internal(format!("gene lookup: {e}")))?
            .ok_or_else(|| AppError::NotFound(format!("gene '{gene_name}' not found")))?;

        // 2. Build query region with flank
        let region_start = (gene.start as i64 - flank as i64).max(0) as u64;
        let region_end = gene.end as u64 + flank;

        // 3. Find reference SID
        let (ref_sid, ref_len) =
            find_ref_sid(&seq_db, &gene.chrom, &ref_hint).ok_or_else(|| {
                AppError::NotFound(format!(
                    "reference sequence '{}' not found (hint='{ref_hint}')",
                    gene.chrom
                ))
            })?;

        let clamped_start = region_start.min(ref_len as u64) as usize;
        let clamped_end = region_end.min(ref_len as u64) as usize;

        // 4. Fetch query subsequence
        let query_seq = seq_db
            .0
            .get_sub_seq_by_id(ref_sid, clamped_start, clamped_end)
            .map_err(|e| AppError::Internal(format!("get_sub_seq_by_id: {e}")))?;

        // 5. Run shimmer query
        let hits = run_query(
            Arc::clone(&seq_db),
            query_seq,
            max_count,
            max_query_count,
            max_target_count,
            merge_range_tol,
            min_anchor_count,
        )
        .map_err(|e| AppError::Internal(e.to_string()))?;

        // 6. Optionally fetch sequences
        let sequences_fasta = if include_sequences {
            Some(
                hits_to_fasta(&seq_db, &hits)
                    .map_err(|e| AppError::Internal(e.to_string()))?,
            )
        } else {
            None
        };

        Ok::<_, AppError>(GeneQueryResponse {
            gene,
            query_region: QueryRegionInfo {
                src: ref_hint.clone(),
                ctg: gene_name.clone(), // will be overwritten below
                bgn: clamped_start as u64,
                end: clamped_end as u64,
            },
            hits,
            sequences_fasta,
        })
    })
    .await
    .map_err(AppError::from)?;

    // Fix: set query_region.ctg to the actual chrom
    let mut resp = result?;
    resp.query_region.ctg = resp.gene.chrom.clone();

    Ok(Json(resp))
}
