use serde::{Deserialize, Serialize};

// ── /sequences ───────────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct SeqEntry {
    pub sid: u32,
    pub src: String,
    pub ctg: String,
    pub length: u32,
}

// ── /query/region ─────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub struct RegionQueryRequest {
    pub db: Option<String>,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub ref_sample_hint: Option<String>,
    pub max_count: Option<u32>,
    pub max_query_count: Option<u32>,
    pub max_target_count: Option<u32>,
    pub merge_range_tol: Option<usize>,
    pub min_anchor_count: Option<usize>,
    #[serde(default = "default_true")]
    pub include_sequences: bool,
}

fn default_true() -> bool { true }

#[derive(Debug, Serialize)]
pub struct QueryRegionInfo {
    pub src: String,
    pub ctg: String,
    pub bgn: u64,
    pub end: u64,
}

#[derive(Debug, Serialize)]
pub struct HitRecord {
    pub src: String,
    pub ctg: String,
    pub bgn: u32,
    pub end: u32,
    pub orientation: u32,
    pub anchor_count: usize,
    pub name: String,
}

#[derive(Debug, Serialize)]
pub struct RegionQueryResponse {
    pub query_region: QueryRegionInfo,
    pub hits: Vec<HitRecord>,
    pub sequences_fasta: Option<String>,
}

// ── /query/gene ───────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub struct GeneQueryRequest {
    pub db: Option<String>,
    pub gene_name: String,
    pub flank: Option<u64>,
    pub ref_sample_hint: Option<String>,
    pub max_count: Option<u32>,
    pub max_query_count: Option<u32>,
    pub max_target_count: Option<u32>,
    pub merge_range_tol: Option<usize>,
    pub min_anchor_count: Option<usize>,
    #[serde(default = "default_true")]
    pub include_sequences: bool,
}

#[derive(Debug, Serialize)]
pub struct GeneInfo {
    pub gene_name: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
}

#[derive(Debug, Serialize)]
pub struct GeneQueryResponse {
    pub gene: GeneInfo,
    pub query_region: QueryRegionInfo,
    pub hits: Vec<HitRecord>,
    pub sequences_fasta: Option<String>,
}

// ── /bundle ───────────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize)]
pub struct BundleRequest {
    pub sequences_fasta: String,
    #[serde(default)]
    pub generate_html: bool,
}

#[derive(Debug, Serialize)]
pub struct BundleResponse {
    pub bed: String,
    pub html: Option<String>,
}

// ── /info ─────────────────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct DbInfoEntry {
    pub name: String,
    pub db_prefix: String,
    pub memory_mode: String,
    pub n_sequences: usize,
    pub gene_db: Option<String>,
    pub ref_sample_hint: String,
}

#[derive(Debug, Serialize)]
pub struct InfoResponse {
    pub default_db: String,
    pub databases: Vec<DbInfoEntry>,
}
