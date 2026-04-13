use serde::{Deserialize, Serialize};
use utoipa::ToSchema;

// ── /sequences ───────────────────────────────────────────────────────────────

#[derive(Debug, Serialize, ToSchema)]
pub struct SeqEntry {
    /// Internal sequence ID
    pub sid: u32,
    /// Source sample / file path
    pub src: String,
    /// Contig name (e.g. `"GRCh38#0#chr1"`)
    pub ctg: String,
    /// Sequence length in bp
    pub length: u32,
}

// ── /query/region ─────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize, ToSchema)]
pub struct RegionQueryRequest {
    /// Database name (defaults to first configured DB)
    pub db: Option<String>,
    /// Contig / chromosome name as stored in the index (e.g. `"chr6"`)
    pub chrom: String,
    /// Start coordinate (0-based, inclusive)
    pub start: u64,
    /// End coordinate (0-based, exclusive)
    pub end: u64,
    /// Substring identifying the reference sample (e.g. `"GRCh38"`)
    pub ref_sample_hint: Option<String>,
    /// Max shimmer-pair occurrence for alignment filtering
    pub max_count: Option<u32>,
    pub max_query_count: Option<u32>,
    pub max_target_count: Option<u32>,
    /// Merge hits within this many bp
    pub merge_range_tol: Option<usize>,
    /// Minimum alignment anchor count to report a hit
    pub min_anchor_count: Option<usize>,
    /// When `true` (default), return aligned sequences as inline FASTA
    #[serde(default = "default_true")]
    pub include_sequences: bool,
}

fn default_true() -> bool { true }

#[derive(Debug, Serialize, ToSchema)]
pub struct QueryRegionInfo {
    /// Reference sample hint used for the query
    pub src: String,
    /// Contig name queried
    pub ctg: String,
    /// Actual start used (clamped to contig length)
    pub bgn: u64,
    /// Actual end used (clamped to contig length)
    pub end: u64,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct HitRecord {
    /// Source sample / file path of the hit sequence
    pub src: String,
    /// Contig name of the hit sequence
    pub ctg: String,
    /// Hit start on the target contig (0-based)
    pub bgn: u32,
    /// Hit end on the target contig (0-based)
    pub end: u32,
    /// 0 = forward, 1 = reverse-complement
    pub orientation: u32,
    /// Number of shimmer alignment anchors supporting this hit
    pub anchor_count: usize,
    /// Display name composed from sample, contig, coordinates and orientation
    pub name: String,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct RegionQueryResponse {
    pub query_region: QueryRegionInfo,
    pub hits: Vec<HitRecord>,
    /// FASTA text of all hit sequences; `null` when `include_sequences` is `false`
    pub sequences_fasta: Option<String>,
}

// ── /query/gene ───────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize, ToSchema)]
pub struct GeneQueryRequest {
    /// Database name (defaults to first configured DB)
    pub db: Option<String>,
    /// Gene symbol to look up (e.g. `"HLA-A"`)
    pub gene_name: String,
    /// Flanking bp added on each side of the gene span (default 100 000)
    pub flank: Option<u64>,
    /// Substring identifying the reference sample (e.g. `"GRCh38"`)
    pub ref_sample_hint: Option<String>,
    pub max_count: Option<u32>,
    pub max_query_count: Option<u32>,
    pub max_target_count: Option<u32>,
    pub merge_range_tol: Option<usize>,
    pub min_anchor_count: Option<usize>,
    /// When `true` (default), return aligned sequences as inline FASTA
    #[serde(default = "default_true")]
    pub include_sequences: bool,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct GeneInfo {
    pub gene_name: String,
    /// Primary-assembly chromosome (e.g. `"chr6"`)
    pub chrom: String,
    /// Gene start (0-based, from annotation DB)
    pub start: i64,
    /// Gene end (0-based exclusive, from annotation DB)
    pub end: i64,
    /// Strand: `"+"` or `"-"`
    pub strand: String,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct GeneQueryResponse {
    pub gene: GeneInfo,
    pub query_region: QueryRegionInfo,
    pub hits: Vec<HitRecord>,
    pub sequences_fasta: Option<String>,
}

// ── /bundle ───────────────────────────────────────────────────────────────────

#[derive(Debug, Deserialize, ToSchema)]
pub struct BundleRequest {
    /// Multi-sequence FASTA text to decompose
    pub sequences_fasta: String,
    /// When `true`, also run `pgr bundle svg --html` and return HTML
    #[serde(default)]
    pub generate_html: bool,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct BundleResponse {
    /// BED-format output from `pgr bundle decomp`
    pub bed: String,
    /// HTML/SVG visualisation; `null` when `generate_html` is `false`
    pub html: Option<String>,
}

// ── /info ─────────────────────────────────────────────────────────────────────

#[derive(Debug, Serialize, ToSchema)]
pub struct DbInfoEntry {
    pub name: String,
    pub db_prefix: String,
    pub memory_mode: String,
    pub n_sequences: usize,
    pub gene_db: Option<String>,
    pub ref_sample_hint: String,
}

#[derive(Debug, Serialize, ToSchema)]
pub struct InfoResponse {
    pub default_db: String,
    pub databases: Vec<DbInfoEntry>,
}

// ── error body ────────────────────────────────────────────────────────────────

/// Standard error response body returned for 4xx / 5xx responses.
#[derive(Debug, Serialize, ToSchema)]
pub struct ErrorBody {
    pub error: String,
}
