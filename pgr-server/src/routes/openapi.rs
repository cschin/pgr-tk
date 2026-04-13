use axum::Json;
use utoipa::OpenApi;

use crate::models::{
    BundleRequest, BundleResponse, DbInfoEntry, ErrorBody, GeneInfo, GeneQueryRequest,
    GeneQueryResponse, HitRecord, InfoResponse, QueryRegionInfo, RegionQueryRequest,
    RegionQueryResponse, SeqEntry,
};

use super::{bundle, health, query, sequences};

/// Top-level OpenAPI document for pgr-server.
#[derive(OpenApi)]
#[openapi(
    info(
        title = "pgr-server",
        version = "0.1.0",
        description = "PGR-TK pangenome HTTP query server. \
                        Loads a shimmer-based pangenome index once at startup \
                        and serves fast JSON queries over all configured haplotype databases."
    ),
    paths(
        health::health,
        health::info,
        sequences::list_sequences,
        query::query_region,
        query::query_gene,
        bundle::run_bundle,
    ),
    components(schemas(
        SeqEntry,
        RegionQueryRequest,
        RegionQueryResponse,
        QueryRegionInfo,
        HitRecord,
        GeneQueryRequest,
        GeneQueryResponse,
        GeneInfo,
        BundleRequest,
        BundleResponse,
        DbInfoEntry,
        InfoResponse,
        ErrorBody,
    )),
    tags(
        (name = "health",    description = "Liveness probe and server metadata"),
        (name = "sequences", description = "List sequences indexed in the database"),
        (name = "query",     description = "Shimmer-based pangenome alignment queries"),
        (name = "bundle",    description = "Principal bundle decomposition via pgr binary"),
    )
)]
pub struct ApiDoc;

/// Return the OpenAPI 3.0 specification as JSON.
#[utoipa::path(
    get,
    path = "/api/v1/openapi.json",
    responses(
        (status = 200, description = "OpenAPI 3.0 specification in JSON format")
    ),
    tag = "health"
)]
pub async fn openapi_json() -> Json<utoipa::openapi::OpenApi> {
    Json(ApiDoc::openapi())
}
