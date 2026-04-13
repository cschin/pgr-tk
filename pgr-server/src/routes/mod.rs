pub mod bundle;
pub mod health;
pub mod query;
pub mod sequences;

use axum::{
    routing::{get, post},
    Router,
};

use crate::state::AppState;

pub fn router(state: AppState) -> Router {
    Router::new()
        .route("/health", get(health::health))
        .route("/info", get(health::info))
        .route("/sequences", get(sequences::list_sequences))
        .route("/query/region", post(query::query_region))
        .route("/query/gene", post(query::query_gene))
        .route("/bundle", post(bundle::run_bundle))
        .with_state(state)
}
