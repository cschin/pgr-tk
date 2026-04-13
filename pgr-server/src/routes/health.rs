use axum::{extract::State, Json};
use serde_json::{json, Value};

use crate::models::{DbInfoEntry, InfoResponse};
use crate::state::AppState;

pub async fn health() -> Json<Value> {
    Json(json!({ "status": "ok" }))
}

pub async fn info(State(state): State<AppState>) -> Json<InfoResponse> {
    let mut databases: Vec<DbInfoEntry> = state
        .databases
        .iter()
        .map(|(_, db)| DbInfoEntry {
            name: db.name.clone(),
            db_prefix: db.db_prefix.clone(),
            memory_mode: db.memory_mode.clone(),
            n_sequences: db.n_sequences,
            gene_db: db.gene_db_path.clone(),
            ref_sample_hint: db.ref_sample_hint.clone(),
        })
        .collect();
    databases.sort_by(|a, b| a.name.cmp(&b.name));

    Json(InfoResponse {
        default_db: state.default_db.clone(),
        databases,
    })
}
