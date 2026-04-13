use axum::{
    extract::{Query, State},
    Json,
};
use serde::Deserialize;
use utoipa::IntoParams;

use crate::error::AppError;
use crate::models::SeqEntry;
use crate::state::AppState;

#[derive(Debug, Deserialize, IntoParams)]
pub struct SeqsQuery {
    /// Database name (defaults to first configured DB)
    pub db: Option<String>,
}

/// List every sequence indexed in the database, sorted by internal sequence ID.
#[utoipa::path(
    get,
    path = "/api/v1/sequences",
    params(SeqsQuery),
    responses(
        (status = 200, description = "Sequence list sorted by sid", body = Vec<SeqEntry>),
        (status = 404, description = "Named database not found",    body = ErrorBody),
    ),
    tag = "sequences"
)]
pub async fn list_sequences(
    State(state): State<AppState>,
    Query(params): Query<SeqsQuery>,
) -> Result<Json<Vec<SeqEntry>>, AppError> {
    let db_name = params.db.as_deref().unwrap_or(&state.default_db).to_string();

    let db = state
        .databases
        .get(&db_name)
        .ok_or_else(|| AppError::NotFound(format!("database '{db_name}' not found")))?;

    let seq_info = db
        .seq_db
        .0
        .seq_info
        .as_ref()
        .ok_or_else(|| AppError::Internal("seq_info not loaded".into()))?;

    let mut entries: Vec<SeqEntry> = seq_info
        .iter()
        .map(|(&sid, (ctg, src, length))| SeqEntry {
            sid,
            src: src.as_deref().unwrap_or("").to_string(),
            ctg: ctg.clone(),
            length: *length,
        })
        .collect();
    entries.sort_by_key(|e| e.sid);

    Ok(Json(entries))
}
