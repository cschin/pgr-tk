use axum::{
    extract::{Query, State},
    Json,
};
use serde::Deserialize;

use crate::error::AppError;
use crate::models::SeqEntry;
use crate::state::AppState;

#[derive(Debug, Deserialize)]
pub struct SeqsQuery {
    pub db: Option<String>,
}

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
