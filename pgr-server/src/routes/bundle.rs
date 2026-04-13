use axum::{extract::State, Json};
use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

use crate::error::AppError;
use crate::models::{BundleRequest, BundleResponse};
use crate::state::AppState;

/// Run principal bundle decomposition on caller-supplied FASTA sequences.
///
/// Writes the FASTA to a temporary file and invokes `pgr bundle decomp`.
/// If `generate_html` is `true`, also runs `pgr bundle svg --html`.
#[utoipa::path(
    post,
    path = "/api/v1/bundle",
    request_body = BundleRequest,
    responses(
        (status = 200, description = "BED decomposition (and optional HTML)",
         body = BundleResponse),
        (status = 503, description = "pgr binary not found or not executable",
         body = ErrorBody),
        (status = 500, description = "pgr bundle command exited with an error",
         body = ErrorBody),
    ),
    tag = "bundle"
)]
pub async fn run_bundle(
    State(state): State<AppState>,
    Json(req): Json<BundleRequest>,
) -> Result<Json<BundleResponse>, AppError> {
    let pgr_binary = state.pgr_binary.clone();

    let result = tokio::task::spawn_blocking(move || {
        // Write FASTA to a temp file
        let mut fasta_tmp = NamedTempFile::new()
            .map_err(|e| AppError::Internal(format!("tempfile: {e}")))?;
        fasta_tmp
            .write_all(req.sequences_fasta.as_bytes())
            .map_err(|e| AppError::Internal(format!("write fasta tmp: {e}")))?;
        let fasta_path = fasta_tmp.path().to_path_buf();

        // Run: pgr bundle decomp --fasta <path>
        let decomp_out = Command::new(&pgr_binary)
            .args(["bundle", "decomp", "--fasta"])
            .arg(&fasta_path)
            .output()
            .map_err(|e| {
                AppError::ServiceUnavailable(format!("pgr bundle decomp failed to start: {e}"))
            })?;

        if !decomp_out.status.success() {
            let stderr = String::from_utf8_lossy(&decomp_out.stderr);
            return Err(AppError::Internal(format!(
                "pgr bundle decomp exited with error: {stderr}"
            )));
        }

        let bed = String::from_utf8_lossy(&decomp_out.stdout).into_owned();

        // Optionally run: pgr bundle svg --html <path>
        let html = if req.generate_html {
            let svg_out = Command::new(&pgr_binary)
                .args(["bundle", "svg", "--html", "--fasta"])
                .arg(&fasta_path)
                .output()
                .map_err(|e| {
                    AppError::ServiceUnavailable(format!(
                        "pgr bundle svg failed to start: {e}"
                    ))
                })?;

            if !svg_out.status.success() {
                let stderr = String::from_utf8_lossy(&svg_out.stderr);
                return Err(AppError::Internal(format!(
                    "pgr bundle svg exited with error: {stderr}"
                )));
            }

            Some(String::from_utf8_lossy(&svg_out.stdout).into_owned())
        } else {
            None
        };

        Ok::<_, AppError>(BundleResponse { bed, html })
    })
    .await
    .map_err(AppError::from)?;

    Ok(Json(result?))
}
