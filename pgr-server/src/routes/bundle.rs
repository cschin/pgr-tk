use axum::{extract::State, Json};
use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

use crate::error::AppError;
use crate::models::{BundleRequest, BundleResponse};
use crate::state::AppState;

/// Run principal bundle decomposition on caller-supplied FASTA sequences.
///
/// Writes the FASTA to a temporary file, invokes `pgr bundle decomp` to
/// produce a BED file, and optionally runs `pgr bundle svg --html` to produce
/// an interactive HTML visualisation.  All intermediate files are written to
/// a temporary directory that is deleted after the response is built.
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
        // Write FASTA to a temp file; keep it alive for the duration of the task.
        let mut fasta_tmp = NamedTempFile::new()
            .map_err(|e| AppError::Internal(format!("tempfile: {e}")))?;
        fasta_tmp
            .write_all(req.sequences_fasta.as_bytes())
            .map_err(|e| AppError::Internal(format!("write fasta tmp: {e}")))?;
        let fasta_path = fasta_tmp.path().to_path_buf();

        // Use a temp dir for output files so they are cleaned up automatically.
        let work_dir = tempfile::tempdir()
            .map_err(|e| AppError::Internal(format!("tempdir: {e}")))?;
        let prefix = work_dir.path().join("bundle");
        let prefix_str = prefix.to_string_lossy();

        // Run: pgr bundle decomp --fastx-path <fasta> --output-prefix <prefix>
        let decomp_out = Command::new(&pgr_binary)
            .args(["bundle", "decomp", "--fastx-path"])
            .arg(&fasta_path)
            .arg("--output-prefix")
            .arg(&*prefix_str)
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

        let bed_path = work_dir.path().join("bundle.bed");
        let bed = std::fs::read_to_string(&bed_path)
            .map_err(|e| AppError::Internal(format!("reading bundle.bed: {e}")))?;

        // Optionally run: pgr bundle svg --html --bed-file-path <bed> --output-prefix <prefix>
        let html = if req.generate_html {
            let svg_out = Command::new(&pgr_binary)
                .args(["bundle", "svg", "--html", "--bed-file-path"])
                .arg(&bed_path)
                .arg("--output-prefix")
                .arg(&*prefix_str)
                .output()
                .map_err(|e| {
                    AppError::ServiceUnavailable(format!("pgr bundle svg failed to start: {e}"))
                })?;

            if !svg_out.status.success() {
                let stderr = String::from_utf8_lossy(&svg_out.stderr);
                return Err(AppError::Internal(format!(
                    "pgr bundle svg exited with error: {stderr}"
                )));
            }

            let html_path = work_dir.path().join("bundle.html");
            let content = std::fs::read_to_string(&html_path)
                .map_err(|e| AppError::Internal(format!("reading bundle.html: {e}")))?;
            Some(content)
        } else {
            None
        };

        Ok::<_, AppError>(BundleResponse { bed, html })
    })
    .await
    .map_err(AppError::from)?;

    Ok(Json(result?))
}
