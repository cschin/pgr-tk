use std::collections::HashMap;
use std::sync::Arc;

use pgr_db::ext::SeqIndexDB;
use tracing::info;

use crate::config::{Config, DbConfig, QueryDefaults};

/// Newtype that makes `SeqIndexDB` usable as shared Axum state.
///
/// Safety: After initial loading the database is used read-only.  The only
/// interior mutability is `AGCFile.inner: Mutex<AgcFile>`, which already
/// provides thread-safe synchronisation for sequence retrieval.
pub struct SyncSeqIndexDB(pub SeqIndexDB);
unsafe impl Send for SyncSeqIndexDB {}
unsafe impl Sync for SyncSeqIndexDB {}

#[derive(Clone)]
pub struct PangenomeDb {
    pub name: String,
    pub db_prefix: String,
    pub memory_mode: String,
    pub ref_sample_hint: String,
    pub gene_db_path: Option<String>,
    pub n_sequences: usize,
    /// The loaded index, shared across all handler tasks.
    pub seq_db: Arc<SyncSeqIndexDB>,
}

/// Shared state injected into every Axum handler.
#[derive(Clone)]
pub struct AppState {
    pub databases: Arc<HashMap<String, PangenomeDb>>,
    pub default_db: String,
    pub query_defaults: QueryDefaults,
    pub pgr_binary: String,
}

impl AppState {
    /// Load all databases from `cfg`, blocking the current thread for each load.
    ///
    /// Call this once at startup before spawning the Axum server.
    pub fn build(cfg: Config) -> anyhow::Result<Self> {
        let default_db = cfg
            .databases
            .first()
            .map(|d| d.name.clone())
            .unwrap_or_default();

        let mut databases = HashMap::new();
        for db_cfg in cfg.databases {
            info!(name = %db_cfg.name, prefix = %db_cfg.db_prefix, mode = %db_cfg.memory_mode, "Loading database");
            let db = load_db(db_cfg)?;
            info!(name = %db.name, n_sequences = db.n_sequences, "Database loaded");
            databases.insert(db.name.clone(), db);
        }

        Ok(AppState {
            databases: Arc::new(databases),
            default_db,
            query_defaults: cfg.query_defaults,
            pgr_binary: cfg.server.pgr_binary,
        })
    }
}

fn load_db(db_cfg: DbConfig) -> anyhow::Result<PangenomeDb> {
    let mut seq_db = SeqIndexDB::new();

    match db_cfg.memory_mode.as_str() {
        "high" => {
            seq_db
                .load_from_agc_index(db_cfg.db_prefix.clone())
                .map_err(|e| anyhow::anyhow!("load_from_agc_index '{}': {e}", db_cfg.db_prefix))?;
            seq_db
                .load_index_to_memory()
                .map_err(|e| anyhow::anyhow!("load_index_to_memory '{}': {e}", db_cfg.db_prefix))?;
        }
        _ => {
            // "moderate" (default) — HashMap in RAM + mmap'd .mdbv
            seq_db
                .load_from_agc_index(db_cfg.db_prefix.clone())
                .map_err(|e| anyhow::anyhow!("load_from_agc_index '{}': {e}", db_cfg.db_prefix))?;
        }
    }

    let n_sequences = seq_db
        .seq_info
        .as_ref()
        .map(|m| m.len())
        .unwrap_or(0);

    Ok(PangenomeDb {
        name: db_cfg.name,
        db_prefix: db_cfg.db_prefix,
        memory_mode: db_cfg.memory_mode,
        ref_sample_hint: db_cfg.ref_sample_hint,
        gene_db_path: db_cfg.gene_db,
        n_sequences,
        seq_db: Arc::new(SyncSeqIndexDB(seq_db)),
    })
}
