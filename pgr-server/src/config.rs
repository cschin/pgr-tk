use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Config {
    #[serde(default)]
    pub server: ServerSettings,

    #[serde(default)]
    pub databases: Vec<DbConfig>,

    #[serde(default)]
    pub query_defaults: QueryDefaults,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            server: ServerSettings::default(),
            databases: vec![],
            query_defaults: QueryDefaults::default(),
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ServerSettings {
    #[serde(default = "default_host")]
    pub host: String,
    #[serde(default = "default_port")]
    pub port: u16,
    #[serde(default = "default_cors")]
    pub cors_origin: String,
    #[serde(default = "default_log_level")]
    pub log_level: String,
    #[serde(default = "default_pgr_binary")]
    pub pgr_binary: String,
}

impl Default for ServerSettings {
    fn default() -> Self {
        Self {
            host: default_host(),
            port: default_port(),
            cors_origin: default_cors(),
            log_level: default_log_level(),
            pgr_binary: default_pgr_binary(),
        }
    }
}

fn default_host() -> String { "127.0.0.1".into() }
fn default_port() -> u16 { 3000 }
fn default_cors() -> String { "*".into() }
fn default_log_level() -> String { "info".into() }
fn default_pgr_binary() -> String { "pgr".into() }

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct DbConfig {
    pub name: String,
    pub db_prefix: String,
    #[serde(default)]
    pub gene_db: Option<String>,
    #[serde(default = "default_memory_mode")]
    pub memory_mode: String,
    #[serde(default = "default_ref_hint")]
    pub ref_sample_hint: String,
}

fn default_memory_mode() -> String { "moderate".into() }
fn default_ref_hint() -> String { "GRCh38".into() }

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct QueryDefaults {
    #[serde(default = "default_flank")]
    pub flank: u64,
    #[serde(default = "default_128")]
    pub max_count: u32,
    #[serde(default = "default_128")]
    pub max_query_count: u32,
    #[serde(default = "default_128")]
    pub max_target_count: u32,
    #[serde(default = "default_merge_tol")]
    pub merge_range_tol: usize,
    #[serde(default = "default_min_anchor")]
    pub min_anchor_count: usize,
}

impl Default for QueryDefaults {
    fn default() -> Self {
        Self {
            flank: default_flank(),
            max_count: default_128(),
            max_query_count: default_128(),
            max_target_count: default_128(),
            merge_range_tol: default_merge_tol(),
            min_anchor_count: default_min_anchor(),
        }
    }
}

fn default_flank() -> u64 { 100_000 }
fn default_128() -> u32 { 128 }
fn default_merge_tol() -> usize { 100_000 }
fn default_min_anchor() -> usize { 10 }

pub fn load_yaml(path: &str) -> Result<Config> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot read config file: {path}"))?;
    let cfg: Config = serde_yaml::from_str(&text)
        .with_context(|| format!("cannot parse YAML config: {path}"))?;
    Ok(cfg)
}
