mod config;
mod error;
mod models;
mod routes;
mod state;

use anyhow::Context;
use axum::http::Method;
use clap::Parser;
use tower_http::{
    cors::{Any, CorsLayer},
    trace::TraceLayer,
};
use tracing::info;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};
use utoipa::OpenApi as _;
use utoipa_swagger_ui::SwaggerUi;

use config::Config;
use routes::ApiDoc;
use state::AppState;

#[derive(Parser, Debug)]
#[clap(name = "pgr-server", about = "PGR-TK pangenome HTTP query server")]
struct Cli {
    /// Path to YAML configuration file
    #[clap(long, short, env = "PGR_SERVER_CONFIG")]
    config: Option<String>,

    /// Bind address (overrides config)
    #[clap(long, env = "PGR_SERVER_HOST")]
    host: Option<String>,

    /// TCP port (overrides config)
    #[clap(long, env = "PGR_SERVER_PORT")]
    port: Option<u16>,

    /// CORS origin (overrides config)
    #[clap(long, env = "PGR_SERVER_CORS_ORIGIN")]
    cors_origin: Option<String>,

    /// Log level: trace|debug|info|warn|error (overrides config)
    #[clap(long, env = "PGR_SERVER_LOG_LEVEL")]
    log_level: Option<String>,

    /// Path to pgr binary used by /bundle (overrides config)
    #[clap(long, env = "PGR_BINARY")]
    pgr_binary: Option<String>,

    /// Database prefix (shorthand for a single-db setup)
    #[clap(long)]
    db_prefix: Option<String>,

    /// Database name to use with --db-prefix
    #[clap(long, default_value = "default")]
    db_name: String,
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Load base config from YAML (or use defaults)
    let mut cfg: Config = if let Some(path) = &cli.config {
        config::load_yaml(path).with_context(|| format!("loading config '{path}'"))?
    } else {
        Config::default()
    };

    // Apply CLI overrides
    if let Some(host) = cli.host {
        cfg.server.host = host;
    }
    if let Some(port) = cli.port {
        cfg.server.port = port;
    }
    if let Some(cors) = cli.cors_origin {
        cfg.server.cors_origin = cors;
    }
    if let Some(ll) = cli.log_level {
        cfg.server.log_level = ll;
    }
    if let Some(bin) = cli.pgr_binary {
        cfg.server.pgr_binary = bin;
    }
    if let Some(prefix) = cli.db_prefix {
        cfg.databases.push(config::DbConfig {
            name: cli.db_name,
            db_prefix: prefix,
            gene_db: None,
            memory_mode: "moderate".into(),
            ref_sample_hint: "GRCh38".into(),
        });
    }

    if cfg.databases.is_empty() {
        anyhow::bail!(
            "No databases configured. Supply --db-prefix on the CLI or add a \
             `databases:` section to your YAML config."
        );
    }

    // Initialise tracing before anything else logs
    let log_filter = EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| EnvFilter::new(&cfg.server.log_level));
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(log_filter)
        .init();

    // Save bind parameters before cfg is moved into the blocking task
    let bind_addr: std::net::SocketAddr = format!("{}:{}", cfg.server.host, cfg.server.port)
        .parse()
        .context("invalid bind address")?;

    info!(
        %bind_addr,
        n_databases = cfg.databases.len(),
        "Starting pgr-server"
    );

    // Load all databases (blocking work; done before accepting connections)
    let state = tokio::task::spawn_blocking(move || AppState::build(cfg))
        .await
        .context("database load task panicked")??;

    info!(default_db = %state.default_db, "All databases loaded");

    // CORS layer
    let cors = CorsLayer::new()
        .allow_methods([Method::GET, Method::POST])
        .allow_headers(Any)
        .allow_origin(Any);

    // Build the Axum router.
    // SwaggerUi is merged at the top level so it can serve its own static
    // assets; the actual API routes are nested under /api/v1.
    let app = axum::Router::new()
        .merge(
            SwaggerUi::new("/api/v1/docs")
                .url("/api/v1/openapi.json", ApiDoc::openapi()),
        )
        .nest("/api/v1", routes::router(state))
        .layer(TraceLayer::new_for_http())
        .layer(cors);

    let listener = tokio::net::TcpListener::bind(bind_addr)
        .await
        .with_context(|| format!("binding to {bind_addr}"))?;
    info!(%bind_addr, "Listening");
    axum::serve(listener, app)
        .await
        .context("axum server error")?;

    Ok(())
}
