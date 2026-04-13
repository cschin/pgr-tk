# pgr-server — Design Document

## 1. Motivation

`pgr query seqs` with `--memory-mode moderate` (HashMap in RAM) gives fast
per-query latency but pays a large startup cost every invocation: the
`frag_location_map` HashMap is rebuilt from the `.mdbi` file each time the
binary is executed. For a 466-haplotype HPRC r2 index this load can take
minutes, making interactive or scripted multi-gene workflows slow.

`pgr-server` solves this by loading the index **once** at startup and
exposing it over HTTP. Subsequent queries incur only the actual
shimmer-lookup + alignment cost, not the index-load overhead. The server
also eliminates the per-query filesystem round-trip for gene coordinate
lookup and reference sequence extraction.

---

## 2. Architecture Overview

```
                        ┌──────────────────────────────────────┐
                        │              pgr-server               │
                        │                                       │
  HTTP client ──────────►  Axum router                         │
  (browser / script)   │      │                                │
                        │      ├─► /health, /info              │
                        │      ├─► /sequences          ┌───────┤
                        │      ├─► /query/region  ─────► Shared│
                        │      ├─► /query/gene    ─────► AppState
                        │      └─► /bundle        ─────►       │
                        │                          │   ├ SeqIndexDB
                        │  spawn_blocking ─────────┘   │  (moderate/high)
                        │  (CPU-heavy ops)              ├ gene SQLite conn
                        │                               └ config
                        └──────────────────────────────────────┘
                                        │
                        ┌───────────────┼──────────────────────┐
                        │    pgr-db     │                       │
                        │  .agcrs  .mdbi  .mdbv  .midx         │
                        └──────────────────────────────────────┘
```

### Key design decisions

| Decision | Rationale |
|---|---|
| Axum + Tokio | Mature async ecosystem; Axum's `State` extractor maps cleanly to the shared DB handle |
| `Arc<RwLock<SeqIndexDB>>` | Index loaded once; queries are read-only → many concurrent readers, zero writers after init |
| `spawn_blocking` for shimmer ops | `query_fragment_to_hps_from_mmap_file` and bundle decomp are CPU-bound; must not block the Tokio executor |
| Stateless request model | Each request is self-contained; no server-side session or job queue needed for the initial version |
| JSON over REST | Simple; easy to consume from Python scripts, shell (`curl`), or a future browser UI |

---

## 3. Crate Structure

```
pgr-server/
├── Cargo.toml
├── pgr-server.yaml.example   # annotated example config
└── src/
    ├── main.rs          # CLI args (clap), config loading, Tokio runtime, server startup
    ├── config.rs        # ServerConfig / DbConfig structs; serde_yaml deserialisation; CLI merge
    ├── state.rs         # AppState, PangenomeDb; DB loading logic
    ├── error.rs         # AppError → axum IntoResponse
    ├── models.rs        # Serde request / response types
    └── routes/
        ├── mod.rs       # Router assembly
        ├── health.rs    # GET /health, GET /info
        ├── sequences.rs # GET /sequences
        ├── query.rs     # POST /query/region, POST /query/gene
        └── bundle.rs    # POST /bundle
```

Key dependencies (`Cargo.toml`):

```toml
[dependencies]
pgr-db      = { path = "../pgr-db" }
axum        = { version = "0.7", features = ["macros"] }
tokio       = { version = "1",   features = ["full"] }
tower-http  = { version = "0.5", features = ["cors", "trace"] }
serde       = { version = "1",   features = ["derive"] }
serde_json  = "1"
serde_yaml  = "0.9"      # YAML config deserialisation
clap        = { version = "4",   features = ["derive", "env"] }
tracing     = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
rusqlite    = { version = "0.39", features = ["bundled"] }
rayon       = "1"
```

---

## 4. Shared State

```rust
/// One loaded pangenome database.
pub struct PangenomeDb {
    /// Pangenome index loaded in moderate or high memory mode.
    /// RwLock: many concurrent read-queries, no writes after init.
    pub seq_db: Arc<RwLock<SeqIndexDB>>,

    /// Path to the SQLite gene annotation DB.
    /// Each request opens a short-lived read-only rusqlite connection.
    pub gene_db_path: Option<PathBuf>,

    /// Logical name from config (e.g. "hprc_r2").
    pub name: String,

    /// Filesystem prefix (informational).
    pub db_prefix: String,

    /// "moderate" or "high" (informational).
    pub memory_mode: String,

    /// Default reference sample hint (e.g. "GRCh38").
    pub ref_sample_hint: String,
}

pub struct AppState {
    /// All configured databases, keyed by logical name.
    /// Populated at startup; never mutated afterward.
    pub databases: HashMap<String, PangenomeDb>,

    /// Name of the default database (first entry in config).
    pub default_db: String,

    /// Server-wide query defaults from config.
    pub query_defaults: QueryDefaults,
}
```

`SeqIndexDB` internally holds:
- `frag_location_map: FxHashMap<(u64,u64), (usize,usize)>` — shimmer → `.mdbv` offset (moderate mode)
- `frag_map_file: Mmap` — mmap'd `.mdbv`
- `agc_db: AGCSeqDB` — sequence fetcher

All fields are read-only after construction, so `RwLock` read-locks
will never contend in practice.

---

## 5. Configuration

Configuration is layered in priority order (highest wins):

```
CLI flags  >  environment variables  >  YAML config file  >  built-in defaults
```

### 5.1 YAML config file

Specified via `--config <path>` (or `PGR_SERVER_CONFIG` env var).
All fields are optional; unset fields fall back to defaults.

```yaml
# pgr-server.yaml

server:
  host: "127.0.0.1"    # bind address
  port: 3000           # TCP port
  cors_origin: "*"     # allowed CORS origin (* = any)
  log_level: "info"    # tracing log level: trace | debug | info | warn | error

databases:
    # At least one entry required.
    # The first entry is the default when the client does not specify a db name.
  - name: "hprc_r2"                        # logical name used in API requests
    db_prefix: "/data/hprc_r2/hprc_r2"     # path prefix for .agcrs/.mdbi/.mdbv/.midx
    gene_db: "/data/hprc_r2/hg38.ncbiRefSeq.db"  # optional; omit to disable /query/gene
    memory_mode: "moderate"                # moderate | high
    ref_sample_hint: "GRCh38"             # substring to identify the reference sample

  - name: "hprc_r1"
    db_prefix: "/data/hprc_r1/hprc_r1"
    gene_db: "/data/hprc_r1/hg38.ncbiRefSeq.db"
    memory_mode: "high"
    ref_sample_hint: "GRCh38"

query_defaults:
  flank:            100000
  max_count:        128
  max_query_count:  128
  max_target_count: 128
  merge_range_tol:  100000
  min_anchor_count: 10
```

### 5.2 CLI flags

CLI flags override the YAML config for single-database quick-start use.

```
pgr-server [OPTIONS]

Options:
  -c, --config <PATH>            YAML config file  [env: PGR_SERVER_CONFIG]
  -p, --db-prefix <PREFIX>       DB prefix (overrides config databases[0].db_prefix)
  -g, --gene-db <PATH>           Gene annotation DB (overrides config databases[0].gene_db)
  -m, --memory-mode <MODE>       moderate | high  (overrides config databases[0].memory_mode)
  -H, --host <HOST>              Bind address  [default: 127.0.0.1]
  -P, --port <PORT>              TCP port  [default: 3000]
      --ref-sample-hint <HINT>   Reference sample hint  [default: GRCh38]
      --cors-origin <ORIGIN>     Allowed CORS origin  [default: *]
      --log-level <LEVEL>        Log level  [default: info]
```

### 5.3 Environment variables

Each YAML key maps to an env var with the `PGR_SERVER_` prefix:

| Env var | YAML equivalent |
|---|---|
| `PGR_SERVER_CONFIG` | — (path to config file) |
| `PGR_SERVER_HOST` | `server.host` |
| `PGR_SERVER_PORT` | `server.port` |
| `PGR_SERVER_LOG_LEVEL` | `server.log_level` |
| `PGR_SERVER_CORS_ORIGIN` | `server.cors_origin` |

Individual database entries cannot be overridden via env vars; use the
CLI flags for single-database quick-start or the YAML file for
multi-database deployments.

### 5.4 Startup behaviour

1. Parse CLI flags (including `--config` path).
2. If a config file is present, deserialize it with `serde_yaml`.
3. Merge: CLI flags overwrite the corresponding config-file fields.
4. Apply env-var overrides on top.
5. For each database entry: load `SeqIndexDB` in the specified memory mode.
   Log progress and memory mode for each.
6. Start the Axum server.

---

## 6. API Specification

### Base path: `/api/v1`

All responses are `Content-Type: application/json` unless noted.
Errors return `{ "error": "<message>" }` with an appropriate HTTP status.

When multiple databases are configured, requests that touch the index
accept an optional `"db"` field naming the target database.
If omitted, the first entry in `databases` is used.

---

### 6.1  `GET /health`

Liveness probe.

**Response 200**
```json
{ "status": "ok" }
```

---

### 6.2  `GET /info`

Server and database metadata.

**Response 200**
```json
{
  "default_db": "hprc_r2",
  "databases": [
    {
      "name":            "hprc_r2",
      "db_prefix":       "/data/hprc_r2/hprc_r2",
      "memory_mode":     "moderate",
      "n_sequences":     466,
      "gene_db":         "hg38.ncbiRefSeq.db",
      "ref_sample_hint": "GRCh38"
    },
    {
      "name":            "hprc_r1",
      "db_prefix":       "/data/hprc_r1/hprc_r1",
      "memory_mode":     "high",
      "n_sequences":     94,
      "gene_db":         null,
      "ref_sample_hint": "GRCh38"
    }
  ]
}
```

---

### 6.3  `GET /sequences`

List all sequences indexed in the database.
Optional query parameters filter the result.

**Query parameters**

| Parameter | Type   | Description |
|-----------|--------|-------------|
| `hint`    | string | Keep rows whose `src` contains this substring (e.g. `GRCh38`) |
| `chrom`   | string | Keep rows whose `ctg` ends with this string (e.g. `chr6`) |

**Query parameters**

| Parameter | Type   | Description |
|-----------|--------|-------------|
| `db`      | string | Database name from config (default: first db) |
| `hint`    | string | Keep rows whose `src` contains this substring (e.g. `GRCh38`) |
| `chrom`   | string | Keep rows whose `ctg` ends with this string (e.g. `chr6`) |

**Response 200**
```json
[
  { "sid": 0, "src": "GRCh38", "ctg": "GRCh38#0#chr1",  "length": 248956422 },
  { "sid": 1, "src": "GRCh38", "ctg": "GRCh38#0#chr2",  "length": 242193529 },
  ...
]
```

---

### 6.4  `POST /query/region`

Query all haplotype sequences covering an explicit genomic region.

**Request body**
```json
{
  "db":               "hprc_r2",
  "chrom":            "chr6",
  "start":            29842531,
  "end":              30045870,
  "ref_sample_hint":  "GRCh38",
  "max_count":        128,
  "max_query_count":  128,
  "max_target_count": 128,
  "merge_range_tol":  100000,
  "min_anchor_count": 10,
  "include_sequences": true
}
```

All fields except `chrom`, `start`, `end` are optional (server defaults apply).

**Response 200**
```json
{
  "query_region": {
    "src": "GRCh38",
    "ctg": "GRCh38#0#chr6",
    "bgn": 29842531,
    "end": 30045870
  },
  "hits": [
    {
      "src":          "HG002.mat",
      "ctg":          "HG002.mat#1#h1tg000004l",
      "bgn":          12345678,
      "end":          12456789,
      "orientation":  0,
      "anchor_count": 42,
      "name":         "HG002.mat::HG002.mat#1#h1tg000004l_12345678_12456789_0"
    },
    ...
  ],
  "sequences_fasta": ">HG002.mat::...\nACGT...\n..."  // null if include_sequences=false
}
```

---

### 6.5  `POST /query/gene`

Look up a gene by name in the gene annotation DB, then run `/query/region`
with the gene span ± flank. Requires `--gene-db` at server startup.

**Request body**
```json
{
  "db":               "hprc_r2",
  "gene_name":        "HLA-A",
  "flank":            100000,
  "ref_sample_hint":  "GRCh38",
  "max_count":        128,
  "max_query_count":  128,
  "max_target_count": 128,
  "merge_range_tol":  100000,
  "min_anchor_count": 10,
  "include_sequences": true
}
```

**Response 200** — same schema as `/query/region`, with an added `gene` field:
```json
{
  "gene": {
    "gene_name": "HLA-A",
    "chrom":     "chr6",
    "start":     29942531,
    "end":       29945870,
    "strand":    "+"
  },
  "query_region": { ... },
  "hits":           [ ... ],
  "sequences_fasta": "..."
}
```

**Response 404** — gene not found
```json
{ "error": "gene 'HLA-X' not found in gene db", "similar": ["HLA-A", "HLA-B"] }
```

**Response 503** — gene db not loaded
```json
{ "error": "gene db not configured; restart server with --gene-db" }
```

---

### 6.6  `POST /bundle`

Run principal bundle decomposition (and optionally SVG/HTML generation) on
a set of sequences provided as FASTA text. Intended to be called with
`sequences_fasta` from a `/query/*` response.

**Request body**
```json
{
  "sequences_fasta": ">seq1\nACGT...\n>seq2\nACGT...",
  "generate_html":   true
}
```

**Response 200**
```json
{
  "bed":  "seq1\t0\t1234\t42:0:0:100\n...",
  "html": "<!DOCTYPE html>..."   // null if generate_html=false
}
```

---

## 7. Concurrency Model

```
Tokio async thread pool
  │
  ├─ I/O tasks (connection accept, JSON (de)serialise)  → Tokio threads
  │
  └─ CPU tasks (shimmer lookup, bundle decomp)
       └─ spawn_blocking() → blocking thread pool
            └─ rayon par_iter() inside pgr-db
                 (rayon manages its own thread pool, sized to available CPUs)
```

`spawn_blocking` prevents the shimmer HashMap lookup (which calls rayon
internally) from starving the Tokio executor. The `RwLock` read-lock is
acquired inside the blocking closure so it does not cross an `.await`.

---

## 8. Error Handling

All handler functions return `Result<Json<T>, AppError>`. `AppError`
implements `IntoResponse`:

| Variant | HTTP status |
|---|---|
| `NotFound(msg)` | 404 |
| `BadRequest(msg)` | 400 |
| `ServiceUnavailable(msg)` | 503 |
| `Internal(msg)` | 500 |

---

## 9. Logging and Observability

- `tracing` + `tracing-subscriber` for structured logging
- Each request logs: method, path, gene/region, hit count, wall time
- `/metrics` endpoint (future): Prometheus-compatible counters for
  request count, error count, query latency histogram

---

## 10. Future Work

| Feature | Notes |
|---|---|
| Multi-DB support | Load several DB prefixes; client selects by name |
| Async gene-db pool | `deadpool-sqlite` to avoid per-request open/close |
| Streaming FASTA response | `StreamBody` for large hit sets |
| `/metrics` Prometheus endpoint | `axum-prometheus` crate |
| Auth / rate limiting | `tower` middleware layers |
| Bundle result caching | Cache decomp results keyed by query region hash |
| WebAssembly frontend | Replace `pgr-web` WASM viewer with a React+Fetch UI calling this API |
