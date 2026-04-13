# pgr-server — Implementation Reference

> Status: **v0.1.0** — initial working implementation.
> For the original architectural rationale see `DESIGN.md`.
> This document describes what is **actually implemented** and is the
> authoritative reference for future contributors.

---

## 1. What It Does

`pgr-server` is a long-running Axum HTTP server that loads a PGR-TK
pangenome index once at startup and answers shimmer-based sequence queries
over a JSON/REST API.  The core motivation is to eliminate the multi-minute
startup cost of `pgr query seqs` (which rebuilds the shimmer HashMap every
invocation) for interactive or scripted multi-gene workflows.

---

## 2. File Layout

```
pgr-server/
├── Cargo.toml                  dependencies
├── pgr-server.yaml.example     annotated example config
├── DESIGN.md                   architectural rationale (written before implementation)
├── IMPLEMENTATION.md           this file
└── src/
    ├── main.rs                 CLI parsing, config loading, server startup; mounts SwaggerUi
    ├── config.rs               Config structs + serde_yaml loader
    ├── state.rs                AppState, PangenomeDb, SyncSeqIndexDB, DB loader
    ├── error.rs                AppError → axum IntoResponse
    ├── models.rs               serde request/response types
    └── routes/
        ├── mod.rs              router wiring (all endpoints under /api/v1); exports ApiDoc
        ├── health.rs           GET /health, GET /info
        ├── sequences.rs        GET /sequences
        ├── query.rs            POST /query/region, POST /query/gene
        ├── bundle.rs           POST /bundle
        └── openapi.rs          GET /openapi.json handler + ApiDoc (#[derive(OpenApi)])
```

Key Cargo dependencies:

| Crate | Version | Role |
|---|---|---|
| `pgr-db` | workspace | shimmer index, AGC sequence I/O |
| `axum` | 0.7 | HTTP framework |
| `tokio` | 1 (full) | async runtime |
| `tower-http` | 0.5 | CORS + tracing middleware |
| `serde` / `serde_json` / `serde_yaml` | 1 / 1 / 0.9 | (de)serialisation |
| `clap` | 4 (derive + env) | CLI + env-var parsing |
| `tracing` + `tracing-subscriber` | 0.1 / 0.3 | structured logging |
| `rusqlite` | 0.39 (bundled) | gene annotation SQLite queries |
| `rayon` | 1 | inner parallelism (via pgr-db) |
| `anyhow` | 1 | error plumbing in non-handler code |
| `tempfile` | 3 | temp FASTA for the /bundle subprocess |
| `utoipa` | 4 (axum_extras) | `#[derive(ToSchema)]`, `#[derive(OpenApi)]`, `#[utoipa::path]` annotations |
| `utoipa-swagger-ui` | 7 (axum) | Serves Swagger UI HTML + static assets |

---

## 3. Startup Sequence

```
main()
  ├─ Cli::parse()                    clap derives CLI from struct fields + env vars
  ├─ config::load_yaml(path)?        serde_yaml; absent fields → built-in defaults
  ├─ apply CLI overrides             host, port, cors, log_level, pgr_binary, db_prefix
  ├─ tracing_subscriber::init()      EnvFilter: RUST_LOG > config log_level
  ├─ spawn_blocking → AppState::build(cfg)
  │    └─ for each DbConfig:
  │         SeqIndexDB::load_from_agc_index(prefix)   builds frag_location_map HashMap
  │         [if memory_mode="high"]
  │           SeqIndexDB::load_index_to_memory()       copies mdbv into RAM
  ├─ Router::new().nest("/api/v1", routes::router(state))
  │    .layer(TraceLayer)
  │    .layer(CorsLayer { GET, POST, Any origin })
  └─ axum::serve(TcpListener::bind(host:port), app).await
```

If no databases are configured (neither via YAML nor `--db-prefix`) the
process exits with a clear error before binding any port.

---

## 4. Shared State

### `SyncSeqIndexDB`

```rust
pub struct SyncSeqIndexDB(pub SeqIndexDB);
unsafe impl Send for SyncSeqIndexDB {}
unsafe impl Sync for SyncSeqIndexDB {}
```

`SeqIndexDB` wraps raw mmap handles and an `AGCFile` that internally uses
`Mutex<rusqlite::Connection>`.  After loading it is **read-only**; the
`unsafe impl` is sound because:

- The shimmer HashMap / mmap reads are immutable (`&self` only).
- Sequence retrieval (`AGCFile.get_seq`) is serialised through the inner
  `Mutex<AgcFile>`, so concurrent fetches are safe.

### `PangenomeDb`

```rust
#[derive(Clone)]
pub struct PangenomeDb {
    pub name:             String,
    pub db_prefix:        String,
    pub memory_mode:      String,   // "moderate" | "high"
    pub ref_sample_hint:  String,   // e.g. "GRCh38"
    pub gene_db_path:     Option<String>,
    pub n_sequences:      usize,
    pub seq_db:           Arc<SyncSeqIndexDB>,
}
```

### `AppState`

```rust
#[derive(Clone)]
pub struct AppState {
    pub databases:      Arc<HashMap<String, PangenomeDb>>,
    pub default_db:     String,
    pub query_defaults: QueryDefaults,
    pub pgr_binary:     String,    // path to the `pgr` binary; used by /bundle
}
```

`AppState` is `Clone` (all fields are either `Arc`, `String`, or primitive)
and is injected into every handler via Axum's `State` extractor.

---

## 5. Configuration

### Priority order (highest wins)

```
CLI flags  >  environment variables  >  YAML file  >  built-in defaults
```

### YAML schema (`pgr-server.yaml.example`)

```yaml
server:
  host:        "127.0.0.1"   # env: PGR_SERVER_HOST
  port:        3000           # env: PGR_SERVER_PORT
  cors_origin: "*"            # env: PGR_SERVER_CORS_ORIGIN
  log_level:   "info"         # env: PGR_SERVER_LOG_LEVEL
  pgr_binary:  "pgr"          # path to pgr binary (for /bundle)

databases:
  - name:            "hprc_r2"
    db_prefix:       "/data/hprc_r2/hprc_r2"
    gene_db:         "/data/hprc_r2/hg38.ncbiRefSeq.db"   # optional
    memory_mode:     "moderate"    # "moderate" | "high"
    ref_sample_hint: "GRCh38"

query_defaults:
  flank:            100000
  max_count:        128
  max_query_count:  128
  max_target_count: 128
  merge_range_tol:  100000
  min_anchor_count: 10
```

The **first** entry in `databases` is the default; requests that omit `"db"`
are routed to it.

### CLI flags

| Flag | Env var | Default |
|---|---|---|
| `--config <PATH>` | `PGR_SERVER_CONFIG` | _(none)_ |
| `--host <HOST>` | `PGR_SERVER_HOST` | `127.0.0.1` |
| `--port <PORT>` | `PGR_SERVER_PORT` | `3000` |
| `--cors-origin <O>` | `PGR_SERVER_CORS_ORIGIN` | `*` |
| `--log-level <L>` | `PGR_SERVER_LOG_LEVEL` | `info` |
| `--pgr-binary <B>` | `PGR_BINARY` | `pgr` |
| `--db-prefix <P>` | _(none)_ | _(none)_ |
| `--db-name <N>` | _(none)_ | `default` |

`--db-prefix` / `--db-name` are a quick-start shorthand that append one
`DbConfig` entry with `memory_mode=moderate` and `ref_sample_hint=GRCh38`.
For anything beyond a single database, use the YAML config.

### Memory modes

| Mode | What is loaded | Lookup cost |
|---|---|---|
| `moderate` | `frag_location_map` HashMap + mmap `.mdbv` | O(1) HashMap + mmap page fault |
| `high` | Above **plus** full `.mdbv` copied into RAM | O(1) HashMap + RAM read (no page fault) |

`low` memory mode is **not** supported by `pgr-server` (it is too slow for
a server workload and does not benefit from staying resident).

---

## 6. API Reference

### Base path: `/api/v1`

All endpoints return `Content-Type: application/json`.
Errors are `{ "error": "<message>" }` with the appropriate HTTP status.

---

### `GET /docs`  *(Swagger UI)*

Interactive API documentation rendered by Swagger UI.
Opens a browser UI that lets you explore and call every endpoint.

```
GET /api/v1/docs
```

Served by `utoipa-swagger-ui`; assets are bundled inside the binary.

---

### `GET /openapi.json`

The raw OpenAPI 3.0 specification as JSON.
Use this to generate client SDKs or import into Postman / Insomnia.

```
GET /api/v1/openapi.json
200 OK  Content-Type: application/json
{ "openapi": "3.0.3", "info": { ... }, "paths": { ... }, ... }
```

Generated by `utoipa` from `routes/openapi.rs::ApiDoc`.  Also served by the
SwaggerUi mount point so the UI can load it without a separate route.

---

### `GET /health`

Liveness probe.

```
200 OK
{ "status": "ok" }
```

---

### `GET /info`

Returns metadata for all loaded databases.

```
200 OK
{
  "default_db": "hprc_r2",
  "databases": [
    {
      "name":            "hprc_r2",
      "db_prefix":       "/data/hprc_r2/hprc_r2",
      "memory_mode":     "moderate",
      "n_sequences":     932,
      "gene_db":         "/data/hprc_r2/hg38.ncbiRefSeq.db",
      "ref_sample_hint": "GRCh38"
    }
  ]
}
```

`n_sequences` is the count of rows in the `.midx` `seq_index` table.

---

### `GET /sequences?db=<name>`

List all sequences in the database, sorted by internal sequence ID (`sid`).

```
GET /api/v1/sequences?db=hprc_r2
200 OK
[
  { "sid": 0, "src": "GRCh38", "ctg": "GRCh38#0#chr1",  "length": 248956422 },
  { "sid": 1, "src": "GRCh38", "ctg": "GRCh38#0#chr2",  "length": 242193529 },
  ...
]
```

`db` defaults to the first configured database.

---

### `POST /query/region`

Query all haplotype sequences that align to a reference genomic region.

**Request body** (all fields optional except `chrom`, `start`, `end`):

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

**Processing pipeline** (executed inside `spawn_blocking`):

1. Locate the reference sequence SID: scan `seq_info` for entries where
   `src` contains `ref_sample_hint` and `ctg == chrom`.
2. Clamp `[start, end]` to `[0, ref_len]`.
3. Fetch the reference subsequence with `get_sub_seq_by_id(sid, start, end)`.
4. Call `query_fragment_to_hps_from_mmap_file` with `gap_penalty=0.025`,
   `max_aln_chain_span=8`.
5. Filter alignments with `anchor_count < min_anchor_count`.
6. Determine orientation per alignment by majority vote of forward vs
   reverse hit-pair orientations.
7. Compute target coordinate range (min bgn, max end across the chain).
8. Separately merge forward and reverse ranges within `merge_range_tol` bp.
9. Build `HitRecord` for each merged range.
10. If `include_sequences=true`, fetch each hit's subsequence from AGC and
    return as inline FASTA (orientation=1 hits are reverse-complemented).

**Response**:

```json
{
  "query_region": {
    "src": "GRCh38",
    "ctg": "chr6",
    "bgn": 29842531,
    "end": 30045870
  },
  "hits": [
    {
      "src":          "/path/to/HG002.mat.agcrs",
      "ctg":          "HG002.mat#1#h1tg000004l",
      "bgn":          12345678,
      "end":          12456789,
      "orientation":  0,
      "anchor_count": 42,
      "name":         "HG002.mat::HG002.mat#1#h1tg000004l_12345678_12456789_0"
    }
  ],
  "sequences_fasta": ">HG002.mat::...\nACGT...\n"
}
```

`sequences_fasta` is `null` when `include_sequences` is `false`.

**Error codes**:

| HTTP | Condition |
|---|---|
| 404 | `db` not found |
| 404 | reference `chrom` not found for the given `ref_sample_hint` |
| 500 | shimmer query failed (index not loaded) |

---

### `POST /query/gene`

Look up a gene by name in the gene annotation SQLite database, then perform
the same query as `/query/region` over the gene span ± flank.

Requires `gene_db` to be configured for the selected database.

**Request body** (all fields optional except `gene_name`):

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
  "include_sequences": false
}
```

**Gene lookup query** (executed against the gene annotation SQLite DB):

```sql
SELECT gene_name, chrom, start, end, strand
FROM genes
WHERE gene_name = ?
ORDER BY
  CASE WHEN chrom NOT GLOB '*_*' THEN 0 ELSE 1 END,
  end - start DESC
LIMIT 1
```

The ordering prefers primary assembly chromosomes (no `_` in name) over
alt/patch contigs, and the longest span among ties.

**Response** — identical to `/query/region` plus a `gene` field:

```json
{
  "gene": {
    "gene_name": "HLA-A",
    "chrom":     "chr6",
    "start":     29942531,
    "end":       29945870,
    "strand":    "+"
  },
  "query_region": { "src": "GRCh38", "ctg": "chr6", "bgn": 29842531, "end": 30045870 },
  "hits":           [ ... ],
  "sequences_fasta": null
}
```

**Error codes**:

| HTTP | Condition |
|---|---|
| 400 | `db` has no `gene_db` configured |
| 404 | `db` not found |
| 404 | gene name not found |
| 500 | SQLite open/query error |

---

### `POST /bundle`

Run principal bundle decomposition on caller-supplied FASTA sequences.
Internally writes the FASTA to a temp file and calls the `pgr` binary.

**Request body**:

```json
{
  "sequences_fasta": ">seq1\nACGT...\n>seq2\nACGT...",
  "generate_html":   false
}
```

**Processing** (inside `spawn_blocking`):

1. Write `sequences_fasta` to a `NamedTempFile`.
2. Execute: `pgr bundle decomp --fasta <tmp>` → capture stdout as `bed`.
3. If `generate_html=true`: execute `pgr bundle svg --html --fasta <tmp>`
   → capture stdout as `html`.

**Response**:

```json
{
  "bed":  "seq1\t0\t1234\t42:0:0:100\n...",
  "html": null
}
```

`html` is `null` when `generate_html` is `false`.

**Error codes**:

| HTTP | Condition |
|---|---|
| 503 | `pgr` binary not found / not executable |
| 500 | `pgr bundle decomp` or `pgr bundle svg` exited with non-zero status |

---

## 7. Concurrency Model

```
Tokio executor (async)
  │
  ├─ accept connection, parse JSON   — async, non-blocking
  │
  └─ spawn_blocking(|| { ... })      — offloads CPU work to blocking pool
       │
       ├─ find_ref_sid               — scans seq_info HashMap (read-only)
       ├─ get_sub_seq_by_id          — hits AGCFile Mutex → rusqlite read
       ├─ query_fragment_to_hps_from_mmap_file
       │    └─ raw_query_fragment_from_mmap_midx
       │         (no rayon; single-threaded per query; multiple
       │          queries run concurrently via the blocking pool)
       └─ hits_to_fasta              — repeated get_sub_seq_by_id calls
```

Concurrent requests each get their own `spawn_blocking` task and their own
clone of `Arc<SyncSeqIndexDB>`.  There is no per-request lock — the
`SyncSeqIndexDB` is read-only after startup.  The only serialisation point
is inside `AGCFile` where sequence retrieval is gated by a `Mutex`.

---

## 8. Error Handling

Every handler returns `Result<Json<T>, AppError>`.

```rust
pub enum AppError {
    NotFound(String),           // 404
    BadRequest(String),         // 400
    ServiceUnavailable(String), // 503
    Internal(String),           // 500
}
```

`AppError` implements `axum::response::IntoResponse`, serialising as
`{ "error": "<message>" }`.

`From<tokio::task::JoinError>` is also implemented so that a panicking
`spawn_blocking` task produces a 500 rather than an unhandled error.

---

## 9. Gene Annotation Database

The gene SQLite database is built by `examples/hprc_r2/fetch_refseq_gtf_db.sh`
(or an equivalent script).  The schema that `pgr-server` expects:

```sql
CREATE TABLE genes (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id   TEXT NOT NULL,      -- Ensembl/NCBI gene ID (e.g. "LOC401127")
    gene_name TEXT NOT NULL,      -- display name queried by /query/gene
    chrom     TEXT NOT NULL,      -- chromosome / contig name
    start     INTEGER NOT NULL,   -- 0-based inclusive
    end       INTEGER NOT NULL,   -- 0-based exclusive
    strand    TEXT NOT NULL       -- "+" or "-"
);
CREATE UNIQUE INDEX idx_genes_id_chrom ON genes(gene_id, chrom);
CREATE INDEX idx_genes_name ON genes(gene_name);
```

The unique index on `(gene_id, chrom)` keeps primary-assembly and
alt-contig rows separate.  The `/query/gene` lookup uses only `gene_name`
and returns the single best row (primary assembly, longest span).

---

## 10. Known Limitations and Future Work

| Item | Notes |
|---|---|
| `chrom` must be an exact contig name | The query currently does `ctg == chrom` with no alias resolution (e.g. `"chr6"` will not match `"GRCh38#0#chr6"` unless that is the exact stored name) |
| No gene-name fuzzy matching | 404 on unknown gene; no "did you mean?" suggestions |
| Sequence fetch serialised through AGCFile Mutex | Concurrent hit-sequence retrieval is single-threaded per DB; large `include_sequences` responses may be slow |
| No connection pooling for gene SQLite | Each `/query/gene` request opens and closes a `rusqlite::Connection`; acceptable for moderate QPS |
| Bundle uses subprocess | `pgr bundle decomp` is called as a child process; tight integration via library API is a future improvement |
| No streaming | Large FASTA responses are buffered in RAM before being sent |
| No authentication or rate limiting | Suitable for trusted internal networks only |
| No `/metrics` endpoint | Prometheus-style counters/histograms not yet implemented |
| `low` memory mode not supported | The binary-search path is too slow for a shared server; use `moderate` or `high` |
