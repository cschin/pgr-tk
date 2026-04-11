# pgr-tk Codebase Improvement Plan

This document surveys the current state of the codebase and proposes concrete improvements
in three areas: error handling, architecture, and code quality. It is intended as a
prioritised roadmap, not a requirement to fix everything at once.

---

## 1. Error Handling — Replacing `unwrap()` and `panic!`

### Current state

The codebase contains approximately **857** panic-prone calls:

| Crate | `unwrap()` | `expect()` | Total |
|-------|:----------:|:----------:|------:|
| pgr-db | 241 | 30 | 271 |
| pgr-bin | 285 | 260 | 545 |
| pgr-tk | 40 | 1 | 41 |
| **Total** | **566** | **291** | **857** |

Most of these are in production code paths, not tests. The three highest-risk hotspots are:

1. **`pgr-db/src/ext.rs`** — cascading `.as_ref().unwrap()` chains in every `Backend` match
   arm of `get_seq_by_id()` and `get_sub_seq_by_id()` (lines 363–499). Three of these even
   carry the comment `// TODO: handle Option unwrap properly`.
2. **`pgr-db/src/gff_db.rs`** — `.unwrap()` on every parsed field; a malformed GFF line
   causes an immediate panic with no diagnostic.
3. **`pgr-db/src/agc_io.rs`** — `CString::new(...).unwrap()` at the C FFI boundary; a
   null byte in a sequence name crashes the process.
4. **`pgr-tk/src/lib.rs`** — `.unwrap()` at the PyO3 boundary; a Rust panic tears down
   the Python interpreter with a raw backtrace rather than a Python exception.

### What to do

#### 1a. Define a crate-level error type in `pgr-db`

Replace ad-hoc panics with a typed error enum that all functions can return.

```rust
// pgr-db/src/error.rs  (new file)
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PgrError {
    #[error("sequence id {0} not found")]
    SeqIdNotFound(u32),

    #[error("contig name '{0}' not found in backend")]
    ContigNotFound(String),

    #[error("backend not initialised")]
    BackendNotInitialised,

    #[error("FFI error: {0}")]
    FfiError(String),

    #[error("parse error in {file} line {line}: {msg}")]
    ParseError { file: String, line: usize, msg: String },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("alignment failed: {0}")]
    AlignmentFailed(String),
}

pub type Result<T> = std::result::Result<T, PgrError>;
```

Then change function signatures from `fn get_seq_by_id(...) -> Vec<u8>` to
`fn get_seq_by_id(...) -> Result<Vec<u8>>` and propagate with `?`.

#### 1b. Harden the backend match arms in `ext.rs`

The current pattern:

```rust
let (ctg_name, sample_name, _) = self.seq_info.as_ref().unwrap().get(&sid).unwrap();
```

Should become:

```rust
let seq_info = self.seq_info.as_ref().ok_or(PgrError::BackendNotInitialised)?;
let (ctg_name, sample_name, _) = seq_info.get(&sid)
    .ok_or(PgrError::SeqIdNotFound(sid))?;
```

#### 1c. Harden the GFF/annotation parsers

`gff_db.rs` and `pgr variant annotate-vcf.rs` should skip or log malformed lines rather than
panic. The existing `TODO: need a proper parser` comment (annotate-vcf-file.rs lines 89–92)
should be resolved.

#### 1d. Fix the AGC FFI boundary

`CString::new(name)` returns `Result` because of potential null bytes. Replace:

```rust
CString::new(name).unwrap().into_raw()
```

with:

```rust
CString::new(name).map_err(|_| PgrError::FfiError(format!("null byte in name: {name}")))?
    .into_raw()
```

#### 1e. Convert PyO3 boundary panics to Python exceptions

At the `pgr-tk` boundary, map `PgrError` to a Python `RuntimeError` using PyO3's
`#[pyclass]` error support or `PyErr::new::<pyo3::exceptions::PyRuntimeError, _>`.

```rust
// Instead of:
self.0.get_seq_by_id(sid).unwrap()

// Use:
self.0.get_seq_by_id(sid)
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e.to_string()))?
```

#### Priority order

| Priority | Target | Effort |
|----------|--------|--------|
| 1 | Define `PgrError` + `Result` in pgr-db | 1–2 days |
| 2 | Fix `ext.rs` backend match arms | 1 day |
| 3 | PyO3 boundary in `pgr-tk/src/lib.rs` | 1 day |
| 4 | GFF / annotation parsers | 0.5 day |
| 5 | AGC FFI boundary | 0.5 day |
| 6 | Remaining `unwrap()` sweep in pgr-bin | 2–3 days |

---

## 2. Architecture Improvements

### 2a. Unified output format (Parquet / Arrow IPC)

The `.alnmap` text format (see `docs_src/alnmap_formap.md`) has variable column counts per
record type, no random access, and no compression. At whole-genome scale this is a
significant bottleneck.

**Proposal:** Add an optional binary output mode to `pgr align alnmap` using Apache Arrow IPC
(via the `arrow2` or `arrow-rs` crate). Three separate Arrow record-batch streams would
replace the interleaved text:

- `chains` — one batch per contig mapping (B/E pairs)
- `blocks` — M/V/S records with a `record_type` discriminant column
- `variants` — V records with positional columns + `ref_seq`/`alt_seq` as `LargeBinary`

Each stream would be accompanied by a lightweight positional index (chromosome → byte offset),
enabling random-access region queries without scanning the full file.

The text `.alnmap` format should be retained as a fallback for human inspection.

### 2b. Caching layer for fragment reconstruction

`seq_db.get_seq()` reconstructs sequences from compressed fragments on every call. For
workflows that repeatedly query the same contigs (e.g. SV analysis iterating over candidate
regions), this is redundant work.

**Proposal:** Add an LRU cache (via the `lru` crate) keyed by `(sid, start, end)` inside
`SeqIndexDB`. Cache size would be user-configurable (default: 256 MB). This is particularly
important for the Python API where users call `get_seq_by_id()` in tight loops.

```rust
use lru::LruCache;

pub struct SeqIndexDB {
    // existing fields ...
    seq_cache: Option<LruCache<(u32, u32, u32), Vec<u8>>>,
}
```

### 2c. Shimmer index sparsification for repeat-rich regions

`frag_map: FxHashMap<(u64, u64), Vec<FragmentSignature>>` accumulates all hits for every
shimmer pair. In repeat-rich regions (centromeres, segmental duplications) a single shimmer
pair can match thousands of locations, making queries slow and memory usage unpredictable.

**Proposal:** Add a configurable hit-count cap (`max_occurrences`, analogous to `minimap2`'s
`-f` flag). Shimmer pairs with hit counts above the cap are excluded from the index or
down-sampled. This would improve both query speed and memory usage for large pangenomes.

```rust
pub struct ShmmrSpec {
    pub w: u32,
    pub k: u32,
    pub r: u32,
    pub min_span: u32,
    pub sketch: bool,
    pub max_occ: Option<u32>,  // new field; None = no cap
}
```

### 2d. Resolve the `Mmap` clone limitation

`pgr bundle decomp.rs` line 259 has a `TODO: fix this` comment explaining that the
database must be rebuilt because `Mmap` is not `Clone`. The workaround is a full re-read
of the file, which is expensive.

**Proposal:** Wrap the `Mmap` in an `Arc<Mmap>`. `Arc<T>` is `Clone` regardless of whether
`T` is, so the same mapped region can be shared across threads without duplication.

```rust
use std::sync::Arc;
use memmap2::Mmap;

pub struct FragFileDB {
    mmap: Arc<Mmap>,
    // ...
}
```

### 2e. Trait-based backend abstraction

The current `Backend` enum in `ext.rs` requires every new storage format to add a new match
arm in every function. This creates scattered code and breaks when compiled without the
`with_agc` feature flag.

**Proposal:** Extract a `SequenceStore` trait and implement it for each backend. Functions
that currently match on `Backend` would accept `&dyn SequenceStore` (or `&impl SequenceStore`
for static dispatch).

```rust
pub trait SequenceStore {
    fn get_seq_by_id(&self, sid: u32) -> Result<Vec<u8>>;
    fn get_sub_seq(&self, sid: u32, start: u32, end: u32) -> Result<Vec<u8>>;
    fn seq_info(&self) -> &SeqInfoMap;
    // ...
}

impl SequenceStore for AgcDB { ... }
impl SequenceStore for FastxDB { ... }
impl SequenceStore for FrgDB { ... }
```

### 2f. Parallelize principal bundle computation

The Dijkstra + DFS in `graph_utils.rs` runs sequentially. For large pangenomes the bundle
decomposition step (`pgr bundle decomp`) becomes the bottleneck.

**Proposal:** Partition the MAP-graph by connected component and process each component in
parallel via `rayon::par_iter()`. Components are independent by definition, so no
synchronisation is needed beyond collecting results.

---

## 3. Code Quality Improvements

### 3a. Centralise magic numbers as named constants

Replace scattered literals with a constants module:

```rust
// pgr-db/src/constants.rs  (new file)

/// Default SHIMMER parameters
pub const DEFAULT_W: u32 = 80;
pub const DEFAULT_K: u32 = 56;
pub const DEFAULT_R: u32 = 4;
pub const DEFAULT_MIN_SPAN: u32 = 64;

/// Hash seed for shimmer pair hashing
pub const SHIMMER_HASH_SEED: u64 = 0xAD12CF59;

/// Gzip magic bytes
pub const GZIP_MAGIC: [u8; 2] = [0x1F, 0x8B];

/// I/O buffer sizes
pub const IO_BUF_SMALL: usize = 1 << 12;   //  4 KB
pub const IO_BUF_DEFAULT: usize = 1 << 14; // 16 KB
pub const IO_BUF_LARGE: usize = 1 << 16;   // 64 KB

/// Fragment reconstruction threshold
pub const FRAG_MIN_LEN: u32 = 128;
```

### 3b. Eliminate code duplication in backend match arms

The four-arm match on `Backend` in `get_seq_by_id`, `get_sub_seq_by_id`, and related
functions (ext.rs lines 363–499) is nearly identical across all functions. Once the
`SequenceStore` trait is in place (see §2e), these match arms collapse to a single
trait-method call.

Until the trait refactor is done, extract a helper:

```rust
fn resolve_seq_info<'a>(
    seq_info: &'a Option<SeqInfoMap>,
    sid: u32,
) -> Result<&'a (String, Option<String>, u32)> {
    seq_info
        .as_ref()
        .ok_or(PgrError::BackendNotInitialised)?
        .get(&sid)
        .ok_or(PgrError::SeqIdNotFound(sid))
}
```

### 3c. Fix or delete the commented-out modules

`gff_db.rs` and `seqs2variants.rs` are commented out in `lib.rs`. They contain
real logic but are untested and unmaintained:

- If the functionality is needed: uncomment, add tests, and wire into the CLI.
- If it is not needed: delete the files to avoid confusion.

### 3d. Resolve all `TODO: Test the output properly` stubs in `aln.rs`

Four test functions (lines 670, 739, 759, 782) contain no assertions. Each should either
be given real input/output test cases or removed.

### 3e. Add integration tests for CLI tools

There are no end-to-end tests for `pgr align alnmap`, `pgr variant diploid-vcf`, or any other
binary. A minimal integration test suite should:

1. Use a small synthetic reference + query FASTA (< 10 kb each).
2. Run each binary and check return code and output file existence.
3. Spot-check key output values (e.g. number of variant records, VCF record count).

These tests can live in `pgr-bin/tests/` and run via `cargo test --package pgr-bin`.

### 3f. Make `ShmmrSpec` implement `Copy`

`ShmmrSpec` is a small struct of five primitive fields. It is currently cloned in multiple
places (e.g. `pgr-tk/src/lib.rs` lines 261, 274, 1188). Adding `#[derive(Copy, Clone)]`
removes the explicit `.clone()` calls and makes the intent clearer.

### 3g. Address `clippy::type_complexity` suppressions

Two `#[allow(clippy::type_complexity)]` annotations in `ext.rs` (lines 561, 985) are
accompanied by `TODO: Define the type for readability`. These should be resolved with
named type aliases:

```rust
// Instead of suppressing clippy, name the type:
pub type BundleIdOrientationMap = FxHashMap<u32, Vec<(u32, u8, u32, u32)>>;
pub type AnchorHitList = Vec<(u32, u32, u32, u32, u8)>;
```

---

## Prioritised Roadmap

### Phase 1 — Error handling foundation (1–2 weeks)

1. Add `thiserror` dependency to `pgr-db`
2. Define `PgrError` and `pgr_db::Result<T>`
3. Harden `ext.rs` backend match arms
4. Fix PyO3 boundary in `pgr-tk`
5. Fix AGC FFI boundary in `agc_io.rs`

**Outcome:** No more panics at the library/API boundary. Users get descriptive errors
instead of backtraces.

### Phase 2 — Code quality (2–3 weeks, can overlap with Phase 1)

1. Add `constants.rs` and replace all magic numbers
2. Resolve or delete `gff_db.rs` and `seqs2variants.rs`
3. Add `Copy` to `ShmmrSpec`
4. Fix the four test stubs in `aln.rs`
5. Resolve `type_complexity` suppressions with named aliases
6. Add a minimal integration test suite for `pgr align alnmap` and `pgr variant diploid-vcf`

**Outcome:** Fewer silent surprises, more maintainable codebase.

### Phase 3 — Architecture (longer term, design before implementing)

1. `Arc<Mmap>` fix for `pgr bundle decomp` (quick win, 1 day)
2. Shimmer hit-count cap (`max_occ` in `ShmmrSpec`) (1 week)
3. Sequence cache in `SeqIndexDB` (1 week)
4. Trait-based backend abstraction (`SequenceStore`) (2 weeks)
5. Parallel principal bundle computation (1 week)
6. Arrow IPC binary output for `pgr align alnmap` (2–3 weeks)

**Outcome:** Better performance at large scale, cleaner extension points for new storage
backends and output formats.

---

## Quick-win Checklist

These can be done in an hour or less each, independently of the larger phases:

- [ ] `Arc<Mmap>` in `pgr bundle decomp.rs` — removes the `TODO: fix this`
- [ ] `#[derive(Copy)]` on `ShmmrSpec` — removes 3 `.clone()` calls
- [ ] Delete or restore `gff_db.rs` and `seqs2variants.rs`
- [ ] Replace `GZIP_MAGIC` literal `[0x1F, 0x8B]` with a named constant in all three files
- [ ] `clippy::type_complexity` → two named type aliases in `ext.rs`
- [ ] `.partial_cmp(...).unwrap()` in `aln.rs:21` → `f64::total_cmp` (stabilised in Rust 1.62)
