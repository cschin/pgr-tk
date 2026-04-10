# pgr-bin Subcommand Consolidation Plan

## Motivation

The 22 `pgr-*` binaries are currently independent executables with no explicit
relationship visible to the user.  They share common shimmer parameters
(`-w/-k/-r/--min-span`) that must be repeated identically on every invocation.
Consolidating them under a single `pgr` binary with subcommand groups:

- Makes the tool surface self-documenting (`pgr --help` gives a full map).
- Eliminates 22 separate install targets (one `pgr` binary in `$PATH`).
- Allows top-level flags (e.g. `--threads`) to be defined once and
  shared across subcommands.
- Groups related tools so users can discover them without reading docs.

---

## Proposed Command Hierarchy

```
pgr
├── index          — build / inspect minimizer indices
│   ├── mdb            (was pgr-mdb)
│   └── shmmr-count    (was pgr-shmmr-count)
│
├── align          — pairwise and multi-way alignment
│   ├── alnmap         (was pgr-alnmap)
│   ├── map-coord      (was pgr-map-coordinate)
│   └── liftover-gtf   (was pgr-liftover-gtf)
│
├── query          — search, fetch and compare sequences
│   ├── seqs           (was pgr-query)
│   ├── fetch          (was pgr-fetch-seqs)
│   ├── cov            (was pgr-compare-cov)
│   └── cov2           (was pgr-compare-cov2)
│
├── bundle         — principal bundle decomposition and analysis
│   ├── decomp         (was pgr-pbundle-decomp)
│   ├── aln            (was pgr-pbundle-aln)
│   ├── dist           (was pgr-pbundle-bed2dist)
│   ├── offset         (was pgr-pbundle-bed2offset)
│   ├── sort           (was pgr-pbundle-bed2sorted)
│   ├── svg            (was pgr-pbundle-bed2svg)
│   └── shmmr-dist     (was pgr-pbundle-shmmr2dist)
│
├── variant        — SV calling, merging and annotation
│   ├── diploid-vcf    (was pgr-generate-diploid-vcf)
│   ├── sv-analysis    (was pgr-generate-sv-analysis)
│   ├── merge-sv       (was pgr-merge-svcnd-bed)
│   ├── annotate-bed   (was pgr-annotate-bed-file)
│   └── annotate-vcf   (was pgr-annotate-vcf-file)
│
└── plot           — visualization outputs
    └── chr-aln        (was pgr-generate-chr-aln-plot)
```

---

## Implementation Strategy

### 1. Clap nested-subcommands pattern

Use clap's `derive` feature with `#[command(subcommand)]` at two levels:

```rust
// pgr-bin/src/bin/pgr.rs  (new single entry point)

#[derive(Parser)]
#[command(name = "pgr", about = "PGR-TK pangenome toolkit")]
struct Cli {
    #[command(subcommand)]
    group: Group,
}

#[derive(Subcommand)]
enum Group {
    Index(index::IndexCmd),
    Align(align::AlignCmd),
    Query(query::QueryCmd),
    Bundle(bundle::BundleCmd),
    Variant(variant::VariantCmd),
    Plot(plot::PlotCmd),
}
```

Each group module exposes its own `Subcommand` enum:

```rust
// pgr-bin/src/index.rs
#[derive(Parser)]
pub struct IndexCmd {
    #[command(subcommand)]
    pub cmd: IndexSubCmd,
}

#[derive(Subcommand)]
pub enum IndexSubCmd {
    Mdb(MdbArgs),
    ShmmrCount(ShmmrCountArgs),
}
```

### 2. Code movement

The current `pgr-*.rs` files each contain:
- A `CmdOptions` struct (clap `Parser`)
- One or more `fn` helpers
- A `main()` calling `CmdOptions::parse()`

For each old binary:

1. **Rename** `CmdOptions` → the subcommand `Args` struct (e.g. `MdbArgs`).
2. **Drop** `main()` and replace it with a `pub fn run(args: &MdbArgs)`.
3. **Move** the file to a module file, e.g.
   `pgr-bin/src/index/mdb.rs`.
4. The new `pgr.rs` entry point dispatches:
   ```rust
   Group::Index(c) => match c.cmd {
       IndexSubCmd::Mdb(a)        => index::mdb::run(&a),
       IndexSubCmd::ShmmrCount(a) => index::shmmr_count::run(&a),
   },
   Group::Bundle(c) => match c.cmd {
       BundleSubCmd::Decomp(a)    => bundle::decomp::run(&a),
       BundleSubCmd::Aln(a)       => bundle::aln::run(&a),
       BundleSubCmd::Dist(a)      => bundle::dist::run(&a),
       // ...
   },
   ```

### 3. Shared shimmer parameters

Seven tools (`mdb`, `bundle decomp`, `query seqs`, `align alnmap`, `bundle aln`,
`query cov`, `query cov2`) all accept identical `w / k / r / min_span / sketch`
flags.  Extract them into a shared struct:

```rust
// pgr-bin/src/shmmr_args.rs
#[derive(Args, Clone)]
pub struct ShmmrArgs {
    #[arg(long, short, default_value_t = 80)]
    pub w: u32,
    #[arg(long, short, default_value_t = 56)]
    pub k: u32,
    #[arg(long, short, default_value_t = 4)]
    pub r: u32,
    #[arg(long, default_value_t = 64)]
    pub min_span: u32,
    #[arg(long, short)]
    pub sketch: bool,
}
```

Use `#[command(flatten)]` to embed it in any `Args` struct that needs it:

```rust
pub struct MdbArgs {
    #[command(flatten)]
    pub shmmr: ShmmrArgs,
    #[arg(long)]
    pub agcrs_input: String,
    // ...
}
```

### 4. Cargo.toml changes

Remove the 22 individual `[[bin]]` entries and replace with one:

```toml
[[bin]]
name = "pgr"
path = "src/bin/pgr.rs"
```

Optionally keep thin compatibility shim binaries (one-liner `main` that calls
`pgr <group> <subcommand>`) for users with scripts depending on the old names,
but these can be deprecated immediately since the project is still pre-1.0.

### 5. Directory layout after refactor

```
pgr-bin/src/
├── bin/
│   └── pgr.rs             ← new single entry point
├── shmmr_args.rs          ← shared ShmmrArgs struct
├── index/
│   ├── mod.rs
│   ├── mdb.rs
│   └── shmmr_count.rs
├── align/
│   ├── mod.rs
│   ├── alnmap.rs
│   ├── map_coord.rs
│   └── liftover_gtf.rs
├── query/
│   ├── mod.rs
│   ├── seqs.rs
│   ├── fetch.rs
│   ├── cov.rs
│   └── cov2.rs
├── bundle/
│   ├── mod.rs
│   ├── decomp.rs
│   ├── aln.rs
│   ├── dist.rs
│   ├── offset.rs
│   ├── sort.rs
│   ├── svg.rs
│   └── shmmr_dist.rs
├── variant/
│   ├── mod.rs
│   ├── diploid_vcf.rs
│   ├── sv_analysis.rs
│   ├── merge_sv.rs
│   ├── annotate_bed.rs
│   └── annotate_vcf.rs
└── plot/
    ├── mod.rs
    └── chr_aln.rs
```

---

## Migration Path

| Phase | Work | Backward compat |
|---|---|---|
| 1 | Create `pgr.rs` entry + group modules; keep old `[[bin]]` entries compiling | Both `pgr index mdb` and `pgr-mdb` work |
| 2 | Move helper functions from old `.rs` files into module files; convert `main()` → `run()` | Both still work |
| 3 | Extract `ShmmrArgs` and flatten into affected `Args` structs | Both still work |
| 4 | Remove old `[[bin]]` entries; add deprecation notice to README | Only `pgr` works |

A single PR can cover all four phases; the refactor is mechanical enough to
do in one sweep once the design is agreed.

---

## Help output after consolidation

```
$ pgr --help
PGR-TK pangenome toolkit

Usage: pgr <COMMAND>

Commands:
  index    Build and inspect minimizer indices and shimmer databases
  align    Pairwise and multi-way genome alignment
  query    Search, fetch and compare sequences from a PGR-TK database
  bundle   Principal bundle decomposition and analysis
  variant  SV detection, merging and variant annotation
  plot     Generate alignment and SV visualisation plots
  help     Print this message or the help of the given subcommand(s)

$ pgr index --help
Build and inspect minimizer indices and shimmer databases

Usage: pgr index <COMMAND>

Commands:
  mdb          Create a minimizer database (.mdbi/.mdbv/.midx) from an AGC archive
  shmmr-count  Count and compare shimmer occurrences across three sequence sets

$ pgr bundle --help
Principal bundle decomposition and analysis

Usage: pgr bundle <COMMAND>

Commands:
  decomp      Principal bundle decomposition from a FASTA/FASTQ file
  aln         Generate alignment from bundle decomposition
  dist        Compute pairwise distance matrix from bundle BED
  offset      Compute offset alignment anchors from bundle BED
  sort        Sort contigs by bundle composition
  svg         Generate interactive SVG from bundle BED
  shmmr-dist  Compute pairwise similarity from shimmer index
```

---

## Key Design Decisions

- **No behaviour changes**: every `run()` function is the verbatim logic of the
  current `main()`.  The refactor is purely structural.
- **Single `VERSION_STRING`** embedded once in `pgr.rs` via `env!("VERSION_STRING")`.
- **`pgr-compare-cov` and `pgr-compare-cov2`** kept as `cov` and `cov2` under
  `query` for now; they can be merged into a single command with a `--v2` flag
  in a follow-up.
- **`pgr-liftover-gtf`** is logically part of the annotation pipeline but its
  primary input is an alignment database, so it sits under `align`.  This is
  a judgment call that can be revisited.
- **`plot chr-aln`** is a thin group with one subcommand.  If `pgr-generate-sv-analysis`
  also produces SVG output it could move here; for now it stays under `variant`
  because its primary output is analysis BED files, not plots.
