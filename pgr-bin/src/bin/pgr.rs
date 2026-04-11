const VERSION_STRING: &str = env!("VERSION_STRING");

use clap::{Parser, Subcommand};
use pgr_bin::{align, bundle, index, plot, query, variant};

/// PGR-TK pangenome toolkit
#[derive(Parser)]
#[clap(name = "pgr")]
#[clap(version = VERSION_STRING)]
#[clap(about = "PGR-TK pangenome toolkit", long_about = None)]
struct Cli {
    #[command(subcommand)]
    group: Group,
}

#[derive(Subcommand)]
enum Group {
    /// Build and inspect minimizer indices and shimmer databases
    Index(index::IndexCmd),
    /// Pairwise and multi-way genome alignment
    Align(align::AlignCmd),
    /// Search, fetch and compare sequences from a PGR-TK database
    Query(query::QueryCmd),
    /// Principal bundle decomposition and analysis
    Bundle(bundle::BundleCmd),
    /// SV detection, merging and variant annotation
    Variant(variant::VariantCmd),
    /// Generate alignment and SV visualisation plots
    Plot(plot::PlotCmd),
}

fn main() {
    let cli = Cli::parse();
    match cli.group {
        Group::Index(cmd) => index::run(cmd),
        Group::Align(cmd) => align::run(cmd),
        Group::Query(cmd) => query::run(cmd),
        Group::Bundle(cmd) => bundle::run(cmd),
        Group::Variant(cmd) => variant::run(cmd),
        Group::Plot(cmd) => plot::run(cmd),
    }
}
