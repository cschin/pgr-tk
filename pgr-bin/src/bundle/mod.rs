pub mod aln;
pub mod decomp;
pub mod dist;
pub mod offset;
pub mod shmmr_dist;
pub mod sort;
pub mod svg;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct BundleCmd {
    #[command(subcommand)]
    pub cmd: BundleSubCmd,
}

#[derive(Subcommand)]
pub enum BundleSubCmd {
    /// Generate the principal bundle decomposition via MAP Graph from a FASTA file
    Decomp(decomp::Args),
    /// Generate pairwise alignment from principal bundle decomposition
    Aln(aln::Args),
    /// Compute a pairwise distance matrix from a bundle BED file
    Dist(dist::Args),
    /// Compute offset alignment anchors from a bundle BED file
    Offset(offset::Args),
    /// Sort contigs by bundle composition from a bundle BED file
    Sort(sort::Args),
    /// Generate an interactive SVG visualisation from a bundle BED file
    Svg(svg::Args),
    /// Compute pairwise similarity from a shimmer index file
    ShmmrDist(shmmr_dist::Args),
}

pub fn run(cmd: BundleCmd) {
    match cmd.cmd {
        BundleSubCmd::Decomp(args) => decomp::run(args).expect("bundle decomp failed"),
        BundleSubCmd::Aln(args) => aln::run(args).expect("bundle aln failed"),
        BundleSubCmd::Dist(args) => dist::run(args).expect("bundle dist failed"),
        BundleSubCmd::Offset(args) => offset::run(args).expect("bundle offset failed"),
        BundleSubCmd::Sort(args) => sort::run(args).expect("bundle sort failed"),
        BundleSubCmd::Svg(args) => svg::run(args).expect("bundle svg failed"),
        BundleSubCmd::ShmmrDist(args) => shmmr_dist::run(args).expect("bundle shmmr-dist failed"),
    }
}
