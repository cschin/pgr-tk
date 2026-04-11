pub mod alnmap;
pub mod liftover_gtf;
pub mod map_coord;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct AlignCmd {
    #[command(subcommand)]
    pub cmd: AlignSubCmd,
}

#[derive(Subcommand)]
pub enum AlignSubCmd {
    /// Align long contigs to a reference and identify potential SV regions
    Alnmap(alnmap::Args),
    /// Map query coordinates to target coordinates via an alnmap/alndb file
    MapCoord(map_coord::Args),
    /// Lift GTF transcript annotations from a reference to haplotype contigs
    LiftoverGtf(liftover_gtf::Args),
}

pub fn run(cmd: AlignCmd) {
    match cmd.cmd {
        AlignSubCmd::Alnmap(args) => alnmap::run(args).expect("alnmap failed"),
        AlignSubCmd::MapCoord(args) => map_coord::run(args).expect("map-coord failed"),
        AlignSubCmd::LiftoverGtf(args) => liftover_gtf::run(args).expect("liftover-gtf failed"),
    }
}
