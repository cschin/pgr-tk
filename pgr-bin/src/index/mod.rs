pub mod mdb;
pub mod shmmr_count;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct IndexCmd {
    #[command(subcommand)]
    pub cmd: IndexSubCmd,
}

#[derive(Subcommand)]
pub enum IndexSubCmd {
    /// Create a minimizer database (.mdbi/.mdbv/.midx) from an AGC archive
    Mdb(mdb::Args),
    /// Count and compare shimmer occurrences across three sequence sets
    ShmmrCount(shmmr_count::Args),
}

pub fn run(cmd: IndexCmd) {
    match cmd.cmd {
        IndexSubCmd::Mdb(args) => mdb::run(args),
        IndexSubCmd::ShmmrCount(args) => shmmr_count::run(args).expect("shmmr-count failed"),
    }
}
