pub mod cov;
pub mod cov2;
pub mod fetch;
pub mod seqs;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct QueryCmd {
    #[command(subcommand)]
    pub cmd: QuerySubCmd,
}

#[derive(Subcommand)]
pub enum QuerySubCmd {
    /// Query a PGR-TK pangenome database and fetch matching target sequences
    Seqs(seqs::Args),
    /// List or fetch sequences from a PGR-TK database by region
    Fetch(fetch::Args),
    /// Compare shimmer-pair coverage between two sets of sequences
    Cov(cov::Args),
    /// Compare shimmer-pair coverage (mmap-backed, improved version)
    Cov2(cov2::Args),
}

pub fn run(cmd: QueryCmd) {
    match cmd.cmd {
        QuerySubCmd::Seqs(args) => seqs::run(args).expect("query seqs failed"),
        QuerySubCmd::Fetch(args) => fetch::run(args).expect("query fetch failed"),
        QuerySubCmd::Cov(args) => cov::run(args),
        QuerySubCmd::Cov2(args) => cov2::run(args),
    }
}
