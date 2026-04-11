pub mod chr_aln;

use clap::{Parser, Subcommand};

#[derive(Parser)]
pub struct PlotCmd {
    #[command(subcommand)]
    pub cmd: PlotSubCmd,
}

#[derive(Subcommand)]
pub enum PlotSubCmd {
    /// Generate chromosome alignment SVG plots from an alnmap JSON file
    ChrAln(chr_aln::Args),
}

pub fn run(cmd: PlotCmd) {
    match cmd.cmd {
        PlotSubCmd::ChrAln(args) => chr_aln::run(args).expect("plot chr-aln failed"),
    }
}
