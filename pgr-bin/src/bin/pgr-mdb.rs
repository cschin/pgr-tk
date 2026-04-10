const VERSION_STRING: &str = env!("VERSION_STRING");

use clap::{self, CommandFactory, Parser};
use pgr_db::agc_io::AGCFile;
use pgr_db::seq_db;

/// Create a pgr minimizer database (.mdbi/.mdbv/.midx) from a single AGC archive
#[derive(Parser, Debug)]
#[clap(name = "pgr-mdb")]
#[clap(author, version)]
#[clap(about, long_about = None)]
struct CmdOptions {
    /// output file prefix (produces <prefix>.mdbi, <prefix>.mdbv, <prefix>.midx)
    #[clap(long)]
    prefix: String,
    /// path to the AGC archive (.agcrs)
    #[clap(long)]
    agcrs_input: String,
    /// minimizer window size
    #[clap(long, short, default_value_t = 80)]
    w: u32,
    /// minimizer k-mer size
    #[clap(long, short, default_value_t = 56)]
    k: u32,
    /// sparse minimizer (shimmer) reduction factor
    #[clap(long, short, default_value_t = 4)]
    r: u32,
    /// min span for neighboring minimizers
    #[clap(long, short, default_value_t = 64)]
    min_span: u32,
    /// use sketch k-mer instead of minimizer
    #[clap(short, long)]
    sketch: bool,
}

fn main() {
    CmdOptions::command().version(VERSION_STRING).get_matches();
    let args = CmdOptions::parse();

    let shmmr_spec = pgr_db::shmmrutils::ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: args.sketch,
    };

    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec);
    let agcfile = AGCFile::new(args.agcrs_input.clone()).expect("failed to open AGC archive");
    let _ = sdb.load_index_from_agcfile(agcfile);
    sdb.write_shmmr_map_index(args.prefix, Some(&args.agcrs_input))
        .expect("failed to write index");
}
