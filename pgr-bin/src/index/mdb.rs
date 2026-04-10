
use clap::{self, Parser};
use pgr_db::agc_io::AGCFile;
use pgr_db::seq_db;
use std::path::Path;

/// Create a pgr minimizer database (.mdbi/.mdbv/.midx) from a single AGC archive
#[derive(Parser, Debug)]
#[clap(about, long_about = None)]
pub struct Args {
    /// output file prefix (produces <prefix>.mdbi, <prefix>.mdbv, <prefix>.midx);
    /// defaults to the stem of --agcrs-input (e.g. "foo" for "foo.agcrs")
    #[clap(long)]
    pub prefix: Option<String>,
    /// path to the AGC archive (.agcrs)
    #[clap(long)]
    pub agcrs_input: String,
    /// minimizer window size
    #[clap(long, short, default_value_t = 80)]
    pub w: u32,
    /// minimizer k-mer size
    #[clap(long, short, default_value_t = 56)]
    pub k: u32,
    /// sparse minimizer (shimmer) reduction factor
    #[clap(long, short, default_value_t = 4)]
    pub r: u32,
    /// min span for neighboring minimizers
    #[clap(long, short, default_value_t = 64)]
    pub min_span: u32,
    /// use sketch k-mer instead of minimizer
    #[clap(short, long)]
    pub sketch: bool,
    /// number of whole samples (haplotypes) to process per batch.
    ///
    /// Reduces peak memory by decompressing and indexing one group of haplotypes
    /// at a time, writing sorted shard files that are merged at the end.  The
    /// default of 16 limits peak RAM to ~70 GB for a 100-haplotype human pangenome
    /// (vs. >300 GB).  Use 0 to process the entire archive at once (legacy
    /// behaviour, no sharding).
    #[clap(long, default_value_t = 16)]
    pub batch_size: usize,
}

pub fn run(args: Args) {

    let shmmr_spec = pgr_db::shmmrutils::ShmmrSpec {
        w: args.w,
        k: args.k,
        r: args.r,
        min_span: args.min_span,
        sketch: args.sketch,
    };

    let prefix = args.prefix.unwrap_or_else(|| {
        Path::new(&args.agcrs_input)
            .with_extension("")
            .to_string_lossy()
            .into_owned()
    });

    let mut sdb = seq_db::CompactSeqDB::new(shmmr_spec);
    let agcfile = AGCFile::new(args.agcrs_input.clone()).expect("failed to open AGC archive");

    if args.batch_size > 0 {
        sdb.load_index_from_agcfile_batched(
            &agcfile,
            args.batch_size,
            &prefix,
            Some(&args.agcrs_input),
        )
        .expect("failed to build batched index");
    } else {
        sdb.load_index_from_agcfile(agcfile)
            .expect("failed to load AGC index");
        sdb.write_shmmr_map_index(prefix, Some(&args.agcrs_input))
            .expect("failed to write index");
    }
}
