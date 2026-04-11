use agc_rs::compressor::Compressor;
use agc_rs::decompressor::AgcFile;
use agc_rs::merge::merge_archives;
use agc_rs::segment::Params;
use anyhow::{bail, Context};
use clap::{Parser, Subcommand};
use std::io::{self, Write};
use std::path::Path;

#[derive(Parser)]
#[command(name = "agc-rs", about = "Assembled Genomes Compressor (Rust)")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Create a new archive from one or more FASTA files
    Create {
        /// Output archive path
        #[arg(short, long)]
        output: String,
        /// Sample name for the first (or only) input
        #[arg(long, default_value = "sample")]
        sample: String,
        /// Input FASTA files (plain or .gz); one sample per file
        inputs: Vec<String>,
    },
    /// Append a FASTA file to an existing archive
    Append {
        /// Path to the existing archive
        archive: String,
        /// Sample name for the new data
        #[arg(long)]
        sample: String,
        /// Input FASTA file
        input: String,
    },
    /// Append multiple FASTA files to an existing archive, building the delta index only once
    BatchAppend {
        /// Path to the existing archive (must already contain a reference sample)
        archive: String,
        /// Input files in the form `sample_name:fasta_path` (plain or .gz)
        #[arg(required = true)]
        inputs: Vec<String>,
    },
    /// Merge multiple archives that share the same reference into one
    Merge {
        /// Output archive path
        #[arg(short, long)]
        output: String,
        /// Source archives to merge (must share the same reference / splitter set)
        #[arg(required = true)]
        inputs: Vec<String>,
    },
    /// List samples, or contigs within a sample
    List {
        /// Path to the archive
        archive: String,
        /// If given, list contigs for this sample instead of listing samples
        #[arg(long)]
        contigs: Option<String>,
    },
    /// Extract a contig sequence (full or subrange) and write FASTA to stdout
    Get {
        /// Path to the archive
        archive: String,
        /// Query: `sample/contig` or `sample/contig:start-end` (0-based, half-open)
        query: String,
    },
    /// Print archive statistics
    Info {
        /// Path to the archive
        archive: String,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        // -------------------------------------------------------------------
        Commands::Create {
            output,
            sample,
            inputs,
        } => {
            if inputs.is_empty() {
                bail!("at least one input file is required");
            }
            let out_path = Path::new(&output);
            let mut c = Compressor::create(out_path, Params::default())
                .with_context(|| format!("creating archive at {output}"))?;

            for (idx, input) in inputs.iter().enumerate() {
                // Generate a distinct sample name for each additional file.
                let name = if idx == 0 {
                    sample.clone()
                } else {
                    format!("{}_{}", sample, idx + 1)
                };
                c.add_fasta(Path::new(input), &name)
                    .with_context(|| format!("compressing {input} as sample {name}"))?;
            }
            c.finish()?;
            eprintln!("Created archive: {output}");
        }

        // -------------------------------------------------------------------
        Commands::Append {
            archive,
            sample,
            input,
        } => {
            let mut c = Compressor::append(Path::new(&archive))
                .with_context(|| format!("opening archive {archive}"))?;
            c.add_fasta(Path::new(&input), &sample)
                .with_context(|| format!("appending {input} as sample {sample}"))?;
            c.finish()?;
            eprintln!("Appended sample '{sample}' to {archive}");
        }

        // -------------------------------------------------------------------
        Commands::BatchAppend { archive, inputs } => {
            if inputs.is_empty() {
                bail!("at least one input is required");
            }
            // Parse "sample_name:fasta_path" pairs.
            let parsed: Vec<(String, String)> = inputs
                .iter()
                .map(|s| {
                    let colon = s.find(':').with_context(|| {
                        format!("input '{}' is not in 'sample_name:fasta_path' format", s)
                    })?;
                    Ok((s[..colon].to_string(), s[colon + 1..].to_string()))
                })
                .collect::<anyhow::Result<Vec<_>>>()?;

            let mut c = Compressor::append(Path::new(&archive))
                .with_context(|| format!("opening archive {archive}"))?;

            let path_name_pairs: Vec<(&Path, &str)> = parsed
                .iter()
                .map(|(name, path)| (Path::new(path.as_str()), name.as_str()))
                .collect();

            c.add_fasta_batch(&path_name_pairs)
                .with_context(|| format!("batch-appending to {archive}"))?;
            c.finish()?;
            eprintln!(
                "Batch-appended {} samples to {archive}",
                path_name_pairs.len()
            );
        }

        // -------------------------------------------------------------------
        Commands::Merge { output, inputs } => {
            if inputs.is_empty() {
                bail!("at least one source archive is required");
            }
            let source_paths: Vec<&Path> = inputs.iter().map(|s| Path::new(s.as_str())).collect();
            merge_archives(Path::new(&output), &source_paths)
                .with_context(|| format!("merging {} archives into {output}", inputs.len()))?;
            eprintln!("Merged {} archives into {output}", inputs.len());
        }

        // -------------------------------------------------------------------
        Commands::List { archive, contigs } => {
            let agc = AgcFile::open(Path::new(&archive))
                .with_context(|| format!("opening archive {archive}"))?;

            match contigs {
                None => {
                    for name in agc.list_samples()? {
                        println!("{}", name);
                    }
                }
                Some(sample) => {
                    for name in agc.list_contigs(&sample)? {
                        println!("{}", name);
                    }
                }
            }
        }

        // -------------------------------------------------------------------
        Commands::Get { archive, query } => {
            let agc = AgcFile::open(Path::new(&archive))
                .with_context(|| format!("opening archive {archive}"))?;

            // Parse `sample/contig` or `sample/contig:start-end`.
            let (sample, contig, range) = parse_query(&query)?;
            let seq = match range {
                None => agc.full_contig(&sample, &contig)?,
                Some((s, e)) => agc.contig_seq(&sample, &contig, s, e)?,
            };

            let stdout = io::stdout();
            let mut out = stdout.lock();
            writeln!(out, ">{}/{}", sample, contig)?;
            for chunk in seq.chunks(60) {
                out.write_all(chunk)?;
                writeln!(out)?;
            }
        }

        // -------------------------------------------------------------------
        Commands::Info { archive } => {
            let agc = AgcFile::open(Path::new(&archive))
                .with_context(|| format!("opening archive {archive}"))?;

            let n_samples = agc.n_samples()?;
            let samples = agc.list_samples()?;

            let mut total_contigs: usize = 0;
            let mut total_bases: u64 = 0;

            for sample in &samples {
                let contigs = agc.list_contigs(sample)?;
                for contig in &contigs {
                    total_bases += agc.contig_len(sample, contig)?;
                }
                total_contigs += contigs.len();
            }

            println!("Archive:       {}", archive);
            println!("Samples:       {}", n_samples);
            println!("Total contigs: {}", total_contigs);
            println!("Total bases:   {}", total_bases);
        }
    }

    Ok(())
}

/// Parse a query string of the form `sample/contig` or `sample/contig:start-end`.
///
/// Returns `(sample, contig, Option<(start, end)>)`.
fn parse_query(query: &str) -> anyhow::Result<(String, String, Option<(u64, u64)>)> {
    // Split on the first '/'.
    let slash = query
        .find('/')
        .context("query must be 'sample/contig' or 'sample/contig:start-end'")?;
    let sample = query[..slash].to_string();
    let rest = &query[slash + 1..];

    // Optionally split on ':' for range.
    if let Some(colon) = rest.find(':') {
        let contig = rest[..colon].to_string();
        let range_str = &rest[colon + 1..];
        let dash = range_str.find('-').context("range must be 'start-end'")?;
        let start: u64 = range_str[..dash].parse().context("invalid start")?;
        let end: u64 = range_str[dash + 1..].parse().context("invalid end")?;
        Ok((sample, contig, Some((start, end))))
    } else {
        Ok((sample, rest.to_string(), None))
    }
}
