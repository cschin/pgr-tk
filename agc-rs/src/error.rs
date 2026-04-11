use thiserror::Error;

/// All errors that can occur in the agc-rs crate.
#[derive(Debug, Error)]
pub enum AgcError {
    /// A SQLite database error.
    #[error("database error: {0}")]
    Db(#[from] rusqlite::Error),

    /// An I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// A zstd compression/decompression error.
    #[error("zstd error: {0}")]
    Zstd(String),

    /// The requested sample was not found in the archive.
    #[error("sample not found: {0}")]
    SampleNotFound(String),

    /// The requested contig was not found for the given sample.
    #[error("contig not found: sample={sample}, contig={contig}")]
    ContigNotFound { sample: String, contig: String },

    /// A requested range exceeds the actual sequence length.
    #[error("range out of bounds: start={start}, end={end}, sequence length={len}")]
    RangeOutOfBounds {
        start: usize,
        end: usize,
        len: usize,
    },

    /// The on-disk schema version is not supported by this library.
    #[error("unsupported schema version: {0}")]
    UnsupportedVersion(u32),

    /// An error from the LZ-diff / delta encoding layer.
    #[error("lz-diff error: {0}")]
    LzDiff(String),

    /// A FASTA parsing error.
    #[error("FASTA parse error: {0}")]
    FastaParse(String),
}

/// Convenience alias used throughout the crate.
pub type Result<T> = std::result::Result<T, AgcError>;
