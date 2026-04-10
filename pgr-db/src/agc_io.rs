use crate::fasta_io::SeqRec;
use agc_rs::decompressor::AgcFile as InnerAgcFile;
use memmap2::Mmap;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

pub type ShmmrToFragMapLocation = FxHashMap<(u64, u64), (usize, usize)>;
use std::io;
use std::path::Path;
use std::sync::Mutex;

#[derive(Debug, Clone)]
pub struct AGCSample {
    pub name: String,
    pub contigs: Vec<(String, usize)>, // (contig_name, length)
}

/// An open agc-rs archive backed by SQLite.
///
/// The inner `AgcFile` (which wraps a `rusqlite::Connection`) is wrapped in a
/// `Mutex` so that `AGCFile` can be shared across threads by Rayon iterators.
pub struct AGCFile {
    pub filepath: String,
    inner: Mutex<InnerAgcFile>,
    pub samples: Vec<AGCSample>,
    pub ctg_lens: FxHashMap<(String, String), usize>,
    sample_ctg: Vec<(String, String)>,
    pub prefetching: bool,       // API-compat no-op
    pub number_iter_thread: usize, // API-compat no-op
}

impl std::fmt::Debug for AGCFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AGCFile")
            .field("filepath", &self.filepath)
            .field("n_samples", &self.samples.len())
            .finish()
    }
}

pub struct AGCSeqDB {
    pub agc_file: AGCFile,
    pub frag_location_map: ShmmrToFragMapLocation,
    pub frag_map_file: Mmap,
}

pub struct AGCFileIter<'a> {
    agc_file: &'a AGCFile,
    current_ctg: usize,
}

impl AGCFile {
    /// Open an `.agcrs` (SQLite) archive at `filepath`.
    pub fn new(filepath: String) -> Result<Self, io::Error> {
        if !Path::new(&filepath).exists() {
            return Err(io::Error::new(io::ErrorKind::NotFound, filepath));
        }

        let inner = InnerAgcFile::open_readonly(Path::new(&filepath))
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

        let sample_names = inner
            .list_samples()
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

        let mut samples = Vec::new();
        let mut ctg_lens_vec = Vec::new();
        let mut sample_ctg = Vec::new();

        for sample_name in &sample_names {
            let contig_names = inner
                .list_contigs(sample_name)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

            let mut ctgs: Vec<(String, usize)> = Vec::new();
            for contig_name in &contig_names {
                let len = inner
                    .contig_len(sample_name, contig_name)
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?
                    as usize;
                ctg_lens_vec.push(((sample_name.clone(), contig_name.clone()), len));
                sample_ctg.push((sample_name.clone(), contig_name.clone()));
                ctgs.push((contig_name.clone(), len));
            }
            samples.push(AGCSample {
                name: sample_name.clone(),
                contigs: ctgs,
            });
        }

        let ctg_lens = FxHashMap::from_iter(ctg_lens_vec);

        Ok(Self {
            filepath,
            inner: Mutex::new(inner),
            samples,
            ctg_lens,
            sample_ctg,
            prefetching: true,
            number_iter_thread: 8,
        })
    }

    pub fn set_iter_thread(&mut self, number_iter_thread: usize) {
        self.number_iter_thread = number_iter_thread;
    }

    pub fn set_prefetching(&mut self, prefetching: bool) {
        self.prefetching = prefetching;
    }

    /// Return the sub-sequence `[bgn, end)` (0-based, half-open) as ASCII bytes.
    pub fn get_sub_seq(
        &self,
        sample_name: String,
        ctg_name: String,
        bgn: usize,
        end: usize,
    ) -> Vec<u8> {
        self.inner
            .lock()
            .expect("agc_io mutex poisoned")
            .contig_seq(&sample_name, &ctg_name, bgn as u64, end as u64)
            .unwrap_or_else(|e| {
                panic!("get_sub_seq({sample_name}/{ctg_name} {bgn}..{end}): {e}")
            })
    }

    /// Return the full sequence of a contig as ASCII bytes.
    pub fn get_seq(&self, sample_name: String, ctg_name: String) -> Vec<u8> {
        self.inner
            .lock()
            .expect("agc_io mutex poisoned")
            .full_contig(&sample_name, &ctg_name)
            .unwrap_or_else(|e| panic!("get_seq({sample_name}/{ctg_name}): {e}"))
    }

    /// Return the full list of `(sample_name, contig_name)` pairs in archive order.
    ///
    /// Used by batch processing to slice the list into chunks without loading
    /// any sequence data.
    pub fn sample_ctg_list(&self) -> &[(String, String)] {
        &self.sample_ctg
    }

    /// Fetch all contigs in parallel using one read-only SQLite connection per
    /// rayon task.
    ///
    /// SQLite WAL mode supports unlimited concurrent readers, so opening
    /// multiple `open_readonly` connections to the same file is safe and
    /// avoids the single-Mutex serialization of `AGCFileIter`.  For a human
    /// genome (50 chromosomes) this turns ~50 sequential decompression calls
    /// into a single parallel batch.
    pub fn par_fetch_seqs(&self) -> Vec<SeqRec> {
        self.par_fetch_seqs_batch(&self.sample_ctg)
    }

    /// Decompress only the given `(sample_name, contig_name)` pairs in parallel.
    ///
    /// This is the core primitive for batched memory-bounded indexing: callers
    /// slice `sample_ctg_list()` into chunks and call this function once per
    /// chunk, keeping peak decompression memory proportional to the chunk size
    /// rather than the total archive size.
    pub fn par_fetch_seqs_batch(&self, batch: &[(String, String)]) -> Vec<SeqRec> {
        let filepath = self.filepath.clone();
        batch
            .par_iter()
            .map(|(sample_name, ctg_name)| {
                // Each rayon task opens its own read-only connection.
                let conn = InnerAgcFile::open_readonly(Path::new(&filepath))
                    .unwrap_or_else(|e| panic!("par_fetch_seqs open {filepath}: {e}"));
                let seq = conn
                    .full_contig(sample_name, ctg_name)
                    .unwrap_or_else(|e| {
                        panic!("par_fetch_seqs {sample_name}/{ctg_name}: {e}")
                    });
                SeqRec {
                    source: Some(sample_name.clone()),
                    id: ctg_name.as_bytes().to_vec(),
                    seq,
                }
            })
            .collect()
    }
}

impl<'a> IntoIterator for &'a AGCFile {
    type Item = io::Result<SeqRec>;
    type IntoIter = AGCFileIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        AGCFileIter::new(self)
    }
}

impl<'a> AGCFileIter<'a> {
    pub fn new(agc_file: &'a AGCFile) -> Self {
        AGCFileIter {
            agc_file,
            current_ctg: 0,
        }
    }
}

impl<'a> Iterator for AGCFileIter<'a> {
    type Item = io::Result<SeqRec>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_ctg >= self.agc_file.sample_ctg.len() {
            return None;
        }
        let (sample_name, ctg_name) = &self.agc_file.sample_ctg[self.current_ctg];
        self.current_ctg += 1;

        let seq = match self
            .agc_file
            .inner
            .lock()
            .expect("agc_io mutex poisoned")
            .full_contig(sample_name, ctg_name)
        {
            Ok(s) => s,
            Err(e) => {
                return Some(Err(io::Error::new(io::ErrorKind::Other, e.to_string())))
            }
        };

        Some(Ok(SeqRec {
            source: Some(sample_name.clone()),
            id: ctg_name.as_bytes().to_vec(),
            seq,
        }))
    }
}
