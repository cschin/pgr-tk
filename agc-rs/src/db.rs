use std::path::Path;

use rusqlite::{params, Connection, OpenFlags};

use crate::error::{AgcError, Result};

/// The schema version stored in the `meta` table that this library expects.
const SCHEMA_VERSION: u32 = 3;

/// SQL statements that create the full agc-rs schema.
const DDL: &str = "
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS sample (
    id   INTEGER PRIMARY KEY,
    name TEXT UNIQUE NOT NULL
);

CREATE TABLE IF NOT EXISTS contig (
    id        INTEGER PRIMARY KEY,
    sample_id INTEGER NOT NULL REFERENCES sample(id),
    name      TEXT    NOT NULL,
    length    INTEGER NOT NULL,
    UNIQUE(sample_id, name)
);

-- kmer_front / kmer_back: canonical k-mer values at the segment boundaries used
-- to split contigs (AGC splitter strategy).  NULL when no splitter is present
-- (e.g. first / last segment of a contig, or fallback-matched segments).
CREATE TABLE IF NOT EXISTS segment_group (
    id         INTEGER PRIMARY KEY,
    ref_data   BLOB NOT NULL,
    params     TEXT NOT NULL,
    kmer_front INTEGER,
    kmer_back  INTEGER,
    delta_blob BLOB
);

CREATE INDEX IF NOT EXISTS idx_seg_group_kmers
    ON segment_group(kmer_front, kmer_back)
    WHERE kmer_front IS NOT NULL AND kmer_back IS NOT NULL;

-- splitter: canonical singleton k-mers chosen from the reference genome used
-- to determine segment boundaries.  Persisted so that subsequent append calls
-- can split new samples at the same positions.
CREATE TABLE IF NOT EXISTS splitter (
    kmer INTEGER PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS segment (
    id          INTEGER PRIMARY KEY,
    contig_id   INTEGER NOT NULL REFERENCES contig(id),
    seg_order   INTEGER NOT NULL,
    group_id    INTEGER NOT NULL REFERENCES segment_group(id),
    in_group_id INTEGER NOT NULL,
    is_rev_comp BOOLEAN NOT NULL DEFAULT 0,
    raw_length  INTEGER NOT NULL,
    delta_data  BLOB,
    UNIQUE(contig_id, seg_order)
);

CREATE INDEX IF NOT EXISTS idx_seg ON segment(contig_id, seg_order);
";

/// Thin wrapper around a [`rusqlite::Connection`] that owns the agc-rs
/// database and enforces the expected schema version.
#[derive(Debug)]
pub struct AgcDb {
    conn: Connection,
}

impl AgcDb {
    // -----------------------------------------------------------------------
    // Construction helpers
    // -----------------------------------------------------------------------

    /// Create a brand-new database at `path`.
    ///
    /// Any existing file at `path` is removed first so that `create` always
    /// produces a fresh archive.  The full schema is applied and the current
    /// `schema_version` is written to the `meta` table.
    /// WAL mode is enabled so that concurrent readers do not block writers.
    pub fn create(path: &Path) -> Result<Self> {
        // Remove any stale archive (old schema, partial write, etc.).
        if path.exists() {
            std::fs::remove_file(path)?;
            // Also remove the WAL and SHM sidecar files if present.
            let _ = std::fs::remove_file(path.with_extension("agcrs-wal"));
            let _ = std::fs::remove_file(path.with_extension("agcrs-shm"));
            // Generic SQLite WAL/SHM sidecars.
            let wal = format!("{}-wal", path.display());
            let shm = format!("{}-shm", path.display());
            let _ = std::fs::remove_file(&wal);
            let _ = std::fs::remove_file(&shm);
        }
        let conn = Connection::open(path)?;

        // Enable WAL and performance pragmas.
        // synchronous=NORMAL: skip fsync on every WAL commit; only sync at
        // checkpoints — safe for non-power-fail scenarios and dramatically
        // reduces system-call overhead on large ingests.
        // cache_size=-131072: 128 MB page cache so recently written blobs are
        // not immediately evicted before the read-back during delta merge.
        conn.execute_batch(
            "PRAGMA journal_mode = WAL;
             PRAGMA synchronous = NORMAL;
             PRAGMA cache_size = -131072;
             PRAGMA temp_store = MEMORY;",
        )?;

        // Apply schema.
        conn.execute_batch(DDL)?;

        // Record schema version.
        conn.execute(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?1, ?2)",
            params!["schema_version", SCHEMA_VERSION.to_string()],
        )?;

        Ok(Self { conn })
    }

    /// Open an existing database at `path` for read-write access.
    ///
    /// Returns [`AgcError::UnsupportedVersion`] if the stored schema version
    /// does not match [`SCHEMA_VERSION`].
    pub fn open(path: &Path) -> Result<Self> {
        let conn = Connection::open(path)?;
        let db = Self { conn };
        db.check_schema_version()?;
        Ok(db)
    }

    /// Open an existing database at `path` in read-only mode.
    ///
    /// Returns [`AgcError::UnsupportedVersion`] if the stored schema version
    /// does not match [`SCHEMA_VERSION`].
    pub fn open_readonly(path: &Path) -> Result<Self> {
        let conn = Connection::open_with_flags(
            path,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        )?;
        let db = Self { conn };
        db.check_schema_version()?;
        Ok(db)
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Return a reference to the underlying [`rusqlite::Connection`].
    pub fn conn(&self) -> &Connection {
        &self.conn
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Read `schema_version` from the `meta` table and verify it equals
    /// [`SCHEMA_VERSION`].
    fn check_schema_version(&self) -> Result<()> {
        let version_str: String = self.conn.query_row(
            "SELECT value FROM meta WHERE key = 'schema_version'",
            [],
            |row| row.get(0),
        )?;

        let version: u32 = version_str
            .parse()
            .map_err(|_| AgcError::UnsupportedVersion(0))?;

        if version != SCHEMA_VERSION {
            return Err(AgcError::UnsupportedVersion(version));
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    /// Helper: create a temporary file path that does not yet exist on disk.
    ///
    /// We keep the [`NamedTempFile`] handle alive so the OS does not reclaim
    /// the path; we delete the file so SQLite creates a fresh database.
    fn temp_db_path() -> (NamedTempFile, std::path::PathBuf) {
        let f = NamedTempFile::new().expect("tempfile");
        let path = f.path().to_path_buf();
        // Remove the placeholder file — AgcDb::create opens a new DB.
        std::fs::remove_file(&path).ok();
        (f, path)
    }

    #[test]
    fn create_and_reopen() {
        let (_guard, path) = temp_db_path();

        // Create a new database.
        {
            let db = AgcDb::create(&path).expect("create");
            // Insert a sample to prove the schema is usable.
            db.conn()
                .execute("INSERT INTO sample (name) VALUES (?1)", params!["sample_A"])
                .expect("insert sample");
        }

        // Reopen read-write and verify the row survived.
        {
            let db = AgcDb::open(&path).expect("open");
            let name: String = db
                .conn()
                .query_row("SELECT name FROM sample WHERE id = 1", [], |r| r.get(0))
                .expect("query");
            assert_eq!(name, "sample_A");
        }

        // Reopen read-only — should also succeed.
        {
            let db = AgcDb::open_readonly(&path).expect("open_readonly");
            let count: i64 = db
                .conn()
                .query_row("SELECT COUNT(*) FROM sample", [], |r| r.get(0))
                .expect("count");
            assert_eq!(count, 1);
        }
    }

    #[test]
    fn wrong_schema_version_rejected() {
        let (_guard, path) = temp_db_path();

        // Create the database normally, then tamper with the version.
        {
            let db = AgcDb::create(&path).expect("create");
            db.conn()
                .execute(
                    "UPDATE meta SET value = '99' WHERE key = 'schema_version'",
                    [],
                )
                .expect("update");
        }

        // Re-opening must fail with UnsupportedVersion.
        match AgcDb::open(&path) {
            Err(AgcError::UnsupportedVersion(99)) => { /* expected */ }
            other => panic!("expected UnsupportedVersion(99), got {:?}", other),
        }
    }
}
