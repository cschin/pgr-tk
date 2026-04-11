# agc-rs: Replace `par_sort_unstable` with Parallel Radix Sort

## Background

Profiling GRCh38 compression shows the `determine_splitters` phase takes **48 s**
wall time — sorting ~3.1 B `u64` k-mers before the singleton scan.

The C++ AGC uses **raduls 2.0** (`RadixSortMSD`), a parallel MSD (Most Significant
Digit) radix sort on raw bytes:

```cpp
raduls::RadixSortMSD(data, tmp=nullptr, n_recs,
                     rec_size=8, key_size=8, n_threads);
```

agc-rs currently uses rayon's `par_sort_unstable()` (parallel pdqsort).

### Algorithmic complexity comparison

| Algorithm | Complexity | Cache behaviour |
|---|---|---|
| `par_sort_unstable` (pdqsort) | O(n log n) | poor for n > L3 cache |
| Radix sort (8 passes over u64) | O(8n) | sequential byte sweeps, cache-friendly |

For n = 3.1 B:
- pdqsort: ~3.1B × 32 ≈ 99 B comparisons
- Radix sort: ~3.1B × 8 = 24.8 B byte-level passes — ~4× fewer operations

---

## Chosen Crate: `voracious_radix_sort`

Two strong candidates exist on crates.io:

| Crate | Version | Notes |
|---|---|---|
| `voracious_radix_sort` | 1.2.0 | Multi-thread Voracious sort; native `u64` support; described as "state of the art" |
| `rdst` | 0.20.14 | Parallel unstable radix sort; supports arbitrary key functions |

**Choice: `voracious_radix_sort`** — it directly sorts `u64` with a multi-threaded
implementation matching raduls' usage pattern most closely. `rdst` is also viable
as a fallback.

---

## Implementation Plan

### Step 1 — Add dependency

In `agc-rs/Cargo.toml`:

```toml
voracious_radix_sort = "1.2.0"
```

### Step 2 — Replace sort in `collect_sorted_singletons`

Location: `agc-rs/src/compressor.rs`, function `collect_sorted_singletons`.

Current code (line ~63):
```rust
all_kmers.par_sort_unstable();
```

Replacement:
```rust
use voracious_radix_sort::{RadixSort, Radixable};
all_kmers.voracious_mt_sort(rayon::current_num_threads());
```

The `u64` type implements `Radixable` in `voracious_radix_sort`, so no wrapper
is needed. `voracious_mt_sort` uses the thread count to drive internal parallelism.

### Step 3 — Verify correctness

The singleton compaction loop (`remove_non_singletons` equivalent) requires the
Vec to be **sorted in ascending order**. `voracious_mt_sort` for `u64` sorts
ascending by value — same contract as `par_sort_unstable`. Verify with:

```
cargo test -p agc-rs --release
```

All four existing round-trip tests must pass (E. coli append, three-genome, etc.).

### Step 4 — Benchmark

Re-run the GRCh38 create command with `\time -v` and compare the
`determine_splitters` phase time from the progress output:

```
[  X.XXs]   splitters: 49223 found (Y.YYs)
```

Expected: the sort phase (currently ~48 s) should drop to ~10–15 s on a
10-thread machine (4× algorithmic improvement × parallel execution).

### Step 5 — Fallback if voracious is not a win

If `voracious_mt_sort` underperforms (e.g. due to memory bandwidth saturation
rather than comparison cost), try `rdst`:

```toml
rdst = "0.20.14"
```

```rust
use rdst::RadixSort;
all_kmers.radix_sort_unstable();  // single-threaded; use par_radix_sort for MT
```

---

## Files Changed

| File | Change |
|---|---|
| `agc-rs/Cargo.toml` | Add `voracious_radix_sort = "1.2.0"` |
| `agc-rs/src/compressor.rs` | Replace `par_sort_unstable()` with `voracious_mt_sort(n_threads)` in `collect_sorted_singletons` |

No other files need to change. The rest of the pipeline (splitter finding,
compression, DB writes) is unaffected.

---

## Expected Outcome

| Phase | Before | Expected after |
|---|---|---|
| Singleton k-mer sort (48s in profile) | `par_sort_unstable` O(n log n) | radix sort O(8n), ~10–15 s |
| Total `determine_splitters` | ~48 s | ~15–20 s |

The DB write system time (589 s in profile) is a separate issue driven by
page-fault pressure from Vec allocations during compression, not the sort.
