/// Convert an ASCII nucleotide byte to its 2-bit encoding.
///
/// | Base | Bits |
/// |------|------|
/// | A    |  0   |
/// | C    |  1   |
/// | G    |  2   |
/// | T    |  3   |
///
/// Returns `None` for any byte that is not one of `A`, `C`, `G`, `T`
/// (case-sensitive).
#[inline]
pub fn base_to_bits(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Compute the 2-bit reverse complement of a `k`-mer stored in the low
/// `2*k` bits of a `u64`.
#[cfg(test)]
///
/// The algorithm:
/// 1. Complement all bits (A↔T, C↔G).
/// 2. Reverse the order of the `k` 2-bit pairs.
#[inline]
fn rev_comp_kmer(kmer: u64, k: u8) -> u64 {
    // Complement: flip all bits, then mask to the relevant 2*k bits.
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let comp = (!kmer) & mask;

    // Reverse the order of the k 2-bit pairs.
    let mut rev: u64 = 0;
    let mut tmp = comp;
    for _ in 0..k {
        rev = (rev << 2) | (tmp & 0x3);
        tmp >>= 2;
    }
    rev
}

/// A sliding-window k-mer that simultaneously tracks the forward encoding
/// and its reverse complement.
///
/// Call [`Kmer::push`] for each incoming base. Once [`Kmer::full`] returns
/// `true` the canonical (lexicographically smaller) encoding is available via
/// [`Kmer::canonical`].
#[derive(Debug, Clone)]
pub struct Kmer {
    forward: u64,
    revcomp: u64,
    len: u8,
    max_len: u8,
}

impl Kmer {
    /// Create a new, empty k-mer builder for k-mers of length `max_len`.
    ///
    /// # Panics
    ///
    /// Panics if `max_len == 0` or `max_len > 32` (the latter would overflow
    /// a `u64`).
    pub fn new(max_len: u8) -> Self {
        assert!(max_len > 0 && max_len <= 32, "max_len must be in 1..=32");
        Self {
            forward: 0,
            revcomp: 0,
            len: 0,
            max_len,
        }
    }

    /// Slide in a new base, updating both the forward and reverse-complement
    /// encodings.
    ///
    /// Returns `true` on success. Returns `false` (and resets the internal
    /// state to empty) if `base` is not a valid ACGT character.
    pub fn push(&mut self, base: u8) -> bool {
        let bits = match base_to_bits(base) {
            Some(b) => b,
            None => {
                self.reset();
                return false;
            }
        };

        let k = self.max_len as u64;
        let mask = if self.max_len == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };

        // Forward: shift left by 2, add new base, keep only 2*k bits.
        self.forward = ((self.forward << 2) | bits) & mask;

        // Reverse complement: the new base's complement is prepended at the
        // most-significant end of the RC word.
        //   complement of bits: (bits ^ 3)  (A↔T: 0↔3, C↔G: 1↔2)
        let comp_bits = bits ^ 3;
        self.revcomp = (self.revcomp >> 2) | (comp_bits << (2 * (k - 1)));

        if self.len < self.max_len {
            self.len += 1;
        }

        true
    }

    /// Return the canonical encoding: the numerically smaller of the forward
    /// and reverse-complement k-mers.
    ///
    /// Only meaningful once [`Kmer::full`] returns `true`.
    #[inline]
    pub fn canonical(&self) -> u64 {
        self.forward.min(self.revcomp)
    }

    /// Return `true` if the forward encoding is the canonical one.
    #[inline]
    pub fn is_forward_canonical(&self) -> bool {
        self.forward <= self.revcomp
    }

    /// Return the raw forward encoding.
    #[inline]
    pub fn forward(&self) -> u64 {
        self.forward
    }

    /// Return the raw reverse-complement encoding.
    #[inline]
    pub fn revcomp(&self) -> u64 {
        self.revcomp
    }

    /// Return `true` once `max_len` bases have been consumed (i.e. the
    /// sliding window is full).
    #[inline]
    pub fn full(&self) -> bool {
        self.len == self.max_len
    }

    /// Push a 2-bit encoded base (0=A, 1=C, 2=G, 3=T).
    ///
    /// Returns `true` on success.  Returns `false` (and resets) if `bits >= 4`
    /// (i.e. the base is N or unknown).
    pub fn push_bits(&mut self, bits: u8) -> bool {
        if bits >= 4 {
            self.reset();
            return false;
        }
        let k = self.max_len as u64;
        let mask = if self.max_len == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * k)) - 1
        };
        self.forward = ((self.forward << 2) | bits as u64) & mask;
        let comp = bits as u64 ^ 3; // complement: A↔T (0↔3), C↔G (1↔2)
        self.revcomp = (self.revcomp >> 2) | (comp << (2 * (k - 1)));
        if self.len < self.max_len {
            self.len += 1;
        }
        true
    }

    /// Reset to the empty state, as if no bases had been pushed.
    pub fn reset(&mut self) {
        self.forward = 0;
        self.revcomp = 0;
        self.len = 0;
    }
}

/// Reverse complement of a 2-bit encoded sequence (A=0, C=1, G=2, T=3, N=4).
///
/// Reverses the order and complements each base (A↔T, C↔G). N (value ≥ 4)
/// is preserved as-is.
pub fn rev_comp_2bit_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| if b < 4 { 3 - b } else { b })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// The canonical k-mer must always equal `min(forward, revcomp)`.
    #[test]
    fn canonical_is_min_of_fwd_and_rc() {
        let seq = b"ACGTACGT";
        let k = 4u8;
        let mut km = Kmer::new(k);

        for &b in seq {
            km.push(b);
            if km.full() {
                let canon = km.canonical();
                assert_eq!(
                    canon,
                    km.forward().min(km.revcomp()),
                    "canonical should be min(fwd, rc)"
                );
            }
        }
    }

    /// After pushing exactly `max_len` valid bases the window must be full,
    /// and `len` must stay capped at `max_len` after additional pushes.
    #[test]
    fn push_fills_correctly() {
        let mut km = Kmer::new(4);
        assert!(!km.full());

        for &b in b"ACGT" {
            let ok = km.push(b);
            assert!(ok, "valid base should return true");
        }
        assert!(km.full(), "should be full after 4 bases");

        // One more push — still full, len does not overflow.
        km.push(b'A');
        assert!(km.full());
        assert_eq!(km.len, 4);
    }

    /// Verify the reverse complement of `ACGT` (k=4).
    ///
    /// ACGT  → 2-bit: 00 01 10 11  (forward)
    /// RC: complement of T,G,C,A in that order = A,C,G,T → same sequence.
    /// So for "ACGT" the forward and RC encodings must be equal.
    #[test]
    fn reverse_complement_acgt() {
        let mut km = Kmer::new(4);
        for &b in b"ACGT" {
            km.push(b);
        }
        assert!(km.full());

        // ACGT is its own reverse complement.
        assert_eq!(
            km.forward(),
            km.revcomp(),
            "ACGT should be its own reverse complement"
        );

        // Cross-check with the standalone helper.
        let rc = rev_comp_kmer(km.forward(), 4);
        assert_eq!(rc, km.forward());
    }

    /// Pushing an invalid base (`N`) must return `false` and reset the
    /// window so subsequent valid bases start fresh.
    #[test]
    fn invalid_base_returns_false() {
        let mut km = Kmer::new(4);

        // Push two valid bases first.
        assert!(km.push(b'A'));
        assert!(km.push(b'C'));
        assert_eq!(km.len, 2);

        // Invalid base.
        assert!(!km.push(b'N'), "N should return false");
        assert_eq!(km.len, 0, "state should be reset after invalid base");
        assert!(!km.full());

        // After reset, four more valid bases should fill the window again.
        for &b in b"TTGG" {
            assert!(km.push(b));
        }
        assert!(km.full());
    }

    /// Verify `base_to_bits` for all valid and a selection of invalid inputs.
    #[test]
    fn base_to_bits_mapping() {
        assert_eq!(base_to_bits(b'A'), Some(0));
        assert_eq!(base_to_bits(b'C'), Some(1));
        assert_eq!(base_to_bits(b'G'), Some(2));
        assert_eq!(base_to_bits(b'T'), Some(3));

        assert_eq!(base_to_bits(b'N'), None);
        assert_eq!(base_to_bits(b'a'), None, "lowercase not supported");
        assert_eq!(base_to_bits(b' '), None);
    }
}
