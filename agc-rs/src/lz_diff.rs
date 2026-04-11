// *******************************************************************************************
// Rust port of AGC lz_diff (V1 and V2) from:
//   https://github.com/refresh-bio/agc  (MIT license)
//   Copyright(C) 2021-2022, S.Deorowicz, A.Danek, H.Li
// *******************************************************************************************

use crate::error::{AgcError, Result};

// ---------------------------------------------------------------------------
// Constants (USE_SPARSE_HT disabled: dense HT matches AGC's default build)
// ---------------------------------------------------------------------------
const HASHING_STEP: u32 = 4;
#[allow(dead_code)]
const DEFAULT_MIN_MATCH_LEN: u32 = 18;
const MIN_NRUN_LEN: u32 = 4;
const N_CODE: u8 = 4;
const N_RUN_STARTER: u8 = 30;
const INVALID_SYMBOL: u8 = 31;
const MAX_NO_TRIES: u32 = 64;
const MAX_LOAD_FACTOR: f64 = 0.7;
const EMPTY_KEY16: u16 = 0xFFFF;
const EMPTY_KEY32: u32 = 0xFFFF_FFFF;

// ---------------------------------------------------------------------------
// MurMur64 finaliser (MurMur3 variant)
// ---------------------------------------------------------------------------
#[inline]
fn murmur64(mut h: u64) -> u64 {
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    h
}

// ---------------------------------------------------------------------------
// Version selector
// ---------------------------------------------------------------------------
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LzVersion {
    V1,
    V2,
}

// ---------------------------------------------------------------------------
// Main struct
// ---------------------------------------------------------------------------
pub struct LzDiff {
    /// Padded reference (original bytes + key_len × INVALID_SYMBOL at end).
    reference: Vec<u8>,
    ht32: Vec<u32>,
    ht16: Vec<u16>,
    ht_size: u64,
    ht_mask: u64,
    key_len: u32,
    min_match_len: u32,
    short_ht_ver: bool,
    index_ready: bool,
    version: LzVersion,
}

impl LzDiff {
    // -----------------------------------------------------------------------
    pub fn new(version: LzVersion, min_match_len: u32) -> Self {
        let key_len = min_match_len - HASHING_STEP + 1;
        LzDiff {
            reference: Vec::new(),
            ht32: Vec::new(),
            ht16: Vec::new(),
            ht_size: 0,
            ht_mask: 0,
            key_len,
            min_match_len,
            short_ht_ver: false,
            index_ready: false,
            version,
        }
    }

    // -----------------------------------------------------------------------
    /// Build the reference and its hash index from 2-bit encoded DNA.
    pub fn prepare(&mut self, reference: &[u8]) {
        self.short_ht_ver = (reference.len() as u64) / (HASHING_STEP as u64) < 65535;
        self.prepare_gen(reference);
        self.prepare_index();
    }

    // -----------------------------------------------------------------------
    /// Encode `text` (2-bit DNA) against the prepared reference.
    ///
    /// Requires that `prepare()` has already been called (as done by
    /// `segment::lz_from_ref_blob`).  The method takes `&self` because
    /// neither `encode_v1` nor `encode_v2` mutates any field of `LzDiff`.
    pub fn encode(&self, text: &[u8]) -> Vec<u8> {
        debug_assert!(self.index_ready, "LzDiff::encode called before prepare()");
        match self.version {
            LzVersion::V1 => self.encode_v1(text),
            LzVersion::V2 => self.encode_v2(text),
        }
    }

    // -----------------------------------------------------------------------
    /// Decode an `encoded` buffer, using the reference stored in `self`.
    /// The caller must supply the *unpadded* reference; here we derive it
    /// from our padded copy (strip the key_len padding bytes at the end).
    pub fn decode(&self, encoded: &[u8]) -> Result<Vec<u8>> {
        // The unpadded reference length.
        let ref_len = if self.reference.len() >= self.key_len as usize {
            self.reference.len() - self.key_len as usize
        } else {
            self.reference.len()
        };
        let reference = &self.reference[..ref_len];
        match self.version {
            LzVersion::V1 => self.decode_v1(reference, encoded),
            LzVersion::V2 => self.decode_v2(reference, encoded),
        }
    }

    // -----------------------------------------------------------------------
    // prepare_gen: copy reference, pad with INVALID_SYMBOL × key_len
    // -----------------------------------------------------------------------
    fn prepare_gen(&mut self, src: &[u8]) {
        self.reference = src.to_vec();
        let kl = self.key_len as usize;
        self.reference
            .resize(self.reference.len() + kl, INVALID_SYMBOL);
    }

    // -----------------------------------------------------------------------
    // prepare_index: count eligible positions, size the HT, fill it.
    //
    // With USE_SPARSE_HT the eligible positions are those i (stepped by
    // hashing_step) where cnt_mod == key_len_mod AND no_prev_valid >= key_len.
    // -----------------------------------------------------------------------
    fn prepare_index(&mut self) {
        let key_len = self.key_len as usize;
        let hashing_step = HASHING_STEP as usize;
        let key_len_mod = key_len % hashing_step;

        // Count eligible entries
        let mut ht_count: u64 = 0;
        let mut no_prev_valid: u32 = 0;
        let mut cnt_mod: usize = 0;

        for &c in &self.reference {
            if c < 4 {
                no_prev_valid += 1;
            } else {
                no_prev_valid = 0;
            }
            cnt_mod += 1;
            if cnt_mod == hashing_step {
                cnt_mod = 0;
            }
            if cnt_mod == key_len_mod && no_prev_valid >= self.key_len {
                ht_count += 1;
            }
        }

        // Round up to next power-of-2 after dividing by load factor
        let mut ht_size = ((ht_count as f64) / MAX_LOAD_FACTOR) as u64;
        if ht_size < 2 {
            ht_size = 2;
        }
        // Strip all but the highest set bit (floor to power-of-2)
        while ht_size & (ht_size - 1) != 0 {
            ht_size &= ht_size - 1;
        }
        // Now ht_size is the largest power-of-2 <= original; multiply by 2
        ht_size <<= 1;

        self.ht_size = ht_size;
        self.ht_mask = ht_size - 1;

        if self.short_ht_ver {
            self.ht16 = vec![EMPTY_KEY16; ht_size as usize];
            self.make_index16();
        } else {
            self.ht32 = vec![EMPTY_KEY32; ht_size as usize];
            self.make_index32();
        }
        self.index_ready = true;
    }

    // -----------------------------------------------------------------------
    fn make_index16(&mut self) {
        let ref_size = self.reference.len();
        let key_len = self.key_len as usize;
        let hashing_step = HASHING_STEP as usize;
        let ht_mask = self.ht_mask;
        let max_no_tries = MAX_NO_TRIES as usize;

        let mut i = 0usize;
        while i + key_len < ref_size {
            let x = Self::get_code_static(&self.reference[i..], key_len);
            if x != !0u64 {
                let pos = (murmur64(x) & ht_mask) as usize;
                for j in 0..max_no_tries {
                    let slot = (pos + j) & (ht_mask as usize);
                    if self.ht16[slot] == EMPTY_KEY16 {
                        self.ht16[slot] = (i / hashing_step) as u16;
                        break;
                    }
                }
            }
            i += hashing_step;
        }
    }

    // -----------------------------------------------------------------------
    fn make_index32(&mut self) {
        let ref_size = self.reference.len();
        let key_len = self.key_len as usize;
        let hashing_step = HASHING_STEP as usize;
        let ht_mask = self.ht_mask;
        let max_no_tries = MAX_NO_TRIES as usize;

        let mut i = 0usize;
        while i + key_len < ref_size {
            let x = Self::get_code_static(&self.reference[i..], key_len);
            if x != !0u64 {
                let pos = (murmur64(x) & ht_mask) as usize;
                for j in 0..max_no_tries {
                    let slot = (pos + j) & (ht_mask as usize);
                    if self.ht32[slot] == EMPTY_KEY32 {
                        self.ht32[slot] = (i / hashing_step) as u32;
                        break;
                    }
                }
            }
            i += hashing_step;
        }
    }

    // -----------------------------------------------------------------------
    // get_code: read key_len bytes, pack as 2-bit per base.
    // Returns !0u64 if any byte > 3.
    // -----------------------------------------------------------------------
    #[inline]
    fn get_code_static(s: &[u8], key_len: usize) -> u64 {
        let mut x: u64 = 0;
        for i in 0..key_len {
            let b = s[i];
            if b > 3 {
                return !0u64;
            }
            x = (x << 2) | (b as u64);
        }
        x
    }

    #[inline]
    fn get_code(&self, s: &[u8]) -> u64 {
        Self::get_code_static(s, self.key_len as usize)
    }

    // -----------------------------------------------------------------------
    // get_Nrun_len: how many consecutive N_CODE bytes starting at s?
    // Returns 0 if fewer than 3 consecutive N_CODE at the start.
    // -----------------------------------------------------------------------
    #[inline]
    fn get_nrun_len(s: &[u8]) -> u32 {
        if s.len() < 3 || s[0] != N_CODE || s[1] != N_CODE || s[2] != N_CODE {
            return 0;
        }
        let mut len = 3usize;
        while len < s.len() && s[len] == N_CODE {
            len += 1;
        }
        len as u32
    }

    // -----------------------------------------------------------------------
    // find_best_match (u16 table)
    // Returns (matched, ref_pos, len_bck, len_fwd)
    // -----------------------------------------------------------------------
    fn find_best_match16(
        &self,
        ht_pos: usize,
        text: &[u8],
        text_offset: usize, // position in text (for back-extension limit)
        max_len: usize,
        no_prev_literals: u32,
    ) -> (bool, u32, u32, u32) {
        let mut len_fwd: u32 = 0;
        let mut len_bck: u32 = 0;
        let mut ref_pos: u32 = 0;
        let ht_mask = self.ht_mask as usize;
        let hashing_step = HASHING_STEP as usize;
        let ref_data = &self.reference;
        let s = &text[text_offset..];

        for i in 0..MAX_NO_TRIES as usize {
            let slot = (ht_pos + i) & ht_mask;
            if self.ht16[slot] == EMPTY_KEY16 {
                break;
            }
            let h_pos = (self.ht16[slot] as usize) * hashing_step;
            let p = &ref_data[h_pos..];

            // Forward extension
            let mut f_len = 0usize;
            while f_len < max_len && f_len < p.len() && s[f_len] == p[f_len] {
                f_len += 1;
            }

            // Backward extension
            let back_limit = no_prev_literals.min(h_pos as u32) as usize;
            let mut b_len = 0usize;
            // s[-b_len-1] vs p[-b_len-1]
            while b_len < back_limit {
                let si = text_offset.wrapping_sub(b_len + 1);
                let pi = h_pos.wrapping_sub(b_len + 1);
                if text[si] != ref_data[pi] {
                    break;
                }
                b_len += 1;
            }

            if b_len + f_len > (len_bck + len_fwd) as usize {
                len_bck = b_len as u32;
                len_fwd = f_len as u32;
                ref_pos = h_pos as u32;
            }
        }

        (
            len_bck + len_fwd >= self.min_match_len,
            ref_pos,
            len_bck,
            len_fwd,
        )
    }

    // -----------------------------------------------------------------------
    // find_best_match (u32 table)
    // -----------------------------------------------------------------------
    fn find_best_match32(
        &self,
        ht_pos: usize,
        text: &[u8],
        text_offset: usize,
        max_len: usize,
        no_prev_literals: u32,
    ) -> (bool, u32, u32, u32) {
        let mut len_fwd: u32 = 0;
        let mut len_bck: u32 = 0;
        let mut ref_pos: u32 = 0;
        let ht_mask = self.ht_mask as usize;
        let hashing_step = HASHING_STEP as usize;
        let ref_data = &self.reference;
        let s = &text[text_offset..];

        for i in 0..MAX_NO_TRIES as usize {
            let slot = (ht_pos + i) & ht_mask;
            if self.ht32[slot] == EMPTY_KEY32 {
                break;
            }
            let h_pos = (self.ht32[slot] as usize) * hashing_step;
            let p = &ref_data[h_pos..];

            let mut f_len = 0usize;
            while f_len < max_len && f_len < p.len() && s[f_len] == p[f_len] {
                f_len += 1;
            }

            let back_limit = no_prev_literals.min(h_pos as u32) as usize;
            let mut b_len = 0usize;
            while b_len < back_limit {
                let si = text_offset.wrapping_sub(b_len + 1);
                let pi = h_pos.wrapping_sub(b_len + 1);
                if text[si] != ref_data[pi] {
                    break;
                }
                b_len += 1;
            }

            if b_len + f_len > (len_bck + len_fwd) as usize {
                len_bck = b_len as u32;
                len_fwd = f_len as u32;
                ref_pos = h_pos as u32;
            }
        }

        (
            len_bck + len_fwd >= self.min_match_len,
            ref_pos,
            len_bck,
            len_fwd,
        )
    }

    // -----------------------------------------------------------------------
    // find_best_match: dispatch to u16 or u32 table
    // -----------------------------------------------------------------------
    fn find_best_match(
        &self,
        ht_pos: usize,
        text: &[u8],
        text_offset: usize,
        max_len: usize,
        no_prev_literals: u32,
    ) -> (bool, u32, u32, u32) {
        if self.short_ht_ver {
            self.find_best_match16(ht_pos, text, text_offset, max_len, no_prev_literals)
        } else {
            self.find_best_match32(ht_pos, text, text_offset, max_len, no_prev_literals)
        }
    }

    // -----------------------------------------------------------------------
    // append_int: write decimal (with optional '-') into encoded
    // -----------------------------------------------------------------------
    #[inline]
    fn append_int(encoded: &mut Vec<u8>, x: i64) {
        if x == 0 {
            encoded.push(b'0');
            return;
        }
        let (neg, mut v) = if x < 0 {
            encoded.push(b'-');
            (true, (-(x as i128)) as u64)
        } else {
            (false, x as u64)
        };
        let _ = neg;
        // build digits in reverse
        let start = encoded.len();
        while v > 0 {
            encoded.push(b'0' + (v % 10) as u8);
            v /= 10;
        }
        encoded[start..].reverse();
    }

    // -----------------------------------------------------------------------
    // read_int: parse decimal integer from encoded slice, return new offset
    // -----------------------------------------------------------------------
    fn read_int(encoded: &[u8], pos: usize) -> Result<(i64, usize)> {
        let mut p = pos;
        if p >= encoded.len() {
            return Err(AgcError::LzDiff("read_int: unexpected end".into()));
        }
        let neg = encoded[p] == b'-';
        if neg {
            p += 1;
        }
        let mut x: i64 = 0;
        let start = p;
        while p < encoded.len() && encoded[p] >= b'0' && encoded[p] <= b'9' {
            x = x * 10 + (encoded[p] - b'0') as i64;
            p += 1;
        }
        if p == start {
            return Err(AgcError::LzDiff(format!(
                "read_int: no digits at pos {pos} (byte={})",
                encoded[pos]
            )));
        }
        if neg {
            x = -x;
        }
        Ok((x, p))
    }

    // -----------------------------------------------------------------------
    // encode_literal: push 'A' + c
    // -----------------------------------------------------------------------
    #[inline]
    fn encode_literal(c: u8, encoded: &mut Vec<u8>) {
        encoded.push(b'A' + c);
    }

    // -----------------------------------------------------------------------
    // encode_Nrun: push [N_RUN_STARTER, decimal(len-MIN_NRUN_LEN), N_CODE]
    // -----------------------------------------------------------------------
    #[inline]
    fn encode_nrun(len: u32, encoded: &mut Vec<u8>) {
        encoded.push(N_RUN_STARTER);
        Self::append_int(encoded, (len - MIN_NRUN_LEN) as i64);
        encoded.push(N_CODE);
    }

    // -----------------------------------------------------------------------
    // is_literal: *p in 'A'..'A'+20 or *p == '!'
    // -----------------------------------------------------------------------
    #[inline]
    fn is_literal(b: u8) -> bool {
        (b >= b'A' && b <= b'A' + 20) || b == b'!'
    }

    // -----------------------------------------------------------------------
    // is_Nrun
    // -----------------------------------------------------------------------
    #[inline]
    fn is_nrun(b: u8) -> bool {
        b == N_RUN_STARTER
    }

    // -----------------------------------------------------------------------
    // decode_literal: returns (c, new_pos)
    // If b == '!' then c = b'!' (special V2 marker).
    // -----------------------------------------------------------------------
    fn decode_literal(encoded: &[u8], pos: usize) -> Result<(u8, usize)> {
        if pos >= encoded.len() {
            return Err(AgcError::LzDiff("decode_literal: unexpected end".into()));
        }
        let b = encoded[pos];
        if b == b'!' {
            Ok((b'!', pos + 1))
        } else {
            Ok((b - b'A', pos + 1))
        }
    }

    // -----------------------------------------------------------------------
    // decode_Nrun: returns (len, new_pos)
    // -----------------------------------------------------------------------
    fn decode_nrun(encoded: &[u8], pos: usize) -> Result<(u32, usize)> {
        // skip N_RUN_STARTER
        let p = pos + 1;
        let (raw_len, p) = Self::read_int(encoded, p)?;
        // skip N_CODE suffix
        if p >= encoded.len() || encoded[p] != N_CODE {
            return Err(AgcError::LzDiff(
                "decode_Nrun: missing N_CODE suffix".into(),
            ));
        }
        let p = p + 1;
        let len = (raw_len as u32) + MIN_NRUN_LEN;
        Ok((len, p))
    }

    // -----------------------------------------------------------------------
    // decode_match_v1: read `delta,raw_len.`
    // -----------------------------------------------------------------------
    fn decode_match_v1(
        encoded: &[u8],
        pos: usize,
        pred_pos: u32,
        min_match_len: u32,
    ) -> Result<(u32, u32, usize)> {
        let (delta, p) = Self::read_int(encoded, pos)?;
        // skip ','
        if p >= encoded.len() || encoded[p] != b',' {
            return Err(AgcError::LzDiff("decode_match_v1: expected ','".into()));
        }
        let p = p + 1;
        let ref_pos = ((pred_pos as i64) + delta) as u32;
        let (len, p) = if p < encoded.len() && encoded[p] == b'.' {
            (!0u32, p)
        } else {
            let (raw_len, p2) = Self::read_int(encoded, p)?;
            ((raw_len as u32) + min_match_len, p2)
        };
        // skip '.'
        if p >= encoded.len() || encoded[p] != b'.' {
            return Err(AgcError::LzDiff("decode_match_v1: expected '.'".into()));
        }
        let p = p + 1;
        Ok((ref_pos, len, p))
    }

    // -----------------------------------------------------------------------
    // decode_match_v2: read `delta[,raw_len].`
    // -----------------------------------------------------------------------
    fn decode_match_v2(
        encoded: &[u8],
        pos: usize,
        pred_pos: u32,
        min_match_len: u32,
    ) -> Result<(u32, u32, usize)> {
        let (delta, p) = Self::read_int(encoded, pos)?;
        let ref_pos = ((pred_pos as i64) + delta) as u32;
        let (len, p) = if p < encoded.len() && encoded[p] == b',' {
            let p = p + 1;
            let (raw_len, p2) = Self::read_int(encoded, p)?;
            // skip '.'
            if p2 >= encoded.len() || encoded[p2] != b'.' {
                return Err(AgcError::LzDiff("decode_match_v2: expected '.'".into()));
            }
            ((raw_len as u32) + min_match_len, p2 + 1)
        } else {
            // next byte must be '.'
            if p >= encoded.len() || encoded[p] != b'.' {
                return Err(AgcError::LzDiff("decode_match_v2: expected '.'".into()));
            }
            (!0u32, p + 1)
        };
        Ok((ref_pos, len, p))
    }

    // -----------------------------------------------------------------------
    // V1 encode
    // -----------------------------------------------------------------------
    fn encode_v1(&self, text: &[u8]) -> Vec<u8> {
        let text_size = text.len();
        let key_len = self.key_len as usize;
        let mut encoded: Vec<u8> = Vec::new();
        let mut i = 0usize;
        let mut pred_pos: u32 = 0;
        let mut no_prev_literals: u32 = 0;

        while i + key_len < text_size {
            let x = self.get_code(&text[i..]);

            if x == !0u64 {
                // Non-ACGT byte: check for N-run
                let nrun = Self::get_nrun_len(&text[i..]);
                if nrun >= MIN_NRUN_LEN {
                    Self::encode_nrun(nrun, &mut encoded);
                    i += nrun as usize;
                    no_prev_literals = 0;
                } else {
                    Self::encode_literal(text[i], &mut encoded);
                    i += 1;
                    pred_pos += 1;
                    no_prev_literals += 1;
                }
                continue;
            }

            let ht_pos = (murmur64(x) & self.ht_mask) as usize;
            let max_len = text_size - i;
            let (found, match_pos, len_bck, len_fwd) =
                self.find_best_match(ht_pos, text, i, max_len, no_prev_literals);

            if !found {
                Self::encode_literal(text[i], &mut encoded);
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
            } else {
                // Backward extension: undo previously emitted literals
                if len_bck > 0 {
                    let lb = len_bck as usize;
                    encoded.truncate(encoded.len() - lb);
                    i -= lb;
                    pred_pos -= len_bck;
                }
                let total_len = len_bck + len_fwd;
                let adjusted_match_pos = match_pos - len_bck;
                // encode_match_v1: delta,raw_len.
                let delta: i64 = (adjusted_match_pos as i64) - (pred_pos as i64);
                Self::append_int(&mut encoded, delta);
                encoded.push(b',');
                Self::append_int(&mut encoded, (total_len - self.min_match_len) as i64);
                encoded.push(b'.');

                pred_pos = adjusted_match_pos + total_len;
                i += (len_bck + len_fwd) as usize;
                no_prev_literals = 0;
            }
        }

        // Emit remaining literals
        while i < text_size {
            Self::encode_literal(text[i], &mut encoded);
            i += 1;
        }

        encoded
    }

    // -----------------------------------------------------------------------
    // V2 encode
    // -----------------------------------------------------------------------
    fn encode_v2(&self, text: &[u8]) -> Vec<u8> {
        let text_size = text.len();
        let key_len = self.key_len as usize;
        let ref_len = self.reference.len() - key_len; // unpadded reference length

        // V2: if text == reference (unpadded), return empty
        if text_size == ref_len && text == &self.reference[..ref_len] {
            return Vec::new();
        }

        let mut encoded: Vec<u8> = Vec::new();
        let mut i = 0usize;
        let mut pred_pos: u32 = 0;
        let mut no_prev_literals: u32 = 0;

        while i + key_len < text_size {
            let x = self.get_code(&text[i..]);

            if x == !0u64 {
                let nrun = Self::get_nrun_len(&text[i..]);
                if nrun >= MIN_NRUN_LEN {
                    Self::encode_nrun(nrun, &mut encoded);
                    i += nrun as usize;
                    no_prev_literals = 0;
                } else {
                    Self::encode_literal(text[i], &mut encoded);
                    i += 1;
                    pred_pos += 1;
                    no_prev_literals += 1;
                }
                continue;
            }

            let ht_pos = (murmur64(x) & self.ht_mask) as usize;
            let max_len = text_size - i;
            let (found, match_pos, len_bck, len_fwd) =
                self.find_best_match(ht_pos, text, i, max_len, no_prev_literals);

            if !found {
                Self::encode_literal(text[i], &mut encoded);
                i += 1;
                pred_pos += 1;
                no_prev_literals += 1;
            } else {
                // Backward extension: undo previously emitted literals
                if len_bck > 0 {
                    let lb = len_bck as usize;
                    encoded.truncate(encoded.len() - lb);
                    i -= lb;
                    pred_pos -= len_bck;
                }
                let total_len = len_bck + len_fwd;
                let adjusted_match_pos = match_pos - len_bck;

                // V2 special: if match_pos == pred_pos, look back in encoded
                // and replace literals that equal reference[match_pos - i] with '!'
                if adjusted_match_pos == pred_pos {
                    let e_size = encoded.len();
                    let mp = adjusted_match_pos as usize;
                    let mut k = 1usize;
                    while k < e_size && k < mp {
                        let eb = encoded[e_size - k];
                        if eb < b'A' || eb > b'Z' {
                            break;
                        }
                        if (eb - b'A') == self.reference[mp - k] {
                            encoded[e_size - k] = b'!';
                        }
                        k += 1;
                    }
                }

                // Check for match-to-end
                let encode_len = if i + total_len as usize == text_size
                    && adjusted_match_pos as usize + total_len as usize == ref_len
                {
                    !0u32
                } else {
                    total_len
                };

                // encode_match_v2
                let delta: i64 = (adjusted_match_pos as i64) - (pred_pos as i64);
                Self::append_int(&mut encoded, delta);
                if encode_len != !0u32 {
                    encoded.push(b',');
                    Self::append_int(&mut encoded, (encode_len - self.min_match_len) as i64);
                }
                encoded.push(b'.');

                pred_pos = adjusted_match_pos + total_len;
                i += total_len as usize;
                no_prev_literals = 0;
            }
        }

        // Emit remaining literals
        while i < text_size {
            Self::encode_literal(text[i], &mut encoded);
            i += 1;
        }

        encoded
    }

    // -----------------------------------------------------------------------
    // V1 decode
    // -----------------------------------------------------------------------
    fn decode_v1(&self, reference: &[u8], encoded: &[u8]) -> Result<Vec<u8>> {
        let mut decoded: Vec<u8> = Vec::new();
        let mut pred_pos: u32 = 0;
        let mut pos = 0usize;

        while pos < encoded.len() {
            let b = encoded[pos];
            if Self::is_literal(b) {
                let (c, new_pos) = Self::decode_literal(encoded, pos)?;
                decoded.push(c);
                pred_pos += 1;
                pos = new_pos;
            } else if Self::is_nrun(b) {
                let (len, new_pos) = Self::decode_nrun(encoded, pos)?;
                decoded.extend(std::iter::repeat(N_CODE).take(len as usize));
                pos = new_pos;
            } else {
                // match token
                let (ref_pos, len, new_pos) =
                    Self::decode_match_v1(encoded, pos, pred_pos, self.min_match_len)?;
                let actual_len = if len == !0u32 {
                    reference.len().saturating_sub(ref_pos as usize)
                } else {
                    len as usize
                };
                let rp = ref_pos as usize;
                if rp + actual_len > reference.len() {
                    return Err(AgcError::LzDiff(format!(
                        "decode_v1: match {rp}..{} out of reference len {}",
                        rp + actual_len,
                        reference.len()
                    )));
                }
                decoded.extend_from_slice(&reference[rp..rp + actual_len]);
                pred_pos = ref_pos + actual_len as u32;
                pos = new_pos;
            }
        }

        Ok(decoded)
    }

    // -----------------------------------------------------------------------
    // V2 decode
    // -----------------------------------------------------------------------
    fn decode_v2(&self, reference: &[u8], encoded: &[u8]) -> Result<Vec<u8>> {
        // V2 empty encoding means identical to reference
        if encoded.is_empty() {
            return Ok(reference.to_vec());
        }

        let mut decoded: Vec<u8> = Vec::new();
        let mut pred_pos: u32 = 0;
        let mut pos = 0usize;

        while pos < encoded.len() {
            let b = encoded[pos];
            if Self::is_literal(b) {
                let (c, new_pos) = Self::decode_literal(encoded, pos)?;
                let actual_c = if c == b'!' {
                    let pp = pred_pos as usize;
                    if pp >= reference.len() {
                        return Err(AgcError::LzDiff(format!(
                            "decode_v2: '!' pred_pos {pp} out of reference len {}",
                            reference.len()
                        )));
                    }
                    reference[pp]
                } else {
                    c
                };
                decoded.push(actual_c);
                pred_pos += 1;
                pos = new_pos;
            } else if Self::is_nrun(b) {
                let (len, new_pos) = Self::decode_nrun(encoded, pos)?;
                decoded.extend(std::iter::repeat(N_CODE).take(len as usize));
                pos = new_pos;
            } else {
                let (ref_pos, len, new_pos) =
                    Self::decode_match_v2(encoded, pos, pred_pos, self.min_match_len)?;
                let actual_len = if len == !0u32 {
                    reference.len().saturating_sub(ref_pos as usize)
                } else {
                    len as usize
                };
                let rp = ref_pos as usize;
                if rp + actual_len > reference.len() {
                    return Err(AgcError::LzDiff(format!(
                        "decode_v2: match {rp}..{} out of reference len {}",
                        rp + actual_len,
                        reference.len()
                    )));
                }
                decoded.extend_from_slice(&reference[rp..rp + actual_len]);
                pred_pos = ref_pos + actual_len as u32;
                pos = new_pos;
            }
        }

        Ok(decoded)
    }
}

// LzDiff holds only Vec<u8>/Vec<u16>/Vec<u32> and plain scalars.
// After prepare() is called all fields are effectively read-only during
// encode(), so it is safe to share across threads.
unsafe impl Send for LzDiff {}
unsafe impl Sync for LzDiff {}

// ===========================================================================
// Tests
// ===========================================================================
#[cfg(test)]
mod tests {
    use super::*;

    fn dna_encode(s: &[u8]) -> Vec<u8> {
        s.iter()
            .map(|&b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            })
            .collect()
    }

    fn dna_decode(v: &[u8]) -> Vec<u8> {
        v.iter()
            .map(|&b| match b {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            })
            .collect()
    }

    /// Build a moderately long test sequence so we exceed min_match_len (18).
    fn make_seq(seed: &[u8], repeat: usize) -> Vec<u8> {
        seed.iter()
            .cycle()
            .take(seed.len() * repeat)
            .copied()
            .collect()
    }

    fn round_trip(version: LzVersion, reference: &[u8], text: &[u8]) {
        let ref_enc = dna_encode(reference);
        let txt_enc = dna_encode(text);

        let mut lz = LzDiff::new(version, DEFAULT_MIN_MATCH_LEN);
        lz.prepare(&ref_enc);
        let encoded = lz.encode(&txt_enc);
        let decoded = lz.decode(&encoded).expect("decode failed");
        assert_eq!(
            dna_decode(&decoded),
            text.to_vec(),
            "round-trip mismatch for version {:?}",
            version
        );
    }

    #[test]
    fn round_trip_identical() {
        let seq = make_seq(b"ACGTACGTACGTACGT", 20);
        round_trip(LzVersion::V1, &seq, &seq);
        round_trip(LzVersion::V2, &seq, &seq);
    }

    #[test]
    fn round_trip_single_snp() {
        let reference = make_seq(b"ACGTACGTACGTACGT", 20);
        let mut text = reference.clone();
        text[100] = b'N'; // introduce a SNP in the middle
        round_trip(LzVersion::V1, &reference, &text);
        round_trip(LzVersion::V2, &reference, &text);
    }

    #[test]
    fn round_trip_insertion() {
        let reference = make_seq(b"ACGTACGTACGTACGT", 20);
        let mut text = reference.clone();
        // Insert 5 bytes at position 50
        let insert: Vec<u8> = b"TTTTT".to_vec();
        text.splice(50..50, insert);
        round_trip(LzVersion::V1, &reference, &text);
        round_trip(LzVersion::V2, &reference, &text);
    }

    #[test]
    fn round_trip_deletion() {
        let reference = make_seq(b"ACGTACGTACGTACGT", 20);
        let mut text = reference.clone();
        // Delete 10 bytes starting at position 60
        text.drain(60..70);
        round_trip(LzVersion::V1, &reference, &text);
        round_trip(LzVersion::V2, &reference, &text);
    }

    #[test]
    fn round_trip_no_shared_kmers() {
        // Reference is all-A; text is all-C  — no shared k-mers at all.
        let reference: Vec<u8> = vec![b'A'; 200];
        let text: Vec<u8> = vec![b'C'; 200];
        round_trip(LzVersion::V1, &reference, &text);
        round_trip(LzVersion::V2, &reference, &text);
    }

    #[test]
    fn round_trip_n_run() {
        // Sequence with 10 N's in the middle
        let mut reference = make_seq(b"ACGTACGTACGTACGT", 20);
        let mut text = reference.clone();
        for b in text[100..110].iter_mut() {
            *b = b'N';
        }
        // Also put N-run in reference to exercise that path
        for b in reference[50..58].iter_mut() {
            *b = b'N';
        }
        round_trip(LzVersion::V1, &reference, &text);
        round_trip(LzVersion::V2, &reference, &text);
    }

    #[test]
    fn v1_and_v2_both_decode_correctly() {
        let reference = make_seq(b"ACGTACGTACGTACGT", 20);
        let mut text = reference.clone();
        // A few scattered edits
        text[30] = b'T';
        text[80] = b'G';
        text[150] = b'A';

        let ref_enc = dna_encode(&reference);
        let txt_enc = dna_encode(&text);

        for version in [LzVersion::V1, LzVersion::V2] {
            let mut lz = LzDiff::new(version, DEFAULT_MIN_MATCH_LEN);
            lz.prepare(&ref_enc);
            let encoded = lz.encode(&txt_enc);
            let decoded = lz.decode(&encoded).expect("decode failed");
            assert_eq!(dna_decode(&decoded), text, "mismatch for {:?}", version);
        }
    }
}
