//! [`MultiwordInt`]: a signed two's-complement integer stored as a **runtime**
//! number of little-endian `u64` limbs, holding the value in
//! `(−2^(64·len−1), 2^(64·len−1))`.
//!
//! This is the runtime-length sibling of [`Packed<L>`](crate::packed::Packed).
//! `Packed` fixes its width `L` at compile time, one width per recursion depth;
//! threading that const generic through the whole NTRU-solve recursion is
//! awkward, so the deep recursion uses this `Vec`-backed width instead. The value
//! semantics are identical: two's complement, base `2^64`, `limbs[0]` least
//! significant, and every binary operation requires both operands to share the
//! same limb count (the arena sizes one width per depth).
//!
//! The hot-loop arithmetic is exposed as free functions over `&[u64]` /
//! `&mut [u64]` slices (the `*_into` family) so that [`MultiwordPoly`] can run
//! them straight on its contiguous word-plane buffer with no per-operation
//! allocation. [`MultiwordInt`] is a thin owned wrapper over the same primitives,
//! used for unit tests and the `BigInt` boundary conversions.
//!
//! The word-level-CRT-into-limbs reconstruction follows fn-dsa's `zint31`.

use std::cmp::Ordering;

use num::{BigInt, Zero};

use crate::rns::{NttPrimeList, Rns};

// ---------------------------------------------------------------------------
// Slice primitives (allocation-free; operands share `out.len()` limbs).
// ---------------------------------------------------------------------------

/// True when the two's-complement value in `a` is negative.
pub(crate) fn is_negative(a: &[u64]) -> bool {
    matches!(a.last(), Some(top) if top >> 63 == 1)
}

/// `out = !a + 1` (two's-complement negation). `out` and `a` share a length.
pub(crate) fn neg_into(out: &mut [u64], a: &[u64]) {
    debug_assert_eq!(out.len(), a.len());
    let mut carry = 1u128;
    for i in 0..out.len() {
        let v = (!a[i]) as u128 + carry;
        out[i] = v as u64;
        carry = v >> 64;
    }
}

/// `a = !a + 1` (two's-complement negation, in place).
pub(crate) fn neg_in_place(a: &mut [u64]) {
    let mut carry = 1u128;
    for w in a.iter_mut() {
        let v = (!*w) as u128 + carry;
        *w = v as u64;
        carry = v >> 64;
    }
}

/// `out = a + b`, wrapping in two's complement.
pub(crate) fn add_into(out: &mut [u64], a: &[u64], b: &[u64]) {
    debug_assert!(out.len() == a.len() && a.len() == b.len());
    let mut carry = 0u128;
    for i in 0..out.len() {
        let s = a[i] as u128 + b[i] as u128 + carry;
        out[i] = s as u64;
        carry = s >> 64;
    }
}

/// `dst += src`, wrapping in two's complement.
pub(crate) fn add_into_self(dst: &mut [u64], src: &[u64]) {
    debug_assert_eq!(dst.len(), src.len());
    let mut carry = 0u128;
    for i in 0..dst.len() {
        let s = dst[i] as u128 + src[i] as u128 + carry;
        dst[i] = s as u64;
        carry = s >> 64;
    }
}

/// `dst -= src`, wrapping in two's complement.
pub(crate) fn sub_into_self(dst: &mut [u64], src: &[u64]) {
    debug_assert_eq!(dst.len(), src.len());
    let mut borrow = 0i128;
    for i in 0..dst.len() {
        let v = dst[i] as i128 - src[i] as i128 - borrow;
        dst[i] = v as u64;
        borrow = if v < 0 { 1 } else { 0 };
    }
}

/// `out = a - b`, wrapping in two's complement.
pub(crate) fn sub_into(out: &mut [u64], a: &[u64], b: &[u64]) {
    debug_assert!(out.len() == a.len() && a.len() == b.len());
    let mut borrow = 0i128;
    for i in 0..out.len() {
        let v = a[i] as i128 - b[i] as i128 - borrow;
        out[i] = v as u64;
        borrow = if v < 0 { 1 } else { 0 };
    }
}

/// `out = a * v`, truncated to `out.len()` limbs.
pub(crate) fn mul_u32_into(out: &mut [u64], a: &[u64], v: u32) {
    debug_assert_eq!(out.len(), a.len());
    let mut carry = 0u128;
    for i in 0..out.len() {
        let prod = a[i] as u128 * v as u128 + carry;
        out[i] = prod as u64;
        carry = prod >> 64;
    }
}

/// Unsigned schoolbook multiply `out = a * b`, treating all three slices as
/// non-negative magnitudes, truncated to `out.len()` limbs. No allocation. `a`,
/// `b`, and `out` may have independent lengths.
pub(crate) fn umul_into(out: &mut [u64], a: &[u64], b: &[u64]) {
    for w in out.iter_mut() {
        *w = 0;
    }
    let lo = out.len();
    for i in 0..a.len() {
        if i >= lo {
            break;
        }
        let ai = a[i] as u128;
        if ai == 0 {
            continue;
        }
        let mut carry = 0u128;
        let mut idx = i;
        for &bj in b.iter() {
            if idx >= lo {
                carry = 0;
                break;
            }
            let prod = ai * bj as u128 + out[idx] as u128 + carry;
            out[idx] = prod as u64;
            carry = prod >> 64;
            idx += 1;
        }
        // Propagate the leftover carry into the higher limbs.
        while carry != 0 && idx < lo {
            let v = out[idx] as u128 + carry;
            out[idx] = v as u64;
            carry = v >> 64;
            idx += 1;
        }
    }
}

/// Signed two's-complement multiply `out = a * b`, truncated to `out.len()`
/// limbs. Allocates two magnitude scratch buffers (sized to `a`/`b`); the
/// allocation-free hot path is [`MultiwordPoly`](crate::multiword_poly)'s
/// schoolbook multiply, which hoists the magnitudes out of the inner loop. Used
/// only by the differential test for the unsigned core ([`umul_into`]).
#[cfg(test)]
pub(crate) fn mul_into(out: &mut [u64], a: &[u64], b: &[u64]) {
    let na = is_negative(a);
    let nb = is_negative(b);
    let mut ma = a.to_vec();
    let mut mb = b.to_vec();
    if na {
        neg_in_place(&mut ma);
    }
    if nb {
        neg_in_place(&mut mb);
    }
    umul_into(out, &ma, &mb);
    if na ^ nb {
        neg_in_place(out);
    }
}

/// Logical left shift `out = a << bits`; bits beyond `out.len()·64` are dropped.
pub(crate) fn shl_into(out: &mut [u64], a: &[u64], bits: u32) {
    debug_assert_eq!(out.len(), a.len());
    let len = out.len();
    let limb = (bits / 64) as usize;
    let bit = bits % 64;
    for i in 0..len {
        let src = i as isize - limb as isize;
        if src < 0 {
            out[i] = 0;
            continue;
        }
        let src = src as usize;
        let mut word = a[src] as u128;
        if bit != 0 {
            let lower = if src >= 1 { a[src - 1] as u128 } else { 0 };
            word = (word << bit) | (lower >> (64 - bit));
        }
        out[i] = word as u64;
    }
}

/// Logical (unsigned) right shift `out = a >> bits`. Used on positive values
/// (e.g. forming `M/2` during CRT centering).
pub(crate) fn shr_unsigned_into(out: &mut [u64], a: &[u64], bits: u32) {
    debug_assert_eq!(out.len(), a.len());
    let len = out.len();
    let limb = (bits / 64) as usize;
    let bit = bits % 64;
    for i in 0..len {
        let src = i + limb;
        let mut word = if src < len { a[src] as u128 } else { 0 };
        if bit != 0 {
            let hi = if src + 1 < len { a[src + 1] as u128 } else { 0 };
            word = (word >> bit) | (hi << (64 - bit));
            word &= u64::MAX as u128;
        }
        out[i] = word as u64;
    }
}

/// Arithmetic right shift by `shift` bits (floor toward −∞), returning the low
/// 128 bits of the result as `i128`. Bit-for-bit identical to `BigInt >> shift`
/// whenever the true result fits `i128`; the windowed input to the `f64`
/// k-estimation FFT. No allocation.
pub(crate) fn shr_to_i128(a: &[u64], shift: u32) -> i128 {
    let len = a.len();
    let sign = if is_negative(a) { u64::MAX } else { 0 };
    let limb = (shift / 64) as usize;
    let bit = shift % 64;
    let get = |idx: usize| -> u64 { if idx < len { a[idx] } else { sign } };
    let mut lo = 0u128;
    for out_idx in 0..2 {
        let src = limb + out_idx;
        let mut word = get(src) as u128;
        if bit != 0 {
            let hi = get(src + 1) as u128;
            word = (word >> bit) | (hi << (64 - bit));
            word &= u64::MAX as u128;
        }
        lo |= word << (64 * out_idx);
    }
    lo as i128
}

/// Unsigned magnitude compare (treats both slices as non-negative magnitudes).
pub(crate) fn ucmp(a: &[u64], b: &[u64]) -> Ordering {
    debug_assert_eq!(a.len(), b.len());
    for i in (0..a.len()).rev() {
        match a[i].cmp(&b[i]) {
            Ordering::Equal => continue,
            ord => return ord,
        }
    }
    Ordering::Equal
}

/// Number of bits in the absolute value, matching [`BigInt::bits`]: `0` for the
/// value 0, otherwise `floor(log2(|v|)) + 1`. Allocates one scratch limb buffer
/// only when `a` is negative.
pub(crate) fn bit_length(a: &[u64]) -> u64 {
    if is_negative(a) {
        let mut mag = vec![0u64; a.len()];
        neg_into(&mut mag, a);
        bit_length_unsigned(&mag)
    } else {
        bit_length_unsigned(a)
    }
}

fn bit_length_unsigned(m: &[u64]) -> u64 {
    for i in (0..m.len()).rev() {
        if m[i] != 0 {
            return (i as u64) * 64 + (64 - m[i].leading_zeros() as u64);
        }
    }
    0
}

/// Write the value of `x` into `out` (two's complement, `out.len()` limbs).
/// `BigInt` boundary only; the magnitude must fit in `out.len()` limbs.
pub(crate) fn from_bigint_into(out: &mut [u64], x: &BigInt) {
    let neg = x.sign() == num::bigint::Sign::Minus;
    let mag = if neg { -x } else { x.clone() };
    let (_, words) = mag.to_u64_digits();
    debug_assert!(
        words.len() <= out.len(),
        "MultiwordInt too narrow: value needs {} limbs ({} bits), have {}",
        words.len(),
        x.bits(),
        out.len()
    );
    for w in out.iter_mut() {
        *w = 0;
    }
    for (i, w) in words.iter().take(out.len()).enumerate() {
        out[i] = *w;
    }
    if neg {
        neg_in_place(out);
    }
}

/// Reconstruct a [`BigInt`] from a two's-complement limb slice (`BigInt`
/// boundary only).
pub(crate) fn to_bigint(a: &[u64]) -> BigInt {
    let neg = is_negative(a);
    let mut mag = a.to_vec();
    if neg {
        neg_in_place(&mut mag);
    }
    let mut acc = BigInt::zero();
    for i in (0..a.len()).rev() {
        acc <<= 64;
        acc += BigInt::from(mag[i]);
    }
    if neg {
        -acc
    } else {
        acc
    }
}

// ---------------------------------------------------------------------------
// Owned wrapper.
// ---------------------------------------------------------------------------

/// A signed integer in `(−2^(64·len−1), 2^(64·len−1))`, two's complement, base
/// `2^64`, little-endian. `len ≥ 2` is assumed (so `from_i128` always fits).
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) struct MultiwordInt {
    limbs: Vec<u64>,
}

impl MultiwordInt {
    /// Zero in `len` limbs.
    pub(crate) fn zero(len: usize) -> Self {
        Self { limbs: vec![0; len] }
    }

    pub(crate) fn len(&self) -> usize {
        self.limbs.len()
    }

    pub(crate) fn limbs(&self) -> &[u64] {
        &self.limbs
    }

    /// Sign-extend a signed 128-bit value into `len` limbs (`len ≥ 2`).
    pub(crate) fn from_i128(v: i128, len: usize) -> Self {
        debug_assert!(len >= 2, "MultiwordInt needs ≥ 2 limbs for a 128-bit input");
        let lo = v as u128;
        let ext = if v < 0 { u64::MAX } else { 0 };
        let mut limbs = vec![ext; len];
        limbs[0] = lo as u64;
        limbs[1] = (lo >> 64) as u64;
        Self { limbs }
    }

    pub(crate) fn is_negative(&self) -> bool {
        is_negative(&self.limbs)
    }

    pub(crate) fn bit_length(&self) -> u64 {
        bit_length(&self.limbs)
    }

    pub(crate) fn shr_to_i128(&self, shift: u32) -> i128 {
        shr_to_i128(&self.limbs, shift)
    }

    pub(crate) fn shl(&self, bits: u32) -> Self {
        debug_assert!(
            self.bit_length() + bits as u64 <= 64 * self.len() as u64,
            "MultiwordInt shl overflow: {}-bit value << {bits} in {} limbs",
            self.bit_length(),
            self.len()
        );
        let mut out = vec![0u64; self.len()];
        shl_into(&mut out, &self.limbs, bits);
        Self { limbs: out }
    }

    pub(crate) fn sub(&self, other: &Self) -> Self {
        let mut out = vec![0u64; self.len()];
        sub_into(&mut out, &self.limbs, &other.limbs);
        Self { limbs: out }
    }

    fn neg(&self) -> Self {
        let mut out = vec![0u64; self.len()];
        neg_into(&mut out, &self.limbs);
        Self { limbs: out }
    }

    fn add(&self, other: &Self) -> Self {
        let mut out = vec![0u64; self.len()];
        add_into(&mut out, &self.limbs, &other.limbs);
        Self { limbs: out }
    }

    fn mul_u32(&self, v: u32) -> Self {
        let mut out = vec![0u64; self.len()];
        mul_u32_into(&mut out, &self.limbs, v);
        Self { limbs: out }
    }

    /// Signed product in `out_len` limbs (truncated). Test/helper wrapper over
    /// [`mul_into`].
    #[cfg(test)]
    pub(crate) fn mul(&self, other: &Self, out_len: usize) -> Self {
        let mut out = vec![0u64; out_len];
        mul_into(&mut out, &self.limbs, &other.limbs);
        Self { limbs: out }
    }

    fn shr_unsigned(&self, bits: u32) -> Self {
        let mut out = vec![0u64; self.len()];
        shr_unsigned_into(&mut out, &self.limbs, bits);
        Self { limbs: out }
    }

    /// Reconstruct a signed value in `len` limbs from mixed-radix Garner digits
    /// and their prime radices via word-level Horner (no `BigInt`), returning the
    /// symmetric representative. `M = ∏ primes` must fit in `64·len` bits and
    /// exceed twice the represented magnitude.
    pub(crate) fn from_garner_digits(digits: &[u32], primes: &[u32], len: usize) -> Self {
        debug_assert_eq!(digits.len(), primes.len());
        // value = Σ_i digit_i · (∏_{j<i} p_j) — Horner over the mixed-radix
        // digits, accumulating the running radix product as the modulus.
        let mut acc = Self::zero(len);
        let mut modulus = {
            let mut l = vec![0u64; len];
            l[0] = 1;
            Self { limbs: l }
        };
        for i in 0..digits.len() {
            acc = acc.add(&modulus.mul_u32(digits[i]));
            modulus = modulus.mul_u32(primes[i]);
        }
        // modulus == M; center to (−M/2, M/2].
        let half = modulus.shr_unsigned(1);
        if ucmp(&acc.limbs, &half.limbs) == Ordering::Greater {
            acc.sub(&modulus)
        } else {
            acc
        }
    }

    /// Const-K convenience wrapper over [`from_garner_digits`](Self::from_garner_digits).
    pub(crate) fn from_rns<const K: usize, P: NttPrimeList<K>>(r: &Rns<K, P>, len: usize) -> Self {
        Self::from_garner_digits(&r.to_garner(), &P::PRIMES, len)
    }

    /// Construct from a [`BigInt`] in `len` limbs (BigInt boundary only). The
    /// magnitude must fit; a too-narrow `len` is caught in debug builds.
    pub(crate) fn from_bigint(x: &BigInt, len: usize) -> Self {
        let mut limbs = vec![0u64; len];
        from_bigint_into(&mut limbs, x);
        Self { limbs }
    }

    /// Reconstruct a [`BigInt`] (BigInt boundary only).
    pub(crate) fn to_bigint(&self) -> BigInt {
        to_bigint(&self.limbs)
    }
}

impl std::ops::SubAssign<&MultiwordInt> for MultiwordInt {
    fn sub_assign(&mut self, rhs: &MultiwordInt) {
        sub_into_self(&mut self.limbs, &rhs.limbs);
    }
}

#[cfg(test)]
mod tests {
    use super::MultiwordInt;
    use crate::rns::{NttPrimes24Bit8, Rns};
    use num::{BigInt, Zero};
    use rand::{rngs::StdRng, RngExt, SeedableRng};

    const LEN: usize = 4; // 256-bit width for the unit tests

    fn rand_value(rng: &mut StdRng) -> (MultiwordInt, BigInt) {
        // Random ~200-bit signed value: hi<<96 + lo (lo in low 96 bits).
        let hi = rng.random::<i64>() as i128;
        let lo = rng.random::<i128>() & ((1i128 << 96) - 1);
        let p = MultiwordInt::from_i128(hi, LEN)
            .shl(96)
            .add(&MultiwordInt::from_i128(lo, LEN));
        let b = p.to_bigint();
        (p, b)
    }

    #[test]
    fn from_to_i128_roundtrip() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..1000 {
            let v = rng.random::<i128>();
            assert_eq!(MultiwordInt::from_i128(v, LEN).to_bigint(), BigInt::from(v));
        }
    }

    #[test]
    fn bigint_roundtrip() {
        let mut rng = StdRng::seed_from_u64(2);
        for _ in 0..1000 {
            let (p, b) = rand_value(&mut rng);
            assert_eq!(MultiwordInt::from_bigint(&b, LEN), p);
            assert_eq!(p.to_bigint(), b);
        }
    }

    #[test]
    fn bit_length_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..1000 {
            let (p, b) = rand_value(&mut rng);
            assert_eq!(p.bit_length(), b.bits(), "value {b}");
        }
        assert_eq!(MultiwordInt::zero(LEN).bit_length(), BigInt::zero().bits());
    }

    #[test]
    fn sub_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(4);
        for _ in 0..1000 {
            let (pa, ba) = rand_value(&mut rng);
            let (pb, bb) = rand_value(&mut rng);
            assert_eq!(pa.sub(&pb).to_bigint(), &ba - &bb);
            let mut acc = pa.clone();
            acc -= &pb;
            assert_eq!(acc.to_bigint(), &ba - &bb);
        }
    }

    #[test]
    fn shl_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(5);
        for _ in 0..500 {
            let lo = rng.random::<i64>() as i128;
            let p = MultiwordInt::from_i128(lo, LEN);
            let b = BigInt::from(lo);
            for &s in &[0u32, 1, 7, 31, 63, 64, 65, 100, 130] {
                assert_eq!(p.shl(s).to_bigint(), &b << s, "v={b} s={s}");
            }
        }
    }

    #[test]
    fn shr_to_i128_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(6);
        for _ in 0..1000 {
            let (p, b) = rand_value(&mut rng);
            for &s in &[0u32, 1, 33, 64, 90, 128] {
                let want = &b >> s;
                if want.bits() <= 126 {
                    assert_eq!(BigInt::from(p.shr_to_i128(s)), want, "v={b} s={s}");
                }
            }
        }
    }

    #[test]
    fn mul_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(8);
        // Wide output so products never truncate.
        const OUT: usize = 8;
        for _ in 0..2000 {
            let (pa, ba) = rand_value(&mut rng);
            let (pb, bb) = rand_value(&mut rng);
            assert_eq!(pa.mul(&pb, OUT).to_bigint(), &ba * &bb);
        }
    }

    #[test]
    fn from_rns_reconstructs_signed() {
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..2000 {
            let v = BigInt::from(rng.random::<i128>()) >> 1; // ≤126 bits, both signs
            let r = Rns::<8, NttPrimes24Bit8>::from_bigint(&v);
            assert_eq!(MultiwordInt::from_rns(&r, LEN).to_bigint(), v);
        }
    }
}
