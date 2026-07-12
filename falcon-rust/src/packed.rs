//! [`Packed`]: a signed two's-complement integer stored as `L` little-endian
//! `u64` limbs, holding the value in `(−2^(64L−1), 2^(64L−1))`.
//!
//! It carries the capital coefficients (`F`, `G`) through RNS Babai reduction.
//! That loop straddles two number representations: the `k·f` correction is
//! formed in the RNS/NTT domain (where multiplication is cheap), but estimating
//! the quotient `k` and testing convergence both need the *positional* value of
//! each coefficient. So the capital has to live in a form that is cheap to read
//! positionally and cheap to update from an RNS product, every iteration. A
//! heap big-integer answers neither well — each iteration would pay a full
//! Garner reconstruction and fresh allocations. Fixed-width limbs make all three
//! hot operations allocation-free:
//!
//! * the quotient FFT only needs the top bits of each coefficient, which
//!   [`shr_to_i128`](Packed::shr_to_i128) supplies as a constant number of limb
//!   reads;
//! * convergence is tested against the magnitude, and [`bit_length`](Packed::bit_length)
//!   reads it from the leading limb in O(1);
//! * the RNS product is reassembled straight into the limbs by a word-level CRT
//!   ([`from_rns`](Packed::from_rns)) that writes one fixed-width buffer rather
//!   than allocating per coefficient.
//!
//! `L` is fixed per recursion depth, just large enough to hold the capital
//! without overflow. The capital coefficients roughly triple in bit length each
//! level: at depth 3 they (together with the shifted partial products formed
//! mid-reduction) peak near 165 bits, so `L = 4` (256-bit) suffices; at depth 4
//! they reach ~303 bits, needing `L = 5` (320-bit); deeper levels scale up
//! likewise.
//!
//! The word-level-CRT-into-fixed-limbs layout follows the approach of fn-dsa's
//! `zint31`.

use num::{BigInt, Zero};

use crate::rns::{NttPrimeList, Rns};

/// A signed integer in `(−2^(64L−1), 2^(64L−1))`, two's complement, base `2^64`,
/// little-endian (`limbs[0]` is least significant).  `L ≥ 2` is assumed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Packed<const L: usize> {
    limbs: [u64; L],
}

impl<const L: usize> Packed<L> {
    pub(crate) const ZERO: Self = Self { limbs: [0; L] };

    /// Sign-extend a signed 128-bit value into the `L`-limb representation.
    /// (Test-only constructor; production builds `Packed` via `from_bigint` and
    /// `from_rns`.)
    #[cfg(test)]
    pub(crate) fn from_i128(v: i128) -> Self {
        const { assert!(L >= 2, "Packed requires at least 2 limbs (128-bit input)") };
        let lo = v as u128;
        let ext = if v < 0 { u64::MAX } else { 0 };
        let mut limbs = [ext; L];
        limbs[0] = lo as u64;
        limbs[1] = (lo >> 64) as u64;
        Self { limbs }
    }

    /// True when the value is negative (top bit of the most significant limb).
    fn is_negative(&self) -> bool {
        self.limbs[L - 1] >> 63 == 1
    }

    /// Two's-complement negation (`!x + 1`).
    fn neg(&self) -> Self {
        let mut out = [0u64; L];
        let mut carry = 1u128;
        for (i, out_i) in out.iter_mut().enumerate().take(L) {
            let v = (!self.limbs[i]) as u128 + carry;
            *out_i = v as u64;
            carry = v >> 64;
        }
        Self { limbs: out }
    }

    /// Magnitude (absolute value) limbs.
    fn magnitude(&self) -> [u64; L] {
        if self.is_negative() {
            self.neg().limbs
        } else {
            self.limbs
        }
    }

    /// Number of bits in the absolute value, matching [`BigInt::bits`]:
    /// zero for the value 0, otherwise `floor(log2(|v|)) + 1`.
    pub(crate) fn bit_length(&self) -> u64 {
        let m = self.magnitude();
        for i in (0..L).rev() {
            if m[i] != 0 {
                return (i as u64) * 64 + (64 - m[i].leading_zeros() as u64);
            }
        }
        0
    }

    /// Arithmetic right shift by `shift` bits (floor division by `2^shift`),
    /// returning the low 128 bits of the result as `i128`.
    ///
    /// Used to feed the `f64` FFT: callers pick `shift` so the result lands in
    /// `i128` range, bit-for-bit identical to the `BigInt` path's
    /// `bi >> capital_shift` (a true arithmetic shift, flooring toward −∞).
    pub(crate) fn shr_to_i128(&self, shift: u32) -> i128 {
        // Shift the two's-complement representation directly, sign-extending the
        // vacated high bits.  This is an arithmetic shift, so it floors toward
        // −∞ for negative values — matching `BigInt::shr` exactly, including the
        // low retained bit (the previous magnitude-based version truncated
        // toward zero and could differ by 1).
        let sign = if self.is_negative() { u64::MAX } else { 0 };
        let limb = (shift / 64) as usize;
        let bit = shift % 64;
        let get = |idx: usize| -> u64 {
            if idx < L {
                self.limbs[idx]
            } else {
                sign
            }
        };
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

    /// Logical left shift by `bits` (multiply by `2^bits`); bits beyond the
    /// `64·L`-bit width are dropped.  The debug assertion guards callers against
    /// a too-small `L`: the shifted magnitude must stay within `64·L` bits.
    pub(crate) fn shl(&self, bits: u32) -> Self {
        debug_assert!(
            self.bit_length() + bits as u64 <= 64 * L as u64,
            "Packed<{L}> shl overflow: {}-bit value << {bits}",
            self.bit_length()
        );
        let limb = (bits / 64) as usize;
        let bit = bits % 64;
        let mut out = [0u64; L];
        for (i, out_i) in out.iter_mut().enumerate().take(L) {
            let src = i as isize - limb as isize;
            if src < 0 {
                continue;
            }
            let src = src as usize;
            let mut word = self.limbs[src] as u128;
            if bit != 0 {
                let lower = if src >= 1 {
                    self.limbs[src - 1] as u128
                } else {
                    0
                };
                word = (word << bit) | (lower >> (64 - bit));
            }
            *out_i = word as u64;
        }
        Self { limbs: out }
    }

    /// `self - other`, wrapping in two's complement.
    pub(crate) fn sub(&self, other: &Self) -> Self {
        let mut out = [0u64; L];
        let mut borrow = 0i128;
        for (i, out_i) in out.iter_mut().enumerate().take(L) {
            let v = self.limbs[i] as i128 - other.limbs[i] as i128 - borrow;
            *out_i = v as u64;
            borrow = if v < 0 { 1 } else { 0 };
        }
        Self { limbs: out }
    }

    fn sub_assign(&mut self, other: &Self) {
        *self = self.sub(other);
    }

    /// `self * (u32)`, truncated to `64·L` bits.
    fn mul_u32(&self, v: u32) -> Self {
        let mut out = [0u64; L];
        let mut carry = 0u128;
        for (i, out_i) in out.iter_mut().enumerate().take(L) {
            let prod = self.limbs[i] as u128 * v as u128 + carry;
            *out_i = prod as u64;
            carry = prod >> 64;
        }
        Self { limbs: out }
    }

    /// Unsigned compare (treats both as non-negative `64·L`-bit magnitudes).
    fn ucmp(&self, other: &Self) -> core::cmp::Ordering {
        for i in (0..L).rev() {
            match self.limbs[i].cmp(&other.limbs[i]) {
                core::cmp::Ordering::Equal => continue,
                ord => return ord,
            }
        }
        core::cmp::Ordering::Equal
    }

    /// Reconstruct a signed value from its RNS residues via **word-level**
    /// Garner CRT (no `BigInt`), returning the symmetric representative.
    ///
    /// The modulus `M = ∏ primes` must fit in `64·L` bits and exceed twice the
    /// represented magnitude (callers size `L` and the prime set accordingly).
    pub(crate) fn from_rns<const K: usize, P: NttPrimeList<K>>(r: &Rns<K, P>) -> Self {
        let digits = r.to_garner();
        // value = ((…((d_{K-1})·p_{K-2} + d_{K-2})…)·p_0 + d_0) — Horner over
        // the mixed-radix digits, most significant first.
        let mut acc = Self::ZERO;
        let mut modulus = {
            let mut l = [0u64; L];
            l[0] = 1;
            Self { limbs: l }
        };
        for (i, &digits_i) in digits.iter().enumerate().take(K) {
            // acc += digit_i · (∏_{j<i} p_j)
            acc = acc.add_packed(&modulus.mul_u32(digits_i));
            modulus = modulus.mul_u32(P::PRIMES[i]);
        }
        // modulus == M; center to (−M/2, M/2].
        let half = modulus.shr_unsigned(1);
        if acc.ucmp(&half) == core::cmp::Ordering::Greater {
            acc.sub(&modulus)
        } else {
            acc
        }
    }

    /// `self + other` (full-width, wrapping).  Helper for CRT accumulation.
    fn add_packed(&self, other: &Self) -> Self {
        let mut out = [0u64; L];
        let mut carry = 0u128;
        for (i, out_i) in out.iter_mut().enumerate().take(L) {
            let s = self.limbs[i] as u128 + other.limbs[i] as u128 + carry;
            *out_i = s as u64;
            carry = s >> 64;
        }
        Self { limbs: out }
    }

    /// Logical (unsigned) right shift by `bits`, used only on the positive
    /// modulus to form `M/2`.
    fn shr_unsigned(&self, bits: u32) -> Self {
        let limb = (bits / 64) as usize;
        let bit = bits % 64;
        let mut out = [0u64; L];
        for (i, out_i) in out.iter_mut().enumerate().take(L) {
            let src = i + limb;
            let mut word = if src < L { self.limbs[src] as u128 } else { 0 };
            if bit != 0 {
                let hi = if src + 1 < L {
                    self.limbs[src + 1] as u128
                } else {
                    0
                };
                word = (word >> bit) | (hi << (64 - bit));
                word &= u64::MAX as u128;
            }
            *out_i = word as u64;
        }
        Self { limbs: out }
    }

    /// Construct from a [`BigInt`] (used at the reduction entry point only).
    ///
    /// The magnitude must fit in `L` limbs; an oversized value would have its
    /// high limbs silently dropped, so a too-small `L` for the depth's capital
    /// is caught here in debug builds rather than corrupting the reduction.
    pub(crate) fn from_bigint(x: &BigInt) -> Self {
        let neg = x.sign() == num::bigint::Sign::Minus;
        let mag = if neg { -x } else { x.clone() };
        let (_, words) = mag.to_u64_digits();
        debug_assert!(
            words.len() <= L,
            "Packed<{L}> too narrow: value needs {} limbs ({} bits)",
            words.len(),
            x.bits()
        );
        let mut limbs = [0u64; L];
        for (i, w) in words.iter().take(L).enumerate() {
            limbs[i] = *w;
        }
        let p = Self { limbs };
        if neg {
            p.neg()
        } else {
            p
        }
    }

    /// Reconstruct a [`BigInt`] (used at the reduction exit point only).
    pub(crate) fn to_bigint(self) -> BigInt {
        let neg = self.is_negative();
        let m = self.magnitude();
        let mut acc = BigInt::zero();
        for i in (0..L).rev() {
            acc <<= 64;
            acc += BigInt::from(m[i]);
        }
        if neg {
            -acc
        } else {
            acc
        }
    }
}

impl<const L: usize> std::ops::SubAssign<&Packed<L>> for Packed<L> {
    fn sub_assign(&mut self, rhs: &Packed<L>) {
        Packed::sub_assign(self, rhs)
    }
}

#[cfg(test)]
mod tests {
    use crate::rns::{NttPrimes24Bit8, Rns};
    use num::{BigInt, Zero};
    use rand::{rngs::StdRng, RngExt, SeedableRng};

    // The unit tests exercise the 4-limb (256-bit) width.
    type Packed = super::Packed<4>;

    fn rand_packed(rng: &mut StdRng) -> (Packed, BigInt) {
        // Random ~200-bit signed value built from i128 + shifted i128.
        let hi = rng.random::<i64>() as i128; // ~64 bits signed (top)
        let lo = rng.random::<i128>();
        let p = Packed::from_i128(hi)
            .shl(96)
            .add_packed(&Packed::from_i128(lo & ((1i128 << 96) - 1)));
        let b = p.to_bigint();
        (p, b)
    }

    #[test]
    fn from_to_i128_roundtrip() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..1000 {
            let v = rng.random::<i128>();
            assert_eq!(Packed::from_i128(v).to_bigint(), BigInt::from(v));
        }
    }

    #[test]
    fn bigint_roundtrip() {
        let mut rng = StdRng::seed_from_u64(2);
        for _ in 0..1000 {
            let (p, b) = rand_packed(&mut rng);
            assert_eq!(Packed::from_bigint(&b), p);
            assert_eq!(p.to_bigint(), b);
        }
    }

    #[test]
    fn bit_length_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..1000 {
            let (p, b) = rand_packed(&mut rng);
            assert_eq!(p.bit_length(), b.bits(), "value {b}");
        }
        assert_eq!(Packed::ZERO.bit_length(), BigInt::zero().bits());
    }

    #[test]
    fn sub_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(4);
        for _ in 0..1000 {
            let (pa, ba) = rand_packed(&mut rng);
            let (pb, bb) = rand_packed(&mut rng);
            assert_eq!(pa.sub(&pb).to_bigint(), &ba - &bb);
        }
    }

    #[test]
    fn shl_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(5);
        for _ in 0..500 {
            // Keep magnitude small enough that the shift stays in range.
            let lo = rng.random::<i64>() as i128;
            let p = Packed::from_i128(lo);
            let b = BigInt::from(lo);
            for &s in &[0u32, 1, 7, 31, 63, 64, 65, 100, 130] {
                assert_eq!(p.shl(s).to_bigint(), &b << s, "v={b} s={s}");
            }
        }
    }

    #[test]
    fn from_rns_reconstructs_signed() {
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..2000 {
            // Values well within the 8-prime (~183-bit) signed capacity: build
            // a ~150-bit signed integer (i128 shifted) and round-trip it.
            let v = BigInt::from(rng.random::<i128>()) >> 1; // ≤126 bits, both signs
            let r = Rns::<8, NttPrimes24Bit8>::from_bigint(&v);
            assert_eq!(Packed::from_rns(&r).to_bigint(), v);
        }
    }

    #[test]
    fn shr_to_i128_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(6);
        for _ in 0..1000 {
            let (p, b) = rand_packed(&mut rng);
            for &s in &[0u32, 1, 33, 64, 90, 128] {
                // shr_to_i128 is a true arithmetic shift (floor toward −∞),
                // bit-for-bit equal to BigInt `>>` whenever the result fits i128.
                let want = &b >> s;
                let got = BigInt::from(p.shr_to_i128(s));
                if want.bits() <= 126 {
                    assert_eq!(got, want, "v={b} s={s}");
                }
            }
        }
    }
}
