//! [`MultiwordPoly`]: a polynomial whose coefficients are signed multi-word
//! integers, stored in one contiguous word-plane buffer so the recursion's
//! ring operations run with no per-coefficient allocation.
//!
//! Coefficient `i` occupies limbs `data[i*w .. (i+1)*w]` (`w` = limbs per
//! coefficient, fixed per depth). This is the allocation-free replacement for
//! `Polynomial<BigInt>` in the deep NTRU-solve recursion: the generic
//! [`Polynomial<F>`](crate::polynomial::Polynomial) allocates a fresh `Vec` (and,
//! for `F = BigInt`, a heap integer per coefficient) in every `Add`/`Sub`/`Mul`,
//! which is the dominant keygen cost.
//!
//! This module holds the *non-multiplying* ring operations (`lift`,
//! `galois_adjoint`, `reduce_by_cyclotomic`, `sub`, per-coefficient shift) plus
//! the `BigInt` boundary conversions. Multiplication (the `field_norm` squarings,
//! the F,G construction product, and the babai `k·f`) goes through the RNS/NTT
//! path; see [`crate::math`] / the RNS multiply added in the next stage.

use num::BigInt;

use crate::multiword_int as mw;
use crate::polynomial::Polynomial;

/// At or below this degree, the flat-word schoolbook negacyclic multiply beats
/// the runtime-`K` RNS/NTT multiply's per-prime overhead. Tuned against the
/// keygen(1024) benchmark with every consumer (field_norm squarings, the
/// asymmetric F,G construction multiplies, and the deep babai `k·f`) routed
/// through [`MultiwordPoly::negacyclic_mul`]: a sweep put the crossover plateau at
/// 64–128 (min ~68 ms), regressing again by 256 as schoolbook's O(n²) loses at
/// large `n`. Shared with the deep babai reduction so it skips building the
/// expensive `RuntimeNtt` in the schoolbook regime.
pub(crate) const SCHOOLBOOK_MAX_N: usize = 64;

/// A degree-`< n` polynomial with signed multi-word coefficients (`w` limbs
/// each), laid out as one contiguous `n·w` buffer.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) struct MultiwordPoly {
    n: usize,
    w: usize,
    data: Vec<u64>,
}

impl MultiwordPoly {
    /// All-zero polynomial with `n` coefficients of `w` limbs each.
    pub(crate) fn zeros(n: usize, w: usize) -> Self {
        Self {
            n,
            w,
            data: vec![0u64; n * w],
        }
    }

    pub(crate) fn n(&self) -> usize {
        self.n
    }

    /// Limb slice of coefficient `i`.
    pub(crate) fn coeff(&self, i: usize) -> &[u64] {
        &self.data[i * self.w..(i + 1) * self.w]
    }

    /// Mutable limb slice of coefficient `i`.
    pub(crate) fn coeff_mut(&mut self, i: usize) -> &mut [u64] {
        &mut self.data[i * self.w..(i + 1) * self.w]
    }

    /// Build from a `BigInt` polynomial, each coefficient widened to `w` limbs.
    /// (`BigInt` boundary — recursion entry.)
    pub(crate) fn from_bigint_poly(p: &Polynomial<BigInt>, w: usize) -> Self {
        let n = p.coefficients.len();
        let mut out = Self::zeros(n, w);
        for (i, c) in p.coefficients.iter().enumerate() {
            mw::from_bigint_into(out.coeff_mut(i), c);
        }
        out
    }

    /// Reconstruct a `BigInt` polynomial. (`BigInt` boundary — recursion exit.)
    pub(crate) fn to_bigint_poly(&self) -> Polynomial<BigInt> {
        Polynomial::new((0..self.n).map(|i| mw::to_bigint(self.coeff(i))).collect())
    }

    /// Multiply each coefficient by `2^bits` (per-coefficient left shift).
    /// Coefficients must stay within `w` limbs (debug-checked via the slice op).
    pub(crate) fn shl_coeffs(&self, bits: u32) -> Self {
        let mut out = Self::zeros(self.n, self.w);
        for i in 0..self.n {
            mw::shl_into(out.coeff_mut(i), self.coeff(i), bits);
        }
        out
    }

    /// `self - other` (same shape). Mirrors `Polynomial - Polynomial`.
    pub(crate) fn sub(&self, other: &Self) -> Self {
        debug_assert!(self.n == other.n && self.w == other.w);
        let mut out = Self::zeros(self.n, self.w);
        for i in 0..self.n {
            mw::sub_into(out.coeff_mut(i), self.coeff(i), other.coeff(i));
        }
        out
    }

    /// In-place `self -= other` (same shape).
    pub(crate) fn sub_assign(&mut self, other: &Self) {
        debug_assert!(self.n == other.n && self.w == other.w);
        for i in 0..self.n {
            let (dst, src) = (
                &mut self.data[i * self.w..(i + 1) * self.w],
                &other.data[i * self.w..(i + 1) * self.w],
            );
            mw::sub_into_self(dst, src);
        }
    }

    /// Bit length of the largest-magnitude coefficient (0 for the zero poly).
    pub(crate) fn max_coeff_bits(&self) -> u64 {
        (0..self.n)
            .map(|i| mw::bit_length(self.coeff(i)))
            .max()
            .unwrap_or(0)
    }

    /// Split into (magnitudes, signs): a same-shape polynomial whose coefficients
    /// are `|c|`, plus the sign of each. Hoists sign handling out of the schoolbook
    /// inner loop so it does only unsigned limb multiplies.
    fn magnitudes(&self) -> (Self, Vec<bool>) {
        let mut mag = self.clone();
        let mut signs = vec![false; self.n];
        for (i, signs_i) in signs.iter_mut().enumerate().take(self.n) {
            let c = mag.coeff_mut(i);
            if mw::is_negative(c) {
                mw::neg_in_place(c);
                *signs_i = true;
            }
        }
        (mag, signs)
    }

    /// `self * other mod (X^n + 1)` by O(n²) schoolbook over flat-word
    /// coefficients, each product accumulated into an `out_w`-limb result. The
    /// allocation-free deep-recursion counterpart to the RNS path: at tiny `n`
    /// (huge coefficients) the RNS multiply's per-prime O(K²) Garner overhead
    /// dominates, and this beats it. Magnitudes/signs are hoisted out of the inner
    /// loop, so the hot path is a single reusable `prod` buffer.
    pub(crate) fn schoolbook_negacyclic_mul(&self, other: &Self, out_w: usize) -> Self {
        debug_assert_eq!(self.n, other.n);
        let n = self.n;
        let (amag, asign) = self.magnitudes();
        let (bmag, bsign) = other.magnitudes();
        let mut out = Self::zeros(n, out_w);
        let mut prod = vec![0u64; out_w];
        for (i, &asign_i) in asign.iter().enumerate().take(n) {
            for (j, &bsign_j) in bsign.iter().enumerate().take(n) {
                mw::umul_into(&mut prod, amag.coeff(i), bmag.coeff(j));
                let k = i + j;
                let (slot, wrap) = if k < n { (k, false) } else { (k - n, true) };
                // Negacyclic wrap flips sign once; operand signs flip it too.
                if asign_i ^ bsign_j ^ wrap {
                    mw::sub_into_self(out.coeff_mut(slot), &prod);
                } else {
                    mw::add_into_self(out.coeff_mut(slot), &prod);
                }
            }
        }
        out
    }

    /// Size-dispatched negacyclic product `self * other mod (X^n + 1)` into
    /// `out_w` limbs. Schoolbook for tiny `n` (deep recursion); runtime-`K` RNS/NTT
    /// otherwise, with the prime count sized from the operands' bit lengths.
    pub(crate) fn negacyclic_mul(&self, other: &Self, out_w: usize) -> Self {
        if self.n <= SCHOOLBOOK_MAX_N {
            self.schoolbook_negacyclic_mul(other, out_w)
        } else {
            // |product coeff| < n · 2^(ba+bb)  ->  ~ ba+bb+log2(n) bits; cover
            // M > 2·that with 24-bit primes (~23 usable bits each), plus a safety
            // prime. Using ba+bb (not 2·max) keeps K tight for asymmetric operands
            // like the F,G construction multiply (large capital × small g_minx).
            let prod_bits =
                self.max_coeff_bits() + other.max_coeff_bits() + self.n.ilog2() as u64 + 2;
            let k = (prod_bits as usize) / 23 + 2;
            let ctx = crate::rns_runtime::RuntimeNtt::cached(self.n, k);
            ctx.negacyclic_mul(self, other, out_w)
        }
    }

    /// Field norm relative to the half-size cyclotomic ring (mirrors
    /// [`Polynomial::field_norm`]): deinterleave into even/odd parts `f0`, `f1`,
    /// then return `f0² - x·f1²` in `Z[X]/(X^{n/2}+1)`. The two squarings go
    /// through [`negacyclic_mul`](Self::negacyclic_mul); the rest is flat-word
    /// shuffles with no per-coefficient allocation.
    pub(crate) fn field_norm(&self) -> Self {
        let n = self.n;
        debug_assert!(n >= 2 && n.is_multiple_of(2));
        let half = n / 2;

        let mut f0 = Self::zeros(half, self.w);
        let mut f1 = Self::zeros(half, self.w);
        for i in 0..half {
            f0.coeff_mut(i).copy_from_slice(self.coeff(2 * i));
            f1.coeff_mut(i).copy_from_slice(self.coeff(2 * i + 1));
        }

        // Squared coefficients: |c|² summed over `half` terms.
        let b = f0.max_coeff_bits().max(f1.max_coeff_bits());
        let prod_bits = 2 * b + (half.max(1)).ilog2() as u64 + 2;
        let out_w = (prod_bits / 64 + 1) as usize;

        let f0_sq = f0.negacyclic_mul(&f0, out_w);
        let f1_sq = f1.negacyclic_mul(&f1, out_w);

        // result = f0² - (x·f1²) reduced by X^half + 1.  x·f1² shifts each
        // coefficient up one slot; the top wraps to slot 0 negated, so
        //   result[0]   = f0²[0]   + f1²[half-1]
        //   result[i>0] = f0²[i]   - f1²[i-1]
        let mut out = f0_sq;
        mw::add_into_self(out.coeff_mut(0), f1_sq.coeff(half - 1));
        for i in 1..half {
            let (dst, src) = (
                &mut out.data[i * out_w..(i + 1) * out_w],
                f1_sq.coeff(i - 1),
            );
            mw::sub_into_self(dst, src);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::Polynomial;
    use num::BigInt;
    use rand::{rngs::StdRng, RngExt, SeedableRng};

    const W: usize = 6; // 384-bit coefficients, room for shifts

    impl MultiwordPoly {
        /// Lift to the next cyclotomic ring by interleaving zeros:
        /// `[c0, c1, …]  ->  [c0, 0, c1, 0, …]` (length doubles). Mirrors
        /// [`Polynomial::lift_next_cyclotomic`].
        pub(crate) fn lift_next_cyclotomic(&self) -> Self {
            let mut out = Self::zeros(self.n * 2, self.w);
            for i in 0..self.n {
                out.coeff_mut(2 * i).copy_from_slice(self.coeff(i));
            }
            out
        }

        /// Galois adjoint: negate the odd-indexed coefficients. Mirrors
        /// [`Polynomial::galois_adjoint`].
        pub(crate) fn galois_adjoint(&self) -> Self {
            let mut out = self.clone();
            for i in (1..self.n).step_by(2) {
                mw::neg_in_place(out.coeff_mut(i));
            }
            out
        }

        /// Reduce by `X^n_target + 1`, folding this polynomial's `self.n`
        /// coefficients down to `n_target`: coefficient block `b = i / n_target`
        /// contributes with sign `(-1)^b` into slot `i % n_target`. Mirrors
        /// [`Polynomial::reduce_by_cyclotomic`].
        pub(crate) fn reduce_by_cyclotomic(&self, n_target: usize) -> Self {
            let mut out = Self::zeros(n_target, self.w);
            for i in 0..self.n {
                let slot = i % n_target;
                let src = self.coeff(i);
                let dst = out.coeff_mut(slot);
                if (i / n_target).is_multiple_of(2) {
                    mw::add_into_self(dst, src);
                } else {
                    mw::sub_into_self(dst, src);
                }
            }
            out
        }
    }

    fn rand_bigint(rng: &mut StdRng) -> BigInt {
        // ~200-bit signed.
        let hi = BigInt::from(rng.random::<i128>());
        let lo = BigInt::from(rng.random::<u64>());
        (hi << 64) + lo
    }

    fn rand_poly(rng: &mut StdRng, n: usize) -> Polynomial<BigInt> {
        Polynomial::new((0..n).map(|_| rand_bigint(rng)).collect())
    }

    fn same(mw: &MultiwordPoly, bi: &Polynomial<BigInt>) {
        assert_eq!(mw.to_bigint_poly().coefficients, bi.coefficients);
    }

    #[test]
    fn roundtrip_bigint_poly() {
        let mut rng = StdRng::seed_from_u64(10);
        for _ in 0..200 {
            let p = rand_poly(&mut rng, 16);
            let m = MultiwordPoly::from_bigint_poly(&p, W);
            same(&m, &p);
        }
    }

    #[test]
    fn lift_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(11);
        for &n in &[2usize, 4, 8, 16] {
            let p = rand_poly(&mut rng, n);
            let m = MultiwordPoly::from_bigint_poly(&p, W);
            same(&m.lift_next_cyclotomic(), &p.lift_next_cyclotomic());
        }
    }

    #[test]
    fn galois_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(12);
        for &n in &[2usize, 4, 8, 16] {
            let p = rand_poly(&mut rng, n);
            let m = MultiwordPoly::from_bigint_poly(&p, W);
            same(&m.galois_adjoint(), &p.galois_adjoint());
        }
    }

    #[test]
    fn sub_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(13);
        for _ in 0..200 {
            let a = rand_poly(&mut rng, 16);
            let b = rand_poly(&mut rng, 16);
            let ma = MultiwordPoly::from_bigint_poly(&a, W);
            let mb = MultiwordPoly::from_bigint_poly(&b, W);
            same(&ma.sub(&mb), &(a.clone() - b.clone()));
            let mut acc = ma.clone();
            acc.sub_assign(&mb);
            same(&acc, &(a - b));
        }
    }

    #[test]
    fn shl_coeffs_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(14);
        for _ in 0..100 {
            // Keep coefficients small so the shift stays within W limbs.
            let p = Polynomial::new((0..16).map(|_| BigInt::from(rng.random::<i64>())).collect());
            let m = MultiwordPoly::from_bigint_poly(&p, W);
            for &s in &[0u32, 1, 13, 64, 130] {
                let want = Polynomial::new(p.coefficients.iter().map(|c| c << s).collect());
                same(&m.shl_coeffs(s), &want);
            }
        }
    }

    /// Schoolbook negacyclic product oracle.
    fn negacyclic_bigint(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>) -> Polynomial<BigInt> {
        let n = a.coefficients.len();
        let mut out = vec![BigInt::from(0); n];
        for i in 0..n {
            for j in 0..n {
                let prod = &a.coefficients[i] * &b.coefficients[j];
                let k = i + j;
                if k < n {
                    out[k] += &prod;
                } else {
                    out[k - n] -= &prod;
                }
            }
        }
        Polynomial::new(out)
    }

    #[test]
    fn schoolbook_negacyclic_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(20);
        const OUT_W: usize = 12;
        for &n in &[1usize, 2, 4, 8] {
            for _ in 0..50 {
                let a = rand_poly(&mut rng, n);
                let b = rand_poly(&mut rng, n);
                let ma = MultiwordPoly::from_bigint_poly(&a, W);
                let mb = MultiwordPoly::from_bigint_poly(&b, W);
                let got = ma.schoolbook_negacyclic_mul(&mb, OUT_W);
                same(&got, &negacyclic_bigint(&a, &b));
            }
        }
    }

    #[test]
    fn negacyclic_mul_dispatch_matches_bigint() {
        // n >= 16 exercises the RNS path; n <= 8 the schoolbook path.
        let mut rng = StdRng::seed_from_u64(21);
        const OUT_W: usize = 12;
        for &n in &[4usize, 8, 16, 32, 64] {
            for _ in 0..10 {
                let a = rand_poly(&mut rng, n);
                let b = rand_poly(&mut rng, n);
                let ma = MultiwordPoly::from_bigint_poly(&a, W);
                let mb = MultiwordPoly::from_bigint_poly(&b, W);
                same(&ma.negacyclic_mul(&mb, OUT_W), &negacyclic_bigint(&a, &b));
            }
        }
    }

    #[test]
    fn field_norm_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(22);
        // Even n; half >= 16 (n >= 32) hits the RNS squaring path.
        for &n in &[2usize, 4, 8, 16, 32, 64] {
            for _ in 0..10 {
                let p = rand_poly(&mut rng, n);
                let mp = MultiwordPoly::from_bigint_poly(&p, W);
                same(&mp.field_norm(), &p.field_norm());
            }
        }
    }

    #[test]
    fn reduce_by_cyclotomic_matches_bigint() {
        let mut rng = StdRng::seed_from_u64(15);
        for &n_target in &[2usize, 4, 8] {
            for &m_len in &[n_target, 2 * n_target - 1, 2 * n_target] {
                let p = rand_poly(&mut rng, m_len);
                let mp = MultiwordPoly::from_bigint_poly(&p, W);
                same(
                    &mp.reduce_by_cyclotomic(n_target),
                    &p.reduce_by_cyclotomic(n_target),
                );
            }
        }
    }
}
