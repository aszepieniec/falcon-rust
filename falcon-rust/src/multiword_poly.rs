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

    pub(crate) fn wordlen(&self) -> usize {
        self.w
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
            if (i / n_target) % 2 == 0 {
                mw::add_into_self(dst, src);
            } else {
                mw::sub_into_self(dst, src);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::MultiwordPoly;
    use crate::polynomial::Polynomial;
    use num::BigInt;
    use rand::{rngs::StdRng, RngExt, SeedableRng};

    const W: usize = 6; // 384-bit coefficients, room for shifts

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
            let p = Polynomial::new(
                (0..16)
                    .map(|_| BigInt::from(rng.random::<i64>()))
                    .collect(),
            );
            let m = MultiwordPoly::from_bigint_poly(&p, W);
            for &s in &[0u32, 1, 13, 64, 130] {
                let want = Polynomial::new(p.coefficients.iter().map(|c| c << s).collect());
                same(&m.shl_coeffs(s), &want);
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
