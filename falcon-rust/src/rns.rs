use std::fmt;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{One, Zero};

use crate::fp_field::{montyred, negqinv_modr, r_sq_modq, FpField};

/// A set of odd prime moduli known at compile time.
pub(crate) trait PrimeList<const N: usize> {
    const PRIMES: [u32; N];
}

/// Dispatch a Montgomery multiplication to the right `montyred<LOG2R>`
/// specialisation based on a runtime LOG2R value.  Each match arm uses a
/// literal const generic so the dead branch inside `montyred` is still
/// eliminated at compile time for that arm.
macro_rules! dispatch_log2r {
    ($log2r:expr, $a:expr, $b:expr, $p:expr, $neg_inv:expr) => {{
        let a = $a;
        let b = $b;
        let p = $p;
        let neg_inv = $neg_inv;
        match $log2r {
            1  => montyred::<1> (a, b, p, neg_inv),
            2  => montyred::<2> (a, b, p, neg_inv),
            3  => montyred::<3> (a, b, p, neg_inv),
            4  => montyred::<4> (a, b, p, neg_inv),
            5  => montyred::<5> (a, b, p, neg_inv),
            6  => montyred::<6> (a, b, p, neg_inv),
            7  => montyred::<7> (a, b, p, neg_inv),
            8  => montyred::<8> (a, b, p, neg_inv),
            9  => montyred::<9> (a, b, p, neg_inv),
            10 => montyred::<10>(a, b, p, neg_inv),
            11 => montyred::<11>(a, b, p, neg_inv),
            12 => montyred::<12>(a, b, p, neg_inv),
            13 => montyred::<13>(a, b, p, neg_inv),
            14 => montyred::<14>(a, b, p, neg_inv),
            15 => montyred::<15>(a, b, p, neg_inv),
            16 => montyred::<16>(a, b, p, neg_inv),
            17 => montyred::<17>(a, b, p, neg_inv),
            18 => montyred::<18>(a, b, p, neg_inv),
            19 => montyred::<19>(a, b, p, neg_inv),
            20 => montyred::<20>(a, b, p, neg_inv),
            21 => montyred::<21>(a, b, p, neg_inv),
            22 => montyred::<22>(a, b, p, neg_inv),
            23 => montyred::<23>(a, b, p, neg_inv),
            24 => montyred::<24>(a, b, p, neg_inv),
            25 => montyred::<25>(a, b, p, neg_inv),
            26 => montyred::<26>(a, b, p, neg_inv),
            27 => montyred::<27>(a, b, p, neg_inv),
            28 => montyred::<28>(a, b, p, neg_inv),
            29 => montyred::<29>(a, b, p, neg_inv),
            30 => montyred::<30>(a, b, p, neg_inv),
            31 => montyred::<31>(a, b, p, neg_inv),
            32 => montyred::<32>(a, b, p, neg_inv),
            _  => unreachable!(),
        }
    }};
}

/// An integer represented as its residues modulo a list of primes.
///
/// Each residue is stored in Montgomery form: `residues[i] = value · R_i mod
/// PRIMES[i]` where `R_i = 2^(PRIMES[i].ilog2()+1)`.  This matches the
/// convention used by `FpField<Q>`.
///
/// Arithmetic (add, sub, mul) is component-wise and stays in Montgomery form
/// throughout, so no extra conversions are needed between operations.
pub(crate) struct Rns<const N: usize, P: PrimeList<N>> {
    residues: [u32; N],
    _phantom: PhantomData<P>,
}

// Manual trait impls so that P need not satisfy Clone/Copy/Debug/PartialEq/Eq.
// PhantomData<P> is always Clone+Copy regardless of P, and residues is [u32;N]
// which is always Clone+Copy+PartialEq+Eq.

impl<const N: usize, P: PrimeList<N>> Clone for Rns<N, P> {
    fn clone(&self) -> Self {
        Rns { residues: self.residues, _phantom: PhantomData }
    }
}

impl<const N: usize, P: PrimeList<N>> Copy for Rns<N, P> {}

impl<const N: usize, P: PrimeList<N>> PartialEq for Rns<N, P> {
    fn eq(&self, other: &Self) -> bool {
        self.residues == other.residues
    }
}

impl<const N: usize, P: PrimeList<N>> Eq for Rns<N, P> {}

impl<const N: usize, P: PrimeList<N>> fmt::Debug for Rns<N, P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Rns").field("residues", &self.residues).finish()
    }
}

impl<const N: usize, P: PrimeList<N>> Rns<N, P> {
    /// `R_i = 2^LOG2R[i]` is the Montgomery base for `PRIMES[i]`.
    const LOG2R: [u32; N] = {
        let primes = P::PRIMES;
        let mut a = [0u32; N];
        let mut i = 0;
        while i < N {
            a[i] = primes[i].ilog2() + 1;
            i += 1;
        }
        a
    };

    /// `-PRIMES[i]^{-1} mod R_i`, the Montgomery reduction constant.
    const NEG_INV: [u32; N] = {
        let primes = P::PRIMES;
        let mut a = [0u32; N];
        let mut i = 0;
        while i < N {
            a[i] = negqinv_modr(primes[i]);
            i += 1;
        }
        a
    };

    /// `R_i^2 mod PRIMES[i]`, used to convert canonical values to Montgomery form.
    const R_SQ: [u32; N] = {
        let primes = P::PRIMES;
        let mut a = [0u32; N];
        let mut i = 0;
        while i < N {
            a[i] = r_sq_modq(primes[i]);
            i += 1;
        }
        a
    };

    /// Convert canonical value `v` (already reduced mod `PRIMES[i]`) to
    /// Montgomery form for prime index `i`.
    fn to_mont_at(v: u32, i: usize) -> u32 {
        dispatch_log2r!(Self::LOG2R[i], v, Self::R_SQ[i], P::PRIMES[i], Self::NEG_INV[i])
    }

    /// Convert Montgomery-form value `a` back to canonical form for prime index `i`.
    fn from_mont_at(a: u32, i: usize) -> u32 {
        dispatch_log2r!(Self::LOG2R[i], a, 1u32, P::PRIMES[i], Self::NEG_INV[i])
    }

    /// Construct from a `u32`, reducing modulo each prime.
    pub(crate) fn from_u32(v: u32) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            residues[i] = Self::to_mont_at(v % P::PRIMES[i], i);
        }
        Rns { residues, _phantom: PhantomData }
    }

    /// Construct from a `u64`, reducing modulo each prime.
    pub(crate) fn from_u64(v: u64) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            residues[i] = Self::to_mont_at((v % P::PRIMES[i] as u64) as u32, i);
        }
        Rns { residues, _phantom: PhantomData }
    }

    /// Return the canonical representative of the `i`-th residue in `[0, PRIMES[i])`.
    pub(crate) fn residue(&self, i: usize) -> u32 {
        Self::from_mont_at(self.residues[i], i)
    }

    /// Return the Garner mixed-radix coefficients `[a_0, …, a_{N-1}]` where
    /// `a_i < PRIMES[i]` and `x = a_0 + PRIMES[0]·(a_1 + PRIMES[1]·(…))`.
    ///
    /// This is the canonical lossless output: the full RNS modulus M = ∏ PRIMES[i]
    /// does not fit in any primitive integer type, so no primitive reconstruction
    /// is provided here.
    pub(crate) fn to_garner(&self) -> [u32; N] {
        // Start with canonical residues.
        let mut u = [0u32; N];
        for i in 0..N {
            u[i] = self.residue(i);
        }
        // Garner's algorithm: for each i, eliminate u[i] from all later entries.
        for i in 0..N {
            for j in i + 1..N {
                let pi = P::PRIMES[i] as u64;
                let pj = P::PRIMES[j] as u64;
                // inv = PRIMES[i]^{-1} mod PRIMES[j] via Fermat's little theorem
                let inv = mod_pow(pi % pj, pj - 2, pj);
                let ui_mod_pj = u[i] as u64 % pj;
                let diff = (u[j] as u64 + pj - ui_mod_pj) % pj;
                u[j] = (diff * inv % pj) as u32;
            }
        }
        u
    }

    /// Reconstruct as a `u64` via Garner's algorithm.
    ///
    /// Only correct when the true value is less than 2^64; wraps silently
    /// otherwise.
    pub(crate) fn to_u64(&self) -> u64 {
        let a = self.to_garner();
        let mut result = 0u64;
        let mut base = 1u64;
        for i in 0..N {
            result = result.wrapping_add(base.wrapping_mul(a[i] as u64));
            base = base.wrapping_mul(P::PRIMES[i] as u64);
        }
        result
    }

    /// Construct from a signed `i32`, reducing modulo each prime with correct
    /// sign handling.
    pub(crate) fn from_i32(x: i32) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            let p = P::PRIMES[i];
            let r = x.rem_euclid(p as i32) as u32;
            residues[i] = Self::to_mont_at(r, i);
        }
        Rns { residues, _phantom: PhantomData }
    }

    /// Construct from a signed `i128`, reducing modulo each prime with correct
    /// sign handling.
    pub(crate) fn from_i128(x: i128) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            let r = x.rem_euclid(P::PRIMES[i] as i128) as u32;
            residues[i] = Self::to_mont_at(r, i);
        }
        Rns { residues, _phantom: PhantomData }
    }

    /// Reconstruct as a signed `i64` via Garner's algorithm.
    ///
    /// Uses the symmetric representative: values above M/2 are interpreted as
    /// negative (analogous to two's complement but with modulus M).  Correct
    /// when the true value fits in i64.  For K ≥ 3 primes of ≥ 30 bits the
    /// modulus M exceeds 2^64 and values outside i64 range wrap silently.
    pub(crate) fn to_i64(&self) -> i64 {
        let a = self.to_garner();
        let mut result = 0u128;
        let mut base = 1u128;
        for i in 0..N {
            result = result.wrapping_add(base.wrapping_mul(a[i] as u128));
            base = base.wrapping_mul(P::PRIMES[i] as u128);
        }
        // base == M; symmetric range: values in (M/2, M) are negative.
        if result > base / 2 {
            result.wrapping_sub(base) as i64
        } else {
            result as i64
        }
    }

    /// Reconstruct as a signed `i128` via Garner's algorithm.
    ///
    /// Correct when the true value fits in i128.  For K ≤ 5 primes of ≤ 24 bits
    /// the modulus M < 2^120 and the u128 accumulator never overflows.
    pub(crate) fn to_i128(&self) -> i128 {
        let a = self.to_garner();
        let mut result = 0u128;
        let mut base = 1u128;
        for i in 0..N {
            result = result.wrapping_add(base.wrapping_mul(a[i] as u128));
            base = base.wrapping_mul(P::PRIMES[i] as u128);
        }
        if result > base / 2 {
            result.wrapping_sub(base) as i128
        } else {
            result as i128
        }
    }
}

/// Binary exponentiation: base^exp mod modulus.  All three values are u64;
/// intermediate products use u128 to avoid overflow.
fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result as u128 * base as u128 % modulus as u128) as u64;
        }
        base = (base as u128 * base as u128 % modulus as u128) as u64;
        exp >>= 1;
    }
    result
}

impl<const N: usize, P: PrimeList<N>> Add for Rns<N, P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            let p = P::PRIMES[i] as u64;
            let s = self.residues[i] as u64 + rhs.residues[i] as u64;
            residues[i] = (if s >= p { s - p } else { s }) as u32;
        }
        Rns { residues, _phantom: PhantomData }
    }
}

impl<const N: usize, P: PrimeList<N>> AddAssign for Rns<N, P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const N: usize, P: PrimeList<N>> Neg for Rns<N, P> {
    type Output = Self;

    fn neg(self) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            let r = self.residues[i];
            residues[i] = if r == 0 { 0 } else { P::PRIMES[i] - r };
        }
        Rns { residues, _phantom: PhantomData }
    }
}

impl<const N: usize, P: PrimeList<N>> Sub for Rns<N, P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}

impl<const N: usize, P: PrimeList<N>> SubAssign for Rns<N, P> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const N: usize, P: PrimeList<N>> Mul for Rns<N, P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut residues = [0u32; N];
        for i in 0..N {
            residues[i] = dispatch_log2r!(
                Self::LOG2R[i],
                self.residues[i],
                rhs.residues[i],
                P::PRIMES[i],
                Self::NEG_INV[i]
            );
        }
        Rns { residues, _phantom: PhantomData }
    }
}

impl<const N: usize, P: PrimeList<N>> MulAssign for Rns<N, P> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const N: usize, P: PrimeList<N>> Zero for Rns<N, P> {
    fn zero() -> Self {
        Rns { residues: [0u32; N], _phantom: PhantomData }
    }

    fn is_zero(&self) -> bool {
        self.residues.iter().all(|&r| r == 0)
    }
}

impl<const N: usize, P: PrimeList<N>> One for Rns<N, P> {
    fn one() -> Self {
        Self::from_u32(1)
    }
}

// ---------------------------------------------------------------------------
// NTT support
// ---------------------------------------------------------------------------

/// A prime list whose primes admit NTT-based polynomial multiplication mod
/// X^n+1 for n ≤ 1024.
///
/// Every prime must satisfy `p ≡ 1 (mod 2048)`, and
/// `ROOTS_OF_UNITY_2048[i]` must be a primitive 2048th root of unity for
/// `PRIMES[i]` in canonical (non-Montgomery) form.
pub(crate) trait NttPrimeList<const N: usize>: PrimeList<N> {
    const ROOTS_OF_UNITY_2048: [u32; N];
}

/// Given a primitive 2048th root of unity (canonical) for prime `p`, return
/// the primitive 2n-th root needed for an NTT of length `n` (n ≤ 1024, power
/// of two) by squaring `log2(1024/n)` times.
fn primitive_root_2n(root_2048: u32, n: usize, p: u32) -> u32 {
    debug_assert!(n >= 1 && n <= 1024 && n.is_power_of_two());
    let squarings = 10u32 - n.ilog2(); // 2^10 = 1024
    let mut root = root_2048 as u64;
    for _ in 0..squarings {
        root = root * root % p as u64;
    }
    root as u32
}

/// Reverse the `bits` least-significant bits of `x`.
fn bitrev(mut x: usize, bits: usize) -> usize {
    let mut r = 0usize;
    for _ in 0..bits {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    r
}

/// Compute bit-reversed powers `[ω^0, ω^1, …, ω^{n-1}]` with `ω` given in
/// Montgomery form, placing them in bit-reversed index order.  The returned
/// slice has length `n`.
fn bitrev_powers_mont(
    root_mont: u32,
    n: usize,
    p: u32,
    log2r: u32,
    neg_inv: u32,
    r_sq: u32,
) -> Vec<u32> {
    // one_mont = R mod p (Montgomery representation of 1)
    let one_mont = dispatch_log2r!(log2r, 1u32, r_sq, p, neg_inv);
    let mut array = vec![0u32; n];
    let mut alpha = one_mont;
    for a in array.iter_mut() {
        *a = alpha;
        alpha = dispatch_log2r!(log2r, alpha, root_mont, p, neg_inv);
    }
    let log2n = n.ilog2() as usize;
    for i in 0..n {
        let j = bitrev(i, log2n);
        if i < j {
            array.swap(i, j);
        }
    }
    array
}

/// Cooley-Tukey NTT butterfly over Z/pZ with Montgomery arithmetic.  All
/// values in `a` and `psi_rev` are in Montgomery form.
fn ntt_u32(a: &mut [u32], psi_rev: &[u32], p: u32, log2r: u32, neg_inv: u32) {
    let n = a.len();
    let mut t = n;
    let mut m = 1;
    while m < n {
        t >>= 1;
        for i in 0..m {
            let j1 = 2 * i * t;
            let s = psi_rev[m + i];
            for j in j1..j1 + t {
                let u = a[j];
                let v = dispatch_log2r!(log2r, a[j + t], s, p, neg_inv);
                let sum = u as u64 + v as u64;
                a[j] = if sum >= p as u64 { (sum - p as u64) as u32 } else { sum as u32 };
                a[j + t] = if u >= v { u - v } else { u + p - v };
            }
        }
        m <<= 1;
    }
}

/// Gentleman-Sande INTT butterfly over Z/pZ with Montgomery arithmetic.
/// Includes the final scaling by `ninv_mont` (= n⁻¹ in Montgomery form).
fn intt_u32(
    a: &mut [u32],
    psi_inv_rev: &[u32],
    ninv_mont: u32,
    p: u32,
    log2r: u32,
    neg_inv: u32,
) {
    let n = a.len();
    let mut t = 1;
    let mut m = n;
    while m > 1 {
        let h = m / 2;
        let mut j1 = 0;
        for i in 0..h {
            let s = psi_inv_rev[h + i];
            for j in j1..j1 + t {
                let u = a[j];
                let v = a[j + t];
                let sum = u as u64 + v as u64;
                a[j] = if sum >= p as u64 { (sum - p as u64) as u32 } else { sum as u32 };
                let sub = if u >= v { u - v } else { u + p - v };
                a[j + t] = dispatch_log2r!(log2r, sub, s, p, neg_inv);
            }
            j1 += 2 * t;
        }
        t <<= 1;
        m >>= 1;
    }
    for ai in a.iter_mut() {
        *ai = dispatch_log2r!(log2r, *ai, ninv_mont, p, neg_inv);
    }
}

/// In-place negacyclic NTT of a coefficient-ordered slice of RNS elements.
///
/// The length must be a power of two ≤ 1024.  Each prime in `P` must satisfy
/// `p ≡ 1 (mod 2n)`.  After this call the slice is in NTT evaluation domain:
/// element-wise multiplication of two such slices corresponds to polynomial
/// multiplication mod X^n + 1 in RNS.
pub(crate) fn ntt_inplace<const N: usize, P: NttPrimeList<N>>(coeffs: &mut [Rns<N, P>]) {
    let n = coeffs.len();
    debug_assert!(n >= 1 && n <= 1024 && n.is_power_of_two());
    for prime_idx in 0..N {
        let p = P::PRIMES[prime_idx];
        let log2r = p.ilog2() + 1;
        let neg_inv = negqinv_modr(p);
        let r_sq = r_sq_modq(p);

        let root_2n = primitive_root_2n(P::ROOTS_OF_UNITY_2048[prime_idx], n, p);
        let root_2n_mont = dispatch_log2r!(log2r, root_2n, r_sq, p, neg_inv);
        let psi_rev = bitrev_powers_mont(root_2n_mont, n, p, log2r, neg_inv, r_sq);

        let mut slice: Vec<u32> = coeffs.iter().map(|r| r.residues[prime_idx]).collect();
        ntt_u32(&mut slice, &psi_rev, p, log2r, neg_inv);
        for (coeff, val) in coeffs.iter_mut().zip(slice) {
            coeff.residues[prime_idx] = val;
        }
    }
}

/// In-place negacyclic INTT: inverse of [`ntt_inplace`].
pub(crate) fn intt_inplace<const N: usize, P: NttPrimeList<N>>(coeffs: &mut [Rns<N, P>]) {
    let n = coeffs.len();
    debug_assert!(n >= 1 && n <= 1024 && n.is_power_of_two());
    for prime_idx in 0..N {
        let p = P::PRIMES[prime_idx];
        let log2r = p.ilog2() + 1;
        let neg_inv = negqinv_modr(p);
        let r_sq = r_sq_modq(p);

        let root_2n = primitive_root_2n(P::ROOTS_OF_UNITY_2048[prime_idx], n, p);
        let root_inv = mod_pow(root_2n as u64, (p - 2) as u64, p as u64) as u32;
        let root_inv_mont = dispatch_log2r!(log2r, root_inv, r_sq, p, neg_inv);
        let psi_inv_rev = bitrev_powers_mont(root_inv_mont, n, p, log2r, neg_inv, r_sq);

        let n_inv = mod_pow(n as u64, (p - 2) as u64, p as u64) as u32;
        let ninv_mont = dispatch_log2r!(log2r, n_inv, r_sq, p, neg_inv);

        let mut slice: Vec<u32> = coeffs.iter().map(|r| r.residues[prime_idx]).collect();
        intt_u32(&mut slice, &psi_inv_rev, ninv_mont, p, log2r, neg_inv);
        for (coeff, val) in coeffs.iter_mut().zip(slice) {
            coeff.residues[prime_idx] = val;
        }
    }
}

/// Two 24-bit NTT-friendly primes (p ≡ 1 mod 2048, 2²³ ≤ p < 2²⁴).
/// Signed capacity ≈ 45 bits; covers `babai_reduce_rns` at recursion depth 1.
pub(crate) struct NttPrimes24Bit2;
impl PrimeList<2> for NttPrimes24Bit2 {
    const PRIMES: [u32; 2] = [8_404_993, 8_427_521];
}
impl NttPrimeList<2> for NttPrimes24Bit2 {
    const ROOTS_OF_UNITY_2048: [u32; 2] = [
        FpField::<8_404_993>::primitive_nth_root_of_unity(2048).value(),
        FpField::<8_427_521>::primitive_nth_root_of_unity(2048).value(),
    ];
}

/// Four 24-bit NTT-friendly primes (p ≡ 1 mod 2048, 2²³ ≤ p < 2²⁴).
/// Signed capacity ≈ 91 bits; covers `babai_reduce_rns` at recursion depth 2.
pub(crate) struct NttPrimes24Bit4;
impl PrimeList<4> for NttPrimes24Bit4 {
    const PRIMES: [u32; 4] = [8_404_993, 8_427_521, 8_441_857, 8_452_097];
}
impl NttPrimeList<4> for NttPrimes24Bit4 {
    const ROOTS_OF_UNITY_2048: [u32; 4] = [
        FpField::<8_404_993>::primitive_nth_root_of_unity(2048).value(),
        FpField::<8_427_521>::primitive_nth_root_of_unity(2048).value(),
        FpField::<8_441_857>::primitive_nth_root_of_unity(2048).value(),
        FpField::<8_452_097>::primitive_nth_root_of_unity(2048).value(),
    ];
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::fp_field::FpField;

    struct SmallPrimes;
    impl PrimeList<3> for SmallPrimes {
        // Three NTT-friendly primes satisfying p ≡ 1 (mod 2048); LOG2R = 20, 30, 31.
        const PRIMES: [u32; 3] = [786_433, 998_244_353, 1_073_754_113];
    }
    impl NttPrimeList<3> for SmallPrimes {
        const ROOTS_OF_UNITY_2048: [u32; 3] = [
            FpField::<786_433>::primitive_nth_root_of_unity(2048).value(),
            FpField::<998_244_353>::primitive_nth_root_of_unity(2048).value(),
            FpField::<1_073_754_113>::primitive_nth_root_of_unity(2048).value(),
        ];
    }
    type Rns3 = Rns<3, SmallPrimes>;

    // A two-prime list that exercises the u128 branch (prime > 2^31).
    struct LargePrimes;
    impl PrimeList<2> for LargePrimes {
        const PRIMES: [u32; 2] = [1_073_754_113, 4_294_967_291];
    }
    type Rns2L = Rns<2, LargePrimes>;

    #[test]
    fn residues_equal_naive_reduction() {
        let primes = SmallPrimes::PRIMES;
        for v in [0u32, 1, 12345, 786_432, 998_244_352, 1_073_754_112] {
            let r = Rns3::from_u32(v);
            for i in 0..3 {
                assert_eq!(r.residue(i), v % primes[i], "v={v}, i={i}");
            }
        }
    }

    #[test]
    fn from_u64_residues_correct() {
        let primes = SmallPrimes::PRIMES;
        let v: u64 = 0x_DEAD_BEEF_1234_5678;
        let r = Rns3::from_u64(v);
        for i in 0..3 {
            assert_eq!(r.residue(i), (v % primes[i] as u64) as u32, "i={i}");
        }
    }

    #[test]
    fn to_u64_roundtrip() {
        for v in [0u64, 1, 999_999_999, 0xFFFF_FFFF, 0x1_0000_0000, 0xDEAD_BEEF_1234] {
            assert_eq!(Rns3::from_u64(v).to_u64(), v, "v={v}");
        }
    }

    #[test]
    fn add_agrees_with_naive() {
        let primes = SmallPrimes::PRIMES;
        let pairs = [(0u32, 1u32), (999, 1_073_754_112), (786_432, 786_432), (12345, 67890)];
        for (a, b) in pairs {
            let ra = Rns3::from_u32(a);
            let rb = Rns3::from_u32(b);
            let rc = ra + rb;
            for i in 0..3 {
                assert_eq!(rc.residue(i), (a as u64 + b as u64) as u32 % primes[i], "a={a}, b={b}, i={i}");
            }
        }
    }

    #[test]
    fn sub_agrees_with_naive() {
        let primes = SmallPrimes::PRIMES;
        let pairs = [(10u32, 3u32), (0, 1), (1_073_754_112, 1_073_754_112)];
        for (a, b) in pairs {
            let ra = Rns3::from_u32(a);
            let rb = Rns3::from_u32(b);
            let rc = ra - rb;
            for i in 0..3 {
                let expected = ((a as i64 - b as i64).rem_euclid(primes[i] as i64)) as u32;
                assert_eq!(rc.residue(i), expected, "a={a}, b={b}, i={i}");
            }
        }
    }

    #[test]
    fn mul_agrees_with_naive() {
        let primes = SmallPrimes::PRIMES;
        let pairs = [(0u32, 0u32), (1, 1), (2, 3), (786_432, 2), (999, 12345)];
        for (a, b) in pairs {
            let ra = Rns3::from_u32(a);
            let rb = Rns3::from_u32(b);
            let rc = ra * rb;
            for i in 0..3 {
                let expected = (a as u64 * b as u64 % primes[i] as u64) as u32;
                assert_eq!(rc.residue(i), expected, "a={a}, b={b}, i={i}");
            }
        }
    }

    #[test]
    fn zero_and_one() {
        let z = Rns3::zero();
        assert!(z.is_zero());
        for i in 0..3 {
            assert_eq!(z.residue(i), 0, "i={i}");
        }

        let o = Rns3::one();
        for i in 0..3 {
            assert_eq!(o.residue(i), 1, "i={i}");
        }
    }

    #[test]
    fn garner_reconstruction_matches_to_u64() {
        for v in [0u64, 1, 42, 0xDEAD_BEEF, 0xFFFF_FFFF_FF] {
            let r = Rns3::from_u64(v);
            // Manually apply Horner from the Garner coefficients.
            let a = r.to_garner();
            let primes = SmallPrimes::PRIMES;
            let mut horner = a[2] as u64;
            horner = horner * primes[1] as u64 + a[1] as u64;
            horner = horner * primes[0] as u64 + a[0] as u64;
            assert_eq!(horner, v, "v={v}");
        }
    }

    #[test]
    fn from_i32_positive() {
        let primes = SmallPrimes::PRIMES;
        for v in [0i32, 1, 12345, 786_432] {
            let r = Rns3::from_i32(v);
            for i in 0..3 {
                assert_eq!(r.residue(i), v as u32 % primes[i], "v={v}, i={i}");
            }
        }
    }

    #[test]
    fn from_i32_negative() {
        let primes = SmallPrimes::PRIMES;
        for v in [-1i32, -12345, -786_432] {
            let r = Rns3::from_i32(v);
            for i in 0..3 {
                let expected = v.rem_euclid(primes[i] as i32) as u32;
                assert_eq!(r.residue(i), expected, "v={v}, i={i}");
            }
        }
    }

    #[test]
    fn to_i64_roundtrip() {
        for v in [-1000i64, -1, 0, 1, 999, 786_432, -786_433] {
            let r = Rns3::from_i32(v as i32);
            assert_eq!(r.to_i64(), v, "v={v}");
        }
    }

    #[test]
    fn ntt_intt_roundtrip() {
        use super::{intt_inplace, ntt_inplace};
        let n = 16;
        let orig: Vec<Rns3> = (0..n as i32).map(Rns3::from_i32).collect();
        let mut buf = orig.clone();
        ntt_inplace::<3, SmallPrimes>(&mut buf);
        // After NTT the values should be different (unless poly is zero)
        assert_ne!(buf[0].residue(0), orig[0].residue(0));
        intt_inplace::<3, SmallPrimes>(&mut buf);
        // After INTT we should recover the original
        for (i, (a, b)) in orig.iter().zip(buf.iter()).enumerate() {
            for pi in 0..3 {
                assert_eq!(a.residue(pi), b.residue(pi), "coeff={i}, prime={pi}");
            }
        }
    }

    #[test]
    fn ntt_multiplication_mod_cyclotomic() {
        use super::{intt_inplace, ntt_inplace};
        // Multiply [1, 1, 0, ...] * [1, 1, 0, ...] mod X^4+1
        // = [1, 2, 1, 0] mod X^4+1 = [1, 2, 1, 0] (degree < 4, no reduction)
        let n = 4;
        let a_coeffs = [1i32, 1, 0, 0];
        let b_coeffs = [1i32, 1, 0, 0];
        let mut a: Vec<Rns3> = a_coeffs.iter().copied().map(Rns3::from_i32).collect();
        let mut b: Vec<Rns3> = b_coeffs.iter().copied().map(Rns3::from_i32).collect();
        ntt_inplace::<3, SmallPrimes>(&mut a);
        ntt_inplace::<3, SmallPrimes>(&mut b);
        let mut c: Vec<Rns3> = a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect();
        intt_inplace::<3, SmallPrimes>(&mut c);
        let expected = [1i32, 2, 1, 0];
        for (i, (e, r)) in expected.iter().zip(c.iter()).enumerate() {
            assert_eq!(*e, r.to_i64() as i32, "coeff={i}");
        }
    }

    #[test]
    fn large_prime_u128_branch() {
        // Exercise the LOG2R=32 branch (4_294_967_291 > 2^31).
        let v: u64 = 3_000_000_000;
        let r = Rns2L::from_u64(v);
        assert_eq!(r.residue(0), v as u32 % 1_073_754_113);
        assert_eq!(r.residue(1), (v % 4_294_967_291) as u32);
        assert_eq!(r.to_u64(), v);

        let a = Rns2L::from_u64(1_000_000_007);
        let b = Rns2L::from_u64(2_000_000_003);
        let c = a * b;
        for i in 0..2 {
            let p = LargePrimes::PRIMES[i] as u64;
            let expected = (1_000_000_007u64 * 2_000_000_003u64 % p) as u32;
            assert_eq!(c.residue(i), expected, "i={i}");
        }
    }
}
