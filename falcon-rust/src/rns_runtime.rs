//! Runtime-`K` RNS negacyclic multiply.
//!
//! The const-generic [`Rns<N, P>`](crate::rns::Rns) machinery bakes the prime
//! count `N` and the prime list into the type. The deep NTRU-solve recursion
//! needs a prime count that grows with depth (into the hundreds at the bottom —
//! the construction products reach ~11000 bits), which cannot be a compile-time
//! constant per call site. This module provides the same negacyclic NTT multiply
//! with the prime list chosen at run time, reusing the low-level per-prime
//! primitives from [`crate::rns`] (`ntt_u32`, `intt_u32`, `bitrev_powers_mont`,
//! `montmul_dyn`, `mod_pow`).
//!
//! A [`RuntimeNtt`] caches everything that depends only on `(n, primes)` — the
//! per-prime Montgomery constants, twiddle tables, and the Garner inverse table
//! — so it is built once per `(depth)` and reused across every multiply / babai
//! pass at that depth. The product of two [`MultiwordPoly`] operands is
//! reconstructed straight into `MultiwordInt` limbs by a runtime word-level
//! Garner CRT, with no `BigInt`.

use std::collections::HashMap;
use std::sync::{Arc, LazyLock, Mutex};

use crate::fp_field::{negqinv_modr, r_sq_modq};
use crate::multiword_int::MultiwordInt;
use crate::multiword_poly::MultiwordPoly;
use crate::rns::{bitrev_powers_mont, intt_u32, mod_pow, montmul_dyn, ntt_u32, primitive_root_2n};

/// A pool of 24-bit NTT-friendly primes (`p ≡ 1 mod 2048`, descending from just
/// below `2^24`), generated once and shared. There are ~1047 such primes; the
/// deepest construction product needs ~478, so the pool is grown on demand up to
/// whatever a run requires.
static PRIME_POOL: LazyLock<Mutex<Vec<u32>>> = LazyLock::new(|| Mutex::new(Vec::new()));

fn is_prime_u32(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n % 2 == 0 {
        return n == 2;
    }
    let mut i = 3u64;
    while i * i <= n as u64 {
        if n as u64 % i == 0 {
            return false;
        }
        i += 2;
    }
    true
}

/// Return the first `k` 24-bit NTT primes (`≡ 1 mod 2048`), descending from
/// `2^24`. Cached and grown on demand.
fn ntt_primes(k: usize) -> Vec<u32> {
    let mut pool = PRIME_POOL.lock().unwrap();
    if pool.len() < k {
        // Largest 24-bit value ≡ 1 mod 2048 is 2^24 - 2048 + 1.
        let mut cand: i64 = ((1i64 << 24) - 2048) + 1;
        // Resume below the smallest prime already found.
        if let Some(&last) = pool.last() {
            cand = last as i64 - 2048;
        }
        while pool.len() < k && cand > 1 {
            if is_prime_u32(cand as u32) {
                pool.push(cand as u32);
            }
            cand -= 2048;
        }
        assert!(
            pool.len() >= k,
            "exhausted 24-bit NTT primes: needed {k}, found {}",
            pool.len()
        );
    }
    pool[..k].to_vec()
}

/// Find a primitive 2048th root of unity modulo `p` (`p ≡ 1 mod 2048`) at run
/// time: take a generator candidate `g` and raise it to `(p-1)/2048`, retrying
/// until the result has full order 2048 (its 1024th power ≠ 1).
fn primitive_root_2048(p: u32) -> u32 {
    let order = (p - 1) as u64;
    let exp = order / 2048;
    let mut g = 2u64;
    loop {
        let root = mod_pow(g, exp, p as u64);
        // Primitive 2048th root iff root^1024 ≠ 1 (so order does not divide 1024).
        if root != 1 && mod_pow(root, 1024, p as u64) != 1 {
            return root as u32;
        }
        g += 1;
    }
}

/// An operand forward-transformed once for reuse across multiplies: per-prime
/// NTT residues in Montgomery form, indexed `[prime][coefficient]`.
pub(crate) struct Transformed {
    ntt: Vec<Vec<u32>>,
}

/// Process-wide cache of [`RuntimeNtt`] contexts keyed by `(n, k)`. Building one
/// costs `k` root-finds + twiddle tables + an O(k²) Garner inverse table, so the
/// same depth's context is built once and shared across keygen attempts.
static NTT_CACHE: LazyLock<Mutex<HashMap<(usize, usize), Arc<RuntimeNtt>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

/// Per-`(n, primes)` precomputed RNS multiply context.
pub(crate) struct RuntimeNtt {
    n: usize,
    primes: Vec<u32>,
    log2r: Vec<u32>,
    neg_inv: Vec<u32>,
    r_sq: Vec<u32>,            // R^2 mod p (to convert into Montgomery form)
    r64: Vec<u32>,            // 2^64 mod p (for multiword reduction)
    psi_rev: Vec<Vec<u32>>,    // forward twiddles, Montgomery, bit-reversed
    psi_inv_rev: Vec<Vec<u32>>, // inverse twiddles
    ninv_mont: Vec<u32>,       // n^{-1}, Montgomery
    /// Garner inverses: `garner_inv[i][j] = primes[i]^{-1} mod primes[j]` for
    /// `j > i` (lower triangle unused). Flattened row-major.
    garner_inv: Vec<u32>,
}

impl RuntimeNtt {
    /// Build the context for an NTT of length `n` over the first `k` 24-bit NTT
    /// primes.
    pub(crate) fn new(n: usize, k: usize) -> Self {
        debug_assert!(n.is_power_of_two());
        let primes = ntt_primes(k);
        let log2r: Vec<u32> = primes.iter().map(|&p| p.ilog2() + 1).collect();
        let neg_inv: Vec<u32> = primes.iter().map(|&p| negqinv_modr(p)).collect();
        let r_sq: Vec<u32> = primes.iter().map(|&p| r_sq_modq(p)).collect();
        let r64: Vec<u32> = primes
            .iter()
            .map(|&p| ((1u128 << 64) % p as u128) as u32)
            .collect();

        let mut psi_rev = Vec::with_capacity(k);
        let mut psi_inv_rev = Vec::with_capacity(k);
        let mut ninv_mont = Vec::with_capacity(k);
        for i in 0..k {
            let p = primes[i];
            let root_2n = primitive_root_2n(primitive_root_2048(p), n, p);
            let root_2n_mont = montmul_dyn(root_2n, r_sq[i], p, log2r[i], neg_inv[i]);
            psi_rev.push(bitrev_powers_mont(
                root_2n_mont,
                n,
                p,
                log2r[i],
                neg_inv[i],
                r_sq[i],
            ));
            let root_inv = mod_pow(root_2n as u64, (p - 2) as u64, p as u64) as u32;
            let root_inv_mont = montmul_dyn(root_inv, r_sq[i], p, log2r[i], neg_inv[i]);
            psi_inv_rev.push(bitrev_powers_mont(
                root_inv_mont,
                n,
                p,
                log2r[i],
                neg_inv[i],
                r_sq[i],
            ));
            let n_inv = mod_pow(n as u64, (p - 2) as u64, p as u64) as u32;
            ninv_mont.push(montmul_dyn(n_inv, r_sq[i], p, log2r[i], neg_inv[i]));
        }

        // Garner inverse table: primes[i]^{-1} mod primes[j], j > i.
        let mut garner_inv = vec![0u32; k * k];
        for i in 0..k {
            for j in (i + 1)..k {
                let pj = primes[j] as u64;
                garner_inv[i * k + j] = mod_pow(primes[i] as u64 % pj, pj - 2, pj) as u32;
            }
        }

        Self {
            n,
            primes,
            log2r,
            neg_inv,
            r_sq,
            r64,
            psi_rev,
            psi_inv_rev,
            ninv_mont,
            garner_inv,
        }
    }

    /// Cached [`RuntimeNtt`] for `(n, k)`, built once per process.
    pub(crate) fn cached(n: usize, k: usize) -> Arc<RuntimeNtt> {
        let mut cache = NTT_CACHE.lock().unwrap();
        cache
            .entry((n, k))
            .or_insert_with(|| Arc::new(RuntimeNtt::new(n, k)))
            .clone()
    }

    pub(crate) fn num_primes(&self) -> usize {
        self.primes.len()
    }

    pub(crate) fn primes(&self) -> &[u32] {
        &self.primes
    }

    /// Reduce a signed two's-complement limb slice modulo prime index `pi`,
    /// returning the canonical residue in `[0, p)`. No allocation.
    ///
    /// Treating the limbs as an unsigned `64·len`-bit integer `U` (Horner over
    /// base `2^64`, most-significant limb first) gives `U mod p`. For a negative
    /// value the true value is `U − 2^(64·len)`, so subtract `2^(64·len) mod p`.
    fn reduce_signed(&self, limbs: &[u64], pi: usize) -> u32 {
        let p = self.primes[pi] as u128;
        let r64 = self.r64[pi] as u128;
        let mut acc = 0u128;
        for &limb in limbs.iter().rev() {
            acc = (acc * r64 + (limb as u128 % p)) % p;
        }
        if crate::multiword_int::is_negative(limbs) {
            let two_pow = mod_pow(2, 64 * limbs.len() as u64, self.primes[pi] as u64) as u128;
            acc = (acc + p - two_pow) % p;
        }
        acc as u32
    }

    /// Forward-transform an operand once: per-prime NTT residues (Montgomery
    /// form), reusable across many multiplies. Hoisting this out of the babai
    /// loop (where the same `f`, `g` multiply every `k`) avoids re-reducing and
    /// re-transforming them on every pass.
    pub(crate) fn forward(&self, a: &MultiwordPoly) -> Transformed {
        debug_assert_eq!(a.n(), self.n);
        let n = self.n;
        let mut ntt = Vec::with_capacity(self.primes.len());
        for pi in 0..self.primes.len() {
            let (p, log2r, neg_inv) = (self.primes[pi], self.log2r[pi], self.neg_inv[pi]);
            let mut ar = vec![0u32; n];
            for i in 0..n {
                ar[i] = montmul_dyn(self.reduce_signed(a.coeff(i), pi), self.r_sq[pi], p, log2r, neg_inv);
            }
            ntt_u32(&mut ar, &self.psi_rev[pi], p, log2r, neg_inv);
            ntt.push(ar);
        }
        Transformed { ntt }
    }

    /// Negacyclic product `k * h mod (X^n + 1)` where `h` is already
    /// forward-transformed (see [`forward`](Self::forward)); only `k` is
    /// transformed here. Each product coefficient is reconstructed into an
    /// `out_w`-limb `MultiwordInt`.
    pub(crate) fn mul_transformed(
        &self,
        kp: &MultiwordPoly,
        h: &Transformed,
        out_w: usize,
    ) -> MultiwordPoly {
        debug_assert!(kp.n() == self.n && h.ntt.len() == self.primes.len());
        let n = self.n;
        let k = self.primes.len();
        let mut residues = vec![0u32; n * k];
        let mut kr = vec![0u32; n];
        for pi in 0..k {
            let (p, log2r, neg_inv) = (self.primes[pi], self.log2r[pi], self.neg_inv[pi]);
            for i in 0..n {
                kr[i] = montmul_dyn(self.reduce_signed(kp.coeff(i), pi), self.r_sq[pi], p, log2r, neg_inv);
            }
            ntt_u32(&mut kr, &self.psi_rev[pi], p, log2r, neg_inv);
            for i in 0..n {
                kr[i] = montmul_dyn(kr[i], h.ntt[pi][i], p, log2r, neg_inv);
            }
            intt_u32(&mut kr, &self.psi_inv_rev[pi], self.ninv_mont[pi], p, log2r, neg_inv);
            for i in 0..n {
                residues[i * k + pi] = montmul_dyn(kr[i], 1, p, log2r, neg_inv);
            }
        }

        // Reconstruct each coefficient: runtime Garner -> mixed-radix digits ->
        // signed word-level CRT into the limbs.
        let mut out = MultiwordPoly::zeros(n, out_w);
        let mut digits = vec![0u32; k];
        for i in 0..n {
            digits.copy_from_slice(&residues[i * k..(i + 1) * k]);
            self.garner_in_place(&mut digits);
            let value = MultiwordInt::from_garner_digits(&digits, &self.primes, out_w);
            out.coeff_mut(i).copy_from_slice(value.limbs());
        }
        out
    }

    /// Negacyclic product `a * b mod (X^n + 1)` (transforms both operands).
    /// Convenience over [`forward`](Self::forward) + [`mul_transformed`](Self::mul_transformed)
    /// for one-shot multiplies and tests.
    pub(crate) fn negacyclic_mul(
        &self,
        a: &MultiwordPoly,
        b: &MultiwordPoly,
        out_w: usize,
    ) -> MultiwordPoly {
        self.mul_transformed(a, &self.forward(b), out_w)
    }

    /// In-place Garner: convert canonical residues `u[j] = x mod primes[j]` to
    /// mixed-radix digits `a[j]` with `x = a_0 + p_0(a_1 + p_1(…))`.
    fn garner_in_place(&self, u: &mut [u32]) {
        let k = self.primes.len();
        for i in 0..k {
            for j in (i + 1)..k {
                let pj = self.primes[j] as u64;
                let inv = self.garner_inv[i * k + j] as u64;
                let diff = (u[j] as u64 + pj - (u[i] as u64 % pj)) % pj;
                u[j] = (diff * inv % pj) as u32;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::RuntimeNtt;
    use crate::multiword_poly::MultiwordPoly;
    use crate::polynomial::Polynomial;
    use num::BigInt;
    use rand::{rngs::StdRng, RngExt, SeedableRng};

    /// Schoolbook negacyclic product of two BigInt polynomials (the oracle).
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

    fn rand_poly(rng: &mut StdRng, n: usize, bits: u32) -> Polynomial<BigInt> {
        let words = bits.div_ceil(32);
        let excess = words * 32 - bits; // trim the top word to exactly `bits` bits
        Polynomial::new(
            (0..n)
                .map(|_| {
                    let mut v = BigInt::from(0);
                    for _ in 0..words {
                        v = (v << 32) + BigInt::from(rng.random::<u32>());
                    }
                    v >>= excess;
                    if rng.random::<bool>() {
                        -v
                    } else {
                        v
                    }
                })
                .collect(),
        )
    }

    fn check(n: usize, bits: u32, k: usize, out_w: usize, seed: u64) {
        let mut rng = StdRng::seed_from_u64(seed);
        let ctx = RuntimeNtt::new(n, k);
        for _ in 0..20 {
            let a = rand_poly(&mut rng, n, bits);
            let b = rand_poly(&mut rng, n, bits);
            let want = negacyclic_bigint(&a, &b);
            let ma = MultiwordPoly::from_bigint_poly(&a, out_w);
            let mb = MultiwordPoly::from_bigint_poly(&b, out_w);
            let got = ctx.negacyclic_mul(&ma, &mb, out_w);
            assert_eq!(
                got.to_bigint_poly().coefficients,
                want.coefficients,
                "n={n} bits={bits} k={k}"
            );
        }
    }

    #[test]
    fn small_operands_few_primes() {
        // ~50-bit operands, product ~106 bits -> 5 primes, 4 limbs out.
        check(128, 50, 5, 4, 1);
    }

    #[test]
    fn medium_operands() {
        // ~200-bit operands, product ~410 bits -> 19 primes, 8 limbs.
        check(32, 200, 19, 8, 2);
    }

    #[test]
    fn deep_small_n_many_primes() {
        // Tiny n, large operands: ~600-bit, product ~1210 bits -> 54 primes.
        check(4, 600, 54, 20, 3);
    }

    #[test]
    fn asymmetric_widths() {
        // One large, one small operand (like F' x g_minx construction).
        let mut rng = StdRng::seed_from_u64(4);
        let n = 8;
        let ctx = RuntimeNtt::new(n, 40);
        for _ in 0..20 {
            let a = rand_poly(&mut rng, n, 400);
            let b = rand_poly(&mut rng, n, 60);
            let want = negacyclic_bigint(&a, &b);
            let ma = MultiwordPoly::from_bigint_poly(&a, 16);
            let mb = MultiwordPoly::from_bigint_poly(&b, 16);
            let got = ctx.negacyclic_mul(&ma, &mb, 16);
            assert_eq!(got.to_bigint_poly().coefficients, want.coefficients);
        }
    }
}
