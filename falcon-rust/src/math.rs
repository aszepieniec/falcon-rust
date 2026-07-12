use std::vec::IntoIter;

use falcon_profiler::profiling;
use itertools::Itertools;
use num::{BigInt, FromPrimitive, One, Zero};
use num_complex::Complex64;
use rand::Rng;

use crate::{
    falcon_field::{Felt, Q},
    fast_fft::FastFft,
    fixed_point::FixedPoint128,
    multiword_int::MultiwordInt,
    multiword_poly::MultiwordPoly,
    packed::Packed,
    polynomial::Polynomial,
    rns::{
        intt_inplace_cached, ntt_inplace_cached, NttPrimeList, NttPrimes24Bit2, NttPrimes24Bit4,
        NttPrimes24Bit5, NttPrimes24Bit8, Rns,
    },
    rns_runtime::RuntimeNtt,
    samplerz::SamplerZCtx,
    U32Field,
};

/// Coefficient sizes (log₂ of max absolute value, in bits) at each recursion
/// depth in [`ntru_solve`], from the Falcon specification
/// (<https://falcon-sign.info/falcon.pdf>, Table on p. 59), gathered
/// experimentally over many Falcon-1024 key-generation runs.
///
/// Depth 0 = outermost call (n = 1024); depth 10 = base case (n = 1, xgcd
/// only, no Babai step).  Columns: `(avg log₂ max|f,g|, σ, avg log₂
/// max|F,G|, σ)`, where F, G are measured before `babai_reduce` runs.
#[rustfmt::skip]
pub const NTRU_SOLVE_BABAI_COEFF_BITS: [(f64, f64, f64, f64); 11] = [
    ( 4.00,  0.00,  19.61,  0.49), // depth  0  n = 1024 (babai_reduce_i32)
    (10.99,  0.08,  39.82,  0.41), // depth  1  n =  512
    (24.07,  0.25,  78.20,  0.73), // depth  2  n =  256
    (50.37,  0.53, 153.65,  1.39), // depth  3  n =  128
    (101.62, 1.02, 303.49,  2.38), // depth  4  n =   64
    (202.22, 1.87, 599.81,  3.87), // depth  5  n =   32
    (400.67, 3.10,1188.68,  6.04), // depth  6  n =   16
    (794.17, 4.98,2361.84,  9.31), // depth  7  n =    8
    (1576.87,7.49,4703.30, 14.77), // depth  8  n =    4
    (3138.35,12.25,9403.29,27.55), // depth  9  n =    2
    (6307.52,24.48,6319.66,24.51), // depth 10  n =    1 (xgcd, no Babai)
];

/// Window / back-shift schedule for one Babai reduction pass, shared by all
/// three estimator loops ([`babai_reduce_bigint`], [`babai_reduce_rns_generic`],
/// [`babai_reduce_rns_runtime`]) so they stay bit-identical (the differential
/// tests depend on it).
///
/// The per-coefficient quotient `k_true = C / f` has `d = capital_size − size`
/// bits, but the f64 `round()` resolves only ~53 bits at once:
/// - `d > 53`: window the top 106 bits (`capital_size − 106`) so the quotient is
///   ≈ 2^53 — k_true's top 53 bits — and `back_shift = d − 53` re-aligns them.
/// - `d ≤ 53`: k_true already fits in 53 bits, so window at the *denominator's*
///   scale (`size − 53`), making the quotient equal k_true exactly and reducing
///   all `d` bits in ONE pass, with `back_shift = 0`.
///
/// The earlier `d ≤ 53` form `(capital_size − 53, d)` scaled the quotient down to
/// ≈ 2^(53 − size), i.e. only a few bits whenever `size` sat near its 53-bit
/// floor (e.g. depth 3, where `max|f,g| ≈ 50`), so it crawled ~1 bit per pass —
/// 24 wasted passes of an expensive n=128 multiply. Windowing at `size − 53`
/// removes the crawl; f64 precision was never the limit (k_true is only ~28 bits
/// there).
fn babai_shifts(capital_size: u64, size: u64, d: u64) -> (u32, u32) {
    if d > 53 {
        ((capital_size - 106) as u32, (d - 53) as u32)
    } else {
        ((size - 53) as u32, 0)
    }
}

/// Reduce the vector (F,G) relative to (f,g). This method follows the python
/// implementation [1].
///
/// Algorithm 7 in the spec [2, p.35]
///
/// [1]: https://github.com/tprest/falcon.py
///
/// [2]: https://falcon-sign.info/falcon.pdf
///
/// This function is marked pub for the purpose of benchmarking; it is not
/// considered part of the public API.
#[doc(hidden)]
#[profiling]
pub fn babai_reduce_bigint(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    let bitsize = |bi: &BigInt| bi.bits();
    let n = f.coefficients.len();
    let size = [
        f.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
        g.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
        53,
    ]
    .into_iter()
    .max()
    .unwrap();
    let shift = (size as i64) - 53;
    let f_adjusted = f
        .map(|bi| Complex64::new(i64::try_from(bi >> shift).unwrap() as f64, 0.0))
        .fft();
    let g_adjusted = g
        .map(|bi| Complex64::new(i64::try_from(bi >> shift).unwrap() as f64, 0.0))
        .fft();

    let f_star_adjusted = f_adjusted.map(|c| c.conj());
    let g_star_adjusted = g_adjusted.map(|c| c.conj());
    let denominator_fft =
        f_adjusted.hadamard_mul(&f_star_adjusted) + g_adjusted.hadamard_mul(&g_star_adjusted);

    let mut prev_capital_size = u64::MAX;
    loop {
        let capital_size = [
            capital_f.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
            capital_g.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
            53,
        ]
        .into_iter()
        .max()
        .unwrap();

        // Stop when we've reached the target size, or when capital_size stopped
        // strictly decreasing (floating-point precision limit reached).
        if capital_size < size {
            break;
        }
        if capital_size >= prev_capital_size {
            break;
        }
        prev_capital_size = capital_size;

        // See `babai_shifts` for the windowing schedule (two-step for d > 53,
        // full-k single pass for d <= 53).
        let d = capital_size - size;
        let (capital_shift, back_shift) = babai_shifts(capital_size, size, d);
        let (capital_shift, back_shift) = (capital_shift as i64, back_shift as i64);

        let capital_f_adjusted = capital_f
            .map(|bi| Complex64::new(i128::try_from(bi >> capital_shift).unwrap() as f64, 0.0))
            .fft();
        let capital_g_adjusted = capital_g
            .map(|bi| Complex64::new(i128::try_from(bi >> capital_shift).unwrap() as f64, 0.0))
            .fft();

        let numerator = capital_f_adjusted.hadamard_mul(&f_star_adjusted)
            + capital_g_adjusted.hadamard_mul(&g_star_adjusted);
        let quotient = numerator.hadamard_div(&denominator_fft).ifft();

        // Use i128 to avoid i64 saturation when the FFT quotient is large
        // (can happen for small n when |f_fft[i]|^2 + |g_fft[i]|^2 is small at
        // some frequency).
        let k = quotient.map(|f| BigInt::from(f.re.round() as i128));

        if k.is_zero() {
            break;
        }
        let kf = (k.clone().karatsuba(f)).reduce_by_cyclotomic(n);
        let kg = (k.clone().karatsuba(g)).reduce_by_cyclotomic(n);
        let shifted_kf = kf.map(|bi| bi << back_shift);
        let shifted_kg = kg.map(|bi| bi << back_shift);

        // Tentative check: if applying shifted_kf would make capital_F grow,
        // the step overshot (common when |f_fft|^2+|g_fft|^2 is very small at
        // one frequency for small n).  In that case fall back to the old
        // single-bit formula for this iteration.
        if d > 53 {
            let new_cs_f = capital_f
                .coefficients
                .iter()
                .zip(shifted_kf.coefficients.iter())
                .map(|(a, b)| (a - b).bits())
                .max()
                .unwrap_or(0);
            let new_cs_g = capital_g
                .coefficients
                .iter()
                .zip(shifted_kg.coefficients.iter())
                .map(|(a, b)| (a - b).bits())
                .max()
                .unwrap_or(0);
            if u64::max(new_cs_f, new_cs_g) >= capital_size {
                // Recompute with old formula (capital_shift = capital_size - 53)
                let cs_old = (capital_size as i64) - 53;
                let cf_old = capital_f
                    .map(|bi| Complex64::new(i64::try_from(bi >> cs_old).unwrap() as f64, 0.0))
                    .fft();
                let cg_old = capital_g
                    .map(|bi| Complex64::new(i64::try_from(bi >> cs_old).unwrap() as f64, 0.0))
                    .fft();
                let num_old =
                    cf_old.hadamard_mul(&f_star_adjusted) + cg_old.hadamard_mul(&g_star_adjusted);
                let quot_old = num_old.hadamard_div(&denominator_fft).ifft();
                let k_old = quot_old.map(|f| BigInt::from(f.re.round() as i64));
                if k_old.is_zero() {
                    break;
                }
                let kf_old = (k_old.clone().karatsuba(f)).reduce_by_cyclotomic(n);
                let kg_old = (k_old.karatsuba(g)).reduce_by_cyclotomic(n);
                *capital_f -= kf_old.map(|bi| bi << d);
                *capital_g -= kg_old.map(|bi| bi << d);
                continue;
            }
        }

        *capital_f -= shifted_kf;
        *capital_g -= shifted_kg;
    }
    Ok(())
}

/// Capital-coefficient representation for the generic RNS Babai reduction
/// [`babai_reduce_rns_generic`].  Implemented for `BigInt` (uncapped — the
/// `babai_reduce_rns_bigint` path) and [`Packed<L>`] (fixed-width limbs — the
/// faster `babai_reduce_rns_packed` path).  The reference `babai_reduce_bigint`
/// (which multiplies via karatsuba, not the NTT) is intentionally kept separate
/// as the correctness oracle.
trait ReductionCapital: Clone {
    fn from_bigint(x: &BigInt) -> Self;
    fn to_bigint(&self) -> BigInt;
    /// Bit-length of the absolute value (matches `BigInt::bits`).
    fn bit_length(&self) -> u64;
    /// Arithmetic right shift by `shift`, low 128 bits as `i128` — the windowed
    /// top-word input to the `f64` k-estimation FFT.
    fn shr_window(&self, shift: u32) -> i128;
    /// Left shift by `bits` (multiply by `2^bits`).
    fn shl(&self, bits: u32) -> Self;
    /// `self - other`.
    fn sub(&self, other: &Self) -> Self;
    /// Reconstruct one `k·h` product coefficient from its RNS residues.
    fn from_rns<const K: usize, P: NttPrimeList<K>>(r: &Rns<K, P>) -> Self;
}

impl ReductionCapital for BigInt {
    fn from_bigint(x: &BigInt) -> Self {
        x.clone()
    }
    fn to_bigint(&self) -> BigInt {
        self.clone()
    }
    fn bit_length(&self) -> u64 {
        self.bits()
    }
    fn shr_window(&self, shift: u32) -> i128 {
        i128::try_from(self >> (shift as i64)).unwrap()
    }
    fn shl(&self, bits: u32) -> Self {
        self << (bits as u64)
    }
    fn sub(&self, other: &Self) -> Self {
        self - other
    }
    fn from_rns<const K: usize, P: NttPrimeList<K>>(r: &Rns<K, P>) -> Self {
        r.to_bigint()
    }
}

impl<const L: usize> ReductionCapital for Packed<L> {
    fn from_bigint(x: &BigInt) -> Self {
        Packed::from_bigint(x)
    }
    fn to_bigint(&self) -> BigInt {
        Packed::to_bigint(self)
    }
    fn bit_length(&self) -> u64 {
        Packed::bit_length(self)
    }
    fn shr_window(&self, shift: u32) -> i128 {
        self.shr_to_i128(shift)
    }
    fn shl(&self, bits: u32) -> Self {
        Packed::shl(self, bits)
    }
    fn sub(&self, other: &Self) -> Self {
        Packed::sub(self, other)
    }
    fn from_rns<const K: usize, P: NttPrimeList<K>>(r: &Rns<K, P>) -> Self {
        Packed::from_rns(r)
    }
}

/// Babai reduction where the capital coefficients live in `C` (a positional
/// multi-word representation) and the `k·f`/`k·g` products are evaluated in RNS
/// via the NTT instead of `BigInt` karatsuba.
///
/// This is the fn-dsa-style split, generic over the capital representation:
/// the capital stays positional so the reduction coefficient `k` can be
/// estimated from its top words (the `f64` FFT windowing, shared verbatim with
/// [`babai_reduce_bigint`]) and is not capped at 127 bits, while the one
/// expensive operation — the polynomial product — runs in the RNS/NTT domain
/// (`f`, `g`, `k` only).  `k` is transformed once per iteration and reused for
/// both products.  Unlike [`babai_reduce_rns`] (capital in RNS, i128-capped to
/// depths 1–2), this works at any depth for which `P` covers the `k·f` product.
fn babai_reduce_rns_generic<C, const K: usize, P>(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String>
where
    C: ReductionCapital,
    P: NttPrimeList<K>,
{
    let bitsize = |bi: &BigInt| bi.bits();
    let n = f.coefficients.len();

    // Precompute the per-prime twiddle tables once; reused by every transform.
    // Production prime lists return compile-time tables; see `NttPrimeList::ntt_tables`.
    let tables = P::ntt_tables(n);

    // Forward NTT of f and g once (reused every iteration).  `from_bigint`
    // works at any depth, including where the coefficients exceed i128.
    let mut f_ntt: Vec<Rns<K, P>> = f
        .coefficients
        .iter()
        .map(Rns::<K, P>::from_bigint)
        .collect();
    let mut g_ntt: Vec<Rns<K, P>> = g
        .coefficients
        .iter()
        .map(Rns::<K, P>::from_bigint)
        .collect();
    ntt_inplace_cached::<K, P>(&mut f_ntt, &tables);
    ntt_inplace_cached::<K, P>(&mut g_ntt, &tables);

    let size = [
        f.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
        g.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
        53,
    ]
    .into_iter()
    .max()
    .unwrap();

    // Guard: the modulus must cover the `k·f` product.  With `k` windowed to
    // ~53 bits the product is `≈ max|f,g| + 54` bits (the measured model — see
    // `product_size_model_matches_measurements`); a too-small prime list `P`
    // for this depth would silently wrap the centered CRT, so catch it here in
    // debug builds instead of corrupting the reduction in release.
    let modulus_signed_bits: f64 = P::PRIMES.iter().map(|&p| (p as f64).log2()).sum::<f64>() - 1.0;
    debug_assert!(
        modulus_signed_bits >= (size + 54) as f64,
        "RNS prime list too small: modulus {modulus_signed_bits:.0} signed bits < \
         product ~{} bits (max|f,g|={size})",
        size + 54
    );

    let shift = (size as i64) - 53;
    let f_adjusted = f
        .map(|bi| Complex64::new(i64::try_from(bi >> shift).unwrap() as f64, 0.0))
        .fft();
    let g_adjusted = g
        .map(|bi| Complex64::new(i64::try_from(bi >> shift).unwrap() as f64, 0.0))
        .fft();

    let f_star_adjusted = f_adjusted.map(|c| c.conj());
    let g_star_adjusted = g_adjusted.map(|c| c.conj());
    let denominator_fft =
        f_adjusted.hadamard_mul(&f_star_adjusted) + g_adjusted.hadamard_mul(&g_star_adjusted);

    // Capital coefficients move into the positional representation `C` for the
    // duration of the loop; `BigInt` is touched only at entry and exit.
    let mut cf: Vec<C> = capital_f.coefficients.iter().map(C::from_bigint).collect();
    let mut cg: Vec<C> = capital_g.coefficients.iter().map(C::from_bigint).collect();

    // Pointwise-multiply the (already forward-transformed) `k_ntt` by `h_ntt`,
    // inverse-transform, and reconstruct each coefficient into `C`.
    let mul = |k_ntt: &[Rns<K, P>], h_ntt: &[Rns<K, P>]| -> Vec<C> {
        let mut prod: Vec<Rns<K, P>> = k_ntt
            .iter()
            .zip(h_ntt.iter())
            .map(|(&a, &b)| a * b)
            .collect();
        intt_inplace_cached::<K, P>(&mut prod, &tables);
        prod.iter().map(C::from_rns).collect()
    };
    // Windowed top-word read of a capital vector, transformed for k-estimation.
    let window = |cs: &[C], sh: u32| -> Polynomial<Complex64> {
        Polynomial::new(
            cs.iter()
                .map(|c| Complex64::new(c.shr_window(sh) as f64, 0.0))
                .collect::<Vec<_>>(),
        )
        .fft()
    };
    let cap_size = |cf: &[C], cg: &[C]| -> u64 {
        cf.iter()
            .chain(cg.iter())
            .map(C::bit_length)
            .fold(53, u64::max)
    };

    let mut prev_capital_size = u64::MAX;
    loop {
        let capital_size = cap_size(&cf, &cg);
        if capital_size < size || capital_size >= prev_capital_size {
            break;
        }
        prev_capital_size = capital_size;

        let d = capital_size - size;
        let (capital_shift, back_shift) = babai_shifts(capital_size, size, d);

        let numerator = window(&cf, capital_shift).hadamard_mul(&f_star_adjusted)
            + window(&cg, capital_shift).hadamard_mul(&g_star_adjusted);
        let quotient = numerator.hadamard_div(&denominator_fft).ifft();

        let k: Vec<i128> = quotient
            .coefficients
            .iter()
            .map(|c| c.re.round() as i128)
            .collect();
        if k.iter().all(|&x| x == 0) {
            break;
        }

        // Transform `k` once, reuse for both the k·f and k·g products.
        let mut k_ntt: Vec<Rns<K, P>> = k.iter().map(|&v| Rns::from_i128(v)).collect();
        ntt_inplace_cached::<K, P>(&mut k_ntt, &tables);
        let kf = mul(&k_ntt, &f_ntt);
        let kg = mul(&k_ntt, &g_ntt);
        let shifted_kf: Vec<C> = kf.iter().map(|p| p.shl(back_shift)).collect();
        let shifted_kg: Vec<C> = kg.iter().map(|p| p.shl(back_shift)).collect();

        if d > 53 {
            let new_cs_f = cf
                .iter()
                .zip(shifted_kf.iter())
                .map(|(a, b)| a.sub(b).bit_length())
                .fold(0, u64::max);
            let new_cs_g = cg
                .iter()
                .zip(shifted_kg.iter())
                .map(|(a, b)| a.sub(b).bit_length())
                .fold(0, u64::max);
            if u64::max(new_cs_f, new_cs_g) >= capital_size {
                let cs_old = (capital_size - 53) as u32;
                let num_old = window(&cf, cs_old).hadamard_mul(&f_star_adjusted)
                    + window(&cg, cs_old).hadamard_mul(&g_star_adjusted);
                let quot_old = num_old.hadamard_div(&denominator_fft).ifft();
                let k_old: Vec<i128> = quot_old
                    .coefficients
                    .iter()
                    .map(|c| c.re.round() as i128)
                    .collect();
                if k_old.iter().all(|&x| x == 0) {
                    break;
                }
                let mut k_old_ntt: Vec<Rns<K, P>> =
                    k_old.iter().map(|&v| Rns::from_i128(v)).collect();
                ntt_inplace_cached::<K, P>(&mut k_old_ntt, &tables);
                let kf_old = mul(&k_old_ntt, &f_ntt);
                let kg_old = mul(&k_old_ntt, &g_ntt);
                for i in 0..n {
                    cf[i] = cf[i].sub(&kf_old[i].shl(d as u32));
                    cg[i] = cg[i].sub(&kg_old[i].shl(d as u32));
                }
                continue;
            }
        }

        for i in 0..n {
            cf[i] = cf[i].sub(&shifted_kf[i]);
            cg[i] = cg[i].sub(&shifted_kg[i]);
        }
    }

    // Move the reduced capital coefficients back out to BigInt.
    for (c, p) in capital_f.coefficients.iter_mut().zip(cf.iter()) {
        *c = p.to_bigint();
    }
    for (c, p) in capital_g.coefficients.iter_mut().zip(cg.iter()) {
        *c = p.to_bigint();
    }
    Ok(())
}

/// Runtime-`K` Babai reduction: capital lives in [`MultiwordPoly`] limbs and the
/// `k·f`/`k·g` corrections are formed by the runtime-prime-count RNS multiply
/// ([`RuntimeNtt`]), so it works at any recursion depth without a compile-time
/// prime list.  Algorithmically identical to [`babai_reduce_rns_generic`] (same
/// windowed k-estimation, same `d > 53` two-step), differing only in the
/// representation: capital and operands are flat words, no `BigInt` in the loop.
///
/// `k_primes` must cover the `k·f` product (`≈ max|f,g| + 54` bits); `cap_w` is
/// the capital limb width (must hold the capital plus the mid-reduction shifts).
/// `BigInt` is touched only at entry (operands, capital) and exit (capital).
#[profiling]
pub(crate) fn babai_reduce_rns_runtime(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
    k_primes: usize,
    cap_w: usize,
) -> Result<(), String> {
    use crate::multiword_int as mw;
    use crate::multiword_poly::SCHOOLBOOK_MAX_N;
    let n = f.coefficients.len();
    // At deep levels (tiny `n`, huge coefficients) the flat-word schoolbook
    // multiply beats the RNS/NTT one, and building the `RuntimeNtt` — `k_primes`
    // root-finds plus an O(k²) Garner table, with `k_primes` reaching ~149 at
    // depth 9 — is pure waste.  Only build it in the large-`n` (RNS) regime.
    let use_schoolbook = n <= SCHOOLBOOK_MAX_N;
    let ntt = if use_schoolbook {
        None
    } else {
        Some(RuntimeNtt::cached(n, k_primes))
    };

    let bitsize = |bi: &BigInt| bi.bits();
    let size = [
        f.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
        g.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
        53,
    ]
    .into_iter()
    .max()
    .unwrap();

    // f64 FFT denominator from windowed f,g (identical to the BigInt path).
    let shift = (size as i64) - 53;
    let f_adjusted = f
        .map(|bi| Complex64::new(i64::try_from(bi >> shift).unwrap() as f64, 0.0))
        .fft();
    let g_adjusted = g
        .map(|bi| Complex64::new(i64::try_from(bi >> shift).unwrap() as f64, 0.0))
        .fft();
    let f_star = f_adjusted.map(|c| c.conj());
    let g_star = g_adjusted.map(|c| c.conj());
    let denominator_fft = f_adjusted.hadamard_mul(&f_star) + g_adjusted.hadamard_mul(&g_star);

    // Operands as flat-word polynomials (multiply inputs); capital likewise.
    let fw = (size as usize / 64) + 2;
    let f_mw = MultiwordPoly::from_bigint_poly(f, fw);
    let g_mw = MultiwordPoly::from_bigint_poly(g, fw);
    // RNS regime only: transform f, g once (only k changes per reduction pass).
    let (f_ntt, g_ntt) = match &ntt {
        Some(ctx) => (Some(ctx.forward(&f_mw)), Some(ctx.forward(&g_mw))),
        None => (None, None),
    };
    // `k·f` / `k·g` per pass: pre-transformed RNS multiply for large `n`, flat-word
    // schoolbook for tiny `n` (where the RNS overhead dominates).
    let mul_kf = |kp: &MultiwordPoly| -> MultiwordPoly {
        match &ntt {
            Some(ctx) => ctx.mul_transformed(kp, f_ntt.as_ref().unwrap(), cap_w),
            None => kp.schoolbook_negacyclic_mul(&f_mw, cap_w),
        }
    };
    let mul_kg = |kp: &MultiwordPoly| -> MultiwordPoly {
        match &ntt {
            Some(ctx) => ctx.mul_transformed(kp, g_ntt.as_ref().unwrap(), cap_w),
            None => kp.schoolbook_negacyclic_mul(&g_mw, cap_w),
        }
    };
    let mut cf = MultiwordPoly::from_bigint_poly(capital_f, cap_w);
    let mut cg = MultiwordPoly::from_bigint_poly(capital_g, cap_w);

    let window = |cs: &MultiwordPoly, sh: u32| -> Polynomial<Complex64> {
        Polynomial::new(
            (0..n)
                .map(|i| Complex64::new(mw::shr_to_i128(cs.coeff(i), sh) as f64, 0.0))
                .collect::<Vec<_>>(),
        )
        .fft()
    };
    let cap_size = |cf: &MultiwordPoly, cg: &MultiwordPoly| -> u64 {
        (0..n)
            .flat_map(|i| [mw::bit_length(cf.coeff(i)), mw::bit_length(cg.coeff(i))])
            .fold(53, u64::max)
    };
    // k (i128 per coefficient) -> width-2 flat-word polynomial.
    let k_poly = |k: &[i128]| -> MultiwordPoly {
        let mut kp = MultiwordPoly::zeros(n, 2);
        for (i, &v) in k.iter().enumerate() {
            kp.coeff_mut(i)
                .copy_from_slice(MultiwordInt::from_i128(v, 2).limbs());
        }
        kp
    };
    let estimate_k = |cf: &MultiwordPoly, cg: &MultiwordPoly, sh: u32| -> Vec<i128> {
        let numerator = window(cf, sh).hadamard_mul(&f_star) + window(cg, sh).hadamard_mul(&g_star);
        let quotient = numerator.hadamard_div(&denominator_fft).ifft();
        quotient
            .coefficients
            .iter()
            .map(|c| c.re.round() as i128)
            .collect()
    };

    let mut prev_capital_size = u64::MAX;
    loop {
        let capital_size = cap_size(&cf, &cg);
        if capital_size < size || capital_size >= prev_capital_size {
            break;
        }
        prev_capital_size = capital_size;

        let d = capital_size - size;
        let (capital_shift, back_shift) = babai_shifts(capital_size, size, d);

        let k = estimate_k(&cf, &cg, capital_shift);
        if k.iter().all(|&x| x == 0) {
            break;
        }
        let kp = k_poly(&k);
        let shifted_kf = mul_kf(&kp).shl_coeffs(back_shift);
        let shifted_kg = mul_kg(&kp).shl_coeffs(back_shift);

        if d > 53 {
            let new_cs = cap_size(&cf.sub(&shifted_kf), &cg.sub(&shifted_kg));
            if new_cs >= capital_size {
                // Overshoot: fall back to the single-bit (d) formula this pass.
                let k_old = estimate_k(&cf, &cg, (capital_size - 53) as u32);
                if k_old.iter().all(|&x| x == 0) {
                    break;
                }
                let kpo = k_poly(&k_old);
                let kf_old = mul_kf(&kpo).shl_coeffs(d as u32);
                let kg_old = mul_kg(&kpo).shl_coeffs(d as u32);
                cf.sub_assign(&kf_old);
                cg.sub_assign(&kg_old);
                continue;
            }
        }

        cf.sub_assign(&shifted_kf);
        cg.sub_assign(&shifted_kg);
    }

    for (c, i) in capital_f.coefficients.iter_mut().zip(0..n) {
        *c = mw::to_bigint(cf.coeff(i));
    }
    for (c, i) in capital_g.coefficients.iter_mut().zip(0..n) {
        *c = mw::to_bigint(cg.coeff(i));
    }
    Ok(())
}

/// Babai reduction with `BigInt` (uncapped) capital and an RNS/NTT `k·f`
/// product.  Thin wrapper over [`babai_reduce_rns_generic`]; see it for the
/// algorithm.
pub(crate) fn babai_reduce_rns_bigint<const K: usize, P: NttPrimeList<K>>(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    babai_reduce_rns_generic::<BigInt, K, P>(f, g, capital_f, capital_g)
}

/// Reduce the vector (F,G) relative to (f,g). This method follows the python
/// implementation [1] but uses multimodular arithmetic for fast operations on
/// big integer polynomials in a cyclotomic ring.
///
/// Algorithm 7 in the spec [2, p.35]
///
/// [1]: https://github.com/tprest/falcon.py
///
/// [2]: https://falcon-sign.info/falcon.pdf
///
///
/// This function is marked pub for the purpose of benchmarking; it is not
/// considered part of the public API.
#[doc(hidden)]
#[profiling]
pub fn babai_reduce_i32(
    f: &Polynomial<i32>,
    g: &Polynomial<i32>,
    capital_f: &mut Polynomial<i32>,
    capital_g: &mut Polynomial<i32>,
) -> Result<(), String> {
    let f_ntt: Polynomial<U32Field> = f.map(|&i| U32Field::from(i)).fft();
    let g_ntt: Polynomial<U32Field> = g.map(|&i| U32Field::from(i)).fft();

    let bitsize = |itr: IntoIter<i32>| {
        // All-zero input has no meaningful bit size; `ilog2(0)` would panic.
        match itr.map(|i| i.abs()).max().unwrap() {
            0 => 0,
            m => (m * 2).ilog2().next_multiple_of(8) as usize,
        }
    };
    let size = usize::max(
        bitsize(
            f.coefficients
                .iter()
                .chain(g.coefficients.iter())
                .cloned()
                .collect_vec()
                .into_iter(),
        ),
        53,
    );

    let shift = (size as i64) - 53;
    let f_adjusted = f
        .map(|i| Complex64::new(i64::from(*i >> shift) as f64, 0.0))
        .fft();
    let g_adjusted = g
        .map(|i| Complex64::new(i64::from(*i >> shift) as f64, 0.0))
        .fft();

    let f_star_adjusted = f_adjusted.map(|c| c.conj());
    let g_star_adjusted = g_adjusted.map(|c| c.conj());
    let denominator_fft =
        f_adjusted.hadamard_mul(&f_star_adjusted) + g_adjusted.hadamard_mul(&g_star_adjusted);

    let mut prev_capital_size = usize::MAX;
    loop {
        let capital_size = [
            bitsize(
                capital_f
                    .coefficients
                    .iter()
                    .chain(capital_g.coefficients.iter())
                    .copied()
                    .collect_vec()
                    .into_iter(),
            ),
            53,
        ]
        .into_iter()
        .max()
        .unwrap();

        if capital_size < size || capital_size >= prev_capital_size {
            break;
        }
        prev_capital_size = capital_size;

        let capital_shift = (capital_size as i64) - 53;
        let capital_f_adjusted = capital_f
            .map(|bi| Complex64::new(i64::from(*bi >> capital_shift) as f64, 0.0))
            .fft();
        let capital_g_adjusted = capital_g
            .map(|bi| Complex64::new(i64::from(*bi >> capital_shift) as f64, 0.0))
            .fft();

        let numerator = capital_f_adjusted.hadamard_mul(&f_star_adjusted)
            + capital_g_adjusted.hadamard_mul(&g_star_adjusted);
        let quotient = numerator.hadamard_div(&denominator_fft).ifft();

        let k_ntt = quotient.map(|f| U32Field::from(f.re.round() as i32)).fft();

        if k_ntt.is_zero() {
            break;
        }

        let kf_ntt = k_ntt.hadamard_mul(&f_ntt).ifft();
        let kg_ntt = k_ntt.hadamard_mul(&g_ntt).ifft();

        let kf = kf_ntt.map(|p| p.balanced_value() as i32);
        let kg = kg_ntt.map(|p| p.balanced_value() as i32);

        *capital_f -= kf;
        *capital_g -= kg;
    }
    Ok(())
}

/// Reduce the vector (F, G) relative to (f, g) using f64 FFT for the quotient
/// step and RNS NTT for polynomial multiplication.
///
/// This mirrors `babai_reduce_i32`: Complex64 FFT computes the rounded quotient
/// k, which is then applied via RNS NTT multiplication.  capital_F,G live in
/// RNS throughout; they are decoded to i128 (via Garner) once per iteration
/// only for the FFT input and convergence check.
///
/// Algorithm 7 in the spec [1, p.35].
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[profiling]
pub(crate) fn babai_reduce_rns<const K: usize, P: NttPrimeList<K>>(
    f: &Polynomial<i32>,
    g: &Polynomial<i32>,
    capital_f: &mut Vec<Rns<K, P>>,
    capital_g: &mut Vec<Rns<K, P>>,
) -> Result<(), String> {
    let n = f.coefficients.len();

    // Precompute the per-prime twiddle tables once; reused by every transform
    // (the one-off f/g forward NTTs and the per-iteration k transform).
    // Production prime lists return compile-time tables; see `NttPrimeList::ntt_tables`.
    let tables = P::ntt_tables(n);

    // Precompute RNS NTT of f,g for polynomial multiplication.
    let mut f_rns_ntt: Vec<Rns<K, P>> = f.coefficients.iter().map(|&i| Rns::from_i32(i)).collect();
    let mut g_rns_ntt: Vec<Rns<K, P>> = g.coefficients.iter().map(|&i| Rns::from_i32(i)).collect();
    ntt_inplace_cached::<K, P>(&mut f_rns_ntt, &tables);
    ntt_inplace_cached::<K, P>(&mut g_rns_ntt, &tables);

    // Precompute Complex64 FFT of f,g. f,g are i32, so values always fit in
    // f64 with no shifting required.
    let f_fft = f.map(|&i| Complex64::new(i as f64, 0.0)).fft();
    let g_fft = g.map(|&i| Complex64::new(i as f64, 0.0)).fft();
    let f_adj = f_fft.map(|c| c.conj());
    let g_adj = g_fft.map(|c| c.conj());
    let denom_fft = f_fft.hadamard_mul(&f_adj) + g_fft.hadamard_mul(&g_adj);

    // "size" for f,g: bit-width rounded to multiple of 8, at least 53.
    // Since f,g are i32, size == 53 always.
    let max_fg = f
        .coefficients
        .iter()
        .chain(g.coefficients.iter())
        .map(|&i| i.unsigned_abs())
        .max()
        .unwrap_or(1);
    let fg_bits = if max_fg == 0 {
        0u32
    } else {
        (max_fg.saturating_mul(2)).ilog2().next_multiple_of(8)
    };
    let size = fg_bits.max(53);

    let mut prev_capital_size = u32::MAX;
    loop {
        // Decode capital_F,G from RNS to i128. Needed both for the convergence
        // check and for the f64 FFT input.
        let cf_i128: Vec<i128> = capital_f.iter().map(|r| r.to_i128()).collect();
        let cg_i128: Vec<i128> = capital_g.iter().map(|r| r.to_i128()).collect();

        let max_cap: u128 = cf_i128
            .iter()
            .chain(cg_i128.iter())
            .map(|&v| v.unsigned_abs())
            .max()
            .unwrap_or(0);

        let cap_bits = if max_cap == 0 {
            0u32
        } else {
            (max_cap.saturating_mul(2)).ilog2().next_multiple_of(8)
        };
        let capital_size = cap_bits.max(53);

        if capital_size < size || capital_size >= prev_capital_size {
            break;
        }
        prev_capital_size = capital_size;

        // Shift capital into f64 range: values become ~2^53 at most.
        let capital_shift = capital_size - 53;

        let capital_f_fft = Polynomial::new(
            cf_i128
                .iter()
                .map(|&v| Complex64::new((v >> capital_shift) as f64, 0.0))
                .collect::<Vec<_>>(),
        )
        .fft();
        let capital_g_fft = Polynomial::new(
            cg_i128
                .iter()
                .map(|&v| Complex64::new((v >> capital_shift) as f64, 0.0))
                .collect::<Vec<_>>(),
        )
        .fft();

        let numerator = capital_f_fft.hadamard_mul(&f_adj) + capital_g_fft.hadamard_mul(&g_adj);
        let quotient = numerator.hadamard_div(&denom_fft).ifft();

        // k_adj = round(k_true / 2^capital_shift).
        let k: Vec<i64> = quotient
            .coefficients
            .iter()
            .map(|c| c.re.round() as i64)
            .collect();

        if k.iter().all(|&x| x == 0) {
            break;
        }

        // k_true = k_adj << capital_shift. Encode in RNS and apply via NTT.
        let mut k_rns_ntt: Vec<Rns<K, P>> = k
            .iter()
            .map(|&k_adj| Rns::<K, P>::from_i128((k_adj as i128) << capital_shift))
            .collect();
        ntt_inplace_cached::<K, P>(&mut k_rns_ntt, &tables);

        let mut kf: Vec<Rns<K, P>> = k_rns_ntt
            .iter()
            .zip(f_rns_ntt.iter())
            .map(|(&a, &b)| a * b)
            .collect();
        let mut kg: Vec<Rns<K, P>> = k_rns_ntt
            .iter()
            .zip(g_rns_ntt.iter())
            .map(|(&a, &b)| a * b)
            .collect();
        intt_inplace_cached::<K, P>(&mut kf, &tables);
        intt_inplace_cached::<K, P>(&mut kg, &tables);

        for i in 0..n {
            capital_f[i] -= kf[i];
            capital_g[i] -= kg[i];
        }
    }
    Ok(())
}

/// Extended Euclidean algorithm for computing the greatest common divisor (g) and
/// Bézout coefficients (u, v) for the relation
///
///  u a + v b = g .
///
/// Implementation adapted from Wikipedia [1].
///
/// [1]: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
#[profiling]
fn xgcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    let (mut old_r, mut r) = (a.clone(), b.clone());
    let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());
    let (mut old_t, mut t) = (BigInt::zero(), BigInt::one());

    while r != BigInt::zero() {
        let quotient = old_r.clone() / r.clone();
        (old_r, r) = (r.clone(), old_r.clone() - quotient.clone() * r);
        (old_s, s) = (s.clone(), old_s.clone() - quotient.clone() * s);
        (old_t, t) = (t.clone(), old_t.clone() - quotient * t);
    }

    (old_r, old_s, old_t)
}

fn try_bigint_poly_to_i32(p: &Polynomial<BigInt>) -> Option<Polynomial<i32>> {
    let coeffs: Option<Vec<i32>> = p
        .coefficients
        .iter()
        .map(|c| i32::try_from(c.clone()).ok())
        .collect();
    coeffs.map(Polynomial::new)
}

fn babai_rns_with_fallback<const K: usize, P: NttPrimeList<K>>(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    let to_rns_vec = |poly: &Polynomial<BigInt>| -> Option<Vec<Rns<K, P>>> {
        poly.coefficients
            .iter()
            .map(|c| i128::try_from(c.clone()).ok().map(Rns::<K, P>::from_i128))
            .collect()
    };
    if let (Some(f_i32), Some(g_i32)) = (try_bigint_poly_to_i32(f), try_bigint_poly_to_i32(g)) {
        if let (Some(mut cf_rns), Some(mut cg_rns)) = (to_rns_vec(capital_f), to_rns_vec(capital_g))
        {
            if babai_reduce_rns::<K, P>(&f_i32, &g_i32, &mut cf_rns, &mut cg_rns).is_ok() {
                for (c, r) in capital_f.coefficients.iter_mut().zip(cf_rns.iter()) {
                    *c = BigInt::from(r.to_i64());
                }
                for (c, r) in capital_g.coefficients.iter_mut().zip(cg_rns.iter()) {
                    *c = BigInt::from(r.to_i64());
                }
                return Ok(());
            }
        }
    }
    babai_reduce_bigint(f, g, capital_f, capital_g)
}

/// Run `babai_reduce_rns` at recursion depth 1 (n = 512, K = 2 primes).
///
/// `capital_f` and `capital_g` are passed and returned as `i128` slices.
/// The conversion to/from `Rns` is included in the timed region because it
/// is part of the hot path in `babai_rns_with_fallback`.
#[doc(hidden)]
pub fn babai_reduce_rns_depth1(
    f: &Polynomial<i32>,
    g: &Polynomial<i32>,
    capital_f: &mut Vec<i128>,
    capital_g: &mut Vec<i128>,
) -> Result<(), String> {
    let mut cf_rns: Vec<Rns<2, NttPrimes24Bit2>> =
        capital_f.iter().map(|&v| Rns::from_i128(v)).collect();
    let mut cg_rns: Vec<Rns<2, NttPrimes24Bit2>> =
        capital_g.iter().map(|&v| Rns::from_i128(v)).collect();
    babai_reduce_rns::<2, NttPrimes24Bit2>(f, g, &mut cf_rns, &mut cg_rns)?;
    for (out, r) in capital_f.iter_mut().zip(cf_rns.iter()) {
        *out = r.to_i128();
    }
    for (out, r) in capital_g.iter_mut().zip(cg_rns.iter()) {
        *out = r.to_i128();
    }
    Ok(())
}

/// Run `babai_reduce_rns` at recursion depth 2 (n = 256, K = 4 primes).
///
/// See [`babai_reduce_rns_depth1`] for calling convention.
#[doc(hidden)]
pub fn babai_reduce_rns_depth2(
    f: &Polynomial<i32>,
    g: &Polynomial<i32>,
    capital_f: &mut Vec<i128>,
    capital_g: &mut Vec<i128>,
) -> Result<(), String> {
    let mut cf_rns: Vec<Rns<4, NttPrimes24Bit4>> =
        capital_f.iter().map(|&v| Rns::from_i128(v)).collect();
    let mut cg_rns: Vec<Rns<4, NttPrimes24Bit4>> =
        capital_g.iter().map(|&v| Rns::from_i128(v)).collect();
    babai_reduce_rns::<4, NttPrimes24Bit4>(f, g, &mut cf_rns, &mut cg_rns)?;
    for (out, r) in capital_f.iter_mut().zip(cf_rns.iter()) {
        *out = r.to_i128();
    }
    for (out, r) in capital_g.iter_mut().zip(cg_rns.iter()) {
        *out = r.to_i128();
    }
    Ok(())
}

/// Babai reduction with **packed** fixed-width (`L`-limb) capital coefficients.
///
/// [`babai_reduce_rns_generic`] specialized to [`Packed<L>`]: the capital lives
/// in `64·L`-bit two's-complement limbs instead of `BigInt`, so the
/// per-iteration hot path is allocation-free and word-level:
///
/// * `capital_size` is a limb scan ([`Packed::bit_length`]);
/// * the FFT input is a top-word read ([`Packed::shr_to_i128`]);
/// * the reduction coefficient `k` stays `i128` (it never exceeds ~2^53), so
///   feeding it into RNS uses the cheap `Rns::from_i128`, not a `BigInt` CRT;
/// * the `k·f` product is reconstructed straight into limbs with word-level
///   CRT ([`Packed::from_rns`]), then shifted and subtracted in place.
///
/// `BigInt` is touched only once at entry and once at exit.  This is the
/// representation fn-dsa uses (base 2^31 there, base 2^64 here); it lets the
/// NTT multiply beat `BigInt` karatsuba at depth 3 (see the per-depth bench).
pub(crate) fn babai_reduce_rns_packed<const LP: usize, const K: usize, P: NttPrimeList<K>>(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    babai_reduce_rns_generic::<Packed<LP>, K, P>(f, g, capital_f, capital_g)
}

/// Run [`babai_reduce_rns_packed`] at recursion depth 3 (n = 128, K = 5 primes,
/// L = 4 limbs / 256-bit capital).
///
/// Only 5 primes are needed: the path pushes the `k·f` product (a hard ≈107-bit
/// quantity, see `ntt_primes24_5_covers_depth3_product`) — not the ≈154-bit
/// capital — through RNS.  The ≈154-bit capital fits in 4 limbs.
#[doc(hidden)]
pub fn babai_reduce_rns_packed_depth3(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    babai_reduce_rns_packed::<4, 5, NttPrimes24Bit5>(f, g, capital_f, capital_g)
}

/// Run [`babai_reduce_rns_packed`] at recursion depth 4 (n = 64, K = 8 primes,
/// L = 5 limbs / 320-bit capital).
///
/// At depth 4 the `k·f` product is ≈155 bits (`bits(f)≈101` plus the ≈53-bit
/// reduction coefficient `k`; see `product_size_model_matches_measurements`),
/// so the 8-prime ≈183-bit list is required — the 5-prime list (≈114 bits)
/// would wrap.  The ≈303-bit capital needs 5 limbs.
///
/// `#[doc(hidden)]` and NOT wired into [`ntru_solve`]: at depth 4 this is
/// **slower** than `BigInt` karatsuba (n=64 bench: ≈1.6 ms vs ≈1.0 ms).  As `n`
/// halves with depth, karatsuba on a small polynomial gets cheap while the RNS
/// per-iteration overhead (k→RNS conversion, the now 8-prime NTT, word-level
/// CRT reconstruction) and the extra reduction iterations dominate — the RNS
/// crossover sits at ~depth 3.  This entry point exists only so the per-depth
/// benchmark and equivalence test can measure that, and to keep the const-`L`
/// machinery exercised; it is not part of the production reduction path.
#[doc(hidden)]
pub fn babai_reduce_rns_packed_depth4(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    babai_reduce_rns_packed::<5, 8, NttPrimes24Bit8>(f, g, capital_f, capital_g)
}

/// Run [`babai_reduce_rns_bigint`] at recursion depth 3 (n = 128, K = 5 primes).
///
/// Capital coefficients at this depth exceed 127 bits, so they are carried as
/// `BigInt` rather than `i128`; the `k·f` product still runs in RNS via the
/// NTT.  Provided for the per-depth benchmark and the depth-3 correctness test.
#[doc(hidden)]
pub fn babai_reduce_rns_bigint_depth3(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    babai_reduce_rns_bigint::<5, NttPrimes24Bit5>(f, g, capital_f, capital_g)
}

/// Run [`babai_reduce_rns_bigint`] at recursion depth 4 (n = 64, K = 8 primes).
/// Capital carried as `BigInt`; the ≈155-bit `k·f` product needs the 8-prime
/// list.
///
/// `#[doc(hidden)]` and NOT wired into [`ntru_solve`]: like its packed sibling
/// [`babai_reduce_rns_packed_depth4`], this loses to `BigInt` karatsuba at
/// depth 4 — see that function's note for why the RNS crossover ends at
/// ~depth 3.  Provided only for the per-depth benchmark and the depth-4
/// equivalence test.
#[doc(hidden)]
pub fn babai_reduce_rns_bigint_depth4(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) -> Result<(), String> {
    babai_reduce_rns_bigint::<8, NttPrimes24Bit8>(f, g, capital_f, capital_g)
}

/// Per-depth parameters for the runtime-`K` flat-word babai reduction
/// ([`babai_reduce_rns_runtime`]): `(k_primes, cap_w)`.  `k_primes` covers the
/// `≈ bits(f)+54`-bit `k·f` product with a +6σ tail; `cap_w` (64-bit limbs)
/// holds the depth's capital (`max|F,G|`, +6σ) plus the mid-reduction shifts.
/// Derived from [`NTRU_SOLVE_BABAI_COEFF_BITS`].
const fn rns_runtime_babai_params(depth: usize) -> (usize, usize) {
    match depth {
        2 => (5, 4),
        3 => (6, 5),
        4 => (8, 7),
        5 => (13, 12),
        6 => (22, 22),
        7 => (41, 40),
        8 => (77, 77),
        9 => (149, 152),
        _ => (0, 0),
    }
}

/// Deepest recursion depth at which the allocation-free flat-word babai reduction
/// ([`babai_reduce_rns_runtime`]) is used instead of `babai_reduce_bigint`.
/// Depths 3..=this use the flat-word path; deeper levels fall back to BigInt.
/// Overridable at runtime via the `RNS_RUNTIME_MAX_DEPTH` env var.
///
/// Default 9 (all depths): the flat-word babai now dispatches its `k·f` multiply
/// internally — RNS/NTT at large `n` (depth 3, n=128) and flat-word schoolbook at
/// tiny `n` (depths 4–9, n≤64), the latter skipping the expensive `RuntimeNtt`
/// entirely.  This replaces the old RNS-only crossover (which had to fall back to
/// BigInt at depth 5 because the RNS multiply lost at small `n`).
fn rns_runtime_max_depth() -> usize {
    use std::sync::LazyLock;
    static D: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("RNS_RUNTIME_MAX_DEPTH")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(9)
    });
    *D
}

/// Field norm via the allocation-free flat-word [`MultiwordPoly`] path instead of
/// `Polynomial<BigInt>::field_norm`.  The two squarings inside dominate the deep
/// NTRU-solve recursion; routing them through the size-dispatched negacyclic
/// multiply (RNS/NTT at shallow large `n`, flat-word schoolbook at deep tiny `n`)
/// removes the per-`BigInt`-operation heap allocation.  Bit-exact with the
/// `BigInt` field norm (see `field_norm_matches_bigint`).
#[profiling]
fn field_norm_flatword(p: &Polynomial<BigInt>) -> Polynomial<BigInt> {
    let max_bits = p.coefficients.iter().map(|c| c.bits()).max().unwrap_or(0);
    // ceil(max_bits/64) limbs for the magnitude, +1 for the sign bit / safety.
    let w_in = (max_bits / 64 + 1) as usize;
    MultiwordPoly::from_bigint_poly(p, w_in)
        .field_norm()
        .to_bigint_poly()
}

/// Negacyclic product `a * b mod (X^n + 1)` via the allocation-free flat-word
/// multiply, replacing `a.karatsuba(&b).reduce_by_cyclotomic(n)` on
/// `Polynomial<BigInt>`.  Used for the two F,G construction multiplies in
/// `ntru_solve` (`capital_*_prime_xsq × *_minx`), which are asymmetric — a large
/// lifted capital times a small galois-adjoint operand.  Bit-exact with
/// karatsuba-then-reduce (`negacyclic_mul_dispatch_matches_bigint`).
#[profiling]
fn construction_mul_flatword(a: &Polynomial<BigInt>, b: &Polynomial<BigInt>) -> Polynomial<BigInt> {
    let n = a.coefficients.len();
    let bits_a = a.coefficients.iter().map(|c| c.bits()).max().unwrap_or(0);
    let bits_b = b.coefficients.iter().map(|c| c.bits()).max().unwrap_or(0);
    // ceil(bits/64) magnitude limbs + 1 sign/safety limb per operand.
    let wa = (bits_a / 64 + 1) as usize;
    let wb = (bits_b / 64 + 1) as usize;
    // |product coeff| < n · 2^(bits_a+bits_b); +1 sign limb.
    let prod_bits = bits_a + bits_b + (n as u64).max(1).ilog2() as u64 + 2;
    let out_w = (prod_bits / 64 + 1) as usize;

    let ma = MultiwordPoly::from_bigint_poly(a, wa);
    let mb = MultiwordPoly::from_bigint_poly(b, wb);
    ma.negacyclic_mul(&mb, out_w).to_bigint_poly()
}

/// Solve the NTRU equation. Given f, g in ZZ[X], find F, G in ZZ[X].
/// such that
///
///    f G - g F = q  mod (X^n + 1)
///
/// Algorithm 6 of the specification [1, p.35].
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[profiling]
fn ntru_solve(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    depth: usize,
    max_rns_depth: usize,
) -> Option<(Polynomial<BigInt>, Polynomial<BigInt>)> {
    let n = f.coefficients.len();
    if n == 1 {
        let (gcd, u, v) = xgcd(&f.coefficients[0], &g.coefficients[0]);
        if gcd != BigInt::one() {
            return None;
        }
        return Some((
            (Polynomial::new(vec![-v * BigInt::from_u32(Q as u32).unwrap()])),
            Polynomial::new(vec![u * BigInt::from_u32(Q as u32).unwrap()]),
        ));
    }

    let f_prime = field_norm_flatword(f);
    let g_prime = field_norm_flatword(g);
    let (capital_f_prime, capital_g_prime) =
        ntru_solve(&f_prime, &g_prime, depth + 1, max_rns_depth)?;

    let capital_f_prime_xsq = capital_f_prime.lift_next_cyclotomic();
    let capital_g_prime_xsq = capital_g_prime.lift_next_cyclotomic();
    let f_minx = f.galois_adjoint();
    let g_minx = g.galois_adjoint();

    let mut capital_f = construction_mul_flatword(&capital_f_prime_xsq, &g_minx);
    let mut capital_g = construction_mul_flatword(&capital_g_prime_xsq, &f_minx);

    let babai_result = match depth {
        1 if max_rns_depth >= 1 => {
            babai_rns_with_fallback::<2, NttPrimes24Bit2>(f, g, &mut capital_f, &mut capital_g)
        }
        2 if max_rns_depth >= 2 => {
            babai_rns_with_fallback::<4, NttPrimes24Bit4>(f, g, &mut capital_f, &mut capital_g)
        }
        // Depth 2 (production, where max_rns_depth=1 so the i128 arm above is
        // skipped) and the deep levels: allocation-free flat-word reduction. At
        // depth 2 (n=256 > schoolbook cutoff) this is the transform-once RNS
        // multiply; deeper levels dispatch to the flat-word schoolbook internally.
        d if (2..=rns_runtime_max_depth()).contains(&d) => {
            let (k_primes, cap_w) = rns_runtime_babai_params(d);
            babai_reduce_rns_runtime(f, g, &mut capital_f, &mut capital_g, k_primes, cap_w)
        }
        _ => babai_reduce_bigint(f, g, &mut capital_f, &mut capital_g),
    };
    match babai_result {
        Ok(_) => Some((capital_f, capital_g)),
        Err(_e) => {
            #[cfg(test)]
            {
                panic!("{}", _e);
            }
            #[cfg(not(test))]
            {
                None
            }
        }
    }
}

#[profiling]
fn ntru_solve_entrypoint(
    f: Polynomial<i32>,
    g: Polynomial<i32>,
    max_rns_depth: usize,
) -> Option<(Polynomial<i32>, Polynomial<i32>)> {
    let g_prime = g.field_norm().map(|c| BigInt::from(*c));
    let f_prime = f.field_norm().map(|c| BigInt::from(*c));
    let (capital_f_prime_bi, capital_g_prime_bi) =
        ntru_solve(&f_prime, &g_prime, 1, max_rns_depth)?;

    let capital_f_prime_coefficients = capital_f_prime_bi
        .coefficients
        .into_iter()
        .map(i32::try_from)
        .collect_vec();
    let capital_g_prime_coefficients = capital_g_prime_bi
        .coefficients
        .into_iter()
        .map(i32::try_from)
        .collect_vec();

    if !capital_f_prime_coefficients
        .iter()
        .chain(capital_g_prime_coefficients.iter())
        .all(|c| c.is_ok())
    {
        return None;
    }
    let capital_f_prime = Polynomial::new(
        capital_f_prime_coefficients
            .into_iter()
            .map(|c| c.unwrap())
            .collect_vec(),
    );
    let capital_g_prime = Polynomial::new(
        capital_g_prime_coefficients
            .into_iter()
            .map(|c| c.unwrap())
            .collect_vec(),
    );

    let capital_f_prime_xsq = capital_f_prime.lift_next_cyclotomic();
    let capital_g_prime_xsq = capital_g_prime.lift_next_cyclotomic();
    let f_minx = f.galois_adjoint();
    let g_minx = g.galois_adjoint();

    // Multiply via the single-prime NTT, reusing the compile-time twiddle
    // tables through `FastFft` (no per-call `bitreversed_powers` recompute).
    let cfp_ntt = capital_f_prime_xsq.map(|c| U32Field::from(*c)).fft();
    let cgp_ntt = capital_g_prime_xsq.map(|c| U32Field::from(*c)).fft();
    let gm_ntt = g_minx.map(|c| U32Field::from(*c)).fft();
    let fm_ntt = f_minx.map(|c| U32Field::from(*c)).fft();
    let cf_ntt = cfp_ntt.hadamard_mul(&gm_ntt).ifft();
    let cg_ntt = cgp_ntt.hadamard_mul(&fm_ntt).ifft();

    let mut capital_f = cf_ntt.map(|c| c.balanced_value() as i32);
    let mut capital_g = cg_ntt.map(|c| c.balanced_value() as i32);

    match babai_reduce_i32(&f, &g, &mut capital_f, &mut capital_g) {
        Ok(_) => Some((capital_f, capital_g)),
        Err(_e) => {
            #[cfg(test)]
            {
                panic!("{}", _e);
            }
            #[cfg(not(test))]
            {
                None
            }
        }
    }
}

/// Sample 4 small polynomials f, g, F, G such that f * G - g * F = q mod (X^n + 1).
/// Algorithm 5 (NTRUgen) of the documentation [1, p.34].
///
/// This function is marked pub for benchmarking purposes only. Not considered part
/// of the public API.
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[doc(hidden)]
#[profiling]
pub fn ntru_gen<R: Rng + ?Sized>(
    n: usize,
    rng: &mut R,
) -> (
    Polynomial<i16>,
    Polynomial<i16>,
    Polynomial<i16>,
    Polynomial<i16>,
) {
    // let mut rng: StdRng = SeedableRng::from_seed(seed);

    loop {
        let f = gen_poly(n, rng);
        let g = gen_poly(n, rng);

        let f_ntt = f.map(|&i| Felt::from(i)).fft();
        if f_ntt.coefficients.iter().any(|e| e.is_zero()) {
            continue;
        }
        let gamma = gram_schmidt_norm_squared(&f, &g);
        if gamma > 1.3689f64 * (Q as f64) {
            continue;
        }

        // max_rns_depth=1: use the i128 RNS babai only at depth 1 (n=512, K=2),
        // where it beats BigInt.  Depth 2 (K=4) is a measured net loss, so it
        // stays on BigInt.  The deeper allocation-free runtime-K RNS path is
        // gated separately by `rns_runtime_max_depth()`, not by this argument.
        if let Some((capital_f, capital_g)) =
            ntru_solve_entrypoint(f.map(|&i| i as i32), g.map(|&i| i as i32), 1)
        {
            // Verify the NTRU equation fG − gF = q.  The depth-0 reduction NTT
            // uses a single 24-bit prime whose ~2^22 reconstruction range only
            // narrowly covers the products; this check turns any rare overflow
            // into a retry instead of a silently invalid key.  Computed in i64
            // so the unreduced products cannot overflow.
            let f_i64 = f.map(|&i| i as i64);
            let g_i64 = g.map(|&i| i as i64);
            let cf_i64 = capital_f.map(|&i| i as i64);
            let cg_i64 = capital_g.map(|&i| i as i64);
            let ntru_eq =
                (f_i64 * cg_i64).reduce_by_cyclotomic(n) - (g_i64 * cf_i64).reduce_by_cyclotomic(n);
            if ntru_eq != Polynomial::constant(Q as i64) {
                continue;
            }
            return (
                f,
                g,
                capital_f.map(|&i| i as i16),
                capital_g.map(|&i| i as i16),
            );
        }
    }
}

/// Like [`ntru_gen`] but with an explicit RNS depth threshold for benchmarking.
///
/// `max_rns_depth` controls how deep into the NTRU recursion `babai_reduce_rns`
/// is used instead of `babai_reduce_bigint`:
/// - 0: always use BigInt (baseline)
/// - 1: RNS at depth 1 only (n=512, K=2 primes)
/// - 2: RNS at depths 1–2 (n=512 and n=256, K=2 and K=4 primes)
///
/// This function is marked pub for benchmarking purposes only.
#[doc(hidden)]
#[profiling]
pub fn ntru_gen_with_rns_depth<R: Rng + ?Sized>(
    n: usize,
    rng: &mut R,
    max_rns_depth: usize,
) -> (
    Polynomial<i16>,
    Polynomial<i16>,
    Polynomial<i16>,
    Polynomial<i16>,
) {
    loop {
        let f = gen_poly(n, rng);
        let g = gen_poly(n, rng);

        let f_ntt = f.map(|&i| Felt::from(i)).fft();
        if f_ntt.coefficients.iter().any(|e| e.is_zero()) {
            continue;
        }
        let gamma = gram_schmidt_norm_squared(&f, &g);
        if gamma > 1.3689f64 * (Q as f64) {
            continue;
        }

        if let Some((capital_f, capital_g)) =
            ntru_solve_entrypoint(f.map(|&i| i as i32), g.map(|&i| i as i32), max_rns_depth)
        {
            return (
                f,
                g,
                capital_f.map(|&i| i as i16),
                capital_g.map(|&i| i as i16),
            );
        }
    }
}

/// Generate a polynomial of degree at most n-1 whose coefficients are
/// distributed according to a discrete Gaussian with mu = 0 and
/// sigma = 1.17 * sqrt(Q / (2n)).
// fn gen_poly(n: usize, rng: &mut dyn Rng) -> Polynomial<i16> {
#[profiling]
fn gen_poly<R: Rng + ?Sized>(n: usize, rng: &mut R) -> Polynomial<i16> {
    let mu = FixedPoint128::ZERO;
    let sigma_star = FixedPoint128::from(1.43300980528773f64);
    // Build the sigma-dependent sampler constants once, not per sample.
    let ctx = SamplerZCtx::new(sigma_star, sigma_star - FixedPoint128::from(0.001f64));
    const NUM_COEFFICIENTS: usize = 4096;
    Polynomial {
        coefficients: (0..NUM_COEFFICIENTS)
            .map(|_| ctx.sample(mu, rng))
            .collect_vec()
            .chunks(NUM_COEFFICIENTS / n)
            .map(|ch| ch.iter().sum())
            .collect_vec(),
    }
}

/// Compute the Gram-Schmidt norm of B = [[g, -f], [G, -F]] from f and g.
/// Corresponds to line 9 in algorithm 5 of the spec [1, p.34]
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[profiling]
fn gram_schmidt_norm_squared(f: &Polynomial<i16>, g: &Polynomial<i16>) -> f64 {
    let n = f.coefficients.len();
    let gamma1 = f64::from(f.l2_norm_squared() + g.l2_norm_squared());

    let fp = |i: &i16| Complex64::new(*i as f64, 0.0);
    let q_fp = Q as f64;
    let n_fp = n as f64;

    let f_fft = f.map(fp).fft();
    let g_fft = g.map(fp).fft();
    let f_adj_fft = f_fft.map(|c| c.conj());
    let g_adj_fft = g_fft.map(|c| c.conj());
    let ffgg_fft = f_fft.hadamard_mul(&f_adj_fft) + g_fft.hadamard_mul(&g_adj_fft);
    let ffgg_fft_inverse = ffgg_fft.hadamard_inv();
    let qf_over_ffgg_fft = f_adj_fft.map(|c| c * q_fp).hadamard_mul(&ffgg_fft_inverse);
    let qg_over_ffgg_fft = g_adj_fft.map(|c| c * q_fp).hadamard_mul(&ffgg_fft_inverse);
    let norm_f_over_ffgg_squared = qf_over_ffgg_fft
        .coefficients
        .iter()
        .map(|c| (c * c.conj()).re)
        .sum::<f64>()
        / n_fp;
    let norm_g_over_ffgg_squared = qg_over_ffgg_fft
        .coefficients
        .iter()
        .map(|c| (c * c.conj()).re)
        .sum::<f64>()
        / n_fp;

    let gamma2 = norm_f_over_ffgg_squared + norm_g_over_ffgg_squared;

    f64::max(gamma1, gamma2)
}

#[cfg(test)]
mod test {

    use std::str::FromStr;

    use itertools::Itertools;
    use num::BigInt;
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest::strategy::Just;
    use rand::{rngs::StdRng, SeedableRng};
    use test_strategy::proptest as strategy_proptest;

    use crate::math::{
        babai_reduce_i32, babai_reduce_rns, gram_schmidt_norm_squared, ntru_gen, ntru_solve,
    };
    use crate::{
        fp_field::FpField,
        polynomial::Polynomial,
        rns::{NttPrimeList, NttPrimes24Bit2, NttPrimes24Bit4, PrimeList, Rns},
    };

    use super::babai_reduce_bigint;

    // Prime list for babai_reduce_rns tests: three primes ≡ 1 (mod 2048).
    struct NttPrimes3;
    impl PrimeList<3> for NttPrimes3 {
        const PRIMES: [u32; 3] = [786_433, 998_244_353, 1_073_754_113];
    }
    impl NttPrimeList<3> for NttPrimes3 {
        const ROOTS_OF_UNITY_2048: [u32; 3] = [
            FpField::<786_433>::primitive_nth_root_of_unity(2048).value(),
            FpField::<998_244_353>::primitive_nth_root_of_unity(2048).value(),
            FpField::<1_073_754_113>::primitive_nth_root_of_unity(2048).value(),
        ];
    }
    type Rns3 = Rns<3, NttPrimes3>;

    #[test]
    fn ntt_primes3_covers_babai_depths_0_and_1() {
        // Signed bit-capacity of the NttPrimes3 configuration.
        let capacity: f64 = NttPrimes3::PRIMES
            .iter()
            .map(|&p| f64::log2(p as f64))
            .sum::<f64>()
            - 1.0;

        // Required bits at depth d: avg + 6·σ + 2 (sign bit + one factor-of-2
        // slack because F − k·f can transiently reach 2·|F| before shrinking).
        let required = |d: usize| {
            let (_, _, avg_f, std_f) = super::NTRU_SOLVE_BABAI_COEFF_BITS[d];
            avg_f + 6.0 * std_f + 2.0
        };

        assert!(
            capacity >= required(0),
            "depth 0: {capacity:.1} < {:.1}",
            required(0)
        );
        assert!(
            capacity >= required(1),
            "depth 1: {capacity:.1} < {:.1}",
            required(1)
        );
        // Depth 2 must exceed our capacity — if this assertion ever fails, the
        // prime list is large enough to extend babai_reduce_rns to depth 2.
        assert!(
            capacity < required(2),
            "depth 2: {capacity:.1} >= {:.1}",
            required(2)
        );
    }

    #[test]
    fn ntt_primes24_2_covers_depth1() {
        let capacity: f64 = NttPrimes24Bit2::PRIMES
            .iter()
            .map(|&p| f64::log2(p as f64))
            .sum::<f64>()
            - 1.0;
        let required = |d: usize| {
            let (_, _, avg_f, std_f) = super::NTRU_SOLVE_BABAI_COEFF_BITS[d];
            avg_f + 6.0 * std_f + 2.0
        };
        assert!(
            capacity >= required(1),
            "depth 1: {capacity:.1} < {:.1}",
            required(1)
        );
        assert!(
            capacity < required(2),
            "depth 2: {capacity:.1} >= {:.1}",
            required(2)
        );
    }

    #[test]
    fn ntt_primes24_4_covers_depth2() {
        let capacity: f64 = NttPrimes24Bit4::PRIMES
            .iter()
            .map(|&p| f64::log2(p as f64))
            .sum::<f64>()
            - 1.0;
        let required = |d: usize| {
            let (_, _, avg_f, std_f) = super::NTRU_SOLVE_BABAI_COEFF_BITS[d];
            avg_f + 6.0 * std_f + 2.0
        };
        assert!(
            capacity >= required(2),
            "depth 2: {capacity:.1} < {:.1}",
            required(2)
        );
        assert!(
            capacity < required(3),
            "depth 3: {capacity:.1} >= {:.1}",
            required(3)
        );
    }

    fn babai_infinite_loop_polynomials() -> (
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
    ) {
        let f = Polynomial::new(
            [
                BigInt::from_str("6426042728002").unwrap(),
                BigInt::from_str("-20675284604736").unwrap(),
                BigInt::from_str("-12121913318466").unwrap(),
                BigInt::from_str("-27836101162563").unwrap(),
            ]
            .to_vec(),
        );

        let g = Polynomial::new(
            [
                BigInt::from_str("-1001246212").unwrap(),
                BigInt::from_str("-1347303037").unwrap(),
                BigInt::from_str("987026048").unwrap(),
                BigInt::from_str("-1001311747").unwrap(),
            ]
            .to_vec(),
        );

        let capital_f = Polynomial::new(
            [
                BigInt::from_str(
                    "563985131491945032326798334533872091781886676547754689048287010878681928",
                )
                .unwrap(),
                BigInt::from_str(
                    "-348444005402208553421931883447687919671423051554816023996113866522386058",
                )
                .unwrap(),
                BigInt::from_str(
                    "-85657170778585026649528684432821341936755757853602491207147473952485632",
                )
                .unwrap(),
                BigInt::from_str(
                    "135623655239747178410899900677875843487151183900794566193191499131611018",
                )
                .unwrap(),
            ]
            .to_vec(),
        );

        let capital_g = Polynomial::new(
            [
                BigInt::from_str(
                    "49040356584788663746447138446729467702643846166576265941049418069366",
                )
                .unwrap(),
                BigInt::from_str(
                    "-57075549200927059197269512430308877512934179841045274854176350745681",
                )
                .unwrap(),
                BigInt::from_str(
                    "18442173959410247991253446345066800376513088376845717824090327663990",
                )
                .unwrap(),
                BigInt::from_str(
                    "19528334302175388221061434098432127604592213845277598673231565264960",
                )
                .unwrap(),
            ]
            .to_vec(),
        );

        (f, g, capital_f, capital_g)
    }

    #[test]
    fn babai_oscillation_terminates() {
        let (f, g, mut capital_f, mut capital_g) = babai_infinite_loop_polynomials();
        let _ = babai_reduce_bigint(&f, &g, &mut capital_f, &mut capital_g);
    }

    /// Generate a random signed `BigInt` with up to `bits` magnitude bits.
    fn rand_signed_bigint(rng: &mut StdRng, bits: u32) -> BigInt {
        use rand::RngExt;
        let mut v = BigInt::from(0);
        for _ in 0..bits {
            v = (v << 1) | BigInt::from(rng.random::<bool>() as u8);
        }
        if rng.random::<bool>() {
            -v
        } else {
            v
        }
    }

    /// Realistic depth-2 inputs: n = 256, max|f,g| ≈ 24 bits, max|F,G| ≈ 78
    /// bits (from `NTRU_SOLVE_BABAI_COEFF_BITS`).  Exercises the ≈80-bit `k·f`
    /// product through the runtime path (K = 5 primes, 4-limb capital) — the
    /// CRT-wrap detector for the depth-2 params now that production routes
    /// depth 2 through `babai_reduce_rns_runtime`.
    fn depth2_inputs(
        rng: &mut StdRng,
    ) -> (
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
    ) {
        let n = 256;
        let f = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 24))
                .collect::<Vec<_>>(),
        );
        let g = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 24))
                .collect::<Vec<_>>(),
        );
        let cf = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 78))
                .collect::<Vec<_>>(),
        );
        let cg = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 78))
                .collect::<Vec<_>>(),
        );
        (f, g, cf, cg)
    }

    /// Realistic depth-3 inputs: n = 128, max|f,g| ≈ 50 bits, max|F,G| ≈ 154
    /// bits (the actual Falcon-1024 depth-3 magnitudes, from
    /// `NTRU_SOLVE_BABAI_COEFF_BITS`).  These sizes exercise the full ≈107-bit
    /// `k·f` product, so the equivalence tests double as a CRT-wrap detector
    /// for the right-sized (K = 5) prime list.
    fn depth3_inputs(
        rng: &mut StdRng,
    ) -> (
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
    ) {
        let n = 128;
        let f = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 50))
                .collect::<Vec<_>>(),
        );
        let g = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 50))
                .collect::<Vec<_>>(),
        );
        let cf = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 154))
                .collect::<Vec<_>>(),
        );
        let cg = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 154))
                .collect::<Vec<_>>(),
        );
        (f, g, cf, cg)
    }

    /// Realistic depth-4 inputs: n = 64, max|f,g| ≈ 101 bits, max|F,G| ≈ 303
    /// bits.  Exercises the full ≈155-bit `k·f` product (CRT-wrap detector for
    /// the K = 8 prime list) and the 5-limb (320-bit) capital.
    fn depth4_inputs(
        rng: &mut StdRng,
    ) -> (
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
    ) {
        let n = 64;
        let f = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 101))
                .collect::<Vec<_>>(),
        );
        let g = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 101))
                .collect::<Vec<_>>(),
        );
        let cf = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 303))
                .collect::<Vec<_>>(),
        );
        let cg = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 303))
                .collect::<Vec<_>>(),
        );
        (f, g, cf, cg)
    }

    /// Sizing model for the multiword RNS reduction paths, established by
    /// measuring the max `k·f` product bit-width on realistic inputs (with the
    /// 183-bit `NttPrimes24Bit8` so the product never wraps during measurement):
    ///
    ///   depth 3 (n=128, bits(f)≈50):  product ≈ 107 bits,  k ≈ 54–55 bits
    ///   depth 4 (n=64,  bits(f)≈101): product ≈ 155 bits,  k ≈ 51–52 bits
    ///
    /// The product does NOT self-normalize: the reduction coefficient `k` stays
    /// ≈53 bits at every depth (it is estimated from windowed top words), but
    /// the product is `k · f_full` with the *full* `f`, so it grows with depth
    /// as `bits(f) + ~54`.  Crucially this is still far below the capital
    /// (`bits(F) ≈ 3·bits(f)`): the modulus need only cover the product, so the
    /// prime count is roughly halved versus sizing for the capital.  Sizing for
    /// the +6σ tail, `ceil((avg_f + 6·σ_f + 54 + 2) / 23)` primes: 5 at depth 3,
    /// 8 at depth 4 (the +6σ tail pushes the product to ≈162 bits, past K=7's
    /// ≈160-bit capacity).
    #[test]
    fn product_size_model_matches_measurements() {
        // (bits(f), measured max product) data points.
        for &(f_bits, measured) in &[(50.0_f64, 107.0_f64), (101.0, 155.0)] {
            let model = f_bits + 54.0;
            assert!(
                (model - measured).abs() <= 8.0,
                "product-size model {model:.0} far from measured {measured:.0} (bits(f)={f_bits})"
            );
        }
    }

    /// The RNS-NTT multiply backend must produce exactly the same reduction as
    /// the BigInt-karatsuba backend, on realistic depth-3-sized inputs (capital
    /// coefficients > 127 bits, exercising the `BigInt` capital path that the
    /// i128-bound `babai_reduce_rns` cannot represent).
    #[test]
    fn babai_reduce_rns_bigint_depth3_matches_bigint() {
        use super::babai_reduce_rns_bigint_depth3;

        let mut rng = StdRng::seed_from_u64(0xba_ba_13_d3);
        for _ in 0..32 {
            let (f, g, cap_f, cap_g) = depth3_inputs(&mut rng);

            let (mut bf, mut bg) = (cap_f.clone(), cap_g.clone());
            babai_reduce_bigint(&f, &g, &mut bf, &mut bg).unwrap();

            let (mut rf, mut rg) = (cap_f, cap_g);
            babai_reduce_rns_bigint_depth3(&f, &g, &mut rf, &mut rg).unwrap();

            assert_eq!(bf, rf, "capital_F mismatch between bigint and rns-bigint");
            assert_eq!(bg, rg, "capital_G mismatch between bigint and rns-bigint");
        }
    }

    /// The packed-limb backend must produce exactly the same reduction as the
    /// BigInt backend (same realistic depth-3-sized inputs).
    #[test]
    fn babai_reduce_rns_packed_depth3_matches_bigint() {
        use super::babai_reduce_rns_packed_depth3;

        let mut rng = StdRng::seed_from_u64(0x9ac_ed_d3);
        for _ in 0..32 {
            let (f, g, cap_f, cap_g) = depth3_inputs(&mut rng);

            let (mut bf, mut bg) = (cap_f.clone(), cap_g.clone());
            babai_reduce_bigint(&f, &g, &mut bf, &mut bg).unwrap();

            let (mut pf, mut pg) = (cap_f, cap_g);
            babai_reduce_rns_packed_depth3(&f, &g, &mut pf, &mut pg).unwrap();

            assert_eq!(bf, pf, "capital_F mismatch between bigint and packed");
            assert_eq!(bg, pg, "capital_G mismatch between bigint and packed");
        }
    }

    /// Realistic depth-5 inputs: n = 32, max|f,g| ≈ 202 bits, max|F,G| ≈ 600
    /// bits.  Exercises a ≈256-bit `k·f` product that needs ~12 runtime primes.
    fn depth5_inputs(
        rng: &mut StdRng,
    ) -> (
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
    ) {
        let n = 32;
        let f = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 202))
                .collect::<Vec<_>>(),
        );
        let g = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 202))
                .collect::<Vec<_>>(),
        );
        let cf = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 600))
                .collect::<Vec<_>>(),
        );
        let cg = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 600))
                .collect::<Vec<_>>(),
        );
        (f, g, cf, cg)
    }

    /// Realistic depth-7 inputs: n = 8, max|f,g| ≈ 850 bits, max|F,G| ≈ 2000
    /// bits.  Deep in the schoolbook regime (n ≤ 64), where the flat-word babai
    /// uses `schoolbook_negacyclic_mul` and skips the `RuntimeNtt` entirely.
    fn depth7_inputs(
        rng: &mut StdRng,
    ) -> (
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
        Polynomial<BigInt>,
    ) {
        let n = 8;
        let f = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 850))
                .collect::<Vec<_>>(),
        );
        let g = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 850))
                .collect::<Vec<_>>(),
        );
        let cf = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 2000))
                .collect::<Vec<_>>(),
        );
        let cg = Polynomial::new(
            (0..n)
                .map(|_| rand_signed_bigint(rng, 2000))
                .collect::<Vec<_>>(),
        );
        (f, g, cf, cg)
    }

    /// The runtime-`K` flat-word backend ([`babai_reduce_rns_runtime`]) must
    /// produce exactly the same reduction as the BigInt-karatsuba oracle, across
    /// depths 2–7 (n = 256, 128, 64, 32, 8): depths 2–3 (n>64) take the RNS
    /// multiply, depths 4–7 the flat-word schoolbook multiply.
    #[test]
    fn babai_reduce_rns_runtime_matches_bigint() {
        use super::babai_reduce_rns_runtime;
        type Inputs = (
            Polynomial<BigInt>,
            Polynomial<BigInt>,
            Polynomial<BigInt>,
            Polynomial<BigInt>,
        );
        type Gen = fn(&mut StdRng) -> Inputs;
        // (inputs, k_primes, cap_w) per depth.  k_primes covers ~bits(f)+54;
        // cap_w (limbs) holds the capital plus the mid-reduction shifts.
        let cases: &[(&str, Gen, usize, usize, u64)] = &[
            ("depth2", depth2_inputs as Gen, 5, 4, 0xd2),
            ("depth3", depth3_inputs as Gen, 5, 6, 0xd3),
            ("depth4", depth4_inputs as Gen, 8, 8, 0xd4),
            ("depth5", depth5_inputs as Gen, 12, 12, 0xd5),
            ("depth7", depth7_inputs as Gen, 41, 40, 0xd7),
        ];
        for &(name, gen, k_primes, cap_w, seed) in cases {
            let mut rng = StdRng::seed_from_u64(0x5117_0000 + seed);
            for _ in 0..16 {
                let (f, g, cap_f, cap_g) = gen(&mut rng);

                let (mut bf, mut bg) = (cap_f.clone(), cap_g.clone());
                babai_reduce_bigint(&f, &g, &mut bf, &mut bg).unwrap();

                let (mut rf, mut rg) = (cap_f, cap_g);
                babai_reduce_rns_runtime(&f, &g, &mut rf, &mut rg, k_primes, cap_w).unwrap();

                assert_eq!(bf, rf, "{name}: capital_F mismatch (runtime vs bigint)");
                assert_eq!(bg, rg, "{name}: capital_G mismatch (runtime vs bigint)");
            }
        }
    }

    /// The multiword `babai_reduce_rns_{packed,bigint}` paths only push the
    /// `k·f` *product* through RNS, not the capital.  At a fixed depth the
    /// product does not grow with the capital's magnitude (the `d>53` windowing
    /// keeps `k ≈ 53` bits), so it is `bits(f) + ~54` — measured to peak at a
    /// hard 107 bits at depth 3.  So the prime list need only cover ~107 bits,
    /// not the ≈154-bit capital — hence `NttPrimes24Bit5` (≈114-bit capacity)
    /// rather than `NttPrimes24Bit8`.  (See `product_size_model_matches_…` for
    /// why this grows with depth — at depth 4 the product is ≈155 bits.)
    #[test]
    fn ntt_primes24_5_covers_depth3_product() {
        use crate::rns::NttPrimes24Bit5;
        const MEASURED_PRODUCT_BITS: f64 = 107.0;
        let capacity: f64 = NttPrimes24Bit5::PRIMES
            .iter()
            .map(|&p| f64::log2(p as f64))
            .sum::<f64>()
            - 1.0; // sign bit for centered reconstruction
        assert!(
            capacity >= MEASURED_PRODUCT_BITS,
            "NttPrimes24Bit5 capacity {capacity:.1} < depth-3 product {MEASURED_PRODUCT_BITS:.1}"
        );
    }

    /// Depth 4: the ≈155-bit `k·f` product (`bits(f)≈101` + ≈54) needs the
    /// 8-prime ≈183-bit list; the 5-prime ≈114-bit list would wrap.  Assert K=8
    /// covers the +6σ tail (≈162 bits) and that K=7 would not.
    #[test]
    fn ntt_primes24_8_covers_depth4_product() {
        use crate::rns::NttPrimes24Bit8;
        // avg_f + 6σ_f + 54 + 2, from NTRU_SOLVE_BABAI_COEFF_BITS[4] = (101.62, 1.02, …).
        let required = 101.62 + 6.0 * 1.02 + 54.0 + 2.0; // ≈ 163.8 bits
        let cap = |primes: &[u32]| -> f64 {
            primes.iter().map(|&p| f64::log2(p as f64)).sum::<f64>() - 1.0
        };
        let cap8 = cap(&NttPrimes24Bit8::PRIMES);
        let cap7 = cap(&NttPrimes24Bit8::PRIMES[..7]);
        assert!(
            cap8 >= required,
            "K=8 capacity {cap8:.1} < depth-4 product {required:.1}"
        );
        assert!(
            cap7 < required,
            "K=7 capacity {cap7:.1} unexpectedly covers {required:.1}"
        );
    }

    /// The capacity guard in babai_reduce_rns_generic must fire (debug builds)
    /// when the prime list is too small for the product, rather than silently
    /// wrapping the CRT.  Depth-3 inputs need ~107 product bits; the 2-prime
    /// list has only ~45-bit capacity.
    #[test]
    #[should_panic(expected = "prime list too small")]
    fn rns_reduction_guards_against_undersized_primes() {
        use crate::rns::NttPrimes24Bit2;
        let mut rng = StdRng::seed_from_u64(0xdead_beef);
        let (f, g, mut cf, mut cg) = depth3_inputs(&mut rng);
        let _ = super::babai_reduce_rns_bigint::<2, NttPrimes24Bit2>(&f, &g, &mut cf, &mut cg);
    }

    /// rns-bigint backend must match BigInt karatsuba on realistic depth-4
    /// inputs (n=64, capital ≈303 bits). Wrap-detector for the K=8 prime list.
    #[test]
    fn babai_reduce_rns_bigint_depth4_matches_bigint() {
        use super::babai_reduce_rns_bigint_depth4;

        let mut rng = StdRng::seed_from_u64(0xba_ba_14_d4);
        for _ in 0..16 {
            let (f, g, cap_f, cap_g) = depth4_inputs(&mut rng);

            let (mut bf, mut bg) = (cap_f.clone(), cap_g.clone());
            babai_reduce_bigint(&f, &g, &mut bf, &mut bg).unwrap();

            let (mut rf, mut rg) = (cap_f, cap_g);
            babai_reduce_rns_bigint_depth4(&f, &g, &mut rf, &mut rg).unwrap();

            assert_eq!(bf, rf, "capital_F mismatch (depth 4, rns-bigint)");
            assert_eq!(bg, rg, "capital_G mismatch (depth 4, rns-bigint)");
        }
    }

    /// packed backend (L=5 / 320-bit capital) must match BigInt karatsuba on
    /// realistic depth-4 inputs.  Exercises both the K=8 product sizing and the
    /// widened limb count.
    #[test]
    fn babai_reduce_rns_packed_depth4_matches_bigint() {
        use super::babai_reduce_rns_packed_depth4;

        let mut rng = StdRng::seed_from_u64(0x9ac_ed_d4);
        for _ in 0..16 {
            let (f, g, cap_f, cap_g) = depth4_inputs(&mut rng);

            let (mut bf, mut bg) = (cap_f.clone(), cap_g.clone());
            babai_reduce_bigint(&f, &g, &mut bf, &mut bg).unwrap();

            let (mut pf, mut pg) = (cap_f, cap_g);
            babai_reduce_rns_packed_depth4(&f, &g, &mut pf, &mut pg).unwrap();

            assert_eq!(bf, pf, "capital_F mismatch (depth 4, packed)");
            assert_eq!(bg, pg, "capital_G mismatch (depth 4, packed)");
        }
    }

    // #[test]
    // fn test_gen_poly() {
    //     let mut rng = rng();
    //     let n = 1024;
    //     let mut sum_norms = 0.0;
    //     let num_iterations = 100;
    //     for _ in 0..num_iterations {
    //         let f = gen_poly(n, &mut rng);
    //         sum_norms += f.l2_norm();
    //     }
    //     let average = sum_norms / (num_iterations as f64);
    //     assert!(90.0 < average);
    //     assert!(average < 94.0);
    // }

    #[test]
    fn test_gs_norm() {
        let n = 512;
        let f = (0..n).map(|i| i % 5).collect_vec();
        let g = (0..n).map(|i| (i % 7) - 4).collect_vec();
        let norm_squared = gram_schmidt_norm_squared(&Polynomial::new(f), &Polynomial::new(g));
        let expected = 5992556.183229722f64;
        let difference = (norm_squared - expected).abs();
        assert!(
            difference < 1.0,
            "norm squared was {norm_squared} =/= {expected} (expected)",
        );
    }

    #[test]
    fn test_ntru_solve() {
        let n = 64;
        let f_coefficients = (0..n).map(|i| ((i % 7) as i32) - 4).collect_vec();
        let f = Polynomial::new(f_coefficients).map(|&i| i.into());
        let g_coefficients = (0..n).map(|i| ((i % 5) as i32) - 3).collect_vec();
        let g = Polynomial::new(g_coefficients).map(|&i| i.into());
        let (capital_f, capital_g) = ntru_solve(&f, &g, 1, 2).unwrap();

        let ntru = (f * capital_g - g * capital_f).reduce_by_cyclotomic(n);
        assert_eq!(Polynomial::constant(12289.into()), ntru);
    }

    #[strategy_proptest]
    fn rns_and_i32_babai_reduce_agree(
        #[strategy(1usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(vec(-5i32..5, #_n))] f_coefficients: Vec<i32>,
        #[strategy(vec(-5i32..5, #_n))] g_coefficients: Vec<i32>,
        #[strategy(vec(-115i32..115, #_n))] capital_f_coefficients: Vec<i32>,
        #[strategy(vec(-115i32..115, #_n))] capital_g_coefficients: Vec<i32>,
    ) {
        let n = f_coefficients.len();
        let f = Polynomial::new(f_coefficients);
        let g = Polynomial::new(g_coefficients);

        if g.coefficients.iter().all(|&x| x == 0) {
            return Ok(());
        }

        // Compute NTRU invariant f·G − g·F on the original inputs.
        // Babai reduction preserves this exactly (each step subtracts k·f from F
        // and k·g from G, so f·(G−k·g) − g·(F−k·f) = f·G − g·F).
        let f_bi = f.map(|&x| BigInt::from(x));
        let g_bi = g.map(|&x| BigInt::from(x));
        let cap_f_bi_orig = Polynomial::new(
            capital_f_coefficients
                .iter()
                .map(|&x| BigInt::from(x))
                .collect::<Vec<_>>(),
        );
        let cap_g_bi_orig = Polynomial::new(
            capital_g_coefficients
                .iter()
                .map(|&x| BigInt::from(x))
                .collect::<Vec<_>>(),
        );
        let invariant =
            (f_bi.clone() * cap_g_bi_orig - g_bi.clone() * cap_f_bi_orig).reduce_by_cyclotomic(n);

        let mut capital_f_i32 = Polynomial::new(capital_f_coefficients.clone());
        let mut capital_g_i32 = Polynomial::new(capital_g_coefficients.clone());
        let mut capital_f_rns: Vec<Rns3> = capital_f_coefficients
            .iter()
            .map(|&i| Rns3::from_i32(i))
            .collect();
        let mut capital_g_rns: Vec<Rns3> = capital_g_coefficients
            .iter()
            .map(|&i| Rns3::from_i32(i))
            .collect();

        let _ = babai_reduce_i32(&f, &g, &mut capital_f_i32, &mut capital_g_i32);
        let _ = babai_reduce_rns::<3, NttPrimes3>(&f, &g, &mut capital_f_rns, &mut capital_g_rns);

        // Verify the invariant is preserved by babai_reduce_rns.
        let cap_f_rns_bi = Polynomial::new(
            capital_f_rns
                .iter()
                .map(|r| BigInt::from(r.to_i64()))
                .collect::<Vec<_>>(),
        );
        let cap_g_rns_bi = Polynomial::new(
            capital_g_rns
                .iter()
                .map(|r| BigInt::from(r.to_i64()))
                .collect::<Vec<_>>(),
        );
        let invariant_after = (f_bi * cap_g_rns_bi - g_bi * cap_f_rns_bi).reduce_by_cyclotomic(n);
        prop_assert_eq!(
            invariant,
            invariant_after,
            "rns babai did not preserve NTRU invariant f·G'−g·F'"
        );
    }

    #[strategy_proptest]
    fn bigint_and_smallint_babai_reduce_agree(
        #[strategy(1usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(vec(-5..5, #_n))] f_coefficients: Vec<i32>,
        #[strategy(vec(-5..5, #_n))] g_coefficients: Vec<i32>,
        #[strategy(vec(-115..115, #_n))] capital_f_coefficients: Vec<i32>,
        #[strategy(vec(-115..115, #_n))] capital_g_coefficients: Vec<i32>,
    ) {
        let f_i32 = Polynomial::new(f_coefficients);
        let g_i32 = Polynomial::new(g_coefficients);
        let mut capital_f_i32 = Polynomial::new(capital_f_coefficients);
        let mut capital_g_i32 = Polynomial::new(capital_g_coefficients);
        let f_bigint = f_i32.map(|i| BigInt::from(*i));
        let g_bigint = g_i32.map(|i| BigInt::from(*i));
        let mut capital_f_bigint = capital_f_i32.map(|i| BigInt::from(*i));
        let mut capital_g_bigint = capital_g_i32.map(|i| BigInt::from(*i));

        let _ = babai_reduce_i32(&f_i32, &g_i32, &mut capital_f_i32, &mut capital_g_i32);
        let _ = babai_reduce_bigint(
            &f_bigint,
            &g_bigint,
            &mut capital_f_bigint,
            &mut capital_g_bigint,
        );

        prop_assert_eq!(capital_f_i32.map(|c| BigInt::from(*c)), capital_f_bigint);
        prop_assert_eq!(capital_g_i32.map(|c| BigInt::from(*c)), capital_g_bigint);
    }

    #[test]
    fn test_ntru_gen() {
        let n = 512;
        let seed: [u8; 32] =
            hex::decode("deadbeef00000000deadbeef00000000deadbeef00000000deadbeef00000000")
                .unwrap()
                .try_into()
                .unwrap();
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let (f, g, capital_f, capital_g) = ntru_gen(n, &mut rng);

        println!("f: {}", f);
        println!("g: {}", g);
        println!("capital f: {}", capital_f);
        println!("capital g: {}", capital_g);
        let f_times_capital_g = (f * capital_g).reduce_by_cyclotomic(n);
        let g_times_capital_f = (g * capital_f).reduce_by_cyclotomic(n);
        let difference = f_times_capital_g - g_times_capital_f;
        assert_eq!(Polynomial::constant(12289), difference);
    }

    /// Full n=1024 keygen exercises the runtime-K flat-word babai reduction down
    /// to depth 9 (n=2, ~149 primes).  Validates `fG − gF = q` end-to-end.
    #[test]
    fn test_ntru_gen_1024() {
        let n = 1024;
        let seed: [u8; 32] = *b"\xc0ffee_multiword_alloc_free_2026!";
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let (f, g, capital_f, capital_g) = ntru_gen(n, &mut rng);
        let f_times_capital_g = (f * capital_g).reduce_by_cyclotomic(n);
        let g_times_capital_f = (g * capital_f).reduce_by_cyclotomic(n);
        assert_eq!(
            Polynomial::constant(12289),
            f_times_capital_g - g_times_capital_f
        );
    }
}
