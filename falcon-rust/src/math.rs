use std::vec::IntoIter;

use falcon_profiler::profiling;
use itertools::Itertools;
use num::{BigInt, FromPrimitive, One, Zero};
use num_complex::Complex64;
use rand::RngCore;

use crate::{
    cyclotomic_fourier::CyclotomicFourier,
    falcon_field::{Felt, Q},
    fast_fft::FastFft,
    fixed_point::FixedPoint64,
    inverse::Inverse,
    polynomial::Polynomial,
    rns::{intt_inplace, ntt_inplace, NttPrimeList, NttPrimes24Bit2, NttPrimes24Bit4, Rns},
    samplerz::sampler_z,
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
    let bitsize = |bi: &BigInt| (bi.bits() + 7) & (u64::MAX ^ 7);
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
        if capital_size < size || capital_size >= prev_capital_size {
            break;
        }
        prev_capital_size = capital_size;

        let capital_shift = (capital_size as i64) - 53;
        let capital_f_adjusted = capital_f
            .map(|bi| Complex64::new(i64::try_from(bi >> capital_shift).unwrap() as f64, 0.0))
            .fft();
        let capital_g_adjusted = capital_g
            .map(|bi| Complex64::new(i64::try_from(bi >> capital_shift).unwrap() as f64, 0.0))
            .fft();

        let numerator = capital_f_adjusted.hadamard_mul(&f_star_adjusted)
            + capital_g_adjusted.hadamard_mul(&g_star_adjusted);
        let quotient = numerator.hadamard_div(&denominator_fft).ifft();

        let k = quotient.map(|f| Into::<BigInt>::into(f.re.round() as i64));

        if k.is_zero() {
            break;
        }
        let kf = (k.clone().karatsuba(f)).reduce_by_cyclotomic(n);
        let shifted_kf = kf.map(|bi| bi << (capital_size - size));
        let kg = (k.clone().karatsuba(g)).reduce_by_cyclotomic(n);
        let shifted_kg = kg.map(|bi| bi << (capital_size - size));

        *capital_f -= shifted_kf;
        *capital_g -= shifted_kg;
    }
    Ok(())
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
        (itr.map(|i| i.abs()).max().unwrap() * 2)
            .ilog2()
            .next_multiple_of(8) as usize
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

    // Precompute RNS NTT of f,g for polynomial multiplication.
    let mut f_rns_ntt: Vec<Rns<K, P>> = f.coefficients.iter().map(|&i| Rns::from_i32(i)).collect();
    let mut g_rns_ntt: Vec<Rns<K, P>> = g.coefficients.iter().map(|&i| Rns::from_i32(i)).collect();
    ntt_inplace::<K, P>(&mut f_rns_ntt);
    ntt_inplace::<K, P>(&mut g_rns_ntt);

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
        ntt_inplace::<K, P>(&mut k_rns_ntt);

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
        intt_inplace::<K, P>(&mut kf);
        intt_inplace::<K, P>(&mut kg);

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

    let f_prime = f.field_norm();
    let g_prime = g.field_norm();
    let (capital_f_prime, capital_g_prime) =
        ntru_solve(&f_prime, &g_prime, depth + 1, max_rns_depth)?;

    let capital_f_prime_xsq = capital_f_prime.lift_next_cyclotomic();
    let capital_g_prime_xsq = capital_g_prime.lift_next_cyclotomic();
    let f_minx = f.galois_adjoint();
    let g_minx = g.galois_adjoint();

    let mut capital_f = (capital_f_prime_xsq.karatsuba(&g_minx)).reduce_by_cyclotomic(n);
    let mut capital_g = (capital_g_prime_xsq.karatsuba(&f_minx)).reduce_by_cyclotomic(n);

    let babai_result = if depth <= max_rns_depth {
        match depth {
            1 => {
                babai_rns_with_fallback::<2, NttPrimes24Bit2>(f, g, &mut capital_f, &mut capital_g)
            }
            2 => {
                babai_rns_with_fallback::<4, NttPrimes24Bit4>(f, g, &mut capital_f, &mut capital_g)
            }
            _ => babai_reduce_bigint(f, g, &mut capital_f, &mut capital_g),
        }
    } else {
        babai_reduce_bigint(f, g, &mut capital_f, &mut capital_g)
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
    let n = f.coefficients.len();

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

    let psi_rev = U32Field::bitreversed_powers(n);
    let psi_rev_inv = U32Field::bitreversed_powers_inverse(n);
    let ninv = U32Field::new(n as u32).inverse_or_zero();
    let mut cfp_ntt = capital_f_prime_xsq.map(|c| U32Field::from(*c));
    let mut cgp_ntt = capital_g_prime_xsq.map(|c| U32Field::from(*c));
    let mut gm_ntt = g_minx.map(|c| U32Field::from(*c));
    let mut fm_ntt = f_minx.map(|c| U32Field::from(*c));
    U32Field::fft(&mut cfp_ntt.coefficients, &psi_rev);
    U32Field::fft(&mut cgp_ntt.coefficients, &psi_rev);
    U32Field::fft(&mut gm_ntt.coefficients, &psi_rev);
    U32Field::fft(&mut fm_ntt.coefficients, &psi_rev);
    let mut cf_ntt = cfp_ntt.hadamard_mul(&gm_ntt);
    let mut cg_ntt = cgp_ntt.hadamard_mul(&fm_ntt);
    U32Field::ifft(&mut cf_ntt.coefficients, &psi_rev_inv, ninv);
    U32Field::ifft(&mut cg_ntt.coefficients, &psi_rev_inv, ninv);

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
pub fn ntru_gen(
    n: usize,
    rng: &mut dyn RngCore,
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

        if let Some((capital_f, capital_g)) =
            ntru_solve_entrypoint(f.map(|&i| i as i32), g.map(|&i| i as i32), 2)
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
pub fn ntru_gen_with_rns_depth(
    n: usize,
    rng: &mut dyn RngCore,
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
// fn gen_poly(n: usize, rng: &mut dyn RngCore) -> Polynomial<i16> {
#[profiling]
fn gen_poly(n: usize, rng: &mut dyn RngCore) -> Polynomial<i16> {
    let mu = FixedPoint64::ZERO;
    let sigma_star = FixedPoint64::from(1.43300980528773f64);
    const NUM_COEFFICIENTS: usize = 4096;
    Polynomial {
        coefficients: (0..NUM_COEFFICIENTS)
            .map(|_| {
                sampler_z(
                    mu,
                    sigma_star,
                    sigma_star - FixedPoint64::from(0.001f64),
                    rng,
                )
            })
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
    use proptest::strategy::Just;
    use proptest::prop_assert_eq;
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
}
