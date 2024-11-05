use std::vec::IntoIter;

use itertools::Itertools;
use num::{BigInt, FromPrimitive, One, Zero};
use num_complex::Complex64;
use rand::RngCore;

use crate::{
    cyclotomic_fourier::CyclotomicFourier,
    falcon_field::{Felt, Q},
    fast_fft::FastFft,
    inverse::Inverse,
    polynomial::Polynomial,
    samplerz::sampler_z,
    u32_field::U32Field,
};

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

    let mut counter = 0;
    loop {
        let capital_size = [
            capital_f.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
            capital_g.map(bitsize).fold(0, |a, &b| u64::max(a, b)),
            53,
        ]
        .into_iter()
        .max()
        .unwrap();

        if capital_size < size {
            break;
        }
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

        counter += 1;
        if counter > 1000 {
            // If we get here, that means that (with high likelihood) we are in an
            // infinite loop. We know it happens from time to time -- seldomly, but it
            // does. It would be nice to fix that! But in order to fix it we need to be
            // able to reproduce it, and for that we need test vectors. So print them
            // and hope that one day they circle back to the implementor.
            return Err(format!("Encountered infinite loop in babai_reduce of falcon-rust.\n\
            Please help the developer(s) fix it! You can do this by sending them the inputs to the function that caused the behavior:\n\
            f: {:?}\n\
            g: {:?}\n\
            capital_f: {:?}\n\
            capital_g: {:?}\n", f.coefficients, g.coefficients, capital_f.coefficients, capital_g.coefficients));
        }
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
pub fn babai_reduce_i32(
    f: &Polynomial<i32>,
    g: &Polynomial<i32>,
    capital_f: &mut Polynomial<i32>,
    capital_g: &mut Polynomial<i32>,
) -> Result<(), String> {
    let f_ntt: Polynomial<U32Field> = f.map(|&i| U32Field::new(i)).fft();
    let g_ntt: Polynomial<U32Field> = g.map(|&i| U32Field::new(i)).fft();

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
        .map(|i| Complex64::new(i64::from(i >> shift) as f64, 0.0))
        .fft();
    let g_adjusted = g
        .map(|i| Complex64::new(i64::from(i >> shift) as f64, 0.0))
        .fft();

    let f_star_adjusted = f_adjusted.map(|c| c.conj());
    let g_star_adjusted = g_adjusted.map(|c| c.conj());
    let denominator_fft =
        f_adjusted.hadamard_mul(&f_star_adjusted) + g_adjusted.hadamard_mul(&g_star_adjusted);

    let mut counter = 0;
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

        if capital_size < size {
            break;
        }
        let capital_shift = (capital_size as i64) - 53;
        let capital_f_adjusted = capital_f
            .map(|bi| Complex64::new(i64::from(bi >> capital_shift) as f64, 0.0))
            .fft();
        let capital_g_adjusted = capital_g
            .map(|bi| Complex64::new(i64::from(bi >> capital_shift) as f64, 0.0))
            .fft();

        let numerator = capital_f_adjusted.hadamard_mul(&f_star_adjusted)
            + capital_g_adjusted.hadamard_mul(&g_star_adjusted);
        let quotient = numerator.hadamard_div(&denominator_fft).ifft();

        let k_ntt = quotient.map(|f| U32Field::new(f.re.round() as i32)).fft();

        if k_ntt.is_zero() {
            break;
        }

        let kf_ntt = k_ntt.hadamard_mul(&f_ntt).ifft();
        let kg_ntt = k_ntt.hadamard_mul(&g_ntt).ifft();

        let kf = kf_ntt.map(|p| p.balanced_value());
        let kg = kg_ntt.map(|p| p.balanced_value());

        *capital_f -= kf;
        *capital_g -= kg;

        counter += 1;
        if counter > 1000 {
            // If we get here, that means that (with high likelihood) we are in an
            // infinite loop. We know it happens from time to time -- seldomly, but it
            // does. It would be nice to fix that! But in order to fix it we need to be
            // able to reproduce it, and for that we need test vectors. So print them
            // and hope that one day they circle back to the implementor.
            return Err(format!("Encountered infinite loop in babai_reduce of falcon-rust.\n\\
            Please help the developer(s) fix it! You can do this by sending them the inputs to the function that caused the behavior:\n\\
            f: {:?}\n\\
            g: {:?}\n\\
            capital_f: {:?}\n\\
            capital_g: {:?}\n", f.coefficients, g.coefficients, capital_f.coefficients, capital_g.coefficients));
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

/// Solve the NTRU equation. Given f, g in ZZ[X], find F, G in ZZ[X].
/// such that
///
///    f G - g F = q  mod (X^n + 1)
///
/// Algorithm 6 of the specification [1, p.35].
///
/// [1]: https://falcon-sign.info/falcon.pdf
fn ntru_solve(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
) -> Option<(Polynomial<BigInt>, Polynomial<BigInt>)> {
    let n = f.coefficients.len();
    if n == 1 {
        let (gcd, u, v) = xgcd(&f.coefficients[0], &g.coefficients[0]);
        if gcd != BigInt::one() {
            return None;
        }
        return Some((
            (Polynomial::new(vec![-v * BigInt::from_u32(Q).unwrap()])),
            Polynomial::new(vec![u * BigInt::from_u32(Q).unwrap()]),
        ));
    }

    let f_prime = f.field_norm();
    let g_prime = g.field_norm();
    let (capital_f_prime, capital_g_prime) = ntru_solve(&f_prime, &g_prime)?;

    let capital_f_prime_xsq = capital_f_prime.lift_next_cyclotomic();
    let capital_g_prime_xsq = capital_g_prime.lift_next_cyclotomic();
    let f_minx = f.galois_adjoint();
    let g_minx = g.galois_adjoint();

    let mut capital_f = (capital_f_prime_xsq.karatsuba(&g_minx)).reduce_by_cyclotomic(n);
    let mut capital_g = (capital_g_prime_xsq.karatsuba(&f_minx)).reduce_by_cyclotomic(n);

    match babai_reduce_bigint(f, g, &mut capital_f, &mut capital_g) {
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

fn ntru_solve_entrypoint(
    f: Polynomial<i32>,
    g: Polynomial<i32>,
) -> Option<(Polynomial<i32>, Polynomial<i32>)> {
    let n = f.coefficients.len();

    let g_prime = g.field_norm().map(|c| BigInt::from(*c));
    let f_prime = f.field_norm().map(|c| BigInt::from(*c));
    let (capital_f_prime_bi, capital_g_prime_bi) = ntru_solve(&f_prime, &g_prime)?;

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
    let ninv = U32Field::new(n as i32).inverse_or_zero();
    let mut cfp_ntt = capital_f_prime_xsq.map(|c| U32Field::new(*c));
    let mut cgp_ntt = capital_g_prime_xsq.map(|c| U32Field::new(*c));
    let mut gm_ntt = g_minx.map(|c| U32Field::new(*c));
    let mut fm_ntt = f_minx.map(|c| U32Field::new(*c));
    U32Field::fft(&mut cfp_ntt.coefficients, &psi_rev);
    U32Field::fft(&mut cgp_ntt.coefficients, &psi_rev);
    U32Field::fft(&mut gm_ntt.coefficients, &psi_rev);
    U32Field::fft(&mut fm_ntt.coefficients, &psi_rev);
    let mut cf_ntt = cfp_ntt.hadamard_mul(&gm_ntt);
    let mut cg_ntt = cgp_ntt.hadamard_mul(&fm_ntt);
    U32Field::ifft(&mut cf_ntt.coefficients, &psi_rev_inv, ninv);
    U32Field::ifft(&mut cg_ntt.coefficients, &psi_rev_inv, ninv);

    let mut capital_f = cf_ntt.map(|c| c.balanced_value());
    let mut capital_g = cg_ntt.map(|c| c.balanced_value());

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

        let f_ntt = f.map(|&i| Felt::new(i)).fft();
        if f_ntt.coefficients.iter().any(|e| e.is_zero()) {
            continue;
        }
        let gamma = gram_schmidt_norm_squared(&f, &g);
        if gamma > 1.3689f64 * (Q as f64) {
            continue;
        }

        if let Some((capital_f, capital_g)) =
            ntru_solve_entrypoint(f.map(|&i| i as i32), g.map(|&i| i as i32))
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
fn gen_poly(n: usize, rng: &mut dyn RngCore) -> Polynomial<i16> {
    let mu = 0.0;
    let sigma_star = 1.43300980528773;
    const NUM_COEFFICIENTS: usize = 4096;
    Polynomial {
        coefficients: (0..NUM_COEFFICIENTS)
            .map(|_| sampler_z(mu, sigma_star, sigma_star - 0.001, rng))
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
fn gram_schmidt_norm_squared(f: &Polynomial<i16>, g: &Polynomial<i16>) -> f64 {
    let n = f.coefficients.len();
    let norm_f_squared = f.l2_norm_squared();
    let norm_g_squared = g.l2_norm_squared();
    let gamma1 = norm_f_squared + norm_g_squared;

    let f_fft = f.map(|i| Complex64::new(*i as f64, 0.0)).fft();
    let g_fft = g.map(|i| Complex64::new(*i as f64, 0.0)).fft();
    let f_adj_fft = f_fft.map(|c| c.conj());
    let g_adj_fft = g_fft.map(|c| c.conj());
    let ffgg_fft = f_fft.hadamard_mul(&f_adj_fft) + g_fft.hadamard_mul(&g_adj_fft);
    let ffgg_fft_inverse = ffgg_fft.hadamard_inv();
    let qf_over_ffgg_fft = f_adj_fft
        .map(|c| c * (Q as f64))
        .hadamard_mul(&ffgg_fft_inverse);
    let qg_over_ffgg_fft = g_adj_fft
        .map(|c| c * (Q as f64))
        .hadamard_mul(&ffgg_fft_inverse);
    let norm_f_over_ffgg_squared = qf_over_ffgg_fft
        .coefficients
        .iter()
        .map(|c| (c * c.conj()).re)
        .sum::<f64>()
        / (n as f64);
    let norm_g_over_ffgg_squared = qg_over_ffgg_fft
        .coefficients
        .iter()
        .map(|c| (c * c.conj()).re)
        .sum::<f64>()
        / (n as f64);

    let gamma2 = norm_f_over_ffgg_squared + norm_g_over_ffgg_squared;

    f64::max(gamma1, gamma2)
}

#[cfg(test)]
mod test {

    use std::str::FromStr;

    use itertools::Itertools;
    use num::{BigInt, FromPrimitive};
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest::strategy::Just;
    use rand::{rngs::StdRng, SeedableRng};
    use test_strategy::proptest as strategy_proptest;

    use crate::{
        math::{babai_reduce_i32, gram_schmidt_norm_squared, ntru_gen, ntru_solve},
        polynomial::Polynomial,
    };

    use super::babai_reduce_bigint;
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
    fn babai_infinite_loop_bigint() {
        let (f, g, mut capital_f, mut capital_g) = babai_infinite_loop_polynomials();
        assert!(babai_reduce_bigint(&f, &g, &mut capital_f, &mut capital_g).is_ok())
    }

    // #[test]
    // fn test_gen_poly() {
    //     let mut rng = thread_rng();
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
        let difference = norm_squared - 5992556.183229722;
        assert!(
            difference * difference < 0.00001,
            "norm squared was {} =/= 5992556.183229722",
            norm_squared,
        );
    }

    #[test]
    fn test_ntru_solve() {
        let n = 64;
        let f_coefficients = (0..n).map(|i| ((i % 7) as i32) - 4).collect_vec();
        let f = Polynomial::new(f_coefficients).map(|&i| i.into());
        let g_coefficients = (0..n).map(|i| ((i % 5) as i32) - 3).collect_vec();
        let g = Polynomial::new(g_coefficients).map(|&i| i.into());
        let (capital_f, capital_g) = ntru_solve(&f, &g).unwrap();

        let expected_capital_f: [i16; 64] = [
            -221, -19, 133, 81, -488, -112, 189, -75, -112, -223, 143, 241, -249, 33, 47, -16, 32,
            -145, 183, -57, -99, 104, -44, 78, -129, 26, 77, -88, 52, -36, 69, -66, -37, 80, -45,
            32, -67, 93, -24, -79, 87, -49, 68, -116, 60, 108, -158, 68, -52, 87, -32, -116, 233,
            -120, -111, 65, 119, 144, -307, -98, 295, -163, -194, -325,
        ];
        let expected_capital_g: [i16; 64] = [
            -861, 625, -531, 151, 80, 11, 132, 547, -308, 4, 184, -134, -74, -61, 215, -2, -188,
            40, 104, -38, -59, 21, 51, -12, -101, 86, 12, -40, 0, -31, 86, -72, 7, 24, -32, 46,
            -71, 53, 0, -21, 23, -49, 60, -16, -38, 30, 18, 3, -41, -42, 114, 2, -119, 80, -64, 95,
            -37, -18, 238, -429, 87, 193, -3, -111,
        ];
        assert_eq!(
            expected_capital_f
                .map(|i| BigInt::from_i16(i).unwrap())
                .to_vec(),
            capital_f.coefficients
        );
        assert_eq!(
            expected_capital_g
                .map(|i| BigInt::from_i16(i).unwrap())
                .to_vec(),
            capital_g.coefficients
        );

        let ntru = (f * capital_g - g * capital_f).reduce_by_cyclotomic(n);
        assert_eq!(Polynomial::constant(12289.into()), ntru);
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

        let small_int_result =
            babai_reduce_i32(&f_i32, &g_i32, &mut capital_f_i32, &mut capital_g_i32);
        let big_int_result = babai_reduce_bigint(
            &f_bigint,
            &g_bigint,
            &mut capital_f_bigint,
            &mut capital_g_bigint,
        );

        prop_assert_eq!(small_int_result.is_err(), big_int_result.is_err());
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
