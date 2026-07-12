use rand::{Rng, RngExt};

use crate::fixed_point::FixedPoint128;

/// Sample an integer from {0, ..., 18} according to the distribution χ, which
/// is close to the half-Gaussian distribution on the natural numbers with mean
/// 0 and standard deviation equal to sigma_max.
fn base_sampler(bytes: [u8; 9]) -> i16 {
    const RCDT: [u128; 18] = [
        3024686241123004913666,
        1564742784480091954050,
        636254429462080897535,
        199560484645026482916,
        47667343854657281903,
        8595902006365044063,
        1163297957344668388,
        117656387352093658,
        8867391802663976,
        496969357462633,
        20680885154299,
        638331848991,
        14602316184,
        247426747,
        3104126,
        28824,
        198,
        1,
    ];
    // Interpret the 9 bytes as the low 72 bits of a big-endian u128 (top 7 bytes
    // zero). Stack buffer — no per-sample allocation.
    let mut buf = [0u8; 16];
    buf[7..16].copy_from_slice(&bytes);
    let u = u128::from_be_bytes(buf);
    RCDT.into_iter().filter(|r| u < *r).count() as i16
}

/// Compute an integer approximation of 2^63 * ccs * exp(-x).
///
/// Evaluated in FixedPoint128 (64 fractional bits) so the argument carries the full 63 bits
/// of precision the FACCT polynomial's 63-bit constants require. Evaluating this at the old
/// FixedPoint64 (32 fractional bits) quantized the argument to 2^-32 and left approx_exp with a
/// ~2^-33 maximum relative error — below Falcon's mandated sampler precision (GHSA-25rm-9wvm-m38v).
fn approx_exp(x: FixedPoint128, ccs: FixedPoint128) -> u64 {
    // The constants C are used to approximate exp(-x); these
    // constants are taken from FACCT (up to a scaling factor
    // of 2^63):
    //   https://eprint.iacr.org/2018/1234
    //   https://github.com/raykzhao/gaussian
    const C: [u64; 13] = [
        0x00000004741183A3u64,
        0x00000036548CFC06u64,
        0x0000024FDCBF140Au64,
        0x0000171D939DE045u64,
        0x0000D00CF58F6F84u64,
        0x000680681CF796E3u64,
        0x002D82D8305B0FEAu64,
        0x011111110E066FD0u64,
        0x0555555555070F00u64,
        0x155555555581FF00u64,
        0x400000000002B400u64,
        0x7FFFFFFFFFFF4800u64,
        0x8000000000000000u64,
    ];

    let mut y: u64 = C[0];
    // x is in [0, ln(2)]; x.0 = x_real * 2^64; we want floor(x_real * 2^63) = x.0 >> 1.
    // The shift KEEPS the low fractional bits (unlike the old 32-bit `<< 31`, which zero-filled
    // them and quantized the argument to 2^-32). Clamp to 0 if somehow negative.
    let z: u64 = if x.0 < 0 { 0 } else { (x.0 >> 1) as u64 };
    for cu in C.iter().skip(1) {
        let zy = (z as u128) * (y as u128);
        y = cu - ((zy >> 63) as u64);
    }

    // ccs is in [0, 1]; ccs.0 = ccs_real * 2^64; we want floor(ccs_real * 2^63) = ccs.0 >> 1.
    let z2: u64 = if ccs.0 < 0 { 0 } else { (ccs.0 >> 1) as u64 };

    (((z2 as u128) * (y as u128)) >> 63) as u64
}

/// A random bool that is true with probability ≈ ccs · exp(−x).
fn ber_exp(x: FixedPoint128, ccs: FixedPoint128, random_bytes: [u8; 7]) -> bool {
    let s = (x / FixedPoint128::LN_2).trunc() as usize;
    let r = x - FixedPoint128::LN_2 * FixedPoint128::from(s as i32);
    let shamt = usize::min(s, 63);
    let z = ((((approx_exp(r, ccs) as u128) << 1) - 1) >> shamt) as u64;
    let mut w = 0i16;
    for (index, i) in (0..64).step_by(8).rev().enumerate() {
        let byte = random_bytes[index];
        w = (byte as i16) - (((z >> i) & 0xff) as i16);
        if w != 0 {
            break;
        }
    }
    w < 0
}

/// Sample an integer from the Gaussian distribution with given mean (mu) and
/// standard deviation (sigma).
/// Precomputed, sigma-dependent constants for [`SamplerZCtx::sample`].  Building
/// these once and sampling many times (as `gen_poly` does — 4096 samples at one
/// fixed sigma) avoids redoing the setup — notably the `FixedPoint64` division
/// `1/sigma` — on every sample.
pub(crate) struct SamplerZCtx {
    inv_2sigma_max_sq: FixedPoint128,
    dss: FixedPoint128,
    ccs: FixedPoint128,
}

impl SamplerZCtx {
    pub(crate) fn new(sigma: FixedPoint128, sigma_min: FixedPoint128) -> Self {
        let sigma_max = FixedPoint128::from(1.8205f64);
        let inv_2sigma_max_sq =
            FixedPoint128::ONE / (FixedPoint128::from(2.0f64) * sigma_max * sigma_max);
        let isigma = FixedPoint128::ONE / sigma;
        let dss = FixedPoint128::from(0.5f64) * isigma * isigma;
        let ccs = sigma_min * isigma;
        Self {
            inv_2sigma_max_sq,
            dss,
            ccs,
        }
    }

    /// Sample one integer centred at `mu` (the only per-sample-varying input).
    pub(crate) fn sample<R: Rng + ?Sized>(&self, mu: FixedPoint128, rng: &mut R) -> i16 {
        let s = mu.floor().trunc();
        let r = mu - FixedPoint128::from(s);
        loop {
            let z0 = base_sampler(rng.random());
            let random_byte: u8 = rng.random();
            let b = (random_byte & 1) as i16;
            let z = b + ((b << 1) - 1) * z0;
            let zf_min_r = FixedPoint128::from(z as i32) - r;
            let x = zf_min_r * zf_min_r * self.dss
                - FixedPoint128::from(z0 as i32 * z0 as i32) * self.inv_2sigma_max_sq;
            if ber_exp(x, self.ccs, rng.random()) {
                return z + (s as i16);
            }
        }
    }
}

pub(crate) fn sampler_z<FP: Into<FixedPoint128>>(
    mu: FP,
    sigma: FP,
    sigma_min: FP,
    rng: &mut dyn Rng,
) -> i16 {
    SamplerZCtx::new(sigma.into(), sigma_min.into()).sample(mu.into(), rng)
}

#[cfg(test)]
mod test {
    use core::convert::Infallible;
    use itertools::Itertools;
    use rand::rand_core::TryRng;
    use rand::{rng, RngExt};
    use std::{thread::sleep, time::Duration};

    use crate::fixed_point::{FixedPoint128, FixedPoint64};
    use crate::samplerz::{approx_exp, ber_exp, sampler_z};
    use rand::{rngs::StdRng, SeedableRng};

    /// Verbatim copy of `base_sampler` *before* the allocation optimization
    /// (commit 892304c): builds the 16-byte buffer via `Vec::concat`.
    fn old_base_sampler(bytes: [u8; 9]) -> i16 {
        const RCDT: [u128; 18] = [
            3024686241123004913666,
            1564742784480091954050,
            636254429462080897535,
            199560484645026482916,
            47667343854657281903,
            8595902006365044063,
            1163297957344668388,
            117656387352093658,
            8867391802663976,
            496969357462633,
            20680885154299,
            638331848991,
            14602316184,
            247426747,
            3104126,
            28824,
            198,
            1,
        ];
        let u = u128::from_be_bytes([vec![0u8; 7], bytes.to_vec()].concat().try_into().unwrap());
        RCDT.into_iter().filter(|r| u < *r).count() as i16
    }

    /// Verbatim copy of `sampler_z` *before* the constant-hoisting
    /// optimization: recomputes the sigma constants and uses `old_base_sampler`.
    fn old_sampler_z(
        mu: FixedPoint128,
        sigma: FixedPoint128,
        sigma_min: FixedPoint128,
        rng: &mut dyn rand::Rng,
    ) -> i16 {
        let sigma_max = FixedPoint128::from(1.8205f64);
        let inv_2sigma_max_sq =
            FixedPoint128::ONE / (FixedPoint128::from(2.0f64) * sigma_max * sigma_max);
        let isigma = FixedPoint128::ONE / sigma;
        let dss = FixedPoint128::from(0.5f64) * isigma * isigma;
        let s = mu.floor().trunc();
        let r = mu - FixedPoint128::from(s);
        let ccs = sigma_min * isigma;
        loop {
            let z0 = old_base_sampler(rng.random());
            let random_byte: u8 = rng.random();
            let b = (random_byte & 1) as i16;
            let z = b + ((b << 1) - 1) * z0;
            let zf_min_r = FixedPoint128::from(z as i32) - r;
            let x = zf_min_r * zf_min_r * dss
                - FixedPoint128::from(z0 as i32 * z0 as i32) * inv_2sigma_max_sq;
            if ber_exp(x, ccs, rng.random()) {
                return z + (s as i16);
            }
        }
    }

    /// The optimized `sampler_z`/`base_sampler` must produce byte-for-byte the
    /// same output sequence as the pre-optimization code on the same RNG stream —
    /// across a range of (sigma, sigma_min, mu). Identical outputs over a long run
    /// also prove the randomness is consumed in the same order/amount (any
    /// mismatch would desync the two RNGs and surface as a divergence).
    #[test]
    fn sampler_z_matches_pre_optimization() {
        let cases = [
            (1.43300980528773f64, 1.43200980528773f64), // keygen gen_poly
            (1.2778336969128337, 1.2778336969128337),   // ~signer floor
            (1.7794, 1.32),
            (1.8205, 1.17),
        ];
        for &(sig, smin) in &cases {
            let sigma = FixedPoint128::from(sig);
            let sigma_min = FixedPoint128::from(smin);
            for &mu_v in &[0.0f64, 0.3, -0.7, 5.5, 12.25] {
                let mu = FixedPoint128::from(mu_v);
                let mut r_old = StdRng::seed_from_u64(0x5A3D_0000 ^ sig.to_bits());
                let mut r_new = StdRng::seed_from_u64(0x5A3D_0000 ^ sig.to_bits());
                for i in 0..20_000 {
                    let a = old_sampler_z(mu, sigma, sigma_min, &mut r_old);
                    let b = sampler_z(mu, sigma, sigma_min, &mut r_new);
                    assert_eq!(a, b, "divergence at sample {i} (sigma={sig}, mu={mu_v})");
                }
            }
        }
    }

    /// RNG used only for testing purposes, whereby the produced
    /// string of random bytes is equal to the one it is initialized
    /// with. Whatever you do, do not use this RNG in production.
    struct UnsafeBufferRng {
        buffer: Vec<u8>,
        index: usize,
    }

    impl UnsafeBufferRng {
        fn new(buffer: &[u8]) -> Self {
            Self {
                buffer: buffer.to_vec(),
                index: 0,
            }
        }

        fn next(&mut self) -> u8 {
            if self.buffer.len() <= self.index {
                // panic!("Ran out of buffer.");
                sleep(Duration::from_millis(10));
                0u8
            } else {
                let return_value = self.buffer[self.index];
                self.index += 1;
                return_value
            }
        }
    }

    impl TryRng for UnsafeBufferRng {
        type Error = Infallible;

        fn try_next_u32(&mut self) -> Result<u32, Self::Error> {
            Ok(u32::from_le_bytes([self.next(), 0, 0, 0]))
        }

        fn try_next_u64(&mut self) -> Result<u64, Self::Error> {
            Ok(u64::from_le_bytes([self.next(), 0, 0, 0, 0, 0, 0, 0]))
        }

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), Self::Error> {
            for d in dest.iter_mut() {
                *d = self.next();
            }
            Ok(())
        }
    }

    #[test]
    fn test_unsafe_buffer_rng() {
        let seed_bytes = hex::decode("7FFECD162AE2").unwrap();
        let mut rng = UnsafeBufferRng::new(&seed_bytes);
        let generated_bytes = (0..seed_bytes.len()).map(|_| rng.next()).collect_vec();
        assert_eq!(seed_bytes, generated_bytes);
    }

    #[test]
    fn test_approx_exp() {
        // Known answers were generated with the following sage script (high precision):
        //   https://eprint.iacr.org/2016/1055 table 3.2
        // With the FixedPoint128 (64-bit fractional) acceptance path the only error is the f64
        // rounding of the inputs (~2^11 on this 2^63 scale). Falcon mandates a sampler relative
        // precision of about 2^-40..-46; 2^-40 on the 2^63 output scale is 2^23 absolute, so we
        // assert the tight 2^23 bound. (The previous FixedPoint64 sampler was ~2^30 here — 2^-33
        // relative — and this loosened tolerance had been raised to 2^40 to hide it; GHSA-25rm-9wvm-m38v.)
        let precision = 1u64 << 23;
        let kats: [(f64, f64, u64); 10] = [
            (0.2314993926072656, 0.8148006314615972, 5962140072160879737),
            (0.2648875572812225, 0.12769669655309035, 903712282351034505),
            (0.11251957513682391, 0.9264611470305881, 7635725498677341553),
            (0.04353439307256617, 0.5306497137523327, 4685877322232397936),
            (0.41834495299784347, 0.879438856118578, 5338392138535350986),
            (
                0.32579398973228557,
                0.16513412873289002,
                1099603299296456803,
            ),
            (0.5939508073919817, 0.029776019144967303, 151637565622779016),
            (0.2932367999399056, 0.37123847662857923, 2553827649386670452),
            (0.5005699297417507, 0.31447208863888976, 1758235618083658825),
            (0.4876437338498085, 0.6159515298936868, 3488632981903743976),
        ];
        for (x, ccs, answer) in kats {
            let result = approx_exp(FixedPoint128::from(x), FixedPoint128::from(ccs));
            let difference = (answer as i128) - (result as i128);
            assert!(
                (difference * difference) as u128 <= (precision as u128) * (precision as u128),
                "answer: {answer} versus approximation: {result}\ndifference: {difference} whereas precision: {precision}"
            );
        }
    }

    /// Regression guard for GHSA-25rm-9wvm-m38v: approx_exp must meet Falcon's mandated sampler
    /// precision (about 2^-40..-46 relative). Measures the maximum relative error of approx_exp
    /// against true exp() over the sampler's argument domain [0, ln 2]. The FixedPoint64 sampler
    /// this replaced was 2^-33 (about 128x too coarse); FixedPoint128 restores about 2^-52 or better.
    #[test]
    fn approx_exp_meets_falcon_precision() {
        let n: u64 = 500_000;
        let two63 = 9_223_372_036_854_775_808.0f64; // 2^63
        let mut maxrel = 0.0f64;
        for i in 0..=n {
            let x = (i as f64 / n as f64) * std::f64::consts::LN_2;
            let approx = approx_exp(FixedPoint128::from(x), FixedPoint128::ONE) as f64 / two63;
            let truev = (-x).exp();
            let rel = ((approx - truev) / truev).abs();
            if rel > maxrel {
                maxrel = rel;
            }
        }
        let requirement = 2f64.powi(-40); // loosest Falcon precision bound in the literature
        assert!(
            maxrel < requirement,
            "approx_exp max relative error {maxrel:.3e} (2^{:.2}) exceeds Falcon's 2^-40 requirement",
            maxrel.log2()
        );
    }

    #[test]
    fn test_ber_exp() {
        let kats = [
            (
                1.268_314_048_020_498_4,
                0.749_990_853_267_664_9,
                hex::decode("ea000000000000").unwrap(),
                false,
            ),
            (
                0.001_563_917_959_143_409_6,
                0.749_990_853_267_664_9,
                hex::decode("6c000000000000").unwrap(),
                true,
            ),
            (
                0.017_921_215_753_999_235,
                0.749_990_853_267_664_9,
                hex::decode("c2000000000000").unwrap(),
                false,
            ),
            (
                0.776_117_648_844_980_6,
                0.751_181_554_542_520_8,
                hex::decode("58000000000000").unwrap(),
                true,
            ),
        ];
        for (x, ccs, bytes, answer) in kats {
            assert_eq!(
                answer,
                ber_exp(
                    FixedPoint128::from(x),
                    FixedPoint128::from(ccs),
                    bytes.try_into().unwrap()
                )
            );
        }
    }

    #[test]
    fn test_sampler_z() {
        let sigma_min = 1.277833697f64;
        let kats = [
            (-91.90471153063714,1.7037990414754918,hex::decode("0fc5442ff043d66e91d1ea000000000000cac64ea5450a22941edc6c").unwrap(),-92i16),
            (-8.322564895434937,1.7037990414754918,hex::decode("f4da0f8d8444d1a77265c2000000000000ef6f98bbbb4bee7db8d9b3").unwrap(),-8),
            (-19.096516109216804,1.7035823083824078,hex::decode("db47f6d7fb9b19f25c36d6000000000000b9334d477a8bc0be68145d").unwrap(),-20),
            (-11.335543982423326, 1.7035823083824078, hex::decode("ae41b4f5209665c74d00dc000000000000c1a8168a7bb516b3190cb42c1ded26cd52000000000000aed770eca7dd334e0547bcc3c163ce0b").unwrap(), -12),
            (7.9386734193997555, 1.6984647769450156, hex::decode("31054166c1012780c603ae0000000000009b833cec73f2f41ca5807c000000000000c89c92158834632f9b1555").unwrap(), 8),
            (-28.990850086867255, 1.6984647769450156, hex::decode("737e9d68a50a06dbbc6477").unwrap(), -30),
            (-9.071257914091655, 1.6980782114808988, hex::decode("a98ddd14bf0bf22061d632").unwrap(), -10),
            (-43.88754568839566, 1.6980782114808988, hex::decode("3cbf6818a68f7ab9991514").unwrap(), -41),
            (-58.17435547946095,1.7010983419195522,hex::decode("6f8633f5bfa5d26848668e0000000000003d5ddd46958e97630410587c").unwrap(),-61),
            (-43.58664906684732, 1.7010983419195522, hex::decode("272bc6c25f5c5ee53f83c40000000000003a361fbc7cc91dc783e20a").unwrap(), -46),
            (-34.70565203313315, 1.7009387219711465, hex::decode("45443c59574c2c3b07e2e1000000000000d9071e6d133dbe32754b0a").unwrap(), -34),
            (-44.36009577368896, 1.7009387219711465, hex::decode("6ac116ed60c258e2cbaeab000000000000728c4823e6da36e18d08da0000000000005d0cc104e21cc7fd1f5ca8000000000000d9dbb675266c928448059e").unwrap(), -44),
            (-21.783037079346236, 1.6958406126012802, hex::decode("68163bc1e2cbf3e18e7426").unwrap(), -23),
            (-39.68827784633828, 1.6958406126012802, hex::decode("d6a1b51d76222a705a0259").unwrap(), -40),
            (-18.488607061056847, 1.6955259305261838, hex::decode("f0523bfaa8a394bf4ea5c10000000000000f842366fde286d6a30803").unwrap(), -22),
            (-48.39610939101591, 1.6955259305261838, hex::decode("87bd87e63374cee62127fc0000000000006931104aab64f136a0485b").unwrap(), -50),
            // Regression vector for GHSA-25rm-9wvm-m38v. Adversarially constructed to land on a
            // rejection boundary the pre-fix ~2^-33 acceptance-probability error flips: the old
            // FixedPoint64 acceptance path returns 102, the current FixedPoint128 path returns 100.
            // Only reproduces on the FixedPoint64 input grid the real signer uses (raw f64 -> 102),
            // which is why this table quantizes its inputs to FixedPoint64 before sampling.
            (100.000_000_5, 1.703_799_041_475_491_8, hex::decode("7895f43f559b370df80170198e14b8b5ffffffffffffffffffff0000000000000000").unwrap(), 100),
        ];
        for (i, (mu, sigma, random_bytes, answer)) in kats.into_iter().enumerate() {
            assert_eq!(
                sampler_z(
                    FixedPoint64::from(mu),
                    FixedPoint64::from(sigma),
                    FixedPoint64::from(sigma_min),
                    &mut UnsafeBufferRng::new(&random_bytes)
                ),
                answer,
                "error in kat {i}"
            );
        }
    }

    /// Regression vector for GHSA-25rm-9wvm-m38v. Commit bf2cf00 evaluated the Gaussian
    /// acceptance probability at ~2^-33 relative precision (FixedPoint64), below the bound
    /// Falcon's Rényi-divergence security analysis requires. This is an adversarially
    /// constructed input that lands on a rejection boundary the 2^-33 error flips: on the
    /// pre-fix (FixedPoint64) code it produces false; the current FixedPoint128 code produces
    /// true. It guards against a re-regression that a random test vector would miss (the
    /// discrepancy has natural frequency ~2^-33). The sampler_z-level counterpart lives in
    /// `test_sampler_z`'s KAT table.
    #[test]
    fn ber_exp_regression_vector_ghsa_25rm_9wvm_m38v() {
        // Pre-fix (FixedPoint64) returned false; the correct-precision path returns true.
        let x = FixedPoint128::from(0.200_000_167_618_381_2_f64);
        let ccs = FixedPoint128::from(0.749_990_853_318_824_9_f64);
        let bytes: [u8; 7] = hex::decode("9d31c1a6ca10d7").unwrap().try_into().unwrap();
        assert!(
            ber_exp(x, ccs, bytes),
            "ber_exp regression: acceptance decision flipped"
        );
    }

    #[test]
    fn endianness() {
        let bytes: [u8; 9] = rng().random();
        let u0 = u128::from_le_bytes(
            [bytes.into_iter().rev().collect_vec(), vec![0u8; 7]]
                .concat()
                .try_into()
                .unwrap(),
        );
        let u1 = u128::from_be_bytes([vec![0u8; 7], bytes.to_vec()].concat().try_into().unwrap());
        assert_eq!(u0, u1);

        assert!(u0 % (1u128 << 64) != 0); // vanishingly small false positive prob
    }
}
