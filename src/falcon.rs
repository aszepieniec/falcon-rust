use std::fmt::Display;

use itertools::Itertools;
use rand::{rngs::StdRng, thread_rng, Rng, RngCore, SeedableRng};
use rand_distr::num_traits::Zero;

use crate::{
    field::{Felt, Q},
    ntt::{intt, ntt},
    polynomial::Polynomial,
    samplerz::sampler_z,
};

pub const fn logn(mut n: u32) -> usize {
    let mut ctr = 0;
    while n != 0 {
        ctr += 1;
        n = n >> 1;
    }
    ctr
}

/// Bytelength of the signing salt and header
const HEAD_LEN: usize = 1;
const SALT_LEN: usize = 40;
const SEED_LEN: usize = 56;

pub struct Params {
    pub n: usize,
    pub sigma: f64,
    pub sigmin: f64,
    pub sig_bound: usize,
    pub sig_bytelen: usize,
}

pub enum FalconVariant {
    Falcon512,
    Falcon1024,
}

impl Params {
    pub const fn new(variant: FalconVariant) -> Self {
        match variant {
            FalconVariant::Falcon512 => Self {
                n: 512,
                sigma: 165.7366171829776,
                sigmin: 1.2778336969128337,
                sig_bound: 34034726,
                sig_bytelen: 666,
            },
            FalconVariant::Falcon1024 => Self {
                n: 1024,
                sigma: 168.38857144654395,
                sigmin: 1.298280334344292,
                sig_bound: 70265242,
                sig_bytelen: 1280,
            },
        }
    }
}

fn gen_poly(n: usize, sigma_min: f64, rng: &mut dyn RngCore) -> Vec<Felt> {
    let mu = 0.0;
    let sigma = 1.43300980528773;
    (0..4096)
        .map(|_| sampler_z(mu, sigma, sigma_min, rng))
        .collect_vec()
        .chunks(4096 / n)
        .map(|ch| Felt(ch.iter().sum()))
        .collect_vec()
}

/// Sample 4 small polynomials f, g, F, G such that
/// f * G - g * F = q mod (X^n + 1). This function implements
/// algorithm 5 (NTRUgen) on page 33 of the documentation [1].
/// [1]: https://falcon-sign.info/falcon.pdf
fn ntru_gen(
    n: usize,
    sigma_min: f64,
    seed: [u8; 32],
) -> (Polynomial, Polynomial, Polynomial, Polynomial) {
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    loop {
        let f = gen_poly(n, sigma_min, &mut rng);
        let g = gen_poly(n, sigma_min, &mut rng);
        // if gs_norm(f, g, Q) > 1.3689 * Q {
        //     continue;
        // }
        let f_ntt = ntt(&f);
        if f_ntt.iter().any(|e| e.is_zero()) {
            continue;
        }
        // let (capital_f, capital_g) = ntru_solve(f, g);
        // return (
        //     Polynomial::new(&f),
        //     Polynomial::new(&g),
        // Polynomial::new(&capital_f),
        // Polynomial::new(&capital_g),
        // );
        todo!()
    }
}

#[derive(Debug, Clone)]
pub struct SecretKey {
    f: Polynomial,
    g: Polynomial,
    capital_f: Polynomial,
    capital_g: Polynomial,
    b0_ntt: Polynomial,
    g0_ntt: Polynomial,
    t_ntt: Polynomial,
}

impl Display for SecretKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "Falcon-{} secret key:\n",
            self.f.coefficients.len()
        ))?;
        f.write_fmt(format_args!(
            "f = {}\n",
            self.f.coefficients.iter().join(", ")
        ))?;
        f.write_fmt(format_args!(
            "g = {}\n",
            self.g.coefficients.iter().join(", ")
        ))?;
        f.write_fmt(format_args!(
            "F = {}\n",
            self.capital_f.coefficients.iter().join(", ")
        ))?;
        f.write_fmt(format_args!(
            "G = {}\n",
            self.capital_g.coefficients.iter().join(", ")
        ))?;
        Ok(())
    }
}

impl SecretKey {
    pub fn generate(variant: FalconVariant) -> Self {
        // According to the docs [1], `thread_rng` uses entropy supplied
        // by the operating system and ChaCha12 to extend it. So it is
        // cryptographically secure, afaict.
        // [1]: https://rust-random.github.io/rand/rand/rngs/struct.ThreadRng.html
        Self::generate_from_seed(variant, thread_rng().gen())
    }

    pub fn generate_from_seed(variant: FalconVariant, seed: [u8; 32]) -> Self {
        let params = Params::new(variant);
        // let (f, g, capital_f, capital_g) = ntru_gen(params.n, seed);
        // let b0 = [g, -f, capital_g, -capital_f];
        // let b0_ntt = b0.iter().map(|b| ntt(&b.coefficients)).collect();
        // let g0 = gram(b0);
        // let g0_ntt = g0.iter().map(|g| ntt(g)).collect();
        // let mut t_ntt = ffldl(g0_ntt);
        // normalize_tree(&mut t_ntt, params.sigma);

        // SecretKey {
        //     f,
        //     g,
        //     capital_f,
        //     capital_g,
        //     b0_ntt,
        //     g0_ntt,
        //     t_ntt,
        // }
        todo!()
    }
}

#[derive(Debug, Clone)]
pub struct PublicKey {
    h: Vec<Felt>,
}

impl PublicKey {
    pub fn from_secret_key(sk: &SecretKey) -> Self {
        let f_ntt = ntt(&sk.f.coefficients);
        let g_ntt = ntt(&sk.g.coefficients);
        let g_inv = Felt::batch_inverse_or_zero(&g_ntt);
        let h_ntt = f_ntt
            .into_iter()
            .zip_eq(g_inv.into_iter())
            .map(|(a, b)| a * b)
            .collect_vec();
        let h = intt(&h_ntt);
        Self {
            h: h.try_into().unwrap(),
        }
    }
}
