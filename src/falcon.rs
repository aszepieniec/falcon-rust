use std::fmt::Display;

use itertools::Itertools;
use num_complex::Complex64;
use rand::{rngs::StdRng, thread_rng, Rng, RngCore, SeedableRng};
use rand_distr::num_traits::Zero;

use crate::{
    ffsampling::{ffldl, gram, normalize_tree, LdlTree},
    fft::fft,
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

fn gen_poly(n: usize, sigma_min: f64, rng: &mut dyn RngCore) -> Polynomial<i16> {
    let mu = 0.0;
    let sigma_star = 1.43300980528773;
    Polynomial {
        coefficients: (0..4096)
            .map(|_| sampler_z(mu, sigma_star, sigma_min, rng))
            .collect_vec()
            .chunks(4096 / n)
            .map(|ch| ch.iter().sum())
            .collect_vec(),
    }
}

/// Sample 4 small polynomials f, g, F, G such that
/// f * G - g * F = q mod (X^n + 1). This function implements
/// algorithm 5 (NTRUgen) on page 33 of the documentation [1].
/// [1]: https://falcon-sign.info/falcon.pdf
fn ntru_gen(
    n: usize,
    seed: [u8; 32],
) -> (
    Polynomial<i16>,
    Polynomial<i16>,
    Polynomial<i16>,
    Polynomial<i16>,
) {
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let sigma_fg = 1.17 * ((Q as f64) / (n as f64)).sqrt();
    loop {
        let f = gen_poly(n, sigma_fg, &mut rng);
        let g = gen_poly(n, sigma_fg, &mut rng);
        let f_ntt = ntt(&f
            .coefficients
            .iter()
            .cloned()
            .map(Felt::new)
            .collect::<Vec<_>>());
        if f_ntt.iter().any(|e| e.is_zero()) {
            continue;
        }
        let norm_f_g: i32 = f
            .coefficients
            .iter()
            .chain(f.coefficients.iter())
            .cloned()
            .map(|c| c as i32)
            .map(|c| c * c)
            .sum();
        let f_star = f.hermitian_adjoint();
        let g_star = g.hermitian_adjoint();
        // if gs_norm(f, g, Q) > 1.3689 * Q {
        //     continue;
        // }
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
    b0_fft: [Vec<Complex64>; 4],
    tree: LdlTree,
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
        let (f, g, capital_f, capital_g) = ntru_gen(params.n, seed);
        let b0 = [g, -f, capital_g, -capital_f];
        let b0_fft = b0
            .map(|v| v.coefficients)
            .map(|c| {
                c.iter()
                    .map(|cc| Complex64::new(*cc as f64, 0.0))
                    .collect_vec()
            })
            .map(|c| fft(&c));
        let g0_fft = gram(b0_fft.clone());
        let mut tree = ffldl(g0_fft);
        normalize_tree(&mut tree, params.sigma);

        SecretKey { b0_fft, tree }
    }
}

#[derive(Debug, Clone)]
pub struct PublicKey {
    h: Vec<Felt>,
}

impl PublicKey {
    pub fn from_secret_key(sk: &SecretKey) -> Self {
        let f_ntt = sk.b0_fft[0]
            .iter()
            .map(|c| Felt::new(c.re as i16))
            .collect_vec();
        let g_ntt = sk.b0_fft[1]
            .iter()
            .map(|c| Felt::new(-c.re as i16))
            .collect_vec();
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
