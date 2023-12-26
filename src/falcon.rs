use itertools::Itertools;
use num_complex::{Complex, Complex64};
use rand::{rngs::StdRng, thread_rng, Rng, RngCore, SeedableRng};
use rand_distr::num_traits::Zero;

use crate::{
    common::split,
    encoding::{compress, decompress},
    ffsampling::{ffldl, ffsampling, gram, normalize_tree, LdlTree},
    fft::{fft, ifft},
    field::{Felt, Q},
    ntt::{intt, ntt},
    polynomial::{hash_to_point, Polynomial},
    samplerz::sampler_z,
};

pub const fn logn(mut n: u32) -> usize {
    let mut ctr = 0;
    while n != 0 {
        ctr += 1;
        n >>= 1;
    }
    ctr
}

/// Bytelength of the signing salt and header
const HEAD_LEN: usize = 1;
const SALT_LEN: usize = 40;
const SEED_LEN: usize = 56;

pub struct SignatureScheme {
    pub n: usize,
    pub sigma: f64,
    pub sigmin: f64,
    pub sig_bound: u64,
    pub sig_bytelen: usize,
}
pub enum FalconVariant {
    Falcon512,
    Falcon1024,
}

impl SignatureScheme {
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

/// Sample 4 small polynomials f, g, F, G such that f * G - g * F = q mod (X^n + 1).
/// Algorithm 5 (NTRUgen) of the documentation [1, p.34].
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
        let norm_f_g = f
            .coefficients
            .iter()
            .chain(g.coefficients.iter())
            .cloned()
            .map(|c| c as i32)
            .map(|c| (c * c) as u64)
            .sum::<u64>();

        let f_fft = fft(&f
            .clone()
            .coefficients
            .clone()
            .iter()
            .map(|a| Complex64::new(*a as f64, 0.0))
            .collect_vec());
        let g_fft = fft(&g
            .clone()
            .coefficients
            .clone()
            .iter()
            .map(|a| Complex64::new(*a as f64, 0.0))
            .collect_vec());
        let f_star_fft = f_fft.iter().map(|&a| a.conj()).collect_vec();
        let g_star_fft = g_fft.iter().map(|&a| a.conj()).collect_vec();
        let ffstar_plus_ggstar_fft = f_fft
            .iter()
            .zip(f_star_fft.iter().zip(g_fft.iter().zip(g_star_fft.iter())))
            .map(|(&a, (&b, (&c, &d)))| a * b + c * d)
            .collect_vec();
        let f_over_d_fft = f_star_fft
            .iter()
            .zip(ffstar_plus_ggstar_fft.iter())
            .map(|(&a, &b)| a / b)
            .collect_vec();
        let g_over_d_fft = g_star_fft
            .iter()
            .zip(ffstar_plus_ggstar_fft.iter())
            .map(|(&a, &b)| a / b)
            .collect_vec();
        let gamma2 = f_over_d_fft
            .iter()
            .chain(g_over_d_fft.iter())
            .map(|&a| (Q as f64) * (a * a.conj()).re)
            .sum::<f64>();
        let gamma = f64::max(norm_f_g as f64, gamma2);
        if gamma > 1.3689f64 * (Q as f64) {
            continue;
        }

        if let Some((capital_f_coefficients, capital_g_coefficients)) = ntru_solve(
            &f.coefficients.iter().map(|&i| i as i32).collect_vec(),
            &g.coefficients.iter().map(|&i| i as i32).collect_vec(),
        ) {
            return (
                f,
                g,
                Polynomial::new(
                    &capital_f_coefficients
                        .iter()
                        .map(|&i| i as i16)
                        .collect_vec(),
                ),
                Polynomial::new(
                    &capital_g_coefficients
                        .iter()
                        .map(|&i| i as i16)
                        .collect_vec(),
                ),
            );
        }
    }
}

/// Extended Euclidean algorithm for computing the greatest common divisor (gcd) and
/// Bézout coefficients (u, v) for the relation
///
///     ua + vb = gcd .   (<-- Bézout relation)
///
/// Implementation adapted from Wikipedia [1].
///
/// [1]: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
fn xgcd(a: i32, b: i32) -> (i32, i32, i32) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);

    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }

    (old_r, old_s, old_t)
}

/// Compute the product of f(X) and g(X) modulo X^n + 1
fn reduced_product(f: &[i32], g: &[i32]) -> Vec<i32> {
    let mut p = vec![0i32; f.len()];
    for i in 0..f.len() {
        for j in 0..g.len() {
            if i + j >= f.len() {
                p[i + j - f.len()] -= f[i] * g[j];
            } else {
                p[i + j] += f[i] * g[j];
            }
        }
    }
    p
}

/// Compute the field norm as required by ntru_solve. Mathematically, let
/// L = QQ[ X ] / < X^n + 1 > and K = QQ[ X ] / < X^{n/2} + 1 > then field_norm maps
///
///     L       -> K
///
///     (a,b)   -> a^2 + x*b^2 .
fn field_norm(f: &[i32]) -> Vec<i32> {
    let (f0, f1) = split(f);
    let mut f0_squared = reduced_product(&f0, &f0);
    let n = f0_squared.len();
    let f1_squared = reduced_product(&f1, &f1);
    for (i, xb2) in f1_squared.into_iter().enumerate() {
        f0_squared[(i + 1) % n] += xb2;
    }
    f0_squared
}

/// Gram-Schmidt orthogonalize the vector (F,G) relative to (f,g).
///
/// Algorithm 7 in the spec [1, p.35]
///
/// [1]: https://falcon-sign.info/falcon.pdf
fn reduce(f: &[i32], g: &[i32], capital_f: &mut [i32], capital_g: &mut [i32]) {
    loop {
        let f_star = Polynomial::new(f).hermitian_adjoint().coefficients;
        let g_star = Polynomial::new(g).hermitian_adjoint().coefficients;
        let ffstar = reduced_product(f, &f_star);
        let ggstar = reduced_product(g, &g_star);
        let ffstar_plus_ggstar = ffstar
            .iter()
            .zip(ggstar.iter())
            .map(|(&a, &b)| a + b)
            .collect_vec();
        let capital_ffstar = reduced_product(capital_f, &f_star);
        let capital_ggstar = reduced_product(capital_g, &g_star);
        let capital_ffstar_plus_capital_ggstar = capital_ffstar
            .into_iter()
            .zip(capital_ggstar)
            .map(|(a, b)| a + b)
            .collect_vec();
        let k = capital_ffstar_plus_capital_ggstar
            .into_iter()
            .zip(ffstar_plus_ggstar)
            .map(|(n, d)| f64::round((n as f64) / (d as f64)) as i32)
            .collect_vec();
        if k == vec![0i32; k.len()] {
            break;
        }
        let k_had_f = k.iter().zip(f.iter()).map(|(&a, &b)| a * b).collect_vec();
        let k_had_g = k.iter().zip(g.iter()).map(|(&a, &b)| a * b).collect_vec();
        for i in 0..k_had_f.len() {
            capital_f[i] -= k_had_f[i];
            capital_g[i] -= k_had_g[i];
        }
    }
}

/// Solve the NTRU equation. Given f, g in ZZ[ X ], find F, G in ZZ[ X ] such that
///
///     fG - gF = q (mod X^n + 1)  (<-- NTRU equation)
///
/// Algorithm 6 of the specification [1, p.35].
///
/// [1]: https://falcon-sign.info/falcon.pdf
fn ntru_solve(f: &[i32], g: &[i32]) -> Option<(Vec<i32>, Vec<i32>)> {
    if f.len() == 1 {
        let (gcd, u, v) = xgcd(f[0], g[0]);
        if gcd != 1 {
            return None;
        }
        return Some(((vec![u * Q]), vec![v * Q]));
    }

    let f_prime = field_norm(f);
    let g_prime = field_norm(g);
    let Some((capital_f_prime, capital_g_prime)) = ntru_solve(&f_prime, &g_prime) else {
        return None;
    };
    let mut capital_f_prime_xsq = vec![0i32; f.len()];
    let mut capital_g_prime_xsq = vec![0i32; g.len()];
    for i in 0..g.len() / 2 {
        capital_f_prime_xsq[2 * i] = capital_f_prime[i];
        capital_g_prime_xsq[2 * i] = capital_g_prime[i];
    }
    let mut f_minx = f.to_vec();
    let mut g_minx = g.to_vec();
    for i in 0..g.len() {
        if i % 2 == 1 {
            f_minx[i] *= -1;
            g_minx[i] *= -1;
        }
    }

    let mut capital_f = reduced_product(&capital_f_prime_xsq, &g_minx);
    let mut capital_g = reduced_product(&capital_g_prime_xsq, &f_minx);
    reduce(f, g, &mut capital_f, &mut capital_g);

    Some((capital_f, capital_g))
}

#[derive(Debug, Clone)]
pub struct SecretKey {
    b0_fft: [Vec<Complex64>; 4],
    tree: LdlTree,
}

impl SecretKey {
    pub fn generate(scheme: &SignatureScheme) -> Self {
        // According to the docs [1], `thread_rng` uses entropy supplied
        // by the operating system and ChaCha12 to extend it. So it is
        // cryptographically secure, afaict.
        // [1]: https://rust-random.github.io/rand/rand/rngs/struct.ThreadRng.html
        Self::generate_from_seed(scheme, thread_rng().gen())
    }

    pub fn generate_from_seed(params: &SignatureScheme, seed: [u8; 32]) -> Self {
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
            .zip_eq(g_inv)
            .map(|(a, b)| a * b)
            .collect_vec();
        let h = intt(&h_ntt);
        Self { h }
    }
}

pub struct Signature {
    r: [u8; 40],
    s: Vec<u8>,
}

impl SignatureScheme {
    // Generate a key pair from a seed.
    pub fn keygen(&self, seed: [u8; 32]) -> (SecretKey, PublicKey) {
        let sk = SecretKey::generate_from_seed(self, seed);
        let pk = PublicKey::from_secret_key(&sk);
        (sk, pk)
    }

    /// Sign a message with the secret key. Algorithm 10 of the specification [1, p.39].
    ///
    /// [1]: https://falcon-sign.info/falcon.pdf
    pub fn sign(&self, m: &[u8], sk: &SecretKey) -> Signature {
        let bound = self.sig_bound;
        let n = self.n;

        let mut rng = thread_rng();
        let mut r = [0u8; 40];
        rng.fill_bytes(&mut r);
        let r_cat_m = [r.to_vec(), m.to_vec()].concat();

        let c = hash_to_point(&r_cat_m, n);
        let one_over_q = 1.0 / (Q as f64);
        let c_over_q_fft = fft(&c
            .coefficients
            .iter()
            .map(|cc| Complex::new(one_over_q * cc.0 as f64, 0.0))
            .collect_vec());

        // B = [[FFT(g), -FFT(f)], [FFT(G), -FFT(F)]]
        let t0 = c_over_q_fft
            .iter()
            .zip(sk.b0_fft[3].iter())
            .map(|(a, b)| a * b)
            .collect_vec();
        let t1 = c_over_q_fft
            .iter()
            .zip(sk.b0_fft[1].iter())
            .map(|(a, b)| -a * b)
            .collect_vec();

        let s = loop {
            let bold_s = loop {
                let z = ffsampling(&(t0.clone(), t1.clone()), &sk.tree, self, &mut rng);
                let t0_min_z0 = t0.iter().zip(z.0).map(|(t, z)| t - z).collect_vec();
                let t1_min_z1 = t1.iter().zip(z.1).map(|(t, z)| t - z).collect_vec();

                // s = (t-z) * B
                let s0 = t0_min_z0
                    .iter()
                    .zip(
                        sk.b0_fft[0]
                            .iter()
                            .zip(t1_min_z1.iter().zip(sk.b0_fft[2].iter())),
                    )
                    .map(|(a, (b, (c, d)))| a * b + c * d)
                    .collect_vec();
                let s1 = t0_min_z0
                    .iter()
                    .zip(
                        sk.b0_fft[1]
                            .iter()
                            .zip(t1_min_z1.iter().zip(sk.b0_fft[3].iter())),
                    )
                    .map(|(a, (b, (c, d)))| a * b + c * d)
                    .collect_vec();

                let length_squared: f64 = s0.iter().map(|a| (a * a.conj()).re).sum::<f64>()
                    + s1.iter().map(|a| (a * a.conj()).re).sum::<f64>();

                if length_squared > (n as u64 * bound) as f64 {
                    continue;
                }

                break [s0, s1].concat();
            };
            let s2: Vec<Complex64> = ifft(&bold_s)[n / 2..].to_vec();
            let maybe_s = compress(
                &s2.iter().map(|a| a.re as i16).collect_vec(),
                8 * self.sig_bytelen - 328,
            );

            match maybe_s {
                Some(s) => {
                    break s;
                }
                None => {
                    continue;
                }
            };
        };

        Signature { r, s }
    }

    /// Verify a signature. Algorithm 16 in the spec [1, p.45].
    /// [1]: https://falcon-sign.info/falcon.pdf
    pub fn verify(&self, m: &[u8], sig: &Signature, pk: &PublicKey) -> bool {
        let n = self.n;
        let r_cat_m = [sig.r.to_vec(), m.to_vec()].concat();
        let c = hash_to_point(&r_cat_m, n);

        if 8 * self.sig_bytelen - 328 != sig.s.len() {
            return false;
        }
        let s2 = match decompress(&sig.s, n) {
            Some(success) => success,
            None => {
                return false;
            }
        };
        let s2_ntt = ntt(&s2.iter().map(|a| Felt(*a)).collect_vec());
        let h_ntt = ntt(&pk.h);
        let c_ntt = ntt(&c.coefficients);

        // s1 = c - s2 * pk.h;
        let s1_ntt = c_ntt
            .iter()
            .zip(s2_ntt.iter().zip(h_ntt.iter()))
            .map(|(&a, (&b, &c))| a - b * c)
            .collect_vec();
        let s1 = intt(&s1_ntt);

        let length_squared = s1.iter().map(|&i| (i.0 * i.0) as u64).sum::<u64>()
            + s2.iter().map(|&i| (i * i) as u64).sum::<u64>();
        length_squared < self.sig_bound
    }
}

#[cfg(test)]
mod test {
    use rand::{thread_rng, Rng, RngCore};

    use super::SignatureScheme;

    #[test]
    fn test_operation() {
        let mut rng = thread_rng();
        let mut msg = [0u8; 5];
        rng.fill_bytes(&mut msg);

        let small_scheme = SignatureScheme::new(super::FalconVariant::Falcon512);
        let (sk, pk) = small_scheme.keygen(rng.gen());
        let sig = small_scheme.sign(&msg, &sk);
        assert!(small_scheme.verify(&msg, &sig, &pk));

        let big_scheme = SignatureScheme::new(super::FalconVariant::Falcon1024);
        let (sk, pk) = big_scheme.keygen(rng.gen());
        let sig = big_scheme.sign(&msg, &sk);
        assert!(big_scheme.verify(&msg, &sig, &pk));
    }
}
