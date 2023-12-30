use bit_vec::BitVec;
use itertools::Itertools;
use num::{BigInt, FromPrimitive, One};
use num_complex::{Complex, Complex64};
use rand::{rngs::StdRng, thread_rng, Rng, RngCore, SeedableRng};
use rand_distr::num_traits::Zero;

use crate::{
    encoding::{compress, decompress},
    ffsampling::{ffldl, ffsampling, gram, normalize_tree, LdlTree},
    fft::{fft, ifft},
    field::{Felt, Q},
    ntt::{intt, ntt},
    polynomial::{hash_to_point, Polynomial},
    samplerz::sampler_z,
};

#[derive(Copy, Clone, Debug)]
pub struct FalconParameters {
    pub(crate) n: usize,
    pub(crate) sigma: f64,
    pub(crate) sigmin: f64,
    pub(crate) sig_bound: i64,
    pub(crate) sig_bytelen: usize,
}

pub enum FalconVariant {
    Falcon512,
    Falcon1024,
}

impl FalconVariant {
    const fn from_n(n: usize) -> Self {
        match n {
            512 => Self::Falcon512,
            1024 => Self::Falcon1024,
            _ => unreachable!(),
        }
    }
    pub(crate) const fn parameters(&self) -> FalconParameters {
        match self {
            FalconVariant::Falcon512 => FalconParameters {
                n: 512,
                sigma: 165.7366171829776,
                sigmin: 1.2778336969128337,
                sig_bound: 34034726,
                sig_bytelen: 666,
            },
            FalconVariant::Falcon1024 => FalconParameters {
                n: 1024,
                sigma: 168.38857144654395,
                sigmin: 1.298280334344292,
                sig_bound: 70265242,
                sig_bytelen: 1280,
            },
        }
    }
}

/// Generate a polynomial of degree at most n-1 whose coefficients are
/// distributed according to a discrete Gaussian with mu = 0 and
/// sigma = 1.17 * sqrt(Q / (2n)).
fn gen_poly(n: usize, rng: &mut dyn RngCore) -> Polynomial<i16> {
    let mu = 0.0;
    let sigma_star = 1.43300980528773;
    Polynomial {
        coefficients: (0..4096)
            .map(|_| sampler_z(mu, sigma_star, sigma_star - 0.001, rng))
            .collect_vec()
            .chunks(4096 / n)
            .map(|ch| ch.iter().sum())
            .collect_vec(),
    }
}

/// Compute the Gram-Schmidt norm of B = [[g, -f], [G, -F]] from f and g.
/// Corresponds to line 9 in algorithm 5 of the spec [1, p.34]
///
/// [1]: https://falcon-sign.info/falcon.pdf
fn gram_schmidt_norm(f: &Polynomial<i16>, g: &Polynomial<i16>) -> f64 {
    let n = f.coefficients.len();
    let norm_f = f.l2_norm();
    let norm_g = g.l2_norm();
    let sqnorm = norm_f * norm_f + norm_g * norm_g;
    let gamma1 = f64::sqrt(sqnorm);

    let f_adj = f.hermitian_adjoint();
    let g_adj = g.hermitian_adjoint();
    let ffgg = (f.clone() * f_adj.clone() + g.clone() * g_adj.clone()).reduce_by_cyclotomic(n);
    let ffgg_float = ffgg.map(|c| *c as f64);
    let ffgginv = ffgg_float.approximate_cyclotomic_ring_inverse(n);
    let qf_over_ffgg =
        (f_adj.map(|c| (*c as f64)) * (Q as f64) * ffgginv.clone()).reduce_by_cyclotomic(n);
    let qg_over_ffgg = (g_adj.map(|c| (*c as f64)) * (Q as f64) * ffgginv).reduce_by_cyclotomic(n);
    let norm_f_over_ffgg = qf_over_ffgg.l2_norm();
    let norm_g_over_ffgg = qg_over_ffgg.l2_norm();

    let gamma2 =
        f64::sqrt(norm_f_over_ffgg * norm_f_over_ffgg + norm_g_over_ffgg * norm_g_over_ffgg);

    f64::max(gamma1, gamma2)
}

/// Sample 4 small polynomials f, g, F, G such that f * G - g * F = q mod (X^n + 1).
/// Algorithm 5 (NTRUgen) of the documentation [1, p.34].
///
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
    loop {
        let f = gen_poly(n, &mut rng);
        let g = gen_poly(n, &mut rng);
        let f_ntt = ntt(&f
            .coefficients
            .iter()
            .cloned()
            .map(Felt::new)
            .collect::<Vec<_>>());
        if f_ntt.iter().any(|e| e.is_zero()) {
            continue;
        }
        let gamma = gram_schmidt_norm(&f, &g);
        if gamma * gamma > 1.3689f64 * (Q as f64) {
            continue;
        }

        if let Some((capital_f, capital_g)) =
            ntru_solve(&f.map(|&i| i.into()), &g.map(|&i| i.into()))
        {
            return (
                f,
                g,
                capital_f.map(|i| i.try_into().unwrap()),
                capital_g.map(|i| i.try_into().unwrap()),
            );
        }
    }
}

/// Extended Euclidean algorithm for computing the greatest common divisor (g) and
/// Bézout coefficients (u, v) for the relation
///
/// $$ u a + v b = g . $$
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

/// Reduce the vector (F,G) relative to (f,g). This method
/// follows the python implementation [1].
///
/// Algorithm 7 in the spec [2, p.35]
///
/// [1]: https://github.com/tprest/falcon.py
///
/// [2]: https://falcon-sign.info/falcon.pdf
fn babai_reduce(
    f: &Polynomial<BigInt>,
    g: &Polynomial<BigInt>,
    capital_f: &mut Polynomial<BigInt>,
    capital_g: &mut Polynomial<BigInt>,
) {
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
    let f_adjusted = f.map(|bi| i64::try_from(bi >> shift).unwrap() as f64);
    let g_adjusted = g.map(|bi| i64::try_from(bi >> shift).unwrap() as f64);
    let f_adjusted_fft = f_adjusted.fft();
    let g_adjusted_fft = g_adjusted.fft();

    let f_star_adjusted_fft = f_adjusted_fft.map(|c| c.conj());
    let g_star_adjusted_fft = g_adjusted_fft.map(|c| c.conj());
    let denominator_fft = f_adjusted_fft.hadamard_mul(&f_star_adjusted_fft)
        + g_adjusted_fft.hadamard_mul(&g_star_adjusted_fft);

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
        let capital_f_adjusted =
            capital_f.map(|bi| i64::try_from(bi >> capital_shift).unwrap() as f64);
        let capital_g_adjusted =
            capital_g.map(|bi| i64::try_from(bi >> capital_shift).unwrap() as f64);

        let capital_f_adjusted_fft = capital_f_adjusted.fft();
        let capital_g_adjusted_fft = capital_g_adjusted.fft();
        let numerator_fft = capital_f_adjusted_fft.hadamard_mul(&f_star_adjusted_fft)
            + capital_g_adjusted_fft.hadamard_mul(&g_star_adjusted_fft);
        let quotient_fft = numerator_fft.hadamard_div(&denominator_fft);
        let quotient = quotient_fft.ifft();

        let k = quotient.map(|f| Into::<BigInt>::into(f.round() as i64));

        if k.is_zero() {
            break;
        }
        let kf = (k.clone() * f.clone())
            .reduce_by_cyclotomic(n)
            .map(|bi| bi << (capital_size - size));
        let kg = (k.clone() * g.clone())
            .reduce_by_cyclotomic(n)
            .map(|bi| bi << (capital_size - size));
        *capital_f -= kf;
        *capital_g -= kg;

        counter += 1;
        if counter > 10000 {
            panic!("Should not have to do more than 10000 iterations! This probably indicates an infinite loop. Please file a bug report, for instance by opening an issue on https://github.com/aszepieniec/falcon-rust/. Debug information:\nf: {:?}\ng: {:?}\ncapital_f: {:?}\ncapital_g: {:?}", f.coefficients, g.coefficients, capital_f.coefficients, capital_g.coefficients);
        }
    }
}

/// Solve the NTRU equation. Given $f, g \in \mathbb{Z}[X]$, find $F, G \in \mathbb{Z}[X]$
/// such that
///
/// $$    f G - g F = q \mod \langle X^n + 1 \rangle $$
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
            (Polynomial::new(vec![-v * BigInt::from_i32(Q).unwrap()])),
            Polynomial::new(vec![u * BigInt::from_i32(Q).unwrap()]),
        ));
    }

    let f_prime = f.field_norm();
    let g_prime = g.field_norm();
    let Some((capital_f_prime, capital_g_prime)) = ntru_solve(&f_prime, &g_prime) else {
        return None;
    };
    let capital_f_prime_xsq = capital_f_prime.lift_next_cyclotomic();
    let capital_g_prime_xsq = capital_g_prime.lift_next_cyclotomic();
    let f_minx = f.galois_adjoint();
    let g_minx = g.galois_adjoint();

    let mut capital_f = (capital_f_prime_xsq * g_minx).reduce_by_cyclotomic(n);
    let mut capital_g = (capital_g_prime_xsq * f_minx).reduce_by_cyclotomic(n);

    babai_reduce(f, g, &mut capital_f, &mut capital_g);

    Some((capital_f, capital_g))
}

#[derive(Debug)]
pub enum FalconDeserializationError {
    CannotDetermineFieldElementEncodingMethod,
    CannotInferFalconVariant,
    InvalidHeaderFormat,
    InvalidLogN,
    BadEncodingLength,
    BadFieldElementEncoding,
    WrongVariant,
}

#[derive(Debug, Clone)]
pub struct SecretKey<const N: usize> {
    /// b0_fft = [[g_fft, -f_fft], [G_fft, -F_fft]]
    b0_fft: [Vec<Complex64>; 4],
    tree: LdlTree,
}

impl<const N: usize> SecretKey<N> {
    /// Generate a secret key using randomness supplied by the operating system.
    pub fn generate() -> Self {
        // According to the docs [1], `thread_rng` uses entropy supplied
        // by the operating system and ChaCha12 to extend it. So it is
        // cryptographically secure, afaict.
        // [1]: https://rust-random.github.io/rand/rand/rngs/struct.ThreadRng.html
        Self::generate_from_seed(thread_rng().gen())
    }

    /// Generate a secret key pseudorandomly by expanding a given seed.
    pub fn generate_from_seed(seed: [u8; 32]) -> Self {
        // separate sk gen for testing purposes
        let b0 = Self::gen_b0(seed);
        Self::from_b0(b0)
    }

    pub(crate) fn gen_b0(seed: [u8; 32]) -> [Polynomial<i16>; 4] {
        let (f, g, capital_f, capital_g) = ntru_gen(N, seed);
        [g, -f, capital_g, -capital_f]
    }

    pub(crate) fn from_b0(b0: [Polynomial<i16>; 4]) -> Self {
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
        let sigma = FalconVariant::from_n(N).parameters().sigma;
        normalize_tree(&mut tree, sigma);

        SecretKey { b0_fft, tree }
    }

    fn field_element_width(n: usize, polynomial_index: usize) -> usize {
        if polynomial_index == 2 {
            8
        } else {
            match n {
                1024 => 5,
                512 => 6,
                _ => unreachable!(),
            }
        }
    }

    fn serialize_field_element(element_width: usize, element: Felt) -> BitVec {
        let mut bits = BitVec::new();
        let int = element.balanced_value();
        for i in (0..element_width).rev() {
            bits.push(int & (1i16 << i) != 0);
        }
        bits
    }

    fn deserialize_field_element(bits: &BitVec) -> Result<Felt, FalconDeserializationError> {
        if bits[0] && bits.iter().skip(1).all(|b| !b) {
            return Err(FalconDeserializationError::BadFieldElementEncoding);
        }

        let mut uint = 0;
        for bit in bits {
            uint = (uint << 1) | (bit as i16);
        }
        if bits[0] {
            uint = (uint << (16 - bits.len())) >> (16 - bits.len());
        }
        Ok(Felt::new(uint))
    }

    /// Serialize the secret key to a vector of bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        // header
        let n = self.b0_fft[0].len();
        let l = n.checked_ilog2().unwrap() as u8;
        let header: u8 = (5 << 4) // fixed bits
                        | l;

        let f = ifft(&self.b0_fft[1])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();
        let g = ifft(&self.b0_fft[0])
            .into_iter()
            .map(|c| Felt::new(c.re.round() as i16))
            .collect_vec();
        let capital_f = ifft(&self.b0_fft[3])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();

        let mut bits = BitVec::from_bytes(&[header]);
        // f
        let width = Self::field_element_width(n, 0);
        for fi in f.into_iter() {
            let mut substring = Self::serialize_field_element(width, fi);
            bits.append(&mut substring);
        }
        // g
        let width = Self::field_element_width(n, 1);
        for fi in g.into_iter() {
            let mut substring = Self::serialize_field_element(width, fi);
            bits.append(&mut substring);
        }
        // capital_f
        let width = Self::field_element_width(n, 2);
        for fi in capital_f.into_iter() {
            let mut substring = Self::serialize_field_element(width, fi);
            bits.append(&mut substring);
        }

        bits.to_bytes()
    }

    /// Deserialize a secret key from a slice of bytes.
    pub fn from_bytes(byte_vector: &[u8]) -> Result<Self, FalconDeserializationError> {
        // check length
        if byte_vector.len() < 2 {
            return Err(FalconDeserializationError::BadEncodingLength);
        }

        // read fields
        let header = byte_vector[0];
        let bit_buffer = BitVec::from_bytes(&byte_vector[1..]);

        // check fixed bits in header
        if (header >> 4) != 5 {
            return Err(FalconDeserializationError::InvalidHeaderFormat);
        }

        // check log n
        let logn = (header & 15) as usize;
        let n = match logn {
            9 => 512,
            10 => 1024,
            _ => return Err(FalconDeserializationError::InvalidLogN),
        };

        // match against const variant generic parameter
        if n != FalconVariant::from_n(N).parameters().n {
            return Err(FalconDeserializationError::WrongVariant);
        }

        // decode integer polynomials using BitVec as intermediate representation

        // f
        let width_f = Self::field_element_width(n, 0);
        let f = bit_buffer
            .iter()
            .take(n * width_f)
            .chunks(width_f)
            .into_iter()
            .map(BitVec::from_iter)
            .map(|subs| Self::deserialize_field_element(&subs))
            .collect::<Result<Vec<Felt>, _>>()?;

        // g
        let width_g = Self::field_element_width(n, 1);
        let g = bit_buffer
            .iter()
            .skip(n * width_f)
            .take(n * width_g)
            .chunks(width_g)
            .into_iter()
            .map(BitVec::from_iter)
            .map(|subs| Self::deserialize_field_element(&subs))
            .collect::<Result<Vec<Felt>, _>>()?;

        // capital_f
        let width_capital_f = Self::field_element_width(n, 2);
        let capital_f = bit_buffer
            .iter()
            .skip(n * width_g + n * width_f)
            .take(n * width_capital_f)
            .chunks(width_capital_f)
            .into_iter()
            .map(BitVec::from_iter)
            .map(|subs| Self::deserialize_field_element(&subs))
            .collect::<Result<Vec<Felt>, _>>()?;

        // all bits in the bit buffer should have been read at this point
        if bit_buffer.len() != n * width_f + n * width_g + n * width_capital_f {
            return Err(FalconDeserializationError::BadEncodingLength);
        }

        // compute capital_g from f, g, capital_f
        let f_ntt = ntt(&f);
        let g_ntt = ntt(&g);
        let capital_f_ntt = ntt(&capital_f);
        let capital_g_ntt = g_ntt
            .into_iter()
            .zip(capital_f_ntt)
            .zip(f_ntt)
            .map(|((g, cf), f)| g * cf / f)
            .collect_vec();
        let capital_g = intt(&capital_g_ntt);

        Ok(Self::from_b0([
            Polynomial::new(g.to_vec()).map(|f| f.balanced_value()),
            -Polynomial::new(f.to_vec()).map(|f| f.balanced_value()),
            Polynomial::new(capital_g.to_vec()).map(|f| f.balanced_value()),
            -Polynomial::new(capital_f.to_vec()).map(|f| f.balanced_value()),
        ]))
    }
}

impl<const N: usize> PartialEq for SecretKey<N> {
    fn eq(&self, other: &Self) -> bool {
        let own_f = ifft(&self.b0_fft[1])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();
        let own_g = ifft(&self.b0_fft[0])
            .into_iter()
            .map(|c| Felt::new(c.re.round() as i16))
            .collect_vec();
        let own_capital_f = ifft(&self.b0_fft[3])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();
        let own_capital_g = ifft(&self.b0_fft[2])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();

        let other_f = ifft(&other.b0_fft[1])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();
        let other_g = ifft(&other.b0_fft[0])
            .into_iter()
            .map(|c| Felt::new(c.re.round() as i16))
            .collect_vec();
        let other_capital_f = ifft(&other.b0_fft[3])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();
        let other_capital_g = ifft(&other.b0_fft[2])
            .into_iter()
            .map(|c| Felt::new(-c.re.round() as i16))
            .collect_vec();

        own_f == other_f
            && own_g == other_g
            && own_capital_f == other_capital_f
            && own_capital_g == other_capital_g
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PublicKey<const N: usize> {
    h: Vec<Felt>,
}

impl<const N: usize> PublicKey<N> {
    /// Compute the public key that matches with this secret key.
    pub fn from_secret_key(sk: &SecretKey<N>) -> Self {
        let f = ifft(&sk.b0_fft[1])
            .iter()
            .map(|c| -Felt::new(c.re.round() as i16))
            .collect_vec();
        let f_ntt = ntt(&f);
        let g = ifft(&sk.b0_fft[0])
            .iter()
            .map(|c| Felt::new(c.re.round() as i16))
            .collect_vec();
        let g_ntt = ntt(&g);
        let f_inv = Felt::batch_inverse_or_zero(&f_ntt);
        let h_ntt = g_ntt
            .clone()
            .into_iter()
            .zip_eq(f_inv)
            .map(|(a, b)| a * b)
            .collect_vec();
        let h = intt(&h_ntt);
        Self { h }
    }

    /// Deserialize the given slice of bytes into a public key.
    pub fn from_bytes(byte_array: &[u8]) -> Result<Self, FalconDeserializationError> {
        let n: usize = match byte_array.len() {
            897 => 512,
            1793 => 1024,
            _ => return Err(FalconDeserializationError::BadEncodingLength),
        };

        // match against variant generic type parameter
        if n != N {
            return Err(FalconDeserializationError::WrongVariant);
        }

        // parse header
        let header = byte_array[0];

        if header >> 4 != 0 {
            return Err(FalconDeserializationError::InvalidHeaderFormat);
        }

        let l = n.ilog2();
        if header != l as u8 {
            return Err(FalconDeserializationError::InvalidLogN);
        }

        // parse h
        let bit_buffer = BitVec::from_bytes(&byte_array[1..]);
        let h = bit_buffer
            .iter()
            .chunks(14)
            .into_iter()
            .map(|ch| {
                let mut int = 0;
                for b in ch {
                    int = (int << 1) | (b as i16);
                }
                int
            })
            .map(Felt::new)
            .collect_vec();

        Ok(PublicKey { h })
    }

    // Serialize the public key as a list of bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let header = self.h.len().ilog2() as u8;
        let mut bit_buffer = BitVec::from_bytes(&[header]);

        for hi in self.h.iter() {
            for i in (0..14).rev() {
                bit_buffer.push(hi.value() & (1 << i) != 0);
            }
        }

        bit_buffer.to_bytes()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Signature<const N: usize> {
    r: [u8; 40],
    s: Vec<u8>,
}

impl<const N: usize> Signature<N> {
    /// Serialize the signature to a vector of bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        // header
        let felt_encoding = 2; // standard (compressed)
        let n = self.s.len();
        let l = n.checked_ilog2().unwrap() as u8;
        let header: u8 = (felt_encoding << 5)
                        | (1 << 4) // fixed bit
                        | l;

        vec![header]
            .into_iter()
            .chain(self.r)
            .chain(self.s.iter().cloned())
            .collect_vec()
    }

    /// Deserialize a signature from a slice of bytes.
    pub fn from_bytes(byte_vector: &[u8]) -> Result<Self, FalconDeserializationError> {
        // check signature length; infer variant
        let n = if byte_vector.len() == FalconVariant::Falcon512.parameters().sig_bytelen {
            512
        } else if byte_vector.len() == FalconVariant::Falcon1024.parameters().sig_bytelen {
            1024
        } else {
            return Err(FalconDeserializationError::CannotInferFalconVariant);
        };

        // match n against const type parameter
        if n != N {
            return Err(FalconDeserializationError::WrongVariant);
        }

        // read fields
        let header = byte_vector[0];
        let salt: [u8; 40] = byte_vector[1..=40].try_into().unwrap();
        let signature_vector = &byte_vector[41..];

        // check encoding and reject if not standard
        let felt_encoding: u8 = 2; // standard
        if (header >> 5) & 3 != felt_encoding {
            return Err(FalconDeserializationError::CannotDetermineFieldElementEncodingMethod);
        }

        // check fixed bits in header
        if (header >> 7) != 0 || ((header >> 4) & 1) == 0 {
            return Err(FalconDeserializationError::InvalidHeaderFormat);
        }

        // check log n
        let logn = (header & 15) as usize;
        if n != (1 << logn) {
            return Err(FalconDeserializationError::InvalidLogN);
        }

        // tests pass; assemble object
        Ok(Signature::<N> {
            r: salt,
            s: signature_vector.to_vec(),
        })
    }
}

// Generate a key pair pseudorandomly by expanding a seed.
pub fn keygen<const N: usize>(seed: [u8; 32]) -> (SecretKey<N>, PublicKey<N>) {
    let sk = SecretKey::generate_from_seed(seed);
    let pk = PublicKey::from_secret_key(&sk);
    (sk, pk)
}

/// Sign a message with the secret key.
///
/// Algorithm 10 of the specification [1, p.39].
///
/// [1]: https://falcon-sign.info/falcon.pdf
pub fn sign<const N: usize>(m: &[u8], sk: &SecretKey<N>) -> Signature<N> {
    let mut rng = thread_rng();
    let mut r = [0u8; 40];
    rng.fill_bytes(&mut r);

    let params = FalconVariant::from_n(N).parameters();
    let bound = params.sig_bound;
    let n = params.n;

    let r_cat_m = [r.to_vec(), m.to_vec()].concat();

    let c = hash_to_point(&r_cat_m, n);
    let one_over_q = 1.0 / (Q as f64);
    let c_over_q_fft = fft(&c
        .coefficients
        .iter()
        .map(|cc| Complex::new(one_over_q * cc.value() as f64, 0.0))
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
            let z = ffsampling(&(t0.clone(), t1.clone()), &sk.tree, &params, &mut rng);
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

            // compute the norm of (s0||s1) and note that they are in FFT representation
            let length_squared: f64 = (s0.iter().map(|a| (a * a.conj()).re).sum::<f64>()
                + s1.iter().map(|a| (a * a.conj()).re).sum::<f64>())
                / (n as f64);

            if length_squared > (bound as f64) {
                continue;
            }

            break [s0, s1];
        };
        let s2: Vec<Complex64> = ifft(&bold_s[1]).to_vec();
        let maybe_s = compress(
            &s2.iter().map(|a| a.re.round() as i16).collect_vec(),
            8 * (params.sig_bytelen - 41),
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
///
/// [1]: https://falcon-sign.info/falcon.pdf
pub fn verify<const N: usize>(m: &[u8], sig: &Signature<N>, pk: &PublicKey<N>) -> bool {
    let n = N;
    let params = FalconVariant::from_n(N).parameters();
    let r_cat_m = [sig.r.to_vec(), m.to_vec()].concat();
    let c = hash_to_point(&r_cat_m, n);

    let s2 = match decompress(&sig.s, (params.sig_bytelen - 41) * 8, n) {
        Some(success) => success,
        None => {
            return false;
        }
    };
    let s2_ntt = ntt(&s2.iter().map(|a| Felt::new(*a)).collect_vec());
    let h_ntt = ntt(&pk.h);
    let c_ntt = ntt(&c.coefficients);

    // s1 = c - s2 * pk.h;
    let s1_ntt = c_ntt
        .iter()
        .zip(s2_ntt.iter().zip(h_ntt.iter()))
        .map(|(&a, (&b, &c))| a - b * c)
        .collect_vec();
    let s1 = intt(&s1_ntt);

    let length_squared = s1
        .iter()
        .map(|i| i.balanced_value() as i64)
        .map(|i| (i * i))
        .sum::<i64>()
        + s2.iter().map(|&i| i as i64).map(|i| (i * i)).sum::<i64>();
    length_squared < params.sig_bound
}

#[cfg(test)]
mod test {
    use std::str::FromStr;

    use itertools::Itertools;
    use num::{BigInt, FromPrimitive};
    use rand::{rngs::StdRng, thread_rng, Rng, RngCore, SeedableRng};

    use crate::{
        encoding::compress,
        falcon::{gram_schmidt_norm, keygen, sign, verify, FalconVariant, Signature},
        field::Felt,
        polynomial::{hash_to_point, Polynomial},
    };

    use super::{babai_reduce, gen_poly, ntru_gen, ntru_solve, PublicKey, SecretKey};

    #[test]
    fn test_operation_falcon_512() {
        let mut rng = thread_rng();
        let mut msg = [0u8; 5];
        rng.fill_bytes(&mut msg);

        println!("testing small scheme ...");
        const N: usize = 512;
        println!("-> keygen ...");
        let (sk, pk) = keygen::<N>(rng.gen());
        println!("-> sign ...");
        let sig = sign::<N>(&msg, &sk);
        println!("-> verify ...");
        assert!(verify::<N>(&msg, &sig, &pk));
        println!("-> ok.");
    }

    #[test]
    fn test_stalling_operation_falcon_1024() {
        // let seed: [u8; 32] = [
        //     119, 186, 1, 120, 0, 255, 165, 121, 56, 149, 105, 255, 53, 63, 192, 102, 231, 197, 233,
        //     249, 212, 179, 1, 18, 33, 42, 137, 10, 172, 179, 168, 35,
        // ];
        let seed: [u8; 32] = thread_rng().gen();
        println!("seed: {:2?}", seed);
        let mut rng: StdRng = SeedableRng::from_seed(seed);

        let mut msg = [0u8; 5];
        rng.fill_bytes(&mut msg);
        println!("testing big scheme ...");
        const N: usize = 1024;
        println!("-> keygen ...");
        let (sk, pk) = keygen::<N>(rng.gen());
        println!("-> sign ...");
        let sig = sign::<N>(&msg, &sk);
        println!("-> verify ...");
        assert!(verify::<N>(&msg, &sig, &pk));
        println!("-> ok.");
    }

    #[test]
    fn test_falcon512_test_vector() {
        let nonce = hex::decode(
            "16c12515258093799956368cdfc182c1ca4a34f077e9244416a8c4c13fb0ca241e8b7ac1712d28eb",
        )
        .unwrap();
        let data = hex::decode("6461746131").unwrap();
        let f = vec![
            -4, -2, -5, -1, 4, -2, 0, -3, -1, 1, -2, -2, -6, -3, 3, -5, -1, 4, -3, -8, 4, -1, 2,
            -1, -8, 5, -6, -3, 6, 0, -2, 4, 5, -6, 2, 3, 6, 4, 2, 3, 3, 7, 0, 1, 5, -3, -1, -9, -1,
            6, -2, -5, 4, 0, 4, -2, 10, -4, -3, 4, -7, -1, -7, -2, -1, -6, 5, -1, -9, 3, 2, -5, 4,
            -2, 2, -4, 4, -3, -1, 0, 5, 2, 2, -1, -9, -7, -2, -1, 0, 3, 1, 0, -1, -2, -5, 4, -1,
            -1, 3, -1, 1, 4, -3, 2, -5, -2, 2, -4, 3, 6, 3, 9, 1, -2, 4, -1, -1, -6, -2, -2, 4, 5,
            -1, 0, 10, -2, 1, -2, -3, 0, -4, -4, -1, 0, 1, -5, -3, -7, -2, -1, 2, -6, 3, 0, 0, 4,
            -4, 0, 0, -5, -2, 5, -8, 8, 5, 4, 10, -4, 3, 8, 5, 1, -7, 0, -5, 0, -4, 3, -4, -2, 2,
            -2, 6, 8, 2, -1, 4, -4, -2, 1, 0, 3, 7, 0, 9, -3, 1, 4, -3, 2, -1, 5, -8, 4, -1, 1, -8,
            2, 4, -9, -3, 1, 3, -1, -7, 5, 5, 4, -3, 0, -7, -3, -1, -6, -7, 0, -3, 0, 3, -3, 0, -3,
            1, 3, 4, -6, -6, -3, 6, 0, 2, -5, 1, -3, -6, -6, -1, -7, -2, -4, 3, 0, -4, -1, 2, 7,
            -7, -2, 4, 2, 0, 1, -1, -3, 2, 1, 8, -1, 1, -2, 1, -1, 1, 4, 0, -4, 4, 3, -2, 6, -3,
            -2, 1, 2, 3, 6, 5, -4, -7, -6, 4, 3, -4, 3, -3, 3, -3, 2, -1, 1, 5, -2, 2, 1, 0, -7, 0,
            0, -1, 4, -3, 2, 1, -3, 5, 4, -6, -1, -3, 2, -1, -8, 4, 2, 4, 0, 1, -5, 8, 5, 4, -3,
            -1, -2, 4, 0, 2, -2, 0, -2, -1, -7, 5, 0, 1, 2, 1, -2, 2, -1, 1, -4, 1, 0, 4, -4, 0, 5,
            1, 4, -5, -2, -3, -2, 1, 3, 1, 2, 5, 12, 0, -1, 4, -6, 1, -4, 3, -5, -4, 4, 2, -2, -6,
            1, 1, 3, -1, 0, -4, -4, -4, 6, -2, 4, -3, 0, -2, -1, 0, -6, -3, -2, 0, 6, 5, -5, -5, 3,
            0, 3, -3, -2, 5, 7, -3, 1, -1, 0, 3, 0, 3, -7, 2, -4, -4, 1, 1, 1, 0, -3, -8, 3, 6, 1,
            -2, -7, 3, 3, 4, -1, -2, -5, 9, 7, 1, 2, -4, 4, 0, -11, 3, 0, -3, -5, 5, -1, -1, 7, 6,
            -1, 6, 3, 9, 5, -2, -3, -3, 1, -2, 0, -1, 1, -2, 2, 0, -5, -1, -4, -2, 2, -1, -3, 0,
            -3, 0, 1, 3, -3, 2, 5, 8, -2, 3, -4, -7, 0, 4, -8, 1, 8, -2, 1, -1, 2, 0, -2, 1, 3, 3,
            4, -2, -4, 3, -4, 2, 3, -2, -4, 1, -4, 10, 2,
        ];
        let g = vec![
            -1, 5, -7, -1, -4, 6, 4, -1, -4, -13, -1, -5, -2, -8, 2, 1, 4, 2, 0, 0, 2, 0, -1, 2, 5,
            -5, -8, 8, 1, 11, 0, -8, -4, 1, 1, -6, -4, 1, -3, 0, -10, -4, -6, -3, -2, 1, 6, 2, 8,
            -2, 2, -2, 1, 3, -4, 2, -1, -1, -2, -2, -3, 0, -3, 2, -3, 2, -3, -4, 2, 3, 4, -5, 6,
            -3, -2, -1, -1, -6, -2, 1, -4, -7, 8, 0, 2, -2, 2, 0, 1, 0, 4, 9, 7, 0, -1, -1, 4, -3,
            -2, 6, 6, 0, 1, 7, -6, -5, 5, 1, 4, -1, 0, -2, 3, -4, 1, -1, -3, -2, 0, -1, -7, -8, -1,
            2, 0, -5, 0, 1, -4, 6, -5, 6, 4, 1, -4, -5, 8, -1, 1, -2, 1, 1, 1, 3, 0, -1, 1, 1, -4,
            -5, -4, 2, -3, 2, -2, 3, 7, -4, 4, -1, -2, 4, -4, -5, 2, 6, -7, 5, -1, 1, 3, 0, -5, -5,
            3, -2, -3, -1, -6, 0, 2, 3, 2, 7, -3, -2, -2, 1, -5, 3, 3, -7, 0, 4, 4, -1, 2, -3, 1,
            3, -1, -1, 0, -7, -6, -3, 7, -3, 5, -5, 1, -2, 0, 9, -2, 3, -1, -5, -3, -5, 3, 1, -4,
            -3, 2, -2, 2, 8, -1, 0, 5, -3, -2, -6, 4, 0, 3, -3, -3, 4, -1, 0, 0, -2, -1, 3, 7, 4,
            5, -1, 8, 0, -1, -6, -3, 4, 3, -3, 5, 2, -1, -2, 1, -1, 3, -2, -6, 4, 0, 0, -4, 1, 6,
            2, 0, 10, 9, 2, -2, 0, 2, 1, -3, -1, -1, 3, 2, 1, 1, -3, -2, 7, 2, -1, 5, -3, -2, 1,
            -2, 2, -2, -4, 3, 2, 1, -4, 1, 4, 3, -7, -4, 2, -5, -2, 5, -3, 1, -4, -5, 1, 0, 0, 0,
            7, -5, -1, 2, 2, -3, 6, -6, 4, -3, -5, -6, -7, -4, 3, -2, -2, -10, -3, 2, -1, -6, -4,
            1, 2, 2, 1, 4, 1, -5, -10, -2, 2, -4, 4, 4, -2, 1, 4, -3, 0, -6, -3, 1, 5, -7, -6, -4,
            8, -1, 0, -1, 6, -3, -2, -2, 6, 2, 3, -3, -3, 5, -2, 1, 1, -4, -4, 8, 0, 3, 2, 3, 7, 4,
            3, 2, -6, -9, 0, -8, 11, -2, 2, -2, -2, 3, 0, -6, 2, -1, 4, 2, -2, 0, -3, -7, -1, -1,
            0, -1, -4, -2, -5, 3, -4, 2, 2, -1, -1, 7, -1, 3, 6, -7, 1, -5, 0, -7, 4, 3, -5, -1, 0,
            3, -4, 1, 2, -7, 1, -2, -8, -2, -5, -5, 1, -4, -4, 4, -3, -2, 2, -4, -8, -1, 0, -9, 5,
            -1, -2, 3, 2, 6, -1, 1, -1, -5, 5, 9, 3, -6, -5, 1, -6, 0, 2, -4, 6, 2, 7, 2, 15, 0,
            -2, 9, 0, 1, 6, 4, -1, -1, -6, -3, 3, 1, -6, -3, 2, 2, -2,
        ];
        let capital_f = vec![
            0, -25, -39, 21, 7, -5, -10, 4, -1, -38, -9, -1, 4, -23, 15, -1, 8, 1, -38, 41, 29, 22,
            9, 12, -46, 0, 9, -17, -19, 32, 38, -3, 14, 6, 2, -6, -18, -1, 23, 80, -12, -20, 24,
            22, -31, -38, -11, 8, 17, 18, 19, -10, 0, -1, 28, -5, -28, -33, 4, -31, -33, -8, -9,
            -44, 46, -11, -5, -21, -22, -7, 1, -11, 33, -8, 12, -7, -6, 63, 17, 12, -49, -11, -31,
            -8, 7, -28, 33, -28, -19, 8, 46, -73, 9, 32, 18, 7, -43, 0, -6, -4, 8, -39, -17, 11,
            15, -25, -9, -28, -2, 24, -23, 10, -15, 4, 41, 46, 18, 2, -3, -29, 11, -3, 20, 35, 21,
            23, 5, -8, -3, -27, -69, 0, 26, -29, -24, 8, 19, 6, -14, -18, 47, 5, 21, -50, 17, -44,
            -36, 24, 9, 16, -38, -5, -54, 34, 13, 31, -2, 9, 8, -12, -14, -17, 28, -59, -20, 19,
            31, 14, 14, 7, -32, 37, 5, -3, -7, -6, 21, -29, -33, 23, -25, -23, 14, 38, -29, -33,
            -9, 23, -43, 18, -12, 2, 30, 32, -28, -21, 42, 1, 6, -6, 58, 34, -22, 1, 5, -2, -8, 14,
            -19, -4, -6, 10, -3, -3, 32, 18, -19, -12, 49, 13, 4, -18, 57, 37, -19, 25, 14, 18,
            -51, 13, 4, 4, 17, -37, -2, 1, 41, -36, -8, -13, 49, -6, 9, 46, -36, -6, -20, -18, -6,
            -29, -42, -21, -25, -29, 5, -41, 51, 49, -20, -22, -9, 3, -6, -52, 10, 41, 12, -27,
            -20, 31, -17, -23, -16, 3, 44, -3, -5, -2, 0, -22, 14, -30, -41, 3, -27, 3, 18, 38, 10,
            49, 45, -13, -27, -4, -10, -67, -1, -17, -2, 72, 46, 20, 24, 22, 16, 25, 6, -6, -31, 2,
            0, -13, -14, 9, 4, 31, 18, 22, 12, 59, -1, -3, -24, -47, -10, 48, 37, -34, -32, -4, 18,
            -2, 52, -8, -7, 34, -44, -14, -21, -49, -35, 41, -4, 31, 3, 23, 9, 8, 0, -24, 38, -9,
            -9, 4, -10, -55, -19, 21, 27, 22, 41, 6, -23, 41, -2, 28, -46, 20, 52, 16, 20, 32, 18,
            2, -3, 9, 16, 33, -18, 12, 6, -9, -19, 1, -5, -15, -17, 6, -3, 4, -22, 30, -34, 43, -4,
            9, -3, -33, -43, -5, -13, -56, 38, 16, 11, -36, 11, -4, -56, 2, 0, -19, -45, -8, -34,
            16, 31, -3, 16, 27, -16, -9, 8, 45, -51, -20, 62, -17, -4, 4, 17, -45, 4, -15, -19, 39,
            39, 15, 17, -19, 2, 45, 36, -22, 16, -23, 28, 34, 12, 5, 10, -7, 28, -35, 17, -37, -50,
            -28, 19, -25, 9, 45, -6, -7, -16, 57, 27, 50, -30, 2, -10, -1, -57, -49, -23, 0, -9,
            -36, -4, -3, 32, -6, -25, 67, -27, -19, 25, -6, 1, -17, -14, 0, 29, 26, -12, -20, 44,
            14, 10, 8, -11, -18, -53, 22, 25, 27, 35, 6, -16, 12, 71, -8,
        ];
        let capital_g = vec![
            27, 6, 12, -3, -31, -42, 27, 17, 11, 8, 34, 6, -3, 2, 11, -11, 18, 48, 1, 21, -7, -6,
            9, 33, -18, -40, -55, -9, -71, -50, 32, -36, 11, 4, 29, 33, 10, -19, -43, -10, 22, -36,
            -23, -21, -14, -47, 25, -4, -14, 30, 16, -18, -11, 6, -37, -27, -12, 6, 7, 33, -36, 33,
            -2, 12, -21, 1, 16, 49, -11, -16, -41, 15, 11, 8, 20, -15, 26, -8, 11, -43, -36, 28, 2,
            -47, -30, -47, -1, 1, 48, -6, -22, 24, -20, -3, -1, -15, -12, 62, 12, 7, -9, 15, -71,
            49, 22, 27, 20, -8, -28, -13, -31, 18, 28, 54, 29, 5, 0, 33, -5, -22, -21, -12, -14,
            -2, 11, -24, 32, -26, -71, 21, -15, -20, -12, 36, -5, 35, 46, 13, -34, -8, 10, -10, 10,
            40, -52, 8, 0, 18, -33, -10, 8, 43, -8, -6, -31, -17, 19, 30, 12, -9, 8, -19, -32, -18,
            -1, -37, 4, 43, 27, 14, -6, -14, -44, -34, -8, 16, -39, 13, 6, -32, 8, 17, -12, 23,
            -44, -25, -66, -12, -31, 30, 14, -9, -5, -10, 44, -12, -2, -43, -22, -18, -7, -9, -15,
            -7, -21, -27, -5, 1, -13, -10, 8, -8, 29, 21, 64, 47, -28, -9, -28, 25, -47, -34, -3,
            -14, -26, -12, -5, -10, -27, -9, -14, -23, -2, -31, 28, 17, -4, -30, 31, 3, -15, 25, 9,
            -32, 0, -6, -22, 20, -37, 3, 12, -19, -17, 13, 30, 11, -15, 15, 50, 66, -31, -31, 16,
            2, 3, -8, 40, -21, -31, -2, 41, -29, -12, 9, 14, -4, 9, 8, -20, 28, 12, 20, -10, 5, -6,
            -33, 6, 21, 51, 30, 9, 3, 8, 7, 19, -53, 19, 15, 4, -38, 19, 29, 18, 6, 19, 3, -17,
            -32, 16, 3, 46, -6, -3, 47, 3, -66, 3, 25, -6, -6, 21, -24, -9, 28, -39, -42, 42, -6,
            -19, -14, 6, -8, 9, 28, -4, 23, 12, -17, -13, 3, 3, 6, 44, 6, -5, 38, -4, -16, 12, -15,
            8, -11, 45, 1, -16, 37, -35, 20, 26, 9, 13, 34, 25, -3, -10, -2, -42, -23, -22, -56,
            -56, 6, 17, -9, 0, 36, 20, 6, -58, 12, 0, -3, -29, -49, -24, -12, -13, 5, -39, -8, 36,
            -9, 44, 35, -64, -22, -12, 26, -15, 41, 36, -19, -37, -20, 46, 35, 9, 32, -5, 27, 21,
            -36, -51, 19, 10, -23, 28, 46, 28, 8, 22, -31, 18, 2, -16, -9, 1, -22, -22, 31, 14, 5,
            44, -3, 38, 0, -12, 50, -23, -19, 1, 42, 15, 1, 13, 32, 45, 37, 15, 11, -9, -23, -6,
            -23, 36, 4, -34, -14, -14, -37, -28, 19, 20, 14, 24, -48, -34, -27, -34, -12, 9, -20,
            -30, 25, 28, -51, -13, 11, -20, -1, -3, 6, -38, -46, -15, 28, 10, -4, 3, -1, 4, -40,
            16, 61, 31, 28, 8, -2, 21, -3, -25, -12, -32, -15, -38, 20, -7, -35, 28, 29, 9, -27,
        ];
        let b0 = [
            Polynomial::new(g),
            Polynomial::new(f.into_iter().map(|i| -i).collect_vec()),
            Polynomial::new(capital_g),
            Polynomial::new(capital_f.into_iter().map(|i| -i).collect_vec()),
        ];
        let sk = SecretKey::<512>::from_b0(b0);

        let expected_signature_vector = vec![
            11, 201, 176, -24, -141, -151, -63, -323, 154, -363, 168, -173, -29, -184, -142, 419,
            -48, 104, 103, -245, -374, 252, -59, 32, 77, -237, 182, -9, 181, -54, -47, 52, -6, 81,
            147, 113, -36, 28, -156, -261, -277, -431, 175, -182, 115, -273, 33, -76, -270, -124,
            -25, -61, -166, 65, -9, 34, 52, -104, 240, -81, 120, 55, 9, 273, -13, -1, -193, 442,
            -43, -58, -86, -100, -14, -96, 245, -120, 10, 2, -40, 341, 8, 112, -260, 100, -24, -22,
            -181, -207, -123, -6, 108, -271, 194, 131, -60, 87, -66, 173, 44, 133, -270, -182, 176,
            59, 289, 25, 98, -47, 153, -257, 160, -21, 73, 58, -4, -39, 79, -124, 31, 119, -175,
            -125, -222, -36, 71, 3, -153, -101, 20, 234, 235, 162, -147, -18, 155, -11, -90, -157,
            -18, -408, -18, -53, -16, 169, 104, -135, 303, -219, 572, 109, -235, -478, 114, 66,
            -17, 186, -13, -57, 31, -132, 73, 134, 35, -165, -279, 27, -360, -3, 44, -40, -262, 60,
            100, 35, 78, -102, -281, -189, -66, 122, -65, -73, -287, -236, -131, -121, -24, 72, 68,
            -156, -69, 54, -127, -185, 154, 60, 144, -99, -81, 139, 80, 98, -93, 227, 170, -338,
            -15, 162, 149, -247, -89, 290, 36, -231, -77, 121, 205, -45, 140, 6, 45, -134, 248,
            -252, 58, 210, 204, 272, 205, 282, 19, -15, 327, 70, 102, -36, 93, 67, -42, -243, 106,
            104, 47, -333, -139, 195, 49, -22, -138, 166, 308, 143, 57, -305, -26, -176, -46, -243,
            -130, 134, -176, -131, -277, 240, -228, -177, 142, -51, 84, 44, 187, 213, 24, 83, -134,
            -202, 286, 48, 58, -199, 7, -18, 173, 113, 52, -190, 1, -117, -177, 122, -229, 83, -90,
            46, 115, 63, -33, -4, 23, -51, 148, 97, 169, -183, -128, 37, 80, 61, 102, -28, 75, 142,
            292, -89, -260, -47, 62, 86, 184, 15, -258, -48, -47, -29, 211, -357, 228, -133, -144,
            275, -110, -127, -83, -74, -89, 149, 9, -44, -208, -46, 121, -157, 147, 216, 133, -96,
            12, 247, 189, 100, -93, 135, -14, 105, 175, -202, 37, 178, 141, 142, -140, -174, -60,
            -13, 95, -208, -84, -52, -144, -125, -2, 63, -436, -273, 47, 106, 122, -221, -180, 104,
            -4, -163, -121, 87, 405, 107, -229, 259, 118, -136, -313, -35, -84, 208, 128, -4, 13,
            304, -40, 75, 165, 183, -196, 7, -48, -21, -250, 160, -280, 370, 91, 198, -228, -70,
            30, -54, -263, -10, -125, -18, -231, -3, 287, -388, -10, 208, -358, -107, 148, -154,
            31, -6, -119, -206, -37, -59, -30, -285, -13, 69, -57, 153, -113, -108, 100, 58, -91,
            -239, -68, -181, 81, 43, 18, -110, -59, -18, 97, -96, 27, 181, -62, -156, -19, -204,
            343, 66, -110, -52, 28, -188, -35, 49, -59, 38, -43, 64, -177, 171, 132, -38, -120,
            214, -42, 110, -324, -34, 158, -102, -4, -61, -117, -134, -310, -99, 79, -308, -306,
            -199, -126, -190, 27, -43, 120, 94, 340, -435, -99, 167, 210, -70, -84, 199,
        ];
        let sig = Signature {
            r: nonce.try_into().unwrap(),
            s: compress(
                &expected_signature_vector,
                (FalconVariant::from_n(512).parameters().sig_bytelen - 41) * 8,
            )
            .unwrap(),
        };

        // We can't recreate this signature because we do not have the seed that generated
        // it.
        // let obtained_signature = SignatureScheme::new(FalconVariant::Falcon512).sign(&data, &sk);

        let pk = PublicKey::from_secret_key(&sk);
        assert!(verify::<512>(&data, &sig, &pk));
    }

    #[test]
    fn test_falcon512_hash_to_point() {
        let nonce = hex::decode(
            "16c12515258093799956368cdfc182c1ca4a34f077e9244416a8c4c13fb0ca241e8b7ac1712d28eb",
        )
        .unwrap();
        let data = hex::decode("6461746131").unwrap();

        let expected_hashed_message = vec![
            977, 11612, 3879, 10128, 2643, 9689, 2895, 3592, 4764, 8492, 968, 407, 11398, 10527,
            923, 5167, 494, 8730, 3485, 6817, 4384, 11423, 11943, 1750, 8829, 10154, 3518, 4588,
            7584, 156, 4405, 9395, 1883, 12126, 12150, 7547, 3120, 8963, 9497, 5096, 5924, 4718,
            1328, 11255, 8140, 1377, 8027, 24, 10527, 1668, 3720, 3820, 5208, 6072, 2256, 741,
            1156, 7665, 1064, 5373, 4650, 10410, 11134, 10688, 10785, 760, 1487, 6035, 716, 10795,
            6615, 1445, 7660, 1637, 12112, 9136, 8753, 8081, 7723, 5783, 11980, 1656, 8283, 12077,
            6793, 1332, 8227, 9045, 1860, 8568, 7432, 9598, 2084, 11042, 3331, 3048, 7274, 7065,
            3761, 11233, 12073, 4560, 3477, 10847, 10512, 4639, 5374, 5082, 5054, 12251, 3088,
            2977, 6803, 2956, 4816, 9634, 11751, 3437, 12106, 8290, 8498, 7231, 6913, 8775, 6571,
            944, 1032, 11603, 4302, 4380, 5407, 10770, 6671, 4000, 11059, 5714, 8030, 9352, 3340,
            6423, 7067, 6530, 5156, 2006, 8675, 10974, 6729, 1761, 9762, 2740, 11483, 1904, 8598,
            3360, 6599, 3538, 4020, 6396, 12226, 11975, 9426, 7804, 343, 7290, 11788, 11834, 11846,
            11894, 6373, 11229, 5011, 11232, 10169, 2464, 6097, 289, 1840, 9180, 10229, 7043,
            10333, 8201, 6892, 8934, 5746, 6782, 3368, 6631, 11854, 8109, 7960, 8229, 2745, 5938,
            4355, 132, 5094, 8153, 3331, 4309, 5245, 2892, 2755, 11201, 5536, 3497, 12073, 4933,
            8447, 11076, 2840, 2497, 2081, 4252, 10941, 10822, 10427, 7986, 6037, 8131, 10966,
            9627, 2332, 11303, 11452, 5431, 3597, 6932, 2712, 2603, 11281, 1012, 5040, 4090, 2542,
            8708, 9038, 457, 6992, 11781, 6356, 4082, 3863, 11526, 7129, 10981, 9562, 4342, 3604,
            3685, 12033, 7627, 8533, 5874, 3955, 4701, 9266, 10778, 2656, 6205, 4977, 588, 7670,
            2065, 2238, 11224, 5351, 1037, 9819, 7915, 7871, 574, 10950, 6989, 2093, 4438, 2690,
            1589, 4797, 2751, 9437, 7643, 11481, 4645, 3034, 11297, 8861, 2951, 10350, 11820, 4896,
            8326, 7779, 11297, 6753, 6525, 8499, 733, 7880, 4909, 11502, 6823, 3170, 11364, 5546,
            11293, 1990, 6678, 7149, 1556, 1401, 2716, 812, 6568, 8520, 418, 982, 12154, 4464,
            7402, 10773, 8443, 4828, 4312, 10068, 11424, 8, 8066, 10838, 5792, 8558, 4009, 1741,
            9716, 3523, 11119, 6253, 2483, 10166, 7762, 5445, 4895, 361, 10352, 11872, 3923, 1005,
            1561, 2553, 3363, 7255, 10705, 6487, 4582, 3882, 618, 97, 12198, 1759, 764, 8882,
            10360, 3752, 3063, 3508, 7641, 4798, 670, 513, 1705, 7189, 8218, 7803, 11232, 3775,
            6056, 11092, 11703, 1881, 11120, 1640, 8222, 11309, 8521, 2462, 5381, 8995, 7441, 9064,
            1727, 7701, 2641, 1050, 9742, 6919, 8989, 10107, 10889, 81, 2397, 1077, 7790, 5486,
            5794, 5217, 7668, 353, 11924, 10005, 11648, 5042, 10776, 1548, 386, 10107, 11119, 322,
            842, 3726, 1678, 4303, 210, 3328, 5753, 9479, 7092, 928, 6163, 7554, 9848, 5259, 3821,
            681, 2527, 10132, 12212, 3163, 9699, 5026, 3727, 1442, 1504, 11759, 9288, 8203, 4091,
            851, 4612, 6287, 10109, 7232, 6913, 7903, 11592, 12135, 5432, 3895, 1597, 11587, 2977,
            3447, 1840, 5445, 11077, 7999, 11472, 10726, 404, 3708, 9221, 9366, 11591, 2898, 1014,
            919, 9524, 7885, 2737, 11699, 8864, 12218, 12243, 4911, 949, 4041, 2898, 6787, 4742,
            3991, 9470, 9737, 968, 7995, 10912, 9080, 9857, 11818, 10201, 8498, 4370, 3341, 2012,
            11164, 11901, 7971, 3049, 10352, 6376, 1011, 1646, 1917, 11359,
        ];
        let r_cat_m = [nonce.to_vec(), data.to_vec()].concat();
        let obtained_hashed_message = hash_to_point(&r_cat_m, 512);
        assert_eq!(
            expected_hashed_message,
            obtained_hashed_message
                .coefficients
                .iter()
                .map(|i| i.value())
                .collect_vec()
        );
    }

    #[test]
    fn test_falcon_1024_test_vector() {
        let nonce = hex::decode(
            "a0b3070ef4bab0ef18c1ebbc79c3014a8611b452d52c14cb5c41087dcc5ee6717e4492cdfe2b8507",
        )
        .unwrap();
        let data = hex::decode("6461746131").unwrap();

        let f = vec![
            3, 5, 1, 2, 0, 4, 2, 2, 2, 3, 2, -4, -2, 6, 2, -4, -3, 3, -1, -3, 3, -3, 3, -4, 2, 3,
            1, 1, 3, -2, -4, -1, 5, 2, 2, -2, 0, 0, -2, 0, -4, 3, 0, -2, -1, 0, -4, 1, -1, -3, 3,
            0, 2, -1, 3, -2, 0, 1, -1, -3, -1, -1, -3, -1, 0, -4, 2, -1, 1, 0, -1, -3, -6, -2, 2,
            -1, -1, -2, 2, 1, -3, 2, 1, -2, 1, -4, -2, -3, -3, -1, -1, 3, -3, -5, -4, -3, 2, 0, -3,
            -4, 2, -3, -1, 3, 0, -1, 1, 0, -3, -3, 3, 3, 2, 5, 3, 0, 1, -2, 0, 5, 6, 5, 0, 5, -1,
            3, 5, -1, 3, -6, -1, 6, 2, -4, 0, 0, -1, -2, -1, -6, -1, 0, 3, -1, -1, 3, -2, -2, -4,
            -1, -2, -3, -8, 4, -4, 5, 1, 1, -1, -1, 2, 0, 3, 7, 0, -2, 1, 0, -2, 0, -3, 0, -7, 4,
            3, 2, 2, 5, -2, -3, -4, 1, -5, 1, 0, -1, -2, -1, -2, 3, 3, 2, -6, -2, 3, -2, -4, 0, -2,
            -2, 1, 2, 2, -3, -1, -1, 1, 5, 0, 3, 4, -1, -6, 6, -3, 0, 0, -1, -3, 0, -1, 2, 0, -1,
            0, -2, 2, 5, 1, -7, -5, 4, 0, 0, -1, 5, 1, 1, 1, 3, -3, 2, 7, -1, 3, -1, 1, 1, 2, -6,
            -1, 2, 2, 3, 1, 0, 1, 1, 8, -1, 0, -6, 2, -1, 1, -2, -4, -6, 2, -3, 1, 0, -1, 2, 6, -4,
            -2, -2, -2, 4, 1, 2, -8, 2, -6, -1, 3, 5, 4, -6, 2, 0, 0, 0, -9, -5, 1, 9, -3, 2, 3, 2,
            2, 3, -1, -1, -2, 3, 1, 0, 3, 2, -2, -2, 1, 2, 2, 1, 0, 1, 3, -2, -2, 0, 2, 0, 2, 1, 0,
            -2, -1, 4, 2, 3, 0, -1, -1, 1, -4, -2, -6, 3, -1, 1, 3, -5, 0, 2, 4, -3, 1, 2, -1, 1,
            1, -4, 1, 3, -6, -2, 5, 0, 4, -2, -1, -4, 2, -1, 0, -1, 2, 2, -2, -2, 1, 0, 1, 0, 0, 5,
            -2, -2, 2, -5, 4, 2, 1, -1, 3, 0, -6, -1, 2, 1, 0, 6, -1, 2, -1, -1, 0, -3, -1, -2, -2,
            2, 8, 0, 0, -2, -1, -1, 0, 3, -1, 3, 2, 1, -1, -2, -3, 0, -2, 4, 0, -2, 3, -1, -1, 1,
            -5, 2, -4, -1, -1, 1, 2, -2, 2, -1, 1, 0, 6, -1, 1, -5, 3, 4, 3, -5, -5, -3, 0, -5, 0,
            -1, 0, 1, 4, 3, 2, -5, 1, -5, -2, 3, 0, 0, 1, -5, -2, -1, -2, 2, -4, -4, 2, -1, 4, 2,
            -2, -3, 0, -3, -1, -5, 1, -2, -1, 3, 2, -2, -3, -4, 4, -1, 1, 2, 3, -4, 1, 0, -5, 0, 0,
            6, 1, -4, 2, -1, -5, 0, 0, -1, 1, 2, 0, -4, -4, -6, 3, 6, 1, -2, 0, -3, 0, 3, 1, 2, -2,
            -1, -1, -1, 2, -1, 3, -4, 0, -3, -4, 0, 3, 2, 2, 3, 0, -1, 5, 1, 1, 2, 4, -8, 3, -4,
            -1, 2, 1, -2, 1, -3, 2, 2, 2, -8, -3, 2, -1, 5, -2, 1, 0, -4, 4, -4, -1, 1, -2, 0, 1,
            -6, -3, 2, 2, -1, 1, 2, -3, -1, -4, -1, -5, 3, -1, -2, -3, -1, -3, 2, 5, 4, 0, -4, -3,
            -3, -4, 0, 5, 1, -6, 0, 1, -7, 1, 0, -1, 0, 5, -1, -2, 3, 0, 3, 1, -4, 0, -1, -2, -3,
            -1, 3, 1, 8, 6, -1, -1, -1, 1, 2, -2, 2, -2, 2, 1, -3, 1, 1, -2, -3, -1, -6, 2, 3, -1,
            1, -1, -1, -1, 2, 0, 0, 0, -4, -7, -4, -1, 2, 0, -3, 8, 6, 4, 1, -4, 2, -3, 3, 3, 7, 0,
            0, -1, -2, 3, 1, 1, 2, -2, 0, 1, -2, 1, 3, 1, 0, 5, -1, -3, 4, 2, 0, 0, 2, 4, -3, 2, 1,
            -1, -1, -1, -2, 1, 1, 3, 3, -5, 3, 0, 2, 1, 0, -3, 1, -5, -1, 0, 2, 0, 1, -7, -1, 5, 6,
            2, 5, 7, 0, 5, -1, -1, -1, 1, -4, -1, -3, -2, -2, 0, 2, 3, 0, -2, 0, -2, 5, -6, 1, 1,
            -4, 3, 3, -1, 1, 6, 2, 1, -2, 1, 2, -3, -3, 7, 0, 1, 1, -2, 0, 2, 0, -4, -4, -4, -4, 5,
            -2, -1, 5, -1, -1, -1, -3, -4, 0, 4, 4, -1, 2, -2, 2, 1, 1, -2, 4, 2, -4, 3, 3, -2, -5,
            3, -4, -2, -1, 0, 0, 0, 4, -4, 7, -2, 5, 5, 5, -3, 0, 0, 1, 3, -1, -1, -2, 0, 0, -4,
            -2, 0, 3, -4, 0, -1, 7, 0, 0, 3, 2, -1, -4, 0, 5, 0, 1, 3, 3, -5, 5, 0, -3, -4, -4, -3,
            -5, -1, -1, 2, 4, 0, -6, 3, 3, -1, 3, -3, -1, 2, 0, -2, -8, 2, -3, 2, 3, 1, 0, -1, -1,
            -6, -4, 2, -5, -1, -5, 2, 2, -1, -1, 5, 6, 3, -3, 0, -1, -2, 2, -4, -1, 0, 1, 2, 1, 2,
            -2, 4, 5, -3, -1, -1, 1, 1, 0, 2, -2, 2, -1, 3, 2, -4, 2, 3, 2, 1, -5, 0, -2, -4, 4,
            -3, 0, 0, 0, -6, 0, -7, -4, -2, 3, -5, -1, -1, -5, 6, -2, 1, -1, 0, 5, 2, 4, -1, -1, 5,
            -4, -2, -6, 6, 1, -2, 0, 5, 5, 4, -1, 1, -2, 2, 3, -4, -1, 1, 3, -6, 0, -3, 0, 3, -5,
            3, -1, -4, 3, -1, 5, -3, 2, -3, 6, 2, -8, 0, 2, 3, -5, 0, 3, -4, 7, 2, -1, -4, 2, 3,
            -8, -3, -2, -5, 1, 1, -3, -2, -2, 0, 0, 2, 2,
        ];
        let g = vec![
            3, -4, -1, -3, 0, -1, 1, 0, -3, 1, 5, 4, -5, 5, 3, 0, -4, 3, -1, -3, -3, 1, 0, 1, -1,
            -5, 1, -2, 2, -1, 3, -1, -3, 0, 4, -2, 2, 2, 0, -2, 2, -1, -3, -1, -2, 3, 1, 0, 2, 7,
            1, 1, 1, 1, -2, -3, -3, 1, 0, -2, -1, -1, -4, 1, 5, -2, 3, 0, -1, -1, 1, -3, 2, 7, 2,
            4, 1, -1, 3, -1, -2, -6, -4, -8, -2, 4, -3, 0, 0, -4, -1, 4, 0, -2, 2, 3, -1, 0, -1, 1,
            0, 1, 1, -2, 3, 2, 3, -3, 1, 7, 3, 0, 3, 1, 1, -2, 4, 0, 4, -2, 0, 2, 1, 0, -1, -3, 1,
            -1, 1, 0, 0, -1, -1, -1, 7, 0, -3, 1, -4, 6, 3, -2, 1, -1, 0, 5, -2, 0, 0, -3, -2, -4,
            1, -1, 0, -2, 5, -3, 4, -5, -2, 1, 1, 4, 0, 0, 1, 0, 3, -2, -4, 3, 1, 1, 0, 3, -4, 1,
            3, -1, -2, 2, -2, 1, 1, -4, 1, 0, 1, 4, 3, 1, 3, -3, -2, 1, 0, -2, -1, -2, 2, 1, -3,
            -3, 2, 4, 5, -3, 0, -2, -2, -3, 0, 2, -3, 0, -9, 0, -5, -3, 2, 0, -2, 0, -2, -2, -1, 0,
            -2, 1, 3, 2, 0, -1, 2, -1, 3, -2, -1, 0, 2, 5, 4, 4, 6, -3, 4, -1, 0, 0, -4, -1, 2, 3,
            1, 1, 0, 0, -4, -3, -1, -1, 1, 0, 2, 0, 4, -5, 1, 5, -5, 1, -4, -2, 2, 4, -2, 1, -3, 3,
            1, 3, 4, 0, 3, -1, -3, 0, 4, 1, 2, 1, 4, 0, -2, -1, 0, 1, 2, 0, 0, -2, 2, 0, 2, -2, -4,
            -2, 2, 0, -1, 1, 2, 2, 1, 5, 0, -1, 6, 1, -5, 5, 2, -1, -1, -3, -1, -3, 3, 0, 0, -1,
            -3, -4, 0, 0, 2, 0, -4, 4, -3, 3, 1, -3, 0, 3, -2, -1, -1, 1, 1, 0, 0, -2, 0, 0, 3, 1,
            2, -2, 3, -1, -3, 1, -1, 3, 0, 1, -7, 0, 5, 4, -1, 0, 4, 3, 1, 0, 1, -4, 4, -2, 2, 2,
            -2, -1, 1, -2, 2, -3, 2, -2, 3, 4, -1, -4, 0, -3, -2, 2, -4, -3, -4, 1, -3, -1, -1, -5,
            0, 0, 1, 2, -5, -4, -2, 2, 0, -2, -2, -4, -2, -3, -4, 4, 4, -5, -1, -2, 1, -6, 2, -2,
            4, 0, -1, -3, 3, 3, 2, -1, 4, 3, 6, 1, 0, 1, 1, 4, -3, 0, 0, 1, 4, 1, 3, -3, 2, -5, 0,
            1, -1, 0, 0, -1, 1, 0, 2, -6, 0, 0, -1, 3, 2, 1, 0, 0, 1, -1, 4, 1, 0, -3, 4, -4, 1, 1,
            1, 1, 2, 0, 7, -2, 0, -5, 5, -1, -3, -2, -4, 5, -1, 1, 1, 0, 3, 0, 1, 3, 6, 3, -2, 1,
            0, -1, -1, 4, -3, -1, 3, 1, -5, 1, 2, -8, 0, -3, -3, -1, -1, 0, -6, -2, 1, -3, 1, -3,
            -4, 1, -1, 4, 5, 0, -3, -3, 3, 1, -2, -2, 0, 1, 3, -2, 1, -4, -2, 2, 0, -4, -4, -1, -1,
            4, 3, 1, 1, -3, 3, 0, 1, 0, -3, 0, -5, -7, 4, 0, 1, 2, -1, -1, 3, 3, 0, -6, 1, 0, 6, 2,
            -5, 3, -5, -4, -4, -1, -1, -3, -3, 0, -2, -2, -1, -3, -2, 2, -1, 0, 5, -3, -2, 2, 2, 2,
            1, -2, 6, 0, -1, 2, 6, 3, -1, 5, -5, -1, -1, -4, 4, -1, -1, -1, 1, -1, 4, 3, 2, -2, -2,
            -2, -5, -2, 0, -2, 2, 2, -3, 2, 0, -4, 3, 3, 0, 1, -2, -6, -2, -1, -2, -1, -2, 5, 3,
            -3, 4, 5, 3, -2, -6, -3, 0, 0, 5, 0, 3, 1, 1, 5, -3, 2, 3, 3, 3, 7, 3, 2, -2, -2, -3,
            2, 4, 0, -4, 0, -2, 0, 0, -2, 1, 3, 4, 1, 2, 1, 4, 1, -3, -5, 3, 1, 6, 4, -1, 0, 3, 3,
            -2, -1, -2, 1, -2, 2, 4, 3, -2, -2, -2, -4, 1, -1, 5, 2, 2, 0, 2, 2, -5, -6, 3, 0, 2,
            -1, 2, -2, 6, 2, -2, 1, -3, 0, 1, 4, 2, -2, 1, 0, 4, -2, 3, 0, 0, 0, 1, 3, -1, -1, -3,
            -2, 2, 2, 3, -2, -2, 0, 3, 5, 4, -4, 1, 0, -1, -1, 2, -3, 8, 0, 1, 2, 2, 3, 2, -4, 2,
            -2, -1, 3, 2, 4, 6, 3, 0, 4, 2, 4, 0, 4, -4, 0, 2, 1, 3, -2, -2, 2, 6, -2, 0, 2, 6, 2,
            -6, -3, 2, -3, 0, 3, -3, -1, 2, -4, -3, -7, 0, 4, 0, 2, -6, -1, 7, 2, 2, 5, -6, -1, 2,
            -1, 1, -4, -1, 2, -2, -1, 3, -9, -2, 1, 3, 3, -1, 2, -3, -1, -2, 3, 3, -5, -1, -3, -2,
            0, -1, -3, 0, 0, 1, 1, 0, -3, -2, 2, 0, -5, 4, 2, 3, 1, 0, -2, -4, -1, -3, 3, -4, 0,
            -1, -6, 2, -2, 4, 1, -1, 3, -3, -4, 5, -2, -1, 2, 0, -3, 3, -2, 3, -4, -1, -4, 4, 0,
            -4, -1, 5, -1, -1, -1, -2, 0, -3, 3, 0, 0, 0, -1, 1, 0, 1, 0, -2, -5, 3, -5, 2, -4, -4,
            -4, -2, -3, -1, 1, 0, 8, -3, -1, -1, -1, -2, -4, -1, 4, 2, 7, -3, 1, 1, 4, 2, -5, -2,
            3, 4, 1, 4, 1, 2, -4, -1, 0, 4, 6, 3, 1, -6, 2, -3, 2, -1, 0, 1, -3, -2, 0, -4, -2, -3,
            2, -1, 2, 1, 1, 0, -1, -1, -4, 1, 4, -2, -3, -3, -1, -5, 0, -2, 2, 0, -2, -2, 0, -2,
            -3, 3, -2, 3, -2, 1, 0, 2, 1,
        ];
        let capital_f = vec![
            -7, 9, 28, 1, 12, 14, -6, -61, 29, -4, 16, -18, -22, 5, 6, -5, 19, 21, -36, -9, 25,
            -47, 50, 40, 23, -71, 25, -58, -51, 2, -26, 12, -5, 16, -46, 13, 17, 36, -53, -45, -34,
            -49, -1, -12, -5, -59, -50, -3, -23, 40, -22, -17, 35, 5, -2, -11, -61, -20, 14, 12,
            -15, -2, -3, 22, 10, 0, 13, 34, 2, 49, -9, -25, -27, -34, -33, 36, 18, 10, -9, 61, -21,
            -10, 41, -27, 47, 52, -17, -47, -30, 69, -1, 6, -8, -10, 17, -10, -13, 37, 0, -23, 27,
            26, 17, -2, 48, -17, -6, 5, 14, -24, 30, 13, 16, 28, 3, 53, -11, -36, -20, 20, -21, 10,
            59, 3, 53, -29, 39, -48, 11, 61, 4, -4, 15, 47, 9, 12, -4, 31, -7, -3, 40, -75, 12, 6,
            10, -34, -15, 5, -11, -46, 31, -17, -28, -18, -15, 18, -10, 12, -4, -47, 0, 39, 7, 13,
            -27, 26, -22, -15, -8, -8, -18, 35, -23, 5, -31, 42, -2, 26, 46, 6, 3, 3, -28, -12, -6,
            -12, -37, 22, -2, -16, -33, 4, -80, 74, 16, 20, -9, 45, -51, 6, -43, 6, 8, -30, 45, 11,
            -14, 44, -9, -16, 16, -25, -42, 10, 32, -8, 0, 44, 26, -49, 19, 17, 5, -11, 36, 1, 23,
            42, -20, -5, 54, -26, 3, 8, 27, 3, 70, 53, -50, 44, -4, -45, -5, -15, -24, -19, 45, 1,
            26, 1, 17, 33, -6, 4, -33, -14, -47, 25, 61, -26, 13, -15, -4, 24, 54, -1, -16, 27, 52,
            12, -2, 0, -25, -26, -1, -15, -47, -14, -16, 34, 37, 25, -62, 9, -6, 7, -15, 14, 13, 6,
            36, 16, -33, 23, 20, 11, 6, -70, -1, -64, 6, 11, -21, 30, -1, 2, -14, 11, 22, -8, 4,
            -27, -1, 25, -25, 40, 26, -30, 27, -13, 19, -10, 15, 21, -8, -6, 12, -56, -3, 15, 31,
            -1, 38, 14, 9, -58, 1, 46, -1, -31, -72, -24, -2, 25, -7, 27, -53, -34, 14, -11, 43,
            -20, 18, -5, 38, 12, -13, -48, 25, 0, 22, -27, -37, -3, 5, 19, -4, -3, 17, 11, 5, -21,
            -48, -18, 15, 0, -6, -4, 15, -1, 8, 49, -2, 59, -33, -22, -24, 2, 13, -31, 48, -31, 8,
            -31, -4, -25, 19, 26, 8, -32, 22, -2, 14, 11, -5, 7, -43, 18, -41, 8, 19, 1, -4, 36,
            -4, 10, -21, -31, 24, -15, 34, 1, -25, 30, -19, 14, -23, -7, 52, -18, 5, -70, 7, -53,
            -20, 22, -8, -25, -27, 33, 9, -57, -7, -56, -24, 17, -28, 47, -16, 2, 61, -28, -86,
            -28, 31, -54, 0, -44, -43, -37, -57, 61, -20, 38, -25, 26, 16, 17, 23, 20, -3, 5, 15,
            -41, 20, 4, -40, -18, -1, 19, 30, 11, 12, -13, 23, 40, 15, 7, -33, -4, -3, 0, 4, -14,
            7, 34, -2, 20, 29, 40, -17, 3, -16, 11, -17, 10, 0, -2, -19, -2, 7, -9, -6, 24, 12, 22,
            51, 17, 23, 8, 14, 16, 3, -50, -12, -12, 2, 19, -17, 41, 19, 70, -36, -15, -3, -6, 33,
            -16, 19, -27, -12, -14, 27, -43, 30, -8, -15, 5, -33, 35, 17, -26, 2, 57, 16, 10, 38,
            5, -5, 2, 3, 12, 21, -4, -29, -13, 2, -4, 16, -25, -23, -13, -39, -15, 38, -36, -23,
            -40, 13, -7, 25, -12, -9, 24, -26, -7, 12, -39, -44, 11, 5, 17, -22, -49, 8, -58, 7,
            17, -23, 0, 19, 32, 54, 11, -5, 19, 23, -19, 52, -36, -9, 38, 9, -12, 44, 33, -3, -29,
            16, 21, 1, -4, -42, 30, -20, 12, 7, 23, 9, 0, -9, 7, 5, 17, -1, -6, -33, -14, -5, 2, 6,
            23, -11, -5, 22, -17, -10, -9, 1, -10, -39, 23, 23, 24, 12, 13, 21, -76, -16, 12, 0,
            12, 31, -6, 50, 9, -11, -38, 8, -25, 37, 16, 1, 31, 47, -15, -2, -9, 9, 9, 54, 18, 11,
            26, 17, 57, -31, 36, -9, 25, -13, 2, 11, -5, 23, -8, -27, -4, 9, -14, -50, -55, 46, 13,
            24, -15, -10, -30, -49, -9, 16, -19, -38, 26, 38, -31, 14, 19, 19, -22, -8, -18, 7,
            -37, 10, -37, -34, -8, 25, 15, -53, 13, 22, 26, 27, -8, -32, 5, -43, 24, -26, 20, -29,
            31, 42, -21, 28, 10, 1, 39, -1, 32, -60, 6, -24, -3, -5, 15, -19, -24, -5, -39, 7, 16,
            17, -40, -6, 6, 30, 14, 10, -12, 1, -1, -26, -6, 38, -81, 45, -20, -14, -22, -7, -14,
            25, -12, 16, -25, 84, -26, 37, 55, -38, -12, 0, 6, -2, 2, 10, 10, -2, 40, 24, -36, -18,
            -26, -12, 63, 40, 14, 15, -40, -8, -6, 11, 5, 11, 39, -44, -18, -15, -30, 10, 2, -50,
            9, -8, 14, -12, -6, -28, 23, 2, -45, -23, -17, -48, -11, 19, -47, -30, 8, -18, -30, 17,
            57, -7, 11, 11, 8, -36, 9, 2, 25, 11, -4, -12, -11, -23, 69, -32, 16, -2, -6, 48, -27,
            2, -42, -11, -10, -37, 29, -4, 8, 41, 15, 4, -23, 39, 9, 32, -8, -1, -30, 16, -1, -16,
            12, -24, 0, -15, -19, -9, -41, -20, -67, 19, -1, 12, 14, 0, 6, 0, 8, 19, -25, 26, -31,
            12, -16, 22, -7, -20, 21, -9, 7, 8, 6, 71, 33, 28, -51, -11, 30, -3, 20, -17, 2, 22,
            12, -39, 7, -24, 48, 18, 34, -10, 29, 23, 52, 11, 46, 29, 54, -55, -13, -30, -25, -24,
            -11, -12, 7, 3, 17, -15, 15, -59, -7, 32, -36, 7, -19, 19, -18, -16, -9, -31, 4, -2, 6,
            13, -69, -14, 15, -1, 44, -11, -15, -19, 24, -3, -34, -29, -20, 17, 49, 17, 42, 50, -9,
            -12, 34, 6, -60, -4, 29, -13, -6, 34, 51, -37, -8, 31, -21, 26, 12, -21, 52, -19, 1,
            -33, -1, 31, -25, 29, 20, -28, 1, -49, -20, -13, -54, 26, -12,
        ];
        let capital_g = vec![
            25, -7, 1, -20, -7, -18, -21, -16, 20, 26, -2, 60, -15, -17, -30, -19, -22, 26, -51, 2,
            -21, 10, 5, 17, 34, -56, -7, -39, -59, -26, -27, -37, -4, -13, -16, 0, -13, -14, -30,
            -36, 10, 1, 5, 36, -6, -26, -38, -9, -18, 8, -17, 11, -32, 4, -24, 24, 4, 28, 9, 82,
            -7, 27, -7, -5, 4, 34, 50, -34, 32, 8, 8, 31, 25, -67, -16, 61, -10, 8, 46, -16, -5,
            56, -25, 44, 2, 42, 19, -1, -21, -9, -26, -41, 5, 16, 20, 9, -18, -27, 5, -48, 18, -11,
            -2, 20, -16, 47, 14, 21, 42, -7, -4, 17, -33, -8, -21, -43, -20, 24, -20, -10, 12, 29,
            41, 46, 70, 7, 9, -28, 30, -8, -9, -3, 52, 3, -4, 45, 2, 2, 4, 3, 36, 6, -22, -12, -13,
            26, -4, -2, -9, 6, 50, -2, 31, 9, -36, -7, 13, 4, -19, 0, -14, 35, -14, 25, -15, -25,
            -44, 21, -9, -1, -18, -2, 22, 28, -1, -32, -40, 0, 7, -12, 8, 12, 2, -17, 9, -56, 26,
            39, 5, 17, 31, -11, 27, -28, 27, 4, -26, 23, 29, -6, 17, -40, -21, 52, -74, 39, 19,
            -16, -30, -28, 30, 58, -4, -5, -18, 14, 30, 7, 19, 31, -7, 9, 17, 25, 26, -57, -53, -7,
            8, -5, -9, -55, -15, 33, 36, 36, 13, 7, -11, -8, 35, 26, -36, -4, -28, 12, 18, 2, 34,
            -35, -27, 20, 26, -12, -5, -38, 23, 38, 29, -10, -42, -10, -16, 17, 11, 8, -22, -35,
            55, -11, -23, -17, -18, -25, 36, -32, 31, 4, -34, -6, 18, 3, 5, 14, -6, 21, -31, 2,
            -17, -39, -18, -17, -42, -16, 22, -58, -13, -8, 6, 29, -17, -36, -5, -16, 18, -1, -42,
            -6, 19, -71, 18, 24, 28, -6, -23, 10, -19, 8, 28, -7, -44, -23, 13, 21, -25, -13, -27,
            -2, -3, -1, 12, -19, 4, 14, -4, -51, -36, -8, 17, -34, -5, 15, 2, 24, -24, -10, 6, -3,
            16, -46, -5, 6, -8, 31, -8, -16, -13, -22, 5, 38, -10, 21, -16, -13, -6, 2, 14, 13,
            -39, 43, 13, -21, 35, -13, -33, 2, 61, 42, 15, -35, 25, -33, 39, -8, -46, -25, 15, -6,
            24, 8, -16, 1, 23, 38, -34, 5, 10, 20, 35, 9, -27, -31, 3, -25, -32, 6, -26, 31, 18, 3,
            16, 8, -11, -3, -11, 0, -10, -4, 17, -17, -12, -26, -18, 5, -16, 12, 12, 58, 12, 34,
            -5, -15, 24, -36, 24, -10, -5, 20, 4, -38, -2, -21, 18, 3, -21, 29, 21, -2, -17, -32,
            12, -4, 4, -15, 29, -12, 35, 1, 31, -29, 16, -17, 6, 20, -15, 41, -7, 29, 6, -26, -43,
            27, 10, 11, -22, -22, 29, -33, 3, -61, -64, 30, 43, -17, 39, -16, 13, 16, 16, 29, 35,
            27, 16, -8, -14, 20, -22, -40, 24, -21, 6, 24, -2, -43, -26, 22, -9, -8, -4, 0, 21, 3,
            34, -77, -80, 10, -38, 13, -50, -48, -25, 25, 47, -3, -9, -38, -35, 2, 17, -26, 6, 38,
            -45, -5, -32, -21, -16, -32, 23, -3, 6, 20, -19, -4, -8, 25, -25, 64, 4, 49, 10, 9, 23,
            10, -10, 15, 33, 7, 25, -17, -15, -17, 41, -36, -5, 6, 30, -19, -53, 37, -8, 24, -45,
            -46, -26, 53, -20, -14, 54, -22, 1, -30, -24, 18, 20, -3, 27, -64, -25, 40, 3, -30, 25,
            3, -12, -5, 19, 50, -11, -9, -23, -17, 19, -3, -5, 1, 37, 19, 5, -14, -9, 26, 24, 1,
            -5, 20, 24, 42, 14, 13, 34, -1, -1, 22, 14, -31, 24, -13, -23, 6, -53, -30, -31, -20,
            -10, 13, 59, -1, -50, 25, -29, 26, 16, -40, 16, -13, 28, 7, -9, -48, 40, 22, -7, -14,
            41, -37, 16, -8, 1, 16, 1, 24, -5, 16, 25, 68, -4, -17, -23, 7, -22, 47, -13, 27, 22,
            27, 40, 21, 64, -42, 14, 13, 35, 24, 27, -26, 17, -29, -18, 6, 9, 54, 8, 3, -27, 29,
            -74, -48, -37, -38, 17, 42, -25, -37, -35, -24, 8, -52, -19, -5, 2, 18, 6, -13, 25,
            -38, -1, 0, -7, -8, 19, 4, 27, 15, -10, 29, 20, -14, 40, -8, 7, 22, 68, -43, 41, 52,
            22, -20, 12, 2, 25, 89, -24, -3, 30, -28, -23, -6, -3, 69, -24, 32, 13, 26, 8, 3, 32,
            12, 34, -6, 15, -25, -10, -33, 4, -11, 12, -66, 2, -40, -4, -10, -7, -23, 38, 9, 35,
            -27, -39, -24, 7, -59, -1, 28, 3, -34, 10, 1, 42, 91, -18, -1, -12, 26, 68, 22, 21,
            -27, -31, -6, 10, -7, 82, -19, 11, 20, 40, -5, 26, -21, 61, -37, 1, 8, -34, 19, -40,
            14, 11, 2, -20, -35, 1, 21, -13, -8, 25, -27, -2, -14, -21, -51, -9, 13, -16, -29, -3,
            5, 7, -53, -27, 1, -53, 42, 40, 61, -23, -20, -14, -1, -25, 4, 42, -10, -8, 6, 27, -23,
            -2, -32, -16, 2, -9, 22, -1, -15, 1, -2, 0, -9, 27, 48, -8, 77, -56, 19, -8, -5, 3,
            -39, -17, -23, 1, -18, 29, -14, -17, 39, 44, 16, -16, -11, -35, -56, 22, -15, 8, -48,
            -47, -46, -5, -44, 5, -11, -30, 2, -8, 6, 11, -25, 32, -19, 0, 4, -2, 9, -2, 49, 38,
            30, -16, 15, -4, -32, -42, -11, -2, 32, 22, 0, -10, 15, -17, -21, 21, 21, 15, 4, -11,
            -12, 53, 39, -32, 19, -18, -8, 8, -12, -23, 22, -2, -65, 10, -8, 6, 0, 26, 2, -8, 6,
            -22, -24, 30, -26, -47, 0, -2, 19, 15, 29, -25, -4, -47, 30, -2, 20, 29, 21, -2, 11,
            56, -25, 45, 37, 44, 58, 5, 8, -20, 24, 23, -39, 43, 0, -67, 7, -5, -1, 6, 46, 14, 32,
            47, -24, 53, -25, 10, 21, -33, 36, 2, -50, 20, -2, -16, 38, 2, 11, -37, 48, 20, -19,
            -14, 16,
        ];
        let b0 = [
            Polynomial::new(g),
            Polynomial::new(f.into_iter().map(|i| -i).collect_vec()),
            Polynomial::new(capital_g),
            Polynomial::new(capital_f.into_iter().map(|i| -i).collect_vec()),
        ];
        let sk = SecretKey::<1024>::from_b0(b0);

        let signature_vector = vec![
            -65, 348, 265, 166, -45, 9, 28, 84, 68, 20, -184, 212, -363, -20, -176, -33, 210, 165,
            228, -47, -68, 225, -173, -222, -235, 47, 206, -105, -391, -71, 99, -175, -161, -72,
            175, 14, 165, -74, -159, 139, -114, 75, 270, -85, -143, 254, 43, 68, 18, -337, 35, 162,
            61, -89, -57, 186, -221, 167, -21, 15, 5, -124, 18, 170, -339, 113, -165, 401, 66, 122,
            67, 93, -43, 107, 5, 24, -3, 232, 309, 108, -164, -65, -193, -143, 246, 317, 67, -80,
            -135, -148, 110, 36, -101, 185, 104, 184, -247, 38, 169, -103, 396, -44, 36, -232,
            -163, 174, 9, 3, 382, 79, 87, 18, -275, 130, -11, 46, 408, 66, -102, 199, 206, 164,
            -77, -27, -167, -246, -36, 102, 179, -169, -103, -44, -160, 292, 53, -57, -128, 24,
            -85, 200, 33, 150, 73, -32, 174, 182, -96, -97, -243, 207, 41, 48, 250, 54, -3, 59, 92,
            90, -153, -273, -130, -200, 61, -143, 197, -120, 125, -102, 282, 104, 87, 12, -196,
            209, 25, 166, -115, -114, 158, 153, 328, 171, -178, -327, 182, 99, -493, -40, -99, 333,
            3, 25, 80, 30, 74, -28, 143, -312, 289, -5, 244, -51, -4, 220, 129, 55, 9, 93, 60, 225,
            158, -156, -406, -133, -401, 162, -15, -6, -310, -294, -70, -28, -91, 200, 219, 126,
            156, 115, 45, -8, 109, 132, 254, -168, 196, 73, -186, -27, -108, -189, 173, -110, 195,
            204, -109, -55, 17, 134, 264, -34, 21, 173, 78, 301, -148, 139, -394, -23, 291, -13,
            40, 193, -41, 395, -83, -77, 5, 153, -282, 18, 34, 43, 52, 104, -234, -93, -64, -100,
            8, -33, -96, -197, -471, 129, 107, -187, -237, -385, 11, 74, 46, -280, 108, -191, 45,
            -112, 261, -273, -203, 233, 353, -56, -6, 6, 30, -120, 205, -158, -168, -166, -151, 56,
            99, 89, -91, 113, 47, 56, 3, 205, 356, -1, 113, -223, -119, 252, -21, -336, 120, 286,
            165, -101, -206, -77, 175, 50, 53, 98, 123, -339, 15, 356, -88, 123, 108, 65, -47, 401,
            -32, 319, 345, -13, -89, 153, 258, 159, 30, 385, 423, -38, 313, 27, -11, -122, -93,
            144, 43, -238, 42, -9, -75, 102, 76, -42, 28, -130, 126, 33, 90, 82, 120, 43, 344,
            -182, 34, -2, 59, 90, 45, -85, 176, 159, 182, -111, 355, 376, 89, -7, -54, 76, 205, 19,
            -117, -57, 74, 14, 221, -96, 222, -86, 111, 33, 66, -166, 267, 60, -170, 92, 196, 221,
            416, -1, -242, -243, -387, 78, 147, -36, 370, -52, 74, -126, -61, 366, -97, 112, 369,
            5, 165, 92, 74, -10, 10, -100, -70, -140, 77, 29, -62, 53, -26, 121, 69, -38, -205,
            186, 98, 81, 71, 89, -313, -335, -195, 56, 22, -193, 42, 273, -67, -3, -236, 59, 225,
            -82, -124, -114, -102, -64, -46, -93, -117, 8, 332, 221, 83, -432, -52, 64, -63, -72,
            74, -214, -135, -64, -316, -336, -212, 240, 28, -164, -60, -80, 96, -13, -11, -151,
            -90, 199, 0, 220, 162, -94, 5, 43, 276, -12, -334, -85, 234, 11, 42, -363, 73, 194, 59,
            253, -100, -124, -92, -22, -260, -198, 84, -231, 134, 106, 184, 72, -2, 7, -11, 147,
            -370, -46, 26, -401, 253, -86, -78, 7, -89, -263, 86, -131, 295, 76, -39, -162, -176,
            -15, 24, 226, -3, 161, 320, -232, 22, -125, -27, 153, 160, -25, -183, 55, 190, -53,
            -68, 142, 4, -14, 165, -30, 104, -143, -48, -72, -316, -27, -111, 98, 195, -63, -17,
            122, 152, 69, 65, -134, 74, 15, 45, -93, -231, 104, 51, 154, -128, 282, 29, -134, 55,
            -107, -161, -181, 140, 45, -271, 74, -67, -148, -5, 209, 62, -50, -180, -1, 32, 53,
            -316, 13, 80, 69, -72, -180, 51, 325, -99, -97, -50, -41, 145, 298, 413, 113, -113,
            102, -56, -47, 147, 152, -106, -95, -97, -297, -29, 192, 269, -71, 298, 73, 25, 194, 0,
            -36, -105, -119, 205, -215, -191, -5, -48, 47, 67, 38, -96, 63, 455, 60, 172, -183,
            -249, 10, 291, -213, -220, -161, -42, -183, 23, -245, -122, 26, -41, -354, 112, -154,
            -39, 56, -199, 61, 232, -49, 30, 50, 175, -26, 100, -3, -301, 292, 102, 198, -122, 123,
            -303, 246, -9, -70, -131, -58, -150, -131, -197, 344, -321, -87, 110, 227, -87, -157,
            30, -58, -2, -162, -107, -61, 33, 65, -429, -177, -160, -341, -200, -165, 261, 100,
            383, 83, -232, -173, 122, 94, -294, -268, -324, -237, 229, -71, 329, 194, 66, 78, 4,
            -255, -140, 8, 133, 328, -50, 325, 72, -92, -16, -84, 101, -23, 381, 45, -159, 202, 7,
            118, 121, -94, -4, 96, 27, 178, -194, -185, 173, 40, 102, -194, 103, 270, 180, 8, 87,
            -99, -19, 130, 58, 138, -63, 357, -218, 112, 157, 23, 58, -231, -49, -94, 333, -34, 5,
            -325, 203, 148, 331, 196, 79, -139, -369, 84, -271, 46, 155, 103, -70, 176, -96, -78,
            72, 52, -116, -64, -34, 151, 291, 246, 4, 148, 106, 126, 148, 244, -293, -160, 112,
            112, 188, 124, 17, -82, 212, -168, -16, 78, -82, -89, -29, 79, 20, -94, 41, -214, -171,
            83, -255, -19, 336, 58, -28, -53, 124, -80, 91, -302, 19, -116, 28, 319, -299, 27, 162,
            -336, -51, -147, 161, -174, 2, -185, -38, -285, 80, -90, 173, 265, 71, 4, 114, -27,
            179, 102, 180, -70, -78, -209, 176, -10, -63, -168, -25, 164, -68, 223, -371, -128, 37,
            233, 273, -259, 192, -132, 131, 144, -52, -5, -60, -71, -318, -78, -221, 82, -262, 66,
            -135, 13, 72, -94, -309, -14, -25, 1, 245, -436, -11, 170, -78, 190, -203, 78, -107,
            57, -185, -136, 100, -125, -105, 300, 82, 110, 65, 16, 11, 443, -138, -146, -67, -76,
            106, -347, -50, -217, -159, 164, -74, -74, 115, -264, 95, -228, 220, 205, 97, -58, 25,
            -149, -101, 267, -299, 90, -185, 247, -82, 62, -100, -63, 182, -178, 159, 66, 121, -21,
            -258, 55, -130, 190, -133, -34, 121, -293, -124, -130, -98, 20, -56, -9, 21, -266, -12,
            -59,
        ];
        let sig = Signature::<1024> {
            r: nonce.try_into().unwrap(),
            s: compress(
                &signature_vector,
                (FalconVariant::Falcon1024.parameters().sig_bytelen - 41) * 8,
            )
            .unwrap(),
        };

        let pk = PublicKey::from_secret_key(&sk);
        assert!(verify::<1024>(&data, &sig, &pk));
    }

    fn signature_vector(n: usize) -> Vec<i16> {
        match n {
            512 => vec![
                11, 201, 176, -24, -141, -151, -63, -323, 154, -363, 168, -173, -29, -184, -142,
                419, -48, 104, 103, -245, -374, 252, -59, 32, 77, -237, 182, -9, 181, -54, -47, 52,
                -6, 81, 147, 113, -36, 28, -156, -261, -277, -431, 175, -182, 115, -273, 33, -76,
                -270, -124, -25, -61, -166, 65, -9, 34, 52, -104, 240, -81, 120, 55, 9, 273, -13,
                -1, -193, 442, -43, -58, -86, -100, -14, -96, 245, -120, 10, 2, -40, 341, 8, 112,
                -260, 100, -24, -22, -181, -207, -123, -6, 108, -271, 194, 131, -60, 87, -66, 173,
                44, 133, -270, -182, 176, 59, 289, 25, 98, -47, 153, -257, 160, -21, 73, 58, -4,
                -39, 79, -124, 31, 119, -175, -125, -222, -36, 71, 3, -153, -101, 20, 234, 235,
                162, -147, -18, 155, -11, -90, -157, -18, -408, -18, -53, -16, 169, 104, -135, 303,
                -219, 572, 109, -235, -478, 114, 66, -17, 186, -13, -57, 31, -132, 73, 134, 35,
                -165, -279, 27, -360, -3, 44, -40, -262, 60, 100, 35, 78, -102, -281, -189, -66,
                122, -65, -73, -287, -236, -131, -121, -24, 72, 68, -156, -69, 54, -127, -185, 154,
                60, 144, -99, -81, 139, 80, 98, -93, 227, 170, -338, -15, 162, 149, -247, -89, 290,
                36, -231, -77, 121, 205, -45, 140, 6, 45, -134, 248, -252, 58, 210, 204, 272, 205,
                282, 19, -15, 327, 70, 102, -36, 93, 67, -42, -243, 106, 104, 47, -333, -139, 195,
                49, -22, -138, 166, 308, 143, 57, -305, -26, -176, -46, -243, -130, 134, -176,
                -131, -277, 240, -228, -177, 142, -51, 84, 44, 187, 213, 24, 83, -134, -202, 286,
                48, 58, -199, 7, -18, 173, 113, 52, -190, 1, -117, -177, 122, -229, 83, -90, 46,
                115, 63, -33, -4, 23, -51, 148, 97, 169, -183, -128, 37, 80, 61, 102, -28, 75, 142,
                292, -89, -260, -47, 62, 86, 184, 15, -258, -48, -47, -29, 211, -357, 228, -133,
                -144, 275, -110, -127, -83, -74, -89, 149, 9, -44, -208, -46, 121, -157, 147, 216,
                133, -96, 12, 247, 189, 100, -93, 135, -14, 105, 175, -202, 37, 178, 141, 142,
                -140, -174, -60, -13, 95, -208, -84, -52, -144, -125, -2, 63, -436, -273, 47, 106,
                122, -221, -180, 104, -4, -163, -121, 87, 405, 107, -229, 259, 118, -136, -313,
                -35, -84, 208, 128, -4, 13, 304, -40, 75, 165, 183, -196, 7, -48, -21, -250, 160,
                -280, 370, 91, 198, -228, -70, 30, -54, -263, -10, -125, -18, -231, -3, 287, -388,
                -10, 208, -358, -107, 148, -154, 31, -6, -119, -206, -37, -59, -30, -285, -13, 69,
                -57, 153, -113, -108, 100, 58, -91, -239, -68, -181, 81, 43, 18, -110, -59, -18,
                97, -96, 27, 181, -62, -156, -19, -204, 343, 66, -110, -52, 28, -188, -35, 49, -59,
                38, -43, 64, -177, 171, 132, -38, -120, 214, -42, 110, -324, -34, 158, -102, -4,
                -61, -117, -134, -310, -99, 79, -308, -306, -199, -126, -190, 27, -43, 120, 94,
                340, -435, -99, 167, 210, -70, -84, 199,
            ],
            1024 => vec![
                -65, 348, 265, 166, -45, 9, 28, 84, 68, 20, -184, 212, -363, -20, -176, -33, 210,
                165, 228, -47, -68, 225, -173, -222, -235, 47, 206, -105, -391, -71, 99, -175,
                -161, -72, 175, 14, 165, -74, -159, 139, -114, 75, 270, -85, -143, 254, 43, 68, 18,
                -337, 35, 162, 61, -89, -57, 186, -221, 167, -21, 15, 5, -124, 18, 170, -339, 113,
                -165, 401, 66, 122, 67, 93, -43, 107, 5, 24, -3, 232, 309, 108, -164, -65, -193,
                -143, 246, 317, 67, -80, -135, -148, 110, 36, -101, 185, 104, 184, -247, 38, 169,
                -103, 396, -44, 36, -232, -163, 174, 9, 3, 382, 79, 87, 18, -275, 130, -11, 46,
                408, 66, -102, 199, 206, 164, -77, -27, -167, -246, -36, 102, 179, -169, -103, -44,
                -160, 292, 53, -57, -128, 24, -85, 200, 33, 150, 73, -32, 174, 182, -96, -97, -243,
                207, 41, 48, 250, 54, -3, 59, 92, 90, -153, -273, -130, -200, 61, -143, 197, -120,
                125, -102, 282, 104, 87, 12, -196, 209, 25, 166, -115, -114, 158, 153, 328, 171,
                -178, -327, 182, 99, -493, -40, -99, 333, 3, 25, 80, 30, 74, -28, 143, -312, 289,
                -5, 244, -51, -4, 220, 129, 55, 9, 93, 60, 225, 158, -156, -406, -133, -401, 162,
                -15, -6, -310, -294, -70, -28, -91, 200, 219, 126, 156, 115, 45, -8, 109, 132, 254,
                -168, 196, 73, -186, -27, -108, -189, 173, -110, 195, 204, -109, -55, 17, 134, 264,
                -34, 21, 173, 78, 301, -148, 139, -394, -23, 291, -13, 40, 193, -41, 395, -83, -77,
                5, 153, -282, 18, 34, 43, 52, 104, -234, -93, -64, -100, 8, -33, -96, -197, -471,
                129, 107, -187, -237, -385, 11, 74, 46, -280, 108, -191, 45, -112, 261, -273, -203,
                233, 353, -56, -6, 6, 30, -120, 205, -158, -168, -166, -151, 56, 99, 89, -91, 113,
                47, 56, 3, 205, 356, -1, 113, -223, -119, 252, -21, -336, 120, 286, 165, -101,
                -206, -77, 175, 50, 53, 98, 123, -339, 15, 356, -88, 123, 108, 65, -47, 401, -32,
                319, 345, -13, -89, 153, 258, 159, 30, 385, 423, -38, 313, 27, -11, -122, -93, 144,
                43, -238, 42, -9, -75, 102, 76, -42, 28, -130, 126, 33, 90, 82, 120, 43, 344, -182,
                34, -2, 59, 90, 45, -85, 176, 159, 182, -111, 355, 376, 89, -7, -54, 76, 205, 19,
                -117, -57, 74, 14, 221, -96, 222, -86, 111, 33, 66, -166, 267, 60, -170, 92, 196,
                221, 416, -1, -242, -243, -387, 78, 147, -36, 370, -52, 74, -126, -61, 366, -97,
                112, 369, 5, 165, 92, 74, -10, 10, -100, -70, -140, 77, 29, -62, 53, -26, 121, 69,
                -38, -205, 186, 98, 81, 71, 89, -313, -335, -195, 56, 22, -193, 42, 273, -67, -3,
                -236, 59, 225, -82, -124, -114, -102, -64, -46, -93, -117, 8, 332, 221, 83, -432,
                -52, 64, -63, -72, 74, -214, -135, -64, -316, -336, -212, 240, 28, -164, -60, -80,
                96, -13, -11, -151, -90, 199, 0, 220, 162, -94, 5, 43, 276, -12, -334, -85, 234,
                11, 42, -363, 73, 194, 59, 253, -100, -124, -92, -22, -260, -198, 84, -231, 134,
                106, 184, 72, -2, 7, -11, 147, -370, -46, 26, -401, 253, -86, -78, 7, -89, -263,
                86, -131, 295, 76, -39, -162, -176, -15, 24, 226, -3, 161, 320, -232, 22, -125,
                -27, 153, 160, -25, -183, 55, 190, -53, -68, 142, 4, -14, 165, -30, 104, -143, -48,
                -72, -316, -27, -111, 98, 195, -63, -17, 122, 152, 69, 65, -134, 74, 15, 45, -93,
                -231, 104, 51, 154, -128, 282, 29, -134, 55, -107, -161, -181, 140, 45, -271, 74,
                -67, -148, -5, 209, 62, -50, -180, -1, 32, 53, -316, 13, 80, 69, -72, -180, 51,
                325, -99, -97, -50, -41, 145, 298, 413, 113, -113, 102, -56, -47, 147, 152, -106,
                -95, -97, -297, -29, 192, 269, -71, 298, 73, 25, 194, 0, -36, -105, -119, 205,
                -215, -191, -5, -48, 47, 67, 38, -96, 63, 455, 60, 172, -183, -249, 10, 291, -213,
                -220, -161, -42, -183, 23, -245, -122, 26, -41, -354, 112, -154, -39, 56, -199, 61,
                232, -49, 30, 50, 175, -26, 100, -3, -301, 292, 102, 198, -122, 123, -303, 246, -9,
                -70, -131, -58, -150, -131, -197, 344, -321, -87, 110, 227, -87, -157, 30, -58, -2,
                -162, -107, -61, 33, 65, -429, -177, -160, -341, -200, -165, 261, 100, 383, 83,
                -232, -173, 122, 94, -294, -268, -324, -237, 229, -71, 329, 194, 66, 78, 4, -255,
                -140, 8, 133, 328, -50, 325, 72, -92, -16, -84, 101, -23, 381, 45, -159, 202, 7,
                118, 121, -94, -4, 96, 27, 178, -194, -185, 173, 40, 102, -194, 103, 270, 180, 8,
                87, -99, -19, 130, 58, 138, -63, 357, -218, 112, 157, 23, 58, -231, -49, -94, 333,
                -34, 5, -325, 203, 148, 331, 196, 79, -139, -369, 84, -271, 46, 155, 103, -70, 176,
                -96, -78, 72, 52, -116, -64, -34, 151, 291, 246, 4, 148, 106, 126, 148, 244, -293,
                -160, 112, 112, 188, 124, 17, -82, 212, -168, -16, 78, -82, -89, -29, 79, 20, -94,
                41, -214, -171, 83, -255, -19, 336, 58, -28, -53, 124, -80, 91, -302, 19, -116, 28,
                319, -299, 27, 162, -336, -51, -147, 161, -174, 2, -185, -38, -285, 80, -90, 173,
                265, 71, 4, 114, -27, 179, 102, 180, -70, -78, -209, 176, -10, -63, -168, -25, 164,
                -68, 223, -371, -128, 37, 233, 273, -259, 192, -132, 131, 144, -52, -5, -60, -71,
                -318, -78, -221, 82, -262, 66, -135, 13, 72, -94, -309, -14, -25, 1, 245, -436,
                -11, 170, -78, 190, -203, 78, -107, 57, -185, -136, 100, -125, -105, 300, 82, 110,
                65, 16, 11, 443, -138, -146, -67, -76, 106, -347, -50, -217, -159, 164, -74, -74,
                115, -264, 95, -228, 220, 205, 97, -58, 25, -149, -101, 267, -299, 90, -185, 247,
                -82, 62, -100, -63, 182, -178, 159, 66, 121, -21, -258, 55, -130, 190, -133, -34,
                121, -293, -124, -130, -98, 20, -56, -9, 21, -266, -12, -59,
            ],
            _ => vec![],
        }
    }

    #[test]
    fn test_signature_deserialize() {
        let n = 1024;
        let sigvec = signature_vector(n);
        let nonce = [0u8; 40];
        let original_signature = Signature::<1024> {
            r: nonce,
            s: compress(
                &sigvec,
                (FalconVariant::Falcon1024.parameters().sig_bytelen - 41) * 8,
            )
            .unwrap(),
        };

        let serialized = original_signature.to_bytes();
        let deserialized = Signature::from_bytes(&serialized).unwrap();

        assert_eq!(original_signature, deserialized);

        let n = 512;
        let sigvec = signature_vector(n);
        let nonce = [0u8; 40];
        let original_signature = Signature::<512> {
            r: nonce,
            s: compress(
                &sigvec,
                (FalconVariant::Falcon512.parameters().sig_bytelen - 41) * 8,
            )
            .unwrap(),
        };

        let serialized = original_signature.to_bytes();
        let deserialized = Signature::from_bytes(&serialized).unwrap();
        let reserialized = deserialized.to_bytes();

        assert_eq!(original_signature, deserialized);
        assert_eq!(serialized, reserialized);
    }

    #[test]
    fn test_signature_deserialize_fail() {
        let n = 512;
        let sigvec = signature_vector(n);
        let nonce = [0u8; 40];
        let original_signature = Signature::<512> {
            r: nonce,
            s: compress(
                &sigvec,
                (FalconVariant::Falcon512.parameters().sig_bytelen - 41) * 8,
            )
            .unwrap(),
        };
        let mut serialized = original_signature.to_bytes();

        // try every byte of header
        for i in 0..8 {
            serialized[0] ^= 1 << i;
            assert!(Signature::<512>::from_bytes(&serialized).is_err());
            serialized[0] ^= 1 << i;
        }

        // try change length
        let longer = [serialized.clone(), vec![0u8]].concat();
        assert!(Signature::<512>::from_bytes(&longer).is_err());

        let len = serialized.len();
        let shorter = serialized[0..(len - 1)].to_vec();
        assert!(Signature::<512>::from_bytes(&shorter).is_err());
    }

    #[test]
    fn test_secret_key_serialization() {
        let sk = SecretKey::<512>::generate();
        let serialized = sk.to_bytes();
        let deserialized = SecretKey::from_bytes(&serialized).unwrap();
        let reserialized = deserialized.to_bytes();

        assert_eq!(sk, deserialized);
        assert_eq!(serialized, reserialized);
    }

    #[test]
    fn test_secret_key_serialization_fail() {
        let sk = SecretKey::<512>::generate();
        let mut serialized = sk.to_bytes();
        let len = serialized.len();

        // test every bit in header
        for i in 0..8 {
            serialized[0] ^= 1 << i;
            assert!(SecretKey::<512>::from_bytes(&serialized).is_err());
            serialized[0] ^= 1 << i;
        }

        // change length
        let longer = [serialized, vec![0]].concat();
        assert!(SecretKey::<512>::from_bytes(&longer).is_err());

        let shorter = &longer[0..len - 1];
        assert!(SecretKey::<512>::from_bytes(shorter).is_err());
    }

    #[test]
    fn test_public_key_serialization() {
        let pk = PublicKey::from_secret_key(&SecretKey::<512>::generate());
        let serialized = pk.to_bytes();
        let deserialized = PublicKey::from_bytes(&serialized).unwrap();
        let reserialized = deserialized.to_bytes();

        assert_eq!(pk, deserialized);
        assert_eq!(serialized, reserialized);
    }

    #[test]
    fn test_public_key_serialization_fail() {
        let pk = PublicKey::from_secret_key(&SecretKey::<512>::generate());
        let mut serialized = pk.to_bytes();
        let len = serialized.len();

        // test every bit in header
        for i in 0..8 {
            serialized[0] ^= 1 << i;
            assert!(SecretKey::<512>::from_bytes(&serialized).is_err());
            serialized[0] ^= 1 << i;
        }

        // change length
        let longer = [serialized, vec![0]].concat();
        assert!(SecretKey::<512>::from_bytes(&longer).is_err());

        let shorter = &longer[0..len - 1];
        assert!(SecretKey::<512>::from_bytes(shorter).is_err());
    }

    #[test]
    fn test_secret_key_field_element_serialization() {
        let mut rng = thread_rng();
        for polynomial_index in [0, 1, 2] {
            let width = SecretKey::<512>::field_element_width(512, polynomial_index);
            for _ in 0..100000 {
                let int = rng.gen_range(-(1 << (width - 1))..(1 << (width - 1)));
                if int == -(1i16 << (width - 1)) {
                    continue;
                }
                let felt = Felt::new(int);
                let ser = SecretKey::<512>::serialize_field_element(width, felt);
                let des = SecretKey::<512>::deserialize_field_element(&ser).unwrap();
                let res = SecretKey::<512>::serialize_field_element(width, des);

                assert_eq!(felt, des);
                assert_eq!(ser, res);
            }
        }
        for polynomial_index in [0, 1, 2] {
            let width = SecretKey::<1024>::field_element_width(1024, polynomial_index);
            for _ in 0..100000 {
                let int = rng.gen_range(-(1 << (width - 1))..(1 << (width - 1)));
                if int == -(1i16 << (width - 1)) {
                    continue;
                }
                let felt = Felt::new(int);
                let ser = SecretKey::<1024>::serialize_field_element(width, felt);
                let des = SecretKey::<1024>::deserialize_field_element(&ser).unwrap();
                let res = SecretKey::<1024>::serialize_field_element(width, des);

                assert_eq!(felt, des);
                assert_eq!(ser, res);
            }
        }
    }

    #[test]
    fn test_falcon1024_hash_to_point() {
        let nonce = hex::decode(
            "a0b3070ef4bab0ef18c1ebbc79c3014a8611b452d52c14cb5c41087dcc5ee6717e4492cdfe2b8507",
        )
        .unwrap();
        let data = hex::decode("6461746131").unwrap();
        let expected_hashed_message = vec![
            12042, 5076, 11082, 8475, 3594, 7967, 1459, 3308, 9106, 2961, 5661, 4394, 6913, 11078,
            9700, 10288, 10693, 3926, 9889, 8313, 9967, 8728, 3097, 8136, 3305, 337, 10006, 12170,
            10438, 11329, 1483, 5761, 5304, 7680, 11509, 11139, 5785, 6618, 10233, 7542, 9200,
            10505, 8737, 12042, 9584, 3965, 9792, 1817, 165, 9292, 6367, 7361, 11940, 6321, 11842,
            9970, 7587, 9759, 9726, 5456, 3229, 8437, 2163, 6712, 9134, 1953, 1439, 6899, 10561,
            11108, 7594, 1191, 8813, 2952, 628, 4942, 7652, 6377, 4763, 6059, 5586, 10931, 7232,
            7963, 8218, 10987, 5760, 2256, 5456, 151, 6569, 7624, 9475, 2943, 2494, 2910, 743, 163,
            5642, 6907, 10323, 10001, 4936, 9805, 2455, 2096, 6324, 4329, 7266, 403, 109, 485,
            6900, 6542, 1932, 616, 11953, 7380, 11304, 7771, 2951, 7506, 2414, 6190, 9796, 4803,
            5528, 2771, 8655, 7443, 7649, 4919, 6315, 12109, 3147, 9968, 5586, 4517, 5179, 2889,
            10052, 9251, 4719, 314, 7170, 4086, 2186, 7566, 5064, 5098, 4547, 7946, 12031, 2812,
            2528, 6580, 7313, 3746, 8823, 273, 4160, 10534, 2407, 2446, 7660, 2075, 10552, 6032,
            10197, 5842, 8159, 11483, 4516, 5693, 4843, 6230, 6228, 203, 1535, 4201, 11654, 8810,
            3432, 1937, 7857, 3853, 3953, 9340, 4289, 3376, 11020, 11257, 8017, 3995, 1247, 9647,
            8109, 3834, 5544, 9000, 8051, 3402, 5173, 3980, 652, 5003, 5861, 1018, 10718, 5525,
            2866, 9027, 6847, 6914, 8004, 6964, 1371, 2226, 2922, 6057, 4646, 3632, 8078, 10268,
            10048, 8842, 12083, 11771, 1286, 1704, 1237, 12286, 5369, 1450, 4608, 7428, 7035, 8368,
            7496, 4735, 6586, 11208, 7320, 1174, 564, 2825, 2232, 302, 5156, 11432, 8660, 2572,
            3791, 10895, 7516, 11939, 5794, 7814, 6927, 6504, 6086, 7562, 7971, 10107, 10457,
            11964, 5883, 3284, 1116, 11177, 8762, 6079, 3895, 8161, 2702, 7060, 4272, 2700, 4248,
            6134, 10226, 12209, 1269, 2126, 529, 3078, 9552, 3957, 2120, 10151, 2884, 2736, 8936,
            2312, 3322, 10434, 1234, 2160, 4277, 768, 10500, 9874, 8286, 1523, 9444, 2688, 5335,
            4868, 6978, 10610, 10616, 3203, 12206, 12212, 3301, 5133, 6444, 8391, 3501, 4433, 3729,
            10363, 6330, 4224, 915, 5011, 9818, 4083, 3989, 2508, 7207, 11530, 4513, 5942, 4891,
            10458, 3901, 9573, 8618, 9230, 9853, 11954, 5680, 9698, 11651, 4759, 3788, 5387, 1101,
            10188, 12202, 7374, 7066, 1777, 6089, 4765, 7611, 731, 5120, 6256, 5506, 9227, 2895,
            7445, 10474, 10470, 3873, 4321, 6060, 6358, 7806, 10506, 6414, 12023, 915, 10709, 6873,
            3832, 11555, 3734, 4527, 5165, 11674, 1328, 988, 8202, 9751, 5616, 5713, 5605, 2773,
            8131, 841, 10294, 3730, 7608, 515, 3689, 2166, 213, 6639, 10963, 4787, 9330, 11288,
            11108, 8784, 1279, 793, 2341, 4702, 12094, 4077, 1924, 2189, 8874, 2528, 3790, 4985,
            9071, 1481, 5837, 11417, 8028, 9227, 5277, 1560, 12113, 7583, 2273, 7697, 4328, 5885,
            11860, 3628, 11975, 10076, 5184, 5337, 5784, 2309, 11642, 8660, 7889, 4098, 2307,
            10728, 1680, 7287, 5601, 11645, 4185, 2749, 4049, 2775, 7546, 6204, 2954, 589, 2279,
            490, 11384, 1761, 605, 8650, 1694, 573, 11168, 7461, 5205, 12114, 2, 3080, 8261, 11312,
            9389, 233, 9778, 9878, 8795, 526, 5805, 3226, 1689, 5878, 110, 7408, 2095, 758, 6010,
            3550, 1857, 4012, 1699, 12211, 1418, 12087, 5659, 623, 1877, 4697, 1584, 7817, 8646,
            4542, 3610, 4898, 4938, 4505, 3836, 7850, 208, 3900, 6176, 11947, 1156, 5150, 6369,
            2610, 886, 9620, 3668, 10611, 4798, 11361, 4423, 2564, 5400, 10993, 7458, 6242, 5031,
            2568, 1279, 10373, 279, 4778, 1904, 4931, 61, 4807, 2625, 2368, 11733, 3252, 9425,
            11160, 1727, 7400, 2859, 1193, 279, 6601, 9726, 8171, 2546, 2079, 4441, 7369, 5458,
            7415, 5616, 3429, 4590, 4055, 1119, 830, 4204, 11149, 6252, 2869, 2182, 5549, 5964,
            9341, 3851, 6609, 1849, 8078, 7959, 4594, 10012, 6102, 4732, 12131, 10072, 10627,
            10673, 10236, 3320, 6131, 3679, 8863, 959, 6178, 2506, 5240, 205, 7666, 12220, 3444,
            657, 2309, 3873, 2589, 10726, 9827, 12093, 3127, 370, 768, 7370, 6384, 11705, 7531,
            8695, 7415, 12280, 9054, 10585, 4879, 2250, 11232, 11133, 7672, 8720, 12066, 11147,
            3448, 7096, 6097, 11820, 743, 6487, 10545, 7861, 8088, 7468, 554, 575, 10024, 8289,
            10596, 4577, 3822, 9952, 3707, 6706, 5550, 6464, 4637, 7717, 4981, 3679, 8900, 7132,
            5091, 581, 8776, 2646, 11821, 7177, 2806, 11060, 9809, 9063, 7579, 4329, 4084, 4934,
            8258, 9963, 12016, 7161, 9105, 6067, 7586, 8607, 10612, 4832, 8241, 10787, 4792, 10490,
            7021, 2327, 10768, 6852, 5759, 7865, 7177, 9124, 2886, 10540, 11967, 3598, 7405, 3706,
            12190, 3848, 5200, 3320, 5299, 617, 5473, 8672, 999, 3787, 9074, 6648, 1399, 9297,
            3071, 5090, 11027, 5282, 4484, 3390, 6580, 5377, 9291, 9371, 5033, 4405, 5086, 1070,
            1430, 7299, 5766, 9464, 9440, 9429, 3459, 8929, 5626, 9763, 692, 10653, 10024, 10787,
            8526, 4629, 8932, 503, 6422, 3272, 9605, 3365, 9151, 829, 732, 4840, 9450, 7131, 3521,
            5093, 5156, 10509, 6342, 8922, 5333, 5084, 9460, 1534, 8823, 3706, 1826, 9785, 7836,
            4784, 4167, 1764, 10536, 11561, 9957, 5348, 11868, 6921, 11652, 11133, 9700, 9272,
            2498, 5606, 10459, 3738, 6157, 7448, 8902, 4434, 4118, 834, 3797, 3586, 7668, 1515,
            7386, 4828, 6483, 11041, 7334, 7454, 705, 168, 4871, 10551, 8622, 10939, 11270, 7069,
            5960, 4824, 2689, 329, 2641, 12100, 8266, 11064, 11568, 11423, 135, 255, 11381, 10108,
            5798, 11123, 11000, 5300, 11283, 4991, 9041, 5313, 5963, 3156, 5638, 2124, 6087, 6655,
            9695, 4987, 3814, 9654, 4969, 10874, 6534, 3705, 12129, 4702, 11588, 2956, 1882, 2773,
            3761, 11458, 2258, 9549, 8678, 693, 2762, 8009, 917, 6170, 7562, 11706, 4929, 1480,
            701, 4648, 6941, 3141, 6031, 618, 7268, 970, 1552, 11065, 2914, 2784, 3364, 3538, 2218,
            1269, 3529, 11971, 3479, 5599, 86, 7026, 10803, 8841, 6293, 11862, 6620, 3545, 2713,
            10159, 9929, 5709, 8793, 1736, 10394, 3301, 1018, 4265, 5440, 11899, 10479, 3225, 6105,
            1895, 573, 7100, 9852, 4823, 7027, 362, 11344, 10640, 7808, 704, 3004, 8235, 10740,
            11496, 9564, 2141, 9886, 3754, 11351, 3749, 3393, 8290, 7064, 538, 1041, 4861, 7990,
            184, 2755, 514, 4967, 11789, 8526, 5013, 8862, 588, 10128, 8448, 10872, 11226, 9916,
            7299, 5737, 7815, 597, 10154, 2355, 8051, 6831, 10056, 11402, 787, 8145, 11003, 11017,
            2964, 3279, 11577, 8816, 12122, 7004, 2260, 8300, 436, 888, 8473, 9315, 3909, 2491,
            4001, 2426, 7211, 8185, 9733, 7803, 11572, 4917, 8084, 7833, 6142, 11806, 4804, 9493,
            7805, 8672, 12026, 4720, 4376, 4891, 11288, 3256, 3684, 790, 3484, 8302, 6340, 10754,
            1643, 1381, 4765, 8063, 7944, 4823, 8810, 4037, 10826, 7747, 9843, 7991, 11933, 11183,
            7149, 4495, 3763,
        ];
        let r_cat_m = [nonce.to_vec(), data.to_vec()].concat();
        let obtained_hashed_message = hash_to_point(&r_cat_m, 1024);
        assert_eq!(
            expected_hashed_message,
            obtained_hashed_message
                .coefficients
                .iter()
                .map(|i| i.value())
                .collect_vec()
        );
    }

    #[test]
    fn test_ntru_gen() {
        let n = 512;
        let mut rng = thread_rng();
        let (f, g, capital_f, capital_g) = ntru_gen(n, rng.gen());

        let f_times_capital_g = (f * capital_g).reduce_by_cyclotomic(n);
        let g_times_capital_f = (g * capital_f).reduce_by_cyclotomic(n);
        let difference = f_times_capital_g - g_times_capital_f;
        assert_eq!(Polynomial::constant(12289), difference);
    }

    #[test]
    fn test_gen_poly() {
        let mut rng = thread_rng();
        let n = 1024;
        let mut sum_norms = 0.0;
        let num_iterations = 100;
        for _ in 0..num_iterations {
            let f = gen_poly(n, &mut rng);
            sum_norms += f.l2_norm();
        }
        let average = sum_norms / (num_iterations as f64);
        assert!(90.0 < average);
        assert!(average < 94.0);
    }

    #[test]
    fn test_gs_norm() {
        let n = 512;
        let f = (0..n).map(|i| i % 5).collect_vec();
        let g = (0..n).map(|i| (i % 7) - 4).collect_vec();
        let norm = gram_schmidt_norm(&Polynomial::new(f), &Polynomial::new(g));
        let difference = (norm * norm) - 5992556.183229722;
        assert!(
            difference * difference < 0.00001,
            "norm was {} with square {} =/= 5992556.183229722",
            norm,
            norm * norm
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

    #[test]
    fn test_stalling_babai_reduce() {
        let f : Polynomial<BigInt> = Polynomial::new(vec![
        BigInt::from_str("608952554397054557983876280254475782156978778966514236337232157903009112871321999683643497269694189442905906970303786045645797783553603819271234233389723961161399883913029875233697100089950361781984561312574797799746866190588557496146920").unwrap(),
        BigInt::from_str("-321198091895342626623905018079452846563740549837372660468562894512341046154970222611540742593012146136893774207013217934448033084417751742329986039507508103824063527509404104910887557717069726121009664272991060311076488206587527713367969").unwrap(),
        BigInt::from_str("-511852458099561624239961196282485646021269605948503734916463216271741105669020643689485962072068026680773633158378214687137638656952667582984049436841521038937832328482550655620883633594646529561294167120997868894250092097125289410469269").unwrap(),
        BigInt::from_str("753499692151616415251662722445327731480539623323000605380895648603405371851338467508575040732846737683390835245839956333926493086479631700998116314658893765004936549965922012760065249547741180589285017455538048794978342726037447740353274").unwrap(),
        BigInt::from_str("-183951725423062894154142247006228629686028553344272108728906669710737867431299866959261567846158043338962173741060227623640390657223232373297828195388077477002393007770978315191221057900275286132154955867517828152370806350497427231775784").unwrap(),
        BigInt::from_str("-601002944822358753730908286404268611783579577362298145349220422219337518495823156994812530727209827139319904644420562078413537022624501170363278280532068568605266646229790495883366738116427379176196534468686595436737518807388742965393276").unwrap(),
        BigInt::from_str("624238805077664493322115890964268570726973934215001879901044583238147554868339541992616428233500727841616571552309882882912193503023373275313004172069274706976281040733564722134422543928965613410098136046852457807013063428923886900325416").unwrap(),
        BigInt::from_str("99247320941428668744094858367782360465579561977614882360836919110717822389271909915623399678999810331939700577612687349369204993498887095160873936787628703559789561872359568953256723455070178252408841389920850163746998930970690582316813").unwrap(),
         ] );
        let g : Polynomial<BigInt> = Polynomial::new(vec![BigInt::from_str("-276748395986364908221094085589576475403380860729035685618554719249731859131042127274390301424851882058378070539404235907974793069695572633477961322056365258258614312832511203946728028995719601790538904487894924461813173730797217250818249982").unwrap(),
         BigInt::from_str("86812596127498460519207154109787274799874383708114200980528940455100071427257285893913175247659497014884186080805798356767162662189621835520839194646265741830798779072646145965638756187598754792152809555119081509898641068179709957911369677").unwrap(),
         BigInt::from_str("116354193347135036416083176315819941089583030678185384093353987930114617677847583589518228207710036313227903774513592022562863039840725428131307408866838133762498621270855327564276484292086769092056758574610298629140002234098489707575021867").unwrap(),
         BigInt::from_str("-301807527485830260261344620490924956518725241514546231783820178240968597102415190293928299169893150661326842964193337260684088468944872626419776842954670641099171003070476911114184625914895635366259003019569137386041706705684810645763462742").unwrap(),
         BigInt::from_str("441299161379624032332866662487752791916685084733535111485266772869967918640142735209019134432254486721439741343615198249730448067610849591335249866703346142573143583974028002224021261123474948313561150275893970137607100222815689139896755492").unwrap(),
         BigInt::from_str("-513595683023472695273563978709461994092425319855638422141058717955228319050940927987007122010294647136937120546974596907309751767539707235746382781851494092323493370150359286202256873867086472582425432902385470952671810828641054498438294804").unwrap(),
         BigInt::from_str("507707498023164261062090937388345987181212863525863034372974903786763499191338982755471619717822792036379197240942561142666158588559211642577155844311933842930176434464164432666289434811839369974005888248331254542016292048307196783213633199").unwrap(),
         BigInt::from_str("-424541034576343804926544811528101218199324642128826281424964488598284815916882148212061156604886226511172215884154920113476057900609708969461302301719829549643270043850114557922078292518535617173173882627190400853927778419866400469640897429").unwrap(),
         ]);
        let mut capital_f : Polynomial<BigInt> = Polynomial::new(vec![BigInt::from_str("-97342744417091948277714584678980786347602848582001232696524345330053898461890905857357345522647395659791576336845694789964269088697943709632494103896150736987754856526288094356135043891715189900121635577150577171090069901904263254875521245315624957844848330067137803380242041663554024842364726847685705570657587737754069873830392344722883838185861856446427166745578895270214744011866852187668634399033333451988880963585072356564797880689406987307640156667060527684509114385593018612483494846977442876341441944273782084465101544305758331261519583332902607660135697475945108048943510260664248615671261521617250026778238255092974321549535707819637475555019670302268006897334037979565876100193827454386193593685346312255337401796479452419197877142101504566717981325166522192756168315219265835427100299487530858064429038294757391177110208660077768189200749278140542062111429636064358919395097447394238420136195092771198853073739528903982988240008769562060567298067006792354872401507668213840301663945476355828676076101348728787027731015289571543503516089252201650204507663487520270900747978472974378427525903102649318505202161064519237006859806152025919985674367021457032169619736078294842549345017292617772475933642016372517247977090292975644484936372855790236039389976377765185897735098315436552073999607483933071419474062067463853536456018502551717098578382604494344063626706791507541341347487247867708099201162394005495917990982824707093353836073985815839390962805421392351458137584205998343007060837735754607584839334483144667130735751831675732263736304562698362103000094056319122933731205822012313010623458271394051212191085033411035440730373545146217479865997198791400437831778193860459516502960270705959988852657161734374152500367040188203072291162413933721047973451881316482230645091519465510001328221204793454622422350467067865130194751282550433054905097011999725696067680468793414211611597394914303831646028924353339697533205884028844101579382205693609648645085691508254038293561864328207767643154085188892009511864883404965252411892034469089299192482787687447903300463868960833041894121824084367351872912172685026256443424998425907314626905144270523502119073066867541925703998772990378281760150482586964383887893172193821027330990438751440849297600603971254412263352398712123032888149167481748357202785915857221768499097028754865069831312645072064132026802150074048326372013245265961126894713312744612852767551595332083618109485960584996854212853771237433546199787147105038322994240082113194596907377332647923923013106078408382535441821706989156149157879244081692025585515685036929313372440240193549740784254890140567535650731900603775908057768046185165661840469421173265812206192326715324961335990477471241601687984503782713417943708664170860307353021708347365930541057868814674955038655657187294680000415392282320921673561445790478111637816487201128034009802124787806543122290829714020004886681948602456883884331810036777197637129173011966379731265688206419148915337562116254731483269286934954561089519320420638066716990363773997578485016737952129851819034690928729075461879175601612691673152541358883546153670560669125838788495232337085711101151160397208791808765076700221211931294779777730263960135578429297516680746776083492672023603438079137531239673668747022696815693083255798284885057597351047083494819001890834737515050552537443111771857261485182291641281398677347510155101389687726666681268449079828201683600477862547960089731272684888577897702311102143550714608").unwrap(),
         BigInt::from_str("-40199423429375026434844861150776825132708397471941903364250719848280546169166117958111017009175151167135469508778313485099735253746988614555786590975620248088232346032255396603818698791407356171329592015782016522793812508618161984917153810102300972936054262023589555344644461502298960465755330339301753863334069805947495472939309203329299885807544382465678599711307460582979914099982811319409019464275835851130444271066974448960933649796098545566395401966091225964101636679728947093535304842903618071677818990096753011219479402822807885922670100651676830633539641562788681575785085184488466094435249018361651025285790177762185206581251329529695490141642523158759393224281463607616751530369066377702086463134779273307005152385834610903775675258022731197412842053435160425541511593823330715410635615493928366057384279371642304514948231382156777014904879359162307741680642773603701437040864313806281995321605083744076110071915490331148945276719089083641908872803547770269145860479379769946701868177761406001749927952691059918639426860482433238887260846650962548730049237077207091946658991527692280043915937217585062228794891360760446206520531652284415950868117954886621995866560047853533992319994761046989818226445958196840865453202365857550219784459402462821039582408991491339421474620200554243268596369172787124471783816445001699869322474723523295474643859272430275065568834497817149253582264034458440158650024275060464891714820380583451412292336703943828142139846147472734036444887603898311491918552998321552961129513176983567870563820696754428047879331162763709888020132290242729685343257062253348359994172944021257361578933060346041870403365208921653519122992339248369042568324814644107431690184556965851809758033897188519770352266312145484794588114117492275941267018473387891864568446708532533939153494280453316934655582178774552892782696683647219992682372679884375191951725297121552715778836596259952271239013765468225503702857720910030791704341210175154649870139932295712822259453268337612007768642331379540377429301626068551619091838696233935694440616663333969592995234005408856402030268774111311054221621206869469819084547434315713570201373321860455180298308686104445906669578923745475111787696133407367846669659776417716369238436856397798463665857708021677619587503918223326278366133835060908049491931494148134136206280775349517094131632148882444878857752093577677840553991393757255010916804751002148847265848603442085000279575294445536411175031965294404920405810810385464494938484898923746399812965718602409254227583454806056421900854353304141437044134489534168624971717278698486022538592087484084213875097363874481924689371350310923800410448063791315058169774653891020232949583406535994806918719303855873715917495572486371828035218353212623405946162884229511067720276291392442580632126915389600748286630377230539928176951037681981216911716279761190396683774126451934035680528588715886786691332565394312800154847212692687897609019110267613258152666046552294439205268289479776181166366213745644470858323461983304370009064232095992284964020480794039597662678078593522548889541753455575970587878832444555502322469328408156815103105435559876477012158451480253124667051810550075753709027819128298567508604365668768218884070792730578870958060020865129811628192423129510671810012796863333773303827340154481430178470030734526187102075826446553595281410037896163317291422088736710507398696911194435746504050742960606866560145839069915552558081388925997812824945124363582900528796").unwrap(),
         BigInt::from_str("23063894192859591173030835913873782486164361699836036563852352800904352214397244642449406941844833955453088146133969028899491573888771404667829410533968057458986302277840054127644309065093153629136866995983617702649682422536220209833284628212041883289242737650074790038163958655782891383765579048762016695048175830659443175703907602717296808669093247287824826703186293786566236927614533458440662604916994162544507850417211506838215396669298763544030593137467624291704183954573330426106628947107215847536539055115455437031168600481887588878137723813822961857343945902798309577408308409376911354703533548368865076755124644243070915347903854647092366943251173142886048587420335600323692198378079624068441218676722589469958222644559900555331556141266416289586239514163236014064582292151198297679629943210642907164554434936406704592665657249953837549508317537401569684044816621909572953969591884017355648119185606383787266158261461395720332072889242890611339752996782965227159833164750535968584248400584151623710915436054877286596982194382854444299762625019662802564599907099118023006838457663621657820335804271818703149830084934745918411934941907348010299193223828181977532597471269968460875691279585628010869614079030676302139541324243371905216233688275573425607082958432402500448466076296274598080371236924598034066604518351908180644873458231640103040865529911997961398881105674354434047869986145587512145741533137910526639528155279993413571441914876524582759664215279275364639385564377049902001122280809809358394261493669116771402266727817983300368955649981421824707262708866774813800927652072330147301152669701938555414906054913315751151066223395431595760680614879151131953733108435265446067707797467337661698589617084410957374951248239605676490890877178340397301260397326504645385287008053951205864383904414574729835729440019467929025494101386006227670381845533411297310207152996191236223109558515001850499450049427966882147830838829580971740268898957062389401141070711329500314214861712465048667616464176963786778618207238659950074857484851194218988055415210091191364776171732456297616045584248314825380521083038133570218852423732245146582946364053858192184869507939163740534742990365563810662846579157962433626712748338932771850732248710673052768931351562677967520449310517876740835385916982391519374343121848112693365193853066122346896968171696311685195592309912573197061514780894425161074249436542895843176006436614304052124141230870036296073241706232253572599189350478916786287583184801807676170625367317695552376436321228407099658312951605123892140091969829104597460929357087498783344596851375823989596649212886288274416051085869022071626755791458701218955261229387820858894479502466834634757278629005218799822616695871021729280465046245445802265549409440197956265308417584483779998097300706486045396835911192235992938461047479065758908243609143953349066093839674857683356988127553638589370490637101085132575936927219244446403074067109978947589235920721881162527387106991190762475138396968554404665117514231954387031051353381638321864113069020148411767989579802124954103045098175890233490311197167687684755104533055330281250230237216482919974114037780933592076307955624901929983083710851599667771011031265286819513332016006290315397734474968044946301491882628437620892012813859211678585260984315208702548378992041235818814968204642191558732426499304704854813907797574961688617722695100276023766096388513790347413858267472533940058171955610336424646782276723254586570447208").unwrap(),
         BigInt::from_str("82815943283266288721606087193975283987940740937438952778835382312767930042874638313350527346912444241126774234786291648146343862022529767795106555516971118710199001491236288946574166947109422940513175508363623204961455912283221036353409348667988864796952708769611037738714480879623768724793885092519328242966714657341083172269602470615169396974943054322116684460143610979367903861798830981742543316478077058252686060004485560516565657801909302687298252787063690129161318147468937655070594529788166559020093818257718957500577015228188889847361759159853186061269075273415855615431364030641334560480786871354836242426796209377447494788719299768111406967605804723241308486960656347443158526075092845678727420646260267609619717995172211336723249712132416073863407083988007791603937121351334481838776444338921000781249139620510104537219764853164090642530072553514517086285028511721197216956906277923750922427989972202115461710872404361958988634339705079769469191221660260376664320681165403537042472544717977874670286354892612107059920839752880333981247337250949226467817620021523728997775718246812634627693115450810586516737228766244318711421936391593321330315409510908880437503322254898191906327782353409631487268159755507442773207840764385578724296063332573061534475377956647178951485153597527330102533131347773078270646224423038192403619412669605229343381780501617858433517168987186036423436875043029260232944406501555403779724925105673738692159327988026814565803753165634345757516843782961431228203695823144650512165168574038171987079687122809424060905757404774952501849267792145645958654493366870533240746079508765793198119822852387387362665541001421842926648444078587865210049172794619757272945543889658321507769846279141231835358234287035598267418567575811749075007031096201884936747899130123416992099598396453644484407991063264110509459943196507906816053356373100655787143255036032942617507381315950239472678934126570731886323118346586064749570439870994086272015616418026701677223939198177953195468527774038806248149829264813312833877166465942542656943096181623547351471330298022813982539265677208925282736724512608895988485204262827558656206154703112195671994523380441968930456364084990832810637585126796610132777596868888125936477445869661091091417215049168675336581086782290237998011412349073503961531997022872318709141713781764830046070375011253962417951819008993837066044969752839361625406121643780943247058102292975259123331481844858624887960371220382563662709298452578463308632421001923938722855000498873848864415289498517868243879151292098503541020426612835514583629827129269299676888121021505149203965358819744210570968714464419943542523574708103432264285136708637940632482018474697111207320138637971422564248722964776282135683194498438713085566697174420071785799333749822576799757257509464953417392898952092061048176169137821838632455615550603173658036474577830025417702866884684819114135026595543482559656603890716133199938755935120373662369868822406922819952941319439716765567977644764416657145382259873821562475059052420511843784691672284391660470727396590790503734432467833217972748757122811380956628970243274937022902511855451979471787364547673906565999829341570488417895974974740278627497988478186895521796201064299126534316554770168777413277256571295997982002070361714634534680918891188157729342281936254433048068397233714110487401523135708272185671225743329140159638070677306450388684582460495547059014682189410752880982555292985685908666048247952352321419436").unwrap(),
         BigInt::from_str("129960017263994841577803130018060127323391520569268206584664330259024493784137578699761135346920720747650961052435817136818605064863743818168054215287099273457302382431496390700791464697162018046377624203872220828112251819523049476537089349938220705905825224582917492902634814244939930604304980085442860491872004227046983722673061256229509910663626942912227586041620007761307481902052658045447981296728700350570132301495130700072393602531615152157637989193766132397458032356087454369477504724248657352413632931054079823577342336966863915409072448838961417766406843903246090505000244239218255479253214141914882145891226576154225310208959529777095309983436325365121352717878321623550961110765921541169001531466622990360645334351852796843986214179047115742057829297115931133452318846745818086981978067522240681138060659389103442170819599671342700918254093726076525397892610643599097464172479650467384412979133023874606230205400002879286871974024069049543194651279213126285354559093564471386753127728347728874639726887316640798575790414522553926803735375455521119449169225770467654038711632205392177931184561058211957075971258640053611431653934728810273724069202427237333929082932595592254670641580312223678918425512357022957462228302058186240509307061369722798068999381597288873789922739495498248436627382184596938790426041251575455346826259046597082940260036660554883416366738367700034028903793598267626464344348164661553083803652588842006633208299741155854010979754156693125279065294331057760897930759110403174913033671489851143625703785711739094879943682043738866464724479359280778843818547506731210381818324193926085290785990846071785502586162660022378245729320538004642101718210300547975520508869189774528086930446092447166871552097322731201728585467032778054628821968830165633814484456490345217652880485490442890528076518041061597000659791876929015232032860262364932447515052222552836896147957307906727339134887620731035814362222501475087123629210341035914730111359550516436143694127525666683329349571492134452560336450015220415659521665139357852440433562307966163046168978922323995586634235260919312146410184816217745114632631060196507544529439209581547571893206651734416361722499235806423525857469766526348355041610945857122084189768174541857249998468546328169181340219156133086449598529937538118491097454405537811478320029028402573659347081037035373985272310229717250375943518504517482736608927209749080909631307118242911197525232048182292615869733852784245326109450849801940833592576911796208475731012819746539526760980423329429610326861853479482405164892539689893697128721734451699260340366360294572568824557838510430374456336260332276512664869673648352919942666340843686795027763265374664288802856938797445195215642092679793751511518376834633997290212858269386997026306283023385261479968495350373871050975411655505876290572949107619627921434920811660751914072447371686801796188509141527669238695021260901817593990704925771542448806400190541530412859072090407214724468341387777287092817623594321341310057739075723141512255631085127180469482895771641184984150957542287577470947292180714811644190786654041673184516603462582715492334809407627254818089277514331452620823124261365687308271393301355584638923818607356356807027393304019042406438680580777432789996883951562077168319543938980636034744455305808983341535870224500902258875354396977756042727047749469177917812246592515448606772789931131381904184797280125883024296920099512948815744510822258517137699988822051003495272").unwrap(),
         BigInt::from_str("157318857411811049942358978851785103736005792046407585858594324999860822018489888486734547061457315358247155658205156535226413839294230770381949835430002112741760699488552794234973193368326485770884744029284554868511973324669248251550757503281990221564610793194100852003352707013081317928532703447807410668178468097189789827849608112035937370672875275243971186910530678879150736406531535971011596617049053959208135122007631371539607413506755202083816718282579974910845135926867323187079757088793763564813242869661109478805699748518409564320123726958684348727365142071054494832063285106059329021111821303189589925760021297372257052343071982483012684722957014016689622160060781153208310802141684729974802725988988342634641713401829699610592790372908044608148455432623228411412239786889159986061871587146533240682182883676518473172381949434633933596835067988620275366112986912547292253916248195613442204653262722804581776208793674038204428362491736544969338457952823731520831129144795379532199324546403701161570323962148029242523225033473073522692677393471608056329540073901805557626169122620185702432220733776078461735565834900905535279781016355327579471238337360760650543268415013715685425519183718875509972736807648408086522924272942735811806893272581203266980848773942045467216015937710012126422807288787056040064145667008441087614139631266402305495207891588253122701386203730477886575225224164764534102702242100892178474178401414229330625918344405663340592217440497270797446598536815482175341364061667444237139447771306893788189710451174613936886219006288165329777927024934889628048898499582149619193286753594661045713496450899632241476034330039584160602757967324641491787936947775734615313765992572168972961137064866292606465821693931593950093884017744449112130521553493849566587684493372305652127255194246468380270033905737170240139352066666222979874910050849487125628537938507193351777047029552754535809632615176726576143325220940990225091947505970490200601723970770125987691089449174008474247081643825132516455773110304118975356845785866493639394369714487412375863686596757959258823298906324713468810330886602102409583631432653193959422730146098433214232939972604917545063005395851076267072021343399163774067419977165166899272804453432833119595164796554626201670778611205218587100052347431458709073828363812013509683548372235054182030735764115990879731470354223691533793986195523447733441955009134227654422708080943643283855524202933017201872575152056180037995990372616116431301629342671086597126650034032367357158311514404302401209188412705773887825066409325557479231305259526221955063494195367004619690543151864566917523658356582125570841223699740580031106478348813718415943429126163189919565599569442904038698052010073189072096438121191653710851614902141002695528783465592216914087127271354902386199156560900416567434166896783382514707075768120661602175826524648238997985093766174858817549638179804511648935105352138175495426277423345761868789931398243582368288857330703001630089086307187209608644203950454208895287545748767700937116345488177636753149289636168934636155938271456657709808357282489907649718636519155401443712154932563619147917812527287851072069421830264227219243473137137745032591132212523383376400464630709553347941589804274261765147655843345460475023406874278851727998354103597255141821957554172767796144619691251549952540608707002955089015885570380386529187684009501163512639571912739163958847537656670265255653146375115855344732608885202213543864925532").unwrap(),
         BigInt::from_str("160727326826345673265933494622052646366281871082182011207978950454349983621713919490680797239917947780771646150986988413900897769303331870818162369811460475433808641155168164101460761450288800085366329802935994383257605710256975923810488599400346050855437632920300012267869300913970025298890470294449738310979263436603013147524992157222775646291490592773330403163567018383583000538304275710085425797260498586836071428260339601220030494860624140022374322197054335196942720266202078718473708128194762917309339151867901629436866997465348074713899889027798338933296089777551714965809038663685044426625622890428335206108044225771345738185973876011650306201833851248946849150860979497241090019717725677372372811511916083343251676939006360268797078970032233781647725474385351844745402845230533790336040533482374579088081542701117030278331292534922050953963304054046377112435817468455641192832086089846916401859263792783527900956134293062949656560280585706011410902942289020049217235722031557167668631346692405048218708066074570204067462546377138815191460128636927231108474939136727222642409523529949332579525277508743310487758597596798673522833317062459757410387713954229739694532179950094467422408240048278964454872297739433181597331850419271610763147172130203768939855471091227316488740920395030822902034332745122520132076489340502616531735689232411577716317283644937140618958799076517285087513856245224536342668004132508982743234806336805961294847256169687901027180192522786685328941600890434527198149817397893853226238561059367254512635170204919370869813476647712548508059275431631000956346817068822870198458950850718595789764286657516106060043731264566387015378644581104274000881570922812657036743836365615628210078167260056928473593289063351372285634610236296278842475148835425612028344981642892765207868339736634747473946728389557159585282998412037124926808408953056899451612920591173191220139887219506744120009346250472718357258320664596831653467928417981328115879930144243858142897255970598142848456631459329778560199588204438759702343895448524332146431156128166613463763201510714813197263529510964984854150064088221454370995749762551399552708792936776845984997511297380744975071330223560747441389307206544127566876137367550228760794153647249641637306304627212842299215844633692592549261768556593025492044983972474236196267163748323134834979184114286157412616718562523501309765658160841104884707164369214065568919724119207695602585749218051276677541900045541216702255712457387394121310078152155313053789483382984805628075990523343282679580639171037565333569393990721007937326671724033160215135982750382683992061084148896329873728967697961323302113277960360779749863820901468211248360193681504648430248195244684299877489488077813689014069256307450180098132812561342717883162757166250157325761337381612841480189792352117144673959623643326363640476591869065758278127265749091498218354910304822110738490815082541531964148951163852875383893376065826182217176401523933605713278684046958725759731777639079336087751638243960315872637465453478713573107275766475318290040216976678996985205561827302851576861947925405235142455245496149424949125526431590363262803858222467943818536741676795578719821649119701390613942969928500464626886044705261865440703733567422980465316862025774423326944572873651137817212254368279866674299541948824010118150710313929225789956752870880879103997596609122915151185685773493082120237330688480975934690950122384371214468993434777480379333928329930708828658040").unwrap(),
         BigInt::from_str("139666516235023288211350018491681865223660616577095431445971152484805090227709511355527731874693054881128496655364769749307770309515844990573794143405139381693840678632201945948985533183937572450121871549805666934528004057139966486219265893620843918046639166875166497259940176562093555797611501983070018166047827677018417739557730136707351371031897827786572959417833125963087475917056269995939303439221755089613396743102665053178388913276481198818564760639045722481399129978336445274126456124963743251395836802657022642049858479782499747258394786870159907333584626791576685484102087102705449791310187459209736274126075681696787460390372577963799469729455239496294946878300115518275771042453542353419931214421232108343288976956375633245443484711817492513868866884205160318586710056572435551897990591098075170557998036706918030424754549453675656908929649299041404784777568142986219993337941660573772956023166231336191835675599087979834849088830250726985327254390100435935451457315979007061545840796177567515183013778053991492621194647666568242086410939793913271969764967175418295713489276793415071842076201984751674326468368093918950318948191458858080084961923280892391924085974134542460846709434740342814243047473303598975784631915255297395606304283195838651399772334619779001873411475435243706782567485343108617873693695716569738650118970538057705218536331697030081693344919243234925525385331390664626653997528502188131646089853796302887130728794127462884333293486060397375314061510212521557156009764054034786333004253081905763380739085900228196747064784100661898262775893499086778400245463993350070241633717811568441956268424914989444274393704192873000511125579423546272706847351354123100950773151058891512739771625434124870386445709485571116426901352466980232554480616927547858035791391893280200201021760182870529282964505030887242559237362652401949102288568356554830577293556245605050579235708611554683197158981732682985277472263386906592118467615795888920573350022624533324937468764991785842110638526949984319626372722770560944498997611332620670909377740891859318638547835818474021604238756054015598650513075163677563385360194305651619168467271821714827448846764477584770124583125884678038994369954884239642467754459962920323577858921965372419356583409344281842414494326970311105865958527470344409435212084991142628555545602720959639330632464043780719809702055478736494573247336709416559478593381304035501795477960787976701730003756239133038818597276328516042547896573192123699289305962196574541495629792666647509383661599613400826251754552013840796151037465806497934539238326234232661212864664308120231073692827141601610317063511602829412956485392343755162736584006469271014622922106835327098319574653404032384227870920686743811291905000530069323569788300744001418071152737435032806915031382406962902189268561319264073838913707782938231178031717944610815900574287231758022335433340511195328967010260938193497975170854870764894704570939323073113307161536334187592884539029100330128763930396825521979361798264670173491612906499960576097692833851912201345742816828484345284418540328409412856228959629358948837029490932312926664065141761935613595530356485628951625968950854116321141506549784004455168351744560449673667651273591536152691821550198588393307363885868994876513807738240426852909856042218829419743545501486813017912360254473631418498414506267868727927044551235164391826286921801984899801243650467215779739586182580619210627585638837167187417646755772284749687114923044").unwrap(),
         ]);
        let mut capital_g : Polynomial<BigInt> = Polynomial::new(vec![BigInt::from_str("638708239932366711425188882448605599695932059889512809584527020036603771959331792003292637552839931599360703148456219559613343411460921520950241436164722495225868320954524607204338441976762485072295766235505921887200681246324080388009931693795331832702110062333676316878296360283279257058355719138909847747096791092810293141458755916539153881039045261921320710001169751106374000284255971432089482683327747247388495482144058065207406866775986174240703086664772482315293658553706378207190934203875550207408578489078410663016576426076397573484331019222754745451632385332250132454461027447104946227091252435986554164860178425358833094324473227639629480260505390171333481957258054356030217580462501177398166722499772740903311935669827107779891546011184596426728054754423879890928796819102045844247470280946061885404057394791079441008577476913195763617757563348655455439317666324833674462652313704189957682817698092791191452939014290791345406545956533556640438916607471535950421595840538595301548756859874720348935974392016223052859055826469967346247505614407462049621160079913083898577625454846522226081226259905628436054927267490239742569783094016971962157434338552385618193185818081300068100598240121601898464807729539609615424134952118743845352768167668589765718136270988531034718343381468107141090640169081358274470989879516172777189051811363439942317667762134030399644604969191148479019274770006936271875451319622203134494096594388538966227297573332763899524480602630220820763786475707421858697006918142649313194142961703608626113653613827439427298221994326515297572617100411200143180371885312911622495191290924986205704183674446966780236860578143813699424035682982902545721797415169313269076101687289760048613145874075380890474624733418633137124187947049972569269725648283403330254746342819824295168946983378145124288035528600144958740980746034935573221483491805559714476849518487781271394389644618009451746151264027704208906537995878422334668176271868833853315627987941738357068726141456727068088650575230043029697488898775340030704202691414530874550706414851659838037618269961334593641359158539503774940691478732853825678728393130426531736937746958199577536374615579208202770502444256567755637648462985678873869274259785254357319772239737130092042134973897538757033466599823333916999877467862127589360535185731906938909056730212681612264635450180454286157715984114175899613115849288227142173157597042854935501039927274019876878384248627803309218062233683245071631088627312226968120122019245603385469539245333348172099541807680513026470812746446471380858803450176868971613424252354243853792944368109744354297615538555573646510003799849724770422653356634120055980997192773671758359685033893387665464353498231183629480224704782234215746506732227390378485603640060399699554902203602965619998501873568213579227551635590103495858587446934298642209510703887141893904189668857899610007754771792138811684570100657421161556898377737821273899304447606372065752838583165396113784157551097776971185019012719291966116321186129320021265087811429003744914078057246944253416282322950315084442059474381557024635136852518045576044921848404917874572206682611155973712665061755081741584651431800522608594751800350463584585404515384775131945169689598298853868927440321580123278394842913962102802952598120728921745055082841269584404787196614471401262465291040473822280371715846424914434703438804146224980254245426518777179772555466424237615552397976636288950986958502105029239858466516410310179").unwrap(),
         BigInt::from_str("357957442106052591755302684406728607618437781196580059495033146167705961351211718090673547738807632845671409478001748875470714263782671542453395869006723361103993874133624948498371241166877129629926248938972378724663704821249175615167884012987350540900916850767606147804182458328724271384688661961528126323036215754115690393116934381144717566248920146581161376762987454240082211508344171815911226167694693248440995173109257873157442808944800681113868110158550535760302755855574666947446268942745847453403324149684919987393810749302507773221647928036454905978977042119070054527205378155043537683340733853039832988828517173743081779588726677217383455497490248394815929699091937467642502037688894319419258291082314251551467451289843777783807389183956201963541047882356951898338525300414370130062485077540896390966566736655633397943684638449796084217350546004248995560277878952690834259602266145583329736045791215452932053774031225442595052741515383728458077815223133437177321401469394849185084031616684608304500983665597377702505754288499368923386083055608659865722662073426589535027338142830225381999195968801573281644215957144109519265793543934616601415232654968165109648029496402512493449486975077615184318486546680702779579848029747871780894907599307757330206541303599028331608981429758218315333223210899660299326543460118353628057529167865105283475340907140963949792212172678319931168333440727895415965474032910804023038208071469984394522205041000300350737228079640415288946599518174810257631463687172903822236967771381501594595570849421766226581628410575395519153141863021330118673355892630189659147635145848041716807323752927174056377974561697130571782638523022196397538033871345124673273544487729847435746628262166465035616528586258299552837434044578712986178521463590938959071665769388236101449679282341549790034614101041618075103292506476525328609272590287017266633459519010275731122696495978625971539478567240868333364444447258743061211562343509928922777905977320287119764360544846631592175122735542902697474560484317665227047955959987656578408148577467178058508355632373656476615430009284591690933440649316484785013226160782730077772346363894821401625047120839738104415464324626771291634052315710768430759085389199790177838489968838122031503744061023715573729093877211729422723733718711171533987992949621396671343036906627882085224308369284506945903723046882519399751178732343151681574428086776923262638875409493827116591223619294227889059079556514338435869239491080208855827324574251710045315284005517216720259715150841633445421660206650249032259462381772238234788166949614289592959509387399821219693823576111171007832711135268566863745313258423464937764239449032413740723774428803423362865766862231970799768906729805941912659693284580983103089941876373208070311170522793483028326514548599472992011315348341703218646956779056861252091598965128757969681547094399323655411785182160256204920463881824927462746700942141217060969769109209995781354742205733653519520752666978553219562335229673432184543005041848241135177600661034965364379117764691139963553198564240602996011344186835764579421734766127713181372906231857591547449215883381001133249837342079472887256324016563455973508323205789207502539096141744106991683578769699855100998640256422978939548980755412867330378512423190310618530188137101997394741478480890568433825234465125456091343864710262313379176415216287921505896848821549499656343895247373878812602604378271874838828047587604635293171716364535530763529").unwrap(),
         BigInt::from_str("438966594880305812061951835079983359125571298624460297232882710573068046060576128889209771856682690173269213302744521384093809866628694950279506051184378492584058580354434416613397236404764452934765167152966474388990086773269415527306472917202332461331246902771773555667719321737334931232070739123457767996738716261737568564051248824236622195569275790514752558049435118328576177589139122644892136282373759843895776321675685188410523493356479464482004680110983115834961885574363776067457590212527509004790477246166539479030878634710652663566998538376545239583118991495701909370447563776979674443359950116857194515035118665134552622636637336923503908265432959544283408682786644943658286443531720727161773287970440806480422482894798936857110179278030221874209597135188716430761932041089377260407321508665222526929483432450686575218704193378772287628441015082169139739308526831190325825561488670125641179641327640960530368205591210431986729841409234893045929015820756420442465474324490735932691857821222934197760124445582736870618421714918501391343416798065526674697272968912583767016141437053911439302080199923339801334941438110874584218788533354731754981102073549494930057213871003255363235942922643233727539863340805151637801912196675173939370353208072760445909753699006537668645876186117808008412820960988846738833676262787589375037983820860431599997311426012790614225274185978598766105608842924218392159322060197474058768487336528655409930591337229770172214096215757414230277270946894910430712343396515145011543923289478720247370769064071175509416128768359902107157793602267671203453442649673162523102177468027546384030806319011999051846955087940037576783720175306686725409883857667378646139295677661883567308179248167816708079279928984463624339549816205119403573576090719275109140882265078549022955098576557301429846714296015465798529062040296570141656125511944343661928577867822813402484477636511526523306199262832344171075496584783575025004155480827380280098557438584899627556186836690821180922030503612257314467459442535296329568442060678226153147294641887886521644434038713412372468581198955398547307462927734705393709939379593886402189928911608705510184561364611674815794788821062953206936254368780970813691576393277745769654382026981974838403289458627134428014092070615325846597175637748114808510596232878788688251292035677973954500895598275934714193310503033303371499448759604601846558566922828523637352553674984832227896552947207624954159592765500230526353056169395831051020831559917168213775129672447534791490159240451640459173512113077308790976393853926542294558261607771989828151989356687760671347403948389645190560855782001680661327821146507009524232648907512478094820938937044934486301161459987129429552261787158650478670270529927663564216444402058774512158159712788658234193739436489634328517310356661126379736758010655080921724372280726677379206956267935217156202481653431647953853144349902438884807782204345070734979754035979732599136952433987615539639936042958668632952103886589400576176039006102487405308595657371027930174147172448096410229057954560497306580272635366973365910696178668636289477641036960488854926405223727012755341527996104911606756756556020428736181472599842739613358254964824040455306884052397670202305024775530608059263096003762417411423640277594826091454453091313505478936319483424221969958276938474179382525427581138562020437427982502937397842987732982808671938365724086003524082795239293915453124018838028365738815673340289083439454").unwrap(),
         BigInt::from_str("-17913490430078572971550255173536329671834618912592071490717636664315169008729929305781269342550284640832934412569874333523260791344719747782379926489547489594575817831323584784272986240330194706393886171964014414363109484096940507730447564628407698446372537903582819632787715488515270222201581353461060517906835397520426857618418516414885589011476771616450578281396494495046185591863457456143701304160169390638953044823039522748172031339203334002079997600049335769085380195585402241953968546503898298852296284539355218059652613164827346022100654758467451435187792439737995599561929910946626374799600131667623117664281974580135215475771126823915274260797408596881467026012321601924885972112107552146597350653352338539280034126428824239605528525576684088079507790028045014976536399243540169338882320554847800743533900727521619034213639900128245967689129721796348233592318112995324751412566770896726719602531905915935137196179910858465231309054923523252919612518119213728916755502883749501451958794211305714467092211309261767303282748268060592712374356344386476128153421123475902011008653430121954857223142139565509614650510165214044456741903816804988548248391678431732224670165131221371280284893361150966010704562380203256104778339126090740172087536976994173978843982299113309931740831553519204272729604163333351611939861703112983998921648468673210303415589793639767594829777716654462229055762450779403293910903176958977536917020357307908715108672920077824819549183779164955938235508720742277356741316155562728055857005970665319797747077169363571558743540470395570661528799177357579344779179408350448841953702034481419220846059778552028455750750368298731341790634032164656299685499193516530513024094317627210341902449115043792203954765178857993166458669009426485114995671486169653001566644900150152422647875178131661038230442363284323127063387353859291425792074828572518570493022358116324855634502124934464887952951503907547131955359156703806309265810794853526654940604750592289234521551618886404720010352620767386471820907962824811847963429755229732559351238612218711160541318596345351033485975465808339242607019645982160294849862707940311465715576418965795364570686160102227563018600128689532457286221243869524888998886380709357634170715796655909225301091965556425031785473751689292042061504313138578464874735942970480139364820585594890470549914188298653562297008886783869097537376377424473946882039160581842521770174611320626415116324967726657212878500552805089902408210951224949351104012255550400268931170066612406797593237193600182813612351872724925921804187735278769494897588557540254671739653034022104063758479558907884792580439054809716294174185505896383471282436842288390059365677679994357388780666937918886200611299825064335524903335372930894170537366661513369484945563457258533321312033105366195017205796359050696579484240694520958260591681914616238354575304030361709940088121034475763682337033068288869006615230482235690303338841569387689207445606740242962953979680983162717319383909247196215234894054588625043871795966241765555529408408692744832374392562268712006534903699733572177338984344685686562208086320953703280520172149832100087822223099133623757719679615724209824894048860563137665427288101307870031310166638514581766212328494464241987927454925724297526353683894053045855406617925231818132187997187144131594024338988651191426164609599927787389301998976313932790232890241559433512950222303636513281937042332744040681115769690896358870623678035367114872245").unwrap(),
         BigInt::from_str("112450568190774082130023493844664098515209968205629836716873566696068142926785237489087624404616463819926968966356619611440677243147431874931855555908759978413789033153905699485459058893888669534097025076320677988056121706014726730938243598133211269605811375047088088497528385994152343239496232667090588576686863640797318029750969163470286027094363809127852820468187907591394220225593925499242095427523015506538266141062120597718065364813458375184502582116988367873230539878027421393831292360006174574451998285968069169597498160776234976747110865463917748628211148920652955203652194292922826988978280913989590395791773995312082410833156335565264366322059962316739283082022318428210046104504995625740339861010359220383950106740838588669858886408143186190246425810674614395428807700397926577634790810841252195384739974227265004040332004717007007584262613649296211193504780118417689603084315987407954906239674385341779048634351821074882416987962498537804743357506191923546751842860398381072432370730357062258440689342001094686593463148933258112536939612729731271855830002530521122729425622459601978020415381598856019906111982711835296327970143670338414502365455027589714676700371801180326155935080200957406240996484745727261800966789525308168134318022442887371569115519627022474872317768901052549748124628020019648377938296433398195348164117871265757145353918920675193392687733855523884151276300732230302922812320604285703204757058504885469371514650597477241503092400929761232984348413547530856665384531740021134440352926429752989393005827093086193049289071358601966844386306891516444587139737051011866124943280010081436488239931650111980953222714794795335355454573014791652425693766275243578970341750617411191061749154796575407175948714104205531649414390996821411359248634762764300149959203268581064898926412374780067537106146924451610976508704385438588605959489349499705018431236161120264435082744857921375639800450369340214529704874353707624070101691149175718023447616010152859979351060468493571839721160217209648863876580295712216622922866284605321598929024514109420270649225862926516249158801430779928507065089978740640078041637004157345684783586104097482287528783402163751490536175341522757790471173596091432396867366625400358383610752106110468533404038697984355097404664703878168202659873231422399524384510884434416252644955401483067819272508153799654183869328121008475827637493913422539492053981691259407699157241223591983984501856342297141356509624863560547689564231309265484638607001315226835110636168663765528461257672726903142017546646610239567961814178106451571074904209798471677559819572608691964601539120943027430839498337956084301736575371053612451558536228430999544761988752602927288256415831040972684547051372499734407799766567667679599304423471590169602856767515556654816815196518328708499418286068214254704397353021712556930444147752568105585348737266247739140385387652436485837246315623378585875425214673156060818269967348035167554397208133491893555769604750157311980729785301727671924712069801064170301172164007360788919538138284466023058636706958360674959381632907636545179967967641733808586845290026106842101908154704234105434209084262497332421209580956197877712031388905001824602547042530843400865703177157924448255297562081975468336776717270759946923509644784731382530193302215378457146793896941425938025955011921129717593073433021933786152912134467729576526029875053663655131019188989841976243345185956677950356691049505118574220198289595304477731778").unwrap(),
         BigInt::from_str("-198931678970070611264094888386680144000292585721266297795539724943360125923791646271487564454481577713986764179615266590005505358868130008076143280330530809600593282556715217648117345661712870397115053355020513708946852887875519472059279327466385216680695182760922122899079117860050440450662733437203435893319626338196629319243814093214207986660832833727633529455288570818748362117281397796638260938923706160060881099774964321084301784246483940290039638295789104206619267655171955284294489404759419708863365051666773466393292796409462108075339691288607253673823448579449752766772078785713766872988462017739097507892505832145907022582960361343938284748541750575718041227580181202330788375329394962414535261633099405412100779914803019690773333530732786784082931596292485824026325065011893736103249457734479995711936168269159191203809617828240164158257417030924237578267886164756042810621317678173741735715265209828664880585054742197249471567655347997855204399524563124793723855390113325316925950913098747111369506448674407717474440323993587845865089537232570213257301126181245321848358229330808750717960491214756199146701406466324996143333387469088239737090788756417159312796129674215075064002529708192317003051920911859192664281957566105301118482693635051612578766958705541817724210584575028019881674966765684441599595105063666357208682714895347386656170868090620976179648173377468027747871022256354435526224047123242029843192335199823851083999590158619687024604659016432414174984655811587289461317429583657346708907263235247242665513108560081264276209969981628802336818612836509955839527433733371585041782744845943836757064494544894645938065087882491936902541255028556569589072569890639433566971236711593934066242537110536310354937064543106469617114066176457279433913344177264769430213663127907189540465017110632491224778814206784773885413603958066608446649225566294989043971180164460752891889489340271889321088495861543925234207792760833195740118046140958499039902365503078868841495047882546754703516171734939715612072125438797106844304720927894850994862442099591620463606280842604591452243058047723119689620396453701482715939148915345012291849880682571142412794914747242396487032185858161769516408926009381067598978070476569935259335548569229137875834423025902002708187284331041540717368085962863173412396201786991860249175683488425313484130740831921138867620740068452041897647085770180538691235181138471366757725198595889857154225047798088085488378056258837899728588202669901860979337921165710061407625450207291062654783758698830811346660138508036001217141811224777875248395712183782256590595789232467768565976567875009595793612276271024789708603197079761933661880553434772814674665376932980978843671085247824559213654704309843600022575961510359376741792136168039133071314223282218602029911921893086341179651075772027624116691836760286540567079303547241066414276589349681470970846898731545273210386571659640502035501988627181931040755763406297201702121075452870373239187418228067098064272128154711601276097168414425469196250106545275662375269281542463747103910997824610706644230339084754947046257869373552532480567521403841202529659871788642806666514495557400286551074946890961784655438142246039835361677602169773318980050139513902767836051279244757243697561078132148330786566490811699946345949439653672187471062920711132485486994748199355501817481047429463851029198996169625156049288026606006702778142436369736281904478119027545291927486836318614706637091668362185618282").unwrap(),
         BigInt::from_str("-464304901754039625155901972963401415991075690213387745367578442323958553270831602973660106216711579367880501808444172391636268021249592403846837899798265351001376084875628741496647734430672363496570619901347766601498772819572773308657480362843661741001626693443322711853735421336881902729795707897937125682563650018574673707617917744477462515892439325954056969321976073734460721630626205005637731435070982767922402039954609284591216222240052008300600743018778601347059171123191191332555147054651547849343444545559044655017175130100020964419157941094901179295588036198202849396998816575249919048910213744761081491387936814935089396315724729667663574551953040703900075110143858478167314625380023940531565892944407152983496811659625714209620132029148767164178806677226331220433688875388855933815914939151075166851178254704897373602399049941679058364476627828633100437332095314120822909101469019990999547805453743407955496772808282818503658617200892883848456110632942536921456692426361485337428946509407726612163144337936909030237131990168718682424829681417631252522171653069918454636051056202781190645874555415278216184641657090131336452027008427208332068008982380686117775146536481573954447160586077898540977112417222382865869586844700560477015396704936943840477842243863540079398212989364562814370661019101206551604517590832469974904682424370553967751286986877276987497837783437732579325735600184100509175928425561548010588710634541146663118840282355885665875979675412303254838625432556755065075192893088712058496751274564155149858191588195281775597980817576426473144225469982467785298672230009046980255471814529384986499609487469476846936009115911346728060405103027607906303988864591802710234953089303853421483209119206382397346324309779722374696597206514251800389525710464501563709346477113149098275348042633036474896527759426350325965376422625878578508136496818278605862810036975435006879488835785159427115010102347114354303385340080198160520944384608752234812720615309221701916777853258528090525962146760636407822408913684689854643169552406817690929996740940963461388977988371622438294003019404779576055559543378018148899425155910283837891383415246028284715182727887592052990625125841553588922401484974791958019208821922142402542627786582214320549041151306260981998548281458049665961725074718795626871038175182986434862235523007368207429113581473788781333345190037565372522900866784407406322008404246306423855055522559855064950789398252178867100706702624375120546869150855695781374069141689752897619970784466797029834136892815988370911618252180933157441877940867475784554117533571419964395506748021824472601418293756353490352329659648340216219987285727035972890989338425882956295271627624190356621745747863210734414785024321625070281801716347183062093717325071437290923147477305546999914592999648619855011819299931311105912824032413171789047139574079631139682270263333166063364394947865163467123501583313625992059945269698660570808712550165349378644408625019418926633309136794600483679502062093384273871870592094488691701522907587110289103619968232547287574554995619195260220919752203712962274577634724525489526113581554302047958329243357423742251920833221480383867499947926526884472867638517066572063525224726226586032195235337806021210093794104273605022292987288110882574295677128379285965311574469179006702193839742104104685907139917059633738687673032284053712208701144388257571197553966170743397245707073775581528110405108857971857916569728090608706891836566704375233").unwrap(),
         BigInt::from_str("-393782769248744330715547506056841795904575235201972281078825887775054994830333286693764480737209760725659601785542885680429102013119208829396171022436781080452235488177109524593837722925712025282843580307132291475672053274107638901943019848327877432566010646757232518649619815648978529284127318971323145951266940420789386323070111229023999445446429197207548199480739768773999832545611888451260787930085167997486022254023651821662359848178785358020521474562150010575218174660986063317498065193845059308647269655066944500920649703408966506715330295333088829339085403085114771042760659598264187783709982871144631315722265435921933299525866220730475614832175460401413076703408419554493197960627212462628131037965859969174492415694274451570485269460437969124029643208800267404175257771126214888575962678256422837272441266629343589561513756466857000243984869035008050958771882997173875952188110659687954976551423601701789282468579910792896008067730059334088093287515852522958603563541367443100930919677122091153227560749598742453976310577039634708824755244542416066238977934883741347790734753742554369877718110512364172810075922956910334173822275001402227260108353265427976049893024315007027833182903463120267242035128855970762939622147757972468428632499552994995982575812708503520348469067792946083046622322269871541018576971052155506647423817557017182062286443979327614427617732184953962274878006827319683298102979478633891854701273115469119883874002644876942965243100499144434201263786097610195050863839920713647626036145178529393530817458909038778142010836762879969450227744127647790458291807936996336527241262999565490624019695661300434247395187625984237411489969398587300083205348615323900818343737712794070584840839088771117369483944610173282994391531194691934146301152603466029106549491127829075822447712672368563082951818983629744499599177258405384904896725089318063682691036899279918389464810056594831537669906039641742921416285123482089131722631769142098209443189748267352529003507030229385518903334534588592644788841857654843250356849409549635309889327766898652380436904201806692197510873007832377460749973213619697944481020409800286008575794132311901890597574673457579859026370532052102960134716511319056916621556595530403992621660070662572134740776789731337449818747551536040387813568240253963377273141938056357776324050523575641180152681955771911758829051380904169583998187328290620242153737676581620705447600135041149204477398730546130052798070621839137148857920388517953516356514187196328945815022007067881313642461594574534679233600219984486392996303082319926819255068576408842890863697479965253864363078584315553802246342180998341850284652612365084446956695415293408941522604673221639841845308296764524009683582241669567951425329501410690729081164163531646445116805306086278198619731135756944206439201557841553179114003670534197788264558397202658507057017427956817311338975039153596918348103328693171560909310473931749676495685427290253927310601010756689983230959102428713384658593341799778164039390611045732412156593246188213109205124814207799508527293872230114057490322264713450433894433051413114961400316210430168611251475313875204619182514980946047877853784855081232010341109846656809556450135890562938637533714153139918845342050243487240425285832686907937143202561280541082111341061680674892065271064845315065586653209792785388614987418798262518429931444250933025891473523971134157989250567514369727217584190945611862567451026192380278431223328710288361177").unwrap(),
         ]);
        babai_reduce(&f, &g, &mut capital_f, &mut capital_g);
    }
}
