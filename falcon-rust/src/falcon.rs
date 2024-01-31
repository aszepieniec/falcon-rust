use bit_vec::BitVec;
use itertools::Itertools;
use num_complex::{Complex, Complex64};
use rand::{thread_rng, Rng, RngCore};

use crate::{
    encoding::{compress, decompress},
    fast_fft::FastFft,
    ffsampling::{ffldl, ffsampling, gram, normalize_tree, LdlTree},
    field::{Felt, Q},
    math::ntru_gen,
    polynomial::{hash_to_point, Polynomial},
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
    /// b0 = [[g, -f], [G, -F]]
    b0: [Polynomial<i16>; 4],
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
            .clone()
            .map(|c| c.map(|cc| Complex64::new(*cc as f64, 0.0)).fft());

        let g0_fft = gram(b0_fft);
        let mut tree = ffldl(g0_fft);
        let sigma = FalconVariant::from_n(N).parameters().sigma;
        normalize_tree(&mut tree, sigma);

        SecretKey { b0, tree }
    }

    /// Determine how many bits to use for each field element of a given polynomial.
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
        let n = self.b0[0].coefficients.len();
        let l = n.checked_ilog2().unwrap() as u8;
        let header: u8 = (5 << 4) // fixed bits
                        | l;

        let f = &self.b0[1];
        let g = &self.b0[0];
        let capital_f = &self.b0[3];

        let mut bits = BitVec::from_bytes(&[header]);
        // f
        let width = Self::field_element_width(n, 0);
        for &fi in f.coefficients.iter() {
            let mut substring = Self::serialize_field_element(width, Felt::new(-fi));
            bits.append(&mut substring);
        }
        // g
        let width = Self::field_element_width(n, 1);
        for &fi in g.coefficients.iter() {
            let mut substring = Self::serialize_field_element(width, Felt::new(fi));
            bits.append(&mut substring);
        }
        // capital_f
        let width = Self::field_element_width(n, 2);
        for &fi in capital_f.coefficients.iter() {
            let mut substring = Self::serialize_field_element(width, Felt::new(-fi));
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
        let f = Polynomial::new(
            bit_buffer
                .iter()
                .take(n * width_f)
                .chunks(width_f)
                .into_iter()
                .map(BitVec::from_iter)
                .map(|subs| Self::deserialize_field_element(&subs))
                .collect::<Result<Vec<Felt>, _>>()?,
        );

        // g
        let width_g = Self::field_element_width(n, 1);
        let g = Polynomial::new(
            bit_buffer
                .iter()
                .skip(n * width_f)
                .take(n * width_g)
                .chunks(width_g)
                .into_iter()
                .map(BitVec::from_iter)
                .map(|subs| Self::deserialize_field_element(&subs))
                .collect::<Result<Vec<Felt>, _>>()?,
        );

        // capital_f
        let width_capital_f = Self::field_element_width(n, 2);
        let capital_f = Polynomial::new(
            bit_buffer
                .iter()
                .skip(n * width_g + n * width_f)
                .take(n * width_capital_f)
                .chunks(width_capital_f)
                .into_iter()
                .map(BitVec::from_iter)
                .map(|subs| Self::deserialize_field_element(&subs))
                .collect::<Result<Vec<Felt>, _>>()?,
        );

        // all bits in the bit buffer should have been read at this point
        if bit_buffer.len() != n * width_f + n * width_g + n * width_capital_f {
            return Err(FalconDeserializationError::BadEncodingLength);
        }

        // capital_g * f - g * capital_f = Q (mod X^n + 1)
        let capital_g = g
            .fft()
            .hadamard_div(&f.fft())
            .hadamard_mul(&capital_f.fft())
            .ifft();

        Ok(Self::from_b0([
            g.map(|f| f.balanced_value()),
            -f.map(|f| f.balanced_value()),
            capital_g.map(|f| f.balanced_value()),
            -capital_f.map(|f| f.balanced_value()),
        ]))
    }
}

impl<const N: usize> PartialEq for SecretKey<N> {
    fn eq(&self, other: &Self) -> bool {
        let own_f = &self.b0[1];
        let own_g = &self.b0[0];
        let own_capital_f = &self.b0[3];
        let own_capital_g = &self.b0[2];

        let other_f = &other.b0[1];
        let other_g = &other.b0[0];
        let other_capital_f = &other.b0[3];
        let other_capital_g = &other.b0[2];

        own_f == other_f
            && own_g == other_g
            && own_capital_f == other_capital_f
            && own_capital_g == other_capital_g
    }
}

impl<const N: usize> Eq for SecretKey<N> {}

#[derive(Debug, Clone, PartialEq)]
pub struct PublicKey<const N: usize> {
    h: Polynomial<Felt>,
}

impl<const N: usize> PublicKey<N> {
    /// Compute the public key that matches with this secret key.
    pub fn from_secret_key(sk: &SecretKey<N>) -> Self {
        let f = sk.b0[1].map(|&c| -Felt::new(c));
        let f_ntt = f.fft();
        let g = sk.b0[0].map(|&c| Felt::new(c));
        let g_ntt = g.fft();
        let h_ntt = g_ntt.hadamard_div(&f_ntt);
        let h = h_ntt.ifft();
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
        let h = Polynomial::new(
            bit_buffer
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
                .collect_vec(),
        );

        Ok(PublicKey { h })
    }

    // Serialize the public key as a list of bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let header = self.h.coefficients.len().ilog2() as u8;
        let mut bit_buffer = BitVec::from_bytes(&[header]);

        for hi in self.h.coefficients.iter() {
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
    let c_over_q_fft = c
        .map(|cc| Complex::new(one_over_q * cc.value() as f64, 0.0))
        .fft();

    // B = [[FFT(g), -FFT(f)], [FFT(G), -FFT(F)]]
    let capital_f_fft = sk.b0[3].map(|&i| Complex64::new(-i as f64, 0.0)).fft();
    let f_fft = sk.b0[1].map(|&i| Complex64::new(-i as f64, 0.0)).fft();
    let capital_g_fft = sk.b0[2].map(|&i| Complex64::new(i as f64, 0.0)).fft();
    let g_fft = sk.b0[0].map(|&i| Complex64::new(i as f64, 0.0)).fft();
    let t0 = c_over_q_fft.hadamard_mul(&capital_f_fft);
    let t1 = -c_over_q_fft.hadamard_mul(&f_fft);

    let s = loop {
        let bold_s = loop {
            let z = ffsampling(&(t0.clone(), t1.clone()), &sk.tree, &params, &mut rng);
            let t0_min_z0 = t0.clone() - z.0;
            let t1_min_z1 = t1.clone() - z.1;

            // s = (t-z) * B
            let s0 = t0_min_z0.hadamard_mul(&g_fft) + t1_min_z1.hadamard_mul(&capital_g_fft);
            let s1 = t0_min_z0.hadamard_mul(&f_fft) + t1_min_z1.hadamard_mul(&capital_f_fft);

            // compute the norm of (s0||s1) and note that they are in FFT representation
            let length_squared: f64 = (s0
                .coefficients
                .iter()
                .map(|a| (a * a.conj()).re)
                .sum::<f64>()
                + s1.coefficients
                    .iter()
                    .map(|a| (a * a.conj()).re)
                    .sum::<f64>())
                / (n as f64);

            if length_squared > (bound as f64) {
                continue;
            }

            break [s0, s1];
        };
        let s2 = bold_s[1].ifft();
        let maybe_s = compress(
            &s2.coefficients
                .iter()
                .map(|a| a.re.round() as i16)
                .collect_vec(),
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
    let s2_ntt = Polynomial::new(s2.iter().map(|a| Felt::new(*a)).collect_vec()).fft();
    let h_ntt = pk.h.fft();
    let c_ntt = c.fft();

    // s1 = c - s2 * pk.h;
    let s1_ntt = c_ntt - s2_ntt.hadamard_mul(&h_ntt);
    let s1 = s1_ntt.ifft();

    let length_squared = s1
        .coefficients
        .iter()
        .map(|i| i.balanced_value() as i64)
        .map(|i| (i * i))
        .sum::<i64>()
        + s2.iter().map(|&i| i as i64).map(|i| (i * i)).sum::<i64>();
    length_squared < params.sig_bound
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use rand::{rngs::StdRng, thread_rng, Rng, RngCore, SeedableRng};

    use crate::{
        encoding::compress,
        falcon::{keygen, sign, verify, FalconVariant, Signature},
        field::Felt,
        polynomial::{hash_to_point, Polynomial},
    };

    use super::{PublicKey, SecretKey};

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
}
