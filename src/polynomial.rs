use num_complex::Complex64;
use rand_distr::num_traits::{One, Zero};
use sha3::{digest::*, Shake256};
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use itertools::Itertools;

use crate::fft::{fft, ifft};
use crate::field::{Felt, Q};

#[derive(Debug, Clone)]
pub struct Polynomial<F> {
    pub coefficients: Vec<F>,
}

impl<F> Polynomial<F>
where
    F: Clone,
{
    pub fn new(coefficients: &[F]) -> Self {
        Self {
            coefficients: coefficients.to_vec(),
        }
    }
}

impl<F> Polynomial<F>
where
    F: Clone + Neg<Output = F>,
{
    /// Compute the Hermitian adjoint of the polynomial f in the
    /// cyclotomic ring Q[ X ] / < Phi_n(X) > where n >= deg(f)+1.
    /// In this structure, the Hermitian adjoint is given by
    ///
    ///     f*(X) = f[0] + sum_{i=1}^{n-1} f[i] * X({n-i}) .
    pub fn hermitian_adjoint(&self) -> Polynomial<F> {
        let coefficients = [
            vec![self.coefficients[0].clone()],
            self.coefficients
                .iter()
                .skip(1)
                .cloned()
                .map(|c| -c)
                .rev()
                .collect_vec(),
        ]
        .concat();
        Polynomial { coefficients }
    }
}

impl<F: Zero + PartialEq + Clone> Polynomial<F> {
    pub fn degree(&self) -> Option<usize> {
        if self.coefficients.is_empty() {
            return None;
        }
        let mut max_index = self.coefficients.len() - 1;
        while self.coefficients[max_index] == F::zero() {
            if let Some(new_index) = max_index.checked_sub(1) {
                max_index = new_index;
            } else {
                return None;
            }
        }
        Some(max_index)
    }
    pub fn lc(&self) -> F {
        match self.degree() {
            Some(non_negative_degree) => self.coefficients[non_negative_degree].clone(),
            None => F::zero(),
        }
    }
}

impl<F: Zero + Clone> Polynomial<F> {
    pub fn shift(&self, shamt: usize) -> Self {
        Self {
            coefficients: [vec![F::zero(); shamt], self.coefficients.clone()].concat(),
        }
    }

    pub fn constant(f: F) -> Self {
        Self {
            coefficients: vec![f],
        }
    }

    pub fn map<G: Clone, C: FnMut(F) -> G>(&self, closure: C) -> Polynomial<G> {
        Polynomial::<G>::new(&self.coefficients.iter().cloned().map(closure).collect_vec())
    }

    pub fn fold<G, C: FnMut(G, F) -> G + Clone>(&self, mut initial_value: G, closure: C) -> G {
        for c in self.coefficients.iter().cloned() {
            initial_value = (closure.clone())(initial_value, c);
        }
        initial_value
    }
}

/// The following implementations are specific to cyclotomic polynomial rings,
/// i.e., F[ X ] / <X^n + 1>, and are used extensively in Falcon.
impl<
        F: One
            + Zero
            + Clone
            + Neg<Output = F>
            + MulAssign
            + AddAssign
            + Div<Output = F>
            + Sub<Output = F>
            + PartialEq,
    > Polynomial<F>
{
    /// Reduce the polynomial by X^n + 1.
    pub fn reduce_by_cyclotomic(&self, n: usize) -> Self {
        let mut coefficients = vec![F::zero(); n];
        let mut sign = -F::one();
        for (i, c) in self.coefficients.iter().cloned().enumerate() {
            if i % n == 0 {
                sign *= -F::one();
            }
            coefficients[i % n] += sign.clone() * c;
        }
        Polynomial::new(&coefficients)
    }

    /// Compute the multiplicative inverse of the polynomial in the ring
    /// F[ X ] / <X^n + 1>
    ///
    /// This function assumes that F is a field; otherwise the gcd will never end.
    pub fn cyclotomic_ring_inverse(&self, n: usize) -> Self {
        let mut cyclotomic_coefficients = vec![F::zero(); n + 1];
        cyclotomic_coefficients[0] = F::one();
        cyclotomic_coefficients[n] = F::one();
        let (_, a, _) = Polynomial::xgcd(self, &Polynomial::new(&cyclotomic_coefficients));
        a
    }

    /// Compute the field norm of the polynomial as an element of the cyclotomic
    /// ring  F[ X ] / <X^n + 1 > relative to one of half the size, i.e.,
    ///  F[ X ] / <X^(n/2) + 1> .
    ///
    /// Corresponds to formula 3.25 in the spec [1, p.30].
    ///
    /// [1]: https://falcon-sign.info/falcon.pdf
    pub fn field_norm(&self) -> Self {
        let n = self.coefficients.len();
        let mut f0_coefficients = vec![F::zero(); n / 2];
        let mut f1_coefficients = vec![F::zero(); n / 2];
        for i in 0..n / 2 {
            f0_coefficients[i] = self.coefficients[2 * i].clone();
            f1_coefficients[i] = self.coefficients[2 * i + 1].clone();
        }
        let f0 = Polynomial::new(&f0_coefficients);
        let f1 = Polynomial::new(&f1_coefficients);
        let f0_squared = (f0.clone() * f0).reduce_by_cyclotomic(n / 2);
        let f1_squared = (f1.clone() * f1).reduce_by_cyclotomic(n / 2);
        let x = Polynomial::new(&[F::zero(), F::one()]);
        f0_squared - (x * f1_squared).reduce_by_cyclotomic(n / 2)
    }

    /// Lift an element from a cyclotomic polynomial ring to one of double the
    /// size.
    pub fn lift_next_cyclotomic(&self) -> Self {
        let n = self.coefficients.len();
        let mut coefficients = vec![F::zero(); n * 2];
        for i in 0..n {
            coefficients[2 * i] = self.coefficients[i].clone();
        }
        Self::new(&coefficients)
    }

    /// Compute the galois adjoint of the polynomial in the cyclotomic ring
    /// F[ X ] / < X^n + 1 > , which corresponds to f(x^2).
    pub fn galois_adjoint(&self) -> Self {
        Self::new(
            &self
                .coefficients
                .iter()
                .cloned()
                .enumerate()
                .map(|(i, c)| if i % 2 == 0 { c } else { -c })
                .collect_vec(),
        )
    }
}

impl<
        F: One
            + Zero
            + Clone
            + Neg<Output = F>
            + MulAssign
            + AddAssign
            + Div<Output = F>
            + Sub<Output = F>
            + PartialEq,
    > Polynomial<F>
{
    /// Extended Euclidean algorithm for polynomials. Uses the EEA to compute
    /// the greatest common divisor g and Bezout coefficients u, v such that
    ///     u * a + v * b = 1
    pub fn xgcd(a: &Self, b: &Self) -> (Self, Self, Self) {
        if a.is_zero() || b.is_zero() {
            return (Self::zero(), Self::zero(), Self::zero());
        }
        let (mut old_r, mut r) = (a.clone(), b.clone());
        let (mut old_s, mut s) = (Self::one(), Self::zero());
        let (mut old_t, mut t) = (Self::zero(), Self::one());

        while !r.is_zero() {
            let quotient = old_r.clone() / r.clone();
            (old_r, r) = (r.clone(), old_r.clone() - quotient.clone() * r.clone());
            (old_s, s) = (s.clone(), old_s.clone() - quotient.clone() * s.clone());
            (old_t, t) = (t.clone(), old_t.clone() - quotient.clone() * t.clone());
        }

        (old_r, old_s, old_t)
    }
}

impl Polynomial<f64> {
    /// Compute a floating-point approximation of the multiplicative inverse of the
    /// polynomial in the ring  F[ X ] / <X^n + 1>
    pub fn approximate_cyclotomic_ring_inverse(&self, n: usize) -> Self {
        let coefficients = [
            self.coefficients.clone(),
            vec![0.0; n - self.coefficients.len()],
        ]
        .concat()
        .iter()
        .map(|r| Complex64::new(*r, 0.0))
        .collect_vec();
        let fft_coefficients = fft(&coefficients);
        let inverse_coefficients =
            ifft(&fft_coefficients.into_iter().map(|c| 1.0 / c).collect_vec());
        Self {
            coefficients: inverse_coefficients.iter().map(|c| c.re).collect_vec(),
        }
    }
    pub fn approximate_cyclotomic_ring_divide(&self, other: &Self, n: usize) -> Self {
        let self_coefficients = [
            self.coefficients.clone(),
            vec![0.0; n - self.coefficients.len()],
        ]
        .concat()
        .iter()
        .map(|r| Complex64::new(*r, 0.0))
        .collect_vec();
        let other_coefficients = [
            other.coefficients.clone(),
            vec![0.0; n - other.coefficients.len()],
        ]
        .concat()
        .iter()
        .map(|r| Complex64::new(*r, 0.0))
        .collect_vec();
        let self_fft = fft(&self_coefficients);
        let other_fft = fft(&other_coefficients);
        let quotient_fft = self_fft
            .into_iter()
            .zip(other_fft)
            .map(|(a, b)| a / b)
            .collect_vec();
        let quotient_coefficients = ifft(&quotient_fft);
        Polynomial::new(
            &quotient_coefficients
                .into_iter()
                .map(|c| c.re)
                .collect_vec(),
        )
    }
}

impl<F: Clone + Into<f64>> Polynomial<F> {
    pub fn l2_norm(&self) -> f64 {
        self.coefficients
            .iter()
            .map(|i| Into::<f64>::into(i.clone()))
            .map(|i| i * i)
            .sum::<f64>()
            .sqrt()
    }
}

impl<F> PartialEq for Polynomial<F>
where
    F: Zero
        + One
        + PartialEq
        + Neg<Output = F>
        + Clone
        + AddAssign
        + MulAssign
        + Div<Output = F>
        + Sub<Output = F>,
{
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() && other.is_zero() {
            true
        } else if self.is_zero() || other.is_zero() {
            false
        } else {
            let self_degree = self.degree().unwrap();
            let other_degree = other.degree().unwrap();
            self.coefficients[0..=self_degree] == other.coefficients[0..=other_degree]
        }
    }
}

impl<F> Eq for Polynomial<F> where
    F: Zero
        + One
        + PartialEq
        + Neg<Output = F>
        + Clone
        + AddAssign
        + MulAssign
        + Div<Output = F>
        + Sub<Output = F>
{
}

impl<F> Add for Polynomial<F>
where
    F: Add<Output = F> + AddAssign + Clone,
{
    type Output = Polynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        let coefficients = if self.coefficients.len() >= rhs.coefficients.len() {
            let mut coefficients = self.coefficients.clone();
            for (i, c) in rhs.coefficients.into_iter().enumerate() {
                coefficients[i] += c;
            }
            coefficients
        } else {
            let mut coefficients = rhs.coefficients.clone();
            for (i, c) in self.coefficients.into_iter().enumerate() {
                coefficients[i] += c;
            }
            coefficients
        };
        Self::Output { coefficients }
    }
}

impl<F> AddAssign for Polynomial<F>
where
    F: Add<Output = F> + AddAssign + Clone,
{
    fn add_assign(&mut self, rhs: Self) {
        if self.coefficients.len() >= rhs.coefficients.len() {
            for (i, c) in rhs.coefficients.into_iter().enumerate() {
                self.coefficients[i] += c;
            }
        } else {
            let mut coefficients = rhs.coefficients.clone();
            for (i, c) in self.coefficients.iter().enumerate() {
                coefficients[i] += c.clone();
            }
            self.coefficients = coefficients;
        }
    }
}

impl<F> Sub for Polynomial<F>
where
    F: Sub<Output = F> + Clone + Neg<Output = F> + Add<Output = F> + AddAssign,
{
    type Output = Polynomial<F>;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<F> SubAssign for Polynomial<F>
where
    F: Add<Output = F> + Neg<Output = F> + AddAssign + Clone + Sub<Output = F>,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.coefficients = self.clone().sub(rhs).coefficients;
    }
}

impl<F: Neg<Output = F> + Clone> Neg for Polynomial<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Output {
            coefficients: self.coefficients.iter().cloned().map(|a| -a).collect(),
        }
    }
}

impl<F> Mul for Polynomial<F>
where
    F: Add
        + AddAssign
        + Mul<Output = F>
        + MulAssign
        + Div<Output = F>
        + Neg<Output = F>
        + Sub<Output = F>
        + Zero
        + One
        + PartialEq
        + Clone,
{
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut coefficients =
            vec![F::zero(); self.coefficients.len() + other.coefficients.len() - 1];
        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                coefficients[i + j] += self.coefficients[i].clone() * other.coefficients[j].clone();
            }
        }
        Self { coefficients }
    }
}

impl<F: Add + Mul<Output = F> + Zero + Clone> Mul<F> for Polynomial<F> {
    type Output = Self;

    fn mul(self, other: F) -> Self::Output {
        Self {
            coefficients: self
                .coefficients
                .iter()
                .cloned()
                .map(|i| i * other.clone())
                .collect_vec(),
        }
    }
}
impl<F> One for Polynomial<F>
where
    F: Zero
        + One
        + Neg<Output = F>
        + PartialEq
        + AddAssign
        + Clone
        + Add<Output = F>
        + Sub<Output = F>
        + Mul<Output = F>
        + Div<Output = F>
        + MulAssign,
{
    fn one() -> Self {
        Self {
            coefficients: vec![F::one()],
        }
    }
}

impl<
        F: Zero
            + One
            + Neg<Output = F>
            + PartialEq
            + AddAssign
            + Clone
            + Add<Output = F>
            + Mul<Output = F>
            + MulAssign
            + Div<Output = F>
            + Sub<Output = F>,
    > Zero for Polynomial<F>
{
    fn zero() -> Self {
        Self {
            coefficients: vec![],
        }
    }

    fn is_zero(&self) -> bool {
        self.degree().is_none()
    }
}

impl<F> Div<Polynomial<F>> for Polynomial<F>
where
    F: Zero
        + One
        + PartialEq
        + AddAssign
        + Clone
        + Mul<Output = F>
        + MulAssign
        + Div<Output = F>
        + Neg<Output = F>
        + Sub<Output = F>,
{
    type Output = Polynomial<F>;

    fn div(self, denominator: Self) -> Self::Output {
        if denominator.is_zero() {
            panic!();
        }
        if self.is_zero() {
            Self::zero();
        }
        let mut remainder = self.clone();
        let mut quotient = Polynomial::<F>::zero();
        while remainder.degree().unwrap() >= denominator.degree().unwrap() {
            let shift = remainder.degree().unwrap() - denominator.degree().unwrap();
            let quotient_coefficient = remainder.lc() / denominator.lc();
            let monomial = Self::constant(quotient_coefficient).shift(shift);
            quotient += monomial.clone();
            remainder -= monomial * denominator.clone();
            if remainder.is_zero() {
                break;
            }
        }
        quotient
    }
}

/// Hash a string to a random polynomial in ZZ[ X ] mod <Phi(X), q>.
/// Algorithm 3, "HashToPoint" in the spec (page 31).
pub fn hash_to_point(string: &[u8], n: usize) -> Polynomial<Felt> {
    const K: u32 = (1u32 << 16) / (Q as u32);

    let mut hasher = Shake256::default();
    hasher.update(string);
    let mut reader = hasher.finalize_xof();

    let mut coefficients: Vec<Felt> = vec![];
    while coefficients.len() != n {
        let mut randomness = [0u8; 2];
        reader.read(&mut randomness);
        // Arabic endianness but so be it
        let t = ((randomness[0] as i32) << 8) | (randomness[1] as i32);
        if (t as u32) < K * (Q as u32) {
            coefficients.push(Felt::new((t % Q) as i16));
        }
    }

    Polynomial { coefficients }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num_complex::Complex64;
    use rand::thread_rng;
    use rand::Rng;
    use rand_distr::num_traits::One;
    use sha3::digest::ExtendableOutput;
    use sha3::digest::Update;
    use sha3::digest::XofReader;
    use sha3::Shake256;

    use crate::fft::fft;
    use crate::fft::ifft;
    use crate::field::Felt;
    use crate::polynomial::hash_to_point;
    use crate::polynomial::Polynomial;

    #[test]
    fn test_shake256() {
        let input = b"\x21\xF1\x34\xAC\x57";
        let kat_output : [u8;64] = *b"\xBB\x8A\x84\x47\x51\x7B\xA9\xCA\x7F\xA3\x4E\xC9\x9A\x80\x00\x4F\x22\x8A\xB2\x82\x47\x28\x41\xEB\x3D\x3A\x76\x22\x5C\x9D\xBE\x77\xF7\xE4\x0A\x06\x67\x76\xD3\x2C\x74\x94\x12\x02\xF9\xF4\xAA\x43\xD1\x2C\x62\x64\xAF\xA5\x96\x39\xC4\x4E\x11\xF5\xE1\x4F\x1E\x56";
        let mut observed_output = [0u8; 64];

        let mut hasher = Shake256::default();
        hasher.update(input);
        let mut reader = hasher.finalize_xof();
        reader.read(&mut observed_output);

        assert_eq!(
            kat_output, observed_output,
            "SHAKE256 from crate sha3 does not match KAT"
        )
    }

    #[test]
    fn test_hash_to_point_sanity() {
        // KAT obtained from python script, run in the same directory
        // as a copy of the official python implementation.
        // ```
        // from falcon import SecretKey
        // sk = SecretKey(512)
        // sk.hash_to_point(b"", b"")
        // ```
        let hash_of_empty = Polynomial {
            coefficients: [
                5816, 7463, 2984, 11537, 9019, 4074, 5180, 11040, 4044, 8937, 694, 7042, 9481,
                10084, 3795, 5677, 5977, 1241, 6332, 2817, 413, 1971, 755, 7241, 6041, 9347, 4136,
                11948, 9140, 1210, 5150, 1630, 4015, 2390, 2346, 2025, 4272, 10978, 7171, 8764,
                11920, 888, 12160, 11275, 7043, 10323, 1181, 1873, 5876, 12213, 627, 11319, 5675,
                8207, 6210, 385, 333, 4581, 1359, 10859, 3346, 3970, 8720, 3640, 8157, 1080, 2794,
                5769, 11618, 6780, 1734, 6484, 1575, 9433, 10353, 2004, 5921, 5013, 4753, 9865,
                10931, 6621, 1417, 9804, 12027, 9437, 10657, 3260, 9541, 4967, 12124, 6827, 333,
                6404, 3498, 6920, 3979, 14, 440, 1293, 8011, 7567, 3899, 3252, 4023, 10727, 11938,
                957, 2412, 9552, 10409, 8063, 9131, 9835, 10305, 3124, 6303, 12241, 6354, 2540,
                10113, 10777, 6803, 4879, 10952, 10503, 1728, 5067, 3339, 7045, 11333, 5469, 11062,
                11666, 5235, 2314, 3345, 2224, 2274, 8060, 4304, 6716, 11595, 1541, 996, 6983, 36,
                449, 7401, 4987, 9177, 810, 1908, 8650, 7646, 6893, 4919, 1971, 4930, 11763, 201,
                12223, 9234, 4081, 6199, 12047, 7646, 9639, 6814, 6739, 5279, 4012, 2101, 10707,
                4241, 12146, 3779, 3999, 3176, 1699, 10294, 5168, 5590, 457, 9709, 6450, 442, 8884,
                6755, 10995, 10923, 3935, 8499, 3508, 12088, 1115, 11336, 1379, 7557, 4954, 7639,
                2514, 8672, 6686, 98, 5676, 8800, 5712, 4724, 7724, 3202, 12128, 10940, 10177,
                9421, 11013, 7372, 8546, 441, 6261, 8779, 2453, 12082, 7922, 5307, 6920, 7726, 823,
                10561, 1251, 10358, 8383, 10905, 8145, 1733, 1718, 3105, 10756, 6798, 10209, 7976,
                11148, 9353, 4746, 1089, 11444, 6571, 409, 8381, 10325, 7649, 10042, 5587, 3625,
                10182, 10494, 228, 4687, 5949, 7995, 12092, 3312, 5339, 5920, 8145, 6796, 1992,
                3205, 2761, 12199, 11342, 9695, 390, 252, 989, 1385, 12148, 8324, 10694, 3690,
                3440, 8888, 12238, 9018, 3354, 5859, 6298, 8098, 4388, 3788, 3045, 11095, 2372,
                10036, 9233, 168, 8500, 3604, 2494, 9854, 5679, 2182, 3350, 7798, 8310, 3544, 948,
                7646, 7235, 2650, 6008, 4610, 2159, 6884, 10545, 688, 4115, 10312, 4408, 4951,
                2891, 9791, 1377, 5909, 11147, 11139, 4969, 5158, 350, 1067, 4242, 10820, 1818,
                6473, 105, 2919, 10892, 7116, 850, 11409, 2652, 6392, 2540, 6892, 8372, 11975,
                4994, 2621, 2763, 11837, 6132, 11293, 9138, 8769, 10964, 9826, 601, 7007, 9078,
                10340, 9410, 8746, 10835, 9053, 11010, 5308, 8851, 1976, 11016, 599, 8348, 9876,
                7100, 1333, 4550, 1735, 4598, 9970, 525, 8320, 1609, 9213, 4178, 484, 10814, 1760,
                9667, 8369, 2286, 10384, 12139, 24, 1178, 5682, 7074, 3676, 3661, 3322, 1831, 5562,
                734, 8059, 8750, 6951, 4760, 10395, 1019, 9404, 2923, 6715, 123, 10157, 4892, 7667,
                1677, 4175, 3455, 12123, 10730, 2000, 8212, 2665, 7088, 8741, 10936, 3172, 225,
                3867, 5140, 2310, 6453, 2898, 3637, 4580, 113, 5991, 3532, 3363, 11457, 11601,
                7280, 6792, 11872, 8127, 2192, 10761, 9019, 8197, 8965, 6061, 10799, 988, 10522,
                1281, 1965, 2716, 9841, 7496, 8456, 5192, 825, 3727, 4664, 7906, 8521, 5901, 10200,
                5167, 9451, 10825, 12011, 2272, 8698, 8174, 11973, 5155, 6890, 9999, 4391, 12044,
                1620, 8310, 111, 4481, 9650, 2077, 7691, 7531, 1956, 494, 3297, 1623, 3266, 7018,
                2031, 6317, 4657, 5206, 2581, 11227, 10508, 4567, 8892, 1363, 6790, 6180, 1588,
                9776, 11998, 10689, 492, 331,
            ]
            .map(Felt::new)
            .to_vec(),
        };
        assert_eq!(hash_of_empty, hash_to_point(&[], 512));
    }

    #[test]
    fn test_div() {
        let mut rng = thread_rng();
        let n = rng.gen_range(1..100);
        let m = rng.gen_range(1..100);
        let expected_coefficients = (0..n).map(|_| rng.gen_range(-5..5)).collect_vec();
        let cofactor_coefficients = (0..m).map(|_| rng.gen_range(-5..5)).collect_vec();
        let cofactor_polynomial = Polynomial::new(&cofactor_coefficients);
        let product = Polynomial::new(&expected_coefficients) * cofactor_polynomial.clone();
        let quotient = product / cofactor_polynomial;
        assert_eq!(Polynomial::new(&expected_coefficients), quotient);
        println!("quotient: {:?}", quotient);
    }

    #[test]
    fn test_approximate_cyclotomic_ring_inverse() {
        let mut rng = thread_rng();
        let n = 64;
        let coefficients = (0..n).map(|_| rng.gen_range(-5..5) as f64).collect_vec();
        let polynomial = Polynomial::new(&coefficients);
        let inverse = polynomial.approximate_cyclotomic_ring_inverse(n);
        let product = polynomial * inverse;
        let difference = product.reduce_by_cyclotomic(n) - Polynomial::<f64>::one();
        assert!(difference.l2_norm() < f64::EPSILON * 100.0);
    }

    #[test]
    fn test_cyclotomic_multiplication() {
        let mut rng = thread_rng();
        let n = 64;
        let coefficients_a = (0..n).map(|_| rng.gen_range(-5..5) as f64).collect_vec();
        let coefficients_b = (0..n).map(|_| rng.gen_range(-5..5) as f64).collect_vec();
        let polynomial_a = Polynomial::new(&coefficients_a);
        let polynomial_b = Polynomial::new(&coefficients_b);
        let reduced = (polynomial_a * polynomial_b).reduce_by_cyclotomic(n);

        let complex_coefficients_a = coefficients_a
            .into_iter()
            .map(|r| Complex64::new(r, 0.0))
            .collect_vec();
        let complex_coefficients_b = coefficients_b
            .into_iter()
            .map(|r| Complex64::new(r, 0.0))
            .collect_vec();
        let fft_a = fft(&complex_coefficients_a);
        let fft_b = fft(&complex_coefficients_b);
        let had = fft_a
            .into_iter()
            .zip(fft_b)
            .map(|(a, b)| a * b)
            .collect_vec();
        let product = ifft(&had);
        let fft_polynomial = Polynomial::new(&product.iter().map(|&c| c.re).collect_vec());

        let difference = reduced.clone() - fft_polynomial.clone();
        assert!(
            difference.l2_norm() <= f64::EPSILON * 10000.0,
            "reduced: {:?}\nfft: {:?}",
            reduced,
            fft_polynomial
        );
    }

    #[test]
    fn test_adjoint() {
        let mut rng = thread_rng();
        let n = 64;
        let coefficients = (0..n).map(|_| rng.gen_range(-5..5) as f64).collect_vec();
        let polynomial = Polynomial::new(&coefficients);
        let hermitian_adjoint = polynomial.hermitian_adjoint();
        let coefficients_fft = fft(&coefficients
            .into_iter()
            .map(|c| Complex64::new(c, 0.0))
            .collect_vec());
        let coefficients_conjugate =
            ifft(&coefficients_fft.into_iter().map(|c| c.conj()).collect_vec());
        let conjugate_polynomial = Polynomial::new(
            &coefficients_conjugate
                .into_iter()
                .map(|c| c.re)
                .collect_vec(),
        );

        let difference = hermitian_adjoint.clone() - conjugate_polynomial.clone();
        assert!(
            difference.l2_norm() <= f64::EPSILON * 100.0,
            "hermitian_adjoint: {:?}\nconjugate: {:?}",
            hermitian_adjoint.coefficients,
            conjugate_polynomial.coefficients
        );
    }

    #[test]
    fn test_approximate_cyclotomic_ring_divide() {
        let mut rng = thread_rng();
        let n = 64;
        let coefficients_lhs = (0..n).map(|_| rng.gen_range(-5..5) as f64).collect_vec();
        let lhs = Polynomial::new(&coefficients_lhs);
        let coefficients_rhs = (0..n).map(|_| rng.gen_range(-5..5) as f64).collect_vec();
        let rhs = Polynomial::new(&coefficients_rhs);

        let inverse = rhs.approximate_cyclotomic_ring_inverse(n);
        let product1 = (lhs.clone() * inverse).reduce_by_cyclotomic(n);

        let product2 = lhs.approximate_cyclotomic_ring_divide(&rhs, n);

        let difference = product1.clone() - product2.clone();
        assert!(
            difference.l2_norm() <= f64::EPSILON * 100.0,
            "lhs*inverse: {:?}\nlhs/rhs: {:?}",
            product1.coefficients,
            product2.coefficients
        );
    }
}
