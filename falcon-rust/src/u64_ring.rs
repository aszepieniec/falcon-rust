use std::fmt::Display;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use rand_distr::{
    num_traits::{One, Zero},
    Distribution, Standard,
};

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;

/// Square of [`u32_field::Q`].
const Q2: u64 = 1152947895184416769_u64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct U64Ring(pub(crate) u64);

impl U64Ring {
    pub const fn new(value: i64) -> Self {
        let gtz_bool = value >= 0;
        let gtz_int = gtz_bool as i64;
        let gtz_sign = gtz_int - ((!gtz_bool) as i64);
        let reduced = gtz_sign * ((gtz_sign * value) % (Q2 as i64));
        let canonical_representative = (reduced + (Q2 as i64) * (1 - gtz_int)) as u64;
        U64Ring(canonical_representative)
    }

    pub const fn value(&self) -> i64 {
        self.0 as i64
    }
}

impl From<usize> for U64Ring {
    fn from(value: usize) -> Self {
        U64Ring::new(value as i64)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for U64Ring {
    fn add(self, rhs: Self) -> Self::Output {
        let (s, _) = self.0.overflowing_add(rhs.0);
        let (d, n) = s.overflowing_sub(Q2);
        let (r, _) = d.overflowing_add(Q2 * (n as u64));
        U64Ring(r)
    }

    type Output = Self;
}

impl AddAssign for U64Ring {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for U64Ring {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl SubAssign for U64Ring {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for U64Ring {
    type Output = U64Ring;

    fn neg(self) -> Self::Output {
        let is_nonzero = self.0 != 0;
        let r = Q2 - self.0;
        U64Ring(r * (is_nonzero as u64))
    }
}

impl Mul for U64Ring {
    fn mul(self, rhs: Self) -> Self::Output {
        U64Ring((((self.0 as u128) * (rhs.0 as u128)) % (Q2 as u128)) as u64)
    }

    type Output = Self;
}

impl MulAssign for U64Ring {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Zero for U64Ring {
    fn zero() -> Self {
        U64Ring::new(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}
impl One for U64Ring {
    fn one() -> Self {
        U64Ring::new(1)
    }
}

impl Distribution<U64Ring> for Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> U64Ring {
        U64Ring::new(((rng.next_u64() >> 1) % Q2) as i64)
    }
}

impl Display for U64Ring {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("{}", self.value()))
    }
}

fn xgcd(a: i64, b: i64) -> (i64, i64, i64) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);

    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }

    (old_s, old_t, old_r)
}

impl Inverse for U64Ring {
    fn inverse_or_zero(self) -> Self {
        let (a, _b, _g) = xgcd(self.0.try_into().unwrap(), Q2.try_into().unwrap());
        Self::new(a)
    }
}

impl CyclotomicFourier for U64Ring {
    fn primitive_root_of_unity(n: usize) -> Self {
        let log2n = n.ilog2();
        assert!(log2n <= 12);
        // and 765470292790715721_i64 is a twelfth root of unity
        let mut a = U64Ring::new(765470292790715721_i64);
        let num_squarings = 12 - n.ilog2();
        for _ in 0..num_squarings {
            a *= a;
        }
        a
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num::One;
    use rand::{thread_rng, Rng, RngCore};

    use crate::{
        cyclotomic_fourier::CyclotomicFourier,
        inverse::Inverse,
        polynomial::Polynomial,
        u64_ring::{U64Ring, Q2},
    };
    use num::Zero;

    #[test]
    fn test_value() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let mut value = (rng.next_u64() & 0x3fffffffffff) as i64;
            if rng.next_u32() % 2 == 1 {
                value *= -1;
            }
            let uf = U64Ring::new(value);
            assert_eq!(
                0,
                (uf.value() - value) % (Q2 as i64),
                "value: {value} but got {}",
                uf.value()
            );
        }
    }

    #[test]
    fn test_add() {
        let mut rng = thread_rng();
        let a_value = (rng.next_u64() % 0x0fffffffffff) as i64;
        let b_value = (rng.next_u64() % 0x0fffffffffff) as i64;
        let a = U64Ring::new(a_value);
        let b = U64Ring::new(b_value);
        assert_eq!(
            a + b,
            U64Ring::new(a.value() + b.value()),
            "a: {a_value}, b: {b_value}, c: {}",
            ((a_value + b_value) as u64) % Q2
        );
    }

    #[test]
    fn test_mul() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let a_value = (rng.next_u64() % 0x3fffffffffff) as i64;
            let b_value = (rng.next_u64() % 0x3fffffffffff) as i64;
            let product = (((a_value as u128) * (b_value as u128)) % (Q2 as u128)) as i64;
            let a = U64Ring::new(a_value);
            let b = U64Ring::new(b_value);
            assert_eq!(
                a * b,
                U64Ring::new(product),
                "{} =/= {}",
                a * b,
                U64Ring::new(product)
            );
        }
    }

    #[test]
    fn test_batch_inverse() {
        let mut rng = thread_rng();
        let a: [U64Ring; 64] = (0..64).map(|_| rng.gen()).collect_vec().try_into().unwrap();
        let b_batch = U64Ring::batch_inverse_or_zero(&a);
        let b_regular = a.iter().map(|e| e.inverse_or_zero()).collect_vec();
        assert_eq!(b_batch.to_vec(), b_regular);
    }

    #[test]
    fn test_inverse() {
        let mut rng = thread_rng();
        let a: U64Ring = rng.gen();
        let b = a.inverse_or_zero();

        assert_eq!(a * b * a, a);
        assert_eq!(a * b * b, b);
    }

    #[test]
    fn test_primitive_nth_root_of_unity() {
        for log2n in 0..=12 {
            let n = 1 << log2n;
            let mut root = U64Ring::primitive_root_of_unity(n);
            for i in 0..log2n {
                assert_ne!(root, U64Ring::one(), "log2n: {log2n} and i: {i}");
                root *= root;
            }
            assert_eq!(root, U64Ring::one());
        }
    }

    #[test]
    fn test_bitreverse() {
        let test_vectors = [
            vec![(0, 0), (1, 1)],
            vec![(2, 1), (3, 3), (0, 0)],
            vec![(4, 1), (5, 5), (6, 3)],
        ];
        for (i, vector) in test_vectors.into_iter().enumerate() {
            let n = 1 << (i + 1);
            for (a, b) in vector.into_iter() {
                assert_eq!(U64Ring::bitreverse_index(a, n), b);
                assert_eq!(U64Ring::bitreverse_index(b, n), a);
            }
        }
    }

    #[test]
    fn test_ntt() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| rng.next_u64() as i64)
            .map(U64Ring::new)
            .collect_vec();
        let mut b = a.clone();
        assert_eq!(a, b);

        let psi_rev = U64Ring::bitreversed_powers(n);
        let psi_inv_rev = U64Ring::bitreversed_powers_inverse(n);
        let ninv = U64Ring::inverse_or_zero(U64Ring::new(n as i64));
        U64Ring::fft(&mut a, &psi_rev);
        U64Ring::ifft(&mut a, &psi_inv_rev, ninv);
        assert_eq!(a, b);

        let x = U64Ring::new(rng.next_u32() as i64);
        let y = U64Ring::new(rng.next_u32() as i64);
        let mut c = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        U64Ring::fft(&mut a, &psi_rev);
        U64Ring::fft(&mut b, &psi_rev);
        U64Ring::fft(&mut c, &psi_rev);

        let c_alt = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        assert_eq!(c, c_alt);
    }

    #[test]
    fn test_multiply_reduce() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| U64Ring::new(rng.gen_range(-20..20)))
            .collect_vec();
        let mut b = (0..n)
            .map(|_| U64Ring::new(rng.gen_range(-20..20)))
            .collect_vec();

        let c = (Polynomial::new(a.clone()) * Polynomial::new(b.clone()))
            .reduce_by_cyclotomic(n)
            .coefficients;

        let psi_rev = U64Ring::bitreversed_powers(n);
        U64Ring::fft(&mut a, &psi_rev);
        U64Ring::fft(&mut b, &psi_rev);
        let mut d = a.iter().zip(b.iter()).map(|(&l, &r)| l * r).collect_vec();
        let psi_inv_rev = U64Ring::bitreversed_powers_inverse(n);
        let ninv = U64Ring::new(n as i64).inverse_or_zero();
        U64Ring::ifft(&mut d, &psi_inv_rev, ninv);

        let diff = |u: &[U64Ring], v: &[U64Ring]| {
            u.iter().zip(v.iter()).map(|(&l, &r)| l - r).collect_vec()
        };

        assert_eq!(diff(&c, &d), vec![U64Ring::zero(); n]);
    }
}
