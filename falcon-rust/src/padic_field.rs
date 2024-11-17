use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num::traits::ConstZero;
use rand_distr::{
    num_traits::{One, Zero},
    Distribution, Standard,
};

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;

// 557057 = 136 * 2^12 + 1
pub(crate) const P: u32 = 557057;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PadicField(pub(crate) u32);

impl PadicField {
    pub const fn new(value: i32) -> Self {
        let gtz_bool = value >= 0;
        let gtz_int = gtz_bool as i32;
        let gtz_sign = gtz_int - ((!gtz_bool) as i32);
        let reduced = gtz_sign * ((gtz_sign * value) % (P as i32));
        let canonical_representative = (reduced + (P as i32) * (1 - gtz_int)) as u32;
        PadicField(canonical_representative)
    }

    pub const fn value(&self) -> i32 {
        self.0 as i32
    }

    pub fn balanced_value(&self) -> i32 {
        let value = self.value();
        let g = (value > ((P as i32) / 2)) as i32;
        value - (P as i32) * g
    }

    pub const fn multiply(&self, other: Self) -> Self {
        PadicField((((self.0 as u64) * (other.0 as u64)) % (P as u64)) as u32)
    }
}

impl ConstZero for PadicField {
    const ZERO: Self = Self::new(0);
}

impl From<usize> for PadicField {
    fn from(value: usize) -> Self {
        PadicField::new(value as i32)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for PadicField {
    fn add(self, rhs: Self) -> Self::Output {
        let (s, _) = self.0.overflowing_add(rhs.0);
        let (d, n) = s.overflowing_sub(P);
        let (r, _) = d.overflowing_add(P * (n as u32));
        PadicField(r)
    }

    type Output = Self;
}

impl AddAssign for PadicField {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for PadicField {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl SubAssign for PadicField {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for PadicField {
    type Output = PadicField;

    fn neg(self) -> Self::Output {
        let is_nonzero = self.0 != 0;
        let r = P - self.0;
        PadicField(r * (is_nonzero as u32))
    }
}

impl Mul for PadicField {
    fn mul(self, rhs: Self) -> Self::Output {
        PadicField((((self.0 as u64) * (rhs.0 as u64)) % (P as u64)) as u32)
    }

    type Output = Self;
}

impl MulAssign for PadicField {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Zero for PadicField {
    fn zero() -> Self {
        PadicField::new(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}
impl One for PadicField {
    fn one() -> Self {
        PadicField::new(1)
    }
}

impl Distribution<PadicField> for Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> PadicField {
        PadicField::new(((rng.next_u32() >> 1) % P) as i32)
    }
}

impl Display for PadicField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("{}", self.value()))
    }
}

impl Inverse for PadicField {
    fn inverse_or_zero(self) -> Self {
        let q_minus_two = P - 2;
        let mut acc = PadicField(1);
        let mut mask = u32::MAX - (u32::MAX >> 1);
        for _ in 0..32 {
            acc = acc * acc;
            if mask & q_minus_two != 0 {
                acc *= self;
            }
            mask >>= 1;
        }
        acc
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for PadicField {
    type Output = PadicField;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            panic!("Cannot divide by zero");
        } else {
            self.multiply(rhs.inverse_or_zero())
        }
    }
}

impl CyclotomicFourier for PadicField {
    fn primitive_root_of_unity(n: usize) -> Self {
        let log2n = n.ilog2();
        assert!(log2n <= 12);
        // and 310570 is a (2^12)th root of unity
        let mut a = PadicField::new(310570);
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
        padic_field::{PadicField, P},
        polynomial::Polynomial,
    };
    use num::Zero;

    #[test]
    fn test_value() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let mut value = (rng.next_u32() & 0x3fff) as i32;
            if rng.next_u32() % 2 == 1 {
                value *= -1;
            }
            let uf = PadicField::new(value);
            assert_eq!(
                0,
                (uf.value() - value) % (P as i32),
                "value: {value} but got {}",
                uf.value()
            );
        }
    }

    #[test]
    fn test_add() {
        let mut rng = thread_rng();
        let a_value = (rng.next_u32() % 0x0fff) as i32;
        let b_value = (rng.next_u32() % 0x0fff) as i32;
        let a = PadicField::new(a_value);
        let b = PadicField::new(b_value);
        assert_eq!(
            a + b,
            PadicField::new(a.value() + b.value()),
            "a: {a_value}, b: {b_value}, c: {}",
            ((a_value + b_value) as u32) % P
        );
    }

    #[test]
    fn test_mul() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let a_value = (rng.next_u32() % 0x3fff) as i32;
            let b_value = (rng.next_u32() % 0x3fff) as i32;
            let product = (((a_value as u32) * (b_value as u32)) % P) as i32;
            let a = PadicField::new(a_value);
            let b = PadicField::new(b_value);
            assert_eq!(
                a * b,
                PadicField::new(product),
                "{} =/= {}",
                a * b,
                PadicField::new(product)
            );
        }
    }

    #[test]
    fn test_batch_inverse() {
        let mut rng = thread_rng();
        let a: [PadicField; 64] = (0..64).map(|_| rng.gen()).collect_vec().try_into().unwrap();
        let b_batch = PadicField::batch_inverse_or_zero(&a);
        let b_regular = a.iter().map(|e| e.inverse_or_zero()).collect_vec();
        assert_eq!(b_batch.to_vec(), b_regular);
    }

    #[test]
    fn test_inverse() {
        let mut rng = thread_rng();
        let a: PadicField = rng.gen();
        let b = a.inverse_or_zero();

        assert_eq!(a * b * a, a);
        assert_eq!(a * b * b, b);
    }

    #[test]
    fn test_primitive_nth_root_of_unity() {
        for log2n in 0..=12 {
            let n = 1 << log2n;
            let mut root = PadicField::primitive_root_of_unity(n);
            for i in 0..log2n {
                assert_ne!(root, PadicField::one(), "log2n: {log2n} and i: {i}");
                root *= root;
            }
            assert_eq!(root, PadicField::one());
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
                assert_eq!(PadicField::bitreverse_index(a, n), b);
                assert_eq!(PadicField::bitreverse_index(b, n), a);
            }
        }
    }

    #[test]
    fn test_ntt() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| rng.next_u32() as i32)
            .map(PadicField::new)
            .collect_vec();
        let mut b = a.clone();
        assert_eq!(a, b);

        let psi_rev = PadicField::bitreversed_powers(n);
        let psi_inv_rev = PadicField::bitreversed_powers_inverse(n);
        let ninv = PadicField::inverse_or_zero(PadicField::new(n as i32));
        PadicField::fft(&mut a, &psi_rev);
        PadicField::ifft(&mut a, &psi_inv_rev, ninv);
        assert_eq!(a, b);

        let x = PadicField::new(rng.next_u32() as i32);
        let y = PadicField::new(rng.next_u32() as i32);
        let mut c = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        PadicField::fft(&mut a, &psi_rev);
        PadicField::fft(&mut b, &psi_rev);
        PadicField::fft(&mut c, &psi_rev);

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
            .map(|_| PadicField::new(rng.gen_range(-20..20)))
            .collect_vec();
        let mut b = (0..n)
            .map(|_| PadicField::new(rng.gen_range(-20..20)))
            .collect_vec();

        let c = (Polynomial::new(a.clone()) * Polynomial::new(b.clone()))
            .reduce_by_cyclotomic(n)
            .coefficients;

        let psi_rev = PadicField::bitreversed_powers(n);
        PadicField::fft(&mut a, &psi_rev);
        PadicField::fft(&mut b, &psi_rev);
        let mut d = a.iter().zip(b.iter()).map(|(&l, &r)| l * r).collect_vec();
        let psi_inv_rev = PadicField::bitreversed_powers_inverse(n);
        let ninv = PadicField::new(n as i32).inverse_or_zero();
        PadicField::ifft(&mut d, &psi_inv_rev, ninv);

        let diff = |u: &[PadicField], v: &[PadicField]| {
            u.iter().zip(v.iter()).map(|(&l, &r)| l - r).collect_vec()
        };

        assert_eq!(diff(&c, &d), vec![PadicField::zero(); n]);
    }
}
