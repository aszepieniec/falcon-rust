use std::fmt::Display;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num::traits::ConstZero;
use rand_distr::{
    num_traits::{One, Zero},
    Distribution, Standard,
};

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;
use crate::padic_field::{PadicField, P};

/// Cube of [`padic_field::P`].
pub(crate) const P3: u64 = 172861751008264193_u64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PadicExtensionRing(pub(crate) u64);

impl PadicExtensionRing {
    pub(crate) const fn new(value: i64) -> Self {
        let gtz_bool = value >= 0;
        let gtz_int = gtz_bool as i64;
        let gtz_sign = gtz_int - ((!gtz_bool) as i64);
        let reduced = gtz_sign * ((gtz_sign * value) % (P3 as i64));
        let canonical_representative = (reduced + (P3 as i64) * (1 - gtz_int)) as u64;
        PadicExtensionRing(canonical_representative)
    }

    pub(crate) const fn value(&self) -> i64 {
        self.0 as i64
    }

    pub(crate) fn lo(&self) -> PadicField {
        PadicField::new((self.0 % (P as u64)).try_into().unwrap())
    }

    pub(crate) fn mid(&self) -> PadicField {
        let lo = self.lo().into();
        PadicField::new(((*self - lo).0 / (P as u64)).try_into().unwrap())
    }

    pub(crate) fn hi(&self) -> PadicField {
        let mid = self.mid().into();
        PadicField::new(
            ((*self - mid).0 / ((P as u64) * (P as u64)))
                .try_into()
                .unwrap(),
        )
    }
}

impl From<bool> for PadicExtensionRing {
    fn from(value: bool) -> Self {
        PadicExtensionRing::from(usize::from(value))
    }
}

impl From<usize> for PadicExtensionRing {
    fn from(value: usize) -> Self {
        PadicExtensionRing::new(value as i64)
    }
}

impl From<PadicField> for PadicExtensionRing {
    fn from(value: PadicField) -> Self {
        Self::new(value.value().into())
    }
}

impl ConstZero for PadicExtensionRing {
    const ZERO: Self = Self::new(0);
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for PadicExtensionRing {
    fn add(self, rhs: Self) -> Self::Output {
        let (s, _) = self.0.overflowing_add(rhs.0);
        let (d, n) = s.overflowing_sub(P3);
        let (r, _) = d.overflowing_add(P3 * (n as u64));
        PadicExtensionRing(r)
    }

    type Output = Self;
}

impl Add<PadicField> for PadicExtensionRing {
    type Output = Self;
    fn add(self, rhs: PadicField) -> Self::Output {
        self + Self::from(rhs)
    }
}

impl AddAssign for PadicExtensionRing {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for PadicExtensionRing {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl SubAssign for PadicExtensionRing {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for PadicExtensionRing {
    type Output = PadicExtensionRing;

    fn neg(self) -> Self::Output {
        let is_nonzero = self.0 != 0;
        let r = P3 - self.0;
        PadicExtensionRing(r * (is_nonzero as u64))
    }
}

impl Mul for PadicExtensionRing {
    fn mul(self, rhs: Self) -> Self::Output {
        PadicExtensionRing((((self.0 as u128) * (rhs.0 as u128)) % (P3 as u128)) as u64)
    }

    type Output = Self;
}

impl MulAssign for PadicExtensionRing {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Zero for PadicExtensionRing {
    fn zero() -> Self {
        PadicExtensionRing::new(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}
impl One for PadicExtensionRing {
    fn one() -> Self {
        PadicExtensionRing::new(1)
    }
}

impl Distribution<PadicExtensionRing> for Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> PadicExtensionRing {
        PadicExtensionRing::new(((rng.next_u64() >> 1) % P3) as i64)
    }
}

impl Display for PadicExtensionRing {
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

impl Inverse for PadicExtensionRing {
    fn inverse_or_zero(self) -> Self {
        let (a, _b, _g) = xgcd(self.0.try_into().unwrap(), P3.try_into().unwrap());
        Self::new(a)
    }
}

impl CyclotomicFourier for PadicExtensionRing {
    fn primitive_root_of_unity(n: usize) -> Self {
        let log2n = n.ilog2();
        assert!(log2n <= 12);
        // and 6408494027267741 is a (2^12)th root of unity
        let mut a = PadicExtensionRing::new(6408494027267741_i64);
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
        padic_extension_ring::{PadicExtensionRing, P3},
        polynomial::Polynomial,
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
            let uf = PadicExtensionRing::new(value);
            assert_eq!(
                0,
                (uf.value() - value) % (P3 as i64),
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
        let a = PadicExtensionRing::new(a_value);
        let b = PadicExtensionRing::new(b_value);
        assert_eq!(
            a + b,
            PadicExtensionRing::new(a.value() + b.value()),
            "a: {a_value}, b: {b_value}, c: {}",
            ((a_value + b_value) as u64) % P3
        );
    }

    #[test]
    fn test_mul() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let a_value = (rng.next_u64() % 0x3fffffffffff) as i64;
            let b_value = (rng.next_u64() % 0x3fffffffffff) as i64;
            let product = (((a_value as u128) * (b_value as u128)) % (P3 as u128)) as i64;
            let a = PadicExtensionRing::new(a_value);
            let b = PadicExtensionRing::new(b_value);
            assert_eq!(
                a * b,
                PadicExtensionRing::new(product),
                "{} =/= {}",
                a * b,
                PadicExtensionRing::new(product)
            );
        }
    }

    #[test]
    fn test_batch_inverse() {
        let mut rng = thread_rng();
        let a: [PadicExtensionRing; 64] =
            (0..64).map(|_| rng.gen()).collect_vec().try_into().unwrap();
        let b_batch = PadicExtensionRing::batch_inverse_or_zero(&a);
        let b_regular = a.iter().map(|e| e.inverse_or_zero()).collect_vec();
        assert_eq!(b_batch.to_vec(), b_regular);
    }

    #[test]
    fn test_inverse() {
        let mut rng = thread_rng();
        let a: PadicExtensionRing = rng.gen();
        let b = a.inverse_or_zero();

        assert_eq!(a * b * a, a);
        assert_eq!(a * b * b, b);
    }

    #[test]
    fn test_primitive_nth_root_of_unity() {
        for log2n in 0..=12 {
            let n = 1 << log2n;
            let mut root = PadicExtensionRing::primitive_root_of_unity(n);
            for i in 0..log2n {
                assert_ne!(root, PadicExtensionRing::one(), "log2n: {log2n} and i: {i}");
                root *= root;
            }
            assert_eq!(root, PadicExtensionRing::one());
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
                assert_eq!(PadicExtensionRing::bitreverse_index(a, n), b);
                assert_eq!(PadicExtensionRing::bitreverse_index(b, n), a);
            }
        }
    }

    #[test]
    fn test_ntt() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| rng.next_u64() as i64)
            .map(PadicExtensionRing::new)
            .collect_vec();
        let mut b = a.clone();
        assert_eq!(a, b);

        let psi_rev = PadicExtensionRing::bitreversed_powers(n);
        let psi_inv_rev = PadicExtensionRing::bitreversed_powers_inverse(n);
        let ninv = PadicExtensionRing::inverse_or_zero(PadicExtensionRing::new(n as i64));
        PadicExtensionRing::fft(&mut a, &psi_rev);
        PadicExtensionRing::ifft(&mut a, &psi_inv_rev, ninv);
        assert_eq!(a, b);

        let x = PadicExtensionRing::new(rng.next_u32() as i64);
        let y = PadicExtensionRing::new(rng.next_u32() as i64);
        let mut c = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        PadicExtensionRing::fft(&mut a, &psi_rev);
        PadicExtensionRing::fft(&mut b, &psi_rev);
        PadicExtensionRing::fft(&mut c, &psi_rev);

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
            .map(|_| PadicExtensionRing::new(rng.gen_range(-20..20)))
            .collect_vec();
        let mut b = (0..n)
            .map(|_| PadicExtensionRing::new(rng.gen_range(-20..20)))
            .collect_vec();

        let c = (Polynomial::new(a.clone()) * Polynomial::new(b.clone()))
            .reduce_by_cyclotomic(n)
            .coefficients;

        let psi_rev = PadicExtensionRing::bitreversed_powers(n);
        PadicExtensionRing::fft(&mut a, &psi_rev);
        PadicExtensionRing::fft(&mut b, &psi_rev);
        let mut d = a.iter().zip(b.iter()).map(|(&l, &r)| l * r).collect_vec();
        let psi_inv_rev = PadicExtensionRing::bitreversed_powers_inverse(n);
        let ninv = PadicExtensionRing::new(n as i64).inverse_or_zero();
        PadicExtensionRing::ifft(&mut d, &psi_inv_rev, ninv);

        let diff = |u: &[PadicExtensionRing], v: &[PadicExtensionRing]| {
            u.iter().zip(v.iter()).map(|(&l, &r)| l - r).collect_vec()
        };

        assert_eq!(diff(&c, &d), vec![PadicExtensionRing::zero(); n]);
    }
}
