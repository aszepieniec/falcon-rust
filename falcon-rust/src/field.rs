use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use rand_distr::{
    num_traits::{One, Zero},
    Distribution, Standard,
};

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;

/// q is the integer modulus which is used in Falcon.
pub(crate) const Q: u32 = 12 * 1024 + 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Felt(u32);

impl Felt {
    pub const fn new(value: i16) -> Self {
        let gtz_bool = value >= 0;
        let gtz_int = gtz_bool as i16;
        let gtz_sign = gtz_int - ((!gtz_bool) as i16);
        let reduced = gtz_sign * ((gtz_sign * value) % (Q as i16));
        let canonical_representative = (reduced + (Q as i16) * (1 - gtz_int)) as u32;
        Felt(canonical_representative)
    }

    pub const fn value(&self) -> i16 {
        self.0 as i16
    }

    pub fn balanced_value(&self) -> i16 {
        let value = self.value();
        let g = (value > ((Q as i16) / 2)) as i16;
        value - (Q as i16) * g
    }

    pub const fn multiply(&self, other: Self) -> Self {
        Felt((self.0 * other.0) % Q)
    }
}

impl From<usize> for Felt {
    fn from(value: usize) -> Self {
        Felt::new(value as i16)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for Felt {
    fn add(self, rhs: Self) -> Self::Output {
        let (s, _) = self.0.overflowing_add(rhs.0);
        let (d, n) = s.overflowing_sub(Q);
        let (r, _) = d.overflowing_add(Q * (n as u32));
        Felt(r)
    }

    type Output = Self;
}

impl AddAssign for Felt {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Felt {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl SubAssign for Felt {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for Felt {
    type Output = Felt;

    fn neg(self) -> Self::Output {
        let is_nonzero = self.0 != 0;
        let r = Q - self.0;
        Felt(r * (is_nonzero as u32))
    }
}

impl Mul for Felt {
    fn mul(self, rhs: Self) -> Self::Output {
        Felt((self.0 * rhs.0) % Q)
    }

    type Output = Self;
}

impl MulAssign for Felt {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Zero for Felt {
    fn zero() -> Self {
        Felt::new(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}
impl One for Felt {
    fn one() -> Self {
        Felt::new(1)
    }
}

impl Distribution<Felt> for Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Felt {
        Felt::new(((rng.next_u32() >> 1) % Q) as i16)
    }
}

impl Display for Felt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("{}", self.value()))
    }
}

impl Inverse for Felt {
    fn inverse_or_zero(self) -> Self {
        // q-2 = 0b10 11 11 11  11 11 11
        let two = self.multiply(self);
        let three = two.multiply(self);
        let six = three.multiply(three);
        let twelve = six.multiply(six);
        let fifteen = twelve.multiply(three);
        let thirty = fifteen.multiply(fifteen);
        let sixty = thirty.multiply(thirty);
        let sixty_three = sixty.multiply(three);

        let sixty_three_sq = sixty_three.multiply(sixty_three);
        let sixty_three_qu = sixty_three_sq.multiply(sixty_three_sq);
        let sixty_three_oc = sixty_three_qu.multiply(sixty_three_qu);
        let sixty_three_hx = sixty_three_oc.multiply(sixty_three_oc);
        let sixty_three_tt = sixty_three_hx.multiply(sixty_three_hx);
        let sixty_three_sf = sixty_three_tt.multiply(sixty_three_tt);

        let all_ones = sixty_three_sf.multiply(sixty_three);
        let two_e_twelve = all_ones.multiply(self);
        let two_e_thirteen = two_e_twelve.multiply(two_e_twelve);

        two_e_thirteen.multiply(all_ones)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Div for Felt {
    type Output = Felt;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            panic!("Cannot divide by zero");
        } else {
            self.multiply(rhs.inverse_or_zero())
        }
    }
}

impl CyclotomicFourier for Felt {
    fn primitive_root_of_unity(n: usize) -> Self {
        let log2n = n.ilog2();
        assert!(log2n <= 12);
        // and 1331 is a twelfth root of unity
        let mut a = Felt::new(1331);
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
        field::{Felt, Q},
        inverse::Inverse,
        polynomial::Polynomial,
    };
    use num::Zero;

    #[test]
    fn test_value() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let mut value = (rng.next_u32() & 0x3fff) as i16;
            if rng.next_u32() % 2 == 1 {
                value *= -1;
            }
            let felt = Felt::new(value);
            assert_eq!(
                0,
                (felt.value() - value) % (Q as i16),
                "value: {value} but got {}",
                felt.value()
            );
        }
    }

    #[test]
    fn test_add() {
        let mut rng = thread_rng();
        let a_value = (rng.next_u32() % 0x0fff) as i16;
        let b_value = (rng.next_u32() % 0x0fff) as i16;
        let a = Felt::new(a_value);
        let b = Felt::new(b_value);
        assert_eq!(
            a + b,
            Felt::new(a.value() + b.value()),
            "a: {a_value}, b: {b_value}, c: {}",
            ((a_value + b_value) as u32) % Q
        );
    }

    #[test]
    fn test_specific_neg() {
        let br = Felt::new(2958);
        let min_br = Felt::new(9331);
        assert_eq!(-br, min_br, "-{} == {} =/= {}", br, -br, min_br);
    }

    #[test]
    fn test_mul() {
        let mut rng = thread_rng();
        for _ in 0..1000 {
            let a_value = (rng.next_u32() % 0x3fff) as i16;
            let b_value = (rng.next_u32() % 0x3fff) as i16;
            let product = (((a_value as u32) * (b_value as u32)) % Q) as i16;
            let a = Felt::new(a_value);
            let b = Felt::new(b_value);
            assert_eq!(
                a * b,
                Felt::new(product),
                "{} =/= {}",
                a * b,
                Felt::new(product)
            );
        }
    }

    #[test]
    fn test_batch_inverse() {
        let mut rng = thread_rng();
        let a: [Felt; 64] = (0..64).map(|_| rng.gen()).collect_vec().try_into().unwrap();
        let b_batch = Felt::batch_inverse_or_zero(&a);
        let b_regular = a.iter().map(|e| e.inverse_or_zero()).collect_vec();
        assert_eq!(b_batch.to_vec(), b_regular);
    }

    #[test]
    fn test_inverse() {
        let mut rng = thread_rng();
        let a: Felt = rng.gen();
        let b = a.inverse_or_zero();

        assert_eq!(a * b * a, a);
        assert_eq!(a * b * b, b);
    }

    #[test]
    fn test_primitive_nth_root_of_unity() {
        for log2n in 0..=12 {
            let n = 1 << log2n;
            let mut root = Felt::primitive_root_of_unity(n);
            for i in 0..log2n {
                assert_ne!(root, Felt::one(), "log2n: {log2n} and i: {i}");
                root *= root;
            }
            assert_eq!(root, Felt::one());
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
                assert_eq!(Felt::bitreverse_index(a, n), b);
                assert_eq!(Felt::bitreverse_index(b, n), a);
            }
        }
    }

    #[test]
    fn test_ntt() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| rng.next_u32() as i16)
            .map(Felt::new)
            .collect_vec();
        let mut b = a.clone();
        assert_eq!(a, b);

        let psi_rev = Felt::bitreversed_powers(n);
        let psi_inv_rev = Felt::bitreversed_powers_inverse(n);
        let ninv = Felt::inverse_or_zero(Felt::new(n as i16));
        Felt::fft(&mut a, &psi_rev);
        Felt::ifft(&mut a, &psi_inv_rev, ninv);
        assert_eq!(a, b);

        let x = Felt::new(rng.next_u32() as i16);
        let y = Felt::new(rng.next_u32() as i16);
        let mut c = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        Felt::fft(&mut a, &psi_rev);
        Felt::fft(&mut b, &psi_rev);
        Felt::fft(&mut c, &psi_rev);

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
            .map(|_| Felt::new(rng.gen_range(-20..20)))
            .collect_vec();
        let mut b = (0..n)
            .map(|_| Felt::new(rng.gen_range(-20..20)))
            .collect_vec();

        let c = (Polynomial::new(a.clone()) * Polynomial::new(b.clone()))
            .reduce_by_cyclotomic(n)
            .coefficients;

        let psi_rev = Felt::bitreversed_powers(n);
        Felt::fft(&mut a, &psi_rev);
        Felt::fft(&mut b, &psi_rev);
        let mut d = a.iter().zip(b.iter()).map(|(&l, &r)| l * r).collect_vec();
        let psi_inv_rev = Felt::bitreversed_powers_inverse(n);
        let ninv = Felt::new(n as i16).inverse_or_zero();
        Felt::ifft(&mut d, &psi_inv_rev, ninv);

        let diff =
            |u: &[Felt], v: &[Felt]| u.iter().zip(v.iter()).map(|(&l, &r)| l - r).collect_vec();

        assert_eq!(diff(&c, &d), vec![Felt::zero(); n]);
    }
}
