use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{One, Zero};
use rand::distr::{Distribution, StandardUniform};

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;

/// q is the integer modulus which is used in Falcon.
#[doc(hidden)]
pub const Q: u16 = 12 * 1024 + 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[doc(hidden)]
pub struct Felt(u16);

impl Felt {
    pub(crate) const Q: u16 = super::falcon_field::Q;
    pub(crate) const LOG2R: usize = 16;
    pub(crate) const R_MINUS_1: u32 = 0xffff_u32;
    pub(crate) const NEGQINV_MODR: u32 = 12287;
    pub(crate) const R_SQ_MODQ: u16 = 10952;

    pub const fn new(value: u16) -> Self {
        let reduced = if value >= Self::Q {
            value - Self::Q
        } else {
            value
        };

        Felt(Self::montyred(reduced, Self::R_SQ_MODQ))
    }

    pub const fn value(&self) -> u16 {
        Self::montyred(self.0, 1)
    }

    pub fn balanced_value(&self) -> i16 {
        let value = self.value() as i16;
        let g = (value > ((Self::Q as i16) / 2)) as i16;
        value - (Q as i16) * g
    }

    pub const fn multiply(&self, other: Self) -> Self {
        Felt(Self::montyred(self.0, other.0))
    }

    /// Multiply by 2⁻¹ mod Q without a full Montgomery reduction.
    ///
    /// The stored value x is an integer in {0, …, Q−1}.  Whether it is in
    /// Montgomery form or not, the following identity holds:
    ///
    ///   (x + Q·(x & 1)) >> 1  ≡  x · 2⁻¹  (mod Q)
    ///
    /// Proof: if x is even, the result is x/2 and 2·(x/2) = x ≡ x (mod Q).
    /// If x is odd, Q is also odd so x+Q is even, and 2·(x+Q)/2 = x+Q ≡ x (mod Q).
    /// The result is always in {0, …, Q−1}, so no further reduction is needed.
    ///
    /// Cost: one AND, one negate-mask, one ADD, one shift — much cheaper than a
    /// full Montgomery multiplication.
    pub const fn half(self) -> Self {
        let x = self.0 as u32;
        // If x is odd, add Q to make it even before halving.
        let r = (x + ((x & 1).wrapping_neg() & Self::Q as u32)) >> 1;
        Felt(r as u16)
    }

    /// Compute the product and montgomery-reduce it.
    ///
    /// Given two field elements 0 ≤ a, b < Q, this function computes a number
    /// 0 ≤ c < Q such that c R = ab mod Q, where R is 2^16.
    const fn montyred(a: u16, b: u16) -> u16 {
        debug_assert!(a < Self::Q);
        debug_assert!(b < Self::Q);

        let product: u32 = (a as u32) * (b as u32);

        // product ≤ (Q-1)²

        // We want to add a magic number to this product such that
        //  - the magic number is a multiple of q, so the residue class does not
        //    change;
        //  - the sum becomes divisible by r (rightmost logr bits are zero) so
        //    that in a later step we can shift (cheap) instead of divide
        //    (expensive).
        //
        // These requirements translate to
        //  - the rightmost logr bits of a·b equals the rightmost logr bits of
        //    -magic_number;
        //  - magic_number = q · something.
        //
        // So:
        //  something is congruent to magic_number · q^-1 modulo r
        //                            a·b · (-q^-1) modulo r

        let tail = product & Self::R_MINUS_1;

        let cofactor = (Self::NEGQINV_MODR * tail) & Self::R_MINUS_1;

        // cofactor ≤ R-1

        let magic_sum = product + cofactor * (Self::Q as u32);

        // magic_sum ≤ (Q-1)² + (R-1)·Q
        //           = Q² -2·Q + 1 + R·Q - Q  (expand)
        //           = Q·R + Q·(Q-3) + 1      (group)
        //           < 2·Q·R                  (because Q-3 < R)

        let mut g = (magic_sum >> Self::LOG2R) as u16;

        // g < 2·Q
        // So one conditional subtraction suffices.
        // The compiler is smart enough to turn this if-statement into
        // branch-free code.
        if g >= Self::Q {
            g -= Self::Q;
        }

        g
    }
}

impl From<usize> for Felt {
    fn from(value: usize) -> Self {
        Felt::new(value as u16)
    }
}

impl From<i16> for Felt {
    fn from(value: i16) -> Self {
        if value >= 0 {
            Self::new(value as u16)
        } else {
            -Self::new((-value) as u16)
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Add for Felt {
    fn add(self, rhs: Self) -> Self::Output {
        let (s, _) = self.0.overflowing_add(rhs.0);
        let (d, n) = s.overflowing_sub(Q);
        let (r, _) = d.overflowing_add(Q * (n as u16));
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
        Felt(r * (is_nonzero as u16))
    }
}

impl Mul for Felt {
    fn mul(self, rhs: Self) -> Self::Output {
        Felt(Self::montyred(self.0, rhs.0))
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

impl Distribution<Felt> for StandardUniform {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Felt {
        Felt::new((rng.next_u32() as u16 >> 1) % Q)
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
    use rand::{rng, Rng, RngExt};

    use crate::{
        cyclotomic_fourier::CyclotomicFourier,
        falcon_field::{Felt, Q},
        inverse::Inverse,
        polynomial::Polynomial,
    };
    use num::Zero;

    impl Felt {
        pub(crate) const RINV_MODQ: u32 = 2304;
    }

    #[test]
    fn test_montyred() {
        let mut rng = rng();
        for _ in 0..1000 {
            let a = rng.random_range(0..Felt::Q);
            let b = rng.random_range(0..Felt::Q);
            let c =
                ((((a as u64) * (b as u64)) * (Felt::RINV_MODQ as u64)) % (Felt::Q as u64)) as u16;
            let cc = Felt::montyred(a, b);
            assert_eq!(c, cc);
        }
    }

    #[test]
    fn test_value() {
        let mut rng = rng();
        for _ in 0..1000 {
            let value = rng.random_range(0..Felt::Q);
            let felt = Felt::new(value);
            assert_eq!(
                value,
                felt.value(),
                "value: {value} but got {}",
                felt.value()
            );
        }
    }

    #[test]
    fn test_add() {
        let mut rng = rng();
        let a_value = (rng.next_u32() % 0x0fff) as u16;
        let b_value = (rng.next_u32() % 0x0fff) as u16;
        let a = Felt::new(a_value);
        let b = Felt::new(b_value);
        assert_eq!(
            a + b,
            Felt::new(a.value() + b.value()),
            "a: {a_value}, b: {b_value}, c: {}",
            ((a_value + b_value) as u16) % Q
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
        let mut rng = rng();
        for _ in 0..1000 {
            let a_value = (rng.next_u32() % 0x3fff) as u16;
            let b_value = (rng.next_u32() % 0x3fff) as u16;
            let product = (((a_value as u32) * (b_value as u32)) % (Felt::Q as u32)) as u16;
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
        let mut rng = rng();
        let a: [Felt; 64] = (0..64)
            .map(|_| rng.random())
            .collect_vec()
            .try_into()
            .unwrap();
        let b_batch = Felt::batch_inverse_or_zero(&a);
        let b_regular = a.iter().map(|e| e.inverse_or_zero()).collect_vec();
        assert_eq!(b_batch.to_vec(), b_regular);
    }

    #[test]
    fn test_inverse() {
        let mut rng = rng();
        let a: Felt = rng.random();
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
        let mut rng = rng();
        let mut a = (0..n)
            .map(|_| rng.random_range(0..Felt::Q))
            .map(Felt::new)
            .collect_vec();
        let mut b = a.clone();
        assert_eq!(a, b);

        let psi_rev = Felt::bitreversed_powers(n);
        let psi_inv_rev = Felt::bitreversed_powers_inverse(n);
        let ninv = Felt::inverse_or_zero(Felt::new(n as u16));
        Felt::fft(&mut a, &psi_rev);
        Felt::ifft(&mut a, &psi_inv_rev, ninv);
        assert_eq!(a, b);

        let x = Felt::new(rng.random_range(0..Felt::Q));
        let y = Felt::new(rng.random_range(0..Felt::Q));
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
        let mut rng = rng();
        let mut a = (0..n)
            .map(|_| Felt::new(rng.random_range(0..20)))
            .collect_vec();
        let mut b = (0..n)
            .map(|_| Felt::new(rng.random_range(0..20)))
            .collect_vec();

        let c = (Polynomial::new(a.clone()) * Polynomial::new(b.clone()))
            .reduce_by_cyclotomic(n)
            .coefficients;

        let psi_rev = Felt::bitreversed_powers(n);
        Felt::fft(&mut a, &psi_rev);
        Felt::fft(&mut b, &psi_rev);
        let mut d = a.iter().zip(b.iter()).map(|(&l, &r)| l * r).collect_vec();
        let psi_inv_rev = Felt::bitreversed_powers_inverse(n);
        let ninv = Felt::new(n as u16).inverse_or_zero();
        Felt::ifft(&mut d, &psi_inv_rev, ninv);

        let diff =
            |u: &[Felt], v: &[Felt]| u.iter().zip(v.iter()).map(|(&l, &r)| l - r).collect_vec();

        assert_eq!(diff(&c, &d), vec![Felt::zero(); n]);
    }
}
