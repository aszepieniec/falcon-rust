use itertools::Itertools;
use num::traits::ConstZero;
use num::BigInt;
use num::BigUint;
use num::One;
use num::Signed;
use num::Zero;
use std::fmt::Debug;
use std::ops::Add;
use std::ops::Rem;
use std::ops::Sub;

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;
use crate::polynomial::Polynomial;
use crate::u32_field::U32Field;
use crate::u32_field::P;
use crate::u64_ring::U64Ring;
use crate::u64_ring::P2;

#[derive(Debug, Clone)]
pub(crate) struct PadicInteger {
    pub(crate) limbs: Vec<U32Field>,
}

impl From<BigUint> for PadicInteger {
    fn from(mut value: BigUint) -> Self {
        let mut limbs = vec![];
        while !value.is_zero() {
            let limb = value.clone().rem(&P);
            limbs.push(U32Field::new(
                limb.to_u32_digits()
                    .first()
                    .copied()
                    .unwrap_or(0)
                    .try_into()
                    .unwrap(),
            ));
            value = (value - limb) / P;
        }
        Self { limbs }
    }
}

impl PadicInteger {
    pub(crate) fn from_bigint(value: BigInt, num_limbs: usize) -> Self {
        let mut borrow = value.is_negative();
        let mut padic_integer = PadicInteger::from(BigUint::try_from(value.abs()).unwrap());
        while padic_integer.limbs.len() < num_limbs {
            padic_integer.limbs.push(U32Field::ZERO);
        }
        if borrow {
            for limb in &mut padic_integer.limbs {
                (*limb, borrow) = (-*limb + U32Field::new(borrow.into()), limb.is_zero());
            }
        }
        padic_integer
    }

    /// Negate p-adic integer, represented as a slice of limbs, given a
    /// modulus P.
    ///
    /// When the limb type is a finite ring, the modulus should be zero.
    /// Otherwise it should be P.
    fn negate<T>(limbs: &mut [T], p: T)
    where
        T: Debug + Copy + From<bool> + Add<T, Output = T> + Zero + Sub<T, Output = T>,
    {
        let mut borrow = false;
        for limb in limbs.iter_mut() {
            let subtrahend = *limb + T::from(borrow);
            if !subtrahend.is_zero() {
                *limb = p - subtrahend;
                borrow = true;
            } else {
                *limb = T::zero();
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PadicPolynomial {
    pub(crate) coefficient_limbs: Vec<U64Ring>,
    pub(crate) num_coeffs: usize,
    pub(crate) limbs_per_coeff: usize,
}

impl PadicPolynomial {
    pub(crate) fn new<const N: usize, const M: usize, T: Into<U64Ring> + Copy>(
        coeff_limbs: [[T; N]; M],
    ) -> Self {
        PadicPolynomial {
            coefficient_limbs: coeff_limbs
                .as_flattened()
                .iter()
                .map(|t| std::convert::Into::<U64Ring>::into(*t))
                .collect_vec(),
            num_coeffs: N,
            limbs_per_coeff: M,
        }
    }
    fn neg(mut self) -> Self {
        self.coefficient_limbs.iter_mut().for_each(|c| {
            *c = -*c;
        });
        self
    }
    fn add_assign(mut self, other: &Self) -> Self {
        self.coefficient_limbs
            .iter_mut()
            .zip(other.coefficient_limbs.iter())
            .for_each(|(a, b)| {
                *a += *b;
            });
        self
    }

    fn hadamard_mul_assign(mut self, other: &Self) -> Self {
        self.coefficient_limbs
            .iter_mut()
            .zip(other.coefficient_limbs.iter())
            .for_each(|(a, b)| {
                *a *= *b;
            });
        self
    }

    fn normalize(&mut self) {
        let mut limbs = vec![U64Ring::ZERO; self.limbs_per_coeff];
        let psi_inverse_num_limbs = U64Ring::bitreversed_powers_inverse(self.limbs_per_coeff);
        let num_limbs_inverse = U64Ring::new(self.limbs_per_coeff.try_into().unwrap());
        let psi_rev_num_limbs = U64Ring::bitreversed_powers(self.limbs_per_coeff);
        for i in 0..self.num_coeffs {
            for (j, limb) in limbs.iter_mut().enumerate() {
                *limb = self.coefficient_limbs[j * self.num_coeffs + i];
            }

            U64Ring::ifft(&mut limbs, &psi_inverse_num_limbs, num_limbs_inverse);

            let mut carry = U32Field::zero();
            for limb in limbs.iter_mut() {
                let temp = *limb + carry;
                carry = temp.hi();
                *limb = temp.lo().into();
            }

            U64Ring::fft(&mut limbs, &psi_rev_num_limbs);

            for (j, limb) in limbs.iter().enumerate() {
                self.coefficient_limbs[j * self.num_coeffs + i] = *limb;
            }
        }
    }

    fn normalized(mut self) -> Self {
        self.normalize();
        self
    }
}

impl PartialEq for PadicPolynomial {
    fn eq(&self, other: &Self) -> bool {
        let self_coefficient_limbs = self.clone().normalized().coefficient_limbs;
        let other_coefficient_limbs = other.clone().normalized().coefficient_limbs;
        self_coefficient_limbs == other_coefficient_limbs
            && self.num_coeffs == other.num_coeffs
            && self.limbs_per_coeff == other.limbs_per_coeff
    }
}

impl From<Polynomial<BigInt>> for PadicPolynomial {
    fn from(value: Polynomial<BigInt>) -> Self {
        let num_coefficients = value.coefficients.len();
        let max_bitsize = value
            .coefficients
            .iter()
            .map(|c| c.bits())
            .max()
            .unwrap_or(0);
        let num_limbs = 2 * ((max_bitsize + 30) / 31).next_power_of_two() as usize;

        let psi_rev_num_limbs = U64Ring::bitreversed_powers(num_limbs);
        let mut coefficient_limbs = vec![U64Ring::ZERO; num_coefficients * num_limbs];
        for (i, coefficient) in value.coefficients.into_iter().enumerate() {
            let was_negative = coefficient.is_negative();
            let expansion = PadicInteger::from_bigint(coefficient.abs(), num_limbs);
            println!("expansion: {:?}", expansion);
            let mut limbs = expansion.limbs.into_iter().map(|l| l.into()).collect_vec();
            println!("limbs: {:?}", limbs);
            if was_negative {
                PadicInteger::negate(&mut limbs, U64Ring::new(P as i64));
            }
            println!("negated: {:?}", limbs);
            U64Ring::fft(&mut limbs, &psi_rev_num_limbs);
            println!("after fft: {:?}", limbs);
            for (j, limb) in limbs.into_iter().enumerate() {
                coefficient_limbs[j * num_coefficients + i] = limb;
            }
        }

        let psi_rev_num_coefficients = U64Ring::bitreversed_powers(num_coefficients);
        for i in 0..num_limbs {
            U64Ring::fft(
                &mut coefficient_limbs[i * num_coefficients..(i + 1) * num_coefficients],
                &psi_rev_num_coefficients,
            );
        }

        Self {
            coefficient_limbs,
            num_coeffs: num_coefficients,
            limbs_per_coeff: num_limbs,
        }
    }
}

impl From<PadicPolynomial> for Polynomial<BigInt> {
    fn from(mut value: PadicPolynomial) -> Self {
        println!("---");
        let num_coefficients = value.num_coeffs;
        let num_limbs = value.limbs_per_coeff;

        let psi_inv_rev_num_coefficients = U64Ring::bitreversed_powers_inverse(num_coefficients);
        let num_coefficients_inv = U64Ring::new(num_coefficients as i64).inverse_or_zero();
        for i in 0..num_limbs {
            U64Ring::ifft(
                &mut value.coefficient_limbs[i * num_coefficients..(i + 1) * num_coefficients],
                &psi_inv_rev_num_coefficients,
                num_coefficients_inv,
            );
        }

        let psi_inv_rev_num_limbs = U64Ring::bitreversed_powers_inverse(num_limbs);
        let num_limbs_inv = U64Ring::new(num_limbs as i64).inverse_or_zero();
        let mut coefficients = vec![];
        for i in 0..num_coefficients {
            let mut limbs = vec![U64Ring::ZERO; num_limbs];
            for (j, limb) in limbs.iter_mut().enumerate() {
                *limb = value.coefficient_limbs[j * num_coefficients + i];
            }
            println!("before ifft: {:?}", limbs);
            U64Ring::ifft(&mut limbs, &psi_inv_rev_num_limbs, num_limbs_inv);
            println!("after ifft: {:?}", limbs);

            let mut carry = U32Field::ZERO;
            for limb in limbs.iter_mut() {
                let temp = *limb + carry;
                carry = temp.hi();
                *limb = temp.lo().into();
            }
            println!("after carry: {:?}", limbs);

            let was_negative = *limbs.last().unwrap() == U64Ring::from(-U32Field::one());
            if was_negative {
                PadicInteger::negate(&mut limbs, U64Ring::new(P as i64));
            }
            println!("after negating: {:?}", limbs);

            let slice = limbs.into_iter().map(|limb| limb.lo().0).collect_vec();
            let mut bi = BigInt::from(BigUint::from_slice(&slice));

            if was_negative {
                bi = -bi;
            }

            coefficients.push(bi);
        }

        Polynomial { coefficients }
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num::BigInt;
    use num::Zero;
    use proptest::collection::vec;
    use test_strategy::proptest;

    use crate::polynomial::Polynomial;
    use crate::u64_ring::P2;

    use super::{PadicInteger, PadicPolynomial, U64Ring};

    fn can_convert_to_and_fro(big_integer_polynomial: Polynomial<BigInt>) -> bool {
        println!("big integer poly: {:?}", big_integer_polynomial);
        let padic_polynomial = PadicPolynomial::from(big_integer_polynomial.clone());
        println!("padic polynomial: {:?}", padic_polynomial);
        let big_integer_polynomial_again = Polynomial::<BigInt>::from(padic_polynomial);
        println!("big int poly again: {:?}", big_integer_polynomial_again);
        big_integer_polynomial == big_integer_polynomial_again
    }

    fn can_convert_fro_and_to(padic_polynomial: PadicPolynomial) -> bool {
        let big_integer_polynomial = Polynomial::<BigInt>::from(padic_polynomial.clone());
        let padic_polynomial_again = PadicPolynomial::from(big_integer_polynomial);
        padic_polynomial == padic_polynomial_again
    }

    #[test]
    fn test_conversion() {
        let bip = |x: Vec<i64>| Polynomial::new(x.into_iter().map(BigInt::from).collect_vec());
        assert!(can_convert_to_and_fro(bip(vec![0])));
        assert!(can_convert_to_and_fro(bip(vec![5])));
        assert!(can_convert_to_and_fro(bip(vec![-5])));
        assert!(can_convert_to_and_fro(bip(vec![0, 5])));
        assert!(can_convert_to_and_fro(bip(vec![0, -5])));
        assert!(can_convert_to_and_fro(bip(vec![5, 0])));
        assert!(can_convert_to_and_fro(bip(vec![-5, 0])));
        assert!(can_convert_to_and_fro(bip(vec![5, 11])));
        assert!(can_convert_to_and_fro(bip(vec![-5, 11])));
        assert!(can_convert_to_and_fro(bip(vec![5, -11])));
        assert!(can_convert_to_and_fro(bip(vec![-5, -11])));
    }

    #[proptest(cases = 10000)]
    fn test_negate_padic_integer(#[strategy(vec(-10_i64..10, 1..10))] limbs: Vec<i64>) {
        let pos_limbs = limbs.into_iter().map(U64Ring::new).collect_vec();
        let mut neg_limbs = pos_limbs.clone();
        PadicInteger::negate(&mut neg_limbs, U64Ring::zero());
        let mut carry = false;
        for (i, (p, n)) in pos_limbs.into_iter().zip(neg_limbs.into_iter()).enumerate() {
            assert_eq!(
                p + n + U64Ring::from(carry),
                U64Ring::zero(),
                "test failed for limb {i}.\np: {:?}\nn: {:?}\nexpected: 0\nobserved: {}\ncarry: {carry}",
                p, n, (p + n + U64Ring::from(carry)).0,
            );
            carry = (p.0 + n.0 + u64::from(carry)) == P2;
        }
    }
}
