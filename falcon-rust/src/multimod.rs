use std::{
    cmp::min,
    fmt::Display,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Shl, Sub, SubAssign},
};

use itertools::Itertools;
use num::{
    bigint::{Sign, ToBigInt},
    BigInt, BigUint, FromPrimitive, Integer, One, ToPrimitive, Zero,
};

use crate::{
    cyclotomic_fourier::CyclotomicFourier,
    inverse::Inverse,
    polynomial::Polynomial,
    product_tree::ProductTree,
    residue_number_system::{BEZOUT_COEFFICIENTS, PARTIAL_PRODUCTS},
};
use crate::{
    product_tree::MASTER_TREE,
    residue_number_system::{PrimeField, MODULI},
};

pub(crate) const MULTIMOD_MAX_CAPACITY: usize = MODULI.len() * 30;

/// Product of all moduli up to modulus n.
///
/// This function is not const because of dependencies, but we would like it to be.
fn product(n: usize) -> BigUint {
    let mut prod = BigUint::one();
    let mut i = 0;
    while i < n {
        prod *= BigUint::from_u32(MODULI[i]).unwrap();
        i += 1;
    }
    prod
}

/// An integer, smaller than the product of all moduli, represented
/// as the list of of residues modulo moduli {0, ..., num-1} for num <= N.
/// Every limb stores at least 30 bits, so the total capacity is 15360 bits.
#[derive(Debug, Clone)]
pub struct MultiModInt {
    /// List of residues modulo p_i
    limbs: Vec<u32>,

    /// Bound on the number of bits
    bitsize_bound: usize,
}

impl MultiModInt {
    pub fn from_biguint_with_num_limbs(int: BigUint, num: usize) -> MultiModInt {
        assert!(
            num <= MODULI.len(),
            "MultiModInt does not support more than {} limbs (given: {num}).",
            MODULI.len()
        );
        let bits = int.bits() as usize;
        let mut tree: &ProductTree = &MASTER_TREE;
        let num_halvings = MODULI.len().ilog2() - num.next_power_of_two().ilog2();
        for _ in 0..num_halvings {
            tree = tree.as_branch().left.as_ref();
        }
        let limbs = tree.reduce(&int)[0..num].to_vec();
        MultiModInt {
            limbs,
            bitsize_bound: bits,
        }
    }

    fn get_more_limbs(&self, total_count: usize) -> Vec<u32> {
        assert!(total_count <= MODULI.len());
        let int = BigInt::from(self.clone());
        let product = product(total_count);
        let uint = if int < BigInt::zero() {
            (product.to_bigint().unwrap() + int).to_biguint().unwrap()
        } else {
            int.to_biguint().unwrap()
        };
        let mut limbs = self.limbs.clone();
        for p in MODULI.into_iter().take(total_count).skip(self.limbs.len()) {
            limbs.push(
                uint.mod_floor(&BigUint::from_u32(p).unwrap())
                    .to_u32()
                    .unwrap(),
            );
        }
        limbs
    }

    pub fn num_limbs_needed(bits: usize) -> usize {
        let sign_margin = 1;
        let bits_per_limb = 30;
        (bits + (bits_per_limb - 1) + sign_margin) / bits_per_limb
    }

    /// Compute the integer's sign by first constructing it as a BigInt and then
    /// testing it against product/2.
    pub fn sign(&self) -> Sign {
        BigInt::from(self.clone()).sign()
    }

    pub fn zero_with_capacity(bits: usize) -> Self {
        let num = Self::num_limbs_needed(bits);
        Self {
            limbs: vec![0u32; num],
            bitsize_bound: 1,
        }
    }

    pub fn expand_capacity(&mut self, bits: usize) {
        let integer = BigInt::from(self.clone());
        for &p in MODULI
            .iter()
            .take(Self::num_limbs_needed(bits))
            .skip(self.limbs.len())
        {
            self.limbs
                .push(integer.mod_floor(&BigInt::from(p)).to_u32().unwrap());
        }
        self.bitsize_bound = bits
    }

    pub fn bin(&self) -> String {
        BigInt::from(self.clone()).to_str_radix(2)
    }

    /// Returns the bitsize of the integer, ignoring its sign.
    pub fn bits(&self) -> usize {
        // BigInts of more than 2^64 bits are never needed in this repo
        BigInt::from(self.clone())
            .bits()
            .try_into()
            .expect("cannot cast bitsize as usize")
    }

    fn to_biguint(&self) -> BigUint {
        self.limbs
            .iter()
            .enumerate()
            .map(|(i, x)| BEZOUT_COEFFICIENTS.get(i).unwrap().mul(*x))
            .sum::<BigUint>()
            .mod_floor(&PARTIAL_PRODUCTS.get(self.limbs.len() - 1).unwrap())
    }
}

impl PartialEq for MultiModInt {
    fn eq(&self, other: &Self) -> bool {
        BigInt::from(self.clone()) == BigInt::from(other.clone())
    }
}

impl Add for MultiModInt {
    type Output = MultiModInt;

    fn add(self, rhs: Self) -> Self::Output {
        let bits = 1 + usize::max(self.bitsize_bound, rhs.bitsize_bound);
        let num = usize::min(MODULI.len(), Self::num_limbs_needed(bits));
        assert!(
            Self::num_limbs_needed(bits) <= num,
            "Not enough limbs for multimod integers of {} bits",
            bits
        );
        let self_limbs = self.get_more_limbs(num);
        let other_limbs = rhs.get_more_limbs(num);
        let mut limbs = Vec::<u32>::with_capacity(num);
        for i in 0..num {
            limbs.push((self_limbs[i] + other_limbs[i]) % MODULI[i]);
        }
        let output = Self::Output {
            limbs,
            bitsize_bound: bits,
        };
        if self.sign() == rhs.sign() {
            assert_eq!(output.sign(), self.sign());
        }
        output
    }
}

impl AddAssign for MultiModInt {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone().add(rhs)
    }
}

impl Sum for MultiModInt {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), |a, b| a + b)
    }
}

impl Sub for MultiModInt {
    type Output = MultiModInt;

    fn sub(self, rhs: Self) -> Self::Output {
        let bits = 1 + usize::max(self.bitsize_bound, rhs.bitsize_bound);
        let num = usize::min(MODULI.len(), Self::num_limbs_needed(bits));
        assert!(
            Self::num_limbs_needed(bits) <= num,
            "Not enough limbs to hold multimod integers of {} bits",
            bits
        );
        let self_limbs = self.get_more_limbs(num);
        let other_limbs = rhs.get_more_limbs(num);
        let mut limbs = Vec::<u32>::with_capacity(num);
        for i in 0..num {
            limbs.push(
                (((self_limbs[i] as u64) + ((MODULI[i] as u64) - (other_limbs[i] as u64)))
                    % (MODULI[i] as u64)) as u32,
            );
        }
        Self::Output {
            limbs,
            bitsize_bound: bits,
        }
    }
}

impl SubAssign for MultiModInt {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone().sub(rhs)
    }
}

impl Mul for MultiModInt {
    type Output = MultiModInt;

    fn mul(self, rhs: Self) -> Self::Output {
        let bits = self.bitsize_bound + rhs.bitsize_bound;
        assert!(
            Self::num_limbs_needed(bits) <= MODULI.len(),
            "Not enough moduli available for {} bits",
            bits
        );
        let num = usize::min(MODULI.len(), Self::num_limbs_needed(bits));
        let self_limbs = self.get_more_limbs(num);
        let other_limbs = rhs.get_more_limbs(num);
        let mut limbs = Vec::<u32>::with_capacity(num);
        for i in 0..num {
            limbs.push(
                (((self_limbs[i] as u64) * (other_limbs[i] as u64)) % (MODULI[i] as u64)) as u32,
            );
        }
        let output = Self {
            limbs,
            bitsize_bound: bits,
        };
        if self.sign() == rhs.sign() && self.sign() != Sign::NoSign {
            assert_eq!(output.sign(), Sign::Plus);
        } else if self.sign() != rhs.sign()
            && self.sign() != Sign::NoSign
            && rhs.sign() != Sign::NoSign
        {
            assert_eq!(output.sign(), Sign::Minus);
        }
        output
    }
}

impl MulAssign for MultiModInt {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone().mul(rhs)
    }
}

impl Neg for MultiModInt {
    type Output = MultiModInt;

    fn neg(self) -> Self::Output {
        let mut limbs = vec![0u32; self.limbs.len()];
        for (i, p) in MODULI.into_iter().take(self.limbs.len()).enumerate() {
            limbs[i] = (p - self.limbs[i]) % p;
        }
        Self::Output {
            limbs,
            bitsize_bound: self.bitsize_bound,
        }
    }
}

impl Zero for MultiModInt {
    fn zero() -> Self {
        Self {
            limbs: vec![0u32; 1],
            bitsize_bound: 1,
        }
    }

    fn is_zero(&self) -> bool {
        self.limbs.iter().all(|l| l.is_zero())
    }
}

impl One for MultiModInt {
    fn one() -> Self {
        Self {
            limbs: vec![1u32; 1],
            bitsize_bound: 1,
        }
    }
}

impl Shl<usize> for MultiModInt {
    type Output = MultiModInt;

    fn shl(self, mut rhs: usize) -> Self::Output {
        let bits = self.bitsize_bound + rhs;
        let num_limbs_needed = Self::num_limbs_needed(bits);
        let mut limbs = if num_limbs_needed > self.limbs.len() {
            let num = usize::min(
                MODULI.len(),
                usize::max(num_limbs_needed, 2 * self.limbs.len()),
            );
            assert!(
                num >= num_limbs_needed,
                "requested bit width ({}) exceeds capacity ({})",
                bits,
                MODULI.len() * 30,
            );
            self.get_more_limbs(num)
        } else {
            self.limbs.clone()
        };
        while rhs > 0 {
            let shamt = min(32, rhs);
            rhs -= shamt;
            for (i, p) in MODULI.into_iter().take(limbs.len()).enumerate() {
                limbs[i] = (((limbs[i] as u64) << shamt) % (p as u64)) as u32;
            }
        }
        Self {
            limbs,
            bitsize_bound: bits,
        }
    }
}

impl Display for MultiModInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", BigUint::from(self.clone()))
    }
}

#[derive(Debug, Clone)]
pub struct MultiModIntConversionError<T: Display + Clone> {
    original: T,
}

impl<T: Display + Clone> Display for MultiModIntConversionError<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "MultiModIntConversionError({})",
            self.original
        ))
    }
}

impl TryFrom<BigUint> for MultiModInt {
    fn try_from(val: BigUint) -> Result<Self, Self::Error> {
        let bits = val.bits() as usize;
        if bits > MULTIMOD_MAX_CAPACITY {
            return Err(MultiModIntConversionError { original: val });
        }
        let num = Self::num_limbs_needed(bits);
        Ok(MultiModInt::from_biguint_with_num_limbs(val, num))
    }

    type Error = MultiModIntConversionError<BigUint>;
}

impl TryFrom<BigInt> for MultiModInt {
    type Error = MultiModIntConversionError<BigInt>;
    fn try_from(value: BigInt) -> Result<Self, Self::Error> {
        if value < BigInt::zero() {
            let Ok(mmi) = MultiModInt::try_from((-value.clone()).to_biguint().unwrap()) else {
                return Err(MultiModIntConversionError { original: value });
            };
            Ok(-mmi)
        } else {
            let Ok(mmi) = MultiModInt::try_from(value.to_biguint().unwrap()) else {
                return Err(MultiModIntConversionError { original: value });
            };
            Ok(mmi)
        }
    }
}

impl From<MultiModInt> for BigUint {
    fn from(value: MultiModInt) -> BigUint {
        value.to_biguint()
    }
}

impl From<MultiModInt> for BigInt {
    fn from(value: MultiModInt) -> Self {
        let num = value.limbs.len();
        let product = product(num);
        let biguint = value
            .limbs
            .iter()
            .enumerate()
            .map(|(i, x)| BEZOUT_COEFFICIENTS.get(i).unwrap().mul(*x))
            .sum::<BigUint>()
            .mod_floor(&product);
        // let biguint = value.to_biguint();
        if biguint > (product.clone() >> 1) {
            biguint.to_bigint().unwrap() - product.to_bigint().unwrap()
        } else {
            biguint.to_bigint().unwrap()
        }
    }
}

impl Display for Polynomial<MultiModInt> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.coefficients.len() {
            match i {
                0 => {
                    write!(f, "{}", self.coefficients[i])?;
                }
                1 => {
                    write!(f, " + {} * x", self.coefficients[i])?;
                }
                _ => {
                    write!(f, " + {} * x^{}", self.coefficients[i], i)?;
                }
            }
        }
        Ok(())
    }
}

impl Polynomial<MultiModInt> {
    fn cyclotomic_mul_ith_limb<const P: u32>(&self, other: &Self, result: &mut Self, i: usize) {
        let n = self.coefficients.len();

        let mut self_vector = self
            .coefficients
            .iter()
            .map(|c| PrimeField::<P>(c.limbs[i]))
            .collect_vec();
        let mut other_vector = other
            .coefficients
            .iter()
            .map(|c| PrimeField::<P>(c.limbs[i]))
            .collect_vec();

        let psi_rev = PrimeField::<P>::bitreversed_powers(n);
        PrimeField::<P>::fft(&mut self_vector, &psi_rev);
        PrimeField::<P>::fft(&mut other_vector, &psi_rev);

        let mut result_vector = self_vector
            .into_iter()
            .zip(other_vector)
            .map(|(a, b)| a * b)
            .collect_vec();

        let psi_rev_inv = PrimeField::<P>::bitreversed_powers_inverse(n);
        let ninv = PrimeField::<P>(n as u32).inverse_or_zero();
        PrimeField::<P>::ifft(&mut result_vector, &psi_rev_inv, ninv);

        for (j, rj) in result_vector.into_iter().enumerate() {
            result.coefficients[j].limbs[i] = rj.0;
        }
    }

    pub(crate) fn cyclotomic_mul(&self, other: &Self) -> Self {
        // calculate new bound on number of bits
        let lhs_max_bits = self
            .coefficients
            .iter()
            .map(|mmi| mmi.bitsize_bound)
            .max()
            .unwrap_or_default();
        let rhs_max_bits = other
            .coefficients
            .iter()
            .map(|mmi| mmi.bitsize_bound)
            .max()
            .unwrap_or_default();
        assert!(
            lhs_max_bits
                >= self
                    .coefficients
                    .iter()
                    .map(|mmi| BigInt::from(mmi.clone()).bits())
                    .max()
                    .unwrap() as usize
        );
        assert!(
            rhs_max_bits
                >= other
                    .coefficients
                    .iter()
                    .map(|mmi| BigInt::from(mmi.clone()).bits())
                    .max()
                    .unwrap() as usize
        );
        let bits =
            1 + self.coefficients.len() + other.coefficients.len() + lhs_max_bits + rhs_max_bits;

        // if not enough capacity for this many bits, activate limbs and recurse once
        if self
            .coefficients
            .iter()
            .chain(other.coefficients.iter())
            .map(|mmi| mmi.bitsize_bound)
            .min()
            .unwrap_or_default()
            < bits
        {
            let mut self_with_activated_limbs = self.clone();
            self_with_activated_limbs
                .coefficients
                .iter_mut()
                .for_each(|mmi| mmi.expand_capacity(bits));
            let mut other_with_activated_limbs = other.clone();
            other_with_activated_limbs
                .coefficients
                .iter_mut()
                .for_each(|mmi| mmi.expand_capacity(bits));
            Self::cyclotomic_mul_with_n_bits(
                &self_with_activated_limbs,
                &other_with_activated_limbs,
                bits,
            )
        } else {
            Self::cyclotomic_mul_with_n_bits(self, other, bits)
        }
    }

    pub fn cyclotomic_mul_with_n_bits(lhs: &Self, rhs: &Self, bits: usize) -> Self {
        let n = lhs.coefficients.len();
        let num_limbs = MultiModInt::num_limbs_needed(bits);
        // create object to be populated and returned
        let mut result =
            Polynomial::<MultiModInt>::new(vec![MultiModInt::zero_with_capacity(bits); n]);

        // populate with cyclotomic products
        let cyclotomic_mul_ith_limbs_p = [
            Self::cyclotomic_mul_ith_limb::<1073754113>,
            Self::cyclotomic_mul_ith_limb::<1073950721>,
            Self::cyclotomic_mul_ith_limb::<1073958913>,
            Self::cyclotomic_mul_ith_limb::<1073983489>,
            Self::cyclotomic_mul_ith_limb::<1074196481>,
            Self::cyclotomic_mul_ith_limb::<1074343937>,
            Self::cyclotomic_mul_ith_limb::<1074442241>,
            Self::cyclotomic_mul_ith_limb::<1074475009>,
            Self::cyclotomic_mul_ith_limb::<1074515969>,
            Self::cyclotomic_mul_ith_limb::<1074524161>,
            Self::cyclotomic_mul_ith_limb::<1074548737>,
            Self::cyclotomic_mul_ith_limb::<1074589697>,
            Self::cyclotomic_mul_ith_limb::<1074597889>,
            Self::cyclotomic_mul_ith_limb::<1074688001>,
            Self::cyclotomic_mul_ith_limb::<1074696193>,
            Self::cyclotomic_mul_ith_limb::<1074810881>,
            Self::cyclotomic_mul_ith_limb::<1074835457>,
            Self::cyclotomic_mul_ith_limb::<1074941953>,
            Self::cyclotomic_mul_ith_limb::<1075007489>,
            Self::cyclotomic_mul_ith_limb::<1075064833>,
            Self::cyclotomic_mul_ith_limb::<1075105793>,
            Self::cyclotomic_mul_ith_limb::<1075351553>,
            Self::cyclotomic_mul_ith_limb::<1075376129>,
            Self::cyclotomic_mul_ith_limb::<1075449857>,
            Self::cyclotomic_mul_ith_limb::<1075507201>,
            Self::cyclotomic_mul_ith_limb::<1075621889>,
            Self::cyclotomic_mul_ith_limb::<1075695617>,
            Self::cyclotomic_mul_ith_limb::<1075720193>,
            Self::cyclotomic_mul_ith_limb::<1075752961>,
            Self::cyclotomic_mul_ith_limb::<1076039681>,
            Self::cyclotomic_mul_ith_limb::<1076064257>,
            Self::cyclotomic_mul_ith_limb::<1076162561>,
            Self::cyclotomic_mul_ith_limb::<1076187137>,
            Self::cyclotomic_mul_ith_limb::<1076211713>,
            Self::cyclotomic_mul_ith_limb::<1076269057>,
            Self::cyclotomic_mul_ith_limb::<1076334593>,
            Self::cyclotomic_mul_ith_limb::<1076391937>,
            Self::cyclotomic_mul_ith_limb::<1076531201>,
            Self::cyclotomic_mul_ith_limb::<1076613121>,
            Self::cyclotomic_mul_ith_limb::<1076908033>,
            Self::cyclotomic_mul_ith_limb::<1076948993>,
            Self::cyclotomic_mul_ith_limb::<1077047297>,
            Self::cyclotomic_mul_ith_limb::<1077178369>,
            Self::cyclotomic_mul_ith_limb::<1077252097>,
            Self::cyclotomic_mul_ith_limb::<1077391361>,
            Self::cyclotomic_mul_ith_limb::<1077399553>,
            Self::cyclotomic_mul_ith_limb::<1077620737>,
            Self::cyclotomic_mul_ith_limb::<1077686273>,
            Self::cyclotomic_mul_ith_limb::<1077792769>,
            Self::cyclotomic_mul_ith_limb::<1077882881>,
            Self::cyclotomic_mul_ith_limb::<1078038529>,
            Self::cyclotomic_mul_ith_limb::<1078079489>,
            Self::cyclotomic_mul_ith_limb::<1078112257>,
            Self::cyclotomic_mul_ith_limb::<1078153217>,
            Self::cyclotomic_mul_ith_limb::<1078210561>,
            Self::cyclotomic_mul_ith_limb::<1078333441>,
            Self::cyclotomic_mul_ith_limb::<1078382593>,
            Self::cyclotomic_mul_ith_limb::<1078398977>,
            Self::cyclotomic_mul_ith_limb::<1078497281>,
            Self::cyclotomic_mul_ith_limb::<1078505473>,
            Self::cyclotomic_mul_ith_limb::<1078546433>,
            Self::cyclotomic_mul_ith_limb::<1078571009>,
            Self::cyclotomic_mul_ith_limb::<1078620161>,
            Self::cyclotomic_mul_ith_limb::<1078775809>,
            Self::cyclotomic_mul_ith_limb::<1078849537>,
            Self::cyclotomic_mul_ith_limb::<1078939649>,
            Self::cyclotomic_mul_ith_limb::<1079185409>,
            Self::cyclotomic_mul_ith_limb::<1079193601>,
            Self::cyclotomic_mul_ith_limb::<1079357441>,
            Self::cyclotomic_mul_ith_limb::<1079431169>,
            Self::cyclotomic_mul_ith_limb::<1079603201>,
            Self::cyclotomic_mul_ith_limb::<1079627777>,
            Self::cyclotomic_mul_ith_limb::<1079635969>,
            Self::cyclotomic_mul_ith_limb::<1079676929>,
            Self::cyclotomic_mul_ith_limb::<1079685121>,
            Self::cyclotomic_mul_ith_limb::<1079734273>,
            Self::cyclotomic_mul_ith_limb::<1079832577>,
            Self::cyclotomic_mul_ith_limb::<1079848961>,
            Self::cyclotomic_mul_ith_limb::<1079873537>,
            Self::cyclotomic_mul_ith_limb::<1079922689>,
            Self::cyclotomic_mul_ith_limb::<1080020993>,
            Self::cyclotomic_mul_ith_limb::<1080348673>,
            Self::cyclotomic_mul_ith_limb::<1080365057>,
            Self::cyclotomic_mul_ith_limb::<1080446977>,
            Self::cyclotomic_mul_ith_limb::<1080496129>,
            Self::cyclotomic_mul_ith_limb::<1080512513>,
            Self::cyclotomic_mul_ith_limb::<1080741889>,
            Self::cyclotomic_mul_ith_limb::<1080864769>,
            Self::cyclotomic_mul_ith_limb::<1080938497>,
            Self::cyclotomic_mul_ith_limb::<1080954881>,
            Self::cyclotomic_mul_ith_limb::<1081077761>,
            Self::cyclotomic_mul_ith_limb::<1081225217>,
            Self::cyclotomic_mul_ith_limb::<1081323521>,
            Self::cyclotomic_mul_ith_limb::<1081348097>,
            Self::cyclotomic_mul_ith_limb::<1081397249>,
            Self::cyclotomic_mul_ith_limb::<1081430017>,
            Self::cyclotomic_mul_ith_limb::<1081495553>,
            Self::cyclotomic_mul_ith_limb::<1081577473>,
            Self::cyclotomic_mul_ith_limb::<1081643009>,
            Self::cyclotomic_mul_ith_limb::<1081700353>,
            Self::cyclotomic_mul_ith_limb::<1081774081>,
            Self::cyclotomic_mul_ith_limb::<1081815041>,
            Self::cyclotomic_mul_ith_limb::<1081839617>,
            Self::cyclotomic_mul_ith_limb::<1081921537>,
            Self::cyclotomic_mul_ith_limb::<1082060801>,
            Self::cyclotomic_mul_ith_limb::<1082093569>,
            Self::cyclotomic_mul_ith_limb::<1082331137>,
            Self::cyclotomic_mul_ith_limb::<1082462209>,
            Self::cyclotomic_mul_ith_limb::<1082552321>,
            Self::cyclotomic_mul_ith_limb::<1082601473>,
            Self::cyclotomic_mul_ith_limb::<1082699777>,
            Self::cyclotomic_mul_ith_limb::<1082724353>,
            Self::cyclotomic_mul_ith_limb::<1082806273>,
            Self::cyclotomic_mul_ith_limb::<1082904577>,
            Self::cyclotomic_mul_ith_limb::<1082929153>,
            Self::cyclotomic_mul_ith_limb::<1083002881>,
            Self::cyclotomic_mul_ith_limb::<1083076609>,
            Self::cyclotomic_mul_ith_limb::<1083117569>,
            Self::cyclotomic_mul_ith_limb::<1083215873>,
            Self::cyclotomic_mul_ith_limb::<1083289601>,
            Self::cyclotomic_mul_ith_limb::<1083371521>,
            Self::cyclotomic_mul_ith_limb::<1083518977>,
            Self::cyclotomic_mul_ith_limb::<1083535361>,
            Self::cyclotomic_mul_ith_limb::<1083928577>,
            Self::cyclotomic_mul_ith_limb::<1083953153>,
            Self::cyclotomic_mul_ith_limb::<1084076033>,
            Self::cyclotomic_mul_ith_limb::<1084133377>,
            Self::cyclotomic_mul_ith_limb::<1084149761>,
            Self::cyclotomic_mul_ith_limb::<1084297217>,
            Self::cyclotomic_mul_ith_limb::<1084477441>,
            Self::cyclotomic_mul_ith_limb::<1084518401>,
            Self::cyclotomic_mul_ith_limb::<1084641281>,
            Self::cyclotomic_mul_ith_limb::<1084649473>,
            Self::cyclotomic_mul_ith_limb::<1084674049>,
            Self::cyclotomic_mul_ith_limb::<1084690433>,
            Self::cyclotomic_mul_ith_limb::<1084813313>,
            Self::cyclotomic_mul_ith_limb::<1084911617>,
            Self::cyclotomic_mul_ith_limb::<1084968961>,
            Self::cyclotomic_mul_ith_limb::<1085140993>,
            Self::cyclotomic_mul_ith_limb::<1085181953>,
            Self::cyclotomic_mul_ith_limb::<1085255681>,
            Self::cyclotomic_mul_ith_limb::<1085411329>,
            Self::cyclotomic_mul_ith_limb::<1085501441>,
            Self::cyclotomic_mul_ith_limb::<1085534209>,
            Self::cyclotomic_mul_ith_limb::<1085550593>,
            Self::cyclotomic_mul_ith_limb::<1085607937>,
            Self::cyclotomic_mul_ith_limb::<1085648897>,
            Self::cyclotomic_mul_ith_limb::<1085673473>,
            Self::cyclotomic_mul_ith_limb::<1085779969>,
            Self::cyclotomic_mul_ith_limb::<1085796353>,
            Self::cyclotomic_mul_ith_limb::<1085870081>,
            Self::cyclotomic_mul_ith_limb::<1085894657>,
            Self::cyclotomic_mul_ith_limb::<1085952001>,
            Self::cyclotomic_mul_ith_limb::<1086066689>,
            Self::cyclotomic_mul_ith_limb::<1086115841>,
            Self::cyclotomic_mul_ith_limb::<1086345217>,
            Self::cyclotomic_mul_ith_limb::<1086410753>,
            Self::cyclotomic_mul_ith_limb::<1086443521>,
            Self::cyclotomic_mul_ith_limb::<1086566401>,
            Self::cyclotomic_mul_ith_limb::<1086713857>,
            Self::cyclotomic_mul_ith_limb::<1086763009>,
            Self::cyclotomic_mul_ith_limb::<1086853121>,
            Self::cyclotomic_mul_ith_limb::<1086902273>,
            Self::cyclotomic_mul_ith_limb::<1086935041>,
            Self::cyclotomic_mul_ith_limb::<1087025153>,
            Self::cyclotomic_mul_ith_limb::<1087148033>,
            Self::cyclotomic_mul_ith_limb::<1087254529>,
            Self::cyclotomic_mul_ith_limb::<1087303681>,
            Self::cyclotomic_mul_ith_limb::<1087418369>,
            Self::cyclotomic_mul_ith_limb::<1087451137>,
            Self::cyclotomic_mul_ith_limb::<1087475713>,
            Self::cyclotomic_mul_ith_limb::<1087492097>,
            Self::cyclotomic_mul_ith_limb::<1087762433>,
            Self::cyclotomic_mul_ith_limb::<1087836161>,
            Self::cyclotomic_mul_ith_limb::<1087860737>,
            Self::cyclotomic_mul_ith_limb::<1087885313>,
            Self::cyclotomic_mul_ith_limb::<1088311297>,
            Self::cyclotomic_mul_ith_limb::<1088376833>,
            Self::cyclotomic_mul_ith_limb::<1088401409>,
            Self::cyclotomic_mul_ith_limb::<1088409601>,
            Self::cyclotomic_mul_ith_limb::<1088458753>,
            Self::cyclotomic_mul_ith_limb::<1088606209>,
            Self::cyclotomic_mul_ith_limb::<1088647169>,
            Self::cyclotomic_mul_ith_limb::<1088679937>,
            Self::cyclotomic_mul_ith_limb::<1088778241>,
            Self::cyclotomic_mul_ith_limb::<1088802817>,
            Self::cyclotomic_mul_ith_limb::<1088851969>,
            Self::cyclotomic_mul_ith_limb::<1088868353>,
            Self::cyclotomic_mul_ith_limb::<1088966657>,
            Self::cyclotomic_mul_ith_limb::<1089220609>,
            Self::cyclotomic_mul_ith_limb::<1089458177>,
            Self::cyclotomic_mul_ith_limb::<1089540097>,
            Self::cyclotomic_mul_ith_limb::<1089835009>,
            Self::cyclotomic_mul_ith_limb::<1089949697>,
            Self::cyclotomic_mul_ith_limb::<1090170881>,
            Self::cyclotomic_mul_ith_limb::<1090179073>,
            Self::cyclotomic_mul_ith_limb::<1090203649>,
            Self::cyclotomic_mul_ith_limb::<1090293761>,
            Self::cyclotomic_mul_ith_limb::<1090400257>,
            Self::cyclotomic_mul_ith_limb::<1090490369>,
            Self::cyclotomic_mul_ith_limb::<1090498561>,
            Self::cyclotomic_mul_ith_limb::<1090646017>,
            Self::cyclotomic_mul_ith_limb::<1090768897>,
            Self::cyclotomic_mul_ith_limb::<1090867201>,
            Self::cyclotomic_mul_ith_limb::<1091014657>,
            Self::cyclotomic_mul_ith_limb::<1091104769>,
            Self::cyclotomic_mul_ith_limb::<1091178497>,
            Self::cyclotomic_mul_ith_limb::<1091399681>,
            Self::cyclotomic_mul_ith_limb::<1091571713>,
            Self::cyclotomic_mul_ith_limb::<1091850241>,
            Self::cyclotomic_mul_ith_limb::<1092022273>,
            Self::cyclotomic_mul_ith_limb::<1092038657>,
            Self::cyclotomic_mul_ith_limb::<1092210689>,
            Self::cyclotomic_mul_ith_limb::<1092268033>,
            Self::cyclotomic_mul_ith_limb::<1092333569>,
            Self::cyclotomic_mul_ith_limb::<1092636673>,
            Self::cyclotomic_mul_ith_limb::<1092653057>,
            Self::cyclotomic_mul_ith_limb::<1092661249>,
            Self::cyclotomic_mul_ith_limb::<1092734977>,
            Self::cyclotomic_mul_ith_limb::<1092882433>,
            Self::cyclotomic_mul_ith_limb::<1092923393>,
            Self::cyclotomic_mul_ith_limb::<1092997121>,
            Self::cyclotomic_mul_ith_limb::<1093005313>,
            Self::cyclotomic_mul_ith_limb::<1093152769>,
            Self::cyclotomic_mul_ith_limb::<1093292033>,
            Self::cyclotomic_mul_ith_limb::<1093414913>,
            Self::cyclotomic_mul_ith_limb::<1093439489>,
            Self::cyclotomic_mul_ith_limb::<1093513217>,
            Self::cyclotomic_mul_ith_limb::<1093537793>,
            Self::cyclotomic_mul_ith_limb::<1093570561>,
            Self::cyclotomic_mul_ith_limb::<1093939201>,
            Self::cyclotomic_mul_ith_limb::<1094184961>,
            Self::cyclotomic_mul_ith_limb::<1094250497>,
            Self::cyclotomic_mul_ith_limb::<1094258689>,
            Self::cyclotomic_mul_ith_limb::<1094299649>,
            Self::cyclotomic_mul_ith_limb::<1094356993>,
            Self::cyclotomic_mul_ith_limb::<1094381569>,
            Self::cyclotomic_mul_ith_limb::<1094430721>,
            Self::cyclotomic_mul_ith_limb::<1094471681>,
            Self::cyclotomic_mul_ith_limb::<1094504449>,
            Self::cyclotomic_mul_ith_limb::<1094725633>,
            Self::cyclotomic_mul_ith_limb::<1094889473>,
            Self::cyclotomic_mul_ith_limb::<1094914049>,
            Self::cyclotomic_mul_ith_limb::<1094946817>,
            Self::cyclotomic_mul_ith_limb::<1094963201>,
            Self::cyclotomic_mul_ith_limb::<1094995969>,
            Self::cyclotomic_mul_ith_limb::<1095012353>,
            Self::cyclotomic_mul_ith_limb::<1095045121>,
            Self::cyclotomic_mul_ith_limb::<1095331841>,
            Self::cyclotomic_mul_ith_limb::<1095462913>,
            Self::cyclotomic_mul_ith_limb::<1095536641>,
            Self::cyclotomic_mul_ith_limb::<1095585793>,
            Self::cyclotomic_mul_ith_limb::<1095700481>,
            Self::cyclotomic_mul_ith_limb::<1095708673>,
            Self::cyclotomic_mul_ith_limb::<1095774209>,
            Self::cyclotomic_mul_ith_limb::<1095847937>,
            Self::cyclotomic_mul_ith_limb::<1096101889>,
            Self::cyclotomic_mul_ith_limb::<1096142849>,
            Self::cyclotomic_mul_ith_limb::<1096175617>,
            Self::cyclotomic_mul_ith_limb::<1096216577>,
            Self::cyclotomic_mul_ith_limb::<1096224769>,
            Self::cyclotomic_mul_ith_limb::<1096339457>,
            Self::cyclotomic_mul_ith_limb::<1096388609>,
            Self::cyclotomic_mul_ith_limb::<1096396801>,
            Self::cyclotomic_mul_ith_limb::<1096486913>,
            Self::cyclotomic_mul_ith_limb::<1096544257>,
            Self::cyclotomic_mul_ith_limb::<1096609793>,
            Self::cyclotomic_mul_ith_limb::<1096716289>,
            Self::cyclotomic_mul_ith_limb::<1096765441>,
            Self::cyclotomic_mul_ith_limb::<1096880129>,
            Self::cyclotomic_mul_ith_limb::<1096937473>,
            Self::cyclotomic_mul_ith_limb::<1096953857>,
            Self::cyclotomic_mul_ith_limb::<1096962049>,
            Self::cyclotomic_mul_ith_limb::<1097076737>,
            Self::cyclotomic_mul_ith_limb::<1097175041>,
            Self::cyclotomic_mul_ith_limb::<1097183233>,
            Self::cyclotomic_mul_ith_limb::<1097199617>,
            Self::cyclotomic_mul_ith_limb::<1097256961>,
            Self::cyclotomic_mul_ith_limb::<1097306113>,
            Self::cyclotomic_mul_ith_limb::<1097322497>,
            Self::cyclotomic_mul_ith_limb::<1097453569>,
            Self::cyclotomic_mul_ith_limb::<1097469953>,
            Self::cyclotomic_mul_ith_limb::<1097527297>,
            Self::cyclotomic_mul_ith_limb::<1097691137>,
            Self::cyclotomic_mul_ith_limb::<1097715713>,
            Self::cyclotomic_mul_ith_limb::<1097748481>,
            Self::cyclotomic_mul_ith_limb::<1097895937>,
            Self::cyclotomic_mul_ith_limb::<1098035201>,
            Self::cyclotomic_mul_ith_limb::<1098084353>,
            Self::cyclotomic_mul_ith_limb::<1098231809>,
            Self::cyclotomic_mul_ith_limb::<1098280961>,
            Self::cyclotomic_mul_ith_limb::<1098313729>,
            Self::cyclotomic_mul_ith_limb::<1098387457>,
            Self::cyclotomic_mul_ith_limb::<1098403841>,
            Self::cyclotomic_mul_ith_limb::<1098485761>,
            Self::cyclotomic_mul_ith_limb::<1098731521>,
            Self::cyclotomic_mul_ith_limb::<1098756097>,
            Self::cyclotomic_mul_ith_limb::<1098780673>,
            Self::cyclotomic_mul_ith_limb::<1098805249>,
            Self::cyclotomic_mul_ith_limb::<1098854401>,
            Self::cyclotomic_mul_ith_limb::<1098919937>,
            Self::cyclotomic_mul_ith_limb::<1099042817>,
            Self::cyclotomic_mul_ith_limb::<1099067393>,
            Self::cyclotomic_mul_ith_limb::<1099272193>,
            Self::cyclotomic_mul_ith_limb::<1099296769>,
            Self::cyclotomic_mul_ith_limb::<1099370497>,
            Self::cyclotomic_mul_ith_limb::<1099411457>,
            Self::cyclotomic_mul_ith_limb::<1099591681>,
            Self::cyclotomic_mul_ith_limb::<1099665409>,
            Self::cyclotomic_mul_ith_limb::<1099755521>,
            Self::cyclotomic_mul_ith_limb::<1099763713>,
            Self::cyclotomic_mul_ith_limb::<1100001281>,
            Self::cyclotomic_mul_ith_limb::<1100132353>,
            Self::cyclotomic_mul_ith_limb::<1100173313>,
            Self::cyclotomic_mul_ith_limb::<1100296193>,
            Self::cyclotomic_mul_ith_limb::<1100451841>,
            Self::cyclotomic_mul_ith_limb::<1100492801>,
            Self::cyclotomic_mul_ith_limb::<1100500993>,
            Self::cyclotomic_mul_ith_limb::<1100574721>,
            Self::cyclotomic_mul_ith_limb::<1100623873>,
            Self::cyclotomic_mul_ith_limb::<1100689409>,
            Self::cyclotomic_mul_ith_limb::<1100697601>,
            Self::cyclotomic_mul_ith_limb::<1100812289>,
            Self::cyclotomic_mul_ith_limb::<1100820481>,
            Self::cyclotomic_mul_ith_limb::<1100861441>,
            Self::cyclotomic_mul_ith_limb::<1100894209>,
            Self::cyclotomic_mul_ith_limb::<1101107201>,
            Self::cyclotomic_mul_ith_limb::<1101139969>,
            Self::cyclotomic_mul_ith_limb::<1101189121>,
            Self::cyclotomic_mul_ith_limb::<1101213697>,
            Self::cyclotomic_mul_ith_limb::<1101352961>,
            Self::cyclotomic_mul_ith_limb::<1101361153>,
            Self::cyclotomic_mul_ith_limb::<1101385729>,
            Self::cyclotomic_mul_ith_limb::<1101434881>,
            Self::cyclotomic_mul_ith_limb::<1101549569>,
            Self::cyclotomic_mul_ith_limb::<1101557761>,
            Self::cyclotomic_mul_ith_limb::<1101598721>,
            Self::cyclotomic_mul_ith_limb::<1101623297>,
            Self::cyclotomic_mul_ith_limb::<1101672449>,
            Self::cyclotomic_mul_ith_limb::<1101729793>,
            Self::cyclotomic_mul_ith_limb::<1101852673>,
            Self::cyclotomic_mul_ith_limb::<1101926401>,
            Self::cyclotomic_mul_ith_limb::<1102041089>,
            Self::cyclotomic_mul_ith_limb::<1102163969>,
            Self::cyclotomic_mul_ith_limb::<1102295041>,
            Self::cyclotomic_mul_ith_limb::<1102360577>,
            Self::cyclotomic_mul_ith_limb::<1102483457>,
            Self::cyclotomic_mul_ith_limb::<1102532609>,
            Self::cyclotomic_mul_ith_limb::<1102565377>,
            Self::cyclotomic_mul_ith_limb::<1102712833>,
            Self::cyclotomic_mul_ith_limb::<1102860289>,
            Self::cyclotomic_mul_ith_limb::<1102934017>,
            Self::cyclotomic_mul_ith_limb::<1103147009>,
            Self::cyclotomic_mul_ith_limb::<1103196161>,
            Self::cyclotomic_mul_ith_limb::<1103245313>,
            Self::cyclotomic_mul_ith_limb::<1103343617>,
            Self::cyclotomic_mul_ith_limb::<1103548417>,
            Self::cyclotomic_mul_ith_limb::<1103835137>,
            Self::cyclotomic_mul_ith_limb::<1103843329>,
            Self::cyclotomic_mul_ith_limb::<1103917057>,
            Self::cyclotomic_mul_ith_limb::<1103966209>,
            Self::cyclotomic_mul_ith_limb::<1104007169>,
            Self::cyclotomic_mul_ith_limb::<1104261121>,
            Self::cyclotomic_mul_ith_limb::<1104408577>,
            Self::cyclotomic_mul_ith_limb::<1104457729>,
            Self::cyclotomic_mul_ith_limb::<1104654337>,
            Self::cyclotomic_mul_ith_limb::<1104695297>,
            Self::cyclotomic_mul_ith_limb::<1104703489>,
            Self::cyclotomic_mul_ith_limb::<1104719873>,
            Self::cyclotomic_mul_ith_limb::<1104867329>,
            Self::cyclotomic_mul_ith_limb::<1105195009>,
            Self::cyclotomic_mul_ith_limb::<1105309697>,
            Self::cyclotomic_mul_ith_limb::<1105367041>,
            Self::cyclotomic_mul_ith_limb::<1105408001>,
            Self::cyclotomic_mul_ith_limb::<1105457153>,
            Self::cyclotomic_mul_ith_limb::<1105489921>,
            Self::cyclotomic_mul_ith_limb::<1105514497>,
            Self::cyclotomic_mul_ith_limb::<1105555457>,
            Self::cyclotomic_mul_ith_limb::<1105580033>,
            Self::cyclotomic_mul_ith_limb::<1105604609>,
            Self::cyclotomic_mul_ith_limb::<1105686529>,
            Self::cyclotomic_mul_ith_limb::<1105784833>,
            Self::cyclotomic_mul_ith_limb::<1105825793>,
            Self::cyclotomic_mul_ith_limb::<1105883137>,
            Self::cyclotomic_mul_ith_limb::<1105907713>,
            Self::cyclotomic_mul_ith_limb::<1105981441>,
            Self::cyclotomic_mul_ith_limb::<1106030593>,
            Self::cyclotomic_mul_ith_limb::<1106071553>,
            Self::cyclotomic_mul_ith_limb::<1106227201>,
            Self::cyclotomic_mul_ith_limb::<1106300929>,
            Self::cyclotomic_mul_ith_limb::<1106350081>,
            Self::cyclotomic_mul_ith_limb::<1106374657>,
            Self::cyclotomic_mul_ith_limb::<1106399233>,
            Self::cyclotomic_mul_ith_limb::<1106464769>,
            Self::cyclotomic_mul_ith_limb::<1106513921>,
            Self::cyclotomic_mul_ith_limb::<1106538497>,
            Self::cyclotomic_mul_ith_limb::<1106546689>,
            Self::cyclotomic_mul_ith_limb::<1106710529>,
            Self::cyclotomic_mul_ith_limb::<1106931713>,
            Self::cyclotomic_mul_ith_limb::<1107030017>,
            Self::cyclotomic_mul_ith_limb::<1107054593>,
            Self::cyclotomic_mul_ith_limb::<1107521537>,
            Self::cyclotomic_mul_ith_limb::<1107619841>,
            Self::cyclotomic_mul_ith_limb::<1107668993>,
            Self::cyclotomic_mul_ith_limb::<1107816449>,
            Self::cyclotomic_mul_ith_limb::<1108013057>,
            Self::cyclotomic_mul_ith_limb::<1108021249>,
            Self::cyclotomic_mul_ith_limb::<1108062209>,
            Self::cyclotomic_mul_ith_limb::<1108135937>,
            Self::cyclotomic_mul_ith_limb::<1108144129>,
            Self::cyclotomic_mul_ith_limb::<1108217857>,
            Self::cyclotomic_mul_ith_limb::<1108381697>,
            Self::cyclotomic_mul_ith_limb::<1108430849>,
            Self::cyclotomic_mul_ith_limb::<1108439041>,
            Self::cyclotomic_mul_ith_limb::<1108463617>,
            Self::cyclotomic_mul_ith_limb::<1108488193>,
            Self::cyclotomic_mul_ith_limb::<1108553729>,
            Self::cyclotomic_mul_ith_limb::<1108586497>,
            Self::cyclotomic_mul_ith_limb::<1108611073>,
            Self::cyclotomic_mul_ith_limb::<1108750337>,
            Self::cyclotomic_mul_ith_limb::<1108774913>,
            Self::cyclotomic_mul_ith_limb::<1108922369>,
            Self::cyclotomic_mul_ith_limb::<1108979713>,
            Self::cyclotomic_mul_ith_limb::<1109004289>,
            Self::cyclotomic_mul_ith_limb::<1109094401>,
            Self::cyclotomic_mul_ith_limb::<1109127169>,
            Self::cyclotomic_mul_ith_limb::<1109217281>,
            Self::cyclotomic_mul_ith_limb::<1109299201>,
            Self::cyclotomic_mul_ith_limb::<1109323777>,
            Self::cyclotomic_mul_ith_limb::<1109618689>,
            Self::cyclotomic_mul_ith_limb::<1109692417>,
            Self::cyclotomic_mul_ith_limb::<1109708801>,
            Self::cyclotomic_mul_ith_limb::<1109880833>,
            Self::cyclotomic_mul_ith_limb::<1109962753>,
            Self::cyclotomic_mul_ith_limb::<1110085633>,
            Self::cyclotomic_mul_ith_limb::<1110151169>,
            Self::cyclotomic_mul_ith_limb::<1110224897>,
            Self::cyclotomic_mul_ith_limb::<1110249473>,
            Self::cyclotomic_mul_ith_limb::<1110282241>,
            Self::cyclotomic_mul_ith_limb::<1110306817>,
            Self::cyclotomic_mul_ith_limb::<1110454273>,
            Self::cyclotomic_mul_ith_limb::<1110716417>,
            Self::cyclotomic_mul_ith_limb::<1111109633>,
            Self::cyclotomic_mul_ith_limb::<1111166977>,
            Self::cyclotomic_mul_ith_limb::<1111183361>,
            Self::cyclotomic_mul_ith_limb::<1111216129>,
            Self::cyclotomic_mul_ith_limb::<1111232513>,
            Self::cyclotomic_mul_ith_limb::<1111257089>,
            Self::cyclotomic_mul_ith_limb::<1111412737>,
            Self::cyclotomic_mul_ith_limb::<1111461889>,
            Self::cyclotomic_mul_ith_limb::<1111502849>,
            Self::cyclotomic_mul_ith_limb::<1111560193>,
            Self::cyclotomic_mul_ith_limb::<1111584769>,
            Self::cyclotomic_mul_ith_limb::<1111625729>,
            Self::cyclotomic_mul_ith_limb::<1111658497>,
            Self::cyclotomic_mul_ith_limb::<1111683073>,
            Self::cyclotomic_mul_ith_limb::<1111756801>,
            Self::cyclotomic_mul_ith_limb::<1111797761>,
            Self::cyclotomic_mul_ith_limb::<1111822337>,
            Self::cyclotomic_mul_ith_limb::<1111928833>,
            Self::cyclotomic_mul_ith_limb::<1111994369>,
            Self::cyclotomic_mul_ith_limb::<1112043521>,
            Self::cyclotomic_mul_ith_limb::<1112174593>,
            Self::cyclotomic_mul_ith_limb::<1112289281>,
            Self::cyclotomic_mul_ith_limb::<1112338433>,
            Self::cyclotomic_mul_ith_limb::<1112363009>,
            Self::cyclotomic_mul_ith_limb::<1112371201>,
            Self::cyclotomic_mul_ith_limb::<1112494081>,
            Self::cyclotomic_mul_ith_limb::<1112567809>,
            Self::cyclotomic_mul_ith_limb::<1112682497>,
            Self::cyclotomic_mul_ith_limb::<1112739841>,
            Self::cyclotomic_mul_ith_limb::<1112887297>,
            Self::cyclotomic_mul_ith_limb::<1112952833>,
            Self::cyclotomic_mul_ith_limb::<1112977409>,
            Self::cyclotomic_mul_ith_limb::<1113059329>,
            Self::cyclotomic_mul_ith_limb::<1113346049>,
            Self::cyclotomic_mul_ith_limb::<1113403393>,
            Self::cyclotomic_mul_ith_limb::<1113567233>,
            Self::cyclotomic_mul_ith_limb::<1113600001>,
            Self::cyclotomic_mul_ith_limb::<1113747457>,
            Self::cyclotomic_mul_ith_limb::<1113886721>,
            Self::cyclotomic_mul_ith_limb::<1113894913>,
            Self::cyclotomic_mul_ith_limb::<1113935873>,
            Self::cyclotomic_mul_ith_limb::<1114034177>,
            Self::cyclotomic_mul_ith_limb::<1114165249>,
            Self::cyclotomic_mul_ith_limb::<1114238977>,
            Self::cyclotomic_mul_ith_limb::<1114451969>,
            Self::cyclotomic_mul_ith_limb::<1114632193>,
            Self::cyclotomic_mul_ith_limb::<1114746881>,
            Self::cyclotomic_mul_ith_limb::<1114771457>,
            Self::cyclotomic_mul_ith_limb::<1114894337>,
            Self::cyclotomic_mul_ith_limb::<1115025409>,
            Self::cyclotomic_mul_ith_limb::<1115074561>,
            Self::cyclotomic_mul_ith_limb::<1115148289>,
            Self::cyclotomic_mul_ith_limb::<1115410433>,
            Self::cyclotomic_mul_ith_limb::<1115492353>,
            Self::cyclotomic_mul_ith_limb::<1115590657>,
            Self::cyclotomic_mul_ith_limb::<1115615233>,
            Self::cyclotomic_mul_ith_limb::<1115631617>,
            Self::cyclotomic_mul_ith_limb::<1115713537>,
            Self::cyclotomic_mul_ith_limb::<1115729921>,
            Self::cyclotomic_mul_ith_limb::<1115762689>,
            Self::cyclotomic_mul_ith_limb::<1115803649>,
            Self::cyclotomic_mul_ith_limb::<1115901953>,
            Self::cyclotomic_mul_ith_limb::<1116131329>,
            Self::cyclotomic_mul_ith_limb::<1116270593>,
            Self::cyclotomic_mul_ith_limb::<1116418049>,
            Self::cyclotomic_mul_ith_limb::<1116549121>,
            Self::cyclotomic_mul_ith_limb::<1116639233>,
            Self::cyclotomic_mul_ith_limb::<1116794881>,
            Self::cyclotomic_mul_ith_limb::<1116917761>,
            Self::cyclotomic_mul_ith_limb::<1116991489>,
        ];

        for (i, function) in cyclotomic_mul_ith_limbs_p
            .iter()
            .enumerate()
            .take(num_limbs)
        {
            function(lhs, rhs, &mut result, i);
        }

        result
    }
}

impl From<Polynomial<MultiModInt>> for Polynomial<BigInt> {
    fn from(value: Polynomial<MultiModInt>) -> Self {
        Polynomial::new(
            value
                .coefficients
                .iter()
                .map(|mmi| BigInt::from(mmi.clone()))
                .collect_vec(),
        )
    }
}

impl From<Polynomial<MultiModInt>> for Polynomial<BigUint> {
    fn from(value: Polynomial<MultiModInt>) -> Self {
        Polynomial::new(
            value
                .coefficients
                .iter()
                .map(|mmi| BigUint::from(mmi.clone()))
                .collect_vec(),
        )
    }
}

#[derive(Debug, Clone)]
pub struct MultiModIntPolynomialConversionError<T: Display + Clone> {
    original: Polynomial<T>,
}

impl<T: Display + Clone> Display for MultiModIntPolynomialConversionError<T>
where
    MultiModInt: TryFrom<T>,
    <MultiModInt as TryFrom<T>>::Error: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "MultiModIntPolynomialConversionError({})",
            self.original
                .coefficients
                .iter()
                .cloned()
                .map(|c| MultiModInt::try_from(c))
                .map(|r| match r {
                    Ok(_) => "ok".to_string(),
                    Err(e) => format!("{}", e),
                })
                .join(", ")
        ))
    }
}

impl TryFrom<Polynomial<BigUint>> for Polynomial<MultiModInt> {
    type Error = MultiModIntPolynomialConversionError<BigUint>;
    fn try_from(val: Polynomial<BigUint>) -> Result<Self, Self::Error> {
        let coefficients = val.coefficients.iter().cloned().map(MultiModInt::try_from);
        if coefficients.clone().all(|c| c.is_ok()) {
            Ok(Polynomial::new(
                coefficients.map(|c| c.unwrap()).collect_vec(),
            ))
        } else {
            Err(MultiModIntPolynomialConversionError { original: val })
        }
    }
}

impl TryFrom<Polynomial<BigInt>> for Polynomial<MultiModInt> {
    type Error = MultiModIntPolynomialConversionError<BigInt>;
    fn try_from(val: Polynomial<BigInt>) -> Result<Self, Self::Error> {
        let coefficients = val.coefficients.iter().cloned().map(MultiModInt::try_from);
        if coefficients.clone().all(|c| c.is_ok()) {
            Ok(Polynomial::new(
                coefficients.map(|c| c.unwrap()).collect_vec(),
            ))
        } else {
            Err(MultiModIntPolynomialConversionError { original: val })
        }
    }
}

#[cfg(test)]
mod test {

    use std::str::FromStr;

    use itertools::Itertools;
    use num::{BigInt, BigUint, Integer, One};
    use proptest::arbitrary::Arbitrary;
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest::proptest;
    use proptest::strategy::BoxedStrategy;
    use proptest::strategy::Just;
    use proptest::strategy::Strategy;
    use test_strategy::proptest as strategy_proptest;

    use crate::polynomial::Polynomial;
    use crate::residue_number_system::BEZOUT_COEFFICIENTS;

    use super::{MultiModInt, MODULI};

    fn arbitrary_bigint(bitlen: usize) -> BoxedStrategy<BigInt> {
        let limbs = vec(u32::arbitrary(), bitlen.div_ceil(32));
        let sign = bool::arbitrary();
        (limbs, sign)
            .prop_map(move |(bb, ss)| {
                let bigint = BigInt::from(BigUint::from_slice(&bb) >> (bb.len() * 32 - bitlen));
                if ss {
                    -bigint
                } else {
                    bigint
                }
            })
            .boxed()
    }

    fn arbitrary_multimodint(bitlen: usize) -> BoxedStrategy<MultiModInt> {
        arbitrary_bigint(bitlen)
            .prop_map(|bi| MultiModInt::try_from(bi).unwrap())
            .boxed()
    }

    fn arbitrary_multimodint_polynomial(
        bitlen: usize,
        num_coefficients: usize,
    ) -> BoxedStrategy<Polynomial<MultiModInt>> {
        vec(arbitrary_multimodint(bitlen), num_coefficients)
            .prop_map(Polynomial::new)
            .boxed()
    }

    fn arbitrary_bigint_polynomial(
        bitlen: usize,
        num_coefficients: usize,
    ) -> BoxedStrategy<Polynomial<BigInt>> {
        vec(arbitrary_bigint(bitlen), num_coefficients)
            .prop_map(Polynomial::new)
            .boxed()
    }

    #[test]
    fn coefficients_are_1_mod_p() {
        for (i, p, c) in MODULI
            .into_iter()
            .enumerate()
            .map(|(i, p)| (i, p, BEZOUT_COEFFICIENTS.get(i).unwrap()))
        {
            assert_eq!(
                c.mod_floor(&Into::<BigUint>::into(p)),
                BigUint::one(),
                "mismatch in position {i}\nmodulus is {p}\ncoefficient is {}",
                c
            );
        }
    }

    #[strategy_proptest]
    fn casting_preserves_bitsize(
        #[strategy(1usize..1000)] bits: usize,
        #[strategy(arbitrary_bigint(#bits))] bigint: BigInt,
    ) {
        let mmi = MultiModInt::try_from(bigint).unwrap();

        assert!(
            BigInt::from(mmi.clone()).bits() <= bits as u64,
            "{} > {}",
            BigInt::from(mmi).bits(),
            bits
        );
    }

    #[strategy_proptest]
    fn multimodint_arithmetic(
        #[strategy(arbitrary_bigint(1000))] a: BigInt,
        #[strategy(arbitrary_bigint(1000))] b: BigInt,
        #[strategy(arbitrary_bigint(1000))] c: BigInt,
    ) {
        let d = a.clone() * b.clone() - c.clone();
        let d_mmi = MultiModInt::try_from(d).unwrap();
        let mmi_d = MultiModInt::try_from(a).unwrap() * MultiModInt::try_from(b).unwrap()
            - MultiModInt::try_from(c).unwrap();
        prop_assert_eq!(d_mmi, mmi_d);
    }

    proptest! {
        #[test]
        fn to_and_fro_proptest(a in "[+-]?[0-9]{1,2980}") {
            let b = BigInt::from_str(&a).unwrap();
            let c = MultiModInt::try_from(b.clone()).unwrap();
            let d = BigInt::from(c.clone());
            let e = MultiModInt::try_from(d.clone()).unwrap();
            prop_assert_eq!(b, d);
            prop_assert_eq!(c, e);
        }
    }

    #[test]
    fn to_and_fro_fixed() {
        let b = BigInt::from_str("-2098365342651410771540682191176265762931285").unwrap();
        let c = MultiModInt::try_from(b.clone()).unwrap();
        let d = BigInt::from(c);
        assert_eq!(b, d);
    }

    #[test]
    fn multimodint_fixed() {
        let a = BigInt::from_str("-322777953413029095759719085619503741066179").unwrap();
        println!("a has {} bits", a.bits());
        let b = MultiModInt::try_from(a.clone()).unwrap();
        println!("b has {} limbs", b.limbs.len());
        let c = BigInt::from(b);
        assert_eq!(c, a);
    }

    #[strategy_proptest]
    fn multimodint_polynomial_arithmetic_a(
        #[strategy(0usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(Just(70))] _bitlen: usize,
        #[strategy(arbitrary_bigint_polynomial(#_bitlen, #_n))] a: Polynomial<BigInt>,
        #[strategy(arbitrary_bigint_polynomial(#_bitlen, #_n))] b: Polynomial<BigInt>,
    ) {
        let mut bigint_terms = vec![];
        let mut multimod_terms = vec![];
        for (ai, ac) in a.coefficients.iter().enumerate() {
            for (bi, bc) in b.coefficients.iter().enumerate() {
                if ai + bi == 11 {
                    bigint_terms.push(ac.clone() * bc);
                    multimod_terms.push(
                        MultiModInt::try_from(ac.clone()).unwrap()
                            * MultiModInt::try_from(bc.clone()).unwrap(),
                    );
                }
            }
        }
        prop_assert_eq!(
            bigint_terms
                .iter()
                .map(|bi| MultiModInt::try_from(bi.clone()).unwrap())
                .collect_vec(),
            multimod_terms.clone()
        );
        let bigint_sum = bigint_terms.into_iter().sum::<BigInt>();
        let multimod_sum = multimod_terms.into_iter().sum();
        prop_assert_eq!(MultiModInt::try_from(bigint_sum).unwrap(), multimod_sum);
        let d = a.clone() * b.clone();
        let mmi_d = Polynomial::<MultiModInt>::try_from(d.clone()).unwrap();

        let d_mmi = Polynomial::<MultiModInt>::try_from(a.clone()).unwrap()
            * Polynomial::<MultiModInt>::try_from(b.clone()).unwrap();

        assert_eq!(mmi_d, d_mmi,);
    }

    #[strategy_proptest]
    fn multimodint_polynomial_arithmetic_b(
        #[strategy(0usize..4)] _logn: usize,
        #[strategy(Just(1 << #_logn))] _n: usize,
        #[strategy(Just(70))] _bitlen: usize,
        #[strategy(arbitrary_bigint_polynomial(#_bitlen, #_n))] a: Polynomial<BigInt>,
        #[strategy(arbitrary_bigint_polynomial(#_bitlen, #_n))] b: Polynomial<BigInt>,
        #[strategy(arbitrary_bigint_polynomial(#_bitlen, #_n))] c: Polynomial<BigInt>,
    ) {
        let d = (a.clone() * b.clone()) - c.clone();
        let mmi_d = Polynomial::<MultiModInt>::try_from(d.clone()).unwrap();

        let d_mmi = Polynomial::<MultiModInt>::try_from(a.clone()).unwrap()
            * Polynomial::<MultiModInt>::try_from(b.clone()).unwrap()
            - Polynomial::<MultiModInt>::try_from(c).unwrap();

        prop_assert_eq!(mmi_d, d_mmi);
    }

    #[strategy_proptest]
    fn cyclotomic_multiplication(
        #[strategy(0usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] n: usize,
        #[strategy(0usize..567)] _bitlen: usize,
        #[strategy(arbitrary_multimodint_polynomial(#_bitlen, #n))] a: Polynomial<MultiModInt>,
        #[strategy(arbitrary_multimodint_polynomial(#_bitlen, #n))] b: Polynomial<MultiModInt>,
    ) {
        let product = a.clone() * b.clone();
        let c_trad = product.reduce_by_cyclotomic(n);

        let c_fast = a.cyclotomic_mul(&b);
        assert_eq!(c_trad, c_fast);
    }

    const MULTIMOD_CAPACITY: usize = 1001;
    #[strategy_proptest(cases = 50)]
    fn shift_left(
        #[filter(|x| *x < MULTIMOD_CAPACITY)]
        #[strategy(1usize..(MULTIMOD_CAPACITY-1))]
        _total_shift: usize,
        #[strategy(0usize..#_total_shift)] shift1: usize,
        #[strategy(Just(#_total_shift - #shift1))] shift2: usize,
        #[strategy(0usize..(MULTIMOD_CAPACITY-(#shift1+#shift2)))] _bitlen: usize,
        #[strategy(0usize..=0)] _logn: usize,
        #[strategy(arbitrary_multimodint_polynomial(#_bitlen, 1 << #_logn))] a: Polynomial<
            MultiModInt,
        >,
    ) {
        let b = a.map(|mmi| mmi.clone() << shift1);
        let c = b.map(|mmi| mmi.clone() << shift2);

        let d = a.map(|mmi| mmi.clone() << (shift1 + shift2));

        prop_assert_eq!(c, d);
    }

    #[strategy_proptest(cases = 10)]
    fn conversion_between_biguint_and_multimodint(
        #[strategy(vec(0..u32::MAX, 6))] limbs: Vec<u32>,
    ) {
        let biguint = BigUint::from_slice(&limbs);
        let multimodint = MultiModInt::try_from(biguint.clone()).unwrap();
        let biguint_again = BigUint::from(multimodint);
        prop_assert_eq!(biguint, biguint_again);
    }
}
