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

use crate::{cyclotomic_fourier::CyclotomicFourier, inverse::Inverse, polynomial::Polynomial};

const N: usize = 330;
pub(crate) const MULTIMOD_MAX_CAPACITY: usize = N * 30;

/// The vector of moduli by which we will be reducing every limb.
/// These were chosen as p where p-1 = 2^12 * cofactor, where
/// cofactor is some 18-bit integer that makes p prime. In fact,
/// this list corresponds to the N smallest such cofactors.
pub const MODULI: [u32; N] = [
    1073754113u32,
    1073950721u32,
    1073958913u32,
    1073983489u32,
    1074196481u32,
    1074343937u32,
    1074442241u32,
    1074475009u32,
    1074515969u32,
    1074524161u32,
    1074548737u32,
    1074589697u32,
    1074597889u32,
    1074688001u32,
    1074696193u32,
    1074810881u32,
    1074835457u32,
    1074941953u32,
    1075007489u32,
    1075064833u32,
    1075105793u32,
    1075351553u32,
    1075376129u32,
    1075449857u32,
    1075507201u32,
    1075621889u32,
    1075695617u32,
    1075720193u32,
    1075752961u32,
    1076039681u32,
    1076064257u32,
    1076162561u32,
    1076187137u32,
    1076211713u32,
    1076269057u32,
    1076334593u32,
    1076391937u32,
    1076531201u32,
    1076613121u32,
    1076908033u32,
    1076948993u32,
    1077047297u32,
    1077178369u32,
    1077252097u32,
    1077391361u32,
    1077399553u32,
    1077620737u32,
    1077686273u32,
    1077792769u32,
    1077882881u32,
    1078038529u32,
    1078079489u32,
    1078112257u32,
    1078153217u32,
    1078210561u32,
    1078333441u32,
    1078382593u32,
    1078398977u32,
    1078497281u32,
    1078505473u32,
    1078546433u32,
    1078571009u32,
    1078620161u32,
    1078775809u32,
    1078849537u32,
    1078939649u32,
    1079185409u32,
    1079193601u32,
    1079357441u32,
    1079431169u32,
    1079603201u32,
    1079627777u32,
    1079635969u32,
    1079676929u32,
    1079685121u32,
    1079734273u32,
    1079832577u32,
    1079848961u32,
    1079873537u32,
    1079922689u32,
    1080020993u32,
    1080348673u32,
    1080365057u32,
    1080446977u32,
    1080496129u32,
    1080512513u32,
    1080741889u32,
    1080864769u32,
    1080938497u32,
    1080954881u32,
    1081077761u32,
    1081225217u32,
    1081323521u32,
    1081348097u32,
    1081397249u32,
    1081430017u32,
    1081495553u32,
    1081577473u32,
    1081643009u32,
    1081700353u32,
    1081774081u32,
    1081815041u32,
    1081839617u32,
    1081921537u32,
    1082060801u32,
    1082093569u32,
    1082331137u32,
    1082462209u32,
    1082552321u32,
    1082601473u32,
    1082699777u32,
    1082724353u32,
    1082806273u32,
    1082904577u32,
    1082929153u32,
    1083002881u32,
    1083076609u32,
    1083117569u32,
    1083215873u32,
    1083289601u32,
    1083371521u32,
    1083518977u32,
    1083535361u32,
    1083928577u32,
    1083953153u32,
    1084076033u32,
    1084133377u32,
    1084149761u32,
    1084297217u32,
    1084477441u32,
    1084518401u32,
    1084641281u32,
    1084649473u32,
    1084674049u32,
    1084690433u32,
    1084813313u32,
    1084911617u32,
    1084968961u32,
    1085140993u32,
    1085181953u32,
    1085255681u32,
    1085411329u32,
    1085501441u32,
    1085534209u32,
    1085550593u32,
    1085607937u32,
    1085648897u32,
    1085673473u32,
    1085779969u32,
    1085796353u32,
    1085870081u32,
    1085894657u32,
    1085952001u32,
    1086066689u32,
    1086115841u32,
    1086345217u32,
    1086410753u32,
    1086443521u32,
    1086566401u32,
    1086713857u32,
    1086763009u32,
    1086853121u32,
    1086902273u32,
    1086935041u32,
    1087025153u32,
    1087148033u32,
    1087254529u32,
    1087303681u32,
    1087418369u32,
    1087451137u32,
    1087475713u32,
    1087492097u32,
    1087762433u32,
    1087836161u32,
    1087860737u32,
    1087885313u32,
    1088311297u32,
    1088376833u32,
    1088401409u32,
    1088409601u32,
    1088458753u32,
    1088606209u32,
    1088647169u32,
    1088679937u32,
    1088778241u32,
    1088802817u32,
    1088851969u32,
    1088868353u32,
    1088966657u32,
    1089220609u32,
    1089458177u32,
    1089540097u32,
    1089835009u32,
    1089949697u32,
    1090170881u32,
    1090179073u32,
    1090203649u32,
    1090293761u32,
    1090400257u32,
    1090490369u32,
    1090498561u32,
    1090646017u32,
    1090768897u32,
    1090867201u32,
    1091014657u32,
    1091104769u32,
    1091178497u32,
    1091399681u32,
    1091571713u32,
    1091850241u32,
    1092022273u32,
    1092038657u32,
    1092210689u32,
    1092268033u32,
    1092333569u32,
    1092636673u32,
    1092653057u32,
    1092661249u32,
    1092734977u32,
    1092882433u32,
    1092923393u32,
    1092997121u32,
    1093005313u32,
    1093152769u32,
    1093292033u32,
    1093414913u32,
    1093439489u32,
    1093513217u32,
    1093537793u32,
    1093570561u32,
    1093939201u32,
    1094184961u32,
    1094250497u32,
    1094258689u32,
    1094299649u32,
    1094356993u32,
    1094381569u32,
    1094430721u32,
    1094471681u32,
    1094504449u32,
    1094725633u32,
    1094889473u32,
    1094914049u32,
    1094946817u32,
    1094963201u32,
    1094995969u32,
    1095012353u32,
    1095045121u32,
    1095331841u32,
    1095462913u32,
    1095536641u32,
    1095585793u32,
    1095700481u32,
    1095708673u32,
    1095774209u32,
    1095847937u32,
    1096101889u32,
    1096142849u32,
    1096175617u32,
    1096216577u32,
    1096224769u32,
    1096339457u32,
    1096388609u32,
    1096396801u32,
    1096486913u32,
    1096544257u32,
    1096609793u32,
    1096716289u32,
    1096765441u32,
    1096880129u32,
    1096937473u32,
    1096953857u32,
    1096962049u32,
    1097076737u32,
    1097175041u32,
    1097183233u32,
    1097199617u32,
    1097256961u32,
    1097306113u32,
    1097322497u32,
    1097453569u32,
    1097469953u32,
    1097527297u32,
    1097691137u32,
    1097715713u32,
    1097748481u32,
    1097895937u32,
    1098035201u32,
    1098084353u32,
    1098231809u32,
    1098280961u32,
    1098313729u32,
    1098387457u32,
    1098403841u32,
    1098485761u32,
    1098731521u32,
    1098756097u32,
    1098780673u32,
    1098805249u32,
    1098854401u32,
    1098919937u32,
    1099042817u32,
    1099067393u32,
    1099272193u32,
    1099296769u32,
    1099370497u32,
    1099411457u32,
    1099591681u32,
    1099665409u32,
    1099755521u32,
    1099763713u32,
    1100001281u32,
    1100132353u32,
    1100173313u32,
    1100296193u32,
    1100451841u32,
    1100492801u32,
    1100500993u32,
    1100574721u32,
    1100623873u32,
    1100689409u32,
    1100697601u32,
    1100812289u32,
    1100820481u32,
    1100861441u32,
    1100894209u32,
    1101107201u32,
    1101139969u32,
    1101189121u32,
    1101213697u32,
];

/// Coefficient c_i in the relation
/// integer = x_0 * c_0  +  x_1 * c_1  +  ...  + x_n-1 * c_n-1  mod prod_i p_i
/// where {x_i} are the integer's limbs.
///
/// The ith coefficient reduces to 1 modulo p_i, and to 0 modulo any other modulus.
/// It can be found as
///
/// c_i = [(prod_j p_j) / p_i]^-1  *  (prod_j p_j) / p_i
///
/// where the inverse is taken modulo p_i.
///
/// This function is not const because of dependencies, but we would like it to be.
fn coefficient(i: usize, n: usize) -> BigUint {
    // aggregate factors
    let mut reduced = 1u32;
    let mut aggregate_big = BigUint::one();
    for (j, &p) in MODULI.iter().enumerate().take(n) {
        if j != i {
            aggregate_big *= BigUint::from_u32(p).unwrap();
            reduced = (((reduced as u64) * (p as u64)) % (MODULI[i] as u64)) as u32;
        }
    }

    // invert reduced
    let mut inverse = 1u32;
    for j in (0..32).rev() {
        inverse = (((inverse as u64) * (inverse as u64)) % (MODULI[i] as u64)) as u32;
        if (MODULI[i] - 2) & (1 << j) != 0 {
            inverse = (((inverse as u64) * (reduced as u64)) % (MODULI[i] as u64)) as u32;
        }
    }

    // combine factors and return
    let inverse_big = BigUint::from_u32(inverse).unwrap();
    inverse_big * aggregate_big
}

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
/// Every limb stores at least 30 bits, so the total capacity is 9900 bits.
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
            num <= 330,
            "MultiModInt does not support more than 330 limbs (given: {num})."
        );
        let mut limbs = Vec::<u32>::with_capacity(num);
        let bits = int.bits() as usize;
        for p in MODULI.into_iter().take(num) {
            limbs.push(
                int.mod_floor(&BigUint::from_u32(p).unwrap())
                    .to_u32()
                    .unwrap(),
            );
        }
        MultiModInt {
            limbs,
            bitsize_bound: bits,
        }
    }

    fn get_more_limbs(&self, total_count: usize) -> Vec<u32> {
        assert!(total_count <= N);
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
        let num = usize::min(N, Self::num_limbs_needed(bits));
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
        let num = usize::min(N, Self::num_limbs_needed(bits));
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
            Self::num_limbs_needed(bits) <= N,
            "Not enough moduli available for {} bits",
            bits
        );
        let num = usize::min(N, Self::num_limbs_needed(bits));
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
            let num = usize::min(N, usize::max(num_limbs_needed, 2 * self.limbs.len()));
            assert!(
                num >= num_limbs_needed,
                "requested bit width ({}) exceeds capacity ({})",
                bits,
                N * 30,
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

impl From<BigUint> for MultiModInt {
    fn from(val: BigUint) -> Self {
        let bits = val.bits() as usize;
        let num = Self::num_limbs_needed(bits);
        MultiModInt::from_biguint_with_num_limbs(val, num)
    }
}

impl From<BigInt> for MultiModInt {
    fn from(value: BigInt) -> Self {
        if value < BigInt::zero() {
            -MultiModInt::from((-value).to_biguint().unwrap())
        } else {
            MultiModInt::from(value.to_biguint().unwrap())
        }
    }
}

impl From<MultiModInt> for BigUint {
    fn from(value: MultiModInt) -> BigUint {
        let num = value.limbs.len();
        value
            .limbs
            .iter()
            .enumerate()
            .map(|(i, x)| x * coefficient(i, num))
            .sum::<BigUint>()
            .mod_floor(&product(num))
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
            .map(|(i, x)| x * coefficient(i, num))
            .sum::<BigUint>()
            .mod_floor(&product);
        if biguint > (product.clone() >> 1) {
            biguint.to_bigint().unwrap() - product.to_bigint().unwrap()
        } else {
            biguint.to_bigint().unwrap()
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct PrimeField<const P: u32>(u32);

impl<const P: u32> Add for PrimeField<P> {
    type Output = PrimeField<P>;

    fn add(self, rhs: Self) -> Self::Output {
        PrimeField::<P>((((self.0 as u64) + (rhs.0 as u64)) % (P as u64)) as u32)
    }
}

impl<const P: u32> Mul for PrimeField<P> {
    type Output = PrimeField<P>;

    fn mul(self, rhs: Self) -> Self::Output {
        PrimeField::<P>((((self.0 as u64) * (rhs.0 as u64)) % (P as u64)) as u32)
    }
}

impl<const P: u32> Sub for PrimeField<P> {
    type Output = PrimeField<P>;

    fn sub(self, rhs: Self) -> Self::Output {
        PrimeField::<P>((((self.0 as u64) + (P as u64) - (rhs.0 as u64)) % (P as u64)) as u32)
    }
}

impl<const P: u32> AddAssign for PrimeField<P> {
    fn add_assign(&mut self, rhs: Self) {
        self.0 = (((self.0 as u64) + (rhs.0 as u64)) % (P as u64)) as u32;
    }
}

impl<const P: u32> SubAssign for PrimeField<P> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 = (((self.0 as u64) + ((P as u64) - (rhs.0 as u64))) % (P as u64)) as u32;
    }
}

impl<const P: u32> MulAssign for PrimeField<P> {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 = (((self.0 as u64) * (rhs.0 as u64)) % (P as u64)) as u32;
    }
}

impl<const P: u32> One for PrimeField<P> {
    fn one() -> Self {
        Self(1u32)
    }
}

impl<const P: u32> Zero for PrimeField<P> {
    fn zero() -> Self {
        Self(0u32)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0u32
    }
}

impl<const P: u32> Inverse for PrimeField<P> {
    fn inverse_or_zero(self) -> Self {
        let mut acc = Self::one();
        for i in (0..32).rev() {
            acc *= acc;
            if (1 << i) & (P - 2) != 0 {
                acc *= self;
            }
        }
        acc
    }
}

impl<const P: u32> CyclotomicFourier for PrimeField<P> {
    fn primitive_root_of_unity(n: usize) -> Self {
        let mut root = Self(root_of_unity_4096(P));
        for _ in n.ilog2()..12 {
            root = root * root;
        }
        root
    }
}

fn root_of_unity_4096(p: u32) -> u32 {
    match p {
        1073754113u32 => 872548469u32,
        1073950721u32 => 820856268u32,
        1073958913u32 => 217067019u32,
        1073983489u32 => 1028957174u32,
        1074196481u32 => 3359408u32,
        1074343937u32 => 295162511u32,
        1074442241u32 => 667034521u32,
        1074475009u32 => 847030914u32,
        1074515969u32 => 916967690u32,
        1074524161u32 => 860257261u32,
        1074548737u32 => 147401299u32,
        1074589697u32 => 104582941u32,
        1074597889u32 => 423456772u32,
        1074688001u32 => 547674691u32,
        1074696193u32 => 775215302u32,
        1074810881u32 => 288166217u32,
        1074835457u32 => 16045962u32,
        1074941953u32 => 14501489u32,
        1075007489u32 => 437673404u32,
        1075064833u32 => 568459387u32,
        1075105793u32 => 995532132u32,
        1075351553u32 => 604979534u32,
        1075376129u32 => 895473478u32,
        1075449857u32 => 334404174u32,
        1075507201u32 => 1019045579u32,
        1075621889u32 => 773225363u32,
        1075695617u32 => 203157826u32,
        1075720193u32 => 542463786u32,
        1075752961u32 => 675625838u32,
        1076039681u32 => 465105033u32,
        1076064257u32 => 288735288u32,
        1076162561u32 => 963146029u32,
        1076187137u32 => 293189229u32,
        1076211713u32 => 528452662u32,
        1076269057u32 => 188392172u32,
        1076334593u32 => 875372414u32,
        1076391937u32 => 614525675u32,
        1076531201u32 => 435656755u32,
        1076613121u32 => 1073414130u32,
        1076908033u32 => 594633712u32,
        1076948993u32 => 1074876355u32,
        1077047297u32 => 500612139u32,
        1077178369u32 => 673610734u32,
        1077252097u32 => 403145659u32,
        1077391361u32 => 913856820u32,
        1077399553u32 => 297993061u32,
        1077620737u32 => 460650866u32,
        1077686273u32 => 892567372u32,
        1077792769u32 => 671889178u32,
        1077882881u32 => 807326881u32,
        1078038529u32 => 1067485582u32,
        1078079489u32 => 350770841u32,
        1078112257u32 => 760387925u32,
        1078153217u32 => 4000735u32,
        1078210561u32 => 752107495u32,
        1078333441u32 => 278553008u32,
        1078382593u32 => 318418273u32,
        1078398977u32 => 76817184u32,
        1078497281u32 => 564767901u32,
        1078505473u32 => 903586484u32,
        1078546433u32 => 892495829u32,
        1078571009u32 => 477493167u32,
        1078620161u32 => 1012982146u32,
        1078775809u32 => 1009032663u32,
        1078849537u32 => 458098999u32,
        1078939649u32 => 709331998u32,
        1079185409u32 => 920135442u32,
        1079193601u32 => 448081243u32,
        1079357441u32 => 736735747u32,
        1079431169u32 => 854874008u32,
        1079603201u32 => 773294351u32,
        1079627777u32 => 675655161u32,
        1079635969u32 => 1058459870u32,
        1079676929u32 => 101995690u32,
        1079685121u32 => 32173423u32,
        1079734273u32 => 636580581u32,
        1079832577u32 => 604476873u32,
        1079848961u32 => 451833406u32,
        1079873537u32 => 551428435u32,
        1079922689u32 => 239732173u32,
        1080020993u32 => 577564692u32,
        1080348673u32 => 60273127u32,
        1080365057u32 => 187456814u32,
        1080446977u32 => 165559212u32,
        1080496129u32 => 420867667u32,
        1080512513u32 => 988127443u32,
        1080741889u32 => 614603713u32,
        1080864769u32 => 706481575u32,
        1080938497u32 => 950281766u32,
        1080954881u32 => 904289190u32,
        1081077761u32 => 372675434u32,
        1081225217u32 => 753138445u32,
        1081323521u32 => 54415730u32,
        1081348097u32 => 1062098443u32,
        1081397249u32 => 108806319u32,
        1081430017u32 => 626915142u32,
        1081495553u32 => 805451990u32,
        1081577473u32 => 984110925u32,
        1081643009u32 => 332412100u32,
        1081700353u32 => 774563238u32,
        1081774081u32 => 539468083u32,
        1081815041u32 => 328118502u32,
        1081839617u32 => 254936075u32,
        1081921537u32 => 1077316205u32,
        1082060801u32 => 69423874u32,
        1082093569u32 => 377574571u32,
        1082331137u32 => 221907010u32,
        1082462209u32 => 4550373u32,
        1082552321u32 => 553249034u32,
        1082601473u32 => 160865354u32,
        1082699777u32 => 356063880u32,
        1082724353u32 => 959445664u32,
        1082806273u32 => 939905480u32,
        1082904577u32 => 221888890u32,
        1082929153u32 => 1076011902u32,
        1083002881u32 => 492291370u32,
        1083076609u32 => 852224816u32,
        1083117569u32 => 562586789u32,
        1083215873u32 => 849898814u32,
        1083289601u32 => 371855449u32,
        1083371521u32 => 26051163u32,
        1083518977u32 => 427397952u32,
        1083535361u32 => 428357584u32,
        1083928577u32 => 433352927u32,
        1083953153u32 => 566277417u32,
        1084076033u32 => 330275291u32,
        1084133377u32 => 118956598u32,
        1084149761u32 => 406288370u32,
        1084297217u32 => 59884544u32,
        1084477441u32 => 860189997u32,
        1084518401u32 => 325922673u32,
        1084641281u32 => 552787249u32,
        1084649473u32 => 305293117u32,
        1084674049u32 => 494639359u32,
        1084690433u32 => 822308407u32,
        1084813313u32 => 971212057u32,
        1084911617u32 => 802486144u32,
        1084968961u32 => 686066254u32,
        1085140993u32 => 389974033u32,
        1085181953u32 => 194931248u32,
        1085255681u32 => 910206612u32,
        1085411329u32 => 113971541u32,
        1085501441u32 => 1929545u32,
        1085534209u32 => 638630316u32,
        1085550593u32 => 953103349u32,
        1085607937u32 => 979447493u32,
        1085648897u32 => 77492782u32,
        1085673473u32 => 503885761u32,
        1085779969u32 => 890340600u32,
        1085796353u32 => 490324652u32,
        1085870081u32 => 585095591u32,
        1085894657u32 => 490296731u32,
        1085952001u32 => 5925332u32,
        1086066689u32 => 885416335u32,
        1086115841u32 => 190204103u32,
        1086345217u32 => 715431564u32,
        1086410753u32 => 256172258u32,
        1086443521u32 => 336400907u32,
        1086566401u32 => 1067959688u32,
        1086713857u32 => 147662199u32,
        1086763009u32 => 48119824u32,
        1086853121u32 => 120659900u32,
        1086902273u32 => 1003428496u32,
        1086935041u32 => 242459986u32,
        1087025153u32 => 75501931u32,
        1087148033u32 => 953053731u32,
        1087254529u32 => 1516956u32,
        1087303681u32 => 569239654u32,
        1087418369u32 => 786605394u32,
        1087451137u32 => 729886271u32,
        1087475713u32 => 151137157u32,
        1087492097u32 => 776112340u32,
        1087762433u32 => 388310431u32,
        1087836161u32 => 494416622u32,
        1087860737u32 => 1004610092u32,
        1087885313u32 => 456447914u32,
        1088311297u32 => 287032140u32,
        1088376833u32 => 172248198u32,
        1088401409u32 => 367557804u32,
        1088409601u32 => 543705124u32,
        1088458753u32 => 1076667398u32,
        1088606209u32 => 872062693u32,
        1088647169u32 => 93429298u32,
        1088679937u32 => 554410168u32,
        1088778241u32 => 344087731u32,
        1088802817u32 => 838603145u32,
        1088851969u32 => 42161934u32,
        1088868353u32 => 777817639u32,
        1088966657u32 => 928901386u32,
        1089220609u32 => 272825340u32,
        1089458177u32 => 743527986u32,
        1089540097u32 => 941988729u32,
        1089835009u32 => 526076662u32,
        1089949697u32 => 49250914u32,
        1090170881u32 => 799148519u32,
        1090179073u32 => 465706394u32,
        1090203649u32 => 687726218u32,
        1090293761u32 => 1077146369u32,
        1090400257u32 => 772489995u32,
        1090490369u32 => 507587641u32,
        1090498561u32 => 657951539u32,
        1090646017u32 => 958651655u32,
        1090768897u32 => 413680831u32,
        1090867201u32 => 972888051u32,
        1091014657u32 => 734854902u32,
        1091104769u32 => 872004536u32,
        1091178497u32 => 788488870u32,
        1091399681u32 => 889570566u32,
        1091571713u32 => 266756696u32,
        1091850241u32 => 149846777u32,
        1092022273u32 => 81479489u32,
        1092038657u32 => 439267915u32,
        1092210689u32 => 207395746u32,
        1092268033u32 => 41425309u32,
        1092333569u32 => 373218382u32,
        1092636673u32 => 126926490u32,
        1092653057u32 => 670001121u32,
        1092661249u32 => 839575328u32,
        1092734977u32 => 121024808u32,
        1092882433u32 => 602295797u32,
        1092923393u32 => 1089026827u32,
        1092997121u32 => 550649822u32,
        1093005313u32 => 716338128u32,
        1093152769u32 => 462515386u32,
        1093292033u32 => 702477405u32,
        1093414913u32 => 1081565042u32,
        1093439489u32 => 575620762u32,
        1093513217u32 => 177675631u32,
        1093537793u32 => 357048573u32,
        1093570561u32 => 345140391u32,
        1093939201u32 => 194522895u32,
        1094184961u32 => 881520392u32,
        1094250497u32 => 771171580u32,
        1094258689u32 => 216718420u32,
        1094299649u32 => 570318207u32,
        1094356993u32 => 491343794u32,
        1094381569u32 => 104624202u32,
        1094430721u32 => 76172906u32,
        1094471681u32 => 709255933u32,
        1094504449u32 => 123760995u32,
        1094725633u32 => 906454852u32,
        1094889473u32 => 678882721u32,
        1094914049u32 => 1033994052u32,
        1094946817u32 => 935106589u32,
        1094963201u32 => 302058382u32,
        1094995969u32 => 247660174u32,
        1095012353u32 => 236026821u32,
        1095045121u32 => 340929558u32,
        1095331841u32 => 562566707u32,
        1095462913u32 => 534147584u32,
        1095536641u32 => 915550825u32,
        1095585793u32 => 369873726u32,
        1095700481u32 => 112705433u32,
        1095708673u32 => 245278477u32,
        1095774209u32 => 427460947u32,
        1095847937u32 => 449858445u32,
        1096101889u32 => 747441117u32,
        1096142849u32 => 630098391u32,
        1096175617u32 => 983326237u32,
        1096216577u32 => 423792479u32,
        1096224769u32 => 293500451u32,
        1096339457u32 => 86552891u32,
        1096388609u32 => 37934936u32,
        1096396801u32 => 619985937u32,
        1096486913u32 => 1028728378u32,
        1096544257u32 => 504978889u32,
        1096609793u32 => 962136079u32,
        1096716289u32 => 892281574u32,
        1096765441u32 => 926997888u32,
        1096880129u32 => 762394306u32,
        1096937473u32 => 998742517u32,
        1096953857u32 => 865994416u32,
        1096962049u32 => 770576444u32,
        1097076737u32 => 158777451u32,
        1097175041u32 => 413260505u32,
        1097183233u32 => 438975645u32,
        1097199617u32 => 809749519u32,
        1097256961u32 => 978766530u32,
        1097306113u32 => 693965359u32,
        1097322497u32 => 413131659u32,
        1097453569u32 => 43346590u32,
        1097469953u32 => 94225258u32,
        1097527297u32 => 561046154u32,
        1097691137u32 => 555315574u32,
        1097715713u32 => 656718977u32,
        1097748481u32 => 521726197u32,
        1097895937u32 => 871117632u32,
        1098035201u32 => 564115921u32,
        1098084353u32 => 681692915u32,
        1098231809u32 => 507436u32,
        1098280961u32 => 98337196u32,
        1098313729u32 => 154855986u32,
        1098387457u32 => 521613139u32,
        1098403841u32 => 573872694u32,
        1098485761u32 => 3596132u32,
        1098731521u32 => 922556928u32,
        1098756097u32 => 842535944u32,
        1098780673u32 => 936946126u32,
        1098805249u32 => 813305111u32,
        1098854401u32 => 747181319u32,
        1098919937u32 => 779654180u32,
        1099042817u32 => 221621754u32,
        1099067393u32 => 229374903u32,
        1099272193u32 => 50695454u32,
        1099296769u32 => 629428489u32,
        1099370497u32 => 22597950u32,
        1099411457u32 => 906998873u32,
        1099591681u32 => 133400545u32,
        1099665409u32 => 683014298u32,
        1099755521u32 => 90542198u32,
        1099763713u32 => 348542213u32,
        1100001281u32 => 876342001u32,
        1100132353u32 => 651143689u32,
        1100173313u32 => 1071346568u32,
        1100296193u32 => 540365668u32,
        1100451841u32 => 1001129997u32,
        1100492801u32 => 1060327745u32,
        1100500993u32 => 531594637u32,
        1100574721u32 => 689765852u32,
        1100623873u32 => 802083108u32,
        1100689409u32 => 864475742u32,
        1100697601u32 => 157977060u32,
        1100812289u32 => 377262549u32,
        1100820481u32 => 163179244u32,
        1100861441u32 => 481121493u32,
        1100894209u32 => 61702150u32,
        1101107201u32 => 991255617u32,
        1101139969u32 => 844342533u32,
        1101189121u32 => 1093652715u32,
        1101213697u32 => 434481171u32,
        _ => unimplemented!(),
    }
}

impl<const P: u32> PartialEq for PrimeField<P> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<const P: u32> Neg for PrimeField<P> {
    type Output = PrimeField<P>;

    fn neg(self) -> Self::Output {
        PrimeField::<P>((P - self.0) % P)
    }
}

impl<const P: u32> From<u32> for PrimeField<P> {
    fn from(value: u32) -> Self {
        PrimeField::<P>(value % P)
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
                    .map(|mmi| BigUint::from(mmi.clone()).bits())
                    .max()
                    .unwrap() as usize
        );
        assert!(
            rhs_max_bits
                >= other
                    .coefficients
                    .iter()
                    .map(|mmi| BigUint::from(mmi.clone()).bits())
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

impl From<Polynomial<BigUint>> for Polynomial<MultiModInt> {
    fn from(val: Polynomial<BigUint>) -> Self {
        Polynomial::new(
            val.coefficients
                .into_iter()
                .map(|bi| bi.into())
                .collect_vec(),
        )
    }
}

impl From<Polynomial<BigInt>> for Polynomial<MultiModInt> {
    fn from(val: Polynomial<BigInt>) -> Self {
        Polynomial::new(
            val.coefficients
                .into_iter()
                .map(|bu| bu.into())
                .collect_vec(),
        )
    }
}

#[cfg(test)]
mod test {

    use std::str::FromStr;

    use itertools::Itertools;
    use num::{bigint::ToBigInt, BigInt, BigUint, Integer, One};
    use proptest::arbitrary::Arbitrary;
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest::proptest;
    use proptest::strategy::BoxedStrategy;
    use proptest::strategy::Just;
    use proptest::strategy::Strategy;
    use rand::{rngs::StdRng, thread_rng, Rng, SeedableRng};
    use test_strategy::proptest as strategy_proptest;

    use crate::{multimod::N, polynomial::Polynomial};

    use super::{coefficient, root_of_unity_4096, MultiModInt, MODULI};

    fn random_bit_vector(length: usize, seed: [u8; 32]) -> Vec<u8> {
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        (0..length).map(|_| rng.gen::<u8>() % 2).collect_vec()
    }

    fn random_biguint(bitlen: usize, seed: [u8; 32]) -> BigUint {
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        BigUint::from_radix_be(&random_bit_vector(bitlen, rng.gen()), 2).unwrap()
    }

    fn random_bigint(bitlen: usize, seed: [u8; 32]) -> BigInt {
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let biguint = random_biguint(bitlen, rng.gen()).to_bigint().unwrap();
        if rng.gen() {
            -biguint
        } else {
            biguint
        }
    }

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
        arbitrary_bigint(bitlen).prop_map(MultiModInt::from).boxed()
    }

    fn arbitrary_multimodint_polynomial(
        bitlen: usize,
        num_coefficients: usize,
    ) -> BoxedStrategy<Polynomial<MultiModInt>> {
        vec(arbitrary_multimodint(bitlen), num_coefficients)
            .prop_map(Polynomial::new)
            .boxed()
    }

    #[test]
    fn coefficients_are_1_mod_p() {
        for (i, p, c) in MODULI
            .into_iter()
            .enumerate()
            .map(|(i, p)| (i, p, coefficient(i, N)))
        {
            assert_eq!(
                c.mod_floor(&Into::<BigUint>::into(p)),
                BigUint::one(),
                "mismatch in position {i}\nmodulus is {p}\ncoefficient is {}",
                c
            );
        }
    }

    #[test]
    fn roots_have_right_order() {
        for p in MODULI {
            let mut r = root_of_unity_4096(p);
            for _ in 0..12 {
                assert_ne!(r, 1);
                r = (((r as u64) * (r as u64)) % (p as u64)) as u32;
            }
            assert_eq!(r, 1);
        }
    }

    #[strategy_proptest]
    fn casting_preserves_bitsize(
        #[strategy(1usize..1000)] bits: usize,
        #[strategy(arbitrary_bigint(#bits))] bigint: BigInt,
    ) {
        let mmi = MultiModInt::from(bigint);

        assert!(
            BigInt::from(mmi.clone()).bits() <= bits as u64,
            "{} > {}",
            BigInt::from(mmi).bits(),
            bits
        );
    }

    #[test]
    fn multimodint_arithmetic() {
        let mut rng = thread_rng();
        let bitlen = 1000;

        let a = random_bigint(bitlen, rng.gen());
        let b = random_bigint(bitlen, rng.gen());
        let c = random_bigint(bitlen, rng.gen());

        let d = a.clone() * b.clone() - c.clone();
        let d_mmi = MultiModInt::from(d);
        let mmi_d = MultiModInt::from(a) * MultiModInt::from(b) - MultiModInt::from(c);
        assert_eq!(d_mmi, mmi_d);
    }

    proptest! {
        #[test]
        fn to_and_fro_proptest(a in "[+-]?[0-9]{1,2980}") {
            let b = BigInt::from_str(&a).unwrap();
            let c = MultiModInt::from(b.clone());
            let d = BigInt::from(c.clone());
            let e = MultiModInt::from(d.clone());
            prop_assert_eq!(b, d);
            prop_assert_eq!(c, e);
        }
    }

    #[test]
    fn to_and_fro_fixed() {
        let b = BigInt::from_str("-2098365342651410771540682191176265762931285").unwrap();
        let c = MultiModInt::from(b.clone());
        let d = BigInt::from(c);
        assert_eq!(b, d);
    }

    #[test]
    fn multimodint_fixed() {
        let a = BigInt::from_str("-322777953413029095759719085619503741066179").unwrap();
        println!("a has {} bits", a.bits());
        let b = MultiModInt::from(a.clone());
        println!("b has {} limbs", b.limbs.len());
        let c = BigInt::from(b);
        assert_eq!(c, a);
    }

    #[test]
    fn multimodint_polynomial_arithmetic_fixed() {
        let seed = [
            0x94, 0x88, 0xc2, 0xe5, 0x3d, 0xde, 0xdc, 0xa4, 0xa6, 0x32, 0x7d, 0x7, 0xbd, 0x22,
            0xf4, 0x7a, 0xa7, 0x94, 0x6d, 0xaf, 0x8f, 0xb7, 0x28, 0x22, 0x53, 0x78, 0x80, 0x64,
            0xea, 0xd6, 0xd2, 0xd6,
        ];
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        // let logn = rng.gen_range(1..10);
        // let logn = rng.gen_range(3..10);
        let logn = 4;
        let n = 1 << logn;
        let bitlen = 70;

        let a = Polynomial::new(
            (0..n)
                .map(|_| random_bigint(bitlen, rng.gen()))
                .collect_vec(),
        );
        let b = Polynomial::new(
            (0..n)
                .map(|_| random_bigint(bitlen, rng.gen()))
                .collect_vec(),
        );

        let mut bigint_terms = vec![];
        let mut multimod_terms = vec![];
        for (ai, ac) in a.coefficients.iter().enumerate() {
            for (bi, bc) in b.coefficients.iter().enumerate() {
                if ai + bi == 11 {
                    bigint_terms.push(ac.clone() * bc);
                    multimod_terms
                        .push(MultiModInt::from(ac.clone()) * MultiModInt::from(bc.clone()));
                }
            }
        }
        assert_eq!(
            bigint_terms
                .iter()
                .map(|bi| MultiModInt::from(bi.clone()))
                .collect_vec(),
            multimod_terms
        );
        let bigint_sum = bigint_terms.into_iter().sum::<BigInt>();
        let multimod_sum = multimod_terms.into_iter().sum();
        assert_eq!(MultiModInt::from(bigint_sum), multimod_sum);
        let d = a.clone() * b.clone();
        let mmi_d = Polynomial::<MultiModInt>::from(d.clone());

        let d_mmi =
            Polynomial::<MultiModInt>::from(a.clone()) * Polynomial::<MultiModInt>::from(b.clone());

        assert_eq!(
            mmi_d,
            d_mmi,
            "\n\nmmi_d: {}\nd_mmi: {} of bitcaps {}\n\nd was: {:?} of bitlens {}\n\na: {:?}\n\nb: {:?}",
            mmi_d,
            d_mmi,
            d_mmi.coefficients.iter().map(|mmi| mmi.bitsize_bound).join(","),
            d,
            d.coefficients.iter().map(|c| c.bits()).join(","),
            a,
            b
        );
    }

    #[test]
    fn multimodint_polynomial_arithmetic() {
        let mut rng = thread_rng();
        let logn = rng.gen_range(3..10);
        let n = 1 << logn;
        let bitlen = 70;

        let a = Polynomial::new(
            (0..n)
                .map(|_| random_bigint(bitlen, rng.gen()))
                .collect_vec(),
        );
        let b = Polynomial::new(
            (0..n)
                .map(|_| random_bigint(bitlen, rng.gen()))
                .collect_vec(),
        );
        let c = Polynomial::new(
            (0..n)
                .map(|_| random_bigint(bitlen, rng.gen()))
                .collect_vec(),
        );

        let d = (a.clone() * b.clone()) - c.clone();
        let mmi_d = Polynomial::<MultiModInt>::from(d.clone());

        let d_mmi = Polynomial::<MultiModInt>::from(a.clone())
            * Polynomial::<MultiModInt>::from(b.clone())
            - Into::<Polynomial<MultiModInt>>::into(c);

        assert_eq!(
            mmi_d,
            d_mmi,
            "\n\nmmi_d: {}\nd_mmi: {} of bitcaps {}\n\nd was: {:?} of bitlens {}\n\na: {:?}\n\nb: {:?}",
            mmi_d,
            d_mmi,
            d_mmi.coefficients.iter().map(|mmi| mmi.bitsize_bound).join(","),
            d,
            d.coefficients.iter().map(|c| c.bits()).join(","),
            a,
            b
        );
    }

    #[test]
    fn cyclotomic_multiplication() {
        let seed: [u8; 32] = thread_rng().gen();
        println!(
            "seed = [{}];",
            seed.iter().map(|c| format!("0x{:x}", c)).join(",")
        );
        let mut rng: StdRng = SeedableRng::from_seed(seed);
        let logn = rng.gen_range(0..10);
        let n = 1 << logn;
        let bitlen = rng.gen_range(0..567);
        let a = Polynomial::new(
            (0..n)
                .map(|_| random_biguint(bitlen, rng.gen()))
                .map(MultiModInt::from)
                .collect_vec(),
        );
        let b = Polynomial::new(
            (0..n)
                .map(|_| random_biguint(bitlen, rng.gen()))
                .map(MultiModInt::from)
                .collect_vec(),
        );

        println!("about to compute old-fashioned product over multimod ints ...");
        let product = a.clone() * b.clone();
        let c_trad = product.reduce_by_cyclotomic(n);

        println!("about to compute cyclotomic product ...");
        let c_fast = a.cyclotomic_mul(&b);
        assert_eq!(c_trad, c_fast, "\nleft: {}\n\nright: {}\n", c_trad, c_fast);
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
}
