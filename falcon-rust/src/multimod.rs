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

const N: usize = 512;
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
    1101352961u32,
    1101361153u32,
    1101385729u32,
    1101434881u32,
    1101549569u32,
    1101557761u32,
    1101598721u32,
    1101623297u32,
    1101672449u32,
    1101729793u32,
    1101852673u32,
    1101926401u32,
    1102041089u32,
    1102163969u32,
    1102295041u32,
    1102360577u32,
    1102483457u32,
    1102532609u32,
    1102565377u32,
    1102712833u32,
    1102860289u32,
    1102934017u32,
    1103147009u32,
    1103196161u32,
    1103245313u32,
    1103343617u32,
    1103548417u32,
    1103835137u32,
    1103843329u32,
    1103917057u32,
    1103966209u32,
    1104007169u32,
    1104261121u32,
    1104408577u32,
    1104457729u32,
    1104654337u32,
    1104695297u32,
    1104703489u32,
    1104719873u32,
    1104867329u32,
    1105195009u32,
    1105309697u32,
    1105367041u32,
    1105408001u32,
    1105457153u32,
    1105489921u32,
    1105514497u32,
    1105555457u32,
    1105580033u32,
    1105604609u32,
    1105686529u32,
    1105784833u32,
    1105825793u32,
    1105883137u32,
    1105907713u32,
    1105981441u32,
    1106030593u32,
    1106071553u32,
    1106227201u32,
    1106300929u32,
    1106350081u32,
    1106374657u32,
    1106399233u32,
    1106464769u32,
    1106513921u32,
    1106538497u32,
    1106546689u32,
    1106710529u32,
    1106931713u32,
    1107030017u32,
    1107054593u32,
    1107521537u32,
    1107619841u32,
    1107668993u32,
    1107816449u32,
    1108013057u32,
    1108021249u32,
    1108062209u32,
    1108135937u32,
    1108144129u32,
    1108217857u32,
    1108381697u32,
    1108430849u32,
    1108439041u32,
    1108463617u32,
    1108488193u32,
    1108553729u32,
    1108586497u32,
    1108611073u32,
    1108750337u32,
    1108774913u32,
    1108922369u32,
    1108979713u32,
    1109004289u32,
    1109094401u32,
    1109127169u32,
    1109217281u32,
    1109299201u32,
    1109323777u32,
    1109618689u32,
    1109692417u32,
    1109708801u32,
    1109880833u32,
    1109962753u32,
    1110085633u32,
    1110151169u32,
    1110224897u32,
    1110249473u32,
    1110282241u32,
    1110306817u32,
    1110454273u32,
    1110716417u32,
    1111109633u32,
    1111166977u32,
    1111183361u32,
    1111216129u32,
    1111232513u32,
    1111257089u32,
    1111412737u32,
    1111461889u32,
    1111502849u32,
    1111560193u32,
    1111584769u32,
    1111625729u32,
    1111658497u32,
    1111683073u32,
    1111756801u32,
    1111797761u32,
    1111822337u32,
    1111928833u32,
    1111994369u32,
    1112043521u32,
    1112174593u32,
    1112289281u32,
    1112338433u32,
    1112363009u32,
    1112371201u32,
    1112494081u32,
    1112567809u32,
    1112682497u32,
    1112739841u32,
    1112887297u32,
    1112952833u32,
    1112977409u32,
    1113059329u32,
    1113346049u32,
    1113403393u32,
    1113567233u32,
    1113600001u32,
    1113747457u32,
    1113886721u32,
    1113894913u32,
    1113935873u32,
    1114034177u32,
    1114165249u32,
    1114238977u32,
    1114451969u32,
    1114632193u32,
    1114746881u32,
    1114771457u32,
    1114894337u32,
    1115025409u32,
    1115074561u32,
    1115148289u32,
    1115410433u32,
    1115492353u32,
    1115590657u32,
    1115615233u32,
    1115631617u32,
    1115713537u32,
    1115729921u32,
    1115762689u32,
    1115803649u32,
    1115901953u32,
    1116131329u32,
    1116270593u32,
    1116418049u32,
    1116549121u32,
    1116639233u32,
    1116794881u32,
    1116917761u32,
    1116991489u32,
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
        1073754113u32 => 48440u32,
        1073950721u32 => 32033u32,
        1073958913u32 => 105315u32,
        1073983489u32 => 661217u32,
        1074196481u32 => 2872u32,
        1074343937u32 => 312461u32,
        1074442241u32 => 135843u32,
        1074475009u32 => 433069u32,
        1074515969u32 => 407521u32,
        1074524161u32 => 988646u32,
        1074548737u32 => 211582u32,
        1074589697u32 => 292708u32,
        1074597889u32 => 660222u32,
        1074688001u32 => 416412u32,
        1074696193u32 => 361109u32,
        1074810881u32 => 1501888u32,
        1074835457u32 => 1002812u32,
        1074941953u32 => 70152u32,
        1075007489u32 => 72221u32,
        1075064833u32 => 91935u32,
        1075105793u32 => 96986u32,
        1075351553u32 => 491583u32,
        1075376129u32 => 23197u32,
        1075449857u32 => 442926u32,
        1075507201u32 => 65212u32,
        1075621889u32 => 28560u32,
        1075695617u32 => 1384994u32,
        1075720193u32 => 2355626u32,
        1075752961u32 => 171856u32,
        1076039681u32 => 582542u32,
        1076064257u32 => 60608u32,
        1076162561u32 => 79179u32,
        1076187137u32 => 11582u32,
        1076211713u32 => 475437u32,
        1076269057u32 => 1506429u32,
        1076334593u32 => 170255u32,
        1076391937u32 => 116594u32,
        1076531201u32 => 457251u32,
        1076613121u32 => 23490u32,
        1076908033u32 => 888888u32,
        1076948993u32 => 175992u32,
        1077047297u32 => 100385u32,
        1077178369u32 => 361242u32,
        1077252097u32 => 878194u32,
        1077391361u32 => 65577u32,
        1077399553u32 => 828959u32,
        1077620737u32 => 139933u32,
        1077686273u32 => 139655u32,
        1077792769u32 => 622031u32,
        1077882881u32 => 794707u32,
        1078038529u32 => 613293u32,
        1078079489u32 => 155863u32,
        1078112257u32 => 75172u32,
        1078153217u32 => 461903u32,
        1078210561u32 => 932090u32,
        1078333441u32 => 1439401u32,
        1078382593u32 => 197705u32,
        1078398977u32 => 68737u32,
        1078497281u32 => 651960u32,
        1078505473u32 => 703312u32,
        1078546433u32 => 284662u32,
        1078571009u32 => 250295u32,
        1078620161u32 => 2982783u32,
        1078775809u32 => 899145u32,
        1078849537u32 => 84761u32,
        1078939649u32 => 22517u32,
        1079185409u32 => 360819u32,
        1079193601u32 => 119679u32,
        1079357441u32 => 55440u32,
        1079431169u32 => 624161u32,
        1079603201u32 => 563187u32,
        1079627777u32 => 1256363u32,
        1079635969u32 => 7546u32,
        1079676929u32 => 587544u32,
        1079685121u32 => 43544u32,
        1079734273u32 => 774321u32,
        1079832577u32 => 720864u32,
        1079848961u32 => 425219u32,
        1079873537u32 => 1471550u32,
        1079922689u32 => 65948u32,
        1080020993u32 => 667559u32,
        1080348673u32 => 507532u32,
        1080365057u32 => 251673u32,
        1080446977u32 => 608032u32,
        1080496129u32 => 218661u32,
        1080512513u32 => 543601u32,
        1080741889u32 => 99383u32,
        1080864769u32 => 297866u32,
        1080938497u32 => 423622u32,
        1080954881u32 => 1590842u32,
        1081077761u32 => 130880u32,
        1081225217u32 => 433113u32,
        1081323521u32 => 64610u32,
        1081348097u32 => 411475u32,
        1081397249u32 => 882483u32,
        1081430017u32 => 110379u32,
        1081495553u32 => 1261186u32,
        1081577473u32 => 1587852u32,
        1081643009u32 => 1187753u32,
        1081700353u32 => 1423968u32,
        1081774081u32 => 94011u32,
        1081815041u32 => 184865u32,
        1081839617u32 => 161130u32,
        1081921537u32 => 202318u32,
        1082060801u32 => 33230u32,
        1082093569u32 => 365317u32,
        1082331137u32 => 432020u32,
        1082462209u32 => 396831u32,
        1082552321u32 => 291384u32,
        1082601473u32 => 1434334u32,
        1082699777u32 => 372839u32,
        1082724353u32 => 550796u32,
        1082806273u32 => 1136235u32,
        1082904577u32 => 1094060u32,
        1082929153u32 => 169112u32,
        1083002881u32 => 523372u32,
        1083076609u32 => 567555u32,
        1083117569u32 => 157143u32,
        1083215873u32 => 550395u32,
        1083289601u32 => 715848u32,
        1083371521u32 => 122597u32,
        1083518977u32 => 1917704u32,
        1083535361u32 => 210077u32,
        1083928577u32 => 214585u32,
        1083953153u32 => 272789u32,
        1084076033u32 => 928745u32,
        1084133377u32 => 1145577u32,
        1084149761u32 => 157065u32,
        1084297217u32 => 261225u32,
        1084477441u32 => 798309u32,
        1084518401u32 => 1109373u32,
        1084641281u32 => 282397u32,
        1084649473u32 => 2195937u32,
        1084674049u32 => 57476u32,
        1084690433u32 => 1162223u32,
        1084813313u32 => 1463065u32,
        1084911617u32 => 955149u32,
        1084968961u32 => 423861u32,
        1085140993u32 => 209024u32,
        1085181953u32 => 1251435u32,
        1085255681u32 => 3229560u32,
        1085411329u32 => 29553u32,
        1085501441u32 => 1192651u32,
        1085534209u32 => 601864u32,
        1085550593u32 => 15692u32,
        1085607937u32 => 78870u32,
        1085648897u32 => 570856u32,
        1085673473u32 => 1745024u32,
        1085779969u32 => 95199u32,
        1085796353u32 => 243260u32,
        1085870081u32 => 1380116u32,
        1085894657u32 => 717175u32,
        1085952001u32 => 501293u32,
        1086066689u32 => 327442u32,
        1086115841u32 => 148163u32,
        1086345217u32 => 20315u32,
        1086410753u32 => 308768u32,
        1086443521u32 => 951292u32,
        1086566401u32 => 668235u32,
        1086713857u32 => 308667u32,
        1086763009u32 => 1770764u32,
        1086853121u32 => 910129u32,
        1086902273u32 => 1269422u32,
        1086935041u32 => 200643u32,
        1087025153u32 => 472048u32,
        1087148033u32 => 1169282u32,
        1087254529u32 => 1472309u32,
        1087303681u32 => 1737639u32,
        1087418369u32 => 453419u32,
        1087451137u32 => 181615u32,
        1087475713u32 => 207058u32,
        1087492097u32 => 31765u32,
        1087762433u32 => 67206u32,
        1087836161u32 => 699743u32,
        1087860737u32 => 668582u32,
        1087885313u32 => 368647u32,
        1088311297u32 => 327190u32,
        1088376833u32 => 194045u32,
        1088401409u32 => 121812u32,
        1088409601u32 => 6713u32,
        1088458753u32 => 1905892u32,
        1088606209u32 => 203u32,
        1088647169u32 => 529868u32,
        1088679937u32 => 1882420u32,
        1088778241u32 => 480252u32,
        1088802817u32 => 406893u32,
        1088851969u32 => 646348u32,
        1088868353u32 => 537738u32,
        1088966657u32 => 76592u32,
        1089220609u32 => 426417u32,
        1089458177u32 => 258908u32,
        1089540097u32 => 183616u32,
        1089835009u32 => 1037166u32,
        1089949697u32 => 785216u32,
        1090170881u32 => 737986u32,
        1090179073u32 => 1601904u32,
        1090203649u32 => 940119u32,
        1090293761u32 => 183634u32,
        1090400257u32 => 652515u32,
        1090490369u32 => 1171802u32,
        1090498561u32 => 621812u32,
        1090646017u32 => 252683u32,
        1090768897u32 => 461198u32,
        1090867201u32 => 1111450u32,
        1091014657u32 => 70051u32,
        1091104769u32 => 578202u32,
        1091178497u32 => 222666u32,
        1091399681u32 => 532774u32,
        1091571713u32 => 615078u32,
        1091850241u32 => 336481u32,
        1092022273u32 => 341326u32,
        1092038657u32 => 184273u32,
        1092210689u32 => 94284u32,
        1092268033u32 => 288226u32,
        1092333569u32 => 32650u32,
        1092636673u32 => 193918u32,
        1092653057u32 => 285361u32,
        1092661249u32 => 157630u32,
        1092734977u32 => 1255483u32,
        1092882433u32 => 331961u32,
        1092923393u32 => 161485u32,
        1092997121u32 => 366401u32,
        1093005313u32 => 195943u32,
        1093152769u32 => 416123u32,
        1093292033u32 => 27677u32,
        1093414913u32 => 4288672u32,
        1093439489u32 => 531933u32,
        1093513217u32 => 2083229u32,
        1093537793u32 => 641327u32,
        1093570561u32 => 847817u32,
        1093939201u32 => 273033u32,
        1094184961u32 => 479761u32,
        1094250497u32 => 200174u32,
        1094258689u32 => 423493u32,
        1094299649u32 => 178872u32,
        1094356993u32 => 357711u32,
        1094381569u32 => 1630200u32,
        1094430721u32 => 141878u32,
        1094471681u32 => 882065u32,
        1094504449u32 => 344233u32,
        1094725633u32 => 401775u32,
        1094889473u32 => 764843u32,
        1094914049u32 => 58256u32,
        1094946817u32 => 587187u32,
        1094963201u32 => 136207u32,
        1094995969u32 => 248448u32,
        1095012353u32 => 102402u32,
        1095045121u32 => 200987u32,
        1095331841u32 => 68193u32,
        1095462913u32 => 2040467u32,
        1095536641u32 => 290154u32,
        1095585793u32 => 390837u32,
        1095700481u32 => 790689u32,
        1095708673u32 => 222115u32,
        1095774209u32 => 1155165u32,
        1095847937u32 => 41642u32,
        1096101889u32 => 142172u32,
        1096142849u32 => 449999u32,
        1096175617u32 => 505600u32,
        1096216577u32 => 471942u32,
        1096224769u32 => 582635u32,
        1096339457u32 => 1853282u32,
        1096388609u32 => 722966u32,
        1096396801u32 => 432655u32,
        1096486913u32 => 556090u32,
        1096544257u32 => 169596u32,
        1096609793u32 => 1047385u32,
        1096716289u32 => 74554u32,
        1096765441u32 => 195662u32,
        1096880129u32 => 431375u32,
        1096937473u32 => 492649u32,
        1096953857u32 => 915628u32,
        1096962049u32 => 1351301u32,
        1097076737u32 => 156621u32,
        1097175041u32 => 254886u32,
        1097183233u32 => 1155302u32,
        1097199617u32 => 281432u32,
        1097256961u32 => 325587u32,
        1097306113u32 => 362611u32,
        1097322497u32 => 301402u32,
        1097453569u32 => 255417u32,
        1097469953u32 => 93403u32,
        1097527297u32 => 153882u32,
        1097691137u32 => 73050u32,
        1097715713u32 => 595907u32,
        1097748481u32 => 132928u32,
        1097895937u32 => 537145u32,
        1098035201u32 => 46695u32,
        1098084353u32 => 610267u32,
        1098231809u32 => 507436u32,
        1098280961u32 => 1746911u32,
        1098313729u32 => 337552u32,
        1098387457u32 => 461008u32,
        1098403841u32 => 1086191u32,
        1098485761u32 => 225278u32,
        1098731521u32 => 530762u32,
        1098756097u32 => 299596u32,
        1098780673u32 => 20725u32,
        1098805249u32 => 448996u32,
        1098854401u32 => 646034u32,
        1098919937u32 => 391747u32,
        1099042817u32 => 111404u32,
        1099067393u32 => 98150u32,
        1099272193u32 => 416107u32,
        1099296769u32 => 92904u32,
        1099370497u32 => 3472254u32,
        1099411457u32 => 818191u32,
        1099591681u32 => 345675u32,
        1099665409u32 => 33u32,
        1099755521u32 => 268338u32,
        1099763713u32 => 1805u32,
        1100001281u32 => 94581u32,
        1100132353u32 => 1290372u32,
        1100173313u32 => 333023u32,
        1100296193u32 => 153236u32,
        1100451841u32 => 76757u32,
        1100492801u32 => 445592u32,
        1100500993u32 => 225357u32,
        1100574721u32 => 37002u32,
        1100623873u32 => 666268u32,
        1100689409u32 => 1135312u32,
        1100697601u32 => 112749u32,
        1100812289u32 => 415011u32,
        1100820481u32 => 543808u32,
        1100861441u32 => 1128678u32,
        1100894209u32 => 825338u32,
        1101107201u32 => 733694u32,
        1101139969u32 => 424569u32,
        1101189121u32 => 47273u32,
        1101213697u32 => 1016392u32,
        1101352961u32 => 309493u32,
        1101361153u32 => 947246u32,
        1101385729u32 => 38614u32,
        1101434881u32 => 136004u32,
        1101549569u32 => 2653628u32,
        1101557761u32 => 1177114u32,
        1101598721u32 => 314884u32,
        1101623297u32 => 45923u32,
        1101672449u32 => 643183u32,
        1101729793u32 => 809336u32,
        1101852673u32 => 69172u32,
        1101926401u32 => 668032u32,
        1102041089u32 => 189480u32,
        1102163969u32 => 310669u32,
        1102295041u32 => 528599u32,
        1102360577u32 => 290563u32,
        1102483457u32 => 652348u32,
        1102532609u32 => 40118u32,
        1102565377u32 => 270152u32,
        1102712833u32 => 107086u32,
        1102860289u32 => 255194u32,
        1102934017u32 => 836951u32,
        1103147009u32 => 129893u32,
        1103196161u32 => 614338u32,
        1103245313u32 => 101343u32,
        1103343617u32 => 223151u32,
        1103548417u32 => 158015u32,
        1103835137u32 => 815159u32,
        1103843329u32 => 367117u32,
        1103917057u32 => 376614u32,
        1103966209u32 => 1136651u32,
        1104007169u32 => 963351u32,
        1104261121u32 => 262452u32,
        1104408577u32 => 359731u32,
        1104457729u32 => 162881u32,
        1104654337u32 => 73618u32,
        1104695297u32 => 741467u32,
        1104703489u32 => 633459u32,
        1104719873u32 => 50337u32,
        1104867329u32 => 222671u32,
        1105195009u32 => 73074u32,
        1105309697u32 => 21330u32,
        1105367041u32 => 72795u32,
        1105408001u32 => 37382u32,
        1105457153u32 => 154285u32,
        1105489921u32 => 152634u32,
        1105514497u32 => 2322063u32,
        1105555457u32 => 189954u32,
        1105580033u32 => 1038715u32,
        1105604609u32 => 665371u32,
        1105686529u32 => 813748u32,
        1105784833u32 => 176627u32,
        1105825793u32 => 282061u32,
        1105883137u32 => 279707u32,
        1105907713u32 => 437100u32,
        1105981441u32 => 33665u32,
        1106030593u32 => 573284u32,
        1106071553u32 => 87231u32,
        1106227201u32 => 688621u32,
        1106300929u32 => 437450u32,
        1106350081u32 => 73065u32,
        1106374657u32 => 799346u32,
        1106399233u32 => 82293u32,
        1106464769u32 => 445923u32,
        1106513921u32 => 1530u32,
        1106538497u32 => 12166u32,
        1106546689u32 => 473383u32,
        1106710529u32 => 112922u32,
        1106931713u32 => 156875u32,
        1107030017u32 => 37145u32,
        1107054593u32 => 1963483u32,
        1107521537u32 => 95438u32,
        1107619841u32 => 896311u32,
        1107668993u32 => 332237u32,
        1107816449u32 => 158636u32,
        1108013057u32 => 693973u32,
        1108021249u32 => 1075931u32,
        1108062209u32 => 327178u32,
        1108135937u32 => 844631u32,
        1108144129u32 => 218304u32,
        1108217857u32 => 592477u32,
        1108381697u32 => 66589u32,
        1108430849u32 => 85392u32,
        1108439041u32 => 5971u32,
        1108463617u32 => 197613u32,
        1108488193u32 => 26650u32,
        1108553729u32 => 1096406u32,
        1108586497u32 => 313332u32,
        1108611073u32 => 216332u32,
        1108750337u32 => 539851u32,
        1108774913u32 => 1300220u32,
        1108922369u32 => 76278u32,
        1108979713u32 => 85299u32,
        1109004289u32 => 2540847u32,
        1109094401u32 => 1626228u32,
        1109127169u32 => 310544u32,
        1109217281u32 => 1167592u32,
        1109299201u32 => 164117u32,
        1109323777u32 => 1415165u32,
        1109618689u32 => 36038u32,
        1109692417u32 => 307845u32,
        1109708801u32 => 361394u32,
        1109880833u32 => 2428759u32,
        1109962753u32 => 1038372u32,
        1110085633u32 => 198062u32,
        1110151169u32 => 479710u32,
        1110224897u32 => 309058u32,
        1110249473u32 => 413300u32,
        1110282241u32 => 658764u32,
        1110306817u32 => 706845u32,
        1110454273u32 => 679209u32,
        1110716417u32 => 514121u32,
        1111109633u32 => 1304575u32,
        1111166977u32 => 927733u32,
        1111183361u32 => 618347u32,
        1111216129u32 => 805101u32,
        1111232513u32 => 918401u32,
        1111257089u32 => 508099u32,
        1111412737u32 => 194754u32,
        1111461889u32 => 428964u32,
        1111502849u32 => 1282672u32,
        1111560193u32 => 26161u32,
        1111584769u32 => 248012u32,
        1111625729u32 => 70070u32,
        1111658497u32 => 277583u32,
        1111683073u32 => 261933u32,
        1111756801u32 => 396371u32,
        1111797761u32 => 1039613u32,
        1111822337u32 => 981652u32,
        1111928833u32 => 920863u32,
        1111994369u32 => 1256515u32,
        1112043521u32 => 685883u32,
        1112174593u32 => 191497u32,
        1112289281u32 => 1387848u32,
        1112338433u32 => 173634u32,
        1112363009u32 => 55556u32,
        1112371201u32 => 364633u32,
        1112494081u32 => 736125u32,
        1112567809u32 => 629534u32,
        1112682497u32 => 658427u32,
        1112739841u32 => 3645846u32,
        1112887297u32 => 561352u32,
        1112952833u32 => 30723u32,
        1112977409u32 => 86423u32,
        1113059329u32 => 279245u32,
        1113346049u32 => 12378u32,
        1113403393u32 => 654172u32,
        1113567233u32 => 1353365u32,
        1113600001u32 => 582870u32,
        1113747457u32 => 410900u32,
        1113886721u32 => 404209u32,
        1113894913u32 => 451817u32,
        1113935873u32 => 281179u32,
        1114034177u32 => 80885u32,
        1114165249u32 => 83446u32,
        1114238977u32 => 126077u32,
        1114451969u32 => 1011248u32,
        1114632193u32 => 6536u32,
        1114746881u32 => 40979u32,
        1114771457u32 => 75118u32,
        1114894337u32 => 273185u32,
        1115025409u32 => 891615u32,
        1115074561u32 => 975014u32,
        1115148289u32 => 683562u32,
        1115410433u32 => 401772u32,
        1115492353u32 => 484440u32,
        1115590657u32 => 413755u32,
        1115615233u32 => 663447u32,
        1115631617u32 => 552041u32,
        1115713537u32 => 18764u32,
        1115729921u32 => 183180u32,
        1115762689u32 => 27878u32,
        1115803649u32 => 1094433u32,
        1115901953u32 => 297757u32,
        1116131329u32 => 774640u32,
        1116270593u32 => 537048u32,
        1116418049u32 => 841858u32,
        1116549121u32 => 905085u32,
        1116639233u32 => 98888u32,
        1116794881u32 => 1214904u32,
        1116917761u32 => 71501u32,
        1116991489u32 => 271036u32,

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

    use crate::{multimod::N, polynomial::Polynomial};

    use super::{coefficient, root_of_unity_4096, MultiModInt, MODULI};

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
}
