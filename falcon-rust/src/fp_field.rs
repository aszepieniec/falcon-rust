use std::fmt::Display;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{One, Zero};
use rand::distr::{Distribution, StandardUniform};

use crate::cyclotomic_fourier::CyclotomicFourier;
use crate::inverse::Inverse;

/// Compute -q^{-1} mod R where R = 2^(q.ilog2()+1).
///
/// Uses five Hensel-lifting steps to get q^{-1} mod 2^32, then negates
/// and masks down to mod R.  Panics at compile time if q is even.
const fn negqinv_modr(q: u32) -> u32 {
    assert!(q & 1 == 1, "FpField requires an odd modulus");
    let mut x = 1u32; // q*x ≡ 1 (mod 2)
    x = x.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(x))); // mod 4
    x = x.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(x))); // mod 16
    x = x.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(x))); // mod 256
    x = x.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(x))); // mod 65536
    x = x.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(x))); // mod 2^32
    let log2r = q.ilog2() + 1;
    let mask = if log2r == 32 { u32::MAX } else { (1u32 << log2r) - 1 };
    x.wrapping_neg() & mask
}

/// Compute R² mod q where R = 2^(q.ilog2()+1).
const fn r_sq_modq(q: u32) -> u32 {
    let log2r = q.ilog2() + 1;
    let r_sq = 1u128 << (2 * log2r);
    (r_sq % q as u128) as u32
}

/// A field element of Z/QZ stored in Montgomery form.
///
/// Q must be an odd prime.  The Montgomery base is R = 2^LOG2R where
/// LOG2R = Q.ilog2() + 1 (the bit length of Q), so R is the smallest
/// power of two strictly greater than Q.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FpField<const Q: u32>(u32);

impl<const Q: u32> FpField<Q> {
    /// log₂ R, where R = 2^LOG2R is the Montgomery base.
    pub(crate) const LOG2R: u32 = Q.ilog2() + 1;
    /// -Q⁻¹ mod R, used in Montgomery reduction.
    pub(crate) const NEGQINV_MODR: u32 = negqinv_modr(Q);
    /// R² mod Q, used to convert canonical values into Montgomery form.
    pub(crate) const R_SQ_MODQ: u32 = r_sq_modq(Q);

    /// Construct a field element from any u32; reduces mod Q automatically.
    pub const fn new(value: u32) -> Self {
        FpField(Self::montyred(value % Q, Self::R_SQ_MODQ))
    }

    /// Return the canonical representative in [0, Q).
    pub const fn value(&self) -> u32 {
        Self::montyred(self.0, 1)
    }

    /// Return the balanced representative in (-Q/2, Q/2].
    #[allow(dead_code)]
    pub fn balanced_value(&self) -> i64 {
        let v = self.value() as i64;
        let g = (v > (Q as i64 / 2)) as i64;
        v - Q as i64 * g
    }

    pub const fn multiply(&self, other: Self) -> Self {
        FpField(Self::montyred(self.0, other.0))
    }

    /// Multiply by 2⁻¹ mod Q without a full Montgomery multiplication.
    ///
    /// (x + Q·(x & 1)) >> 1 ≡ x · 2⁻¹ (mod Q) for any x in [0, Q).
    #[allow(dead_code)]
    pub const fn half(self) -> Self {
        let x = self.0 as u64;
        let r = (x + ((x & 1).wrapping_neg() & Q as u64)) >> 1;
        FpField(r as u32)
    }

    /// Montgomery multiplication: given a, b in [0, Q),
    /// compute a·b·R⁻¹ mod Q where R = 2^LOG2R.
    ///
    /// When LOG2R ≤ 31 the intermediate magic_sum fits in u64, so the
    /// compiler emits only 64-bit arithmetic.  When LOG2R = 32 (Q ≥ 2^31)
    /// it falls back to u128.  The branch is on a compile-time constant so
    /// the dead arm is eliminated by the optimizer.
    const fn montyred(a: u32, b: u32) -> u32 {
        let product: u64 = (a as u64) * (b as u64);

        if Self::LOG2R <= 31 {
            let r_mask = (1u32 << Self::LOG2R) - 1;
            let tail = (product as u32) & r_mask;
            let cofactor = Self::NEGQINV_MODR.wrapping_mul(tail) & r_mask;
            // magic_sum < 2·Q·R < 2^(2·LOG2R+1) ≤ 2^63 — fits in u64
            let magic_sum = product + (cofactor as u64) * (Q as u64);
            let g = magic_sum >> Self::LOG2R;
            (if g >= Q as u64 { g - Q as u64 } else { g }) as u32
        } else {
            // LOG2R == 32; magic_sum can reach ~2^65, use u128
            let cofactor = Self::NEGQINV_MODR.wrapping_mul(product as u32);
            let magic_sum: u128 = product as u128 + cofactor as u128 * Q as u128;
            let g = (magic_sum >> 32) as u64;
            (if g >= Q as u64 { g - Q as u64 } else { g }) as u32
        }
    }
}

impl<const Q: u32> From<usize> for FpField<Q> {
    fn from(value: usize) -> Self {
        FpField::new((value % Q as usize) as u32)
    }
}

impl<const Q: u32> From<i32> for FpField<Q> {
    fn from(value: i32) -> Self {
        if value >= 0 {
            FpField::new((value as u32) % Q)
        } else {
            -FpField::new(value.unsigned_abs() % Q)
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const Q: u32> Add for FpField<Q> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let s = self.0 as u64 + rhs.0 as u64;
        let r = if s >= Q as u64 { s - Q as u64 } else { s };
        FpField(r as u32)
    }
}

impl<const Q: u32> AddAssign for FpField<Q> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<const Q: u32> Sub for FpField<Q> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl<const Q: u32> SubAssign for FpField<Q> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<const Q: u32> Neg for FpField<Q> {
    type Output = FpField<Q>;

    fn neg(self) -> Self::Output {
        let is_nonzero = self.0 != 0;
        let r = Q - self.0;
        FpField(r * (is_nonzero as u32))
    }
}

impl<const Q: u32> Mul for FpField<Q> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        FpField(Self::montyred(self.0, rhs.0))
    }
}

impl<const Q: u32> MulAssign for FpField<Q> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<const Q: u32> Zero for FpField<Q> {
    fn zero() -> Self {
        FpField::new(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl<const Q: u32> One for FpField<Q> {
    fn one() -> Self {
        FpField::new(1)
    }
}

impl<const Q: u32> Distribution<FpField<Q>> for StandardUniform {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> FpField<Q> {
        FpField::new((rng.next_u32() >> 1) % Q)
    }
}

impl<const Q: u32> Display for FpField<Q> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value())
    }
}

impl<const Q: u32> Inverse for FpField<Q> {
    fn inverse_or_zero(self) -> Self {
        // Fermat: self^{Q-2} mod Q. Q-2 has at most LOG2R bits.
        let q_minus_two = Q - 2;
        let mut acc = Self::new(1);
        let mut mask = 1u32 << (Self::LOG2R - 1);
        for _ in 0..Self::LOG2R {
            acc = acc.multiply(acc);
            if mask & q_minus_two != 0 {
                acc = acc.multiply(self);
            }
            mask >>= 1;
        }
        acc
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const Q: u32> Div for FpField<Q> {
    type Output = FpField<Q>;

    fn div(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            panic!("Cannot divide by zero");
        } else {
            self.multiply(rhs.inverse_or_zero())
        }
    }
}

impl CyclotomicFourier for FpField<1073754113> {
    fn primitive_root_of_unity(n: usize) -> Self {
        let log2n = n.ilog2();
        assert!(log2n <= 12);
        // 48440 is a primitive 12th root of unity mod 1073754113
        let mut a = Self::new(48440);
        for _ in 0..(12 - log2n) {
            a *= a;
        }
        a
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num::Zero;
    use rand::{rng, Rng, RngCore};

    use crate::{fp_field::FpField, inverse::Inverse};

    // Small prime — LOG2R = 14, uses u64 branch in montyred.
    type SF = FpField<12289>;
    // Larger prime — LOG2R = 31, still uses u64 branch.
    type BF = FpField<1073754113>;
    // Prime > 2^31 — LOG2R = 32, exercises the u128 branch.
    type LF = FpField<4294967291>;

    /// Compute R⁻¹ mod q by Fermat: (R mod q)^{q-2} mod q.
    fn rinv_modq(q: u64) -> u64 {
        let log2r = (q as u32).ilog2() + 1;
        let r_mod_q = if log2r == 32 {
            (1u64 << 32) % q
        } else {
            (1u64 << log2r) % q
        };
        // Binary exponentiation for r_mod_q^{q-2} mod q
        let mut base = r_mod_q;
        let mut exp = q - 2;
        let mut acc = 1u64;
        while exp > 0 {
            if exp & 1 == 1 {
                acc = (acc as u128 * base as u128 % q as u128) as u64;
            }
            base = (base as u128 * base as u128 % q as u128) as u64;
            exp >>= 1;
        }
        acc
    }

    fn test_montyred_for<const Q: u32>() {
        let mut rng = rng();
        let rinv = rinv_modq(Q as u64);
        for _ in 0..1000 {
            let a = rng.random_range(0..Q);
            let b = rng.random_range(0..Q);
            let expected =
                ((a as u128 * b as u128 % Q as u128) * rinv as u128 % Q as u128) as u32;
            let got = FpField::<Q>::montyred(a, b);
            // montyred(a, b) = a*b*R^{-1} mod Q — but a,b here are canonical,
            // not Montgomery form, so use the direct definition.
            let direct = (a as u128 * b as u128 * rinv as u128 % Q as u128) as u32;
            assert_eq!(got, direct, "Q={Q}, a={a}, b={b}");
            assert_eq!(got, expected);
        }
    }

    #[test]
    fn test_montyred() {
        test_montyred_for::<12289>();
        test_montyred_for::<1073754113>();
        test_montyred_for::<4294967291>();
    }

    fn test_value_for<const Q: u32>() {
        let mut rng = rng();
        for _ in 0..1000 {
            let v = rng.random_range(0..Q);
            let elem = FpField::<Q>::new(v);
            assert_eq!(elem.value(), v, "Q={Q}, v={v}");
        }
    }

    #[test]
    fn test_value() {
        test_value_for::<12289>();
        test_value_for::<1073754113>();
        test_value_for::<4294967291>();
    }

    #[test]
    fn test_add() {
        let mut rng = rng();
        for (a_val, b_val) in (0..100).map(|_| (rng.next_u32() % 12289, rng.next_u32() % 12289)) {
            let a = SF::new(a_val);
            let b = SF::new(b_val);
            let expected = SF::new((a_val + b_val) % 12289);
            assert_eq!(a + b, expected, "a={a_val}, b={b_val}");
        }
    }

    #[test]
    fn test_specific_neg() {
        // -2958 ≡ 9331 (mod 12289) since 2958 + 9331 = 12289
        assert_eq!(-SF::new(2958), SF::new(9331));
        assert_eq!(-SF::new(0), SF::new(0));
    }

    fn test_mul_for<const Q: u32>() {
        let mut rng = rng();
        for _ in 0..1000 {
            let a_val = rng.random_range(0..Q);
            let b_val = rng.random_range(0..Q);
            let expected = (a_val as u64 * b_val as u64 % Q as u64) as u32;
            let a = FpField::<Q>::new(a_val);
            let b = FpField::<Q>::new(b_val);
            assert_eq!((a * b).value(), expected, "Q={Q}, a={a_val}, b={b_val}");
        }
    }

    #[test]
    fn test_mul() {
        test_mul_for::<12289>();
        test_mul_for::<1073754113>();
        test_mul_for::<4294967291>();
    }

    fn test_inverse_for<const Q: u32>() {
        let mut rng = rng();
        let a: FpField<Q> = rng.random();
        let b = a.inverse_or_zero();
        assert_eq!(a * b * a, a);
        assert_eq!(a * b * b, b);
    }

    #[test]
    fn test_inverse() {
        test_inverse_for::<12289>();
        test_inverse_for::<1073754113>();
        test_inverse_for::<4294967291>();
    }

    fn test_batch_inverse_for<const Q: u32>() {
        let mut rng = rng();
        let a: [FpField<Q>; 64] = (0..64)
            .map(|_| rng.random())
            .collect_vec()
            .try_into()
            .unwrap();
        let b_batch = FpField::<Q>::batch_inverse_or_zero(&a);
        let b_regular: Vec<_> = a.iter().map(|e| e.inverse_or_zero()).collect();
        assert_eq!(b_batch, b_regular);
    }

    #[test]
    fn test_batch_inverse() {
        test_batch_inverse_for::<12289>();
        test_batch_inverse_for::<1073754113>();
        test_batch_inverse_for::<4294967291>();
    }

    #[test]
    fn test_zero_inverse() {
        assert_eq!(SF::zero().inverse_or_zero(), SF::zero());
        assert_eq!(BF::zero().inverse_or_zero(), BF::zero());
        assert_eq!(LF::zero().inverse_or_zero(), LF::zero());
    }

    fn test_half_for<const Q: u32>() {
        let two_inv = FpField::<Q>::new(2).inverse_or_zero();
        let mut rng = rng();
        for _ in 0..100 {
            let x: FpField<Q> = rng.random();
            assert_eq!(x.half(), x * two_inv, "Q={Q}");
        }
    }

    #[test]
    fn test_half() {
        test_half_for::<12289>();
        test_half_for::<1073754113>();
        test_half_for::<4294967291>();
    }
}
