use std::fmt::{self, Display, Formatter};
use std::iter::Sum;
use std::mem::size_of;
use std::ops::{
    Add, AddAssign, BitAnd, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Shl, Shr, Sub,
    SubAssign,
};

use num::traits::{ConstOne, ConstZero};
use num::{One, Zero};

/// Constraint trait for signed integer backing types used in `FixedPoint<T>`.
/// Provides type-specific constants and I/O conversions; all arithmetic
/// algorithms are implemented generically in `FixedPoint<T>`.
pub(crate) trait FixedInt:
    Copy
    + Clone
    + PartialEq
    + Eq
    + PartialOrd
    + Ord
    + fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + BitAnd<Output = Self>
    + ConstZero
    + ConstOne
{
    /// Raw representation of 1.0: `1 << FRAC_BITS` where `FRAC_BITS = 4 * size_of::<Self>()`.
    const FRAC_UNIT: Self;
    /// Raw representation of ln(2): `round(ln(2) × 2^FRAC_BITS)`.
    const LN_2_RAW: Self;

    fn from_i64(v: i64) -> Self;
    fn from_f64(v: f64) -> Self;
    fn to_f64(self) -> f64;
    /// Reinterpret self (which must be non-negative) as u128 for sqrt.
    fn raw_to_u128(self) -> u128;
    fn raw_from_u128(v: u128) -> Self;
    fn wrapping_abs(self) -> Self;
    fn trunc_to_i64(self) -> i64;

    /// Clear the lower FRAC_BITS bits (floor toward −∞).
    /// Compute `a_raw * b_raw / 2^FRAC_BITS` using whatever widening is available.
    fn widen_mul(a: Self, b: Self) -> Self;
    /// Compute `a_raw * 2^FRAC_BITS / b_raw` using whatever widening is available.
    fn widen_div(a: Self, b: Self) -> Self;

    fn floor_bits(self) -> Self {
        self & (Self::ZERO - Self::FRAC_UNIT)
    }

    /// Add half the scale then floor (round half toward +∞).
    fn round_bits(self) -> Self {
        (self + (Self::FRAC_UNIT >> 1)).floor_bits()
    }

    fn is_nonneg(self) -> bool {
        self >= Self::ZERO
    }
}

impl FixedInt for i64 {
    const FRAC_UNIT: i64 = 1i64 << 32;
    const LN_2_RAW: i64 = 2_977_044_472i64; // round(ln(2) × 2³²)

    fn from_i64(v: i64) -> i64 {
        v << 32
    }
    fn from_f64(v: f64) -> i64 {
        (v * 4_294_967_296.0).round() as i64
    }
    fn to_f64(self) -> f64 {
        self as f64 / 4_294_967_296.0
    }
    fn raw_to_u128(self) -> u128 {
        self as u64 as u128
    }
    fn raw_from_u128(v: u128) -> i64 {
        v as i64
    }
    fn wrapping_abs(self) -> i64 {
        i64::wrapping_abs(self)
    }
    fn trunc_to_i64(self) -> i64 {
        self / (1i64 << 32)
    }
    fn widen_mul(a: i64, b: i64) -> i64 {
        (((a as i128) * (b as i128)) >> 32) as i64
    }
    fn widen_div(a: i64, b: i64) -> i64 {
        (((a as i128) << 32) / (b as i128)) as i64
    }
}

impl FixedInt for i128 {
    const FRAC_UNIT: i128 = 1i128 << 64;
    const LN_2_RAW: i128 = 12_786_308_645_202_655_660i128; // round(ln(2) × 2⁶⁴)

    fn from_i64(v: i64) -> i128 {
        (v as i128) << 64
    }
    fn from_f64(v: f64) -> i128 {
        (v * (1u128 << 64) as f64).round() as i128
    }
    fn to_f64(self) -> f64 {
        self as f64 / (1u128 << 64) as f64
    }
    fn raw_to_u128(self) -> u128 {
        self as u128
    }
    fn raw_from_u128(v: u128) -> i128 {
        v as i128
    }
    fn wrapping_abs(self) -> i128 {
        i128::wrapping_abs(self)
    }
    fn trunc_to_i64(self) -> i64 {
        (self / (1i128 << 64)) as i64
    }
    fn widen_mul(a: i128, b: i128) -> i128 {
        // Full schoolbook with sign handling. No native i256, so we compute
        // |a| * |b| / 2^64 in u128 (four terms including lo*lo), then fix sign.
        let a_neg = a < 0;
        let b_neg = b < 0;
        let abs_a: u128 = a.unsigned_abs();
        let abs_b: u128 = b.unsigned_abs();
        let half = 64u32;
        let mask: u128 = u64::MAX as u128;
        let a_hi = abs_a >> half;
        let a_lo = abs_a & mask;
        let b_hi = abs_b >> half;
        let b_lo = abs_b & mask;
        // Wrapping arithmetic: intermediates may wrap in u128 but the final
        // result is correct modulo 2^128 — valid as long as the true result
        // fits in i128 (always true for Falcon values bounded by ~2^28).
        let result_u = (a_hi.wrapping_mul(b_hi) << half)
            .wrapping_add(a_hi.wrapping_mul(b_lo))
            .wrapping_add(a_lo.wrapping_mul(b_hi))
            .wrapping_add(a_lo.wrapping_mul(b_lo) >> half);
        let result_s = result_u as i128;
        if a_neg ^ b_neg { -result_s } else { result_s }
    }
    fn widen_div(a: i128, b: i128) -> i128 {
        let q = 32usize;
        let a1 = a << q;
        let d1 = a1 / b;
        let r1 = a1 % b;
        (d1 << q) + (r1 << q) / b
    }
}

/// Fixed-point number: `T` holds the raw integer, scale is `2^FRAC_BITS`
/// where `FRAC_BITS = 4 * size_of::<T>()` (half the type width).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct FixedPoint<T>(pub(crate) T);

pub(crate) type FixedPoint64 = FixedPoint<i64>;
pub(crate) type FixedPoint128 = FixedPoint<i128>;

fn isqrt_u128(n: u128) -> u128 {
    if n == 0 {
        return 0;
    }
    let bits = 128 - n.leading_zeros();
    let mut x = 1u128 << (bits / 2 + 1);
    loop {
        let x1 = (x + n / x) / 2;
        if x1 >= x {
            break;
        }
        x = x1;
    }
    x
}

impl<T: FixedInt> FixedPoint<T> {
    pub const ZERO: Self = FixedPoint(T::ZERO);
    pub const ONE: Self = FixedPoint(T::FRAC_UNIT);
    pub const LN_2: Self = FixedPoint(T::LN_2_RAW);

    pub fn floor(self) -> Self {
        FixedPoint(self.0.floor_bits())
    }

    pub fn round(self) -> Self {
        FixedPoint(self.0.round_bits())
    }

    #[allow(dead_code)]
    pub fn abs(self) -> Self {
        FixedPoint(self.0.wrapping_abs())
    }

    pub fn trunc(self) -> i64 {
        self.0.trunc_to_i64()
    }

    /// Integer square root via Newton-Raphson on u128, then fractional correction.
    pub fn sqrt(self) -> Self {
        assert!(self.0.is_nonneg(), "sqrt of negative FixedPoint");
        let n = T::raw_to_u128(self.0);
        if n == 0 {
            return Self::ZERO;
        }
        let frac_half = 2 * size_of::<T>();
        let r = isqrt_u128(n);
        let e = n - r * r;
        let delta = if r > 0 { (e << frac_half) / (2 * r) } else { 0 };
        FixedPoint(T::raw_from_u128((r << frac_half) + delta))
    }
}

/// Convert FixedPoint128 (scale 2⁶⁴) to FixedPoint64 (scale 2³²).
impl From<FixedPoint128> for FixedPoint64 {
    fn from(v: FixedPoint128) -> FixedPoint64 {
        FixedPoint((v.0 >> 32) as i64)
    }
}

// ---------------------------------------------------------------------------
// Arithmetic operators
// ---------------------------------------------------------------------------

impl<T: FixedInt> Add for FixedPoint<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FixedPoint(self.0 + rhs.0)
    }
}

impl<T: FixedInt> AddAssign for FixedPoint<T> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<T: FixedInt> Sub for FixedPoint<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FixedPoint(self.0 - rhs.0)
    }
}

impl<T: FixedInt> SubAssign for FixedPoint<T> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<T: FixedInt> Neg for FixedPoint<T> {
    type Output = Self;
    fn neg(self) -> Self {
        FixedPoint(-self.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<T: FixedInt> Mul for FixedPoint<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FixedPoint(T::widen_mul(self.0, rhs.0))
    }
}

impl<T: FixedInt> MulAssign for FixedPoint<T> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<T: FixedInt> Div for FixedPoint<T> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        FixedPoint(T::widen_div(self.0, rhs.0))
    }
}

impl<T: FixedInt> DivAssign for FixedPoint<T> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

/// Floor-division remainder: result has the same sign as the divisor.
#[allow(clippy::suspicious_arithmetic_impl)]
impl<T: FixedInt> Rem for FixedPoint<T> {
    type Output = Self;
    fn rem(self, rhs: Self) -> Self {
        self - (self / rhs).floor() * rhs
    }
}

impl<T: FixedInt> RemAssign for FixedPoint<T> {
    fn rem_assign(&mut self, rhs: Self) {
        *self = *self % rhs;
    }
}

// ---------------------------------------------------------------------------
// num traits
// ---------------------------------------------------------------------------

impl<T: FixedInt> Zero for FixedPoint<T> {
    fn zero() -> Self {
        Self::ZERO
    }
    fn is_zero(&self) -> bool {
        self.0 == T::ZERO
    }
}

impl<T: FixedInt> One for FixedPoint<T> {
    fn one() -> Self {
        Self::ONE
    }
}

impl<T: FixedInt> num::Num for FixedPoint<T> {
    type FromStrRadixErr = &'static str;
    fn from_str_radix(_str: &str, _radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Err("FixedPoint::from_str_radix not supported")
    }
}

impl<T: FixedInt> Display for FixedPoint<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0.to_f64())
    }
}

// ---------------------------------------------------------------------------
// Conversions
// ---------------------------------------------------------------------------

impl<T: FixedInt> From<i16> for FixedPoint<T> {
    fn from(v: i16) -> Self {
        FixedPoint(T::from_i64(v as i64))
    }
}

impl<T: FixedInt> From<i32> for FixedPoint<T> {
    fn from(v: i32) -> Self {
        FixedPoint(T::from_i64(v as i64))
    }
}

impl<T: FixedInt> From<i64> for FixedPoint<T> {
    fn from(v: i64) -> Self {
        FixedPoint(T::from_i64(v))
    }
}

impl<T: FixedInt> From<f64> for FixedPoint<T> {
    fn from(v: f64) -> Self {
        FixedPoint(T::from_f64(v))
    }
}

impl<T: FixedInt> From<FixedPoint<T>> for f64 {
    fn from(v: FixedPoint<T>) -> f64 {
        v.0.to_f64()
    }
}

impl<T: FixedInt> Sum for FixedPoint<T> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::ZERO, Add::add)
    }
}

impl<'a, T: FixedInt> Sum<&'a FixedPoint<T>> for FixedPoint<T> {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        iter.copied().sum()
    }
}

// ---------------------------------------------------------------------------
// Tests (all use FixedPoint64 = FixedPoint<i64>)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod test {
    use proptest::prop_assert;
    use proptest::prop_assert_eq;
    use test_strategy::proptest as strategy_proptest;

    use super::{FixedPoint128, FixedPoint64};

    const ZERO: FixedPoint64 = FixedPoint64::ZERO;
    const ONE: FixedPoint64 = FixedPoint64::ONE;

    fn fp(v: i32) -> FixedPoint64 {
        FixedPoint64::from(v)
    }

    // ── deterministic ────────────────────────────────────────────────────────

    #[test]
    fn constants() {
        assert_eq!(f64::from(ZERO), 0.0);
        assert_eq!(f64::from(ONE), 1.0);
        assert!((f64::from(FixedPoint64::LN_2) - std::f64::consts::LN_2).abs() < 1e-9);
    }

    #[test]
    fn constants_fp128() {
        assert_eq!(f64::from(FixedPoint128::ZERO), 0.0);
        assert_eq!(f64::from(FixedPoint128::ONE), 1.0);
        assert!((f64::from(FixedPoint128::LN_2) - std::f64::consts::LN_2).abs() < 1e-9);
    }

    #[test]
    #[should_panic]
    fn div_by_zero() {
        let _ = ONE / ZERO;
    }

    #[test]
    #[should_panic]
    fn sqrt_negative() {
        let _ = fp(-1).sqrt();
    }

    #[test]
    fn sqrt_exact_cases() {
        assert_eq!(fp(4).sqrt(), fp(2));
        assert_eq!(fp(9).sqrt(), fp(3));
        assert_eq!(
            FixedPoint64::from(0.25f64).sqrt(),
            FixedPoint64::from(0.5f64)
        );
        assert_eq!(ZERO.sqrt(), ZERO);
        assert_eq!(ONE.sqrt(), ONE);
    }

    #[test]
    fn floor_negative_edge_cases() {
        assert_eq!(FixedPoint64::from(-1.5f64).floor(), fp(-2));
        assert_eq!(FixedPoint64::from(-0.0001f64).floor(), fp(-1));
        assert_eq!(FixedPoint64::from(-1.0f64).floor(), fp(-1));
    }

    #[test]
    fn round_half_up() {
        assert_eq!(FixedPoint64::from(1.5f64).round(), fp(2));
        assert_eq!(FixedPoint64::from(-1.5f64).round(), fp(-1));
    }

    #[test]
    fn fp128_basic_arithmetic() {
        let a = FixedPoint128::from(3i32);
        let b = FixedPoint128::from(4i32);
        assert_eq!(f64::from(a + b), 7.0);
        assert_eq!(f64::from(a * b), 12.0);
        let inv = FixedPoint128::ONE / b;
        let diff = (f64::from(inv) - 0.25).abs();
        assert!(diff < 1e-9, "1/4 = {}", f64::from(inv));
    }

    #[test]
    fn fp128_from_fp64() {
        let a = FixedPoint128::from(42i32);
        let b = FixedPoint64::from(a);
        assert_eq!(b, FixedPoint64::from(42i32));
    }

    // ── proptests ────────────────────────────────────────────────────────────

    #[strategy_proptest]
    fn add_commutative(
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] a: i32,
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] b: i32,
    ) {
        prop_assert_eq!(fp(a) + fp(b), fp(b) + fp(a));
    }

    #[strategy_proptest]
    fn add_associative(
        #[strategy(-0x1fff_ffffi32..=0x1fff_ffffi32)] a: i32,
        #[strategy(-0x1fff_ffffi32..=0x1fff_ffffi32)] b: i32,
        #[strategy(-0x1fff_ffffi32..=0x1fff_ffffi32)] c: i32,
    ) {
        prop_assert_eq!((fp(a) + fp(b)) + fp(c), fp(a) + (fp(b) + fp(c)));
    }

    #[strategy_proptest]
    fn add_identity(#[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] a: i32) {
        prop_assert_eq!(fp(a) + ZERO, fp(a));
    }

    #[strategy_proptest]
    fn sub_is_neg_add(
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] a: i32,
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] b: i32,
    ) {
        prop_assert_eq!(fp(a) - fp(b), fp(a) + (-fp(b)));
    }

    #[strategy_proptest]
    fn neg_involution(#[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] a: i32) {
        prop_assert_eq!(-(-fp(a)), fp(a));
    }

    #[strategy_proptest]
    fn mul_commutative(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(-0x7fffi32..=0x7fffi32)] b: i32,
    ) {
        prop_assert_eq!(fp(a) * fp(b), fp(b) * fp(a));
    }

    #[strategy_proptest]
    fn mul_identity(#[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] a: i32) {
        prop_assert_eq!(fp(a) * ONE, fp(a));
    }

    #[strategy_proptest]
    fn mul_integer_exact(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(-0x7fffi32..=0x7fffi32)] b: i32,
    ) {
        prop_assert_eq!(fp(a) * fp(b), fp(a * b));
    }

    #[strategy_proptest]
    fn div_inverse_of_mul(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(1i32..=0x7fffi32)] b: i32,
    ) {
        prop_assert_eq!((fp(a) * fp(b)) / fp(b), fp(a));
    }

    #[strategy_proptest]
    fn floor_matches_f64(#[strategy(-1_000_000_000f64..1_000_000_000f64)] v: f64) {
        prop_assert_eq!(FixedPoint64::from(v).floor(), FixedPoint64::from(v.floor()));
    }

    #[strategy_proptest]
    fn floor_idempotent(#[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] a: i32) {
        let x = fp(a);
        prop_assert_eq!(x.floor(), x);
    }

    #[strategy_proptest]
    fn round_result_is_integer(#[strategy(-1_000_000_000f64..1_000_000_000f64)] v: f64) {
        let r = FixedPoint64::from(v).round();
        prop_assert_eq!(r.0 & 0xFFFF_FFFF, 0);
    }

    #[strategy_proptest]
    fn abs_non_negative(#[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] a: i32) {
        prop_assert!(fp(a).abs() >= ZERO);
    }

    #[strategy_proptest]
    fn abs_equals_neg_for_negative(#[strategy(-0x7fff_ffffi32..=-1i32)] a: i32) {
        prop_assert_eq!(fp(a).abs(), -fp(a));
    }

    #[strategy_proptest]
    fn sqrt_vs_f64(#[strategy(0.0f64..1_000_000_000f64)] v: f64) {
        let result = f64::from(FixedPoint64::from(v).sqrt());
        prop_assert!((result - v.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn trunc_exact_cases() {
        assert_eq!(fp(3).trunc(), 3);
        assert_eq!(fp(-3).trunc(), -3);
        assert_eq!(FixedPoint64::from(1.9f64).trunc(), 1);
        assert_eq!(FixedPoint64::from(-1.9f64).trunc(), -1);
        assert_eq!(FixedPoint64::from(-1.0f64).trunc(), -1);
        assert_eq!(ZERO.trunc(), 0);
    }

    #[strategy_proptest]
    fn trunc_matches_floor_for_nonneg(#[strategy(0f64..1_000_000_000f64)] v: f64) {
        prop_assert_eq!(FixedPoint64::from(v).trunc(), v.trunc() as i64);
    }

    #[strategy_proptest]
    fn trunc_toward_zero(#[strategy(-1_000_000_000f64..1_000_000_000f64)] v: f64) {
        let t = FixedPoint64::from(v).trunc();
        prop_assert_eq!(t, v.trunc() as i64);
    }

    #[test]
    fn rem_positive() {
        assert_eq!(
            FixedPoint64::from(7.0f64) % FixedPoint64::from(3.0f64),
            FixedPoint64::from(1.0f64)
        );
    }

    #[test]
    fn rem_negative_dividend() {
        // -7 mod 3 with floor-division semantics = 2
        let result = FixedPoint64::from(-7.0f64) % FixedPoint64::from(3.0f64);
        let expected = FixedPoint64::from(2.0f64);
        let diff = (f64::from(result) - f64::from(expected)).abs();
        assert!(diff < 1e-9, "got {result}, expected {expected}");
    }

    #[strategy_proptest]
    fn rem_sign_matches_divisor(
        #[strategy(-100.0f64..100.0f64)] a: f64,
        #[strategy(1.0f64..100.0f64)] b: f64,
    ) {
        let r = f64::from(FixedPoint64::from(a) % FixedPoint64::from(b));
        prop_assert!(r >= -1e-9 && r < b + 1e-9);
    }

    #[strategy_proptest]
    fn rem_assign_matches_rem(
        #[strategy(-100.0f64..100.0f64)] a: f64,
        #[strategy(1.0f64..100.0f64)] b: f64,
    ) {
        let mut x = FixedPoint64::from(a);
        x %= FixedPoint64::from(b);
        prop_assert_eq!(x, FixedPoint64::from(a) % FixedPoint64::from(b));
    }

    #[strategy_proptest]
    fn sum_equals_fold_add(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(-0x7fffi32..=0x7fffi32)] b: i32,
        #[strategy(-0x7fffi32..=0x7fffi32)] c: i32,
    ) {
        let v = [fp(a), fp(b), fp(c)];
        let s: FixedPoint64 = v.iter().copied().sum();
        prop_assert_eq!(s, fp(a) + fp(b) + fp(c));
    }

    #[strategy_proptest]
    fn sum_ref_equals_sum(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(-0x7fffi32..=0x7fffi32)] b: i32,
    ) {
        let v = vec![fp(a), fp(b)];
        let by_ref: FixedPoint64 = v.iter().sum();
        let by_val: FixedPoint64 = v.into_iter().sum();
        prop_assert_eq!(by_ref, by_val);
    }

    #[strategy_proptest]
    fn ordering_consistent_with_integers(
        #[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] a: i32,
        #[strategy(-0x7fff_ffffi32..=0x7fff_ffffi32)] b: i32,
    ) {
        prop_assert_eq!(fp(a) < fp(b), a < b);
        prop_assert_eq!(fp(a) == fp(b), a == b);
    }

    #[strategy_proptest]
    fn add_assign_matches_add(
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] a: i32,
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] b: i32,
    ) {
        let mut x = fp(a);
        x += fp(b);
        prop_assert_eq!(x, fp(a) + fp(b));
    }

    #[strategy_proptest]
    fn sub_assign_matches_sub(
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] a: i32,
        #[strategy(-0x3fff_ffffi32..=0x3fff_ffffi32)] b: i32,
    ) {
        let mut x = fp(a);
        x -= fp(b);
        prop_assert_eq!(x, fp(a) - fp(b));
    }

    #[strategy_proptest]
    fn mul_assign_matches_mul(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(-0x7fffi32..=0x7fffi32)] b: i32,
    ) {
        let mut x = fp(a);
        x *= fp(b);
        prop_assert_eq!(x, fp(a) * fp(b));
    }

    #[strategy_proptest]
    fn div_assign_matches_div(
        #[strategy(-0x7fffi32..=0x7fffi32)] a: i32,
        #[strategy(1i32..=0x7fffi32)] b: i32,
    ) {
        let mut x = fp(a);
        x /= fp(b);
        prop_assert_eq!(x, fp(a) / fp(b));
    }
}
