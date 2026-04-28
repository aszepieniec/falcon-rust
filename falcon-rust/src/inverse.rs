use std::ops::MulAssign;

use itertools::Itertools;
use num::{One, Zero};
use num_complex::{Complex, Complex64};

use crate::fixed_point::{FixedInt, FixedPoint};

pub(crate) trait Inverse: Copy + Zero + MulAssign + One {
    /// Get the inverse of a, or zero if it is zero.
    fn inverse_or_zero(self) -> Self;

    /// Get the inverses of a batch of elements, and skip over any
    /// that are zero.
    fn batch_inverse_or_zero(batch: &[Self]) -> Vec<Self> {
        let mut acc = Self::one();
        let mut rp: Vec<Self> = Vec::with_capacity(batch.len());
        for batch_item in batch {
            if !batch_item.is_zero() {
                rp.push(acc);
                acc = *batch_item * acc;
            } else {
                rp.push(Self::zero());
            }
        }
        let mut inv = Self::inverse_or_zero(acc);
        for i in (0..batch.len()).rev() {
            if !batch[i].is_zero() {
                rp[i] *= inv;
                inv *= batch[i];
            }
        }
        rp
    }
}

impl Inverse for Complex64 {
    fn inverse_or_zero(self) -> Self {
        let modulus = self.re * self.re + self.im * self.im;
        Complex64::new(self.re / modulus, -self.im / modulus)
    }
    fn batch_inverse_or_zero(batch: &[Self]) -> Vec<Self> {
        batch
            .iter()
            .map(|&c| Complex64::new(1.0, 0.0) / c)
            .collect_vec()
    }
}

impl Inverse for f64 {
    fn inverse_or_zero(self) -> Self {
        1.0 / self
    }
    fn batch_inverse_or_zero(batch: &[Self]) -> Vec<Self> {
        batch.iter().map(|&c| 1.0 / c).collect_vec()
    }
}

impl<T: FixedInt> Inverse for FixedPoint<T> {
    fn inverse_or_zero(self) -> Self {
        if self.is_zero() {
            Self::ZERO
        } else {
            Self::ONE / self
        }
    }
}

impl<T: FixedInt> Inverse for Complex<FixedPoint<T>>
where
    Complex<FixedPoint<T>>: Copy + Zero + MulAssign + One,
{
    fn inverse_or_zero(self) -> Self {
        if self.im.is_zero() {
            // Fast path for purely real values: avoids computing re²
            // whose denominator would overflow for large Gram elements.
            if self.re.is_zero() {
                Complex::zero()
            } else {
                Complex::new(FixedPoint::ONE / self.re, FixedPoint::ZERO)
            }
        } else {
            let m = self.re * self.re + self.im * self.im;
            if m.is_zero() {
                Complex::zero()
            } else {
                Complex::new(self.re / m, -(self.im / m))
            }
        }
    }

    /// Element-wise inversion — the accumulator method in the default impl
    /// overflows when batch values are large (e.g. Gram matrix elements).
    fn batch_inverse_or_zero(batch: &[Self]) -> Vec<Self> {
        batch.iter().map(|&c| c.inverse_or_zero()).collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{rng, RngCore};

    #[test]
    fn test_complex_inverse() {
        let mut rng = rng();
        let c = Complex64::new(rng.next_u32() as f64, rng.next_u32() as f64);
        let i = c.inverse_or_zero();
        let diff = c * i - Complex64::one();
        let norm = diff.re * diff.re + diff.im * diff.im;
        assert!(norm < f64::EPSILON * 100.0, "norm: {norm}");
    }

    #[test]
    fn test_complex_fp_inverse() {
        // c = (3 + 4i): |c|^2 = 25, c^-1 = (3/25 - 4i/25)
        let c = Complex::new(FixedPoint::<i64>::from(3i32), FixedPoint::from(4i32));
        let inv = c.inverse_or_zero();
        let prod = c * inv;
        let re_err = (f64::from(prod.re) - 1.0).abs();
        let im_err = f64::from(prod.im).abs();
        assert!(re_err < 1e-9 && im_err < 1e-9, "prod={prod:?}");
    }

    #[test]
    fn test_fp_inverse_zero() {
        assert_eq!(FixedPoint::<i64>::ZERO.inverse_or_zero(), FixedPoint::ZERO);
    }

    #[test]
    fn test_fp_inverse() {
        let x = FixedPoint::<i64>::from(4i32);
        let inv = x.inverse_or_zero();
        let prod = x * inv;
        let err = (f64::from(prod) - 1.0).abs();
        assert!(err < 1e-9, "prod={prod}");
    }
}
