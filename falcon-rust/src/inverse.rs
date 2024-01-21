use std::ops::MulAssign;

use itertools::Itertools;
use num::{One, Zero};
use num_complex::Complex64;

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

#[cfg(test)]
mod test {
    use super::*;
    use rand::{thread_rng, RngCore};

    #[test]
    fn test_complex_inverse() {
        let mut rng = thread_rng();
        let c = Complex64::new(rng.next_u32() as f64, rng.next_u32() as f64);
        let i = c.inverse_or_zero();
        let diff = c * i - Complex64::one();
        let norm = diff.re * diff.re + diff.im * diff.im;
        assert!(norm < f64::EPSILON * 100.0, "norm: {norm}");
    }
}
