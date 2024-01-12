use std::ops::MulAssign;

use num::{One, Zero};

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
