use rand_distr::num_traits::Zero;

/// This file contains methods and objects which are reused through multiple files.
///

/// Split a polynomial f into two polynomials using even-odd split.
pub(crate) fn split<T: Zero + Copy>(f: &[T]) -> (Vec<T>, Vec<T>) {
    let n = f.len();
    let mut f0 = vec![T::zero(); n / 2];
    let mut f1 = vec![T::zero(); n / 2];
    for i in 0..n / 2 {
        f0[i] = f[2 * i];
        f1[i] = f[2 * i + 1];
    }
    (f0, f1)
}

/// Merge two polynomials into a single polynomial f, using even-odd
/// split.
pub(crate) fn merge<T: Zero + Copy>(f0: &[T], f1: &[T]) -> Vec<T> {
    let n = 2 * f0.len();
    let mut f = vec![T::zero(); n];
    for i in 0..n / 2 {
        f[2 * i] = f0[i];
        f[2 * i + 1] = f1[i];
    }
    f
}
