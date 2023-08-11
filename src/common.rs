use rand_distr::num_traits::Zero;

use crate::field::Felt;

/// This file contains methods and objects which are reused through multiple files.
///

/// Split a polynomial f into two polynomials
/// of N coefficients, using even-off split.
///
///  Args:
///      f: a polynomial
///  Format: coefficients
pub(crate) fn split(f: &[Felt]) -> (Vec<Felt>, Vec<Felt>) {
    let n = f.len();
    let mut f0 = vec![Felt::zero(); n / 2];
    let mut f1 = vec![Felt::zero(); n / 2];
    for i in 0..n / 2 {
        f0[i] = f[2 * i];
        f1[i] = f[2 * i + 1];
    }
    (f0, f1)
}

/// Merge two polynomials into a single
/// polynomial f, using even-odd split.
/// Args:
///     f_list: a pair of polynomials
/// Format: coefficients
pub(crate) fn merge(f0: &[Felt], f1: &[Felt]) -> Vec<Felt> {
    let n = 2 * f0.len();
    let mut f = vec![Felt::zero(); n];
    for i in 0..n / 2 {
        f[2 * i + 0] = f0[i];
        f[2 * i + 1] = f1[i];
    }
    f
}

/// Compute the square euclidean norm of the vector v.
pub(crate) fn sqnorm(v: &[Vec<Felt>]) -> u32 {
    let mut res = 0u32;
    for elt in v {
        for coeff in elt {
            res += (coeff.0 as u32) * (coeff.0 as u32);
        }
    }
    res
}
