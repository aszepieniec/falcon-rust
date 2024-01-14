use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{One, Zero};
use num_complex::Complex64;

use crate::{cyclotomic_fourier::CyclotomicFourier, inverse::Inverse, polynomial::Polynomial};

pub trait FftVector {
    type Field: Add + Mul + AddAssign + MulAssign + Neg + Sub + SubAssign + One + Zero;
    fn fft(&mut self);
    fn ifft(&mut self);
}

impl FftVector for Polynomial<Complex64> {
    type Field = Complex64;
    fn fft(&mut self) {
        let n = self.coefficients.len();
        let psi_rev = Complex64::bitreversed_powers(n);
        Complex64::fft(&mut self.coefficients, &psi_rev);
    }

    fn ifft(&mut self) {
        let n = self.coefficients.len();
        let psi_inv_rev = Complex64::bitreversed_powers_inverse(n);
        let ninv = Complex64::inverse_or_zero(Complex64::new(n as f64, 0.0));
        Complex64::ifft(&mut self.coefficients, &psi_inv_rev, ninv);
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use num_complex::Complex64;
    use rand::{thread_rng, Rng};

    use crate::polynomial::Polynomial;

    use super::FftVector;

    #[test]
    fn test_fft() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| Complex64::new(rng.gen_range(0..2) as f64, rng.gen_range(0..2) as f64))
            .collect_vec();
        // let mut a = vec![Complex64::zero(), Complex64::one()];
        let mut b = crate::fft::fft(&a);
        let mut a_fft = Polynomial::new(a);
        a_fft.fft();
        a = a_fft.coefficients;

        // a and b should be the same but in different orders
        // so check their sums of 2^nth roots instead (more stable)

        let mut sums_a = vec![];
        let mut sums_b = vec![];
        for _ in 0..n {
            sums_a.push(a.iter().cloned().sum::<Complex64>());
            sums_b.push(b.iter().cloned().sum::<Complex64>());

            a.iter_mut().for_each(|e| {
                e.sqrt();
            });
            b.iter_mut().for_each(|e| {
                e.sqrt();
            });
        }

        assert!(sums_a
            .into_iter()
            .zip(sums_b.into_iter())
            .all(|(sum_a, sum_b)| (sum_a - sum_b).norm() < 100.0 * f32::EPSILON as f64));
    }
}
