use std::{
    f64::consts::PI,
    ops::{Add, Mul, MulAssign, Sub},
};

use num::{One, Zero};
use num_complex::Complex64;

use crate::inverse::Inverse;

pub(crate) trait CyclotomicFourier
where
    Self: Sized
        + Copy
        + One
        + Zero
        + Add<Output = Self>
        + Sub<Output = Self>
        + Mul<Output = Self>
        + MulAssign
        + Inverse,
{
    /// Get the inverse of 2^n.
    fn power_of_two_inverse(n: usize) -> Self {
        let mut a = Self::one() + Self::one();
        for _ in 0..n {
            a *= a;
        }
        Self::inverse_or_zero(a)
    }

    /// Get a primitive nth (with n a power of 2) root of unity.
    fn primitive_root_of_unity(n: usize) -> Self;

    /// Compute the integer whose n-bit binary expansion is the reverse of
    /// that of the argument.
    fn bitreverse_index(arg: usize, n: usize) -> usize {
        assert!(n > 0);
        assert_eq!(n & (n - 1), 0);
        let mut rev = 0;
        let mut m = n >> 1;
        let mut k = 1;
        while m > 0 {
            rev |= (((arg & m) != 0) as usize) * k;
            k <<= 1;
            m >>= 1;
        }
        rev
    }

    /// Compute the first n powers of the 2nth root of unity, and put them in
    /// bit-reversed order.
    fn bitreversed_powers(n: usize) -> Vec<Self> {
        let psi = Self::primitive_root_of_unity(2 * n);
        let mut array = vec![Self::zero(); n];
        let mut alpha = Self::one();
        for a in array.iter_mut() {
            *a = alpha;
            alpha *= psi;
        }
        Self::bitreverse_array(&mut array);
        array
    }

    /// Compute the first n powers of the 2nth root of unity, invert them, and
    /// put them in bit-reversed order.
    fn bitreversed_powers_inverse(n: usize) -> Vec<Self> {
        let psi = Self::primitive_root_of_unity(2 * n).inverse_or_zero();
        let mut array = vec![Self::zero(); n];
        let mut alpha = Self::one();
        for a in array.iter_mut() {
            *a = alpha;
            alpha *= psi;
        }
        Self::bitreverse_array(&mut array);
        array
    }

    /// Reorder the given elements in the array by reversing the binary
    /// expansions of their indices.
    fn bitreverse_array<T>(array: &mut [T]) {
        let n = array.len();
        for i in 0..n {
            let j = Self::bitreverse_index(i, n);
            if i < j {
                array.swap(i, j);
            }
        }
    }

    /// Compute the evaluations of the polynomial on the roots of the
    /// polynomial X^n + 1 using a fast Fourier transform. Algorithm 1
    /// from https://eprint.iacr.org/2016/504.pdf.
    ///
    /// Arguments:
    ///
    ///  - a : &mut [Self]
    ///    (a reference to) a mutable array of field elements which is to
    ///    be transformed under the FFT. The transformation happens in-
    ///    place.
    ///
    ///  - psi_rev: &[Self]
    ///    (a reference to) an array of powers of psi, from 0 to n-1,
    ///    but ordered by bit-reversed index. Here psi is a primitive root
    ///    of order 2n. You can use
    ///    `Self::bitreversed_powers(psi, n)` for this purpose, but this
    ///    trait implementation is not const. For the performance benefit
    ///    you want a precompiled array, which you can get if you can get
    ///    by implementing the same method and marking it "const".
    fn fft(a: &mut [Self], psi_rev: &[Self]) {
        let n = a.len();
        let mut t = n;
        let mut m = 1;
        while m < n {
            t >>= 1;
            for i in 0..m {
                let j1 = 2 * i * t;
                let j2 = j1 + t - 1;
                let s = psi_rev[m + i];
                for j in j1..=j2 {
                    let u = a[j];
                    let v = a[j + t] * s;
                    a[j] = u + v;
                    a[j + t] = u - v;
                }
            }
            m <<= 1;
        }
    }

    /// Compute the coefficients of the polynomial with the given evaluations
    /// on the roots of X^n + 1 using an inverse fast Fourier transform.
    /// Algorithm 2 from https://eprint.iacr.org/2016/504.pdf.
    ///
    /// Arguments:
    ///
    ///  - a : &mut [Self]
    ///    (a reference to) a mutable array of field elements which is to
    ///    be transformed under the IFFT. The transformation happens in-
    ///    place.
    ///
    ///  - psi_inv_rev: &[Self]
    ///    (a reference to) an array of powers of psi^-1, from 0 to n-1,
    ///    but ordered by bit-reversed index. Here psi is a primitive root of
    ///    order 2n. You can use
    ///    `Self::bitreversed_powers(Self::inverse_or_zero(psi), n)` for
    ///    this purpose, but this trait implementation is not const. For
    ///    the performance benefit you want a precompiled array, which you
    ///    can get if you can get by implementing the same methods and marking
    ///    them "const".
    fn ifft(a: &mut [Self], psi_inv_rev: &[Self], ninv: Self) {
        let n = a.len();
        let mut t = 1;
        let mut m = n;
        while m > 1 {
            let h = m / 2;
            let mut j1 = 0;
            for i in 0..h {
                let j2 = j1 + t - 1;
                let s = psi_inv_rev[h + i];
                for j in j1..=j2 {
                    let u = a[j];
                    let v = a[j + t];
                    a[j] = u + v;
                    a[j + t] = (u - v) * s;
                }
                j1 += 2 * t;
            }
            t <<= 1;
            m >>= 1;
        }
        for ai in a.iter_mut() {
            *ai *= ninv;
        }
    }

    fn split_fft(f: &[Self], psi_inv_rev: &[Self]) -> (Vec<Self>, Vec<Self>) {
        let n_over_2 = f.len() / 2;
        let mut f0 = vec![Self::zero(); n_over_2];
        let mut f1 = vec![Self::zero(); n_over_2];
        let two_inv = (Self::one() + Self::one()).inverse_or_zero();
        for i in 0..n_over_2 {
            let two_i = i * 2;
            let two_zeta_inv = two_inv * psi_inv_rev[n_over_2 + i];
            f0[i] = two_inv * (f[two_i] + f[two_i + 1]);
            f1[i] = two_zeta_inv * (f[two_i] - f[two_i + 1]);
        }
        (f0, f1)
    }

    fn merge_fft(f0: &[Self], f1: &[Self], psi_rev: &[Self]) -> Vec<Self> {
        let n_over_2 = f0.len();
        let n = 2 * n_over_2;
        let mut f = vec![Self::zero(); n];
        for i in 0..n_over_2 {
            let two_i = i * 2;
            f[two_i] = f0[i] + psi_rev[n_over_2 + i] * f1[i];
            f[two_i + 1] = f0[i] - psi_rev[n_over_2 + i] * f1[i];
        }
        f
    }
}

impl CyclotomicFourier for Complex64 {
    fn primitive_root_of_unity(n: usize) -> Self {
        let angle = 2. * PI / (n as f64);
        Complex64::new(f64::cos(angle), f64::sin(angle))
    }

    /// Custom implementation of CyclotomicFourier::bitreversed_powers for
    /// better precision.
    fn bitreversed_powers(n: usize) -> Vec<Self> {
        let mut array = vec![Self::zero(); n];
        let half_circle = PI;
        for (i, a) in array.iter_mut().enumerate() {
            let angle = (i as f64) * half_circle / (n as f64);
            *a = Self::new(f64::cos(angle), f64::sin(angle));
        }
        Self::bitreverse_array(&mut array);
        array
    }

    /// Custom implementation of CyclotomicFourier::bitreversed_powers_inverse
    /// for better precision.
    fn bitreversed_powers_inverse(n: usize) -> Vec<Self> {
        let mut array = vec![Self::zero(); n];
        let half_circle = PI;
        for (i, a) in array.iter_mut().enumerate() {
            let angle = (i as f64) * half_circle / (n as f64);
            *a = Self::new(f64::cos(angle), -f64::sin(angle));
        }
        Self::bitreverse_array(&mut array);
        array
    }
}

#[cfg(test)]
mod test {
    use crate::inverse::Inverse;
    use crate::{cyclotomic_fourier::CyclotomicFourier, polynomial::Polynomial};
    use itertools::Itertools;
    use num::One;
    use num_complex::Complex64;
    use rand::{thread_rng, Rng, RngCore};

    fn diff(u: &[Complex64], v: &[Complex64]) -> f64 {
        u.iter()
            .zip(v.iter())
            .map(|(l, r)| l - r)
            .map(|c| c * c.conj())
            .sum::<Complex64>()
            .re
    }

    #[test]
    fn test_primitive_nth_root_of_unity() {
        for log2n in 0..12 {
            let n = 1 << log2n;
            let mut z = Complex64::primitive_root_of_unity(n);
            for _ in 0..log2n {
                assert!((z - Complex64::one()).norm() > f32::EPSILON as f64);
                z *= z;
            }
            assert!((z - Complex64::one()).norm() < f32::EPSILON as f64);
        }
    }

    #[test]
    fn test_fft() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| Complex64::new(rng.gen_range(0..2) as f64, rng.gen_range(0..2) as f64))
            .collect_vec();
        let mut b = a.clone();

        assert!(diff(&a, &b) < 100.0 * f64::EPSILON);

        let psi_rev = Complex64::bitreversed_powers(2 * n);
        let psi_inv_rev = Complex64::bitreversed_powers_inverse(2 * n);
        let ninv = Complex64::inverse_or_zero(Complex64::new(n as f64, 0.0));
        Complex64::fft(&mut a, &psi_rev);
        Complex64::ifft(&mut a, &psi_inv_rev, ninv);
        assert!(
            diff(&a, &b) < f32::EPSILON as f64,
            "a: {:?}\nb: {:?}\nnorm: {}",
            a,
            b,
            diff(&a, &b)
        );

        let x = Complex64::new(rng.next_u32() as f64, rng.next_u32() as f64);
        let y = Complex64::new(rng.next_u32() as f64, rng.next_u32() as f64);
        let mut c = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        Complex64::fft(&mut a, &psi_rev);
        Complex64::fft(&mut b, &psi_rev);
        Complex64::fft(&mut c, &psi_rev);

        let c_alt = a
            .iter()
            .zip(b.iter())
            .map(|(&l, &r)| x * l + y * r)
            .collect_vec();

        assert!(
            diff(&c, &c_alt) < f32::EPSILON as f64,
            "norm of difference: {}",
            diff(&c, &c_alt)
        );
    }

    #[test]
    fn test_multiply_reduce() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| Complex64::new(rng.gen_range(0..2) as f64, 0.0))
            .collect_vec();
        let mut b = (0..n)
            .map(|_| Complex64::new(rng.gen_range(0..2) as f64, 0.0))
            .collect_vec();

        let c = (Polynomial::new(a.clone()) * Polynomial::new(b.clone()))
            .reduce_by_cyclotomic(n)
            .coefficients;

        let psi_rev = Complex64::bitreversed_powers(n);
        Complex64::fft(&mut a, &psi_rev);
        Complex64::fft(&mut b, &psi_rev);
        let mut d = a.iter().zip(b.iter()).map(|(l, r)| l * r).collect_vec();
        let psi_inv_rev = Complex64::bitreversed_powers_inverse(n);
        let ninv = Complex64::new(1.0 / (n as f64), 0.0);
        Complex64::ifft(&mut d, &psi_inv_rev, ninv);

        assert!(
            diff(&c, &d) < f32::EPSILON as f64,
            "lhs: {:?}\nrhs: {:?}\nnorm: {}",
            c,
            d,
            diff(&c, &d)
        );
    }

    #[test]
    fn test_split_fft() {
        let n = 32;
        let mut rng = thread_rng();
        let mut a = (0..n)
            .map(|_| Complex64::new(rng.gen_range(0..2) as f64, 0.0))
            .collect_vec();

        let mut e = a.chunks(2).map(|ch| ch[0]).collect_vec();
        let mut o = a.chunks(2).map(|ch| ch[1]).collect_vec();

        Complex64::fft(&mut a, &Complex64::bitreversed_powers(n));
        let (f0, f1) = Complex64::split_fft(&a, &Complex64::bitreversed_powers_inverse(n));

        Complex64::fft(&mut e, &Complex64::bitreversed_powers(n));
        Complex64::fft(&mut o, &Complex64::bitreversed_powers(n));

        assert!(
            diff(&e, &f0) <= 100.0 * f64::EPSILON,
            "diff: {}",
            diff(&e, &f0)
        );

        assert!(
            diff(&o, &f1) <= 100.0 * f64::EPSILON,
            "diff: {}",
            diff(&o, &f1)
        );
    }

    #[test]
    fn test_merge_fft() {
        let n = 32;
        let mut rng = thread_rng();
        let mut e = (0..n)
            .map(|_| Complex64::new(rng.gen_range(-50..50) as f64, 0.0))
            .collect_vec();
        let mut o = (0..n)
            .map(|_| Complex64::new(rng.gen_range(-50..50) as f64, 0.0))
            .collect_vec();

        let mut ab = Complex64::merge_fft(&e, &o, &Complex64::bitreversed_powers(2 * n));
        Complex64::ifft(
            &mut ab,
            &Complex64::bitreversed_powers_inverse(2 * n),
            Complex64::new(1.0 / (2.0 * n as f64), 0.0),
        );

        Complex64::ifft(
            &mut e,
            &Complex64::bitreversed_powers_inverse(n),
            Complex64::new(1.0 / (n as f64), 0.0),
        );
        Complex64::ifft(
            &mut o,
            &Complex64::bitreversed_powers_inverse(n),
            Complex64::new(1.0 / (n as f64), 0.0),
        );

        let f = e
            .into_iter()
            .zip(o)
            .flat_map(|(ee, oo)| vec![ee, oo])
            .collect_vec();

        assert!(
            diff(&f, &ab) <= 100.0 * f64::EPSILON,
            "diff: {}",
            diff(&f, &ab)
        );
    }
}
