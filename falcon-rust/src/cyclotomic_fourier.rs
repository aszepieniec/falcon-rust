use std::ops::{Add, Mul, MulAssign, Sub};

use num::{One, Zero};

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
        + Inverse
        + From<usize>,
{
    /// Get the inverse of 2^n.
    fn power_of_two_inverse(n: usize) -> Self {
        let mut a = Self::from(2);
        for _ in 0..n {
            a *= a;
        }
        Self::inverse_or_zero(a)
    }

    /// Get a primitive nth (with n a power of 2) root of unity.
    fn primitive_root_of_unity(n: usize) -> Self;

    /// Compute the integer whose n-bit binary expansion is the reverse of
    /// that of the argument.
    fn bitreverse(arg: usize, n: usize) -> usize {
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

    /// Compute the n powers of the nth root of unity psi, and put them in
    /// bit-reversed order.
    fn bitreversed_powers(psi: Self, n: usize) -> Vec<Self> {
        let mut array = vec![Self::from(0); n];
        let mut alpha = Self::from(1);
        for a in array.iter_mut() {
            *a = alpha;
            alpha *= psi;
        }
        for i in 0..n {
            let j = Self::bitreverse(i, n);
            if i < j {
                array.swap(i, j);
            }
        }
        array
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
    ///    but ordered by bit-reversed index. You can use
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
    ///  - psi_rev: &[Self]
    ///    (a reference to) an array of powers of psi^-1, from 0 to n-1,
    ///    but ordered by bit-reversed index. You can use
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
}
