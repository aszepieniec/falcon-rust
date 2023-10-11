use num_complex::Complex64;
use rand_distr::num_traits::Zero;

use crate::{
    common::{merge, split},
    fft_constants::roots_dict_complex,
};

/// Split a polynomial into even and odd halves, except the polynomial
/// is given in FFT representation.
pub fn split_fft(f_fft: &[Complex64]) -> (Vec<Complex64>, Vec<Complex64>) {
    let n = f_fft.len();
    let w = roots_dict_complex(n);
    let mut f0_fft = vec![Complex64::zero(); n / 2];
    let mut f1_fft = vec![Complex64::zero(); n / 2];
    for i in 0..(n / 2) {
        f0_fft[i] = 0.5 * (f_fft[2 * i] + f_fft[2 * i + 1]);
        f1_fft[i] = 0.5 * (f_fft[2 * i] - f_fft[2 * i + 1]) * w[2 * i].conj();
    }
    (f0_fft, f1_fft)
}

/// Merge two polynomials into one by interleaving their coefficients,
/// except the polynomials are given in FFT representation.
pub fn merge_fft(f0_fft: &[Complex64], f1_fft: &[Complex64]) -> Vec<Complex64> {
    let n = 2 * f0_fft.len();
    let w = roots_dict_complex(n);
    let mut f_fft = vec![Complex64::zero(); n];
    for i in 0..(n / 2) {
        f_fft[2 * i + 0] = f0_fft[i] + w[2 * i] * f1_fft[i];
        f_fft[2 * i + 1] = f0_fft[i] - w[2 * i] * f1_fft[i];
    }
    f_fft
}

/// Compute the discrete Fourier transform of the given polynomial
/// using the FFT algorithm (or rather, a variant of the FFT family).
pub fn fft(f: &[Complex64]) -> Vec<Complex64> {
    let n = f.len();
    if n > 2 {
        let (f0, f1) = split(f);
        let f0_fft = fft(&f0);
        let f1_fft = fft(&f1);
        merge_fft(&f0_fft, &f1_fft)
    } else {
        let mut f_fft = vec![Complex64::zero(); n];
        f_fft[0] = f[0] + Complex64::i() * f[1];
        f_fft[1] = f[0] - Complex64::i() * f[1];
        f_fft
    }
}

/// Compute the inverse discrete Fourier transform of the given
/// polynomial using the FFT algorithm (the logical inverse of
/// the previous one).
fn ifft(f_fft: &[Complex64]) -> Vec<Complex64> {
    let n = f_fft.len();
    if n > 2 {
        let (f0_fft, f1_fft) = split_fft(f_fft);
        let f0 = ifft(&f0_fft);
        let f1 = ifft(&f1_fft);
        merge(&f0, &f1)
    } else {
        let mut f = vec![Complex64::zero(); 2];
        f[0] = 0.5 * (f_fft[0] + f_fft[1]);
        f[1] = 0.5 * (f_fft[0] - f_fft[1]) * Complex64::i().conj();
        f
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num_complex::Complex64;
    use rand::{thread_rng, Rng};

    use crate::fft::{fft, ifft};

    fn assert_approx_equal(a: &[Complex64], b: &[Complex64]) {
        let norm = a
            .iter()
            .zip(b.iter())
            .map(|(x, y)| x - y)
            .map(|x| x * x.conj())
            .sum::<Complex64>();
        assert!(
            norm.re <= 1e-10,
            "squared norm of difference: {}\n\nlhs: {:?}\n\nrhs: {:?}",
            norm.re,
            a,
            b
        );
    }

    #[test]
    fn test_fft_linearity() {
        let mut rng = thread_rng();
        let a: Complex64 = Complex64::new(rng.gen(), rng.gen());
        let b: Complex64 = Complex64::new(rng.gen(), rng.gen());
        const N: usize = 256;
        let v: [Complex64; N] = (0..N)
            .map(|_| Complex64::new(rng.gen(), rng.gen()))
            .collect_vec()
            .try_into()
            .unwrap();
        let w: [Complex64; N] = (0..N)
            .map(|_| Complex64::new(rng.gen(), rng.gen()))
            .collect_vec()
            .try_into()
            .unwrap();

        let v_ntt = fft(&v);
        let w_ntt = fft(&w);
        let lhs = v_ntt
            .into_iter()
            .zip(w_ntt.into_iter())
            .map(|(vv, ww)| a * vv + b * ww)
            .collect_vec();

        let rhs_coeffs = v
            .into_iter()
            .zip(w.into_iter())
            .map(|(vv, ww)| a * vv + b * ww)
            .collect_vec();
        let rhs = fft(&rhs_coeffs);

        assert_approx_equal(&lhs, &rhs);
    }

    #[test]
    fn test_fft_inverse() {
        let mut rng = thread_rng();
        const N: usize = 1024;
        let v: [Complex64; N] = (0..N)
            .map(|_| Complex64::new(rng.gen(), rng.gen()))
            .collect_vec()
            .try_into()
            .unwrap();

        let v_fft = fft(&v);
        let v_again = ifft(&v_fft);
        let v_fft_again = fft(&v_again);

        assert_approx_equal(&v.to_vec(), &v_again);
        assert_approx_equal(&v_fft, &v_fft_again);
    }
}
