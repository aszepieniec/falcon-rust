use num_complex::Complex64;
use rand_distr::num_traits::Zero;

use crate::{
    common::{merge, split},
    fft_constants::roots_dict_complex,
};

/// Split a polynomial into even and odd halves, except the polynomial
/// is given in FFT representation.
pub(crate) fn split_fft(f_fft: &[Complex64]) -> (Vec<Complex64>, Vec<Complex64>) {
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
pub(crate) fn merge_fft(f0_fft: &[Complex64], f1_fft: &[Complex64]) -> Vec<Complex64> {
    let n = 2 * f0_fft.len();
    let w = roots_dict_complex(n);
    let mut f_fft = vec![Complex64::zero(); n];
    for i in 0..(n / 2) {
        f_fft[2 * i] = f0_fft[i] + w[2 * i] * f1_fft[i];
        f_fft[2 * i + 1] = f0_fft[i] - w[2 * i] * f1_fft[i];
    }
    f_fft
}

/// Compute the discrete Fourier transform of the given polynomial
/// using the FFT algorithm (or rather, a variant of the FFT family).
pub(crate) fn fft(f: &[Complex64]) -> Vec<Complex64> {
    let n = f.len();
    match n {
        1 => f.to_vec(),
        2 => {
            let mut f_fft = vec![Complex64::zero(); n];
            f_fft[0] = f[0] + Complex64::i() * f[1];
            f_fft[1] = f[0] - Complex64::i() * f[1];
            f_fft
        }
        _ => {
            let (f0, f1) = split(f);
            let f0_fft = fft(&f0);
            let f1_fft = fft(&f1);
            merge_fft(&f0_fft, &f1_fft)
        }
    }
}

/// Compute the inverse discrete Fourier transform of the given
/// polynomial using the FFT algorithm (the logical inverse of
/// the previous one).
pub(crate) fn ifft(f_fft: &[Complex64]) -> Vec<Complex64> {
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
