use itertools::Itertools;
use num_complex::{Complex, Complex64};
use rand::RngCore;
use rand_distr::num_traits::{One, Zero};

use crate::{
    falcon,
    fft::{merge_fft, split_fft},
    samplerz::sampler_z,
};

/// Computes the Gram matrix. The argument must be a 2x2 matrix
/// whose elements are equal-length vectors of complex numbers,
/// representing polynomials in FFT domain.
pub(crate) fn gram(b: [Vec<Complex64>; 4]) -> [Vec<Complex64>; 4] {
    const N: usize = 2;
    let mut g: [Vec<Complex<f64>>; 4] = (0..4)
        .map(|_| (0..b[0].len()).map(|_| Complex64::zero()).collect_vec())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    for i in 0..N {
        for j in 0..N {
            for k in 0..N {
                g[N * i + j] = g[N * i + j]
                    .iter()
                    .zip(
                        b[N * i + k]
                            .iter()
                            .zip(b[N * j + k].iter().map(|c| c.conj()))
                            .map(|(a, b)| a * b),
                    )
                    .map(|(a, b)| a + b)
                    .collect_vec();
            }
        }
    }
    g
}

/// Compute the LDL decomposition of a 2x2 matrix G such that
///     L D L* = G
/// where D is diagonal, and L is lower-triangular. The elements of matrices
/// are in FFT domain.
pub(crate) fn ldl(g: [Vec<Complex64>; 4]) -> ([Vec<Complex64>; 4], [Vec<Complex64>; 4]) {
    let n = g[0].len();

    let zero = (0..n).map(|_| Complex64::zero()).collect_vec();
    let one = (0..n).map(|_| Complex64::one()).collect_vec();

    let l10 = g[2]
        .iter()
        .zip(g[0].iter())
        .map(|(a, b)| a / b)
        .collect_vec();
    let bc = l10.iter().map(|c| c * c.conj());
    let abc = g[0].iter().zip(bc).map(|(a, bc)| a * bc);
    let d11 = g[3]
        .iter()
        .zip(abc)
        .map(|(g11, abc)| g11 - abc)
        .collect_vec();

    let l = [one.clone(), zero.clone(), l10.clone(), one];
    let d = [g[0].clone(), zero.clone(), zero, d11];
    (l, d)
}

#[derive(Debug, Clone)]
pub(crate) enum LdlTree {
    Branch(Vec<Complex64>, Box<LdlTree>, Box<LdlTree>),
    Leaf([Complex64; 2]),
}

/// Compute the LDL Tree of G. Corresponds to Algorithm 9 of the
/// specification [1, p.37]. The argument is a 2x2 matrix of
/// polynomials, given in FFT form.
///
/// [1]: https://falcon-sign.info/falcon.pdf
pub(crate) fn ffldl(g: [Vec<Complex64>; 4]) -> LdlTree {
    let n = g[0].len();
    let (l, d) = ldl(g);

    if n > 2 {
        let (d00, d01) = split_fft(&d[0]);
        let (d10, d11) = split_fft(&d[3]);
        let g0 = [
            d00.clone(),
            d01.clone(),
            d01.iter().map(|c| c.conj()).collect_vec(),
            d00,
        ];
        let g1 = [
            d10.clone(),
            d11.clone(),
            d11.iter().map(|c| c.conj()).collect_vec(),
            d10,
        ];
        LdlTree::Branch(l[2].clone(), Box::new(ffldl(g0)), Box::new(ffldl(g1)))
    } else {
        LdlTree::Branch(
            l[2].clone(),
            Box::new(LdlTree::Leaf(d[0].clone().try_into().unwrap())),
            Box::new(LdlTree::Leaf(d[3].clone().try_into().unwrap())),
        )
    }
}

pub(crate) fn normalize_tree(tree: &mut LdlTree, sigma: f64) {
    match tree {
        LdlTree::Branch(_ell, left, right) => {
            normalize_tree(left, sigma);
            normalize_tree(right, sigma);
        }
        LdlTree::Leaf(vector) => {
            vector[0] = Complex::new(sigma / vector[0].re.sqrt(), 0.0);
            vector[1] = Complex64::zero();
        }
    }
}

/// Sample short polynomials using a Falcon tree. Algorithm 11 from the spec [1, p.40].
///
/// [1]: https://falcon-sign.info/falcon.pdf
pub(crate) fn ffsampling(
    t: &(Vec<Complex64>, Vec<Complex64>),
    tree: &LdlTree,
    parameters: &falcon::FalconParameters,
    rng: &mut dyn RngCore,
) -> (Vec<Complex64>, Vec<Complex64>) {
    match tree {
        LdlTree::Branch(ell, left, right) => {
            let bold_t1 = split_fft(&t.1);
            let bold_z1 = ffsampling(&bold_t1, right, parameters, rng);
            let z1 = merge_fft(&bold_z1.0, &bold_z1.1);

            // t0' = t0  + (t1 - z1) * l
            let t0_prime =
                t.0.iter()
                    .zip(t.1.iter().zip(z1.iter().zip(ell.iter())))
                    .map(|(t0_, (t1_, (z1_, l_)))| t0_ + (t1_ - z1_) * l_)
                    .collect_vec();

            let bold_t0 = split_fft(&t0_prime);
            let bold_z0 = ffsampling(&bold_t0, left, parameters, rng);
            let z0 = merge_fft(&bold_z0.0, &bold_z0.1);

            (z0, z1)
        }
        LdlTree::Leaf(value) => {
            let z0 = sampler_z(t.0[0].re, value[0].re, parameters.sigmin, rng);
            let z1 = sampler_z(t.1[0].re, value[0].re, parameters.sigmin, rng);
            (
                vec![Complex64::new(z0 as f64, 0.0)],
                vec![Complex64::new(z1 as f64, 0.0)],
            )
        }
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num_complex::{Complex, Complex64};
    use rand::{thread_rng, Rng};
    use rand_distr::num_traits::Zero;

    use crate::ffsampling::gram;

    #[test]
    fn test_gram() {
        let mut rng = thread_rng();
        let n = rng.gen_range(2..10);
        let a: [Vec<Complex64>; 4] = (0..4)
            .map(|_| {
                (0..n)
                    .map(|_| Complex::new(rng.gen(), rng.gen()))
                    .collect_vec()
            })
            .collect_vec()
            .try_into()
            .unwrap();
        let mut b = a.clone();
        b[0] = a[0].iter().map(|c| c.conj()).collect_vec();
        b[2] = a[1].iter().map(|c| c.conj()).collect_vec();
        b[1] = a[2].iter().map(|c| c.conj()).collect_vec();
        b[3] = a[3].iter().map(|c| c.conj()).collect_vec();

        let mut c: [Vec<Complex64>; 4] = (0..4)
            .map(|_| (0..n).map(|_| Complex64::zero()).collect_vec())
            .collect_vec()
            .try_into()
            .unwrap();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c[2 * i + j] = c[2 * i + j]
                        .iter()
                        .zip(
                            a[2 * i + k]
                                .iter()
                                .zip(b[2 * k + j].iter())
                                .map(|(aa, bb)| aa * bb),
                        )
                        .map(|(cc, ab)| cc + ab)
                        .collect_vec();
                }
            }
        }

        let g = gram(a);

        assert_eq!(c, g);
    }
}
