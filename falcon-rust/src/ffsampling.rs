use falcon_profiler::profiling;
use num::{One, Zero};
use num_complex::Complex;
use rand::RngCore;

use crate::{
    falcon, fast_fft::FastFft, fixed_point::{FixedPoint64, FixedPoint128}, polynomial::Polynomial,
    samplerz::sampler_z,
};

type ComplexFP = Complex<FixedPoint64>;
type ComplexFixed128 = Complex<FixedPoint128>;

// ---------------------------------------------------------------------------
// FixedPoint128-based tree building (private)
// ---------------------------------------------------------------------------

fn gram_fixed128(b: [Polynomial<ComplexFixed128>; 4]) -> [Polynomial<ComplexFixed128>; 4] {
    const N: usize = 2;
    let mut g: [Polynomial<ComplexFixed128>; 4] = [
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
        Polynomial::zero(),
    ];
    for i in 0..N {
        for j in 0..N {
            for k in 0..N {
                g[N * i + j] = g[N * i + j].clone()
                    + b[N * i + k].hadamard_mul(&b[N * j + k].map(|c| c.conj()));
            }
        }
    }
    g
}

fn ldl_fixed128(
    g: [Polynomial<ComplexFixed128>; 4],
) -> (
    [Polynomial<ComplexFixed128>; 4],
    [Polynomial<ComplexFixed128>; 4],
) {
    let zero = Polynomial::<ComplexFixed128>::one();
    let one = Polynomial::<ComplexFixed128>::zero();

    let l10 = g[2].hadamard_div(&g[0]);
    let bc = l10.map(|c| c * c.conj());
    let abc = g[0].hadamard_mul(&bc);
    let d11 = g[3].clone() - abc;

    let l = [one.clone(), zero.clone(), l10.clone(), one];
    let d = [g[0].clone(), zero.clone(), zero, d11];
    (l, d)
}

enum LdlTreeFixed128 {
    Branch(
        Polynomial<ComplexFixed128>,
        Box<LdlTreeFixed128>,
        Box<LdlTreeFixed128>,
    ),
    Leaf([ComplexFixed128; 2]),
}

fn ffldl_fixed128(gram_matrix: [Polynomial<ComplexFixed128>; 4]) -> LdlTreeFixed128 {
    let n = gram_matrix[0].coefficients.len();
    let (l, d) = ldl_fixed128(gram_matrix);

    if n > 2 {
        let (d00, d01) = d[0].split_fft();
        let (d10, d11) = d[3].split_fft();
        let g0 = [d00.clone(), d01.clone(), d01.map(|c| c.conj()), d00];
        let g1 = [d10.clone(), d11.clone(), d11.map(|c| c.conj()), d10];
        LdlTreeFixed128::Branch(
            l[2].clone(),
            Box::new(ffldl_fixed128(g0)),
            Box::new(ffldl_fixed128(g1)),
        )
    } else {
        LdlTreeFixed128::Branch(
            l[2].clone(),
            Box::new(LdlTreeFixed128::Leaf(
                d[0].clone().coefficients.try_into().unwrap(),
            )),
            Box::new(LdlTreeFixed128::Leaf(
                d[3].clone().coefficients.try_into().unwrap(),
            )),
        )
    }
}

fn normalize_tree_fixed128(tree: &mut LdlTreeFixed128, sigma: FixedPoint128) {
    match tree {
        LdlTreeFixed128::Branch(_ell, left, right) => {
            normalize_tree_fixed128(left, sigma);
            normalize_tree_fixed128(right, sigma);
        }
        LdlTreeFixed128::Leaf(vector) => {
            vector[0] = Complex::new(sigma / vector[0].re.sqrt(), FixedPoint128::ZERO);
            vector[1] = ComplexFixed128::zero();
        }
    }
}

fn convert_tree_fixed128(tree: LdlTreeFixed128) -> LdlTree {
    match tree {
        LdlTreeFixed128::Branch(ell, left, right) => LdlTree::Branch(
            ell.map(|c| {
                Complex::new(
                    FixedPoint64::from(c.re),
                    FixedPoint64::from(c.im),
                )
            }),
            Box::new(convert_tree_fixed128(*left)),
            Box::new(convert_tree_fixed128(*right)),
        ),
        LdlTreeFixed128::Leaf(leaf) => LdlTree::Leaf([
            Complex::new(FixedPoint64::from(leaf[0].re), FixedPoint64::from(leaf[0].im)),
            Complex::new(FixedPoint64::from(leaf[1].re), FixedPoint64::from(leaf[1].im)),
        ]),
    }
}

/// Build the normalised LDL tree from b0 in FFT domain using FixedPoint128
/// arithmetic, then convert leaf values to FixedPoint64 for use by ffsampling.
#[profiling]
pub(crate) fn build_falcon_tree(
    b0_fft: [Polynomial<ComplexFixed128>; 4],
    sigma: FixedPoint128,
) -> LdlTree {
    let gram = gram_fixed128(b0_fft);
    let mut tree = ffldl_fixed128(gram);
    normalize_tree_fixed128(&mut tree, sigma);
    convert_tree_fixed128(tree)
}


// ---------------------------------------------------------------------------
// LDL tree (FixedPoint64) — stored in SecretKey, used by ffsampling
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub(crate) enum LdlTree {
    Branch(Polynomial<ComplexFP>, Box<LdlTree>, Box<LdlTree>),
    Leaf([ComplexFP; 2]),
}

/// Sample short polynomials using a Falcon tree. Algorithm 11 from the spec [1, p.40].
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[profiling]
pub(crate) fn ffsampling(
    t: &(Polynomial<ComplexFP>, Polynomial<ComplexFP>),
    tree: &LdlTree,
    parameters: &falcon::FalconParameters,
    rng: &mut dyn RngCore,
) -> (Polynomial<ComplexFP>, Polynomial<ComplexFP>) {
    match tree {
        LdlTree::Branch(ell, left, right) => {
            let bold_t1 = t.1.split_fft();
            let bold_z1 = ffsampling(&bold_t1, right, parameters, rng);
            let z1 = Polynomial::<ComplexFP>::merge_fft(&bold_z1.0, &bold_z1.1);

            // t0' = t0  + (t1 - z1) * l
            let t0_prime = t.0.clone() + (t.1.clone() - z1.clone()).hadamard_mul(ell);

            let bold_t0 = t0_prime.split_fft();
            let bold_z0 = ffsampling(&bold_t0, left, parameters, rng);
            let z0 = Polynomial::<ComplexFP>::merge_fft(&bold_z0.0, &bold_z0.1);

            (z0, z1)
        }
        LdlTree::Leaf(value) => {
            let z0 = sampler_z(
                t.0.coefficients[0].re,
                value[0].re,
                parameters.sigmin,
                rng,
            );
            let z1 = sampler_z(
                t.1.coefficients[0].re,
                value[0].re,
                parameters.sigmin,
                rng,
            );
            (
                Polynomial::new(vec![Complex::new(FixedPoint64::from(z0), FixedPoint64::ZERO)]),
                Polynomial::new(vec![Complex::new(FixedPoint64::from(z1), FixedPoint64::ZERO)]),
            )
        }
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num::Zero;
    use num_complex::Complex;
    use rand::{rng, Rng};

    use crate::{fixed_point::FixedPoint64, polynomial::Polynomial};

    type ComplexFP = Complex<FixedPoint64>;

    fn gram(b: [Polynomial<ComplexFP>; 4]) -> [Polynomial<ComplexFP>; 4] {
        const N: usize = 2;
        let mut g: [Polynomial<ComplexFP>; 4] = [
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
            Polynomial::zero(),
        ];
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    g[N * i + j] = g[N * i + j].clone()
                        + b[N * i + k].hadamard_mul(&b[N * j + k].map(|c| c.conj()));
                }
            }
        }
        g
    }

    #[test]
    fn test_gram() {
        let mut rng = rng();
        let n = rng.random_range(2..10);
        let a: [Polynomial<ComplexFP>; 4] = (0..4)
            .map(|_| {
                (0..n)
                    .map(|_| Complex::new(
                        FixedPoint64::from(rng.random::<f64>()),
                        FixedPoint64::from(rng.random::<f64>()),
                    ))
                    .collect_vec()
            })
            .map(Polynomial::new)
            .collect_vec()
            .try_into()
            .unwrap();
        let mut b = a.clone();
        b[0] = a[0].map(|c| c.conj());
        b[2] = a[1].map(|c| c.conj());
        b[1] = a[2].map(|c| c.conj());
        b[3] = a[3].map(|c| c.conj());

        let mut c: [Polynomial<ComplexFP>; 4] = (0..4)
            .map(|_| (0..n).map(|_| ComplexFP::zero()).collect_vec())
            .map(Polynomial::new)
            .collect_vec()
            .try_into()
            .unwrap();
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    c[2 * i + j] = c[2 * i + j].clone() + a[2 * i + k].hadamard_mul(&b[2 * k + j]);
                }
            }
        }

        let g = gram(a);

        assert_eq!(c, g);
    }
}
