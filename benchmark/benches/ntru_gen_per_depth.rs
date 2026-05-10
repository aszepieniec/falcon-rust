use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use falcon_rust::{math::NTRU_SOLVE_BABAI_COEFF_BITS, polynomial::Polynomial};
use itertools::Itertools;
use num::{bigint::Sign, BigInt};
use rand::{rng, rngs::StdRng, Rng, SeedableRng};

fn random_bigint(bitsize: usize, rng: &mut StdRng) -> BigInt {
    let num_bytes = bitsize.div_ceil(8);
    let bytes = (0..num_bytes).map(|_| rng.random::<u8>()).collect_vec();
    let sign = if rng.random() { Sign::Minus } else { Sign::Plus };
    BigInt::from_bytes_be(sign, &bytes)
}

fn random_poly(n: usize, bits: usize, rng: &mut StdRng) -> Polynomial<BigInt> {
    Polynomial::new((0..n).map(|_| random_bigint(bits.max(1), rng)).collect_vec())
}

pub fn ntru_gen_per_depth(c: &mut Criterion) {
    let mut rng: StdRng = SeedableRng::from_seed(rng().random());
    let mut group = c.benchmark_group("ntru-gen per depth");

    // field_norm at each depth (going down the recursion).
    // Uses schoolbook multiplication on two half-sized polynomials (O((n/2)²)).
    // Input: size n = 1024 >> depth, coefficient bit size from the table.
    for depth in 0usize..=9 {
        let n = 1024usize >> depth;
        let fg_bits = NTRU_SOLVE_BABAI_COEFF_BITS[depth].0.ceil() as usize;

        group.bench_function(format!("field-norm/d={depth}"), |b| {
            b.iter_batched(
                || random_poly(n, fg_bits, &mut rng),
                |f| {
                    let _ = f.field_norm();
                },
                BatchSize::SmallInput,
            );
        });
    }

    // karatsuba + reduce_by_cyclotomic at each depth (coming back up the recursion).
    //
    // Depth 0 uses U32Field NTT in ntru_solve_entrypoint rather than BigInt
    // karatsuba, so it is omitted here.
    //
    // poly1 ≈ capital_f_prime_lifted: output of Babai at depth+1, whose
    //   post-reduction bit size ≈ fg_bits at depth+1 (column 0 of the table).
    // poly2 ≈ g_adjoint: same bit size as f,g at this depth.
    //
    // Two such karatsuba+reduce calls happen per ntru_solve invocation (one for
    // capital_F, one for capital_G); this benchmark times a single one.
    for depth in 1usize..=9 {
        let n = 1024usize >> depth;
        let bits1 = NTRU_SOLVE_BABAI_COEFF_BITS[depth + 1].0.ceil() as usize;
        let bits2 = NTRU_SOLVE_BABAI_COEFF_BITS[depth].0.ceil() as usize;

        group.bench_function(format!("karatsuba/d={depth}"), |b| {
            b.iter_batched(
                || (
                    random_poly(n, bits1, &mut rng),
                    random_poly(n, bits2, &mut rng),
                ),
                |(poly1, poly2)| {
                    let _ = poly1.karatsuba(&poly2).reduce_by_cyclotomic(n);
                },
                BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group!(benches, ntru_gen_per_depth);
criterion_main!(benches);
