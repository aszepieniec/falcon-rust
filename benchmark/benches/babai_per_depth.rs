use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use falcon_rust::{math::NTRU_SOLVE_BABAI_COEFF_BITS, polynomial::Polynomial};
use itertools::Itertools;
use num::{bigint::Sign, BigInt};
use rand::{rng, rngs::StdRng, Rng, SeedableRng};

fn random_integer(bitsize: usize, rng: &mut StdRng) -> BigInt {
    let num_bytes = bitsize.div_ceil(8);
    let bytes = (0..num_bytes).map(|_| rng.random::<u8>()).collect_vec();
    let sign = if rng.random() { Sign::Minus } else { Sign::Plus };
    BigInt::from_bytes_be(sign, &bytes)
}

fn random_i32(bitsize: usize, rng: &mut StdRng) -> i32 {
    let bits = bitsize.min(31);
    let mag: i64 = rng.random_range(0..(1i64 << bits));
    if rng.random() { mag as i32 } else { -(mag as i32) }
}

/// Returns `(f_i32, g_i32, capital_f_bi, capital_g_bi)` at the given depth.
fn generate_depth_inputs_bigint(
    depth: usize,
    rng: &mut StdRng,
) -> (Polynomial<BigInt>, Polynomial<BigInt>, Polynomial<BigInt>, Polynomial<BigInt>) {
    let n = 1024usize >> depth;
    let (fg_bits_f, _, cap_bits_f, _) = NTRU_SOLVE_BABAI_COEFF_BITS[depth];
    let fg_bits = fg_bits_f.ceil() as usize;
    let cap_bits = cap_bits_f.ceil() as usize;
    let f = Polynomial::new((0..n).map(|_| random_integer(fg_bits, rng)).collect_vec());
    let g = Polynomial::new((0..n).map(|_| random_integer(fg_bits, rng)).collect_vec());
    let cf = Polynomial::new((0..n).map(|_| random_integer(cap_bits, rng)).collect_vec());
    let cg = Polynomial::new((0..n).map(|_| random_integer(cap_bits, rng)).collect_vec());
    (f, g, cf, cg)
}

/// Returns `(f_i32, g_i32, capital_f_i128, capital_g_i128)` for RNS depths.
/// f,g use i32 (they fit: depth 1 ≈11 bits, depth 2 ≈24 bits).
fn generate_depth_inputs_rns(
    depth: usize,
    rng: &mut StdRng,
) -> (Polynomial<i32>, Polynomial<i32>, Vec<i128>, Vec<i128>) {
    let n = 1024usize >> depth;
    let (fg_bits_f, _, cap_bits_f, _) = NTRU_SOLVE_BABAI_COEFF_BITS[depth];
    let fg_bits = fg_bits_f.ceil() as usize;
    let cap_bits = cap_bits_f.ceil() as usize;
    let f = Polynomial::new((0..n).map(|_| random_i32(fg_bits, rng)).collect_vec());
    let g = Polynomial::new((0..n).map(|_| random_i32(fg_bits, rng)).collect_vec());
    let cf: Vec<i128> = (0..n).map(|_| {
        let bits = cap_bits.min(127);
        let mag: u128 = rng.random_range(0..(1u128 << bits));
        if rng.random() { mag as i128 } else { -(mag as i128) }
    }).collect();
    let cg: Vec<i128> = (0..n).map(|_| {
        let bits = cap_bits.min(127);
        let mag: u128 = rng.random_range(0..(1u128 << bits));
        if rng.random() { mag as i128 } else { -(mag as i128) }
    }).collect();
    (f, g, cf, cg)
}

pub fn babai_per_depth(c: &mut Criterion) {
    let mut rng: StdRng = SeedableRng::from_seed(rng().random());
    let mut group = c.benchmark_group("babai per depth");

    // BigInt variant at every depth 1–9.
    for depth in 1usize..=9 {
        group.bench_function(format!("bigint/d={depth}"), |b| {
            b.iter_batched(
                || generate_depth_inputs_bigint(depth, &mut rng),
                |(f, g, mut cf, mut cg)| {
                    let _ = falcon_rust::math::babai_reduce_bigint(&f, &g, &mut cf, &mut cg);
                },
                BatchSize::SmallInput,
            );
        });
    }

    // RNS variant at depth 1 (K=2, n=512).
    group.bench_function("rns-k2/d=1", |b| {
        b.iter_batched(
            || generate_depth_inputs_rns(1, &mut rng),
            |(f, g, mut cf, mut cg)| {
                let _ = falcon_rust::math::babai_reduce_rns_depth1(&f, &g, &mut cf, &mut cg);
            },
            BatchSize::SmallInput,
        );
    });

    // RNS variant at depth 2 (K=4, n=256).
    group.bench_function("rns-k4/d=2", |b| {
        b.iter_batched(
            || generate_depth_inputs_rns(2, &mut rng),
            |(f, g, mut cf, mut cg)| {
                let _ = falcon_rust::math::babai_reduce_rns_depth2(&f, &g, &mut cf, &mut cg);
            },
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, babai_per_depth);
criterion_main!(benches);
