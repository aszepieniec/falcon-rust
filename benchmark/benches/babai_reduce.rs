use criterion::{criterion_group, criterion_main, Criterion};
use falcon_rust::polynomial::Polynomial;
use itertools::Itertools;
use num::{bigint::Sign, BigInt};
use rand::{rngs::StdRng, thread_rng, Rng, SeedableRng};

const BITSIZE_SMALL: [usize; 10] = [4, 11, 24, 51, 100, 199, 399, 792, 1580, 3151];
const BITSIZE_LARGE: [usize; 10] = [
    19, 39, 12501, 12480, 12428, 12328, 12139, 11745, 10975, 9417,
];
const N: [usize; 10] = [1024, 512, 256, 128, 64, 32, 16, 8, 4, 2];

fn random_integer(bitsize: usize, rng: &mut StdRng) -> BigInt {
    let num_bytes = (bitsize + 7) / 8;
    let bytes = (0..num_bytes).map(|_| rng.gen::<u8>()).collect_vec();
    let sign = match rng.gen() {
        true => Sign::Minus,
        false => Sign::Plus,
    };
    let int = BigInt::from_bytes_be(sign, &bytes);
    int
}

fn generate_reducible_polynomials(
    bitsize_lowercase: usize,
    bitsize_capital: usize,
    n: usize,
    rng: &mut StdRng,
) -> (
    Polynomial<BigInt>,
    Polynomial<BigInt>,
    Polynomial<BigInt>,
    Polynomial<BigInt>,
) {
    let f = Polynomial::new(
        (0..n)
            .map(|_| random_integer(bitsize_lowercase, rng))
            .collect_vec(),
    );
    let g = Polynomial::new(
        (0..n)
            .map(|_| random_integer(bitsize_lowercase, rng))
            .collect_vec(),
    );
    let capital_f = Polynomial::new(
        (0..n)
            .map(|_| random_integer(bitsize_capital, rng))
            .collect_vec(),
    );
    let capital_g = Polynomial::new(
        (0..n)
            .map(|_| random_integer(bitsize_capital, rng))
            .collect_vec(),
    );
    (f, g, capital_f, capital_g)
}

pub fn babai_reduce(c: &mut Criterion) {
    let mut rng: StdRng = SeedableRng::from_seed(thread_rng().gen());
    let mut group = c.benchmark_group("babai reduce bigint");

    for (&n, (&small, large)) in N.iter().zip(BITSIZE_SMALL.iter().zip(BITSIZE_LARGE)) {
        group.bench_function(format!("{small}/{large}/{n}"), |b| {
            b.iter(|| {
                let (f, g, mut capital_f, mut capital_g) =
                    generate_reducible_polynomials(small, large, n, &mut rng);
                let _ =
                    falcon_rust::math::babai_reduce_bigint(&f, &g, &mut capital_f, &mut capital_g);
            })
        });
    }

    group.finish();

    let mut group = c.benchmark_group("babai reduce mmi");

    for (&n, (&small, large)) in N.iter().zip(BITSIZE_SMALL.iter().zip(BITSIZE_LARGE)) {
        group.bench_function(format!("{small}/{large}/{n}"), |b| {
            b.iter(|| {
                let (f, g, mut capital_f, mut capital_g) =
                    generate_reducible_polynomials(small, large, n, &mut rng);
                let _ = falcon_rust::math::babai_reduce_mmi(&f, &g, &mut capital_f, &mut capital_g);
            })
        });
    }

    group.finish();
}

criterion_group!(benches, babai_reduce);
criterion_main!(benches);
