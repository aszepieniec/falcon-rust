use criterion::{criterion_group, criterion_main, Criterion};
use falcon_rust::{multimod::MultiModInt, polynomial::Polynomial};
use itertools::Itertools;
use num::{bigint::Sign, BigInt};
use rand::{rngs::StdRng, thread_rng, Rng, SeedableRng};

const BITSIZE_SMALL: [usize; 10] = [4, 11, 24, 51, 102, 201, 399, 792, 1580, 3151];
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

fn generate_polynomials(
    bitsize_lowercase: usize,
    bitsize_capital: usize,
    n: usize,
    rng: &mut StdRng,
) -> (Polynomial<BigInt>, Polynomial<BigInt>) {
    let lowercase = Polynomial::new(
        (0..n)
            .map(|_| random_integer(bitsize_lowercase, rng))
            .collect_vec(),
    );
    let uppercase = Polynomial::new(
        (0..n)
            .map(|_| random_integer(bitsize_capital, rng))
            .collect_vec(),
    );
    (uppercase, lowercase)
}

fn naive_cyclomul(
    uppercase: &Polynomial<BigInt>,
    lowercase: &Polynomial<BigInt>,
    n: usize,
) -> Polynomial<BigInt> {
    (uppercase * lowercase).reduce_by_cyclotomic(n)
}

fn karatsuba_cyclomul(
    uppercase: &Polynomial<BigInt>,
    lowercase: &Polynomial<BigInt>,
    n: usize,
) -> Polynomial<BigInt> {
    uppercase.karatsuba(lowercase).reduce_by_cyclotomic(n)
}

fn multimodint_cyclomul(
    uppercase: &Polynomial<BigInt>,
    lowercase: &Polynomial<BigInt>,
    capacity: usize,
) -> Polynomial<BigInt> {
    let u = uppercase.map(|c| MultiModInt::from_bigint_with_capacity(c, capacity));
    let l = lowercase.map(|c| MultiModInt::from_bigint_with_capacity(c, capacity));
    l.cyclotomic_mul(&u).map(|m| BigInt::from(m.clone()))
}

pub fn cyclomul(c: &mut Criterion) {
    let mut rng: StdRng = SeedableRng::from_seed(thread_rng().gen());
    let mut group = c.benchmark_group("cyclomul naive");

    for (&n, (&small, large)) in N.iter().zip(BITSIZE_SMALL.iter().zip(BITSIZE_LARGE)) {
        group.bench_function(format!("{small}/{large}/{n}"), |b| {
            b.iter(|| {
                let (uppercase, lowercase) = generate_polynomials(small, large, n, &mut rng);
                let _ = naive_cyclomul(&uppercase, &lowercase, n);
            })
        });
    }

    group.finish();

    let mut group = c.benchmark_group("cyclomul karatsuba");

    for (&n, (&small, large)) in N.iter().zip(BITSIZE_SMALL.iter().zip(BITSIZE_LARGE)) {
        group.bench_function(format!("{small}/{large}/{n}"), |b| {
            b.iter(|| {
                let (uppercase, lowercase) = generate_polynomials(small, large, n, &mut rng);
                let _ = karatsuba_cyclomul(&uppercase, &lowercase, n);
            })
        });
    }

    group.finish();

    let mut group = c.benchmark_group("cyclomul multimod");

    for (&n, (&small, large)) in N.iter().zip(BITSIZE_SMALL.iter().zip(BITSIZE_LARGE)) {
        group.bench_function(format!("{small}/{large}/{n}"), |b| {
            b.iter(|| {
                let (uppercase, lowercase) = generate_polynomials(small, large, n, &mut rng);
                let _ = multimodint_cyclomul(&uppercase, &lowercase, large);
            })
        });
    }

    group.finish();
}

criterion_group!(benches, cyclomul);
criterion_main!(benches);
