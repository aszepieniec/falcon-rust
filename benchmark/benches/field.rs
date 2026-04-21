use criterion::{criterion_group, criterion_main, Criterion};
use falcon_rust::falcon_field::Felt;
use rand::rng;
use rand::Rng;

pub fn field(c: &mut Criterion) {
    let mut rng = rng();
    let mut group = c.benchmark_group("field");
    group.sample_size(10);

    const NUM_OPS: usize = 1000_usize;
    let mut a = [Felt::new(0); NUM_OPS];
    let mut b = [Felt::new(0); NUM_OPS];
    let mut c = [Felt::new(0); NUM_OPS];
    for i in 0..NUM_OPS {
        a[i] = rng.random();
        b[i] = rng.random();
    }

    group.bench_function(format!("{NUM_OPS} additions"), |bencher| {
        bencher.iter(|| {
            for i in 0..NUM_OPS {
                c[i] = a[i] + b[i];
            }
        })
    });

    group.bench_function(format!("{NUM_OPS} multiplications"), |bencher| {
        bencher.iter(|| {
            for i in 0..NUM_OPS {
                c[i] = a[i] * b[i];
            }
        })
    });
    group.finish();
}

criterion_group!(benches, field);
criterion_main!(benches);
