use criterion::{criterion_group, criterion_main, Criterion};
use rand::rng;

pub fn ntru_gen(c: &mut Criterion) {
    let mut rng = rng();
    let mut group = c.benchmark_group("ntru-gen");
    group.sample_size(10);

    group.bench_function("ntru-gen-512", |b| {
        b.iter(|| {
            let _ = falcon_rust::math::ntru_gen(512, &mut rng);
        })
    });

    group.bench_function("ntru-gen-1024", |b| {
        b.iter(|| {
            let _ = falcon_rust::math::ntru_gen(1024, &mut rng);
        })
    });
    group.finish();
}

criterion_group!(benches, ntru_gen);
criterion_main!(benches);
