use criterion::{criterion_group, criterion_main, Criterion};
use rand::{thread_rng, Rng};

pub fn ntru_gen(c: &mut Criterion) {
    let mut rng = thread_rng();
    let mut group = c.benchmark_group("ntru-gen");
    group.sample_size(10);

    group.bench_function("ntru-gen-512", |b| {
        b.iter(|| {
            let _ = falcon_rust::math::ntru_gen(512, rng.gen());
        })
    });

    group.bench_function("ntru-gen-1024", |b| {
        b.iter(|| {
            let _ = falcon_rust::math::ntru_gen(1024, rng.gen());
        })
    });
    group.finish();
}

criterion_group!(benches, ntru_gen);
criterion_main!(benches);
