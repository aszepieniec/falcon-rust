#![allow(unused)]

use std::{
    fs::File,
    io::{Read, Write},
};

use benchmark::{load_helper_data, save_helper_data};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use falcon_rust::falcon::{self, PublicKey, SecretKey, Signature, FALCON_1024, FALCON_512};
use rand::{thread_rng, Rng};

const NUM_KEYS: usize = 10;
const SIGS_PER_KEY: usize = 10;

pub fn falcon_rust_operation(c: &mut Criterion) {
    let mut rng = thread_rng();
    let mut keys512 = vec![];
    let mut keys1024 = vec![];

    let mut group = c.benchmark_group("keygen");
    group.sample_size(NUM_KEYS);
    group.bench_function("keygen 512", |b| {
        b.iter(|| {
            keys512.push(FALCON_512.keygen(rng.gen()));
        })
    });
    group.bench_function("keygen 1024", |b| {
        b.iter(|| {
            keys1024.push(FALCON_1024.keygen(rng.gen()));
        })
    });
    group.finish();

    let mut group = c.benchmark_group("sign");
    let mut sigs512 = vec![];
    let mut msgs512 = vec![];
    let mut sigs1024 = vec![];
    let mut msgs1024 = vec![];
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut i = 0;
    group.bench_function("sign 512", |b| {
        b.iter(|| {
            let msg: [u8; 32] = rng.gen();
            let sig = FALCON_512.sign(&msg, &keys512[i % NUM_KEYS].0);
            msgs512.push(msg);
            sigs512.push(sig);
            i += 1;
        })
    });
    i = 0;
    group.bench_function("sign 1024", |b| {
        b.iter(|| {
            let msg: [u8; 32] = rng.gen();
            let sig = FALCON_1024.sign(&msg, &keys1024[i % NUM_KEYS].0);
            msgs1024.push(msg);
            sigs1024.push(sig);
            i += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("verify");
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    i = 0;
    group.bench_function("verify 512", |b| {
        b.iter(|| {
            let sig = FALCON_512.verify(
                &msgs512[i % msgs512.len()],
                &sigs512[i % sigs512.len()],
                &keys512[i % NUM_KEYS].1,
            );
            i += 1;
        })
    });
    i = 0;
    group.bench_function("verify 1024", |b| {
        b.iter(|| {
            let sig = FALCON_512.verify(
                &msgs1024[i % msgs1024.len()],
                &sigs1024[i % sigs1024.len()],
                &keys1024[i % NUM_KEYS].1,
            );
            i += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("serialize");
    let mut serialized_secret_keys_512 = vec![vec![]; (NUM_KEYS)];
    let mut serialized_public_keys_512 = vec![vec![]; (NUM_KEYS)];
    let mut serialized_signatures_512 = vec![vec![]; (NUM_KEYS * SIGS_PER_KEY)];
    let mut serialized_secret_keys_1024 = vec![vec![]; (NUM_KEYS)];
    let mut serialized_public_keys_1024 = vec![vec![]; (NUM_KEYS)];
    let mut serialized_signatures_1024 = vec![vec![]; (NUM_KEYS * SIGS_PER_KEY)];
    i = 0;
    group.bench_function("serialize secret key 512", |b| {
        b.iter(|| {
            serialized_secret_keys_512[i % (NUM_KEYS)] =
                keys512[i % keys512.len()].0.clone().to_bytes();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("serialize secret key 1024", |b| {
        b.iter(|| {
            serialized_secret_keys_1024[i % (NUM_KEYS)] =
                keys1024[i % keys1024.len()].0.clone().to_bytes();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("serialize public key 512", |b| {
        b.iter(|| {
            serialized_public_keys_512[i % (NUM_KEYS)] =
                keys512[i % keys512.len()].1.clone().to_bytes();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("serialize public key 1024", |b| {
        b.iter(|| {
            serialized_public_keys_1024[i % (NUM_KEYS)] =
                keys1024[i % keys1024.len()].1.clone().to_bytes();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("serialize signature 512", |b| {
        b.iter(|| {
            serialized_signatures_512[i % (NUM_KEYS * SIGS_PER_KEY)] =
                sigs512[i % sigs512.len()].clone().to_bytes();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("serialize signature 1024", |b| {
        b.iter(|| {
            serialized_signatures_1024[i % (NUM_KEYS * SIGS_PER_KEY)] =
                sigs1024[i % sigs1024.len()].clone().to_bytes();
            i += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("deserialize");
    i = 0;
    group.bench_function("deserialize secret key 512", |b| {
        b.iter(|| {
            SecretKey::from_bytes(&serialized_secret_keys_512[i % (NUM_KEYS)]).unwrap();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("deserialize secret key 1024", |b| {
        b.iter(|| {
            SecretKey::from_bytes(&serialized_secret_keys_1024[i % (NUM_KEYS)]).unwrap();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("deserialize public key 512", |b| {
        b.iter(|| {
            PublicKey::from_bytes(&serialized_public_keys_512[i % (NUM_KEYS)]).unwrap();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("deserialize public key 1024", |b| {
        b.iter(|| {
            PublicKey::from_bytes(&serialized_public_keys_1024[i % (NUM_KEYS)]).unwrap();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("deserialize signature 512", |b| {
        b.iter(|| {
            Signature::from_bytes(&serialized_signatures_512[i % (NUM_KEYS * SIGS_PER_KEY)])
                .unwrap();
            i += 1;
        })
    });
    i = 0;
    group.bench_function("deserialize signature 1024", |b| {
        b.iter(|| {
            Signature::from_bytes(&serialized_signatures_1024[i % (NUM_KEYS * SIGS_PER_KEY)])
                .unwrap();
            i += 1;
        })
    });
    group.finish();
}

criterion_group!(benches, falcon_rust_operation);
criterion_main!(benches);
