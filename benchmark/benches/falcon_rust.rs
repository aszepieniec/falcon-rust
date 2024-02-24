#![allow(unused)]

use std::{
    fs::File,
    io::{Read, Write},
};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use pqcrypto_falcon::*;
use rand::{thread_rng, Rng};

const NUM_KEYS: usize = 10;
const SIGS_PER_KEY: usize = 10;

pub fn falcon_rust_operation(c: &mut Criterion) {
    let mut rng = thread_rng();
    let mut keys512 = (0..NUM_KEYS)
        .map(|_| falcon_rust::falcon512::keygen(rng.gen()))
        .collect_vec();
    let mut keys1024 = (0..NUM_KEYS)
        .map(|_| falcon_rust::falcon1024::keygen(rng.gen()))
        .collect_vec();
    let mut msgs512 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.gen::<[u8; 15]>())
        .collect_vec();
    let mut msgs1024 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.gen::<[u8; 15]>())
        .collect_vec();
    let mut sigs512 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|i| falcon_rust::falcon512::sign(&msgs512[i], &keys512[i % NUM_KEYS].0))
        .collect_vec();
    let mut sigs1024 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|i| falcon_rust::falcon1024::sign(&msgs1024[i], &keys1024[i % NUM_KEYS].0))
        .collect_vec();
    let mut serialized_secret_keys_512 = keys512.iter().map(|k| k.0.to_bytes()).collect_vec();
    let mut serialized_public_keys_512 = keys512.iter().map(|k| k.1.to_bytes()).collect_vec();
    let mut serialized_signatures_512 = sigs512.iter().map(|s| s.to_bytes()).collect_vec();
    let mut serialized_secret_keys_1024 = keys1024.iter().map(|k| k.0.to_bytes()).collect_vec();
    let mut serialized_public_keys_1024 = keys1024.iter().map(|k| k.1.to_bytes()).collect_vec();
    let mut serialized_signatures_1024 = sigs1024.iter().map(|s| s.to_bytes()).collect_vec();

    let mut group = c.benchmark_group("falcon-rust");
    group.sample_size(NUM_KEYS);
    group.bench_function("keygen 512", |b| {
        b.iter(|| {
            falcon_rust::falcon512::keygen(rng.gen());
        })
    });
    group.bench_function("keygen 1024", |b| {
        b.iter(|| {
            falcon_rust::falcon1024::keygen(rng.gen());
        })
    });
    group.finish();

    let mut group = c.benchmark_group("falcon-rust");
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_sign_512 = 0;
    group.bench_function("sign 512", |b| {
        b.iter(|| {
            falcon_rust::falcon512::sign(
                &msgs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)],
                &keys512[iterator_sign_512 % NUM_KEYS].0,
            );
            iterator_sign_512 += 1;
        })
    });
    let mut iterator_sign_1024 = 0;
    group.bench_function("sign 1024", |b| {
        b.iter(|| {
            falcon_rust::falcon1024::sign(
                &msgs1024[iterator_sign_1024 % (NUM_KEYS * SIGS_PER_KEY)],
                &keys1024[iterator_sign_1024 % NUM_KEYS].0,
            );
            iterator_sign_1024 += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("falcon-rust");
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_verify_512 = 0;
    group.bench_function("verify 512", |b| {
        b.iter(|| {
            assert!(falcon_rust::falcon512::verify(
                &msgs512[iterator_verify_512 % msgs512.len()],
                &sigs512[iterator_verify_512 % sigs512.len()],
                &keys512[iterator_verify_512 % NUM_KEYS].1,
            ));
            iterator_verify_512 += 1;
        })
    });
    let mut iterator_verify_1024 = 0;
    group.bench_function("verify 1024", |b| {
        b.iter(|| {
            assert!(falcon_rust::falcon1024::verify(
                &msgs1024[iterator_verify_1024 % msgs1024.len()],
                &sigs1024[iterator_verify_1024 % sigs1024.len()],
                &keys1024[iterator_verify_1024 % NUM_KEYS].1,
            ));
            iterator_verify_1024 += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("falcon-rust");
    let mut iterator_serialize_sk_512 = 0;
    group.bench_function("serialize sk 512", |b| {
        b.iter(|| {
            keys512[iterator_serialize_sk_512 % keys512.len()]
                .0
                .clone()
                .to_bytes();
            iterator_serialize_sk_512 += 1;
        })
    });
    let mut iterator_serialize_sk_1024 = 0;
    group.bench_function("serialize sk 1024", |b| {
        b.iter(|| {
            keys1024[iterator_serialize_sk_1024 % keys1024.len()]
                .0
                .clone()
                .to_bytes();
            iterator_serialize_sk_1024 += 1;
        })
    });
    let mut iterator_serialize_pk_512 = 0;
    group.bench_function("serialize pk 512", |b| {
        b.iter(|| {
            keys512[iterator_serialize_pk_512 % keys512.len()]
                .1
                .clone()
                .to_bytes();
            iterator_serialize_pk_512 += 1;
        })
    });
    let mut iterator_serialize_pk_1024 = 0;
    group.bench_function("serialize pk 1024", |b| {
        b.iter(|| {
            keys1024[iterator_serialize_pk_1024 % keys1024.len()]
                .1
                .clone()
                .to_bytes();
            iterator_serialize_pk_1024 += 1;
        })
    });
    let mut iterator_serialize_sig_512 = 0;
    group.bench_function("serialize sig 512", |b| {
        b.iter(|| {
            sigs512[iterator_serialize_sig_512 % sigs512.len()]
                .clone()
                .to_bytes();
            iterator_serialize_sig_512 += 1;
        })
    });
    let mut iterator_serialize_sig_1024 = 0;
    group.bench_function("serialize sig 1024", |b| {
        b.iter(|| {
            sigs1024[iterator_serialize_sig_1024 % sigs1024.len()]
                .clone()
                .to_bytes();
            iterator_serialize_sig_1024 += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("falcon-rust");
    let mut iterator_deserialize_sk_512 = 0;
    group.bench_function("deserialize sk 512", |b| {
        b.iter(|| {
            falcon_rust::falcon512::SecretKey::from_bytes(
                &serialized_secret_keys_512
                    [iterator_deserialize_sk_512 % serialized_secret_keys_512.len()],
            )
            .unwrap();
            iterator_deserialize_sk_512 += 1;
        })
    });
    let mut iterator_deserialize_sk_1024 = 0;
    group.bench_function("deserialize sk 1024", |b| {
        b.iter(|| {
            falcon_rust::falcon1024::SecretKey::from_bytes(
                &serialized_secret_keys_1024
                    [iterator_deserialize_sk_1024 % serialized_secret_keys_1024.len()],
            )
            .unwrap();
            iterator_deserialize_sk_1024 += 1;
        })
    });
    let mut iterator_deserialize_pk_512 = 0;
    group.bench_function("deserialize pk 512", |b| {
        b.iter(|| {
            falcon_rust::falcon512::PublicKey::from_bytes(
                &serialized_public_keys_512
                    [iterator_deserialize_pk_512 % serialized_public_keys_512.len()],
            )
            .unwrap();
            iterator_deserialize_pk_512 += 1;
        })
    });
    let mut iterator_deserialize_pk_1024 = 0;
    group.bench_function("deserialize pk 1024", |b| {
        b.iter(|| {
            falcon_rust::falcon1024::PublicKey::from_bytes(
                &serialized_public_keys_1024
                    [iterator_deserialize_pk_1024 % serialized_public_keys_1024.len()],
            )
            .unwrap();
            iterator_deserialize_pk_1024 += 1;
        })
    });
    let mut iterator_deserialize_sig_512 = 0;
    group.bench_function("deserialize sig 512", |b| {
        b.iter(|| {
            falcon_rust::falcon512::Signature::from_bytes(
                &serialized_signatures_512
                    [iterator_deserialize_sig_512 % serialized_signatures_512.len()],
            )
            .unwrap();
            iterator_deserialize_sig_512 += 1;
        })
    });
    let mut iterator_deserialize_sig_1024 = 0;
    group.bench_function("deserialize sig 1024", |b| {
        b.iter(|| {
            falcon_rust::falcon1024::Signature::from_bytes(
                &serialized_signatures_1024
                    [iterator_deserialize_sig_1024 % serialized_signatures_1024.len()],
            )
            .unwrap();
            iterator_deserialize_sig_1024 += 1;
        })
    });
    group.finish();
}

fn falcon_c_ffi_operation(c: &mut Criterion) {
    let mut rng = thread_rng();
    let mut keys512 = (0..NUM_KEYS).map(|_| falcon512::keypair()).collect_vec();
    let mut keys1024 = (0..NUM_KEYS).map(|_| falcon1024::keypair()).collect_vec();

    let mut group = c.benchmark_group("c ffi");
    group.sample_size(NUM_KEYS);
    group.bench_function("keygen 512", |b| {
        b.iter(|| {
            falcon512::keypair();
        })
    });
    group.bench_function("keygen 1024", |b| {
        b.iter(|| {
            falcon1024::keypair();
        })
    });
    group.finish();

    let mut group = c.benchmark_group("c ffi");
    let mut msgs512 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.gen::<[u8; 15]>())
        .collect_vec();
    let mut sigs512 = msgs512
        .iter()
        .enumerate()
        .map(|(i, msg)| falcon512::detached_sign(msg, &keys512[i % NUM_KEYS].1))
        .collect_vec();
    let mut msgs1024 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.gen::<[u8; 15]>())
        .collect_vec();
    let mut sigs1024 = msgs1024
        .iter()
        .enumerate()
        .map(|(i, msg)| falcon1024::detached_sign(msg, &keys1024[i % NUM_KEYS].1))
        .collect_vec();
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_sign_512 = 0;
    group.bench_function("sign 512", |b| {
        b.iter(|| {
            sigs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)] = falcon512::detached_sign(
                &msgs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)],
                &keys512[iterator_sign_512 % NUM_KEYS].1,
            );
            iterator_sign_512 += 1;
        })
    });
    let mut iterator_sign_1024 = 0;
    group.bench_function("sign 1024", |b| {
        b.iter(|| {
            sigs1024[iterator_sign_1024 % (NUM_KEYS * SIGS_PER_KEY)] = falcon1024::detached_sign(
                &msgs1024[iterator_sign_1024 % (NUM_KEYS * SIGS_PER_KEY)],
                &keys1024[iterator_sign_1024 % NUM_KEYS].1,
            );
            iterator_sign_1024 += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("c ffi");
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_verify_512 = 0;
    group.bench_function("verify 512", |b| {
        b.iter(|| {
            assert!(falcon512::verify_detached_signature(
                &sigs512[iterator_verify_512 % sigs512.len()],
                &msgs512[iterator_verify_512 % msgs512.len()],
                &keys512[iterator_verify_512 % NUM_KEYS].0,
            )
            .is_ok());
            iterator_verify_512 += 1;
        })
    });
    let mut iterator_verify_1024 = 0;
    group.bench_function("verify 1024", |b| {
        b.iter(|| {
            assert!(falcon1024::verify_detached_signature(
                &sigs1024[iterator_verify_1024 % sigs1024.len()],
                &msgs1024[iterator_verify_1024 % msgs1024.len()],
                &keys1024[iterator_verify_1024 % NUM_KEYS].0,
            )
            .is_ok());
            iterator_verify_1024 += 1;
        })
    });
    group.finish();
}

criterion_group!(benches, falcon_rust_operation, falcon_c_ffi_operation);
criterion_main!(benches);
