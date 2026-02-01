#![allow(unused)]

use std::fs::File;
use std::io::{Read, Write};

use criterion::{criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use rand::rngs::{OsRng, StdRng};
use rand::{rng, Rng, SeedableRng};

const NUM_KEYS: usize = 10;
const SIGS_PER_KEY: usize = 10;

pub fn falcon_rust_operation(c: &mut Criterion) {
    let mut rng = rng();
    let mut keys512 = (0..NUM_KEYS)
        .map(|_| falcon_rust::falcon512::keygen(rng.random()))
        .collect_vec();
    let mut keys1024 = (0..NUM_KEYS)
        .map(|_| falcon_rust::falcon1024::keygen(rng.random()))
        .collect_vec();
    let mut msgs512 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.random::<[u8; 15]>())
        .collect_vec();
    let mut msgs1024 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.random::<[u8; 15]>())
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
            falcon_rust::falcon512::keygen(rng.random());
        })
    });
    group.bench_function("keygen 1024", |b| {
        b.iter(|| {
            falcon_rust::falcon1024::keygen(rng.random());
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

fn falcon_cffi_operation(c: &mut Criterion) {
    let mut rng = rng();
    let mut keys512 = (0..NUM_KEYS)
        .map(|_| pqcrypto_falcon::falcon512::keypair())
        .collect_vec();
    let mut keys1024 = (0..NUM_KEYS)
        .map(|_| pqcrypto_falcon::falcon1024::keypair())
        .collect_vec();

    let mut group = c.benchmark_group("c ffi");
    group.sample_size(NUM_KEYS);
    group.bench_function("keygen 512", |b| {
        b.iter(|| {
            pqcrypto_falcon::falcon512::keypair();
        })
    });
    group.bench_function("keygen 1024", |b| {
        b.iter(|| {
            pqcrypto_falcon::falcon1024::keypair();
        })
    });
    group.finish();

    let mut group = c.benchmark_group("c ffi");
    let mut msgs512 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.random::<[u8; 15]>())
        .collect_vec();
    let mut sigs512 = msgs512
        .iter()
        .enumerate()
        .map(|(i, msg)| pqcrypto_falcon::falcon512::detached_sign(msg, &keys512[i % NUM_KEYS].1))
        .collect_vec();
    let mut msgs1024 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.random::<[u8; 15]>())
        .collect_vec();
    let mut sigs1024 = msgs1024
        .iter()
        .enumerate()
        .map(|(i, msg)| pqcrypto_falcon::falcon1024::detached_sign(msg, &keys1024[i % NUM_KEYS].1))
        .collect_vec();
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_sign_512 = 0;
    group.bench_function("sign 512", |b| {
        b.iter(|| {
            sigs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)] =
                pqcrypto_falcon::falcon512::detached_sign(
                    &msgs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)],
                    &keys512[iterator_sign_512 % NUM_KEYS].1,
                );
            iterator_sign_512 += 1;
        })
    });
    let mut iterator_sign_1024 = 0;
    group.bench_function("sign 1024", |b| {
        b.iter(|| {
            sigs1024[iterator_sign_1024 % (NUM_KEYS * SIGS_PER_KEY)] =
                pqcrypto_falcon::falcon1024::detached_sign(
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
            assert!(pqcrypto_falcon::falcon512::verify_detached_signature(
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
            assert!(pqcrypto_falcon::falcon1024::verify_detached_signature(
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

fn falcon_fndsa_operation(c: &mut Criterion) {
    let mut kg = fn_dsa::KeyPairGeneratorStandard::default();
    let mut rng = StdRng::from_os_rng();

    let mut keys512 = (0..NUM_KEYS)
        .map(|_| {
            let mut sign_key = [0u8; fn_dsa::sign_key_size(fn_dsa::FN_DSA_LOGN_512)];
            let mut vrfy_key = [0u8; fn_dsa::vrfy_key_size(fn_dsa::FN_DSA_LOGN_512)];
            fn_dsa::KeyPairGenerator::keygen(
                &mut kg,
                fn_dsa::FN_DSA_LOGN_512,
                &mut rand_old::thread_rng(),
                &mut sign_key,
                &mut vrfy_key,
            );
            (sign_key, vrfy_key)
        })
        .collect_vec();
    let mut keys1024 = (0..NUM_KEYS)
        .map(|_| {
            let mut sign_key = [0u8; fn_dsa::sign_key_size(fn_dsa::FN_DSA_LOGN_1024)];
            let mut vrfy_key = [0u8; fn_dsa::vrfy_key_size(fn_dsa::FN_DSA_LOGN_1024)];
            fn_dsa::KeyPairGenerator::keygen(
                &mut kg,
                fn_dsa::FN_DSA_LOGN_1024,
                &mut rand_old::thread_rng(),
                &mut sign_key,
                &mut vrfy_key,
            );
            (sign_key, vrfy_key)
        })
        .collect_vec();

    let mut group = c.benchmark_group("fn dsa");
    group.sample_size(NUM_KEYS);
    group.bench_function("keygen 512", |b| {
        b.iter(|| {
            let mut sign_key = [0u8; fn_dsa::sign_key_size(fn_dsa::FN_DSA_LOGN_512)];
            let mut vrfy_key = [0u8; fn_dsa::vrfy_key_size(fn_dsa::FN_DSA_LOGN_512)];
            fn_dsa::KeyPairGenerator::keygen(
                &mut kg,
                fn_dsa::FN_DSA_LOGN_512,
                &mut rand_old::thread_rng(),
                &mut sign_key,
                &mut vrfy_key,
            );
        })
    });
    group.bench_function("keygen 1024", |b| {
        b.iter(|| {
            let mut sign_key = [0u8; fn_dsa::sign_key_size(fn_dsa::FN_DSA_LOGN_1024)];
            let mut vrfy_key = [0u8; fn_dsa::vrfy_key_size(fn_dsa::FN_DSA_LOGN_1024)];
            fn_dsa::KeyPairGenerator::keygen(
                &mut kg,
                fn_dsa::FN_DSA_LOGN_1024,
                &mut rand_old::thread_rng(),
                &mut sign_key,
                &mut vrfy_key,
            );
        })
    });
    group.finish();

    let mut group = c.benchmark_group("fn dsa");
    let mut msgs512 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.random::<[u8; 15]>())
        .collect_vec();
    let mut sigs512 = msgs512
        .iter()
        .enumerate()
        .map(|(i, msg)| {
            let mut sk = <fn_dsa::SigningKeyStandard as fn_dsa::SigningKey>::decode(
                &keys512[i % NUM_KEYS].0,
            )
            .unwrap_or_else(|| panic!("cannot decode"));
            let mut sig = vec![0u8; fn_dsa::signature_size(fn_dsa::SigningKey::get_logn(&sk))];
            fn_dsa::SigningKey::sign(
                &mut sk,
                &mut rand_old::thread_rng(),
                &fn_dsa::DOMAIN_NONE,
                &fn_dsa::HASH_ID_RAW,
                msg,
                &mut sig,
            );
            sig
        })
        .collect_vec();
    let mut msgs1024 = (0..NUM_KEYS * SIGS_PER_KEY)
        .map(|_| rng.random::<[u8; 15]>())
        .collect_vec();
    let mut sigs1024 = msgs1024
        .iter()
        .enumerate()
        .map(|(i, msg)| {
            let mut sk = <fn_dsa::SigningKeyStandard as fn_dsa::SigningKey>::decode(
                &keys1024[i % NUM_KEYS].0,
            )
            .unwrap_or_else(|| panic!("cannot decode"));
            let mut sig = vec![0u8; fn_dsa::signature_size(fn_dsa::SigningKey::get_logn(&sk))];
            fn_dsa::SigningKey::sign(
                &mut sk,
                &mut rand_old::thread_rng(),
                &fn_dsa::DOMAIN_NONE,
                &fn_dsa::HASH_ID_RAW,
                msg,
                &mut sig,
            );
            sig
        })
        .collect_vec();
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_sign_512 = 0;
    group.bench_function("sign 512", |b| {
        b.iter(|| {
            let mut sk = <fn_dsa::SigningKeyStandard as fn_dsa::SigningKey>::decode(
                &keys512[iterator_sign_512 % NUM_KEYS].0,
            )
            .unwrap_or_else(|| panic!("cannot decode"));
            let mut sig = vec![0u8; fn_dsa::signature_size(fn_dsa::SigningKey::get_logn(&sk))];
            fn_dsa::SigningKey::sign(
                &mut sk,
                &mut rand_old::thread_rng(),
                &fn_dsa::DOMAIN_NONE,
                &fn_dsa::HASH_ID_RAW,
                &msgs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)],
                &mut sig,
            );
            sigs512[iterator_sign_512 % (NUM_KEYS * SIGS_PER_KEY)] = sig;
            iterator_sign_512 += 1;
        })
    });
    let mut iterator_sign_1024 = 0;
    group.bench_function("sign 1024", |b| {
        b.iter(|| {
            let mut sk = <fn_dsa::SigningKeyStandard as fn_dsa::SigningKey>::decode(
                &keys1024[iterator_sign_1024 % NUM_KEYS].0,
            )
            .unwrap_or_else(|| panic!("cannot decode"));
            let mut sig = vec![0u8; fn_dsa::signature_size(fn_dsa::SigningKey::get_logn(&sk))];
            fn_dsa::SigningKey::sign(
                &mut sk,
                &mut rand_old::thread_rng(),
                &fn_dsa::DOMAIN_NONE,
                &fn_dsa::HASH_ID_RAW,
                &msgs1024[iterator_sign_1024 % (NUM_KEYS * SIGS_PER_KEY)],
                &mut sig,
            );
            sigs1024[iterator_sign_1024 % (NUM_KEYS * SIGS_PER_KEY)] = sig;
            iterator_sign_1024 += 1;
        })
    });
    group.finish();

    let mut group = c.benchmark_group("fn dsa");
    group.sample_size(NUM_KEYS * SIGS_PER_KEY);
    let mut iterator_verify_512 = 0;
    group.bench_function("verify 512", |b| {
        b.iter(|| {
            let vk = <fn_dsa::VerifyingKeyStandard as fn_dsa::VerifyingKey>::decode(
                &keys512[iterator_verify_512 % NUM_KEYS].1,
            )
            .unwrap_or_else(|| panic!("cannot decode"));
            assert!(fn_dsa::VerifyingKey::verify(
                &vk,
                &sigs512[iterator_verify_512 % sigs512.len()],
                &fn_dsa::DOMAIN_NONE,
                &fn_dsa::HASH_ID_RAW,
                &msgs512[iterator_verify_512 % msgs512.len()]
            ));
            iterator_verify_512 += 1;
        })
    });
    let mut iterator_verify_1024 = 0;
    group.bench_function("verify 1024", |b| {
        b.iter(|| {
            let vk = <fn_dsa::VerifyingKeyStandard as fn_dsa::VerifyingKey>::decode(
                &keys1024[iterator_verify_1024 % NUM_KEYS].1,
            )
            .unwrap_or_else(|| panic!("cannot decode"));
            assert!(fn_dsa::VerifyingKey::verify(
                &vk,
                &sigs1024[iterator_verify_1024 % sigs1024.len()],
                &fn_dsa::DOMAIN_NONE,
                &fn_dsa::HASH_ID_RAW,
                &msgs1024[iterator_verify_1024 % msgs1024.len()]
            ));
            iterator_verify_1024 += 1;
        })
    });
    group.finish();
}

criterion_group!(
    benches,
    falcon_rust_operation,
    falcon_cffi_operation,
    falcon_fndsa_operation
);
criterion_main!(benches);
