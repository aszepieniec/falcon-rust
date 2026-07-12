# Falcon-Rust

![GitHub License](https://img.shields.io/github/license/aszepieniec/falcon-rust)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/aszepieniec/falcon-rust/rust.yml)
![Crates.io Version](https://img.shields.io/crates/v/falcon-rust)
![docs.rs](https://img.shields.io/docsrs/falcon-rust)
![Crates.io Downloads (latest version)](https://img.shields.io/crates/dv/falcon-rust)
![Deps.rs Crate Dependencies (latest)](https://img.shields.io/deps-rs/falcon-rust/latest)
![Crates.io Dependents](https://img.shields.io/crates/dependents/falcon-rust)

Unofficial rust implementation of the [Falcon](https://falcon-sign.info/) post-quantum
digital signature scheme.

Falcon was submitted to the [NIST PQC](https://csrc.nist.gov/projects/post-quantum-cryptography)
standardization project and was [selected](https://csrc.nist.gov/Projects/post-quantum-cryptography/selected-algorithms-2022) for 
standardization. The final standard is still outstanding. We do anticipate slight changes
between the standard and the submission, and these changes might break compatibility.

Falcon comes in two variants. Falcon512 claims at least 108 bits of security, and
Falcon1024 claims at least 252 bits of security, both against quantum computers.

This implementation adheres to the [specification](https://falcon-sign.info/falcon.pdf). It was originally written following the the official [python implementation](https://github.com/tprest/falcon.py), but has since deviated.

## Example

```rust
use rand::rng;
use rand::RngExt;
use falcon_rust::falcon512;

let mut rng = rng();
let mut msg : [u8; 5] = rng.random();
let (sk, pk) = falcon512::keygen(rng.random());
let sig = falcon512::sign(&msg, &sk);
assert!(falcon512::verify(&msg, &sig, &pk));
```

## Performance

If you are after performance, you are probably better off with one of the implementations by the inventors, either the foreign function interface (FFI) into the optimized C code ([`pqcrypto-falcon`](https://crates.io/crates/pqcrypto-falcon)), or the optimized rust crate ([`fn-dsa`](https://crates.io/crates/fn-dsa)). The following benchmark was produced by my 12th Gen Intel(R) Core(TM) i9-12900K (which supports AVX2). You can make your own by running `cargo bench`.

|                      | Keygen          | Sign      | Verify    |
|----------------------|-----------------|-----------|-----------|
|      falcon-rust 512 |   14.272 ms     | 456.44 µs | 15.573 µs |
|     falcon-rust 1024 |   30.498 ms     | 941.82 µs | 31.370 µs |
|            C FFI 512 |   3.9177 ms     | 113.53 µs | 24.685 µs |
|           C FFI 1024 |   11.532 ms     | 225.28 µs | 48.562 µs |
|           FN DSA 512 |   2.0290 ms     | 154.70 µs | 10.343 µs |
|          FN DSA 1024 |   9.0410 ms     | 300.76 µs | 20.696 µs |


## Features

 - [x] key generation
 - [x] signature generation
 - [x] signature verification
 - [x] derandomized algorithms
 - [x] (de)serialization
 - [x] Montgomery representation
 - [x] residue number system
 - [ ] uncompressed signature format
 - [ ] signed-message interface
 - [ ] hardware optimizations
 - [ ] message-recovery mode
 - [ ] constant-time (?)

## To-do's

 - [ ] NIST KATs
 - [ ] make LdlTree straightforward
 - [ ] optimize representation of secret key, signature, public key
 - [ ] test interoperability against the reference implementation
 - [ ] negative tests
 - [ ] profile, and fix bottlenecks
 - [ ] streaming (de)serialization
 - [ ] investigate secret-dependent time variability

## Contributing

Contributions are welcome! If accepted, contributions will be released under the same
license.
