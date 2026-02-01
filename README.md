# Falcon-Rust

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
|      falcon-rust 512 |   226.84 ms     | 329.26 µs | 20.277 µs |
|     falcon-rust 1024 |   1.3106 **s**  | 665.00 µs | 41.933 µs |
|            C FFI 512 |   4.0688 ms     | 129.49 µs | 25.555 µs |
|           C FFI 1024 |   12.048 ms     | 255.59 µs | 50.142 µs |
|           FN DSA 512 |   2.0587 ms     | 143.86 µs | 10.486 µs |
|          FN DSA 1024 |   9.0216 ms     | 278.71 µs | 20.793 µs |


## Features

 - [x] key generation
 - [x] signature generation
 - [x] signature verification
 - [x] derandomized algorithms
 - [x] (de)serialization
 - [ ] better algorithms (e.g. RNS)
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
 - [ ] Montgomery representation for field elements
 - [ ] Residue number system (RNS) for big integer arithmetic
 - [ ] streaming (de)serialization
 - [ ] investigate secret-dependent time variability

## Contributing

Contributions are welcome! If accepted, contributions will be released under the same
license.
