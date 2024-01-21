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
let mut rng = thread_rng();
let mut msg : [u8; 5] = rng.gen();
let (sk, pk) = falcon512::keygen(rng.gen());
let sig = falcon512::sign(&msg, &sk);
assert!(falcon512::verify(&msg, &sig, &pk));
```

## Performance

Performance is still inferior to the optimized C code accessible from rust via the [foreign function interface](https://crates.io/crates/pqcrypto-falcon) "`pqcrypto-falcon`". These measurements were taken on my Intel(R) Core(TM) i7-10750H CPU @
2.60GHz (which supports AVX2). You can make your own by running `cargo bench`.

|                      | Keygen      | Sign      | Verify    |
|----------------------|-------------|-----------|-----------|
|      falcon-rust 512 | 419.18 ms   | 692.68 µs | 41.668 µs |
|     falcon-rust 1024 |   2.4038 **s**  | 1.3891 **ms** | 86.385 µs |
|  pqcrypto-falcon 512 |   7.5356 ms | 253.44 µs | 48.065 µs |
| pqcrypto-falcon 1024 |  21.454 ms  | 510.43 µs | 94.669 µs |


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
