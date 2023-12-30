# Falcon-Rust

Unofficial rust implementation of the [Falcon](https://falcon-sign.info/) post-quantum
digital signature scheme.

Falcon was submitted to the [NIST PQC](https://csrc.nist.gov/projects/post-quantum-cryptography)
standardization project and was [selected](https://csrc.nist.gov/Projects/post-quantum-cryptography/selected-algorithms-2022) for 
standardization. The final standard is still outstanding. We do anticipate slight changes
between the standard and the submission, and these changes might break compatibility.

Falcon comes in two variants. Falcon512 claims at least 108 bits of security, and
Falcon1024 claims at least 252 bits of security, both against quantum computers.

This implementation was written following the [specification](https://falcon-sign.info/falcon.pdf)
and the official [python implementation](https://github.com/tprest/falcon.py).

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
|      falcon-rust 512 | 668.55 ms   | 736.07 µs | 149.55 µs |
|     falcon-rust 1025 |   4.0982 s  | 1.4869 ms | 310.58 µs |
|  pqcrypto-falcon 512 |   7.5356 ms | 253.44 µs | 48.065 µs |
| pqcrypto-falcon 1024 |  21.454 ms  | 510.43 µs | 94.669 µs |

Note that performance improvements may be expected as naïve algorithms are replaced by
better ones.

## Features

 - [x] key generation
 - [x] signature generation
 - [x] signature verification
 - [x] derandomized algorithms
 - [x] (de)serialization
 - [ ] optimal algorithms (e.g. Karatsuba)
 - [ ] uncompressed signature format
 - [ ] signed-message interface
 - [ ] hardware optimizations
 - [ ] message-recovery mode
 - [ ] compatible interface with [`pqcrypto-falcon`](https://crates.io/crates/pqcrypto-falcon)
 - [ ] constant-time (?)

## To-do's

 - [ ] NIST KATs
 - [ ] make LdlTree straightforward
 - [ ] optimize representation of secret key, signature, public key
 - [ ] test interoperability against the reference implementation
 - [ ] negative tests
 - [ ] Montgomery representation for field elements
 - [ ] Karatsuba for big integer polynomial multiplication
 - [ ] streaming (de)serialization
 - [ ] investigate secret-dependent time variability

## Contributing

Contributions are welcome! If accepted, contributions will be released under the same
license.
