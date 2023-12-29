# Falcon-Rust

Unofficial rust implementation of the [Falcon](https://falcon-sign.info/) post-quantum
digital signature scheme.

Falcon was submitted to the [NIST PQC](https://csrc.nist.gov/projects/post-quantum-cryptography)
standardization project and was [selected](https://csrc.nist.gov/Projects/post-quantum-cryptography/selected-algorithms-2022) for 
standardization. The final standard is still outstanding. We do anticipate slight changes
between the standard and the submission, and these changes might break compatibility.

This implementation was written following the [specification](https://falcon-sign.info/falcon.pdf)
and the official [python implementation](https://github.com/tprest/falcon.py).

## Features

 - [x] key generation
 - [x] signature generation
 - [x] signature verification
 - [x] derandomized algorithms
 - [x] (de)serialization
 - [ ] uncompressed signature format
 - [ ] signed-message interface
 - [ ] optimal algorithms (e.g. Karatsuba)
 - [ ] hardware optimizations
 - [ ] message-recovery mode

## To-do's

 - [ ] benchmarking suite
 - [ ] optimize representation of secret key, signature, public key
 - [ ] include usage example in readme
 - [ ] make LdlTree straightforward
 - [ ] test interoperability against the reference implementation
 - [ ] negative tests
 - [ ] Montgomery representation for field elements
 - [ ] Karatsuba for big integer polynomial multiplication

## Contributing

Contributions are welcome! If accepted, contributions will be released under the same
license.
