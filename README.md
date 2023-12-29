# Falcon-Rust

Unofficial rust implementation of the [Falcon](https://falcon-sign.info/) post-quantum
digital signature scheme.

Falcon was submitted to the NIST PQC standardization project and was selected for 
standardization. The final standard is still outstanding. We do anticipate slight changes
between the standard and the submission, and these changes might break compatibility.

## Features

 - [x] key generation
 - [x] signature generation
 - [x] signature verification
 - [x] derandomized algorithms
 - [ ] serialization
 - [ ] uncompressed signature format
 - [ ] signed-message interface
 - [ ] optimal algorithms (e.g. Karatsuba)
 - [ ] hardware optimizations

## To-do's

 - [ ] serialization
 - [ ] test against the reference implementation
 - [ ] benchmarking suite
 - [ ] Montgomery representation for field elements
 - [ ] Karatsuba for big integer polynomial multiplication

## Contributing

Contributions are welcome! If accepted, contributions will be released under the same
license.
