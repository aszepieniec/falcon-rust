//! Unofficial rust implementation of the [Falcon] post-quantum
//! digital signature scheme.
//!
//! Falcon was submitted to the [NIST PQC]
//! standardization project and was [selected] for
//! standardization. The final standard is still outstanding. We do anticipate slight changes
//! between the standard and the submission, and these changes might break compatibility.
//!
//! Falcon comes in two variants. Falcon512 claims at least 108 bits of security, and
//! Falcon1024 claims at least 252 bits of security, both against quantum computers.
//!
//! This implementation was written following the [specification]
//! and using the [python implementation] as a guide, although later versions diverge from this
//! reference point.
//!
//! [Falcon]: https://falcon-sign.info/
//! [NIST PQC]: https://csrc.nist.gov/projects/post-quantum-cryptography
//! [selected]: https://csrc.nist.gov/Projects/post-quantum-cryptography/selected-algorithms-2022
//! [specification]: https://falcon-sign.info/falcon.pdf
//! [python implementation]: https://github.com/tprest/falcon.py
//!
//! # Usage
//!
//! First, `falcon-rust = "0.1.3"` to your `Cargo.toml` file.
//!
//! Then to use the interface:
//! ```
//! use falcon_rust::falcon512;
//!
//! use rand::rng;
//! use rand::RngExt;
//!
//! let msg = b"Hello, world!";
//! let (sk, pk) = falcon512::keygen(rng().random());
//! let sig = falcon512::sign(msg, &sk);
//! assert!(falcon512::verify(msg, &sig, &pk));
//! ```
//!
//! For serialization / deserialization:
//! ```
//! use falcon_rust::falcon512;
//!
//! use rand::rng;
//! use rand::RngExt;
//!
//! let msg = b"Hello, world!";
//! let (sk, pk) = falcon512::keygen(rng().random());
//! let sig = falcon512::sign(msg, &sk);
//!
//! let sk_buffer = sk.to_bytes();
//! let pk_buffer = pk.to_bytes();
//! let sig_buffer = sig.to_bytes();
//! falcon512::SecretKey::from_bytes(&sk_buffer);
//! falcon512::PublicKey::from_bytes(&pk_buffer);
//! falcon512::Signature::from_bytes(&sig_buffer);
//! ```

pub(crate) mod cyclotomic_fourier;
pub(crate) mod rns;
pub(crate) mod encoding;
pub(crate) mod falcon;
pub(crate) mod fixed_point;
pub(crate) mod fp_field;
pub(crate) mod packed;
pub mod falcon1024;
pub mod falcon512;
pub mod falcon_field;
pub(crate) mod fast_fft;
pub(crate) mod ffsampling;
pub(crate) mod inverse;
pub mod math; // pub for benching
pub mod polynomial; // pub for benching
pub(crate) mod samplerz;
// The depth-0 NTRU-solve NTT is RNS with a single prime; reuse a 24-bit prime
// from the multi-prime RNS list (`NttPrimes24Bit*`) rather than a bespoke 30-bit
// one.  The depth-0 products measure ~2^19.5, well within this prime's ~2^22
// signed reconstruction range; `ntru_gen` self-checks `fG − gF = q` to catch any
// rare overflow.
pub(crate) type U32Field = fp_field::FpField<8404993>;

#[cfg(feature = "profiling")]
pub mod profiling;
