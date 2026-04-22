use bit_vec::BitVec;
use falcon_profiler::profiling;
use itertools::Itertools;
use num::Integer;

/// Take as input a list of integers v and a byte length `byte_length``, and
/// return a bytestring of length `byte_length` that encode/compress v.
/// If this is not possible, return False.
///
/// For each coefficient of v:
/// - the sign is encoded on 1 bit
/// - the 7 lower bits are encoded naively (binary)
/// - the high bits are encoded in unary encoding
///
/// This method can fail, in which case it returns None. The signature
/// generation algorithm knows this and will re-run the loop.
///
/// Algorithm 17 p. 47 of the specification [1].
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[profiling]
pub(crate) fn compress(v: &[i16], byte_length: usize) -> Option<Vec<u8>> {
    // encode each coefficient separately; join later
    let lengths_and_coefficients = v.iter().map(|c| compress_coefficient(*c)).collect_vec();
    let total_length = lengths_and_coefficients
        .iter()
        .map(|(l, _c)| *l)
        .sum::<usize>();

    // if we can't fit all coefficients in the allotted bytes
    if total_length > byte_length * 8 {
        return None;
    }

    // no coefficients are given
    if v.is_empty() {
        return None;
    }

    // join all but one coefficients assuming enough space
    let mut bytes = vec![0u8; byte_length];
    let mut counter = 0;
    for (length, coefficient) in lengths_and_coefficients.iter().take(v.len() - 1) {
        let (cdiv8, cmod8) = counter.div_mod_floor(&8);
        bytes[cdiv8] |= coefficient >> cmod8;
        bytes[cdiv8 + 1] |= ((*coefficient as u16) << (8 - cmod8)) as u8;
        let (cldiv8, clmod8) = (counter + length - 1).div_mod_floor(&8);
        bytes[cldiv8] |= 128u8 >> clmod8;
        bytes[cldiv8 + 1] |= (128u16 << (8 - clmod8)) as u8;
        counter += length;
    }

    // treat last coefficient special
    let (length, coefficient) = lengths_and_coefficients.last().unwrap();
    {
        let (cdiv8, cmod8) = counter.div_mod_floor(&8);
        bytes[cdiv8] |= coefficient >> cmod8;
        bytes[cdiv8 + 1] |= ((*coefficient as u16) << (8 - cmod8)) as u8;
        let (cldiv8, clmod8) = (counter + length - 1).div_mod_floor(&8);
        bytes[cldiv8] |= 128u8 >> clmod8;
        if cldiv8 + 1 < byte_length {
            bytes[cldiv8 + 1] |= (128u16 << (8 - clmod8)) as u8;
        } else if (128u16 << (8 - clmod8)) as u8 != 0 {
            return None;
        }
        counter += length;
    }
    Some(bytes)
}

/// Helper function for compress; isolate attention to one coefficient.
fn compress_coefficient(coeff: i16) -> (usize, u8) {
    let sign = (coeff < 0) as u8;
    let abs = coeff.unsigned_abs();
    let low = abs as u8 & 127;
    let high = abs >> 7;
    (1 + 7 + high as usize + 1, ((sign << 7) | low))
}

///  This is a deprecated decompress routine used now only for testing
/// compatibility with the new, faster implementation (below).
#[allow(dead_code)]
pub(crate) fn decompress_slow(x: &[u8], n: usize) -> Option<Vec<i16>> {
    let bitvector = BitVec::from_bytes(x);
    let mut index = 0;
    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        // early return if
        if index + 8 >= bitvector.len() {
            return None;
        }

        // read sign
        let sign = if bitvector[index] { -1 } else { 1 };
        index += 1;

        // read low bits
        let mut low_bits = 0i16;
        for _ in 0..7 {
            low_bits = (low_bits << 1) | if bitvector[index] { 1 } else { 0 };
            index += 1;
        }

        // read high bits
        let mut high_bits = 0;
        while !bitvector[index] {
            index += 1;
            high_bits += 1;
        }
        index += 1;

        // compose integer and collect it
        let integer = sign * ((high_bits << 7) | low_bits);
        result.push(integer);
    }
    Some(result)
}

/// Take as input an encoding x, and a length n, and return a list of
/// integers v of length n such that x encode v. If such a list does
/// not exist, the encoding is invalid and we output None.
///
/// Algorithm 18 p. 48 of the specification [1].
///
/// [1]: https://falcon-sign.info/falcon.pdf
#[profiling]
pub(crate) fn decompress(x: &[u8], n: usize) -> Option<Vec<i16>> {
    // Bit-accumulator decoder, structured the same as fn-dsa's comp_decode.
    //
    // `acc` is a 32-bit shift register; `acc_len` counts how many of its
    // least-significant bits are valid.  Invariant: acc_len <= 7 at the top
    // of each coefficient iteration.
    //
    // For each coefficient we load exactly one new byte into acc (covering
    // the sign bit and 7 low bits), then refill byte-by-byte as needed while
    // scanning the unary-encoded high bits.  The key win over a per-bit index
    // approach is that acc_len is updated with a simple decrement — no
    // division or modulo by 8 on every bit.
    let mut result = Vec::with_capacity(n);
    let mut i = 0usize;   // byte cursor into x
    let mut acc: u32 = 0;
    let mut acc_len: u32 = 0;

    for _ in 0..n {
        // Load the next byte; it holds the sign bit (MSB) and the 7 low bits
        // of this coefficient's absolute value.
        if i >= x.len() {
            return None;
        }
        acc = (acc << 8) | x[i] as u32;
        i += 1;
        // After the shift, the fresh byte occupies bits acc_len..acc_len+7.
        let s = (acc >> (acc_len + 7)) & 1;  // sign: 0 = positive, 1 = negative
        let mut m = (acc >> acc_len) & 0x7F; // |coeff| low 7 bits

        // Read the unary-encoded high bits: zero or more 0-bits, terminated
        // by a 1-bit.  Each 0-bit adds 128 to m (one high-bit position).
        loop {
            if acc_len == 0 {
                if i >= x.len() {
                    return None;
                }
                acc = (acc << 8) | x[i] as u32;
                i += 1;
                acc_len = 8;
            }
            acc_len -= 1;
            if (acc >> acc_len) & 1 != 0 {
                break; // found the terminating 1
            }
            m += 0x80;
            if m > 2047 {
                return None;
            }
        }

        // Reject the "-0" encoding (sign=negative, value=0).
        // m.wrapping_sub(1) >> 31 is 1 iff m == 0.
        if s & (m.wrapping_sub(1) >> 31) != 0 {
            return None;
        }

        // Apply sign branchlessly: sw = 0 if positive, 0xFFFF_FFFF if negative.
        let sw = s.wrapping_neg();
        result.push(((m ^ sw).wrapping_sub(sw)) as i16);
    }

    // Verify zero-padding: remaining bits in the accumulator and all
    // unread bytes must be zero.
    if acc_len > 0 && acc & ((1 << acc_len) - 1) != 0 {
        return None;
    }
    for &byte in &x[i..] {
        if byte != 0 {
            return None;
        }
    }

    Some(result)
}

#[cfg(test)]
mod test {

    use crate::encoding::{compress, decompress, decompress_slow};
    use crate::falcon_field::Q;
    use bit_vec::BitVec;
    use itertools::Itertools;
    use rand::distr::Distribution;
    use rand::{rng, Rng};

    use proptest::prelude::*;

    /// Compression and decompression routines for signatures.
    ///
    /// This is a deprecated compress routine used now only for testing
    /// compatibility with the new, faster implementation, `compress`.
    #[allow(dead_code)]
    pub(crate) fn compress_slow(v: &[i16], slen: usize) -> Option<Vec<u8>> {
        let mut bitvector: BitVec = BitVec::with_capacity(slen);
        for coeff in v {
            // encode sign
            bitvector.push(*coeff < 0);
            // encode low bits
            let s = (*coeff).abs();
            for i in (0..7).rev() {
                bitvector.push(((s >> i) & 1) != 0);
            }
            // encode high bits
            for _ in 0..(s >> 7) {
                bitvector.push(false);
            }
            bitvector.push(true);
        }
        // return failure if encoding is too long
        if bitvector.len() > slen {
            return None;
        }
        // pad
        while bitvector.len() < slen {
            bitvector.push(false);
        }
        Some(bitvector.to_bytes())
    }

    fn short_elements(n: usize) -> Vec<i16> {
        let sigma = 1.5 * f64::from(Q).sqrt();
        let distribution = rand_distr::Normal::<f64>::new(0.0, sigma).unwrap();
        let mut rng = rng();
        (0..n)
            .map(|_| (distribution.sample(&mut rng) + 0.5).floor() as i16)
            .collect::<Vec<_>>()
    }
    proptest! {
        #[test]
        fn compress_does_not_crash(v in (0..2000usize).prop_map(short_elements)) {
            compress(&v, 2*v.len());
        }
    }
    proptest! {
        #[test]
        fn decompress_recovers(v in (0..2000usize).prop_map(short_elements)) {
            let slen = 2 * v.len();
            let n = v.len();
            if let Some(compressed) = compress(&v, slen) {
                let recovered = decompress(&compressed, n).unwrap();
                prop_assert_eq!(v, recovered.clone());
                let recompressed = compress(&recovered, slen).unwrap();
                prop_assert_eq!(compressed, recompressed);
            }
        }
    }

    #[test]
    fn compress_empty_vec_does_not_crash() {
        compress(&[], 0);
    }

    #[test]
    fn test_compress_decompress() {
        let num_iterations = 1000;

        let sigma = 1.5 * f64::from(Q).sqrt();
        let distribution = rand_distr::Normal::<f64>::new(0.0, sigma).unwrap();
        let mut rng = rng();

        let mut num_successes_512 = 0;
        let mut num_successes_1024 = 0;
        for _ in 0..num_iterations {
            const SALT_LEN: usize = 40;
            const HEAD_LEN: usize = 1;

            // N = 512
            {
                const N: usize = 512;
                const SIG_BYTELEN: usize = 666;
                let slen = SIG_BYTELEN - SALT_LEN - HEAD_LEN;

                let initial: [i16; N] = (0..N)
                    .map(|_| (distribution.sample(&mut rng) + 0.5).floor() as i16)
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                if let Some(compressed) = compress(&initial, slen * 8) {
                    if let Some(decompressed) = decompress(&compressed, N) {
                        assert_eq!(initial.to_vec(), decompressed);
                        num_successes_512 += 1;
                    }
                }
            }

            // N = 1024
            {
                const N: usize = 1024;
                const SIG_BYTELEN: usize = 1280;
                let slen = SIG_BYTELEN - SALT_LEN - HEAD_LEN;

                let initial: [i16; 1024] = (0..N)
                    .map(|_| (distribution.sample(&mut rng) + 0.5).floor() as i16)
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                if let Some(compressed) = compress(&initial, slen * 8) {
                    if let Some(decompressed) = decompress(&compressed, N) {
                        assert_eq!(initial.to_vec(), decompressed);
                        num_successes_1024 += 1;
                    }
                }
            }
        }
        assert!((num_successes_512 as f64) / (num_iterations as f64) > 0.995);
        assert!((num_successes_1024 as f64) / (num_iterations as f64) > 0.995);
    }

    #[test]
    fn test_compress_equiv() {
        let sigma = 1.5 * f64::from(Q).sqrt();
        let distribution = rand_distr::Normal::<f64>::new(0.0, sigma).unwrap();
        let mut rng = rng();

        let n = 200;
        let initial = (0..n)
            .map(|_| (distribution.sample(&mut rng) + 0.5).floor() as i16)
            .collect::<Vec<_>>();
        let slen = 2 * n * 8;
        let compressed = compress_slow(&initial, slen).unwrap();
        let compressed_fast = compress(&initial, slen / 8).unwrap();
        assert_eq!(
            compressed,
            compressed_fast,
            "\n{:#?}\n{:#?}",
            BitVec::from_bytes(&compressed),
            BitVec::from_bytes(&compressed_fast)
        );
    }

    #[test]
    fn test_decompress_equiv() {
        let sigma = 1.5 * f64::from(Q).sqrt();
        let distribution = rand_distr::Normal::<f64>::new(0.0, sigma).unwrap();
        let mut rng = rng();

        let num_iterations = 1000;

        for _ in 0..num_iterations {
            let n = rng.random_range(1..100);
            let initial = (0..n)
                .map(|_| (distribution.sample(&mut rng) + 0.5).floor() as i16)
                .collect::<Vec<_>>();
            let slen = 2 * n * 8;
            let compressed = compress(&initial, slen).unwrap();

            let decompressed = decompress(&compressed, n);
            let decompressed_fast = decompress_slow(&compressed, n);

            assert_eq!(decompressed, decompressed_fast);
        }
    }

    #[test]
    fn test_decompress_failures() {
        let sigma = 1.5 * f64::from(Q).sqrt();
        let distribution = rand_distr::Normal::<f64>::new(0.0, sigma).unwrap();
        let mut rng = rng();

        let num_iterations = 1000;

        for _ in 0..num_iterations {
            let n = rng.random_range(1..100);
            let initial = (0..n)
                .map(|_| (distribution.sample(&mut rng) + 0.5).floor() as i16)
                .collect::<Vec<_>>();
            let slen = 2 * n * 8;
            let compressed = compress(&initial, slen).unwrap();

            assert!(decompress(&compressed, n + 1).is_none());
            // assert!(decompress(&compressed, n - 1).is_none()); // should work

            // flip last set bit -- should cause failure
            let mut compressed_bitvec = BitVec::from_bytes(&compressed);
            let mut index = compressed_bitvec.len();
            while !compressed_bitvec.get(index - 1).unwrap() {
                index -= 1;
            }
            compressed_bitvec.set(index - 1, false);
            let last_bit_flipped = compressed_bitvec.to_bytes();
            assert!(decompress(&last_bit_flipped, n).is_none());

            // try random string -- might fail, but if not must re-encode to the same
            // let random = (0..compressed.len()).map(|_| rng.gen::<u8>()).collect_vec();
            let mut random = compressed.iter().map(|_| rng.random::<u8>()).collect_vec();
            let num_trailing_zeros = compressed
                .iter()
                .cloned()
                .rev()
                .find_position(|&x| x != 0)
                .map(|(pos, _val)| pos)
                .unwrap_or(0);
            let len = random.len();
            for i in 0..num_trailing_zeros {
                random[len - 1 - i] = 0;
            }
            if let Some(decompressed) = decompress(&random, n) {
                let recompressed = compress(&decompressed, slen).unwrap();
                assert_eq!(
                    random,
                    recompressed,
                    "decompressed: {:?}\ndifference: {:?}",
                    decompressed,
                    random
                        .iter()
                        .enumerate()
                        .zip(recompressed.iter().enumerate())
                        .filter(|((_rai, rav), (_rei, rev))| rav != rev)
                        .map(|((rai, rav), (_rei, rev))| format!("{}. {} vs {}", rai, rav, rev))
                        .join(" ")
                );
            }
        }
    }
}
