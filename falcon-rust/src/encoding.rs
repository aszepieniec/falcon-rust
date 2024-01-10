use bit_vec::BitVec;

/// Compression and decompression routines for signatures.

/// Take as input a list of integers v and a bytelength slen, and
/// return a bytestring of length slen that encode/compress v.
/// If this is not possible, return False.
///
/// For each coefficient of v:
/// - the sign is encoded on 1 bit
/// - the 7 lower bits are encoded naively (binary)
/// - the high bits are encoded in unary encoding
///
/// This method can fail, in which case it returns None. The signature
/// generation algorithm knows this and will re-run the loop.
pub(crate) fn compress(v: &[i16], slen: usize) -> Option<Vec<u8>> {
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

/// Take as input an encoding x, and a length n, and return a list of
/// integers v of length n such that x encode v. If such a list does
/// not exist, the encoding is invalid and we output None.
pub(crate) fn decompress(x: &[u8], bitlength: usize, n: usize) -> Option<Vec<i16>> {
    if x.len() * 8 != bitlength {
        return None;
    }

    let bitvector = BitVec::from_bytes(x);
    let mut index = 0;
    let mut result = Vec::with_capacity(n);
    for _ in 0..n {
        // early return if
        if index + 8 >= bitvector.len() {
            println!(
                "index + 8 = {} but bitvector is only {} long",
                index + 8,
                bitvector.len()
            );
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

#[cfg(test)]
mod test {

    use crate::{
        encoding::{compress, decompress},
        falcon::FalconVariant,
        field::Q,
    };
    use rand::thread_rng;
    use rand_distr::{num_traits::ToPrimitive, Distribution, Normal};

    #[test]
    fn test_compress_decompress() {
        let num_iterations = 200;

        let sigma = 1.5 * (Q.to_f64().unwrap()).sqrt();
        let distribution = Normal::<f64>::new(0.0, sigma).unwrap();
        let mut rng = thread_rng();

        for _ in 0..num_iterations {
            const SALT_LEN: usize = 40;
            const HEAD_LEN: usize = 1;

            // N = 512
            {
                const N: usize = 512;
                const SIG_BYTELEN: usize = 666;
                let slen = SIG_BYTELEN - SALT_LEN - HEAD_LEN;

                let initial: [i16; N] = (0..N)
                    .map(|_| {
                        (distribution.sample(&mut rng) + 0.5)
                            .floor()
                            .to_i32()
                            .unwrap() as i16
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                let compressed = compress(&initial, slen * 8).unwrap();
                let decompressed = decompress(
                    &compressed,
                    (FalconVariant::Falcon512.parameters().sig_bytelen - SALT_LEN - HEAD_LEN) * 8,
                    N,
                )
                .unwrap();
                assert_eq!(initial.to_vec(), decompressed);
            }

            // N = 1024
            {
                const N: usize = 1024;
                const SIG_BYTELEN: usize = 1280;
                let slen = SIG_BYTELEN - SALT_LEN - HEAD_LEN;

                let initial: [i16; 1024] = (0..N)
                    .map(|_| {
                        (distribution.sample(&mut rng) + 0.5)
                            .floor()
                            .to_i32()
                            .unwrap() as i16
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                let compressed = compress(&initial, slen * 8).unwrap();
                let decompressed = decompress(&compressed, slen * 8, N).unwrap();
                assert_eq!(initial.to_vec(), decompressed);
            }
        }
    }
}
