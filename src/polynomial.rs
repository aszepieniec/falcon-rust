use sha3::{digest::*, Shake256};
use std::ops::{Add, Neg, Sub};

use itertools::Itertools;

use crate::field::{Felt, Q};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial {
    pub coefficients: Vec<Felt>,
}

impl Polynomial {
    pub fn new(coefficients: &[Felt]) -> Self {
        Self {
            coefficients: coefficients.to_vec(),
        }
    }
}

impl Add for Polynomial {
    type Output = Polynomial;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            coefficients: self
                .coefficients
                .iter()
                .zip_eq(rhs.coefficients.iter())
                .map(|(a, b)| *a + *b)
                .collect(),
        }
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            coefficients: self
                .coefficients
                .iter()
                .zip_eq(rhs.coefficients.iter())
                .map(|(a, b)| *a - *b)
                .collect(),
        }
    }
}

impl Neg for Polynomial {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::Output {
            coefficients: self.coefficients.iter().cloned().map(|a| -a).collect(),
        }
    }
}

/// Hash a string to a random polynomial in ZZ[X] mod <Phi(X), q>.
/// Algorithm 3, "HashToPoint" in the spec (page 31).
pub fn hash_to_point(string: &[u8], n: usize) -> Polynomial {
    const K: u32 = (1u32 << 16) / (Q as u32);

    let mut hasher = Shake256::default();
    hasher.update(string);
    let mut reader = hasher.finalize_xof();

    let mut coefficients: Vec<Felt> = vec![];
    while coefficients.len() != n {
        let mut randomness = [0u8; 2];
        reader.read(&mut randomness);
        // Arabic endianness but so be it
        let t = ((randomness[0] as i32) << 8) | (randomness[1] as i32);
        if (t as u32) < K * (Q as u32) {
            coefficients.push(Felt::new((t % Q) as i16));
        }
    }

    Polynomial { coefficients }
}

#[cfg(test)]
mod test {
    use sha3::digest::ExtendableOutput;
    use sha3::digest::Update;
    use sha3::digest::XofReader;
    use sha3::Shake256;

    use crate::field::Felt;
    use crate::polynomial::hash_to_point;
    use crate::polynomial::Polynomial;

    #[test]
    fn test_shake256() {
        let input = b"\x21\xF1\x34\xAC\x57";
        let kat_output : [u8;64] = *b"\xBB\x8A\x84\x47\x51\x7B\xA9\xCA\x7F\xA3\x4E\xC9\x9A\x80\x00\x4F\x22\x8A\xB2\x82\x47\x28\x41\xEB\x3D\x3A\x76\x22\x5C\x9D\xBE\x77\xF7\xE4\x0A\x06\x67\x76\xD3\x2C\x74\x94\x12\x02\xF9\xF4\xAA\x43\xD1\x2C\x62\x64\xAF\xA5\x96\x39\xC4\x4E\x11\xF5\xE1\x4F\x1E\x56";
        let mut observed_output = [0u8; 64];

        let mut hasher = Shake256::default();
        hasher.update(input);
        let mut reader = hasher.finalize_xof();
        reader.read(&mut observed_output);

        assert_eq!(
            kat_output, observed_output,
            "SHAKE256 from crate sha3 does not match KAT"
        )
    }

    #[test]
    fn test_hash_to_point_sanity() {
        // KAT obtained from python script, run in the same directory
        // as a copy of the official python implementation.
        // ```
        // from falcon import SecretKey
        // sk = SecretKey(512)
        // sk.hash_to_point(b"", b"")
        // ```
        let hash_of_empty = Polynomial {
            coefficients: [
                5816, 7463, 2984, 11537, 9019, 4074, 5180, 11040, 4044, 8937, 694, 7042, 9481,
                10084, 3795, 5677, 5977, 1241, 6332, 2817, 413, 1971, 755, 7241, 6041, 9347, 4136,
                11948, 9140, 1210, 5150, 1630, 4015, 2390, 2346, 2025, 4272, 10978, 7171, 8764,
                11920, 888, 12160, 11275, 7043, 10323, 1181, 1873, 5876, 12213, 627, 11319, 5675,
                8207, 6210, 385, 333, 4581, 1359, 10859, 3346, 3970, 8720, 3640, 8157, 1080, 2794,
                5769, 11618, 6780, 1734, 6484, 1575, 9433, 10353, 2004, 5921, 5013, 4753, 9865,
                10931, 6621, 1417, 9804, 12027, 9437, 10657, 3260, 9541, 4967, 12124, 6827, 333,
                6404, 3498, 6920, 3979, 14, 440, 1293, 8011, 7567, 3899, 3252, 4023, 10727, 11938,
                957, 2412, 9552, 10409, 8063, 9131, 9835, 10305, 3124, 6303, 12241, 6354, 2540,
                10113, 10777, 6803, 4879, 10952, 10503, 1728, 5067, 3339, 7045, 11333, 5469, 11062,
                11666, 5235, 2314, 3345, 2224, 2274, 8060, 4304, 6716, 11595, 1541, 996, 6983, 36,
                449, 7401, 4987, 9177, 810, 1908, 8650, 7646, 6893, 4919, 1971, 4930, 11763, 201,
                12223, 9234, 4081, 6199, 12047, 7646, 9639, 6814, 6739, 5279, 4012, 2101, 10707,
                4241, 12146, 3779, 3999, 3176, 1699, 10294, 5168, 5590, 457, 9709, 6450, 442, 8884,
                6755, 10995, 10923, 3935, 8499, 3508, 12088, 1115, 11336, 1379, 7557, 4954, 7639,
                2514, 8672, 6686, 98, 5676, 8800, 5712, 4724, 7724, 3202, 12128, 10940, 10177,
                9421, 11013, 7372, 8546, 441, 6261, 8779, 2453, 12082, 7922, 5307, 6920, 7726, 823,
                10561, 1251, 10358, 8383, 10905, 8145, 1733, 1718, 3105, 10756, 6798, 10209, 7976,
                11148, 9353, 4746, 1089, 11444, 6571, 409, 8381, 10325, 7649, 10042, 5587, 3625,
                10182, 10494, 228, 4687, 5949, 7995, 12092, 3312, 5339, 5920, 8145, 6796, 1992,
                3205, 2761, 12199, 11342, 9695, 390, 252, 989, 1385, 12148, 8324, 10694, 3690,
                3440, 8888, 12238, 9018, 3354, 5859, 6298, 8098, 4388, 3788, 3045, 11095, 2372,
                10036, 9233, 168, 8500, 3604, 2494, 9854, 5679, 2182, 3350, 7798, 8310, 3544, 948,
                7646, 7235, 2650, 6008, 4610, 2159, 6884, 10545, 688, 4115, 10312, 4408, 4951,
                2891, 9791, 1377, 5909, 11147, 11139, 4969, 5158, 350, 1067, 4242, 10820, 1818,
                6473, 105, 2919, 10892, 7116, 850, 11409, 2652, 6392, 2540, 6892, 8372, 11975,
                4994, 2621, 2763, 11837, 6132, 11293, 9138, 8769, 10964, 9826, 601, 7007, 9078,
                10340, 9410, 8746, 10835, 9053, 11010, 5308, 8851, 1976, 11016, 599, 8348, 9876,
                7100, 1333, 4550, 1735, 4598, 9970, 525, 8320, 1609, 9213, 4178, 484, 10814, 1760,
                9667, 8369, 2286, 10384, 12139, 24, 1178, 5682, 7074, 3676, 3661, 3322, 1831, 5562,
                734, 8059, 8750, 6951, 4760, 10395, 1019, 9404, 2923, 6715, 123, 10157, 4892, 7667,
                1677, 4175, 3455, 12123, 10730, 2000, 8212, 2665, 7088, 8741, 10936, 3172, 225,
                3867, 5140, 2310, 6453, 2898, 3637, 4580, 113, 5991, 3532, 3363, 11457, 11601,
                7280, 6792, 11872, 8127, 2192, 10761, 9019, 8197, 8965, 6061, 10799, 988, 10522,
                1281, 1965, 2716, 9841, 7496, 8456, 5192, 825, 3727, 4664, 7906, 8521, 5901, 10200,
                5167, 9451, 10825, 12011, 2272, 8698, 8174, 11973, 5155, 6890, 9999, 4391, 12044,
                1620, 8310, 111, 4481, 9650, 2077, 7691, 7531, 1956, 494, 3297, 1623, 3266, 7018,
                2031, 6317, 4657, 5206, 2581, 11227, 10508, 4567, 8892, 1363, 6790, 6180, 1588,
                9776, 11998, 10689, 492, 331,
            ]
            .map(Felt::new)
            .to_vec(),
        };
        assert_eq!(hash_of_empty, hash_to_point(&[], 512));
    }
}
