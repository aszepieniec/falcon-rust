use std::ops::MulAssign;

use itertools::Itertools;
use num::{BigUint, FromPrimitive, One, Zero};

#[derive(Debug, Clone)]
pub(crate) struct ModularInversePair {
    pub(crate) zeros_then_one: BigUint,
    pub(crate) ones_then_zero: BigUint,
}

impl ModularInversePair {
    pub fn rust_expression(&self) -> String {
        format!("ModularInversePair {{\nzeros_then_one: BigUint::from_bytes_le(&[{}]),\nones_then_zeros: BigUint::from_bytes_le(&[{}]),}}", self.zeros_then_one.to_bytes_le().iter().join(", "), self.ones_then_zero.to_bytes_le().iter().join(", "))
    }
}

#[derive(Debug, Clone)]
pub(crate) struct ModularInversesSequence(pub(crate) Vec<ModularInversePair>);

impl ModularInversesSequence {
    pub fn from_moduli(moduli: &[u32]) -> Self {
        let mut partial_products = vec![BigUint::one()];
        for &modulus in moduli.iter() {
            partial_products
                .push(partial_products.last().unwrap() * BigUint::from_u32(modulus).unwrap());
        }

        let mut sequence = vec![];
        for (i, &current_modulus) in moduli.iter().enumerate() {
            let current_modulus_as_biguint = BigUint::from_u32(current_modulus).unwrap();
            let current_modulus_inverse = current_modulus_as_biguint
                .modinv(&partial_products[i])
                .unwrap();
            let running_product_inverse = partial_products[i]
                .modinv(&current_modulus_as_biguint)
                .unwrap();
            sequence.push(ModularInversePair {
                zeros_then_one: running_product_inverse * &partial_products[i],
                ones_then_zero: current_modulus_inverse * current_modulus_as_biguint,
            });
        }
        Self(sequence)
    }

    pub fn construct(&self, modular_residues: &[u32]) -> BigUint {
        let mut acc = BigUint::zero();
        for i in 0..modular_residues.len() {
            let pair = &self.0[i];
            acc = acc * &pair.ones_then_zero + modular_residues[i] * &pair.zeros_then_one;
        }
        acc
    }

    pub fn rust_expression(&self) -> String {
        format!(
            "ModularInversesSequence(vec![{}])",
            self.0.iter().map(|pair| pair.rust_expression()).join(",\n")
        )
    }
}

#[cfg(test)]
mod test {
    use itertools::Itertools;
    use num::{BigUint, FromPrimitive};
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use test_strategy::proptest as strategy_proptest;

    use super::ModularInversesSequence;

    const MODULI: [u32; 6] = [
        1073754113u32,
        1073950721u32,
        1073958913u32,
        1073983489u32,
        1074196481u32,
        1074343937u32,
    ];

    #[strategy_proptest]
    fn constructed_integer_agrees_modulo_primes(
        #[strategy(vec(0u32..MODULI[0], MODULI.len()))] residues: Vec<u32>,
    ) {
        let sequence = ModularInversesSequence::from_moduli(&MODULI);
        let integer = sequence.construct(&residues);

        println!("residues: {}", residues.iter().join(","));
        println!("integer: {}", integer);

        for (expected_residue, &p) in residues.into_iter().zip(MODULI.iter()) {
            let observed_residue = *(integer.clone() % BigUint::from_u32(p).unwrap())
                .to_u32_digits()
                .first()
                .unwrap_or(&0);
            prop_assert_eq!(expected_residue, observed_residue);
        }
    }
}
