use std::ops::{Add, Neg, Sub};

use itertools::Itertools;

use crate::field::Felt;

#[derive(Debug, Clone)]
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
