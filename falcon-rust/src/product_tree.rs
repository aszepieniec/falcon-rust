use core::fmt;
use std::collections::VecDeque;

use itertools::Itertools;
use num::BigUint;
use num::FromPrimitive;
use num::ToPrimitive;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductBranch {
    pub product: BigUint,
    pub left: Box<ProductTree>,
    pub right: Box<ProductTree>,
}

impl ProductBranch {
    fn as_rust_expression(&self) -> String {
        format!(
            "ProductBranch {{ product: BigUint::from_bytes_le(&[{}]), left: Box::new({}), right: Box::new({}) }} ",
            self.product.to_bytes_le().iter().join(", "),
            self.left.as_rust_expression(),
            self.right.as_rust_expression()
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProductTree {
    Leaf(u32),
    Branch(ProductBranch),
}

impl ProductTree {
    pub fn value(&self) -> BigUint {
        match self {
            ProductTree::Leaf(leaf) => BigUint::from_u32(*leaf).unwrap(),
            ProductTree::Branch(branch) => branch.product.clone(),
        }
    }

    pub fn from_leafs(integers: &[u32]) -> ProductTree {
        let mut deque: VecDeque<ProductTree> = integers
            .into_iter()
            .map(|&i| ProductTree::Leaf(i))
            .collect();
        while deque.len() >= 2 {
            let left = deque.pop_front().unwrap();
            let right = deque.pop_front().unwrap();
            let new_node = ProductTree::Branch(ProductBranch {
                product: left.value() * right.value(),
                left: Box::new(left),
                right: Box::new(right),
            });
            deque.push_back(new_node);
        }
        deque[0].clone()
    }

    pub fn reduce(&self, int: &BigUint) -> Vec<u32> {
        match self {
            ProductTree::Leaf(_leaf) => vec![int.to_u32().unwrap()],
            ProductTree::Branch(branch) => {
                let left_remainder = int % branch.left.value();
                let right_remainder = int % branch.right.value();
                [
                    branch.left.reduce(&left_remainder),
                    branch.right.reduce(&right_remainder),
                ]
                .concat()
            }
        }
    }

    pub fn as_branch(&self) -> &ProductBranch {
        match self {
            ProductTree::Leaf(leaf) => {
                panic!("cannot cast product tree as branch when it is a leaf")
            }
            ProductTree::Branch(branch) => branch,
        }
    }

    pub fn as_rust_expression(&self) -> String {
        match self {
            ProductTree::Leaf(leaf) => format!("ProductTree::Leaf({})", *leaf),
            ProductTree::Branch(branch) => {
                format!("ProductTree::Branch({})", branch.as_rust_expression())
            }
        }
    }
}

#[cfg(test)]
mod test {

    use std::ops::Mul;

    use itertools::Itertools;
    use num::BigUint;
    use num::FromPrimitive;
    use num::One;
    use num::ToPrimitive;
    use proptest::arbitrary::Arbitrary;
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest::strategy::BoxedStrategy;
    use proptest::strategy::Just;
    use proptest::strategy::Strategy;
    use test_strategy::proptest as strategy_proptest;

    use crate::product_tree::ProductTree;

    fn arbitrary_biguint(bitlen: usize) -> BoxedStrategy<BigUint> {
        let limbs = vec(u32::arbitrary(), bitlen.div_ceil(32));
        limbs
            .prop_map(move |bb| BigUint::from_slice(&bb) >> (bb.len() * 32 - bitlen))
            .boxed()
    }

    #[strategy_proptest]
    fn cumulative_product(
        #[strategy(1usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(vec(u32::arbitrary(), #_n))] ints: Vec<u32>,
    ) {
        let cumulative_product = ints
            .iter()
            .cloned()
            .map(|i| BigUint::from_u32(i).unwrap())
            .fold(BigUint::one(), BigUint::mul);
        let tree = ProductTree::from_leafs(&ints);
        prop_assert_eq!(cumulative_product, tree.value());
    }

    #[strategy_proptest]
    fn reduce_individually(
        #[strategy(1usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(vec(u32::arbitrary(), #_n))] moduli: Vec<u32>,
        #[strategy(arbitrary_biguint(#_n*30))] integer: BigUint,
    ) {
        let tree = ProductTree::from_leafs(&moduli);
        let leaf_remainders = tree.reduce(&integer);
        let individual_remainders = moduli
            .into_iter()
            .map(|m| BigUint::from_u32(m).unwrap())
            .map(|p| integer.clone() % p)
            .map(|b| b.to_u32().unwrap())
            .collect_vec();
        prop_assert_eq!(individual_remainders, leaf_remainders);
    }
}
