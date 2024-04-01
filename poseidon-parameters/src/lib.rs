#![cfg_attr(not(feature = "std"), no_std)]
#![allow(non_snake_case)]
//! Parameters for Poseidon1 and Poseidon2.
//!
//! The API here is split into [`v1`] and [`v2`] to avoid confusion
//! between the two versions.

use decaf377::Fq;

mod alpha;
mod arc_matrix;
mod error;
mod matrix;
mod matrix_ops;
mod mds_matrix;
mod round_numbers;

#[cfg(test)]
mod tests;

/// Structures related to Poseidon version 1 parameters.
pub mod v1;

/// Structures related to Poseidon version 2 parameters.
pub mod v2;

pub trait StuffThatNeedsToGoInDecaf377 {
    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self;
}

// TEMP
impl StuffThatNeedsToGoInDecaf377 for Fq {
    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Fq::from(1u64);
        let exp_u64 = exp.as_ref();
        for _ in 0..exp_u64[0] {
            res *= self;
        }
        res
    }
}
