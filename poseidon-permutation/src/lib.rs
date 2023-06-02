//! An implemention of the Poseidon permutation for fixed-width
//! hashing.
#![cfg_attr(not(feature = "std"), no_std)]

mod permutation;

pub use permutation::Instance;

#[cfg(feature = "r1cs")]
pub mod r1cs;
