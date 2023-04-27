#![cfg_attr(not(feature = "std"), no_std)]
#![allow(non_snake_case)]
//! Parameters for Poseidon1 and Poseidon2.
//!
//! The API here is split into [`v1`] and [`v2`] to avoid confusion
//! between the two versions.

mod alpha;
mod arc_matrix;
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
