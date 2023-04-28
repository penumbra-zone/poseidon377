#![cfg_attr(not(feature = "std"), no_std)]
#![allow(non_snake_case)]
#![deny(missing_docs)]
//! Module for generating parameters for the Poseidon SNARK-friendly hash function.
//!
//! This crate will, given a choice of:
//!
//! * M, the desired security level (in bits),
//! * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
//! * p, the prime modulus,
//! * `allow_inverse`, whether or not to allow an inverse alpha for the Sbox layer.
//!
//! generate the best choice of parameters, for both the unoptimized version of Poseidon
//! specified in the [Poseidon paper], as well as the optimizations described in Appendix
//! B.
//!
//! [Poseidon paper]: https://eprint.iacr.org/2019/458.pdf

#[cfg(not(feature = "std"))]
extern crate alloc;

/// Logic for generating Poseidon1 parameters.
pub mod v1;

/// Logic for generating Poseidon2 parameters.
pub mod v2;

mod alpha;
mod appendix_g;
mod input;
mod mds;
mod round_constants;
mod rounds;
mod transcript;
mod utils;

/// For generating parameters at build time.
#[cfg(feature = "std")]
mod poseidon_build;

use utils::log2;
