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

mod alpha;
mod appendix_g;
pub(crate) mod input;
mod mds;
mod round_constants;
mod rounds;
mod transcript;
mod utils;

/// For generating parameters at build time.
#[cfg(feature = "std")]
pub mod poseidon_build;

use ark_ff::PrimeField;
use poseidon_parameters::PoseidonParameters;
use utils::log2;

/// Generate a Poseidon instance mapped over Fp given a choice of:
///
/// * M, the desired security level (in bits),
/// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
/// * p, the prime modulus,
/// * `allow_inverse`, whether or not to allow an inverse alpha.
pub fn generate<F: PrimeField>(
    M: usize,
    t: usize,
    p: F::BigInt,
    allow_inverse: bool,
) -> PoseidonParameters<F> {
    let input = input::generate(M, t, p, allow_inverse);
    let alpha = alpha::generate::<F>(p, allow_inverse);
    let rounds = rounds::generate(&input, &alpha);
    let mds = mds::generate(&input);
    let arc = round_constants::generate(&input, rounds, alpha);
    let optimized_mds = mds::generate_optimized(&mds, t, &rounds);
    let optimized_arc = round_constants::generate_optimized(&arc, &mds, &rounds);

    PoseidonParameters::<F> {
        M: input.M,
        t: input.t,
        alpha,
        rounds,
        mds,
        arc,
        optimized_mds,
        optimized_arc,
    }
}
