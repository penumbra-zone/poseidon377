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
mod input;
mod matrix;
mod mds;
mod round_constants;
mod rounds;
mod transcript;
mod utils;

/// For generating parameters at build time.
#[cfg(feature = "std")]
pub mod poseidon_build;

pub use alpha::Alpha;
use ark_ff::PrimeField;
pub use input::InputParameters;
pub use matrix::{
    dot_product, mat_mul, Matrix, MatrixOperations, SquareMatrix, SquareMatrixOperations,
};
pub use mds::{MdsMatrix, OptimizedMdsMatrices};
// use poseidon_parameters::PoseidonParameters;
pub use round_constants::{ArcMatrix, OptimizedArcMatrix};
pub use rounds::RoundNumbers;
pub use utils::log2;

// /// A set of Poseidon parameters for a given set of input parameters.
// #[derive(Clone, Debug)]
// pub struct PoseidonParameters<F: PrimeField> {
//     // Input parameters.
//     /// Security level.
//     pub M: usize,
//     /// Width of desired hash function, e.g. $t=3$ corresponds to a 2-to-1 hash.
//     pub t: usize,

//     // Generated parameters.
//     /// Exponent of the Sbox, i.e. S-box(x) = x^{\alpha} used in the `SubWords` step
//     pub alpha: Alpha,

//     /// Round numbers
//     pub rounds: rounds::RoundNumbers,

//     /// `t x t` MDS matrix used in the `MixLayer` step
//     pub mds: mds::MdsMatrix<F>,

//     /// `num_total_rounds x t` matrix of constants used in the `AddRoundConstant` step
//     pub arc: round_constants::ArcMatrix<F>,

//     /// Optimized round constants.
//     pub optimized_arc: round_constants::OptimizedArcMatrix<F>,

//     /// Optimized MDS matrices.
//     pub optimized_mds: mds::OptimizedMdsMatrices<F>,
// }

pub struct PoseidonParameters<T>(pub T);

impl<F: PrimeField> PoseidonParameters<poseidon_parameters::PoseidonParameters<F>> {
    /// Generate a Poseidon instance mapped over Fp given a choice of:
    ///
    /// * M, the desired security level (in bits),
    /// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    /// * p, the prime modulus,
    /// * `allow_inverse`, whether or not to allow an inverse alpha.
    pub fn new(
        M: usize,
        t: usize,
        p: F::BigInt,
        allow_inverse: bool,
    ) -> poseidon_parameters::PoseidonParameters<F> {
        let input = InputParameters::new(M, t, p, allow_inverse);
        let alpha = alpha::Alpha::generate::<F>(p, allow_inverse);
        let rounds = rounds::RoundNumbers::new(&input, &alpha);
        let mds = mds::MdsMatrix::generate(&input);
        let arc = round_constants::ArcMatrix::generate(&input, rounds, alpha);
        let optimized_mds = mds::OptimizedMdsMatrices::generate(&mds, t, &rounds);
        let optimized_arc = round_constants::OptimizedArcMatrix::generate(&arc, &mds, &rounds);

        poseidon_parameters::PoseidonParameters::<F> {
            M: input.M,
            t: input.t,
            alpha,
            rounds: rounds.0,
            mds,
            arc,
            optimized_mds,
            optimized_arc,
        }

        // Self(instance)
    }
}
