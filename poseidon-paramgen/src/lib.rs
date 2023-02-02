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

use core::ops::Deref;

use alpha::AlphaWrapper;
use ark_ff::PrimeField;
use input::InputParametersWrapper;
use mds::{MdsMatrixWrapper, OptimizedMdsMatricesWrapper};
use poseidon_parameters::PoseidonParameters;
use round_constants::{ArcMatrixWrapper, OptimizedArcMatrixWrapper};
use rounds::RoundNumbersWrapper;
pub use utils::log2;

/// A set of Poseidon parameters for a given set of input parameters.
pub struct PoseidonParametersWrapper<F: PrimeField>(pub PoseidonParameters<F>);

impl<F: PrimeField> From<PoseidonParameters<F>> for PoseidonParametersWrapper<F> {
    fn from(value: PoseidonParameters<F>) -> Self {
        Self(value)
    }
}

impl<F: PrimeField> Deref for PoseidonParametersWrapper<F> {
    type Target = PoseidonParameters<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: PrimeField> PoseidonParametersWrapper<F> {
    /// Generate a Poseidon instance mapped over Fp given a choice of:
    ///
    /// * M, the desired security level (in bits),
    /// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    /// * p, the prime modulus,
    /// * `allow_inverse`, whether or not to allow an inverse alpha.
    #[allow(clippy::new_ret_no_self)]
    pub fn new(M: usize, t: usize, p: F::BigInt, allow_inverse: bool) -> PoseidonParameters<F> {
        let input = InputParametersWrapper::new(M, t, p, allow_inverse);
        let alpha = AlphaWrapper::generate::<F>(p, allow_inverse);
        let rounds = RoundNumbersWrapper::new(&input, &alpha);
        let mds = MdsMatrixWrapper::generate(&input);
        let arc = ArcMatrixWrapper::generate(&input, rounds, alpha);
        let optimized_mds = OptimizedMdsMatricesWrapper::generate(&mds, t, &rounds);
        let optimized_arc = OptimizedArcMatrixWrapper::generate(&arc, &mds, &rounds);

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
}
