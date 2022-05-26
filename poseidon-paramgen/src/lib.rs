#![allow(non_snake_case)]
//! Module for generating Poseidon parameters

mod alpha;
mod matrix;
mod mds;
mod rounds;
mod transcript;
mod utils;

pub use alpha::Alpha;
pub use matrix::{Matrix, SquareMatrix};
pub use mds::MdsMatrix;
pub use rounds::RoundNumbers;
pub use utils::log2;

use ark_ff::{BigInteger, PrimeField};

/// A set of Poseidon parameters for a given set of input parameters.
///
/// This is based on the attacks described in the original Poseidon paper [1].
///
/// References:
/// * Original Poseidon paper: https://eprint.iacr.org/2019/458.pdf
/// * MDS Matrices https://eprint.iacr.org/2020/500/20200702:141143
///
/// TODO(later): Add an Into for converting this to be the ark-sponge Parameter struct.
pub struct PoseidonParameters<F: PrimeField> {
    // Saved input parameters.
    input: InputParameters<F::BigInt>,

    // Generated parameters.
    pub alpha: Alpha,
    pub rounds: rounds::RoundNumbers,
    pub mds: mds::MdsMatrix<F>,
}

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Clone)]
pub struct InputParameters<T: BigInteger> {
    /// Security level in bits.
    pub M: usize,

    /// Width of desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    pub t: usize,

    /// Modulus of the prime field.
    pub p: T, // let modulus = <F as PrimeField>::Params::MODULUS;

    // The below are derived values, stored for convenience.
    /// log_2(p)
    pub log_2_p: f64,
}

impl<T: BigInteger> InputParameters<T> {
    pub fn new(M: usize, t: usize, p: T) -> Self {
        let log_2_p = log2(p);
        InputParameters { M, t, p, log_2_p }
    }
}

impl<F: PrimeField> PoseidonParameters<F> {
    /// Generate a Poseidon instance mapped over Fp given a choice of:
    ///
    /// * M, a desired security level (in bits),
    /// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    /// * p, the prime modulus,
    pub fn new(M: usize, t: usize, p: F::BigInt) -> Self {
        let input = InputParameters::new(M, t, p);
        let alpha = alpha::Alpha::generate::<F>(p);
        let rounds = rounds::RoundNumbers::new(&input, &alpha);
        let mds = mds::MdsMatrix::new(&input);

        Self {
            input,
            alpha,
            rounds,
            mds,
        }
    }
}
