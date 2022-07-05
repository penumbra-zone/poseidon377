#![allow(non_snake_case)]
//! Module for generating Poseidon parameters

mod alpha;
mod input;
mod matrix;
mod mds;
mod round_constants;
mod rounds;
mod transcript;
mod utils;

pub use alpha::Alpha;
pub use input::InputParameters;
pub use matrix::{Matrix, SquareMatrix};
pub use mds::{MdsMatrix, OptimizedMdsMatrices};
pub use round_constants::{ArcMatrix, OptimizedArcMatrix};
pub use rounds::RoundNumbers;
pub use utils::log2;

use ark_ff::PrimeField;
use ark_sponge::poseidon::Parameters as ArkPoseidonParameters;

/// A set of Poseidon parameters for a given set of input parameters.
#[derive(Clone, Debug)]
pub struct PoseidonParameters<F: PrimeField> {
    // Saved input parameters.
    pub input: InputParameters<F::BigInt>,

    // Generated parameters.
    /// Exponent of the Sbox, i.e. S-box(x) = x^{\alpha} used in the `SubWords` step
    pub alpha: Alpha,

    /// Round numbers
    pub rounds: rounds::RoundNumbers,

    /// `t x t` MDS matrix used in the `MixLayer` step
    pub mds: mds::MdsMatrix<F>,

    /// `num_total_rounds x t` matrix of constants used in the `AddRoundConstant` step
    pub arc: round_constants::ArcMatrix<F>,

    /// Optional optimizations
    pub optimized_arc: round_constants::OptimizedArcMatrix<F>,
    pub optimized_mds: mds::OptimizedMdsMatrices<F>,
}

impl<F: PrimeField> PoseidonParameters<F> {
    /// Generate a Poseidon instance mapped over Fp given a choice of:
    ///
    /// * M, the desired security level (in bits),
    /// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    /// * p, the prime modulus,
    /// * `allow_inverse`, whether or not to allow an inverse alpha.
    pub fn new(M: usize, t: usize, p: F::BigInt, allow_inverse: bool) -> Self {
        let input = InputParameters::new(M, t, p, allow_inverse);
        let alpha = alpha::Alpha::generate::<F>(p, allow_inverse);
        let rounds = rounds::RoundNumbers::new(&input, &alpha);
        let mds = mds::MdsMatrix::new(&input);
        let arc = round_constants::ArcMatrix::generate(&input, rounds, alpha);
        let optimized_mds = mds::OptimizedMdsMatrices::generate(&mds, t);
        let optimized_arc =
            round_constants::OptimizedArcMatrix::generate(&arc, &optimized_mds, &rounds, t);

        Self {
            input,
            alpha,
            rounds,
            mds,
            arc,
            optimized_mds,
            optimized_arc,
        }
    }
}

impl<F: PrimeField> Into<ArkPoseidonParameters<F>> for PoseidonParameters<F> {
    fn into(self) -> ArkPoseidonParameters<F> {
        let alpha = match self.alpha {
            Alpha::Exponent(exp) => exp as u64,
            Alpha::Inverse => panic!("ark-sponge does not allow inverse alpha"),
        };
        // TODO: let user specify different capacity choices
        let capacity = 1;
        let rate = self.input.t - capacity;

        ArkPoseidonParameters {
            full_rounds: self.rounds.full(),
            partial_rounds: self.rounds.partial(),
            alpha,
            ark: self.arc.into(),
            mds: self.mds.into(),
            rate,
            capacity,
        }
    }
}
