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
pub use matrix::{
    dot_product, mat_mul, Matrix, MatrixOperations, SquareMatrix, SquareMatrixOperations,
};
pub use mds::{MdsMatrix, OptimizedMdsMatrices};
use num::BigUint;
pub use round_constants::{ArcMatrix, OptimizedArcMatrix};
pub use rounds::RoundNumbers;
pub use utils::log2;

use ark_ff::PrimeField;
use ark_sponge::poseidon::Parameters as ArkPoseidonParameters;

/// A set of Poseidon parameters for a given set of input parameters.
#[derive(Clone, Debug)]
pub struct PoseidonParameters<F: PrimeField> {
    /// Saved input parameters.
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

    /// Optimized round constants.
    pub optimized_arc: round_constants::OptimizedArcMatrix<F>,

    /// Optimized MDS matrices.
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
        let mds = mds::MdsMatrix::generate(&input);
        let arc = round_constants::ArcMatrix::generate(&input, rounds, alpha);
        let optimized_mds = mds::OptimizedMdsMatrices::generate(&mds, t, &rounds);
        let optimized_arc =
            round_constants::OptimizedArcMatrix::generate(&arc, &optimized_mds, &rounds);

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

    /// Instantiate poseidon parameters from existing unoptimized components.
    ///
    /// # Safety
    ///
    /// This method does not check the choices are secure. It is the caller's
    /// responsibility to ensure that the components passed are secure.
    pub fn from_unoptimized_components(
        alpha: u32,
        mds_elements: Vec<Vec<F>>,
        arc_elements: Vec<Vec<F>>,
        M: usize,
        t: usize,
        p: F::BigInt,
        r_F: usize,
        r_P: usize,
    ) -> Self {
        let input = InputParameters::new(M, t, p, false);
        let alpha = alpha::Alpha::Exponent(alpha);
        let rounds = rounds::RoundNumbers::from_rounds(r_F, r_P);
        let mds = mds::MdsMatrix::new(t, t, flatten(mds_elements));
        let arc = round_constants::ArcMatrix::new(r_P + r_F, t, flatten(arc_elements));
        dbg!("MDS!");
        for elem in mds.elements().into_iter() {
            // We use the BigUint type here since the Display of the field element
            // is not in decimal: see https://github.com/arkworks-rs/algebra/issues/320
            let elem_bigint: BigUint = (*elem).into();
            dbg!("{} ", elem_bigint.to_string());
        }
        let optimized_mds = mds::OptimizedMdsMatrices::generate(&mds, t, &rounds);
        let optimized_arc =
            round_constants::OptimizedArcMatrix::generate(&arc, &optimized_mds, &rounds);

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

/// Flatten a vec of vecs
fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
    nested.into_iter().flatten().collect()
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
