use crate::r1cs::vendor::sponge::poseidon::PoseidonParameters;
use ark_ff::{fields::models::*, FpConfig, PrimeField};
use ark_std::{vec, vec::Vec};

/// An entry in the default Poseidon parameters
pub struct PoseidonDefaultParametersEntry {
    /// The rate (in terms of number of field elements).
    pub rate: usize,
    /// Exponent used in S-boxes.
    pub alpha: usize,
    /// Number of rounds in a full-round operation.
    pub full_rounds: usize,
    /// Number of rounds in a partial-round operation.
    pub partial_rounds: usize,
    /// Number of matrices to skip when generating parameters using the Grain LFSR.
    ///
    /// The matrices being skipped are those that do not satisfy all the desired properties.
    /// See the [reference implementation](https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/generate_parameters_grain.sage) for more detail.
    pub skip_matrices: usize,
}

impl PoseidonDefaultParametersEntry {
    /// Create an entry in PoseidonDefaultParameters.
    pub const fn new(
        rate: usize,
        alpha: usize,
        full_rounds: usize,
        partial_rounds: usize,
        skip_matrices: usize,
    ) -> Self {
        Self {
            rate,
            alpha,
            full_rounds,
            partial_rounds,
            skip_matrices,
        }
    }
}

/// A field with Poseidon parameters associated
pub trait PoseidonDefaultParametersField: PrimeField {
    /// Obtain the default Poseidon parameters for this rate and for this prime field,
    /// with a specific optimization goal.
    fn get_default_poseidon_parameters(
        rate: usize,
        optimized_for_weights: bool,
    ) -> Option<PoseidonParameters<Self>>;
}
