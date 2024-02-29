use ark_ff::PrimeField;

use crate::{alpha, input::InputParameters, mds, round_constants, rounds};
use poseidon_parameters::v1::PoseidonParameters;

/// For generating parameters at build time.
pub mod poseidon_build {
    pub use crate::poseidon_build::v1_compile as compile;
}

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
    let input = InputParameters::generate(M, t, p, allow_inverse);
    let alpha = alpha::generate::<F>(p, allow_inverse);
    let rounds = rounds::v1_generate(&input, &alpha);
    let mds = mds::v1_generate(&input);
    let arc = round_constants::v1_generate(&input, rounds, alpha);
    let optimized_mds = mds::generate_optimized(&mds, t, &rounds);
    let optimized_arc = round_constants::generate_optimized(&arc, &mds, &rounds);

    PoseidonParameters::<F> {
        M: input.M,
        alpha,
        rounds,
        mds,
        arc,
        optimized_mds,
        optimized_arc,
    }
}
