use ark_ff::PrimeField;

use crate::{alpha, input::InputParameters, round_constants, rounds};
use poseidon_parameters::v2::PoseidonParameters;

/// For generating parameters at build time.
pub mod poseidon_build {
    pub use crate::poseidon_build::v2_compile as compile;
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
    let rounds = rounds::v2_generate(&input, &alpha);
    let arc = round_constants::v2_generate(&input, rounds, alpha);
    let m_e = todo!();
    let m_i = todo!();

    PoseidonParameters::<F> {
        M: input.M,
        t: input.t,
        alpha,
        rounds,
        arc,
        m_e,
        m_i,
    }
}
