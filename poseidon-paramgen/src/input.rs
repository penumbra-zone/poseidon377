use ark_ff::BigInteger;
use poseidon_parameters::InputParameters;

use crate::log2;

/// Create a new set of input parameters for a new Poseidon instance.
pub fn generate<T: BigInteger>(
    M: usize,
    t: usize,
    p: T,
    allow_inverse: bool,
) -> InputParameters<T> {
    let log_2_p = log2(p);
    InputParameters {
        M,
        t,
        p,
        log_2_p,
        allow_inverse,
    }
}
