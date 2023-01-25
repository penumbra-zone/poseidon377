use ark_ff::BigInteger;

use crate::log2;

/// Input parameters that are used to generate Poseidon parameters.
pub struct InputParameters<T>(pub T);

impl<T: BigInteger> InputParameters<poseidon_parameters::InputParameters<T>> {
    /// Create a new set of input parameters for a new Poseidon instance.
    pub fn new(
        M: usize,
        t: usize,
        p: T,
        allow_inverse: bool,
    ) -> poseidon_parameters::InputParameters<T> {
        let log_2_p = log2(p);
        poseidon_parameters::InputParameters {
            M,
            t,
            p,
            log_2_p,
            allow_inverse,
        }
    }
}
