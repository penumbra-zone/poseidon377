use ark_ff::BigInteger;

use crate::log2;

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Clone, Debug)]
pub struct InputParameters<T: BigInteger> {
    /// Whether or not to allow inverse alpha.
    pub allow_inverse: bool,

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
    /// Create a new set of input parameters for a new Poseidon instance.
    pub fn generate(M: usize, t: usize, p: T, allow_inverse: bool) -> InputParameters<T> {
        let log_2_p = log2(p);
        InputParameters {
            M,
            t,
            p,
            log_2_p,
            allow_inverse,
        }
    }
}
