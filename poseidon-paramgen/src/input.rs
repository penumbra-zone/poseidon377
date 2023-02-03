use core::ops::Deref;

use ark_ff::BigInteger;
use poseidon_parameters::InputParameters;

use crate::log2;

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Debug)]
pub struct InputParametersWrapper<T: BigInteger>(pub InputParameters<T>);

impl<T: BigInteger> From<InputParameters<T>> for InputParametersWrapper<T> {
    fn from(value: InputParameters<T>) -> Self {
        Self(value)
    }
}

impl<T: BigInteger> Deref for InputParametersWrapper<T> {
    type Target = InputParameters<T>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T: BigInteger> InputParametersWrapper<T> {
    /// Create a new set of input parameters for a new Poseidon instance.
    #[allow(clippy::new_ret_no_self)]
    pub fn new(M: usize, t: usize, p: T, allow_inverse: bool) -> InputParameters<T> {
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
