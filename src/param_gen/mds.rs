use ark_ff::BigInteger;
use ark_ff::PrimeField;

use crate::{InputParameters, Matrix};

/// The number of attempts to find a secure MDS matrix.
const NUM_ATTEMPTS: usize = 100;

/// Represents a secure MDS (maximum distance separable) matrix.
pub(super) struct MdsMatrix<F: PrimeField>(Matrix<F>);

impl<F> MdsMatrix<F>
where
    F: PrimeField,
{
    pub fn new<P: BigInteger>(input: &InputParameters<P>) -> Self {
        todo!()
    }
}
