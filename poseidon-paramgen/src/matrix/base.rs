use core::ops::Deref;

use ark_ff::PrimeField;
use poseidon_parameters::Matrix;

/// Represents a matrix over `PrimeField` elements.
#[derive(Debug)]
pub struct MatrixWrapper<F: PrimeField>(pub Matrix<F>);

impl<F: PrimeField> From<Matrix<F>> for MatrixWrapper<F> {
    fn from(value: Matrix<F>) -> Self {
        Self(value)
    }
}

impl<F: PrimeField> Deref for MatrixWrapper<F> {
    type Target = Matrix<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
