use core::ops::Deref;

use ark_ff::PrimeField;
use poseidon_parameters::SquareMatrix;

/// Represents a square matrix over `PrimeField` elements
#[derive(Debug)]
pub struct SquareMatrixWrapper<F: PrimeField>(pub SquareMatrix<F>);

impl<F: PrimeField> From<SquareMatrix<F>> for SquareMatrixWrapper<F> {
    fn from(value: SquareMatrix<F>) -> Self {
        Self(value)
    }
}

impl<F: PrimeField> Deref for SquareMatrixWrapper<F> {
    type Target = SquareMatrix<F>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
