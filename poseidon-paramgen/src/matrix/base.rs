use core::ops::Deref;

use anyhow::{anyhow, Result};
use ark_ff::PrimeField;
use ark_std::vec::Vec;
use poseidon_parameters::Matrix;
use poseidon_parameters::MatrixOperations;

/// Represents a matrix over `PrimeField` elements.
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

impl<F: PrimeField> MatrixWrapper<F> {
    /// Get row vector at a specified row index
    pub fn row_vector(&self, i: usize) -> Matrix<F> {
        let mut row_elements = Vec::with_capacity(self.0.n_cols);
        for j in 0..self.0.n_cols {
            row_elements.push(self.0.get_element(i, j));
        }
        Matrix::new(1, self.0.n_cols, row_elements)
    }
}
