use anyhow::{anyhow, Result};
use ark_ff::PrimeField;
use ark_std::vec::Vec;
use poseidon_parameters::BasicMatrixOperations;

use crate::MatrixOperations;

/// Represents a matrix over `PrimeField` elements.
pub struct Matrix<T>(pub T);

impl<F: PrimeField> MatrixOperations<F> for poseidon_parameters::Matrix<F> {
    /// Take transpose of the matrix
    fn transpose(&self) -> Self {
        let mut transposed_elements = Vec::with_capacity(self.n_rows * self.n_cols);

        for j in 0..self.n_cols {
            for i in 0..self.n_rows {
                transposed_elements.push(self.get_element(i, j))
            }
        }
        Self::new(self.n_cols, self.n_rows, transposed_elements)
    }

    /// Hadamard (element-wise) matrix product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self> {
        if self.n_rows != rhs.n_rows || self.n_cols != rhs.n_cols {
            return Err(anyhow!("Hadamard product requires same shape matrices"));
        }

        let mut new_elements = Vec::with_capacity(self.n_rows * self.n_cols);
        for i in 0..self.n_rows {
            for j in 0..self.n_cols {
                new_elements.push(self.get_element(i, j) * rhs.get_element(i, j));
            }
        }

        Ok(Self::new(self.n_rows, self.n_cols, new_elements))
    }
}

impl<F: PrimeField> Matrix<&poseidon_parameters::Matrix<F>> {
    /// Get row vector at a specified row index
    pub fn row_vector(&self, i: usize) -> poseidon_parameters::Matrix<F> {
        let mut row_elements = Vec::with_capacity(self.0.n_cols);
        for j in 0..self.0.n_cols {
            row_elements.push(self.0.get_element(i, j));
        }
        poseidon_parameters::Matrix::new(1, self.0.n_cols, row_elements)
    }
}
