use anyhow::{anyhow, Result};
use ark_ff::PrimeField;

use crate::MatrixOperations;

/// Represents a matrix over `PrimeField` elements.
///
/// This matrix can be used to represent row or column
/// vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F: PrimeField> {
    /// Elements of the matrix.
    pub elements: Vec<F>,
    /// Number of columns.
    pub n_cols: usize,
    /// Number of rows.
    pub n_rows: usize,
}

impl<F: PrimeField> MatrixOperations<F> for Matrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Matrix<F> {
        if elements.len() != n_rows * n_cols {
            panic!("Matrix has insufficient elements")
        }
        Matrix {
            elements,
            n_cols,
            n_rows,
        }
    }

    fn n_rows(&self) -> usize {
        self.n_rows
    }

    fn n_cols(&self) -> usize {
        self.n_cols
    }

    fn elements(&self) -> &Vec<F> {
        &self.elements
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.elements[i * self.n_cols + j]
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.elements[i * self.n_cols + j] = val
    }

    fn rows(&self) -> Vec<&[F]> {
        self.elements.chunks(self.n_cols).collect()
    }

    /// Take transpose of the matrix
    fn transpose(&self) -> Matrix<F> {
        let mut transposed_elements = Vec::with_capacity(self.n_rows * self.n_cols);

        for j in 0..self.n_cols {
            for i in 0..self.n_rows {
                transposed_elements.push(self.get_element(i, j))
            }
        }
        Matrix::new(self.n_cols, self.n_rows, transposed_elements)
    }

    /// Hadamard (element-wise) matrix product
    fn hadamard_product(&self, rhs: &Matrix<F>) -> Result<Matrix<F>> {
        if self.n_rows != rhs.n_rows || self.n_cols != rhs.n_cols {
            return Err(anyhow!("Hadamard product requires same shape matrices"));
        }

        let mut new_elements = Vec::with_capacity(self.n_rows * self.n_cols);
        for i in 0..self.n_rows {
            for j in 0..self.n_cols {
                new_elements.push(self.get_element(i, j) * rhs.get_element(i, j));
            }
        }

        Ok(Matrix::new(self.n_rows, self.n_cols, new_elements))
    }
}

impl<F: PrimeField> Matrix<F> {
    /// Get row vector at a specified row index
    pub fn row_vector(&self, i: usize) -> Matrix<F> {
        let mut row_elements = Vec::with_capacity(self.n_cols);
        for j in 0..self.n_cols {
            row_elements.push(self.get_element(i, j));
        }
        Matrix::new(1, self.n_cols, row_elements)
    }
}
