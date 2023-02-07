use core::slice::Chunks;

use anyhow::{anyhow, Result};
use ark_ff::{vec::Vec, PrimeField};

pub trait MatrixOperations<F> {
    /// Create a new matrix
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self;
    /// Access elements as a vector
    fn elements(&self) -> &Vec<F>;
    /// Get element[i,j]
    fn get_element(&self, i: usize, j: usize) -> F;
    /// Set element[i,j]
    fn set_element(&mut self, i: usize, j: usize, val: F);
    /// Get rows
    fn rows(&self) -> Vec<&[F]>;
    /// Get rows in chunks
    fn iter_rows(&self) -> Chunks<F> {
        self.elements().chunks(self.n_cols())
    }
    /// Number of rows
    fn n_rows(&self) -> usize;
    /// Number of columns
    fn n_cols(&self) -> usize;
    /// Take transpose of the matrix    
    fn transpose(&self) -> Self;
    /// Compute Hadamard (element-wise) product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized;
}

/// Multiply two matrices
pub fn mat_mul<F: PrimeField, M: MatrixOperations<F>>(lhs: &M, rhs: &M) -> Result<M> {
    if lhs.n_cols() != rhs.n_rows() {
        return Err(anyhow!(
            "matrix dimensions do not allow matrix multiplication"
        ));
    }

    let rhs_T = rhs.transpose();

    Ok(M::new(
        lhs.n_rows(),
        rhs.n_cols(),
        lhs.iter_rows()
            .flat_map(|row| {
                // Rows of the transposed matrix are the columns of the original matrix
                rhs_T
                    .iter_rows()
                    .map(|column| dot_product(row, column))
                    .collect::<Vec<F>>()
            })
            .collect(),
    ))
}

/// Compute vector dot product
pub fn dot_product<F: PrimeField>(a: &[F], b: &[F]) -> F {
    if a.len() != b.len() {
        panic!("vecs not same len")
    }

    a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
}

/// Matrix operations that are defined on square matrices.
pub trait SquareMatrixOperations<F> {
    /// Compute the matrix inverse, if it exists
    fn inverse(&self) -> Result<Self>
    where
        Self: Sized;
    /// Construct an n x n identity matrix
    fn identity(n: usize) -> Self;
    /// Compute the matrix of minors
    fn minors(&self) -> Self;
    /// Compute the matrix of cofactors
    fn cofactors(&self) -> Self;
    /// Compute the matrix determinant
    fn determinant(&self) -> F;
}
