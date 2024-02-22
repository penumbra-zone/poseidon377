use core::slice::Chunks;

use crate::error::PoseidonParameterError;
use decaf377::Fq;

pub trait MatrixOperations {
    /// Create a new matrix
    fn new(elements: &[Fq]) -> Self;
    /// Access elements as an array of arrays
    fn elements(&self) -> &[Fq];
    /// Get element[i,j]
    fn get_element(&self, i: usize, j: usize) -> Fq;
    /// Set element[i,j]
    fn set_element(&mut self, i: usize, j: usize, val: Fq);
    /// Get rows
    fn rows(&self) -> &[&[Fq]];
    /// Get rows in chunks
    fn iter_rows(&self) -> Chunks<Fq> {
        self.elements().chunks(self.n_cols())
    }
    /// Number of rows
    fn n_rows(&self) -> usize;
    /// Number of columns
    fn n_cols(&self) -> usize;
    /// Take transpose of the matrix
    fn transpose(&self) -> Self;
    /// Compute Hadamard (element-wise) product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized;
}

// TODO: mat_mul

/// Compute vector dot product
pub fn dot_product(a: &[Fq], b: &[Fq]) -> Fq {
    if a.len() != b.len() {
        panic!("vecs not same len")
    }

    a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
}

/// Matrix operations that are defined on square matrices.
pub trait SquareMatrixOperations {
    /// Compute the matrix inverse, if it exists
    fn inverse(&self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized;
    /// Construct an identity matrix
    fn identity() -> Self;
    /// Compute the matrix of minors
    fn minors(&self) -> Self;
    /// Compute the matrix of cofactors
    fn cofactors(&self) -> Self;
    /// Compute the matrix determinant
    fn determinant(&self) -> Fq;
}
