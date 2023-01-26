use core::slice::Chunks;

use anyhow::Result;
use ark_std::vec::Vec;

// TODO : remove basic operations covered by the basic trait
/// Basic matrix operations all matrices should implement.
pub trait MatrixOperations<F> {
    // /// Create a new matrix
    // fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self;
    // /// Access elements as a vector
    // fn elements(&self) -> &Vec<F>;
    // /// Get element[i,j]
    // fn get_element(&self, i: usize, j: usize) -> F;
    // /// Set element[i,j]
    // fn set_element(&mut self, i: usize, j: usize, val: F);
    // /// Get rows
    // fn rows(&self) -> Vec<&[F]>;
    // /// Get rows in chunks
    // fn iter_rows(&self) -> Chunks<F> {
    //     self.elements().chunks(self.n_cols())
    // }
    /// Compute matrix transpose
    fn transpose(&self) -> Self;
    /// Compute Hadamard (element-wise) product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized;
    // /// Number of rows
    // fn n_rows(&self) -> usize;
    // /// Number of columns
    // fn n_cols(&self) -> usize;
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
