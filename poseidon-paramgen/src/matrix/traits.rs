use anyhow::Result;

/// Basic matrix operations all matrices should implement.
pub trait MatrixOperations<F> {
    /// Access elements as a vector
    fn elements(&self) -> &Vec<F>;
    /// Get element[i,j]
    fn get_element(&self, i: usize, j: usize) -> F;
    /// Set element[i,j]
    fn set_element(&mut self, i: usize, j: usize, val: F);
    /// Get rows
    fn rows(&self) -> Vec<&[F]>;
    /// Compute matrix transpose
    fn transpose(&self) -> Self;
    /// Compute Hadamard (element-wise) product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized;
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
