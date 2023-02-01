use anyhow::Result;

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
