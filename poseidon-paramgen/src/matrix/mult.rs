use anyhow::{anyhow, Result};
use ark_ff::PrimeField;
use ark_std::{ops::Mul, vec::Vec};

use crate::{Matrix, MatrixOperations, SquareMatrix};

/// Compute vector dot product
pub fn dot_product<F: PrimeField>(a: &[F], b: &[F]) -> F {
    if a.len() != b.len() {
        panic!("vecs not same len")
    }

    a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
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

/// Multiply scalar by Matrix
impl<F: PrimeField> Mul<F> for Matrix<F> {
    type Output = Matrix<F>;

    fn mul(self, rhs: F) -> Self::Output {
        let elements = self.elements();
        let new_elements: Vec<F> = elements.iter().map(|element| *element * rhs).collect();
        Matrix::new(self.n_rows, self.n_cols, new_elements)
    }
}

/// Multiply scalar by SquareMatrix
impl<F: PrimeField> Mul<F> for SquareMatrix<F> {
    type Output = SquareMatrix<F>;

    fn mul(self, rhs: F) -> Self::Output {
        let elements = self.elements();
        let new_elements: Vec<F> = elements.iter().map(|element| *element * rhs).collect();
        SquareMatrix::from_vec(new_elements)
    }
}
