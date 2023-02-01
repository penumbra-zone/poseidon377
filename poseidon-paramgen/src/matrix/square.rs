use core::ops::Deref;

use anyhow::{anyhow, Result};
use ark_ff::PrimeField;
use ark_std::{vec, vec::Vec};
use poseidon_parameters::{Matrix, MatrixOperations, SquareMatrix};

use crate::{mat_mul, SquareMatrixOperations};

use super::MatrixWrapper;

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

impl<F: PrimeField> SquareMatrixOperations<F> for SquareMatrixWrapper<F> {
    /// Construct a dim x dim identity matrix
    fn identity(dim: usize) -> Self {
        let mut m = SquareMatrix::from_vec(vec![F::zero(); dim * dim]);

        // Set diagonals to 1
        for i in 0..dim {
            m.set_element(i, i, F::one());
        }

        m.into()
    }

    /// Compute the inverse of the matrix
    fn inverse(&self) -> Result<Self> {
        let identity = Self::identity(self.n_rows());

        if self.n_rows() == 1 {
            return Ok(Self(SquareMatrix::from_vec(vec![self
                .get_element(0, 0)
                .inverse()
                .expect("inverse of single element must exist for 1x1 matrix")])));
        }

        let determinant = self.determinant();
        if determinant == F::zero() {
            return Err(anyhow!("err: matrix has no inverse"));
        }

        let minors = self.minors();
        let cofactor_matrix = self.cofactors();
        let signed_minors = minors
            .hadamard_product(&cofactor_matrix)
            .expect("minor and cofactor matrix have correct dimensions");
        let adj = signed_minors.transpose();
        let matrix_inverse = Self(adj) * (F::one() / determinant);

        debug_assert_eq!(
            mat_mul(self as &SquareMatrix<F>, &matrix_inverse)
                .expect("matrix and its inverse should have same dimensions"),
            identity.0
        );
        Ok(Self(matrix_inverse))
    }

    /// Compute the (unsigned) minors of this matrix
    fn minors(&self) -> Self {
        match self.n_cols() {
            0 => panic!("matrix has no elements!"),
            1 => SquareMatrix::from_vec(vec![self.get_element(0, 0)]).into(),
            2 => {
                let a = self.get_element(0, 0);
                let b = self.get_element(0, 1);
                let c = self.get_element(1, 0);
                let d = self.get_element(1, 1);
                SquareMatrix::from_vec(vec![d, c, b, a]).into()
            }
            _ => {
                let dim = self.n_rows();
                let mut minor_matrix_elements = Vec::with_capacity(dim * dim);
                for i in 0..dim {
                    for j in 0..dim {
                        let mut elements: Vec<F> = Vec::new();
                        for k in 0..i {
                            for l in 0..j {
                                elements.push(self.get_element(k, l))
                            }
                            for l in (j + 1)..dim {
                                elements.push(self.get_element(k, l))
                            }
                        }
                        for k in i + 1..dim {
                            for l in 0..j {
                                elements.push(self.get_element(k, l))
                            }
                            for l in (j + 1)..dim {
                                elements.push(self.get_element(k, l))
                            }
                        }
                        let minor = SquareMatrixWrapper(SquareMatrix::from_vec(elements));
                        minor_matrix_elements.push(minor.determinant());
                    }
                }
                SquareMatrix::from_vec(minor_matrix_elements).into()
            }
        }
    }

    /// Compute the cofactor matrix, i.e. $C_{ij} = (-1)^{i+j}$
    fn cofactors(&self) -> Self {
        let dim = self.n_rows();
        let mut elements = Vec::with_capacity(dim);
        for i in 0..dim {
            for j in 0..dim {
                elements.push((-F::one()).pow([(i + j) as u64]))
            }
        }
        SquareMatrix::from_vec(elements).into()
    }

    /// Compute the matrix determinant
    fn determinant(&self) -> F {
        match self.n_cols() {
            0 => panic!("matrix has no elements!"),
            1 => self.get_element(0, 0),
            2 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                a11 * a22 - a21 * a12
            }
            3 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);

                a11 * (Self(SquareMatrixWrapper::new_2x2(a22, a23, a32, a33)).determinant())
                    - a12 * (Self(SquareMatrixWrapper::new_2x2(a21, a23, a31, a33)).determinant())
                    + a13 * (Self(SquareMatrixWrapper::new_2x2(a21, a22, a31, a32)).determinant())
            }
            _ => {
                // Unoptimized, but MDS matrices are fairly small, so we do the naive thing
                let mut det = F::zero();
                let mut levi_civita = true;
                let dim = self.n_rows();

                for i in 0..dim {
                    let mut elements: Vec<F> = Vec::new();
                    for k in 0..i {
                        for l in 1..dim {
                            elements.push(self.get_element(k, l))
                        }
                    }
                    for k in i + 1..dim {
                        for l in 1..dim {
                            elements.push(self.get_element(k, l))
                        }
                    }
                    let minor = SquareMatrix::from_vec(elements);
                    if levi_civita {
                        det += self.get_element(i, 0) * SquareMatrixWrapper(minor).determinant();
                    } else {
                        det -= self.get_element(i, 0) * SquareMatrixWrapper(minor).determinant();
                    }
                    levi_civita = !levi_civita;
                }

                det
            }
        }
    }
}

impl<F: PrimeField> SquareMatrixWrapper<F> {
    /// Get row vector at a specified row index.
    pub fn row_vector(&self, i: usize) -> Matrix<F> {
        MatrixWrapper(self.0 .0.clone()).row_vector(i)
    }

    /// Create a 2x2 `SquareMatrix` from four elements.
    pub fn new_2x2(a: F, b: F, c: F, d: F) -> SquareMatrix<F> {
        SquareMatrix::from_vec(vec![a, b, c, d])
    }
}
