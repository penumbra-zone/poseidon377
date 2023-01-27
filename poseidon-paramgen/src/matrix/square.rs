use anyhow::{anyhow, Result};
use ark_ff::PrimeField;
use ark_std::{vec, vec::Vec};
use poseidon_parameters::BasicMatrixOperations;

use crate::{mat_mul, Matrix, MatrixOperations, SquareMatrixOperations};

/// Represents a square matrix over `PrimeField` elements
pub struct SquareMatrix<T>(pub T);

impl<F: PrimeField> MatrixOperations<F> for poseidon_parameters::SquareMatrix<F> {
    /// Take transpose of the matrix
    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    /// Hadamard (element-wise) matrix product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self> {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F: PrimeField> SquareMatrixOperations<F> for poseidon_parameters::SquareMatrix<F> {
    /// Construct a dim x dim identity matrix
    fn identity(dim: usize) -> Self {
        let mut m = poseidon_parameters::SquareMatrix::from_vec(vec![F::zero(); dim * dim]);

        // Set diagonals to 1
        for i in 0..dim {
            m.set_element(i, i, F::one());
        }

        m
    }

    /// Compute the inverse of the matrix
    fn inverse(&self) -> Result<Self> {
        let identity = Self::identity(self.n_rows());

        if self.n_rows() == 1 {
            return Ok(poseidon_parameters::SquareMatrix::from_vec(vec![self
                .get_element(0, 0)
                .inverse()
                .expect("inverse of single element must exist for 1x1 matrix")]));
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
        let matrix_inverse = SquareMatrix(adj) * (F::one() / determinant);

        debug_assert_eq!(
            mat_mul(self, &matrix_inverse)
                .expect("matrix and its inverse should have same dimensions"),
            identity
        );
        Ok(matrix_inverse)
    }

    /// Compute the (unsigned) minors of this matrix
    fn minors(&self) -> Self {
        match self.0.n_cols {
            0 => panic!("matrix has no elements!"),
            1 => poseidon_parameters::SquareMatrix::from_vec(vec![self.get_element(0, 0)]),
            2 => {
                let a = self.get_element(0, 0);
                let b = self.get_element(0, 1);
                let c = self.get_element(1, 0);
                let d = self.get_element(1, 1);
                poseidon_parameters::SquareMatrix::from_vec(vec![d, c, b, a])
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
                        let minor = poseidon_parameters::SquareMatrix::from_vec(elements);
                        minor_matrix_elements.push(minor.determinant());
                    }
                }
                poseidon_parameters::SquareMatrix::from_vec(minor_matrix_elements)
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
        poseidon_parameters::SquareMatrix::from_vec(elements)
    }

    /// Compute the matrix determinant
    fn determinant(&self) -> F {
        match self.0.n_cols {
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

                a11 * (SquareMatrix::new_2x2(a22, a23, a32, a33).determinant())
                    - a12 * (SquareMatrix::new_2x2(a21, a23, a31, a33).determinant())
                    + a13 * (SquareMatrix::new_2x2(a21, a22, a31, a32).determinant())
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
                    let minor = poseidon_parameters::SquareMatrix::from_vec(elements);
                    if levi_civita {
                        det += self.get_element(i, 0) * minor.determinant();
                    } else {
                        det -= self.get_element(i, 0) * minor.determinant();
                    }
                    levi_civita = !levi_civita;
                }

                det
            }
        }
    }
}

impl<F: PrimeField> SquareMatrix<poseidon_parameters::SquareMatrix<F>> {
    /// Get row vector at a specified row index.
    pub fn row_vector(&self, i: usize) -> poseidon_parameters::Matrix<F> {
        let m = Matrix(&self.0 .0);
        m.row_vector(i)
    }

    /// Create a 2x2 `SquareMatrix` from four elements.
    pub fn new_2x2(a: F, b: F, c: F, d: F) -> poseidon_parameters::SquareMatrix<F> {
        poseidon_parameters::SquareMatrix::from_vec(vec![a, b, c, d])
    }
}
