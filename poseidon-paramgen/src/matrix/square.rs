use anyhow::{anyhow, Result};
use ark_ff::PrimeField;

use crate::{Matrix, MatrixOperations};

/// Represents a square matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquareMatrix<F: PrimeField> {
    pub inner: Matrix<F>,
}

impl<F: PrimeField> SquareMatrix<F> {
    pub fn elements(&self) -> &Vec<F> {
        &self.inner.elements
    }

    pub fn rows(&self) -> Vec<&[F]> {
        self.elements().chunks(self.dim()).collect()
    }

    pub fn get_element(&self, i: usize, j: usize) -> F {
        self.inner.get_element(i, j)
    }

    pub fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.inner.set_element(i, j, val)
    }

    pub fn from_vec(elements: Vec<F>) -> Self {
        if (elements.len() as f64).sqrt().fract() != 0.0 {
            panic!("SquareMatrix must be square")
        }

        let dim = (elements.len() as f64).sqrt() as usize;

        SquareMatrix {
            inner: Matrix::new(dim, dim, elements),
        }
    }

    /// Dimension of the dim x dim matrix
    pub fn dim(&self) -> usize {
        self.inner.n_rows
    }

    /// Construct a dim x dim identity matrix
    pub fn identity(dim: usize) -> SquareMatrix<F> {
        let mut m = SquareMatrix::from_vec(vec![F::zero(); dim * dim]);

        // Set diagonals to 1
        for i in 0..dim {
            m.set_element(i, i, F::one());
        }

        m
    }

    /// Take transpose of the matrix
    pub fn transpose(&self) -> SquareMatrix<F> {
        let dim = self.dim();
        let mut transposed_elements = Vec::with_capacity(dim * dim);

        for j in 0..dim {
            for i in 0..dim {
                transposed_elements.push(self.get_element(i, j))
            }
        }
        SquareMatrix::from_vec(transposed_elements)
    }

    pub fn new_2x2(a: F, b: F, c: F, d: F) -> Self {
        SquareMatrix::from_vec(vec![a, b, c, d])
    }

    /// Hadamard (element-wise) matrix product
    pub fn hadamard_product(&self, rhs: &SquareMatrix<F>) -> Result<SquareMatrix<F>> {
        Ok(SquareMatrix {
            inner: self.inner.hadamard_product(&rhs.inner)?,
        })
    }

    /// Compute the inverse of the matrix
    pub fn inverse(&self) -> SquareMatrix<F> {
        let identity: SquareMatrix<F> = SquareMatrix::identity(self.dim());

        if self.dim() == 1 {
            return SquareMatrix::from_vec(vec![self
                .get_element(0, 0)
                .inverse()
                .expect("inverse of single element must exist for 1x1 matrix")]);
        }

        let determinant = self.determinant();
        if determinant == F::zero() {
            panic!("err: matrix has no inverse")
        }

        let minors = self.minors();
        let cofactor_matrix = self.cofactors();
        let signed_minors = minors
            .hadamard_product(&cofactor_matrix)
            .expect("minor and cofactor matrix have correct dimensions");
        let adj = signed_minors.transpose();
        let matrix_inverse = adj * (F::one() / determinant);

        debug_assert_eq!(self * &matrix_inverse, identity);
        matrix_inverse
    }

    /// Compute the (unsigned) minors of this matrix
    pub fn minors(&self) -> SquareMatrix<F> {
        match self.inner.n_cols {
            0 => panic!("matrix has no elements!"),
            1 => SquareMatrix::from_vec(vec![self.get_element(0, 0)]),
            2 => {
                let a = self.get_element(0, 0);
                let b = self.get_element(0, 1);
                let c = self.get_element(1, 0);
                let d = self.get_element(1, 1);
                SquareMatrix::from_vec(vec![d, c, b, a])
            }
            _ => {
                let dim = self.dim();
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
                        let minor = SquareMatrix::from_vec(elements);
                        minor_matrix_elements.push(minor.determinant());
                    }
                }
                SquareMatrix::from_vec(minor_matrix_elements)
            }
        }
    }

    /// Compute the cofactor matrix, i.e. $C_{ij} = (-1)^{i+j}$
    pub fn cofactors(&self) -> SquareMatrix<F> {
        let dim = self.dim();
        let mut elements = Vec::with_capacity(dim);
        for i in 0..dim {
            for j in 0..dim {
                elements.push((-F::one()).pow(&[(i + j) as u64]))
            }
        }
        SquareMatrix::from_vec(elements)
    }

    /// Compute the matrix determinant
    pub fn determinant(&self) -> F {
        match self.inner.n_cols {
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
                let dim = self.dim();

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
