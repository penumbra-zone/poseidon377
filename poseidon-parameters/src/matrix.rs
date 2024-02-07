use core::ops::Mul;

// TODO: This library needs to be replaced
use num_integer::Roots;

use crate::error::PoseidonParameterError;
use crate::matrix_ops::{mat_mul, MatrixOperations, SquareMatrixOperations};
use decaf377::Fq;

/// Represents a matrix over `PrimeField` elements.
///
/// This matrix can be used to represent row or column
/// vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix {
    /// Elements of the matrix.
    pub elements: Vec<Fq>,
    /// Number of columns.
    pub n_cols: usize,
    /// Number of rows.
    pub n_rows: usize,
}

impl MatrixOperations for Matrix {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<Fq>) -> Self {
        if elements.len() != n_rows * n_cols {
            panic!("Matrix has an insufficient number of elements")
        }
        Self {
            elements,
            n_cols,
            n_rows,
        }
    }

    fn elements(&self) -> &Vec<Fq> {
        &self.elements
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.elements[i * self.n_cols + j]
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.elements[i * self.n_cols + j] = val
    }

    fn rows(&self) -> Vec<&[Fq]> {
        self.elements.chunks(self.n_cols).collect()
    }

    fn n_rows(&self) -> usize {
        self.n_rows
    }

    fn n_cols(&self) -> usize {
        self.n_cols
    }

    fn transpose(&self) -> Self {
        let mut transposed_elements = Vec::with_capacity(self.n_rows * self.n_cols);

        for j in 0..self.n_cols {
            for i in 0..self.n_rows {
                transposed_elements.push(self.get_element(i, j))
            }
        }
        Self::new(self.n_cols, self.n_rows, transposed_elements)
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        if self.n_rows != rhs.n_rows || self.n_cols != rhs.n_cols {
            return Err(PoseidonParameterError::InvalidMatrixDimensions);
        }

        let mut new_elements = Vec::with_capacity(self.n_rows * self.n_cols);
        for i in 0..self.n_rows {
            for j in 0..self.n_cols {
                new_elements.push(self.get_element(i, j) * rhs.get_element(i, j));
            }
        }

        Ok(Self::new(self.n_rows, self.n_cols, new_elements))
    }
}

/// Multiply scalar by Matrix
impl Mul<Fq> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Fq) -> Self::Output {
        let elements = self.elements();
        let new_elements: Vec<Fq> = elements.iter().map(|element| *element * rhs).collect();
        Self::new(self.n_rows(), self.n_cols(), new_elements)
    }
}

impl Matrix {
    /// Get row vector at a specified row index
    pub fn row_vector(&self, i: usize) -> Matrix {
        let mut row_elements = Vec::with_capacity(self.n_cols);
        for j in 0..self.n_cols {
            row_elements.push(self.get_element(i, j));
        }
        Matrix::new(1, self.n_cols, row_elements)
    }
}

/// Represents a square matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquareMatrix(pub Matrix);

impl MatrixOperations for SquareMatrix {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<Fq>) -> Self {
        Self(Matrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<Fq> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[Fq]> {
        self.0.rows()
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl SquareMatrixOperations for SquareMatrix {
    /// Compute the inverse of the matrix
    fn inverse(&self) -> Result<Self, PoseidonParameterError> {
        let identity = Self::identity(self.n_rows());

        if self.n_rows() == 1 {
            return Ok(Self::from_vec(vec![self
                .get_element(0, 0)
                .inverse()
                .expect("inverse of single element must exist for 1x1 matrix")]));
        }

        let determinant = self.determinant();
        if determinant == Fq::zero() {
            return Err(PoseidonParameterError::NoMatrixInverse);
        }

        let minors = self.minors();
        let cofactor_matrix = self.cofactors();
        let signed_minors = minors
            .hadamard_product(&cofactor_matrix)
            .expect("minor and cofactor matrix have correct dimensions");
        let adj = signed_minors.transpose();
        let matrix_inverse = adj * (Fq::one() / determinant);

        debug_assert_eq!(
            mat_mul(self, &matrix_inverse)
                .expect("matrix and its inverse should have same dimensions"),
            identity
        );
        Ok(matrix_inverse)
    }

    /// Construct a dim x dim identity matrix
    fn identity(dim: usize) -> Self {
        let mut m = Self::from_vec(vec![Fq::zero(); dim * dim]);

        // Set diagonals to 1
        for i in 0..dim {
            m.set_element(i, i, Fq::one());
        }

        m
    }

    /// Compute the (unsigned) minors of this matrix
    fn minors(&self) -> Self {
        match self.n_cols() {
            0 => panic!("matrix has no elements!"),
            1 => Self::from_vec(vec![self.get_element(0, 0)]),
            2 => {
                let a = self.get_element(0, 0);
                let b = self.get_element(0, 1);
                let c = self.get_element(1, 0);
                let d = self.get_element(1, 1);
                Self::from_vec(vec![d, c, b, a])
            }
            _ => {
                let dim = self.n_rows();
                let mut minor_matrix_elements = Vec::with_capacity(dim * dim);
                for i in 0..dim {
                    for j in 0..dim {
                        let mut elements: Vec<Fq> = Vec::new();
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
                        let minor = Self::from_vec(elements);
                        minor_matrix_elements.push(minor.determinant());
                    }
                }
                Self::from_vec(minor_matrix_elements)
            }
        }
    }

    /// Compute the cofactor matrix, i.e. $C_{ij} = (-1)^{i+j}$
    fn cofactors(&self) -> Self {
        let dim = self.n_rows();
        let mut elements = Vec::with_capacity(dim);
        // TODO: non arkworks Fq::pow
        use crate::StuffThatNeedsToGoInDecaf377;
        for i in 0..dim {
            for j in 0..dim {
                elements.push((-Fq::one()).pow([(i + j) as u64]))
            }
        }
        Self::from_vec(elements)
    }

    /// Compute the matrix determinant
    fn determinant(&self) -> Fq {
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

                a11 * (Self::new_2x2(a22, a23, a32, a33).determinant())
                    - a12 * (Self::new_2x2(a21, a23, a31, a33).determinant())
                    + a13 * (Self::new_2x2(a21, a22, a31, a32).determinant())
            }
            _ => {
                // Unoptimized, but MDS matrices are fairly small, so we do the naive thing
                let mut det = Fq::zero();
                let mut levi_civita = true;
                let dim = self.n_rows();

                for i in 0..dim {
                    let mut elements: Vec<Fq> = Vec::new();
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
                    let minor = Self::from_vec(elements);
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

/// Multiply scalar by SquareMatrix
impl Mul<Fq> for SquareMatrix {
    type Output = SquareMatrix;

    fn mul(self, rhs: Fq) -> Self::Output {
        let elements = self.0.elements();
        let new_elements: Vec<Fq> = elements.iter().map(|element| *element * rhs).collect();
        Self::from_vec(new_elements)
    }
}

impl SquareMatrix {
    /// Create a `SquareMatrix` from a vector of elements.
    pub fn from_vec(elements: Vec<Fq>) -> Self {
        let dim = elements.len().sqrt();
        if dim * dim != elements.len() {
            panic!("SquareMatrix must be square")
        }
        Self(Matrix::new(dim, dim, elements))
    }

    /// Get row vector at a specified row index.
    pub fn row_vector(&self, i: usize) -> Matrix {
        self.0.row_vector(i)
    }

    /// Create a 2x2 `SquareMatrix` from four elements.
    pub fn new_2x2(a: Fq, b: Fq, c: Fq, d: Fq) -> SquareMatrix {
        Self::from_vec(vec![a, b, c, d])
    }
}
