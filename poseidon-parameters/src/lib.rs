#![cfg_attr(not(feature = "std"), no_std)]
#![allow(non_snake_case)]

#[cfg(test)]
mod tests;

use anyhow::anyhow;
use anyhow::Result;
use core::slice::Chunks;
use num_integer::Roots;

use ark_ff::BigInteger;
use ark_ff::PrimeField;
use ark_std::{ops::Mul, vec, vec::Vec};

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Clone, Debug)]
pub struct InputParameters<T: BigInteger> {
    /// Whether or not to allow inverse alpha.
    pub allow_inverse: bool,

    /// Security level in bits.
    pub M: usize,

    /// Width of desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    pub t: usize,

    /// Modulus of the prime field.
    pub p: T, // let modulus = <F as PrimeField>::Params::MODULUS;

    // The below are derived values, stored for convenience.
    /// log_2(p)
    pub log_2_p: f64,
}

/// The exponent in `Sbox(x) = x^\alpha`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Alpha {
    /// A positive exponent $x^{alpha}$.
    Exponent(u32),
    /// 1/x
    Inverse,
}

/// `RoundNumbers` required for security based on known attacks.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct RoundNumbers {
    /// Number of partial rounds.
    pub r_P: usize,
    /// Number of full rounds.
    pub r_F: usize,
}

impl RoundNumbers {
    /// Number of full rounds.    
    pub fn full(&self) -> usize {
        self.r_F
    }

    /// Number of partial rounds.    
    pub fn partial(&self) -> usize {
        self.r_P
    }

    /// Number of total rounds.
    pub fn total(&self) -> usize {
        self.r_P + self.r_F
    }
}

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

/// Represents a matrix over `PrimeField` elements.
///
/// This matrix can be used to represent row or column
/// vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F: PrimeField> {
    /// Elements of the matrix.
    pub elements: Vec<F>,
    /// Number of columns.
    pub n_cols: usize,
    /// Number of rows.
    pub n_rows: usize,
}

impl<F: PrimeField> MatrixOperations<F> for Matrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        if elements.len() != n_rows * n_cols {
            panic!("Matrix has an insufficient number of elements")
        }
        Self {
            elements,
            n_cols,
            n_rows,
        }
    }

    fn elements(&self) -> &Vec<F> {
        &self.elements
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.elements[i * self.n_cols + j]
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.elements[i * self.n_cols + j] = val
    }

    fn rows(&self) -> Vec<&[F]> {
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

    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized,
    {
        if self.n_rows != rhs.n_rows || self.n_cols != rhs.n_cols {
            return Err(anyhow!("Hadamard product requires same shape matrices"));
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
impl<F: PrimeField> Mul<F> for Matrix<F> {
    type Output = Matrix<F>;

    fn mul(self, rhs: F) -> Self::Output {
        let elements = self.elements();
        let new_elements: Vec<F> = elements.iter().map(|element| *element * rhs).collect();
        Self::new(self.n_rows(), self.n_cols(), new_elements)
    }
}

impl<F: PrimeField> Matrix<F> {
    /// Get row vector at a specified row index
    pub fn row_vector(&self, i: usize) -> Matrix<F> {
        let mut row_elements = Vec::with_capacity(self.n_cols);
        for j in 0..self.n_cols {
            row_elements.push(self.get_element(i, j));
        }
        Matrix::new(1, self.n_cols, row_elements)
    }
}

/// Represents a square matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquareMatrix<F: PrimeField>(pub Matrix<F>);

impl<F: PrimeField> MatrixOperations<F> for SquareMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        Self(Matrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<F> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[F]> {
        self.0.rows()
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F: PrimeField> SquareMatrixOperations<F> for SquareMatrix<F> {
    /// Construct a dim x dim identity matrix
    fn identity(dim: usize) -> Self {
        let mut m = Self::from_vec(vec![F::zero(); dim * dim]);

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
            return Ok(Self::from_vec(vec![self
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
        let matrix_inverse = adj * (F::one() / determinant);

        debug_assert_eq!(
            mat_mul(self, &matrix_inverse)
                .expect("matrix and its inverse should have same dimensions"),
            identity
        );
        Ok(matrix_inverse)
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
        for i in 0..dim {
            for j in 0..dim {
                elements.push((-F::one()).pow([(i + j) as u64]))
            }
        }
        Self::from_vec(elements)
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

                a11 * (Self::new_2x2(a22, a23, a32, a33).determinant())
                    - a12 * (Self::new_2x2(a21, a23, a31, a33).determinant())
                    + a13 * (Self::new_2x2(a21, a22, a31, a32).determinant())
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
impl<F: PrimeField> Mul<F> for SquareMatrix<F> {
    type Output = SquareMatrix<F>;

    fn mul(self, rhs: F) -> Self::Output {
        let elements = self.0.elements();
        let new_elements: Vec<F> = elements.iter().map(|element| *element * rhs).collect();
        Self::from_vec(new_elements)
    }
}

impl<F: PrimeField> SquareMatrix<F> {
    /// Create a `SquareMatrix` from a vector of elements.
    pub fn from_vec(elements: Vec<F>) -> Self {
        let dim = elements.len().sqrt();
        if dim * dim != elements.len() {
            panic!("SquareMatrix must be square")
        }
        Self(Matrix::new(dim, dim, elements))
    }

    /// Get row vector at a specified row index.
    pub fn row_vector(&self, i: usize) -> Matrix<F> {
        self.0.row_vector(i)
    }

    /// Create a 2x2 `SquareMatrix` from four elements.
    pub fn new_2x2(a: F, b: F, c: F, d: F) -> SquareMatrix<F> {
        Self::from_vec(vec![a, b, c, d])
    }
}

/// Represents an MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MdsMatrix<F: PrimeField>(pub SquareMatrix<F>);

impl<F: PrimeField> MatrixOperations<F> for MdsMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        Self(SquareMatrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<F> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[F]> {
        self.0.rows()
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F> MdsMatrix<F>
where
    F: PrimeField,
{
    /// Instantiate an MDS matrix from a list of elements.
    ///
    /// # Security
    ///
    /// You must ensure this matrix was generated securely,
    /// using the Cauchy method in `fixed_cauchy_matrix` or
    /// using the random subsampling method described in the original
    /// paper.
    pub fn from_elements(elements: Vec<F>) -> Self {
        Self(SquareMatrix::from_vec(elements))
    }
}

impl<F: PrimeField> From<MdsMatrix<F>> for Vec<Vec<F>> {
    fn from(val: MdsMatrix<F>) -> Self {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..val.0.n_rows() {
            let mut row = Vec::new();
            for j in 0..val.0.n_rows() {
                row.push(val.0 .0.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

/// Represents an matrix of round constants.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArcMatrix<F: PrimeField>(pub Matrix<F>);

impl<F: PrimeField> MatrixOperations<F> for ArcMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        Self(Matrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<F> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[F]> {
        self.0.rows()
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F: PrimeField> From<ArcMatrix<F>> for Vec<Vec<F>> {
    fn from(arc: ArcMatrix<F>) -> Self {
        let mut rows = Vec::<Vec<F>>::new();
        let m = &arc.0;

        for i in 0..arc.n_rows() {
            let mut row = Vec::new();
            for j in 0..arc.n_cols() {
                row.push(m.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

/// Represents an optimized matrix of round constants.
///
/// This modifies the partial rounds in the middle of the permutation,
/// wherein you add constants _first_ before iterating through the partial
/// rounds.
///
/// This method follows `calc_equivalent_constants` from Appendix B's
/// `poseidonperm_x3_64_24_optimized.sage`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedArcMatrix<F: PrimeField>(pub ArcMatrix<F>);

impl<F: PrimeField> MatrixOperations<F> for OptimizedArcMatrix<F> {
    /// Create a `OptimizedArcMatrix` from its elements.
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        Self(ArcMatrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<F> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[F]> {
        self.0.rows()
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

/// Represents an optimized MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedMdsMatrices<F: PrimeField> {
    /// A (t - 1) x (t - 1) MDS submatrix derived from the MDS matrix.
    pub M_hat: SquareMatrix<F>,
    /// A 1 x (t - 1) (row) vector derived from the MDS matrix.
    pub v: Matrix<F>,
    /// A (t - 1) x 1 (column) vector derived from the MDS matrix.
    pub w: Matrix<F>,
    /// A matrix formed from Mhat (an MDS submatrix of the MDS matrix).
    pub M_prime: SquareMatrix<F>,
    /// A sparse matrix formed from M,
    pub M_doubleprime: SquareMatrix<F>,
    /// The inverse of the t x t MDS matrix (needed to compute round constants).
    pub M_inverse: SquareMatrix<F>,
    /// The inverse of the (t - 1) x (t - 1) Mhat matrix.
    pub M_hat_inverse: SquareMatrix<F>,
    /// Element at M00
    pub M_00: F,
    /// M_i
    pub M_i: Matrix<F>,
    /// v_collection: one per round.
    pub v_collection: Vec<Matrix<F>>,
    /// w_hat_collection: one per round
    pub w_hat_collection: Vec<Matrix<F>>,
}

// TODO: arc and mds could be collections
/// A set of Poseidon parameters for a given set of input parameters.
#[derive(Clone, Debug)]
pub struct PoseidonParameters<F: PrimeField> {
    // Input parameters.
    /// Security level.
    pub M: usize,
    /// Width of desired hash function, e.g. $t=3$ corresponds to a 2-to-1 hash.
    pub t: usize,

    // Generated parameters.
    /// Exponent of the Sbox, i.e. S-box(x) = x^{\alpha} used in the `SubWords` step
    pub alpha: Alpha,

    /// Round numbers
    pub rounds: RoundNumbers,

    /// `t x t` MDS matrix used in the `MixLayer` step
    pub mds: MdsMatrix<F>,

    /// `num_total_rounds x t` matrix of constants used in the `AddRoundConstant` step
    pub arc: ArcMatrix<F>,

    /// Optimized round constants.
    pub optimized_arc: OptimizedArcMatrix<F>,

    /// Optimized MDS matrices.
    pub optimized_mds: OptimizedMdsMatrices<F>,
}
