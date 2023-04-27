use anyhow::Result;
use ark_ff::{vec::Vec, PrimeField};

use crate::{
    matrix::{Matrix, SquareMatrix},
    matrix_ops::{MatrixOperations, SquareMatrixOperations},
};

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

impl<F: PrimeField> MdsMatrix<F> {
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

    /// Compute inverse of MDS matrix
    pub fn inverse(&self) -> SquareMatrix<F> {
        self.0
            .inverse()
            .expect("all well-formed MDS matrices should have inverses")
    }

    /// Return the elements M_{0,1} .. M_{0,t} from the first row
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn v(&self) -> Matrix<F> {
        let elements: Vec<F> = self.0 .0.elements()[1..self.0 .0.n_rows()].to_vec();
        Matrix::new(1, self.0.n_rows() - 1, elements)
    }

    /// Return the elements M_{1,0} .. M_{t,0}from the first column
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn w(&self) -> Matrix<F> {
        let mut elements = Vec::with_capacity(self.0.n_rows() - 1);
        for i in 1..self.n_rows() {
            elements.push(self.get_element(i, 0))
        }
        Matrix::new(&self.n_rows() - 1, 1, elements)
    }

    /// Compute the (t - 1) x (t - 1) Mhat matrix from the MDS matrix
    ///
    /// This is simply the MDS matrix with the first row and column removed
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn hat(&self) -> SquareMatrix<F> {
        let dim = self.n_rows();
        let mut mhat_elements = Vec::with_capacity((dim - 1) * (dim - 1));
        for i in 1..dim {
            for j in 1..dim {
                mhat_elements.push(self.0.get_element(i, j))
            }
        }

        SquareMatrix::from_vec(mhat_elements)
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
