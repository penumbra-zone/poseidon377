use crate::{
    error::PoseidonParameterError,
    matrix::{Matrix, SquareMatrix},
    matrix_ops::{MatrixOperations, SquareMatrixOperations},
    MAX_DIMENSION,
};
use decaf377::Fq;
use heapless::Vec;

/// Represents an MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MdsMatrix(pub SquareMatrix);

impl MatrixOperations for MdsMatrix {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<Fq, MAX_DIMENSION>) -> Self {
        Self(SquareMatrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<Fq, MAX_DIMENSION> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[Fq], MAX_DIMENSION> {
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

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl MdsMatrix {
    /// Instantiate an MDS matrix from a list of elements.
    ///
    /// # Security
    ///
    /// You must ensure this matrix was generated securely,
    /// using the Cauchy method in `fixed_cauchy_matrix` or
    /// using the random subsampling method described in the original
    /// paper.
    pub fn from_elements(elements: Vec<Fq, MAX_DIMENSION>) -> Self {
        Self(SquareMatrix::from_vec(elements))
    }

    /// Compute inverse of MDS matrix
    pub fn inverse(&self) -> SquareMatrix {
        self.0
            .inverse()
            .expect("all well-formed MDS matrices should have inverses")
    }

    /// Return the elements M_{0,1} .. M_{0,t} from the first row
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn v(&self) -> Matrix {
        let mut elements = Vec::<Fq, MAX_DIMENSION>::new();
        elements
            .extend_from_slice(&self.0 .0.elements()[1..self.0 .0.n_rows()])
            .expect("capacity should not be exceeded");
        Matrix::new(1, self.0.n_rows() - 1, elements)
    }

    /// Return the elements M_{1,0} .. M_{t,0}from the first column
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn w(&self) -> Matrix {
        let mut elements = Vec::<Fq, MAX_DIMENSION>::new();
        for i in 1..self.n_rows() {
            elements
                .push(self.get_element(i, 0))
                .expect("capacity should not be exceeded");
        }
        Matrix::new(&self.n_rows() - 1, 1, elements)
    }

    /// Compute the (t - 1) x (t - 1) Mhat matrix from the MDS matrix
    ///
    /// This is simply the MDS matrix with the first row and column removed
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn hat(&self) -> SquareMatrix {
        let dim = self.n_rows();
        let mut mhat_elements = Vec::<Fq, MAX_DIMENSION>::new();
        for i in 1..dim {
            for j in 1..dim {
                mhat_elements
                    .push(self.0.get_element(i, j))
                    .expect("capacity should not be exceeded");
            }
        }

        SquareMatrix::from_vec(mhat_elements)
    }
}

impl From<MdsMatrix> for Vec<Vec<Fq, MAX_DIMENSION>, MAX_DIMENSION> {
    fn from(val: MdsMatrix) -> Self {
        let mut rows = Vec::<Vec<Fq, MAX_DIMENSION>, MAX_DIMENSION>::new();
        for i in 0..val.0.n_rows() {
            let mut row = Vec::new();
            for j in 0..val.0.n_rows() {
                row.push(val.0 .0.get_element(i, j))
                    .expect("capacity should not be exceeded");
            }
            rows.push(row).expect("capacity should not be exceeded");
        }
        rows
    }
}

/// Represents an optimized MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedMdsMatrices {
    /// A (t - 1) x (t - 1) MDS submatrix derived from the MDS matrix.
    pub M_hat: SquareMatrix,
    /// A 1 x (t - 1) (row) vector derived from the MDS matrix.
    pub v: Matrix,
    /// A (t - 1) x 1 (column) vector derived from the MDS matrix.
    pub w: Matrix,
    /// A matrix formed from Mhat (an MDS submatrix of the MDS matrix).
    pub M_prime: SquareMatrix,
    /// A sparse matrix formed from M,
    pub M_doubleprime: SquareMatrix,
    /// The inverse of the t x t MDS matrix (needed to compute round constants).
    pub M_inverse: SquareMatrix,
    /// The inverse of the (t - 1) x (t - 1) Mhat matrix.
    pub M_hat_inverse: SquareMatrix,
    /// Element at M00
    pub M_00: Fq,
    /// M_i
    pub M_i: Matrix,
    /// v_collection: one per round.
    pub v_collection: Vec<Matrix, MAX_DIMENSION>,
    /// w_hat_collection: one per round
    pub w_hat_collection: Vec<Matrix, MAX_DIMENSION>,
}
