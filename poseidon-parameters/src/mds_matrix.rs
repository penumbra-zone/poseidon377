use crate::{
    error::PoseidonParameterError,
    matrix::{Matrix, SquareMatrix},
    matrix_ops::{MatrixOperations, SquareMatrixOperations},
};
use decaf377::Fq;

/// Represents an MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MdsMatrix<
    const STATE_SIZE: usize,
    const STATE_SIZE_MINUS_1: usize,
    const NUM_ELEMENTS: usize,
    const NUM_ELEMENTS_STATE_SIZE_MINUS_1_2: usize,
>(pub SquareMatrix<STATE_SIZE, NUM_ELEMENTS>);

impl<
        const STATE_SIZE: usize,
        const STATE_SIZE_MINUS_1: usize,
        const NUM_ELEMENTS: usize,
        const NUM_ELEMENTS_STATE_SIZE_MINUS_1_2: usize,
    > MatrixOperations
    for MdsMatrix<STATE_SIZE, STATE_SIZE_MINUS_1, NUM_ELEMENTS, NUM_ELEMENTS_STATE_SIZE_MINUS_1_2>
{
    fn new(elements: &[Fq]) -> Self {
        assert!(STATE_SIZE == STATE_SIZE_MINUS_1 + 1);
        assert!(STATE_SIZE * STATE_SIZE == NUM_ELEMENTS);
        assert!(STATE_SIZE_MINUS_1 * STATE_SIZE_MINUS_1 == NUM_ELEMENTS_STATE_SIZE_MINUS_1_2);
        Self(SquareMatrix::new(elements))
    }

    fn elements(&self) -> &[Fq] {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.0.set_element(i, j, val)
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<
        const STATE_SIZE: usize,
        const STATE_SIZE_MINUS_1: usize,
        const NUM_ELEMENTS: usize,
        const NUM_ELEMENTS_STATE_SIZE_MINUS_1_2: usize,
    > MdsMatrix<STATE_SIZE, STATE_SIZE_MINUS_1, NUM_ELEMENTS, NUM_ELEMENTS_STATE_SIZE_MINUS_1_2>
{
    /// Instantiate an MDS matrix from a list of elements.
    ///
    /// # Security
    ///
    /// You must ensure this matrix was generated securely,
    /// using the Cauchy method in `fixed_cauchy_matrix` or
    /// using the random subsampling method described in the original
    /// paper.
    pub fn from_elements(elements: &[Fq]) -> Self {
        Self(SquareMatrix::new(elements))
    }

    pub fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    /// Compute inverse of MDS matrix
    pub fn inverse(&self) -> SquareMatrix<STATE_SIZE, NUM_ELEMENTS> {
        self.0
            .inverse()
            .expect("all well-formed MDS matrices should have inverses")
    }

    /// Return the elements M_{0,1} .. M_{0,t} from the first row
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn v(&self) -> Matrix<1, STATE_SIZE_MINUS_1, STATE_SIZE_MINUS_1> {
        let elements = &self.0 .0.elements()[1..self.0 .0.n_rows()];
        Matrix::new(elements)
    }

    /// Return the elements M_{1,0} .. M_{t,0} from the first column
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn w(&self) -> Matrix<STATE_SIZE_MINUS_1, 1, STATE_SIZE_MINUS_1> {
        let mut elements = [Fq::from(0u64); STATE_SIZE_MINUS_1];
        for i in 1..self.n_rows() {
            elements[i - 1] = self.get_element(i, 0);
        }
        Matrix::new(&elements)
    }

    /// Compute the (t - 1) x (t - 1) Mhat matrix from the MDS matrix
    ///
    /// This is simply the MDS matrix with the first row and column removed
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn hat(&self) -> SquareMatrix<STATE_SIZE_MINUS_1, NUM_ELEMENTS_STATE_SIZE_MINUS_1_2> {
        let dim = self.n_rows();
        let mut mhat_elements = [Fq::from(0u64); NUM_ELEMENTS_STATE_SIZE_MINUS_1_2];
        let mut index = 0;
        for i in 1..dim {
            for j in 1..dim {
                mhat_elements[index] = self.0.get_element(i, j);
                index += 1;
            }
        }
        SquareMatrix::new(&mhat_elements)
    }

    /// Create a new matrix from a slice of elements.
    ///
    /// # Security
    ///
    /// You must ensure this matrix was generated securely,
    /// using the Cauchy method in `fixed_cauchy_matrix` or
    /// using the random subsampling method described in the original
    /// paper.
    pub const fn new_from_known(elements: [Fq; NUM_ELEMENTS]) -> Self {
        Self(SquareMatrix::new_from_known(elements))
    }
}

/// Represents an optimized MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedMdsMatrices<
    const N_ROUNDS: usize,
    const N_PARTIAL_ROUNDS: usize,
    const STATE_SIZE: usize,
    const STATE_SIZE_MINUS_1: usize,
    const NUM_MDS_ELEMENTS: usize,
    const NUM_STATE_SIZE_MINUS_1_ELEMENTS: usize,
> {
    /// A (t - 1) x (t - 1) MDS submatrix derived from the MDS matrix.
    pub M_hat: SquareMatrix<STATE_SIZE_MINUS_1, NUM_STATE_SIZE_MINUS_1_ELEMENTS>,
    /// A 1 x (t - 1) (row) vector derived from the MDS matrix.
    pub v: Matrix<1, STATE_SIZE_MINUS_1, STATE_SIZE_MINUS_1>,
    /// A (t - 1) x 1 (column) vector derived from the MDS matrix.
    pub w: Matrix<STATE_SIZE_MINUS_1, 1, STATE_SIZE_MINUS_1>,
    /// A matrix formed from Mhat (an MDS submatrix of the MDS matrix).
    pub M_prime: SquareMatrix<STATE_SIZE, NUM_MDS_ELEMENTS>,
    /// A sparse matrix formed from M,
    pub M_doubleprime: SquareMatrix<STATE_SIZE, NUM_MDS_ELEMENTS>,
    /// The inverse of the t x t MDS matrix (needed to compute round constants).
    pub M_inverse: SquareMatrix<STATE_SIZE, NUM_MDS_ELEMENTS>,
    /// The inverse of the (t - 1) x (t - 1) Mhat matrix.
    pub M_hat_inverse: SquareMatrix<STATE_SIZE_MINUS_1, NUM_STATE_SIZE_MINUS_1_ELEMENTS>,
    /// Element at M00
    pub M_00: Fq,
    /// M_i
    pub M_i: Matrix<STATE_SIZE, STATE_SIZE, NUM_MDS_ELEMENTS>,
    /// v_collection: one per partial round.
    pub v_collection: [Matrix<1, STATE_SIZE_MINUS_1, STATE_SIZE_MINUS_1>; N_PARTIAL_ROUNDS],
    /// w_hat_collection: one per round
    pub w_hat_collection: [Matrix<STATE_SIZE_MINUS_1, 1, STATE_SIZE_MINUS_1>; N_PARTIAL_ROUNDS],
}
