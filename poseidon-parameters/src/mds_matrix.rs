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
>(pub SquareMatrix<STATE_SIZE, NUM_ELEMENTS>);

impl<const STATE_SIZE: usize, const STATE_SIZE_MINUS_1: usize, const NUM_ELEMENTS: usize>
    MatrixOperations for MdsMatrix<STATE_SIZE, STATE_SIZE_MINUS_1, NUM_ELEMENTS>
{
    fn new(elements: &[Fq]) -> Self {
        assert!(STATE_SIZE == STATE_SIZE_MINUS_1 + 1);
        assert!(STATE_SIZE * STATE_SIZE == NUM_ELEMENTS);
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

    fn rows(&self) -> &[&[Fq]] {
        self.0.rows()
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

impl<const STATE_SIZE: usize, const STATE_SIZE_MINUS_1: usize, const NUM_ELEMENTS: usize>
    MdsMatrix<STATE_SIZE, STATE_SIZE_MINUS_1, NUM_ELEMENTS>
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

    fn transpose(&self) -> Self {
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
        Matrix::new(&elements)
    }

    /// Return the elements M_{1,0} .. M_{t,0} from the first column
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn w(&self) -> Matrix<STATE_SIZE_MINUS_1, 1, STATE_SIZE_MINUS_1> {
        let mut elements = [Fq::zero(); STATE_SIZE_MINUS_1];
        for i in 1..self.n_rows() {
            elements[i - 1] = self.get_element(i, 0);
        }
        Matrix::new(&elements)
    }

    // /// Compute the (t - 1) x (t - 1) Mhat matrix from the MDS matrix
    // ///
    // /// This is simply the MDS matrix with the first row and column removed
    // ///
    // /// Ref: p.20 of the Poseidon paper
    // TODO: Need a const generic parameter for the number of elements in the Mhat matrix
    // pub fn hat(&self) -> SquareMatrix<STATE_SIZE_MINUS_1, STATE_SIZE_MINUS_1> {
    //     let dim = self.n_rows();
    //     let mut mhat_elements = Vec::<Fq, MAX_DIMENSION>::new();
    //     for i in 1..dim {
    //         for j in 1..dim {
    //             mhat_elements
    //                 .push(self.0.get_element(i, j))
    //                 .expect("capacity should not be exceeded");
    //         }
    //     }
    //     SquareMatrix::from_vec(mhat_elements)
    // }
}

// impl From<MdsMatrix> for Vec<Vec<Fq, MAX_DIMENSION>, MAX_DIMENSION> {
//     fn from(val: MdsMatrix) -> Self {
//         let mut rows = Vec::<Vec<Fq, MAX_DIMENSION>, MAX_DIMENSION>::new();
//         for i in 0..val.0.n_rows() {
//             let mut row = Vec::new();
//             for j in 0..val.0.n_rows() {
//                 row.push(val.0 .0.get_element(i, j))
//                     .expect("capacity should not be exceeded");
//             }
//             rows.push(row).expect("capacity should not be exceeded");
//         }
//         rows
//     }
// }

/// Represents an optimized MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedMdsMatrices<
    const N_ROUNDS: usize,
    const STATE_SIZE: usize,
    const STATE_SIZE_MINUS_1: usize,
    const NUM_ELEMENTS_STATE_SIZE_2: usize,
    const NUM_ELEMENTS_STATE_SIZE_MINUS_1_2: usize,
> {
    /// A (t - 1) x (t - 1) MDS submatrix derived from the MDS matrix.
    pub M_hat: SquareMatrix<STATE_SIZE_MINUS_1, NUM_ELEMENTS_STATE_SIZE_MINUS_1_2>,
    /// A 1 x (t - 1) (row) vector derived from the MDS matrix.
    pub v: Matrix<1, STATE_SIZE_MINUS_1, STATE_SIZE_MINUS_1>,
    /// A (t - 1) x 1 (column) vector derived from the MDS matrix.
    pub w: Matrix<STATE_SIZE_MINUS_1, 1, STATE_SIZE_MINUS_1>,
    /// A matrix formed from Mhat (an MDS submatrix of the MDS matrix).
    pub M_prime: SquareMatrix<STATE_SIZE, NUM_ELEMENTS_STATE_SIZE_2>,
    /// A sparse matrix formed from M,
    pub M_doubleprime: SquareMatrix<STATE_SIZE, NUM_ELEMENTS_STATE_SIZE_2>,
    /// The inverse of the t x t MDS matrix (needed to compute round constants).
    pub M_inverse: SquareMatrix<STATE_SIZE, NUM_ELEMENTS_STATE_SIZE_2>,
    /// The inverse of the (t - 1) x (t - 1) Mhat matrix.
    pub M_hat_inverse: SquareMatrix<STATE_SIZE_MINUS_1, NUM_ELEMENTS_STATE_SIZE_MINUS_1_2>,
    /// Element at M00
    pub M_00: Fq,
    /// M_i
    pub M_i: Matrix<STATE_SIZE, STATE_SIZE, NUM_ELEMENTS_STATE_SIZE_2>,
    /// v_collection: one per round.
    pub v_collection: [Matrix<1, STATE_SIZE_MINUS_1, STATE_SIZE_MINUS_1>; N_ROUNDS],
    /// w_hat_collection: one per round
    pub w_hat_collection: [Matrix<STATE_SIZE_MINUS_1, 1, STATE_SIZE_MINUS_1>; N_ROUNDS],
}
