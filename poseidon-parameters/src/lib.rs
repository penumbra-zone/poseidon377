#![no_std]
#![allow(non_snake_case)]

use core::slice::Chunks;

use ark_ff::BigInteger;
use ark_ff::PrimeField;
use ark_std::vec::Vec;

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

// impl<F: PrimeField> Matrix<F> {
//     pub fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
//         // if elements.len() != n_rows * n_cols {
//         //     panic!("Matrix has insufficient elements")
//         // }
//         Self {
//             elements,
//             n_cols,
//             n_rows,
//         }
//     }
// }

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

pub trait BasicMatrixOperations<F> {
    /// Create a new matrix
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self;
    /// Access elements as a vector
    fn elements(&self) -> &Vec<F>;
    /// Get element[i,j]
    fn get_element(&self, i: usize, j: usize) -> F;
    /// Get rows in chunks
    fn iter_rows(&self) -> Chunks<F> {
        self.elements().chunks(self.n_cols())
    }
    /// Number of rows
    fn n_rows(&self) -> usize;
    /// Number of columns
    fn n_cols(&self) -> usize;
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

impl<F: PrimeField> BasicMatrixOperations<F> for Matrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        todo!()
    }

    fn elements(&self) -> &Vec<F> {
        todo!()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        todo!()
    }

    fn n_rows(&self) -> usize {
        todo!()
    }

    fn n_cols(&self) -> usize {
        todo!()
    }
}

/// Represents a square matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquareMatrix<F: PrimeField>(pub Matrix<F>);

impl<F: PrimeField> BasicMatrixOperations<F> for SquareMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        todo!()
    }

    fn elements(&self) -> &Vec<F> {
        todo!()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        todo!()
    }

    fn n_rows(&self) -> usize {
        todo!()
    }

    fn n_cols(&self) -> usize {
        todo!()
    }
}

/// Represents an MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MdsMatrix<F: PrimeField>(pub SquareMatrix<F>);

impl<F: PrimeField> BasicMatrixOperations<F> for MdsMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        todo!()
    }

    fn elements(&self) -> &Vec<F> {
        todo!()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        todo!()
    }

    fn n_rows(&self) -> usize {
        todo!()
    }

    fn n_cols(&self) -> usize {
        todo!()
    }
}

/// Represents an matrix of round constants.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArcMatrix<F: PrimeField>(pub Matrix<F>);

impl<F: PrimeField> BasicMatrixOperations<F> for ArcMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
        todo!()
    }

    fn elements(&self) -> &Vec<F> {
        todo!()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        todo!()
    }

    fn n_rows(&self) -> usize {
        todo!()
    }

    fn n_cols(&self) -> usize {
        todo!()
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

// pub fn meh<F: PrimeField>(elements: &Vec<F>) -> Chunks<F> {
//     // F::from_str("aa").map_err(|_| ()).unwrap();
//     // let elements: Vec<F> = Vec::new();
//     elements.chunks(2)
//     let elements: Vec<F> = Vec::new();
//     let len = elements.len();
//     let a = len.sqrt();
// }
