pub use crate::alpha::Alpha;
pub use crate::arc_matrix::ArcMatrix;
pub use crate::matrix::SquareMatrix;
pub use crate::round_numbers::RoundNumbers;

pub use crate::{
    matrix_ops::MatrixOperations, matrix_ops::Polynomial, matrix_ops::SquareMatrixOperations,
};

/// A set of Poseidon2 parameters for a given set of input parameters.
#[derive(Clone, Debug)]
pub struct PoseidonParameters {
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

    /// External matrix
    pub m_e: SquareMatrix,

    /// Internal matrix
    pub m_i: SquareMatrix,

    /// Round constants
    pub arc: ArcMatrix,
}
