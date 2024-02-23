pub use crate::alpha::Alpha;
pub use crate::arc_matrix::ArcMatrix;
pub use crate::matrix::SquareMatrix;
pub use crate::round_numbers::RoundNumbers;

pub use crate::{matrix_ops::MatrixOperations, matrix_ops::SquareMatrixOperations};

/// A set of Poseidon2 parameters for a given set of input parameters.
///
/// The const `STATE_SIZE` corresponds to $t$ in the paper, the width of the hash function,
/// e.g. $t=3$ corresponds to a 2-to-1 hash.
///
/// The const `NUM_MDS_ELEMENTS` corresponds to the number of elements in the MDS matrices, which
/// should equal `STATE_SIZE * STATE_SIZE`.
#[derive(Clone, Debug)]
pub struct PoseidonParameters<
    const STATE_SIZE: usize,
    const NUM_MDS_ELEMENTS: usize,
    const NUM_ROUND_ROWS: usize,
    const NUM_ROUND_COLS: usize,
    const NUM_ROUND_ELEMENTS: usize,
> {
    // Input parameters.
    /// Security level.
    pub M: usize,

    // Generated parameters.
    /// Exponent of the Sbox, i.e. S-box(x) = x^{\alpha} used in the `SubWords` step
    pub alpha: Alpha,

    /// Round numbers
    pub rounds: RoundNumbers,

    /// External matrix
    pub m_e: SquareMatrix<STATE_SIZE, NUM_MDS_ELEMENTS>,

    /// Internal matrix
    pub m_i: SquareMatrix<STATE_SIZE, NUM_MDS_ELEMENTS>,

    /// Round constants
    pub arc: ArcMatrix<NUM_ROUND_ROWS, NUM_ROUND_COLS, NUM_ROUND_ELEMENTS>,
}
