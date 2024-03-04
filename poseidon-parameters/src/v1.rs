pub use crate::alpha::Alpha;
pub use crate::round_numbers::RoundNumbers;

pub use crate::{matrix::Matrix, matrix::SquareMatrix};

// pub use crate::matrix_ops::mat_mul;

pub use crate::{
    arc_matrix::ArcMatrix, arc_matrix::OptimizedArcMatrix, matrix_ops::MatrixOperations,
    matrix_ops::SquareMatrixOperations, mds_matrix::MdsMatrix, mds_matrix::OptimizedMdsMatrices,
};

/// A set of Poseidon1 parameters for a given set of input parameters.
///
/// The const `STATE_SIZE` corresponds to $t$ in the paper, the width of the hash function,
/// e.g. $t=3$ corresponds to a 2-to-1 hash.
#[derive(Clone, Debug)]
pub struct PoseidonParameters<
    const STATE_SIZE: usize,
    const STATE_SIZE_MINUS_1: usize,
    const NUM_MDS_ELEMENTS: usize,
    const NUM_STATE_SIZE_MINUS_1_ELEMENTS: usize,
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

    /// `t x t` MDS matrix used in the `MixLayer` step
    pub mds: MdsMatrix<
        STATE_SIZE,
        STATE_SIZE_MINUS_1,
        NUM_MDS_ELEMENTS,
        NUM_STATE_SIZE_MINUS_1_ELEMENTS,
    >,

    /// `num_total_rounds x t` matrix of constants used in the `AddRoundConstant` step
    pub arc: ArcMatrix<NUM_ROUND_ROWS, NUM_ROUND_COLS, NUM_ROUND_ELEMENTS>,

    /// Optimized round constants.
    pub optimized_arc: OptimizedArcMatrix<NUM_ROUND_ROWS, NUM_ROUND_COLS, NUM_ROUND_ELEMENTS>,
    // TODO:
    // /// Optimized MDS matrices.
    //pub optimized_mds: OptimizedMdsMatrices,
}
