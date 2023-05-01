use ark_ff::PrimeField;

pub use crate::alpha::Alpha;
pub use crate::round_numbers::RoundNumbers;

pub use crate::{matrix::Matrix, matrix::SquareMatrix, matrix_ops::mat_mul};

pub use crate::{
    arc_matrix::ArcMatrix, arc_matrix::OptimizedArcMatrix, matrix_ops::MatrixOperations,
    matrix_ops::SquareMatrixOperations, mds_matrix::MdsMatrix, mds_matrix::OptimizedMdsMatrices,
};

// TODO: arc and mds could be collections
/// A set of Poseidon1 parameters for a given set of input parameters.
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
