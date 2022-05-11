use ark_ff::BigInteger;
use ark_ff::PrimeField;

use crate::{InputParameters, Matrix};

/// The number of attempts to find a secure MDS matrix.
const NUM_ATTEMPTS: usize = 100;

/// Represents an MDS (maximum distance separable) matrix.
pub(super) struct MdsMatrix<F: PrimeField>(Matrix<F>);

impl<F> MdsMatrix<F>
where
    F: PrimeField,
{
    pub fn new<P: BigInteger>(input: &InputParameters<P>) -> Self {
        // A t x t MDS matrix only exists if: 2t + 1 <= p
        let two_times_t_bigint: P = (2 * input.t as u64).into();
        if two_times_t_bigint < input.p {
            panic!("no MDS matrix exists");
        }

        for _ in 0..NUM_ATTEMPTS {
            let candidate = MdsMatrix::cauchy_matrix(input.t);
            if !candidate.is_secure() {
                continue;
            } else {
                return candidate;
            }
        }

        // If we get here, we were not able to find a valid matrix
        panic!("could not find a valid MDS matrix")
    }

    /// Whether this choice of MDS matrix is secure
    fn is_secure(&self) -> bool {
        todo!()
    }

    /// Generate a `t x t` Cauchy matrix
    fn cauchy_matrix(t: usize) -> Self {
        /// For x_i, y_i in Fp where i=[1,t], the entries are
        /// 1 / (x_i + y_i)

        /// The entries of {x_i}_{1<=i<=t} and {y_i}_{1<=i<=t}
        /// are pairwise distinct and x_i + y_i != 0
        /// where i = {1,...,t} and j = {1,...,t}
        ///
        todo!()
    }
}
