use std::collections::HashSet;

use ark_ff::PrimeField;
use merlin::Transcript;

use crate::{transcript::TranscriptProtocol, InputParameters, SquareMatrix};

/// The number of attempts to find a secure MDS matrix.
const NUM_ATTEMPTS: usize = 100;

/// Represents an MDS (maximum distance separable) matrix.
pub struct MdsMatrix<F: PrimeField>(pub SquareMatrix<F>);

impl<F> MdsMatrix<F>
where
    F: PrimeField,
{
    pub fn new(input: &InputParameters<F::BigInt>) -> Self {
        // A t x t MDS matrix only exists if: 2t + 1 <= p
        let two_times_t_bigint: F::BigInt = (2 * input.t as u64).into();
        if two_times_t_bigint > input.p {
            panic!("no MDS matrix exists");
        }

        for _ in 0..NUM_ATTEMPTS {
            let candidate = MdsMatrix::cauchy_matrix(&input);
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

    /// Attempt to generate a `t x t` Cauchy matrix
    ///
    /// For random x_i, y_i in Fp where `i=[0,t)`, the entries are
    /// `1 / (x_i + y_i)`.
    ///
    /// The entries of {x_i}_{0<=i<t} and {y_i}_{0<=i<t}
    /// are pairwise distinct and
    /// where i = {0,...,t-1} and j = {0,...,t-1}.
    fn cauchy_matrix(input: &InputParameters<F::BigInt>) -> Self {
        'attempt_loop: for _ in 0..NUM_ATTEMPTS {
            let mut transcript = Transcript::new(b"cauchy-matrix");
            transcript.domain_sep::<F>(input);
            let mut xs = Vec::with_capacity(input.t);
            let mut ys = Vec::with_capacity(input.t);
            for n in 0..input.t {
                let x_i = transcript.cauchy_coefficient();
                let y_j = transcript.cauchy_coefficient();
                xs.push(x_i);
                ys.push(y_j);
            }

            let mut elements = Vec::<Vec<F>>::with_capacity(input.t);
            for i in 0..input.t {
                let mut row = Vec::<F>::with_capacity(input.t);
                for j in 0..input.t {
                    // Check x_i + y_j != 0
                    if xs[i] + ys[j] == F::zero() {
                        continue 'attempt_loop;
                    }

                    row.push(F::one() / (xs[i] + ys[j]))
                }
                elements.push(row);
            }

            // Check if pairwise distinct.
            let xi_plus_xj: HashSet<F> = xs.into_iter().chain(ys.into_iter()).collect();
            if xi_plus_xj.len() != 2 * input.t {
                continue 'attempt_loop;
            }

            // TODO: Check Cauchy determinant is nonzero

            return Self(SquareMatrix::new(elements));
        }

        panic!("could not find a valid MDS matrix in Cauchy form")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PoseidonParameters;

    use ark_ed_on_bls12_381::FqParameters as Fq381Parameters;
    use ark_ff::fields::FpParameters;

    #[test]
    fn cauchy_method_mds() {
        // let M = 128;
        // let alpha = 5;
        // let t = 3;
        // let params_rate_3 = PoseidonParameters::new(M, alpha, t, Fq381Parameters::MODULUS);
        // params_rate_3.mds
    }
}
