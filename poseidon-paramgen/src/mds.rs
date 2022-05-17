use std::collections::HashSet;

use ark_ff::PrimeField;
use merlin::Transcript;
use nalgebra::DMatrix;

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
        // TODO: run algorithms 1-3 to check if matrix is a secure choice
        true
    }

    /// Attempt to generate a `t x t` Cauchy matrix
    ///
    /// For random x_i, y_i in Fp where `i=[0,t)`, the entries are
    /// `1 / (x_i - y_i)`.
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
            for _ in 0..input.t {
                let x_i = transcript.cauchy_coefficient();
                let y_j = transcript.cauchy_coefficient();
                xs.push(x_i);
                ys.push(y_j);
            }

            let mut elements = Vec::<F>::with_capacity(input.t);
            for i in 0..input.t {
                for j in 0..input.t {
                    // Check x_i - y_j != 0
                    if xs[i] - ys[j] == F::zero() {
                        continue 'attempt_loop;
                    }

                    elements.push(F::one() / (xs[i] - ys[j]))
                }
            }

            // Check if pairwise distinct.
            let xi_plus_xj: HashSet<F> = xs
                .clone()
                .into_iter()
                .chain(ys.clone().into_iter())
                .collect();
            if xi_plus_xj.len() != 2 * input.t {
                continue 'attempt_loop;
            }

            let cauchy_matrix = SquareMatrix::from_vec(elements);

            let computed_determinant = cauchy_matrix.determinant();
            // All Cauchy matrices should be invertible
            assert!(computed_determinant != F::zero());

            // Check the determinant with the explicit formula for a Cauchy determinant
            // from Cauchy 1841 p.154 (Exercices d'analyse et de physique mathematique, Vol 2)
            let mut x_prod = F::one();
            for i in 0..input.t {
                for j in 0..i {
                    x_prod *= xs[i] - xs[j];
                }
            }
            let mut y_prod = F::one();
            for i in 0..input.t {
                for j in 0..i {
                    y_prod *= ys[i] - ys[j];
                }
            }
            let mut xy_prod = F::one();
            for i in 0..input.t {
                for j in 0..input.t {
                    xy_prod *= xs[i] - ys[j];
                }
            }

            // Prefix is given by the expression for k at the bottom of p.154 Cauchy 1841
            let prefix = (-F::one()).pow(&[(input.t * (input.t - 1) / 2) as u64]);
            assert_eq!(computed_determinant, prefix * x_prod * y_prod / xy_prod);

            // TODO: Check the matrix has no eigenvalues in Fp
            let matrix: DMatrix<F> = DMatrix::from_vec(
                cauchy_matrix.dim,
                cauchy_matrix.dim,
                cauchy_matrix.elements().to_vec(),
            );
            // Problem: ComplexField must be implemented on F for this to work
            // let all_eigenvalues = matrix.eigenvalues();
            // let complex_eigenvalues = matrix.complex_eigenvalues();
            //assert!(all_eigenvalues.len() - complex_eigenvalues.len() == 0);

            return Self(cauchy_matrix);
        }

        panic!("could not find a valid MDS matrix in Cauchy form")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ed_on_bls12_381::{Fq, FqParameters as Fq381Parameters};
    use ark_ff::{fields::FpParameters, Zero};

    #[test]
    fn cauchy_method_mds() {
        let M = 128;
        let alpha = 5;
        let t = 3;

        let input = InputParameters::new(alpha, M, 3, Fq381Parameters::MODULUS);
        let MDS_matrix: MdsMatrix<Fq> = MdsMatrix::new(&input);

        assert!(MDS_matrix.0.determinant() != Fq::zero());
        assert_eq!(MDS_matrix.0.nrows(), t);
        assert!(MDS_matrix.0.get_element(0, 0) != Fq::zero());
    }
}
