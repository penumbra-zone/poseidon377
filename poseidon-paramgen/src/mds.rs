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

        MdsMatrix::fixed_cauchy_matrix(&input)
    }

    /// Generate a deterministic Cauchy matrix
    ///
    /// The original Poseidon paper describes a method for constructing MDS matrices
    /// from randomly selecting $x_i$, $y_j$ from the field. The resulting matrix needs
    /// to be passed to algorithms 1-3 described in
    /// [Grassi, Rechberger, Schofnegger 2020](https://eprint.iacr.org/archive/2020/500/20200702:141143)
    /// in order to determine if infinitely long subspace trails can be constructed for
    /// the Cauchy matrix. If yes, then the MDS matrix must be thrown away, and the process
    /// must begin again for another random choice of $x_i$, $y_j$.
    ///
    /// However, Section 5.4 of [Keller and Rosemarin 2020](https://eprint.iacr.org/2020/179.pdf)
    /// describes how the MDS matrix can be constructed in a deterministic fashion
    /// where infinitely long subspace trails cannot be constructed. This method constructs an MDS
    /// matrix using that method.
    pub fn fixed_cauchy_matrix(input: &InputParameters<F::BigInt>) -> Self {
        let xs: Vec<F> = (0..input.t as u64).map(F::from).collect();
        let ys: Vec<F> = (input.t as u64..2 * input.t as u64).map(F::from).collect();

        let mut elements = Vec::<F>::with_capacity(input.t);
        for i in 0..input.t {
            for j in 0..input.t {
                // Check x_i + y_j != 0
                assert_ne!(xs[i] + ys[i], F::zero());
                elements.push(F::one() / (xs[i] + ys[j]))
            }
        }

        let cauchy_matrix = SquareMatrix::from_vec(elements);
        // Sanity check: All Cauchy matrices should be invertible
        assert!(cauchy_matrix.determinant() != F::zero());

        Self(cauchy_matrix)
    }

    /// Attempt to generate a `t x t` Cauchy matrix using random sampling in Fp.
    ///
    /// For random x_i, y_j in Fp where `i,j=[0,t)`, the entries are
    /// `1 / (x_i - y_j)`.
    ///
    /// The entries of {x_i}_{0<=i<t} and {y_j}_{0<=j<t}
    /// are pairwise distinct and
    /// where i = {0,...,t-1} and j = {0,...,t-1}.
    ///
    /// # Security
    ///
    /// As described in the original Poseidon paper, these MDS matrices
    /// must be run through Algorithms 1-3 in
    /// [Grassi, Rechberger, Schofnegger 2020](https://eprint.iacr.org/archive/2020/500/20200702:141143).
    fn random_cauchy_matrix(input: &InputParameters<F::BigInt>) -> Self {
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
