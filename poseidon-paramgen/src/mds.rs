use ark_ff::PrimeField;

use crate::{InputParameters, SquareMatrix};

/// Represents an MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
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

        MdsMatrix::fixed_cauchy_matrix(input)
    }

    /// Dimension of the (square) MDS matrix
    pub fn dim(&self) -> usize {
        self.0.inner.n_rows
    }

    /// Generate a deterministic Cauchy matrix
    ///
    /// The original Poseidon paper describes a method for constructing MDS matrices
    /// from randomly selecting $x_i$, $y_j$ from the field and then constructing each element in the
    /// matrix using $1/(x_i + y_j)$. The resulting Cauchy matrix needs to then be passed to algorithms 1-3
    /// described in
    /// [Grassi, Rechberger, Schofnegger 2020](https://eprint.iacr.org/archive/2020/500/20200702:141143)
    /// in order to determine if infinitely long subspace trails can be constructed for
    /// the Cauchy matrix. If yes, then the MDS matrix must be thrown away, and the process
    /// must begin again for another random choice of $x_i$, $y_j$ until a secure choice is found.
    ///
    /// However, Section 5.4 of [Keller and Rosemarin 2020](https://eprint.iacr.org/2020/179.pdf)
    /// describes how the MDS matrix can be constructed in a deterministic fashion
    /// where infinitely long subspace trails cannot be constructed. This method constructs an MDS
    /// matrix using that deterministic method, avoiding the need to implement Algorithms 1-3.
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
}

impl<F: PrimeField> Into<Vec<Vec<F>>> for MdsMatrix<F> {
    fn into(self) -> Vec<Vec<F>> {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..self.dim() {
            let mut row = Vec::new();
            for j in 0..self.dim() {
                row.push(self.0.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ed_on_bls12_381::{Fq, FqParameters as Fq381Parameters};
    use ark_ff::{fields::FpParameters, Zero};

    #[test]
    fn convert_from_mds_to_vec_of_vecs() {
        let MDS_matrix = MdsMatrix(SquareMatrix::from_vec(vec![
            Fq::from(1u32),
            Fq::from(2u32),
            Fq::from(3u32),
            Fq::from(4u32),
        ]));
        let vec_of_vecs: Vec<Vec<Fq>> = MDS_matrix.into();
        assert_eq!(vec_of_vecs[0][0], Fq::from(1u32));
        assert_eq!(vec_of_vecs[0][1], Fq::from(2u32));
        assert_eq!(vec_of_vecs[1][0], Fq::from(3u32));
        assert_eq!(vec_of_vecs[1][1], Fq::from(4u32));
    }

    #[test]
    fn cauchy_method_mds() {
        let M = 128;
        let t = 3;

        let input = InputParameters::new(M, 3, Fq381Parameters::MODULUS, true);
        let MDS_matrix: MdsMatrix<Fq> = MdsMatrix::new(&input);

        assert!(MDS_matrix.0.determinant() != Fq::zero());
        assert_eq!(MDS_matrix.dim(), t);
        assert!(MDS_matrix.0.get_element(0, 0) != Fq::zero());
    }
}
