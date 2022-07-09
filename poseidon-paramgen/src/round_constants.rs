use anyhow::Result;
use ark_ff::PrimeField;
use merlin::Transcript;

use crate::{
    matrix::mat_mul, transcript::TranscriptProtocol, Alpha, InputParameters, Matrix,
    MatrixOperations, OptimizedMdsMatrices, RoundNumbers,
};

/// Represents an matrix of round constants.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArcMatrix<F: PrimeField>(pub Matrix<F>);

impl<F> ArcMatrix<F>
where
    F: PrimeField,
{
    /// Generate round constants.
    pub fn generate(
        input: &InputParameters<F::BigInt>,
        round_numbers: RoundNumbers,
        alpha: Alpha,
    ) -> ArcMatrix<F> {
        let mut transcript = Transcript::new(b"round-constants");
        transcript.domain_sep::<F>(input, round_numbers, alpha);

        let num_total_rounds = round_numbers.total();
        let elements = (0..num_total_rounds * input.t)
            .map(|_| transcript.round_constant())
            .collect();
        ArcMatrix(Matrix::new(num_total_rounds, input.t, elements))
    }

    /// Get row vector of constants by round
    pub(crate) fn constants_by_round(&self, r: usize) -> Matrix<F> {
        self.0.row_vector(r)
    }
}

impl<F: PrimeField> MatrixOperations<F> for ArcMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> ArcMatrix<F> {
        ArcMatrix(Matrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<F> {
        self.0.elements()
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.0.set_element(i, j, val)
    }
    fn rows(&self) -> Vec<&[F]> {
        self.0.rows()
    }

    fn transpose(&self) -> Self {
        ArcMatrix(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self> {
        Ok(ArcMatrix(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F: PrimeField> Into<Vec<Vec<F>>> for ArcMatrix<F> {
    fn into(self) -> Vec<Vec<F>> {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..self.n_rows() {
            let mut row = Vec::new();
            for j in 0..self.n_cols() {
                row.push(self.0.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

/// Represents an optimized matrix of round constants.
///
/// From the text of Appendix B of the Poseidon paper:
///
/// "it is possible to swap the order of the linear layer and the
/// round constant addition as both operations are linear. The round
/// constant then needs to be exchanged with an equivalent one. For
/// round constant $c^{(i)}$ the equivalent one can be written as:
/// $\hat{c}^{(i)} = MC^{-1}(c^{(i)}$ where MC is the linear layer in the
/// i-th round."
///
/// MC is the linear layer where we multiply the state by the MDS matrix.
/// So to compute the optimized round constants, we need the MDS matrix inverse.
///
/// This modifies the partial rounds in the middle of the permutation,
/// wherein you add constants _first_ before iterating through the partial
/// rounds.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedArcMatrix<F: PrimeField>(pub Matrix<F>);

impl<F> OptimizedArcMatrix<F>
where
    F: PrimeField,
{
    /// Generate the optimized round constants.
    pub fn generate(
        arc: &ArcMatrix<F>,
        mds: &OptimizedMdsMatrices<F>,
        rounds: &RoundNumbers,
        t: usize,
    ) -> OptimizedArcMatrix<F> {
        let n_rows = arc.n_rows();
        let n_cols = arc.n_cols();
        let mut new_constants = Vec::with_capacity(n_rows * n_cols);
        let original_constants = arc.0.clone().elements;

        let r_f = rounds.full() / 2;
        let r_T = rounds.total();
        let r_P = rounds.partial();

        // First round is unchanged
        for i in 0..t {
            new_constants.push(original_constants[i])
        }

        // Next r_f - 1 rounds are multiplied by Minv
        for r in 1..r_f - 1 {
            // Multiply row vector (1 x t) by Minv (t x t)
            let new_round_constants = mat_mul(&arc.constants_by_round(r), &mds.M_inverse.0)
                .expect("parameter matrices have expected dimensions");
            assert_eq!(new_round_constants.n_rows, 1);
            assert_eq!(new_round_constants.n_cols, t);
            new_constants.extend(new_round_constants.elements);
        }

        let mut partial_constants = Vec::with_capacity(r_P);
        let mut acc = arc.constants_by_round(r_f + r_P);

        for r in (r_f..r_f + r_P).rev() {
            let mut acc_prime = mat_mul(&acc, &mds.M_inverse.0)
                .expect("parameter matrices have expected dimensions");
            partial_constants.push(acc_prime.get_element(0, 0));
            acc_prime.set_element(0, 0, F::zero());
            acc = acc_prime
                .hadamard_product(&arc.constants_by_round(r))
                .expect("parameter matrices have expected dimensions");
        }

        let acc_times_m_inv =
            mat_mul(&acc, &mds.M_inverse.0).expect("parameter matrices have expected dimensions");
        new_constants.extend(acc_times_m_inv.elements);
        let final_partial_constants: Vec<F> = partial_constants.into_iter().rev().collect();
        new_constants.extend(final_partial_constants);

        // Final r_f rounds are multiplied by Minv
        for r in r_T - r_f..r_T {
            let new_round_constants = mat_mul(&arc.constants_by_round(r), &mds.M_inverse.0)
                .expect("parameter matrices have expected dimensions");
            new_constants.extend(new_round_constants.elements);
        }

        // The unoptimized constants are a matrix of size r x t.
        // However, the optimized constants are a vector consisting of:
        // * First t * r_f elements: applied to full rounds
        // * Next 1 * r_P elements: applied to partial rounds
        // * Final t * r_f elements: applied to full rounds
        let expected_optimized_constant_length = t * 2 * r_f + r_P;
        assert_eq!(new_constants.len(), expected_optimized_constant_length);
        OptimizedArcMatrix(Matrix::new(
            1,
            expected_optimized_constant_length,
            new_constants,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ed_on_bls12_381::Fq;

    #[test]
    fn convert_from_arc_to_vec_of_vecs() {
        let arc_matrix = ArcMatrix(Matrix::new(
            2,
            3,
            vec![
                Fq::from(1u32),
                Fq::from(2u32),
                Fq::from(0u32),
                Fq::from(4u32),
                Fq::from(5u32),
                Fq::from(6u32),
            ],
        ));
        let vec_of_vecs: Vec<Vec<Fq>> = arc_matrix.into();
        assert_eq!(vec_of_vecs[0][0], Fq::from(1u32));
        assert_eq!(vec_of_vecs[0][1], Fq::from(2u32));
        assert_eq!(vec_of_vecs[0][2], Fq::from(0u32));
        assert_eq!(vec_of_vecs[1][0], Fq::from(4u32));
        assert_eq!(vec_of_vecs[1][1], Fq::from(5u32));
        assert_eq!(vec_of_vecs[1][2], Fq::from(6u32));
    }
}
