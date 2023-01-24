use anyhow::Result;
use ark_ff::PrimeField;
use ark_std::{vec, vec::Vec};
use merlin::Transcript;

use crate::{
    matrix::mat_mul, transcript::TranscriptProtocol, Alpha, InputParameters, Matrix,
    MatrixOperations, MdsMatrix, RoundNumbers,
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
        round_numbers: RoundNumbers<poseidon_parameters::RoundNumbers>,
        alpha: Alpha<poseidon_parameters::Alpha>,
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
    pub fn constants_by_round(&self, r: usize) -> Matrix<F> {
        self.0.row_vector(r)
    }

    /// Set row vector of constants by round
    pub(crate) fn set_constants_by_round(&mut self, r: usize, constants: Vec<F>) {
        assert_eq!(constants.len(), self.n_cols());
        for (j, value) in constants.into_iter().enumerate() {
            self.set_element(r, j, value);
        }
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

impl<F: PrimeField> Into<Vec<Vec<F>>> for OptimizedArcMatrix<F> {
    fn into(self) -> Vec<Vec<F>> {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..self.0.n_rows() {
            let mut row = Vec::new();
            for j in 0..self.0.n_cols() {
                row.push(self.0.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

/// Represents an optimized matrix of round constants.
///
/// This modifies the partial rounds in the middle of the permutation,
/// wherein you add constants _first_ before iterating through the partial
/// rounds.
///
/// This method follows `calc_equivalent_constants` from Appendix B's
/// `poseidonperm_x3_64_24_optimized.sage`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedArcMatrix<F: PrimeField>(pub ArcMatrix<F>);

impl<F> OptimizedArcMatrix<F>
where
    F: PrimeField,
{
    /// Generate the optimized round constants.
    pub fn generate(
        arc: &ArcMatrix<F>,
        mds: &MdsMatrix<F>,
        rounds: &RoundNumbers<poseidon_parameters::RoundNumbers>,
    ) -> OptimizedArcMatrix<F> {
        let n_cols = arc.n_cols();
        let mut constants_temp = arc.clone();
        let r_f = rounds.full() / 2;
        let r_T = rounds.total();
        let mds_T = mds.transpose();
        let mds_inv = &mds_T.inverse();

        // C_i = M^-1 * C_(i+1)
        for r in ((r_f)..(r_T - 1 - r_f)).rev() {
            // inv_cip1 = list(vector(constants_temp[i+1]) * MDS_matrix_field_transpose.inverse())
            let inv_cip1 = mat_mul(&constants_temp.constants_by_round(r + 1), &mds_inv.0).expect(
                "matrix multiplication of row of ARC matrix and MDS must have correct dimensions",
            );
            assert_eq!(inv_cip1.n_cols(), n_cols);
            assert_eq!(inv_cip1.n_rows(), 1);

            // constants_temp[i] = list(vector(constants_temp[i]) + vector([0] + inv_cip1[1:]))
            let mut delta_new_constants_r = Vec::with_capacity(n_cols);
            delta_new_constants_r.push(F::zero());
            for j in 1..n_cols {
                delta_new_constants_r.push(inv_cip1.get_element(0, j));
            }
            for j in 0..n_cols {
                let curr_element = constants_temp.get_element(r, j);
                constants_temp.set_element(r, j, curr_element + delta_new_constants_r[j]);
            }

            // constants_temp[i+1] = [inv_cip1[0]] + [0] * (t-1)
            let mut new_constants_row_i_plus_1 = Vec::with_capacity(n_cols);
            new_constants_row_i_plus_1.push(inv_cip1.get_element(0, 0));
            new_constants_row_i_plus_1.extend(vec![F::zero(); n_cols - 1]);
            constants_temp.set_constants_by_round(r + 1, new_constants_row_i_plus_1);
        }

        OptimizedArcMatrix(constants_temp)
    }

    /// Create a `OptimizedArcMatrix` from its elements.
    pub fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> OptimizedArcMatrix<F> {
        OptimizedArcMatrix(ArcMatrix::new(n_rows, n_cols, elements))
    }
}

#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_377::Fq;

    use super::*;

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
