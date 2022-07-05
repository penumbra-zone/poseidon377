use ark_ff::PrimeField;
use merlin::Transcript;

use crate::{
    transcript::TranscriptProtocol, Alpha, InputParameters, Matrix, OptimizedMdsMatrices,
    RoundNumbers,
};

/// Represents an matrix of round constants.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArcMatrix<F: PrimeField>(pub Matrix<F>);

impl<F> ArcMatrix<F>
where
    F: PrimeField,
{
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

    pub fn n_rows(&self) -> usize {
        self.0.n_rows
    }

    pub fn n_cols(&self) -> usize {
        self.0.n_cols
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
    pub fn generate(
        arc: &ArcMatrix<F>,
        mds: &OptimizedMdsMatrices<F>,
        rounds: &RoundNumbers,
        t: usize,
    ) -> OptimizedArcMatrix<F> {
        todo!()
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
