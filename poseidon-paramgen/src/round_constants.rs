use ark_ff::PrimeField;
use merlin::Transcript;

use crate::{transcript::TranscriptProtocol, Alpha, InputParameters, Matrix, RoundNumbers};

/// Represents an matrix of round constants.
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
        ArcMatrix(Matrix::new(input.t, num_total_rounds, elements))
    }
}
