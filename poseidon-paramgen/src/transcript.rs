use ark_ff::{BigInteger, PrimeField};
use ark_std::vec;
use merlin::Transcript;

use crate::{Alpha, InputParameters, RoundNumbers};

pub(crate) trait TranscriptProtocol {
    fn domain_sep<F: PrimeField>(
        &mut self,
        input: &InputParameters<F::BigInt>,
        round_numbers: RoundNumbers<poseidon_parameters::RoundNumbers>,
        alpha: Alpha<poseidon_parameters::Alpha>,
    );
    fn round_constant<F: PrimeField>(&mut self) -> F;
}

impl TranscriptProtocol for Transcript {
    fn domain_sep<F: PrimeField>(
        &mut self,
        input: &InputParameters<F::BigInt>,
        round_numbers: RoundNumbers<poseidon_parameters::RoundNumbers>,
        alpha: Alpha<poseidon_parameters::Alpha>,
    ) {
        self.append_message(b"dom-sep", b"poseidon-paramgen");
        // Bind transcript to input parameter choices
        self.append_message(b"t", &input.t.to_le_bytes());
        self.append_message(b"M", &input.M.to_le_bytes());
        self.append_message(b"p", &input.p.to_bytes_le());

        // Bind transcript also to specific instance as done with the Grain LFSR
        // in Appendix F of the Poseidon paper.
        self.append_message(b"r_F", &[round_numbers.full() as u8]);
        self.append_message(b"r_P", &[round_numbers.partial() as u8]);
        self.append_message(b"alpha", &alpha.to_bytes_le());
    }

    fn round_constant<F: PrimeField>(&mut self) -> F {
        let size_in_bytes = (F::size_in_bits() + 135) / 8;
        let mut dest = vec![0u8; size_in_bytes];
        self.challenge_bytes(b"round-constant", &mut dest);
        F::from_le_bytes_mod_order(&dest)
    }
}
