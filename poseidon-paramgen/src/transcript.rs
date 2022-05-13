use ark_ff::{BigInteger, PrimeField};
use merlin::Transcript;

use crate::InputParameters;

trait TranscriptProtocol {
    fn domain_sep<F: PrimeField>(&mut self, input: InputParameters<F::BigInt>);
    fn cauchy_coefficient<F: PrimeField>(&mut self) -> F;
}

impl TranscriptProtocol for Transcript {
    fn domain_sep<F: PrimeField>(&mut self, input: InputParameters<F::BigInt>) {
        self.append_message(b"dom-sep", b"poseidon-paramgen");
        // Bind transcript to input parameter choices
        self.append_message(b"t", &input.t.to_le_bytes());
        self.append_message(b"M", &input.M.to_le_bytes());
        self.append_message(b"p", &input.p.to_bytes_le());
        self.append_message(b"alpha", &input.alpha.to_bytes_le());
    }

    fn cauchy_coefficient<F: PrimeField>(&mut self) -> F {
        let size_in_bytes = (F::size_in_bits() + 128) / 8;
        let mut dest = vec![0u8; size_in_bytes];
        self.challenge_bytes(b"cauchy-coefficient", &mut dest);
        F::from_le_bytes_mod_order(&dest)
    }
}
