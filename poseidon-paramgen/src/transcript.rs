use ark_ff::PrimeField;

use merlin::Transcript;

trait TranscriptProtocol {
    fn domain_sep(&mut self);
    fn cauchy_coefficient<F: PrimeField>(&mut self) -> F;
}

impl TranscriptProtocol for Transcript {
    // todo: take inputParameters
    fn domain_sep(&mut self) {
        self.append_message(b"dom-sep", b"poseidon-paramgen");
        //self.append_message(b"prime", xxxx);
    }

    fn cauchy_coefficient<F: PrimeField>(&mut self) -> F {
        let size_in_bytes = (F::size_in_bits() + 128) / 8;
        let mut dest = vec![0u8; size_in_bytes];
        self.challenge_bytes(b"cauchy-coefficient", &mut dest);
        F::from_le_bytes_mod_order(&dest)
    }
}
