use ark_ed_on_bls12_377::Fq;
use ark_ff::One;
use ark_sponge::{poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge};

use crate::rate_2;

// This struct will let us replace the implementation of the inner
// PoseidonSponge if we later choose to do so.
struct Sponge {
    inner: PoseidonSponge<Fq>,
    // TODOs: Domain separation, mode variable, padding fcn, cache?
}

impl Sponge {
    fn new() -> Self {
        let mut sponge = PoseidonSponge::<Fq>::new(&rate_2::params());
        Sponge { inner: sponge }
    }

    /// Take a single field element into the sponge.
    fn absorb(&mut self, element: Fq) {
        self.inner.absorb(&element)
    }

    /// Produce a single field element.
    fn squeeze(&mut self) -> Fq {
        self.inner.squeeze_native_field_elements(1)[0]
    }

    /// Hash variable-length input into a hash.
    pub fn hash(&mut self, message: Vec<Fq>, out_len: usize) -> Vec<Fq> {
        for i in 0..message.len() {
            self.absorb(message[i]);
        }

        // Domain separation
        self.absorb(Fq::one());

        let mut output = Vec::<Fq>::new();
        for _i in 0..out_len {
            output.push(self.squeeze());
        }
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_variable_len_hash() {
        let out_len = 1;
        let mut sponge = Sponge::new();
        let message = vec![Fq::one(), Fq::one()];
        let result = sponge.hash(message, out_len);
        assert_eq!(result.len(), out_len);
    }
}
