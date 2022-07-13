#![allow(non_snake_case)]

use ark_ff::PrimeField;

use poseidon_paramgen::{mat_mul, Alpha, Matrix, MatrixOperations, PoseidonParameters};

/// Represents a generic instance of `Poseidon`.
///
/// Intended for generic fixed-width hashing.
pub struct Instance<F: PrimeField> {
    /// Parameters for this instance of Poseidon.
    parameters: PoseidonParameters<F>,

    /// Inner state.
    state_words: Vec<F>,
}

impl<F: PrimeField> Instance<F> {
    /// Instantiate a new hash function over GF(p) given `Parameters`.
    pub fn new(parameters: PoseidonParameters<F>) -> Self {
        let t = parameters.input.t;
        Self {
            parameters,
            state_words: vec![F::zero(); t],
        }
    }

    /// Fixed width hash from n:1. Outputs a F given `t` input words.
    pub fn n_to_1_fixed_hash(&mut self, input_words: Vec<F>) -> F {
        // Check input words are `t` elements long
        if input_words.len() != self.parameters.input.t {
            panic!("err: input words must be t elements long")
        }

        // Set internal state words.
        for (i, input_word) in input_words.into_iter().enumerate() {
            self.state_words[i] = input_word
        }

        // Apply Poseidon permutation.
        self.unoptimized_permute();

        // Emit a single element since this is a n:1 hash.
        self.state_words[1]
    }

    /// Permutes the internal state.
    ///
    /// This implementation is based on the unoptimized Sage implementation
    /// `poseidonperm_x5_254_3.sage` provided in Appendix B of the Poseidon paper.
    fn unoptimized_permute(&mut self) {
        let R_f = self.parameters.rounds.full() / 2;
        let R_P = self.parameters.rounds.partial();
        let mut round_constants_counter = 0;
        let t = self.parameters.input.t;
        let round_constants = self.parameters.arc.elements().clone();

        // First full rounds
        for _ in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..t {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.full_sub_words();
            self.mix_layer();
        }

        // Partial rounds
        for _ in 0..R_P {
            // Apply `AddRoundConstants` layer
            for i in 0..t {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.partial_sub_words();
            self.mix_layer();
        }

        // Final full rounds
        for _ in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..t {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.full_sub_words();
            self.mix_layer();
        }
    }

    /// Applies the partial `SubWords` layer.
    fn partial_sub_words(&mut self) {
        match self.parameters.alpha {
            Alpha::Exponent(exp) => self.state_words[0] = (self.state_words[0]).pow(&[exp as u64]),
            Alpha::Inverse => self.state_words[0] = F::one() / self.state_words[0],
        }
    }

    /// Applies the full `SubWords` layer.
    fn full_sub_words(&mut self) {
        for i in 0..self.parameters.input.t {
            match self.parameters.alpha {
                Alpha::Exponent(exp) => {
                    self.state_words[i] = (self.state_words[i]).pow(&[exp as u64])
                }
                Alpha::Inverse => self.state_words[i] = F::one() / self.state_words[i],
            }
        }
    }

    /// Applies the `MixLayer`.
    fn mix_layer(&mut self) {
        let mds_matrix = &self.parameters.mds.0;
        let t = self.parameters.input.t;

        let state_words_col_vector = Matrix::new(t, 1, self.state_words.clone());
        let new_state_words = mat_mul(&mds_matrix.0, &state_words_col_vector).expect(
            "MDS matrix and state words column vector should have correct matrix dimensions",
        );
        // Set new state words.
        for (i, new_state_word) in new_state_words.elements().iter().enumerate() {
            self.state_words[i] = *new_state_word
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    use ark_ed_on_bls12_377::{Fq, FqParameters};
    use ark_ff::FpParameters;
    use ark_sponge::poseidon::{Parameters, State};

    fn fq_strategy() -> BoxedStrategy<Fq> {
        any::<[u8; 32]>()
            .prop_map(|bytes| Fq::from_le_bytes_mod_order(&bytes[..]))
            .boxed()
    }

    proptest! {
        #[test]
        fn ark_sponge_and_unoptimized_permutation_consistent(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy()) {
            let params_2_to_1 = PoseidonParameters::<Fq>::new(128, 3, FqParameters::MODULUS, true);

            let params_ark: Parameters<Fq> = params_2_to_1.clone().into();
            let mut ark_state = State::from(params_ark);
            ark_state[0] = elem_1;
            ark_state[1] = elem_2;
            ark_state[2] = elem_3;
            ark_state.permute();
            let ark_result = ark_state[1];

            let mut our_instance = Instance::new(params_2_to_1);
            let our_result = our_instance.n_to_1_fixed_hash(vec![elem_1, elem_2, elem_3]);

            assert_eq!(ark_result, our_result);
        }
    }
}

// TODO: Vectorize
// TODO: Optimized permutation
