#![allow(non_snake_case)]

use ark_ff::PrimeField;

use poseidon_paramgen::{Alpha, MatrixOperations, PoseidonParameters};

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
        self.permute();

        // Emit a single element since this is a n:1 hash.
        self.state_words[1]
    }

    #[cfg(test)]
    pub(crate) fn output_words(&self) -> Vec<F> {
        self.state_words.clone()
    }

    /// Permutes the internal state.
    ///
    /// This implementation is based on the optimized Sage implementation
    /// `poseidonperm_x3_64_optimized.sage` provided in Appendix B of the Poseidon paper.
    fn permute(&mut self) {
        let R_f = self.parameters.rounds.full() / 2;

        // First chunk of full rounds
        for r in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..self.parameters.input.t {
                self.state_words[i] += self.parameters.optimized_arc.0.get_element(r, i);
            }
            self.full_sub_words();
            self.mix_layer_mds();
        }
        let mut round_constants_counter = R_f;

        // Partial rounds
        // First part of `AddRoundConstants` layer
        for i in 0..self.parameters.input.t {
            self.state_words[i] += self
                .parameters
                .optimized_arc
                .0
                .get_element(round_constants_counter, i);
        }
        // First full matrix multiplication.
        self.mix_layer_mi();

        for r in 0..self.parameters.rounds.partial() - 1 {
            self.partial_sub_words();
            // Rest of `AddRoundConstants` layer, moved to after the S-box layer
            round_constants_counter += 1;
            self.state_words[0] += self
                .parameters
                .optimized_arc
                .0
                .get_element(round_constants_counter, 0);
            self.sparse_mat_mul(self.parameters.rounds.partial() - r - 1);
        }
        // Last partial round
        self.partial_sub_words();
        self.sparse_mat_mul(0);
        round_constants_counter += 1;

        // Final full rounds
        for _ in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..self.parameters.input.t {
                self.state_words[i] += self
                    .parameters
                    .optimized_arc
                    .0
                    .get_element(round_constants_counter, i);
            }
            self.full_sub_words();
            self.mix_layer_mds();
            round_constants_counter += 1;
        }
    }

    /// Fixed width hash from n:1. Outputs a F given `t` input words. Unoptimized.
    pub fn unoptimized_n_to_1_fixed_hash(&mut self, input_words: Vec<F>) -> F {
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
            self.mix_layer_mds();
        }

        // Partial rounds
        for _ in 0..R_P {
            // Apply `AddRoundConstants` layer
            for i in 0..t {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.partial_sub_words();
            self.mix_layer_mds();
        }

        // Final full rounds
        for _ in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..t {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.full_sub_words();
            self.mix_layer_mds();
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
        match self.parameters.alpha {
            Alpha::Exponent(exp) => {
                self.state_words = self
                    .state_words
                    .iter()
                    .map(|x| x.pow(&[exp as u64]))
                    .collect()
            }
            Alpha::Inverse => {
                self.state_words = self.state_words.iter().map(|x| F::one() / x).collect()
            }
        }
    }

    /// Applies the `MixLayer` using the M_i matrix.
    fn mix_layer_mi(&mut self) {
        self.state_words = self
            .parameters
            .optimized_mds
            .M_i
            .iter_rows()
            .map(|row| {
                row.iter()
                    .zip(&self.state_words)
                    .map(|(x, y)| *x * *y)
                    .sum()
            })
            .collect();
    }

    /// Applies the `MixLayer` using the MDS matrix.
    fn mix_layer_mds(&mut self) {
        self.state_words = self
            .parameters
            .mds
            .0
             .0
            .iter_rows()
            .map(|row| {
                row.iter()
                    .zip(&self.state_words)
                    .map(|(x, y)| *x * *y)
                    .sum()
            })
            .collect();
    }

    /// This is `cheap_matrix_mul` in the Sage spec
    fn sparse_mat_mul(&mut self, round_number: usize) {
        // mul_row = [(state_words[0] * v[i]) for i in range(0, t-1)]
        // add_row = [(mul_row[i] + state_words[i+1]) for i in range(0, t-1)]
        let add_row: Vec<F> = self.parameters.optimized_mds.v_collection[round_number]
            .elements()
            .iter()
            .enumerate()
            .map(|(i, x)| *x * self.state_words[0] + self.state_words[i + 1])
            .collect();

        // column_1 = [M_0_0] + w_hat
        // state_words_new[0] = sum([column_1[i] * state_words[i] for i in range(0, t)])
        // state_words_new = [state_words_new[0]] + add_row
        self.state_words[0] = self.parameters.optimized_mds.M_00 * self.state_words[0]
            + self.parameters.optimized_mds.w_hat_collection[round_number]
                .elements()
                .iter()
                .zip(self.state_words[1..self.parameters.input.t].iter())
                .map(|(x, y)| *x * *y)
                .sum::<F>();

        self.state_words[1..self.parameters.input.t]
            .copy_from_slice(&add_row[..(self.parameters.input.t - 1)]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    use ark_ed_on_bls12_377::{Fq, FqParameters};
    use ark_ff::FpParameters;
    use ark_ff::{PrimeField, Zero};
    use ark_sponge::poseidon::{Parameters, State};
    use poseidon_paramgen::PoseidonParameters;

    #[test]
    fn check_optimized_impl_vs_sage() {
        let params_2_to_1 = PoseidonParameters::<Fq>::new(128, 3, FqParameters::MODULUS, true);
        let mut our_instance = Instance::new(params_2_to_1);
        let hash_output =
            our_instance.n_to_1_fixed_hash(vec![Fq::zero(), Fq::from(1u64), Fq::from(2u64)]);
        let output_words = our_instance.output_words();
        assert_eq!(hash_output, output_words[1]);
        let expected_output_words = [
            ark_ff::field_new!(
                Fq,
                "1005395416230692022226189338173977461720389945215833518989246316298938732662"
            ),
            ark_ff::field_new!(
                Fq,
                "334624358114973565971415282870268792208067693392598174364241659220644131605"
            ),
            ark_ff::field_new!(
                Fq,
                "4551703942608256274690841944452321045558690573870761888368849297930709301300"
            ),
        ];
        for (a, b) in expected_output_words.iter().zip(output_words.iter()) {
            assert_eq!(*a, *b);
        }
    }

    #[test]
    fn check_unoptimized_impl_vs_sage() {
        let params_2_to_1 = PoseidonParameters::<Fq>::new(128, 3, FqParameters::MODULUS, true);
        let mut our_instance = Instance::new(params_2_to_1);
        let hash_output = our_instance.unoptimized_n_to_1_fixed_hash(vec![
            Fq::zero(),
            Fq::from(1u64),
            Fq::from(2u64),
        ]);
        let output_words = our_instance.output_words();
        assert_eq!(hash_output, output_words[1]);
        let expected_output_words = [
            ark_ff::field_new!(
                Fq,
                "1005395416230692022226189338173977461720389945215833518989246316298938732662"
            ),
            ark_ff::field_new!(
                Fq,
                "334624358114973565971415282870268792208067693392598174364241659220644131605"
            ),
            ark_ff::field_new!(
                Fq,
                "4551703942608256274690841944452321045558690573870761888368849297930709301300"
            ),
        ];
        for (a, b) in expected_output_words.iter().zip(output_words.iter()) {
            assert_eq!(*a, *b);
        }
    }

    fn fq_strategy() -> BoxedStrategy<Fq> {
        any::<[u8; 32]>()
            .prop_map(|bytes| Fq::from_le_bytes_mod_order(&bytes[..]))
            .boxed()
    }

    proptest! {
        #[test]
        fn ark_sponge_and_permutation_consistent(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy()) {
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

        #[test]
        fn optimized_and_unoptimized_permutation_consistent(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy()) {
            let params_2_to_1 = PoseidonParameters::<Fq>::new(128, 3, FqParameters::MODULUS, true);

            let mut our_instance = Instance::new(params_2_to_1.clone());
            let our_result = our_instance.n_to_1_fixed_hash(vec![elem_1, elem_2, elem_3]);

            let mut unoptimized_instance = Instance::new(params_2_to_1);
            let unoptimized_result =
                unoptimized_instance.unoptimized_n_to_1_fixed_hash(vec![elem_1, elem_2, elem_3]);

            assert_eq!(unoptimized_result, our_result);
        }
    }
}
