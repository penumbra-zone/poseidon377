use ark_ff::PrimeField;
use ark_std::{vec, vec::Vec};

use ark_r1cs_std::{fields::fp::FpVar, prelude::*};
use ark_relations::r1cs::ConstraintSystemRef;
use poseidon_parameters::v1::{Alpha, MatrixOperations, PoseidonParameters};

/// Represents a Poseidon permutation instance.
pub struct InstanceVar<F: PrimeField> {
    /// Parameters for this instance of Poseidon.
    pub parameters: PoseidonParameters<F>,

    /// Constraint system
    pub cs: ConstraintSystemRef<F>,

    /// Current state
    pub state_words: Vec<FpVar<F>>,
}

impl<F> InstanceVar<F>
where
    F: PrimeField,
{
    /// Initialize a new Poseidon instance.
    pub fn new(parameters: PoseidonParameters<F>, cs: ConstraintSystemRef<F>) -> Self {
        let zero = FpVar::<F>::zero();
        // t = rate + capacity
        let state_words = vec![zero; parameters.t];

        Self {
            parameters,
            cs,
            state_words,
        }
    }

    /// Fixed width hash from n:1. Outputs a F given `t` input words.
    pub fn n_to_1_fixed_hash(&mut self, input_words: Vec<FpVar<F>>) -> FpVar<F> {
        // Check input words are `t` elements long
        if input_words.len() != self.parameters.t {
            panic!("err: input words must be t elements long")
        }

        // Set internal state words.
        for (i, input_word) in input_words.into_iter().enumerate() {
            self.state_words[i] = input_word
        }

        // Apply Poseidon permutation.
        self.permute();

        // Emit a single element since this is a n:1 hash.
        self.state_words[1].clone()
    }

    /// Fixed width hash from n:1. Outputs a F given `t` input words.
    pub fn unoptimized_n_to_1_fixed_hash(&mut self, input_words: Vec<FpVar<F>>) -> FpVar<F> {
        // Check input words are `t` elements long
        if input_words.len() != self.parameters.t {
            panic!("err: input words must be t elements long")
        }

        // Set internal state words.
        for (i, input_word) in input_words.into_iter().enumerate() {
            self.state_words[i] = input_word
        }

        // Apply Poseidon permutation.
        self.unoptimized_permute();

        // Emit a single element since this is a n:1 hash.
        self.state_words[1].clone()
    }

    /// Unoptimized Poseidon permutation.
    #[allow(non_snake_case)]
    pub fn unoptimized_permute(&mut self) {
        let R_f = self.parameters.rounds.full() / 2;
        let R_P = self.parameters.rounds.partial();
        let mut round_constants_counter = 0;
        let t = self.parameters.t;
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

    /// Permutes the internal state.
    ///
    /// This implementation is based on the optimized Sage implementation
    /// `poseidonperm_x3_64_optimized.sage` provided in Appendix B of the Poseidon paper.
    fn permute(&mut self) {
        let R_f = self.parameters.rounds.full() / 2;

        // First chunk of full rounds
        for r in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..self.parameters.t {
                self.state_words[i] += self.parameters.optimized_arc.0.get_element(r, i);
            }
            self.full_sub_words();
            self.mix_layer_mds();
        }
        let mut round_constants_counter = R_f;

        // Partial rounds
        // First part of `AddRoundConstants` layer
        for i in 0..self.parameters.t {
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
            for i in 0..self.parameters.t {
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

    /// Applies the partial `SubWords` layer.
    fn partial_sub_words(&mut self) {
        match self.parameters.alpha {
            Alpha::Exponent(exp) => {
                self.state_words[0] = (self.state_words[0])
                    .pow_by_constant([exp as u64])
                    .expect("can compute pow")
            }
            Alpha::Inverse => unimplemented!("err: inverse alpha not implemented"),
        }
    }

    /// Applies the full `SubWords` layer.
    fn full_sub_words(&mut self) {
        match self.parameters.alpha {
            Alpha::Exponent(exp) => {
                self.state_words = self
                    .state_words
                    .iter()
                    .map(|x| x.pow_by_constant([exp as u64]).expect("can compute pow"))
                    .collect()
            }
            Alpha::Inverse => {
                unimplemented!("err: inverse alpha not implemented")
            }
        }
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
                let temp_vec: Vec<FpVar<F>> = row
                    .iter()
                    .zip(&self.state_words)
                    .map(|(x, y)| {
                        FpVar::<F>::new_constant(self.cs.clone(), x).expect("can create constant")
                            * y
                    })
                    .collect();
                let result = temp_vec.iter().sum();
                result
            })
            .collect();
    }

    /// Applies the `MixLayer` using the M_i matrix.
    fn mix_layer_mi(&mut self) {
        self.state_words = self
            .parameters
            .optimized_mds
            .M_i
            .iter_rows()
            .map(|row| {
                let temp_vec: Vec<FpVar<F>> = row
                    .iter()
                    .zip(&self.state_words)
                    .map(|(x, y)| {
                        FpVar::<F>::new_constant(self.cs.clone(), x).expect("can create constant")
                            * y
                    })
                    .collect();
                let result = temp_vec.iter().sum();
                result
            })
            .collect();
    }

    /// This is `cheap_matrix_mul` in the Sage spec
    fn sparse_mat_mul(&mut self, round_number: usize) {
        // mul_row = [(state_words[0] * v[i]) for i in range(0, t-1)]
        // add_row = [(mul_row[i] + state_words[i+1]) for i in range(0, t-1)]
        let add_row: Vec<FpVar<F>> = self.parameters.optimized_mds.v_collection[round_number]
            .elements
            .iter()
            .enumerate()
            .map(|(i, x)| {
                FpVar::<F>::new_constant(self.cs.clone(), x).expect("can create constant")
                    * self.state_words[0].clone()
                    + self.state_words[i + 1].clone()
            })
            .collect();

        // column_1 = [M_0_0] + w_hat
        // state_words_new[0] = sum([column_1[i] * state_words[i] for i in range(0, t)])
        // state_words_new = [state_words_new[0]] + add_row
        let temp_vec: Vec<FpVar<F>> = self.parameters.optimized_mds.w_hat_collection[round_number]
            .elements
            .iter()
            .zip(self.state_words[1..self.parameters.t].iter())
            .map(|(x, y)| {
                FpVar::<F>::new_constant(self.cs.clone(), x).expect("can create constant") * y
            })
            .collect();
        let result: FpVar<F> = temp_vec.iter().sum();
        let M00_var = FpVar::<F>::new_constant(self.cs.clone(), self.parameters.optimized_mds.M_00)
            .expect("can create constant");
        self.state_words[0] = M00_var * self.state_words[0].clone() + result;

        // self.state_words[1..self.parameters.t].copy_from_slice(&add_row[..(self.parameters.t - 1)]);
        for index in 1..self.parameters.t {
            self.state_words[index] = add_row[index - 1].clone();
        }
    }
}
