#![allow(non_snake_case)]

use decaf377::Fq;
use poseidon_parameters::v1::{Alpha, MatrixOperations, PoseidonParameters};

//temp
use poseidon_parameters::StuffThatNeedsToGoInDecaf377;

/// Represents a generic instance of `Poseidon`.
///
/// Intended for generic fixed-width hashing.
pub struct Instance<
    'a,
    const STATE_SIZE: usize,
    const STATE_SIZE_MINUS_1: usize,
    const NUM_MDS_ELEMENTS: usize,
    const NUM_STATE_SIZE_MINUS_1_ELEMENTS: usize,
    const NUM_ROUND_ROWS: usize,
    const NUM_ROUND_COLS: usize,
    const NUM_ROUND_ELEMENTS: usize,
    const NUM_PARTIAL_ROUNDS: usize,
> {
    /// Parameters for this instance of Poseidon.
    parameters: &'a PoseidonParameters<
        STATE_SIZE,
        STATE_SIZE_MINUS_1,
        NUM_MDS_ELEMENTS,
        NUM_STATE_SIZE_MINUS_1_ELEMENTS,
        NUM_ROUND_ROWS,
        NUM_ROUND_COLS,
        NUM_ROUND_ELEMENTS,
        NUM_PARTIAL_ROUNDS,
    >,

    /// Inner state.
    state_words: [Fq; STATE_SIZE],
}

impl<
        'a,
        const STATE_SIZE: usize,
        const STATE_SIZE_MINUS_1: usize,
        const NUM_MDS_ELEMENTS: usize,
        const NUM_STATE_SIZE_MINUS_1_ELEMENTS: usize,
        const NUM_ROUND_ROWS: usize,
        const NUM_ROUND_COLS: usize,
        const NUM_ROUND_ELEMENTS: usize,
        const NUM_PARTIAL_ROUNDS: usize,
    >
    Instance<
        'a,
        STATE_SIZE,
        STATE_SIZE_MINUS_1,
        NUM_MDS_ELEMENTS,
        NUM_STATE_SIZE_MINUS_1_ELEMENTS,
        NUM_ROUND_ROWS,
        NUM_ROUND_COLS,
        NUM_ROUND_ELEMENTS,
        NUM_PARTIAL_ROUNDS,
    >
{
    /// Instantiate a new hash function over Fq given `Parameters`.
    pub fn new(
        parameters: &'a PoseidonParameters<
            STATE_SIZE,
            STATE_SIZE_MINUS_1,
            NUM_MDS_ELEMENTS,
            NUM_STATE_SIZE_MINUS_1_ELEMENTS,
            NUM_ROUND_ROWS,
            NUM_ROUND_COLS,
            NUM_ROUND_ELEMENTS,
            NUM_PARTIAL_ROUNDS,
        >,
    ) -> Self {
        Self {
            parameters,
            state_words: [Fq::from(0u64); STATE_SIZE],
        }
    }

    /// Fixed width hash from n:1. Outputs a F given `t` input words.
    pub fn n_to_1_fixed_hash(&mut self, input_words: &[Fq; STATE_SIZE]) -> Fq {
        // Set internal state words.
        for (i, input_word) in input_words.iter().enumerate() {
            self.state_words[i] = *input_word
        }

        // Apply Poseidon permutation.
        self.permute();

        // Emit a single element since this is a n:1 hash.
        self.state_words[1]
    }

    /// Print out internal state.
    pub fn output_words(&self) -> [Fq; STATE_SIZE] {
        self.state_words
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
            for i in 0..STATE_SIZE {
                self.state_words[i] += self.parameters.optimized_arc.0.get_element(r, i);
            }
            self.full_sub_words();
            self.mix_layer_mds();
        }
        let mut round_constants_counter = R_f;

        // Partial rounds
        // First part of `AddRoundConstants` layer
        for i in 0..STATE_SIZE {
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
            for i in 0..STATE_SIZE {
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
    pub fn unoptimized_n_to_1_fixed_hash(&mut self, input_words: [Fq; STATE_SIZE]) -> Fq {
        // Set internal state words.
        for (i, input_word) in input_words.iter().enumerate() {
            self.state_words[i] = *input_word
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
        let round_constants = self.parameters.arc.elements();

        // First full rounds
        for _ in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..STATE_SIZE {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.full_sub_words();
            self.mix_layer_mds();
        }

        // Partial rounds
        for _ in 0..R_P {
            // Apply `AddRoundConstants` layer
            for i in 0..STATE_SIZE {
                self.state_words[i] += round_constants[round_constants_counter];
                round_constants_counter += 1;
            }
            self.partial_sub_words();
            self.mix_layer_mds();
        }

        // Final full rounds
        for _ in 0..R_f {
            // Apply `AddRoundConstants` layer
            for i in 0..STATE_SIZE {
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
            Alpha::Exponent(exp) => self.state_words[0] = (self.state_words[0]).pow([exp as u64]),
            Alpha::Inverse => self.state_words[0] = Fq::from(1u64) / self.state_words[0],
        }
    }

    /// Applies the full `SubWords` layer.
    fn full_sub_words(&mut self) {
        match self.parameters.alpha {
            Alpha::Exponent(exp) => {
                for i in 0..STATE_SIZE {
                    self.state_words[i] = self.state_words[i].pow([exp as u64]);
                }
            }
            Alpha::Inverse => {
                for i in 0..STATE_SIZE {
                    self.state_words[i] = Fq::from(1u64) / self.state_words[i];
                }
            }
        }
    }

    /// Applies the `MixLayer` using the M_i matrix.
    fn mix_layer_mi(&mut self) {
        let mut new_state_words = [Fq::from(0u64); STATE_SIZE];
        for (i, row) in self.parameters.optimized_mds.M_i.iter_rows().enumerate() {
            let sum = row
                .iter()
                .zip(&self.state_words)
                .map(|(x, y)| *x * *y)
                .sum();
            new_state_words[i] = sum;
        }
        self.state_words = new_state_words;
    }

    /// Applies the `MixLayer` using the MDS matrix.
    fn mix_layer_mds(&mut self) {
        let mut new_state_words = [Fq::from(0u64); STATE_SIZE];

        for (i, row) in self.parameters.mds.0 .0.iter_rows().enumerate() {
            let sum = row
                .iter()
                .zip(&self.state_words)
                .map(|(x, y)| *x * *y)
                .sum();
            new_state_words[i] = sum;
        }
        self.state_words = new_state_words;
    }

    /// This is `cheap_matrix_mul` in the Sage spec
    fn sparse_mat_mul(&mut self, round_number: usize) {
        // mul_row = [(state_words[0] * v[i]) for i in range(0, t-1)]
        // add_row = [(mul_row[i] + state_words[i+1]) for i in range(0, t-1)]
        let mut add_row = [Fq::from(0u64); STATE_SIZE_MINUS_1];
        for (i, x) in self.parameters.optimized_mds.v_collection[round_number]
            .elements
            .iter()
            .enumerate()
        {
            add_row[i] = *x * self.state_words[0] + self.state_words[i + 1];
        }

        // column_1 = [M_0_0] + w_hat
        // state_words_new[0] = sum([column_1[i] * state_words[i] for i in range(0, t)])
        // state_words_new = [state_words_new[0]] + add_row
        self.state_words[0] = self.parameters.optimized_mds.M_00 * self.state_words[0]
            + self.parameters.optimized_mds.w_hat_collection[round_number]
                .elements
                .iter()
                .zip(self.state_words[1..STATE_SIZE].iter())
                .map(|(x, y)| *x * *y)
                .sum::<Fq>();

        self.state_words[1..STATE_SIZE].copy_from_slice(&add_row[..(STATE_SIZE - 1)]);
    }
}
