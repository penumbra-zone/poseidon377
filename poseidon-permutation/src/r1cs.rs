#![allow(non_snake_case)]
use ark_std::vec::Vec;

use ark_r1cs_std::{fields::fp::FpVar, prelude::*};
use ark_relations::r1cs::ConstraintSystemRef;
use decaf377::Fq;
use poseidon_parameters::v1::{Alpha, MatrixOperations, PoseidonParameters};

/// Represents a Poseidon permutation instance.
pub struct InstanceVar<
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
    pub parameters: PoseidonParameters<
        STATE_SIZE,
        STATE_SIZE_MINUS_1,
        NUM_MDS_ELEMENTS,
        NUM_STATE_SIZE_MINUS_1_ELEMENTS,
        NUM_ROUND_ROWS,
        NUM_ROUND_COLS,
        NUM_ROUND_ELEMENTS,
        NUM_PARTIAL_ROUNDS,
    >,

    /// Constraint system
    pub cs: ConstraintSystemRef<Fq>,

    /// Current state
    pub state_words: Vec<FpVar<Fq>>,
}

impl<
        const STATE_SIZE: usize,
        const STATE_SIZE_MINUS_1: usize,
        const NUM_MDS_ELEMENTS: usize,
        const NUM_STATE_SIZE_MINUS_1_ELEMENTS: usize,
        const NUM_ROUND_ROWS: usize,
        const NUM_ROUND_COLS: usize,
        const NUM_ROUND_ELEMENTS: usize,
        const NUM_PARTIAL_ROUNDS: usize,
    >
    InstanceVar<
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
    /// Fixed width hash from n:1. Outputs a Fq given `t` input words.
    pub fn n_to_1_fixed_hash(
        parameters: PoseidonParameters<
            STATE_SIZE,
            STATE_SIZE_MINUS_1,
            NUM_MDS_ELEMENTS,
            NUM_STATE_SIZE_MINUS_1_ELEMENTS,
            NUM_ROUND_ROWS,
            NUM_ROUND_COLS,
            NUM_ROUND_ELEMENTS,
            NUM_PARTIAL_ROUNDS,
        >,
        cs: ConstraintSystemRef<Fq>,
        input_words: [FpVar<Fq>; STATE_SIZE],
    ) -> FpVar<Fq> {
        // t = rate + capacity

        let mut instance = InstanceVar {
            parameters,
            cs,
            state_words: input_words.to_vec(),
        };

        // Apply Poseidon permutation.
        instance.permute();

        // Emit a single element since this is a n:1 hash.
        instance.state_words[1].clone()
    }

    /// Poseidon permutation.
    pub fn permute(&mut self) {
        let R_f = self.parameters.rounds.full() / 2;
        let R_P = self.parameters.rounds.partial();
        let mut round_constants_counter = 0;
        let round_constants: [Fq; NUM_ROUND_ELEMENTS] = self.parameters.arc.inner_elements();

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
                for i in 0..STATE_SIZE {
                    self.state_words[i] = (self.state_words[i])
                        .pow_by_constant([exp as u64])
                        .expect("can compute pow");
                }
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
                let temp_vec: Vec<FpVar<Fq>> = row
                    .iter()
                    .zip(&self.state_words)
                    .map(|(x, y)| {
                        FpVar::<Fq>::new_constant(self.cs.clone(), x).expect("can create constant")
                            * y
                    })
                    .collect();
                let result = temp_vec.iter().sum();
                result
            })
            .collect();
    }
}
