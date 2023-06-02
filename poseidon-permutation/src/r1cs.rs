use ark_ff::{PrimeField, Zero};

use poseidon_parameters::v1::PoseidonParameters;

/// Represents a Poseidon permutation instance.
pub struct InstanceVar<F: PrimeField> {
    /// Parameters for this instance of Poseidon.
    pub parameters: PoseidonParameters,

    /// Constraint system
    pub cs: ConstraintSystemRef<F>,

    /// Current state
    pub state_words: Vec<FpVar<F>>,
}

impl InstanceVar {
    /// Initialize a new Poseidon instance.
    pub fn new(parameters: PoseidonParameters, cs: ConstraintSystemRef<F>) -> Self {
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
        self.unoptimized_permute();
        // self.permute();

        // Emit a single element since this is a n:1 hash.
        self.state_words[1]
    }

    /// Unoptimized Poseidon permutation.
    pub fn unoptimized_permute(&mut self) {
        todo!()
    }

    /// Optimized Poseidon permutation.
    pub fn permute(&mut self) {
        todo!()
    }
}
