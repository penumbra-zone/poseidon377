//! Module for generating Poseidon parameters

mod mds;
mod rounds;

use ark_ff::BigInteger256;

/// A set of Poseidon parameters for a given set of input parameters.
struct Instance {
    // Saved input parameters.
    input: InputParameters,

    // Generated parameters.
    rounds: rounds::RoundNumbers,
    //mds: mds::MdsMatrix,
}

/// Input parameters that are used to generate Poseidon parameters.
struct InputParameters {
    alpha: i64,
    M: usize,
    t: usize,
    p: BigInteger256,
}

impl Instance {
    /// Generate a Poseidon instance mapped over Fp given a choice of:
    ///
    /// * $\alpha$,
    /// * M, a desired security level (in bits),
    /// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    /// * p, the prime p,
    pub fn new(M: usize, alpha: i64, t: usize, p: BigInteger256) -> Self {
        // Alpha must be a positive odd integer (p.10), or -1.
        if alpha != -1 || alpha < 1 || alpha % 2 == 0 {
            panic!("invalid value for alpha: {}", alpha);
        }

        let input = InputParameters { alpha, M, t, p };
        let rounds = rounds::RoundNumbers::new(input);

        // TODO: MDS matrix

        Self { input, rounds }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ed_on_bls12_377::Fq;
    use ark_ff::fields::FpParameters;

    #[test]
    fn poseidon_instance() {
        let params_rate_3 = Instance::new(128, 17, 3, Fq::MODULUS);
    }
}
