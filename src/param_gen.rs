//! Module for generating Poseidon parameters

mod mds;
mod rounds;

use ark_ff::BigInteger256;
use num_bigint::BigUint;

/// A set of Poseidon parameters for a given set of input parameters.
///
/// TODO: Modify this to be the ark-sponge Parameter struct.
struct Instance {
    // Saved input parameters.
    input: InputParameters,

    // Generated parameters.
    rounds: rounds::RoundNumbers,
    //mds: mds::MdsMatrix,
}

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Clone)]
struct InputParameters {
    alpha: i64,
    M: usize,
    t: usize,
    p: BigInteger256,
    // The number of bits needed to represent p.
    n: usize,
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

        let p_biguint: BigUint = p.into();
        let n = p_biguint.bits() as usize;
        let input = InputParameters { alpha, M, t, p, n };
        let rounds = rounds::RoundNumbers::new(input.clone());

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
