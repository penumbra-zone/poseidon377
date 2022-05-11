//! Module for generating Poseidon parameters
mod addition_chains;
mod matrix;
mod mds;
mod rounds;
mod utils;

use matrix::Matrix;

use ark_ff::BigInteger;

use crate::utils::log2;

#[derive(Clone)]
pub enum Alpha {
    Exponent(u32),
    Inverse,
}

/// A set of Poseidon parameters for a given set of input parameters.
///
/// This is based on the attacks described in the original Poseidon paper [1].
///
/// References:
/// * Original Poseidon paper: https://eprint.iacr.org/2019/458.pdf
/// * MDS Matrices https://eprint.iacr.org/2020/500/20200702:141143
///
/// TODO(later): Add an Into for converting this to be the ark-sponge Parameter struct.
pub struct PoseidonParameters<P: BigInteger> {
    // Saved input parameters.
    input: InputParameters<P>,

    // Generated parameters.
    pub rounds: rounds::RoundNumbers,
    //mds: mds::MdsMatrix,
}

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Clone)]
pub struct InputParameters<P: BigInteger> {
    alpha: Alpha, // TODO: Choose best alpha based on choice of p.
    M: usize,
    t: usize,
    p: P,

    // The below are derived values, stored for convenience.
    /// log_2(p)
    log_2_p: f64,
}

impl<P> PoseidonParameters<P>
where
    P: BigInteger,
{
    /// Generate a Poseidon instance mapped over Fp given a choice of:
    ///
    /// * $\alpha$,
    /// * M, a desired security level (in bits),
    /// * t, the width of the desired hash function, e.g. $t=3$ corresponds to 2-to-1 hash.
    /// * p, the prime p,
    pub fn new(M: usize, alpha: i64, t: usize, p: P) -> Self {
        // Alpha must be a positive odd integer (p.10), or -1.
        let alpha_var: Alpha;
        if alpha == -1 {
            alpha_var = Alpha::Inverse;
        } else if alpha > 1 && alpha % 2 != 0 {
            alpha_var = Alpha::Exponent(alpha as u32)
        } else {
            panic!("invalid value for alpha: {}", alpha);
        }

        let log_2_p = log2(p);
        let input = InputParameters {
            alpha: alpha_var,
            M,
            t,
            p,
            log_2_p,
        };
        let rounds = rounds::RoundNumbers::new(&input);

        // TODO: MDS matrix

        Self { input, rounds }
    }
}
