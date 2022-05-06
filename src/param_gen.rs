//! Module for generating Poseidon parameters

mod mds;
mod rounds;

use ark_ff::BigInteger;
use num_bigint::BigUint;

#[derive(Clone)]
pub enum Alpha {
    Exponent(u32),
    Inverse,
}

/// A set of Poseidon parameters for a given set of input parameters.
///
/// TODO(later): Add an Into for converting this to be the ark-sponge Parameter struct.
struct PoseidonParameters<P: BigInteger> {
    // Saved input parameters.
    input: InputParameters<P>,

    // Generated parameters.
    pub rounds: rounds::RoundNumbers,
    //mds: mds::MdsMatrix,
}

/// Input parameters that are used to generate Poseidon parameters.
#[derive(Clone)]
struct InputParameters<P: BigInteger> {
    alpha: Alpha, // TODO: Choose best alpha based on choice of p.
    M: usize,
    t: usize,
    p: P,
    // The number of bits needed to represent p.
    n: usize,
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

        let p_biguint: BigUint = p.into();
        let n = p_biguint.bits() as usize;
        let input = InputParameters {
            alpha: alpha_var,
            M,
            t,
            p,
            n,
        };
        let rounds = rounds::RoundNumbers::new(&input);

        // TODO: MDS matrix

        Self { input, rounds }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    //use ark_ed_on_bls12_377::Fq;
    use ark_ed_on_bls12_381::FqParameters as Fq381Parameters;
    use ark_ff::fields::FpParameters;

    #[test]
    fn poseidon_bls12_381_instance() {
        let alpha = 5;
        let security_margin = 128;

        let params_rate_3 =
            PoseidonParameters::new(security_margin, alpha, 3, Fq381Parameters::MODULUS);
        dbg!(params_rate_3.rounds.total());
        dbg!(params_rate_3.rounds.full());
        dbg!(params_rate_3.rounds.partial());
        assert_eq!(params_rate_3.rounds.full(), 8);
        assert_eq!(params_rate_3.rounds.partial(), 57);

        let params_rate_5 =
            PoseidonParameters::new(security_margin, alpha, 5, Fq381Parameters::MODULUS);
        dbg!(params_rate_5.rounds.total());
        dbg!(params_rate_5.rounds.full());
        dbg!(params_rate_5.rounds.partial());
        assert_eq!(params_rate_5.rounds.full(), 8);
        // assert_eq!(params_rate_5.rounds.partial(), 60);
    }

    #[test]
    fn poseidon_bls12_377_instance() {
        //let params_rate_3 = PoseidonParameters::new(128, 17, 3, Fq::MODULUS);
    }
}
