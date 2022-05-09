//! Module for generating Poseidon parameters
use std::convert::From;
use std::convert::TryFrom;
use std::convert::TryInto;
mod addition_chains;
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

/// Take the binary log of a `BigInteger`
pub fn log2<P>(x: P) -> f64
where
    P: BigInteger,
{
    let mut p_biguint: BigUint = x.into();
    let two_to_50: BigUint = 1125899906842624u64.into(); // 2**50
    let mut log_bit_boundaries = 0;
    while p_biguint > two_to_50 {
        p_biguint >>= 48u64;
        log_bit_boundaries += 48;
    }
    let mut p_le_bytes = p_biguint.to_bytes_le();
    // `u64::from_le_bytes` takes `[u8; 8]` but `p_le_bytes` is a Vec<u8>
    // so we must pad the vec to the proper size:
    while p_le_bytes.len() < 8 {
        p_le_bytes.push(0x00);
    }

    let x_u64: u64 = u64::from_le_bytes(
        p_le_bytes
            .try_into()
            .expect("we have padded to the expected length"),
    );
    log_bit_boundaries as f64 + ((x_u64) as f64).log2()
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ff::BigInteger256;

    //use ark_ed_on_bls12_377::Fq;
    use ark_ed_on_bls12_381::FqParameters as Fq381Parameters;
    use ark_ff::fields::FpParameters;

    #[test]
    fn log2_bigint() {
        let test_val: BigInteger256 = 4.into();
        assert_eq!(log2(test_val), 2.0);

        let test_val: BigInteger256 = 257.into();
        assert!(log2(test_val) > 8.005);
        assert!(log2(test_val) < 8.006);

        let test_val: BigInteger256 = 65536.into();
        assert_eq!(log2(test_val), 16.0);

        let test_val = Fq381Parameters::MODULUS;
        assert!(log2(test_val) > 254.856);
        assert!(log2(test_val) < 254.858);
    }

    /// Represents a row in Table 7-9 in Appendix G of the paper.
    struct TableRow {
        // Security margin
        M: usize,
        // Text size: N = n * t
        N: usize,
        // S-box size: n
        n: usize,
        // Number of S-boxes (t)
        t: usize,
        // Number of full rounds
        r_F: usize,
        // Number of partial rounds
        r_P: usize,
        // Cost
        cost: usize,
    }

    #[test]
    #[ignore]
    fn table_7_concrete_instances_alpha_3() {
        todo!()
    }

    #[test]
    #[ignore]
    fn table_8_concrete_instances_alpha_5() {
        let alpha = 5;
        let N = 1536;

        // Appendix G of the paper provides in Table 8 concrete instances of x^{5} Poseidon
        let table_8 = [
            [128, N, 768, 2, 8, 56, 72],
            [128, N, 384, 4, 8, 56, 88],
            [128, N, 256, 6, 8, 57, 105],
            [128, N, 192, 8, 8, 57, 121],
            [128, N, 96, 16, 8, 42, 170],
            [256, N, 768, 2, 8, 116, 132],
            [256, N, 384, 4, 8, 116, 148],
            [256, N, 256, 6, 8, 117, 165],
            [256, N, 192, 8, 8, 86, 150],
            [256, N, 96, 16, 8, 42, 170],
        ];

        for row in table_8 {
            let table_row = TableRow {
                M: row[0],
                N: row[1],
                n: row[2],
                t: row[3],
                r_F: row[4],
                r_P: row[5],
                cost: row[6],
            };
            let computed_instance =
                PoseidonParameters::new(table_row.M, alpha, table_row.t, Fq381Parameters::MODULUS);
            assert_eq!(computed_instance.rounds.full(), table_row.r_F);
            assert_eq!(computed_instance.rounds.partial(), table_row.r_P);
        }
    }

    #[test]
    #[ignore]
    fn table_9_concrete_instances_inverse_alpha() {
        todo!()
    }

    #[test]
    fn poseidon_bls12_377_instance() {
        //let params_rate_3 = PoseidonParameters::new(128, 17, 3, Fq::MODULUS);
    }
}
