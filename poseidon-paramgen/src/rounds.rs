use crate::input::InputParameters;
use ark_ff::BigInteger;
use ark_std::cmp::{Ordering, PartialOrd};
use num_bigint::{BigInt, ToBigInt};
use poseidon_parameters::v1::{Alpha, RoundNumbers};

/// Generate round numbers.
///
/// For Poseidon 1, we currently panic if the security level is not at
/// the 128-bit or 256-bit level. This is because in the work by Ashur,
/// Buschman, and Mahzoun 2023, they find the Grobner basis attacks are
/// stronger than described in the original Poseidon paper, however they
/// only find partial and full attacks beyond the 256-bit level [0]. The
/// parameter generation logic for 128-bit and 256-bit security levels
/// is unchanged.
///
/// [0]: https://eprint.iacr.org/2023/537
pub fn v1_generate<T: BigInteger>(input: &InputParameters<T>, alpha: &Alpha) -> RoundNumbers {
    let mut choice: Option<RoundNumbers> = None;
    let mut cost = usize::MAX;
    let mut cost_rf = usize::MAX;

    if input.M > 256 {
        panic!("This crate does not support security levels beyond 256 bits, use Poseidon2");
    }

    // Loop through choices of r_F, r_P
    for r_P in 1..400 {
        for r_F in 4..100 {
            let mut candidate = RoundNumbersBuilder(RoundNumbers { r_F, r_P });
            if !candidate.is_secure_v1(input, alpha) {
                continue;
            }

            candidate.apply_security_margin();
            let candidate_cost = candidate.sbox_count(input.t);
            // Pick the minimum cost Candidate, and if the cost is tied with another
            // candidate, we switch to the new candidate if the total number of full rounds is lower.
            if (candidate_cost < cost) || ((candidate_cost == cost) && (r_F < cost_rf)) {
                cost = candidate_cost;
                cost_rf = r_F;
                choice = Some(candidate.0);
            }
        }
    }

    choice.unwrap()
}

/// Generate round numbers for Poseidon2.
pub fn v2_generate<T: BigInteger>(input: &InputParameters<T>, alpha: &Alpha) -> RoundNumbers {
    let mut choice: Option<RoundNumbers> = None;
    let mut cost = usize::MAX;
    let mut cost_rf = usize::MAX;

    // Loop through choices of r_F, r_P
    for r_P in 1..400 {
        for r_F in 4..100 {
            let mut candidate = RoundNumbersBuilder(RoundNumbers { r_F, r_P });
            if !candidate.is_secure_v2(input, alpha) {
                continue;
            }

            candidate.apply_security_margin();
            let candidate_cost = candidate.sbox_count(input.t);
            // Pick the minimum cost Candidate, and if the cost is tied with another
            // candidate, we switch to the new candidate if the total number of full rounds is lower.
            if (candidate_cost < cost) || ((candidate_cost == cost) && (r_F < cost_rf)) {
                cost = candidate_cost;
                cost_rf = r_F;
                choice = Some(candidate.0);
            }
        }
    }

    choice.unwrap()
}

struct RoundNumbersBuilder(pub RoundNumbers);

impl RoundNumbersBuilder {
    /// Check if this `RoundNumbers` choice is secure given all known attacks.
    fn is_secure_v1<T: BigInteger>(&self, input: &InputParameters<T>, alpha: &Alpha) -> bool {
        // Check if the number of full rounds are sufficient.
        if self.0.full() < statistical_attack_full_rounds(input, alpha) {
            return false;
        }

        match alpha {
            // For positive alpha, the interpolation and Grobner bounds are on the total
            // number of rounds.
            Alpha::Exponent(_) => {
                if self.0.total() <= algebraic_attack_interpolation(input, alpha) {
                    return false;
                }
                if self.0.total() <= algebraic_attack_grobner_basis_v1(input, alpha) {
                    return false;
                }
            }
            // For inverse alpha, the interpolation and Grobner bounds are on r_F scaled
            // by the binary log of `t` plus r_P. See Eqn 4.
            Alpha::Inverse => {
                if (self.0.full() as f64 * (input.t as f64).log2()).floor() as usize
                    + self.0.partial()
                    <= algebraic_attack_interpolation(input, alpha)
                {
                    return false;
                }
                if (self.0.full() as f64 * (input.t as f64).log2()).floor() as usize
                    + self.0.partial()
                    <= algebraic_attack_grobner_basis_v1(input, alpha)
                {
                    return false;
                }
            }
        }

        true
    }

    /// Check if this `RoundNumbers` choice is secure given all known attacks.
    fn is_secure_v2<T: BigInteger>(&self, input: &InputParameters<T>, alpha: &Alpha) -> bool {
        // Check if the number of full rounds are sufficient.
        if self.0.full() < statistical_attack_full_rounds(input, alpha) {
            return false;
        }

        match alpha {
            // For positive alpha, the interpolation and Grobner bounds are on the total
            // number of rounds.
            Alpha::Exponent(_) => {
                if self.0.total() <= algebraic_attack_interpolation(input, alpha) {
                    return false;
                }
                if self.0.total() <= algebraic_attack_grobner_basis_v1(input, alpha) {
                    return false;
                }
                if algebraic_attack_grobner_basis_v2_possible(input, alpha, &self.0) {
                    return false;
                }
            }
            Alpha::Inverse => unimplemented!("no support for inverse alpha!"),
        }

        true
    }

    /// Get the number of SBoxes for these `RoundNumbers` on a permutation of width `t`.
    pub(crate) fn sbox_count(&self, t: usize) -> usize {
        t * self.0.full() + self.0.partial()
    }

    /// Add suggested security margin of +2 R_F and +7.5% R_P.
    /// Ref: Section 5.4.
    fn apply_security_margin(&mut self) {
        *self.0.full_mut() += 2;
        *self.0.partial_mut() = (1.075 * (self.0.partial() as f64)).ceil() as usize;
    }
}

/// Number of full rounds required to defend against statistical attacks.
///
/// These are the differential/linear distinguisher attacks described
/// in Section 5.5.1 of the paper.
fn statistical_attack_full_rounds<T: BigInteger>(
    input: &InputParameters<T>,
    alpha: &Alpha,
) -> usize {
    // C is defined in Section 5.5.1, p.10.
    let C = match alpha {
        Alpha::Inverse => 2.0,
        Alpha::Exponent(exp) => (*exp as f64 - 1.0).log2(),
    };

    // Statistical attacks require at least 6 full rounds.
    // Differential/Linear Distinguishers.
    // Ref: Section 5.5.1.
    if input.M as f64 <= ((input.log_2_p.floor() - C) * (input.t as f64 + 1.0)) {
        6
    } else {
        10
    }
}

/// Number of total rounds to defend against interpolation attacks for positive alpha.
/// For negative alpha, this bound is $\floor{R_F \log_2(t)} + R_P$ rounds to defend against
/// interpolation attacks.
///
/// These attacks are described in Section 5.5.2 of the paper.
/// For positive alpha, we use Eqn 3.
/// For negative alpha, we use Eqn 4.
fn algebraic_attack_interpolation<T: BigInteger>(
    input: &InputParameters<T>,
    alpha: &Alpha,
) -> usize {
    let min_args = [input.M as f64, input.log_2_p];
    match alpha {
        Alpha::Inverse => {
            return (((input.t as f64).log2()).ceil()
                + (0.5 * min_args.iter().min_by(cmp_f64).expect("no NaNs")).ceil())
                as usize;
        }
        Alpha::Exponent(exp) => {
            return ((2f64.log(*exp as f64) * min_args.iter().min_by(cmp_f64).expect("no NaNs"))
                .ceil()
                + (input.t as f64).log(*exp as f64))
            .ceil() as usize;
        }
    };
}

/// If the improved Grobner basis attacks from https://eprint.iacr.org/2023/537 are possible.
///
/// See: https://github.com/HorizenLabs/poseidon2/commit/44bdcbc37887390442c7e743bad655a7ab8a7b7d
fn algebraic_attack_grobner_basis_v2_possible<T: BigInteger>(
    input: &InputParameters<T>,
    alpha: &Alpha,
    choice: &RoundNumbers,
) -> bool {
    let r = (input.t as f64 / 3.0).floor() as usize;
    let exponent = match alpha {
        Alpha::Inverse => unimplemented!("no support for inverse alpha!"),
        Alpha::Exponent(exp) => *exp as usize,
    };

    // Equation from Item 1, Section 4.3 in https://eprint.iacr.org/2023/537
    let upper = ((choice.full() - 1) * input.t)
        + choice.partial()
        + r
        + r * (choice.full() / 2)
        + choice.partial()
        + exponent;
    let lower = r * (choice.full() / 2) + choice.partial() + exponent;

    // Start computing the binomial coefficient (upper choose lower).
    let mut acc = BigInt::from(1);
    for i in 0..lower {
        acc = (acc
            * (upper - i)
                .to_bigint()
                .expect("always succeeds if operating on unsigned ints"))
            / (i + 1)
                .to_bigint()
                .expect("always succeeds if operating on unsigned ints");

        // 2 in the below expression is the minimum value of the linear
        // algebra constant $\omega$ defined in section 4.1
        // of https://eprint.iacr.org/2023/537. Note that in the equation presented
        // in section 4.3 of the paper, they use the maximum value of the linear
        // algebra constant, however the minimum value is the more conservative choice.
        //
        // If we find a binomial coeff bitsize that is equal to or beyond the security
        // level, we return early to avoid computing the full binomial coefficient.
        if (2.0 * acc.bits() as f64).ceil() as usize >= input.M {
            return false;
        }
    }

    // If we get here, it means the attack is possible.
    true
}

/// Number of total rounds to defend against Grobner basis attacks.
///
/// These are described in Section 5.5.2 of the paper.
///
/// We use the first two conditions described in Section C.2.2,
/// eliding the third since if the first condition is satisfied, then
/// the third will be also.
fn algebraic_attack_grobner_basis_v1<T: BigInteger>(
    input: &InputParameters<T>,
    alpha: &Alpha,
) -> usize {
    let grobner_1: f64;
    let grobner_2: f64;

    match alpha {
        Alpha::Inverse => {
            // First Grobner constraint
            let grobner_1_min_args = [input.M as f64, input.log_2_p];
            grobner_1 = (0.5f64 * grobner_1_min_args.iter().min_by(cmp_f64).expect("no NaNs"))
                .ceil()
                + ((input.t as f64).log2()).ceil();
            // Second Grobner constraint
            let grobner_2_min_args = [
                (input.M as f64 / (input.t as f64 + 1.0)).ceil(),
                (0.5 * input.log_2_p).ceil(),
            ];
            grobner_2 = (input.t - 1) as f64
                + ((input.t as f64).log2()).ceil()
                + grobner_2_min_args.iter().min_by(cmp_f64).expect("no NaNs");
        }
        Alpha::Exponent(exp) => {
            // First Grobner constraint
            let grobner_1_min_args = [(input.M as f64 / 3.0), (input.log_2_p / 2.0)];
            grobner_1 =
                2f64.log(*exp as f64) * grobner_1_min_args.iter().min_by(cmp_f64).expect("no NaNs");
            // Second Grobner constraint
            let grobner_2_min_args = [
                2f64.log(*exp as f64) * input.M as f64 / (input.t as f64 + 1.0),
                2f64.log(*exp as f64) * input.log_2_p / 2.0,
            ];
            grobner_2 =
                (input.t - 1) as f64 + grobner_2_min_args.iter().min_by(cmp_f64).expect("no NaNs");
        }
    };

    // Return the _most_ strict Grobner basis constraint.
    let grobner_values = [grobner_1, grobner_2];
    return grobner_values
        .iter()
        .max_by(cmp_f64)
        .expect("no NaNs")
        .floor() as usize;
}

fn cmp_f64(lhs: &&f64, rhs: &&f64) -> Ordering {
    lhs.partial_cmp(rhs).unwrap()
}
