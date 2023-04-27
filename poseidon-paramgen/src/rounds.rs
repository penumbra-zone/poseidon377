use crate::input::InputParameters;
use ark_ff::BigInteger;
use ark_std::cmp::{Ordering, PartialOrd};
use poseidon_parameters::v1::{Alpha, RoundNumbers};

/// Generate round numbers.
pub fn generate<T: BigInteger>(input: &InputParameters<T>, alpha: &Alpha) -> RoundNumbers {
    let mut choice: Option<RoundNumbers> = None;
    let mut cost = usize::MAX;
    let mut cost_rf = usize::MAX;

    // Loop through choices of r_F, r_P
    for r_P in 1..400 {
        for r_F in 4..100 {
            let mut candidate = RoundNumbersBuilder(RoundNumbers { r_F, r_P });
            if !candidate.is_secure(input, alpha) {
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
    fn is_secure<T: BigInteger>(&self, input: &InputParameters<T>, alpha: &Alpha) -> bool {
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
                if self.0.total() <= algebraic_attack_grobner_basis(input, alpha) {
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
                    <= algebraic_attack_grobner_basis(input, alpha)
                {
                    return false;
                }
            }
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

/// Number of total rounds to defend against Grobner basis attacks.
///
/// These are described in Section 5.5.2 of the paper.
///
/// We use the first two conditions described in Section C.2.2,
/// eliding the third since if the first condition is satisfied, then
/// the third will be also.
fn algebraic_attack_grobner_basis<T: BigInteger>(
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
