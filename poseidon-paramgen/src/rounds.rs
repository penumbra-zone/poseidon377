use std::cmp::{Ordering, PartialOrd};

use ark_ff::BigInteger;

use super::{Alpha, InputParameters};

/// `RoundNumbers` based on the equations from the original Poseidon paper.
#[derive(Clone)]
pub struct RoundNumbers {
    /// Number of partial rounds.
    r_P: usize,
    /// Number of full rounds.
    r_F: usize,
}

impl RoundNumbers {
    pub fn new<T: BigInteger>(input: &InputParameters<T>) -> Self {
        let mut choice: Option<RoundNumbers> = None;
        let mut cost = usize::MAX;
        let mut cost_rf = usize::MAX;

        // Loop through choices of r_F, r_P
        for r_P in 1..400 {
            for r_F in 4..100 {
                let mut candidate = RoundNumbers { r_F, r_P };
                if !candidate.is_secure(input) {
                    continue;
                }

                candidate.apply_security_margin();
                let candidate_cost = candidate.sbox_count(input.t);
                // Pick the minimum cost Candidate, and if the cost is tied with another
                // candidate, we switch to the new candidate if the total number of full rounds is lower.
                if (candidate_cost < cost) || ((candidate_cost == cost) && (r_F < cost_rf)) {
                    cost = candidate_cost;
                    cost_rf = r_F;
                    choice = Some(candidate);
                }
            }
        }

        choice.unwrap()
    }

    /// Determine whether this `RoundNumbers` choice is secure given all known attacks.
    fn is_secure<T: BigInteger>(&self, input: &InputParameters<T>) -> bool {
        // Check if the number of full rounds are sufficient.
        if self.r_F < RoundNumbers::statistical_attack_full_rounds(input) {
            return false;
        }

        match input.alpha {
            // For positive alpha, the interpolation and Grobner bounds are on the total
            // number of rounds.
            Alpha::Exponent(_) => {
                if self.total() <= RoundNumbers::algebraic_attack_interpolation(input) {
                    return false;
                }
                if self.total() <= RoundNumbers::algebraic_attack_grobner_basis(input) {
                    return false;
                }
            }
            // For inverse alpha, the interpolation and Grobner bounds are on r_F scaled
            // by the binary log of `t` plus r_P. See Eqn 4.
            Alpha::Inverse => {
                if (self.r_F as f64 * (input.t as f64).log2()).floor() as usize + self.r_P
                    <= RoundNumbers::algebraic_attack_interpolation(input)
                {
                    return false;
                }
                if (self.r_F as f64 * (input.t as f64).log2()).floor() as usize + self.r_P
                    <= RoundNumbers::algebraic_attack_grobner_basis(input)
                {
                    return false;
                }
            }
        }

        true
    }

    /// Get the number of SBoxes for these `RoundNumbers` on a permutation of width `t`.
    pub(crate) fn sbox_count(&self, t: usize) -> usize {
        t * self.r_F + self.r_P
    }

    /// Add suggested security margin of +2 R_F and +7.5% R_P.
    /// Ref: Section 5.4.
    fn apply_security_margin(&mut self) {
        self.r_F += 2;
        self.r_P = (1.075 * (self.r_P as f64)).ceil() as usize;
    }

    /// Number of full rounds required to defend against statistical attacks.
    ///
    /// These are the differential/linear distinguisher attacks described
    /// in Section 5.5.1 of the paper.
    fn statistical_attack_full_rounds<T: BigInteger>(input: &InputParameters<T>) -> usize {
        // C is defined in Section 5.5.1, p.10.
        let C = match input.alpha {
            Alpha::Inverse => 2.0,
            Alpha::Exponent(exp) => (exp as f64 - 1.0).log2(),
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
    fn algebraic_attack_interpolation<T: BigInteger>(input: &InputParameters<T>) -> usize {
        let min_args = [input.M as f64, input.log_2_p];
        match input.alpha {
            Alpha::Inverse => {
                return (((input.t as f64).log2()).ceil()
                    + (0.5 * min_args.iter().min_by(cmp_f64).expect("no NaNs")).ceil())
                    as usize;
            }
            Alpha::Exponent(exp) => {
                return ((2f64.log(exp as f64) * min_args.iter().min_by(cmp_f64).expect("no NaNs"))
                    .ceil()
                    + (input.t as f64).log(exp as f64))
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
    fn algebraic_attack_grobner_basis<T: BigInteger>(input: &InputParameters<T>) -> usize {
        let grobner_1: f64;
        let grobner_2: f64;

        match input.alpha {
            Alpha::Inverse => {
                // First Grobner constraint
                let grobner_1_min_args = [input.M as f64, (input.log_2_p as f64)];
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
                grobner_1 = 2f64.log(exp as f64)
                    * grobner_1_min_args.iter().min_by(cmp_f64).expect("no NaNs");
                // Second Grobner constraint
                let grobner_2_min_args = [
                    2f64.log(exp as f64) * input.M as f64 / (input.t as f64 + 1.0),
                    2f64.log(exp as f64) * input.log_2_p / 2.0,
                ];
                grobner_2 = (input.t - 1) as f64
                    + grobner_2_min_args.iter().min_by(cmp_f64).expect("no NaNs");
            }
        };

        // Return the _most_ strict Grobner basis constraint.
        let grobner_values = [grobner_1, grobner_2];
        return grobner_values
            .iter()
            .max_by(cmp_f64)
            .expect("no NaNs")
            .ceil() as usize;
    }

    pub fn total(&self) -> usize {
        self.r_P + self.r_F
    }

    pub fn partial(&self) -> usize {
        self.r_P
    }

    pub fn full(&self) -> usize {
        self.r_F
    }
}

fn cmp_f64(lhs: &&f64, rhs: &&f64) -> Ordering {
    lhs.partial_cmp(rhs).unwrap()
}
