#![allow(non_snake_case)]
use std::cmp::{Ordering, PartialOrd};

use ark_ff::BigInteger;

use super::{Alpha, InputParameters};

pub(super) struct RoundNumbers {
    /// Number of partial rounds.
    r_P: usize,
    /// Number of full rounds.
    r_F: usize,
}

impl RoundNumbers {
    pub fn new<P: BigInteger>(input: &InputParameters<P>) -> Self {
        let mut choice: Option<RoundNumbers> = None;
        let mut cost = usize::MAX;

        // Loop through choices of r_F, r_P
        for r_P in 1..500 {
            for r_F in 6..100 {
                let mut candidate = RoundNumbers { r_F, r_P };
                if !candidate.is_secure(input) {
                    continue;
                }

                candidate.apply_security_margin();
                let candidate_cost = candidate.sbox_count(input.t);
                if candidate_cost < cost {
                    cost = candidate_cost;
                    choice = Some(candidate);
                }
            }
        }

        choice.unwrap()
    }

    /// Determine whether this `RoundNumbers` choice is secure given all known attacks.
    fn is_secure<P: BigInteger>(&self, input: &InputParameters<P>) -> bool {
        // Check if the number of full rounds are sufficient.
        if RoundNumbers::statistical_attack_full_rounds(input) > self.r_F {
            return false;
        }

        // Check the total number of rounds is sufficient.
        if RoundNumbers::algebraic_attack_interpolation(input) >= self.total() {
            return false;
        }
        if RoundNumbers::algebraic_attack_grobner_basis(input) >= self.total() {
            return false;
        }

        true
    }

    /// Get the number of SBoxes for these `RoundNumbers` on a permutation of width `t`.
    fn sbox_count(&self, t: usize) -> usize {
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
    fn statistical_attack_full_rounds<P: BigInteger>(input: &InputParameters<P>) -> usize {
        let r_F = 0usize;

        // C is defined in Section 5.5.1, p.10.
        let C: f64;
        match input.alpha {
            Alpha::Inverse => C = 2.0,
            Alpha::Exponent(exp) => C = (exp as f64 - 1.0).log2(),
        };

        // Statistical attacks require at least 6 full rounds.
        // Differential/Linear Distinguishers.
        // Ref: Section 5.5.1.
        let r_F: usize;
        if input.M as f64 <= ((input.n as f64 - C).floor() * (input.t as f64 + 1.0)) {
            r_F = 6;
        } else {
            r_F = 10;
        }

        r_F
    }

    /// Number of total rounds to defend against interpolation attacks.
    ///
    /// These attacks are described in Section 5.5.2 of the paper.
    /// For positive alpha, we use Eqn 3.
    /// For negative alpha, we use Eqn 4.
    fn algebraic_attack_interpolation<P: BigInteger>(input: &InputParameters<P>) -> usize {
        match input.alpha {
            Alpha::Inverse => todo!("havent done alpha=-1 case"),
            Alpha::Exponent(exp) => {
                let min_args = [input.M as f64, input.n as f64];
                return (1f64
                    + (2f64.log(exp as f64) * min_args.iter().min_by(cmp_f64).expect("no NaNs"))
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
    fn algebraic_attack_grobner_basis<P: BigInteger>(input: &InputParameters<P>) -> usize {
        match input.alpha {
            Alpha::Inverse => todo!("havent done alpha=-1 case"),
            Alpha::Exponent(exp) => {
                // First Grobner constraint
                let grobner_1_min_args = [(input.M as f64 / 3.0), (input.n as f64 / 2.0)];
                let grobner_1 = 2f64.log(exp as f64)
                    * grobner_1_min_args.iter().min_by(cmp_f64).expect("no NaNs");
                // Second Grobner constraint
                let grobner_2_min_args = [
                    2f64.log(exp as f64) * input.M as f64 / (input.t as f64 + 1.0),
                    2f64.log(exp as f64) * input.n as f64 / 2.0,
                ];
                let grobner_2 = (input.t - 1) as f64
                    + grobner_2_min_args.iter().min_by(cmp_f64).expect("no NaNs");

                // Return the _most_ strict Grobner basis constraint.
                // This is the total round requirement.
                let grobner_values = [grobner_1, grobner_2];
                return grobner_values
                    .iter()
                    .max_by(cmp_f64)
                    .expect("no NaNs")
                    .ceil() as usize;
            }
        };
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
