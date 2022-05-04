#![allow(non_snake_case)]
use std::cmp::{Ordering, PartialOrd};

use ark_ff::BigInteger;

use super::InputParameters;

pub(super) struct RoundNumbers {
    /// Number of partial rounds.
    r_P: usize,
    /// Number of full rounds.
    r_F: usize,
}

impl RoundNumbers {
    pub fn new<P: BigInteger>(input: InputParameters<P>) -> Self {
        let r_P = 1usize;
        let r_F = 1usize;

        // We compute the round requirements based on the provided input parameters
        // and all known attacks.
        let r_F_stat = RoundNumbers::statistical_attack_full_rounds(input.clone());
        let r_F_interp = RoundNumbers::algebraic_attack_interpolation(input.clone()) - r_P;
        let r_F_grobner = RoundNumbers::algebraic_attack_grobner_basis(input) - r_P;

        // TODO: r_P

        // TODO: Pick highest r_F, r_P to defend against all known attacks.
        let r_F_choices = [r_F_stat, r_F_interp, r_F_grobner];
        let r_F = *(r_F_choices.iter().max().expect("is not None"));

        // Then, apply the suggested security margin.
        let mut rounds = Self { r_P, r_F };
        rounds.apply_security_margin();
        rounds
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
    fn statistical_attack_full_rounds<P: BigInteger>(input: InputParameters<P>) -> usize {
        let r_F = 0usize;

        // C is defined in Section 5.5.1, p.10.
        let C: f64;
        if input.alpha == -1 {
            C = 2.0
        } else {
            C = (input.alpha as f64 - 1.0).log2();
        }

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
    ///
    /// For positive alpha, we use the equation in Appendix C (Section C.2.1)
    /// which states that the constraint is on the number of total rounds
    /// for alpha > 1.
    ///
    /// For negative alpha, we use Appendix D (TODO).
    fn algebraic_attack_interpolation<P: BigInteger>(input: InputParameters<P>) -> usize {
        if input.alpha > 1 {
            let min_args = [input.M as f64, input.n as f64];
            (1f64
                + (2f64.log(input.alpha as f64)
                    * min_args.iter().min_by(cmp_f64).expect("no NaNs"))
                .ceil()
                + (input.t as f64).log(input.alpha as f64)) as usize
        } else {
            todo!("havent done alpha=-1 case")
        }
    }

    /// Number of total rounds to defend against Grobner basis attacks.
    ///
    /// These are described in Section 5.5.2 of the paper.
    ///
    /// We use the first two conditions described in Section C.2.2,
    /// eliding the third since if the first condition is satisfied, then
    /// the third will be also.
    fn algebraic_attack_grobner_basis<P: BigInteger>(input: InputParameters<P>) -> usize {
        if input.alpha > 1 {
            // First Grobner constraint
            let grobner_1_min_args = [(input.M as f64 / 3.0), (input.n as f64 / 2.0)];
            let grobner_1 = 2f64.log(input.alpha as f64)
                * grobner_1_min_args.iter().min_by(cmp_f64).expect("no NaNs");
            // Second Grobner constraint
            let grobner_2_min_args = [
                2f64.log(input.alpha as f64) * input.M as f64 / (input.t as f64 + 1.0),
                2f64.log(input.alpha as f64) * input.n as f64 / 2.0,
            ];
            let grobner_2 =
                (input.t - 1) as f64 + grobner_2_min_args.iter().min_by(cmp_f64).expect("no NaNs");

            // Return the _most_ strict Grobner basis constraint.
            // This is the total round requirement.
            let grobner_values = [grobner_1, grobner_2];
            return grobner_values
                .iter()
                .max_by(cmp_f64)
                .expect("no NaNs")
                .ceil() as usize;
        } else {
            todo!("havent done alpha=-1 case")
        }
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
