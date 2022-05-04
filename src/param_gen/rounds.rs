use super::InputParameters;

pub(super) struct RoundNumbers {
    /// Number of partial rounds.
    r_P: usize,
    /// Number of full rounds.
    r_F: usize,
}

impl RoundNumbers {
    pub fn new(input: InputParameters) -> Self {
        let r_P = 1usize;
        let r_F = 1usize;

        // We compute the round requirements based on the provided input parameters
        // and all known attacks.
        let r_F_stat = RoundNumbers::statistical_attack_full_rounds(input);

        // TODO: Pick highest r_F, r_P to defend against all known attacks.

        // Then, apply the suggested security margin.
        let rounds = Self { r_P, r_F };
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
    fn statistical_attack_full_rounds(input: InputParameters) -> usize {
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
        if input.M <= ((input.p.log2() - C).floor() * (input.t + 1)) {
            r_F = 6;
        } else {
            r_F = 10;
        }

        r_F
    }

    /// Number of rounds to defend against interpolation attacks.
    ///
    /// These are described in Section 5.5.2 of the paper.
    ///
    /// TODO: Is it right that the interpolation attack constraints
    /// are over TOTAL rounds instead of R_F? Eqns imply total rounds for
    /// alpha > 0, but the associated code is for full rounds.
    /// Eqns also imply total rounds for 1/x but the associated code implies
    /// partial rounds.
    ///
    /// Appendix C (p.24) states that the constraint is the total rounds.
    fn algebraic_attack_interpolation(input: InputParameters) -> usize {
        todo!()
    }

    /// Number of rounds to defend against Grobner basis attacks.
    ///
    /// These are described in Section 5.5.2 of the paper.
    fn algebraic_attack_grobner_basis(input: InputParameters) -> usize {
        todo!()
    }

    pub fn total_rounds(&self) -> usize {
        self.r_P + self.r_F
    }

    pub fn partial_rounds(&self) -> usize {
        self.r_P
    }

    pub fn full_rounds(&self) -> usize {
        self.r_F
    }
}

// TODO: tests asserting Table 2 in paper
