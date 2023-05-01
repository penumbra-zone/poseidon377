/// `RoundNumbers` required for security based on known attacks.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct RoundNumbers {
    /// Number of partial rounds.
    pub r_P: usize,
    /// Number of full rounds.
    pub r_F: usize,
}

impl RoundNumbers {
    /// Number of full rounds.    
    pub fn full(&self) -> usize {
        self.r_F
    }

    /// Number of full rounds as mutable reference.
    pub fn full_mut(&mut self) -> &mut usize {
        &mut self.r_F
    }

    /// Number of partial rounds.    
    pub fn partial(&self) -> usize {
        self.r_P
    }

    /// Number of full rounds as mutable reference.
    pub fn partial_mut(&mut self) -> &mut usize {
        &mut self.r_P
    }

    /// Number of total rounds.
    pub fn total(&self) -> usize {
        self.r_P + self.r_F
    }
}
