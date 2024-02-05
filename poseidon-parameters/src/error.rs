#[derive(Debug)]
pub enum PoseidonParameterError {
    InvalidMatrixDimensions,
    NoMatrixInverse,
}

impl core::fmt::Display for PoseidonParameterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let msg = match self {
            Self::InvalidMatrixDimensions => "Invalid matrix dimensions",
            Self::NoMatrixInverse => "No matrix inverse",
        };

        msg.fmt(f)
    }
}
