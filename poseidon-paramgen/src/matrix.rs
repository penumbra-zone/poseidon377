use ark_ff::PrimeField;

/// Represents a matrix over Fp
pub struct Matrix<F: PrimeField>(Vec<Vec<F>>);
