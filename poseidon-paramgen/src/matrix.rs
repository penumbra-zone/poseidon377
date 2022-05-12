use ark_ff::PrimeField;

/// Represents a matrix over Fp
pub struct Matrix<F: PrimeField>(Vec<Vec<F>>);

/// Represents a square matrix over Fp
pub struct SquareMatrix<F: PrimeField> {
    elements: Matrix<F>,
    dim: usize,
}

impl<F: PrimeField> SquareMatrix<F> {
    pub fn new(elements: Vec<Vec<F>>) -> Self {
        let dim = elements.len();

        if dim == 0 {
            panic!("Matrix must have at least one row")
        }

        for row in elements.iter() {
            if row.len() != dim {
                panic!("SquareMatrix must be square!");
            }
        }

        SquareMatrix {
            elements: Matrix(elements),
            dim,
        }
    }
}
