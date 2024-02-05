use core::slice::Chunks;

use ark_ff::{vec::Vec, PrimeField};

use crate::error::PoseidonParameterError;

pub trait MatrixOperations<F> {
    /// Create a new matrix
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self;
    /// Access elements as a vector
    fn elements(&self) -> &Vec<F>;
    /// Get element[i,j]
    fn get_element(&self, i: usize, j: usize) -> F;
    /// Set element[i,j]
    fn set_element(&mut self, i: usize, j: usize, val: F);
    /// Get rows
    fn rows(&self) -> Vec<&[F]>;
    /// Get rows in chunks
    fn iter_rows(&self) -> Chunks<F> {
        self.elements().chunks(self.n_cols())
    }
    /// Number of rows
    fn n_rows(&self) -> usize;
    /// Number of columns
    fn n_cols(&self) -> usize;
    /// Take transpose of the matrix
    fn transpose(&self) -> Self;
    /// Compute Hadamard (element-wise) product
    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized;
}

/// Multiply two matrices
pub fn mat_mul<F: PrimeField, M: MatrixOperations<F>>(
    lhs: &M,
    rhs: &M,
) -> Result<M, PoseidonParameterError> {
    if lhs.n_cols() != rhs.n_rows() {
        return Err(PoseidonParameterError::InvalidMatrixDimensions);
    }

    let rhs_T = rhs.transpose();

    Ok(M::new(
        lhs.n_rows(),
        rhs.n_cols(),
        lhs.iter_rows()
            .flat_map(|row| {
                // Rows of the transposed matrix are the columns of the original matrix
                rhs_T
                    .iter_rows()
                    .map(|column| dot_product(row, column))
                    .collect::<Vec<F>>()
            })
            .collect(),
    ))
}

/// Compute vector dot product
pub fn dot_product<F: PrimeField>(a: &[F], b: &[F]) -> F {
    if a.len() != b.len() {
        panic!("vecs not same len")
    }

    a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
}

pub struct Polynomial<F> {
    /// The coefficients of the polynomial a_0, ..., a_i
    pub coeffs: Vec<F>,
}

impl<F: PrimeField> Polynomial<F> {
    /// Construct a new polynomial
    pub fn new(coeffs: Vec<F>) -> Self {
        Self { coeffs }
    }

    /// Degree of the polynomial
    pub fn max_degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// Evaluate the polynomial at a given point
    pub fn evaluate(&self, x: F) -> F {
        self.coeffs
            .iter()
            .rev()
            .fold(F::zero(), |acc, coeff| acc * x + coeff)
    }

    /// Check if the polynomial is irreducible using Perron's irreducibility criterion.
    pub fn is_irreducible(&self) -> bool {
        // We first need to check the polynomial is monic.
        if self.coeffs.last() != Some(&F::one()) {
            unimplemented!("polynomial is not monic, not sure how to check irreducibility")
        } else {
            // The polynomial is monic, so we can apply Perron's criterion.
            // See https://en.wikipedia.org/wiki/Perron%27s_irreducibility_criterion
            // for more details.
            let n = self.max_degree();
            let mut sum = F::one();
            for i in 0..n - 1 {
                sum += self.coeffs[i];
            }

            match self.coeffs[n - 1] {
                // Condition 1:
                // $\abs{a_{n-1}} > 1 + \abs{a_{n-2}} + ... + \abs{a_0}$
                coeff if coeff > sum => true,
                // Condition 2:
                // $\abs{a_{n-1}} = 1 + \abs{a_{n-2}} + ... + \abs{a_0}$
                // AND
                // f(1) != 0 AND f(-1) != 0
                coeff if coeff == sum => {
                    let f_of_1 = self.evaluate(F::one());
                    let f_of_neg_1 = self.evaluate(-F::one());
                    f_of_1 != F::zero() && f_of_neg_1 != F::zero()
                }
                _ => false,
            }
        }
    }
}

/// Matrix operations that are defined on square matrices.
pub trait SquareMatrixOperations<F> {
    /// Compute the matrix inverse, if it exists
    fn inverse(&self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized;
    /// Construct an n x n identity matrix
    fn identity(n: usize) -> Self;
    /// Compute the matrix of minors
    fn minors(&self) -> Self;
    /// Compute the matrix of cofactors
    fn cofactors(&self) -> Self;
    /// Compute the matrix determinant
    fn determinant(&self) -> F;
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ed_on_bls12_377::Fq;

    #[test]
    fn poly_evaluate() {
        // f(x) = 1 + 2x + 3x^2
        let poly = Polynomial::new(vec![Fq::from(1), Fq::from(2), Fq::from(3)]);
        assert_eq!(poly.max_degree(), 2);
        assert_eq!(poly.evaluate(Fq::from(2)), Fq::from(17));
    }
}
