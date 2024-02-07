use core::slice::Chunks;
use heapless::Vec;

use crate::error::PoseidonParameterError;
use crate::MAX_DIMENSION;
use decaf377::Fq;

pub trait MatrixOperations {
    /// Create a new matrix
    fn new(n_rows: usize, n_cols: usize, elements: Vec<Fq, MAX_DIMENSION>) -> Self;
    /// Access elements as a vector
    fn elements(&self) -> &Vec<Fq, MAX_DIMENSION>;
    /// Get element[i,j]
    fn get_element(&self, i: usize, j: usize) -> Fq;
    /// Set element[i,j]
    fn set_element(&mut self, i: usize, j: usize, val: Fq);
    /// Get rows
    fn rows(&self) -> Vec<&[Fq], MAX_DIMENSION>;
    /// Get rows in chunks
    fn iter_rows(&self) -> Chunks<Fq> {
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
pub fn mat_mul<M: MatrixOperations>(lhs: &M, rhs: &M) -> Result<M, PoseidonParameterError> {
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
                    .collect::<Vec<Fq, MAX_DIMENSION>>()
            })
            .collect(),
    ))
}

/// Compute vector dot product
pub fn dot_product(a: &[Fq], b: &[Fq]) -> Fq {
    if a.len() != b.len() {
        panic!("vecs not same len")
    }

    a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
}

pub struct Polynomial {
    /// The coefficients of the polynomial a_0, ..., a_i
    pub coeffs: Vec<Fq, MAX_DIMENSION>,
}

impl Polynomial {
    /// Construct a new polynomial
    pub fn new(coeffs: Vec<Fq, MAX_DIMENSION>) -> Self {
        Self { coeffs }
    }

    /// Degree of the polynomial
    pub fn max_degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// Evaluate the polynomial at a given point
    pub fn evaluate(&self, x: Fq) -> Fq {
        self.coeffs
            .iter()
            .rev()
            .fold(Fq::zero(), |acc, coeff| acc * x + coeff)
    }

    /// Check if the polynomial is irreducible using Perron's irreducibility criterion.
    pub fn is_irreducible(&self) -> bool {
        // We first need to check the polynomial is monic.
        if self.coeffs.last() != Some(&Fq::one()) {
            unimplemented!("polynomial is not monic, not sure how to check irreducibility")
        } else {
            // The polynomial is monic, so we can apply Perron's criterion.
            // See https://en.wikipedia.org/wiki/Perron%27s_irreducibility_criterion
            // for more details.
            let n = self.max_degree();
            let mut sum = Fq::one();
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
                    let f_of_1 = self.evaluate(Fq::one());
                    let f_of_neg_1 = self.evaluate(-Fq::one());
                    f_of_1 != Fq::zero() && f_of_neg_1 != Fq::zero()
                }
                _ => false,
            }
        }
    }
}

/// Matrix operations that are defined on square matrices.
pub trait SquareMatrixOperations {
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
    fn determinant(&self) -> Fq;
}

#[cfg(test)]
mod tests {
    use super::*;
    use decaf377::Fq;

    #[test]
    fn poly_evaluate() {
        // f(x) = 1 + 2x + 3x^2
        let mut coeffs_vec = Vec::<Fq, MAX_DIMENSION>::new();
        coeffs_vec.push(Fq::from(1u64));
        coeffs_vec.push(Fq::from(2u64));
        coeffs_vec.push(Fq::from(3u64));
        let poly = Polynomial::new(coeffs_vec);
        assert_eq!(poly.max_degree(), 2);
        assert_eq!(poly.evaluate(Fq::from(2u64)), Fq::from(17u64));
    }
}
