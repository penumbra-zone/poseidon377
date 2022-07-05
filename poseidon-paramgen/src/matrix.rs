use std::ops::Mul;

use ark_ff::PrimeField;

/// Represents a matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F: PrimeField> {
    pub elements: Vec<F>,
    pub n_cols: usize,
    pub n_rows: usize,
}

impl<F: PrimeField> Matrix<F> {
    pub fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Matrix<F> {
        if elements.len() != n_rows * n_cols {
            panic!("Matrix has insufficient elements")
        }
        Matrix {
            elements,
            n_cols,
            n_rows,
        }
    }

    pub fn get_element(&self, i: usize, j: usize) -> F {
        self.elements[i * self.n_cols + j]
    }

    pub fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.elements[i * self.n_cols + j] = val
    }
}

/// Represents a square matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquareMatrix<F: PrimeField> {
    pub inner: Matrix<F>,
}

impl<F: PrimeField> SquareMatrix<F> {
    pub fn elements(&self) -> &Vec<F> {
        &self.inner.elements
    }

    pub fn rows(&self) -> Vec<&[F]> {
        self.elements().chunks(self.dim()).collect()
    }

    pub fn get_element(&self, i: usize, j: usize) -> F {
        self.inner.get_element(i, j)
    }

    pub fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.inner.set_element(i, j, val)
    }

    pub fn from_vec(elements: Vec<F>) -> Self {
        if (elements.len() as f64).sqrt().fract() != 0.0 {
            panic!("SquareMatrix must be square")
        }

        let dim = (elements.len() as f64).sqrt() as usize;

        SquareMatrix {
            inner: Matrix::new(dim, dim, elements),
        }
    }

    /// Dimension of the dim x dim matrix
    pub fn dim(&self) -> usize {
        self.inner.n_rows
    }

    /// Construct a dim x dim identity matrix
    pub fn identity(dim: usize) -> SquareMatrix<F> {
        let mut m = SquareMatrix::from_vec(vec![F::zero(); dim * dim]);

        // Set diagonals to 1
        for i in 0..dim {
            m.set_element(i, i, F::one());
        }

        m
    }

    /// Take transpose of the matrix
    pub fn transpose(&self) -> SquareMatrix<F> {
        let dim = self.dim();
        let mut transposed_elements = Vec::with_capacity(dim * dim);

        for j in 0..dim {
            for i in 0..dim {
                transposed_elements.push(self.get_element(i, j))
            }
        }
        SquareMatrix::from_vec(transposed_elements)
    }

    pub fn new_2x2(a: F, b: F, c: F, d: F) -> Self {
        SquareMatrix::from_vec(vec![a, b, c, d])
    }

    /// Compute the inverse of the matrix
    pub fn inverse(&self) -> SquareMatrix<F> {
        let identity: SquareMatrix<F> = SquareMatrix::identity(self.dim());
        let matrix_inverse = todo!();

        debug_assert_eq!(self * matrix_inverse, identity);
    }

    /// Compute the matrix determinant
    pub fn determinant(&self) -> F {
        match self.inner.n_cols {
            0 => panic!("matrix has no elements!"),
            1 => self.get_element(0, 0),
            2 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                a11 * a22 - a21 * a12
            }
            3 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);

                a11 * (SquareMatrix::new_2x2(a22, a23, a32, a33).determinant())
                    - a12 * (SquareMatrix::new_2x2(a21, a23, a31, a33).determinant())
                    + a13 * (SquareMatrix::new_2x2(a21, a22, a31, a32).determinant())
            }
            _ => {
                // Unoptimized, but MDS matrices are fairly small, so we do the naive thing
                let mut det = F::zero();
                let mut levi_civita = true;

                for i in 0..self.inner.n_cols {
                    let mut elements: Vec<F> = Vec::new();
                    for k in 0..i {
                        for l in 1..self.inner.n_cols {
                            elements.push(self.get_element(k, l))
                        }
                    }
                    for k in i + 1..self.inner.n_cols {
                        for l in 1..self.inner.n_cols {
                            elements.push(self.get_element(k, l))
                        }
                    }
                    let minor = SquareMatrix::from_vec(elements);
                    if levi_civita {
                        det += self.get_element(i, 0) * minor.determinant();
                    } else {
                        det -= self.get_element(i, 0) * minor.determinant();
                    }
                    levi_civita = !levi_civita;
                }

                det
            }
        }
    }
}

/// Flatten a vec of vecs
fn flatten<T>(nested: Vec<Vec<T>>) -> Vec<T> {
    nested.into_iter().flatten().collect()
}

/// Compute dot product between two vectors
pub fn dot_product<F: PrimeField>(a: Vec<F>, b: Vec<F>) -> F {
    a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
}

impl<F: PrimeField> Mul for SquareMatrix<F> {
    type Output = SquareMatrix<F>;

    // Only multiplying square matrices is infallible
    // since the number of rows in the LHS must be equal to the
    // number of columns in the RHS.
    fn mul(self, rhs: Self) -> Self::Output {
        let rhs_T = rhs.transpose();

        let res: Vec<Vec<F>> = self
            .rows()
            .into_iter()
            .map(|row| {
                // Rows of the transposed matrix are the columns of the original matrix
                rhs_T
                    .rows()
                    .into_iter()
                    .map(|column| dot_product(row.to_vec(), column.to_vec()))
                    .collect()
            })
            .collect();

        SquareMatrix::from_vec(flatten(res))
    }
}

impl<F: PrimeField> Mul for &SquareMatrix<F> {
    type Output = SquareMatrix<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.clone() * rhs.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ed_on_bls12_381::Fq;
    use ark_ff::{One, Zero};

    #[test]
    fn identity_matrix() {
        let identity: SquareMatrix<Fq> = SquareMatrix::identity(2);
        assert_eq!(identity.get_element(0, 0), Fq::one());
        assert_eq!(identity.get_element(0, 1), Fq::zero());
        assert_eq!(identity.get_element(1, 1), Fq::one());
        assert_eq!(identity.get_element(1, 0), Fq::zero());
    }

    #[test]
    fn matmul() {
        let identity: SquareMatrix<Fq> = SquareMatrix::identity(2);

        let matrix_2x2 = SquareMatrix::from_vec(vec![
            Fq::one(),
            Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(4u64),
        ]);

        let res = matrix_2x2 * identity;
        assert_eq!(res.get_element(0, 0), Fq::one());
        assert_eq!(res.get_element(0, 1), Fq::from(2u64));
        assert_eq!(res.get_element(1, 0), Fq::from(3u64));
        assert_eq!(res.get_element(1, 1), Fq::from(4u64));
    }

    #[test]
    fn transpose() {
        let matrix_2x2 = SquareMatrix::from_vec(vec![
            Fq::one(),
            Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(4u64),
        ]);

        let res = matrix_2x2.transpose();
        assert_eq!(res.get_element(0, 0), Fq::one());
        assert_eq!(res.get_element(0, 1), Fq::from(3u64));
        assert_eq!(res.get_element(1, 0), Fq::from(2u64));
        assert_eq!(res.get_element(1, 1), Fq::from(4u64));
    }

    #[test]
    fn create_matrix_from_vec() {
        let matrix_2x2 = SquareMatrix::from_vec(vec![
            Fq::one(),
            Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(4u64),
        ]);
        assert_eq!(matrix_2x2.get_element(0, 0), Fq::one());
        assert_eq!(matrix_2x2.get_element(0, 1), Fq::from(2u64));
        assert_eq!(matrix_2x2.get_element(1, 0), Fq::from(3u64));
        assert_eq!(matrix_2x2.get_element(1, 1), Fq::from(4u64));

        let matrix_2x3 = Matrix::new(
            2,
            3,
            vec![
                Fq::one(),
                Fq::from(2u64),
                Fq::from(3u64),
                Fq::from(4u64),
                Fq::from(5u64),
                Fq::from(6u64),
            ],
        );
        assert_eq!(matrix_2x3.get_element(0, 0), Fq::one());
        assert_eq!(matrix_2x3.get_element(0, 1), Fq::from(2u64));
        assert_eq!(matrix_2x3.get_element(0, 2), Fq::from(3u64));
        assert_eq!(matrix_2x3.get_element(1, 0), Fq::from(4u64));
        assert_eq!(matrix_2x3.get_element(1, 1), Fq::from(5u64));
        assert_eq!(matrix_2x3.get_element(1, 2), Fq::from(6u64));
    }

    #[test]
    fn determinant() {
        let matrix_1x1 = SquareMatrix::from_vec(vec![Fq::one()]);
        assert_eq!(matrix_1x1.determinant(), Fq::one());

        let a = Fq::one();
        let b = Fq::one() + Fq::one();
        let c = Fq::from(3u64);
        let d = Fq::from(4u64);
        let matrix_2x2 = SquareMatrix::from_vec(vec![a, b, c, d]);
        assert_eq!(matrix_2x2.determinant(), -Fq::from(2u64));

        let e = Fq::from(5u64);
        let f = Fq::from(6u64);
        let g = Fq::from(7u64);
        let h = Fq::from(8u64);
        let i = Fq::from(9u64);
        let matrix_3x3 = SquareMatrix::from_vec(vec![a, b, c, d, e, f, g, h, i]);
        assert_eq!(matrix_3x3.determinant(), Fq::from(0u64));

        let elem = Fq::from(10u64);
        let matrix_4x4 = SquareMatrix::from_vec(vec![
            a, b, c, d, e, f, g, h, i, elem, elem, elem, elem, elem, elem, elem,
        ]);
        assert_eq!(matrix_4x4.determinant(), Fq::from(0u64));
    }
}
