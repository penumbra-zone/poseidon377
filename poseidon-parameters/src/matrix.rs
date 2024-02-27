use core::convert::TryInto;
use core::ops::Mul;

use crate::error::PoseidonParameterError;
use crate::matrix_ops::{dot_product, MatrixOperations, SquareMatrixOperations};
use decaf377::Fq;

/// Represents a matrix over `PrimeField` elements.
///
/// This matrix can be used to represent row or column
/// vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<const N_ROWS: usize, const N_COLS: usize, const N_ELEMENTS: usize> {
    /// Elements of the matrix, stored in a fixed-size array.
    ///
    pub elements: [Fq; N_ELEMENTS],
}

impl<const N_ROWS: usize, const N_COLS: usize, const N_ELEMENTS: usize>
    Matrix<N_ROWS, N_COLS, N_ELEMENTS>
{
    pub fn transpose(&self) -> Matrix<N_COLS, N_ROWS, N_ELEMENTS> {
        let mut transposed_elements = [Fq::default(); N_ELEMENTS];

        let mut index = 0;
        for j in 0..self.n_cols() {
            for i in 0..self.n_rows() {
                transposed_elements[index] = self.get_element(i, j);
                index += 1;
            }
        }
        Matrix::<N_COLS, N_ROWS, N_ELEMENTS>::new(&transposed_elements)
    }
}

impl<const N_ROWS: usize, const N_COLS: usize, const N_ELEMENTS: usize> MatrixOperations
    for Matrix<N_ROWS, N_COLS, N_ELEMENTS>
{
    fn new(elements: &[Fq]) -> Self {
        // Note: We use a third const generic to denote the number of elements in the
        // matrix here due to `generic_const_exprs` being an unstable Rust feature at
        // the time of writing.
        if N_ELEMENTS != N_ROWS * N_COLS {
            panic!("Matrix has an insufficient number of elements")
        }

        let elements: [Fq; N_ELEMENTS] = elements
            .try_into()
            .expect("Matrix has the correct number of elements");

        Self { elements }
    }

    fn elements(&self) -> &[Fq] {
        &self.elements
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.elements[i * N_COLS + j]
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.elements[i * N_COLS + j] = val
    }

    fn rows(&self) -> &[&[Fq]] {
        // self.elements.chunks(self.n_cols()).collect()
        todo!()
    }

    fn n_rows(&self) -> usize {
        N_ROWS
    }

    fn n_cols(&self) -> usize {
        N_COLS
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        let mut new_elements = [Fq::default(); N_ELEMENTS];
        let mut index = 0;
        for i in 0..self.n_rows() {
            for j in 0..self.n_cols() {
                new_elements[index] = self.get_element(i, j) * rhs.get_element(i, j);
                index += 1;
            }
        }

        Ok(Self::new(&new_elements))
    }
}

/// Multiply two matrices
pub fn mat_mul<
    const LHS_N_ROWS: usize,
    const LHS_N_COLS: usize,
    const LHS_N_ELEMENTS: usize,
    const RHS_N_ROWS: usize,
    const RHS_N_COLS: usize,
    const RHS_N_ELEMENTS: usize,
    const RESULT_N_ELEMENTS: usize,
>(
    lhs: &Matrix<LHS_N_ROWS, LHS_N_COLS, LHS_N_ELEMENTS>,
    rhs: &Matrix<RHS_N_ROWS, RHS_N_COLS, RHS_N_ELEMENTS>,
) -> Matrix<LHS_N_ROWS, RHS_N_COLS, RESULT_N_ELEMENTS> {
    let rhs_T = rhs.transpose();

    let mut new_elements = [Fq::default(); RESULT_N_ELEMENTS];

    let mut index = 0;
    for row in lhs.iter_rows() {
        // Rows of the transposed matrix are the columns of the original matrix
        for column in rhs_T.iter_rows() {
            new_elements[index] = dot_product(row, column);
            index += 1;
        }
    }

    Matrix::<LHS_N_ROWS, RHS_N_COLS, RESULT_N_ELEMENTS>::new(&new_elements)
}

/// Multiply scalar by Matrix
impl<const N_ROWS: usize, const N_COLS: usize, const N_ELEMENTS: usize> Mul<Fq>
    for Matrix<N_ROWS, N_COLS, N_ELEMENTS>
{
    type Output = Matrix<N_ROWS, N_COLS, N_ELEMENTS>;

    fn mul(self, rhs: Fq) -> Self::Output {
        let elements = self.elements();
        let mut new_elements = [Fq::default(); N_ELEMENTS];
        for (i, &element) in elements.iter().enumerate() {
            new_elements[i] = element * rhs;
        }
        Self::new(&new_elements)
    }
}

impl<const N_ROWS: usize, const N_COLS: usize, const N_ELEMENTS: usize>
    Matrix<N_ROWS, N_COLS, N_ELEMENTS>
{
    /// Get row vector at a specified row index
    pub fn row_vector(&self, i: usize) -> Matrix<1, N_COLS, N_ELEMENTS> {
        let mut row_elements = [Fq::default(); N_COLS];
        for j in 0..N_COLS {
            row_elements[j] = self.get_element(i, j);
        }
        Matrix::new(&row_elements)
    }
}

impl<const N_ROWS: usize, const N_ELEMENTS: usize> SquareMatrix<N_ROWS, N_ELEMENTS> {
    pub fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }
}

/// Represents a square matrix over `PrimeField` elements
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SquareMatrix<const N_ROWS: usize, const N_ELEMENTS: usize>(
    pub Matrix<N_ROWS, N_ROWS, N_ELEMENTS>,
);

impl<const N_ROWS: usize, const N_ELEMENTS: usize> MatrixOperations
    for SquareMatrix<N_ROWS, N_ELEMENTS>
{
    fn new(elements: &[Fq]) -> Self {
        Self(Matrix::new(elements))
    }

    fn elements(&self) -> &[Fq] {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> &[&[Fq]] {
        todo!()
    }

    fn n_rows(&self) -> usize {
        N_ROWS
    }

    fn n_cols(&self) -> usize {
        // Matrix is square
        N_ROWS
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<const N_ROWS: usize, const N_ELEMENTS: usize> SquareMatrixOperations
    for SquareMatrix<N_ROWS, N_ELEMENTS>
{
    /// Compute the inverse of the matrix
    fn inverse(&self) -> Result<Self, PoseidonParameterError> {
        let identity = Self::identity();

        if self.n_rows() == 1 {
            let elements = [self
                .get_element(0, 0)
                .inverse()
                .expect("inverse of single element must exist for 1x1 matrix")];
            return Ok(Self::new(&elements));
        }

        let determinant = self.determinant();
        if determinant == Fq::zero() {
            return Err(PoseidonParameterError::NoMatrixInverse);
        }

        let minors = self.minors();
        let cofactor_matrix = self.cofactors();
        let signed_minors = minors
            .hadamard_product(&cofactor_matrix)
            .expect("minor and cofactor matrix have correct dimensions");
        let adj = signed_minors.transpose();
        let matrix_inverse = adj * (Fq::one() / determinant);

        debug_assert_eq!(square_mat_mul(self, &matrix_inverse), identity);
        Ok(matrix_inverse)
    }

    /// Construct an identity matrix
    fn identity() -> Self {
        let elements = [Fq::zero(); N_ELEMENTS];
        let mut m = Self::new(&elements);

        // Set diagonals to 1
        for i in 0..N_ROWS {
            m.set_element(i, i, Fq::one());
        }

        m
    }

    /// Compute the (unsigned) minors of this matrix
    fn minors(&self) -> Self {
        match N_ROWS {
            0 => panic!("matrix has no elements!"),
            1 => Self::new(&[self.get_element(0, 0)]),
            2 => {
                let a = self.get_element(0, 0);
                let b = self.get_element(0, 1);
                let c = self.get_element(1, 0);
                let d = self.get_element(1, 1);
                Self::new(&[d, c, b, a])
            }
            3 => minor_matrix::<N_ROWS, 2, N_ELEMENTS, 4>(self),
            4 => minor_matrix::<N_ROWS, 3, N_ELEMENTS, 9>(self),
            5 => minor_matrix::<N_ROWS, 4, N_ELEMENTS, 16>(self),
            6 => minor_matrix::<N_ROWS, 5, N_ELEMENTS, 25>(self),
            7 => minor_matrix::<N_ROWS, 6, N_ELEMENTS, 36>(self),
            8 => minor_matrix::<N_ROWS, 7, N_ELEMENTS, 49>(self),
            _ => {
                unimplemented!("poseidon-parameters only supports square matrices up to 8")
            }
        }
    }

    /// Compute the cofactor matrix, i.e. $C_{ij} = (-1)^{i+j}$
    fn cofactors(&self) -> Self {
        let dim = self.n_rows();
        let mut elements = [Fq::zero(); N_ELEMENTS];

        // TODO: non arkworks Fq::pow
        use crate::StuffThatNeedsToGoInDecaf377;
        let mut index = 0;
        for i in 0..dim {
            for j in 0..dim {
                elements[index] = (-Fq::one()).pow([(i + j) as u64]);
                index += 1;
            }
        }
        Self::new(&elements)
    }

    /// Compute the matrix determinant
    fn determinant(&self) -> Fq {
        match self.n_cols() {
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

                a11 * (Self::new_2x2(a22, a23, a32, a33).determinant())
                    - a12 * (Self::new_2x2(a21, a23, a31, a33).determinant())
                    + a13 * (Self::new_2x2(a21, a22, a31, a32).determinant())
            }
            4 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a14 = self.get_element(0, 3);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a24 = self.get_element(1, 3);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);
                let a34 = self.get_element(2, 3);
                let a41 = self.get_element(3, 0);
                let a42 = self.get_element(3, 1);
                let a43 = self.get_element(3, 2);
                let a44 = self.get_element(3, 3);

                a11 * (SquareMatrix::<3, 9>::new(&[a22, a23, a24, a32, a33, a34, a42, a43, a44])
                    .determinant())
                    - a12
                        * (SquareMatrix::<3, 9>::new(&[
                            a21, a23, a24, a31, a33, a34, a41, a43, a44,
                        ])
                        .determinant())
                    + a13
                        * (SquareMatrix::<3, 9>::new(&[
                            a21, a22, a24, a31, a32, a34, a41, a42, a44,
                        ])
                        .determinant())
                    - a14
                        * (SquareMatrix::<3, 9>::new(&[
                            a21, a22, a23, a31, a32, a33, a41, a42, a43,
                        ])
                        .determinant())
            }
            5 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a14 = self.get_element(0, 3);
                let a15 = self.get_element(0, 4);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a24 = self.get_element(1, 3);
                let a25 = self.get_element(1, 4);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);
                let a34 = self.get_element(2, 3);
                let a35 = self.get_element(2, 4);
                let a41 = self.get_element(3, 0);
                let a42 = self.get_element(3, 1);
                let a43 = self.get_element(3, 2);
                let a44 = self.get_element(3, 3);
                let a45 = self.get_element(3, 4);
                let a51 = self.get_element(4, 0);
                let a52 = self.get_element(4, 1);
                let a53 = self.get_element(4, 2);
                let a54 = self.get_element(4, 3);
                let a55 = self.get_element(4, 4);

                a11 * (SquareMatrix::<4, 16>::new(&[
                    a22, a23, a24, a25, a32, a33, a34, a35, a42, a43, a44, a45, a52, a53, a54, a55,
                ])
                .determinant())
                    - a12
                        * (SquareMatrix::<4, 16>::new(&[
                            a21, a23, a24, a25, a31, a33, a34, a35, a41, a43, a44, a45, a51, a53,
                            a54, a55,
                        ])
                        .determinant())
                    + a13
                        * (SquareMatrix::<4, 16>::new(&[
                            a21, a22, a24, a25, a31, a32, a34, a35, a41, a42, a44, a45, a51, a52,
                            a54, a55,
                        ])
                        .determinant())
                    - a14
                        * (SquareMatrix::<4, 16>::new(&[
                            a21, a22, a23, a25, a31, a32, a33, a35, a41, a42, a43, a45, a51, a52,
                            a53, a55,
                        ])
                        .determinant())
                    + a15
                        * (SquareMatrix::<4, 16>::new(&[
                            a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44, a51, a52,
                            a53, a54,
                        ])
                        .determinant())
            }
            6 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a14 = self.get_element(0, 3);
                let a15 = self.get_element(0, 4);
                let a16 = self.get_element(0, 5);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a24 = self.get_element(1, 3);
                let a25 = self.get_element(1, 4);
                let a26 = self.get_element(1, 5);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);
                let a34 = self.get_element(2, 3);
                let a35 = self.get_element(2, 4);
                let a36 = self.get_element(2, 5);
                let a41 = self.get_element(3, 0);
                let a42 = self.get_element(3, 1);
                let a43 = self.get_element(3, 2);
                let a44 = self.get_element(3, 3);
                let a45 = self.get_element(3, 4);
                let a46 = self.get_element(3, 5);
                let a51 = self.get_element(4, 0);
                let a52 = self.get_element(4, 1);
                let a53 = self.get_element(4, 2);
                let a54 = self.get_element(4, 3);
                let a55 = self.get_element(4, 4);
                let a56 = self.get_element(4, 5);
                let a61 = self.get_element(5, 0);
                let a62 = self.get_element(5, 1);
                let a63 = self.get_element(5, 2);
                let a64 = self.get_element(5, 3);
                let a65 = self.get_element(5, 4);
                let a66 = self.get_element(5, 5);

                a11 * (SquareMatrix::<5, 25>::new(&[
                    a22, a23, a24, a25, a26, a32, a33, a34, a35, a36, a42, a43, a44, a45, a46, a52,
                    a53, a54, a55, a56, a62, a63, a64, a65, a66,
                ])
                .determinant())
                    - a12
                        * (SquareMatrix::<5, 25>::new(&[
                            a21, a23, a24, a25, a26, a31, a33, a34, a35, a36, a41, a43, a44, a45,
                            a46, a51, a53, a54, a55, a56, a61, a63, a64, a65, a66,
                        ])
                        .determinant())
                    + a13
                        * (SquareMatrix::<5, 25>::new(&[
                            a21, a22, a24, a25, a26, a31, a32, a34, a35, a36, a41, a42, a44, a45,
                            a46, a51, a52, a54, a55, a56, a61, a62, a64, a65, a66,
                        ])
                        .determinant())
                    - a14
                        * (SquareMatrix::<5, 25>::new(&[
                            a21, a22, a23, a25, a26, a31, a32, a33, a35, a36, a41, a42, a43, a45,
                            a46, a51, a52, a53, a55, a56, a61, a62, a63, a65, a66,
                        ])
                        .determinant())
                    + a15
                        * (SquareMatrix::<5, 25>::new(&[
                            a21, a22, a23, a24, a26, a31, a32, a33, a34, a36, a41, a42, a43, a44,
                            a46, a51, a52, a53, a54, a56, a61, a62, a63, a64, a66,
                        ])
                        .determinant())
                    - a16
                        * (SquareMatrix::<5, 25>::new(&[
                            a21, a22, a23, a24, a25, a31, a32, a33, a34, a35, a41, a42, a43, a44,
                            a45, a51, a52, a53, a54, a56, a61, a62, a63, a64, a66,
                        ])
                        .determinant())
            }
            7 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a14 = self.get_element(0, 3);
                let a15 = self.get_element(0, 4);
                let a16 = self.get_element(0, 5);
                let a17 = self.get_element(0, 6);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a24 = self.get_element(1, 3);
                let a25 = self.get_element(1, 4);
                let a26 = self.get_element(1, 5);
                let a27 = self.get_element(1, 6);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);
                let a34 = self.get_element(2, 3);
                let a35 = self.get_element(2, 4);
                let a36 = self.get_element(2, 5);
                let a37 = self.get_element(2, 6);
                let a41 = self.get_element(3, 0);
                let a42 = self.get_element(3, 1);
                let a43 = self.get_element(3, 2);
                let a44 = self.get_element(3, 3);
                let a45 = self.get_element(3, 4);
                let a46 = self.get_element(3, 5);
                let a47 = self.get_element(3, 6);
                let a51 = self.get_element(4, 0);
                let a52 = self.get_element(4, 1);
                let a53 = self.get_element(4, 2);
                let a54 = self.get_element(4, 3);
                let a55 = self.get_element(4, 4);
                let a56 = self.get_element(4, 5);
                let a57 = self.get_element(4, 6);
                let a61 = self.get_element(5, 0);
                let a62 = self.get_element(5, 1);
                let a63 = self.get_element(5, 2);
                let a64 = self.get_element(5, 3);
                let a65 = self.get_element(5, 4);
                let a66 = self.get_element(5, 5);
                let a67 = self.get_element(5, 6);
                let a71 = self.get_element(6, 0);
                let a72 = self.get_element(6, 1);
                let a73 = self.get_element(6, 2);
                let a74 = self.get_element(6, 3);
                let a75 = self.get_element(6, 4);
                let a76 = self.get_element(6, 5);
                let a77 = self.get_element(6, 6);

                a11 * (SquareMatrix::<6, 36>::new(&[
                    a22, a23, a24, a25, a26, a27, a32, a33, a34, a35, a36, a37, a42, a43, a44, a45,
                    a46, a47, a52, a53, a54, a55, a56, a57, a62, a63, a64, a65, a66, a67, a72, a73,
                    a74, a75, a76, a77,
                ])
                .determinant())
                    - a12
                        * (SquareMatrix::<6, 36>::new(&[
                            a21, a23, a24, a25, a26, a27, a31, a33, a34, a35, a36, a37, a41, a43,
                            a44, a45, a46, a47, a51, a53, a54, a55, a56, a57, a61, a63, a64, a65,
                            a66, a67, a71, a73, a74, a75, a76, a77,
                        ])
                        .determinant())
                    + a13
                        * (SquareMatrix::<6, 36>::new(&[
                            a21, a22, a24, a25, a26, a27, a31, a32, a34, a35, a36, a37, a41, a42,
                            a44, a45, a46, a47, a51, a52, a54, a55, a56, a57, a61, a62, a64, a65,
                            a66, a67, a71, a72, a74, a75, a76, a77,
                        ])
                        .determinant())
                    - a14
                        * (SquareMatrix::<6, 36>::new(&[
                            a21, a22, a23, a25, a26, a27, a31, a32, a33, a35, a36, a37, a41, a42,
                            a43, a45, a46, a47, a51, a52, a53, a55, a56, a57, a61, a62, a63, a65,
                            a66, a67, a71, a72, a73, a75, a76, a77,
                        ])
                        .determinant())
                    + a15
                        * (SquareMatrix::<6, 36>::new(&[
                            a21, a22, a23, a24, a26, a27, a31, a32, a33, a34, a36, a37, a41, a42,
                            a43, a44, a46, a47, a51, a52, a53, a54, a56, a57, a61, a62, a63, a64,
                            a66, a67, a71, a72, a73, a74, a76, a77,
                        ])
                        .determinant())
                    - a16
                        * (SquareMatrix::<6, 36>::new(&[
                            a21, a22, a23, a24, a25, a27, a31, a32, a33, a34, a35, a37, a41, a42,
                            a43, a44, a45, a47, a51, a52, a53, a54, a55, a57, a61, a62, a63, a64,
                            a65, a67, a71, a72, a73, a74, a75, a77,
                        ])
                        .determinant())
                    + a17
                        * (SquareMatrix::<6, 36>::new(&[
                            a21, a22, a23, a24, a25, a26, a31, a32, a33, a34, a35, a36, a41, a42,
                            a43, a44, a45, a46, a51, a52, a53, a54, a55, a56, a61, a62, a63, a64,
                            a65, a66, a71, a72, a73, a74, a75, a76,
                        ])
                        .determinant())
            }
            8 => {
                let a11 = self.get_element(0, 0);
                let a12 = self.get_element(0, 1);
                let a13 = self.get_element(0, 2);
                let a14 = self.get_element(0, 3);
                let a15 = self.get_element(0, 4);
                let a16 = self.get_element(0, 5);
                let a17 = self.get_element(0, 6);
                let a18 = self.get_element(0, 7);
                let a21 = self.get_element(1, 0);
                let a22 = self.get_element(1, 1);
                let a23 = self.get_element(1, 2);
                let a24 = self.get_element(1, 3);
                let a25 = self.get_element(1, 4);
                let a26 = self.get_element(1, 5);
                let a27 = self.get_element(1, 6);
                let a28 = self.get_element(1, 7);
                let a31 = self.get_element(2, 0);
                let a32 = self.get_element(2, 1);
                let a33 = self.get_element(2, 2);
                let a34 = self.get_element(2, 3);
                let a35 = self.get_element(2, 4);
                let a36 = self.get_element(2, 5);
                let a37 = self.get_element(2, 6);
                let a38 = self.get_element(2, 7);
                let a41 = self.get_element(3, 0);
                let a42 = self.get_element(3, 1);
                let a43 = self.get_element(3, 2);
                let a44 = self.get_element(3, 3);
                let a45 = self.get_element(3, 4);
                let a46 = self.get_element(3, 5);
                let a47 = self.get_element(3, 6);
                let a48 = self.get_element(3, 7);
                let a51 = self.get_element(4, 0);
                let a52 = self.get_element(4, 1);
                let a53 = self.get_element(4, 2);
                let a54 = self.get_element(4, 3);
                let a55 = self.get_element(4, 4);
                let a56 = self.get_element(4, 5);
                let a57 = self.get_element(4, 6);
                let a58 = self.get_element(4, 7);
                let a61 = self.get_element(5, 0);
                let a62 = self.get_element(5, 1);
                let a63 = self.get_element(5, 2);
                let a64 = self.get_element(5, 3);
                let a65 = self.get_element(5, 4);
                let a66 = self.get_element(5, 5);
                let a67 = self.get_element(5, 6);
                let a68 = self.get_element(5, 7);
                let a71 = self.get_element(6, 0);
                let a72 = self.get_element(6, 1);
                let a73 = self.get_element(6, 2);
                let a74 = self.get_element(6, 3);
                let a75 = self.get_element(6, 4);
                let a76 = self.get_element(6, 5);
                let a77 = self.get_element(6, 6);
                let a78 = self.get_element(6, 7);
                let a81 = self.get_element(7, 0);
                let a82 = self.get_element(7, 1);
                let a83 = self.get_element(7, 2);
                let a84 = self.get_element(7, 3);
                let a85 = self.get_element(7, 4);
                let a86 = self.get_element(7, 5);
                let a87 = self.get_element(7, 6);
                let a88 = self.get_element(7, 7);

                a11 * (SquareMatrix::<7, 49>::new(&[
                    a22, a23, a24, a25, a26, a27, a28, a32, a33, a34, a35, a36, a37, a38, a42, a43,
                    a44, a45, a46, a47, a48, a52, a53, a54, a55, a56, a57, a58, a62, a63, a64, a65,
                    a66, a67, a68, a72, a73, a74, a75, a76, a77, a78, a82, a83, a84, a85, a86, a87,
                    a88,
                ])
                .determinant())
                    - a12
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a23, a24, a25, a26, a27, a28, a31, a33, a34, a35, a36, a37, a38,
                            a41, a43, a44, a45, a46, a47, a48, a51, a53, a54, a55, a56, a57, a58,
                            a61, a63, a64, a65, a66, a67, a68, a71, a73, a74, a75, a76, a77, a78,
                            a81, a83, a84, a85, a86, a87, a88,
                        ])
                        .determinant())
                    + a13
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a22, a24, a25, a26, a27, a28, a31, a32, a34, a35, a36, a37, a38,
                            a41, a42, a44, a45, a46, a47, a48, a51, a52, a54, a55, a56, a57, a58,
                            a61, a62, a64, a65, a66, a67, a68, a71, a72, a74, a75, a76, a77, a78,
                            a81, a82, a84, a85, a86, a87, a88,
                        ])
                        .determinant())
                    - a14
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a22, a23, a25, a26, a27, a28, a31, a32, a33, a35, a36, a37, a38,
                            a41, a42, a43, a45, a46, a47, a48, a51, a52, a53, a55, a56, a57, a58,
                            a61, a62, a63, a65, a66, a67, a68, a71, a72, a73, a75, a76, a77, a78,
                            a81, a82, a83, a85, a86, a87, a88,
                        ])
                        .determinant())
                    + a15
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a22, a23, a24, a26, a27, a28, a31, a32, a33, a34, a36, a37, a38,
                            a41, a42, a43, a44, a46, a47, a48, a51, a52, a53, a54, a56, a57, a58,
                            a61, a62, a63, a64, a66, a67, a68, a71, a72, a73, a74, a76, a77, a78,
                            a81, a82, a83, a84, a86, a87, a88,
                        ])
                        .determinant())
                    - a16
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a22, a23, a24, a25, a27, a28, a31, a32, a33, a34, a35, a37, a38,
                            a41, a42, a43, a44, a45, a47, a48, a51, a52, a53, a54, a55, a57, a58,
                            a61, a62, a63, a64, a65, a67, a68, a71, a72, a73, a74, a75, a77, a78,
                            a81, a82, a83, a84, a85, a87, a88,
                        ])
                        .determinant())
                    + a17
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a22, a23, a24, a25, a26, a28, a31, a32, a33, a34, a35, a36, a38,
                            a41, a42, a43, a44, a45, a46, a48, a51, a52, a53, a54, a55, a56, a58,
                            a61, a62, a63, a64, a65, a66, a68, a71, a72, a73, a74, a75, a76, a78,
                            a81, a82, a83, a84, a85, a86, a88,
                        ])
                        .determinant())
                    - a18
                        * (SquareMatrix::<7, 49>::new(&[
                            a21, a22, a23, a24, a25, a26, a27, a31, a32, a33, a34, a35, a36, a37,
                            a41, a42, a43, a44, a45, a46, a47, a51, a52, a53, a54, a55, a56, a57,
                            a61, a62, a63, a64, a65, a66, a67, a71, a72, a73, a74, a75, a76, a77,
                            a81, a82, a83, a84, a85, a86, a87,
                        ])
                        .determinant())
            }
            _ => {
                unimplemented!("poseidon-parameters only supports square matrices up to 8")
            }
        }
    }
}

/// Multiply scalar by SquareMatrix
impl<const N_ROWS: usize, const N_ELEMENTS: usize> Mul<Fq> for SquareMatrix<N_ROWS, N_ELEMENTS> {
    type Output = SquareMatrix<N_ROWS, N_ELEMENTS>;

    fn mul(self, rhs: Fq) -> Self::Output {
        let elements = self.elements();
        let mut new_elements = [Fq::default(); N_ELEMENTS];
        for (i, &element) in elements.iter().enumerate() {
            new_elements[i] = element * rhs;
        }
        Self::new(&new_elements)
    }
}

impl<const N_ROWS: usize, const N_ELEMENTS: usize> SquareMatrix<N_ROWS, N_ELEMENTS> {
    /// Get row vector at a specified row index.
    pub fn row_vector(&self, i: usize) -> Matrix<1, N_ROWS, N_ELEMENTS> {
        self.0.row_vector(i)
    }

    /// Create a 2x2 `SquareMatrix` from four elements.
    pub fn new_2x2(a: Fq, b: Fq, c: Fq, d: Fq) -> SquareMatrix<2, 4> {
        SquareMatrix::<2, 4>::new(&[a, b, c, d])
    }
}

/// Multiply two matrices
pub fn square_mat_mul<
    const LHS_N_ROWS: usize,
    const LHS_N_ELEMENTS: usize,
    const RHS_N_ROWS: usize,
    const RHS_N_ELEMENTS: usize,
    const RESULT_N_ELEMENTS: usize,
>(
    lhs: &SquareMatrix<LHS_N_ROWS, LHS_N_ELEMENTS>,
    rhs: &SquareMatrix<RHS_N_ROWS, RHS_N_ELEMENTS>,
) -> SquareMatrix<LHS_N_ROWS, RESULT_N_ELEMENTS> {
    let rhs_T = rhs.transpose();

    let mut new_elements = [Fq::default(); RESULT_N_ELEMENTS];

    let mut index = 0;
    for row in lhs.iter_rows() {
        // Rows of the transposed matrix are the columns of the original matrix
        for column in rhs_T.iter_rows() {
            new_elements[index] = dot_product(row, column);
            index += 1;
        }
    }

    SquareMatrix::<LHS_N_ROWS, RESULT_N_ELEMENTS>::new(&new_elements)
}

/// Helper function for computing matrix minor
fn minor_matrix<
    const DIM: usize,
    const DIM_MINUS_1: usize,
    const N_ELEMENTS: usize,
    const N_ELEMENTS_DIM_MINUS_1: usize,
>(
    matrix: &SquareMatrix<DIM, N_ELEMENTS>,
) -> SquareMatrix<DIM, N_ELEMENTS> {
    let mut minor_matrix_elements = [Fq::default(); N_ELEMENTS];
    let mut outer_index = 0;
    for i in 0..DIM {
        for j in 0..DIM {
            let mut elements = [Fq::default(); N_ELEMENTS_DIM_MINUS_1];
            let mut index = 0;
            for k in 0..i {
                for l in 0..j {
                    elements[index] = matrix.get_element(k, l);
                    index += 1;
                }
                for l in (j + 1)..DIM {
                    elements[index] = matrix.get_element(k, l);
                    index += 1;
                }
            }
            for k in i + 1..DIM {
                for l in 0..j {
                    elements[index] = matrix.get_element(k, l);
                    index += 1;
                }
                for l in (j + 1)..DIM {
                    elements[index] = matrix.get_element(k, l);
                    index += 1;
                }
            }
            let minor = SquareMatrix::<DIM_MINUS_1, N_ELEMENTS_DIM_MINUS_1>::new(&elements);
            minor_matrix_elements[outer_index] = minor.determinant();
            outer_index += 1;
        }
    }
    SquareMatrix::<DIM, N_ELEMENTS>::new(&minor_matrix_elements)
}
