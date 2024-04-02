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

    /// Create a new matrix from a slice of elements.
    pub const fn new_from_known(elements: [Fq; N_ELEMENTS]) -> Self {
        if N_ELEMENTS != N_ROWS * N_COLS {
            panic!("Matrix has an insufficient number of elements")
        }

        Self { elements }
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

/// Multiply two `Matrix`
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
        if determinant == Fq::from(0u64) {
            return Err(PoseidonParameterError::NoMatrixInverse);
        }

        let minors = self.minors();
        let cofactor_matrix = self.cofactors();
        let signed_minors = minors
            .hadamard_product(&cofactor_matrix)
            .expect("minor and cofactor matrix have correct dimensions");
        let adj = signed_minors.transpose();
        let matrix_inverse = adj * (Fq::from(1u64) / determinant);

        debug_assert_eq!(square_mat_mul(self, &matrix_inverse), identity);
        Ok(matrix_inverse)
    }

    /// Construct an identity matrix
    fn identity() -> Self {
        let elements = [Fq::from(0u64); N_ELEMENTS];
        let mut m = Self::new(&elements);

        // Set diagonals to 1
        for i in 0..N_ROWS {
            m.set_element(i, i, Fq::from(1u64));
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
        let mut elements = [Fq::from(0u64); N_ELEMENTS];

        let mut index = 0;
        for i in 0..dim {
            for j in 0..dim {
                elements[index] = (-Fq::from(1u64)).power([(i + j) as u64]);
                index += 1;
            }
        }
        Self::new(&elements)
    }

    /// Compute the matrix determinant
    fn determinant(&self) -> Fq {
        match N_ROWS {
            0 => panic!("matrix has no elements!"),
            1 => self.get_element(0, 0),
            2 => determinant::<N_ROWS, 1, N_ELEMENTS, 1>(self),
            3 => determinant::<N_ROWS, 2, N_ELEMENTS, 4>(self),
            4 => determinant::<N_ROWS, 3, N_ELEMENTS, 9>(self),
            5 => determinant::<N_ROWS, 4, N_ELEMENTS, 16>(self),
            6 => determinant::<N_ROWS, 5, N_ELEMENTS, 25>(self),
            7 => determinant::<N_ROWS, 6, N_ELEMENTS, 36>(self),
            8 => determinant::<N_ROWS, 7, N_ELEMENTS, 49>(self),
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

    /// Create a new matrix from a slice of elements.
    pub const fn new_from_known(elements: [Fq; N_ELEMENTS]) -> Self {
        Self(Matrix::new_from_known(elements))
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

/// Helper function for computing matrix minors
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

/// Helper function for computing matrix determinant
fn determinant<
    const DIM: usize,
    const DIM_MINUS_1: usize,
    const N_ELEMENTS: usize,
    const N_ELEMENTS_DIM_MINUS_1: usize,
>(
    matrix: &SquareMatrix<DIM, N_ELEMENTS>,
) -> Fq {
    let mut det = Fq::from(0u64);
    let mut levi_civita = true;

    for i in 0..DIM {
        let mut elements = [Fq::default(); N_ELEMENTS_DIM_MINUS_1];
        let mut index = 0;
        for k in 0..i {
            for l in 1..DIM {
                elements[index] = matrix.get_element(k, l);
                index += 1;
            }
        }
        for k in i + 1..DIM {
            for l in 1..DIM {
                elements[index] = matrix.get_element(k, l);
                index += 1;
            }
        }
        let minor = SquareMatrix::<DIM_MINUS_1, N_ELEMENTS_DIM_MINUS_1>::new(&elements);
        if levi_civita {
            det += matrix.get_element(i, 0) * minor.determinant();
        } else {
            det -= matrix.get_element(i, 0) * minor.determinant();
        }
        levi_civita = !levi_civita;
    }
    det
}
