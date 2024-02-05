use crate::{error::PoseidonParameterError, matrix::Matrix, matrix_ops::MatrixOperations};
use decaf377::Fq;

/// Represents an matrix of round constants.
///
/// Arc stands for `AddRoundConstant` which is the
/// step in the permutation where this matrix is used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArcMatrix(pub Matrix);

impl MatrixOperations for ArcMatrix {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<Fq>) -> Self {
        Self(Matrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<Fq> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[Fq]> {
        self.0.rows()
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl From<ArcMatrix> for Vec<Vec<Fq>> {
    fn from(arc: ArcMatrix) -> Self {
        let mut rows = Vec::<Vec<Fq>>::new();
        let m = &arc.0;

        for i in 0..arc.n_rows() {
            let mut row = Vec::new();
            for j in 0..arc.n_cols() {
                row.push(m.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

/// Represents an optimized matrix of round constants.
///
/// This modifies the partial rounds in the middle of the permutation,
/// wherein you add constants _first_ before iterating through the partial
/// rounds.
///
/// This method follows `calc_equivalent_constants` from Appendix B's
/// `poseidonperm_x3_64_24_optimized.sage`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedArcMatrix(pub ArcMatrix);

impl MatrixOperations for OptimizedArcMatrix {
    /// Create a `OptimizedArcMatrix` from its elements.
    fn new(n_rows: usize, n_cols: usize, elements: Vec<Fq>) -> Self {
        Self(ArcMatrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<Fq> {
        self.0.elements()
    }

    fn get_element(&self, i: usize, j: usize) -> Fq {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: Fq) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[Fq]> {
        self.0.rows()
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self, PoseidonParameterError>
    where
        Self: Sized,
    {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}
