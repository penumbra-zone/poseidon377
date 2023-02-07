use std::fmt::Display;

use ark_ff::PrimeField;
use ark_std::vec::Vec;
use num::BigUint;
use poseidon_parameters::{
    Alpha, ArcMatrix, Matrix, MatrixOperations, MdsMatrix, OptimizedArcMatrix,
    OptimizedMdsMatrices, PoseidonParameters, SquareMatrix,
};

use crate::generate;

/// Create parameter code.
pub fn compile<F: PrimeField>(
    M: usize,
    t_values: Vec<usize>,
    p: F::BigInt,
    allow_inverse: bool,
) -> String {
    let mut params_code = "use ark_ff::PrimeField;\n
use poseidon_parameters::{Alpha, ArcMatrix, RoundNumbers, SquareMatrix, Matrix, MdsMatrix, OptimizedArcMatrix, OptimizedMdsMatrices, PoseidonParameters, MatrixOperations};\n\n"
        .to_string();

    for t in t_values {
        let params = generate::<F>(M, t, p, allow_inverse);
        params_code.push_str(&format!("{}", DisplayablePoseidonParameters(&params))[..]);
    }

    params_code
}

struct DisplayablePoseidonParameters<'a, F: PrimeField>(&'a PoseidonParameters<F>);
impl<F: PrimeField> Display for DisplayablePoseidonParameters<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let this = self.0;

        let capacity = 1;
        let rate = this.t - capacity;

        let rounds = this.rounds;

        let r_P = rounds.partial();
        let r_F = rounds.full();

        let arc = &this.arc;
        let mds = &this.mds;
        let alpha = this.alpha;
        let optimized_mds = &this.optimized_mds;
        let optimized_arc = &this.optimized_arc;

        write!(
            f,
            r"/// Parameters for the rate-{rate} instance of Poseidon.
pub fn rate_{rate}<F: PrimeField>() -> PoseidonParameters<F> {{
    PoseidonParameters {{
        M: {},
        t: {},
        arc: {},
        mds: {},
        alpha: {},
        rounds: RoundNumbers {{r_P: {r_P}, r_F: {r_F}}},
        optimized_mds: {},
        optimized_arc: {},
    }}
}}
",
            this.M,
            this.t,
            DisplayableArcMatrix(arc),
            DisplayableMdsMatrix(mds),
            DisplayableAlpha(alpha),
            DisplayableOptimizedMdsMatrices(optimized_mds),
            DisplayableOptimizedArcMatrix(optimized_arc),
        )
    }
}

struct DisplayableAlpha(Alpha);
impl Display for DisplayableAlpha {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Alpha::Exponent(exp) => write!(f, "Alpha::Exponent({exp})"),
            Alpha::Inverse => write!(f, "Alpha::Inverse"),
        }
    }
}

struct DisplayableMatrix<'a, F: PrimeField>(&'a Matrix<F>);
impl<F: PrimeField> Display for DisplayableMatrix<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.0.n_rows();
        let n_cols = self.0.n_cols();
        let elements = serialize_slice_f(&self.0.elements().to_vec());
        write!(f, r"Matrix::new({n_rows}, {n_cols}, {elements})",)
    }
}

struct DisplayableSquareMatrix<'a, F: PrimeField>(&'a SquareMatrix<F>);
impl<F: PrimeField> Display for DisplayableSquareMatrix<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.0.n_rows();
        let n_cols = self.0.n_cols();
        let elements = serialize_slice_f(&self.0.elements().to_vec());
        write!(f, r"SquareMatrix::new({n_rows}, {n_cols}, {elements})",)
    }
}

fn serialize_slice_matrix_f<F: PrimeField>(elements: &[Matrix<F>]) -> String {
    let mut new_str = "vec![".to_string();
    for elem in elements {
        new_str.push_str(&format!("{}, ", DisplayableMatrix(elem)).to_string());
    }
    // Remove the trailing ", "
    new_str.pop();
    new_str.pop();
    new_str.push(']');
    new_str
}

fn serialize_slice_f<F: PrimeField>(elements: &[F]) -> String {
    let mut new_str = "vec![".to_string();
    for elem in elements {
        // We use the BigUint type here since the Display of the field element
        // is not in decimal: see https://github.com/arkworks-rs/algebra/issues/320
        let elem_bigint: BigUint = (*elem).into();
        new_str.push_str("F::from_str(\"");
        new_str.push_str(&format!("{}", elem_bigint).to_string());
        new_str.push_str("\").map_err(|_| ()).unwrap(), ");
    }
    // Remove the trailing ", "
    new_str.pop();
    new_str.pop();
    new_str.push(']');
    new_str
}

fn serialize_slice_of_vecs_f<F: PrimeField>(elements: &[Vec<F>]) -> String {
    let mut new_str = "vec![".to_string();
    for r in elements {
        for c in r {
            let elem_bigint: BigUint = (*c).into();
            new_str.push_str("F::from_str(\"");
            new_str.push_str(&format!("{}", elem_bigint).to_string());
            new_str.push_str("\").map_err(|_| ()).unwrap(), ");
        }
    }
    // Remove the trailing ", "
    new_str.pop();
    new_str.pop();
    new_str.push(']');
    new_str
}

fn serialize_f<F: PrimeField>(single_element: &F) -> String {
    let mut new_str = "F::from_str(\"".to_string();
    let elem_bigint: BigUint = (*single_element).into();
    new_str.push_str(&format!("{}", elem_bigint));
    new_str.push_str("\").map_err(|_| ()).unwrap()");
    new_str
}

struct DisplayableOptimizedMdsMatrices<'a, F: PrimeField>(&'a OptimizedMdsMatrices<F>);
impl<F: PrimeField> Display for DisplayableOptimizedMdsMatrices<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let this = self.0;

        let M_hat = &this.M_hat;
        let v = &this.v;
        let w = &this.w;
        let M_prime = &this.M_prime;
        let M_doubleprime = &this.M_doubleprime;
        let M_inverse = &this.M_inverse;
        let M_hat_inverse = &this.M_hat_inverse;
        let M_i = &this.M_i;
        let v_collection = &this.v_collection;
        let w_hat_collection = &this.w_hat_collection;

        write!(
            f,
            r"OptimizedMdsMatrices {{
                M_hat: {},
                v: {},
                w: {},
                M_prime: {},
                M_doubleprime: {},
                M_inverse: {},
                M_hat_inverse: {},
                M_00: {},
                M_i: {},
                v_collection: {},
                w_hat_collection: {},
            }}",
            DisplayableSquareMatrix(M_hat),
            DisplayableMatrix(v),
            DisplayableMatrix(w),
            DisplayableSquareMatrix(M_prime),
            DisplayableSquareMatrix(M_doubleprime),
            DisplayableSquareMatrix(M_inverse),
            DisplayableSquareMatrix(M_hat_inverse),
            serialize_f(&this.M_00),
            DisplayableMatrix(M_i),
            serialize_slice_matrix_f(v_collection),
            serialize_slice_matrix_f(w_hat_collection),
        )
    }
}

struct DisplayableMdsMatrix<'a, F: PrimeField>(&'a MdsMatrix<F>);
impl<F: PrimeField> Display for DisplayableMdsMatrix<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mds_elements: Vec<F> = self.0.elements().to_vec();

        let mut mds_str = "MdsMatrix::from_elements(".to_string();
        mds_str.push_str(&serialize_slice_f(&mds_elements));
        mds_str.push(')');
        write!(f, "{}", &mds_str[..])
    }
}

struct DisplayableArcMatrix<'a, F: PrimeField>(&'a ArcMatrix<F>);
impl<F: PrimeField> Display for DisplayableArcMatrix<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.0.n_rows();
        let n_cols = self.0.n_cols();

        let arc = self.0.clone();
        let elements: Vec<Vec<F>> = arc.into();

        let mut arc_str = "ArcMatrix::new(".to_string();
        arc_str.push_str(&n_rows.to_string());
        arc_str.push_str(", ");
        arc_str.push_str(&n_cols.to_string());
        arc_str.push_str(r", ");
        arc_str.push_str(&serialize_slice_of_vecs_f(&elements));
        arc_str.push(')');
        write!(f, "{}", &arc_str[..])
    }
}

fn optimized_arc_matrix_to_vec_vec<F: PrimeField>(matrix: &OptimizedArcMatrix<F>) -> Vec<Vec<F>> {
    let mut rows = Vec::<Vec<F>>::new();
    let arc: &ArcMatrix<F> = &matrix.0;

    for i in 0..arc.n_rows() {
        let mut row = Vec::new();
        for j in 0..arc.n_cols() {
            row.push(arc.get_element(i, j));
        }
        rows.push(row);
    }
    rows
}

struct DisplayableOptimizedArcMatrix<'a, F: PrimeField>(&'a OptimizedArcMatrix<F>);
impl<F: PrimeField> Display for DisplayableOptimizedArcMatrix<'_, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.0 .0.n_rows();
        let n_cols = self.0 .0.n_cols();
        let elements = optimized_arc_matrix_to_vec_vec(self.0);

        let mut arc_str = "OptimizedArcMatrix::new(".to_string();
        arc_str.push_str(&n_rows.to_string());
        arc_str.push_str(", ");
        arc_str.push_str(&n_cols.to_string());
        arc_str.push_str(r", ");
        arc_str.push_str(&serialize_slice_of_vecs_f(&elements));
        arc_str.push(')');
        write!(f, "{}", &arc_str[..])
    }
}
