use std::fmt::Display;

use ark_ff::PrimeField;
use ark_std::vec::Vec;
use num::BigUint;
use poseidon_parameters::{Alpha, Matrix, MatrixOperations};

use crate::{
    matrix::{MatrixWrapper, SquareMatrixWrapper},
    mds::MdsMatrixWrapper,
    round_constants::ArcMatrixWrapper,
    AlphaWrapper, OptimizedArcMatrixWrapper, OptimizedMdsMatricesWrapper,
    PoseidonParametersWrapper,
};

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
        let params: PoseidonParametersWrapper<F> =
            PoseidonParametersWrapper::<F>::new(M, t, p, allow_inverse).into();
        params_code.push_str(&format!("{}", params)[..]);
    }

    params_code
}

impl<F: PrimeField> Display for PoseidonParametersWrapper<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let this = &self.0;

        let capacity = 1;
        let rate = this.t - capacity;

        let rounds = this.rounds;

        let r_P = rounds.partial();
        let r_F = rounds.full();

        let arc = ArcMatrixWrapper(this.arc.clone());
        let mds = MdsMatrixWrapper(this.mds.clone());
        let alpha = AlphaWrapper(this.alpha);
        let optimized_mds = OptimizedMdsMatricesWrapper(this.optimized_mds.clone());
        let optimized_arc = OptimizedArcMatrixWrapper(this.optimized_arc.clone());

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
            this.M, this.t, arc, mds, alpha, optimized_mds, optimized_arc
        )
    }
}

impl Display for AlphaWrapper {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            Alpha::Exponent(exp) => write!(f, "Alpha::Exponent({exp})"),
            Alpha::Inverse => write!(f, "Alpha::Inverse"),
        }
    }
}

impl<F: PrimeField> Display for MatrixWrapper<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.n_rows();
        let n_cols = self.n_cols();
        let elements = serialize_vec_f(self.elements().to_vec());
        write!(f, r"Matrix::new({n_rows}, {n_cols}, {elements})",)
    }
}

impl<F: PrimeField> Display for SquareMatrixWrapper<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.0.n_rows();
        let n_cols = self.0.n_cols();
        let elements = serialize_vec_f(self.0.elements().to_vec());
        write!(f, r"SquareMatrix::new({n_rows}, {n_cols}, {elements})",)
    }
}

fn serialize_vec_matrix_f<F: PrimeField>(elements: Vec<Matrix<F>>) -> String {
    let mut new_str = "vec![".to_string();
    for elem in elements {
        let elem: MatrixWrapper<F> = elem.into();
        new_str.push_str(&format!("{}, ", elem).to_string());
    }
    // Remove the trailing ", "
    new_str.pop();
    new_str.pop();
    new_str.push(']');
    new_str
}

fn serialize_vec_f<F: PrimeField>(elements: Vec<F>) -> String {
    let mut new_str = "vec![".to_string();
    for elem in elements {
        // We use the BigUint type here since the Display of the field element
        // is not in decimal: see https://github.com/arkworks-rs/algebra/issues/320
        let elem_bigint: BigUint = elem.into();
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

fn serialize_vec_of_vecs_f<F: PrimeField>(elements: Vec<Vec<F>>) -> String {
    let mut new_str = "vec![".to_string();
    for r in elements {
        for c in r {
            let elem_bigint: BigUint = c.into();
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

fn serialize_f<F: PrimeField>(single_element: F) -> String {
    let mut new_str = "F::from_str(\"".to_string();
    let elem_bigint: BigUint = single_element.into();
    new_str.push_str(&format!("{}", elem_bigint));
    new_str.push_str("\").map_err(|_| ()).unwrap()");
    new_str
}

impl<F: PrimeField> Display for OptimizedMdsMatricesWrapper<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let this = self.0.clone();

        let M_hat = SquareMatrixWrapper(this.M_hat);
        let v: MatrixWrapper<F> = this.v.into();
        let w: MatrixWrapper<F> = this.w.into();
        let M_prime = SquareMatrixWrapper(this.M_prime);
        let M_doubleprime = SquareMatrixWrapper(this.M_doubleprime);
        let M_inverse = SquareMatrixWrapper(this.M_inverse);
        let M_hat_inverse = SquareMatrixWrapper(this.M_hat_inverse);
        let M_i: MatrixWrapper<F> = this.M_i.into();
        let v_collection = this.v_collection.clone();
        let w_hat_collection = this.w_hat_collection.clone();

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
            M_hat,
            v,
            w,
            M_prime,
            M_doubleprime,
            M_inverse,
            M_hat_inverse,
            serialize_f(this.M_00),
            M_i,
            serialize_vec_matrix_f(v_collection),
            serialize_vec_matrix_f(w_hat_collection),
        )
    }
}

impl<F: PrimeField> Display for MdsMatrixWrapper<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mds_elements: Vec<F> = self.into();

        let mut mds_str = "MdsMatrix::from_elements(".to_string();
        mds_str.push_str(&serialize_vec_f(mds_elements));
        mds_str.push(')');
        write!(f, "{}", &mds_str[..])
    }
}

impl<F: PrimeField> Display for ArcMatrixWrapper<F> {
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
        arc_str.push_str(&serialize_vec_of_vecs_f(elements));
        arc_str.push(')');
        write!(f, "{}", &arc_str[..])
    }
}

impl<F: PrimeField> Display for OptimizedArcMatrixWrapper<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_rows = self.0 .0.n_rows();
        let n_cols = self.0 .0.n_cols();
        let elements: Vec<Vec<F>> = self.into();

        let mut arc_str = "OptimizedArcMatrix::new(".to_string();
        arc_str.push_str(&n_rows.to_string());
        arc_str.push_str(", ");
        arc_str.push_str(&n_cols.to_string());
        arc_str.push_str(r", ");
        arc_str.push_str(&serialize_vec_of_vecs_f(elements));
        arc_str.push(')');
        write!(f, "{}", &arc_str[..])
    }
}
