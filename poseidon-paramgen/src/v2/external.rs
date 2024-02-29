use ark_ff::PrimeField;
use poseidon_parameters::{
    v2::MatrixOperations,
    v2::{SquareMatrix, SquareMatrixOperations},
};

/// Generate external matrix
pub fn generate<F: PrimeField>(t: usize) -> SquareMatrix<F> {
    if t < 4 {
        // For t=[2, 3], we don't need to generate an external matrix
        // because we also use the internal matrix in the full rounds.
        panic!("unexpected state size: should use internal matrix for t < 4")
    } else if t % 4 != 0 {
        // The internal matrix is only defined for t = 4t' where t' is an integer.
        panic!("unexpected state size: internal matrix only defined for t % 4 != 0")
    }

    // For t>= 4, we use the following fixed matrix (Section 5.1, Poseidon2 paper).
    let M4 = SquareMatrix::<F>::from_vec(vec![
        F::from(5u64),
        F::from(7u64),
        F::one(),
        F::from(3u64),
        F::from(4u64),
        F::from(6u64),
        F::one(),
        F::one(),
        F::one(),
        F::from(3u64),
        F::from(5u64),
        F::from(7u64),
        F::one(),
        F::one(),
        F::from(4u64),
        F::from(6u64),
    ]);

    if t == 4 {
        M4
    } else {
        let mut matrix = SquareMatrix::identity(t);
        let d = t / 4;
        for i in 0..d {
            for j in 0..d {
                if i == j {
                    for inner_row in 0..4 {
                        for inner_col in 0..4 {
                            matrix.set_element(
                                i * 4 + inner_row,
                                j * 4 + inner_col,
                                F::from(2u64) * M4.get_element(inner_row, inner_col),
                            )
                        }
                    }
                } else {
                    for inner_row in 0..4 {
                        for inner_col in 0..4 {
                            matrix.set_element(
                                i * 4 + inner_row,
                                j * 4 + inner_col,
                                M4.get_element(inner_row, inner_col),
                            )
                        }
                    }
                }
            }
        }
        matrix
    }
}

#[cfg(test)]
mod tests {
    use decaf377::Fq;

    use super::*;

    #[test]
    fn external_matrix_t_equals_4() {
        let matrix: SquareMatrix<Fq> = generate(4);
        // If t=4, the matrix should simply be the fixed M4, unmodified.

        // Row 0
        assert_eq!(Fq::from(5u64), matrix.get_element(0, 0));
        assert_eq!(Fq::from(7u64), matrix.get_element(0, 1));
        assert_eq!(Fq::from(1u64), matrix.get_element(0, 2));
        assert_eq!(Fq::from(3u64), matrix.get_element(0, 3));

        // Row 1
        assert_eq!(Fq::from(4u64), matrix.get_element(1, 0));
        assert_eq!(Fq::from(6u64), matrix.get_element(1, 1));
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 2));
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 3));

        // Row 2
        assert_eq!(Fq::from(1u64), matrix.get_element(2, 0));
        assert_eq!(Fq::from(3u64), matrix.get_element(2, 1));
        assert_eq!(Fq::from(5u64), matrix.get_element(2, 2));
        assert_eq!(Fq::from(7u64), matrix.get_element(2, 3));

        // Row 3
        assert_eq!(Fq::from(1u64), matrix.get_element(3, 0));
        assert_eq!(Fq::from(1u64), matrix.get_element(3, 1));
        assert_eq!(Fq::from(4u64), matrix.get_element(3, 2));
        assert_eq!(Fq::from(6u64), matrix.get_element(3, 3));
    }

    #[test]
    fn external_matrix_t_equals_8() {
        let matrix: SquareMatrix<Fq> = generate(8);

        // Row 0
        assert_eq!(Fq::from(10u64), matrix.get_element(0, 0));
        assert_eq!(Fq::from(14u64), matrix.get_element(0, 1));
        assert_eq!(Fq::from(2u64), matrix.get_element(0, 2));
        assert_eq!(Fq::from(6u64), matrix.get_element(0, 3));
        assert_eq!(Fq::from(5u64), matrix.get_element(0, 4));
        assert_eq!(Fq::from(7u64), matrix.get_element(0, 5));
        assert_eq!(Fq::from(1u64), matrix.get_element(0, 6));
        assert_eq!(Fq::from(3u64), matrix.get_element(0, 7));

        // Row 1
        assert_eq!(Fq::from(8u64), matrix.get_element(1, 0));
        assert_eq!(Fq::from(12u64), matrix.get_element(1, 1));
        assert_eq!(Fq::from(2u64), matrix.get_element(1, 2));
        assert_eq!(Fq::from(2u64), matrix.get_element(1, 3));
        assert_eq!(Fq::from(4u64), matrix.get_element(1, 4));
        assert_eq!(Fq::from(6u64), matrix.get_element(1, 5));
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 6));
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 7));

        // Row 2
        assert_eq!(Fq::from(2u64), matrix.get_element(2, 0));
        assert_eq!(Fq::from(6u64), matrix.get_element(2, 1));
        assert_eq!(Fq::from(10u64), matrix.get_element(2, 2));
        assert_eq!(Fq::from(14u64), matrix.get_element(2, 3));
        assert_eq!(Fq::from(1u64), matrix.get_element(2, 4));
        assert_eq!(Fq::from(3u64), matrix.get_element(2, 5));
        assert_eq!(Fq::from(5u64), matrix.get_element(2, 6));
        assert_eq!(Fq::from(7u64), matrix.get_element(2, 7));

        // Row 3
        assert_eq!(Fq::from(2u64), matrix.get_element(3, 0));
        assert_eq!(Fq::from(2u64), matrix.get_element(3, 1));
        assert_eq!(Fq::from(8u64), matrix.get_element(3, 2));
        assert_eq!(Fq::from(12u64), matrix.get_element(3, 3));
        assert_eq!(Fq::from(1u64), matrix.get_element(3, 4));
        assert_eq!(Fq::from(1u64), matrix.get_element(3, 5));
        assert_eq!(Fq::from(4u64), matrix.get_element(3, 6));
        assert_eq!(Fq::from(6u64), matrix.get_element(3, 7));

        // Row 4
        assert_eq!(Fq::from(5u64), matrix.get_element(4, 0));
        assert_eq!(Fq::from(7u64), matrix.get_element(4, 1));
        assert_eq!(Fq::from(1u64), matrix.get_element(4, 2));
        assert_eq!(Fq::from(3u64), matrix.get_element(4, 3));
        assert_eq!(Fq::from(10u64), matrix.get_element(4, 4));
        assert_eq!(Fq::from(14u64), matrix.get_element(4, 5));
        assert_eq!(Fq::from(2u64), matrix.get_element(4, 6));
        assert_eq!(Fq::from(6u64), matrix.get_element(4, 7));

        // Row 5
        assert_eq!(Fq::from(4u64), matrix.get_element(5, 0));
        assert_eq!(Fq::from(6u64), matrix.get_element(5, 1));
        assert_eq!(Fq::from(1u64), matrix.get_element(5, 2));
        assert_eq!(Fq::from(1u64), matrix.get_element(5, 3));
        assert_eq!(Fq::from(8u64), matrix.get_element(5, 4));
        assert_eq!(Fq::from(12u64), matrix.get_element(5, 5));
        assert_eq!(Fq::from(2u64), matrix.get_element(5, 6));
        assert_eq!(Fq::from(2u64), matrix.get_element(5, 7));

        // Row 6
        assert_eq!(Fq::from(1u64), matrix.get_element(6, 0));
        assert_eq!(Fq::from(3u64), matrix.get_element(6, 1));
        assert_eq!(Fq::from(5u64), matrix.get_element(6, 2));
        assert_eq!(Fq::from(7u64), matrix.get_element(6, 3));
        assert_eq!(Fq::from(2u64), matrix.get_element(6, 4));
        assert_eq!(Fq::from(6u64), matrix.get_element(6, 5));
        assert_eq!(Fq::from(10u64), matrix.get_element(6, 6));
        assert_eq!(Fq::from(14u64), matrix.get_element(6, 7));

        // Row 7
        assert_eq!(Fq::from(1u64), matrix.get_element(7, 0));
        assert_eq!(Fq::from(1u64), matrix.get_element(7, 1));
        assert_eq!(Fq::from(4u64), matrix.get_element(7, 2));
        assert_eq!(Fq::from(6u64), matrix.get_element(7, 3));
        assert_eq!(Fq::from(2u64), matrix.get_element(7, 4));
        assert_eq!(Fq::from(2u64), matrix.get_element(7, 5));
        assert_eq!(Fq::from(8u64), matrix.get_element(7, 6));
        assert_eq!(Fq::from(12u64), matrix.get_element(7, 7));
    }
}
