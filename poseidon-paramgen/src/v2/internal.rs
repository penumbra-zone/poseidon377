use ark_ff::PrimeField;
use poseidon_parameters::v2::{SquareMatrix, SquareMatrixOperations};

/// Generate internal matrix
///
/// This matrix needs to be invertible, and no arbitrarily long
/// subspace trails should exist.
pub fn generate<F: PrimeField>(t: usize) -> SquareMatrix<F> {
    let M_i: SquareMatrix<F>;

    if t == 2 {
        M_i = SquareMatrix::<F>::from_vec(vec![F::from(2u64), F::one(), F::one(), F::from(3u64)]);
    } else if t == 3 {
        M_i = SquareMatrix::<F>::from_vec(vec![
            F::from(2u64),
            F::one(),
            F::one(),
            F::one(),
            F::from(2u64),
            F::one(),
            F::one(),
            F::one(),
            F::from(3u64),
        ]);
    } else {
        // From Section 5.3 of the Poseidon2 paper, in lieu of implementing
        // the three algorithms defined in Grassi et al. 2020 [0] to check
        // for arbitrarily long subspace trails, we can instead check that the
        // minimal polynomials of the matrices M_i, M_i^2, ..., are irreducible
        // and of maximum degree. If that is true, then no arbitrarily long
        // subspace trails exist. See Proposition 12 and its proof in [0].
        //
        // [0] https://eprint.iacr.org/2020/500
        unimplemented!("internal matrix for t >= 4 not yet implemented")
    }

    // Check the matrix is invertible.
    assert!(M_i.inverse().is_ok());

    M_i
}

#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_377::Fq;
    use poseidon_parameters::v2::MatrixOperations;

    use super::*;

    #[test]
    fn internal_matrix_t_equals_2() {
        let matrix: SquareMatrix<Fq> = generate(2);
        // The off-diagonal elements should be 1. The diagonals are non-zero.

        // Row 0
        assert_eq!(Fq::from(2u64), matrix.get_element(0, 0));
        assert_eq!(Fq::from(1u64), matrix.get_element(0, 1));

        // Row 1
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 0));
        assert_eq!(Fq::from(3u64), matrix.get_element(1, 1));
    }

    #[test]
    fn internal_matrix_t_equals_3() {
        let matrix: SquareMatrix<Fq> = generate(3);
        // The off-diagonal elements should be 1. The diagonals are non-zero.

        // Row 0
        assert_eq!(Fq::from(2u64), matrix.get_element(0, 0));
        assert_eq!(Fq::from(1u64), matrix.get_element(0, 1));
        assert_eq!(Fq::from(1u64), matrix.get_element(0, 2));

        // Row 1
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 0));
        assert_eq!(Fq::from(2u64), matrix.get_element(1, 1));
        assert_eq!(Fq::from(1u64), matrix.get_element(1, 2));

        // Row 2
        assert_eq!(Fq::from(1u64), matrix.get_element(2, 0));
        assert_eq!(Fq::from(1u64), matrix.get_element(2, 1));
        assert_eq!(Fq::from(3u64), matrix.get_element(2, 2));
    }
}
