//! Matrix operations over field elements.
//!
//! This module provides the essential linear algebra needed for
//! parameter generation.
//!
//! The [`Matrix`] structure can be used to represent any matrix, row vector,
//! or column vector. The [`SquareMatrix`] should be used to represent any
//! matrix that is square.
//!
//! Matrix multiplication is defined using `mat_mul`.

mod base;
mod mult;
mod square;
mod traits;

pub use base::Matrix;
pub use mult::mat_mul;
pub use square::SquareMatrix;
pub use traits::{MatrixOperations, SquareMatrixOperations};

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
    fn square_matmul() {
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
    fn nonsquare_matmul() {
        let test_elements = vec![
            Fq::one(),
            Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(4u64),
            Fq::from(5u64),
            Fq::from(6u64),
        ];
        let matrix_2x3 = Matrix::new(3, 2, test_elements);

        let res = mat_mul(&matrix_2x3, &matrix_2x3);
        assert!(res.is_err());

        let matrix_3x2 = matrix_2x3.transpose();
        let res = mat_mul(&matrix_2x3, &matrix_3x2).expect("is ok");
        assert_eq!(res.get_element(0, 0), Fq::from(5u64));
        assert_eq!(res.get_element(0, 1), Fq::from(11u64));
        assert_eq!(res.get_element(0, 2), Fq::from(17u64));
        assert_eq!(res.get_element(1, 0), Fq::from(11u64));
        assert_eq!(res.get_element(1, 1), Fq::from(25u64));
        assert_eq!(res.get_element(1, 2), Fq::from(39u64));
        assert_eq!(res.get_element(2, 0), Fq::from(17u64));
        assert_eq!(res.get_element(2, 1), Fq::from(39u64));
        assert_eq!(res.get_element(2, 2), Fq::from(61u64));
    }

    #[test]
    fn hadamard_product() {
        let test_elements = vec![
            Fq::one(),
            Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(4u64),
            Fq::from(5u64),
            Fq::from(6u64),
        ];
        let matrix_2x3 = Matrix::new(3, 2, test_elements);

        let res = matrix_2x3.hadamard_product(&matrix_2x3).expect("is ok");
        assert_eq!(res.get_element(0, 0), Fq::from(1u64));
        assert_eq!(res.get_element(0, 1), Fq::from(4u64));
        assert_eq!(res.get_element(1, 0), Fq::from(9u64));
        assert_eq!(res.get_element(1, 1), Fq::from(16u64));
        assert_eq!(res.get_element(2, 0), Fq::from(25u64));
        assert_eq!(res.get_element(2, 1), Fq::from(36u64));
    }

    #[test]
    fn transpose() {
        let matrix_2x3 = Matrix::new(
            3,
            2,
            vec![
                Fq::one(),
                Fq::from(2u64),
                Fq::from(3u64),
                Fq::from(4u64),
                Fq::from(5u64),
                Fq::from(6u64),
            ],
        );
        assert_eq!(matrix_2x3.get_element(0, 1), Fq::from(2u64));
        assert_eq!(matrix_2x3.get_element(1, 0), Fq::from(3u64));
        assert_eq!(matrix_2x3.get_element(1, 1), Fq::from(4u64));
        assert_eq!(matrix_2x3.get_element(2, 0), Fq::from(5u64));
        assert_eq!(matrix_2x3.get_element(2, 1), Fq::from(6u64));
        let res = matrix_2x3.transpose();
        assert_eq!(res.get_element(1, 0), Fq::from(2u64));
        assert_eq!(res.get_element(0, 1), Fq::from(3u64));
        assert_eq!(res.get_element(1, 1), Fq::from(4u64));
        assert_eq!(res.get_element(0, 2), Fq::from(5u64));
        assert_eq!(res.get_element(1, 2), Fq::from(6u64));

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
    fn cofactors() {
        let identity_1x1 = SquareMatrix::identity(1);
        let expected_res = SquareMatrix::from_vec(vec![Fq::one()]);
        assert_eq!(identity_1x1.cofactors(), expected_res);

        let identity_2x2 = SquareMatrix::identity(2);
        let expected_res =
            SquareMatrix::from_vec(vec![Fq::one(), -Fq::one(), -Fq::one(), Fq::one()]);
        assert_eq!(identity_2x2.cofactors(), expected_res);
    }

    #[test]
    fn inverse() {
        let matrix_1x1 = SquareMatrix::from_vec(vec![Fq::from(2u64)]);
        let res = matrix_1x1.inverse();
        assert_eq!(matrix_1x1 * res, SquareMatrix::identity(1));

        let matrix_2x2 = SquareMatrix::from_vec(vec![
            Fq::one(),
            Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(4u64),
        ]);

        let res = matrix_2x2.inverse();
        assert_eq!(matrix_2x2 * res, SquareMatrix::identity(2));

        let identity_3x3: SquareMatrix<Fq> = SquareMatrix::identity(3);
        assert_eq!(identity_3x3, identity_3x3.inverse());

        let matrix_3x3 = SquareMatrix::from_vec(vec![
            Fq::from(3u64),
            Fq::from(0u64),
            Fq::from(2u64),
            Fq::from(2u64),
            Fq::from(0u64),
            -Fq::from(2u64),
            Fq::from(0u64),
            Fq::from(1u64),
            Fq::from(1u64),
        ]);
        let res = matrix_3x3.inverse();
        assert_eq!(matrix_3x3 * res.clone(), SquareMatrix::identity(3));
        let expected_res = SquareMatrix::from_vec(vec![
            Fq::from(2u64),
            Fq::from(2u64),
            Fq::from(0u64),
            -Fq::from(2u64),
            Fq::from(3u64),
            Fq::from(10u64),
            Fq::from(2u64),
            -Fq::from(3u64),
            Fq::from(0u64),
        ]) * (Fq::one() / Fq::from(10u64));
        assert_eq!(res, expected_res);
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
