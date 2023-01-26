use anyhow::Result;
use ark_ff::{BigInteger, PrimeField};
use ark_std::vec::Vec;
use poseidon_parameters::BasicMatrixOperations;

use crate::{
    matrix::mat_mul, InputParameters, Matrix, MatrixOperations, PoseidonParameters, RoundNumbers,
    SquareMatrix, SquareMatrixOperations,
};

/// Represents an MDS (maximum distance separable) matrix.
pub struct MdsMatrix<T>(pub T);

impl<F> MdsMatrix<&poseidon_parameters::MdsMatrix<F>>
where
    F: PrimeField,
{
    /// Generate the MDS matrix.
    pub fn generate(
        input: &poseidon_parameters::InputParameters<F::BigInt>,
    ) -> poseidon_parameters::MdsMatrix<F> {
        // A t x t MDS matrix only exists if: 2t + 1 <= p
        let two_times_t_bigint: F::BigInt = (2 * input.t as u64).into();
        if two_times_t_bigint > input.p {
            panic!("no MDS matrix exists");
        }

        MdsMatrix::fixed_cauchy_matrix(input)
    }

    /// Instantiate an MDS matrix from a list of elements.
    ///
    /// # Security
    ///
    /// You must ensure this matrix was generated securely,
    /// using the Cauchy method in `fixed_cauchy_matrix` or
    /// using the random subsampling method described in the original
    /// paper.
    pub fn from_elements(elements: Vec<F>) -> poseidon_parameters::MdsMatrix<F> {
        poseidon_parameters::MdsMatrix::<F>(SquareMatrix::from_vec(elements))
    }

    /// Generate a deterministic Cauchy matrix
    ///
    /// The original Poseidon paper describes a method for constructing MDS matrices
    /// from randomly selecting $x_i$, $y_j$ from the field and then constructing each element in the
    /// matrix using $1/(x_i + y_j)$. The resulting Cauchy matrix needs to then be passed to algorithms 1-3
    /// described in
    /// [Grassi, Rechberger, Schofnegger 2020](https://eprint.iacr.org/archive/2020/500/20200702:141143)
    /// in order to determine if infinitely long subspace trails can be constructed for
    /// the Cauchy matrix. If yes, then the MDS matrix must be thrown away, and the process
    /// must begin again for another random choice of $x_i$, $y_j$ until a secure choice is found.
    ///
    /// However, here we use a deterministic method for creating Cauchy matrices that has
    /// been empirically checked to be safe using the three algorithms above over `decaf377` for t=1-100.
    pub fn fixed_cauchy_matrix(
        input: &poseidon_parameters::InputParameters<F::BigInt>,
    ) -> poseidon_parameters::MdsMatrix<F> {
        // We explicitly check for small fields where the deterministic procedure can fail.
        // In these cases, the full algorithms 1-3 should be implemented.
        if input.p.num_bits() < 128 {
            panic!("field too small to use deterministic MDS matrix generation")
        }

        let xs: Vec<F> = (0..input.t as u64).map(F::from).collect();
        let ys: Vec<F> = (input.t as u64..2 * input.t as u64).map(F::from).collect();

        let mut elements = Vec::<F>::with_capacity(input.t);
        for i in 0..input.t {
            for j in 0..input.t {
                // Check x_i + y_j != 0
                assert_ne!(xs[i] + ys[i], F::zero());
                elements.push(F::one() / (xs[i] + ys[j]))
            }
        }

        let cauchy_matrix = SquareMatrix::from_vec(elements);
        // Sanity check: All Cauchy matrices should be invertible
        assert!(cauchy_matrix.determinant() != F::zero());

        poseidon_parameters::MdsMatrix(cauchy_matrix)
    }

    /// Compute inverse of MDS matrix
    pub fn inverse(&self) -> poseidon_parameters::SquareMatrix<F> {
        self.0
             .0
            .inverse()
            .expect("all well-formed MDS matrices should have inverses")
    }

    /// Compute the (t - 1) x (t - 1) Mhat matrix from the MDS matrix
    ///
    /// This is simply the MDS matrix with the first row and column removed
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn hat(&self) -> poseidon_parameters::SquareMatrix<F> {
        let dim = self.0.n_rows();
        let mut mhat_elements = Vec::with_capacity((dim - 1) * (dim - 1));
        for i in 1..dim {
            for j in 1..dim {
                mhat_elements.push(self.0.get_element(i, j))
            }
        }

        SquareMatrix::from_vec(mhat_elements)
    }

    /// Return the elements M_{0,1} .. M_{0,t} from the first row
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn v(&self) -> poseidon_parameters::Matrix<F> {
        let elements: Vec<F> = self.0 .0.elements()[1..self.0 .0.n_rows()].to_vec();
        poseidon_parameters::Matrix::new(1, self.0.n_rows() - 1, elements)
    }

    /// Return the elements M_{1,0} .. M_{t,0}from the first column
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn w(&self) -> poseidon_parameters::Matrix<F> {
        let mut elements = Vec::with_capacity(self.0.n_rows() - 1);
        for i in 1..self.0.n_rows() {
            elements.push(self.0.get_element(i, 0))
        }
        poseidon_parameters::Matrix::new(self.0.n_rows() - 1, 1, elements)
    }
}

impl<F: PrimeField> MatrixOperations<F> for poseidon_parameters::MdsMatrix<F> {
    // fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> Self {
    //     Self(poseidon_parameters::SquareMatrix::new(
    //         n_rows, n_cols, elements,
    //     ))
    // }

    // fn elements(&self) -> &Vec<F> {
    //     self.0.elements()
    // }

    // fn n_rows(&self) -> usize {
    //     self.0.n_rows()
    // }

    // fn n_cols(&self) -> usize {
    //     self.0.n_cols()
    // }

    // fn get_element(&self, i: usize, j: usize) -> F {
    //     self.0.get_element(i, j)
    // }

    // fn set_element(&mut self, i: usize, j: usize, val: F) {
    //     self.0.set_element(i, j, val)
    // }
    // fn rows(&self) -> Vec<&[F]> {
    //     self.0.rows()
    // }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self> {
        Ok(Self(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F: PrimeField> From<MdsMatrix<poseidon_parameters::MdsMatrix<F>>> for Vec<Vec<F>> {
    fn from(val: MdsMatrix<poseidon_parameters::MdsMatrix<F>>) -> Self {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..val.0.n_rows() {
            let mut row = Vec::new();
            for j in 0..val.0.n_rows() {
                row.push(val.0 .0.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

impl<F: PrimeField> From<&MdsMatrix<poseidon_parameters::MdsMatrix<F>>> for Vec<F> {
    fn from(val: &MdsMatrix<poseidon_parameters::MdsMatrix<F>>) -> Self {
        let elements = val.0 .0.elements();
        elements.to_vec()
    }
}

/// Represents an optimized MDS (maximum distance separable) matrix.
pub struct OptimizedMdsMatrices<T>(pub T);

impl<F> OptimizedMdsMatrices<poseidon_parameters::OptimizedMdsMatrices<F>>
where
    F: PrimeField,
{
    /// Generate the optimized MDS matrices.
    pub fn generate(
        mds: &poseidon_parameters::MdsMatrix<F>,
        t: usize,
        rounds: &poseidon_parameters::RoundNumbers,
    ) -> poseidon_parameters::OptimizedMdsMatrices<F> {
        let M_hat = MdsMatrix(mds).hat();
        let M_hat_inverse = M_hat
            .inverse()
            .expect("all well-formed MDS matrices should have inverses");
        let v = MdsMatrix(mds).v();
        let w = MdsMatrix(mds).w();
        let M_prime = OptimizedMdsMatrices::prime(&M_hat);
        let M_00 = mds.get_element(0, 0);
        let M_doubleprime = OptimizedMdsMatrices::doubleprime(&M_hat_inverse, &w, &v, M_00);

        // Sanity checks
        assert_eq!(M_prime.n_cols(), mds.n_cols());
        assert_eq!(M_hat.n_cols() + 1, mds.n_cols());

        // If M' and M'' are well-formed, then M = M' * M'' (Eqn. 7, Appendix B)
        assert_eq!(
            mds.0.clone(),
            mat_mul(&M_prime, &M_doubleprime).expect("M' and M'' must have the same dimensions")
        );

        // If M'' is well-formed, it should be sparse with:
        // (t - 1)^2 - (t - 1) coefficients equal to 0
        // t - 1 coefficients equal to 1
        // (Text under Eqn. 7, Appendix B)
        assert_eq!(
            M_doubleprime
                .elements()
                .iter()
                .filter(|&n| *n == F::zero())
                .count(),
            (t - 1) * (t - 1) - (t - 1)
        );
        assert_eq!(
            M_doubleprime
                .elements()
                .iter()
                .filter(|&n| *n == F::one())
                .count(),
            t - 1
        );

        // From `calc_equivalent_matrices` in `poseidonperm_x3_64_24_optimized.sage`.
        let (M_i, v_collection, w_hat_collection) =
            OptimizedMdsMatrices::calc_equivalent_matrices(mds, rounds);

        poseidon_parameters::OptimizedMdsMatrices {
            M_hat,
            M_hat_inverse,
            v,
            w,
            M_prime,
            M_doubleprime,
            M_inverse: MdsMatrix(mds).inverse(),
            M_i,
            v_collection,
            w_hat_collection,
            M_00,
        }
    }

    pub(crate) fn calc_equivalent_matrices(
        mds: &poseidon_parameters::MdsMatrix<F>,
        rounds: &poseidon_parameters::RoundNumbers,
    ) -> (
        poseidon_parameters::Matrix<F>,
        Vec<poseidon_parameters::Matrix<F>>,
        Vec<poseidon_parameters::Matrix<F>>,
    ) {
        let rounds: RoundNumbers<poseidon_parameters::RoundNumbers> = RoundNumbers(*rounds);
        let r_P = rounds.partial();
        let mut w_hat_collection = Vec::with_capacity(rounds.partial());
        let mut v_collection = Vec::with_capacity(rounds.partial());

        let M_T = mds.transpose();
        let mut M_mul = M_T.clone();
        let mut M_i = OptimizedMdsMatrices::prime(&M_mul.0);

        for _ in (0..r_P).rev() {
            let M_hat = MdsMatrix(&M_mul).hat();
            let w = MdsMatrix(&M_mul).w();

            let v = MdsMatrix(&M_mul).v();
            v_collection.push(v);
            let w_hat = mat_mul(&M_hat.inverse().expect("can invert Mhat").0, &w)
                .expect("can compute w_hat");
            w_hat_collection.push(w_hat);

            // Now we compute M' and M * M' for the previous round
            M_i = OptimizedMdsMatrices::prime(&M_hat);

            M_mul = poseidon_parameters::MdsMatrix(
                mat_mul(&M_T.0, &M_i).expect("mds and M_i have same dimensions"),
            );
        }

        (M_i.0.transpose(), v_collection, w_hat_collection)
    }

    fn prime(M_hat: &poseidon_parameters::SquareMatrix<F>) -> poseidon_parameters::SquareMatrix<F> {
        let dim = M_hat.n_cols() + 1;
        let mut new_elements = Vec::with_capacity(dim * dim);

        for i in 0..dim {
            for j in 0..dim {
                if i == 0 && j == 0 {
                    // M_00 = 1
                    new_elements.push(F::one());
                } else if i == 0 || j == 0 {
                    new_elements.push(F::zero());
                } else {
                    new_elements.push(M_hat.get_element(i - 1, j - 1))
                }
            }
        }

        SquareMatrix::from_vec(new_elements)
    }

    fn doubleprime(
        M_hat_inverse: &poseidon_parameters::SquareMatrix<F>,
        w: &poseidon_parameters::Matrix<F>,
        v: &poseidon_parameters::Matrix<F>,
        M_00: F,
    ) -> poseidon_parameters::SquareMatrix<F> {
        let dim = M_hat_inverse.n_cols() + 1;
        let mut new_elements = Vec::with_capacity(dim * dim);
        let identity = poseidon_parameters::SquareMatrix::identity(dim - 1);
        let w_hat =
            mat_mul(&M_hat_inverse.0, w).expect("matrix multiplication should always exist");

        for i in 0..dim {
            for j in 0..dim {
                if i == 0 && j == 0 {
                    new_elements.push(M_00);
                } else if i == 0 {
                    new_elements.push(v.get_element(0, j - 1));
                } else if j == 0 {
                    new_elements.push(w_hat.get_element(i - 1, 0));
                } else {
                    new_elements.push(identity.get_element(i - 1, j - 1))
                }
            }
        }

        SquareMatrix::from_vec(new_elements)
    }
}

#[cfg(test)]
mod tests {
    use ark_ed_on_bls12_377::{Fq as Fq377, FqParameters as Fq377Parameters};
    use ark_ed_on_bls12_381::{Fq, FqParameters as Fq381Parameters};
    use ark_ff::{fields::FpParameters, One, Zero};

    use super::*;

    #[test]
    fn convert_from_mds_to_vec_of_vecs() {
        let MDS_matrix = poseidon_parameters::MdsMatrix(SquareMatrix::from_vec(vec![
            Fq::from(1u32),
            Fq::from(2u32),
            Fq::from(3u32),
            Fq::from(4u32),
        ]));
        let vec_of_vecs: Vec<Vec<Fq>> = MdsMatrix(MDS_matrix).into();
        assert_eq!(vec_of_vecs[0][0], Fq::from(1u32));
        assert_eq!(vec_of_vecs[0][1], Fq::from(2u32));
        assert_eq!(vec_of_vecs[1][0], Fq::from(3u32));
        assert_eq!(vec_of_vecs[1][1], Fq::from(4u32));
    }

    #[test]
    fn cauchy_method_mds() {
        let M = 128;
        let t = 3;

        let input = InputParameters::new(M, 3, Fq381Parameters::MODULUS, true);
        let MDS_matrix: poseidon_parameters::MdsMatrix<Fq> = MdsMatrix::generate(&input);

        assert!(MDS_matrix.0.determinant() != Fq::zero());
        assert_eq!(MDS_matrix.n_rows(), t);
        assert!(MDS_matrix.0.get_element(0, 0) != Fq::zero());
    }

    #[test]
    fn check_calc_equivalent_matrices_vs_sage() {
        let M = 128;

        let input = InputParameters::new(M, 3, Fq377Parameters::MODULUS, true);
        let rounds = RoundNumbers::new(&input, &poseidon_parameters::Alpha::Exponent(17));
        let mds: poseidon_parameters::MdsMatrix<Fq377> = MdsMatrix::generate(&input);
        let M_00 = mds.get_element(0, 0);
        // Sanity check
        assert_eq!(
            M_00,
            ark_ff::field_new!(
                Fq377,
                "5629641166285580282832549959187697687583932890102709218623488970611606159361"
            ),
        );

        let (M_i, v_collection, w_hat_collection) =
            OptimizedMdsMatrices::calc_equivalent_matrices(&mds, &rounds);

        // There are 31 (number of partial rounds) of these, we check the first 2 since it's the same method.
        let v_collection_expected = [
            [
                ark_ff::field_new!(
                    Fq377,
                    "6333346312071277818186618704086159898531924501365547870951425091938056929281"
                ),
                ark_ff::field_new!(
                    Fq377,
                    "6755569399542696339399059951025237225100719468123251062348186764733927391233"
                ),
            ],
            [
                ark_ff::field_new!(
                    Fq377,
                    "7740756603642672888894756193883084320427907723891225175607297334590958469121"
                ),
                ark_ff::field_new!(
                    Fq377,
                    "7851338840837568215878966996652842667862592119946814106687401582227972161537"
                ),
            ],
        ];
        for i in 0..v_collection_expected.len() {
            for (j, v_entry_computed) in v_collection[i].elements().iter().enumerate() {
                assert_eq!(*v_entry_computed, v_collection_expected[i][j]);
            }
        }

        let w_hat_collection_expected = [
            [
                ark_ff::field_new!(Fq377, "3"),
                ark_ff::field_new!(
                    Fq377,
                    "844446174942837042424882493878154653137589933515406382793523345591740923902"
                ),
            ],
            [
                ark_ff::field_new!(Fq377, "981"),
                ark_ff::field_new!(
                    Fq377,
                    "1688892349885674084849764987756309306275179867030812765587046691183481846649"
                ),
            ],
        ];
        for i in 0..w_hat_collection_expected.len() {
            for (j, w_hat_entry_computed) in w_hat_collection[i].elements().iter().enumerate() {
                assert_eq!(*w_hat_entry_computed, w_hat_collection_expected[i][j]);
            }
        }

        let M_i_expected = vec![
            Fq377::one(),
            Fq377::zero(),
            Fq377::zero(),
            Fq377::zero(),
            ark_ff::field_new!(
                Fq377,
                "1949629285152675843545617098663080067734218406516000484720630379218497119024"
            ),
            ark_ff::field_new!(
                Fq377,
                "6804287869450188502728877251894011667833647269738979685488937504164506768586"
            ),
            Fq377::zero(),
            ark_ff::field_new!(
                Fq377,
                "6804287869450188502728877251894011667833647269738979685488937504164506768586"
            ),
            ark_ff::field_new!(
                Fq377,
                "4924677972410444052137834859533533887056104638988047570112284264367323462906"
            ),
        ];
        for (a, b) in M_i.elements().iter().zip(M_i_expected.iter()) {
            assert_eq!(*a, *b);
        }
    }
}
