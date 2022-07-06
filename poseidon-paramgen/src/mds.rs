use ark_ff::PrimeField;

use crate::{matrix::mat_mul, InputParameters, Matrix, SquareMatrix};

/// Represents an MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MdsMatrix<F: PrimeField>(pub SquareMatrix<F>);

impl<F> MdsMatrix<F>
where
    F: PrimeField,
{
    pub fn new(input: &InputParameters<F::BigInt>) -> Self {
        // A t x t MDS matrix only exists if: 2t + 1 <= p
        let two_times_t_bigint: F::BigInt = (2 * input.t as u64).into();
        if two_times_t_bigint > input.p {
            panic!("no MDS matrix exists");
        }

        MdsMatrix::fixed_cauchy_matrix(input)
    }

    /// Dimension of the (square) MDS matrix
    pub fn dim(&self) -> usize {
        self.0.inner.n_rows
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
    /// However, Section 5.4 of [Keller and Rosemarin 2020](https://eprint.iacr.org/2020/179.pdf)
    /// describes how the MDS matrix can be constructed in a deterministic fashion
    /// where infinitely long subspace trails cannot be constructed. This method constructs an MDS
    /// matrix using that deterministic method, avoiding the need to implement Algorithms 1-3.
    pub fn fixed_cauchy_matrix(input: &InputParameters<F::BigInt>) -> Self {
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

        Self(cauchy_matrix)
    }

    /// Compute inverse of MDS matrix
    pub fn inverse(&self) -> SquareMatrix<F> {
        self.0.inverse()
    }

    pub fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    /// Compute the (t - 1) x (t - 1) Mhat matrix from the MDS matrix
    ///
    /// This is simply the MDS matrix with the first row and column removed
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn hat(&self) -> SquareMatrix<F> {
        let dim = self.dim();
        let mut mhat_elements = Vec::with_capacity((dim - 1) * (dim - 1));
        for i in 1..dim {
            for j in 1..dim {
                mhat_elements.push(self.get_element(i, j))
            }
        }

        SquareMatrix::from_vec(mhat_elements)
    }

    /// Return the elements M_{0,1} .. M_{0,t} from the first row
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn v(&self) -> Matrix<F> {
        let elements: Vec<F> = self.0.elements()[1..self.dim()].to_vec();
        Matrix::new(1, self.dim() - 1, elements)
    }

    /// Return the elements M_{1,0} .. M_{t,0} from the first column
    ///
    /// Ref: p.20 of the Poseidon paper
    pub fn w(&self) -> Matrix<F> {
        let mut elements = Vec::with_capacity(self.dim() - 1);
        for i in 1..self.dim() {
            elements.push(self.get_element(i, 0))
        }
        Matrix::new(self.dim() - 1, 1, elements)
    }
}

impl<F: PrimeField> Into<Vec<Vec<F>>> for MdsMatrix<F> {
    fn into(self) -> Vec<Vec<F>> {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..self.dim() {
            let mut row = Vec::new();
            for j in 0..self.dim() {
                row.push(self.0.get_element(i, j));
            }
            rows.push(row);
        }
        rows
    }
}

/// Represents an optimized MDS (maximum distance separable) matrix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OptimizedMdsMatrices<F: PrimeField> {
    /// The t x t MDS matrix of the linear layer.
    pub M: MdsMatrix<F>,
    /// A (t - 1) x (t - 1) MDS submatrix derived from the MDS matrix.
    pub M_hat: SquareMatrix<F>,
    /// A 1 x (t - 1) (row) vector derived from the MDS matrix.
    pub v: Matrix<F>,
    /// A (t - 1) x 1 (column) vector derived from the MDS matrix.
    pub w: Matrix<F>,
    /// A matrix formed from Mhat (an MDS submatrix of the MDS matrix).
    pub M_prime: SquareMatrix<F>,
    /// A sparse matrix formed from M,
    pub M_doubleprime: SquareMatrix<F>,
    /// The inverse of the t x t MDS matrix (needed to compute round constants).
    pub M_inverse: SquareMatrix<F>,
    /// The inverse of the (t - 1) x (t - 1) Mhat matrix.
    pub M_hat_inverse: SquareMatrix<F>,
}

impl<F> OptimizedMdsMatrices<F>
where
    F: PrimeField,
{
    pub fn generate(mds: &MdsMatrix<F>, t: usize) -> OptimizedMdsMatrices<F> {
        let M_hat = mds.hat();
        let M_hat_inverse = M_hat.inverse();
        let v = mds.v();
        let w = mds.w();
        let M_prime = OptimizedMdsMatrices::prime(&M_hat);
        let M_00 = mds.get_element(0, 0);
        let M_doubleprime = OptimizedMdsMatrices::doubleprime(&M_hat_inverse, &w, &v, M_00);

        // Sanity checks
        assert_eq!(M_prime.dim(), mds.dim());
        assert_eq!(M_hat.dim() + 1, mds.dim());

        // If M' and M'' are well-formed, then M = M' * M'' (Eqn. 7, Appendix B)
        assert_eq!(mds.0.clone(), &M_prime * &M_doubleprime);

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

        OptimizedMdsMatrices {
            M: mds.clone(),
            M_hat,
            M_hat_inverse,
            v,
            w,
            M_prime,
            M_doubleprime,
            M_inverse: mds.inverse(),
        }
    }

    fn prime(M_hat: &SquareMatrix<F>) -> SquareMatrix<F> {
        let dim = M_hat.dim() + 1;
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
        M_hat_inverse: &SquareMatrix<F>,
        w: &Matrix<F>,
        v: &Matrix<F>,
        M_00: F,
    ) -> SquareMatrix<F> {
        let dim = M_hat_inverse.dim() + 1;
        let mut new_elements = Vec::with_capacity(dim * dim);
        let identity = SquareMatrix::identity(dim - 1);
        let w_hat =
            mat_mul(&M_hat_inverse.inner, w).expect("matrix multiplication should always exist");

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
    use super::*;

    use ark_ed_on_bls12_381::{Fq, FqParameters as Fq381Parameters};
    use ark_ff::{fields::FpParameters, Zero};

    #[test]
    fn convert_from_mds_to_vec_of_vecs() {
        let MDS_matrix = MdsMatrix(SquareMatrix::from_vec(vec![
            Fq::from(1u32),
            Fq::from(2u32),
            Fq::from(3u32),
            Fq::from(4u32),
        ]));
        let vec_of_vecs: Vec<Vec<Fq>> = MDS_matrix.into();
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
        let MDS_matrix: MdsMatrix<Fq> = MdsMatrix::new(&input);

        assert!(MDS_matrix.0.determinant() != Fq::zero());
        assert_eq!(MDS_matrix.dim(), t);
        assert!(MDS_matrix.0.get_element(0, 0) != Fq::zero());
    }
}
