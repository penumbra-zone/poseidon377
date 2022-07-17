use anyhow::Result;
use ark_ff::PrimeField;
use merlin::Transcript;

use crate::{
    matrix::mat_mul, transcript::TranscriptProtocol, Alpha, InputParameters, Matrix,
    MatrixOperations, OptimizedMdsMatrices, RoundNumbers,
};

/// Represents an matrix of round constants.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArcMatrix<F: PrimeField>(pub Matrix<F>);

impl<F> ArcMatrix<F>
where
    F: PrimeField,
{
    /// Generate round constants.
    pub fn generate(
        input: &InputParameters<F::BigInt>,
        round_numbers: RoundNumbers,
        alpha: Alpha,
    ) -> ArcMatrix<F> {
        let mut transcript = Transcript::new(b"round-constants");
        transcript.domain_sep::<F>(input, round_numbers, alpha);

        let num_total_rounds = round_numbers.total();
        let elements = (0..num_total_rounds * input.t)
            .map(|_| transcript.round_constant())
            .collect();
        ArcMatrix(Matrix::new(num_total_rounds, input.t, elements))
    }

    /// Get row vector of constants by round
    pub fn constants_by_round(&self, r: usize) -> Matrix<F> {
        self.0.row_vector(r)
    }

    /// Set row vector of constants by round
    pub(crate) fn set_constants_by_round(&mut self, r: usize, constants: Vec<F>) {
        assert_eq!(constants.len(), self.n_cols());
        for (j, value) in constants.into_iter().enumerate() {
            self.set_element(r, j, value);
        }
    }
}

impl<F: PrimeField> MatrixOperations<F> for ArcMatrix<F> {
    fn new(n_rows: usize, n_cols: usize, elements: Vec<F>) -> ArcMatrix<F> {
        ArcMatrix(Matrix::new(n_rows, n_cols, elements))
    }

    fn elements(&self) -> &Vec<F> {
        self.0.elements()
    }

    fn n_rows(&self) -> usize {
        self.0.n_rows()
    }

    fn n_cols(&self) -> usize {
        self.0.n_cols()
    }

    fn get_element(&self, i: usize, j: usize) -> F {
        self.0.get_element(i, j)
    }

    fn set_element(&mut self, i: usize, j: usize, val: F) {
        self.0.set_element(i, j, val)
    }

    fn rows(&self) -> Vec<&[F]> {
        self.0.rows()
    }

    fn transpose(&self) -> Self {
        ArcMatrix(self.0.transpose())
    }

    fn hadamard_product(&self, rhs: &Self) -> Result<Self> {
        Ok(ArcMatrix(self.0.hadamard_product(&rhs.0)?))
    }
}

impl<F: PrimeField> Into<Vec<Vec<F>>> for ArcMatrix<F> {
    fn into(self) -> Vec<Vec<F>> {
        let mut rows = Vec::<Vec<F>>::new();
        for i in 0..self.n_rows() {
            let mut row = Vec::new();
            for j in 0..self.n_cols() {
                row.push(self.0.get_element(i, j));
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
pub struct OptimizedArcMatrix<F: PrimeField>(pub ArcMatrix<F>);

impl<F> OptimizedArcMatrix<F>
where
    F: PrimeField,
{
    /// Generate the optimized round constants.
    pub fn generate(
        arc: &ArcMatrix<F>,
        mds: &OptimizedMdsMatrices<F>,
        rounds: &RoundNumbers,
    ) -> OptimizedArcMatrix<F> {
        let n_cols = arc.n_cols();
        let mut constants_temp = arc.clone();
        let r_f = rounds.full() / 2;
        let r_T = rounds.total();
        let mds_T = mds.M.transpose();
        let mds_inv = &mds_T.inverse();

        // C_i = M^-1 * C_(i+1)
        for r in ((r_f)..(r_T - 1 - r_f)).rev() {
            // inv_cip1 = list(vector(constants_temp[i+1]) * MDS_matrix_field_transpose.inverse())
            let inv_cip1 = mat_mul(&constants_temp.constants_by_round(r + 1), &mds_inv.0).expect(
                "matrix multiplication of row of ARC matrix and MDS must have correct dimensions",
            );
            assert_eq!(inv_cip1.n_cols(), n_cols);
            assert_eq!(inv_cip1.n_rows(), 1);

            // constants_temp[i] = list(vector(constants_temp[i]) + vector([0] + inv_cip1[1:]))
            let mut delta_new_constants_r = Vec::with_capacity(n_cols);
            delta_new_constants_r.push(F::zero());
            for j in 1..n_cols {
                delta_new_constants_r.push(inv_cip1.get_element(0, j));
            }
            for j in 0..n_cols {
                let curr_element = constants_temp.get_element(r, j);
                constants_temp.set_element(r, j, curr_element + delta_new_constants_r[j]);
            }

            // constants_temp[i+1] = [inv_cip1[0]] + [0] * (t-1)
            let mut new_constants_row_i_plus_1 = Vec::with_capacity(n_cols);
            new_constants_row_i_plus_1.push(inv_cip1.get_element(0, 0));
            new_constants_row_i_plus_1.extend(vec![F::zero(); n_cols - 1]);
            constants_temp.set_constants_by_round(r + 1, new_constants_row_i_plus_1);
        }

        OptimizedArcMatrix(constants_temp)
    }
}

#[cfg(test)]
mod tests {
    use crate::MdsMatrix;

    use super::*;

    use ark_ed_on_bls12_377::{Fq, FqParameters};
    use ark_ff::FpParameters;

    #[test]
    fn convert_from_arc_to_vec_of_vecs() {
        let arc_matrix = ArcMatrix(Matrix::new(
            2,
            3,
            vec![
                Fq::from(1u32),
                Fq::from(2u32),
                Fq::from(0u32),
                Fq::from(4u32),
                Fq::from(5u32),
                Fq::from(6u32),
            ],
        ));
        let vec_of_vecs: Vec<Vec<Fq>> = arc_matrix.into();
        assert_eq!(vec_of_vecs[0][0], Fq::from(1u32));
        assert_eq!(vec_of_vecs[0][1], Fq::from(2u32));
        assert_eq!(vec_of_vecs[0][2], Fq::from(0u32));
        assert_eq!(vec_of_vecs[1][0], Fq::from(4u32));
        assert_eq!(vec_of_vecs[1][1], Fq::from(5u32));
        assert_eq!(vec_of_vecs[1][2], Fq::from(6u32));
    }

    #[test]
    fn check_cal_equivalent_constants_vs_sage() {
        let M = 128;
        let t = 3;
        let alpha = Alpha::Exponent(17);

        let input = InputParameters::new(M, 3, FqParameters::MODULUS, true);
        let rounds = RoundNumbers::new(&input, &alpha);
        let mds: MdsMatrix<Fq> = MdsMatrix::generate(&input);
        let optimized_mds = OptimizedMdsMatrices::generate(&mds, t, &rounds);
        let arc = ArcMatrix::generate(&input, rounds, alpha);
        let computed_constants = OptimizedArcMatrix::generate(&arc, &optimized_mds, &rounds);

        // Check of the first three rows, TODO: the rest
        let constants_expected = [
            // i = 0
            ark_ff::field_new!(
                Fq,
                "6988674540205220184161726804233283990758776299863834762040727861674579978543"
            ),
            ark_ff::field_new!(
                Fq,
                "7709355290596865338819285880138672457952222650642807855871317708209074795474"
            ),
            ark_ff::field_new!(
                Fq,
                "5201859429164217460310636175253377554540105593885866289756080180570368021951"
            ),
            // i = 1
            ark_ff::field_new!(
                Fq,
                "3049272212137844517764629946141402720245409346788255860423236916151652752101"
            ),
            ark_ff::field_new!(
                Fq,
                "731916334316210902202615695074797267498900349008237468358490947620590797518"
            ),
            ark_ff::field_new!(
                Fq,
                "6948583625231201615547345467810234311420481469329722516969162382644726485729"
            ),
            // i = 2
            ark_ff::field_new!(
                Fq,
                "4490233413574409524500190972213779994057100500913568615764229347275941966621"
            ),
            ark_ff::field_new!(
                Fq,
                "7396725314962185388189273989459617572911903025336675876756586397590970741108"
            ),
            ark_ff::field_new!(
                Fq,
                "2921080621736132447186011889175418237339394704314886136898495644039034276407"
            ),
            // i = 3
            ark_ff::field_new!(
                Fq,
                "6913381781295218724328701128337217653607031228663005211451263519103715771612"
            ),
            ark_ff::field_new!(
                Fq,
                "6019749413834726725125133350812094485497924094481882289540053799149578344435"
            ),
            ark_ff::field_new!(
                Fq,
                "4995954379914570178339745532281620079808702064916164361475898096048200561372"
            ),
            // i = 4
            ark_ff::field_new!(
                Fq,
                "2124186257537670932969443833628233916733000677872857619541493107812792994523"
            ),
            ark_ff::field_new!(
                Fq,
                "8046680783590808119458575993353618411934235335611096278481414605732256966748"
            ),
            ark_ff::field_new!(
                Fq,
                "7728864115296354631254506702385454419649201463379437684964755668131665373770"
            ),
            // i = 5
            ark_ff::field_new!(
                Fq,
                "1587876443634129168020854493450277603551026930347072939013547982925007616083"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 6
            ark_ff::field_new!(
                Fq,
                "1872517857299829216110854288967924097774976347427133022280076211211193761008"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 7
            ark_ff::field_new!(
                Fq,
                "2865968282298738946916952328725344253782807421536873039355328169023888530245"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 8
            ark_ff::field_new!(
                Fq,
                "7455327835606292788170559178990685552956519725450625518991150460337898502981"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 9
            ark_ff::field_new!(
                Fq,
                "6329682633220718897156926033656135046726106667999290597037736437537723647602"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 10
            ark_ff::field_new!(
                Fq,
                "3613027807341394963755611118675678841206509001829959688448877468113780695724"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 11
            ark_ff::field_new!(
                Fq,
                "4788537307194901360939650381012807804357630544329484629714912311521437518228"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 12
            ark_ff::field_new!(
                Fq,
                "5047006190758477600618328833686225035382965846715389589078702055169801277678"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 12
            ark_ff::field_new!(
                Fq,
                "5256063589523754991030306148490153819472126518330311867790764914234366347354"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 13
            ark_ff::field_new!(
                Fq,
                "7285844395274516280140535132510134860589365682920403294051692871565318200362"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 14
            ark_ff::field_new!(
                Fq,
                "3600450171718300138992516292595059700619341198177343728601472554774500889476"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 15
            ark_ff::field_new!(
                Fq,
                "7787989087100846143340988293215203273810426179739496802830516609907837611081"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 16
            ark_ff::field_new!(
                Fq,
                "7301882523327699565270753526439954192354475444570480963263305003830476657786"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 17
            ark_ff::field_new!(
                Fq,
                "763594137157752635985768278684803492830821444914620994297553472211315661849"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 18
            ark_ff::field_new!(
                Fq,
                "6046498698570211031132668535520621848965660746784736190544802469261553385823"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 19
            ark_ff::field_new!(
                Fq,
                "3150441867260293943485507376579844708517672260313707707401124636715975577108"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 20
            ark_ff::field_new!(
                Fq,
                "4769790453160546731149476332055196771343686765835180717452851622105817348420"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 21
            ark_ff::field_new!(
                Fq,
                "2730101715176213557952241542583779359812962074117477792886312349093815716063"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 22
            ark_ff::field_new!(
                Fq,
                "4038712952177210776672167805667007025979746203446904732638144157814162874974"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 23
            ark_ff::field_new!(
                Fq,
                "8255884447375225368974508845230187073113405931036329251808757909462690117877"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 24
            ark_ff::field_new!(
                Fq,
                "6885357291866532881065325687296979182338400685466768253292476430882961049206"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 25
            ark_ff::field_new!(
                Fq,
                "5926450599079603851023669220114471084498684957196677061516431568252479679758"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 26
            ark_ff::field_new!(
                Fq,
                "3377356905616162065433605641434977254712955961607034684650444069886455658547"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 27
            ark_ff::field_new!(
                Fq,
                "7087123987541754752647886830141457307401543688775490970072745121708862300128"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 28
            ark_ff::field_new!(
                Fq,
                "552337688682599454634788361316947247608548473790723890908543562593850437500"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 29
            ark_ff::field_new!(
                Fq,
                "8316842625469401119320722072820434224919330111623009975390943860375758190013"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 30
            ark_ff::field_new!(
                Fq,
                "7470799711651047173375742255027801566960556043062715772396462568348766982673"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 31
            ark_ff::field_new!(
                Fq,
                "1103088797790769386248055205940160107160509043268595290127658486406970058937"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 32
            ark_ff::field_new!(
                Fq,
                "2066764216865210667046992475368855181241427871253646185209411533419820909441"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 33
            ark_ff::field_new!(
                Fq,
                "3419605300888091486969655282591785054139800483892141159786787616140421786221"
            ),
            ark_ff::field_new!(Fq, "0"),
            ark_ff::field_new!(Fq, "0"),
            // i = 34
            ark_ff::field_new!(
                Fq,
                "3224051331877783069937246154883192018725119865352760355089784651290920906906"
            ),
            ark_ff::field_new!(
                Fq,
                "6344299245364913590728081469479468396346568591287952250450228130036748752624"
            ),
            ark_ff::field_new!(
                Fq,
                "7163298144698758147180460135771716060107526826354172769438296931285647166270"
            ),
            // i = 34
            ark_ff::field_new!(
                Fq,
                "5439669811527007432129035770771078324065108376477325977772842645864184474263"
            ),
            ark_ff::field_new!(
                Fq,
                "2388540132512935462910769151388489826759515381039083966425687579783397147680"
            ),
            ark_ff::field_new!(
                Fq,
                "5552242451981337643030834536688771349936013591654366196694153732254259211846"
            ),
            // i = 35
            ark_ff::field_new!(
                Fq,
                "2293718155480151242211995443457360332861427853767599092412745658932866271280"
            ),
            ark_ff::field_new!(
                Fq,
                "1114026302922263706493655417918518257260056401972366120358912283469308695784"
            ),
            ark_ff::field_new!(
                Fq,
                "8382504281541902573191937024713512277106734793833129272685742500601751387802"
            ),
            // i = 36
            ark_ff::field_new!(
                Fq,
                "7918246066354990409761145553536472732928632746154569124091286249273194443468"
            ),
            ark_ff::field_new!(
                Fq,
                "1937283435333591890382964471447412728772736446554658078729798048971164286629"
            ),
            ark_ff::field_new!(
                Fq,
                "2653607876727675876560671684767499894868406286850781422102673956317680932197"
            ),
        ];
        for (a, b) in constants_expected
            .iter()
            .zip(computed_constants.0.elements().iter())
        {
            assert_eq!(*a, *b);
        }
    }
}
