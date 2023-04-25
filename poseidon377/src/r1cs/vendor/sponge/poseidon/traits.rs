use crate::poseidon::grain_lfsr::PoseidonGrainLFSR;
use crate::poseidon::PoseidonParameters;
use ark_ff::{fields::models::*, FpParameters, PrimeField};
use ark_std::{vec, vec::Vec};

/// An entry in the default Poseidon parameters
pub struct PoseidonDefaultParametersEntry {
    /// The rate (in terms of number of field elements).
    pub rate: usize,
    /// Exponent used in S-boxes.
    pub alpha: usize,
    /// Number of rounds in a full-round operation.
    pub full_rounds: usize,
    /// Number of rounds in a partial-round operation.
    pub partial_rounds: usize,
    /// Number of matrices to skip when generating parameters using the Grain LFSR.
    ///
    /// The matrices being skipped are those that do not satisfy all the desired properties.
    /// See the [reference implementation](https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/generate_parameters_grain.sage) for more detail.
    pub skip_matrices: usize,
}

impl PoseidonDefaultParametersEntry {
    /// Create an entry in PoseidonDefaultParameters.
    pub const fn new(
        rate: usize,
        alpha: usize,
        full_rounds: usize,
        partial_rounds: usize,
        skip_matrices: usize,
    ) -> Self {
        Self {
            rate,
            alpha,
            full_rounds,
            partial_rounds,
            skip_matrices,
        }
    }
}

/// A trait for default Poseidon parameters associated with a prime field
pub trait PoseidonDefaultParameters: FpParameters {
    /// An array of the parameters optimized for constraints
    /// (rate, alpha, full_rounds, partial_rounds, skip_matrices)
    /// for rate = 2, 3, 4, 5, 6, 7, 8
    ///
    /// Here, `skip_matrices` denote how many matrices to skip before
    /// finding one that satisfy all the requirements.
    const PARAMS_OPT_FOR_CONSTRAINTS: [PoseidonDefaultParametersEntry; 7];

    /// An array of the parameters optimized for weights
    /// (rate, alpha, full_rounds, partial_rounds, skip_matrices)
    /// for rate = 2, 3, 4, 5, 6, 7, 8
    const PARAMS_OPT_FOR_WEIGHTS: [PoseidonDefaultParametersEntry; 7];
}

/// A field with Poseidon parameters associated
pub trait PoseidonDefaultParametersField: PrimeField {
    /// Obtain the default Poseidon parameters for this rate and for this prime field,
    /// with a specific optimization goal.
    fn get_default_poseidon_parameters(
        rate: usize,
        optimized_for_weights: bool,
    ) -> Option<PoseidonParameters<Self>>;
}

/// Internal function that uses the `PoseidonDefaultParameters` to compute the Poseidon parameters.
pub fn get_default_poseidon_parameters_internal<F: PrimeField, P: PoseidonDefaultParameters>(
    rate: usize,
    optimized_for_weights: bool,
) -> Option<PoseidonParameters<F>> {
    let params_set = if !optimized_for_weights {
        P::PARAMS_OPT_FOR_CONSTRAINTS
    } else {
        P::PARAMS_OPT_FOR_WEIGHTS
    };

    for param in params_set.iter() {
        if param.rate == rate {
            let (ark, mds) = find_poseidon_ark_and_mds::<F>(
                P::MODULUS_BITS as u64,
                rate,
                param.full_rounds as u64,
                param.partial_rounds as u64,
                param.skip_matrices as u64,
            );

            return Some(PoseidonParameters {
                full_rounds: param.full_rounds,
                partial_rounds: param.partial_rounds,
                alpha: param.alpha as u64,
                ark,
                mds,
                rate: param.rate,
                capacity: 1,
            });
        }
    }

    None
}

/// Internal function that computes the ark and mds from the Poseidon Grain LFSR.
pub fn find_poseidon_ark_and_mds<F: PrimeField>(
    prime_bits: u64,
    rate: usize,
    full_rounds: u64,
    partial_rounds: u64,
    skip_matrices: u64,
) -> (Vec<Vec<F>>, Vec<Vec<F>>) {
    let mut lfsr = PoseidonGrainLFSR::new(
        false,
        prime_bits,
        (rate + 1) as u64,
        full_rounds,
        partial_rounds,
    );

    let mut ark = Vec::<Vec<F>>::new();
    for _ in 0..(full_rounds + partial_rounds) {
        ark.push(lfsr.get_field_elements_rejection_sampling(rate + 1));
    }

    let mut mds = Vec::<Vec<F>>::new();
    mds.resize(rate + 1, vec![F::zero(); rate + 1]);
    for _ in 0..skip_matrices {
        let _ = lfsr.get_field_elements_mod_p::<F>(2 * (rate + 1));
    }

    // a qualifying matrix must satisfy the following requirements
    // - there is no duplication among the elements in x or y
    // - there is no i and j such that x[i] + y[j] = p
    // - the resultant MDS passes all the three tests

    let xs = lfsr.get_field_elements_mod_p::<F>(rate + 1);
    let ys = lfsr.get_field_elements_mod_p::<F>(rate + 1);

    for i in 0..(rate + 1) {
        for j in 0..(rate + 1) {
            mds[i][j] = (xs[i] + &ys[j]).inverse().unwrap();
        }
    }

    (ark, mds)
}

macro_rules! impl_poseidon_default_parameters_field {
    ($field: ident, $params: ident) => {
        impl<P: $params + PoseidonDefaultParameters> PoseidonDefaultParametersField for $field<P> {
            fn get_default_poseidon_parameters(
                rate: usize,
                optimized_for_weights: bool,
            ) -> Option<PoseidonParameters<Self>> {
                get_default_poseidon_parameters_internal::<Self, P>(rate, optimized_for_weights)
            }
        }
    };
}

impl_poseidon_default_parameters_field!(Fp64, Fp64Parameters);
impl_poseidon_default_parameters_field!(Fp256, Fp256Parameters);
impl_poseidon_default_parameters_field!(Fp320, Fp320Parameters);
impl_poseidon_default_parameters_field!(Fp384, Fp384Parameters);
impl_poseidon_default_parameters_field!(Fp448, Fp448Parameters);
impl_poseidon_default_parameters_field!(Fp768, Fp768Parameters);
impl_poseidon_default_parameters_field!(Fp832, Fp832Parameters);

#[cfg(test)]
mod test {
    use crate::poseidon::{
        PoseidonDefaultParameters, PoseidonDefaultParametersEntry, PoseidonDefaultParametersField,
    };
    use ark_ff::{field_new, fields::Fp256};
    use ark_ff::{BigInteger256, FftParameters, Fp256Parameters, FpParameters};
    use ark_test_curves::bls12_381::FrParameters;

    pub struct TestFrParameters;

    impl Fp256Parameters for TestFrParameters {}
    impl FftParameters for TestFrParameters {
        type BigInt = <FrParameters as FftParameters>::BigInt;
        const TWO_ADICITY: u32 = FrParameters::TWO_ADICITY;
        const TWO_ADIC_ROOT_OF_UNITY: Self::BigInt = FrParameters::TWO_ADIC_ROOT_OF_UNITY;
    }

    // This TestFrParameters is the same as the BLS12-381's Fr.
    // MODULUS = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    impl FpParameters for TestFrParameters {
        const MODULUS: BigInteger256 = FrParameters::MODULUS;
        const MODULUS_BITS: u32 = FrParameters::MODULUS_BITS;
        const CAPACITY: u32 = FrParameters::CAPACITY;
        const REPR_SHAVE_BITS: u32 = FrParameters::REPR_SHAVE_BITS;
        const R: BigInteger256 = FrParameters::R;
        const R2: BigInteger256 = FrParameters::R2;
        const INV: u64 = FrParameters::INV;
        const GENERATOR: BigInteger256 = FrParameters::GENERATOR;
        const MODULUS_MINUS_ONE_DIV_TWO: BigInteger256 = FrParameters::MODULUS_MINUS_ONE_DIV_TWO;
        const T: BigInteger256 = FrParameters::T;
        const T_MINUS_ONE_DIV_TWO: BigInteger256 = FrParameters::T_MINUS_ONE_DIV_TWO;
    }

    impl PoseidonDefaultParameters for TestFrParameters {
        const PARAMS_OPT_FOR_CONSTRAINTS: [PoseidonDefaultParametersEntry; 7] = [
            PoseidonDefaultParametersEntry::new(2, 17, 8, 31, 0),
            PoseidonDefaultParametersEntry::new(3, 5, 8, 56, 0),
            PoseidonDefaultParametersEntry::new(4, 5, 8, 56, 0),
            PoseidonDefaultParametersEntry::new(5, 5, 8, 57, 0),
            PoseidonDefaultParametersEntry::new(6, 5, 8, 57, 0),
            PoseidonDefaultParametersEntry::new(7, 5, 8, 57, 0),
            PoseidonDefaultParametersEntry::new(8, 5, 8, 57, 0),
        ];
        const PARAMS_OPT_FOR_WEIGHTS: [PoseidonDefaultParametersEntry; 7] = [
            PoseidonDefaultParametersEntry::new(2, 257, 8, 13, 0),
            PoseidonDefaultParametersEntry::new(3, 257, 8, 13, 0),
            PoseidonDefaultParametersEntry::new(4, 257, 8, 13, 0),
            PoseidonDefaultParametersEntry::new(5, 257, 8, 13, 0),
            PoseidonDefaultParametersEntry::new(6, 257, 8, 13, 0),
            PoseidonDefaultParametersEntry::new(7, 257, 8, 13, 0),
            PoseidonDefaultParametersEntry::new(8, 257, 8, 13, 0),
        ];
    }

    pub type TestFr = Fp256<TestFrParameters>;

    #[test]
    fn bls12_381_fr_poseidon_default_parameters_test() {
        // constraints
        let constraints_rate_2 = TestFr::get_default_poseidon_parameters(2, false).unwrap();
        assert_eq!(
            constraints_rate_2.ark[0][0],
            field_new!(
                TestFr,
                "27117311055620256798560880810000042840428971800021819916023577129547249660720"
            )
        );
        assert_eq!(
            constraints_rate_2.mds[0][0],
            field_new!(
                TestFr,
                "26017457457808754696901916760153646963713419596921330311675236858336250747575"
            )
        );

        let constraints_rate_3 = TestFr::get_default_poseidon_parameters(3, false).unwrap();
        assert_eq!(
            constraints_rate_3.ark[0][0],
            field_new!(
                TestFr,
                "11865901593870436687704696210307853465124332568266803587887584059192277437537"
            )
        );
        assert_eq!(
            constraints_rate_3.mds[0][0],
            field_new!(
                TestFr,
                "18791275321793747281053101601584820964683215017313972132092847596434094368732"
            )
        );

        let constraints_rate_4 = TestFr::get_default_poseidon_parameters(4, false).unwrap();
        assert_eq!(
            constraints_rate_4.ark[0][0],
            field_new!(
                TestFr,
                "41775194144383840477168997387904574072980173775424253289429546852163474914621"
            )
        );
        assert_eq!(
            constraints_rate_4.mds[0][0],
            field_new!(
                TestFr,
                "42906651709148432559075674119637355642263148226238482628104108168707874713729"
            )
        );

        let constraints_rate_5 = TestFr::get_default_poseidon_parameters(5, false).unwrap();
        assert_eq!(
            constraints_rate_5.ark[0][0],
            field_new!(
                TestFr,
                "24877380261526996562448766783081897666376381975344509826094208368479247894723"
            )
        );
        assert_eq!(
            constraints_rate_5.mds[0][0],
            field_new!(
                TestFr,
                "30022080821787948421423927053079656488514459012053372877891553084525866347732"
            )
        );

        let constraints_rate_6 = TestFr::get_default_poseidon_parameters(6, false).unwrap();
        assert_eq!(
            constraints_rate_6.ark[0][0],
            field_new!(
                TestFr,
                "37928506567864057383105673253383925733025682403141583234734361541053005808936"
            )
        );
        assert_eq!(
            constraints_rate_6.mds[0][0],
            field_new!(
                TestFr,
                "49124738641420159156404016903087065194698370461819821829905285681776084204443"
            )
        );

        let constraints_rate_7 = TestFr::get_default_poseidon_parameters(7, false).unwrap();
        assert_eq!(
            constraints_rate_7.ark[0][0],
            field_new!(
                TestFr,
                "37848764121158464546907147011864524711588624175161409526679215525602690343051"
            )
        );
        assert_eq!(
            constraints_rate_7.mds[0][0],
            field_new!(
                TestFr,
                "28113878661515342855868752866874334649815072505130059513989633785080391114646"
            )
        );

        let constraints_rate_8 = TestFr::get_default_poseidon_parameters(8, false).unwrap();
        assert_eq!(
            constraints_rate_8.ark[0][0],
            field_new!(
                TestFr,
                "51456871630395278065627483917901523970718884366549119139144234240744684354360"
            )
        );
        assert_eq!(
            constraints_rate_8.mds[0][0],
            field_new!(
                TestFr,
                "12929023787467701044434927689422385731071756681420195282613396560814280256210"
            )
        );

        // weights
        let weights_rate_2 = TestFr::get_default_poseidon_parameters(2, true).unwrap();
        assert_eq!(
            weights_rate_2.ark[0][0],
            field_new!(
                TestFr,
                "25126470399169474618535500283750950727260324358529540538588217772729895991183"
            )
        );
        assert_eq!(
            weights_rate_2.mds[0][0],
            field_new!(
                TestFr,
                "46350838805835525240431215868760423854112287760212339623795708191499274188615"
            )
        );

        let weights_rate_3 = TestFr::get_default_poseidon_parameters(3, true).unwrap();
        assert_eq!(
            weights_rate_3.ark[0][0],
            field_new!(
                TestFr,
                "16345358380711600255519479157621098002794924491287389755192263320486827897573"
            )
        );
        assert_eq!(
            weights_rate_3.mds[0][0],
            field_new!(
                TestFr,
                "37432344439659887296708509941462699942272362339508052702346957525719991245918"
            )
        );

        let weights_rate_4 = TestFr::get_default_poseidon_parameters(4, true).unwrap();
        assert_eq!(
            weights_rate_4.ark[0][0],
            field_new!(
                TestFr,
                "2997721997773001075802235431463112417440167809433966871891875582435098138600"
            )
        );
        assert_eq!(
            weights_rate_4.mds[0][0],
            field_new!(
                TestFr,
                "43959024692079347032841256941012668338943730711936867712802582656046301966186"
            )
        );

        let weights_rate_5 = TestFr::get_default_poseidon_parameters(5, true).unwrap();
        assert_eq!(
            weights_rate_5.ark[0][0],
            field_new!(
                TestFr,
                "28142027771717376151411984909531650866105717069245696861966432993496676054077"
            )
        );
        assert_eq!(
            weights_rate_5.mds[0][0],
            field_new!(
                TestFr,
                "13157425078305676755394500322568002504776463228389342308130514165393397413991"
            )
        );

        let weights_rate_6 = TestFr::get_default_poseidon_parameters(6, true).unwrap();
        assert_eq!(
            weights_rate_6.ark[0][0],
            field_new!(
                TestFr,
                "7417004907071346600696060525974582183666365156576759507353305331252133694222"
            )
        );
        assert_eq!(
            weights_rate_6.mds[0][0],
            field_new!(
                TestFr,
                "51393878771453405560681338747290999206747890655420330824736778052231938173954"
            )
        );

        let weights_rate_7 = TestFr::get_default_poseidon_parameters(7, true).unwrap();
        assert_eq!(
            weights_rate_7.ark[0][0],
            field_new!(
                TestFr,
                "47093173418416013663709314805327945458844779999893881721688570889452680883650"
            )
        );
        assert_eq!(
            weights_rate_7.mds[0][0],
            field_new!(
                TestFr,
                "51455917624412053400160569105425532358410121118308957353565646758865245830775"
            )
        );

        let weights_rate_8 = TestFr::get_default_poseidon_parameters(8, true).unwrap();
        assert_eq!(
            weights_rate_8.ark[0][0],
            field_new!(
                TestFr,
                "16478680729975035007348178961232525927424769683353433314299437589237598655079"
            )
        );
        assert_eq!(
            weights_rate_8.mds[0][0],
            field_new!(
                TestFr,
                "39160448583049384229582837387246752222769278402304070376350288593586064961857"
            )
        );
    }
}
