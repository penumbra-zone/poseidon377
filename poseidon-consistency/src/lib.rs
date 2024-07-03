#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(all(test, feature = "arkworks"))]
mod tests {
    use ark_ed_on_bls12_377::Fq as ArkFq;
    use ark_ff::BigInteger256;

    use decaf377::Fq;
    use poseidon377::{RATE_2_PARAMS, RATE_4_PARAMS};
    use poseidon_permutation::Instance;
    use proptest::prelude::*;

    pub(crate) fn from_ark_fq(x: ArkFq) -> Fq {
        BigInteger256::from(x).into()
    }

    #[test]
    fn check_optimized_impl_vs_sage() {
        let params_2_to_1 = RATE_2_PARAMS;
        let mut our_instance = Instance::new(&params_2_to_1);
        let hash_output =
            our_instance.n_to_1_fixed_hash(&[Fq::from(0u64), Fq::from(1u64), Fq::from(2u64)]);
        let output_words = our_instance.output_words();
        assert_eq!(hash_output, output_words[1]);
        let expected_output_words = [
            from_ark_fq(ark_ff::MontFp!(
                "6368779772888548211318735707249600947486536081021109980085678920065117075165"
            )),
            from_ark_fq(ark_ff::MontFp!(
                "546637332213889354237126997303352456465330747031466737868777261691321847212"
            )),
            from_ark_fq(ark_ff::MontFp!(
                "1488544471679348337017344865262529731114801536476862121626711131361325263279"
            )),
        ];
        for (a, b) in expected_output_words.iter().zip(output_words.iter()) {
            assert_eq!(*a, *b);
        }
    }

    #[test]
    fn check_unoptimized_impl_vs_sage() {
        let params_2_to_1 = RATE_2_PARAMS;
        let mut our_instance = Instance::new(&params_2_to_1);
        let hash_output = our_instance.unoptimized_n_to_1_fixed_hash([
            Fq::from(0u64),
            Fq::from(1u64),
            Fq::from(2u64),
        ]);
        let output_words = our_instance.output_words();
        assert_eq!(hash_output, output_words[1]);
        let expected_output_words = [
            from_ark_fq(ark_ff::MontFp!(
                "6368779772888548211318735707249600947486536081021109980085678920065117075165"
            )),
            from_ark_fq(ark_ff::MontFp!(
                "546637332213889354237126997303352456465330747031466737868777261691321847212"
            )),
            from_ark_fq(ark_ff::MontFp!(
                "1488544471679348337017344865262529731114801536476862121626711131361325263279"
            )),
        ];
        for (a, b) in expected_output_words.iter().zip(output_words.iter()) {
            assert_eq!(*a, *b);
        }
    }

    fn fq_strategy() -> BoxedStrategy<Fq> {
        any::<[u8; 32]>()
            .prop_map(|bytes| Fq::from_le_bytes_mod_order(&bytes[..]))
            .boxed()
    }

    proptest! {
        #[test]
        fn optimized_and_unoptimized_permutation_consistent(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy(), elem_4 in fq_strategy(), elem_5 in fq_strategy()) {
            let params_4_to_1 = RATE_4_PARAMS;

            let mut our_instance = Instance::new(&params_4_to_1);
            let our_result = our_instance.n_to_1_fixed_hash(&[elem_1, elem_2, elem_3, elem_4, elem_5]);

            let mut unoptimized_instance = Instance::new(&params_4_to_1);
            let unoptimized_result =
                unoptimized_instance.unoptimized_n_to_1_fixed_hash([elem_1, elem_2, elem_3, elem_4, elem_5]);

            assert_eq!(unoptimized_result, our_result);
        }
    }
}
