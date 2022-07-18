use crate::{Fq, Instance};

/// Hash a single [`Fq`] element with the provided `domain_separator`.
pub fn hash_1(domain_separator: &Fq, value: Fq) -> Fq {
    let params = &crate::RATE_1_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(vec![domain_separator.clone(), value])
}

/// Hash two [`Fq`] elements with the provided `domain_separator`.
pub fn hash_2(domain_separator: &Fq, value: (Fq, Fq)) -> Fq {
    let params = &crate::RATE_2_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(vec![domain_separator.clone(), value.0, value.1])
}

/// Hash three [`Fq`] elements with the provided `domain_separator`.
pub fn hash_3(domain_separator: &Fq, value: (Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_3_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(vec![domain_separator.clone(), value.0, value.1, value.2])
}

/// Hash four [`Fq`] elements with the provided `domain_separator`.
pub fn hash_4(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_4_PARAMS;
    let mut state = Instance::new(&params);
    state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
    ])
}

/// Hash five [`Fq`] elements with the provided `domain_separator`.
pub fn hash_5(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_5_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
        value.4,
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    use ark_ed_on_bls12_377::Fq;
    use ark_ff::PrimeField;
    use ark_sponge::poseidon::State;

    use crate::params;

    fn fq_strategy() -> BoxedStrategy<Fq> {
        any::<[u8; 32]>()
            .prop_map(|bytes| Fq::from_le_bytes_mod_order(&bytes[..]))
            .boxed()
    }

    fn old_hash_1(domain_separator: &Fq, value: Fq) -> Fq {
        let mut state = State::from(params::rate_1());

        // Use the domain separator as the sponge's capacity element
        state[0] = domain_separator.clone();
        state[1] = value;

        state.permute();
        state[1]
    }

    /// Hash two [`Fq`] elements with the provided `domain_separator`.
    fn old_hash_2(domain_separator: &Fq, value: (Fq, Fq)) -> Fq {
        let mut state = State::from(params::rate_2());

        // Use the domain separator as the sponge's capacity element
        state[0] = domain_separator.clone();
        state[1] = value.0;
        state[2] = value.1;

        state.permute();
        state[1]
    }

    /// Hash three [`Fq`] elements with the provided `domain_separator`.
    fn old_hash_3(domain_separator: &Fq, value: (Fq, Fq, Fq)) -> Fq {
        let mut state = State::from(params::rate_3());

        // Use the domain separator as the sponge's capacity element
        state[0] = domain_separator.clone();
        state[1] = value.0;
        state[2] = value.1;
        state[3] = value.2;

        state.permute();
        state[1]
    }

    /// Hash four [`Fq`] elements with the provided `domain_separator`.
    fn old_hash_4(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq)) -> Fq {
        let mut state = State::from(params::rate_4());

        // Use the domain separator as the sponge's capacity element
        state[0] = domain_separator.clone();
        state[1] = value.0;
        state[2] = value.1;
        state[3] = value.2;
        state[4] = value.3;

        state.permute();
        state[1]
    }

    /// Hash five [`Fq`] elements with the provided `domain_separator`.
    fn old_hash_5(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq)) -> Fq {
        let mut state = State::from(params::rate_5());

        // Use the domain separator as the sponge's capacity element
        state[0] = domain_separator.clone();
        state[1] = value.0;
        state[2] = value.1;
        state[3] = value.2;
        state[4] = value.3;
        state[5] = value.4;

        state.permute();
        state[1]
    }

    // These tests check we're not breaking the hashes in the wild that are using the old parameter set
    proptest! {
        #[test]
        fn consistent_hash_1(elem_1 in fq_strategy(), elem_2 in fq_strategy()) {
            let ark_result = old_hash_1(&elem_1, elem_2);
            let our_result = hash_1(&elem_1, elem_2);
            assert_eq!(ark_result, our_result);
        }

        #[test]
        fn consistent_hash_2(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy()) {
            let ark_result = old_hash_2(&elem_1, (elem_2, elem_3));
            let our_result = hash_2(&elem_1, (elem_2, elem_3));
            assert_eq!(ark_result, our_result);
        }

        #[test]
        fn consistent_hash_3(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy(),  elem_4 in fq_strategy()) {
            let ark_result = old_hash_3(&elem_1, (elem_2, elem_3, elem_4));
            let our_result = hash_3(&elem_1, (elem_2, elem_3, elem_4));
            assert_eq!(ark_result, our_result);
        }

        #[test]
        fn consistent_hash_4(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy(),  elem_4 in fq_strategy(), elem_5 in fq_strategy()) {
            let ark_result = old_hash_4(&elem_1, (elem_2, elem_3, elem_4, elem_5));
            let our_result = hash_4(&elem_1, (elem_2, elem_3, elem_4, elem_5));
            assert_eq!(ark_result, our_result);
        }

        #[test]
        fn consistent_hash_5(elem_1 in fq_strategy(), elem_2 in fq_strategy(), elem_3 in fq_strategy(),  elem_4 in fq_strategy(), elem_5 in fq_strategy(), elem_6 in fq_strategy()) {
            let ark_result = old_hash_5(&elem_1, (elem_2, elem_3, elem_4, elem_5, elem_6));
            let our_result = hash_5(&elem_1, (elem_2, elem_3, elem_4, elem_5, elem_6));
            assert_eq!(ark_result, our_result);
        }
    }
}
