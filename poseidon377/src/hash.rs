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

/// Hash six [`Fq`] elements with the provided `domain_separator`.
pub fn hash_6(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq, Fq)) -> Fq {
    let params = &crate::RATE_6_PARAMS;
    let mut state = Instance::new(params);
    state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
        value.4,
        value.5,
    ])
}
