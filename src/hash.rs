use crate::{Fq, State};

/// Hash a single [`Fq`] element with the provided `domain_separator`.
pub fn hash_1(domain_separator: &Fq, value: Fq) -> Fq {
    let mut state = State::from(crate::RATE_1_PARAMS.clone());

    // Use the domain separator as the sponge's capacity element
    state[0] = domain_separator.clone();
    state[1] = value;

    state.permute();
    state[1]
}

/// Hash two [`Fq`] elements with the provided `domain_separator`.
pub fn hash_2(domain_separator: &Fq, value: (Fq, Fq)) -> Fq {
    let mut state = State::from(crate::RATE_2_PARAMS.clone());

    // Use the domain separator as the sponge's capacity element
    state[0] = domain_separator.clone();
    state[1] = value.0;
    state[2] = value.1;

    state.permute();
    state[1]
}

/// Hash three [`Fq`] elements with the provided `domain_separator`.
pub fn hash_3(domain_separator: &Fq, value: (Fq, Fq, Fq)) -> Fq {
    let mut state = State::from(crate::RATE_3_PARAMS.clone());

    // Use the domain separator as the sponge's capacity element
    state[0] = domain_separator.clone();
    state[1] = value.0;
    state[2] = value.1;
    state[3] = value.2;

    state.permute();
    state[1]
}

/// Hash four [`Fq`] elements with the provided `domain_separator`.
pub fn hash_4(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq)) -> Fq {
    let mut state = State::from(crate::RATE_4_PARAMS.clone());

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
pub fn hash_5(domain_separator: &Fq, value: (Fq, Fq, Fq, Fq, Fq)) -> Fq {
    let mut state = State::from(crate::RATE_5_PARAMS.clone());

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
