use std::convert::TryInto;

use ark_ed_on_bls12_377::FqParameters;
use ark_ff::FpParameters;

use ark_sponge::poseidon::Parameters;
use poseidon_paramgen::PoseidonParameters;

use crate::{params, Fq};

/// Security level
pub const M: usize = 128;

/// Optimized parameters for the rate-1 instance of Poseidon.
pub fn rate_1() -> PoseidonParameters<Fq> {
    let ark_params: Parameters<Fq> = params::rate_1();
    PoseidonParameters::from_unoptimized_components(
        ark_params.alpha.try_into().unwrap(),
        ark_params.mds.clone(),
        ark_params.ark.clone(),
        M,
        2,
        FqParameters::MODULUS,
        ark_params.full_rounds,
        ark_params.partial_rounds,
    )
}

/// Optimized parameters for the rate-2 instance of Poseidon.
pub fn rate_2() -> PoseidonParameters<Fq> {
    let ark_params: Parameters<Fq> = params::rate_2();
    PoseidonParameters::from_unoptimized_components(
        ark_params.alpha.try_into().unwrap(),
        ark_params.mds,
        ark_params.ark,
        M,
        3,
        FqParameters::MODULUS,
        ark_params.full_rounds,
        ark_params.partial_rounds,
    )
}

/// Optimized parameters for the rate-3 instance of Poseidon.
pub fn rate_3() -> PoseidonParameters<Fq> {
    let ark_params: Parameters<Fq> = params::rate_3();
    PoseidonParameters::from_unoptimized_components(
        ark_params.alpha.try_into().unwrap(),
        ark_params.mds,
        ark_params.ark,
        M,
        4,
        FqParameters::MODULUS,
        ark_params.full_rounds,
        ark_params.partial_rounds,
    )
}

/// Optimized parameters for the rate-4 instance of Poseidon.
pub fn rate_4() -> PoseidonParameters<Fq> {
    let ark_params: Parameters<Fq> = params::rate_4();
    PoseidonParameters::from_unoptimized_components(
        ark_params.alpha.try_into().unwrap(),
        ark_params.mds,
        ark_params.ark,
        M,
        5,
        FqParameters::MODULUS,
        ark_params.full_rounds,
        ark_params.partial_rounds,
    )
}

/// Optimized parameters for the rate-5 instance of Poseidon.
pub fn rate_5() -> PoseidonParameters<Fq> {
    let ark_params: Parameters<Fq> = params::rate_5();
    PoseidonParameters::from_unoptimized_components(
        ark_params.alpha.try_into().unwrap(),
        ark_params.mds,
        ark_params.ark,
        M,
        6,
        FqParameters::MODULUS,
        ark_params.full_rounds,
        ark_params.partial_rounds,
    )
}
