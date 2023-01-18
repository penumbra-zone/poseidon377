//! An instantiation of Poseidon for the BLS12-377 scalar field.
#![cfg_attr(not(feature = "std"), no_std)]

use once_cell::sync::Lazy;

mod hash;

// We want the full parameters to be available in the documentation.
#[cfg(doc)]
pub mod params {
    use ark_ff::vec;
    include!(concat!(env!("OUT_DIR"), "/params.rs"));
}
#[cfg(not(doc))]
mod params {
    use ark_ff::vec;
    include!(concat!(env!("OUT_DIR"), "/params.rs"));
}

pub use hash::{hash_1, hash_2, hash_3, hash_4, hash_5, hash_6, hash_7};

/// Parameters for the rate-1 instance of Poseidon.
pub static RATE_1_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_1);
/// Parameters for the rate-2 instance of Poseidon.
pub static RATE_2_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_2);
/// Parameters for the rate-3 instance of Poseidon.
pub static RATE_3_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_3);
/// Parameters for the rate-4 instance of Poseidon.
pub static RATE_4_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_4);
/// Parameters for the rate-5 instance of Poseidon.
pub static RATE_5_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_5);
/// Parameters for the rate-6 instance of Poseidon.
pub static RATE_6_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_6);
/// Parameters for the rate-7 instance of Poseidon.
pub static RATE_7_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_7);

pub use decaf377::Fq;
pub use poseidon_paramgen::PoseidonParameters;
pub use poseidon_permutation::Instance;

#[cfg(feature = "r1cs")]
pub mod r1cs;
