//! An instantiation of Poseidon for the BLS12-377 scalar field.
#![cfg_attr(not(feature = "std"), no_std)]

mod hash;
mod params;

// pub use hash::{hash_1, hash_2, hash_3, hash_4, hash_5, hash_6, hash_7};
pub use hash::{hash_1, hash_2, hash_3};

/// Parameters for the rate-1 instance of Poseidon.
pub const RATE_1_PARAMS: PoseidonParameters<2, 1, 4, 1, 39, 2, 78, 31> = params::rate_1::rate_1();

/// Parameters for the rate-2 instance of Poseidon.
pub const RATE_2_PARAMS: PoseidonParameters<3, 2, 9, 4, 39, 3, 117, 31> = params::rate_2::rate_2();

/// Parameters for the rate-3 instance of Poseidon.
pub const RATE_3_PARAMS: PoseidonParameters<4, 3, 16, 9, 39, 4, 156, 31> = params::rate_3::rate_3();

// /// Parameters for the rate-4 instance of Poseidon.
// pub static RATE_4_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_4);
// /// Parameters for the rate-5 instance of Poseidon.
// pub static RATE_5_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_5);
// /// Parameters for the rate-6 instance of Poseidon.
// pub static RATE_6_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_6);
// /// Parameters for the rate-7 instance of Poseidon.
// pub static RATE_7_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(params::rate_7);
pub use decaf377::Fq;
pub use poseidon_parameters::v1::PoseidonParameters;
pub use poseidon_permutation::Instance;

// #[cfg(feature = "r1cs")]
// pub mod r1cs;
