//! An instantiation of Poseidon for the BLS12-377 scalar field.

use once_cell::sync::Lazy;

mod hash;
mod params;
/// Module containing optimized parameters. Temporary until we are ready
/// to make the next round of address-breaking changes.
mod temp_params;

pub use hash::{hash_1, hash_2, hash_3, hash_4, hash_5};

/// Parameters for the rate-1 instance of Poseidon.
pub const RATE_1_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(temp_params::rate_1);
/// Parameters for the rate-2 instance of Poseidon.
pub const RATE_2_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(temp_params::rate_2);
/// Parameters for the rate-3 instance of Poseidon.
pub const RATE_3_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(temp_params::rate_3);
/// Parameters for the rate-4 instance of Poseidon.
pub const RATE_4_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(temp_params::rate_4);
/// Parameters for the rate-5 instance of Poseidon.
pub const RATE_5_PARAMS: Lazy<PoseidonParameters<Fq>> = Lazy::new(temp_params::rate_5);

pub use ark_ed_on_bls12_377::Fq;
pub use poseidon_paramgen::PoseidonParameters;
pub use poseidon_permutation::Instance;
