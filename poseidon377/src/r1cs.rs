use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use ark_sponge::{constraints::CryptographicSpongeVar, poseidon::constraints::PoseidonSpongeVar};

use ark_ed_on_bls12_377::constraints::FqVar;

use crate::Fq;

pub fn hash_1(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: FqVar,
) -> Result<FqVar, SynthesisError> {
    todo!()
}
