use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use ark_sponge::{constraints::CryptographicSpongeVar, poseidon::constraints::PoseidonSpongeVar};

use decaf377::r1cs::FqVar;

use crate::Fq;

pub fn hash_1(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: FqVar,
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_1_PARAMS).clone();
    let ark_params = (params).into();

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![domain_separator, &value])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}
