use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use decaf377::r1cs::FqVar;
use poseidon_parameters::{Alpha, PoseidonParameters};

use crate::Fq;

pub mod vendor;

use vendor::sponge::{
    constraints::CryptographicSpongeVar,
    poseidon::{constraints::PoseidonSpongeVar, PoseidonParameters as ArkPoseidonParameters},
};

fn convert_to_ark_sponge_parameters(params: PoseidonParameters<Fq>) -> ArkPoseidonParameters<Fq> {
    let alpha = match params.alpha {
        Alpha::Exponent(exp) => exp as u64,
        Alpha::Inverse => panic!("ark-sponge does not allow inverse alpha"),
    };
    // TODO: let user specify different capacity choices
    let capacity = 1;
    let rate = params.t - capacity;
    let full_rounds = params.rounds.full();
    let partial_rounds = params.rounds.partial();

    ArkPoseidonParameters {
        full_rounds,
        partial_rounds,
        alpha,
        ark: params.arc.into(),
        mds: params.mds.into(),
        rate,
        capacity,
    }
}

pub fn hash_1(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: FqVar,
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_1_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![domain_separator, &value])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}

pub fn hash_2(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_2_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![domain_separator, &value.0, &value.1])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}

pub fn hash_3(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_3_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![domain_separator, &value.0, &value.1, &value.2])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}

pub fn hash_4(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_4_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![
        domain_separator,
        &value.0,
        &value.1,
        &value.2,
        &value.3,
    ])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}

pub fn hash_5(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_5_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![
        domain_separator,
        &value.0,
        &value.1,
        &value.2,
        &value.3,
        &value.4,
    ])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}

pub fn hash_6(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_6_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![
        domain_separator,
        &value.0,
        &value.1,
        &value.2,
        &value.3,
        &value.4,
        &value.5,
    ])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}

pub fn hash_7(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let params = (*crate::RATE_7_PARAMS).clone();
    let ark_params = convert_to_ark_sponge_parameters(params);

    let mut poseidon_instance: PoseidonSpongeVar<Fq> = PoseidonSpongeVar::new(cs, &ark_params);
    poseidon_instance.absorb(&vec![
        domain_separator,
        &value.0,
        &value.1,
        &value.2,
        &value.3,
        &value.4,
        &value.5,
        &value.6,
    ])?;
    let output = poseidon_instance.squeeze_field_elements(1)?;
    Ok(output[0].clone())
}
