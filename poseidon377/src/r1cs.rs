use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use decaf377::r1cs::FqVar;
use poseidon_permutation::r1cs::InstanceVar;

use crate::Fq;

pub fn hash_1(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: FqVar,
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_1_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![domain_separator.clone(), value]))
}

pub fn hash_2(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_2_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![domain_separator.clone(), value.0, value.1]))
}

pub fn hash_3(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_3_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![domain_separator.clone(), value.0, value.1, value.2]))
}

pub fn hash_4(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_4_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
    ]))
}

pub fn hash_5(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_5_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
        value.4,
    ]))
}

pub fn hash_6(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_6_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
        value.4,
        value.5,
    ]))
}

pub fn hash_7(
    cs: ConstraintSystemRef<Fq>,
    domain_separator: &FqVar,
    value: (FqVar, FqVar, FqVar, FqVar, FqVar, FqVar, FqVar),
) -> Result<FqVar, SynthesisError> {
    let mut state = InstanceVar::new(crate::RATE_7_PARAMS.clone(), cs.clone());
    Ok(state.n_to_1_fixed_hash(vec![
        domain_separator.clone(),
        value.0,
        value.1,
        value.2,
        value.3,
        value.4,
        value.5,
        value.6,
    ]))
}
