use ark_ff::{One, PrimeField};
use ark_groth16::{Groth16, ProvingKey, VerifyingKey};
use once_cell::sync::Lazy;
use proptest::prelude::*;

use ark_r1cs_std::prelude::{AllocVar, EqGadget};
use ark_relations::r1cs::{ConstraintSynthesizer, ToConstraintField};
use ark_snark::SNARK;
use decaf377::{
    r1cs::{CountConstraints, FqVar},
    Bls12_377, Fq,
};
use rand_core::OsRng;

// This is a domain separator we'll use as a constant in our circuit below.
static DOMAIN_SEP: Lazy<Fq> = Lazy::new(|| Fq::from(666));

#[derive(Clone)]
struct PreimageCircuit {
    // Witness
    preimage: Fq,

    // Public input
    pub hash_output: Fq,
}

impl ConstraintSynthesizer<Fq> for PreimageCircuit {
    fn generate_constraints(
        self,
        cs: ark_relations::r1cs::ConstraintSystemRef<Fq>,
    ) -> ark_relations::r1cs::Result<()> {
        let preimage_var = FqVar::new_witness(cs.clone(), || Ok(&self.preimage))?;

        let domain_separator_var = FqVar::new_constant(cs.clone(), *DOMAIN_SEP)?;
        let hash_output_var = FqVar::new_input(cs.clone(), || Ok(self.hash_output))?;

        let test_hash_output = poseidon377::r1cs::hash_1(cs, &domain_separator_var, preimage_var)?;
        hash_output_var.enforce_equal(&test_hash_output)?;

        Ok(())
    }
}

impl PreimageCircuit {
    fn generate_test_parameters() -> (ProvingKey<Bls12_377>, VerifyingKey<Bls12_377>) {
        let circuit = PreimageCircuit {
            preimage: Fq::from(2),
            hash_output: Fq::from(2),
        };
        let (pk, vk) = Groth16::circuit_specific_setup(circuit, &mut OsRng)
            .expect("can perform circuit specific setup");
        (pk, vk)
    }
}

fn fq_strategy() -> BoxedStrategy<Fq> {
    any::<[u8; 32]>()
        .prop_map(|bytes| Fq::from_le_bytes_mod_order(&bytes[..]))
        .boxed()
}

proptest! {
#![proptest_config(ProptestConfig::with_cases(5))]
#[test]
fn groth16_hash_proof_happy_path(preimage in fq_strategy()) {
        let (pk, vk) = PreimageCircuit::generate_test_parameters();
        let mut rng = OsRng;

        let hash_output = poseidon377::hash_1(&DOMAIN_SEP, preimage);

        // Prover POV
        let circuit = PreimageCircuit { preimage,hash_output };
        let proof = Groth16::prove(&pk, circuit.clone(), &mut rng)
            .expect("can generate proof");
        dbg!(circuit.clone().num_constraints_and_instance_variables());

        // Verifier POV
        let processed_pvk = Groth16::process_vk(&vk).expect("can process verifying key");
        let public_inputs = hash_output.to_field_elements().unwrap();
        let proof_result =
            Groth16::verify_with_processed_vk(&processed_pvk, &public_inputs, &proof).unwrap();

        assert!(proof_result);
    }
}

proptest! {
#![proptest_config(ProptestConfig::with_cases(5))]
#[test]
fn groth16_hash_proof_unhappy_path(preimage in fq_strategy()) {
        let (pk, vk) = PreimageCircuit::generate_test_parameters();
        let mut rng = OsRng;

        let hash_output = poseidon377::hash_1(&DOMAIN_SEP, preimage);

        // Prover POV
        let circuit = PreimageCircuit { preimage,hash_output };
        let proof = Groth16::prove(&pk, circuit.clone(), &mut rng)
            .expect("can generate proof");
        dbg!(circuit.clone().num_constraints_and_instance_variables());

        // Verifier POV
        let processed_pvk = Groth16::process_vk(&vk).expect("can process verifying key");
        let public_inputs = (hash_output + Fq::one()).to_field_elements().unwrap();
        let proof_result =
            Groth16::verify_with_processed_vk(&processed_pvk, &public_inputs, &proof).unwrap();

        assert!(!proof_result);
    }
}
