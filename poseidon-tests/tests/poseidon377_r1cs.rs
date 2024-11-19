use std::{fs, io::BufWriter, path::PathBuf};

use ark_groth16::{r1cs_to_qap::LibsnarkReduction, Groth16, ProvingKey, VerifyingKey};
use ark_r1cs_std::prelude::{AllocVar, EqGadget};
use ark_relations::r1cs::{ConstraintSynthesizer, ToConstraintField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_snark::SNARK;
use decaf377::{
    r1cs::{CountConstraints, FqVar},
    Bls12_377, Fq,
};
use once_cell::sync::Lazy;
use proptest::prelude::*;
use rand_core::OsRng;

// This is a domain separator we'll use as a constant in our circuits below.
const DOMAIN_SEP: Fq = Fq::from_montgomery_limbs([
    15545862963729521748,
    8301038308012579564,
    6940338389606830037,
    431207966161902270,
]);

/// The maximum fixed-width Poseidon hash exposed to downstream users of this crate.
const MAX_WIDTH_POSEIDON_HASH: usize = 7;

static PREIMAGE_PK: Lazy<ProvingKey<Bls12_377>> = Lazy::new(|| {
    let pk_bytes = include_bytes!("../../poseidon-tests/test_vectors/preimage_pk.bin");
    ProvingKey::deserialize_uncompressed(&pk_bytes[..])
        .expect("can parse public element input proving key")
});

static PREIMAGE_VK: Lazy<VerifyingKey<Bls12_377>> = Lazy::new(|| {
    let vk_bytes = include_bytes!("../../poseidon-tests/test_vectors/preimage_vk.param");
    VerifyingKey::deserialize_uncompressed(&vk_bytes[..])
        .expect("can parse public element input verifying key")
});

#[derive(Clone)]
struct PreimageCircuit {
    // Witnesses
    preimages: [Fq; MAX_WIDTH_POSEIDON_HASH],

    // Public inputs
    //
    // The expected public inputs are the result of
    // hashing the above witness values, i.e.:
    //
    // hash_output[0] = hash_1(preimage[0])
    // hash_output[1] = hash_2(preimage[0], preimage[1])
    // ..
    // hash_output[6] = hash_7(preimage[0], preimage[1], .., preimage[6])
    //
    // This ensures that this test circuit covers all fixed-width hashes we expose
    // for use in Penumbra circuits.
    pub hash_outputs: [Fq; MAX_WIDTH_POSEIDON_HASH],
}

impl ConstraintSynthesizer<Fq> for PreimageCircuit {
    fn generate_constraints(
        self,
        cs: ark_relations::r1cs::ConstraintSystemRef<Fq>,
    ) -> ark_relations::r1cs::Result<()> {
        // Add all witnesses
        let mut preimage_vars = Vec::new();
        for value in self.preimages {
            let preimage_var = FqVar::new_witness(cs.clone(), || Ok(&value))?;
            preimage_vars.push(preimage_var);
        }

        let domain_separator_var = FqVar::new_constant(cs.clone(), DOMAIN_SEP)?;
        // Add all public inputs
        let mut hash_output_vars = Vec::new();
        for hash_output in self.hash_outputs {
            let hash_output_var = FqVar::new_input(cs.clone(), || Ok(hash_output))?;
            hash_output_vars.push(hash_output_var);
        }

        // Add constraints to check the hash outputs are the result of computing the fixed-width
        // hashes of the corresponding inputs.
        let test_hash1_output =
            poseidon377::r1cs::hash_1(cs.clone(), &domain_separator_var, preimage_vars[0].clone())?;
        hash_output_vars[0].enforce_equal(&test_hash1_output)?;

        let test_hash2_output = poseidon377::r1cs::hash_2(
            cs.clone(),
            &domain_separator_var,
            (preimage_vars[0].clone(), preimage_vars[1].clone()),
        )?;
        hash_output_vars[1].enforce_equal(&test_hash2_output)?;

        let test_hash3_output = poseidon377::r1cs::hash_3(
            cs.clone(),
            &domain_separator_var,
            (
                preimage_vars[0].clone(),
                preimage_vars[1].clone(),
                preimage_vars[2].clone(),
            ),
        )?;
        hash_output_vars[2].enforce_equal(&test_hash3_output)?;

        let test_hash4_output = poseidon377::r1cs::hash_4(
            cs.clone(),
            &domain_separator_var,
            (
                preimage_vars[0].clone(),
                preimage_vars[1].clone(),
                preimage_vars[2].clone(),
                preimage_vars[3].clone(),
            ),
        )?;
        hash_output_vars[3].enforce_equal(&test_hash4_output)?;

        let test_hash5_output = poseidon377::r1cs::hash_5(
            cs.clone(),
            &domain_separator_var,
            (
                preimage_vars[0].clone(),
                preimage_vars[1].clone(),
                preimage_vars[2].clone(),
                preimage_vars[3].clone(),
                preimage_vars[4].clone(),
            ),
        )?;
        hash_output_vars[4].enforce_equal(&test_hash5_output)?;

        let test_hash6_output = poseidon377::r1cs::hash_6(
            cs.clone(),
            &domain_separator_var,
            (
                preimage_vars[0].clone(),
                preimage_vars[1].clone(),
                preimage_vars[2].clone(),
                preimage_vars[3].clone(),
                preimage_vars[4].clone(),
                preimage_vars[5].clone(),
            ),
        )?;
        hash_output_vars[5].enforce_equal(&test_hash6_output)?;

        let test_hash7_output = poseidon377::r1cs::hash_7(
            cs,
            &domain_separator_var,
            (
                preimage_vars[0].clone(),
                preimage_vars[1].clone(),
                preimage_vars[2].clone(),
                preimage_vars[3].clone(),
                preimage_vars[4].clone(),
                preimage_vars[5].clone(),
                preimage_vars[6].clone(),
            ),
        )?;
        hash_output_vars[6].enforce_equal(&test_hash7_output)?;
        Ok(())
    }
}

impl PreimageCircuit {
    fn generate_test_parameters() -> (ProvingKey<Bls12_377>, VerifyingKey<Bls12_377>) {
        let circuit = PreimageCircuit {
            preimages: [Fq::from(2u64); MAX_WIDTH_POSEIDON_HASH],
            hash_outputs: [Fq::from(2u64); MAX_WIDTH_POSEIDON_HASH],
        };
        let (pk, vk) =
            Groth16::<Bls12_377, LibsnarkReduction>::circuit_specific_setup(circuit, &mut OsRng)
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
#![proptest_config(ProptestConfig::with_cases(1))]
#[test]
fn groth16_hash_proof_happy_path(v1 in fq_strategy(), v2 in fq_strategy(), v3 in fq_strategy(), v4 in fq_strategy(), v5 in fq_strategy(), v6 in fq_strategy(), v7 in fq_strategy()) {
        let pk = PREIMAGE_PK.clone();
        let vk = PREIMAGE_VK.clone();
        let mut rng = OsRng;

        let preimages = [v1, v2, v3, v4, v5, v6, v7];
        let mut hash_outputs = [Fq::from(0u64); MAX_WIDTH_POSEIDON_HASH];
        hash_outputs[0] = poseidon377::hash_1(&DOMAIN_SEP, preimages[0]);
        hash_outputs[1] = poseidon377::hash_2(&DOMAIN_SEP, (preimages[0], preimages[1]));
        hash_outputs[2] = poseidon377::hash_3(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2]));
        hash_outputs[3] = poseidon377::hash_4(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3]));
        hash_outputs[4] = poseidon377::hash_5(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3], preimages[4]));
        hash_outputs[5] = poseidon377::hash_6(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3], preimages[4], preimages[5]));
        hash_outputs[6] = poseidon377::hash_7(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3], preimages[4], preimages[5], preimages[6]));

        // Prover POV
        let circuit = PreimageCircuit { preimages, hash_outputs };
        let proof = Groth16::<Bls12_377, LibsnarkReduction>::prove(&pk, circuit.clone(), &mut rng)
            .expect("can generate proof");
        dbg!(circuit.clone().num_constraints_and_instance_variables());

        // Verifier POV
        let processed_pvk = Groth16::<Bls12_377, LibsnarkReduction>::process_vk(&vk).expect("can process verifying key");
        let public_inputs = hash_outputs.to_field_elements().unwrap();
        let proof_result =
            Groth16::<Bls12_377, LibsnarkReduction>::verify_with_processed_vk(&processed_pvk, &public_inputs, &proof).unwrap();

        assert!(proof_result);
    }
}

proptest! {
#![proptest_config(ProptestConfig::with_cases(5))]
#[test]
fn groth16_hash_proof_unhappy_path(v1 in fq_strategy(), v2 in fq_strategy(), v3 in fq_strategy(), v4 in fq_strategy(), v5 in fq_strategy(), v6 in fq_strategy(), v7 in fq_strategy()) {
        let pk = PREIMAGE_PK.clone();
        let vk = PREIMAGE_VK.clone();
        let mut rng = OsRng;

        let preimages = [v1, v2, v3, v4, v5, v6, v7];
        let mut hash_outputs = [Fq::from(0u64); MAX_WIDTH_POSEIDON_HASH];
        hash_outputs[0] = poseidon377::hash_1(&DOMAIN_SEP, preimages[0]);
        hash_outputs[1] = poseidon377::hash_2(&DOMAIN_SEP, (preimages[0], preimages[1]));
        hash_outputs[2] = poseidon377::hash_3(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2]));
        hash_outputs[3] = poseidon377::hash_4(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3]));
        hash_outputs[4] = poseidon377::hash_5(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3], preimages[4]));
        hash_outputs[5] = poseidon377::hash_6(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3], preimages[4], preimages[5]));
        hash_outputs[6] = poseidon377::hash_7(&DOMAIN_SEP, (preimages[0], preimages[1], preimages[2], preimages[3], preimages[4], preimages[5], preimages[6]));

        // Prover POV
        let circuit = PreimageCircuit { preimages, hash_outputs };
        let proof = Groth16::<Bls12_377, LibsnarkReduction>::prove(&pk, circuit.clone(), &mut rng)
            .expect("can generate proof");
        dbg!(circuit.clone().num_constraints_and_instance_variables());

        // Verifier POV
        let processed_pvk = Groth16::<Bls12_377, LibsnarkReduction>::process_vk(&vk).expect("can process verifying key");
        let public_inputs = [(hash_outputs[0] + Fq::from(1u64)); MAX_WIDTH_POSEIDON_HASH].to_field_elements().unwrap();
        let proof_result =
            Groth16::<Bls12_377, LibsnarkReduction>::verify_with_processed_vk(&processed_pvk, &public_inputs, &proof).unwrap();

        assert!(!proof_result);
    }
}

fn write_params(
    target_dir: &PathBuf,
    name: &str,
    pk: &ProvingKey<Bls12_377>,
    vk: &VerifyingKey<Bls12_377>,
) -> anyhow::Result<()> {
    let pk_location = target_dir.join(format!("{}_pk.bin", name));
    let vk_location = target_dir.join(format!("{}_vk.param", name));

    let pk_file = fs::File::create(&pk_location)?;
    let vk_file = fs::File::create(&vk_location)?;

    let pk_writer = BufWriter::new(pk_file);
    let vk_writer = BufWriter::new(vk_file);

    ProvingKey::serialize_uncompressed(pk, pk_writer).expect("can serialize ProvingKey");
    VerifyingKey::serialize_uncompressed(vk, vk_writer).expect("can serialize VerifyingKey");

    Ok(())
}

#[ignore]
#[test]
fn generate_test_vectors() {
    let (pk, vk) = PreimageCircuit::generate_test_parameters();
    write_params(&PathBuf::from("test_vectors"), "preimage", &pk, &vk)
        .expect("can write test vectors");
}
