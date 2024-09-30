#![no_main]
sp1_zkvm::entrypoint!(main);

use decaf377::Fq;

use poseidon377::{RATE_2_PARAMS, RATE_4_PARAMS};
use poseidon_permutation::Instance;
use rand_chacha::ChaChaRng;
use rand_core::{RngCore, SeedableRng};

#[sp1_derive::cycle_tracker]
fn hash_4_1_our_impl(i: &Fq, j: &Fq, k: &Fq, l: &Fq, m: &Fq) -> Fq {
    let mut our_instance = Instance::new(&RATE_4_PARAMS);
    our_instance.n_to_1_fixed_hash(&[*i, *j, *k, *l, *m])
}

#[sp1_derive::cycle_tracker]
fn hash_2_1_our_impl(i: &Fq, j: &Fq, k: &Fq) -> Fq {
    let mut our_instance = Instance::new(&RATE_2_PARAMS);
    our_instance.n_to_1_fixed_hash(&[*i, *j, *k])
}

pub fn main() {
    let n = 10;

    println!("cycle-tracker-start: setup");
    let mut rng = ChaChaRng::seed_from_u64(666);
    let mut test_field_elements = Vec::with_capacity(n);

    for _ in 0..n {
        let mut i_bytes = [0u8; 32];
        let mut j_bytes = [0u8; 32];
        let mut k_bytes = [0u8; 32];
        let mut l_bytes = [0u8; 32];
        let mut m_bytes = [0u8; 32];
        rng.fill_bytes(&mut i_bytes);
        rng.fill_bytes(&mut j_bytes);
        rng.fill_bytes(&mut k_bytes);
        rng.fill_bytes(&mut l_bytes);
        rng.fill_bytes(&mut m_bytes);
        test_field_elements.push((
            Fq::from_le_bytes_mod_order(&i_bytes[..]),
            Fq::from_le_bytes_mod_order(&j_bytes[..]),
            Fq::from_le_bytes_mod_order(&k_bytes[..]),
            Fq::from_le_bytes_mod_order(&l_bytes[..]),
            Fq::from_le_bytes_mod_order(&m_bytes[..]),
        ))
    }
    println!("cycle-tracker-end: setup");

    println!("cycle-tracker-start: main-body");
    for (i, j, k, l, m) in test_field_elements.iter() {
        let result = hash_4_1_our_impl(i, j, k, l, m);
        println!("result: {}", result);
    }
    println!("cycle-tracker-end: main-body");

    println!("cycle-tracker-start: main-body");
    for (i, j, k, _, _) in test_field_elements.iter() {
        let result = hash_2_1_our_impl(i, j, k);
        println!("result: {}", result);
    }
    println!("cycle-tracker-end: main-body");
}
