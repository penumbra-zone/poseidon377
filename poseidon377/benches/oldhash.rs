use ark_ed_on_bls12_377::Fq;
use ark_ff::PrimeField;
use ark_sponge::poseidon::State;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand_chacha::ChaChaRng;
use rand_core::{RngCore, SeedableRng};

use poseidon377::{hash_4, params};

fn hash_4_1_ark_sponge(i: &Fq, j: &Fq, k: &Fq, l: &Fq, m: &Fq) -> Fq {
    let mut state = State::from(params::rate_4());

    // Use the domain separator as the sponge's capacity element
    state[0] = *i;
    state[1] = *j;
    state[2] = *k;
    state[3] = *l;
    state[4] = *m;

    state.permute();
    state[1]
}

fn hash_4_1_our_impl(i: &Fq, j: &Fq, k: &Fq, l: &Fq, m: &Fq) -> Fq {
    hash_4(i, (*j, *k, *l, *m))
}

pub fn poseidon377_ark_sponge_vs_optimized(c: &mut Criterion) {
    let mut group = c.benchmark_group("poseidon377_ark_sponge_vs_optimized");
    let n = 10;
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

    for (index, (i, j, k, l, m)) in test_field_elements.iter().enumerate() {
        group.bench_with_input(
            BenchmarkId::new("ark-sponge", format!("ark-sponge: {}", index)),
            &(&i, &j, &k, &l, &m),
            |b, (i, j, k, l, m)| b.iter(|| hash_4_1_ark_sponge(i, j, k, l, m)),
        );
        group.bench_with_input(
            BenchmarkId::new(
                "poseidon-permutation",
                format!("poseidon-permutation: {}", index),
            ),
            &(&i, &j, &k, &l, &m),
            |b, (i, j, k, l, m)| b.iter(|| hash_4_1_our_impl(i, j, k, l, m)),
        );
    }
    group.finish();
}

criterion_group!(benches, poseidon377_ark_sponge_vs_optimized);
criterion_main!(benches);
