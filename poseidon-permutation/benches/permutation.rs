use ark_ed_on_bls12_377::{Fq, FqParameters};
use ark_ff::{FpParameters, PrimeField};
use ark_sponge::poseidon::{Parameters, State};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use once_cell::sync::Lazy;
use rand_chacha::ChaChaRng;
use rand_core::{RngCore, SeedableRng};

use poseidon_paramgen::PoseidonParameters;
use poseidon_permutation::Instance;

static PARAMS_4_TO_1: Lazy<PoseidonParameters<Fq>> =
    Lazy::new(|| PoseidonParameters::<Fq>::new(128, 5, FqParameters::MODULUS, true));

fn hash_4_1_ark_sponge(i: &Fq, j: &Fq, k: &Fq, l: &Fq, m: &Fq) -> Fq {
    let params_ark: Parameters<Fq> = (PARAMS_4_TO_1.clone()).into();
    let mut ark_state = State::from(params_ark);
    ark_state[0] = *i;
    ark_state[1] = *j;
    ark_state[2] = *k;
    ark_state[3] = *l;
    ark_state[4] = *m;
    ark_state.permute();

    ark_state[1]
}

fn hash_4_1_our_impl(i: &Fq, j: &Fq, k: &Fq, l: &Fq, m: &Fq) -> Fq {
    let mut our_instance = Instance::new(&PARAMS_4_TO_1);
    our_instance.n_to_1_fixed_hash(vec![*i, *j, *k, *l, *m])
}

fn hash_4_1_our_impl_unoptimized(i: &Fq, j: &Fq, k: &Fq, l: &Fq, m: &Fq) -> Fq {
    let mut our_instance = Instance::new(&PARAMS_4_TO_1);
    our_instance.unoptimized_n_to_1_fixed_hash(vec![*i, *j, *k, *l, *m])
}

pub fn bench_ark_sponge_vs_optimized(c: &mut Criterion) {
    let mut group = c.benchmark_group("ark_sponge_vs_optimized");
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

pub fn bench_ark_sponge_vs_unoptimized(c: &mut Criterion) {
    let mut group = c.benchmark_group("ark_sponge_vs_unoptimized");
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
            BenchmarkId::new("unoptimized", format!("unoptimized: {}", index)),
            &(&i, &j, &k, &l, &m),
            |b, (i, j, k, l, m)| b.iter(|| hash_4_1_our_impl_unoptimized(i, j, k, l, m)),
        );
    }
    group.finish();
}

pub fn bench_unoptimized_vs_optimized(c: &mut Criterion) {
    let mut group = c.benchmark_group("unoptimized_vs_optimized");
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
            BenchmarkId::new("base_alg", format!("base_alg: {}", index)),
            &(&i, &j, &k, &l, &m),
            |b, (i, j, k, l, m)| b.iter(|| hash_4_1_our_impl_unoptimized(i, j, k, l, m)),
        );
        group.bench_with_input(
            BenchmarkId::new("optimized_alg", format!("optimized_alg: {}", index)),
            &(&i, &j, &k, &l, &m),
            |b, (i, j, k, l, m)| b.iter(|| hash_4_1_our_impl(i, j, k, l, m)),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_ark_sponge_vs_optimized);
//criterion_group!(benches, bench_ark_sponge_vs_unoptimized);
//criterion_group!(benches, bench_unoptimized_vs_optimized);
criterion_main!(benches);
