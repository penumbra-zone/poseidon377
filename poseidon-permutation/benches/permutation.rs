use ark_ed_on_bls12_377::{Fq, FqParameters};
use ark_ff::{FpParameters, PrimeField};
use ark_sponge::poseidon::{Parameters, State};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand_chacha::ChaChaRng;
use rand_core::{RngCore, SeedableRng};

use poseidon_paramgen::PoseidonParameters;
use poseidon_permutation::Instance;

fn hash_2_1_ark_sponge(i: &Fq, j: &Fq, k: &Fq) -> Fq {
    let params_2_to_1 = PoseidonParameters::<Fq>::new(128, 3, FqParameters::MODULUS, true);

    let params_ark: Parameters<Fq> = params_2_to_1.into();
    let mut ark_state = State::from(params_ark);
    ark_state[0] = *i;
    ark_state[1] = *j;
    ark_state[2] = *k;
    ark_state.permute();

    ark_state[1]
}

fn hash_2_1_our_impl(i: &Fq, j: &Fq, k: &Fq) -> Fq {
    let params_2_to_1 = PoseidonParameters::<Fq>::new(128, 3, FqParameters::MODULUS, true);

    let mut our_instance = Instance::new(params_2_to_1);
    our_instance.n_to_1_fixed_hash(vec![*i, *j, *k])
}

pub fn bench_ark_sponge_vs_optimized(c: &mut Criterion) {
    let mut group = c.benchmark_group("ark_sponge_vs_unoptimized");
    let n = 10;
    let mut rng = ChaChaRng::seed_from_u64(666);
    let mut test_field_elements = Vec::with_capacity(n);

    for _ in 0..n {
        let mut i_bytes = [0u8; 32];
        let mut j_bytes = [0u8; 32];
        let mut k_bytes = [0u8; 32];
        rng.fill_bytes(&mut i_bytes);
        rng.fill_bytes(&mut j_bytes);
        rng.fill_bytes(&mut k_bytes);
        test_field_elements.push((
            Fq::from_le_bytes_mod_order(&i_bytes[..]),
            Fq::from_le_bytes_mod_order(&j_bytes[..]),
            Fq::from_le_bytes_mod_order(&k_bytes[..]),
        ))
    }

    for (i, j, k) in test_field_elements {
        group.bench_with_input(
            BenchmarkId::new("ark-sponge", format!("{}/{}/{}", i, j, k)),
            &(&i, &j, &k),
            |b, (i, j, k)| b.iter(|| hash_2_1_ark_sponge(i, j, k)),
        );
        group.bench_with_input(
            BenchmarkId::new("poseidon-permutation", format!("{}/{}", i, j)),
            &(&i, &j, &k),
            |b, (i, j, k)| b.iter(|| hash_2_1_our_impl(i, j, k)),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_ark_sponge_vs_optimized);
criterion_main!(benches);
