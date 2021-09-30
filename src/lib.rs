mod rate_2;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        use ark_ed_on_bls12_377::Fq; // lazy import, fix
        use ark_ff::One;
        use ark_sponge::{
            poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge,
        };

        // how do we set the initial state ??
        // how do we set rate + capacity? settable in PoseidonSpongeVar but not PoseidonSponge ?
        // hardcoded here https://docs.rs/ark-sponge/0.3.0/src/ark_sponge/poseidon/mod.rs.html#239-240
        let mut sponge = PoseidonSponge::<Fq>::new(&rate_2::params());

        // comment suggests the impl hardcodes a different field ?
        // https://docs.rs/ark-sponge/0.3.0/src/ark_sponge/poseidon/mod.rs.html#230

        sponge.absorb(&Fq::one());
        sponge.absorb(&Fq::one());

        let output = sponge.squeeze_native_field_elements(1);
        dbg!(output);
        panic!();
    }
}
