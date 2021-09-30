mod rate_2;
mod sponge;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        use ark_ed_on_bls12_377::Fq; // lazy import, fix
        use ark_ff::One;
        use ark_sponge::{
            poseidon::PoseidonSponge, CryptographicSponge, FieldBasedCryptographicSponge, SpongeExt,
        };

        // how do we set rate + capacity? settable in PoseidonSpongeVar but not PoseidonSponge ?
        // hardcoded here https://docs.rs/ark-sponge/0.3.0/src/ark_sponge/poseidon/mod.rs.html#239-240
        let mut sponge = PoseidonSponge::<Fq>::new(&rate_2::params());

        // How do we set the initial state ??
        // Hmm.. seems like the `SpongeExt` trait does not let us directly set the state:
        let state = vec![Fq::one()]; // Can't instantiate from this.
        let poseidon_state = sponge.into_state();

        // At this point one cannot directly modify the state... but can create a new sponge from this `PoseidonSpongeState` type:
        let sponge_from_state = PoseidonSponge::from_state(poseidon_state, &rate_2::params());

        let mut sponge = PoseidonSponge::<Fq>::new(&rate_2::params());

        // comment suggests the impl hardcodes a different field ?
        // https://docs.rs/ark-sponge/0.3.0/src/ark_sponge/poseidon/mod.rs.html#230

        sponge.absorb(&Fq::one());
        sponge.absorb(&Fq::one());

        let output = sponge.squeeze_native_field_elements(1);
        dbg!(output);
    }
}
