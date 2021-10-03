mod sponge;

pub mod params;

// Since we depend on a git version of ark-sponge, re-exporting it here means
// our deps can access it without having to keep git revisions in sync.
//
// Going forward, this re-export should be removed and the functionality our
// deps need from direct use of ark-sponge should be folded into this crate.
// However, it's faster to iterate on required functionality without imposing
// hard compartmentalization boundaries from the start.
pub use ark_sponge;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        use ark_ed_on_bls12_377::Fq; // lazy import, fix
        use ark_ff::{One, Zero};
        use ark_sponge::{
            poseidon::PoseidonSponge, CryptographicSponge, DuplexSpongeMode,
            FieldBasedCryptographicSponge,
        };

        // Current API has a `new()` method as part of the `CryptographicSponge`
        // trait, but this method doesn't allow setting the initial state
        // manually.  Instead, the fields can be set manually.
        // Slightly inconvenient that we have to initialize the mode.
        let mut sponge = PoseidonSponge {
            parameters: params::rate_2(),
            state: vec![Fq::zero(); 3],
            mode: DuplexSpongeMode::Absorbing {
                next_absorb_index: 0,
            },
        };

        sponge.absorb(&Fq::one());
        sponge.absorb(&Fq::one());

        let output = sponge.squeeze_native_field_elements(1);
        dbg!(output);
    }
}
