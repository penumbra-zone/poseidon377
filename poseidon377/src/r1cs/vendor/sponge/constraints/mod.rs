use crate::r1cs::vendor::sponge::{Absorb, CryptographicSponge, FieldElementSize};
use ark_ff::PrimeField;
// use ark_nonnative_field::params::{get_params, OptimizationType};
// use ark_nonnative_field::{AllocatedNonNativeFieldVar, NonNativeFieldVar};
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::bits::uint8::UInt8;
use ark_r1cs_std::fields::fp::{AllocatedFp, FpVar};
use ark_r1cs_std::R1CSVar;
use ark_relations::lc;
use ark_relations::r1cs::{ConstraintSystemRef, LinearCombination, SynthesisError};
use ark_std::vec;
use ark_std::vec::Vec;

mod absorb;
pub use absorb::*;

/// Enables simple access to the "gadget" version of the sponge.
/// Simplifies trait bounds in downstream generic code.
pub trait SpongeWithGadget<CF: PrimeField>: CryptographicSponge {
    /// The gadget version of `Self`.
    type Var: CryptographicSpongeVar<CF, Self>;
}

/// The interface for a cryptographic sponge constraints on field `CF`.
/// A sponge can `absorb` or take in inputs and later `squeeze` or output bytes or field elements.
/// The outputs are dependent on previous `absorb` and `squeeze` calls.
pub trait CryptographicSpongeVar<CF: PrimeField, S: CryptographicSponge>: Clone {
    /// Parameters used by the sponge.
    type Parameters;

    /// Initialize a new instance of the sponge.
    fn new(cs: ConstraintSystemRef<CF>, params: &Self::Parameters) -> Self;

    /// Returns a ref to the underlying constraint system the sponge is operating in.
    fn cs(&self) -> ConstraintSystemRef<CF>;

    /// Absorb an input into the sponge.
    fn absorb(&mut self, input: &impl AbsorbGadget<CF>) -> Result<(), SynthesisError>;

    /// Squeeze `num_bytes` bytes from the sponge.
    fn squeeze_bytes(&mut self, num_bytes: usize) -> Result<Vec<UInt8<CF>>, SynthesisError>;

    /// Squeeze `num_bit` bits from the sponge.
    fn squeeze_bits(&mut self, num_bits: usize) -> Result<Vec<Boolean<CF>>, SynthesisError>;

    /// Creates a new sponge with applied domain separation.
    fn fork(&self, domain: &[u8]) -> Result<Self, SynthesisError> {
        let mut new_sponge = self.clone();

        let mut input = Absorb::to_sponge_bytes_as_vec(&domain.len());
        input.extend_from_slice(domain);

        let elems: Vec<CF> = input.to_sponge_field_elements_as_vec();
        let elem_vars = elems
            .into_iter()
            .map(|elem| FpVar::Constant(elem))
            .collect::<Vec<_>>();

        new_sponge.absorb(&elem_vars)?;

        Ok(new_sponge)
    }

    /// Squeeze `num_elements` field elements from the sponge.
    fn squeeze_field_elements(
        &mut self,
        num_elements: usize,
    ) -> Result<Vec<FpVar<CF>>, SynthesisError>;
}
