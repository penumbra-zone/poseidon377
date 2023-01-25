use ark_ff::PrimeField;
use ark_std::vec::Vec;
use num::integer::gcd;
use num_bigint::BigUint;

/// Shortest addition chains for small numbers.
///
/// The numbers at each depth are sorted in descending order.
///
/// Used for selecting the optimal choice of alpha.
///
/// Courtesy of:
/// https://wwwhomes.uni-bielefeld.de/achim/addition_chain.html
pub struct ShortestAdditionChains {
    depth_2: [u32; 2],
    depth_3: [u32; 3],
    depth_4: [u32; 5],
    depth_5: [u32; 9],
    depth_6: [u32; 15],
    depth_7: [u32; 26],
}

const SHORTEST_ADDITION_CHAINS: ShortestAdditionChains = ShortestAdditionChains {
    depth_2: [4, 3],
    depth_3: [8, 6, 5],
    depth_4: [16, 12, 10, 9, 7],
    depth_5: [32, 24, 20, 18, 17, 15, 14, 13, 11],
    depth_6: [64, 48, 40, 36, 34, 33, 30, 28, 27, 26, 25, 23, 22, 21, 19],
    depth_7: [
        128, 96, 80, 72, 68, 66, 65, 60, 56, 54, 52, 51, 50, 49, 46, 45, 44, 43, 42, 41, 39, 38,
        37, 35, 31, 29,
    ],
};

impl ShortestAdditionChains {
    fn depths_in_order(&self) -> Vec<&u32> {
        let mut depths: Vec<&u32> = self.depth_2.iter().chain(self.depth_3.iter()).collect();
        let depths_4_5: Vec<&u32> = self.depth_4.iter().chain(self.depth_5.iter()).collect();
        let depths_6_7: Vec<&u32> = self.depth_6.iter().chain(self.depth_7.iter()).collect();
        depths.extend(depths_4_5);
        depths.extend(depths_6_7);
        depths
    }
}

/// The exponent in `Sbox(x) = x^\alpha`.
pub struct Alpha<T>(pub T);

impl Alpha<poseidon_parameters::Alpha> {
    /// Select the best choice of `Alpha` given the parameters.
    pub fn generate<F: PrimeField>(
        p: F::BigInt,
        allow_inverse: bool,
    ) -> poseidon_parameters::Alpha {
        // Move through the addition chains in increasing depth,
        // picking the leftmost choice that meets the coprime requirement.
        for candidate in SHORTEST_ADDITION_CHAINS.depths_in_order() {
            if candidate % 2 != 0 && Self::alpha_coprime_to_p_minus_one::<F>(*candidate, p) {
                return poseidon_parameters::Alpha::Exponent(*candidate);
            }
        }

        if allow_inverse {
            poseidon_parameters::Alpha::Inverse
        } else {
            panic!("could not find a small positive exponent and allow_inverse was not enabled")
        }
    }

    fn alpha_coprime_to_p_minus_one<F: PrimeField>(alpha: u32, p: F::BigInt) -> bool {
        let one: BigUint = F::one().into();
        let p_minus_one: BigUint = p.into() - one;
        let alpha_bigint: BigUint = F::from(alpha).into();
        let computed_gcd = gcd(alpha_bigint, p_minus_one);
        F::from(computed_gcd) == F::one()
    }

    /// Return the memory representation of alpha as a byte array in little-endian byte order.
    pub fn to_bytes_le(&self) -> [u8; 4] {
        match self.0 {
            poseidon_parameters::Alpha::Exponent(exp) => exp.to_le_bytes(),
            poseidon_parameters::Alpha::Inverse => (-1i32).to_le_bytes(),
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::{Fq as Fq254, FqParameters as FqParameters254};
    use ark_ed_on_bls12_377::{Fq as Fq377, FqParameters as FqParameters377};
    use ark_ed_on_bls12_381::{Fq as Fq381, FqParameters as FqParameters381};
    use ark_ff::FpParameters;
    use num_bigint::BigInt;

    use super::*;

    #[test]
    fn test_gcd() {
        let computed_gcd = gcd(BigInt::from(16u32), BigInt::from(24u32));
        assert_eq!(computed_gcd, BigInt::from(8u32));
        let computed_gcd = gcd(BigInt::from(33u32), BigInt::from(11u32));
        assert_eq!(computed_gcd, BigInt::from(11u32));
        let computed_gcd = gcd(BigInt::from(17u32), BigInt::from(11u32));
        assert_eq!(computed_gcd, BigInt::from(1u32));
        let computed_gcd = gcd(BigInt::from(11586173u32), BigInt::from(63141853u32));
        assert_eq!(computed_gcd, BigInt::from(1u32));
        let p: BigUint = FqParameters377::MODULUS.into();
        let p_minus_one: BigUint = p - BigUint::from(1u32);
        let computed_gcd = gcd(BigInt::from(5u32), p_minus_one.into());
        assert_eq!(computed_gcd, BigInt::from(5u32))
    }

    #[test]
    fn check_alpha_5() {
        // We know from the Poseidon paper that we should get an alpha of 5 for
        // BLS12-381 and BN254 (see Table 2)
        let p = FqParameters381::MODULUS;
        assert_eq!(
            Alpha::generate::<Fq381>(p, true),
            poseidon_parameters::Alpha::Exponent(5)
        );

        let p = FqParameters254::MODULUS;
        assert_eq!(
            Alpha::generate::<Fq254>(p, true),
            poseidon_parameters::Alpha::Exponent(5)
        );
    }

    #[test]
    fn check_alpha_17() {
        // For Poseidon377, we should get an alpha of 17 (from our own work).
        let p = FqParameters377::MODULUS;
        assert_eq!(
            Alpha::generate::<Fq377>(p, true),
            poseidon_parameters::Alpha::Exponent(17)
        );
    }
}
