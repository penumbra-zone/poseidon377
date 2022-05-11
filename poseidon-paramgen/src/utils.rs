use std::convert::TryInto;

use ark_ff::BigInteger;
use num_bigint::BigUint;

/// Take the binary log of a `BigInteger`
pub fn log2<P>(x: P) -> f64
where
    P: BigInteger,
{
    let mut p_biguint: BigUint = x.into();
    let two_to_50: BigUint = 1125899906842624u64.into(); // 2**50
    let mut log_bit_boundaries = 0;
    while p_biguint > two_to_50 {
        p_biguint >>= 48u64;
        log_bit_boundaries += 48;
    }
    let mut p_le_bytes = p_biguint.to_bytes_le();
    // `u64::from_le_bytes` takes `[u8; 8]` but `p_le_bytes` is a Vec<u8>
    // so we must pad the vec to the proper size:
    while p_le_bytes.len() < 8 {
        p_le_bytes.push(0x00);
    }

    let x_u64: u64 = u64::from_le_bytes(
        p_le_bytes
            .try_into()
            .expect("we have padded to the expected length"),
    );
    log_bit_boundaries as f64 + ((x_u64) as f64).log2()
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ed_on_bls12_381::FqParameters as Fq381Parameters;
    use ark_ff::{BigInteger256, FpParameters};

    #[test]
    fn log2_bigint() {
        let test_val: BigInteger256 = 4.into();
        assert_eq!(log2(test_val), 2.0);

        let test_val: BigInteger256 = 257.into();
        assert!(log2(test_val) > 8.005);
        assert!(log2(test_val) < 8.006);

        let test_val: BigInteger256 = 65536.into();
        assert_eq!(log2(test_val), 16.0);

        let test_val = Fq381Parameters::MODULUS;
        assert!(log2(test_val) > 254.856);
        assert!(log2(test_val) < 254.858);
    }
}
