/// The exponent in `Sbox(x) = x^\alpha`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Alpha {
    /// A positive exponent $x^{alpha}$.
    Exponent(u32),
    /// 1/x
    Inverse,
}

impl Alpha {
    /// Return the memory representation of alpha as a byte array in little-endian byte order.
    pub fn to_bytes_le(&self) -> [u8; 4] {
        match self {
            Alpha::Exponent(exp) => exp.to_le_bytes(),
            Alpha::Inverse => (-1i32).to_le_bytes(),
        }
    }
}
