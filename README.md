# Reference Implementations for Various Versions of HadesHASH
This repository contains the source code of reference implementations for various versions of HadesHASH [1]. Source code is available in both Sage and C++. Moreover, scripts to calculate the round numbers, the round constants, and the MDS matrices are also included.

## *sage_1512_24_x3_bf*
This is an implementation of the Starkad permutation with f(x) = x^3, n = 63, t = 24, N = 1512, and p(x) = x^63 + x + 1.

## *sage_1518_6_inv_pf*
This is an implementation of the Poseidon permutation with f(x) = x^(-1), n = 253, t = 6, N = 1518, and p = 2^252 + 27742317777372353535851937790883648493.

## *sage_1536_24_x3_pf*
This is an implementation of the Poseidon implementation with f(x) = x^3, n = 64, t = 24, N = 1536, and p = 2^64 - 2^8 - 1.

## *cpp_1536_24_x3_pf*
This is an implementation of the Sponge function using the Poseidon permutation with f(x) = x^3, n = 64, t = 24, N = 1536, and p = 2^64 - 2^8 - 1.

## *scripts*
- `calc_round_numbers.py`: Calculate round numbers given N, t, M, the field, and the S-box.
- `create_mds_bf.sage`: Create an MDS matrix for Starkad given N and t.
- `create_mds_pf.sage`: Create an MDS matrix for Poseidon given N, t, and a prime number.
- `create_rcs_grain.sage`: Create round constants given the field, the S-box, n, t, R_F, R_P, and (potentially) a prime number in hexadecimal form.

[1] *Starkad and Poseidon: New Hash Functions for Zero Knowledge Proof Systems*. Cryptology ePrint Archive, Report 2019/458. https://eprint.iacr.org/2019/458.