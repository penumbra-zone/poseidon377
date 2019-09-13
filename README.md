# Reference Implementations for Various Versions of HadesHASH
This repository contains the source code of reference implementations for various versions of HadesHASH [1]. Source code is available in both Sage and C++. Moreover, scripts to calculate the round numbers, the round constants, and the MDS matrices are also included.

## *sage_1512_24_x3_bf*
This is an implementation of the Starkad permutation with f(x) = x^3, n = 63, t = 24, N = 1512, and p(x) = x^63 + x + 1.
- Input:
`0x2e000000000000005800000000000000a8000000000000014000000000000002600000000000000480000000000000088000000000000010000000000000001e0000000000000038000000000000006800000000000000c000000000000001600000000000000280000000000000048000000000000008000000000000000e00000000000000180000000000000028000000000000004000000000000000600000000000000080000000000000008000000000000000`
- Output:
`0xd000b39309d3817f002dcb7bad9f2cd9f11e5032542171508f2cff192df03986010e5e7fd650c328155751eae822bac00146420a51d0aec2f368cb8c259f93a1688926778fdf6220387650ff6ecc0805010aaa4bcca523f35c63a79ea8f5717becb60b16159c7f238573aaad5f2302f18779a144205027cd52ff1cfdd84ea58671894a9301ef01a91d9fcc9370f0e804b14f5874686cf30121bbb7abd0bf00afa4e74dc7bc573bdc64be06f5935233e09cd3e71c851609ba8cc6c194d2`

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