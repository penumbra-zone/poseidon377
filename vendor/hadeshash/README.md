# Reference Implementations for Various Versions of Starkad and Poseidon
This repository contains the source code of reference implementations for various versions of Poseidon [1]. The source code is available in Sage. Moreover, scripts to calculate the round numbers, the round constants, and the MDS matrices are also included.

### Update from 07/03/2021
We fixed several bugs in the implementation. First, the linear layer was computed as `state = state * M` instead of `state = M * state`, and secondly the final matrix multiplication was missing. The test vectors were also changed accordingly.

[1] *Poseidon: A New Hash Function for Zero-Knowledge Proof Systems*. Cryptology ePrint Archive, Report 2019/458. https://eprint.iacr.org/2019/458. Accepted at USENIX'21.
