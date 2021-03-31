# Reference Implementations for Various Versions of Starkad and Poseidon
This repository contains the source code of reference implementations for various versions of Poseidon [1]. The source code is available in Sage. Moreover, scripts to calculate the round numbers, the round constants, and the MDS matrices are also included.

### Update from 31/03/2021
The round numbers were adjusted in order to match the results from `code/calc_round_numbers.py`. The test vectors were also changed accordingly.

### Update from 07/03/2021
We fixed several bugs in the implementation. First, the linear layer was computed differently, and secondly the final matrix multiplication was missing. The test vectors were also changed accordingly.

In more detail, the linear layer was computed as `state = state * M` instead of `state = M * state` in the previous version, hence essentially the transpose `M^T` of `M` was used. Security-wise, this was no problem, since `M^T` is also secure w.r.t. subspace trails in the partial rounds. However, the current version is now exactly following the specification given in the paper.

Regarding the second change, no linear layer operation was done in the last round. This was against our spec and may indeed lead to less security in the specific hash setting we are considering. It is now fixed (i.e., it now follows the specification given in the paper), and the last round includes a linear layer operation in the current version.


<br>
<br>

[1] *Poseidon: A New Hash Function for Zero-Knowledge Proof Systems*. Cryptology ePrint Archive, Report 2019/458. https://eprint.iacr.org/2019/458. Accepted at USENIX'21.
