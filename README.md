# Poseidon

This repository contains:

* [`poseidon377`](../main/poseidon377): an instantiation of the Poseidon hash function for [`decaf377`](https://github.com/penumbra-zone/decaf377)
* [`poseidon-paramgen`](../main/poseidon-paramgen): an independent implementation of Poseidon parameter generation
* [`poseidon-parameters`](../main/poseidon-parameters): types that represent Poseidon
parameters
* [`poseidon-permutation`](../main/poseidon-permutation): an independent implementation of the Poseidon permutation
* [`poseidon-consistency`](../main/poseidon-consistency): property-based tests for consistency between Poseidon implementations
* [`poseidon-tests`](../main/poseidon-tests): test vectors for `poseidon-parameters` and `poseidon377`

## Audits

`poseidon-paramgen` [was audited](https://research.nccgroup.com/2022/09/12/public-report-penumbra-labs-decaf377-implementation-and-poseidon-parameter-selection-review/) by NCC Group in Summer 2022.
Note that the audit covered only the parameter generation described in the original Poseidon paper. The Poseidon2 parameter
generation is not yet audited.

## Benchmarks

Run `criterion` benchmarks using:

```
cargo bench
```

This will generate a report at `target/criterion/report/index.html`.

Performance benchmarked on commit 9750f5ff01d11f158f111a1a75401901049e5575 on
a 2023 Macbook Pro M2 (12 core CPU) with 32 GB memory using our 4-to-1 optimized
poseidon hash takes 47.3µs, or ~21,141 hashes/second.
