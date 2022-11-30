# Poseidon

This repository contains:

* [`poseidon377`](../main/poseidon377): an instantiation of the Poseidon hash function for [`decaf377`](https://github.com/penumbra-zone/decaf377)
* [`poseidon-paramgen`](../main/poseidon-paramgen): an independent implementation of Poseidon parameter generation
* [`poseidon-permutation`](../main/poseidon-permutation): an independent implementation of the Poseidon permutation

## Audits

`poseidon-paramgen` [was audited](https://research.nccgroup.com/2022/09/12/public-report-penumbra-labs-decaf377-implementation-and-poseidon-parameter-selection-review/) by NCC Group in Summer 2022.

## Benchmarks

Run `criterion` benchmarks using:

```
cargo bench
```

This will generate a report at `target/criterion/report/index.html`.
