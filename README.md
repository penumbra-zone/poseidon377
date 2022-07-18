# Poseidon

This repository contains:

* [`poseidon377`](../main/poseidon377): an instantiation of the Poseidon hash function for [`decaf377`](https://github.com/penumbra-zone/decaf377)
* [`poseidon-paramgen`](../main/poseidon-paramgen): an independent implementation of Poseidon parameter generation
* [`poseidon-permutation`](../main/poseidon-permutation): an independent implementation of the Poseidon permutation

*Warning:* These are in active development and have not yet been security audited.

## Benchmarks

Run `criterion` benchmarks using:

```
cargo bench
```

This will generate a report at `target/criterion/report/index.html`.
