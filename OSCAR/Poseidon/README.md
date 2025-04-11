# Poseidon
This directory contains the [OSCAR](https://www.oscar-system.org/) implementation of the polynomial model for [Poseidon](https://www.usenix.org/conference/usenixsecurity21/presentation/grassi) in permutation truncation mode, and scripts F4 experiments.

For various Poseidon instances, the script
- [Poseidon_Compression_Partial_Guessing_Experiment.jl](./Poseidon_Compression_Partial_Guessing_Experiment.jl) guesses the last $\frac{n}{2}$ S-boxes in the partial rounds, and computes the DRL Gröbner basis with F4.
- [Poseidon_Compression_Second_Last_Guessing_Experiment.jl](./Poseidon_Compression_Second_Last_Guessing_Experiment.jl) guesses $\frac{n}{2}$ S-boxes in the second to last full round, and computes the DRL Gröbner basis with F4.

## Requirements
- [OSCAR](https://www.oscar-system.org/) 1.2.2.

## Usage
For a usage example see [Poseidon_Demo.ipynb](./Poseidon_Demo.ipynb).

The scripts for small scale experiments can be executed with
```Shell
PATH_TO_JULIA\julia Poseidon_Compression_Partial_Guessing_Experiment.jl

PATH_TO_JULIA\julia Poseidon_Compression_Second_Last_Guessing_Experiment.jl
```
The Poseidon parameters (finite field, state sizes, round numbers interval) for the experiments can be modified within the respective sript.
