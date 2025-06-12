# Ciminion2
This directory contains the [OSCAR](https://www.oscar-system.org/) implementation of the polynomial model for the [Ciminion2](https://doi.org/10.46586/tosc.v2025.i1.240-275) Pseudo-Random Function, and scripts for small scale characteristic polynomial and FGLM experiments.

The files
 - [Ciminion_2.jl](./Ciminion_2.jl),
 - [Ciminion_2_polynomial_model.jl](./Ciminion_2_polynomial_model.jl),
 
 are forked from [Groebner-Basis-Cryptanalysis-of-Ciminion-and-Hydra](https://github.com/sca-research/Groebner-Basis-Cryptanalysis-of-Ciminion-and-Hydra.git) (commit `18e719c04a646ade8a699ef91423d467a79be9ee`) and slightly modified to adhere to the present directory structure.

In addition, it contains the scripts:
- [Ciminion_2_characteristic_polynomial_experiments.jl](./Ciminion_2_characteristic_polynomial_experiments.jl) which constructs the Ciminion2 characteristic polynomial with respect to the last variable for small scale instances.
- [Ciminion_2_fglm_experiments.jl](./Ciminion_2_fglm_experiments.jl) which performs term order conversion to LEX for the Ciminion2 DRL Gröbner basis for small scale instances.

For each round number, $10$ trials are performed and average running times are reported.

## Requirements
- [OSCAR](https://www.oscar-system.org/) 1.4.1.

## Usage
For a usage example see [Ciminion_2_Demo.ipynb](./Ciminion_2_Demo.ipynb).

The scripts for small scale experiments can be executed with
```Shell
PATH_TO_JULIA\julia Ciminion_2_characteristic_polynomial_experiments.jl

PATH_TO_JULIA\julia Ciminion_2_fglm_experiments.jl
```
The Ciminion2 parameters (finite field, round numbers interval) for the experiments can be modified within the respective sript.
