# Hydra
This directory contains the [OSCAR](https://www.oscar-system.org/) implementation of the polynomial model for the [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9) Pseudo-Random Function, and scripts for small scale characteristic polynomial and F4 experiments.

The files
 - [Hydra.jl](./Hydra.jl),
 - [Hydra_polynomial_model.jl](./Hydra_polynomial_model.jl),
 - [utilities.jl](./utilities.jl)
 
 are forked from [Groebner-Basis-Cryptanalysis-of-Ciminion-and-Hydra](https://github.com/sca-research/Groebner-Basis-Cryptanalysis-of-Ciminion-and-Hydra.git) (commit `18e719c04a646ade8a699ef91423d467a79be9ee`) and slightly modified to adhere to the present directory structure.

In addition, it contains the scripts:
- [Hydra_characteristic_polynomial_experiments.jl](./Hydra_characteristic_polynomial_experiments.jl) which constructs the Hydra characteristic polynomial with respect to the last variable for small scale instances.
- [Hydra_step_degrees_experiments.jl](./Hydra_step_degrees_experiments.jl) which performs term order conversion to LEX for the Hydra DRL Gröbner basis for small scale instances.

For each round number, $10$ trials are performed and average running times are reported.

## Requirements
- [OSCAR](https://www.oscar-system.org/) 1.2.2.

## Usage
For a usage example see [Hydra_Demo.ipynb](./Hydra_Demo.ipynb).

The scripts for small scale experiments can be executed with
```Shell
PATH_TO_JULIA\julia Hydra_characteristic_polynomial_experiments.jl

PATH_TO_JULIA\julia Hydra_fglm_experiments.jl
```
The Hydra parameters (finite field, round numbers interval) for the experiments can be modified within the respective sript.
