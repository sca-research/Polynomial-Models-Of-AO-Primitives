# Hydra
This directory contains the [SageMath](https://www.sagemath.org/) implementation of the polynomial model for the [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9) Pseudo-Random Function, and a script for generic coordinates verification.

The files
 - [Hydra.sage](./Hydra.sage),
 - [Hydra_polynomial_model.sage](./Hydra_polynomial_model.sage),
 - [utilities.sage](./utilities.sage)
 
 are forked from [Groebner-Basis-Cryptanalysis-of-Ciminion-and-Hydra](https://github.com/sca-research/Groebner-Basis-Cryptanalysis-of-Ciminion-and-Hydra.git) (commit `18e719c04a646ade8a699ef91423d467a79be9ee`) and slightly modified to adhere to the present directory structure.

In addition, it contains the script [Hydra_generic_coordinates.sage](./Hydra_generic_coordinates.sage) which verifies generic coordinates for a Hydra instance in a fixed round numbers interval $r_{min} \leq r_\mathcal{H} \leq r_{max}$.

## Requirements
- [SageMath](https://www.sagemath.org/) 10.5.

## Usage
For a usage example see [Hydra_Demo.sage](./Hydra_Demo.ipynb).

The generic coordinates verification scripts can be executed with
```Shell
PATH_TO_SAGE\sage Hydra_generic_coordinates.sage
```
The Hydra parameters (finite field, matrices, round numbers interval) for the verification can be modified within the respective sript.
