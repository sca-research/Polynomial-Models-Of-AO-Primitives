# GMiMC
This directory contains [SageMath](https://www.sagemath.org/) implementations of the polynomial models for [GMiMC](https://doi.org/10.1007/978-3-030-29962-0_8) with:
- contracting/expanding round function and [multivariate keys](./GMiMC.sage).
- contracting/expanding round function and [univariate keys](./GMiMC_univariate.sage).
- [GMiMCHash](GMiMC_erf_sponge.sage), i.e., GMiMC with expanding round function in sponge mode.

In addition, it contains scripts
- [GMiMC_generic_coordinates.sage](./GMiMC_generic_coordinates.sage),
- [GMiMC_univariate_generic_coordinates.sage](./GMiMC_univariate_generic_coordinates.sage),
- [GMiMC_erf_sponge_generic_coordinates.sage](./GMiMC_erf_sponge_generic_coordinates.sage),

which verify generic coordinates for a GMiMC instance in a fixed round numbers interval $r_{min} \leq r \leq r_{max}$.

## Requirements
- [SageMath](https://www.sagemath.org/) 10.6.

## Usage
For usage examples of the various instances see the usage example sage scripts or the Demo notebooks.

The generic coordinates verification scripts can be executed with
```Shell
PATH_TO_SAGE\sage GMiMC_generic_coordinates.sage
```
The GMiMC parameters (finite field, state size, key schedule matrix, round numbers interval) for the verification can be modified within the respective sript.
