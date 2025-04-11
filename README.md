This repository contains polynomial models of ciphers, hash functions and pseudo-random functions analyzed in the PhD thesis [Algebraic Aspects of Symmetric Cryptography]().

In [SageMath](https://www.sagemath.org/) the following primitives and experiments are implemented:
- [GMiMC](https://doi.org/10.1007/978-3-030-29962-0_8).
    - Polynomial models for contracting/expanding round functions, univariate keys, and sponge functions.
    - Experiments for generic coordinates.
- [MiMC](https://doi.org/10.1007/978-3-662-53887-6_7).
    - Polynomial models for univariate and Feistel-MiMC.
    - Experiments for the empirical solving degree and Castelnuovo-Mumford regularity.
- [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9).
    - Script for generic coordinates verification.
- Toy Feistel ciphers with expanding and contracting round functions.
    - Polynomial models for contracting/expanding round function.
    - Experiments for generic coordinates

In [OSCAR](https://www.oscar-system.org/) the following primitives and experiments are implemented:
- [Ciminion2](https://doi.org/10.46586/tosc.v2025.i1.240-275).
    - Small scale experiment for the computation of the characteristic polynomial.
    - Small scale experiment for term order conversion with FGLM.
- [Hydra](https://doi.org/10.1007/978-3-031-30634-1_9).
    - Small scale experiment for the computation of the characteristic polynomial.
    - Small scale experiment for term order conversion with FGLM.
    - Small scale experiment for the step degrees in a Gröbner basis computation with F4.
- [Poseidon](https://www.usenix.org/conference/usenixsecurity21/presentation/grassi).
    - Polynomial model for permutation truncation mode.
    - Guessing experiments for Gröbner basis computations.


## Requirements
- [SageMath](https://www.sagemath.org/) 10.5.
- [OSCAR](https://www.oscar-system.org/) 1.2.2.
