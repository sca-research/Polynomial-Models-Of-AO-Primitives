# MiMC
This directory contains [SageMath](https://www.sagemath.org/) implementations of the polynomial models for [MiMC](https://doi.org/10.1007/978-3-662-53887-6_7) with:
- A single plain/ciphertext sample.
- Two plain/ciphertext sample.
- A single plain/ciphertext for MiMC.

The file [lazard_gb_algorithm.py](./lazard_gb_algorithm.py) contains a plain SageMath implementation of Gaussian elimination on the Macaulay matrix.
The code is based on a blog post from Jan Ferdinand Sauer.
The original post cannot be accessed anymore, but it is still available through the [internet archive](https://web.archive.org/web/20241016004144/https://asdm.gmbh/2021/03/15/d_reg/).

In addition, it contains Notebooks which compute the empirical solving degree and Castelnuovo-Mumford regularity for increasing round numbers:
- [MiMC_field_equation_solving_degree.ipynb](./MiMC_field_equation_solving_degree.ipynb).
- [MiMC_field_equation_remainder_solving_degree.ipynb](./MiMC_field_equation_remainder_solving_degree.ipynb).
- [MiMC_all_field_equations_solving_degree.ipynb](./MiMC_all_field_equations_solving_degree.ipynb).
- [MiMC_two_plain_text_solving_degree.ipynb](./MiMC_two_plain_text_solving_degree.ipynb)


## Requirements
- [SageMath](https://www.sagemath.org/) 10.5.

## Usage
For usage examples of the various instances see [MiMC_Demo.ipynb](./MiMC_Demo.ipynb).
