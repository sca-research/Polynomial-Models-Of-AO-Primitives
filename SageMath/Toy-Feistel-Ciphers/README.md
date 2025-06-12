# Toy Feistel Ciphers
This directory contains [SageMath](https://www.sagemath.org/) implementations of the polynomial models for toy Feistel ciphers with:
- Expanding round function
```math
    \begin{pmatrix}{c}
        x_1 \\ 
        \vdots \\ 
        x_n
    \end{pmatrix}
    \mapsto
    \begin{pmatrix}
        x_n \\
        x_1 + x_n^d \\
        \vdots \\
        x_{n - 1} + x_n^d 
    \end{pmatrix}
    .
```
- Contracting round function
```math
    \begin{pmatrix}
        x_1 \\ 
        \vdots \\ 
        x_n
    \end{pmatrix}
    \mapsto
    \begin{pmatrix}
        x_n + \left( \sum_{i = 1}^{n - 1} x_i \right)^d \\
        x_1 \\
        \vdots \\
        x_{n - 1}
    \end{pmatrix}
    .
```

Both variants do not deploy a key schedule, i.e., the master key is added after every round.

In addition, it contains the Notebook [Feistel_Generic_Coordinates.ipynb](./Feistel_Generic_Coordinates.ipynb) which verifies generic coordinates for various instances of contracting/expanding Feistel ciphers.

## Requirements
- [SageMath](https://www.sagemath.org/) 10.6.

## Usage
For usage examples of the various instances see [Feistel_usage_example.sage](./Feistel_usage_example.sage) and [Feistel_Demo.ipynb](./Feistel_Demo.ipynb).
