# Multivariate Gaussian linear model evaluation

Generates samples from a multivariate Gaussian distribution and
evaluates a simple linear transformation model.

## Usage

``` r
gaussian_fun(N)
```

## Arguments

- N:

  Number of input samples to generate.

## Value

A list with two elements:

- `x`: a numeric matrix of size `N x 8` containing the input samples.

- `y`: a numeric vector of length `N` with the corresponding function
  outputs.

## Details

Inputs `x` are sampled from: \$\$ \mathbf{X} \sim
\mathcal{N}(\boldsymbol{\mu}, \Sigma), \quad \boldsymbol{\mu} = \[1, 1,
1\], \quad \Sigma = \begin{bmatrix} 1 & 0.5 & 0.5 \\ 0.5 & 1 & 0.5 \\
0.5 & 0.5 & 1 \end{bmatrix} \$\$

The output is given by: \$\$ \mathbf{Y} = A \mathbf{X}^{\top}, \quad A =
\begin{bmatrix} 4 & -2 & 1 \\ 2 & 5 & -1 \end{bmatrix} \$\$

## See also

[`sobol_fun`](https://pietrocipolla.github.io/gsaot/reference/sobol_fun.md),
[`ishi_homma_fun`](https://pietrocipolla.github.io/gsaot/reference/ishi_homma_fun.md)

## Examples

``` r
result <- gaussian_fun(1000)
head(result$x)
#>            X1        X2         X3
#> [1,] 1.210646 0.6066786 0.48065877
#> [2,] 1.717017 0.5524071 1.66414399
#> [3,] 2.504612 2.3095378 1.74635838
#> [4,] 1.466895 0.6142316 0.01331935
#> [5,] 1.524946 2.1775382 1.30649161
#> [6,] 1.268600 2.8835283 1.93367697
head(result$y)
#>            Y1        Y2
#> [1,] 4.109884  4.974026
#> [2,] 7.427396  4.531925
#> [3,] 7.145730 14.810554
#> [4,] 4.652436  5.991629
#> [5,] 3.051201 12.631092
#> [6,] 1.241019 15.021164
```
