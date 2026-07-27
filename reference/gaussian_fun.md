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
#>             X1          X2         X3
#> [1,] 0.5245388 -0.14806128  0.6375819
#> [2,] 1.8116193 -0.08498721 -0.1162312
#> [3,] 1.7390710  3.41798653  1.7744252
#> [4,] 0.6432006  1.06648598 -0.2259466
#> [5,] 1.4622031  0.20912384 -1.5212449
#> [6,] 1.9270348  1.20723027  2.6087168
head(result$y)
#>             Y1         Y2
#> [1,] 3.0318598 -0.3288106
#> [2,] 7.3002203  3.3145337
#> [3,] 1.8947362 18.7936495
#> [4,] 0.2138838  6.8447777
#> [5,] 3.9093198  5.4912704
#> [6,] 7.9023955  7.2815041
```
