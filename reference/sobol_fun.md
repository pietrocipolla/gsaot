# Sobol G function evaluation

This function evaluates the Sobol G function on a set of input samples
generated via crude Monte Carlos. It returns both the sampled inputs and
the corresponding function outputs.

## Usage

``` r
sobol_fun(N, a = c(0, 1, 4.5, 9, 99, 99, 99, 99))
```

## Arguments

- N:

  Integer. Number of input samples to generate.

- a:

  (default: `c(0, 1, 4.5, 9, 99, 99, 99, 99)`) Numeric vector of
  non-negative parameters of length 8. These parameters control the
  sensitivity of each input dimension.

## Value

A list with two elements:

- `x`: a numeric matrix of size `N x 8` containing the input samples.

- `y`: a numeric vector of length `N` with the corresponding function
  outputs.

## Details

The Sobol G function is defined as: \$\$ Y = \prod\_{j=1}^{8} \frac{\|\\
4 X_j - 2\\ \| + a_j}{1 + a_j} \$\$ where \\X_j \sim \mathcal{U}(0, 1)\\
independently.

## See also

[`ishi_homma_fun`](https://pietrocipolla.github.io/gsaot/reference/ishi_homma_fun.md),
[`gaussian_fun`](https://pietrocipolla.github.io/gsaot/reference/gaussian_fun.md)

## Examples

``` r
result <- sobol_fun(1000)
head(result$x)
#>             X1        X2         X3        X4        X5          X6         X7
#> [1,] 0.8336769 0.3482587 0.68735827 0.3343985 0.1063165 0.648367345 0.17086912
#> [2,] 0.6894955 0.4466470 0.29258650 0.7959399 0.2118618 0.007476662 0.73724916
#> [3,] 0.5746480 0.8284600 0.59923234 0.7537201 0.9551472 0.225094068 0.30623252
#> [4,] 0.1020558 0.7184446 0.90163744 0.4497658 0.1815488 0.496164552 0.04597023
#> [5,] 0.4036010 0.4013458 0.47526941 0.6295620 0.2172066 0.313823511 0.22045919
#> [6,] 0.9683687 0.9211852 0.02474348 0.7432475 0.2371152 0.690876706 0.90556190
#>             X8
#> [1,] 0.7974711
#> [2,] 0.4254239
#> [3,] 0.6467902
#> [4,] 0.3462589
#> [5,] 0.6898496
#> [6,] 0.9185341
head(result$y)
#> [1] 0.9956661 0.4554667 0.3088842 1.5190787 0.2134684 2.9514805
```
