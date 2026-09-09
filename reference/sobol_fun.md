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
#>                X1        X2        X3         X4        X5          X6
#> [1,] 0.5824343285 0.2604401 0.2035265 0.83027657 0.6944382 0.334624904
#> [2,] 0.0009700714 0.6029773 0.7669474 0.87926387 0.2442127 0.661772836
#> [3,] 0.9628692346 0.8949505 0.7793144 0.09655633 0.3632611 0.055355085
#> [4,] 0.5684291318 0.3948880 0.4655722 0.25423834 0.9001638 0.006870335
#> [5,] 0.6584779941 0.4305764 0.5843230 0.78229013 0.6681610 0.807332404
#> [6,] 0.8872893630 0.1149982 0.5116943 0.92812649 0.9221846 0.939284217
#>             X7         X8
#> [1,] 0.4541136 0.03535645
#> [2,] 0.5816327 0.76951503
#> [3,] 0.4659482 0.23788377
#> [4,] 0.3117278 0.07702918
#> [5,] 0.9331192 0.59350335
#> [6,] 0.2589919 0.37212269
head(result$y)
#> [1] 0.3426730 1.4864385 2.5759957 0.1669649 0.3607881 1.7583542
```
