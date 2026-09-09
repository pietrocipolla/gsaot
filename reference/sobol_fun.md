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
#>              X1        X2         X3        X4        X5         X6        X7
#> [1,] 0.93817631 0.7343972 0.15627652 0.1825949 0.8387762 0.93928958 0.5887905
#> [2,] 0.08487855 0.5164544 0.07002028 0.4265761 0.6151018 0.49714912 0.6191580
#> [3,] 0.14492355 0.4649305 0.48249120 0.4128417 0.5108307 0.64075592 0.9841738
#> [4,] 0.44146936 0.4112679 0.05258835 0.7331824 0.6561861 0.74096539 0.1297382
#> [5,] 0.21849819 0.4464265 0.55658482 0.5272283 0.6089121 0.80222051 0.9864185
#> [6,] 0.31135777 0.2964584 0.09442593 0.5175281 0.6441465 0.09744557 0.5959897
#>             X8
#> [1,] 0.9684624
#> [2,] 0.8194618
#> [3,] 0.5271801
#> [4,] 0.3035369
#> [5,] 0.4066131
#> [6,] 0.5798095
head(result$y)
#> [1] 1.8876291 0.9136084 0.6205143 0.1798975 0.5348975 0.6833754
```
