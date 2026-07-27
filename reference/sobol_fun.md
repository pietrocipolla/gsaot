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
#>              X1        X2         X3         X4        X5        X6        X7
#> [1,] 0.04552900 0.7795283 0.06714148 0.01327198 0.6054234 0.5360592 0.3676227
#> [2,] 0.82367214 0.6335010 0.08543624 0.03543032 0.9162302 0.1351674 0.3837757
#> [3,] 0.05253056 0.5010970 0.62526495 0.24791247 0.3098489 0.1531811 0.4615969
#> [4,] 0.36573959 0.8052777 0.59763193 0.44752727 0.5344002 0.2286109 0.2566831
#> [5,] 0.03408520 0.9418859 0.50184952 0.95659238 0.5010350 0.7913539 0.3236823
#> [6,] 0.83065819 0.3581014 0.04889978 0.27704814 0.6782374 0.5854085 0.4236421
#>              X8
#> [1,] 0.41263543
#> [2,] 0.29086120
#> [3,] 0.57768952
#> [4,] 0.31725012
#> [5,] 0.06745867
#> [6,] 0.08828817
head(result$y)
#> [1] 2.3273854 1.2123990 0.8066764 0.4831918 2.2789030 1.1637396
```
