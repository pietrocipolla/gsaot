# Compute confidence intervals for sensitivity indices

Computes confidence intervals for a `gsaot_indices` object using
bootstrap results.

## Usage

``` r
# S3 method for class 'gsaot_indices'
confint(object, parm = NULL, level = 0.95, type = "norm", ...)
```

## Arguments

- object:

  An object of class `gsaot_indices`, with bootstrap results included.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  (default is 0.95) Confidence level for the interval.

- type:

  (default is `"norm"`) Method to compute the confidence interval. For
  more information, check the `type` option of
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- ...:

  Additional arguments (currently unused).

## Value

A data frame with the following columns:

- `input`: Name of the input variable.

- `component`: The index component for Wasserstein-Bures.

- `index`: Estimated indices

- `original`: Original estimates.

- `bias`: Bootstrap bias estimate.

- `low.ci`: Lower bound of the confidence interval.

- `high.ci`: Upper bound of the confidence interval.

## Examples

``` r
N <- 1000

mx <- c(1, 1, 1)
Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)

x1 <- rnorm(N)
x2 <- rnorm(N)
x3 <- rnorm(N)

x <- cbind(x1, x2, x3)
x <- mx + x %*% chol(Sigmax)

A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
y <- t(A %*% t(x))

x <- data.frame(x)
y <- y

res <- ot_indices_wb(x, y, 10, boot = TRUE, R = 100)
confint(res, parm = c(1,3), level = 0.9)
#>   input  component      index   original        bias      low.ci    high.ci
#> 1    X1 wass-bures 0.46760919 0.47199925 0.004390061 0.452283910 0.48293447
#> 2    X3 wass-bures 0.12797107 0.13377761 0.005806539 0.111097756 0.14484438
#> 3    X1  advective 0.29460519 0.29686317 0.002257976 0.285055250 0.30415514
#> 4    X3  advective 0.11674128 0.11975069 0.003009411 0.102611338 0.13087122
#> 5    X1  diffusive 0.17300400 0.17513608 0.002132084 0.166122576 0.17988542
#> 6    X3  diffusive 0.01122979 0.01402692 0.002797127 0.008019041 0.01444054
```
