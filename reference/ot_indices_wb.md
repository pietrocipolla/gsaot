# Estimate Wasserstein-Bures approximation of the optimal transport solution

Estimate Wasserstein-Bures approximation of the optimal transport
solution

## Usage

``` r
ot_indices_wb(
  x,
  y,
  M,
  boot = FALSE,
  R = NULL,
  parallel = "no",
  ncpus = 1,
  conf = 0.95,
  type = "norm"
)
```

## Arguments

- x:

  A matrix or data.frame containing the input(s) values. The values can
  be numeric, factors, or strings. The type of data changes the
  partitioning. If the values are continuous (double), the function
  partitions the data into `M` sets. If the values are discrete
  (integers, strings, factors), the number of partitioning sets is
  data-driven.

- y:

  A matrix containing the output values. Each column represents a
  different output variable, and each row represents a different
  observation. Only numeric values are allowed.

- M:

  A scalar representing the number of partitions for continuous inputs.

- boot:

  (default `FALSE`) Logical that sets whether or not to perform
  bootstrapping of the OT indices.

- R:

  (default `NULL`) Positive integer, number of bootstrap replicas.

- parallel:

  (default `"no"`) The type of parallel operation to be used (if any).
  If missing, the default is taken from the option `boot.parallel` (and
  if that is not set, `"no"`). Only considered if `boot = TRUE`. For
  more information, check the
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) function.

- ncpus:

  (default `1`) Positive integer: number of processes to be used in
  parallel operation: typically one would chose this to the number of
  available CPUs. Check the `ncpus` option in the
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) function of
  the boot package.

- conf:

  (default `0.95`) Number between `0` and `1` representing the default
  confidence level. Only considered if `boot = TRUE`. Different
  confidence levels can be computed as a postprocessing using
  [`confint.gsaot_indices()`](https://pietrocipolla.github.io/gsaot/reference/confint.gsaot_indices.md).

- type:

  (default `"norm"`) Method to compute the default confidence interval.
  Only considered if `boot = TRUE`. For more information, check the
  `type` argument of
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html).
  Different confidence intervals can be computed as a postprocessing
  using
  [`confint.gsaot_indices()`](https://pietrocipolla.github.io/gsaot/reference/confint.gsaot_indices.md).

## Value

A `gsaot_indices` object containing:

- `method`: a string that identifies the type of indices computed.

- `indices`: a names array containing the sensitivity indices between 0
  and 1 for each column in x, indicating the influence of each input
  variable on the output variables.

- `bound`: a double representing the upper bound of the separation
  measure or an array representing the mean of the separation for each
  input according to the bootstrap replicas.

- `x`, `y`: input and output data provided as arguments of the function.

- `inner_statistic`: a list of matrices containing the values of the
  inner statistics for the partitions defined by `partitions`. If
  `method = wasserstein-bures`, each matrix has three rows containing
  the Wasserstein-Bures indices, the Advective, and the Diffusive
  components.

- `partitions`: a matrix containing the partitions built to calculate
  the sensitivity indices. Each column contains the partition associated
  to the same column in `x`.

If `boot = TRUE`, the object contains also:

- `indices_ci`: a `data.frame` with first column the input, second and
  third columns the lower and upper bound of the confidence interval.

- `inner_statistic_ci`: a list of matrices. Each element of the list
  contains the lower and upper confidence bounds for the partition
  defined by the row.

- `bound_ci`: a list containing the lower and upper bounds of the
  confidence intervals of the separation measure bound.

- `type`, `conf`: type of confidence interval and confidence level,
  provided as arguments.

- `W_boot`: list of bootstrap objects, one for each input.

## See also

[`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md),
[`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md)

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

ot_indices_wb(x, y, 10)
#> Method: wass-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4558156 0.4650694 0.1258777 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2862736 0.3033568 0.1100902 
#> 
#> Diffusive component:
#>         X1         X2         X3 
#> 0.16954200 0.16171258 0.01578748 
```
