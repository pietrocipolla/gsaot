# Estimate sensitivity maps using optimal transport indices

Estimate sensitivity maps using optimal transport indices

## Usage

``` r
ot_indices_smap(x, y, M)
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

  A matrix containing the output values. Each column is interpreted as a
  different output.

- M:

  A scalar representing the number of partitions for continuous inputs.

## Value

A matrix where each column represents an input and each row represents
an output. The values are indices between 0 and 1 computed using
[`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md).

## Examples

``` r
N <- 1000

x1 <- rnorm(N)
x2 <- rnorm(N)
x <- cbind(x1, x2)

y1 <- 10 * x1
y2 <- x1 + x2
y <- cbind(y1, y2)

ot_indices_smap(data.frame(x), y, 30)
#>           x1         x2
#> y1 0.9395661 0.04029454
#> y2 0.3241271 0.35022339
```
