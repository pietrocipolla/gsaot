# Entropic lower bounds for entropic optimal transport sensitivity indices

Calculate entropic lower bounds for entropic Optimal Transport
sensitivity indices

## Usage

``` r
entropic_bound(
  y,
  M,
  cost = "L2",
  discrete_out = FALSE,
  solver = "sinkhorn",
  solver_optns = NULL,
  scaling = TRUE
)
```

## Arguments

- y:

  An array or a matrix containing the output values.

- M:

  A scalar representing the number of partitions for continuous inputs.

- cost:

  (default `"L2"`) A string or function defining the cost function of
  the Optimal Transport problem. It should be "L2" or a function taking
  as input y and returning a cost matrix. If `cost="L2"`, `ot_indices`
  uses the squared Euclidean metric.

- discrete_out:

  (default `FALSE`) Logical, by default the output sample in `y` are
  equally weighted. If `discrete_out=TRUE`, the function tries to create
  an histogram of the realizations and to use the histogram as weights.
  It works if the output is discrete or mixed and the number of
  realizations is large. The advantage of this option is to reduce the
  dimension of the cost matrix.

- solver:

  Solver for the Optimal Transport problem. Currently supported options
  are:

  - `"sinkhorn"` (default), the Sinkhorn's solver (Cuturi 2013) .

  - `"sinkhorn_log"`, the Sinkhorn's solver in log scale (Peyré et
    al. 2019) .

- solver_optns:

  (optional) A list containing the options for the Optimal Transport
  solver. See details for allowed options and default ones.

- scaling:

  (default `TRUE`) Logical that sets whether or not to scale the cost
  matrix.

## Value

A scalar representing the entropic lower bound.

## Details

The function allows the computation of the entropic lower bounds.
`solver` should be either `"sinkhorn"` or `"sinkhorn_log"`.

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

M <- 25

sink_lb <- entropic_bound(y, M)
```
