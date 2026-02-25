# Irrelevance threshold for optimal transport sensitivity indices

Calculate irrelevance threshold using dummy variable for Optimal
Transport sensitivity indices

## Usage

``` r
irrelevance_threshold(
  y,
  M,
  dummy_optns = NULL,
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

- dummy_optns:

  (default `NULL`) A list containing the options on the distribution of
  the dummy variable. See `details` for more information.

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

  - `"1d"`, the one-dimensional analytic solution.

  - `"wasserstein-bures"`, the Wasserstein-Bures solution.

  - `"sinkhorn"` (default), the Sinkhorn's solver (Cuturi 2013) .

  - `"sinkhorn_log"`, the Sinkhorn's solver in log scale (Peyré et
    al. 2019) .

  - `"transport"`, a solver of the non regularized OT problem using
    [`transport::transport()`](https://rdrr.io/pkg/transport/man/transport.html).

- solver_optns:

  (optional) A list containing the options for the Optimal Transport
  solver. See details for allowed options and default ones.

- scaling:

  (default `TRUE`) Logical that sets whether or not to scale the cost
  matrix.

## Value

An object of class `gsaot_indices`.

## Details

The function allows the computation of irrelevance threshold. The
function samples from a distribution defined in `dummy_optns` (by
default a standard normal), independent from the output `y` and then
computes the indices using the algorithm specified in `solver`. Under
the hood, `lower_bound` calls the other available functions in the
package:

- [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md)
  (for `solver="1d"`)

- [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
  (for `solver="wasserstein-bures"`)

- [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
  (for `solver %in% c("sinkhorn", "sinkhorn_log", "wasserstein")`) The
  user can choose the distribution of the dummy variable using the
  argument `dummy_optns`. `dummy_optns` should be a named list with at
  least a term called `"distr"` defining the sampling function. The
  other terms in the list are used as arguments to the sampling
  function.

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

dummy_lb <- irrelevance_threshold(y, M)

# Custom sampling funtion and network simplex solver
dummy_optns <- list(distr = "rgamma", shape = 3)
dummy_lb_cust <- irrelevance_threshold(y, M,
                                      dummy_optns = dummy_optns,
                                      solver = "transport")
```
