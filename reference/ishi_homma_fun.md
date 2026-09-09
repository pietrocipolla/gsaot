# Ishigami-Homma function evaluation

Evaluates the Ishigami-Homma function. Input samples are drawn from a
uniform distribution over \\\[-\pi, \pi\]^3\\

## Usage

``` r
ishi_homma_fun(N, A = 2, B = 1)
```

## Arguments

- N:

  Number of input samples to generate.

- A:

  (default: `2`) Numeric, amplitude of the second sine component .

- B:

  (default: `1`) Numeric, coefficient of the interaction term.

## Value

A list with two elements:

- `x`: a numeric matrix of size `N x 8` containing the input samples.

- `y`: a numeric vector of length `N` with the corresponding function
  outputs.

## Details

The Ishigami-Homma function is defined as: \$\$Y = \sin(X_1) + A \cdot
\sin^2(X_2) + B \cdot X_3^4 \cdot \sin(X_1)\$\$ where \\X_i \sim
\mathcal{U}(-\pi, \pi)\\.

## See also

[`sobol_fun`](https://pietrocipolla.github.io/gsaot/reference/sobol_fun.md),
[`gaussian_fun`](https://pietrocipolla.github.io/gsaot/reference/gaussian_fun.md)

## Examples

``` r
result <- ishi_homma_fun(1000)
head(result$x)
#>             X1         X2       X3
#> [1,]  1.101371 -0.3673841 2.398859
#> [2,] -2.367151  0.4339278 2.468320
#> [3,]  3.047575 -2.4072754 2.112918
#> [4,]  2.807415 -0.2482628 1.169294
#> [5,] -1.078276 -1.9243950 1.279449
#> [6,] -2.124207  0.2548824 1.800859
head(result$y)
#> [1]  30.682321 -26.304282   2.863012   1.061890  -1.482189  -9.671349
```
