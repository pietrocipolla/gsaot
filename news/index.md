# Changelog

## gsaot 1.1.1

CRAN release: 2025-09-17

- Solved compatibility issue with `ggplot2` 4.0.0

- Solved a bug for `solver="sinkhorn_log"`

## gsaot 1.1.0

CRAN release: 2025-07-18

- Added a [`confint()`](https://rdrr.io/r/stats/confint.html) method for
  class `gsaot_indices`

- Added bootstrap data to `gsaot_indices` objects

## gsaot 1.0.1

- Added a [`summary()`](https://rdrr.io/r/base/summary.html) method for
  class `gsaot_indices`

- In
  [`plot_comparison()`](https://pietrocipolla.github.io/gsaot/reference/plot_comparison.md),
  added checks for inputs

- All titles in title style

## gsaot 1.0.0

CRAN release: 2025-06-17

- In `plot_functions()`, changed the y axis label for consistency with
  theory

- Split `lower_bound()` into
  [`entropic_bound()`](https://pietrocipolla.github.io/gsaot/reference/entropic_bound.md)
  and
  [`irrelevance_threshold()`](https://pietrocipolla.github.io/gsaot/reference/irrelevance_threshold.md)
  for the two use cases previously identified by `"bound"`.

- In
  [`plot_separations()`](https://pietrocipolla.github.io/gsaot/reference/plot_separations.md),
  changed the label on the y-axis to `"Separation Measurements"`

- In
  [`print.gsaot_indices()`](https://pietrocipolla.github.io/gsaot/reference/print.gsaot_indices.md),
  removed the upper bound from display

- In `ot_indices_bw()`, changed names of the indices

- In
  [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md),
  removed checks on `y` is `cost!="L2"`

- Added
  [`plot_comparison()`](https://pietrocipolla.github.io/gsaot/reference/plot_comparison.md)
  to graphically compare different indices

- [`print.gsaot_indices()`](https://pietrocipolla.github.io/gsaot/reference/print.gsaot_indices.md)
  output simplified

- [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md),
  [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md),
  and
  [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
  now work correctly with NAs in input

## gsaot 0.2.0

CRAN release: 2025-04-21

- In
  [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md),
  improved solver `sinkhorn_log`.

- `lower_bound()` has an improved documentation.

- [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md)
  now allows the use of different L_p^p norms.

- [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md),
  [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md),
  and
  [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
  have an extended bootstrap output.

- [`plot_separations()`](https://pietrocipolla.github.io/gsaot/reference/plot_separations.md)
  is the new name of `plot_inner_stats()`.

## gsaot 0.1.0

CRAN release: 2025-01-10

- Added a NEWS.md file to track changes to the package.
- Initial CRAN submission.
