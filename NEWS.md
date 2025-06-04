* Split `lower_bound()` into `entropic_bound()` and `irrelevance_threshold()` for 
  the two use cases previously identified by `"bound"`.

* In `plot_separations()`, changed the label on the y-axis to `"Separation Measurements"`

* In `print.gsaot_indices()`, removed the upper bound from display

* In `ot_indices_bw()`, changed names of the indices

* In `ot_indices()`, removed checks on `y` is `cost!="L2"`

* Added `plot_comparison()` to graphically compare different indices

* `print.gsaot_indices()` output simplified 

* `ot_indices()`, `ot_indices_1d()`, and `ot_indices_wb()` now work correctly with NAs in input

# gsaot 0.2.0

* In `ot_indices()`, improved solver `sinkhorn_log`.

* `lower_bound()` has an improved documentation.

* `ot_indices_1d()` now allows the use of different L_p^p norms.

* `ot_indices()`, `ot_indices_1d()`, and `ot_indices_wb()` have an extended bootstrap output.

* `plot_separations()` is the new name of `plot_inner_stats()`.

# gsaot 0.1.0

* Added a NEWS.md file to track changes to the package.
* Initial CRAN submission.
