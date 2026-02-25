# Summary method for `gsaot_indices` objects

Summary method for `gsaot_indices` objects

## Usage

``` r
# S3 method for class 'gsaot_indices'
summary(object, digits = 3, ranking = NULL, ...)
```

## Arguments

- object:

  An object of class `"gsaot_indices"`.

- digits:

  (default: `3`) Number of significant digits to print for numeric
  values.

- ranking:

  An integer with absolute value less or equal than the number of
  inputs. If positive, select the first `ranking` inputs per importance.

- ...:

  Further arguments (currently ignored).

## Value

(Invisibly) a named `list` containing the main elements summarised on
screen.
