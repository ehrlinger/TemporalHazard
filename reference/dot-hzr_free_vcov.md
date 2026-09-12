# Free-parameter vcov submatrix for the delta-method sandwich

Fixed parameters (e.g. `fixed = "shapes"`) leave NA rows/cols in the
expanded vcov. Treat them as known-with-zero-variance: restrict the
sandwich to the free submatrix. (The CoE-conserved `log_mu` normally
participates, because CoE fits use the full-information vcov, but it
leaves an NA row, like a fixed parameter, when that recomputation was
unavailable.) Returns `NULL` (with a warning) when CLs cannot be
computed. Shared by the aggregate and decomposed se.fit paths.

## Usage

``` r
.hzr_free_vcov(vcov_mat, p)
```

## Arguments

- vcov_mat:

  The fitted vcov (or NULL / wrong shape).

- p:

  Length of the parameter vector.

## Value

`list(vcov_use, free_idx)`, or `NULL` if unusable.
