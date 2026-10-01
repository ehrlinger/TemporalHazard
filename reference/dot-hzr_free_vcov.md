# Vcov submatrix for the delta-method sandwich

Restricts the sandwich to the parameters the prediction uses. A
parameter held fixed (`fixed`, from the fit's `fixed_mask`; e.g.
`fixed = "shapes"`) with an NA row is dropped as known. (The
CoE-conserved `log_mu` is one when its full-information recompute
failed; the caller warns.) A fixed one with a variance, such as a shape
a g3 constraint derives, stays in. An NA variance on a parameter that is
NOT fixed was masked by `.hzr_safe_solve()` as non-positive: it
withholds the standard error if the prediction uses it (#586). Returns
`NULL` (with a warning) when CLs cannot be computed. Shared by the
aggregate and decomposed se.fit paths.

## Usage

``` r
.hzr_free_vcov(
  vcov_mat,
  p,
  unused = integer(0),
  fixed = NULL,
  param_names = NULL
)
```

## Arguments

- vcov_mat:

  The fitted vcov (or NULL / wrong shape).

- p:

  Length of the parameter vector.

- unused:

  Indices of parameters the prediction does not depend on, by the
  prediction's form (`.hzr_weibull_used()`, `.hzr_multiphase_used()`),
  never read off numeric zeros in the Jacobian. They are dropped
  exactly, whatever their variance.

- fixed:

  Logical, `TRUE` for a parameter held fixed in the fit, or `NULL` when
  nothing is fixed (every single-distribution fit).

- param_names:

  Names for the warning, or `NULL`.

## Value

`list(vcov_use, free_idx)`, or `NULL` if unusable.
