# Extract variance-covariance matrix from hazard model

Returns the estimated variance-covariance matrix of the fitted
coefficients.

## Usage

``` r
# S3 method for class 'hazard'
vcov(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Value

A numeric matrix containing the estimated variance-covariance matrix of
the fitted coefficients, with rows and columns named by the coefficient
labels (phase-prefixed for multiphase models, e.g. `early.x`). Rows and
columns for parameters held fixed (e.g. fixed shape parameters) are `NA`
because they carry no variance; the finite free-parameter block is still
usable. For Conservation-of-Events fits the conserved phase `log_mu`
*normally* carries a variance: it is removed from the optimizer search
but the vcov is the full-information matrix at the optimum (the CoE
solution is the unconstrained MLE). That recomputation requires numDeriv
and an invertible Hessian; if either is unavailable the fit emits a
warning and the conserved `log_mu` stays `NA` (the rest of the matrix is
unaffected). Returns a scalar `NA` only when the model has not been
fitted or no covariance matrix is available.

## Examples

``` r
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
vcov(fit)
#>              mu           nu
#> mu  0.011755419 -0.004708797
#> nu -0.004708797  0.018917131
```
