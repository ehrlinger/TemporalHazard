# Read a SAS `outhaz` estimate dataset

`PROC HAZARD`'s `outhaz=` dataset stores the converged estimates and the
asymptotic variance-covariance matrix at full double precision, where
the printed `.lst` carries about seven significant figures. For any
quantity the dataset holds, it is the better parity reference: print
precision stops being the binding constraint and optimizer convergence
takes over.

## Usage

``` r
hzr_read_outhaz(path)
```

## Arguments

- path:

  Path to a `.sas7bdat` written by `outhaz=`, or to an `.rds` holding
  the same data frame.

## Value

An `hzr_outhaz` object: a list with `estimates` (named numeric),
`status` (named integer, 1 free / 0 fixed), `vcov` (matrix over free
parameters, dimnames set) and `flags` (named numeric of model-structure
flags). When no parameter is free, `vcov` is a 0x0 matrix rather than
`NULL`, so check its dimensions rather than
[`is.null()`](https://rdrr.io/r/base/NULL.html).

## Details

The log-likelihood is *not* stored here; take it from the `.lst`.

## Experimental

This function is experimental and its return shape is expected to
change. The result carries the fitted model, so
[`predict.hzr_outhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hzr_outhaz.md)
predicts from it – that is what the `hzr_translate_sas(librefs = )` path
emits – but it is not a `hazard` object and none of the other `hazard`
methods apply to it. The `_STATUS_` coding is asserted against a
synthetic fixture, so a real `OUTHAZ=` file using a different convention
would yield an empty `vcov` alongside a fully populated `estimates`.

## Examples

``` r
f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
if (nzchar(f)) str(hzr_read_outhaz(f))
#> List of 4
#>  $ estimates: Named num [1:11] 0 0.0314 1.4142 0 0 ...
#>   ..- attr(*, "names")= chr [1:11] "DELTA" "THALF" "NU" "M" ...
#>  $ status   : Named int [1:11] 0 1 1 0 0 0 0 0 1 1 ...
#>   ..- attr(*, "names")= chr [1:11] "DELTA" "THALF" "NU" "M" ...
#>  $ vcov     : num [1:4, 1:4] 0.1013 0.0245 -0.0187 0.0138 0.0245 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "THALF" "NU" "E0" "C0"
#>   .. ..$ : chr [1:4] "THALF" "NU" "E0" "C0"
#>  $ flags    : Named num [1:6] 2 1 0 0 0 0
#>   ..- attr(*, "names")= chr [1:6] "G1FLAG" "FIXDEL0" "FIXMNU1" "G3FLAG" ...
#>  - attr(*, "class")= chr "hzr_outhaz"
```
