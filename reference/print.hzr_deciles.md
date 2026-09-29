# Print method for hzr_deciles

Print method for hzr_deciles

## Usage

``` r
# S3 method for class 'hzr_deciles'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An `hzr_deciles` object.

- digits:

  Number of decimal places for formatting.

- ...:

  Additional arguments (ignored).

## Value

The object `x` of class `c("hzr_deciles", "data.frame")`, invisibly. The
data frame has one row per risk group and columns: `group` (integer
group index, 1 = lowest risk), `n` (group size), `events` (observed
event count; weighted for a weighted fit), `expected` (expected event
count from model predictions, net of any entry time, and weighted as
`events` is), `observed_rate`, `expected_rate` (events / n; per unit of
case weight for a weighted fit, and `NA` for a group whose weights sum
to 0), `chi_sq` (per-group (O-E)^2/V contribution, where V is the
Poisson variance of the weighted count, equal to E for an unweighted
fit), `p_value` (1-df chi-square upper-tail p), `mean_survival`,
`mean_cumhaz` (mean predicted values in group). An `"overall"` attribute
contains the omnibus chi-square test (fields: `chi_sq`, `df`, `p_value`,
`time`, `groups`, `total_events`, `total_expected`, `n_included`,
`n_excluded`).
