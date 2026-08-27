# Legacy HAZARD to TemporalHazard argument mapping

Returns a formal mapping table that defines how legacy SAS
HAZARD/C-style inputs map to `hazard(...)` arguments in this package.

## Usage

``` r
hzr_argument_mapping(include_planned = TRUE)
```

## Arguments

- include_planned:

  Logical; if `FALSE`, only rows marked as implemented are returned.

## Value

A data frame with one row per mapping rule.

## See also

[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
for the target arguments,
[`hzr_phase()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_phase.md)
for phase construction, and
[`hzr_decompos()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos.md)/[`hzr_decompos_g3()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_decompos_g3.md)
for the parameter-level early/late translations.

## Examples

``` r
hzr_argument_mapping()
#>    sas_statement               legacy_input                   c_concept
#> 1         HAZARD              TIME variable              obs time array
#> 2         HAZARD      EVENT/censor variable       event indicator array
#> 3         HAZARD          X covariate block               design matrix
#> 4         HAZARD         initial parameters            parameter vector
#> 5         HAZARD      baseline distribution phase distribution selector
#> 6         HAZARD            control options    optimizer/control struct
#> 7         HAZARD  additional legacy options        misc legacy switches
#> 8           TIME                          t                 time vector
#> 9          EVENT                     status                event vector
#> 10         PARMS                     theta0               starting coef
#> 11          DIST                       dist               dist selector
#> 12        HAZARD phases (3-phase structure)    3-phase Early/Const/Late
#> 13        HAZARD           MU_1, MU_2, MU_3     per-phase scale factors
#> 14            G1        THALF / RHO (early)             early half-life
#> 15            G1                 NU (early)         early time exponent
#> 16            G1                  M (early)                 early shape
#> 17            G1              DELTA (early)        early time transform
#> 18            G2          G2 constant phase  constant hazard rate phase
#> 19            G3                 TAU (late)               late G3 scale
#> 20            G3               GAMMA (late)       late G3 time exponent
#> 21            G3               ALPHA (late)               late G3 shape
#> 22            G3                 ETA (late)      late G3 outer exponent
#>                      r_parameter required                expected_type
#> 1                           time     TRUE               numeric vector
#> 2                         status     TRUE       numeric/logical vector
#> 3                              x    FALSE numeric matrix or data.frame
#> 4                          theta    FALSE               numeric vector
#> 5                           dist    FALSE             character scalar
#> 6                        control    FALSE                   named list
#> 7                            ...    FALSE              named arguments
#> 8                           time     TRUE               numeric vector
#> 9                         status     TRUE       numeric/logical vector
#> 10                         theta    FALSE               numeric vector
#> 11                          dist    FALSE             character scalar
#> 12  phases (list of hzr_phase())    FALSE    list of hzr_phase objects
#> 13 mu (via exp(log_mu) in theta)    FALSE          numeric (per-phase)
#> 14            hzr_phase(t_half=)    FALSE              positive scalar
#> 15                hzr_phase(nu=)    FALSE               numeric scalar
#> 16                 hzr_phase(m=)    FALSE               numeric scalar
#> 17        (absorbed by decompos)    FALSE               numeric scalar
#> 18         hzr_phase('constant')    FALSE        hzr_phase('constant')
#> 19         hzr_phase('g3', tau=)    FALSE              positive scalar
#> 20       hzr_phase('g3', gamma=)    FALSE               numeric scalar
#> 21       hzr_phase('g3', alpha=)    FALSE               numeric scalar
#> 22         hzr_phase('g3', eta=)    FALSE               numeric scalar
#>                                                                                transform_rule
#> 1                                                                     pass through as numeric
#> 2                                                                       coerce to numeric 0/1
#> 3                                                                   data.frame -> data.matrix
#> 4                                                 length must equal ncol(x) when x is present
#> 5                                                                 normalized lower-case label
#> 6                                                                      stored in spec$control
#> 7                                                            stored in legacy_args for parity
#> 8                                                                                pass through
#> 9                                                                           coerce to numeric
#> 10                                                                 map PARMS/INITIAL to theta
#> 11                                                                          map DIST= to dist
#> 12 list(early=hzr_phase('cdf',...), constant=hzr_phase('constant'), late=hzr_phase('g3',...))
#> 13                          exp(alpha_j) in internal parameterization; estimated on log scale
#> 14                                         maps directly to hzr_phase(t_half=) starting value
#> 15                                             maps directly to hzr_phase(nu=) starting value
#> 16                                              maps directly to hzr_phase(m=) starting value
#> 17                  time transform B(t) = (exp(delta*t)-1)/delta absorbed into decompos shape
#> 18                                             hzr_phase('constant') with no shape parameters
#> 19                                      maps directly to hzr_phase('g3', tau=) for late phase
#> 20                                    maps directly to hzr_phase('g3', gamma=) for late phase
#> 21                                    maps directly to hzr_phase('g3', alpha=) for late phase
#> 22                                      maps directly to hzr_phase('g3', eta=) for late phase
#>    implementation_status
#> 1            implemented
#> 2            implemented
#> 3            implemented
#> 4            implemented
#> 5            implemented
#> 6            implemented
#> 7            implemented
#> 8            implemented
#> 9            implemented
#> 10               planned
#> 11           implemented
#> 12           implemented
#> 13           implemented
#> 14           implemented
#> 15           implemented
#> 16           implemented
#> 17           implemented
#> 18           implemented
#> 19           implemented
#> 20           implemented
#> 21           implemented
#> 22           implemented
#>                                                                                                       notes
#> 1                                                                              Core observation time input.
#> 2                                      Event indicator currently retained as numeric in object$data$status.
#> 3                                              Future versions will support richer design encoding helpers.
#> 4                                                             Used by predict.hazard as coefficient vector.
#> 5                                                       Current default is 'weibull'; more options planned.
#> 6                                                 Control list is stored and reserved for optimizer parity.
#> 7                                              Supports legacy-style pass-through options during migration.
#> 8                                                               Canonical SAS migration uses TIME= mapping.
#> 9                                                              Canonical SAS migration uses EVENT= mapping.
#> 10                                                             SAS PARMS syntax parser not yet implemented.
#> 11                                                                  SAS DIST keyword maps directly to dist.
#> 12              Use dist='multiphase' with phases argument. N-phase generalization of legacy 3-phase model.
#> 13          Each phase has its own scale mu_j(x) = exp(alpha_j + x*beta_j). Starting value via hzr_phase().
#> 14                                 Half-life: time at which G(t_half) = 0.5. Same concept as SAS RHO/THALF.
#> 15                            Time exponent controlling rate dynamics. Same parameter name as SAS early NU.
#> 16                      Shape exponent controlling distributional form. Same parameter name as SAS early M.
#> 17          The C DELTA controlled B(t) = (exp(delta*t)-1)/delta. This transform is absorbed by decompos().
#> 18                                  Flat background rate. No shape parameters estimated. SAS G2 equivalent.
#> 19                                   Late-phase G3 scale parameter. Maps directly to hzr_phase('g3', tau=).
#> 20                                   Late-phase G3 time exponent. Maps directly to hzr_phase('g3', gamma=).
#> 21 Late-phase G3 shape parameter. alpha=0 gives exponential case. Maps directly to hzr_phase('g3', alpha=).
#> 22                                    Late-phase G3 outer exponent. Maps directly to hzr_phase('g3', eta=).
hzr_argument_mapping(include_planned = FALSE)
#>    sas_statement               legacy_input                   c_concept
#> 1         HAZARD              TIME variable              obs time array
#> 2         HAZARD      EVENT/censor variable       event indicator array
#> 3         HAZARD          X covariate block               design matrix
#> 4         HAZARD         initial parameters            parameter vector
#> 5         HAZARD      baseline distribution phase distribution selector
#> 6         HAZARD            control options    optimizer/control struct
#> 7         HAZARD  additional legacy options        misc legacy switches
#> 8           TIME                          t                 time vector
#> 9          EVENT                     status                event vector
#> 10          DIST                       dist               dist selector
#> 11        HAZARD phases (3-phase structure)    3-phase Early/Const/Late
#> 12        HAZARD           MU_1, MU_2, MU_3     per-phase scale factors
#> 13            G1        THALF / RHO (early)             early half-life
#> 14            G1                 NU (early)         early time exponent
#> 15            G1                  M (early)                 early shape
#> 16            G1              DELTA (early)        early time transform
#> 17            G2          G2 constant phase  constant hazard rate phase
#> 18            G3                 TAU (late)               late G3 scale
#> 19            G3               GAMMA (late)       late G3 time exponent
#> 20            G3               ALPHA (late)               late G3 shape
#> 21            G3                 ETA (late)      late G3 outer exponent
#>                      r_parameter required                expected_type
#> 1                           time     TRUE               numeric vector
#> 2                         status     TRUE       numeric/logical vector
#> 3                              x    FALSE numeric matrix or data.frame
#> 4                          theta    FALSE               numeric vector
#> 5                           dist    FALSE             character scalar
#> 6                        control    FALSE                   named list
#> 7                            ...    FALSE              named arguments
#> 8                           time     TRUE               numeric vector
#> 9                         status     TRUE       numeric/logical vector
#> 10                          dist    FALSE             character scalar
#> 11  phases (list of hzr_phase())    FALSE    list of hzr_phase objects
#> 12 mu (via exp(log_mu) in theta)    FALSE          numeric (per-phase)
#> 13            hzr_phase(t_half=)    FALSE              positive scalar
#> 14                hzr_phase(nu=)    FALSE               numeric scalar
#> 15                 hzr_phase(m=)    FALSE               numeric scalar
#> 16        (absorbed by decompos)    FALSE               numeric scalar
#> 17         hzr_phase('constant')    FALSE        hzr_phase('constant')
#> 18         hzr_phase('g3', tau=)    FALSE              positive scalar
#> 19       hzr_phase('g3', gamma=)    FALSE               numeric scalar
#> 20       hzr_phase('g3', alpha=)    FALSE               numeric scalar
#> 21         hzr_phase('g3', eta=)    FALSE               numeric scalar
#>                                                                                transform_rule
#> 1                                                                     pass through as numeric
#> 2                                                                       coerce to numeric 0/1
#> 3                                                                   data.frame -> data.matrix
#> 4                                                 length must equal ncol(x) when x is present
#> 5                                                                 normalized lower-case label
#> 6                                                                      stored in spec$control
#> 7                                                            stored in legacy_args for parity
#> 8                                                                                pass through
#> 9                                                                           coerce to numeric
#> 10                                                                          map DIST= to dist
#> 11 list(early=hzr_phase('cdf',...), constant=hzr_phase('constant'), late=hzr_phase('g3',...))
#> 12                          exp(alpha_j) in internal parameterization; estimated on log scale
#> 13                                         maps directly to hzr_phase(t_half=) starting value
#> 14                                             maps directly to hzr_phase(nu=) starting value
#> 15                                              maps directly to hzr_phase(m=) starting value
#> 16                  time transform B(t) = (exp(delta*t)-1)/delta absorbed into decompos shape
#> 17                                             hzr_phase('constant') with no shape parameters
#> 18                                      maps directly to hzr_phase('g3', tau=) for late phase
#> 19                                    maps directly to hzr_phase('g3', gamma=) for late phase
#> 20                                    maps directly to hzr_phase('g3', alpha=) for late phase
#> 21                                      maps directly to hzr_phase('g3', eta=) for late phase
#>    implementation_status
#> 1            implemented
#> 2            implemented
#> 3            implemented
#> 4            implemented
#> 5            implemented
#> 6            implemented
#> 7            implemented
#> 8            implemented
#> 9            implemented
#> 10           implemented
#> 11           implemented
#> 12           implemented
#> 13           implemented
#> 14           implemented
#> 15           implemented
#> 16           implemented
#> 17           implemented
#> 18           implemented
#> 19           implemented
#> 20           implemented
#> 21           implemented
#>                                                                                                       notes
#> 1                                                                              Core observation time input.
#> 2                                      Event indicator currently retained as numeric in object$data$status.
#> 3                                              Future versions will support richer design encoding helpers.
#> 4                                                             Used by predict.hazard as coefficient vector.
#> 5                                                       Current default is 'weibull'; more options planned.
#> 6                                                 Control list is stored and reserved for optimizer parity.
#> 7                                              Supports legacy-style pass-through options during migration.
#> 8                                                               Canonical SAS migration uses TIME= mapping.
#> 9                                                              Canonical SAS migration uses EVENT= mapping.
#> 10                                                                  SAS DIST keyword maps directly to dist.
#> 11              Use dist='multiphase' with phases argument. N-phase generalization of legacy 3-phase model.
#> 12          Each phase has its own scale mu_j(x) = exp(alpha_j + x*beta_j). Starting value via hzr_phase().
#> 13                                 Half-life: time at which G(t_half) = 0.5. Same concept as SAS RHO/THALF.
#> 14                            Time exponent controlling rate dynamics. Same parameter name as SAS early NU.
#> 15                      Shape exponent controlling distributional form. Same parameter name as SAS early M.
#> 16          The C DELTA controlled B(t) = (exp(delta*t)-1)/delta. This transform is absorbed by decompos().
#> 17                                  Flat background rate. No shape parameters estimated. SAS G2 equivalent.
#> 18                                   Late-phase G3 scale parameter. Maps directly to hzr_phase('g3', tau=).
#> 19                                   Late-phase G3 time exponent. Maps directly to hzr_phase('g3', gamma=).
#> 20 Late-phase G3 shape parameter. alpha=0 gives exponential case. Maps directly to hzr_phase('g3', alpha=).
#> 21                                    Late-phase G3 outer exponent. Maps directly to hzr_phase('g3', eta=).
```
