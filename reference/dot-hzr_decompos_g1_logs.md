# Log-scale terms shared by the m \> 0 branches of [`hzr_decompos()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos.md)

Cases 1 (`nu > 0`) and 3 (`nu < 0`) evaluate the same intermediate
quantities; only the sign of `rho` and the final assembly of `G` differ.
Forming them directly overflows well inside the range an optimizer can
reach: `2^m` is `Inf` from `m = 1024`, and `bt^(-1/nu)` overflows sooner
still (from `m ~ 750` at `t/t_half = 0.5`, `m*nu = 3`). `btnu` then
becomes `Inf` and `btnu^(-1/m)` collapses silently to 0.

## Usage

``` r
.hzr_decompos_g1_logs(time, t_half, nu, m)
```

## Arguments

- time:

  Numeric vector of positive times.

- t_half:

  Positive scalar.

- nu:

  Nonzero scalar.

- m:

  Positive scalar.

## Value

List with `log_S` (\\= -\log(\mathrm{btnu})/m\\) and `log_g`.

## Details

The `m` in `rho`'s `(2^m - 1)/m` factor cancels the explicit `m`
multiplier, so exactly \$\$m \\ b(t)^{-1/\nu} = (t\_{1/2}/t)^{1/\nu}
(2^m - 1),\$\$ and `log(btnu)` follows from a softplus of the log of
that product. Every intermediate below stays finite for any `m > 0` a
double can hold.
