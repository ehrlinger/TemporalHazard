# Compute point predictions + delta-method CLs

Dispatched from predict.hazard() when `se.fit = TRUE`.

## Usage

``` r
.hzr_predict_with_se(
  object,
  type,
  time = NULL,
  x = NULL,
  x_list = NULL,
  cov_counts = NULL,
  phases = NULL,
  level = 0.95,
  diff_fn,
  conf_type = c("log-log", "logit")
)
```

## Arguments

- object:

  Fitted `hazard` object.

- type:

  One of "hazard", "linear_predictor", "survival", "cumulative_hazard".

- time:

  Numeric prediction times or NULL (for type = "hazard" /
  "linear_predictor").

- x:

  Design matrix (single-distribution) or NULL. Ignored for multiphase;
  use `x_list` instead.

- x_list:

  Named list of per-phase design matrices (multiphase only).

- cov_counts:

  Named integer covariate counts (multiphase only).

- phases:

  Named list of `hzr_phase` objects (multiphase only).

- level:

  Confidence level (default 0.95).

- diff_fn:

  Required function(theta) -\> numeric vector of length n returning the
  delta-method target (H, exp(eta), or eta depending on `type`). Used
  for the point estimate AND for the numeric jacobian fallback in exp /
  loglogistic / lognormal.

- conf_type:

  Survival CL transform: `"log-log"` (default, on `log(-log S)`) or
  `"logit"` (on `logit(1 - S)`, HAZPRED-compatible). Only consulted when
  `type == "survival"`.

## Value

data.frame with columns `fit`, `se.fit`, `lower`, `upper`.

## Details

### DELTA-METHOD TARGET

The quantity we actually differentiate depends on the prediction type
and the CL scale the type uses:

type = "cumulative_hazard" -\> target = H, log-scale CL type =
"survival" -\> target = H, log-log-survival CL (final `fit` is exp(-H);
reported `se.fit` is S \* se(H), the SE of S) type = "hazard" -\> target
= exp(eta) (single-dist only), log-scale CL type = "linear_predictor"
-\> target = eta, natural-scale CL

The caller supplies `diff_fn(theta)` that returns the target vector of
length n. For Weibull and multiphase we build J analytically; for exp /
loglogistic / lognormal we fall through to a numeric jacobian of
`diff_fn`.
