# TemporalHazard 1.2.0

## Breaking changes

This release contains a breaking change but ships as a minor version. The
`1.x` line is the run-up to a first production release; the major digit is
reserved for that milestone rather than spent on a single changed default.
The change below is also closer to a correction than a redesign — the previous
default deviated from the SAS/C reference this package exists to reproduce.
Read the entry regardless: it can change which variables a stepwise run
selects.

* `hzr_stepwise()` now defaults to `criterion = "score"`, reproducing SAS/C
  HAZARD's `SELECTION` statistic. Previously it defaulted to `"wald"`, which
  refit the model once per candidate and used the refit's Wald chi-square --
  a deviation from the reference implementation this package exists to
  reproduce. **Re-running an existing stepwise analysis can now select a
  different variable set**, because the score and Wald paths take different
  step sequences. Pass `criterion = "wald"` to restore the previous behaviour
  exactly.

  The score criterion also removes the per-candidate refit, which dominated
  runtime: a 92-variable two-phase screen fell from roughly 25 minutes per
  bootstrap replicate to seconds.

  Following SAS, the variance used during *selection* is approximate --
  shaping-parameter covariances are ignored. Final-model standard errors are
  unchanged and still use the full Hessian.

  Score is an *entry* criterion. The drop path never refit per candidate, so
  removals are still tested on the current model's Wald p-value against
  `slstay`, as SAS does; drop rows in `$steps` are labelled `"wald"`
  accordingly.

## New features

* The SAS `.lst` parsers now ship with the installed package, under
  `sas-parity/` (`inst/sas-parity/` in the source tree -- `R CMD INSTALL`
  strips the `inst/` prefix). They previously lived in `tests/testthat/`,
  which `R CMD INSTALL` skips unless `--install-tests` is passed -- so a plain
  `install.packages()` or `remotes::install_github()` left them unreachable,
  and a downstream analysis wanting to check its own SAS output against them
  had to clone the repository. Reach them with:

  ```r
  source(system.file("sas-parity", "helper-sas-parity.R",
                     package = "TemporalHazard"))
  ```

  The parsers themselves are unchanged; only their location is. The package's
  own parity tests load them through a shim at
  `tests/testthat/helper-sas-parity.R`, so testthat's helper auto-sourcing
  still applies and no test file changed.

  These functions remain internal (`.hzr_`-prefixed) and unexported. They
  parse a specific vintage of SAS HAZARD listing output and carry no API
  stability guarantee.

* `hzr_bootstrap(verbose = TRUE)` now shows a text progress bar over the
  bootstrap replicates (via `utils::txtProgressBar()`) instead of an
  every-50-replicates message.

* `hzr_bootstrap()` gains a `scope` argument for embedded stepwise variable
  selection during each bootstrap replicate -- the R equivalent of SAS's
  `%HAZBOOT` procedure. Each replicate runs a fresh `hzr_stepwise()`
  selection (starting from a fixed-shape refit of the base model) instead
  of a plain refit, so `summary$pct` reports the variable's selection
  frequency across resamples and `summary$mean`/`sd`/`ci_*` describe the
  coefficient distribution conditional on selection. `scope = NULL`
  (the default) preserves the original fixed-formula bootstrap unchanged.

## Bug fixes

* `hzr_bootstrap()` now resamples fits built with the **vector interface**
  (`time =` / `status =` rather than a formula plus `data`). Previously it
  resampled `data` only, but a vector-interface call stores `time = d$col` as
  an *expression*, so every replicate re-evaluated it against the original
  data and returned the original fit. The result was `n_success = n_boot`,
  `n_failed = 0`, no warning, and `n_boot` **identical** replicates -- a
  summary table that looked complete and contained nothing, with `sd` exactly
  0 on every parameter. The evaluated `time`, `status`, `time_lower` and
  `time_upper` vectors are already stored on the fitted object, so they are
  now resampled by the same index as the rows and rewired into each
  replicate's call, exactly as `data` and `weights` already were. Both
  interfaces now produce identical bootstrap replicates for the same model,
  data and seed. Found running a 500-replicate production bagging job that
  completed in 9.5 minutes and produced no usable output.

* A fit that cannot compute a Hessian now says so. The analytic Hessian
  declines for left- and interval-censored rows by design, the optimizer falls
  back to `numDeriv::hessian()`, and `numDeriv` is a `Suggests` -- so on a
  machine installed without Suggests, an interval-censored multiphase fit
  produced no standard errors, `rcond = NA`, `pd = NA` and a `vcov()` of bare
  `logical`, with nothing naming the cause. The user-visible symptom was
  `diag(vcov(fit))` reporting an invalid `'nrow'`, which is unrecognisable
  from the cause. Three paths now warn: `numDeriv` absent (naming it and the
  install command), `numDeriv::hessian()` failing (carrying its message), and
  no Hessian available at all. A `hessian_fn` hook that *errors* is also no
  longer swallowed into silence, so a broken analytic hook is distinguishable
  from one that deliberately declines. Behaviour is unchanged -- the
  diagnostics are still `NA` -- but the reason is now stated. Found while
  fitting a production interval-censored study.


* `hzr_bootstrap()` no longer returns a silent `n_success = 0` (and
  `n_failed = n_boot`, with no error and no warning) when the model was fitted
  inside a function. `hazard()` stored its call but not the environment that
  call was written in, so each replicate's refit resolved arguments passed by
  symbol -- `theta`, `phases`, `control` -- against the package namespace and
  `globalenv()` rather than the caller's locals. Fits built at the top level
  appeared to work by falling through to `globalenv()`; fits built inside a
  function failed on every replicate, and the per-replicate `tryCatch()`
  swallowed the error. `hazard()` now records the fitting environment, and each
  replicate is evaluated in a child of it that carries the resampled data and
  weights. Affects both `refit` and `select` modes.

* **`hzr_bootstrap(scope = ...)` selected nothing when the base fit's formula
  was passed by symbol.** `hazard()` records its call with `match.call()`, so a
  formula assigned to a variable first (`f <- Surv(t, d) ~ 1; hazard(f, ...)`)
  is stored as a *symbol* rather than a call. The scope-mutating refit
  recovered it with `as.formula(deparse(...))`, which turns that symbol into
  the string `"f"` and errors with `invalid formula "f": not a call`. Every
  post-entry refit therefore failed, no candidate ever entered, and the run
  reported `n_success = n_boot`, `n_failed = 0`, no error and no warning --
  with a summary holding only the base model's parameters. The stored formula
  is now evaluated in the fit's recorded calling environment, which handles
  the literal and by-symbol forms alike, and a stored formula that fails to
  resolve raises an error naming the problem instead of degrading to an empty
  screen. The same defect affected `hzr_stepwise()` directly. (#114)

* **A select-mode `hzr_bootstrap()` run that selects no covariate now warns.**
  The base model's own parameters appear in every replicate by construction,
  so they fill the summary at `pct = 100` and an empty screen reads as a set
  of perfectly reliable variables; nothing in the output prompted the reader
  to compare the parameter names against `names(coef(object))`. The warning
  names the likely causes: an entry criterion stricter than intended, a
  `scope` naming columns absent from the data, or a base fit whose stored call
  cannot be rewritten. Legitimate empty screens warn too -- an entry criterion
  no candidate can clear is also worth reporting. (#115)

* `hzr_bootstrap()` no longer floods the console with per-replicate numerical
  warnings (e.g. ill-conditioned-Hessian notes from unstable resamples), which
  are not individually actionable when the bootstrap aggregates over replicates.
  Structural problems (a mistyped `scope` column, an invalid scope) still
  surface once, up front.

* `hzr_bootstrap(scope = ..., trace = ...)` no longer errors with "formal
  argument matched by multiple actual arguments". Select-mode forwarded
  `...` to `hzr_stepwise()` alongside an explicit `trace = FALSE`, so any
  caller-supplied `trace=` collided with it.

* Multiphase models with a `"cdf"`/`"hazard"` phase whose shape sits exactly
  at the `m = 0` (Case 3L) or `nu = 0` (Case 2L) limiting-case boundary no
  longer lose their analytic Hessian. The finite-difference second
  derivative used to probe the *other* shape parameter's `-h` side, which
  can cross into the mathematically undefined `m < 0 && nu < 0` region and
  raise an error; this silently fell back to a numerical Hessian (or, if
  that also failed to invert, to `NA` standard errors) for every affected
  fit, not just `hzr_bootstrap()`'s Conservation-of-Events full-information
  recompute. The boundary direction now uses a one-sided finite difference
  instead.

* Multiphase fits with a single free parameter (a two-phase model with all
  shapes fixed, where Conservation of Events fixes one of the two `log_mu`)
  now use the analytic Hessian for standard errors instead of silently
  falling back to a numerical one. Restricting the Hessian to the lone free
  parameter dropped it from a 1x1 matrix to a scalar, which was rejected as
  non-conformant; it is now kept as a matrix (`drop = FALSE`).

# TemporalHazard 1.1.0

## New features

* `predict.hazard(type = "hazard")` now works for **multiphase** models,
  returning the instantaneous additive hazard
  `h(t|x) = sum_j mu_j(x) phi_j'(t)` (previously only single-distribution
  models supported `"hazard"`, via `exp(eta)`). Like `"survival"` /
  `"cumulative_hazard"` it is time-based (requires `newdata$time`), supports
  covariate `newdata`, and `se.fit = TRUE` (delta-method limits on the log
  scale via a numeric Jacobian of the hazard evaluator). `decompose = TRUE` is
  not supported for `"hazard"`. This gives the multiphase instantaneous hazard a
  public route (it was previously reachable only through internal functions).

* `predict.hazard(..., se.fit = TRUE, conf.type = "logit")` selects the survival
  confidence-limit transform. The default `"log-log"` builds limits on
  `log(-log S)` (the `survival::survfit` standard); `"logit"` builds them on
  `logit(1 - S)`, reproducing SAS HAZARD's `HAZPRED` survival limits. With the
  full-information vcov for CoE fits, `conf.type = "logit"` matches the SAS
  `hp.death.AVC` survival CLs to ~1e-5. Hazard / cumulative-hazard limits are
  unaffected (their log scale already matches HAZPRED).

* `predict.hazard(type = "cumulative_hazard", decompose = TRUE, se.fit = TRUE)`
  now returns per-phase **and** total delta-method confidence limits for
  multiphase models, as a long data frame
  (`time`, `component`, `fit`, `se.fit`, `lower`, `upper`). Each phase's CL uses
  only that phase's parameters, so per-phase limits do not sum to the total.
  Previously this combination raised an error.

## Changes

* **`hzr_deciles()` now matches the SAS `deciles.hazard` macro exactly.**
  Previously it excluded subjects censored before the horizon and defined the
  expected count as `sum(1 - S(horizon))`. It now follows the SAS method: **all**
  subjects are ranked into equal-sized risk groups by predicted survival at the
  horizon, and the expected count per group is the **sum of predicted cumulative
  hazard at each subject's own follow-up time** (so group totals sum to the total
  observed events under conservation of events). The `time` argument now only
  stratifies subjects into risk groups; it no longer restricts or excludes any
  subject, and the expected/observed totals are horizon-independent. Verified to
  reproduce the `hm.death.AVC.deciles` SAS decile table (CASES/EXPECTED/ACTUAL)
  to print precision. The output columns are unchanged; their definitions are
  updated in `?hzr_deciles`.

## Bug fixes

* **Conservation-of-Events fits now report the full-information variance.**
  CoE removes one phase's `log_mu` from the optimizer *search* (its score
  equation is the CoE constraint), but the previous code also dropped it from
  the *uncertainty* -- the conserved phase got an `NA` standard error, and
  anything depending on it (other SEs, `se(H)`, prediction confidence limits)
  was understated wherever that phase contributed. At the optimum the CoE
  solution is the unconstrained MLE, so `vcov()` is now recomputed from the
  unconstrained-objective Hessian over the full free set (including the
  conserved `log_mu`), matching an all-`mu`-free (`conserve = FALSE`) fit at the
  same point. On `hz.death.AVC` every parameter SE now matches the SAS HAZARD
  reference (e.g. the conserved early `log_mu`: 0.133 vs the previous ~0.059).
  The recomputation uses `numDeriv` (Suggests) and an invertible Hessian; if
  either is unavailable the fit emits a warning and the conserved `log_mu`
  retains an `NA` standard error (as before).

* **Conservation of Events ignored left-truncation (counting-process entry
  times).** For multiphase fits on `Surv(start, stop, event)` data, the CoE
  reparameterization conserved `Sum H(stop)` while the likelihood scores the
  intercepts on the entry-time scale, `Sum E = Sum [H(stop) - H(start)]`. The
  conserved phase therefore absorbed the spurious `Sum H(start)`, biasing its
  intercept and lowering the attained log-likelihood (the `hz.te123.OMC` fit-1
  parity offset, gap-list P1 #6). `.hzr_conserve_events()` and
  `.hzr_select_fixmu_phase()` now subtract the per-phase entry-time cumulative
  hazard, matching the likelihood and C HAZARD `setcoe` under `LCENSOR`/
  `STARTTME`. Plain right-censored fits (no `start` time) are unaffected.

* **`vcov()` was unusable for multiphase fits and returned an unnamed matrix.**
  `vcov.hazard()` collapsed the entire matrix to a scalar `NA` whenever any cell
  was `NA`. Multiphase fits legitimately have `NA` variance rows -- for
  parameters held fixed (e.g. early shapes) and for the
  Conservation-of-Events-conserved phase `log_mu` -- so the finite
  free-parameter block was discarded for almost every multiphase model. The
  method now returns the full matrix with `NA` rows preserved and labels rows
  and columns with the coefficient names (phase-prefixed for multiphase, e.g.
  `early.x` vs `constant.x`), so a covariate shared across phases resolves to
  distinct, name-addressable slots. A scalar `NA` is returned only when no
  covariance matrix is available.

* **Weibull analytic gradient produced `NaN` for right-censored `time = 0` rows.**
  `.hzr_gradient_weibull()` used an unguarded `log(time)` in the shape (`nu`)
  score; a legal right-censored row at `time = 0` made `0 * -Inf = NaN`, which
  poisoned the entire summed shape-gradient component (then silently zeroed by
  the optimizer, harming convergence). `log(time)` is now guarded with
  `log(pmax(time, .Machine$double.xmin))`, matching the analytic Hessian. The
  other families were audited: exponential (no `log(time)` in the score),
  log-normal (rejects `time = 0`), and multiphase (the decomposition clamps
  `time`) are unaffected.

* **Weibull event hazard was inconsistent with its cumulative hazard.**
  `.hzr_logl_weibull()` defined the event hazard as `mu*nu*t^(nu-1)*exp(eta)`
  while the cumulative hazard was `(mu*t)^nu*exp(eta)`; the former is missing a
  `mu^(nu-1)` factor (the exact derivative is `nu*mu^nu*t^(nu-1)*exp(eta) =
  (nu/t)*H`, Form A as in the C/SAS HAZARD reference). The natural-scale
  log-likelihood and its analytic gradient (`d/dmu`, `d/dnu` event terms) are
  corrected to match. Pure event/right-censored fits were already correct (they
  use the self-consistent internal reparameterization); the visible effect is on
  **mixed event + interval/left-censored Weibull fits**, which delegate to this
  likelihood and previously optimized a slightly mis-specified event term.

* **Weibull gradient attribute ignored observation weights.**
  `.hzr_logl_weibull(..., return_gradient = TRUE)` attached an unweighted
  gradient even when `weights` were supplied (the analytic gradient was off by
  the weight scale, e.g. halved under `weights = 2`). `weights` is now forwarded
  to the score computation. The model-fitting path was unaffected (it uses a
  separate internal weighted gradient); this only changes callers reading the
  `return_gradient = TRUE` attribute on weighted data.

* **`hzr_bootstrap()` was non-functional for weighted fits** (Phase 7c).
  The resample loop rewired only `data` in the refit call, leaving the
  original `weights` argument bound to a symbol in the *caller's* frame.
  The internal `eval()` could not resolve that symbol, so **every** replicate
  of a weighted model errored out (`n_success == 0`) regardless of `fraction`;
  even had it resolved, the un-resampled weights would have been misaligned
  with the bootstrapped rows.  `weights` is now evaluated once and resampled
  in lockstep with the data on each replicate (mirroring how `data` is
  handled).  Unweighted bootstraps are unaffected.  A regression test covers
  both the `fraction < 1` and full-size weighted paths in
  `test-diagnostics.R`.  Follow-up: `hzr_bootstrap()` now resamples the
  weights already stored on the fitted object (`object$data$weights`) rather
  than re-evaluating the call's `weights` expression in `parent.frame()`,
  which fails when the original symbol is no longer in scope (e.g. the fit
  was built inside a helper that has returned).  Caller-frame evaluation
  remains a fallback for objects fitted before weights were stored.  The same
  fragility applied to the call's `data` argument: `hazard()` now stores the
  evaluated `data` argument (the data frame passed to `hazard()`, not a
  `model.frame()` result) on the fitted object (`object$data$frame`), and
  `hzr_bootstrap()` resamples that stored frame instead of re-evaluating
  `cl$data` in `parent.frame()`, so bootstrap succeeds even when the original
  `data` symbol is out of scope.  Caller-frame evaluation remains a fallback
  for objects fitted before the frame was stored.

* **4-phase CoE fixmu-phase selection** (Phase 7d).
  `.hzr_select_fixmu_phase()` used `which.max()` over raw per-phase cumhaz
  at the starting theta.  G3 late phases with typical shape parameters have
  unnormalized cumhaz orders of magnitude larger than other phases, causing
  CoE to pin the G3 `log_mu` away from its true near-zero MLE.  Fixed by
  excluding phases whose cumhaz contribution exceeds 10× the median before
  selecting (falls back to `which.max` when all phases are outliers).  On the
  4-phase CABGKUL fit the CoE vs no-CoE LL gap closes from 6.9 to < 0.1
  units.  Six new tests cover the 4-phase code path in
  `test-conservation-of-events.R`.
* **`time_lower` dual-use bug in Weibull and multiphase likelihoods.**
  When `time_lower` was supplied for a mixed interval-censored + right-censored
  dataset, the Weibull LL interpreted `time_lower` as the counting-process
  *entry time* for right-censored rows, computing H(stop) − H(start) = 0 and
  silently zeroing those rows' likelihood contribution.  Fixed in
  `likelihood-weibull.R` (4 sites: LL, gradient, L-BFGS-B internal LL/gradient)
  and `likelihood-multiphase.R`: `start_vec` is now set from `time_lower` only
  for genuine epoch rows (`status %in% c(0L, 1L)` and `time_lower < time`).
  Two regression tests added to `test-interval-censoring-weibull.R`.

* **`hzr_decompos()` Case 3 corrected and `nu = 0, m >= 0` now fails loud**
  (Phase 7d).  Two issues in the early-phase (G1) sign dispatch:
    - **Case 3 (`m > 0, nu < 0`, "bounded cumulative") carried a spurious
      factor of `m`.**  Its `rho` used a bare `(2^m - 1)^nu` instead of the
      `((2^m - 1)/m)^nu` form used by Case 1, leaving an `m` factor on the
      `bt^(-1/nu)` term.  The CDF diverged from the C HAZARD G1 evaluator
      (`g1flag = 5`) by up to ~0.2 and was discontinuous with its `m -> 0`
      limit (Case 3L).  Adding the `/m` divisor makes the `m` factors cancel,
      reproducing the C evaluator exactly and restoring continuity (verified
      against `src/common/hzd_ln_G1_and_SG1.c`).  No shipped phase uses
      Case 3, so fitted models are unaffected; the synthetic 3-phase golden
      fixture was regenerated because its free-shape optimizer path crosses
      Case 3 territory.
    - **`nu = 0` with `m >= 0`** fell through every dispatch branch, leaving
      the CDF unassigned and raising the cryptic `object 'G' not found`.  The
      `nu -> 0` limit is defined only for `m < 0`; for `m >= 0` it is
      degenerate.  The function now raises a clear, explanatory error.
  New `test-decompos-boundary.R` locks in continuity of all limiting branches
  (Case 1 -> 1L, 2 -> 1L, 2 -> 2L, 3 -> 3L), Case 3 <-> C `g1flag=5` parity,
  `g = dG/dt` internal consistency, CDF sanity, and stability at extreme
  `t_half`.

## Improvements

* **Hardened Hessian inversion for standard errors (Phase 7c).**
  Post-fit variance-covariance estimation now symmetrizes the Hessian,
  checks its reciprocal condition number, inverts via Cholesky with a
  `solve()` fallback for non-positive-definite Hessians, and guards
  non-positive variances instead of silently emitting `NaN` standard
  errors. Ill-conditioned, non-positive-definite, and non-finite Hessians
  now raise specific, named warnings, and fits carry `rcond` / `pd`
  diagnostics that `summary()` surfaces as a note when a fit is flagged.
  This closes the "12+-parameter Hessian stability" hardening item for the
  inversion layer; analytic Hessians (more accurate standard errors) follow
  in subsequent releases.

* **Analytic Hessian for exponential standard errors (Phase 7c, Layer 2).**
  The exponential distribution now computes its post-fit Hessian in closed form
  (`X~' diag(wH) X~` over event + right-censored rows) rather than numerically,
  giving more accurate standard errors. The shared optimizer gained a
  `hessian_fn` hook that analytic Hessians for the remaining families will reuse;
  left/interval-censored exponential fits fall back to the numerical Hessian.
* **Analytic Hessian for Weibull standard errors (Phase 7c, Layer 2).**
  The Weibull distribution now computes its post-fit Hessian in closed form on
  the internal `(alpha, psi, beta)` optimization scale (then mapped to the
  natural scale by the existing delta method) rather than numerically, giving
  more accurate standard errors. Covers event + right-censored data (including
  counting-process start times); left/interval-censored fits fall back to the
  numerical Hessian.
* **Analytic Hessian for log-logistic standard errors (Phase 7c, Layer 2).**
  The log-logistic distribution now computes its post-fit Hessian in closed form
  on the internal `(log alpha, log beta, beta_coef)` scale rather than numerically,
  giving more accurate standard errors. Covers event + right-censored data;
  left/interval-censored fits fall back to the numerical Hessian.

* **Analytic Hessian for log-normal standard errors (Phase 7c, Layer 2).**
  The log-normal distribution now computes its post-fit Hessian in closed form
  on the internal `(mu, log_sigma, beta_coef)` scale rather than numerically,
  giving more accurate standard errors. Covers event + right-censored data;
  left/interval-censored fits fall back to the numerical Hessian.

* **Analytic Hessian for multiphase standard errors (Phase 7c, Layer 2 PR-6).**
  Post-fit standard errors for all multiphase fits now come from a closed-form
  Hessian of the negative log-likelihood rather than a numerical Richardson
  approximation. The Hessian is assembled from three terms: (A) a
  phase-block-diagonal curvature of Σᵢ wᵢ H(tᵢ), (B) a dense Fisher
  information outer product Σₑ (wᵢ/hᵢ²) ∇h ∇hᵀ capturing cross-phase
  parameter interactions, and (C) a phase-block-diagonal curvature of
  −Σₑ wᵢ log h(tᵢ). μ/β parameters use fully closed-form expressions;
  shape parameters (t_half, ν, m, and G3 parameters) use second-order
  central differences. The Conservation-of-Events full-information vcov
  path also switches to the analytic Hessian.
  Left/interval-censored fits fall back to the numerical Hessian.
  Completes the 6-PR analytic-Hessian rollout across all five families.

## Documentation

* `vignette("fitting-hazard-models")` gains an **Interval and left censoring**
  section covering: status coding reference (`-1`/`0`/`1`/`2`), a cardiac
  clinic-visit simulation with right- and interval-censored observations,
  the direct `time_lower`/`time_upper` API, and a comparison showing the
  interval-censored fit recovering `nu` close to 1.0 (true value) while the
  naive exact-at-upper fit incurs a shape bias of ~+0.45.  Includes a callout note on the correct
  use of `time_lower = 0` for right-censored rows.
* `vignette("fitting-hazard-models")` gains a **Convergence troubleshooting**
  section covering: reading the KM cumulative hazard for Weibull starting
  values (log-log plot), when to fix shape parameters vs. estimate freely,
  diagnosing overparameterization via near-zero phase scales and `NA` from
  `vcov()`, and `control` options (`n_starts`, `maxit`).
* Added a package-level overview help page (`?TemporalHazard`) giving the
  additive multiphase model, the phase-type vocabulary, the SAS/C HAZARD
  bridge, and a map of the main entry points.
* Expanded the mathematical content of the core help files in the style of
  `randomForestSRC`: explicit display equations for the generalized temporal
  decomposition `G(t)` (`?hzr_decompos`), the additive cumulative-hazard model
  on `?hzr_phase` and `?hazard`, and defining formulas plus the
  Mächler (2012) reference for the numerical primitives (`?hzr_log1pexp`,
  `?hzr_log1mexp`, `?hzr_clamp_prob`).
* Added methodological references to the nonparametric diagnostics
  (Kaplan-Meier/Greenwood, Nelson-Aalen, Aalen-Johansen) and filled in missing
  cross-references across the exported help pages.
* Explained the remaining enumerated options in the style of the `hzr_phase()`
  phase-type help. `?hazard` gains a **Baseline distributions** section
  describing each `dist` value (`"weibull"`, `"exponential"`, `"loglogistic"`,
  `"lognormal"`, `"multiphase"`) by its hazard shape and when to use it;
  `?hzr_stepwise` gains a **Selection direction and criterion** section
  explaining each `direction` (`"forward"`/`"backward"`/`"both"`) and
  `criterion` (`"wald"`/`"aic"`), including how Wald selection differs from
  C/SAS HAZARD's score-statistic path.

## Testing

* **Patient-specific HAZPRED prediction parity** (Group A fixtures
  `hp.death.AVC.hm1` / `hm2`).  New `test-sas-parity.R` blocks predict survival
  and instantaneous hazard -- with logit survival CLs and log hazard CLs at the
  SAS 1-SD level -- from the saved multivariable both-phase model
  (`hm.death.AVC` final fit, "HMDEATH") for two covariate profiles each
  (hm1: with/without an associated cardiac anomaly; hm2: complete vs partial
  canal by date of repair), matching SAS to ~5e-4 (survival) / ~8e-3 (hazard;
  the looser hazard tolerance reflects the near-singular 9-coefficient fit and
  the steep early-phase times).  Adds a header-driven
  `.hzr_parse_sas_nomogram_mv()` (parses the BY-group "digital nomogram" whose
  rows each carry their own covariate vector) and a shared
  `.hzr_fit_avc_hmdeath()` helper.

* **Stratified HAZPRED calibration parity** (Group A fixture
  `hs.death.AVC.hm1`).  New `test-sas-parity.R` blocks reproduce the
  population-averaged, stratified-by-`COM_IV` outputs from the same HMDEATH
  model: (1) the observed-vs-expected "predict number of deaths" table --
  per stratum, EXPECTED = sum of predicted cumulative hazard at each subject's
  own follow-up, PEXPECT = sum of predicted death probability, ACTUAL =
  observed deaths (totals conserve events, 14.76 + 55.24 = 70), to ~5e-3; and
  (2) the per-stratum mean survival curve (MSURVIV) at the digital time grid, to
  ~5e-4.  Adds `.hzr_parse_sas_calibration()` and
  `.hzr_parse_sas_strata_survival()`.

* **`hm.death.AVC` stepwise documented as a non-parity gap** (Group A).  The
  phase-aware forward `SELECTION SLE=0.2 SLS=0.1` fit's *final* selected model
  is the saved "HMDEATH" fit already verified by the `hm.death.AVC.deciles` /
  `hp.death.AVC.hm1` / `hm2` parity tests; its *selection path* cannot be
  reproduced (SAS uses approximate variances during selection while R's full
  Hessian is near-singular here; SAS's `/I` `/S` flags are phase-level but R's
  `force_in` is phase-blind; R oscillates at p ~ slstay and lands in a worse
  basin -- the same divergence already documented for `hm.deadp.VALVES`).
  `test-sas-parity.R` gains a regression-guard test that exercises the
  multiphase phase-aware stepwise path end-to-end on real data without
  asserting path parity; see `inst/dev/FIXTURE-GAP-LIST.md`.

* **`bs.death.AVC` bootstrap documented as a non-parity gap** (Group A).  SAS
  `%HAZBOOT` runs a fresh stepwise selection on each bootstrap resample and
  reports a variable-selection frequency; R's `hzr_bootstrap()` resamples and
  refits a *fixed* model (no embedded-selection mode), and reimplementing the
  SAS procedure would inherit the documented `hm.death.AVC` stepwise
  divergence.  `test-sas-parity.R` adds `.hzr_parse_sas_bootstrap()` and asserts
  the SAS reference selection frequencies in parseable form (so the parity test
  is half-written for a future bootstrap-with-selection capability), plus a
  regression guard that R's fixed-model bootstrap runs on the cohort; see
  `inst/dev/FIXTURE-GAP-LIST.md`.

* **Phase-specific covariate recovery tests** (Phase 7d).  New
  `test-phase-specific-covariates.R` confirms that `hzr_phase(formula = ~ ...)`
  is correct, not just runnable: simulation-based recovery tests verify that a
  covariate entered into one phase recovers its true coefficient, that the same
  covariate carries independent (here opposite-sign) effects across two phases,
  and that a covariate confined to one phase does not leak into another.  This
  is the honest substitute for a SAS parity fixture and guards against the
  "accepts the formal but never applies it" regression that has surfaced
  before with weights and counting-process times.
* Added fractional (non-integer) weight coverage to close the roadmap 7a gap.
  Prior weight tests verified weighting only via integer row duplication, which
  cannot express fractional (e.g. inverse-probability) weights. The new tests
  assert the two properties that define a correct per-row weighted
  log-likelihood: an **additive split** (a row of weight `a + b` equals two
  identical copies of weights `a` and `b`) and **linear scaling**
  (`L(theta; c*w) = c * L(theta; w)`, gradient likewise, MLE invariant), across
  the Weibull, exponential, and multiphase-with-covariates paths.
* Made the single-distribution weighted-fit tests exercise a real fit. They
  previously omitted `theta` start values, so `hazard(fit = TRUE)` took its
  unfitted branch and the assertions compared `NULL`/`NA` vacuously; they now
  supply starts and genuinely compare the weighted MLE to the duplicated-row
  MLE.
* Added interval-censoring coverage under the multiphase model (roadmap 7c).
  The multiphase likelihood's interval-/left-censored branch had a working code
  path but no isolated test. New R-only self-consistency invariants in
  `test-interval-censoring-multiphase.R` verify the interval contribution equals
  `log(S(lower) - S(upper))`, the left-censored term equals
  `log(1 - exp(-H(u)))`, right-censoring stays `-(H(stop) - H(start))`
  (including left truncation), invalid bounds (`lower > upper`) yield `-Inf`,
  integer weights match row duplication on interval rows, and an
  interval-censored multiphase fit converges.
* Added a SAS fractional-weight parity capture scaffold under
  `inst/extdata/weights-fixtures/` (roadmap 7a / FIXTURE-GAP-LIST B5): a
  `PROC HAZARD ... WEIGHT IPW` template, a deterministic non-integer weight
  dataset, a `.lst` parser, and `test-weights-sas-parity.R`. The parity test
  re-fits the SAS specification in R and compares covariate estimates and
  log-likelihood; it skips when the capture fixture is absent (as it is by
  default), so CI and installation are unaffected until a SAS run is dropped
  in. R-side fractional-weight correctness is already proven by the invariants
  above; this is the drop-in external SAS confirmation.

---

# TemporalHazard 1.0.3

## Bug fixes / CRAN compliance

* `hzr_bootstrap()` no longer touches `.GlobalEnv` directly. The 1.0.2
  `oldseed`/`on.exit()`/`assign(".Random.seed", ...)` save-restore wrapper
  added in 1.0.2 violated CRAN policy on writing to `.GlobalEnv` and has
  been removed. When `seed` is supplied the function simply calls
  `set.seed(seed)` (the documented R API for seeded reproducibility); the
  `@param seed` documentation now notes that the caller's RNG state is not
  restored on exit. With `seed = NULL` (the default) the function does
  not call `set.seed()` at entry, so it starts from the caller's current
  RNG state; the bootstrap still consumes random numbers and advances
  that state in the usual way.

# TemporalHazard 1.0.2

## Bug fixes / CRAN compliance

* The golden-fixture generators (`.hzr_create_*_golden_fixture()`,
  previously `R/golden_fixtures.R`) have been moved out of the package to
  `data-raw/golden_fixtures.R`. They are maintainer-only helpers for
  regenerating the bundled `inst/fixtures/*.rds` reference outputs and are
  not part of the installed package, so they are no longer shipped, checked,
  or user-reachable. This resolves the home-filespace concern at its root:
  the earlier fallback resolved to `system.file("fixtures", ...)` — i.e. the
  installed package directory — whenever the package was installed, so the
  1.0.1 "falls back to `tempdir()`" fix did not actually prevent writing to
  the user library. The bundled `.rds` fixtures still ship and the parity
  tests still read them via `system.file()`.
* `.hzr_generate_golden_fixture()` (the C-binary reference writer in
  `R/parity-helpers.R`, which shares a file with test-time helpers and so
  was kept in the package) now takes a required `output_dir` argument with
  no default path.
* Removed the remaining hardcoded `seed = 42` literals from the relocated
  generators; recorded fixture metadata reflects the actual `seed` argument
  passed (`NULL` by default, so no seed is set inside the function).
* `hzr_bootstrap()` no longer leaves the caller's random-number stream
  altered when `seed` is supplied: the global `.Random.seed` is saved before
  `set.seed()` and restored via `on.exit()`, matching the fixture generators.
  Bootstrap reproducibility under a given `seed` is unchanged.

# TemporalHazard 1.0.1

## Bug fixes / CRAN compliance

* Added `\value` documentation to all exported functions that were missing it:
  `hazard()`, `coef.hazard()`, `vcov.hazard()`, `print.hzr_calibrate()`,
  `print.hzr_deciles()`, `print.hzr_gof()`, and `print.hzr_kaplan()`.
* Internal fixture generators (`R/golden_fixtures.R`) no longer set a specific
  seed unconditionally. Generators now accept an optional `seed` argument;
  when provided, the global RNG state is saved and restored via `on.exit()`.
* Default `output_dir` for fixture generators falls back to `tempdir()` instead
  of the package source directory, keeping the home filespace unmodified.

# TemporalHazard 0.9.8

## New features

* **Delta-method confidence limits on `predict.hazard()`** — Phase 4g of
  the development plan lands. Two new arguments: `se.fit = FALSE` and
  `level = 0.95`. When `se.fit = TRUE`, the return value becomes a
  data frame with columns `fit`, `se.fit`, `lower`, `upper`.
  - **Weibull and multiphase use closed-form Jacobians**
    (`dH/dtheta`, `dexp(eta)/dtheta`, `deta/dtheta`); exponential /
    log-logistic / log-normal fall back to `numDeriv::jacobian` on a
    per-call cumhaz closure.
  - **Transforms match SAS HAZARD** (`hzp_calc_haz_CL.c` /
    `hzp_calc_srv_CL.c`): `hazard` and `cumulative_hazard` use
    log-scale CLs; `survival` uses log(-log S) CLs (equivalent to
    log-cumhaz) so 0 <= lower <= upper <= 1; `linear_predictor` is
    symmetric on the natural scale.
  - **Fixed-shape / CoE multiphase fits produce meaningful CLs** — the
    delta-method sandwich is restricted to the free-parameter submatrix
    of `vcov`, treating fixed parameters as known-with-zero-variance.
  - Backward compatible: `se.fit = FALSE` (default) preserves the
    pre-0.9.8 scalar-vector / decompose-data-frame return shape.

# TemporalHazard 0.9.7

## New features

* **Counting-process / repeating-events likelihood wired up** — Phase 4f
  of the development plan lands. `Surv(start, stop, event)` with any
  `start > 0` is now accepted. The Weibull and multiphase log-likelihoods
  apply `H(stop) - H(start)` to event and right-censored terms; the
  trivial `start = 0` case degenerates to `H(stop)` and recovers the
  plain-Surv fit exactly. Splitting each row into contiguous epochs
  preserves both the log-likelihood and the MLE to optimizer tolerance
  (split-invariance).
* **Weibull + multiphase analytic gradients handle H(start).** The
  closed-form Weibull score adds a `-d H(start)/d theta` term per row
  (guarded at `start = 0`). The multiphase analytic gradient computes
  per-phase `Phi_j(start)` and its shape derivatives, then adds
  `+w_H_start * mu_j * dPhi_j(start)` to each parameter's score; G3
  phase derivatives at `start` use the same finite-difference machinery
  as at `stop`.
* **0.9.5 narrowing removed.** The `hazard()` guard that rejected
  counting-process `Surv(start, stop, event)` with any `start > 0` is
  gone.

# TemporalHazard 0.9.6

## New features

* **`weights` now supported for all distributions** — Phase 4e of the
  development plan lands. The exponential, log-logistic, and
  log-normal likelihoods and their analytic gradients now apply row
  weights to every censoring term (event, right-censored,
  left-censored, interval-censored). The 0.9.5 guard in `hazard()`
  that rejected `weights` for `dist %in% c("exponential",
  "loglogistic", "lognormal")` has been removed. Fits with integer
  weights reproduce the row-duplicated fit to optimizer tolerance
  across all five distributions.
* **Conservation of Events now honours weights.**
  `.hzr_conserve_events()` and `.hzr_select_fixmu_phase()` take an
  optional `weights` argument; the multiphase optimizer threads it
  through so per-phase cumulative hazards are summed on the same
  scale as the (weighted) observed event count. CoE no longer
  auto-disables when weights are non-uniform — the dimension
  reduction stays on and the MLE matches the full-dim path.

## Bug fixes

* **Multiphase analytic gradient now applies `weights`.**
  `.hzr_gradient_multiphase()` accepted neither `weights` nor its
  downstream equivalents: the per-row score weights `w_H` / `inv_h`
  were set to ±1 and the interval-censored finite-difference
  correction summed an unweighted LL. Weighted multiphase fits
  therefore optimised a weighted objective with an unweighted score;
  BFGS line search still converged near the correct MLE but the
  final gradient norm did not go to zero. All three paths now honour
  row weights, and the optimizer's `gradient_fn` wrapper (including
  the all-zero numeric fallback and the CoE wrapper) forwards
  `weights` consistently. Regression test covers weighted analytic
  vs numerical gradient parity. Surfaced by Copilot review on PR #18.

# TemporalHazard 0.9.5

## New features

* **Stepwise covariate selection** — `hzr_stepwise()` runs forward,
  backward, or two-way stepwise selection on an existing `hazard` fit
  using Wald p-values or AIC deltas as the entry / retention criterion.
  Phase-specific entry is supported for multiphase models: a covariate
  can enter one phase and not another. Defaults match SAS `PROC HAZARD`
  (`SLENTRY = 0.30`, `SLSTAY = 0.20`); AIC mode uses `ΔAIC < 0`
  uniformly. SAS-style `MOVE` oscillation guard freezes variables that
  enter + exit more than `max_move` times. Returns an object of class
  `c("hzr_stepwise", "hazard")` with a `$steps` selection trace, scope
  record, and elapsed timer. Implements the core algorithm from C
  HAZARD `stepw.c` / `backw.c`.

## Bug fixes

* **Multiphase convergence after weights/repeating-events merge** —
  restored multiphase optimization that regressed in 0.9.4: three
  interacting defects in the new `weights` threading (dup-arg
  collision in the multiphase / Weibull closures, positional-arg
  corruption in every distribution's gradient call) made every
  optimizer iteration error silently inside `tryCatch`. Diagnosed and
  fixed via commit 73b4657.
* **Weibull analytic gradient now applies `weights`** — both
  `.hzr_gradient_weibull()` and the `grad_internal` closure inside
  `.hzr_optim_weibull()` accepted `weights` as a formal but did not
  apply it to the score vector. The optimizer still converged via
  line search on the (weighted) log-likelihood, but the gradient
  direction was wrong and the final gradient norm did not go to
  zero. Both gradient paths now weight the event indicator and
  cumulative hazard building blocks. Fits with integer weights
  reproduce the equivalent row-duplicated fit to optimizer tolerance.

## Scope change

* `weights` is now only accepted for `dist = "weibull"` and
  `dist = "multiphase"`. The 0.9.4 NEWS claimed weights were
  threaded through all distribution-specific likelihoods; in fact the
  exponential, log-logistic, and log-normal single-distribution paths
  accepted the formal but never applied it, so the fit was silently
  unweighted. `hazard()` now raises an explicit error when `weights`
  is supplied with one of those distributions rather than returning
  an unweighted fit. Full support for the remaining single-dist paths
  is tracked in `inst/dev/DEVELOPMENT-PLAN.md` Phase 4e.
* **Conservation of Events is auto-disabled when weights are not
  all 1.** `.hzr_conserve_events()` receives the weighted event count
  as its target but sums per-phase cumulative hazards across rows
  *without* applying weights, so Turner's adjustment comes out on a
  mismatched scale. The multiphase optimizer now detects non-unit
  weights and skips the CoE dimension reduction, falling through to
  the (correctly weighted) full-dimensional path. Fits are still
  correct; they just don't benefit from the one-parameter
  analytical closed-form solve. Weighted CoE wire-up is tracked
  alongside the other weights completion work in
  `inst/dev/DEVELOPMENT-PLAN.md` Phase 4e.
* **Repeating-events / counting-process notation narrowed.**
  `Surv(start, stop, event)` with `start > 0` is no longer accepted
  by `hazard()`. The 0.9.4 NEWS claimed each epoch contributed
  `H(stop) - H(start)` to the likelihood, but downstream likelihoods
  only read `time_lower` for interval-censored rows (`status == 2`);
  counting-process rows (`status` in `{0, 1}`) were silently scored
  with `H(stop)` alone, so any fit with nonzero entry times was
  silently wrong. `hazard()` now raises an explicit error. The
  trivial case `Surv(0, t, d)` -- equivalent to `Surv(t, d)` --
  continues to work. Full wire-up of `H(stop) - H(start)` for all
  distribution paths is tracked in
  `inst/dev/DEVELOPMENT-PLAN.md` Phase 4f.

# TemporalHazard 0.9.4

## New features

* **Observation weights** — `weights` argument in `hazard()` applies Fisher
  weighting to the log-likelihood for `dist = "weibull"` and
  `dist = "multiphase"`. Each observation's contribution is multiplied
  by its weight, enabling severity-weighted event analyses. Implements
  the SAS `WEIGHT` statement. _The original 0.9.4 entry claimed
  coverage of all distribution paths; the 0.9.5 patch corrected the
  claim and fixed a gradient wire-up bug in the Weibull path._
* **Repeating events** — `Surv(start, stop, event)` start-stop notation
  is parsed. _The original 0.9.4 entry claimed each epoch contributed
  `H(stop) - H(start)` to the likelihood, but the downstream
  likelihoods never applied the lower bound for counting-process rows;
  the 0.9.5 patch narrowed the feature to the trivial `start = 0`
  case and added an explicit error for nonzero starts._

# TemporalHazard 0.9.3

## New features

* `hzr_deciles()` — Decile-of-risk calibration function comparing observed
  vs. expected event counts across risk groups with chi-square GOF testing.
  Implements the SAS `deciles.hazard.sas` macro workflow.
* `hzr_gof()` — Goodness-of-fit function comparing parametric predictions
  against nonparametric (Kaplan-Meier) estimates with observed vs. expected
  event counting. Implements the SAS `hazplot.sas` macro workflow.
* `hzr_kaplan()` — Kaplan-Meier survival estimator with logit-transformed
  confidence limits that respect the [0, 1] boundary, interval hazard rate,
  density, and restricted mean survival time (life integral). Implements the
  SAS `kaplan.sas` macro output structure.
* `hzr_calibrate()` — Variable calibration function for assessing functional
  form before model entry. Groups a continuous covariate into quantile bins
  and applies logit, Gompertz, or Cox link transforms. Supports
  stratification via the `by` parameter. Implements the SAS `logit.sas` and
  `logitgr.sas` macros.
* `hzr_nelson()` — Wayne Nelson cumulative hazard estimator with lognormal
  confidence limits. Supports weighted events for severity-adjusted repeated
  event analyses. Implements the SAS `nelsonl.sas` macro.
* `hzr_bootstrap()` — Bootstrap resampling for hazard model coefficients with
  bagging support (fractional sampling). Returns per-replicate estimates and
  summary statistics (mean, SD, percentile CI). Implements the SAS
  `bootstrap.hazard.sas` macro workflow.
* `hzr_competing_risks()` — Competing risks cumulative incidence using the
  Aalen-Johansen estimator with Greenwood variance. Handles any number of
  competing event types. Implements the SAS `markov.sas` macro.
* **Conservation of Events (CoE)** — Turner's theorem is now integrated into
  the multiphase optimizer. One phase's log_mu scaling parameter is solved
  analytically at each iteration, reducing the optimization dimension by 1
  and improving numerical stability and convergence. Enabled by default;
  disable with `control = list(conserve = FALSE)`. Implements the core
  algorithm from C HAZARD `setcoe.c` / `consrv.c`.
* New vignette: "Complete Clinical Analysis Walkthrough" — end-to-end
  workflow from Kaplan-Meier baseline through validated multivariable model,
  mirroring the SAS HAZARD analytical sequence.

## Improvements

* Multi-start optimizer now respects user-set RNG seeds for reproducibility
  (removed `set.seed(NULL)` that was actively breaking determinism).
* Vignette metadata normalized to YAML `vignette:` key across all 8 files.
* `fit` parameter documentation corrected to state default is FALSE.
* README now includes key capabilities table and development plan link.

# TemporalHazard 0.9.1

## New features

* G3 late-phase decomposition (`hzr_phase("g3", ...)`) now fully integrated
  into the multiphase optimizer, Hessian, and prediction pipeline.
* `fixed = "shapes"` parameter in `hzr_phase()` allows fixing shape parameters
  during estimation (matching C/SAS HAZARD workflow of estimating only log-mu
  scale parameters).

## Bug fixes

* `summary.hazard()` now correctly reports standard errors when some
  parameters are fixed. Previously, `anyNA(vcov)` rejected the entire
  variance-covariance matrix when fixed parameters had NA entries.
* `print.summary.hazard()` coefficient table now shows the correct label
  for G3 phases (was printing empty parentheses).
* `print.summary.hazard()` phase listing now uses the phase name in
  CDF labels (e.g., "cdf (late risk)") instead of hardcoded "early risk".
* SAS missing value markers (`.`) in CSV datasets are now handled via
  `na.strings = c("NA", ".")` in `data-raw/make_data.R`, preventing
  numeric columns from being read as character.

## Documentation

* Seven Quarto vignettes: getting-started, fitting-hazard-models,
  prediction-visualization, inference-diagnostics, mathematical-foundations,
  package-architecture, and sas-to-r-migration.
* Roxygen examples now include both single-phase and multiphase models.
* README switched to self-contained CABGKUL examples with G3 late phase.
* Dataset axis labels corrected to "Months" (not "Years").

## Infrastructure

* CI workflows updated to use `roxygen2::load_pkgload` for lazy data
  compatibility.
* Added lintr CI workflow with `.lintr` configuration.
* pkgdown action bumped to `peaceiris/actions-gh-pages@v4`.
* Added `use-public-rspm: true` to all CI workflows.
* Added `lintr` to Suggests.

# TemporalHazard 0.9.0

## New features

* Multiphase engine: N-phase additive cumulative hazard models via
 `dist = "multiphase"` with `hzr_phase()` specification.
* `hzr_decompos()` parametric family implementing the three-parameter
  temporal decomposition of Blackstone, Naftel, and Turner (1986).
* Multi-start optimizer with Hessian-based variance-covariance estimation.
* C binary parity tests against the KUL CABG reference dataset.
* Five clinical reference datasets: `avc`, `cabgkul`, `omc`, `tga`, `valves`.

# TemporalHazard 0.1.0

## New features

* Single-phase engine: Weibull, exponential, log-logistic, and log-normal
  distributions with formula interface.
* `hazard()` API with `predict()`, `summary()`, `coef()`, `vcov()` S3 methods.
* Golden fixture regression testing system.
* Numerically stable helper primitives (`hzr_log1pexp`, `hzr_log1mexp`,
  `hzr_clamp_prob`).

# TemporalHazard 0.0.0.9000

* Initial package scaffold.
* Added numerically stable helper primitives.
* Added baseline unit tests and CI workflow.
