# TemporalHazard 1.2.11

## Breaking changes

* **`hazard(fit = TRUE)` without `theta` is now an error for the
  single-distribution models.** For `dist = "weibull"`, `"exponential"`,
  `"loglogistic"` and `"lognormal"`, the optimizer ran only when `theta` was
  supplied, so a call that left it out returned an unfitted object -- `NULL`
  coefficients, an `NA` objective -- with no error and no warning. `print()`
  gave no sign of it. The call now stops and asks for starting values, as it
  does for a zero-length `theta`, which used to fail inside `optim()` with a
  message that did not name the cause.
  `dist = "multiphase"` is unaffected: it assembles its own start from
  `phases`. `fit = FALSE` without `theta` still builds an unfitted model.

  Code that relied on the old behaviour was getting no fit. Two tests in this
  package were: they compared `NULL` coefficients with `NULL` coefficients,
  so the counting-process equivalence and epoch-split invariance they claimed
  for the Weibull were never checked. Both now fit, and both pass.

* **`hzr_translate_sas()` now emits a `stop()` in place of the fit when a
  `PARMS` statement builds no phase it could use.** Operands the translator
  could not read (a template's `MUE=?`, or `MUE = 0.2` written with spaces
  around `=`, which `PROC HAZARD` accepts) or could not use (a `MUE` or `MUL`
  with no shape operand) are recorded in `$untranslated`, but the fit chunk
  used to be emitted anyway, as `hazard(fit = TRUE, theta = c())` under the
  default Weibull. That chunk rendered an unfitted object, and would now fail
  on the error above with a message about `theta` that does not name the real
  cause. The emitted `stop()` names it. This is a limit of the translation,
  not a `PROC HAZARD` refusal, so it is kept apart from the existing
  "selects no phase" stop.

* **A multiphase formula that names a phase as a function is now an error**
  (#275). `hazard(Surv(int_dead, dead) ~ constant(age), dist = "multiphase",
  phases = ...)` read `constant(age)` as a phase-scoped term, then replaced
  the whole right-hand side with `~ 1`, and nothing sent the term to its
  phase. The fit converged without an `age` coefficient, with no error and
  no warning. `hazard()` now stops, names the phase or phases it found, and points to
  `hzr_phase(..., formula = ~ var)`, which is where a phase's covariates
  belong. Code that relied on the old behaviour was fitting a model without
  those covariates; drop the terms from the formula to keep that model, or
  move them into `hzr_phase()` to get the one the formula described.
  Formulas that call an ordinary function such as `log()`, and plain global
  covariates, are unaffected.

* **`hzr_bootstrap()` now refuses a fit whose formula uses a per-row
  variable that is not a column of its `data`** (#278). Replicates resample
  the rows of `data`, so such a variable was held fixed while the rows moved
  under it. The interval was wrong, and every replicate still reported
  success, with no warning: on `avc`, a copy of `age` kept outside `data`
  gave an interval that excluded its own estimate. The check covers the
  response, the covariates of the global and phase formulas, and a
  select-mode `scope`, and reads variables from the terms, so `log(age)`
  needs only the column `age`. A constant outside `data`, such as `pi`, a
  cutoff or a knots vector, is still allowed. The error names the
  variables; add them to `data` and refit. A scope variable that only the
  scope formula's own frame can see, which the refits could never test, is
  refused the same way. So is a vector-interface fit whose design matrix
  was passed directly as `x`: it was re-evaluated without resampling in
  every replicate, or in select mode dropped from the candidate refits, so
  all of them succeeded and the interval was wrong.

* **`hzr_stepwise()` now refuses a base fit written as `Surv(...) ~ .`**
  (#279). It read the base model's terms without the data, which cannot
  expand `.`, and treated the failure as a model with no terms: the screen
  reported zero steps, which looks the same as finding nothing to drop. It
  now stops before printing anything and asks for the base model's terms to
  be written out. A `scope` of `~ .`, which used to give a screen with no
  candidates, stops the same way, and so does a multiphase base fit whose
  global formula uses `.` while a phase has no formula of its own and so
  inherits it. A screen whose base model has its terms written out is
  unchanged, as is a multiphase screen in which every phase has its own
  formula.

* **An entry time after the exit time is now an error, and
  `time_lower = time` now means "no entry" in every family** (#253). On a
  row with status 0 or 1, `time_lower` is the counting-process entry time
  when `0 < time_lower < time`. `hazard()` used to warn when
  `time_lower >= time` on such a row and fit anyway, and each family then
  did something different. Multiphase returned a "log-likelihood" of
  +47915.76 with `converged = TRUE` on the AVC data. The Weibull read the
  rows as entering at time 0, and the other three families ignored
  `time_lower` altogether.
  - `time_lower > time` on a status 0/1 row now stops, as SAS HAZARD
    rejects a start time after the exit time (error `SETCOE960`).
  - `time_lower == time` on a status 0/1 row is read as no entry time, with
    no warning, in all five families. This is the mixed-interval layout,
    where exact and right-censored rows carry `time_lower = time` and only
    interval-censored rows (status 2) carry a real lower bound. It was
    already the Weibull rule; multiphase used to degenerate on it.
  - `time_lower = 0` still means no entry time.
  - Rows with `time_lower == time > 0` beside rows with a genuine entry
    time now stop. In counting-process data they are zero-length epochs,
    which `hzr_repeated_events()` can emit; read as "no entry", each would
    be charged its full cumulative hazard from time 0.

## New features

* **Every fit now says what it did not do** (#242, following #197). A
  `hazard` object carries `degraded`, the steps the fit did not perform, and
  `degraded_causes`, the reason for each. `print()` and `summary()` always
  show them as a "Not done in this run" block, and the block reads "none"
  when nothing was lost: a line that appears only on bad news cannot be told
  from one that was never written. Five steps are recorded:
  - fitting itself: `fit = FALSE`, or a fit imported from SAS output;
  - standard errors, naming whether numDeriv was missing, `numDeriv::hessian()`
    failed, or the Hessian was non-finite or singular, and, when a covariance
    was computed, which estimated parameters were left without a usable
    variance;
  - the variance of the phase that Conservation of Events conserves;
  - the weak-direction check;
  - Conservation of Events itself.

  The block replaces two notes that could name the wrong cause: "not examined
  for a weakly identified direction" and "standard errors unavailable; the
  Hessian could not be inverted". `fit$fit$weak` keeps its meaning, and is
  `NA` exactly when `"weak_direction_check"` is listed. An object saved by an
  earlier version prints "not recorded" rather than "none".

* **`hzr_translate_sas()` translates a `%repeat` call** into
  `hzr_repeated_events()` (#241), renaming its outputs to the names the job
  gives them so the job's fit reads them. The macro's input is still the
  reader's to supply. Any step between the macro and the fit that names the
  macro's output, or uses a macro variable that might, stops the document
  with the step quoted, rather than fitting data the job changed; a plain
  `PROC SORT` is the one step let through. A step that changes the output
  without naming it, such as a macro that writes it internally, is not
  detected.

## Bug fixes

* **The multiphase gradient and Hessian are now right when an early phase's
  `m` is near 0.** Both differentiate in `m` by finite differences, and
  their stencils straddled 0: the gradient's (half-width about 6e-6)
  whenever `|m|` was smaller than that, the analytic Hessian's (half-width
  1.2e-4) whenever `nu > 0` and `|m|` was below 1.2e-4. `hzr_decompos()`
  changes formula at `m = 0`, and the `m < 0` family meets the `m >= 0` one
  in a cusp rather than continuing it, so each difference mixed two
  branches and returned neither side's derivative. The gradient gave +8.4
  where the true value was -20.2, on a 13-parameter fit that converged to
  `m = 2.9e-6`; the Hessian put -4.8e5 on the `m` diagonal where the value
  is about 0.7, so the standard errors of such fits were wrong too. Every
  other parameter's gradient was unaffected.

  Both stencils now keep the sign of `m`, one-sided on the `m >= 0` side,
  where the family is smooth: there the gradient and the analytic Hessian
  are now right. Below 0 the gradient's step is at most 1% of `|m|`, because
  the cusp varies on that scale, floored at 1e-10 so rounding stays bounded.
  The Hessian's steps below 0 are one-sided but still 1.2e-4 wide: they stop
  the branches mixing, but for `m` within about 1e-4 below 0 and `nu < 2`
  they understate the curvature, so the standard error of `m` there is too
  large (about five times, on the fit the tests use). At `nu = 0` that path
  used to stop with an error; it now returns these values. Fits whose
  standard errors come from `numDeriv` -- those with left- or
  interval-censored rows -- and the score test's information still
  difference across 0. Likelihood values are unchanged.

  The likelihood itself is still not differentiable at `m = 0`, so a fit
  whose optimum sits there reports a nonzero gradient. SAS/C never meets the
  point: it estimates `log|M|` with the sign fixed by the starting value, so
  `M` cannot reach or cross 0. `hazard()` estimates `m` directly and can.

* **A fit that reports convergence is now checked against SAS/C HAZARD's
  own test for it, and continued with `stats::nlm()` when it fails the
  test.** `hazard()`'s BFGS optimizer stops on the relative change in the
  log-likelihood (`control$reltol`, default 1e-5), which lets a flat ridge
  end short of the maximum with `converged = TRUE`: a 13-parameter
  early-CDF plus late-G3 model stopped 0.013 below the SAS listing's
  log-likelihood, and synthetic fits of the same shape up to 5 units below.
  SAS/C accepts an optimum only when the relative gradient,
  `max |g_i| * max(|x_i|, 1) / max(|f|, 1)`, is at most `eps^(1/3)`, about
  6e-6. When BFGS reports convergence and that test fails, every
  distribution's fit is now continued with `stats::nlm()`, the
  Dennis-Schnabel algorithm SAS/C's optimizer was ported from, at SAS's
  tolerances, and the continued point is kept only if the log-likelihood
  improves. The default `reltol` is unchanged: tightening it instead cost
  30% to 60% more time on the test suite and broke eight or nine tests.

  Every fit records the test in `fit$fit$rel_gradient` (`NA` when the test
  was not applied, because the optimizer did not report convergence, or the
  gradient cannot be evaluated) and, when the continuation improved the fit,
  `nlm()`'s termination code in `fit$fit$polish_code`, and `print()` and
  `summary()` show it. Under Conservation of Events the analytic score omits
  how the conserved scale moves with the other parameters, so there the test
  is computed from finite differences of the log-likelihood, as SAS/C does.
  The continuation keeps the analytic score, so a CoE fit can honestly end
  with the test not met. Only SAS/C's two hard failures warn: code 4, the
  iteration limit, and code 5, where the likelihood kept rising along some
  direction and may have no maximum. Codes 2 and 3, where SAS/C prints a
  caution and retries, are recorded without a warning. The test is relative
  to the size of the log-likelihood, so a fit that meets it is within SAS's
  tolerance of the maximum rather than exactly at it.

  Estimates of fits that used to stop short now change. One test depended
  on a detail of where BFGS stopped: it showed `gamma` and `eta`
  non-identified at `alpha = 1` by a large standard error. On that exactly
  flat ridge the Hessian is singular in theory, so whether a finite standard
  error comes out at all is numerical noise, and at the polished point it
  does not. The test now accepts either a missing or a 100-fold larger
  standard error, and also checks that both fits reach the same
  log-likelihood and the same `gamma * eta`.

* **Exponential, log-logistic and log-normal fits ignored left truncation**
  (#253). On status 0/1 rows these three families used `time_lower` only as
  a censoring bound, which applies to status 2, so a left-truncated fit was
  silently fitted as if every subject had been at risk from time 0. They
  now subtract the cumulative hazard at entry, H(time) - H(time_lower), as
  the Weibull and multiphase likelihoods already did. The log-likelihood,
  its gradient and the closed-form Hessian all carry the entry term.

* **`hzr_gof()` reported a conservation ratio that was not one** (#254).
  For a model with covariates it computed expected events from a single
  curve at the covariate means, then printed the total as the
  "Conservation ratio (E/O)". On the covariate model in the clinical
  walkthrough vignette that printed 0.606, while the fit conserved events exactly (68.000 expected
  against 68 observed). Expected events are now summed per subject, each
  subject's cumulative hazard at exit minus that at entry. For Weibull,
  exponential and multiphase fits with conservation of events, E/O is then
  the conservation-of-events identity; for log-logistic and log-normal fits
  it checks calibration in total. For a weighted fit, both observed and
  expected events now carry the case weights, since that is what a weighted
  fit conserves (the sum of w·H equals the sum of w·d). The `par_surv` and
  `par_cumhaz` columns are still the covariate-mean curve, for plotting
  against Kaplan-Meier, and the risk-set counts and Kaplan-Meier columns stay
  unweighted. Unweighted intercept-only fits without entry times are
  unchanged.

* **`hzr_gof()` drew the mean-patient curve at covariates of 0** for a
  multiphase fit whose covariates enter only through the phase formulas.
  Such a fit has no global design matrix, so the `par_surv` and
  `par_cumhaz` columns were predicted from time alone, which set every
  phase covariate to 0. They now use each phase's design-matrix column
  means, so a factor enters as the proportion of patients in each level.

* **`hzr_gof()` stopped on a multiphase fit with both a global covariate and
  phase-formula covariates, and every multiphase fit carried a
  `par_cumhaz_time` column** (#263, #264). With `Surv(...) ~ age` and a phase
  formula `~ mal`, the mean-patient curve was built from the global
  covariates alone, so `predict()` found no `mal` column and stopped before
  the expected-event tally. The curve is now evaluated at the column means
  of each phase's own design matrix, for global, phase-formula and mixed
  covariates alike. With `time_windows`, a multiphase fit's output had twice
  as many rows as grid times. The mean patient now carries the covariate
  means in the window that contains each time. Separately, the
  `par_cumhaz_<phase>` columns were chosen
  by dropping `total` from the decomposition, which let its `time` column
  through as a phase; they are now chosen by phase name.

* **`hzr_gof()` had five smaller errors**, found by Copilot's and
  r-reviewer's reviews of #285.
  - `seq()` builds grid times that differ from the data times in the last
    binary digits, and they were matched with a fixed tolerance of 100
    machine epsilons. Above 128 that is smaller than the gap between
    adjacent doubles, so at times in days or months events fell off the
    grid: 180 of 197 were counted on a `seq(150, 300, by = 0.1)` grid. The
    tolerance now scales with the time.
  - With a custom `time_grid`, `n_risk` between Kaplan-Meier times carried
    the previous count forward, so it kept subjects who had left and
    missed ones who had entered. It now counts the risk set at each grid
    time.
  - An unsorted grid made the cumulative columns non-cumulative. The grid
    is now sorted, with repeated times dropped.
  - With an event at time 0, `km_surv` at 0 was averaged with 1.
  - For a fit with both `time_windows` and entry times, expected events
    took H(entry) in the entry-time covariate window, while the likelihood
    uses the exit-time window. E/O came out 1.052 on a Weibull fit that
    conserves events.

* **A Weibull fit with one masked variance reported the others on the wrong
  scale.** When the Hessian inverse has a non-positive variance, its row and
  column are set to `NA`. The delta-method transform from the internal
  `(alpha, psi)` scale to the reported `(mu, nu)` scale was then skipped for
  the whole matrix, so the surviving standard errors stayed on the internal
  scale beside `(mu, nu)` estimates, with no error. A masked `alpha` printed
  `SE(psi)` as the standard error of `nu`, off by a factor of `nu`; a masked
  `psi` printed `SE(alpha)` as the standard error of `mu`, which depends on
  `psi` and has no valid standard error there. The transform now runs on the
  finite block and carries the mask through the Jacobian: a reported
  parameter is `NA` exactly when it depends on a masked one. Exponential,
  log-logistic and log-normal fits are unaffected; they report on the scale
  they are optimised on.

* **A `survival::Surv()` object passed as `status` is now translated, as the
  formula interface always did** (#226). `Surv()` codes censoring with
  different integers from this package, and the vector interface took the
  object's second column unchanged. Under `type = "left"` a left-censored row
  was fitted as right-censored; under `"interval"` and `"counting"` the
  second column is not the status at all, so the fit read `time2` or `stop`
  as status codes. There was no error and no warning, and
  `objective = "sas"` could not see a left-censored row to refuse it. Both
  interfaces now read the `Surv` through one internal helper, driven by its
  `type`, so they store the same status and bounds and give the same fit.
  The bounds a `Surv` carries are taken from it; a `time`, `time_lower` or
  `time_upper` that disagrees with them is an error rather than being
  silently replaced. `hzr_bootstrap()` resamples those bounds too, although
  they never appear in the stored call.

* **`.` in a `hazard()` formula no longer puts the response in the design**
  (#273). `Surv(int_dead, dead) ~ .` expanded `.` to every column of `data`,
  including `int_dead` and `dead`, so the outcome was fitted as a predictor.
  With starting values sized for those extra columns the fit converged, with
  no error and a log-likelihood far above the correct model's. With starting
  values sized for the real covariates it stopped with "non-conformable
  arguments", which did not
  name the cause, and `predict(newdata = )` demanded the response columns.
  `.` now means every column the `Surv()` term does not use, as in
  `survival::coxph()`, so a `~ .` fit gives the same design and estimates as
  the formula written out in full. A `data` with no other column gives a
  model with no covariates, and there `.` beside other terms is an error.
  **Estimates from an earlier `~ .` fit change**, and so does the length of
  `theta` it needs. A `.` in `hzr_phase(formula = )` is fixed separately
  (#277). A right-hand-side variable that is not a column of `data` is now
  looked up where the formula was written, so a
  variable local to the calling function resolves instead of failing with
  "object not found". `hzr_bootstrap()`, which resamples only the rows
  of `data`, refuses such a fit (#278).

* **`hzr_phase(formula = ~ .)` no longer puts the response in the phase
  design** (#277). A phase formula's `.` was expanded by `model.frame()`
  against every column of `data`, including the columns of the `Surv()`
  term, so the outcome was fitted as a phase covariate and the fit
  converged with no error. `hazard()` now writes `.` out once, before
  fitting, the same way as for the global formula (#273): every column the
  `Surv()` term does not use. The fitted object stores the written-out
  formula, so `predict(newdata = )` no longer needs the response columns.
  On the vector interface (`time =`, `status =`) no `Surv()` term says which
  columns hold the response, so a phase formula with `.` is now an error
  there; write the phase's terms out. **Estimates from an earlier fit with
  `.` in a phase formula change.**

* **`hzr_argument_mapping()` listed DELTA as implemented.** Its
  `implementation_status` was `"implemented"` and its `r_parameter` read
  "(absorbed by decompos)", while the row's own notes say a non-zero DELTA
  is refused or flagged and never fitted (#181). The row is now
  `"planned"` with `r_parameter` "(not implemented)", so
  `hzr_argument_mapping(include_planned = FALSE)` no longer includes it.

# TemporalHazard 1.2.10

## New features

* New `hzr_repeated_events()` rebuilds the input to a repeated-events hazard
  model. From a long data set with one row per candidate event per subject, it
  returns one row per inter-event segment, with the segment's start time,
  duration and running event count. It reproduces the SAS macro `%repeat`,
  which built this input for the repeated-events `HAZARD` jobs and whose output
  was seldom saved, so those jobs can now be run again in R. It refuses input
  that would otherwise give a plausible but wrong result -- a non-numeric time,
  follow-up or indicator column, a factor `id`, a missing `followup` value, or
  an empty data frame -- and warns, naming the subjects, when an event falls
  after the end of follow-up, when `followup` varies within a subject, or when
  a missing time leaves a segment undefined. As in the macro, `rcensor` and
  `event` can both be 1 on the same row; see `?hzr_repeated_events`.

## Bug fixes

* **`hzr_translate_sas()` now starts an unspecified shape parameter where
  `PROC HAZARD` starts it, not where `hzr_phase()` does.** A `PARMS` statement
  that named only some of a phase's shape operands had the rest filled from
  `hzr_phase()`'s defaults, and three of them disagree with the reference C
  (`src/hazard/stmtprc.c`): `NU` starts at 2 rather than 1, `M` at 1 rather
  than 0, and `ETA` at 2 rather than 1. Since the multiphase likelihood is
  multimodal, a different starting vector can reach a different optimum, so
  this was a fidelity divergence rather than a cosmetic one. The defaults are
  now `PROC HAZARD`'s, and the emitted `hzr_phase()` call names every shape
  argument explicitly so that what is printed and what reaches the optimizer
  cannot disagree. `TAU` is the exception, because `PROC HAZARD` derives it
  from the data: an unspecified `TAU` becomes `0.75 * Tmax`
  (`src/hazard/readobs.c`) and one written as non-positive becomes
  `2 * Tmax / 3` (`SETG3()`). Neither can be reproduced at parse time, so such
  a phase is emitted at `tau = 1` and recorded in `$untranslated`, naming
  whichever rule applies -- unless `SETG3_ignore_tau()` does, which pins `TAU`
  at 1 anyway. No job in the *public corpus* is partially specified, so no corpus
  translation changes; the package's own end-to-end fits test does carry a
  partial block (`MUE THALF NU MUC`, no `M`), and its early phase now starts at
  `m = 1`.

* **`hzr_translate_sas()` now pins `TAU` where `SETG3_ignore_tau()` pins it.**
  When `ALPHA` is fixed at 1, `setg3.c:378-379` sets `TAU` to 1 *and* fixes
  it; the emitted `hzr_phase()` call mirrored neither, leaving `TAU` free. At
  `alpha = 1` the `G3` form collapses to `(t/tau)^(gamma*eta)`, so `log_mu`
  and `log_tau` are exactly aliased: the translated fit converged onto a flat
  ridge and returned no standard errors for either, with nothing in
  `$untranslated` to say so. The emitted call now carries `tau = 1` and
  `"tau"` in `fixed`, which identifies `log_mu` again. Jobs whose `PARMS`
  named a different `TAU` additionally record a row, since that value is used
  by neither `PROC HAZARD` nor the translation. No corpus translation changes:
  every late-phase block in the public corpus already writes `TAU=1 FIXTAU`,
  so the emitted calls are byte-identical before and after (checked by running
  the parser over all of them). What changes there is that those jobs no
  longer need a warning.

  The same branch's `ETA` fix (`setg3.c:405`) is mirrored too: with `GAMMA` and
  `ETA` both free at `alpha = 1` the pair is exactly singular, and
  `hzr_phase()` previously left both free and fitted the ridge. On the
  package's own fixture that moves `gamma`'s standard error from 7.5 to 0.02.
  This replaces an `$untranslated` row with a faithful translation.

* **`hzr_translate_sas()` now reports the whole of what `SETG3()` would do to
  a late phase, and mirrors only the part that has to be mirrored.**
  `PROC HAZARD` does not optimize from the operands `PARMS` supplies: `SETG3()`
  rewrites them first, and refuses some jobs outright. The translator now walks
  that function (`src/model/setg3.c`) and records both, splitting them on a
  single principle:

  - A rewrite that resolves an **exact non-identifiability** is *mirrored*,
    because the degeneracy is algebra and is just as real in R. There are two,
    both in `SETG3_ignore_tau()`: the `TAU` pin and the `ETA` fix, described
    above.
  - Every other rewrite keeps `PROC HAZARD` inside a numerical branch it can
    evaluate -- the role `g3flag` plays, which `hzr_decompos_g3()` does not
    need because it carries the general four-parameter `G3` form. Copying those
    would import a SAS limitation into R, and reaching a late shape the
    reference implementation cannot is a purpose of this package. They are
    *recorded* in `$untranslated`, so a SAS parity run knows why the starting
    values differ.

  In practice this covers `SETG3_verify_ge_2()`'s push of `GAMMA * ETA` clear
  of 2 (which the `PROC HAZARD` defaults `gamma = 1`, `eta = 2` trip exactly,
  so it applied to every defaulted non-`WEIBULL` late phase), the value
  substitutions in all eight sign branches, and `SETG3_alpha_fixup()` /
  `SETG3_alpha_gener()` deriving `ALPHA` from `GAMMA * ETA`.

  Nine refusal codes can now be recorded rather than emitted as runnable fits,
  where previously only `SETG3980` was -- and that one only for `alpha = 0`,
  not for the negative `ALPHA` that raises it too. (All sixteen are mirrored
  from the C, but seven guard conditions the entry checks at `setg3.c:269-284`
  have already refused, so no input reaches them; an exhaustive search in the
  tests pins which nine are live.)

  One correction to a rule this package had recorded wrongly: an **unspecified**
  `TAU` does not reach `SETG3()` as the value `stmtprc.c` starts it at, 0.
  `src/hazard/readobs.c:153-154` replaces it with `0.75 * Tmax` on an active
  late phase, and `readobs()` runs before `SETG3()` (`hazard.c:276` against
  `:292`). So an absent `TAU` starts at `0.75 * Tmax`, not `2 * Tmax / 3` --
  that rule (`setg3.c:317`) governs only a `TAU` the job wrote as non-positive
  -- and a bare `FIXTAU` cannot raise `SETG3900`. `ALPHA = 0` **with**
  `FIXALPHA` is the limiting exponential in both implementations and still
  translates cleanly; left free, `SETG3` derives an `ALPHA` instead, which is
  now reported.

  No corpus job is affected: every late-phase block in the public corpus
  carries `WEIBULL`, which returns at `setg3.c:347` before the sign dispatch,
  and all of them write `TAU=1 FIXTAU ALPHA=1 FIXALPHA` with a positive `GAMMA`
  and `ETA`. Verified by running the parser over all of them before and after.
* **`hzr_translate_sas()` no longer builds a phase that `PROC HAZARD` would
  not.** A `PARMS` statement names its phases with `MUE`, `MUC` and `MUL`; the
  shape operands (`THALF`/`NU`/`M` early, `TAU`/`GAMMA`/`ALPHA`/`ETA` late)
  only shape a phase that already exists. The translator had this the other way
  round and keyed on the shape operands, so
  `PARMS MUE=0.2 THALF=1 NU=1 TAU=2 GAMMA=1.5` emitted a two-phase model
  against `PROC HAZARD`'s one -- carrying an invented `mu` starting value of
  0.1 for the phase that should not have been there, and offering nothing in
  `$untranslated` to say so. In the reference C only `setparmno()` sets
  `C->phase[n]`, and only when that `MU` is greater than zero
  (`src/hazard/setparmno.c`); the seven shape operands are registered by
  `setprmf()`, which never touches it. A phase is now built only when its own
  `MU` was specified and positive -- so a `MUC` of exactly zero no longer
  builds a constant phase either, where before the translator tested only
  whether the keyword was present. Shape operands belonging to a phase that
  never activated are recorded in `$untranslated` rather than dropped, since
  `PROC HAZARD` zeroes them (`src/hazard/stmtprc.c`) and skips their covariates
  (`src/hazard/setstat.c`). A job that activates no phase at all is now
  **refused** rather than translated -- whether its `PARMS` named no positive
  `MU`, or it carried no `PARMS` statement whatsoever. `PROC HAZARD` does not
  run such a job: `src/hazard/modterm.c` raises `ERROR 1001: No phase
  selected` and the procedure exits before computing any results, so there is
  no fit for a translation to be faithful to. The translator previously
  emitted a runnable single-distribution fit and reported full token coverage;
  it now emits a `stop()` in place of the `hazard()` call, as it already did
  for `LCENSOR` combined with `ICENSOR`. The two cases are one state rather
  than two: `src/hazard/stmtprc.c` zeroes all three phases at initialization
  and only `setparmno()` turns one back on, so a job with no `PARMS` has no
  active phase for the same reason a `PARMS` naming `MUE=0` does. `modterm()`
  is reached on every job, not only once a multiphase model has been selected
  -- its one call site in `outmods()` is unconditional in the procedure's main
  sequence.

  The refusal fires only when every `PARMS` operand was understood. A
  statement this parser could not read is recorded, operand by operand, but
  never refused: `PARMS MUE = 0.2 THALF = 1` (spaces around `=`) parses to
  nothing here while `PROC HAZARD`'s own lexer discards whitespace and runs
  the job with an active early phase, so refusing it would stop a job the
  reference accepts. A second `PARMS` statement also now adds to the first
  rather than replacing it, matching the single field table `parmprc()` reads
  once after all statements are processed.

  One job in the public `hazard` corpus changes, and only to drop a claim that
  was wrong: the second `%HAZARD` block of
  `dist/examples/hm.dthar.TGA.sas` is a documentation template carrying
  literal `?` placeholders, which the reference would reject as a syntax
  error rather than as `ERROR 1001`. (That file's first block is a valid
  `PARMS` activating two phases and is unaffected.) Its per-operand rows are
  unchanged. No corpus job is refused.

* **`hzr_translate_sas()` no longer discards the `ALPHA` and `ETA` a `PARMS`
  statement specified alongside `WEIBULL`.** The translator read the bare
  `WEIBULL` keyword as a request to constrain the late phase to
  `alpha = eta = 1` and pinned both. `SETG3_weibull()` in the reference C
  (`src/model/setg3.c`) is the *generalized* Weibull -- it admits all positive
  parameter values, validates `gamma > 0`, `eta > 0` and `alpha >= 0`, and
  assigns nothing. A production job carrying
  `alpha = 2.501719 ... eta = 0.1365255 weibull` therefore translated to a
  seven-parameter fit where `PROC HAZARD` estimated nine, with both starting
  values replaced by 1 and no warning. `WEIBULL` now leaves `ALPHA` and `ETA`
  as `PARMS` gave them, free unless an explicit `FIXALPHA`/`FIXETA` pins them.
  The `G3`-collapses-to-Weibull identity at `alpha = eta = 1` is real and
  unchanged; it was simply not what the keyword requests.

# TemporalHazard 1.2.9

## Bug fixes

* **The phase-identifiability warning no longer names a cause that did not
  occur.** `variation` is the relative range of a phase's contribution across
  the observed times, so it collapses for every phase when those times carry
  nothing that separates them. The saturated scan fired anyway: with
  `time = rep(2, 20)` it told a `constant` phase that its shape parameters were
  unidentified, when a `constant` phase has none, and blamed a short half-life.
  A single-row fit produced the same text.

  That case is now reported in its own words, under a condition that takes two
  things together: the observed times' own relative range is below the
  tolerance, *and* no phase's contribution varies across them. Neither half is
  sufficient. The range alone would condemn times bunched far from the origin,
  where a steep phase can be fully identified; the phase measures alone would
  condemn well-spread times whenever a single phase saturates. Because the
  range is relative, near-ties reach the condition as exact ties do --
  `c(100, 100.000001)` produced the old wording verbatim. Per-phase saturation
  and absence over well-spread times are reported exactly as before, and the
  condition is now caught even when every phase carries covariates, where the
  per-phase measures are withheld and previously nothing was said at all.

  Two related corrections. A phase that has not started is reported as absent
  even in a left-truncated fit, since that test rests on the share and not on
  how the times are spread. And when the counting-process entry times or
  interval bounds supply enough evaluation points, the degenerate-times verdict
  is withheld -- and "enough" is counted rather than assumed. With the event
  times tied, each added point supplies one more functional of the parameters,
  so one more than the number of added points must reach the number of **free**
  parameters. The package's default two-phase model has five, and pinning a
  phase's shapes lowers the bar because those parameters are no longer free.
  Points that add nothing do not count: a zero entry time, a bound equal to an
  event time -- both of which `Surv(type = "interval")` synthesises for every
  row -- or a bound the objective never reads, since `time_upper` is live only
  for left-censored and interval rows. Counting any of those silenced this
  warning on data that needed it while changing no likelihood value.

  The message also no longer names shape parameters for a model that has none.
  A single `constant` phase is now silent, since one functional determines its
  `mu` exactly; two `constant` phases still warn, since only their sum is
  determined, but the wording says the phases cannot be told apart rather than
  naming shapes they do not have.

  Where the times are NOT degenerate, the per-phase saturation message is
  reported regardless of the bounds, so adding a `start` column no longer
  removes a correct diagnostic (#211).

  Known limitation, unchanged by this work and stated rather than implicit: the
  per-phase measures are taken over the event times alone, and the count of
  evaluation points is deliberately conservative -- it does not credit the
  extra functional a zero entry time makes separately observable, so a fit can
  be warned about while being identified. That cuts both ways.
  A phase can be flat across the event times while its shape is identified
  through the bounds, so the saturated message can overstate what is lost; and
  the counting rule above is a *local* identification argument for tied event
  times, so a fit can clear it and still have a flat direction the guard does
  not see, in which case nothing is said at all. Tracked in #228.

* **The SAS nomogram parser no longer deletes a whole-year time column.** The
  `PROC PRINT` observation counter is identified as a leading column running
  exactly `1..n`, which rested on no measurement column ever being a gapless
  1-based integer run -- true of the corpus in hand, not of SAS listings. A
  counter-suppressed nomogram on a whole-year grid prints `YEARS` as `1.0000,
  2.0000, 3.0000`, and the parser deleted the time key. Because the remaining
  columns then matched the header count, the shape check that would have caught
  it passed: the listing parsed, and the comparison ran against a table with no
  time column.

  The header now decides: a leading column the header names as a counter is
  dropped whatever values it carries, which also fixes the reverse case -- page
  two of a paginated listing starts at 41, and requiring a `1..n` run left that
  counter in the data. The labels are matched case-insensitively, after the
  leading underscore the parser strips, so SAS's `_N_` is recognised. A `1..n`
  run under a name the parser does not know is now kept and warned about rather
  than silently dropped: the listing cannot distinguish a counter spelled a new
  way from a measurement that happens to run `1..n`, and an extra column is
  visible where a deleted one is not (#212).

* **A phase named `total` is now rejected rather than silently losing its
  contribution.** `total` is the key the multiphase cumulative-hazard
  accumulator is stored under, so `phases = list(total = hzr_phase("cdf"), ...)`
  overwrote that phase's own contribution vector with the summed hazard. The fit
  returned normally; the phase then reported a contribution share of exactly 1
  and the identifiability diagnostic read it as dominating the model. The name is
  reserved across the decomposition and prediction output as well, where it
  labels the total column. `hazard()` and `hzr_theta_names()` -- which applies
  the same validation -- now stop with an explanatory message, and the
  reservation is documented on the `phases` argument (#214).

* **`objective = "sas"` data defects are reported as data defects.** The two
  conditions the SAS objective imposes -- no left-censored rows, and a positive
  width on every interval-censored row -- are pure functions of the data, but
  were checked inside the objective. The optimizer's per-start handler caught
  them and reported "produced no usable fit from N starts", which reads as a
  convergence problem and invites raising `n_starts`: a remedy that cannot work,
  because the condition is identical at every start. `hazard()` now checks both
  before any optimization, and the interval-width message reports the offending
  row in the data rather than its position among the interval rows. The guards
  inside the objective and gradient are unchanged, since the gradient is
  reachable without `hazard()` (#213).

* **Behaviour change:** because those preconditions are properties of the data
  and not of the fit, they are now checked whenever `objective = "sas"` is
  specified, including `fit = FALSE`. `hazard(..., fit = FALSE, objective =
  "sas")` on left-censored or zero-width-interval data previously returned a
  `hazard` object carrying an objective it could never be fitted with; it now
  stops. An `NA` in `status` under `objective = "sas"` is also reported by
  argument and row, rather than as a bare `missing value where TRUE/FALSE
  needed`.

  Note this reaches only the codes `hazard()` is given. On the **vector**
  interface a `survival::Surv()` object is unclassed without translating its
  codes, so they mean something else on arrival, and what goes wrong depends on
  `type`: under `type = "left"` a left-censored row is coded `0` and is fitted
  as *right-censored*; under `type = "interval"` it is coded `2` and is read as
  an *interval* of zero width, which this guard then rejects; and a genuine
  interval row, coded `3`, matches no branch of the likelihood and contributes
  nothing. The formula interface translates and is guarded. That asymmetry is a
  pre-existing defect of the vector path, tracked in #226.

## Testing

* **The SAS parity tests for `hz.te123.OMC` fit 1 and `hz.tm123.OMC` still
  described the P1 #6 Conservation-of-Events gap that PR #65 closed, and their
  tolerances were set from that gap rather than from the code's behaviour.**
  Fit 1 asserted the log-likelihood within `0.2` where R and SAS now agree to
  1.3e-04, `hz.tm123.OMC` within `0.5` against 4.6e-04, and neither asserted `MUE`
  at all, each carrying a comment quoting a value the fix had already moved. Restoring
  the pre-PR-#65 defect -- dropping the entry-time term from
  `.hzr_conserve_events()` -- left every one of those assertions passing, so
  they could not have caught the regression they were nominally about. The
  log-likelihood tolerances are now `1e-5` (relative), both
  Conservation-of-Events intercepts `MUE` and `MUL` are asserted, and each fit
  first checks `conserve_applied`, since CoE disables itself silently on
  unsupported data. The same mutation now fails six assertions.

  `MUL` is compared as a ratio to the SAS value rather than directly, and that
  distinction is the point rather than a detail. `expect_equal()` divides by
  `mean(abs(expected))` only when that exceeds `tolerance`; `MUL` is 2.1e-04,
  below any tolerance worth setting for it, so a direct comparison silently
  becomes an *absolute* one -- `MUL = 0` passes at `5e-04`, and so does the
  regression above. A first version of this change asserted `MUL` directly and
  reproduced, in the fix, the defect it was removing. Against an expected value
  of `1` the comparison is relative, as the tolerance implies. `hz.te123.OMC`
  fit 2's pre-existing `MUL` assertion sat 7% above that branch point and moves
  to the same form. `hz.te123.OMC` fit 2's log-likelihood tolerance
  goes from `1e-2` to `1e-5` for the same reason; it carried no stale claim, but
  `1e-2` admitted a drift of 3.1 in log-likelihood on a quantity matching to
  1.5e-04. No package code changed.

# TemporalHazard 1.2.8

## New features

* **New `hzr_theta_names()`** returns the names of a multiphase `theta`
  vector, in the order `theta` requires, before any fit runs. Use it to check
  a hand-written starting vector against the specification it belongs to:

  ```r
  phases <- list(early = hzr_phase("cdf"), late = hzr_phase("g3"))
  stopifnot(length(theta0) == length(hzr_theta_names(phases)))
  setNames(theta0, hzr_theta_names(phases))
  ```

  `theta` is positional and its entries are not on a common scale -- the late
  phase logs `mu` and `tau` but carries `gamma`, `alpha` and `eta` naturally --
  and **wrapping the wrong element in `log()` produces a fit, not an error**.
  A comment describing the order is therefore not enough for a template that
  ships to many studies, which is what prompted this: the order is a property
  of the phase specification and changes the moment an author adds, removes or
  retypes a phase.

  The function is a thin wrapper over the naming the optimizer itself uses,
  not a second implementation of it. `hazard()`'s optimizer, the score test's
  re-expansion and this function now all go through one internal helper, so
  the order documented here cannot drift from the order a fit produces --- a
  test asserts the two are identical for three different phase
  specifications. Phase validation is shared too, so unnamed phases get the
  same `phase_1` / `phase_2` labels a fit will give them.

* **`hzr_stepwise()`'s `$steps` frame gains a `stat_type` column**, saying
  what the `stat` on that row is and so which reference distribution
  recomputes its p-value: `"score_q"` (chi-square on `df`), `"wald_z"`
  (standard normal) or `"wald_chisq"` (chi-square on `df`).

  `df` could not tell these apart. A scalar Wald is reported as a *z*, not as
  its square, so it and a score Q are both recorded at `df = 1` while calling
  for different distributions. The gap was widest on a candidate rescued by
  the Wald fallback added in 1.2.7: that row carries a Wald z under
  `criterion = "score"`, where every neighbouring row carries a Q. The
  selection was right and `p_value` was right, but a reader recomputing a
  p-value from `stat` the way the neighbouring rows permit got an answer wrong
  by dozens of orders of magnitude. Under `criterion = "score"` the *entry*
  rows reading `"wald_z"` are the ones the fallback rescued --- filter on
  `action == "enter" & stat_type == "wald_z"`. Drop rows are always
  Wald-tested under that criterion, following SAS, so they read `"wald_z"`
  whether or not the fallback fired.

* **`hzr_bootstrap()` reports Wald fallbacks in select mode**, through
  `$n_wald_fallback_replicates` and `$n_wald_fallbacks`, with a warning when
  either is non-zero. Every replicate runs under `suppressWarnings()`, so a
  run in which the fallback fired throughout previously reported nothing at
  all -- and the bootstrap is where a wholesale substitution matters most,
  since those entries drive the pooled selection frequencies.

* **The nomogram parser accepts `lines =` as well as `path =`.** A multi-fit
  listing needs each nomogram attributed to the fit whose block contains it,
  not to the file. Previously a caller who had already split the listing by
  fit had to write each block back out to a temporary file to parse it, which
  is a workaround rather than an interface. Passing the lines directly now
  works, and is what makes the multiple-nomogram warning actionable.

## Bug fixes

* **`objective` is now recorded on the fit, so refits keep the estimand.**
  `hazard(objective = "sas")` was accepted and acted on by the optimizer, but
  the choice was never stored on the object. Everything that rebuilds a
  `hazard()` call from `fit$spec` therefore reverted to
  `objective = "likelihood"`: every `hzr_stepwise()` candidate refit, the
  Wald-fallback refit and each accepted move, plus the score test's numeric
  Hessian and gradient.

  Selecting on a `"sas"` base fit consequently differenced `delta_logLik`,
  `aic` and `delta_aic` across *two different estimands*, and the score
  statistic was computed against a likelihood the fit had not been fitted to.
  On the esophagectomy reference the two objectives differ by about 22
  log-likelihood units -- larger than most single-variable effects -- and the
  run produced a full `$steps` table with no error and no warning.

  `fit$spec$objective` now records it, and one accessor feeds every consumer
  so they cannot drift apart. A fit that carries no `objective` predates the
  argument and is read as `"likelihood"`, which is what it was.

  Note the asymmetry this had created: `hzr_bootstrap()` *refit* mode
  re-evaluates the stored call, which did carry `objective = "sas"`, so it was
  already consistent; *select* mode goes through the scope refit and was not.

  Following from this, **`objective` cannot be overridden per refit.** Passing it
  through `hzr_stepwise()`'s `...` is accepted when it restates the base fit's
  own objective and refused when it conflicts, with a message saying why:
  a candidate refit under a different objective would make `delta_logLik`,
  `aic` and `delta_aic` differences between two estimands rather than between
  two models. It is refused rather than quietly ignored, because a full
  `$steps` table that silently disregarded an explicit argument is the failure
  this entry is about.

* **A candidate that neither criterion could test is now reported as such.**
  `criterion = "score"` declines a candidate whose observed information is
  indefinite at `beta = 0` -- which happens when the effect is *large* -- and
  then refits it and tests it by Wald instead. That rescue can converge,
  producing a perfectly good point estimate, while its Hessian is singular:
  no standard error, so the Wald test cannot be computed either.

  The fallback dropped that case silently. The row kept the *score's* reason,
  `information_indefinite`, which describes the first of two independent
  failures and says nothing about the second, and no counter moved. A
  strongly predictive variable vanished and the screen rendered as an honest
  "nothing met `slentry`" -- the same failure this package has twice fixed
  elsewhere, where a clean-looking screen and an honest null result cannot be
  told apart.

  Such rows now report `fallback_no_variance`, distinct from
  `information_indefinite` (the rescue errored or did not converge, and is
  listed in `$criteria$refit_failures`). `hzr_stepwise()` and
  `hzr_bootstrap()` both warn on either, saying the candidate was tested by
  **neither** criterion and which mechanism applied. The
  `information_indefinite` prose claimed "that refit also failed", which was
  wrong for the new case and sent readers to an empty `refit_failures`; both
  it and the bootstrap warning now describe the two mechanisms separately.

  No behaviour changed in what gets selected: a candidate that could not be
  tested is still not entered. What changed is that the run says so.

* **`fit$fit$weak` now distinguishes "no ridge" from "not checked".** The
  weak-identification diagnostic introduced in 1.2.7 returned `NULL` both when
  a fit had been examined and found well identified and when it could not be
  examined at all -- most importantly when no Hessian was available, which
  happens on an install without the suggested numDeriv package and on fits whose
  rows are left- or interval-censored, where the analytic Hessian declines by
  design. `NEWS` offered `fit$fit$weak` as the programmatic check, so on those
  installs it certified as well identified a fit nothing had looked at.

  The field now takes three values: a list when a ridge was found, `NULL` when
  the fit was examined and is well identified, and `NA` when the check could
  not run. Test it with `is.list(fit$fit$weak)` rather than `!is.null()`.
  `summary()` prints a note in the `NA` case saying the check did not run, and
  a fit imported with `hzr_read_outhaz()` -- which has no R Hessian to examine
  -- now reports `NA` rather than reading as certified clean.

* **A fit with more than one flat direction says so.** The detector reported a
  single direction, which invited reading every parameter it did not name as
  identified. It now reports `n_directions`, the number of distinct
  near-flat directions found, and the warning says when there is more than
  one. Distinctness is counted over the parameters spanning each direction,
  not over eigenvectors: a ridge between two parameters clears the gate twice,
  once on the flat direction and once on its stiff partner, so counting
  eigenvectors would report one ridge as two.

* **A rescued candidate that goes on to win is no longer fitted twice.** The
  Wald fallback refits each candidate it rescues, and the acceptance step then
  refit the winner again with identical arguments. The two fits were
  bit-identical, so this was cost rather than incorrectness -- but it doubled
  the price of every accepted fallback entry, against a criterion whose whole
  advantage is that it does not refit per candidate. The rescuing fit is now
  kept and reused, as the Wald path already did with its candidate fits.

* **`DELTA` is not implemented, and two comments said it was absorbed.** The
  headers of `R/decomposition.R` and `R/argument_mapping.R` both stated that
  the C `DELTA` parameter's time transformation
  `B(t) = (exp(delta * t) - 1) / delta` is "absorbed by `decompos()`". It is
  not absorbed; it is unimplemented, and `delta = 0` is assumed. `DELTA`
  enters the C reference in three separate places -- it builds `rho` from
  `B(t_half)` rather than `t_half`, it replaces the time argument with
  `B(t)`, and `delta * t` enters the log-density additively so the density
  carries a factor of `exp(delta * t)` -- and R computes the `delta = 0`
  branch of all three.

  The comment was the harmful part. It made the omission look deliberate and
  safe, so a reader looking for exactly this discrepancy was told to stop
  looking, while a `PROC HAZARD` job with `DELTA != 0` was reproduced against a
  different function with no error. Both comments now say what is true.

  The SAS-facing paths now distinguish the two cases rather than treating
  `DELTA` as one unmapped keyword. `hzr_read_outhaz()` already stopped on a
  non-zero `DELTA`; `hzr_translate_sas()` now records `PARMS DELTA = <nonzero>`
  as untranslated with a reason saying the emitted call fits a *different*
  model, and the `.lst` natural-estimates parser warns when a listing carries
  one. `DELTA = 0` and a bare `FIXDELTA` are treated as faithful translations,
  because that is the branch R implements -- previously both the safe and the
  unsafe case produced the identical generic note "PARMS keyword has no phase
  target", which distinguished nothing.

* **A multiphase fit now records whether Conservation of Events was actually
  applied.** CoE counts exact events, so it is disabled whenever any `status`
  falls outside \{0, 1\} -- which interval or left censoring guarantees -- and
  whenever the model has fewer than two phases. That is deliberate and
  correct. What was missing is that **nothing on the returned object said it
  had happened**: a caller who passed `control = list(conserve = TRUE)` got a
  fit carrying `conserve = TRUE` over a computation that did not run.

  It is not an edge case. `ICENSOR` appears on 42 to 74 blocks per production
  study, and `ICENSOR` guarantees `status` leaves \{0, 1\}, so the auto-disable
  fires constantly. Production also writes `NOCONSERVE` on 16 to 68 blocks per
  study, so R and SAS usually agree on the *outcome* -- but for different
  reasons, and a job carrying both `CONSERVE` and `ICENSOR` is exactly where
  the two could diverge unobserved.

  A fitted multiphase object now carries two new fields, both alongside the
  requested `conserve` under `fit$spec$control`:

  * `fit$spec$control$conserve_applied` -- logical, whether CoE was actually
    applied;
  * `fit$spec$control$conserve_disabled_reason` -- one of `"not_requested"`,
    `"unsupported_censoring"`, `"single_phase"`, `"no_events"` or
    `"setup_failed"`, and `NA` when CoE was applied.

  Read `fit$spec$control$conserve_applied`, not `fit$spec$control$conserve`:
  the latter says only what you asked for. `conserve` is a
  `dist = "multiphase"` control; the single-distribution fits do not use it.

  The reason is recorded rather than a bare logical because the causes want
  different responses -- and a bare `FALSE` reads as "you turned it off" to a
  user who did the opposite.

* **The SAS `.lst` nomogram parser no longer returns the first of several
  tables silently.** `.hzr_parse_sas_nomogram()` matched every nomogram header
  in a listing and read only the first, with no warning and nothing in the
  return value to say a second existed. That is the same shape as the three
  layout defects fixed in 1.2.1 -- silent, and discoverable only by pointing
  the parser at a second study.

  It now warns when a listing holds more than one, and attaches `n_found` to
  the returned frame whatever the count, so a caller can distinguish "one
  nomogram" from "the first of several" without reading the source. The
  behaviour is otherwise unchanged: the first table is still what comes back.

  Every file in the corpus this was found against happens to print exactly one
  nomogram, so the parser got the right answer there -- by luck of the corpus
  rather than by construction.

## Documentation

* **Corrected the documented paths for two fit fields.** The weak-identification
  result is at `fit$fit$weak` and the phase shares at `fit$fit$phase_share`;
  `NEWS.md` and, for the shares, `?hazard` had both named them one level too
  high. Code following the documented recipe got `FALSE` from
  `is.list(fit$weak)` for every fit, ridge or not -- the same wrong answer the
  1.2.8 three-value change was made to prevent, reached by a different route.

* **`stat_type` no longer over-claims which rows the Wald fallback rescued.**
  Under `criterion = "score"` drop rows are always Wald-tested, following SAS,
  so they read `"wald_z"` whether or not the fallback fired. The rescued rows
  are the *entry* rows: filter on
  `action == "enter" & stat_type == "wald_z"`.

* **`predict()` and `hzr_nelson()` now say that SAS draws narrower bands.**
  SAS `%KAPLAN`, `%NELSONT` and `PROC HAZPRED` all take their band width from
  `CLEVEL`, whose default is `0.68268948` -- documented in the macro source as
  "(1 sd)". That makes the multiplier `1` to seven decimals -- the literal is
  truncated -- so the band is one standard error, 68.3%, not 95%.

  Nothing here computed the wrong thing: parity is tested and passing, and
  `hzr_kaplan()` already documented the convention. The gap was that the other
  three entry points did not, and they are the ones a reader meets when
  checking an R fit against an existing SAS figure. At the R default of
  `level = 0.95` the reproduced band is about 1.96 times wider than the one
  being checked against, with no error on either side -- so the two look like
  they disagree numerically when they do not.

  The defaults are unchanged. `0.95` is the right R-side default, and adopting
  SAS's silently would make `predict()` disagree with every other R modelling
  function. The help pages now carry the level to pass instead:
  `level = 2 * stats::pnorm(1) - 1`.

* `summary()`'s documentation and the *Inference and diagnostics* vignette
  both listed the notes the method prints and had not been updated for the
  ridge note. Both now include it, and the vignette says plainly that a flat
  direction means the point estimates along it are unreliable, not only their
  standard errors.

* **`hazard()` now states that SAS's `STEEPEST` has no equivalent.**
  `PROC HAZARD` jobs write `STEEPEST QUASI` together -- steepest descent, then
  quasi-Newton -- and `STEEPEST` appears 14 to 109 times per study across the
  corpus. `QUASI`/`QUASINEWTON` maps to `method = "bfgs"`; there is no
  steepest-descent option and no two-stage strategy. Since the multiphase
  likelihood is multimodal, a different descent path can land on a different
  optimum, so a fit translated from such a job may not reproduce SAS's
  estimates. `hzr_translate_sas()` already recorded the keyword as
  untranslated rather than dropping it; the `control$method` documentation now
  says why.

## Internal

* Two tests in `test-score-wald-fallback.R` did not catch the mutations their
  comments named. Widening `.hzr_score_fallback_reasons` to include `constant`
  and `collinear` left both green, because a degenerate candidate still fails
  to enter -- it merely costs a refit on the way out -- and the "noise stays
  out" assertion is guarded by `slentry` rather than by how narrow the
  fallback is (the fixture's own Wald p-values are 0.0997 and 0.149, so a
  fallback that refit everything would still decline both). Both now assert
  `n_wald_fallbacks`, which is the quantity that moves.

* A roxygen block in `R/score-test.R` bound to the character vector declared
  after it rather than to the function it documents. `@noRd`, so no Rd was
  affected; source readability only.

# TemporalHazard 1.2.7

## New features

* **A fit sitting on a likelihood ridge now says so, and names the
  parameters.** An ill-conditioned Hessian already warned that standard
  errors were unreliable. That understates the problem when the
  ill-conditioning is a ridge: the likelihood is near-flat along some
  combination of parameters, and there the individual *point estimates* are
  not determined by the data either -- only the combination is. The fit still
  reports `converged`, and the coefficient table still prints a number for
  every parameter, so nothing on the object signalled it.

  `hazard()` now warns, once per fit, naming the parameters that span the flat
  direction along with their correlation and the Hessian's `rcond`, and
  `summary()` prints the same note. The finding is recorded on the object as
  `fit$fit$weak`, so it can be checked programmatically rather than scraped from a
  warning. See the 1.2.8 notes below for the three values that field takes.

  The direction is read off the *correlation* of the estimates rather than
  their raw covariance. Parameters here sit on very different scales -- an `m`
  of 27 against a `nu` of 0.027 -- and in raw units a direction that moves
  both equally in statistical terms loads almost entirely on the larger one,
  which would report a two-parameter ridge as a single unidentified parameter.
  A parameter that is merely imprecise, without trading off against another,
  is deliberately not reported: that is ordinary low precision, and the
  existing `rcond` warning and the parameter's own standard error already
  cover it.

  The check is generic -- it runs for every distribution and knows nothing
  about phase shapes -- and is gated on the `rcond` threshold the package
  already uses, so it never fires where the ill-conditioning warning stays
  silent.

## Bug fixes

* **The score criterion no longer declines a candidate for being too
  predictive.** `criterion = "score"` computes SAS HAZARD's Q exactly --
  `Q = grad^2 * I22`, the reciprocal Schur complement of the *observed*
  information at `beta = 0` (`src/vars/q1.c`). When a candidate's true effect
  is far from zero the log-likelihood is convex there, the Schur complement
  turns negative, and Q is undefined. The criterion therefore declined
  candidates in proportion to how predictive they were: on a fixture with one
  planted effect (`beta = 0.9`, LR = 178) and two pure-noise columns, the
  screen entered both noise columns and never tested the real one (#130).

  This is not a deviation from the reference -- it is inherited from it. SAS
  documents the same failure in `q1.c` ("IT IS POSSIBLE THAT THE PROGRAM WILL
  RETURN A NEGATIVE Q VALUE ... THE USER SHOULD USE THE MORE EXPENSIVE Q2 AS
  AN ALTERNATIVE"), and `dqstat.c` declines the candidate with `p = 1`. `Q2`
  is named once in the C tree and never implemented.

  So Q itself is unchanged and stays bit-faithful; only the *handling*
  diverges, and only where SAS says its own answer is unusable. A candidate
  the score cannot test is now refit and tested by Wald -- the substitute the
  unbuilt `Q2` was for. This is a deliberate, documented divergence from the
  reference implementation.

  The fallback is deliberately narrow: it applies to `information_indefinite`
  and `coefficient_diverging`, the two causes that mean "the approximation at
  zero broke down". Collinear, constant and non-numeric candidates are still
  declined without a refit, so the screen keeps the speed advantage that the
  score criterion exists for -- the cost is paid only on the few candidates
  that trip it. The returned object's `$criteria` gains `n_wald_fallbacks`,
  so the substitution is reported rather than silent, and a fallback refit
  that fails is recorded in `$criteria$refit_failures` and warned about
  rather than leaving a row indistinguishable from one never refit.


# TemporalHazard 1.2.6

## Bug fixes

* `.hzr_parse_sas_nomogram()` no longer discards a nomogram whose PROC PRINT
  counter column is labelled `OBS` rather than `Obs` (#184). The counter was
  dropped by exact name, so on the other casing it stayed among the header
  names, the row-width guard rejected every data row, and the parser returned
  `NULL` -- indistinguishable from a listing that printed no nomogram at all.
  A corpus sweep therefore reported its own parse failures as gaps in the SAS
  output. The counter is now identified structurally, as a leading column
  running exactly `1..n`, so any label parses.

* `.hzr_parse_sas_nomogram()` now warns rather than returning `NULL` in
  silence when a nomogram header matched but no data rows could be read.
  `NULL` again means "no such table", and only that.

# TemporalHazard 1.2.5

## New features

* A multiphase fit now warns when a phase has effectively left the model. Such
  a phase is silent in every other way: the fit converges, reports no trouble,
  and the affected parameters simply drift.

  Two modes are distinguished, because the consequences differ. A phase that is
  **absent** -- contributing essentially none of the cumulative hazard at any
  observed time -- has not started by the end of follow-up, and neither its
  `mu` nor its shape is identified. A phase that is **saturated** -- one whose
  contribution is constant across the observed times, typically a `cdf` phase
  whose half-life is far shorter than the first observation -- has already
  finished, and then acts as a constant offset: its `mu` stays well identified
  while the shape parameters (`t_half`, `nu`, `m`) go exactly flat. Pinning
  those at any value leaves the log-likelihood unchanged.

  The distinction is the point. It is tempting to describe a phase that
  supplies no late hazard as one whose `mu` has stopped being identified;
  `mu` is in fact the one parameter that survives, through the offset the
  phase already contributed. The share is measured against the cumulative
  hazard rather than the instantaneous hazard for the same reason.

  The shares are recorded on the fit as `fit$fit$phase_share`, so the warning can
  be checked rather than taken on trust, and the threshold is
  `control$phase_share_tol` (default 1e-8).

# TemporalHazard 1.2.4

## New features

* `hazard()` gains an `objective` argument. The default, `"likelihood"`, is
  unchanged: interval-censored rows contribute the interval probability
  `log(S(l) - S(u))`. The new `"sas"` reproduces what `PROC HAZARD` actually
  accumulates for such a row -- the ordinary event-density term with the
  *instantaneous* hazard replaced by the *interval-mean* hazard over (l, u]:

  ```
    d * log[ S(u) * (Lambda(u) - Lambda(l)) / (u - l) ]
  ```

  which makes the three row types one family: right-censored contributes
  `log S(u)`, an exact event `log S(u) + log h(u)`, and an interval-censored
  row `log S(u) + log h_bar(l, u]`. Exact-event and right-censored rows are
  untouched by the switch, and it applies only to `dist = "multiphase"`.

  **This is a different estimator, not a reparameterization, and must not be
  used for new analyses.** It is a density, not a probability, so it is
  inconsistent for wide intervals -- on a 12-year-interval reference fit the
  two forms differ by 22 log-likelihood units. It exists to reproduce legacy
  `PROC HAZARD` runs, and it is deliberately an explicit top-level argument
  rather than a `control` element, because it changes the estimand.

  Interval-censored rows with `u <= l`, and any left-censored row, are errors
  under `"sas"` rather than silently-substituted values: `PROC HAZARD` has no
  left-censoring statement, so no SAS run corresponds to such a result.

* New dataset `uslife2023`: the NCHS United States life table for 2023 on a
  synthetic 100,000 radix, 124 rows, every one interval-censored and exactly
  one year wide. Published aggregate counts only. It is the reference fixture
  for `objective = "sas"`, which reproduces its SAS log-likelihood of -410414
  at the printed estimates and at three off-optimum points of SAS's own
  iteration trace.

## Internal

* The interval-censored contribution was written twice -- once in the
  log-likelihood and again in the finite-difference closure inside the
  gradient. Those copies had to agree or the optimizer would step by the
  gradient of a different objective than it evaluated. Both now delegate to a
  single `.hzr_logl_interval()`. Behavior under the default is unchanged and
  bit-identical, log-likelihood and gradient alike.


# TemporalHazard 1.2.3

## Bug fixes

* A logical column no longer stops a stepwise screen dead. `hzr_stepwise(scope
  = NULL)` enumerates its own candidates and counts logical columns among them,
  on the grounds that a 0/1 field arriving logical rather than numeric is a
  property of the reader that produced the frame, not of the variable. Both
  criteria then refused what the package had offered: the score criterion
  errored with "is not numeric (logical)", and the Wald criterion failed
  looking up `phase.var` when `model.matrix()` had named the column
  `phase.varTRUE`. Either way the screen stopped before its first step, on a
  column nobody had chosen by hand.

  A logical candidate is now modelled as the 0/1 predictor it is, and gives
  the same screen as the identical numeric column: same variables entered, in
  the same order, at the same p-values, with the same coefficients.

  The coefficient-name half of this also fixes a two-level factor, which
  expands the same way (`varb`). A candidate that expands to more than one
  column is still refused, unchanged.

  That makes `criterion = "wald"` a real answer for a two-level factor or
  character column named in an explicit `scope`, which the score criterion
  still cannot expand. Its refusal used to say switching criterion would not
  help, and now points at it instead.



# TemporalHazard 1.2.2

## New features

* `hazard()`'s vector interface now evaluates `time`, `status`, `time_lower`,
  `time_upper` and `weights` in `data`'s scope. `hazard(data = df, time = tt)`
  previously failed with `object 'tt' not found`: `data` was consulted only by
  the formula path, and the vector path accepted it and ignored it. The rule is
  `subset()`'s -- a column of `data` wins, and `df$col`, a local vector or a
  literal falls through to the calling frame unchanged -- so with `data = NULL`
  nothing changes and the formula path is untouched.

  Because a column winning can silently redirect a wrapper that forwards its
  own argument by name, `hazard()` now **warns**, once per call, when a symbol
  is both a column of `data` and visible from the calling frame -- that frame
  or a lexical parent of it, up to and including the global environment --
  naming every such symbol and the argument it appeared in. `data` must now be
  a data frame or a list: `hazard(data = <matrix>, ...)` errors, where a matrix
  was previously accepted and silently ignored along with everything else in
  `data`.

* `hzr_translate_sas()` translates a SAS `PROC HAZARD` / `PROC HAZPRED` job
  into a Quarto document of equivalent R calls. It parses the SAS statements,
  builds the calls, and renders them into `.qmd` chunks -- the model state is
  stored as unevaluated calls, so rendering is `deparse()`, not string
  templating.

  **This function is experimental.** A job that translates now renders: the
  emitted `hazard()` chunk binds its fit to a name and asks for an actual
  fit, and the `predict()` chunks have something to predict from. Measured
  on the public `hazard` corpus of 110 `.sas` files, 57 translate into 22
  distinct documents; the 11 of those that synthetic data can drive end to
  end evaluate every chunk and bind a converged fit, and the other 11 --
  `PROC HAZPRED`-only jobs with no local fit to bind -- are exercised up to
  their fit chunks. Read that as a measurement, not as "the translator
  works": the rest of the corpus is refusals or jobs whose external `INHAZ=`
  could not be resolved. It remains a translation aid rather than a turnkey
  reproduction, and the API, the `hzr_sas_job` field layout and the emitted
  document format may all still change.

  Two SAS constructs are refused outright rather than mistranslated into
  something that computes a wrong answer. Each records an `UNTRANSLATED` row
  and emits a `stop()` in place of the fit, so the document fails where the
  fit would have been:

  - a `SELECTION` statement requesting a stepwise screen. `hzr_stepwise()`'s
    refit path needs a formula-interface base fit and this translator emits
    the vector interface, so every candidate refit would error and the screen
    would report zero steps -- indistinguishable from "nothing met
    `slentry`" (#152, #160; the underlying `hzr_stepwise()` silent no-op is
    #159).
  - `LCENSOR` combined with `ICENSOR`. `hazard()`'s single `time_lower`
    argument carries the entry time for status 0/1 rows and the interval's
    lower bound for status 2 rows, so one column cannot express both (#155).

  Two gaps that made the emitted calls compute a different answer from the
  SAS job are closed. The log prediction grid now takes its step from the
  job's own `INC=` expression rather than a hardcoded one (#153): three
  denominators appear across the public corpus (`/49.9`, `/99.9`, `/999.9`)
  and the denominator sets both the step and the number of points SAS's
  `DO lo TO hi BY INC` lands, so reading every job as `/99.9` gave the
  `/999.9` jobs 100 points on a step ten times too large -- every time
  wrong, over an empty `untranslated` and full coverage. The span is the
  loop's own `log(hi) - lo`, which coincides with `5 + log(hi)` only because
  every corpus job starts at -5, and an `INC=` in a form the parser cannot
  read is refused with an `UNTRANSLATED` row rather than stepped by a guess.
  Every log grid in the public corpus also writes its `DO` statement with an
  explicit trailing element (`DO lo TO hi BY INC, hi`), which SAS's `DO` list
  syntax evaluates as the loop *plus* one final point at exactly `hi` -- the
  translator emitted only the loop, so every translated grid stopped short
  of the time the job actually asked for (about 8% short for a `/99.9`
  step). The trailing element is now read from the job's own `DO` statement
  and emitted as the grid's last point; a trailing element this cannot
  resolve to the loop's own bound is refused with an `UNTRANSLATED` row.
  `EVENT`, `ICENSOR` and `WEIGHT` are **counts** in the reference
  implementation, not flags, and `setlik.c` combines one record's
  contribution as `c1c2c3 = c1w + c2 + c3w` with `c1w = C1 * WT` and
  `c3w = C3 * WT`. All three now reach the fit that way. An `ICENSOR` event
  count is no longer discarded (#154); an `EVENT` count carries into
  `weights` and `status` derives from `EVENT > 0`, where `EVENT = 2` used to
  map straight onto `status = 2` and be fitted as **interval-censored** --
  a different likelihood branch, not an under-count (#157); and a `WEIGHT`
  variable no longer weights right-censored rows, because `c2` is the one
  term entering that sum unweighted and `readc2.c` sets it to `1` on exactly
  those rows -- a `WEIGHT` that was `0` there previously deleted them from
  the fit silently (#158). A row where the `EVENT` and `ICENSOR` counts both
  fire is two contributions at once, which one `status` and one `weight`
  cannot express, so the emitted status chunk now stops before the fit
  rather than picking the event branch and discarding the interval one.

  `RCENSOR` is the third of those counts and was being ignored outright
  (#162). It names `C2` -- "COUNT OF CENSORED INDIVIDUALS AT TIME=T" --
  and when a job names it, `readc2.c` reads the column straight from the
  data and skips the `C2 = 1` derivation that a job without `RCENSOR` gets.
  Four censored individuals were therefore fitted as one observation. The
  censored branch of `weights` is now that variable, still unmultiplied by
  `WEIGHT`, and the both-fire guard now covers `EVENT` + `RCENSOR` and
  `ICENSOR` + `RCENSOR` as well: `readobs.c` deletes an all-zero row only
  when `RCENSOR` is named *and* exactly one of the other two is, so a row
  with two counts positive always survives to be summed.

  A `0/1` `RCENSOR` flag that is exactly `1 - EVENT` -- which is what both
  corpus jobs carry, and what the statement is usually used for -- fits
  exactly as before. A `0/1` flag that is **not** its complement does not,
  and both ways it can differ are deliberate: a row with the event and the
  censoring flag both set now stops the document instead of being fitted as
  an event alone, and a row with neither set now carries weight `0` instead
  of a fabricated weight of `1`, matching the row SAS would have deleted.
  A count column that is **negative or missing** on any row now stops the
  document with a message naming the variable, rather than being folded into
  the censored branch or propagating `NA` into `weights` until `hazard()`
  refused it as "non-negative and finite". `readc1.c`, `readc2.c` and
  `readc3.c` apply the same rule to every count the job names -- a missing
  value sets `mdel`, a negative one sets `del` -- and `readobs.c` then skips
  `setobs()` for that row and subtracts it from `Nobs`. Such a row
  contributes nothing at all, so translating it as a right-censored
  observation of weight `1` adds survival mass at a time SAS had removed.
  This is the one case where a missing count is **not** interchangeable with
  a zero one: a zero count is kept and contributes, a missing one is deleted.
  The translator cannot drop rows without changing `n` behind the reader's
  back, so it stops and says to filter them.

  Every translated job now emits a `status` chunk ahead of its fit, where
  these guards live -- previously only jobs with `ICENSOR` or more than one
  named count did. The emitted document format remains experimental.

  Loading a fit from an external `INHAZ=` dataset returns a classed
  `hzr_outhaz` object with a `predict()` method (#151). That method takes
  the same arguments in the same order as `predict.hazard()` --
  `newdata`, `type`, `decompose`, `se.fit`, `level`, `conf.type` -- so a
  positional call means the same thing for both methods of the generic, and
  `conf.type` (the `PROC HAZPRED` parity switch the translator emits) is a
  real argument rather than one that a misspelling could drop into `...`,
  returning the log-log limits the SAS job did not ask for. Its *value* is
  checked only on the survival standard-error path that reads it, exactly as
  `predict.hazard()` does, so an ignored value does not fail a point or
  hazard prediction; a mistyped argument *name* still errors. `type` defaults to
  `"hazard"`, as in `predict.hazard()`, and `decompose = TRUE` is an error:
  an `OUTHAZ=` dataset carries fitted parameters, not a per-phase
  decomposition. Point predictions work; `se.fit = TRUE` is **refused** whenever the SAS fit *estimated* a
  late shape parameter that `PROC HAZARD` put on a composite scale --
  `log(GAMMA*ETA - 2)` and friends, which is the generic unconstrained
  three-phase case rather than an exotic one -- and likewise under `FIXMNU1`
  or where one late parameter is derived from another. A translated `PROC
  HAZPRED` block requests confidence limits unless the SAS job says `NOCL`,
  so such a job stops at its `predict()` chunks with an explicit message
  naming the parameter and its scale, rather than reporting standard errors
  built on the wrong one.

  Treat `job$coverage` as a measure of *parsing* -- tokens recognised --
  not of whether the result runs. The
  parameter translation itself is verified separately: refitting the
  `hz.death.AVC.sas` job's parameters through `hazard()` directly reproduces
  the SAS log-likelihood to the six significant figures the reference
  listing prints (`-210.501`).

  The keyword grammar behind the parser -- 122 keyword rules, 67 of them
  mapped to an R target -- is **generated from the reference
  `HAZARD`/`HAZPRED` C implementation's own lex sources**
  (`data-raw/hazard-grammar.R`), not hand-written. Only the extracted table
  ships; no GPL-2 source enters the tarball. A hand-written table would
  capture only the spellings a study happened to use, and the grammar has
  real context-dependent collisions --
  `M` means a phase shape parameter inside `PARMS` and `MOVE` inside a
  `PHOP`/`STEP` statement -- that a context-free lookup gets silently wrong.

  Constructs the translator does not cover are recorded on the returned
  `hzr_sas_job` object and rendered as visible `UNTRANSLATED` callouts in
  the `.qmd`, never dropped. Two limits are worth stating plainly:

  - **Prediction grids built from `SET`-derived values, function calls, or
    unknown names are not translated.** The parser resolves a `PROC
    HAZPRED` grid's `DO` loop bounds when they are literal numbers or
    DATA-step constants it can fold (e.g. `DO MONTHS = 1*DTY, 2*DTY, ...;`
    with `DTY` assigned earlier in the same DATA step) -- but a bound
    read from `SET`, computed by a function call, or naming something the
    parser can't resolve is refused whole rather than partially read: a
    partially read grid is a partial `newdata`, which is a hollow result.
    Such grids emit an explicit `UNTRANSLATED` block instead, and the
    `predict()` chunks that would have read the grid become a `stop()`
    naming it: emitting `predict(fit, newdata = <name>)` that nothing
    builds either fails on an unbound name or, if the rendering session
    happens to hold an object of that name, reports predictions over
    unrelated times. On the
    public corpus, grid resolution is 19 of 55 (35%), up from 10 of 55
    (18%) before constant folding.
  - **An unresolved `INHAZ=` fails the render, on purpose.** A `PROC
    HAZPRED` job whose fitted-model dataset can't be located -- neither
    from another translated job's `OUTHAZ=` nor from the `librefs`
    argument -- gets an `inhaz-unresolved` chunk, ahead of the grid and
    `predict()` chunks, whose whole body is a `stop()` naming the
    unresolved libref. The document fails to render rather than reporting
    predictions over a model it never loaded.

## Bug fixes

* `$se` on a fitted object is now one standard error per parameter, whatever
  the variance matrix looks like. A multiphase fit legitimately carries NA
  variance rows for the parameters it holds fixed, and a single NA anywhere in
  the matrix collapsed the whole vector to a length-1 `NA`. A five-parameter
  fit came back with a length-1 `$se`, so naming the standard errors against
  the parameters failed with "'names' attribute [5] must be the same length as
  the vector [1]".

  A scalar `NA` now carries the meaning it already has for `vcov()`, that there
  is no variance matrix at all. Where a matrix is present, `$se` is computed
  element by element: a parameter held fixed carries an NA variance row and
  earns an `NA` standard error, while every parameter whose variance *was*
  computed keeps its own. Sizing the vector to the matrix is not enough on its
  own -- filling it with `NA` throughout would leave `$se` conformable and
  empty, and contradicting `summary()`, which reads the same matrix and reports
  those standard errors. `summary()` was never affected either way, so this was
  a quiet inconsistency on the fit object rather than a visible break.

* `hzr_decompos()` no longer returns a wrong value for large `|m|`. The three
  branches with a nonzero `m` all formed `2^m` and the terms built from it
  all three lost the answer well inside the range a fit can reach.

  For `m > 0` the failure is overflow. `2^m` is `Inf` from `m = 1024`, but
  `bt^(-1/nu)` goes first: at `t/t_half = 0.5` with `m * nu = 3` it overflows by
  `m = 750`, and sooner for `nu > 1`. Either way `btnu` becomes `Inf`, and
  `Inf^(-1/m)` is `0`. So `hzr_decompos(0.5, t_half = 1, nu = 3/1000, m = 1000)`
  reported `G = 0` where the answer is `0.3969`, with `g` and `h` `NaN`. Nothing
  warned. `G = 0` is a perfectly ordinary probability, and a fit whose optimizer
  wandered into large `m` used it. The `nu < 0` branch collapsed the same way,
  to `G = 1`.

  For `m < 0` the failure is cancellation instead. `1 - 2^m` rounds to exactly
  `1` once `2^m` falls below machine epsilon, so `(1 - 2^m)^(-nu) - 1` is `0`,
  `rho` is `Inf`, and `G` is again `0`. That collapse is at `m = -53`, which an
  optimizer reaches much more easily than `m = 750`, and the accuracy decays
  before it: at `m = -20` the old code was already wrong in the tenth digit.

  All three branches now work on the log scale. The `m` in the `(2^m - 1)/m`
  factor of `rho` cancels the explicit multiplier, so `m * bt^(-1/nu)` is
  exactly `(t_half/t)^(1/nu) * (2^m - 1)`, and `log(btnu)` follows from
  `hzr_log1pexp()` applied to the log of that product. The `m < 0` branch takes
  `log(1 - 2^m)` from `hzr_log1mexp()` rather than forming the difference. Both
  primitives were already in the package. Checked against a reference computed
  at 100 or more decimal digits, `G` and `g` are now accurate to machine
  precision from `m = -1000` up to `m = 5000`, they track the analytic
  large-`m` limit at `m = 1e6`, and they are unchanged where the old code was
  already right. Small `m` improves too, by eight orders of magnitude or more:
  `log(2^m - 1)` is taken as `x + hzr_log1mexp(x)` for `x = m * log(2)`, which
  holds at both ends, where the direct `log1p(-2^(-m))` decays from about
  `m = 1e-3` down and reaches `-Inf` once `2^(-m)` rounds to `1`.

  One boundary remains, and it is now visible rather than silent. Below about
  `m = -1074` the term `2^m` underflows outright and no rearrangement recovers
  it in double precision; `hzr_decompos()` returns `NA` there.

  The multiphase log-likelihood is evaluable again over the same range. On a
  two-phase fixture it returned `-Inf` from about `m = 450`, for this same
  reason: the smallest observed time is what makes `log(t_half/t)` largest, so
  the overflow arrives earlier than the `m = 750` above. That is what made the
  likelihood surface along the ridge `m * nu = const` hard to characterize.

* `hzr_stepwise()` can no longer return a zero-step result that is silently
  empty. Every accepted move goes through a refit, and a refit that failed was
  downgraded to a warning and then dropped: the returned object carried no
  record of it, so a screen that could not fit a single candidate looked
  exactly like one that tested them all and liked none. `$criteria` now carries
  `refit_failures`, `n_refit_failures` and `stopped_refit_failed`; a run that
  ends on an iteration with failed refits warns that its candidates were never
  tested; and the trace names the cause instead of claiming "no further
  action". A base fit built with the vector interface (`time =` / `status =`)
  stores no formula for the refit to mutate, so every candidate would fail --
  `hzr_stepwise()` now rejects it up front with one message naming the remedy,
  through the same predicate the refit itself uses.

* A multiphase fit is now reproducible. `hazard(dist = "multiphase")` offsets
  the starting values for every optimization start after the first, and those
  offsets were drawn from the ambient RNG stream. The identical call run twice
  returned a different answer: on a 150-row two-phase fit the estimates moved
  by about 0.3 on the log scale and the objective by about 0.08, which is
  enough to change what the fit says. Fitting also advanced the caller's
  stream, so a later `sample()` or `rnorm()` depended on whether a model had
  been fitted first.

  Where the assembled starting values did not converge on their own, the draw
  decided whether there was a fit at all: the fit succeeded only from a
  perturbed start, and about a quarter of draws stopped outright. The same call
  could raise that error on one run and not the next. The reason those starting
  values failed is fixed below, so the draw no longer decides that; it decides
  only which optimum is reached.

  The offsets now come from an internally seeded stream, and the ambient
  `.Random.seed` is restored afterwards. The same data and the same control
  give the same fit, with no `set.seed()` needed, and fitting leaves the
  caller's stream where it found it. The new `control$start_seed` (default 3)
  selects a different ensemble of starts. That is worth reaching for when a fit
  looks like it settled in a local optimum: fit at a few seeds and compare the
  `objective` values. See the note on `starts` below for how to read two
  objectives that differ, which is not always a pair of rival optima.

  `start_seed` takes any whole number within integer range, negatives included
  -- `set.seed(-1)` is perfectly valid and deterministic, so restricting to
  non-negative values would discard half the seed space for no reason. A
  fractional value is rejected rather than truncated: `set.seed()` truncates,
  so `3.9` and `3` would select the same ensemble, and a sweep over
  `3.1 / 3.5 / 3.9` would report three fits having tried one set of starts.
  Coercing quietly would keep that aliasing and merely move it. A value out of
  integer range is rejected too, because `set.seed()` would otherwise fail with
  "supplied seed is not a valid integer" and name neither the argument nor the
  fit it came from.

  `hzr_bootstrap()` draws its own resample before each refit, so replicates are
  still distinct. Its numbers do shift, because the refits no longer advance
  the stream between resamples, and a run with `seed=` is now reproducible end
  to end.

* A multiphase optimization start no longer dies on an infeasible shape. The
  multiphase cumulative hazard short-circuits to an infinite hazard when a
  phase's `m` and `nu` are both negative, so the optimizer sees a penalty and
  backs out of the region. Asked for a per-phase decomposition it
  short-circuited in the wrong shape -- a bare vector where the caller expects
  a named list -- and the Conservation-of-Events adjustment, which runs inside
  the objective on every evaluation, raised `$ operator is invalid for atomic
  vectors`. BFGS steps into that region routinely, so the error came back out
  of `optim()` and the multi-start loop threw the whole start away.

  A discarded start was reported as a failure to converge, so a crash read as a
  numerical problem. On the two-phase fixture in `test-multiphase-gradient.R`
  it cost the fit its assembled starting values outright: `n_starts = 1`
  stopped with an error, the fit survived only on a perturbed start, and 12 of
  50 `start_seed` values failed. All 50 converge now, `n_starts = 1` converges
  on its own, and across 50 seeds every one of the 250 starts is usable.

* A multiphase fit now says which of its starts survived. `fit$fit$starts`
  gives one row per optimization start: its `status`, its `objective`, its
  `convergence` code from `optim()`, whether it was the `best` one and so the
  fit you are looking at, and the `message` of any error it raised. Worth a
  look when a fit is in doubt: on the fixture above the assembled start reaches
  -159.15 and a perturbed start -158.30, and start 1 wins 5 of 50 seeds at the
  default `n_starts = 5` (17 of 50 at `n_starts = 3` -- the rate depends on how
  many starts there are to lose to, so read it against your own setting).

  Read two such numbers as objectives, not as rival optima. That fixture has no
  interior maximum in `m`. Profiled, its objective climbs past -158.30 toward a
  finite limit of -157.88 that is reached only as `m` grows without bound, so
  the better number is a point on a flat ridge where the optimizer met its
  tolerance, and its standard error on `m` is 42.8 against an estimate of 27.2.
  Starts that disagree like that are telling you the shape is barely
  identified, which is the reading `starts` is there to support. On data that
  does identify the shape the picture is the ordinary one: the `avc`
  early+constant profile has an interior maximum near `m = 1` and falls away on
  either side.

  `status` separates the four ways a start can end, and in particular a start
  that stopped at `maxit` reads as `"nonconverged"`, not `"ok"`. That
  distinction is not cosmetic: `optim()` attaches a perfectly finite objective
  to a run it abandoned at the iteration limit, and such a start can carry a
  better objective than one that genuinely converged and so become the
  reported fit. Which start wins is unchanged -- it is still the best
  objective -- but you can now see whether it converged. `fit$fit$converged`
  continues to report that for the fit as a whole.

  A start that errors now also warns rather than being absorbed, and when every
  start fails the error names what was raised instead of calling it a
  convergence failure. An error thrown inside the objective used to be
  indistinguishable from a start that merely optimized badly, which is how the
  defect above stayed hidden.


# TemporalHazard 1.2.1

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
  step sequences. Pass `criterion = "wald"` to restore the previous behavior
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
  `%HAZBOOT` procedure. **This is experimental**: the selection arguments and
  the shape of what they return may change in a future release, and
  `?hzr_bootstrap` says why under "Selection mode is experimental". The
  fixed-formula bootstrap (`scope = NULL`) is unaffected and unchanged.
  The short version: the design is still being read off production runs, and
  a screen large enough to matter runs for hours while this function writes
  nothing until its last replicate, so splitting a run across processes is
  currently the caller's job. Each replicate runs a fresh `hzr_stepwise()`
  selection (starting from a fixed-shape refit of the base model) instead
  of a plain refit, so `summary$pct` reports the variable's selection
  frequency across resamples and `summary$mean`/`sd`/`ci_*` describe the
  coefficient distribution conditional on selection. `scope = NULL`
  (the default) preserves the original fixed-formula bootstrap unchanged.

* `hzr_read_outhaz()` reads a `PROC HAZARD` `outhaz=` estimate dataset,
  returning the estimates, each parameter's free/fixed status, the
  variance-covariance matrix over the free parameters, and the model-structure
  flags. `outhaz` stores its numbers at full double precision where the
  printed `.lst` carries about seven significant figures, so for any quantity
  it holds it is the better parity reference -- print precision stops being
  the binding constraint and optimizer convergence takes over. The
  log-likelihood is not among them; that still comes from the `.lst`.

## Bug fixes

* **`hzr_stepwise()` never checked that an accepted step improved the fit.**
  A forward step enters a model that *contains* the one it started from, so at
  the optimum the log-likelihood cannot fall. It was written into `$steps` at
  every step and compared at none, so a step whose refit failed to converge
  entered anyway and every later step was then scored against a model that was
  not at its own optimum. In the production screen that surfaced this, three of
  ten steps lowered the log-likelihood and the run still reported convergence,
  ten entries and `p = 0.000` throughout; the final 19-coefficient model fitted
  57 units worse than the nested 16-coefficient model from three steps earlier,
  which cannot happen at a maximum.

  `$steps` now carries `delta_logLik`, a forward step that lowers the objective
  warns and is counted in `$criteria$n_nonmonotone_entries`, and
  `hzr_bootstrap(scope = )` reports `$n_nonmonotone_replicates` — a replicate
  whose path went backwards still contributes its selections to the pooled
  frequencies, and each replicate runs under `suppressWarnings()` so the
  step-level warning cannot reach the user. The comparison carries a small
  tolerance so optimizer noise does not fire it. Reported as issue #134.
* **A score statistic could be finite, enormous and meaningless.** A
  production screen accepted a candidate with `stat` = 92,211 on 1 df and
  `p = 0.000`, after which the refit made the model worse. Neither existing
  guard reached it: the adjusted variance stayed positive and well above the
  collinearity floor, so the statistic was reported as evidence.

  Near-collinearity alone does not do this — as a candidate approaches
  collinearity its score shrinks along with its variance and `Q` stays small.
  `Q` explodes when the model being scored against is *not at its optimum*,
  because the reduced-model score is then no longer zero: the numerator is
  inflated while the denominator stays small. That is the state a failed refit
  leaves behind. `hzr_stepwise()` now declines a candidate whose implied
  coefficient exceeds ±50, reporting `coefficient_diverging`, which is what
  the SAS/C reference has always done (`dqstat.c` rejects `|QBETA| > 50` as
  "the model is going to infinity"). Measured on the bundled `avc` data, a
  model displaced 0.25 from its optimum produced `Q` = 6.5e7 with no reason
  reported at all; a legitimate candidate reaches an implied coefficient of
  about 14, so the threshold has real headroom. Reported as issue #134.

* **The multiphase gradient and Hessian disagreed with the log-likelihood on
  left-truncated data.** For a row with `status` in `{0, 1}` the log-likelihood
  subtracts `H(time_lower)` unconditionally, but the analytic derivatives
  defined the entry time with an extra `time_lower < time` filter. A subject
  entering the risk set at its own event or censoring time was therefore
  differentiated as though it had no entry time, while its weight was still
  applied -- so the derivative was taken of a different function from the one
  being evaluated, and the optimizer left any sensible region immediately.
  Measured on `avc` at fixed parameters, the analytic gradient was out by 382
  where every row entered at its exit time, and by 126 where only *some* did
  -- which is ordinary left-truncated data, not a pathological input. Both
  derivatives now define the entry time exactly as the likelihood does, and
  new tests assert agreement with `numDeriv` across five entry-time layouts
  rather than the one the old filter happened to admit.

* **`hazard()` documented `time_lower` incorrectly, and now warns when it is
  self-defeating.** The argument was described only as the lower bound of a
  censoring interval, "defaulting to `time` if NULL". For `status` in
  `{0, 1}` it is in fact the counting-process **entry time**, and leaving it
  `NULL` means entry at `0`, *not* at `time`. Read literally, the old wording
  said that passing `time_lower = time` changes nothing; it in fact states
  that every subject left the risk set at the instant it entered, which
  removes every such row from the likelihood and leaves the objective
  unbounded above. The documentation now gives both roles, and supplying
  `time_lower >= time` on a `status` 0 or 1 row warns, naming the count and
  the `NULL` default. Reported as issue #136.

* **`hzr_stepwise()` now says *why* a candidate could not be scored, and warns
  when the reason is that the candidate looks strong.** Under
  `criterion = "score"` a candidate whose Q statistic cannot be computed drops
  out of the step, and the run previously reported only a count of them. Two of
  the causes mean opposite things. A collinear column should be dropped. But
  the observed information at `beta = 0` is not positive definite away from a
  maximum, and when a candidate's effect is *large* the log-likelihood curves
  upward there, the adjusted variance goes negative, and the candidate is
  declined -- so the criterion is least able to score exactly the variables a
  screen most wants to find. The old warning attributed both to "a degenerate
  or collinear candidate column", which tells a user to discard their best
  variable.

  `$criteria$uncomputable_reasons` (and `$uncomputable_reasons` on a
  `mode = "select"` bootstrap) now counts the causes by name, `$all_scores`
  carries a `reason` column per candidate, and both warnings name them. The
  reference implementation separates these too, and the R side now matches its
  split: a candidate whose *own* observed information is not positive is
  reported apart from one that is unusable only given what is already in the
  model (`information_nonpositive` against `collinear` and
  `information_indefinite`). The first is reachable on a multiphase fit with a
  large share of interval-censored rows. A run
  that *completed* while declining a candidate for this reason now warns too:
  it previously returned a selection -- sometimes an empty one -- in complete
  silence, which is the case where the omission is least visible. The
  underlying limitation of the score criterion is unchanged and is tracked
  separately; `criterion = "wald"` tests these candidates.

  One behavior change comes with it: the guard on the adjusted variance is now
  a magnitude test rather than a sign test. A variance within rounding distance
  of zero is reported as collinear whichever side of zero it lands on, and only
  a materially negative one is reported as indefinite. The previous floor was
  signed and relative to `I_bb`, so where `I_bb` was itself negative a slightly
  negative variance passed through and produced a negative Q.

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

* **The formula interface mistranslated left- and interval-censored
  `Surv()` objects.** `survival::Surv()` and this package use different
  integer codings for censoring status, and the parser passed `Surv()`'s
  through unchanged. `Surv(time, event, type = "left")` codes a left-censored
  row as `0`, which this package reads as *right*-censored: a wrong answer
  with no error, warning, or other outward sign. Under
  `type = "interval"` / `"interval2"`, `Surv()` codes rows `0`/`1`/`2`/`3`
  for right / event / left / interval against this package's `0`/`1`/`-1`/`2`,
  so left-censored rows were read as interval-censored and interval rows
  carried a status the likelihood does not recognise at all.

  Two related faults in the same branch: `Surv()` stores the status in its
  `time2` column for every non-interval row, and the parser read that
  sentinel as an upper bound; and it set `time_lower` for every row, which
  the likelihood treats as a counting-process *entry* time when status is
  `0` or `1`, cancelling each exact-event and right-censored row out of the
  likelihood. Together these made an interval-censored formula fit return
  the optimizer's failure sentinel rather than a fit.

  Status codes are now translated, an upper bound is taken only from a
  genuine interval row, and `time_lower` left-truncates only interval rows.
  A regression test asserts that a `Surv(type = "interval")` fit reproduces
  the equivalent vector-interface fit to 1e-8 in log-likelihood.
  Found when a production study's three interval-censored records could only
  be expressed through the vector interface.

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
  from one that deliberately declines. Behavior is unchanged -- the
  diagnostics are still `NA` -- but the reason is now stated. Found while
  fitting a production interval-censored study.

* **The score criterion could not test a single candidate on an interval- or
  left-censored multiphase fit.** The analytic multiphase Hessian declines by
  design for `status` in `{-1, 2}`, and the score path had no fallback on that
  branch -- the single-distribution branch has had one all along. The `NULL`
  propagated into the step's reusable nuisance block, every candidate scored
  `NA`, and `hzr_stepwise()` stopped having tested nothing, reporting it in the
  language of a degenerate candidate. Both halves became reachable in this
  release and only together: the `Surv()` translation fix above made left- and
  interval-censored rows expressible through the formula interface, and
  `criterion = "score"` became the default. No test exercised the two at once.

  The observed information is now computed numerically where the analytic form
  declines, as the single-distribution path already did. It agrees with the
  analytic Hessian to 1e-4 on the equivalent right-censored fit, which is what
  licenses using it in place of one. The cost is a numeric Hessian per
  candidate -- the per-candidate work the score criterion exists to avoid --
  but it is paid only where there would otherwise be no information matrix at
  all, and slower is the right trade against selecting nothing. `numDeriv` is a
  `Suggests` here as elsewhere: when it is absent this now stops and names both
  it and `criterion = "wald"`, rather than returning a screen that tested
  nothing.


* **`hzr_stepwise(scope = NULL)` still failed on a formula passed by
  variable.** The fix for that defect reached `.hzr_refit_with_scope()` but
  not three sibling sites, so the default-scope path still raised
  `invalid formula "f": not a call` -- the very string the entry below says
  no longer occurs. All four sites now resolve the stored formula through one
  internal helper, so a fifth cannot drift: `match.call()` records `formula`
  unevaluated, and `deparse(quote(f))` is `"f"`, which `as.formula()` rejects.

  Two consequences of that path becoming reachable, both fixed here.
  `scope = NULL` now skips columns it cannot model instead of erroring on
  them -- numeric and logical columns are kept, since whether a 0/1 field
  arrives logical or numeric depends on the reader that built the frame
  rather than on the variable:
  under an explicit scope the caller named the column, so an error is right,
  but under `scope = NULL` the package enumerates the candidates itself and a
  column it cannot model is its own choice to make better. Any data frame
  carrying a character or factor column -- which is most of them -- was
  otherwise unusable with the default scope.

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

* **A stepwise screen that could not score anything now says so, instead of
  looking like one that finished.** Under `criterion = "score"` a candidate
  whose Q statistic cannot be computed -- a degenerate or collinear column,
  or an information matrix that will not invert on this data -- yields `NA`
  and is dropped from consideration. When that happened to every remaining
  candidate the step returned exactly what a legitimate "no candidate met
  `slentry`" stop returns, so a screen that stopped because it was *unable
  to test* its candidates was indistinguishable from one that tested them
  and found nothing. The per-step diagnostic existed on the returned object
  the whole time and had no readers.

  `hzr_stepwise()` now warns when a run stops this way and reports
  `$criteria$n_uncomputable_scores`. Because `hzr_bootstrap()` runs each
  replicate under `suppressWarnings()` -- deliberately, so per-replicate
  numerical noise does not swamp the console -- that warning cannot surface
  in the mode where it matters most, so the count is aggregated instead:
  `hzr_bootstrap()` gains `$n_uncomputable_replicates` and warns once when it
  is non-zero. A replicate that scored nothing still counts toward
  `n_success` while contributing no selections, so it silently depresses
  every reported selection frequency -- which is the whole deliverable of a
  bootstrap screen.

  Found by a pre-release review pass, not by a failing test: the package's
  own `print.hzr_bootstrap` test runs a five-replicate screen in which four
  replicates cannot score a candidate and none selects anything, and it
  passed throughout because it only ever asserted the printed label.

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
