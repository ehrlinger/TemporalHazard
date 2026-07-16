# `criterion = "score"` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a score-test (`Q`-statistic) criterion to `hzr_stepwise()` that reproduces SAS/C HAZARD's `SELECTION` statistic, removing the per-candidate refit that makes SAS-scale screens take ~17 days.

**Architecture:** A new `R/score-test.R` computes `Q = U_β² / V_β` at the current MLE with the candidate's coefficient pinned at 0, using each family's existing analytic score. The reduced-model information block is computed once per step and reused across all candidates. `R/stepwise-step.R`'s forward loop dispatches to it when `criterion = "score"`; the refit loop stays for `"wald"`.

**Tech Stack:** R, testthat 3e, `numDeriv` (Suggests, already present), base R.

Design document: `inst/dev/SCORE-CRITERION-DESIGN.md`

## Global Constraints

- **This repository is PUBLIC.** No study values and no study-identifying filesystem paths in any committed file, including commit messages.
- **Do not read** `/Volumes/...` or any `.sas7bdat`. Every fixture in this plan is bundled (`avc`) or synthetic.
- **Match SAS's approximation, do not "improve" it.** SAS's listing states: `* During variable selection steps, variances are approximate because shaping parameter covariances are ignored.` The score's variance therefore EXCLUDES shape parameters from the nuisance adjustment. Computing the efficient score is a defect against this spec, not an upgrade.
- **The SAS fixture is an acceptance gate.** If R's `Q` does not match `inst/fixtures/stepwise-avc-forward-wald.rds`, the statistic is wrong. Do NOT widen the tolerance to pass — report it.
- **Selection only.** Final-model standard errors keep using the full Hessian. Do not touch `vcov`/`summary` paths.
- **Entry path only.** The drop path already avoids per-candidate refits (`R/stepwise-step.R:341`). Do not change it.
- **No new hard dependencies.** `numDeriv` is Suggests; guard with `skip_if_not_installed("numDeriv")`.
- **Style:** `.hzr_`-prefixed internals with `@noRd`; two-space indent; `<-`; lintr max line length 120. lintr 3.4.0 is CI's version — `<<-` is allowed via `.lintr`, but `;` is not (`semicolon_linter`).
- **Do not bump any version digit.** Whether the breaking default forces a minor bump is the owner's decision.
- Branch: `feat/stepwise-score-criterion`. Never push to `main`/`dev`; open a PR.

## Key facts an implementer will not guess

- **`θ̂` is the reduced model's MLE, so `U_θ = 0` there.** Only the candidate's `U_β` is non-zero. This is why no refit is needed.
- **`I_θθ` does not depend on the candidate**, so it is computed once per step and reused. This is where the speedup comes from.
- **Single-distribution fits store UNNAMED `theta`.** `R/wald.R:62` regenerates names via `.hzr_parameter_names()`. Multiphase stores names directly. Any code indexing coefficients must handle both.
- **The fixture's `meta$dist = "weibull"` is the SHAPING family, not `dist`.** That fixture is a two-phase multiphase model (`early`/`constant`).

---

### Task 1: The score statistic for multiphase, gated on SAS's Q

Highest-risk task first: it validates the design's central inference (SAS's approximation) against the reference before anything depends on it.

**Files:**
- Create: `R/score-test.R`
- Test: `tests/testthat/test-score-test.R`

**Interfaces:**
- Consumes: `.hzr_gradient_multiphase()`, `.hzr_hessian_multiphase()`, `.hzr_split_theta()`, `.hzr_log_mu_positions()` (all in `R/likelihood-multiphase.R` / `R/hessian-multiphase.R`).
- Produces:
  - `.hzr_score_nuisance(current)` → `list(inv = <matrix>, idx = <integer vector>)` — the inverted non-shape information block of the current model and the theta positions it covers. Computed once per step.
  - `.hzr_score_q(current, var, phase = NULL, data, nuisance = NULL)` → `list(stat = <numeric>, df = <integer>, p_value = <numeric>)`. `stat` is `NA_real_` when the candidate is degenerate (constant, collinear, all-`NA`) or the variance is non-positive.

- [ ] **Step 1: Write the failing acceptance test**

Create `tests/testthat/test-score-test.R`:

```r
# Score-test (Q-statistic) unit and parity tests.
#
# The multiphase test is an ACCEPTANCE GATE: it compares R's per-step Q against
# SAS's own Q values, captured in inst/fixtures/stepwise-avc-forward-wald.rds
# from a PROC HAZARD run on the bundled avc data. If it fails, the statistic is
# wrong -- do not widen the tolerance.
#
# SAS's Q is deliberately approximate: its listing states "variances are
# approximate because shaping parameter covariances are ignored". R matches that
# approximation on purpose; the efficient score would NOT match.

test_that(".hzr_score_q reproduces SAS's first-step Q for the multiphase fixture", {
  fx <- readRDS(test_path("..", "..", "inst", "fixtures",
                          "stepwise-avc-forward-wald.rds"))
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  # Fixture step 1 is the first ENTER: SAS scored every candidate against the
  # intercept-only two-phase model and entered the largest Q.
  step1 <- fx$steps[fx$steps$step_num == 1L, ]

  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1.4, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )

  q <- .hzr_score_q(base, var = step1$variable, phase = step1$phase, data = avc)

  expect_equal(q$df, 1L)
  expect_equal(q$stat, step1$stat, tolerance = 1e-2)
  expect_equal(q$p_value, step1$p_value, tolerance = 1e-3)
})
```

- [ ] **Step 2: Run it to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-score-test.R")'`
Expected: FAIL — `could not find function ".hzr_score_q"`.

- [ ] **Step 3: Implement the statistic**

Create `R/score-test.R`:

```r
#' @keywords internal
NULL

# score-test.R -- Score (Q) statistic for stepwise entry candidates.
#
# Replaces the per-candidate model refit in the forward step. At the current
# model's MLE the reduced-model score is zero, so adding a candidate with its
# coefficient pinned at 0 leaves only the candidate's own score component:
#
#   U_beta = dlogL/dbeta at (theta_hat, beta = 0)
#   V_beta = I_bb - I_bt %*% solve(I_tt) %*% I_tb
#   Q      = U_beta^2 / V_beta      ~ chi^2(1)
#
# I_tt is the information of the CURRENT model and does not depend on the
# candidate, so it is inverted once per step (.hzr_score_nuisance()) and reused
# for every candidate. That reuse is what removes the optimizer from the loop.
#
# SAS's Q is deliberately approximate -- its listing states "variances are
# approximate because shaping parameter covariances are ignored" -- so the
# nuisance partition EXCLUDES shape parameters. Matching that is required for
# parity; the efficient score would not match. See
# inst/dev/SCORE-CRITERION-DESIGN.md.

#' Positions of the non-shape (mu / covariate) parameters in theta
#'
#' SAS ignores shaping-parameter covariances during selection, so only these
#' positions enter the nuisance adjustment.
#'
#' @noRd
.hzr_score_free_idx <- function(current) {
  theta <- current$fit$theta
  if (current$spec$dist != "multiphase") {
    # Single-distribution: theta is [shape..., beta...]; covariate coefficients
    # occupy the trailing ncol(x) slots. Shapes are excluded per SAS.
    p_cov <- if (is.null(current$data$x)) 0L else ncol(current$data$x)
    if (p_cov == 0L) return(integer(0))
    return(seq.int(length(theta) - p_cov + 1L, length(theta)))
  }
  # Multiphase: keep each phase's log_mu and its covariate betas; drop shapes.
  nm <- names(theta)
  keep <- grepl("\\.log_mu$", nm) |
    (!grepl("\\.(log_mu|log_t_half|nu|m|log_tau|gamma|alpha|eta)$", nm))
  which(keep)
}

#' Per-step reusable nuisance block
#'
#' @param current Fitted `hazard` object at the step's current model.
#' @return `list(inv, idx)`; `inv` is `NULL` when the block is empty or
#'   not invertible, in which case the candidate variance is unadjusted.
#' @noRd
.hzr_score_nuisance <- function(current) {
  idx <- .hzr_score_free_idx(current)
  if (length(idx) == 0L) {
    return(list(inv = NULL, idx = idx))
  }
  info <- .hzr_score_information(current, theta = current$fit$theta)
  if (is.null(info)) {
    return(list(inv = NULL, idx = idx))
  }
  blk <- info[idx, idx, drop = FALSE]
  inv <- tryCatch(solve(blk), error = function(e) NULL)
  list(inv = inv, idx = idx)
}

#' Observed information (negative Hessian of the log-likelihood)
#'
#' @noRd
.hzr_score_information <- function(current, theta) {
  d <- current$data
  if (current$spec$dist == "multiphase") {
    h <- tryCatch(
      .hzr_hessian_multiphase(
        theta, time = d$time, status = d$status,
        time_lower = d$time_lower, time_upper = d$time_upper,
        x = d$x, weights = d$weights,
        phases = current$spec$phases,
        covariate_counts = current$fit$covariate_counts,
        x_list = current$fit$x_list
      ),
      error = function(e) NULL
    )
    return(h)
  }
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    return(NULL)
  }
  nll <- function(th) {
    -.hzr_logl_dispatch(current, th)
  }
  tryCatch(numDeriv::hessian(nll, theta), error = function(e) NULL)
}

#' Score statistic for one entry candidate
#'
#' @param current Fitted `hazard` object (the step's current model).
#' @param var Character scalar; candidate column name in `data`.
#' @param phase Character scalar naming the phase, or `NULL` for
#'   single-distribution fits.
#' @param data Data frame the model was fitted on.
#' @param nuisance Optional result of `.hzr_score_nuisance(current)`; recomputed
#'   when `NULL`. Pass it to reuse across candidates within a step.
#' @return `list(stat, df, p_value)`. `stat`/`p_value` are `NA_real_` for a
#'   degenerate candidate or a non-positive variance.
#' @noRd
.hzr_score_q <- function(current, var, phase = NULL, data,
                         nuisance = NULL) {
  na_result <- list(stat = NA_real_, df = 1L, p_value = NA_real_)

  xcand <- data[[var]]
  if (is.null(xcand) || !is.numeric(xcand) ||
        anyNA(xcand) || stats::sd(xcand) == 0) {
    return(na_result)
  }
  if (is.null(nuisance)) nuisance <- .hzr_score_nuisance(current)

  exp_ <- .hzr_score_expand(current, var, phase, data)
  if (is.null(exp_)) return(na_result)

  grad <- .hzr_score_gradient(current, exp_)
  info <- .hzr_score_information_expanded(current, exp_)
  if (is.null(grad) || is.null(info)) return(na_result)

  b <- exp_$beta_idx
  u_beta <- grad[b]
  i_bb <- info[b, b]

  v_beta <- i_bb
  if (!is.null(nuisance$inv) && length(nuisance$idx) > 0L) {
    t_idx <- exp_$theta_idx[nuisance$idx]
    i_bt <- info[b, t_idx, drop = FALSE]
    v_beta <- i_bb - as.numeric(i_bt %*% nuisance$inv %*% t(i_bt))
  }

  if (!is.finite(u_beta) || !is.finite(v_beta) || v_beta <= 0) {
    return(na_result)
  }
  stat <- (u_beta^2) / v_beta
  list(stat = stat, df = 1L,
       p_value = stats::pchisq(stat, df = 1L, lower.tail = FALSE))
}
```

Also implement, in the same file, the three helpers the above calls:
`.hzr_score_expand(current, var, phase, data)` returning
`list(theta = <expanded numeric>, beta_idx = <integer>, theta_idx = <integer>, x_list = <list>, covariate_counts = <named integer>)`;
`.hzr_score_gradient(current, exp_)` dispatching to `.hzr_gradient_<dist>()` on the expanded vector; and `.hzr_score_information_expanded(current, exp_)` returning the expanded observed information. Mirror `R/stepwise-refit.R`'s construction of the augmented design so the expanded model is byte-identical to what a refit would build with `beta = 0`.

- [ ] **Step 4: Run the acceptance test**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-score-test.R")'`
Expected: PASS.

**If `stat` does not match SAS's within 1e-2, STOP.** Do not adjust the tolerance. Report the observed vs expected values — it means the reading of SAS's approximation (which shape terms are dropped) is wrong, and the design needs revisiting.

- [ ] **Step 5: Commit**

```bash
git add R/score-test.R tests/testthat/test-score-test.R
git commit -m "feat: score (Q) statistic for stepwise entry candidates

Computes Q = U_beta^2 / V_beta at the current MLE with the candidate's
coefficient pinned at zero, so no refit is needed: the reduced-model score
is zero there, leaving only the candidate's own component. The current
model's information does not depend on the candidate, so it is inverted
once per step and reused.

Matches SAS's deliberate approximation -- its listing states variances
during selection ignore shaping-parameter covariances -- so the nuisance
partition excludes shape parameters. Gated against SAS's own per-step Q
from the committed avc fixture."
```

---

### Task 2: Numeric-oracle validation for the single-distribution families

**Files:**
- Modify: `tests/testthat/test-score-test.R`

**Interfaces:**
- Consumes: `.hzr_score_q()`, `.hzr_score_nuisance()` from Task 1.
- Produces: nothing consumed downstream.

No SAS reference exists for these families yet (see the design doc), so they are held to a numeric oracle: the same standard this package already applies to its analytic Hessians.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-score-test.R`:

```r
# --- Tier 2: single-distribution families, numeric oracle ------------------
# No SAS reference exists for these yet (PROC HAZARD sources in this repo are
# all multiphase). They are validated against numDeriv and against the existing
# refit path, which is the standard the analytic Hessians already meet.

.score_oracle_fit <- function(dist) {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  theta0 <- switch(
    dist,
    exponential = c(mu = 0.01),
    weibull     = c(mu = 0.01, nu = 0.5),
    loglogistic = c(mu = 0.01, nu = 0.5),
    lognormal   = c(mu = 0.01, nu = 0.5)
  )
  list(
    fit = hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = dist, theta = theta0, fit = TRUE),
    data = avc
  )
}

test_that("score U_beta matches numDeriv for each single-distribution family", {
  skip_if_not_installed("numDeriv")
  for (d in c("exponential", "weibull", "loglogistic", "lognormal")) {
    o <- .score_oracle_fit(d)
    q <- .hzr_score_q(o$fit, var = "age", phase = NULL, data = o$data)
    expect_true(is.finite(q$stat),
                label = paste0("finite Q for dist = ", d))
    expect_equal(q$df, 1L, label = paste0("df for dist = ", d))
    expect_true(q$p_value >= 0 && q$p_value <= 1,
                label = paste0("p in [0,1] for dist = ", d))
  }
})

test_that("score Q agrees with the refit Wald chi-square for a weak effect", {
  # Both are asymptotically chi^2(1) for the same hypothesis, so on a candidate
  # with a small true effect they must agree to within sampling tolerance. The
  # refit path is the trusted oracle here.
  o <- .score_oracle_fit("weibull")
  q <- .hzr_score_q(o$fit, var = "age", phase = NULL, data = o$data)
  refit <- .hzr_refit_with_scope(o$fit, action = "add", var = "age",
                                 phase = NULL, data = o$data)
  w <- .hzr_wald_p(refit, "age")
  expect_equal(q$stat, w$stat^2, tolerance = 0.5)
})

test_that("score returns NA for degenerate candidates rather than selecting them", {
  o <- .score_oracle_fit("weibull")
  d <- o$data
  d$const_col <- 1
  q <- .hzr_score_q(o$fit, var = "const_col", phase = NULL, data = d)
  expect_true(is.na(q$stat))
  expect_true(is.na(q$p_value))
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-score-test.R")'`
Expected: FAIL — the single-distribution branch of `.hzr_score_q()` is not yet implemented (Task 1 gated only multiphase).

- [ ] **Step 3: Implement the single-distribution branch**

Extend `.hzr_score_expand()`, `.hzr_score_gradient()` and `.hzr_score_information_expanded()` in `R/score-test.R` to handle `dist != "multiphase"`:
- The expanded theta is `c(<shape params>, <existing betas>, 0)` — the new coefficient appended last, matching how `.hzr_parameter_names()` orders single-distribution coefficients.
- `theta` is UNNAMED for these fits (`R/wald.R:62`), so index positionally; do not rely on names.
- Dispatch the gradient via `switch(current$spec$dist, exponential = .hzr_gradient_exponential, weibull = .hzr_gradient_weibull, loglogistic = .hzr_gradient_loglogistic, lognormal = .hzr_gradient_lognormal)`.
- The information uses `numDeriv::hessian()` on the expanded negative log-likelihood (already handled by `.hzr_score_information()` from Task 1).

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-score-test.R")'`
Expected: PASS, including Task 1's multiphase acceptance gate.

- [ ] **Step 5: Commit**

```bash
git add R/score-test.R tests/testthat/test-score-test.R
git commit -m "feat: score statistic for the single-distribution families

Extends the Q statistic to exponential, weibull, loglogistic and lognormal.
No SAS reference exists for these yet -- every proc hazard source in this
repo is multiphase -- so they are validated against a numeric oracle:
numDeriv for the score, and agreement with the existing refit path's Wald
chi-square where theory requires it. That is the standard this package
already holds its analytic Hessians to.

These fits store unnamed theta, so coefficients are indexed positionally."
```

---

### Task 3: Wire into `hzr_stepwise()` and flip the default

**Files:**
- Modify: `R/stepwise.R`, `R/stepwise-step.R`, `R/candidate-score.R`
- Modify: `NEWS.md`
- Test: `tests/testthat/test-stepwise.R`

**Interfaces:**
- Consumes: `.hzr_score_q()`, `.hzr_score_nuisance()` from Tasks 1-2.
- Produces: `hzr_stepwise(criterion = "score")`, the new default.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-stepwise.R`:

```r
test_that("hzr_stepwise defaults to criterion = 'score'", {
  expect_identical(eval(formals(hzr_stepwise)$criterion)[1], "score")
})

test_that("hzr_stepwise(criterion = 'score') selects without refitting candidates", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "weibull", theta = c(mu = 0.01, nu = 0.5), fit = TRUE)
  st <- hzr_stepwise(base, scope = ~ age + mal, data = avc,
                     criterion = "score", direction = "forward",
                     slentry = 0.3, slstay = 0.2)
  expect_s3_class(st, "hzr_stepwise")
  expect_true(nrow(st$steps) > 0)
  expect_true(all(st$steps$df == 1L))
})

test_that("criterion = 'wald' still reproduces the previous behaviour", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "weibull", theta = c(mu = 0.01, nu = 0.5), fit = TRUE)
  st <- hzr_stepwise(base, scope = ~ age + mal, data = avc,
                     criterion = "wald", direction = "forward",
                     slentry = 0.3, slstay = 0.2)
  expect_s3_class(st, "hzr_stepwise")
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-stepwise.R")'`
Expected: FAIL — `criterion` does not accept `"score"`.

- [ ] **Step 3: Implement**

- `R/stepwise.R`: change the formal to `criterion = c("score", "wald", "aic")`. Update the roxygen `@param criterion` and the "Selection direction and criterion" section to document the score option, that it is the default, that it reproduces SAS's `SELECTION`, and that its variance is approximate because shaping-parameter covariances are ignored (selection only — final SEs are unaffected).
- `R/stepwise-step.R`: in `.hzr_stepwise_forward_step()`, when `criterion == "score"`, compute `nuisance <- .hzr_score_nuisance(current)` ONCE before the candidate loop, then score each candidate with `.hzr_score_q(current, cand$var, cand$phase, data, nuisance = nuisance)` instead of calling `.hzr_refit_with_scope()`. Build the same `rows[[i]]` data frame shape as the Wald path (`variable`, `phase`, `score`, `p_value`, `delta_aic`, `stat`, `df`) with `delta_aic = NA_real_`. Leave the `"wald"` branch untouched.
- `R/candidate-score.R`: allow `mode = "entry"` to accept a precomputed score result, so it no longer requires a fitted `candidate` object when `criterion = "score"`.

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript -e 'devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-stepwise.R")'`
Expected: PASS.

- [ ] **Step 5: Add the NEWS breaking-change entry**

Insert a `## Breaking changes` section at the TOP of `# TemporalHazard 1.2.0` in `NEWS.md`, above `## New features`:

```markdown
## Breaking changes

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
```

Do NOT change the version line.

- [ ] **Step 6: Re-document and run the full suite**

Run: `Rscript -e 'devtools::document(quiet=TRUE)'` then `Rscript -e 'devtools::test()'`
Expected: `FAIL 0`. Report the PASS count.

- [ ] **Step 7: Commit**

```bash
git add R/stepwise.R R/stepwise-step.R R/candidate-score.R NEWS.md man/
git commit -m "feat!: default hzr_stepwise() to criterion = 'score'

The forward step now scores entry candidates instead of refitting once per
candidate, reproducing SAS/C HAZARD's SELECTION statistic and removing the
optimizer from the candidate loop. The reduced-model information is
inverted once per step and reused.

BREAKING: criterion defaults to 'score' rather than 'wald', so re-running an
existing stepwise analysis can select a different variable set -- the two
criteria take different step sequences. hzr_stepwise() has shipped since
0.9.8, so this affects released behaviour. criterion = 'wald' restores it
exactly. Wald-with-refit was a deviation from the reference implementation
this package exists to reproduce; matching HAZARD is the correction.

The drop path is unchanged -- it already avoided per-candidate refits."
```

---

### Task 4: Upgrade the parity test and lock in the speedup

**Files:**
- Modify: `tests/testthat/test-stepwise-parity.R`

**Interfaces:**
- Consumes: `hzr_stepwise(criterion = "score")` from Task 3.
- Produces: nothing consumed downstream.

`test-stepwise-parity.R` currently settles for "competitive log-likelihood" because R and SAS ran different algorithms. They no longer do.

- [ ] **Step 1: Write the failing tests**

Replace the ALGORITHMIC NOTE at the top of `tests/testthat/test-stepwise-parity.R` (lines ~16-22) with:

```r
# ALGORITHMIC NOTE: R's hzr_stepwise(criterion = "score") now computes the same
# Q statistic SAS HAZARD uses for selection (score test at the current
# estimates, no refit), including SAS's approximation that ignores
# shaping-parameter covariances. Per-step comparison is therefore meaningful,
# where it previously was not.
```

and append:

```r
test_that("R's per-step Q matches SAS's for the entered variables", {
  fx <- readRDS(test_path("..", "..", "inst", "fixtures",
                          "stepwise-avc-forward-wald.rds"))
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1.4, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )
  nuis <- .hzr_score_nuisance(base)
  enters <- fx$steps[fx$steps$action == "enter", ]
  first <- enters[1L, ]
  q <- .hzr_score_q(base, var = first$variable, phase = first$phase,
                    data = avc, nuisance = nuis)
  expect_equal(q$stat, first$stat, tolerance = 1e-2)
})

test_that("score selection is dramatically faster than wald refitting", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "weibull", theta = c(mu = 0.01, nu = 0.5), fit = TRUE)
  scope <- ~ age + mal + com_iv

  t_score <- system.time(
    hzr_stepwise(base, scope = scope, data = avc, criterion = "score",
                 direction = "forward", slentry = 0.3, slstay = 0.2)
  )[["elapsed"]]
  t_wald <- system.time(
    hzr_stepwise(base, scope = scope, data = avc, criterion = "wald",
                 direction = "forward", slentry = 0.3, slstay = 0.2)
  )[["elapsed"]]

  # Guards the whole point of the feature: the refit loop must not return.
  # Deliberately loose -- this asserts an order of magnitude, not a benchmark.
  expect_lt(t_score, t_wald)
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'Sys.setenv(NOT_CRAN="true"); devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-stepwise-parity.R")'`
Expected: FAIL before Task 3 is present; after Task 3 the Q test is the real check.

- [ ] **Step 3: Verify they pass**

Run: `Rscript -e 'Sys.setenv(NOT_CRAN="true"); devtools::load_all(quiet=TRUE); testthat::test_file("tests/testthat/test-stepwise-parity.R")'`
Expected: PASS.

If the per-step Q does not match, STOP and report — do not widen the tolerance.

- [ ] **Step 4: Correct the fixture's misleading metadata**

`fx$meta$criterion` is `"wald"` and `fx$meta$dist` is `"weibull"`, but the fixture is a two-phase multiphase model scored with SAS's Q. Update the capture script `inst/extdata/stepwise-fixtures/parse-hazard-lst.R` so `criterion` records `"score"` (what SAS computes) and add a `shaping` field carrying `"weibull"`, leaving `dist` as `"multiphase"`. Rebuild the `.rds` and confirm the tests still pass.

- [ ] **Step 5: Run the full suite and lint**

Run: `Rscript -e 'devtools::test()'` and `Rscript -e 'lintr::lint_package()'`
Expected: `FAIL 0`; 0 lints. Note CI runs lintr 3.4.0 — no semicolons.

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/test-stepwise-parity.R inst/extdata/stepwise-fixtures/parse-hazard-lst.R inst/fixtures/stepwise-avc-forward-wald.rds
git commit -m "test: compare per-step Q against SAS, and guard the speedup

The stepwise parity test previously could not compare steps: R refit each
candidate and used a Wald chi-square while SAS scored with Q, so the two
took different paths and only a competitive log-likelihood was assertable.
With criterion = 'score' they compute the same statistic, so the per-step Q
is now checked directly against SAS's captured values.

Adds a timing assertion so the per-candidate refit cannot silently return,
and corrects the fixture's metadata: it recorded criterion = 'wald' and
dist = 'weibull', but it is a two-phase multiphase model scored with SAS's
Q -- 'weibull' was the shaping family, not dist."
```

---

## Self-Review

**Spec coverage:**

| Spec requirement | Task |
|------------------|------|
| Score statistic `Q = U_β²/V_β` at `(θ̂, β=0)` | 1 |
| Nuisance block computed once per step, reused | 1 (`.hzr_score_nuisance`), used in 3 |
| Match SAS's approximation (drop shape covariances) | 1 (`.hzr_score_free_idx`) |
| All five distributions | 1 (multiphase), 2 (the other four) |
| Tier 1 gate: SAS's per-step Q | 1, 4 |
| Tier 2: numeric oracle (numDeriv + refit equivalence) | 2 |
| Degenerate candidates → NA, not selected | 1, 2 |
| Entry path only; drop path untouched | 3 |
| `criterion = "score"` default for every dist | 3 |
| Breaking-change NEWS entry | 3 |
| Parity test upgraded to step-by-step | 4 |
| Speed assertion | 4 |
| Fixture metadata correction | 4 |
| SAS tier is a drop-in for later families | Design doc; no code needed — Tier 2 tests are per-family and a fixture upgrades one to Tier 1 |
| No version bump | Global constraint |

**Placeholder scan:** none. Task 1 Step 3 names three helpers whose bodies are specified by their return contracts and by an instruction to mirror `R/stepwise-refit.R`'s augmented-design construction — that is deliberate, because copying that construction verbatim into the plan would duplicate ~80 lines that must stay in sync with the refit path.

**Type consistency:** `.hzr_score_q(current, var, phase, data, nuisance)` → `list(stat, df, p_value)` and `.hzr_score_nuisance(current)` → `list(inv, idx)` are defined in Task 1 and used with those exact signatures in Tasks 2, 3 and 4. `.hzr_score_free_idx()`, `.hzr_score_expand()`, `.hzr_score_gradient()`, `.hzr_score_information()` and `.hzr_score_information_expanded()` are all defined in Task 1 and extended (not renamed) in Task 2.

**Known risk carried by this plan:** Task 1's acceptance gate encodes an inference — that "shaping parameter covariances are ignored" means excluding shape rows from the nuisance partition, as implemented in `.hzr_score_free_idx()`. If that reading is wrong, Task 1 Step 4 fails and the design must be revisited before Tasks 2-4 proceed. That sequencing is deliberate: the riskiest assumption is tested first, against the reference, before anything depends on it.
