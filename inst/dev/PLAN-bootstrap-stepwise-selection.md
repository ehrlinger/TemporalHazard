# `scope=` Selection Mode for `hzr_bootstrap()` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `hzr_bootstrap()` an optional `scope` argument that runs a fresh `hzr_stepwise()` selection on each bootstrap replicate (instead of a fixed-formula refit), so `summary$pct` reports variable **selection frequency** — the R equivalent of SAS's `%HAZBOOT`, and the missing half of the `hz.dead.sas` → `bh.dead.sas` workflow.

**Architecture:** One function (`hzr_bootstrap()` in `R/diagnostics.R`) grows a `scope`-gated branch. `scope = NULL` keeps today's behavior byte-for-byte. `scope` supplied swaps the per-replicate body from "refit `object`'s exact formula" to "refit the shape-fixed base on resampled data, then run `hzr_stepwise()` from there." A small helper (`.hzr_bootstrap_param_names()`) is extracted first so both branches share the same parameter-name-resolution logic.

**Tech Stack:** R (package `TemporalHazard`), testthat 3rd edition, roxygen2 8.0.0.

## Global Constraints

- `scope = NULL` behavior must not change in any observable way (existing tests in `tests/testthat/test-diagnostics.R` and `tests/testthat/test-sas-parity.R` must keep passing unmodified).
- No PHI-adjacent data in this repo — all new tests use the shipped `avc` fixture only.
- No new S3 class; `hzr_bootstrap()` keeps returning class `"hzr_bootstrap"`.
- Do not attempt SAS selection-*path* parity (documented, accepted gap — see `inst/dev/BOOTSTRAP-SELECTION-DESIGN.md`).
- Every step that touches `R/` or `tests/` ends with running the affected tests before moving on; do not batch verification to the end.
- Branch `feature/bootstrap-stepwise-selection` is already checked out with the design spec committed (`inst/dev/BOOTSTRAP-SELECTION-DESIGN.md`). Stay on this branch for all commits in this plan.

---

### Task 1: Extract `.hzr_bootstrap_param_names()` helper

**Files:**
- Modify: `R/diagnostics.R:1293-1317` (inline block inside `hzr_bootstrap()`)
- Test: `tests/testthat/test-diagnostics.R` (no new test file — this is a pure refactor verified by the *existing* `hzr_bootstrap` test block at lines 578-776)

**Interfaces:**
- Produces: `.hzr_bootstrap_param_names(fit_obj)` — internal helper, takes any fitted-`hazard`-like object (must have `$fit$theta` and `$data$x`), returns a character vector of resolved parameter names the same length as `fit_obj$fit$theta`. Task 2 calls this once for the fixed-refit path and once per replicate for the select path.

This is a pure extraction: the existing inline logic in `hzr_bootstrap()` (currently computed once, before the loop, from `object`) moves into a standalone function with no behavior change.

- [ ] **Step 1: Add the helper function directly above `hzr_bootstrap()` in `R/diagnostics.R`**

Insert immediately before the `#' Bootstrap resampling for hazard model coefficients` roxygen block (currently line 1191):

```r
#' Resolve parameter names for a bootstrap replicate's fitted theta
#'
#' Shape parameters are already named in `theta`; covariate betas often
#' come through with empty names. Covariate coefficients occupy the last
#' `ncol(x)` positions of theta -- fill any blanks within that block from
#' the design matrix column names by relative index, so downstream pivots
#' (e.g. `reshape(wide)`) get a distinct column per covariate, even when
#' some betas are already named and others are not.
#'
#' @param fit_obj A fitted `hazard`-like object (has `$fit$theta` and
#'   `$data$x`).
#' @return Character vector of resolved parameter names, same length as
#'   `fit_obj$fit$theta`.
#' @keywords internal
#' @noRd
.hzr_bootstrap_param_names <- function(fit_obj) {
  theta <- fit_obj$fit$theta
  param_names <- names(theta)
  if (is.null(param_names)) {
    param_names <- character(length(theta))
  }
  if (!is.null(fit_obj$data$x)) {
    x_names <- colnames(fit_obj$data$x)
    p <- ncol(fit_obj$data$x)
    n_theta <- length(param_names)
    if (!is.null(x_names) && p > 0L && n_theta >= p) {
      cov_idx <- seq.int(n_theta - p + 1L, n_theta)
      blank_in_block <- !nzchar(param_names[cov_idx])
      param_names[cov_idx[blank_in_block]] <- x_names[blank_in_block]
    }
  }
  still_blank <- !nzchar(param_names)
  if (any(still_blank)) {
    param_names[still_blank] <- paste0("param_", which(still_blank))
  }
  param_names
}

```

- [ ] **Step 2: Replace the inline block inside `hzr_bootstrap()` with a call to the helper**

Find this block (currently `R/diagnostics.R:1293-1317`, inside `hzr_bootstrap()`, right after the `orig_weights` computation):

```r
  # Parameter names from the fitted model. Shape parameters (e.g. mu, nu) are
  # named in theta, but covariate betas often come through with empty names.
  # Covariate coefficients occupy the last ncol(x) positions of theta; fill
  # any blanks within that block from the design matrix column names by
  # relative index, so downstream pivots (e.g. reshape(wide)) get a distinct
  # column per covariate -- even when some betas are already named and
  # others are not.
  param_names <- names(object$fit$theta)
  if (is.null(param_names)) {
    param_names <- character(length(object$fit$theta))
  }
  if (!is.null(object$data$x)) {
    x_names <- colnames(object$data$x)
    p <- ncol(object$data$x)
    n_theta <- length(param_names)
    if (!is.null(x_names) && p > 0L && n_theta >= p) {
      cov_idx <- seq.int(n_theta - p + 1L, n_theta)
      blank_in_block <- !nzchar(param_names[cov_idx])
      param_names[cov_idx[blank_in_block]] <- x_names[blank_in_block]
    }
  }
  still_blank <- !nzchar(param_names)
  if (any(still_blank)) {
    param_names[still_blank] <- paste0("param_", which(still_blank))
  }
```

Replace it with:

```r
  param_names <- .hzr_bootstrap_param_names(object)
```

- [ ] **Step 3: Run the existing `hzr_bootstrap` test block to confirm no regression**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-diagnostics.R')"`
Expected: All tests pass, including the `hzr_bootstrap` block (lines 578-776) — `hzr_bootstrap returns correct structure`, `... seed gives reproducible results`, `... labels covariate rows with design-matrix names`, `... preserves named betas while filling blanks`, `... resamples weights alongside the data`, and the two out-of-scope-symbol tests. Zero failures.

- [ ] **Step 4: Commit**

```bash
git add R/diagnostics.R
git commit -m "refactor: extract .hzr_bootstrap_param_names() helper from hzr_bootstrap()

Pure extraction, no behavior change -- prepares for the scope= select
mode in the next commit, which needs the same name-resolution logic
applied per-replicate instead of once up front."
```

---

### Task 2: Add `scope=` selection mode to `hzr_bootstrap()`

**Files:**
- Modify: `R/diagnostics.R` (the `hzr_bootstrap()` function body and its roxygen block, plus `print.hzr_bootstrap()`)
- Test: `tests/testthat/test-diagnostics.R` (new tests inserted after line 776, before the `hzr_competing_risks tests` section header)

**Interfaces:**
- Consumes: `.hzr_bootstrap_param_names(fit_obj)` from Task 1. `hzr_stepwise(fit, scope, data, direction, criterion, slentry, slstay, max_steps, max_move, force_in, force_out, trace, ...)` from `R/stepwise.R:151` (already exported, signature unchanged).
- Produces: `hzr_bootstrap()`'s new signature (below); result list gains `mode` (`"refit"`/`"select"`) and, in select mode, `scope`. `print.hzr_bootstrap()` gains a "Mode:" line. Task 3 and the neo_therapy worked example both call this new signature directly.

New signature:

```r
hzr_bootstrap(object, n_boot = 200L, fraction = 1.0, seed = NULL,
              verbose = FALSE, scope = NULL,
              direction = c("both", "forward", "backward"),
              criterion = c("wald", "aic"), slentry = 0.30, slstay = 0.20,
              max_steps = 50L, max_move = 4L, force_in = character(),
              force_out = character(), ...)
```

- [ ] **Step 1: Write the failing tests**

Insert into `tests/testthat/test-diagnostics.R`, immediately after the last existing `hzr_bootstrap` test (`hzr_bootstrap uses stored data frame when the original symbol is out of scope`, ending at line 776) and before the `# hzr_competing_risks tests` section header (line 779):

```r
test_that("hzr_bootstrap scope = NULL reports mode = refit and is unaffected", {
  set.seed(42)
  fit <- .fit_avc_weibull()
  bs <- hzr_bootstrap(fit, n_boot = 10, seed = 123)

  expect_identical(bs$mode, "refit")
  expect_null(bs$scope)
  expect_true(all(bs$summary$pct > 0))
})

test_that("hzr_bootstrap scope runs embedded stepwise selection per replicate", {
  set.seed(13)
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  bs <- hzr_bootstrap(base, n_boot = 10, seed = 321,
                       scope = ~ age + mal + com_iv,
                       slentry = 0.3, slstay = 0.2,
                       control = list(n_starts = 1))

  expect_s3_class(bs, "hzr_bootstrap")
  expect_identical(bs$mode, "select")
  expect_true(is.list(bs$scope) || inherits(bs$scope, "formula"))
  expect_gt(bs$n_success, 0)
  expect_true(all(bs$summary$pct > 0 & bs$summary$pct <= 100))
  # Shape parameters are never dropped by stepwise selection (only
  # covariates are), so they must be present -- and selected -- in every
  # successful replicate.
  shape_rows <- bs$summary[bs$summary$parameter %in% c("mu", "nu"), ]
  expect_true(all(shape_rows$pct == 100))
})

test_that("hzr_bootstrap scope raises immediately on a structurally invalid scope", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  # scope must be NULL, a one-sided formula, or a character vector for a
  # single-distribution fit (hzr_stepwise()'s own validation); this must
  # raise via the pre-loop validation call, not after n_boot replicates.
  expect_error(
    hzr_bootstrap(base, n_boot = 10, scope = 42),
    "one-sided formula"
  )
})

test_that("hzr_bootstrap scope with a nonexistent column warns but does not raise", {
  # A typo'd (but syntactically valid) scope column is NOT caught by the
  # pre-loop validation -- hzr_stepwise() itself only warns per candidate
  # refit failure (see stepwise-step.R .hzr_stepwise_forward_step) and
  # never selects it. Confirms hzr_bootstrap() inherits that behavior
  # rather than silently reporting it as n_failed replicates.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )

  expect_warning(
    bs <- hzr_bootstrap(base, n_boot = 3, seed = 1,
                         scope = ~ age + not_a_real_column,
                         control = list(n_starts = 1)),
    "not_a_real_column"
  )
  expect_false("not_a_real_column" %in% bs$summary$parameter)
})

test_that("print.hzr_bootstrap reports the mode", {
  fit <- .fit_avc_weibull()
  bs_refit <- hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_output(print(bs_refit), "fixed refit")

  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5),
    fit   = TRUE
  )
  bs_sel <- hzr_bootstrap(base, n_boot = 5, seed = 42, scope = ~ age + mal,
                           control = list(n_starts = 1))
  expect_output(print(bs_sel), "stepwise selection")
})
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-diagnostics.R')"`
Expected: FAIL — `hzr_bootstrap()` does not accept a `scope` argument yet (`unused argument (scope = ...)` or similar), so all five new tests fail. The pre-existing tests above them still pass.

- [ ] **Step 3: Replace `hzr_bootstrap()` and `print.hzr_bootstrap()` with the extended implementation**

Replace the entire current `hzr_bootstrap()` function (roxygen block + body, `R/diagnostics.R` from `#' Bootstrap resampling for hazard model coefficients` through the closing `}` of the function — after Task 1, roughly lines 1191-1360-ish) with:

```r
#' Bootstrap resampling for hazard model coefficients
#'
#' Resample data with replacement, refit the hazard model on each
#' replicate, and accumulate coefficient distributions. Returns a tidy
#' data frame of per-replicate estimates with summary statistics.
#' This is the R equivalent of the SAS `bootstrap.hazard.sas` macro.
#'
#' When `scope` is supplied, each replicate instead runs a fresh
#' [hzr_stepwise()] selection on the resampled data (starting from a
#' fixed-shape refit of `object`) instead of refitting `object`'s exact
#' formula. This is the R equivalent of the SAS `%HAZBOOT` macro: fit the
#' hazard shape with no covariates (fixing it via `hzr_phase(..., fixed =
#' "shapes")`), then bootstrap-screen candidate covariates for how often
#' they enter the model. `summary$pct` then reports the selection
#' frequency across replicates, and `summary$mean`/`sd`/`ci_*` describe the
#' coefficient distribution conditional on selection.
#'
#' @param object A fitted `hazard` object (with `fit = TRUE`).
#' @param n_boot Integer: number of bootstrap replicates (default 200).
#' @param fraction Numeric in (0, 1]: fraction of data to sample per
#'   replicate (default 1.0 for full bootstrap; < 1 for bagging).
#' @param seed Optional integer random seed for reproducibility. When
#'   supplied, `set.seed(seed)` is called at function entry, jumping the
#'   global RNG to the seeded state; it is not restored on exit. Pass
#'   `NULL` (the default) to skip the `set.seed()` call and start from
#'   the caller's current RNG state. Note that the bootstrap consumes
#'   random numbers either way, so the global RNG state will advance
#'   during the call -- `seed = NULL` avoids the *reset* at entry, not
#'   the advance during resampling.
#' @param verbose Logical; if `TRUE`, print progress every 50 replicates.
#' @param scope Candidate variable scope for embedded stepwise selection
#'   during each bootstrap replicate. `NULL` (default) preserves the
#'   original fixed-formula bootstrap: every replicate refits `object`'s
#'   exact model, and `summary$pct` is always ~100. When supplied (a
#'   one-sided formula, character vector, or -- for multiphase fits -- a
#'   named list of one-sided formulas keyed by phase, matching
#'   [hzr_stepwise()]'s `scope`), each replicate runs a fresh
#'   [hzr_stepwise()] selection instead; see Details.
#' @param direction,criterion,slentry,slstay,max_steps,max_move,force_in,force_out
#'   Passed through to [hzr_stepwise()] on each replicate when `scope` is
#'   supplied; ignored when `scope = NULL`. See [hzr_stepwise()] for
#'   definitions and defaults.
#' @param ... Additional arguments forwarded to [hzr_stepwise()] (e.g.
#'   `control = list(n_starts = 1)`) when `scope` is supplied; ignored
#'   otherwise.
#'
#' @return A list with class `"hzr_bootstrap"` containing:
#' \describe{
#'   \item{replicates}{Data frame with columns `replicate`, `parameter`,
#'     and `estimate` -- one row per parameter per successful replicate.}
#'   \item{summary}{Data frame with columns `parameter`, `n`, `pct`,
#'     `mean`, `sd`, `min`, `max`, `ci_lower`, `ci_upper` -- one row per
#'     parameter. In `mode = "select"`, `pct` is the selection frequency
#'     and the other statistics are conditional on selection.}
#'   \item{n_success}{Number of successfully converged replicates.}
#'   \item{n_failed}{Number of replicates that failed to converge.}
#'   \item{mode}{`"refit"` (fixed-formula bootstrap) or `"select"`
#'     (embedded stepwise selection).}
#'   \item{scope}{Only present when `mode == "select"`: the candidate
#'     scope used.}
#' }
#'
#' @examples
#' \donttest{
#' data(avc)
#' avc <- na.omit(avc)
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ age + mal,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5, 0, 0),
#'   fit   = TRUE
#' )
#' bs <- hzr_bootstrap(fit, n_boot = 50, seed = 123)
#' print(bs)
#'
#' # Embedded stepwise selection: screen candidate covariates for how
#' # often they enter the model across resamples (R equivalent of SAS
#' # %HAZBOOT).
#' base <- hazard(
#'   survival::Surv(int_dead, dead) ~ 1,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5),
#'   fit   = TRUE
#' )
#' bs_sel <- hzr_bootstrap(base, n_boot = 20, seed = 123,
#'                          scope = ~ age + mal,
#'                          slentry = 0.3, slstay = 0.2)
#' print(bs_sel)
#' }
#'
#' @seealso [hazard()] for model fitting, [vcov.hazard()] for
#'   Hessian-based standard errors, [hzr_stepwise()] for the selection
#'   procedure used when `scope` is supplied.
#' @export
hzr_bootstrap <- function(object, n_boot = 200L, fraction = 1.0,
                           seed = NULL, verbose = FALSE,
                           scope = NULL,
                           direction = c("both", "forward", "backward"),
                           criterion = c("wald", "aic"),
                           slentry = 0.30, slstay = 0.20,
                           max_steps = 50L, max_move = 4L,
                           force_in = character(), force_out = character(),
                           ...) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a fitted hazard object.", call. = FALSE)
  }
  if (is.null(object$fit$theta) ||
      (is.logical(object$fit$converged) && is.na(object$fit$converged))) {
    stop("'object' has no fitted parameters. Refit with fit = TRUE.",
         call. = FALSE)
  }

  n_boot <- as.integer(n_boot)
  if (n_boot < 1L) stop("'n_boot' must be at least 1.", call. = FALSE)
  if (fraction <= 0 || fraction > 1) {
    stop("'fraction' must be in (0, 1].", call. = FALSE)
  }

  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  select_mode <- !is.null(scope)

  if (!is.null(seed)) set.seed(seed)

  # Reconstruct the call components
  cl <- object$call
  # Prefer the evaluated `data` argument stored on the fitted object
  # (`object$data$frame` -- the data frame passed to hazard(), not a
  # model.frame() result): it is guaranteed available and needs no caller-frame
  # lookup. Re-evaluating `cl$data` in `parent.frame()` is fragile -- the
  # original `data` symbol may no longer be in scope (e.g. the fit was built
  # inside a helper that has returned) -- so fall back to it only for objects
  # fitted before the frame was stored on `object$data`.
  orig_data <- if (!is.null(object$data$frame)) {
    object$data$frame
  } else {
    eval(cl$data, envir = parent.frame())
  }
  n_obs <- nrow(orig_data)
  sample_size <- max(1L, as.integer(n_obs * fraction))

  # Observation weights, if any, must be resampled in lockstep with the data.
  # Prefer the weights already evaluated and stored on the fitted object: they
  # are guaranteed aligned with the original rows and need no caller-frame
  # lookup. Re-evaluating `cl$weights` in `parent.frame()` is fragile -- the
  # original symbol/expression may no longer be in scope, or may now resolve to
  # a different value -- so fall back to it only for objects fitted before
  # weights were stored on `object$data`. Each replicate is then rewired to a
  # locally bound, resampled copy (mirroring `data`).
  orig_weights <- if (!is.null(object$data$weights)) {
    object$data$weights
  } else if (!is.null(cl$weights)) {
    eval(cl$weights, envir = parent.frame())
  } else {
    NULL
  }

  # Select-mode: fail loud on a structurally invalid scope (wrong type, an
  # unnamed list for a multiphase fit, an unknown phase name) up front, by
  # running one stepwise search against the untouched data, instead of
  # surfacing it n_boot replicates later. This does NOT catch a merely
  # nonexistent column name -- hzr_stepwise() itself converts that into a
  # per-candidate warning() rather than an error (see
  # inst/dev/BOOTSTRAP-SELECTION-DESIGN.md, "Error handling").
  if (select_mode) {
    hzr_stepwise(
      object, scope = scope, data = orig_data,
      direction = direction, criterion = criterion,
      slentry = slentry, slstay = slstay,
      max_steps = max_steps, max_move = max_move,
      force_in = force_in, force_out = force_out,
      trace = FALSE, ...
    )
  }

  # Parameter names from the fitted model. In fixed-refit mode every
  # replicate shares the same theta layout, so names are resolved once, up
  # front. In select-mode, each replicate can select a different variable
  # set, so names are resolved per replicate inside the loop instead.
  param_names <- if (!select_mode) .hzr_bootstrap_param_names(object) else NULL

  # Accumulate results
  rep_list <- vector("list", n_boot)
  n_success <- 0L
  n_failed <- 0L

  for (b in seq_len(n_boot)) {
    if (verbose && b %% 50 == 0) {
      cat("Bootstrap replicate", b, "/", n_boot, "\n")
    }

    # Resample with replacement
    idx <- sample.int(n_obs, size = sample_size, replace = TRUE)
    boot_data <- orig_data[idx, , drop = FALSE] # nolint: object_usage_linter.
    # boot_weights is referenced via quote() inside eval -- lintr cannot trace it
    boot_weights <- if (is.null(orig_weights)) NULL else orig_weights[idx] # nolint: object_usage_linter.

    if (select_mode) {
      # Refit the (shape-fixed) base model on the resampled data first, so
      # the stepwise search's entry/retention tests compare candidates
      # against a base likelihood computed on the SAME resampled data --
      # then run a fresh stepwise selection from that base.
      boot_fit <- tryCatch({
        cl_base <- cl
        cl_base$data <- quote(boot_data)
        if (!is.null(orig_weights)) cl_base$weights <- quote(boot_weights)
        cl_base$fit <- TRUE
        base_boot <- eval(cl_base)
        if (!is.finite(base_boot$fit$objective)) {
          stop("base refit did not converge")
        }
        hzr_stepwise(
          base_boot, scope = scope, data = boot_data,
          direction = direction, criterion = criterion,
          slentry = slentry, slstay = slstay,
          max_steps = max_steps, max_move = max_move,
          force_in = force_in, force_out = force_out,
          trace = FALSE, ...
        )
      }, error = function(e) NULL)
    } else {
      # Refit using the same call but with resampled data (and weights, if any)
      # (boot_data/boot_weights are referenced via quote() inside eval)
      boot_fit <- tryCatch({
        cl_boot <- cl
        cl_boot$data <- quote(boot_data)
        if (!is.null(orig_weights)) cl_boot$weights <- quote(boot_weights)
        cl_boot$fit <- TRUE
        eval(cl_boot)
      }, error = function(e) NULL)
    }

    if (!is.null(boot_fit) && is.finite(boot_fit$fit$objective)) {
      n_success <- n_success + 1L
      theta_b <- boot_fit$fit$theta
      names_b <- if (select_mode) {
        .hzr_bootstrap_param_names(boot_fit)
      } else {
        param_names[seq_along(theta_b)]
      }
      rep_list[[b]] <- data.frame(
        replicate = b,
        parameter = names_b,
        estimate  = as.numeric(theta_b),
        stringsAsFactors = FALSE
      )
    } else {
      n_failed <- n_failed + 1L
    }
  }

  # Combine replicates
  replicates <- do.call(rbind, Filter(Negate(is.null), rep_list))
  if (is.null(replicates)) {
    replicates <- data.frame(replicate = integer(0),
                              parameter = character(0),
                              estimate = numeric(0),
                              stringsAsFactors = FALSE)
  }
  rownames(replicates) <- NULL

  # Summary statistics per parameter
  if (nrow(replicates) > 0) {
    summary_list <- lapply(split(replicates, replicates$parameter), function(d) {
      data.frame(
        parameter = d$parameter[1],
        n         = nrow(d),
        pct       = 100 * nrow(d) / n_boot,
        mean      = mean(d$estimate),
        sd        = stats::sd(d$estimate),
        min       = min(d$estimate),
        max       = max(d$estimate),
        ci_lower  = stats::quantile(d$estimate, 0.025),
        ci_upper  = stats::quantile(d$estimate, 0.975),
        stringsAsFactors = FALSE
      )
    })
    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL
    if (select_mode) {
      # No single canonical variable order across replicates (different
      # replicates select different sets) -- rank by selection frequency,
      # most-selected first, ties broken alphabetically.
      summary_df <- summary_df[order(-summary_df$pct, summary_df$parameter), ]
    } else {
      # Sort by parameter order in the original model
      idx_order <- match(summary_df$parameter, param_names)
      summary_df <- summary_df[order(idx_order), ]
    }
  } else {
    summary_df <- data.frame(parameter = character(0), n = integer(0),
                              pct = numeric(0), mean = numeric(0),
                              sd = numeric(0), min = numeric(0),
                              max = numeric(0), ci_lower = numeric(0),
                              ci_upper = numeric(0),
                              stringsAsFactors = FALSE)
  }

  result <- list(
    replicates = replicates,
    summary    = summary_df,
    n_success  = n_success,
    n_failed   = n_failed,
    mode       = if (select_mode) "select" else "refit"
  )
  if (select_mode) result$scope <- scope
  class(result) <- "hzr_bootstrap"
  result
}

#' @rdname hzr_bootstrap
#' @param x An `hzr_bootstrap` object.
#' @param digits Number of decimal places for formatting.
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_bootstrap <- function(x, digits = 4, ...) {
  cat("Bootstrap inference for hazard model\n")
  mode <- x$mode %||% "refit"
  cat("Mode:", if (identical(mode, "select")) {
    "embedded stepwise selection"
  } else {
    "fixed refit"
  }, "\n")
  cat("Replicates:", x$n_success, "successful,", x$n_failed, "failed\n\n")
  if (nrow(x$summary) > 0) {
    display <- x$summary
    for (col in c("pct", "mean", "sd", "min", "max",
                   "ci_lower", "ci_upper")) {
      display[[col]] <- round(display[[col]], digits)
    }
    print(display, row.names = FALSE)
  }
  invisible(x)
}
```

Note: `%||%` is already defined package-wide in `R/stepwise-step.R:571` (`` `%||%` <- function(a, b) if (is.null(a)) b else a ``) — no new definition needed.

- [ ] **Step 4: Run the new tests to verify they pass**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-diagnostics.R')"`
Expected: PASS — all 5 new tests plus all pre-existing tests in the file (including the Task 1 regression check) pass. Zero failures.

- [ ] **Step 5: Regenerate documentation**

Run: `Rscript -e "devtools::document()"`
Expected: `man/hzr_bootstrap.Rd` is regenerated with the new `@param` entries and the second `\donttest{}` example; `NAMESPACE` is unchanged (no new exports — `print.hzr_bootstrap` was already an S3 method registration that exists). Confirm with `git diff --stat man/ NAMESPACE` that only `man/hzr_bootstrap.Rd` changed.

- [ ] **Step 6: Run the full test suite**

Run: `Rscript -e "devtools::test()"`
Expected: 0 failures, 0 errors. (Warnings are expected and fine: the `hzr_bootstrap scope with a nonexistent column warns but does not raise` test intentionally triggers `warning()`s.)

- [ ] **Step 7: Commit**

```bash
git add R/diagnostics.R man/hzr_bootstrap.Rd tests/testthat/test-diagnostics.R
git commit -m "feat: add scope= embedded stepwise selection to hzr_bootstrap()

Each replicate can now run a fresh hzr_stepwise() selection (starting
from a fixed-shape refit of the base model) instead of a plain refit,
so summary\$pct becomes the SAS %HAZBOOT-style selection frequency.
scope = NULL (default) is unchanged. Closes the gap documented in
inst/dev/FIXTURE-GAP-LIST.md and the bs.death.AVC test block."
```

---

### Task 3: Exercise the multiphase `scope=` path against real AVC data

**Files:**
- Modify: `tests/testthat/test-sas-parity.R` (insert a new `test_that()` block after the existing `bs.death.AVC: R fixed-model bootstrap runs on the AVC cohort` test, which currently ends at line 730)

**Interfaces:**
- Consumes: `hzr_bootstrap(object, scope = ..., ...)` from Task 2.

This upgrades the "half-written" `bs.death.AVC` parity test (documented at `test-sas-parity.R:667-709`, which only parses SAS reference frequencies) so R's new capability is actually exercised on the shipped, non-PHI `avc` dataset — a multiphase fit with `fixed = "shapes"`, mirroring the `hz.death.AVC` → `bh.death.AVC` SAS hand-off. This test does **not** assert R-vs-SAS frequency parity (still a documented, accepted gap) and does **not** require SAS `.lst` fixtures, so it runs in public CI.

- [ ] **Step 1: Write the test**

Insert into `tests/testthat/test-sas-parity.R`, immediately after the `test_that("bs.death.AVC: R fixed-model bootstrap runs on the AVC cohort", { ... })` block (ends at line 730) and before the `# hp.death.AVC.hm1 / hm2` comment block (line 732):

```r
test_that("bs.death.AVC: R scope= embedded stepwise selection runs on the AVC cohort", {
  testthat::skip_on_cran()
  # No SAS .lst needed here -- exercises R's new scope= bootstrap-selection
  # path (Task 2 of inst/dev/PLAN-bootstrap-stepwise-selection.md) on the
  # shipped avc dataset, so it must run in public CI (do NOT gate on
  # fixtures). This intentionally does NOT assert R-vs-SAS selection
  # frequency parity -- see the "bs.death.AVC" SAS-frequency test above and
  # inst/dev/BOOTSTRAP-SELECTION-DESIGN.md for why that remains a
  # documented, accepted gap.
  set.seed(222)

  data(avc, package = "TemporalHazard")
  # Same null 2-phase AVC model as the fixed-refit test above, with shapes
  # fixed -- the R equivalent of hz.death.AVC's output feeding bh.death.AVC.
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512, nu = 1.44, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, conserve = TRUE))

  bs <- hzr_bootstrap(base, n_boot = 10L, seed = 222L,
                       scope = list(early    = ~ age + mal,
                                    constant = ~ status + com_iv),
                       slentry = 0.3, slstay = 0.2,
                       control = list(n_starts = 1, conserve = TRUE))

  expect_s3_class(bs, "hzr_bootstrap")
  expect_identical(bs$mode, "select")
  expect_gt(bs$n_success, 0L)
  # Phase-qualified shape parameters are never dropped by stepwise
  # selection, so both phases must be represented in the summary.
  params <- bs$summary$parameter
  expect_true(any(grepl("^early\\.", params)))
  expect_true(any(grepl("^constant\\.", params)))
})
```

- [ ] **Step 2: Run it**

Run: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-sas-parity.R')"`
Expected: This new test passes. SAS-fixture-gated tests in the same file will `skip()` if you don't have the local SAS fixture directory configured — that's expected and unrelated to this change; confirm no new failures among the tests that do run.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-sas-parity.R
git commit -m "test: exercise hzr_bootstrap(scope=) against real AVC multiphase data

Upgrades the bs.death.AVC block from 'parses SAS reference frequencies
only' to 'also runs the equivalent R scope= path.' Still does not
assert R-vs-SAS frequency parity (documented, accepted gap)."
```

---

### Task 4: NEWS.md entry

**Files:**
- Modify: `NEWS.md:1`

**Interfaces:** None (documentation only).

- [ ] **Step 1: Add a "New features" section under the open dev-cycle header**

Current `NEWS.md` top:

```markdown
# TemporalHazard 1.1.0.9000 (development version)

# TemporalHazard 1.1.0
```

Replace with:

```markdown
# TemporalHazard 1.1.0.9000 (development version)

## New features

* `hzr_bootstrap()` gains a `scope` argument for embedded stepwise variable
  selection during each bootstrap replicate -- the R equivalent of SAS's
  `%HAZBOOT` procedure. Each replicate runs a fresh `hzr_stepwise()`
  selection (starting from a fixed-shape refit of the base model) instead
  of a plain refit, so `summary$pct` reports the variable's selection
  frequency across resamples and `summary$mean`/`sd`/`ci_*` describe the
  coefficient distribution conditional on selection. `scope = NULL`
  (the default) preserves the original fixed-formula bootstrap unchanged.

# TemporalHazard 1.1.0
```

This is a `.9000` dev-cycle addition, not a version bump — per project convention the minor/major digit is reserved for a deliberate, human-reviewed release.

- [ ] **Step 2: Confirm the DESCRIPTION/NEWS version-sync test still passes**

Run: `Rscript -e "devtools::test()"` (or, if there is a dedicated test file for this, e.g. `tests/testthat/test-version-sync.R`: `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-version-sync.R')"`)
Expected: Passes — `DESCRIPTION` (`Version: 1.1.0.9000`) is unchanged, and `NEWS.md`'s top `# TemporalHazard 1.1.0.9000 (development version)` header is unchanged, only content was added beneath it.

- [ ] **Step 3: Commit**

```bash
git add NEWS.md
git commit -m "docs: add NEWS.md entry for hzr_bootstrap(scope=)"
```

---

### Task 5 (separate project, not git-tracked): neo_therapy worked example

**Location:** `/Volumes/qhsstudies/thoracic/esophagus/malignant/neo_therapy/analyses/` (not a git repository — no commit step; the deliverable is the saved script file).

**Files:**
- Create: `/Volumes/qhsstudies/thoracic/esophagus/malignant/neo_therapy/analyses/bh.dead_r_bootstrap.qmd`

**Interfaces:**
- Consumes: `hzr_bootstrap(object, scope = ..., slentry = ..., slstay = ..., ...)` from Task 2, installed from the `feature/bootstrap-stepwise-selection` branch (or `dev`, once this plan's commits are merged).

This mirrors the `hz.dead.sas` → `bh.dead.sas` hand-off using the real `bhblt.sas7bdat` dataset, already read the same way in `bh.dead_r_test.qmd:305`. Per the earlier scoping decision: a **trimmed** ~15-20 variable subset (not the full ~90-variable `bh.dead.sas` candidate list) and a **small** `n_boot` (50-100, not SAS's 1000), so this runs in a couple of minutes as a first pass. Scaling to the full list/`n_boot` is a deliberate follow-up left for you to run once the mechanics are confirmed working.

- [ ] **Step 1: Install the branch locally**

Run: `Rscript -e "devtools::install_github('ehrlinger/temporal_hazard', ref = 'feature/bootstrap-stepwise-selection')"` (or, if working from a local clone: `Rscript -e "devtools::install('/Users/ehrlinj/Documents/GitHub/temporal_hazard')"` with that branch checked out).
Expected: `TemporalHazard` installs with `hzr_bootstrap()` accepting `scope=`. Confirm with `Rscript -e "library(TemporalHazard); 'scope' %in% names(formals(hzr_bootstrap))"` → `TRUE`.

- [ ] **Step 2: Create the worked-example script**

Create `/Volumes/qhsstudies/thoracic/esophagus/malignant/neo_therapy/analyses/bh.dead_r_bootstrap.qmd` with:

````markdown
---
title: "BH Death -- R Bootstrap Variable Screening (scope=)"
date: today
format: html
---

```{r}
#| label: setup
library(TemporalHazard)
library(haven)

caa <- read_sas("../datasets/bhblt.sas7bdat")
colnames(caa) <- tolower(colnames(caa))
```

## Step 1: fit the hazard shape (hz.dead.sas equivalent)

Shape values below match the fitted `hz.dead.sas` output already used in
`bh.dead_r_test.qmd` (`t_half = 1.666871, nu = -.7462, m = 0`), fixed here
via `fixed = "shapes"` so the bootstrap-screening step below only varies
covariates.

```{r}
#| label: hz-shape

base_fit <- hazard(
  survival::Surv(iv_dead, dead) ~ 1,
  data   = caa,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 1.666871, nu = -.7462, m = 0,
                          fixed = "shapes"),
    constant = hzr_phase("constant")
  ),
  fit = TRUE
)
summary(base_fit)
```

## Step 2: bootstrap-screen candidate covariates (bh.dead.sas equivalent)

Trimmed candidate list (representative subset of `bh.dead.sas`'s ~90-variable
list) and `n_boot = 50` for a first pass -- scale both up once this runs
cleanly. `slentry`/`slstay` match `bh.dead.sas`'s actual `SLE = 0.12`,
`SLS = 0.10` call (`analyses/bh.dead.sas:326`).

```{r}
#| label: bh-bootstrap-screen

candidate_scope <- list(
  early = ~ male + race_wh + hx_mi + dysphpri + symp_wl + symp_cp +
    barretts + preadj + ct1 + ct2 + cn0 + cn1 + hist_ad + grade +
    tum_low + resm1 + lvi + age + bmipr,
  constant = ~ male + race_wh + hx_mi + dysphpri + symp_wl + symp_cp +
    barretts + preadj + ct1 + ct2 + cn0 + cn1 + hist_ad + grade +
    tum_low + resm1 + lvi + age + bmipr
)

bs <- hzr_bootstrap(
  base_fit, n_boot = 50, seed = 123,
  scope    = candidate_scope,
  slentry  = 0.12, slstay = 0.10,
  control  = list(n_starts = 1)
)
print(bs)
```

## Step 3: screening output

```{r}
#| label: screening-table

bs$summary[order(-bs$summary$pct), ]
```
````

- [ ] **Step 2: Run it and inspect the output**

Run: `Rscript -e "quarto::quarto_render('/Volumes/qhsstudies/thoracic/esophagus/malignant/neo_therapy/analyses/bh.dead_r_bootstrap.qmd')"`
Expected: Renders without error; `bh-bootstrap-screen` chunk completes in roughly 1-3 minutes (50 replicates x trimmed scope); the `screening-table` chunk prints a `summary` data frame with `pct` (selection frequency) descending. Report the rendered `pct` values back — variables entering in a small fraction of replicates are the ones `bh.dead.sas` would flag for screening out.

No commit step: this directory is not a git repository (per `CLAUDE.md`'s "No PHI in code, notes, or repos" and the earlier scoping decision, this file and the real `bhblt.sas7bdat` data it reads stay local, never pushed to any git remote).

---

## Self-Review

**Spec coverage:** API signature (Task 2) ✓. Accumulation logic / per-replicate param names (Task 1 + Task 2) ✓. Summary/print unchanged mechanism + mode line (Task 2) ✓. Result object `mode`/`scope` fields (Task 2) ✓. Error handling — corrected pre-loop validation + documented warn-not-raise behavior (Task 2, tests 3 & 4) ✓. Out-of-scope items (no SAS path parity, no `scope=NULL` behavior change, no new S3 class) — respected throughout, not implemented ✓. Testing items 1-5 from the spec (regression, toy scope test, multiphase AVC exercise, structural bad-argument, typo'd-column) — Tasks 2 & 3 ✓. Worked example — Task 5 ✓.

**Placeholder scan:** No TBD/TODO; every step has complete, runnable code.

**Type consistency:** `.hzr_bootstrap_param_names(fit_obj)` signature and return type consistent between Task 1 (definition) and Task 2 (two call sites). `hzr_bootstrap()`'s new argument names (`scope`, `direction`, `criterion`, `slentry`, `slstay`, `max_steps`, `max_move`, `force_in`, `force_out`) match `hzr_stepwise()`'s existing parameter names exactly (`R/stepwise.R:151-163`) at every forwarding call site (pre-loop validation call and per-replicate call).
