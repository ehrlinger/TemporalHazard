# Limit restart (#683, #729) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** When a multiphase fit's shape parameter runs to a limit, `hazard()` restarts from a
fixed grid of starts and returns the best interior fit, or refuses with a classed error. It
never returns the drifted iterate.

**Architecture:** `hazard()` gains an outer layer for `fit = TRUE, dist = "multiphase"`.
1. The layer evaluates the caller's own call once, as an "inner" call, with its warnings
   captured.
2. It reads a trigger from the fit object's existing `$weak` and `$boundary`, plus one new
   boundary record.
3. On a trigger, it re-evaluates the same call once per grid start, and keeps the best fit
   that does not trigger, or refuses.

Every restart is a whole `hazard()` run, so CoE, vcov and all the post-fit records are
ordinary. A package-private flag stops recursion.

**Tech Stack:** R (>= 4.1), testthat 3e, roxygen2. No new dependencies.

**Spec:** `inst/dev/LIMIT-RESTART-DESIGN.md` (commit `c5bb5cec`). Read it first.

## Global Constraints

- Branch `design/limit-restart-683`. Never push to `main`. Never merge. The maintainer
  merges.
- **Version.** No version bump unless the maintainer asks. Never `.9000`. Never touch the
  MINOR or MAJOR digit.
- **NEWS.** A NEWS entry goes in `news/design-limit-restart-683.md` (bullets only, no
  heading). Do not edit `NEWS.md`. Spell-check the fragment with
  `spelling::spell_check_files()`.
- **Style.** `stats::` prefix on every stats generic. No `print()`, `library()` or
  `browser()` in `R/`. `cat()` only inside `print.*`.
- **No em-dashes** in anything written. A hook blocks them.
- **Definition of done**, in order: `Rscript -e 'devtools::document()'`, then
  `Rscript -e 'lintr::lint_package()'` with 0 lints, then
  `NOT_CRAN=true Rscript -e 'devtools::test()'` with 0 failures.
- **Once per PR:** commit, then `git archive HEAD | tar -x -C "$TMPDIR/tree"`, then
  `R CMD build`, then `R CMD check --as-cran` with the manual (never `--no-manual`).
- **Slow tests.** A new test file slower than about 1 s with `NOT_CRAN=false` gets a
  file-level `skip_on_cran()` as its first code line.
- **Assert that a computation varied or covered what it claims**, not that it returned.
  Prefer assertions that can fail.
- **Stamp counts.** Write "N at <sha>", `saveRDS()` the full suite result, and report
  per-file pass deltas, explaining every drop.
- **No PHI.** `valves` and `avc` are package datasets. No study-share paths. Run nothing
  over `/Volumes`.

---

## File structure

| File | Responsibility |
|---|---|
| `R/limit-restart.R` (new) | The private state flag, `.hzr_limit_trigger()`, `.hzr_limit_grid()`, `.hzr_limit_capture()`, `.hzr_limit_refusal()`, `.hzr_limit_fit()`, `.hzr_limit_line()` |
| `R/boundary-detect.R` | New `.hzr_phase_beyond_data_record()` beside `.hzr_phase_step_record()` |
| `R/hazard_api.R` | Call the new record in `.hzr_boundary_check_impl()`; a lead-in in `.hzr_boundary_message()`; the outer layer at the top of `hazard()`; `$fit$limit` default; `limit_restart` in `.hzr_control_names`; the print line; roxygen for `control$limit_restart` and `$fit$limit` |
| `R/likelihood-multiphase.R` | Strip `control$limit_restart` before `optim()` |
| `R/stepwise-step.R` | `.hzr_refit_failure_reason()` reads the refusal's short reason |
| `R/diagnostics.R` | `hzr_bootstrap()` failure tally reads the refusal's short reason |
| `tests/testthat/test-limit-trigger.R` (new) | Trigger and grid unit tests on hand-built objects (fast, CRAN-safe) |
| `tests/testthat/test-boundary-beyond-data.R` (new) | The `phase_beyond_data` record |
| `tests/testthat/test-limit-restart.R` (new) | End-to-end: #729 acceptance, refusal, varied grid, print (file-level `skip_on_cran()`) |
| `tests/testthat/test-limit-callers.R` (new) | Stepwise and bootstrap handling of a refusal (file-level `skip_on_cran()`) |
| `tests/testthat/test-boundary-phase-discontinuity.R`, `test-unbounded-phase.R`, `test-boundary-step-detail.R` | Migrated to the new behaviour (Task 8) |
| `tests/testthat/test-sas-parity.R` | `hm.death.patient` warm start becomes a blind start (Task 8) |
| `news/design-limit-restart-683.md` (new) | NEWS fragment |

---

### Task 1: Measure before building (no commit)

This task checks the spec's premises against live code before anything is built. Its
output is a short table for the maintainer. **Stop at the checkpoint.**

**Files:** none committed.
- Create: `$SCRATCH/limit-probe.R`, where `$SCRATCH` is the session scratchpad. Never write
  into the repo or `R_HOME`.

- [ ] **Step 1: Record the baseline suite at the branch base**

```bash
git rev-parse --short HEAD
NOT_CRAN=true Rscript -e 'r <- as.data.frame(devtools::test(reporter = "silent")); saveRDS(r, Sys.getenv("SCRATCH_RDS")); cat(sum(r$passed), "passed,", sum(r$failed), "failed,", sum(r$skipped), "skipped\n")'
```

Run it with `SCRATCH_RDS=$SCRATCH/suite-base.rds`.
- Expected: 0 failed.
- Write down "N passed at <sha>". This is the baseline for Task 9's per-file delta.

- [ ] **Step 2: Write the probe script**

```r
suppressMessages(devtools::load_all(quiet = TRUE, helpers = FALSE))
stopifnot(packageVersion("TemporalHazard") >= "1.3.0")
capture <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  list(value = v, warnings = w)
}
report <- function(label, r) {
  v <- r$value
  if (inherits(v, "error")) { cat(label, ": ERROR", conditionMessage(v), "\n"); return(invisible()) }
  b <- v$fit$boundary
  cat(sprintf("%s: LL %.2f converged %s weak %s boundary %s\n", label, v$fit$objective,
    isTRUE(v$fit$converged),
    if (is.list(v$fit$weak)) paste(v$fit$weak$params, collapse = ",") else format(v$fit$weak),
    if (is.list(b)) paste(vapply(b, `[[`, "", "mechanism"), collapse = ",") else format(b)))
  print(round(v$fit$theta[grepl("log_t_half|\\.nu$|log_mu", names(v$fit$theta))], 3))
}
data(valves, package = "TemporalHazard"); d <- valves
d$nyha[is.na(d$nyha)] <- mean(d$nyha, na.rm = TRUE); d$mv_nyha <- d$mitral * d$nyha
fit729 <- function(n_starts, t_half = 10, nu = 0.2, theta = NULL) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = t_half, nu = nu, m = 1, fixed = "m",
                    formula = ~ age_cop + nyha + mv_nyha + double_ + ao_pinc),
                  constant = hzr_phase("constant", formula = ~ age_cop + black + i_path)),
    fit = TRUE, theta = theta,
    control = list(n_starts = n_starts, maxit = 2000, conserve = TRUE))
}
report("729 n1", capture(fit729(1L)))
report("729 n5", capture(fit729(5L)))
# The grid the spec proposes, run by hand: event-time quartiles x nu in {0.5, 1, 2}.
q <- unname(stats::quantile(d$int_dead[d$dead == 1], c(.25, .5, .75), na.rm = TRUE))
for (th in q) for (nu in c(0.5, 1, 2)) {
  r <- capture(fit729(1L, t_half = th, nu = nu))
  cat(sprintf("grid t_half %.3f nu %.1f -> ", th, nu)); report("", r)
}
```

Then append the #448 fixture (copy `avc_fit_448()` verbatim from
`tests/testthat/test-boundary-phase-discontinuity.R:12-26`) and the #444 fixture (copy the
smallest fitting helper from `tests/testthat/test-unbounded-phase.R`), and `report()` each.

- [ ] **Step 3: Run it**

```bash
Rscript "$SCRATCH/limit-probe.R" 2>&1 | tee "$SCRATCH/limit-probe.log"
```

Expected:
- #729 n1 at LL about -2036.66, n5 at about -1817.64.
- Most grid rows at about -1792.07.
- Record each fixture's `weak` params and `boundary` mechanisms exactly as printed.

- [ ] **Step 4: Checkpoint. Report the table to the maintainer and stop.**

The table has one row per fixture: LL, weak params, boundary mechanisms, and whether the
spec's trigger would fire. It also gives the grid's 9 LLs.

**Stop and ask before Task 2 if any of these hold:**
- (a) #729 n1's weak params name no shape. The trigger would then miss it.
- (b) Fewer than 1 of the 9 grid starts reaches -1792.07 to 0.01.
- (c) The #448 or #444 fixture would trigger. Their fits would then restart or refuse
  instead of returning with a warning, which is a user-visible change to fits SAS
  completes. The maintainer decides whether that is wanted before their tests are
  migrated in Task 8.

---

### Task 2: The `phase_beyond_data` boundary record

**Files:**
- Modify: `R/boundary-detect.R` (append after `.hzr_phase_step_record()`)
- Modify: `R/hazard_api.R:95-104` (inside the loop of `.hzr_boundary_check_impl()`) and
  `R/hazard_api.R:189-192` (the `lead` vector of `.hzr_boundary_message()`)
- Test: `tests/testthat/test-boundary-beyond-data.R`

**Interfaces:**
- Produces: `.hzr_phase_beyond_data_record(name, type, theta, t_obs)`. It returns `NULL`, or
  `list(mechanism = "phase_beyond_data", phase, parameter = "t_half", t_half, t_max, G_max, detail)`.

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-boundary-beyond-data.R
# #683: a G1 phase whose rise begins after the last observed time. Its t_half
# is not identified; the phase is a power law in t over the data (#729 drifted
# to log t_half 275.8 with converged = TRUE).

theta_for <- function(log_t_half, nu = 1, m = 1) {
  c(early.log_mu = 0, early.log_t_half = log_t_half, early.nu = nu, early.m = m)
}

test_that("a cdf phase whose rise lies past the data is recorded", {
  rec <- .hzr_phase_beyond_data_record("early", "cdf", theta_for(275.8),
                                       t_obs = c(0.1, 5, 17))
  expect_identical(rec$mechanism, "phase_beyond_data")
  expect_identical(rec$phase, "early")
  expect_identical(rec$parameter, "t_half")
  expect_lt(rec$G_max, .hzr_step_resolution)
  expect_match(rec$detail, "after the data ends")
})

test_that("a t_half past the data whose rise is visible is NOT recorded", {
  # t_half 2x beyond the last time: G(17) is far from 0, the shape is seen.
  expect_null(.hzr_phase_beyond_data_record("early", "cdf", theta_for(log(34)),
                                            t_obs = c(0.1, 5, 17)))
})

test_that("an interior t_half is not recorded, and nor is a constant phase", {
  expect_null(.hzr_phase_beyond_data_record("early", "cdf", theta_for(log(3)),
                                            t_obs = c(0.1, 5, 17)))
  expect_null(.hzr_phase_beyond_data_record(
    "c", "constant", c(c.log_mu = 0), t_obs = c(0.1, 5, 17)))
})

test_that("a t_half that overflows exp() is still recorded", {
  rec <- .hzr_phase_beyond_data_record("early", "cdf", theta_for(800),
                                       t_obs = c(0.1, 5, 17))
  expect_identical(rec$mechanism, "phase_beyond_data")
})

test_that("the boundary check carries the record into $boundary", {
  out <- .hzr_boundary_check_impl(
    theta = theta_for(275.8),
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1)),
    time = c(0.1, 5, 17), fitted = TRUE)
  mech <- vapply(out$boundary, `[[`, "", "mechanism")
  expect_true("phase_beyond_data" %in% mech)
})
```

- [ ] **Step 2: Run them and verify they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-boundary-beyond-data.R")'`

Expected: FAIL, could not find function `.hzr_phase_beyond_data_record`.

- [ ] **Step 3: Implement the record**

Append to `R/boundary-detect.R`:

```r
#' Build a `$boundary` record for a G1 phase whose rise lies past the data
#'
#' The mirror of #444's `unbounded_phase`: there `t_half` sits below the
#' first observed time; here it sits beyond the last, so far that `G` is still
#' numerically 0 there. Over the observed range the phase is then a power law
#' in t and `t_half` is not identified: `log_mu` and `log_t_half` trade along
#' a ridge to infinity (#683; #729 reported converged = TRUE at log t_half
#' 275.8). Like #448 this asks the data, not a threshold on the parameter:
#' `.hzr_step_resolution` is the "indistinguishable from 0" scale.
#'
#' @param name,type The phase's name and type.
#' @param theta The fitted, named full theta.
#' @param t_obs Every positive, finite observed time the likelihood read.
#' @return A `$boundary` record, or `NULL`.
#' @keywords internal
#' @noRd
.hzr_phase_beyond_data_record <- function(name, type, theta, t_obs) {
  if (!type %in% c("cdf", "hazard") || !length(t_obs)) return(NULL)
  keys <- paste0(name, c(".log_t_half", ".nu", ".m"))
  if (!all(keys %in% names(theta))) return(NULL)
  log_t_half <- unname(theta[[keys[[1L]]]])
  t_max <- max(t_obs)
  if (!is.finite(log_t_half) || log_t_half <= log(t_max)) return(NULL)
  t_half <- exp(log_t_half)
  # A t_half that overflows, or a decomposition that fails, leaves G
  # unknown. That is still a phase past the data, so it is recorded with
  # G "not computable" rather than dropped: dropping it would let the
  # extreme case through and catch only the mild one.
  g_max <- if (is.finite(t_half)) {
    tryCatch(hzr_decompos(t_max, t_half = t_half,
                          nu = unname(theta[[keys[[2L]]]]),
                          m = unname(theta[[keys[[3L]]]]))$G,
             error = function(e) NA_real_)
  } else {
    NA_real_
  }
  if (is.finite(g_max) && g_max >= .hzr_step_resolution) return(NULL)
  list(
    mechanism = "phase_beyond_data", phase = name, parameter = "t_half",
    t_half = t_half, t_max = t_max, G_max = g_max,
    detail = paste0(
      "phase '", name, "' has log t_half ", format(log_t_half, digits = 4),
      " against a last observed time of ", format(t_max, digits = 4),
      ", and G there is ",
      if (is.finite(g_max)) format(g_max, digits = 3) else "not computable",
      ": its rise begins after the data ends, so t_half is not identified ",
      "and the phase is a power law in t over the observed range.")
  )
}
```

- [ ] **Step 4: Call it from the boundary check**

In `R/hazard_api.R`, inside `.hzr_boundary_check_impl()`, directly after the
`if (!is.null(step_rec)) found[[length(found) + 1L]] <- step_rec` line (currently 104) and
before `if (!.hzr_phase_type_unbounded(...)) next`, insert:

```r
    # #683: the mirror of the unbounded record below, a rise past the data.
    beyond_rec <- .hzr_phase_beyond_data_record(nm, phases[[k]]$type, theta,
                                                t_obs = t_ok)
    if (!is.null(beyond_rec)) found[[length(found) + 1L]] <- beyond_rec
```

In `.hzr_boundary_message()`, add to the `lead` vector:

```r
            phase_beyond_data = "phase rises after the data ends: ",
```

- [ ] **Step 5: Run the tests and verify they pass**

Run: the command from Step 2. Expected: all PASS.

- [ ] **Step 6: Run the boundary neighbours for regressions**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in c("test-unbounded-phase.R","test-boundary-phase-discontinuity.R","test-boundary-step-detail.R","test-degraded-record.R")) testthat::test_file(file.path("tests/testthat", f))'`

Expected: 0 failures. A new `phase_beyond_data` record on an existing fixture is a finding:
list it, and do not loosen any assertion.

- [ ] **Step 7: Commit**

```bash
git add R/boundary-detect.R R/hazard_api.R tests/testthat/test-boundary-beyond-data.R
git commit -m "Boundary record for a phase whose rise lies past the data (#683)"
```

---

### Task 3: The trigger

**Files:**
- Create: `R/limit-restart.R`
- Test: `tests/testthat/test-limit-trigger.R`

**Interfaces:**
- Consumes: `$fit$boundary` records, now including `phase_beyond_data` from Task 2;
  `$fit$weak` (`NULL`, `NA`, or a list with `$params`); `$fit$fixed_mask` (logical aligned
  with `$fit$theta`, or `NULL`); `$spec$phases`.
- Produces: `.hzr_limit_trigger(object)`. It returns `NULL`, or a list of
  `list(phase, parameter, limit, source)`, where:
  - `parameter` is the user-facing name: `"t_half"`, `"nu"`, `"m"`, `"tau"`, `"gamma"`,
    `"alpha"` or `"eta"`;
  - `limit` is a string;
  - `source` is a mechanism name or `"weak"`.
- Also produces `.hzr_limit_shape_names(phases)`: a named character vector mapping each
  internal shape name (`"early.log_t_half"`) to its phase.

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-limit-trigger.R
# #683: which fits the restart path acts on. Hand-built objects: fast.

phases2 <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1, fixed = "m"),
                constant = hzr_phase("constant"))
mk <- function(boundary = NULL, weak = NULL, fixed = NULL) {
  th <- c(early.log_mu = 0, early.log_t_half = 0, early.nu = 1, early.m = 1,
          constant.log_mu = 0, early.age = 0.1)
  list(spec = list(dist = "multiphase", phases = phases2),
       fit = list(theta = th, boundary = boundary, weak = weak,
                  fixed_mask = if (is.null(fixed)) NULL else names(th) %in% fixed))
}

test_that("a clean fit does not trigger", {
  expect_null(.hzr_limit_trigger(mk()))
  expect_null(.hzr_limit_trigger(mk(boundary = NA, weak = NA)))
})

test_that("each boundary mechanism maps to its parameter and limit", {
  cases <- list(
    phase_beyond_data   = c("t_half", "Inf (the rise lies past the data)"),
    unbounded_phase     = c("t_half", "0 (below the first observed time)"),
    phase_discontinuity = c("nu", "0 (a step at t_half)"))
  for (mech in names(cases)) {
    tr <- .hzr_limit_trigger(mk(boundary = list(list(mechanism = mech, phase = "early"))))
    expect_length(tr, 1L)
    expect_identical(tr[[1]]$phase, "early")
    expect_identical(tr[[1]]$parameter, cases[[mech]][[1]])
    expect_identical(tr[[1]]$limit, cases[[mech]][[2]])
    expect_identical(tr[[1]]$source, mech)
  }
})

test_that("records that are not limits do not trigger", {
  for (mech in c("g3_alpha_one", "g3_fixge2_alpha_start", "g3_corner_supremum",
                 "coe_no_events_left", "coe_phase_vanished")) {
    expect_null(.hzr_limit_trigger(mk(boundary = list(list(mechanism = mech,
                                                           phase = "early")))))
  }
})

test_that("a weak direction triggers only when it names a shape", {
  tr <- .hzr_limit_trigger(mk(weak = list(params = c("early.log_mu", "early.nu"))))
  expect_identical(vapply(tr, `[[`, "", "parameter"), "nu")
  expect_identical(tr[[1]]$source, "weak")
  # A covariate ridge is collinearity, not a limit.
  expect_null(.hzr_limit_trigger(mk(weak = list(params = c("early.age", "constant.log_mu")))))
})

test_that("a covariate NAMED like a shape is not read as one", {
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1))
  obj <- list(spec = list(dist = "multiphase", phases = ph),
              fit = list(theta = c(early.log_mu = 0, early.log_t_half = 0,
                                   early.nu = 1, early.m = 1, early.nu_1 = 0),
                         weak = list(params = "early.nu_1"), boundary = NULL))
  expect_null(.hzr_limit_trigger(obj))
})

test_that("a FIXED parameter never triggers", {
  expect_null(.hzr_limit_trigger(mk(weak = list(params = "early.m"), fixed = "early.m")))
  expect_null(.hzr_limit_trigger(mk(
    boundary = list(list(mechanism = "phase_discontinuity", phase = "early")),
    fixed = "early.nu")))
})

test_that("a weak reading in a VANISHED phase does not trigger", {
  van <- list(list(mechanism = "coe_phase_vanished", phase = "early"))
  expect_null(.hzr_limit_trigger(mk(boundary = van,
                                    weak = list(params = c("early.log_t_half", "early.nu")))))
  # but a limit record in that phase still does
  tr <- .hzr_limit_trigger(mk(boundary = c(van, list(list(mechanism = "unbounded_phase",
                                                            phase = "early")))))
  expect_identical(tr[[1]]$source, "unbounded_phase")
})

test_that("one parameter reached two ways is reported once", {
  tr <- .hzr_limit_trigger(mk(
    boundary = list(list(mechanism = "phase_discontinuity", phase = "early")),
    weak = list(params = c("early.nu", "early.log_mu"))))
  expect_length(tr, 1L)
})
```

- [ ] **Step 2: Run and verify they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-limit-trigger.R")'`

Expected: FAIL, could not find function `.hzr_limit_trigger`.

- [ ] **Step 3: Implement**

Create `R/limit-restart.R`:

```r
# limit-restart.R -- a phase parameter running to its limit (#683, #729)
#
# Design: inst/dev/LIMIT-RESTART-DESIGN.md. A multiphase fit can walk a
# shape parameter toward a limit of the BNT parameterisation and return the
# last iterate with converged = TRUE. hazard() now restarts from a fixed grid
# and returns the best interior fit, or refuses. It never switches the model
# and never returns the drifted iterate. SAS/C has no such check: it chooses
# a limiting form only from values the analyst fixed (setg1.c:620-636).

# Boundary mechanisms that ARE a shape at a limit, with the parameter and the
# limit each one names. Anything else on $boundary (setup holds, the g3
# corner, CoE) is handled where it is raised and does not trigger.
.hzr_limit_mechanisms <- list(
  phase_beyond_data   = c("t_half", "Inf (the rise lies past the data)"),
  unbounded_phase     = c("t_half", "0 (below the first observed time)"),
  phase_discontinuity = c("nu", "0 (a step at t_half)")
)

# Internal shape name -> the hzr_phase() argument a user would fix.
.hzr_limit_user_name <- c(log_t_half = "t_half", nu = "nu", m = "m",
                          log_tau = "tau", gamma = "gamma", alpha = "alpha",
                          eta = "eta")

#' Every shape parameter's full name, mapped to its phase
#'
#' Built from the phases' own naming, never by suffix: a covariate column
#' called `nu` is `early.nu` too, and would be misread (classify by the
#' fit's structure, not by a string's shape).
#' @param phases Named list of `hzr_phase` objects.
#' @return Named character vector: names are `<phase>.<shape>`, values the
#'   phase name.
#' @keywords internal
#' @noRd
.hzr_limit_shape_names <- function(phases) {
  out <- character(0)
  for (nm in names(phases)) {
    shapes <- .hzr_phase_theta_names(phases[[nm]], nm)[-1L]  # drop log_mu
    out <- c(out, stats::setNames(rep(nm, length(shapes)), shapes))
  }
  out
}

#' Does a fitted multiphase object show a shape parameter at a limit?
#'
#' Reads only what the object carries: `$fit$boundary` records whose
#' mechanism is a limit, and a `$fit$weak` direction that names a shape.
#' A parameter the user fixed never triggers: holding it there was their
#' choice, as FIXNU is in SAS.
#'
#' @param object A `hazard` object, or a list with `$spec$phases` and `$fit`.
#' @return `NULL`, or a list of `list(phase, parameter, limit, source)`.
#' @keywords internal
#' @noRd
.hzr_limit_trigger <- function(object) {
  phases <- object$spec$phases
  fs <- object$fit
  if (!length(phases)) return(NULL)
  shape_phase <- .hzr_limit_shape_names(phases)
  theta_names <- names(fs$theta)
  fixed <- if (length(fs$fixed_mask) == length(theta_names)) {
    theta_names[as.logical(fs$fixed_mask)]
  } else {
    character(0)
  }
  internal_of <- function(phase, user) {
    key <- names(.hzr_limit_user_name)[.hzr_limit_user_name == user]
    paste0(phase, ".", key)
  }
  hits <- list()
  add <- function(phase, parameter, limit, source) {
    if (internal_of(phase, parameter) %in% fixed) return(invisible())
    id <- paste0(phase, ".", parameter)
    if (id %in% names(hits)) return(invisible())
    hits[[id]] <<- list(phase = phase, parameter = parameter, limit = limit,
                        source = source)
  }
  if (is.list(fs$boundary)) {
    for (r in fs$boundary) {
      m <- .hzr_limit_mechanisms[[r$mechanism %||% ""]]
      if (!is.null(m)) add(r$phase, m[[1L]], m[[2L]], r$mechanism)
    }
  }
  # A phase recorded as vanished has unidentified shapes because its mu went
  # to 0: a different limit, already warned (#565), and not a shape at a
  # limit. Its weak reading does not trigger (amended after Task 1: #444's
  # t_half = 3 fit vanished from every grid start).
  vanished <- if (is.list(fs$boundary)) {
    unlist(lapply(fs$boundary, function(r) {
      if (identical(r$mechanism, "coe_phase_vanished")) r$phase
    }))
  }
  if (is.list(fs$weak)) {
    for (p in intersect(fs$weak$params, names(shape_phase))) {
      if (shape_phase[[p]] %in% vanished) next
      short <- substring(p, nchar(shape_phase[[p]]) + 2L)
      add(shape_phase[[p]], .hzr_limit_user_name[[short]],
          "a flat direction (not identified)", "weak")
    }
  }
  if (!length(hits)) NULL else unname(hits)
}
```

- [ ] **Step 4: Run and verify they pass**

Run: the command from Step 2. Expected: all PASS.

- [ ] **Step 5: Planted faults**

Do each of the following in turn, re-run, confirm the named test goes red, then revert:
- delete the `phase_discontinuity` entry from `.hzr_limit_mechanisms`;
- replace `intersect(fs$weak$params, names(shape_phase))` with `fs$weak$params`;
- delete the `if (internal_of(...) %in% fixed) return(invisible())` line.

Record in the PR body which test killed each mutant.

- [ ] **Step 6: Commit**

```bash
git add R/limit-restart.R tests/testthat/test-limit-trigger.R
git commit -m "Limit trigger read from \$weak and \$boundary (#683)"
```

---

### Task 4: The restart grid

**Files:**
- Modify: `R/limit-restart.R`
- Test: `tests/testthat/test-limit-trigger.R` (append)

**Interfaces:**
- Consumes: the trigger from Task 3; `.hzr_limit_shape_names()`.
- Produces: `.hzr_limit_grid(start, trigger, phases, event_times, fixed)`. It returns a list
  of named numeric theta vectors, each the same length and names as `start`, with no two
  identical.
  - `start`: a named full theta, the user's start.
  - `event_times`: the times of rows with `status == 1`.
  - `fixed`: names not to move.

- [ ] **Step 1: Write the failing tests (append)**

```r
test_that("the grid is 3 t_half x 3 nu, covariates zeroed, nothing else moved", {
  start <- c(early.log_mu = -1, early.log_t_half = log(10), early.nu = 0.2,
             early.m = 1, constant.log_mu = -3, early.age = 0.4, constant.age = 0.2)
  tr <- list(list(phase = "early", parameter = "nu", limit = "x", source = "weak"))
  et <- c(1, 2, 3, 4, 5, 6, 7, 8)
  g <- .hzr_limit_grid(start, tr, phases2, event_times = et, fixed = "early.m")
  expect_length(g, 9L)
  th <- vapply(g, `[[`, 0, "early.log_t_half")
  nu <- vapply(g, `[[`, 0, "early.nu")
  expect_setequal(unique(round(exp(th), 6)),
                  round(unname(stats::quantile(et, c(.25, .5, .75))), 6))
  expect_setequal(unique(nu), c(0.5, 1, 2))
  for (s in g) {
    expect_identical(names(s), names(start))
    expect_equal(unname(s[c("early.age", "constant.age")]), c(0, 0))
    expect_equal(unname(s[c("early.log_mu", "early.m", "constant.log_mu")]),
                 unname(start[c("early.log_mu", "early.m", "constant.log_mu")]))
  }
  expect_false(anyDuplicated(lapply(g, unname)) > 0L)
})

test_that("a fixed t_half collapses the grid to nu alone", {
  start <- c(early.log_mu = -1, early.log_t_half = 0, early.nu = 1, early.m = 1,
             constant.log_mu = -3)
  tr <- list(list(phase = "early", parameter = "nu", limit = "x", source = "weak"))
  g <- .hzr_limit_grid(start, tr, phases2, event_times = 1:8,
                       fixed = c("early.m", "early.log_t_half"))
  expect_length(g, 3L)
  expect_true(all(vapply(g, `[[`, 0, "early.log_t_half") == 0))
})

test_that("no events gives no grid", {
  start <- c(early.log_mu = -1, early.log_t_half = 0, early.nu = 1, early.m = 1,
             constant.log_mu = -3)
  tr <- list(list(phase = "early", parameter = "t_half", limit = "x", source = "weak"))
  expect_length(.hzr_limit_grid(start, tr, phases2, numeric(0), "early.m"), 0L)
})

test_that("a g3 trigger varies tau at the quartiles and nothing else", {
  ph <- list(late = hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1))
  start <- c(late.log_mu = 0, late.log_tau = 0, late.gamma = 2, late.alpha = 1,
             late.eta = 3)
  tr <- list(list(phase = "late", parameter = "tau", limit = "x", source = "weak"))
  g <- .hzr_limit_grid(start, tr, ph, event_times = 1:8, fixed = character(0))
  expect_length(g, 3L)
  expect_true(all(vapply(g, `[[`, 0, "late.gamma") == 2))
})
```

- [ ] **Step 2: Run and verify they fail**

Expected: FAIL, could not find function `.hzr_limit_grid`.

- [ ] **Step 3: Implement (append to `R/limit-restart.R`)**

```r
# The grid's nu values. Spread across the sigmoid's usable range on the
# natural scale, away from both limits (0 and large).
.hzr_limit_grid_nu <- c(0.5, 1, 2)

#' Restart starts for a triggered fit
#'
#' Covariate coefficients are zeroed, as SAS's RETRY_T does
#' (hazrd3.c:343-361). For each triggered G1 phase, t_half is set to the
#' 25/50/75% quantiles of the event times and crossed with nu in
#' `.hzr_limit_grid_nu`; for a triggered g3 phase, tau takes the same
#' quantiles. Every other entry keeps `start`, and a fixed entry never moves.
#' Deterministic: no RNG, so a restarted fit is reproducible from its call.
#'
#' @param start Named full theta: the user's start.
#' @param trigger Output of `.hzr_limit_trigger()`.
#' @param phases The fit's phases.
#' @param event_times Times of the rows with an event.
#' @param fixed Names in `start` that must not move.
#' @return A list of named theta vectors, without duplicates.
#' @keywords internal
#' @noRd
.hzr_limit_grid <- function(start, trigger, phases, event_times, fixed) {
  et <- event_times[is.finite(event_times) & event_times > 0]
  if (!length(et)) return(list())
  q <- unique(unname(stats::quantile(et, c(0.25, 0.5, 0.75))))
  shape_phase <- .hzr_limit_shape_names(phases)
  base <- start
  is_cov <- !names(base) %in% c(names(shape_phase), paste0(names(phases), ".log_mu"))
  base[is_cov & !names(base) %in% fixed] <- 0
  set_free <- function(s, key, value) {
    if (key %in% names(s) && !key %in% fixed) s[[key]] <- value
    s
  }
  out <- list()
  for (ph in unique(vapply(trigger, `[[`, "", "phase"))) {
    if (identical(phases[[ph]]$type, "g3")) {
      for (tq in q) out[[length(out) + 1L]] <- set_free(base, paste0(ph, ".log_tau"), log(tq))
    } else {
      for (tq in q) for (nu in .hzr_limit_grid_nu) {
        s <- set_free(base, paste0(ph, ".log_t_half"), log(tq))
        out[[length(out) + 1L]] <- set_free(s, paste0(ph, ".nu"), nu)
      }
    }
  }
  out[!duplicated(lapply(out, unname))]
}
```

- [ ] **Step 4: Run and verify they pass**

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add R/limit-restart.R tests/testthat/test-limit-trigger.R
git commit -m "Deterministic restart grid for a triggered phase (#683)"
```

---

### Task 5: The outer layer in `hazard()`: restart, refuse, `$fit$limit`, `control$limit_restart`

**Files:**
- Modify: `R/limit-restart.R`
- Modify: `R/hazard_api.R`:
  - the first statement of `hazard()` (line 927);
  - `.hzr_control_names` (line 3697);
  - the `$fit$limit` default, set just before `class(obj) <- "hazard"` (around line 2074);
  - the roxygen `@param control` and `@return`.
- Modify: `R/likelihood-multiphase.R:2725` (strip the key)
- Test: `tests/testthat/test-limit-restart.R`

**Interfaces:**
- Consumes: Tasks 2-4.
- Produces:
  - `control$limit_restart`: `TRUE` (default) or `FALSE`. Anything else is an error.
  - `fit$fit$limit`: `NA` when not examined, `NULL` when examined and clean, or
    `list(trigger, original = list(objective, theta, weak, boundary), grid = data.frame(start, objective, triggered, error), winner)`
    when a restart happened.
  - The condition class `hzr_limit_refusal` (also `error`), with fields `phase`,
    `parameter`, `limit`, `best_objective` and `trigger`.
  - The warning class `hzr_limit_restarted`.

- [ ] **Step 1: Write the failing end-to-end tests**

```r
# tests/testthat/test-limit-restart.R
skip_on_cran()
# #683/#729: a shape parameter running to its limit is restarted, or refused;
# the drifted iterate is never returned. hm.death.patient blind start.

valves_729 <- function() {
  data(valves, package = "TemporalHazard", envir = environment())
  d <- valves
  d$nyha[is.na(d$nyha)] <- mean(d$nyha, na.rm = TRUE)
  d$mv_nyha <- d$mitral * d$nyha
  d
}
fit_729 <- function(n_starts = 1L, limit_restart = TRUE, d = valves_729()) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 10, nu = 0.2, m = 1, fixed = "m",
                    formula = ~ age_cop + nyha + mv_nyha + double_ + ao_pinc),
                  constant = hzr_phase("constant", formula = ~ age_cop + black + i_path)),
    fit = TRUE,
    control = list(n_starts = n_starts, maxit = 2000, conserve = TRUE,
                   limit_restart = limit_restart))
}
sas_ll_729 <- -1792.07   # hm.death.patient.lst, asserted in test-sas-parity.R

for (ns in c(1L, 5L)) {
  test_that(paste0("#729 blind start, n_starts = ", ns, ", reaches SAS's optimum"), {
    expect_warning(fit <- fit_729(ns), class = "hzr_limit_restarted")
    expect_equal(fit$fit$objective / sas_ll_729, 1, tolerance = 1e-5)
    lim <- fit$fit$limit
    expect_type(lim, "list")
    # The drift it recovered from is on the record, and is the bad one.
    expect_lt(lim$original$objective, sas_ll_729 - 10)
    expect_true(any(vapply(lim$trigger, `[[`, "", "phase") == "early"))
    # The grid ran, varied, and the winner is a grid row.
    expect_gte(nrow(lim$grid), 2L)
    expect_gt(length(unique(round(lim$grid$objective[is.finite(lim$grid$objective)], 4))), 1L)
    expect_false(lim$grid$triggered[[lim$winner]])
    expect_equal(lim$grid$objective[[lim$winner]], fit$fit$objective)
    # A rescued fit is an ordinary fit, and is not itself at a limit.
    expect_null(.hzr_limit_trigger(fit))
    expect_true(isTRUE(fit$fit$converged))
    expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
    # Its call is the user's, so hzr_bootstrap() refits from the user's start.
    expect_null(fit$call$theta)
  })
}

test_that("limit_restart = FALSE refuses at once, naming the parameter", {
  err <- tryCatch(suppressWarnings(fit_729(1L, limit_restart = FALSE)),
                  hzr_limit_refusal = function(e) e)
  expect_s3_class(err, "hzr_limit_refusal")
  expect_identical(err$phase, "early")
  expect_true(err$parameter %in% c("nu", "t_half"))
  expect_match(conditionMessage(err), "hzr_phase\\(")
  expect_lt(err$best_objective, sas_ll_729 - 10)
})

test_that("#448: when the best objective is at the step, the fit refuses", {
  # Task 1 grid: best -21.66 at the step, best clean start -208.31. Returning
  # the clean start would be 130 units worse than today's -78.15.
  data(avc, package = "TemporalHazard", envir = environment())
  tt <- avc$int_dead; st <- avc$dead
  k <- is.finite(tt) & is.finite(st) & tt > 0
  err <- tryCatch(suppressWarnings(hazard(
    time = tt[k], status = st[k], dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
                  late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                                   constraint = "alpha_gamma_eta")),
    fit = TRUE, control = list(n_starts = 1L))),
    hzr_limit_refusal = function(e) e)
  expect_s3_class(err, "hzr_limit_refusal")
  expect_identical(err$parameter, "nu")
  # The best objective it reports is the step's, above today's returned fit.
  expect_gt(err$best_objective, -78.15)
})

test_that("a clean fit is unchanged and records $limit NULL", {
  data(avc, package = "TemporalHazard", envir = environment())
  k <- is.finite(avc$int_dead) & avc$int_dead > 0
  fit <- hazard(time = avc$int_dead[k], status = avc$dead[k], dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1L))
  expect_true("limit" %in% names(fit$fit))
  expect_null(fit$fit$limit)
})

test_that("an unfitted object records $limit NA", {
  fit <- hazard(time = c(1, 2, 3), status = c(1, 0, 1), dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = FALSE)
  expect_true(is.na(fit$fit$limit))
})

test_that("limit_restart must be TRUE or FALSE", {
  expect_error(fit_729(1L, limit_restart = NA), "limit_restart")
  expect_error(fit_729(1L, limit_restart = "yes"), "limit_restart")
})
```

**On the clean-fit fixture:** if Task 1 showed that this `avc` early+constant fit triggers,
replace it with the `hm.death.patient` fit warm-started at the SAS estimates. Its
`theta_sas` construction is in `test-sas-parity.R:205-300`, and a fit started there must
read `$limit` `NULL`.

- [ ] **Step 2: Run and verify they fail**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-limit-restart.R")'`

Expected: FAIL. At this point `limit_restart` triggers the "not an element any fit reads"
warning, and there is no restart.

- [ ] **Step 3: Accept and strip the control key**

In `R/hazard_api.R`, set the multiphase entry of `.hzr_control_names` to:

```r
  multiphase = c("n_starts", "conserve", "phase_share_tol", "start_seed",
                 "limit_restart")
```

In `R/likelihood-multiphase.R`, directly after `control$start_seed <- NULL` (line 2725),
add:

```r
  control$limit_restart <- NULL  # read by hazard()'s limit layer (#683)
```

- [ ] **Step 4: Implement the layer (append to `R/limit-restart.R`)**

```r
# TRUE while hazard() is evaluating a call on the limit layer's behalf, so
# the inner call fits and returns without checking again. Package-private;
# not a control element, so no user setting can turn the check off.
.hzr_limit_state <- new.env(parent = emptyenv())
.hzr_limit_state$inner <- FALSE

#' Evaluate a hazard() call as the limit layer's inner call
#'
#' Its warnings are held back (a discarded fit's warnings must not reach the
#' user as if they described the result) and its error is returned, not
#' raised.
#' @param cl A `hazard()` call. @param env Where to evaluate it.
#' @return `list(value, warnings)`; `value` is a fit or an error condition.
#' @keywords internal
#' @noRd
.hzr_limit_capture <- function(cl, env) {
  old <- .hzr_limit_state$inner
  .hzr_limit_state$inner <- TRUE
  on.exit(.hzr_limit_state$inner <- old, add = TRUE)
  caught <- list()
  value <- tryCatch(
    withCallingHandlers(eval(cl, env), warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }),
    error = function(e) e)
  list(value = value, warnings = caught)
}

.hzr_limit_describe <- function(trigger) {
  paste(vapply(trigger, function(r) {
    paste0(r$phase, ".", r$parameter, " ran to ", r$limit)
  }, character(1)), collapse = "; ")
}

#' The refusal condition
#' @keywords internal
#' @noRd
.hzr_limit_refusal <- function(trigger, best_objective, n_tried,
                               disabled = FALSE) {
  first <- trigger[[1L]]
  msg <- paste0(
    "hazard(): ", .hzr_limit_describe(trigger), ". ",
    if (disabled) {
      "No restart was tried (control$limit_restart = FALSE)"
    } else if (n_tried > 0L) {
      paste0("None of ", n_tried, " restarts reached a fit without a ",
             "parameter at a limit")
    } else {
      "No restart start could be built (no event times, or every shape fixed)"
    },
    "; the best objective reached was ", format(best_objective, digits = 8),
    ". The drifted estimate is not returned. Fix the parameter at a value ",
    "you choose, as SAS's FIX", toupper(first$parameter), " does: ",
    "hzr_phase(..., ", first$parameter, " = <value>, fixed = \"",
    first$parameter, "\"), or supply other starting values (`theta`).")
  structure(
    class = c("hzr_limit_refusal", "error", "condition"),
    list(message = msg, call = NULL, phase = first$phase,
         parameter = first$parameter, limit = first$limit,
         best_objective = best_objective, trigger = trigger))
}

.hzr_limit_usable <- function(v) {
  inherits(v, "hazard") && !isFALSE(v$fit$converged) &&
    isTRUE(is.finite(v$fit$objective))
}

#' hazard()'s limit layer: fit, check, restart or refuse
#'
#' @param cl `match.call()` of the user's `hazard()` call.
#' @param env Its `parent.frame()`.
#' @param control The call's evaluated `control`.
#' @return A `hazard` object; or an error is raised.
#' @keywords internal
#' @noRd
.hzr_limit_fit <- function(cl, env, control) {
  restart <- control$limit_restart %||% TRUE
  if (!(isTRUE(restart) || isFALSE(restart))) {
    stop("control$limit_restart must be TRUE or FALSE.", call. = FALSE)
  }
  first <- .hzr_limit_capture(cl, env)
  emit <- function(r) for (w in r$warnings) warning(w)
  if (inherits(first$value, "condition")) {
    emit(first)
    stop(first$value)
  }
  trigger <- .hzr_limit_trigger(first$value)
  if (is.null(trigger)) {
    emit(first)
    return(first$value)
  }
  # On a refusal the drifted fit's own warnings (weak, boundary) are raised
  # first: they describe the drift the refusal names.
  if (!restart) {
    emit(first)
    stop(.hzr_limit_refusal(trigger, first$value$fit$objective, 0L,
                            disabled = TRUE))
  }
  fit0 <- first$value
  prep <- .hzr_evaluate_prepare(fit0)
  # The user's start: their theta if they gave one, else the phase specs'
  # default, assembled as the optimizer assembles it.
  start <- if (!is.null(cl$theta)) {
    stats::setNames(eval(cl$theta, env), prep$names)
  } else {
    stats::setNames(unlist(lapply(names(prep$phases), function(nm) {
      .hzr_phase_start(prep$phases[[nm]],
                       n_covariates = prep$covariate_counts[[nm]])
    }), use.names = FALSE), prep$names)
  }
  fixed <- prep$names[as.logical(fit0$fit$fixed_mask %||%
                                   logical(length(prep$names)))]
  grid <- .hzr_limit_grid(start, trigger, prep$phases,
                          event_times = fit0$data$time[fit0$data$status == 1],
                          fixed = fixed)
  ctrl <- control
  ctrl$n_starts <- 1L   # the grid IS the multi-start; 9 x n_starts is waste
  runs <- lapply(grid, function(s) {
    cl_i <- cl
    cl_i$theta <- s
    cl_i$control <- ctrl
    .hzr_limit_capture(cl_i, env)
  })
  rows <- data.frame(
    start = seq_along(runs),
    objective = vapply(runs, function(r) {
      if (.hzr_limit_usable(r$value)) r$value$fit$objective else NA_real_
    }, numeric(1)),
    triggered = vapply(runs, function(r) {
      .hzr_limit_usable(r$value) && !is.null(.hzr_limit_trigger(r$value))
    }, logical(1)),
    error = vapply(runs, function(r) {
      if (inherits(r$value, "condition")) conditionMessage(r$value) else NA_character_
    }, character(1)),
    stringsAsFactors = FALSE)
  # The winner is the best objective OVERALL, the original included. A clean
  # fit that is worse than a fit at a limit is not a rescue: there the data
  # support the limit, and the fit refuses (amended after Task 1: on #448's
  # avc fixture the best clean start was -208.31, against -21.66 at the step
  # and -78.15 returned today).
  best_any <- suppressWarnings(max(c(fit0$fit$objective, rows$objective),
                                   na.rm = TRUE))
  ok <- which(is.finite(rows$objective) & !rows$triggered)
  win <- if (length(ok)) ok[which.max(rows$objective[ok])] else NA_integer_
  # Relative tolerance: grid starts that reach the same optimum differ in
  # the last digits, and a tie must not read as "a limit fit was better".
  if (is.na(win) ||
        best_any - rows$objective[[win]] > 1e-6 * max(1, abs(best_any))) {
    emit(first)
    stop(.hzr_limit_refusal(trigger, best_any, length(runs)))
  }
  out <- runs[[win]]$value
  emit(runs[[win]])
  # The user's call, not the grid start's, so print() shows what they wrote
  # and hzr_bootstrap() refits from their start (through this layer again).
  out$call <- fit0$call
  out$call_env <- fit0$call_env
  out$fit$limit <- list(
    trigger = trigger,
    original = list(objective = fit0$fit$objective, theta = fit0$fit$theta,
                    weak = fit0$fit$weak, boundary = fit0$fit$boundary),
    grid = rows, winner = win)
  warning(structure(
    class = c("hzr_limit_restarted", "warning", "condition"),
    list(message = paste0(
      .hzr_limit_describe(trigger), " (objective ",
      format(fit0$fit$objective, digits = 8), "). Restarted from ",
      length(runs), " starts; returned start ", win, " (objective ",
      format(out$fit$objective, digits = 8), "). See fit$fit$limit."),
      call = NULL)))
  out
}
```

- [ ] **Step 5: Wire the layer into `hazard()`**

In `R/hazard_api.R`, make this the first statement of `hazard()`, before
`objective <- match.arg(objective)`:

```r
  # #683: a multiphase fit is checked for a shape parameter at a limit, and
  # restarted or refused. The layer re-evaluates this very call (inner),
  # so everything below runs unchanged for it.
  if (isTRUE(fit) && identical(unname(dist), "multiphase") &&
        !.hzr_limit_state$inner) {
    return(.hzr_limit_fit(match.call(), parent.frame(), control))
  }
```

Directly before `class(obj) <- "hazard"` (around line 2074), add the following. It's
written as a single-bracket list assignment because `obj$fit$limit <- NULL` would *delete*
the element, and the clean-fit test asserts that the name is present:

```r
  # #683: NA not examined; NULL examined and clean. The limit layer replaces
  # it with a record when it restarts.
  obj$fit["limit"] <- list(
    if (fit_ran && identical(dist, "multiphase")) NULL else NA)
```

- [ ] **Step 6: Document**

In `hazard()`'s roxygen `@param control`, add an item for `limit_restart`:

> `limit_restart` (multiphase, default `TRUE`): when a phase's shape parameter runs to a
> limit of its parameterisation (`t_half` toward 0 or past the data, `nu` toward 0 as a
> step, or a flat direction through a shape), refit from a fixed grid of starts and return
> the best fit without a parameter at a limit. The original drift is recorded in
> `fit$fit$limit`. When no start succeeds, or with `FALSE`, `hazard()` stops with an error
> of class `hzr_limit_refusal` naming the parameter. The drifted estimate is never
> returned.

Add `fit$fit$limit` to the `@return` list, with the three values.

- [ ] **Step 7: Run and verify they pass**

Run: the Step 2 command. Expected: all PASS. If the #729 n_starts = 5 test fails on the
objective, do not loosen the tolerance. Print `fit$fit$limit$grid` and take it to the
maintainer. The grid itself is a design decision.

- [ ] **Step 8: Run the whole suite**

Run: `NOT_CRAN=true Rscript -e 'devtools::document(); r <- as.data.frame(devtools::test(reporter = "silent")); saveRDS(r, Sys.getenv("SCRATCH_RDS")); cat(sum(r$passed), sum(r$failed), "\n")'`

Use `SCRATCH_RDS=$SCRATCH/suite-task5.rds`.
- Expected: the only failures are in the #444 and #448 files that Task 1 flagged, and in
  `test-sas-parity.R` if its warm start now triggers. Every one of those is handled in
  Task 8.
- **Any other failure is a finding.** Stop, and do not edit that test to pass.

- [ ] **Step 9: Commit**

```bash
git add R/limit-restart.R R/hazard_api.R R/likelihood-multiphase.R man/ tests/testthat/test-limit-restart.R
git commit -m "hazard(): restart or refuse a multiphase fit at a parameter limit (#683, #729)"
```

---

### Task 6: `print()` and `summary()` always show the limit check

**Files:**
- Modify: `R/limit-restart.R` (add `.hzr_limit_line()`)
- Modify: `R/hazard_api.R`:
  - `print.hazard()`: before `cat(.hzr_format_not_done(...))`, about line 3054;
  - `print.summary.hazard()`: at the matching place;
  - `summary.hazard()`: carry `limit` into the summary object.
- Test: `tests/testthat/test-limit-trigger.R` (append; fast)

**Interfaces:**
- Produces: `.hzr_limit_line(limit)`, which returns one string beginning
  `"  Limit check: "`.

- [ ] **Step 1: Failing tests**

```r
test_that("the limit line has a sentence for each state", {
  expect_identical(.hzr_limit_line(NULL), "  Limit check: none")
  expect_identical(.hzr_limit_line(NA), "  Limit check: not run")
  lim <- list(trigger = list(list(phase = "early", parameter = "nu",
                                  limit = "0 (a step at t_half)", source = "weak")),
              grid = data.frame(start = 1:2, objective = c(-5, -4)), winner = 2L)
  expect_identical(.hzr_limit_line(lim),
    "  Limit check: restarted (early.nu ran to 0 (a step at t_half)); start 2 of 2 kept")
})

test_that("print() shows the line for a multiphase object", {
  obj <- hazard(time = c(1, 2, 3), status = c(1, 0, 1), dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = FALSE)
  expect_match(paste(utils::capture.output(print(obj)), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
})
```

- [ ] **Step 2: Run and verify they fail.**

- [ ] **Step 3: Implement**

```r
#' The always-printed limit line (a line shown only on bad news cannot be
#' told from one its author forgot to write; R/degraded.R)
#' @keywords internal
#' @noRd
.hzr_limit_line <- function(limit) {
  paste0("  Limit check: ", if (is.null(limit)) {
    "none"
  } else if (!is.list(limit)) {
    "not run"
  } else {
    paste0("restarted (", .hzr_limit_describe(limit$trigger), "); start ",
           limit$winner, " of ", nrow(limit$grid), " kept")
  })
}
```

In `print.hazard()`, before the `.hzr_format_not_done` line, add (multiphase only, since
single-distribution fits have no limit layer):

```r
  if (identical(x$spec$dist, "multiphase")) cat(.hzr_limit_line(x$fit$limit), "\n", sep = "")
```

Do the same in `print.summary.hazard()`, reading the field `summary.hazard()` carries.
Add `limit = object$fit$limit` to the list that `summary.hazard()` builds.

- [ ] **Step 4: Run, verify pass, and run `test-print*.R` and `test-summary*.R` for regressions.**

Print snapshots that change gain exactly one line. Update each by hand and name it in the
PR.

- [ ] **Step 5: Commit**

```bash
git add R/ tests/
git commit -m "print/summary: always show the limit check (#683)"
```

---

### Task 7: Stepwise and bootstrap read a refusal as a counted reason

Both callers already catch refit errors and tally `conditionMessage()`:
- stepwise in `.hzr_refit_failure_reason()`, `R/stepwise-step.R:704`;
- bootstrap at `R/diagnostics.R:3008`.

A refusal's long message would become a unique tally key per replicate. This task gives it
one short, stable reason. It also proves each route by executing it.

**Files:**
- Modify: `R/limit-restart.R` (add `.hzr_limit_reason()`)
- Modify: `R/stepwise-step.R:704-712`
- Modify: `R/diagnostics.R:3008-3013`
- Test: `tests/testthat/test-limit-callers.R` (file-level `skip_on_cran()`)

**Interfaces:**
- Produces: `.hzr_limit_reason(e)`, which returns
  `"limit: <phase>.<parameter> ran to <limit>"`.

- [ ] **Step 1: Failing tests**

```r
# tests/testthat/test-limit-callers.R
skip_on_cran()
# #683: a refusal reached through a refit is a counted, named failure; it never
# reads as "no candidate good enough" or as a clean replicate.

test_that("a refusal has one short reason", {
  e <- .hzr_limit_refusal(list(list(phase = "early", parameter = "nu",
                                    limit = "0 (a step at t_half)", source = "weak")),
                          -10, 9L)
  expect_identical(.hzr_limit_reason(e), "limit: early.nu ran to 0 (a step at t_half)")
  expect_identical(.hzr_refit_failure_reason(e), .hzr_limit_reason(e))
})
```

For each route, add a test that forces a refusal by mocking the layer:

```r
with_refusal <- function(code) {
  testthat::local_mocked_bindings(
    .hzr_limit_trigger = function(object) {
      list(list(phase = names(object$spec$phases)[[1]], parameter = "nu",
                limit = "0 (a step at t_half)", source = "weak"))
    },
    .hzr_limit_grid = function(...) list(),
    .package = "TemporalHazard")
  code
}
```

Every refit then refuses. Using the fixtures in `helper-stepwise.R`, which have a small
multiphase base fit (built **outside** `with_refusal()`, so the base exists):
- **Stepwise forward:** `hzr_stepwise(..., direction = "forward")`. Assert:
  - a warning matching `"candidate refit failed for .*: limit: "`;
  - every candidate in `$failures`;
  - `$failure_reasons` all starting `"limit: "`;
  - no candidate entered.
- **Bootstrap:** `hzr_bootstrap(base, n_boot = 3L)`. Assert:
  - `n_failed == 3`;
  - the failure tally has the single name `"limit: <phase>.nu ran to 0 (a step at t_half)"`
    with count 3.
- **Stepwise drop:** with a base that has one term and `direction = "backward"`, and the
  drop accepted, assert that the step stops with class `hzr_limit_refusal` or records it
  as a failed drop.
  - Read `R/stepwise-step.R:1041` first, and assert what the code does.
  - If it silently keeps the model, that is a finding. Report it, and do not encode it.

- [ ] **Step 2: Run and verify they fail.**

- [ ] **Step 3: Implement**

In `R/limit-restart.R`:

```r
#' One short, stable reason for a refusal, for failure tallies
#' @keywords internal
#' @noRd
.hzr_limit_reason <- function(e) {
  paste0("limit: ", .hzr_limit_describe(e$trigger))
}
```

In `.hzr_refit_failure_reason()` (`R/stepwise-step.R`), make the first branch:

```r
  if (inherits(refit, "hzr_limit_refusal")) {
    return(.hzr_limit_reason(refit))
  }
```

In `hzr_bootstrap()` (`R/diagnostics.R`), in the `failure <- if (inherits(boot_fit, "error"))`
branch, before `msg <- conditionMessage(boot_fit)`:

```r
      if (inherits(boot_fit, "hzr_limit_refusal")) {
        .hzr_limit_reason(boot_fit)
      } else {
        msg <- conditionMessage(boot_fit)
        if (nzchar(msg)) msg else "error with an empty message"
      }
```

Replace the two existing lines with that block.

- [ ] **Step 4: Run, and verify they pass.** Then run the stepwise and bootstrap files.

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all(); for (f in list.files("tests/testthat", "stepwise|bootstrap", full.names = TRUE)) testthat::test_file(f)'`

Expected: 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/ tests/testthat/test-limit-callers.R
git commit -m "Stepwise and bootstrap tally a limit refusal by one reason (#683)"
```

---

### Task 8: Migrate the existing boundary tests, the blind-start parity test, and the false-positive gate

**Files:**
- Modify: `tests/testthat/test-boundary-phase-discontinuity.R`,
  `tests/testthat/test-unbounded-phase.R`, `tests/testthat/test-boundary-step-detail.R`
  (only the tests that Task 5 Step 8 showed failing)
- Modify: `tests/testthat/test-sas-parity.R:205-300`
- Create (not committed): `$SCRATCH/fp-gate.R`

- [ ] **Step 1: Migrate each failing boundary test by what it asserted**

The rule is that **the record is still asserted, now where it lives**. It was on the
returned fit. Now it is on `fit$fit$limit$original$boundary` when restarted, or on the
refusal's `trigger`.

For each failing test:
- **If the fit now restarts:** wrap it in `expect_warning(..., class = "hzr_limit_restarted")`.
  Assert the original record via `fit$fit$limit$original$boundary`, with the same fields
  and the same values the test asserted before.
- **If it now refuses:** `expect_error(..., class = "hzr_limit_refusal")`, and assert the
  `trigger` source equals the mechanism the test was about.
- **If the test is about the record's construction** (detail text, magnitude), call
  `.hzr_boundary_check_impl()` directly on the fitted theta, as the new file in Task 2
  does. Keep every value assertion.

Do not delete an assertion. A pass count that drops in these files must be explained per
file in the PR.

**Stop if Task 1's checkpoint was not cleared for these fixtures.**

- [ ] **Step 2: Blind-start the `hm.death.patient` parity test**

In `test-sas-parity.R`:
- Delete the "warm-started from the SAS estimates ... Blind multi-start lands in a
  different (shrunk-beta) basin" comment paragraph (lines 214-218).
- Keep the existing warm-start test unchanged. It checks likelihood parity at SAS's point.
- Add a new test beside it:

```r
test_that("hm.death.patient: a blind start reaches SAS's optimum (#729)", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  ref <- .hzr_parse_sas_lst(file.path(dir, "hm.death.patient.lst"))$fits[[1]]
  data(valves, package = "TemporalHazard", envir = environment())
  d <- valves
  d$nyha[is.na(d$nyha)] <- mean(d$nyha, na.rm = TRUE)
  d$mv_nyha <- d$mitral * d$nyha
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
    dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 10, nu = 0.2, m = 1, fixed = "m",
                    formula = ~ age_cop + nyha + mv_nyha + double_ + ao_pinc),
                  constant = hzr_phase("constant", formula = ~ age_cop + black + i_path)),
    fit = TRUE, control = list(n_starts = 1L, maxit = 2000, conserve = TRUE)))
  expect_equal(fit$fit$objective / ref$loglik, 1, tolerance = 1e-5)
  expect_true(is.list(fit$fit$limit))  # it got there by restarting, on the record
})
```

- [ ] **Step 3: The false-positive gate**

Write `$SCRATCH/fp-gate.R`. It runs every multiphase parity fit with the layer in place
and counts restarts:

```r
suppressMessages(devtools::load_all(quiet = TRUE, helpers = TRUE))
fired <- list()
trace(".hzr_limit_fit", exit = quote({
  r <- returnValue(default = NULL)
  if (is.list(r$fit$limit)) fired[[length(fired) + 1L]] <<- list(
    call = deparse(cl)[1], trigger = .hzr_limit_describe(r$fit$limit$trigger))
}), where = asNamespace("TemporalHazard"), print = FALSE)
files <- list.files("tests/testthat", "parity", full.names = TRUE)
stopifnot(length(files) > 0)   # a zero needs a known positive
res <- lapply(files, function(f) as.data.frame(testthat::test_file(f, reporter = "silent")))
cat("parity files:", length(files), " tests run:", sum(vapply(res, nrow, 1L)),
    " skipped:", sum(vapply(res, function(r) sum(r$skipped), 1)), "\n")
cat("restarts fired:", length(fired), "\n")
for (x in fired) cat("-", x$trigger, "|", x$call, "\n")
```

Run it with `NOT_CRAN=true` and `HAZARD_REPO` / `HAZARD_EXAMPLES_DIR` pointing at the
local checkouts. **Assert the skip count is 0 for the SAS-fixture files** before reading
the restart count. A skipped parity file reads as "0 fired".

Expected: **0 restarts fired** in the parity files, apart from the new #729 blind-start
test, which must fire exactly once.
- Any other firing is a false positive by the spec's gate. Stop and report it with its
  trigger.
- Then run the same counter over the whole suite, and list every firing with its file in
  the PR.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/
git commit -m "Migrate boundary tests to the limit layer; blind-start hm.death.patient (#683, #729)"
```

---

### Task 9: Definition of done, NEWS, review, PR

- [ ] **Step 1: NEWS fragment.** Create `news/design-limit-restart-683.md`:

```markdown
* `hazard()` no longer returns a multiphase fit whose shape parameter ran to a
  limit of its parameterisation (`t_half` toward 0 or past the last observed
  time, `nu` toward 0 as a step, or a flat direction through a shape). Such a
  fit could report `converged = TRUE` at a drifted estimate (#683). It now
  refits from a fixed grid of starts and returns the best fit with no
  parameter at a limit, recording the drift in `fit$fit$limit` and warning
  (class `hzr_limit_restarted`); when no start succeeds it stops with an error
  of class `hzr_limit_refusal` naming the parameter. `control$limit_restart =
  FALSE` refuses without restarting. `print()` and `summary()` show the
  outcome on a new "Limit check" line. On the `hm.death.patient` model a blind
  start that previously stopped 244 log-likelihood units below SAS now reaches
  SAS's optimum (#729).
```

Run: `Rscript -e 'print(spelling::spell_check_files("news/design-limit-restart-683.md", lang = "en-US"))'`

Expected: no words, or only words already in `inst/WORDLIST`.

- [ ] **Step 2: The local gate, in order**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
NOT_CRAN=true Rscript -e 'r <- as.data.frame(devtools::test(reporter = "silent")); saveRDS(r, "<SCRATCH>/suite-final.rds"); cat(sum(r$passed), "passed", sum(r$failed), "failed", sum(r$skipped), "skipped\n")'
Rscript -e 'spelling::spell_check_package(use_wordlist = TRUE)'
```

Expected: 0 lints and 0 failures.
- Compare `suite-final.rds` with `suite-base.rds` per file. Explain every file whose pass
  count fell.
- Stamp the total as "N passed at <sha>".

- [ ] **Step 3: The check, from a committed tree**

```bash
git status --short   # must be empty
rm -rf "$TMPDIR/tree" && mkdir -p "$TMPDIR/tree"
git archive HEAD | tar -x -C "$TMPDIR/tree"
cd "$TMPDIR" && R CMD build tree && R CMD check --as-cran TemporalHazard_*.tar.gz
```

Expected: 0 errors, 0 warnings, 0 notes. Read the per-step `[Ns]` timings, and record the
overall and tests times with the sha.

- [ ] **Step 4: r-reviewer.** Run the `r-reviewer` agent over `git diff origin/main...HEAD`.
  Verify each finding against the code before acting on it. Stop once a pass finds no
  silent wrong answer.

- [ ] **Step 5: Push the branch and open the PR.** Push `design/limit-restart-683` and open
  the PR with `gh pr create`.
  - The body says `Closes #683` and `Closes #729`.
  - It includes the Task 1 table, the planted-fault kills, the false-positive gate output,
    the per-file pass delta, and the check timings.
  - Then stop. The maintainer merges.
