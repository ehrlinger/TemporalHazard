# SAS Translator Runnable-Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `hzr_translate_sas()` emit Quarto documents that actually render, and stop it reporting success over documents that cannot.

**Architecture:** The root fix is in `hazard()`, not the translator: `hazard()` gains `subset()`-style data-masked evaluation of its censoring arguments, which makes the currently-ignored `data =` argument meaningful. The translator then only needs to bind its fit, request an actual fit, scope the status chunk, and name the grid column what `predict.hazard()` requires. A render-simulation test harness is built **first** and is the acceptance test for everything after it.

**Tech Stack:** R (base `eval`/`substitute` NSE -- no rlang), roxygen2 8.1.0, testthat 3, lintr 3.4.0.

## Global Constraints

- Design spec: `inst/dev/PLAN-sas-translator-fix.md`. Read it before starting.
- Operational contract: `AGENTS.md`. Read it before starting.
- **No new package dependencies.** Data masking uses base `eval(substitute(x), data, parent.frame())`.
- Censoring status coding is `-1` left, `0` right, `1` event, `2` interval -- NOT `survival::Surv()`'s integers.
- `stats::` prefixing is house style for stats generics.
- No `browser()`, no bare `print()`, no `library()` inside `R/`.
- Version stays `1.2.2` in both `DESCRIPTION` and the `NEWS.md` top heading. Do not roll minor or major.
- Definition of done for every task: `Rscript -e 'devtools::document()'`, then `Rscript -e 'lintr::lint_package()'` (0 lints), then `Rscript -e 'devtools::test()'` (0 failures). Run them; do not infer them.
- Tests that take more than a second or two get `skip_on_cran()`.
- Branch is `fix/sas-translator-runnable`, already created off `dev` and carrying the spec commit. Never push to `main` or `dev`.

---

### Task 1: Render-simulation harness

The single assertion that catches every defect in #151 and #152: write the `.qmd`, parse it, evaluate its chunks in a **clean** environment holding only the input data. Build it first; it fails on every job today.

**Files:**
- Create: `tests/testthat/helper-render-sim.R`
- Test: `tests/testthat/test-sas-render-sim.R`

**Interfaces:**
- Produces: `render_sim(job, data = list())` -> `list(ok = logical, results = named character, env = environment)`. `results` holds `"ok"` or `"ERROR: <message>"` per chunk, in emission order. Every later task uses this.

- [ ] **Step 1: Write the harness**

```r
# tests/testthat/helper-render-sim.R
#
# Evaluate a translated job's emitted chunks the way rendering the .qmd
# would: in a clean environment holding only the input data. Asserting the
# SHAPE of emitted calls cannot distinguish a document that runs from one
# that errors -- every defect in #151/#152 was invisible to shape assertions
# and visible to this in one line.

render_sim <- function(job, data = list()) {
  env <- new.env(parent = globalenv())
  for (nm in names(data)) assign(nm, data[[nm]], envir = env)
  res <- character(0)
  for (nm in names(job$calls)) {
    out <- tryCatch({
      eval(job$calls[[nm]], env)
      "ok"
    }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
    res[[nm]] <- out
  }
  list(ok = all(res == "ok"), results = res, env = env)
}
```

- [ ] **Step 2: Write a test that the harness itself works**

```r
# tests/testthat/test-sas-render-sim.R
test_that("render_sim reports per-chunk success and failure", {
  job <- list(calls = list(
    good = quote(x <- 1 + 1),
    bad  = quote(stop("boom")),
    uses = quote(y <- x * 2)
  ))
  got <- render_sim(job)
  expect_false(got$ok)
  expect_equal(unname(got$results[["good"]]), "ok")
  expect_match(got$results[["bad"]], "^ERROR: boom")
  expect_equal(got$env$y, 4)
})

test_that("render_sim binds supplied data and nothing else", {
  job <- list(calls = list(a = quote(n <- nrow(D)), b = quote(z <- MISSING_VAR)))
  got <- render_sim(job, data = list(D = data.frame(x = 1:3)))
  expect_equal(got$env$n, 3L)
  expect_match(got$results[["b"]], "object 'MISSING_VAR' not found")
})
```

- [ ] **Step 3: Run and verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: `FAIL 0 | PASS 6`

- [ ] **Step 4: Add the failing acceptance test for a real job**

```r
test_that("a minimal translated job renders end to end", {
  skip_on_cran()
  set.seed(1)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_s3_class(got$env$fit, "hazard")
  expect_true(isTRUE(got$env$fit$fit$converged))
})
```

- [ ] **Step 5: Run it and confirm it fails for the expected reasons**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: FAIL 1 -- the `info` string shows `fit  ERROR: object 'INT_DEAD' not found`, and `got$env$fit` does not exist.

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/helper-render-sim.R tests/testthat/test-sas-render-sim.R
git commit -m "test: add a render simulation harness for translated jobs

Evaluates emitted chunks in a clean environment holding only the input data,
which is what rendering the .qmd does. Shape assertions on emitted calls
cannot tell a document that runs from one that errors; this can. The
acceptance test fails today, as it should: object 'INT_DEAD' not found.

Refs #151"
```

---

### Task 2: Data masking in `hazard()`

**Files:**
- Modify: `R/hazard_api.R` (insert after the formula-dispatch block ending at ~line 374, before the `is.null(time) || is.null(status)` check)
- Test: `tests/testthat/test-hazard-data-masking.R` (create)

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `hazard(data = <df>, time = <bare column>, status = <bare column>, ...)` resolves columns from `data`. Tasks 3-8 rely on this.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-hazard-data-masking.R
#
# hazard()'s `data` argument was accepted and silently ignored on the vector
# path: only the formula path consulted it. That is why the SAS translator's
# emitted hazard(data = AVCS, time = INT_DEAD, ...) failed with
# "object 'INT_DEAD' not found" (#151).

make_df <- function(n = 100) {
  set.seed(4)
  data.frame(tt = stats::rexp(n, 0.3),
             ev = rep(c(1, 0), length.out = n),
             entry = rep(0, n))
}

test_that("bare column names resolve against `data` on the vector path", {
  df <- make_df()
  fit <- hazard(data = df, time = tt, status = ev, dist = "weibull",
                theta = c(0.5, 1), fit = TRUE)
  expect_s3_class(fit, "hazard")
  expect_equal(fit$data$time, df$tt)
  expect_equal(fit$data$status, df$ev)
})

test_that("df$col, a local vector and a bare name all still work", {
  df <- make_df()
  local_time <- df$tt
  a <- hazard(data = df, time = df$tt, status = df$ev, dist = "weibull",
              theta = c(0.5, 1), fit = TRUE)
  b <- hazard(data = df, time = local_time, status = ev, dist = "weibull",
              theta = c(0.5, 1), fit = TRUE)
  d <- hazard(time = df$tt, status = df$ev, dist = "weibull",
              theta = c(0.5, 1), fit = TRUE)
  expect_equal(a$fit$theta, d$fit$theta)
  expect_equal(b$fit$theta, d$fit$theta)
})

test_that("a column masks a same-named caller variable", {
  df <- make_df()
  tt <- rep(999, nrow(df))
  fit <- hazard(data = df, time = tt, status = ev, dist = "weibull",
                theta = c(0.5, 1), fit = TRUE)
  expect_equal(fit$data$time, df$tt)
  expect_false(any(fit$data$time == 999))
})

test_that("time_lower and weights are masked too", {
  df <- make_df()
  df$w <- rep(c(1, 2), length.out = nrow(df))
  fit <- hazard(data = df, time = tt, status = ev, time_lower = entry,
                weights = w, dist = "weibull", theta = c(0.5, 1), fit = TRUE)
  expect_equal(fit$data$weights, df$w)
})

test_that("the formula path is unaffected", {
  df <- make_df()
  fit <- hazard(survival::Surv(tt, ev) ~ 1, data = df, dist = "weibull",
                theta = c(0.5, 1), fit = TRUE)
  expect_s3_class(fit, "hazard")
  expect_equal(fit$data$time, df$tt)
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hazard-data-masking.R")'`
Expected: FAIL with `object 'tt' not found`.

- [ ] **Step 3: Implement the masking**

Insert immediately after the formula-dispatch `}` and before `# After formula dispatch, require time and status`:

```r
  # Data masking on the vector path. `data` used to be consulted only by the
  # formula path, so hazard(data = df, time = tt) failed with "object 'tt'
  # not found" while looking like it should work -- the defect behind the SAS
  # translator's unrenderable documents (#151). This is the base-R idiom
  # subset()/transform()/with() use: evaluate the argument expression with
  # `data` as the environment and the caller's frame as its parent, so a
  # column wins, and anything that is not a column (df$col, a local vector,
  # a literal) falls through to the caller unchanged.
  if (is.null(formula) && !is.null(data)) {
    if (!is.data.frame(data) && !is.list(data)) {
      stop("'data' must be a data frame or a list.", call. = FALSE)
    }
    mask_env <- parent.frame()
    time <- eval(substitute(time), data, mask_env)
    status <- eval(substitute(status), data, mask_env)
    time_lower <- eval(substitute(time_lower), data, mask_env)
    time_upper <- eval(substitute(time_upper), data, mask_env)
    weights <- eval(substitute(weights), data, mask_env)
  }
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hazard-data-masking.R")'`
Expected: PASS, 0 failures.

- [ ] **Step 5: Verify nothing that re-evaluates a stored call regressed**

`AGENTS.md` warns that anything rewriting or resampling a stored call has to handle both interfaces. Masking changes what a stored call means, so check the resamplers explicitly rather than assuming.

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-bootstrap.R"); testthat::test_file("tests/testthat/test-stepwise.R")'`
Expected: 0 failures. If a bootstrap replicate now fails to resolve a column, the resampler is rewriting `time`/`status` without carrying `data`; fix by having it also substitute the resampled `data`, and record what you changed in the commit message.

- [ ] **Step 6: Update the `@param data` documentation**

In `R/hazard_api.R`, the roxygen for `@param data` must state that it is now consulted on both paths:

```r
#' @param data Optional data frame. On the formula path it supplies the model
#'   frame. On the vector path `time`, `status`, `time_lower`, `time_upper`
#'   and `weights` are evaluated in its scope, the way [base::subset()] and
#'   [base::transform()] do: a bare column name resolves to that column, and
#'   anything that is not a column (`df$col`, a local vector, a literal)
#'   falls through to the calling environment. A column of the same name as
#'   a caller variable wins.
```

- [ ] **Step 7: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/hazard_api.R man/hazard.Rd tests/testthat/test-hazard-data-masking.R
git commit -m "feat: evaluate hazard()'s censoring arguments in data's scope

hazard(data = df, time = tt) failed with 'object tt not found': `data` was
consulted only by the formula path, while the vector path accepted and
ignored it. That is the root cause of the SAS translator emitting documents
that cannot render (#151), and it bites anyone who passes data = by hand.

Uses the base subset()/transform() idiom, eval(substitute(x), data,
parent.frame()) -- no new dependency. Columns win; df\$col, local vectors and
literals fall through to the caller unchanged. Purely additive: with
data = NULL nothing changes, and the formula path is untouched.

Refs #151"
```

---

### Task 3: Bind the fit and request an actual fit

**Files:**
- Modify: `R/translate-sas.R` (the `calls[[fit_slot]] <- r$call` assignment, ~line 147)
- Modify: `R/sas-parse-job.R` (the `args` list built at ~lines 285-292)
- Test: `tests/testthat/test-sas-render-sim.R` (the acceptance test from Task 1)

**Interfaces:**
- Consumes: `render_sim()` from Task 1; data masking from Task 2.
- Produces: every fit chunk is `<name> <- hazard(..., fit = TRUE)` where `<name>` is the resolved slot name (`fit`, `fit_2`, ...). Tasks 5, 8 and 10 reference those names.

- [ ] **Step 1: Confirm the Task 1 acceptance test still fails, and how**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: FAIL. With Task 2 in place the `fit` chunk now evaluates, so the failure moves to `expect_s3_class(got$env$fit, "hazard")` -- `object 'fit' not found` -- because the chunk is a bare call that binds nothing.

- [ ] **Step 2: Emit `fit = TRUE`**

In `R/sas-parse-job.R`, immediately after `if (!is.null(cens$time_lower)) args$time_lower <- cens$time_lower`, add:

```r
  # Without fit = TRUE the emitted call returns an unfitted object: converged
  # is NA, objective is NA, and theta holds the SAS *starting* values, while
  # print.hazard() shows a populated summary that says none of that (#151).
  args$fit <- TRUE
```

- [ ] **Step 3: Bind the fit to its slot name**

In `R/translate-sas.R`, replace `calls[[fit_slot]] <- r$call` with:

```r
      # Bind the fit: predict() chunks reference the fit by its slot name, and
      # a bare hazard(...) call binds nothing, so those chunks failed with
      # "object 'fit' not found" -- or worse, silently used an unrelated
      # object of that name already in the rendering session (#151).
      calls[[fit_slot]] <- call("<-", as.name(fit_slot), r$call)
```

- [ ] **Step 4: Run the acceptance test to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: PASS -- `got$ok` TRUE, `got$env$fit` is a `hazard`, and `converged` is TRUE.

- [ ] **Step 5: Fix the tests that encoded the old shape**

`tests/testthat/test-sas-render-qmd.R` builds fixtures as `calls = list(fit = quote(hazard(...)))`, which encodes the missing assignment. Update those fixtures to `quote(fit <- hazard(...))`. `tests/testthat/test-sas-translate-fits.R` evaluates `job$calls` in a loop and takes the loop's last value as the fit; with the binding in place it should read `env$fit` instead, and its `env <- list2env(as.list(AVCS)); env$AVCS <- AVCS` workaround should collapse to `env$AVCS <- AVCS` alone -- that workaround existed only because the columns were unresolvable.

- [ ] **Step 6: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/translate-sas.R R/sas-parse-job.R tests/testthat/
git commit -m "fix(sas): bind the emitted fit and request an actual fit

Two defects that together made every translated document unrunnable: the fit
chunk was a bare hazard(...) call binding nothing, so predict() chunks failed
with 'object fit not found' -- or silently used an unrelated object of that
name already in the session -- and fit = TRUE was never emitted, so the call
stopped at the SAS starting values with converged = NA.

test-sas-render-qmd.R's fixtures encoded the missing assignment, which is why
no assertion could see it; updated to the bound form.

Refs #151"
```

---

### Task 4: Scope the status chunk to its data

**Files:**
- Modify: `R/sas-parse-job.R` (the `status_call <- call("<-", cens$status_name, cens$status_expr)` line, ~line 289)
- Test: `tests/testthat/test-sas-render-sim.R`

**Interfaces:**
- Consumes: Tasks 1-3.
- Produces: the status chunk is `.hzr_status <- with(<data>, <expr>)`.

- [ ] **Step 1: Write the failing test**

```r
test_that("a censoring job's status chunk renders", {
  skip_on_cran()
  set.seed(2)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_true(any(got$env$.hzr_status == 2))
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: FAIL -- `status  ERROR: object 'DEAD' not found`. Masking covers the inside of `hazard()`, not a free-standing chunk.

- [ ] **Step 3: Wrap the status expression in `with()`**

Replace the `status_call` construction with:

```r
    # The status chunk is evaluated outside hazard(), so hazard()'s data
    # masking does not reach it: wrap it in with() so its bare column names
    # resolve. A plain local binding (not a new column) keeps the caller's
    # data frame unmutated, and hazard()'s mask falls through to the caller's
    # frame to find it.
    status_expr <- if (is.null(data_name)) {
      cens$status_expr
    } else {
      call("with", as.name(data_name), cens$status_expr)
    }
    status_call <- call("<-", cens$status_name, status_expr)
```

- [ ] **Step 4: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: PASS.

- [ ] **Step 5: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-render-sim.R
git commit -m "fix(sas): scope the emitted status chunk to its data frame

hazard()'s data masking covers its own arguments, not a free-standing chunk,
so .hzr_status <- ifelse(DEAD == 0 & ...) still failed on bare column names.
Wrapped in with(<data>, ...). Kept as a local binding rather than a new
column so the reader's data frame is not mutated; hazard()'s mask falls
through to the caller's frame to find it.

Refs #151"
```

---

### Task 5: Name the prediction grid column `time`

**Files:**
- Modify: `R/sas-parse-job.R` (`names(cl) <- c("", time_var)` in the log-grid branch ~line 546; `names(cl) <- c("", var)` in the explicit-DO branch ~line 598)
- Test: `tests/testthat/test-sas-render-sim.R`

**Interfaces:**
- Consumes: Tasks 1-4.
- Produces: emitted grids are `data.frame(time = <expr>)`.

- [ ] **Step 1: Write the failing test**

```r
test_that("a HAZARD + HAZPRED job renders, predictions included", {
  skip_on_cran()
  set.seed(3)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14 OUTHAZ=E.H;",
    "EVENT DEAD; TIME INT_DEAD; PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );",
    "DATA PGRID; DO MONTHS = 1 TO 12 BY 1; OUTPUT; END; RUN;",
    "%HAZPRED( PROC HAZPRED DATA=PGRID INHAZ=E.H OUT=P; TIME MONTHS; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_true("time" %in% names(got$env$PGRID))
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: FAIL -- `pred  ERROR: 'newdata' must contain a 'time' column for 'survival' predictions.`

- [ ] **Step 3: Rename the emitted column in both grid branches**

In the log-grid branch replace `names(cl) <- c("", time_var)` with:

```r
    # predict.hazard() requires a column literally named `time`
    # (R/hazard_api.R:1006). Naming the grid column after the SAS DO variable
    # produced a newdata that predict() rejects outright, so the "grids
    # resolve" coverage figure counted grids that could not be used (#151).
    names(cl) <- c("", "time")
```

and in the explicit-DO branch replace `names(cl) <- c("", var)` with the same `names(cl) <- c("", "time")`.

- [ ] **Step 4: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-render-sim.R")'`
Expected: PASS.

- [ ] **Step 5: Update the grid tests that assert the SAS variable name**

Grep for existing assertions on the grid column name and update them:

Run: `grep -rn 'MONTHS' tests/testthat/ | grep -v '\.sas'`
Update any that expect the grid's column to be `MONTHS`; the SAS `DO` variable name is still what the parser reads, but it is no longer what the emitted column is called.

- [ ] **Step 6: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/
git commit -m "fix(sas): name the emitted prediction grid column 'time'

The grid column was named after the SAS DO variable (MONTHS, ...), but
predict.hazard() requires a column literally named 'time' and rejects
anything else. Every 'resolved' grid therefore produced a newdata that
predict() refused -- so the 19/55 (35%) grid-resolution figure counted grids
that could not be used.

Refs #151"
```

---

### Task 6: Carry ICENSOR's event count as weights

**Files:**
- Modify: `R/sas-parse-job.R` (the `args` build, after `args$time_lower`)
- Test: `tests/testthat/test-sas-censoring.R`

**Interfaces:**
- Consumes: Tasks 1-5.
- Produces: `args$weights` is the ICENSOR count variable when the job has one.

Confirmed against the reference implementation (`~/Documents/GitHub/hazard/src/llike/setlik.c`): `c3w = c3 * weight`, `c1c2c3 = c1w + c2 + c3w`, `llike = -(c1c2c3)*(cumhaz - cumhst)`, `if (c3 > ZERO) llike += (c3w * lct)`. `C3` multiplies the contribution in both terms. `ICENSOR` appears in 64 of 93 production jobs.

- [ ] **Step 1: Write the failing test**

```r
test_that("ICENSOR's event count is emitted as weights", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]], txt)
  # setlik.c: C3 is a COUNT that multiplies the log-likelihood contribution
  # in both terms, so dropping it changes the likelihood (#154).
  expect_equal(got$call[["weights"]], as.name("C3"))
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-censoring.R")'`
Expected: FAIL -- `weights` is NULL.

- [ ] **Step 3: Emit the weights**

In `R/sas-parse-job.R`, after the `args$time_lower` assignment:

```r
  # ICENSOR's first operand is an event COUNT, not a flag. setlik.c multiplies
  # the log-likelihood contribution by it (c3w = c3 * weight, then both
  # -(c1w+c2+c3w)*(cumhaz-cumhst) and += c3w*lct), so discarding the
  # magnitude fits a 3-event row as a single observation and silently changes
  # the likelihood (#154).
  if (!is.null(cens$weights_name)) args$weights <- as.name(cens$weights_name)
```

Then in `.hzr_censor_spec()` set `weights_name` to the ICENSOR count variable's name when an `ICENSOR` statement is present, and `NULL` otherwise. Return it in the spec list alongside `status_name` / `status_expr` / `time_lower`.

- [ ] **Step 4: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-censoring.R")'`
Expected: PASS.

- [ ] **Step 5: Verify the weights actually reach the fit**

Add to `tests/testthat/test-sas-render-sim.R`, extending the Task 4 job:

```r
test_that("the ICENSOR count reaches the fitted object as weights", {
  skip_on_cran()
  set.seed(5)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- ifelse(seq_len(n) %% 5 == 0, 3, 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  # The decisive check: a count of 3 must not arrive as a weight of 1.
  expect_equal(got$env$fit$data$weights, AVCS$C3FLAG)
  expect_true(any(got$env$fit$data$weights == 3))
})
```

- [ ] **Step 6: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/
git commit -m "fix(sas): carry ICENSOR's event count through as weights

ICENSOR's first operand is an event COUNT. setlik.c multiplies the
log-likelihood contribution by it -- c3w = c3 * weight, then both
-(c1w+c2+c3w)*(cumhaz-cumhst) and += c3w*lct -- so discarding the magnitude
fitted a 3-event row as a single interval-censored observation and silently
changed the likelihood. Confirmed against the reference source, not from the
formula transcribed in our own roxygen.

ICENSOR appears in 64 of 93 PROC HAZARD jobs across the production studies,
so this is the common path.

Closes #154"
```

---

### Task 7: Refuse LCENSOR + ICENSOR loudly

**Files:**
- Modify: `R/sas-parse-job.R` (`.hzr_censor_spec()`)
- Test: `tests/testthat/test-sas-censoring.R`

**Interfaces:**
- Consumes: Tasks 1-6.
- Produces: an `untranslated` row with construct `"LCENSOR + ICENSOR"`, and a `stop()` chunk in the emitted document.

`hazard()`'s `time_lower` carries two meanings selected by status -- entry time for status 0/1, interval lower bound for status 2 -- so one row cannot express both. SAS carries three distinct times (`TIME`, `CTIME`, `STIME`) and subtracts `H(STIME)` for all rows. Supporting this needs a new `hazard()` argument, which is out of scope; 0 of 93 production jobs use the combination.

- [ ] **Step 1: Write the failing test**

```r
test_that("LCENSOR combined with ICENSOR is refused, not silently mis-fitted", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; LCENSOR ST; ICENSOR C3 = CT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_true(any(grepl("LCENSOR", got$untranslated$construct)))
  # time_lower cannot carry both an entry time and an interval lower bound,
  # so the job must not translate to a silently wrong fit (#155).
  expect_null(got$call[["time_lower"]])
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-censoring.R")'`
Expected: FAIL -- nothing is recorded and `time_lower` is emitted.

- [ ] **Step 3: Implement the refusal**

In `.hzr_censor_spec()`, when both an `LCENSOR` and an `ICENSOR` statement are present, record and refuse:

```r
  # hazard()'s time_lower has two meanings selected by status: entry time for
  # status 0/1, interval lower bound for status 2. One row cannot carry both,
  # so a left-truncated interval-censored subject would be fitted as at risk
  # from time 0 -- a converging, plausible, wrong answer. SAS carries three
  # distinct times (TIME, CTIME, STIME) and subtracts H(STIME) for every row.
  # Supporting this needs a new hazard() argument; refuse until it exists.
  if (has_lcensor && has_icensor) {
    # NOTE: the helper is .hzr_untranslated_frame(line, construct, detail),
    # called positionally -- see .hzr_selection_spec() for the existing usage.
    untr <- rbind(untr, .hzr_untranslated_frame(
      NA_integer_, "LCENSOR + ICENSOR",
      paste("left truncation combined with interval censoring needs a",
            "separate entry-time argument hazard() does not have (#155);",
            "translate this job by hand")
    ))
    return(list(status_name = NULL, status_expr = NULL, time_lower = NULL,
                weights_name = NULL, untranslated = untr, refused = TRUE))
  }
```

Then in the caller, when `cens$refused` is TRUE, emit a `stop()` chunk instead of a fit, so the document fails loudly rather than rendering over a wrong model.

- [ ] **Step 4: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-censoring.R")'`
Expected: PASS.

- [ ] **Step 5: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-censoring.R
git commit -m "fix(sas): refuse LCENSOR + ICENSOR instead of dropping truncation

time_lower carries two meanings selected by status -- entry time for status
0/1, interval lower bound for status 2 -- so a single row cannot express
both. The translator was overwriting the entry time with the interval bound,
fitting a left-truncated interval-censored subject as at risk from time 0: a
converging, plausible, wrong answer with no warning.

SAS carries three distinct times (TIME, CTIME, STIME) and subtracts H(STIME)
for every row including interval ones, so supporting this needs a new
hazard() argument. 0 of 93 PROC HAZARD jobs across the production studies use
the combination, so refuse loudly now and extend the API separately.

Refs #155"
```

---

### Task 8: Make SELECTION emit a callable `hzr_stepwise()`

**Files:**
- Modify: `R/sas-parse-job.R` (the `head <- quote(hzr_stepwise)` branch, ~lines 319-329)
- Modify: `R/translate-sas.R` (so the stepwise call can reference the bound fit)
- Test: `tests/testthat/test-sas-selection.R`, `tests/testthat/test-sas-render-sim.R`

**Interfaces:**
- Consumes: Tasks 1-3 (the fit binding provides the `fit =` argument).
- Produces: `hzr_stepwise(fit = <fitname>, scope = ~ <candidates>, direction =, slentry =, slstay =)`.

`SELECTION` appears in **36 of 93** production jobs -- the second most common construct after the censoring statements.

- [ ] **Step 1: Write the failing test**

```r
test_that("SELECTION emits a callable hzr_stepwise() with a real scope", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14 SELECTION=STEPWISE",
    "SLENTRY=0.05 SLSTAY=0.1;",
    "EVENT DEAD; TIME T; EARLY X1, X2, X3;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]], txt)
  expect_equal(got$call[[1L]], as.name("hzr_stepwise"))
  expect_false(is.null(got$call[["fit"]]))
  expect_equal(got$call[["scope"]], quote(~X1 + X2 + X3))
  # The candidate pool must NOT be baked into the phase formula as forced-in
  # terms -- that inverts the meaning of the SAS statement (#152).
  ph <- got$call[["phases"]]
  expect_false(grepl("X1", paste(deparse(ph), collapse = " ")))
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-selection.R")'`
Expected: FAIL -- `fit` and `scope` are both absent and `X1` appears in the phase formula.

- [ ] **Step 3: Implement**

In the `isTRUE(sel$stepwise)` branch: build `scope` as a one-sided formula from the phase-covariate names, remove those names from the phase formula so they are candidates rather than forced terms, and add a `fit` placeholder that `R/translate-sas.R` rewrites to the resolved fit slot name (the same placeholder-rewriting `.hzr_point_predict_at()` already does for `predict()`):

```r
      head <- quote(hzr_stepwise)
      # SAS's SELECTION screens the phase-statement variables; they are
      # candidates, not model terms. Baking them into the phase formula made
      # them forced-in, inverting the statement's meaning, and the emitted
      # call had neither `fit` nor `scope` so it could not run at all (#152).
      args$fit <- quote(fit)
      args$scope <- stats::reformulate(cand_names)
      args$direction <- sel$direction
```

`cand_names` and the forced-term removal come from **where the covariates are
gathered**, not from post-processing built phases. In `.hzr_parse_hazard()`,
`covars` is assembled at `R/sas-parse-job.R:217-263` as a list with `$early`,
`$constant` and `$late` elements, then passed to `.hzr_parse_parms(parms_ops,
covars = covars)` at line 272. Split the same way `.hzr_parse_phase_covars()`
does, and when the job is a stepwise SELECTION build the phases with **no**
covariates so they are candidates rather than forced terms:

```r
  # Under SELECTION the phase variables are the candidate pool, so they must
  # not reach .hzr_parse_parms() as phase formulas -- that is what made them
  # forced-in model terms (#152). Collect their names for `scope` instead.
  cand_names <- unique(unlist(lapply(covars, .hzr_parse_phase_covars),
                              use.names = FALSE))
  is_stepwise <- !is.null(sel_ops) && isTRUE(.hzr_selection_spec(sel_ops)$stepwise)
  parms <- .hzr_parse_parms(parms_ops,
                            covars = if (is_stepwise) list() else covars)
```

Order matters: this runs before line 272's existing `parms <-` call, which it
replaces. `.hzr_selection_spec()` is called twice as written -- hoist it to a
single local if the duplicate call is awkward, but do not change its
behaviour.

- [ ] **Step 4: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-selection.R")'`
Expected: PASS.

- [ ] **Step 5: Add a render test that executes the stepwise call**

```r
test_that("a SELECTION job renders end to end", {
  skip_on_cran()
  set.seed(6)
  n <- 300
  A <- data.frame(T = stats::rexp(n, 0.2), DEAD = rep(c(1, 0), length.out = n),
                  X1 = stats::rnorm(n), X2 = stats::rnorm(n),
                  X3 = stats::rnorm(n))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14 SELECTION=STEPWISE",
    "SLENTRY=0.05 SLSTAY=0.1;",
    "EVENT DEAD; TIME T; EARLY X1, X2, X3;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(A = A))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
})
```

- [ ] **Step 6: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-parse-job.R R/translate-sas.R tests/testthat/
git commit -m "fix(sas): emit a callable hzr_stepwise() for SELECTION

Two defects. The emitted call had neither fit nor scope -- hzr_stepwise()'s
signature is (fit, scope, data, ...) -- so it failed with 'argument fit is
missing'. And the EARLY candidate pool was baked into the phase formula as
forced-in model terms, inverting the meaning of the SAS statement: even once
callable it would have screened the wrong thing.

The old test asserted only that call[[1]] was the symbol hzr_stepwise, never
evaluating it. SELECTION appears in 36 of 93 PROC HAZARD jobs across the
production studies.

Closes #152"
```

---

### Task 9: Match SAS's log-grid step

**Files:**
- Modify: `R/sas-parse-job.R` (the log-grid branch, ~lines 528-546)
- Test: `tests/testthat/test-sas-parse-hazpred.R`

**Interfaces:**
- Consumes: Task 5 (the grid column is now `time`).
- Produces: log grids whose points match SAS's `INC = (5 + LN_MAX)/99.9`.

- [ ] **Step 1: Write the failing test**

```r
test_that("the log grid matches SAS's INC = (5 + LN_MAX)/99.9 step", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; MAX = 180; LN_MAX = LOG(MAX);",
    "INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  cl <- .hzr_parse_grid(txt)
  grid <- eval(cl)
  expect_equal(nrow(grid), 100L)
  # SAS's last point is exp(-5 + 99*(log(180)+5)/99.9) = 164.2, not 180.
  expect_equal(max(grid$time), exp(-5 + 99 * (log(180) + 5) / 99.9),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(max(grid$time), 180)))
})

test_that("a directly-assigned log bound is not mistaken for MAX", {
  txt <- .hzr_sas_normalise(paste(
    "DATA PGRID; LN_MAX = 5.2; INC = (5 + LN_MAX)/99.9;",
    "DO LN_TIME = -5 TO LN_MAX BY INC; MONTHS = EXP(LN_TIME);",
    "OUTPUT; END; RUN;"
  ))
  cl <- .hzr_parse_grid(txt)
  grid <- if (is.null(cl)) NULL else eval(cl)
  # Must not silently produce a grid ending at t = 5.2 by reading LN_MAX=5.2
  # as MAX=5.2 and taking its log (#153).
  expect_true(is.null(grid) || max(grid$time) > 100)
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-parse-hazpred.R")'`
Expected: FAIL -- max is 180, and the second test builds a grid ending at 5.2.

- [ ] **Step 3: Fix the step and the regex**

Replace the `inner <- bquote(exp(seq(...length.out = 100)))` construction with the explicit SAS step, and anchor the `MAX=` regex so it cannot match `LN_MAX=`:

```r
    hi_txt <- local({
      # Anchor on a non-word character before MAX so LN_MAX=5.2 is not read
      # as MAX=5.2 -- that produced a grid ending at t = 5.2 instead of 181,
      # with no error (#153).
      m <- regmatches(body, regexec("(^|[^A-Z0-9_])MAX *= *([0-9.]+)", body))[[1L]]
      if (length(m) > 2L) m[[3L]] else NA_character_
    })
```

```r
    # SAS writes INC = (5 + LN_MAX)/99.9 and steps DO LN_TIME = -5 TO LN_MAX
    # BY INC, which lands 100 points whose LAST is exp(-5 + 99*INC), NOT
    # LN_MAX. seq(length.out = 100) implies a /99 step, so only the first
    # point agreed and the divergence grew to ~9.6% by the last (#153).
    inner <- bquote(
      exp(.(str2lang(lo_txt)) +
            seq(0, 99) * ((5 + log(.(str2lang(hi_txt)))) / 99.9))
    )
```

- [ ] **Step 4: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-parse-hazpred.R")'`
Expected: PASS.

- [ ] **Step 5: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/sas-parse-job.R tests/testthat/test-sas-parse-hazpred.R
git commit -m "fix(sas): match SAS's log-grid step instead of seq(length.out)

SAS writes INC = (5 + LN_MAX)/99.9; seq(lo, log(hi), length.out = 100)
implies a /99 step. Both give 100 points and both start in the same place, so
the grids looked equivalent while only the first point agreed -- the last was
180 against SAS's 164.2, a ~9.6% divergence that grows monotonically along
the grid and is attributable to nothing when comparing against a SAS listing.

Also anchored the MAX= regex so LN_MAX=5.2 is no longer read as MAX=5.2,
which produced a grid ending at t = 5.2 instead of 181.

Closes #153"
```

---

### Task 10: Class `hzr_read_outhaz()` and give it a `predict()` method

**Files:**
- Modify: `R/read-outhaz.R`
- Create: `tests/testthat/test-predict-outhaz.R`
- Modify: `NAMESPACE` (via roxygen `@export` on the S3 method)

**Interfaces:**
- Consumes: Tasks 1-5.
- Produces: `hzr_read_outhaz()` returns an object of class `hzr_outhaz`; `predict.hzr_outhaz(object, newdata, type, se.fit, conf.type, ...)` mirrors `predict.hazard()`'s signature.

OUTHAZ datasets carry the phase structure: `_NAME_` rows include `G1FLAG`, `G3FLAG`, `FIXDEL0`, `FIXMNU1`, `FIXGE2`, `FIXGAE2`, then `DELTA/THALF/NU/M`, `TAU/GAMMA/ALPHA/ETA`, `E0/C0/L0`. The phase spec is reconstructible from the file alone.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-predict-outhaz.R
test_that("hzr_read_outhaz() returns a classed object", {
  f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
  skip_if(!nzchar(f), "fixture not installed")
  got <- hzr_read_outhaz(f)
  expect_s3_class(got, "hzr_outhaz")
})

test_that("predict() works on a loaded OUTHAZ fit", {
  f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
  skip_if(!nzchar(f), "fixture not installed")
  obj <- hzr_read_outhaz(f)
  got <- predict(obj, newdata = data.frame(time = c(1, 6, 12)),
                 type = "survival")
  expect_length(got, 3L)
  expect_true(all(got >= 0 & got <= 1))
  expect_true(all(diff(got) <= 0))
})

test_that("an all-fixed parameter set yields a 0x0 vcov, not NULL", {
  f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
  skip_if(!nzchar(f), "fixture not installed")
  obj <- hzr_read_outhaz(f)
  # Documented hollow-object shape: check dimensions, never is.null().
  expect_true(is.matrix(obj$vcov))
})
```

- [ ] **Step 2: Run it and verify it fails**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-predict-outhaz.R")'`
Expected: FAIL -- class is `list`, and `predict()` reports `no applicable method for 'predict' applied to an object of class "list"`.

- [ ] **Step 3: Class the return**

In `R/read-outhaz.R`, replace the final `list(...)` with:

```r
  structure(
    list(estimates = est[is_param], status = status, vcov = vcov,
         flags = flags),
    class = "hzr_outhaz"
  )
```

- [ ] **Step 4: Implement `predict.hzr_outhaz()`**

Reconstruct the phase spec from the flag rows and delegate to the existing machinery:

```r
#' Predictions from a fit loaded out of a SAS `OUTHAZ=` dataset
#'
#' Reconstructs the phase specification from the `G1FLAG` / `G3FLAG` /
#' `FIX*` rows the `OUTHAZ=` dataset carries, then predicts as
#' [predict.hazard()] does.
#'
#' @param object An `hzr_outhaz` object from [hzr_read_outhaz()].
#' @param newdata Data frame with a `time` column.
#' @param type One of `"survival"`, `"hazard"`, `"cumulative_hazard"`.
#' @param ... Passed to [predict.hazard()].
#' @return Numeric vector, or a data frame when `se.fit = TRUE`.
#' @export
predict.hzr_outhaz <- function(object, newdata, type = "survival", ...) {
  spec <- .hzr_outhaz_to_spec(object)
  predict(spec, newdata = newdata, type = type, ...)
}
```

Write `.hzr_outhaz_to_spec()` (internal, `@noRd`) to build a `hazard`-classed object from the estimates, vcov and flag rows. `G1FLAG` non-zero means an early (`cdf`) phase is present; `G3FLAG` non-zero means a late (`g3`) phase; the constant phase is present when `C0` has an estimate row.

- [ ] **Step 5: Run to verify it passes**

Run: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-predict-outhaz.R")'`
Expected: PASS.

- [ ] **Step 6: Add a render test for the librefs path**

```r
test_that("a librefs-resolved INHAZ job renders", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  file.copy(system.file("extdata", "outhaz-fixture.rds",
                        package = "TemporalHazard"),
            file.path(dir, "hzdeath.rds"))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c(
    "DATA PGRID; DO MONTHS = 1 TO 12 BY 1; OUTPUT; END; RUN;",
    "%HAZPRED( PROC HAZPRED DATA=PGRID INHAZ=EX.HZDEATH OUT=P; TIME MONTHS; );"
  ), f)
  job <- suppressWarnings(
    hzr_translate_sas(f, librefs = c(EX = file.path(dir, "hzdeath.rds"))))
  got <- render_sim(job)
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
})
```

- [ ] **Step 7: Run the full gate and commit**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
git add R/read-outhaz.R man/ NAMESPACE tests/testthat/test-predict-outhaz.R tests/testthat/test-sas-render-sim.R
git commit -m "feat: class hzr_read_outhaz() and add predict.hzr_outhaz()

hzr_read_outhaz() returned an unclassed list, so the librefs/INHAZ= path --
documented as the supported way to resolve an external fitted model -- emitted
predict() against something with no method, failing with 'no applicable method
for predict applied to an object of class list'. Nothing recorded it in
\$untranslated.

OUTHAZ datasets carry the phase structure (G1FLAG, G3FLAG, FIX* rows), so the
spec is reconstructible from the file alone.

Refs #151"
```

---

### Task 11: Upgrade the corpus test from "parses" to "renders"

**Files:**
- Modify: `tests/testthat/test-sas-translate-corpus.R`

**Interfaces:**
- Consumes: `render_sim()` from Task 1, and Tasks 2-10.

The current assertions -- `tokens_seen > 0`, `tokens_seen >= tokens_mapped`, `n_ok > 20` -- are all true over a corpus of documents that cannot render. That is what let this ship.

- [ ] **Step 1: Add the render assertion**

```r
test_that("public-corpus jobs that translate also render", {
  skip_on_cran()
  fs <- list.files(Sys.getenv("HZR_SAS_CORPUS", unset = ""),
                   pattern = "[.]sas$", full.names = TRUE, recursive = TRUE)
  skip_if(length(fs) == 0L, "no .sas corpus available")
  n_render <- 0L
  n_trans <- 0L
  for (f in fs) {
    job <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) NULL)
    if (is.null(job)) next
    n_trans <- n_trans + 1L
    # Render with no data bound: a job whose DATA= dataset is absent must
    # fail in its guard chunk and nowhere else. Any other chunk erroring is
    # a translator defect, not a missing-data artefact.
    sim <- render_sim(job)
    bad <- names(sim$results)[sim$results != "ok"]
    bad <- setdiff(bad, grep("^data", names(sim$results), value = TRUE))
    if (length(bad) == 0L) n_render <- n_render + 1L
  }
  # Assert coverage, not just a non-zero count: a threshold that cannot fail
  # is how the previous corpus test passed over unrenderable documents.
  expect_gt(n_trans, 20L)
  expect_gte(n_render, n_trans)
})
```

- [ ] **Step 2: Run it**

Run: `HZR_SAS_CORPUS=$HOME/Documents/GitHub/hazard NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-sas-translate-corpus.R")'`
Expected: PASS. If `n_render < n_trans`, print the failing filenames and fix the translator -- do not lower the threshold.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-sas-translate-corpus.R
git commit -m "test: assert corpus jobs render, not merely parse

The old assertions -- tokens_seen > 0, tokens_seen >= tokens_mapped,
n_ok > 20 -- were all true over a corpus of documents that could not render,
which is how #151 shipped. Now every job that translates must also evaluate
its chunks, with only the missing-data guard allowed to fail.

Refs #151"
```

---

### Task 12: Rewrite the release documentation and re-gate

**Files:**
- Modify: `NEWS.md`, `R/translate-sas.R` (the `@section Experimental`), `R/read-outhaz.R`, `vignettes/sas-to-r-migration.qmd`

- [ ] **Step 1: Rewrite the NEWS entry**

Replace the experimental note's gap list with what is now true. Keep the feature marked experimental -- the API and emitted format may still change, and #155 remains open -- but remove every claim that is no longer accurate. State plainly which SAS constructs translate and which are refused.

- [ ] **Step 2: Rewrite the vignette callout**

The `callout-warning` currently says the emitted document does not render. Replace with an accurate statement of the remaining limits: `LCENSOR` + `ICENSOR` refused (#155), `SET`-derived grids refused, coverage counts parsing.

- [ ] **Step 3: Run the full local gate**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'
Rscript -e 'devtools::test()'
Rscript -e 'spelling::spell_check_package(use_wordlist = TRUE)'
```
Expected: 0 lints, 0 test failures, no spelling findings. Add any new prose words to `inst/WORDLIST` and re-sort with `LC_ALL=C sort -u`.

- [ ] **Step 4: Run the release check from a clean export**

```bash
git status --short          # must be empty; the archive exports HEAD, not the working tree
mkdir -p "$TMPDIR/tree"
git archive HEAD | tar -x -C "$TMPDIR/tree"
R CMD build "$TMPDIR/tree"
R CMD check --as-cran TemporalHazard_1.2.2.tar.gz
tar tzf TemporalHazard_1.2.2.tar.gz | grep -iE 'CLAUDE|AGENTS|data-raw' || echo "clean"
```
Expected: `Status: OK`, 0/0/0, well under the 10-minute budget (3m06s at the last measurement), tarball under 5 MB, no developer files.

- [ ] **Step 5: Commit and open the PR**

```bash
git add NEWS.md R/ man/ vignettes/ inst/WORDLIST
git commit -m "docs: describe what the translator now does

Refs #151, #152, #153, #154, #155"
git push -u origin fix/sas-translator-runnable
gh pr create --base dev --title "fix(sas): make translated documents actually render"
```

---

## Self-review notes

- **Spec coverage:** §1 -> Task 2. §2 -> Tasks 3, 4, 5. §3 #152 -> Task 8, #153 -> Task 9, #155 -> Task 7. §2(e) -> Task 10. #154 -> Task 6. §7 -> Tasks 1 and 11. Release docs -> Task 12. No spec section is unimplemented.
- **Deliberately deferred:** a separate entry-time argument in `hazard()` (the #155 root cause) and the `cran-comments.md` rewrite; both are named in the spec's out-of-scope section.
- **Naming consistency:** `render_sim()` (Task 1) is used verbatim in Tasks 3-11. `hzr_outhaz` is the class in Task 10 and in the Task 10 render test. `weights_name` is introduced in Task 6 and referenced in Task 7's refusal return.
- **Riskiest task:** Task 2, because it changes the core entry point. Its Step 5 exists specifically to catch resampler regressions rather than assume there are none.
