# tests/testthat/test-limit-callers.R
skip_on_cran()
# #683: a refusal reached through a refit is a counted, named failure; it never
# reads as "no candidate good enough" or as a clean replicate.

limit_trigger_hard <- function(phase) {
  list(list(phase = phase, parameter = "nu", limit = "0 (a step at t_half)",
            source = "phase_discontinuity"))
}

# Every refit inside `code` meets a hard trigger and an empty restart grid, so
# the layer refuses. Bases are built OUTSIDE this, so they exist.
with_refusal <- function(code) {
  testthat::local_mocked_bindings(
    .hzr_limit_trigger = function(object) {
      limit_trigger_hard(names(object$spec$phases)[[1]])
    },
    .hzr_limit_grid = function(...) list(),
    .package = "TemporalHazard")
  force(code)
}

# The warnings a call raised, as a character vector, with the call's value.
collect_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

lc_data <- function(n = 300L, seed = 683L) {
  set.seed(seed)
  d <- data.frame(z = stats::rnorm(n), w = stats::rnorm(n), v = stats::rnorm(n))
  d$time <- stats::rexp(n) * exp(-0.6 * d$z)
  d$status <- 1L
  d
}
lc_base <- function(d, formula = ~ 1) {
  hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
         phases = list(constant = hzr_phase("constant", formula = formula)),
         fit = TRUE)
}

test_that("a refusal has one short reason", {
  e <- .hzr_limit_refusal(limit_trigger_hard("early"), -10, 9L)
  expect_identical(.hzr_limit_reason(e),
                   "limit: early.nu ran to 0 (a step at t_half)")
  expect_identical(.hzr_refit_failure_reason(e), .hzr_limit_reason(e))
})

test_that("the mocked layer refuses every refit (the route is reached)", {
  d <- lc_data()
  base <- lc_base(d)
  err <- with_refusal(tryCatch(
    hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
           phases = list(constant = hzr_phase("constant", formula = ~ z)),
           fit = TRUE),
    hzr_limit_refusal = function(e) e))
  expect_s3_class(err, "hzr_limit_refusal")
  expect_s3_class(base, "hazard")
})

test_that("stepwise forward tallies a refusal as a failed candidate", {
  d <- lc_data()
  base <- lc_base(d)
  res <- collect_warnings(
    with_refusal(hzr_stepwise(base, scope = list(constant = ~ z + w),
                              data = d, direction = "forward",
                              criterion = "wald", trace = FALSE)))
  sw <- res$value
  expect_true(any(grepl("candidate refit failed for .*: limit: ", res$warnings)))
  expect_true(any(grepl("Stepwise selection stopped", res$warnings)))
  cr <- sw$criteria
  # Every candidate refit refused, so every candidate is a counted failure.
  expect_setequal(sub("@.*", "", cr$refit_failures), c("z", "w"))
  expect_identical(cr$n_refit_failures, 2L)
  expect_length(cr$refit_failure_reasons, 2L)
  expect_true(all(startsWith(cr$refit_failure_reasons, "limit: ")))
  expect_true(cr$stopped_refit_failed)
  expect_identical(nrow(sw$steps), 0L)
})

test_that("stepwise backward records a refused drop as a failed drop", {
  d <- lc_data()
  full <- lc_base(d, ~ z + v)
  # The drop is accepted by the screen (v is noise); only its refit refuses.
  res <- collect_warnings(
    with_refusal(hzr_stepwise(full, data = d, direction = "backward",
                              criterion = "wald", slstay = 0.9,
                              trace = FALSE)))
  sw <- res$value
  expect_true(any(grepl("post-drop refit failed for .*: limit: ", res$warnings)))
  cr <- sw$criteria
  expect_gte(cr$n_refit_failures, 1L)
  expect_true(all(startsWith(cr$refit_failure_reasons, "limit: ")))
  expect_true(cr$stopped_refit_failed)
  # The model was not changed: the refused drop is not a removal.
  expect_identical(nrow(sw$steps), 0L)
})

test_that("bootstrap tallies every refused replicate under one reason", {
  d <- lc_data()
  base <- lc_base(d)
  boot <- with_refusal(suppressWarnings(hzr_bootstrap(base, n_boot = 3L)))
  expect_identical(boot$n_failed, 3L)
  tally <- boot$failure_reasons
  expect_length(tally, 1L)
  expect_identical(names(tally), "limit: constant.nu ran to 0 (a step at t_half)")
  expect_equal(unname(as.integer(tally)), 3L)
})

# Every refit inside `code` meets a hard trigger on its first (outer) fit and
# a clean one on its single grid start, so the real layer runs its grid and
# returns outcome "restarted". Each hazard() call reads the trigger twice
# (outer, then the one grid run), so the mock alternates; a desync fails the
# n_restarted assertions below loudly.
with_restart <- function(code) {
  state <- new.env()
  state$outer <- TRUE
  testthat::local_mocked_bindings(
    .hzr_limit_trigger = function(object) {
      hit <- state$outer
      state$outer <- !state$outer
      if (hit) limit_trigger_hard(names(object$spec$phases)[[1]])
    },
    .hzr_limit_grid = function(start, ...) list(start),
    .package = "TemporalHazard")
  force(code)
}

test_that("bootstrap counts every restarted replicate, and prints the count", {
  d <- lc_data()
  base <- lc_base(d)
  boot <- with_restart(suppressWarnings(hzr_bootstrap(base, n_boot = 3L)))
  expect_identical(boot$n_failed, 0L)
  expect_identical(boot$n_success, 3L)
  expect_identical(boot$n_restarted, 3L)
  tally <- boot$restart_reasons
  expect_length(tally, 1L)
  expect_identical(names(tally),
                   "limit: constant.nu ran to 0 (a step at t_half)")
  expect_equal(unname(as.integer(tally)), 3L)
  expect_true(all(boot$replicates$restarted))
  expect_match(utils::capture.output(print(boot)),
               "Restarted replicates: 3 of 3 successful", all = FALSE,
               fixed = TRUE)
})

test_that("a clean bootstrap reports zero restarted, and still prints it", {
  d <- lc_data()
  base <- lc_base(d)
  boot <- suppressWarnings(hzr_bootstrap(base, n_boot = 3L))
  expect_identical(boot$n_restarted, 0L)
  expect_length(boot$restart_reasons, 0L)
  expect_false(any(boot$replicates$restarted))
  expect_match(utils::capture.output(print(boot)),
               "Restarted replicates: 0 of 3 successful", all = FALSE,
               fixed = TRUE)
})
