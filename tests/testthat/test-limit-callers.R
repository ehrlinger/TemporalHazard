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
    .hzr_limit_trigger = function(object, unusable = TRUE) {
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
  # Only the accepted drop is refitted, so exactly one refit refused.
  expect_identical(cr$n_refit_failures, 1L)
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
    .hzr_limit_trigger = function(object, unusable = TRUE) {
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

# #736: a refusal from the start that lost a stepwise refit is kept and said.
test_that("a refit's losing start that refused is kept and warned", {
  d <- lc_data()
  base <- lc_base(d)
  refusal <- .hzr_limit_refusal(limit_trigger_hard("constant"), 50, 9L,
                                start_objective = 40)
  warm <- list(value = refusal, warnings = list())
  default <- list(value = base, warnings = list())
  w <- NULL
  fit <- withCallingHandlers(
    .hzr_refit_best_start(warm, default),
    hzr_limit_refit_refused = function(c) {
      w <<- c
      invokeRestart("muffleWarning")
    })
  expect_s3_class(w, "hzr_limit_refit_refused")
  expect_match(conditionMessage(w), paste0(
    "The refit from the warm start refused (constant.nu ran to 0 ",
    "(a step at t_half); objective 40); the default start's fit is kept"),
    fixed = TRUE)
  expect_identical(fit$fit$refit_start, "default")
  expect_identical(fit$fit$refit_refusal, refusal)
  expect_identical(unname(is.na(fit$fit$refit_objectives)), c(TRUE, FALSE))
  # A loser that failed for another reason is not a limit refusal.
  other <- list(value = simpleError("boom"), warnings = list())
  fit2 <- expect_no_warning(.hzr_refit_best_start(other, default))
  expect_null(fit2$fit$refit_refusal)
})

# Only an intercept-only constant phase triggers (theta of length 1), and each
# such hazard() call reads the trigger twice, outer then its one grid start,
# so the base refit of every replicate restarts while every candidate and
# final fit, which carry z, are clean.
with_base_restart <- function(code) {
  state <- new.env()
  state$outer <- TRUE
  testthat::local_mocked_bindings(
    .hzr_limit_trigger = function(object, unusable = TRUE) {
      if (length(object$fit$theta) != 1L) return(NULL)
      hit <- state$outer
      state$outer <- !state$outer
      if (hit) limit_trigger_hard(names(object$spec$phases)[[1]])
    },
    .hzr_limit_grid = function(start, ...) list(start),
    .package = "TemporalHazard")
  force(code)
}

test_that("select-mode bootstrap counts a replicate whose base refit restarted", {
  d <- lc_data()
  base <- lc_base(d)
  boot <- with_base_restart(suppressWarnings(hzr_bootstrap(
    base, n_boot = 3L, scope = list(constant = ~ z),
    criterion = "wald", direction = "forward")))
  expect_identical(boot$n_success, 3L)
  # z entered every replicate, so each final fit is not the restarted base.
  expect_identical(
    sum(grepl("z", boot$replicates$parameter)), 3L)
  expect_identical(boot$n_restarted, 3L)
  expect_true(all(boot$replicates$restarted))
  expect_identical(names(boot$restart_reasons),
                   "limit: constant.nu ran to 0 (a step at t_half)")
})

test_that("a bootstrap with no successful replicate prints 0 of 0 restarted", {
  d <- lc_data()
  base <- lc_base(d)
  boot <- with_refusal(suppressWarnings(hzr_bootstrap(base, n_boot = 2L)))
  expect_identical(boot$n_success, 0L)
  out <- utils::capture.output(print(boot))
  expect_true("Restarted replicates: 0 of 0 successful" %in% out)
})

test_that("control$limit_restart is validated on an unfitted call too", {
  d <- lc_data(50L)
  expect_error(
    hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
           phases = list(constant = hzr_phase("constant")), fit = FALSE,
           control = list(limit_restart = "yes")),
    "control$limit_restart must be TRUE or FALSE.", fixed = TRUE)
  expect_s3_class(
    hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
           phases = list(constant = hzr_phase("constant")), fit = FALSE,
           control = list(limit_restart = FALSE)), "hazard")
})
