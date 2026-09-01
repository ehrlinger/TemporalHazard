# Conservation of Events: recording whether it was actually applied (#146).
#
# R/likelihood-multiphase.R disables CoE whenever any status falls outside
# {0, 1}. That is deliberate and correct -- the CoE identity counts only exact
# events -- but the fact that it happened was not recorded anywhere. A caller
# who passed control = list(conserve = TRUE) got a fit carrying conserve = TRUE
# over a computation that did not run: the package's documented signature
# defect shape.
#
# It is not an edge case. ICENSOR appears on 42-74 blocks per production study,
# and ICENSOR guarantees status leaves {0, 1}, so this fires constantly.

# Explicit starts rather than hzr_phase("cdf")'s defaults: the default start
# does not survive left-censored rows on this fixture ("t_half must be a
# positive scalar" out of the optimizer), which is unrelated to CoE and would
# make the -1 case below fail for the wrong reason. Explicit starts with
# n_starts = 1 is also the house preference -- it removes the optimizer path
# from what the test depends on.
coe_fit <- function(status_pattern, conserve = TRUE,
                    phases = list(
                      early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                      const = hzr_phase("constant")
                    )) {
  set.seed(42)
  n <- 60
  time <- stats::rexp(n, 0.2) + 0.05
  status <- rep(status_pattern, length.out = n)

  # Interval rows (status 2) need bounds; left-censored rows (-1) take their
  # upper bound from `time` and must NOT be given interval-style bounds, or
  # the fit walks onto an invalid shape. Supply the bounds only where the
  # fixture actually has interval rows.
  args <- list(time = time, status = status)
  if (any(status == 2)) {
    args$time_lower <- time
    args$time_upper <- ifelse(status == 2, time * 1.5, time)
  }

  # fit = TRUE is required: hazard() defaults to FALSE, and an unfitted object
  # never reaches the optimizer, so every assertion below would be about a
  # model that was never run. (The test proposed in #146 omitted it, and would
  # have passed against a model that never ran.)
  suppressWarnings(do.call(hazard, c(args, list(
    phases = phases, dist = "multiphase", fit = TRUE,
    control = list(conserve = conserve, n_starts = 1)
  ))))
}

test_that("CoE is disabled, and says so, when interval rows are present", {
  skip_on_cran()
  fit <- coe_fit(c(1, 0, 2))
  # Guard the guard: assert the fit actually ran, or the assertions below are
  # about an unfitted object where every field is absent.
  expect_false(is.null(fit$fit$theta))
  expect_true(fit$spec$control$conserve)          # what was asked for
  expect_false(fit$spec$control$conserve_applied) # what happened
  expect_identical(fit$spec$control$conserve_disabled_reason,
                   "unsupported_censoring")
})

test_that("CoE is applied, and says so, on exact + right-censored data", {
  skip_on_cran()
  # The column has to vary or it reports nothing. Same fixture, same request,
  # only the censoring differs.
  fit <- coe_fit(c(1, 0))
  expect_false(is.null(fit$fit$theta))
  expect_true(fit$spec$control$conserve_applied)
  expect_true(is.na(fit$spec$control$conserve_disabled_reason))
})

test_that("left-censored rows disable CoE too, not only interval ones", {
  skip_on_cran()
  # status -1 is left-censored in this package's coding. The identity counts
  # exact events, so anything outside {0, 1} disables it.
  fit <- coe_fit(c(1, 0, -1))
  expect_false(is.null(fit$fit$theta))
  expect_false(fit$spec$control$conserve_applied)
  expect_identical(fit$spec$control$conserve_disabled_reason,
                   "unsupported_censoring")
})

test_that("a caller who declined CoE is distinguished from one who was refused", {
  skip_on_cran()
  # A bare FALSE reads as "you turned it off" to a user who did the opposite.
  # The reason is what separates the two.
  fit <- coe_fit(c(1, 0), conserve = FALSE)
  expect_false(is.null(fit$fit$theta))
  expect_false(fit$spec$control$conserve_applied)
  expect_identical(fit$spec$control$conserve_disabled_reason, "not_requested")
})

test_that("a single-phase model reports its own reason", {
  skip_on_cran()
  # CoE distributes a conserved event total across phases; with one phase there
  # is nothing to distribute.
  fit <- coe_fit(c(1, 0), phases = list(const = hzr_phase("constant")))
  expect_false(is.null(fit$fit$theta))
  expect_false(fit$spec$control$conserve_applied)
  expect_identical(fit$spec$control$conserve_disabled_reason, "single_phase")
})
