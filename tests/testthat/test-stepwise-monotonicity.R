# A forward step must not lower the log-likelihood --------------------------
#
# The model a forward step enters CONTAINS the model it started from, so at
# the optimum the objective cannot fall. When it does, the refit did not
# converge, and every later step is scored against a model that is not at its
# own optimum.
#
# `hzr_stepwise()` wrote the objective into `$steps` at every step and
# compared it at none: `grep -n objective R/stepwise.R` returned three write
# sites and no read. A production 2-phase screen (issue #134) accepted four
# collinear forms of one covariate, three of which lowered the
# log-likelihood, and finished reporting `converged`, ten entries and
# `p = 0.000` throughout. The final 19-coefficient model fitted 57 units worse
# than the nested 16-coefficient model from three steps earlier -- which is
# arithmetically impossible at a maximum, and the only thing that noticed was
# a monotonicity check written by hand downstream.

.mono_fit <- function() {
  data("avc", package = "TemporalHazard", envir = environment())
  set.seed(1L)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L))
  list(fit = fit, data = avc)
}

test_that("a healthy forward run records delta_logLik and reports none worse", {
  o <- .mono_fit()
  r <- suppressWarnings(hzr_stepwise(
    o$fit, scope = list(early = ~ age, constant = ~ age), data = o$data,
    direction = "forward", criterion = "wald", slentry = 0.5, trace = FALSE,
    max_steps = 2L, control = list(n_starts = 1L, maxit = 500L)))

  expect_true("delta_logLik" %in% names(r$steps))
  entries <- r$steps[r$steps$action == "enter", ]
  expect_gt(nrow(entries), 0L)
  # Every entry gained: this is the invariant, not merely a recorded number.
  expect_true(all(entries$delta_logLik >= 0))
  expect_identical(r$criteria$n_nonmonotone_entries, 0L)
})

test_that("a forward step that lowers the log-likelihood warns and is counted", {
  o <- .mono_fit()
  # Force the reported failure: take the real accepted step and degrade its
  # objective. Reproducing it from data needs a near-singular design and a
  # refit that genuinely fails, which is neither quick nor stable across
  # platforms -- and the branch under test is the comparison, not the cause.
  # SET relative to the model being entered FROM, rather than subtracting from
  # the entered value. The genuine step here gains ~14.6 units, so a fixed
  # subtraction has to beat a number the fixture does not control -- an
  # earlier draft subtracted 5, left the step 9.6 units ahead, and asserted
  # nothing.
  base_obj <- o$fit$fit$objective
  orig <- .hzr_stepwise_forward_step
  local_mocked_bindings(
    .hzr_stepwise_forward_step = function(...) {
      out <- orig(...)
      if (isTRUE(out$accepted)) {
        out$fit$fit$objective <- base_obj - 5
      }
      out
    }
  )

  w <- capture_warnings(
    r <- hzr_stepwise(
      o$fit, scope = list(early = ~ age, constant = ~ age), data = o$data,
      direction = "forward", criterion = "wald", slentry = 0.5, trace = FALSE,
      max_steps = 1L, control = list(n_starts = 1L, maxit = 500L)))

  expect_match(w, "log-likelihood FELL", all = FALSE)
  # The message must say WHY it is impossible, not merely that it happened.
  expect_match(w, "contains the current one", all = FALSE)
  expect_gt(r$criteria$n_nonmonotone_entries, 0L)

  entries <- r$steps[r$steps$action == "enter", ]
  expect_true(any(entries$delta_logLik < 0))
})

test_that("the tolerance does not fire on optimizer noise", {
  # A drop of 1e-12 is not evidence of a failed refit, and a check that
  # cried wolf on rounding would be turned off rather than read.
  o <- .mono_fit()
  base_obj <- o$fit$fit$objective
  orig <- .hzr_stepwise_forward_step
  local_mocked_bindings(
    .hzr_stepwise_forward_step = function(...) {
      out <- orig(...)
      if (isTRUE(out$accepted)) {
        # Just below the starting objective, not below the entered one: the
        # step must land inside the tolerance band to test it at all.
        out$fit$fit$objective <- base_obj - 1e-12
      }
      out
    }
  )

  r <- suppressWarnings(hzr_stepwise(
    o$fit, scope = list(early = ~ age, constant = ~ age), data = o$data,
    direction = "forward", criterion = "wald", slentry = 0.5, trace = FALSE,
    max_steps = 1L, control = list(n_starts = 1L, maxit = 500L)))
  expect_identical(r$criteria$n_nonmonotone_entries, 0L)
})
