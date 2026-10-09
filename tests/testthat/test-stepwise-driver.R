# `.fit_driver_base()` lives in helper-stepwise.R (auto-sourced by
# testthat) so it can be shared with the integration tests.

# Return class and shape ----------------------------------------------------

test_that("hzr_stepwise returns an `hzr_stepwise`-classed `hazard` object", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "forward", trace = FALSE
  )
  expect_s3_class(res, "hzr_stepwise")
  expect_s3_class(res, "hazard")
  # Still a usable hazard fit
  expect_s3_class(summary(res), "summary.hzr_stepwise")
  expect_true(is.data.frame(res$steps))
  expect_true(is.character(res$trace_msg))
})


# Forward-only happy path ---------------------------------------------------

test_that("forward selection picks the strong signal and ignores noise", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "forward", criterion = "wald",
    slentry = 0.05, trace = FALSE
  )
  # x1 should be the only (or first) variable to enter
  entered <- res$steps$variable[res$steps$action == "enter"]
  expect_true("x1" %in% entered)
  # Scores for entered vars all clear slentry
  expect_true(all(res$steps$p_value[res$steps$action == "enter"] < 0.05,
                  na.rm = TRUE))
})


# Backward-only drops noise -------------------------------------------------

test_that("backward selection drops weak covariates from an overfit model", {
  obj <- .fit_driver_base()
  overfit <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "forward", criterion = "wald",
    slentry = 0.99, trace = FALSE   # admit everyone
  )
  # Now strip back
  trimmed <- hzr_stepwise(
    overfit, data = obj$data,
    direction = "backward", criterion = "wald",
    slstay = 0.20, trace = FALSE
  )
  dropped <- trimmed$steps$variable[trimmed$steps$action == "drop"]
  expect_true(any(c("x2", "x3", "x4") %in% dropped))
  expect_false("x1" %in% dropped)
})


# Two-way convergence -------------------------------------------------------

test_that("direction = both converges to a stable model", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "both", criterion = "wald",
    slentry = 0.30, slstay = 0.20, trace = FALSE
  )
  # Final model contains x1 and no strictly-dropped var
  final_vars <- colnames(res$data$x)
  expect_true("x1" %in% final_vars)
  # The elapsed field is populated
  expect_s3_class(res$elapsed, "difftime")
})


# AIC criterion -------------------------------------------------------------

test_that("criterion = aic uses ΔAIC and records it in the trace", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "forward", criterion = "aic", trace = FALSE
  )
  expect_equal(unique(res$steps$criterion), "aic")
  enters <- res$steps[res$steps$action == "enter", ]
  expect_true(all(enters$delta_aic < 0, na.rm = TRUE))
  # ΔAIC formatting appears in the trace for each entry
  expect_true(any(grepl("AIC", res$trace_msg)))
})


# force_out ----------------------------------------------------------------

test_that("force_out keeps a variable from ever entering", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "forward", criterion = "wald",
    slentry = 0.99,                    # admit everyone else
    force_out = "x1",
    trace = FALSE
  )
  expect_false("x1" %in% res$steps$variable[res$steps$action == "enter"])
})


# force_in -----------------------------------------------------------------

test_that("force_in variables stay in even with high p-value", {
  # Start from an overfit model including x2 (pure noise) and force it in
  obj <- .fit_driver_base(signal_beta = 0)
  full <- hazard(
    Surv(time, status) ~ x2,
    data = obj$data,
    theta = c(0.5, 1.0, 0),
    dist = "weibull", fit = TRUE
  )
  res <- hzr_stepwise(
    full, scope = NULL, data = obj$data,
    direction = "backward", criterion = "wald",
    slstay = 0.01,          # very strict — would normally drop x2
    force_in = "x2", trace = FALSE
  )
  # x2 never drops
  expect_false("x2" %in% res$steps$variable[res$steps$action == "drop"])
  # x2 still in final model
  expect_true("x2" %in% colnames(res$data$x))
})


# MOVE oscillation guard ----------------------------------------------------

test_that("max_move = 0 freezes a variable after its first entry", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    direction = "forward", criterion = "wald",
    slentry = 0.99, max_move = 0L, trace = FALSE
  )
  # The first accepted entry triggers a freeze immediately.
  expect_true("frozen" %in% res$steps$action)
  # The frozen variable is in the scope record
  expect_true(length(res$scope$frozen) >= 1L)
})


# max_steps cap -------------------------------------------------------------

test_that("hitting max_steps warns and stops", {
  obj <- .fit_driver_base()
  expect_warning(
    res <- hzr_stepwise(
      obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
      direction = "forward", criterion = "wald",
      slentry = 0.99, max_steps = 1L, trace = FALSE
    ),
    "max_steps"
  )
  expect_true(res$criteria$hit_max_steps)
  expect_true(nrow(res$steps) >= 1L)
})


# Trace output --------------------------------------------------------------

test_that("trace captures header, per-step lines, and a final summary", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1, data = obj$data,
    direction = "forward", criterion = "wald", trace = FALSE
  )
  # Header mentions the criterion
  expect_match(res$trace_msg[1L], "criterion = wald")
  # At least one step line
  expect_true(any(grepl("Step 1", res$trace_msg)))
  # Final line mentions AIC
  expect_true(any(grepl("AIC", res$trace_msg)))
})


# Print + summary methods ---------------------------------------------------

test_that("print.hzr_stepwise emits the captured trace", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1, data = obj$data,
    direction = "forward", criterion = "wald", trace = FALSE
  )
  out <- capture.output(print(res))
  expect_true(any(grepl("Stepwise selection", out)))
  expect_true(any(grepl("Step 1", out)))
})

test_that("summary.hzr_stepwise dispatches to summary.hazard on the final fit", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1, data = obj$data,
    direction = "forward", criterion = "wald", trace = FALSE
  )
  s <- summary(res)
  expect_s3_class(s, "summary.hzr_stepwise")
  expect_s3_class(s, "summary.hazard")
  # Final fit's coef table is accessible
  expect_true(is.data.frame(s$coefficients))
})


# Input validation ---------------------------------------------------------

test_that("hzr_stepwise rejects non-hazard fit", {
  expect_error(
    hzr_stepwise(list(), data = data.frame()),
    "must be a `hazard` object"
  )
})

test_that("hzr_stepwise requires `data`", {
  obj <- .fit_driver_base()
  expect_error(hzr_stepwise(obj$fit), "`data` must be a data frame")
  expect_error(hzr_stepwise(obj$fit, data = "not a df"),
               "`data` must be a data frame")
})


# Polish helpers -----------------------------------------------------------

test_that("as.data.frame returns the $steps trace", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2, data = obj$data,
    direction = "forward", criterion = "wald", trace = FALSE
  )
  df <- as.data.frame(res)
  expect_s3_class(df, "data.frame")
  expect_identical(df, res$steps)
  expect_true(all(c("step_num", "action", "variable", "phase",
                    "score", "p_value", "delta_aic") %in% names(df)))
})

test_that("stepwise_trace returns the captured console lines", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1, data = obj$data,
    direction = "forward", criterion = "wald", trace = FALSE
  )
  tr <- stepwise_trace(res)
  expect_identical(tr, res$trace_msg)
  expect_true(any(grepl("Stepwise selection", tr)))
})

test_that("stepwise_trace rejects non-hzr_stepwise input", {
  expect_error(stepwise_trace(list()),
               "must be an `hzr_stepwise` object")
})


# Score criterion -----------------------------------------------------------

test_that("hzr_stepwise defaults to criterion = 'score'", {
  expect_identical(eval(formals(hzr_stepwise)$criterion)[1], "score")
})

test_that("hzr_stepwise(criterion = 'score') selects without refitting candidates", {
  skip_if_not_installed("numDeriv")
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    criterion = "score", direction = "forward", trace = FALSE
  )
  expect_s3_class(res, "hzr_stepwise")
  expect_true(nrow(res$steps) > 0)
  expect_identical(res$steps$variable[1], "x1")
  expect_true(all(res$steps$df == 1L))
  # The score path reports no dAIC: there is no candidate refit to take it from.
  expect_true(all(is.na(res$steps$delta_aic)))
})

test_that("criterion = 'wald' still reproduces the previous behaviour", {
  obj <- .fit_driver_base()
  res <- hzr_stepwise(
    obj$fit, scope = ~ x1 + x2 + x3 + x4, data = obj$data,
    criterion = "wald", direction = "forward", trace = FALSE
  )
  expect_s3_class(res, "hzr_stepwise")
  expect_identical(res$steps$variable[1], "x1")
})

test_that("a scope naming a nonexistent column warns under both criteria", {
  # NEWS claims a mistyped scope column "still surfaces once, up front" --
  # and hzr_bootstrap()'s pre-loop scope validation (see
  # test-diagnostics.R's "hzr_bootstrap scope with a nonexistent column
  # warns but does not raise") relies on hzr_stepwise() warning here rather
  # than silently dropping the candidate. Confirm the score path matches
  # the pre-existing wald behaviour instead of staying silent.
  skip_if_not_installed("numDeriv")
  obj <- .fit_driver_base()

  expect_warning(
    res_wald <- hzr_stepwise(
      obj$fit, scope = c("x1", "not_a_real_column"), data = obj$data,
      criterion = "wald", direction = "forward", trace = FALSE
    ),
    "not_a_real_column"
  )
  expect_false("not_a_real_column" %in% res_wald$steps$variable)

  expect_warning(
    res_score <- hzr_stepwise(
      obj$fit, scope = c("x1", "not_a_real_column"), data = obj$data,
      criterion = "score", direction = "forward", trace = FALSE
    ),
    "not_a_real_column"
  )
  expect_false("not_a_real_column" %in% res_score$steps$variable)
})

test_that("criterion = 'score' rejects a non-converged base up front", {
  data(avc)
  avc <- na.omit(avc)
  # One iteration from a real start: the optimizer runs and stops short, so
  # the base has coefficients but did not converge. (This used to omit theta,
  # which never ran the optimizer at all -- converged NA, not FALSE -- and
  # hazard() now refuses that call.)
  base_bad <- suppressWarnings(
    hazard(survival::Surv(int_dead, dead) ~ age, data = avc,
           dist = "weibull", fit = TRUE, theta = c(0.5, 1, 0),
           control = list(maxit = 1))
  )
  expect_identical(base_bad$fit$converged, FALSE)
  expect_length(base_bad$fit$theta, 3L)

  expect_error(
    hzr_stepwise(base_bad, scope = ~ age + mal, data = avc,
                 criterion = "score", direction = "forward", trace = FALSE),
    "did not converge"
  )

  # wald refits each candidate from the call and tolerates the bad base.
  # Assert the screen tested something: a screen whose refits all failed
  # also returns without error.
  res <- hzr_stepwise(base_bad, scope = ~ age + mal, data = avc,
                      criterion = "wald", direction = "forward",
                      trace = FALSE)
  expect_length(res$criteria$refit_failures, 0L)
  expect_identical(res$steps$variable[res$steps$action == "enter"], "mal")
})

# delta AIC is the candidate's AIC minus the base's, so a base that stopped
# short of its maximum hands every candidate the shortfall: pure noise
# entered at delta AIC -16.5 from this base (-249 from maxit = 1).
.sw_noise_data <- function() {
  set.seed(5)
  n <- 300
  noise <- stats::rnorm(n)
  noise2 <- stats::rnorm(n)
  t <- stats::rweibull(n, shape = 1.5, scale = 2)
  data.frame(time = pmin(t, 3), dead = as.integer(t <= 3),
             noise = noise, noise2 = noise2)
}

test_that("criterion = 'aic' rejects a non-converged base up front", {
  d <- .sw_noise_data()
  base_bad <- suppressWarnings(
    hazard(survival::Surv(time, dead) ~ 1, data = d, dist = "weibull",
           theta = c(5, 0.5), control = list(maxit = 5), fit = TRUE)
  )
  # Premise: the optimizer ran and stopped short.
  expect_identical(base_bad$fit$converged, FALSE)
  expect_length(base_bad$fit$theta, 2L)
  for (dir in c("forward", "both")) {
    expect_error(
      hzr_stepwise(base_bad, scope = ~ noise + noise2, data = d,
                   criterion = "aic", direction = dir, trace = FALSE),
      "criterion = 'aic' requires a converged base model"
    )
  }
  # A backward screen takes no scope; it drops from the base's own terms.
  full_bad <- suppressWarnings(
    hazard(survival::Surv(time, dead) ~ noise + noise2, data = d,
           dist = "weibull", theta = c(5, 0.5, 0, 0),
           control = list(maxit = 5), fit = TRUE)
  )
  expect_identical(full_bad$fit$converged, FALSE)
  expect_error(
    hzr_stepwise(full_bad, data = d, criterion = "aic",
                 direction = "backward", trace = FALSE),
    "criterion = 'aic' requires a converged base model"
  )

  # Control: the same base run to convergence screens, and enters no noise.
  base_ok <- suppressWarnings(
    hazard(survival::Surv(time, dead) ~ 1, data = d, dist = "weibull",
           theta = c(5, 0.5), fit = TRUE)
  )
  expect_identical(base_ok$fit$converged, TRUE)
  res <- suppressWarnings(
    hzr_stepwise(base_ok, scope = ~ noise + noise2, data = d,
                 criterion = "aic", direction = "forward", trace = FALSE)
  )
  expect_length(res$criteria$refit_failures, 0L)
  expect_identical(sum(res$steps$action == "enter"), 0L)
})

test_that("criterion = 'wald' rejects a non-converged base it would drop from", {
  # A Wald entry is tested at the candidate's own converged refit, so a
  # forward screen tolerates the bad base (see above). A removal is tested on
  # the current model's own estimates and variance, which for the first step
  # are the base's, at a point the optimizer had not finished with.
  data(avc)
  avc <- na.omit(avc)
  base_bad <- suppressWarnings(
    hazard(survival::Surv(int_dead, dead) ~ age + mal, data = avc,
           dist = "weibull", fit = TRUE, theta = c(0.5, 1, 0, 0),
           control = list(maxit = 1))
  )
  expect_identical(base_bad$fit$converged, FALSE)
  for (dir in c("backward", "both")) {
    expect_error(
      hzr_stepwise(base_bad, data = avc, criterion = "wald",
                   direction = dir, trace = FALSE),
      "direction = '.*' requires a converged base model"
    )
  }
})

test_that("a two-level factor candidate is tested by wald and refused by score", {
  skip_if_not_installed("numDeriv")
  obj <- .fit_driver_base()
  obj$data$fac <- factor(ifelse(obj$data$x2 > 0, "hi", "lo"))

  # Wald refits the candidate and tests its one column, `faclo`. It used to
  # look the coefficient up as `fac`, find nothing, and stop -- while the
  # score refusal below names `criterion = "wald"` as the way through.
  direct <- hazard(Surv(time, status) ~ fac, data = obj$data,
                   theta = c(0.5, 1, 0), dist = "weibull", fit = TRUE)
  expect_identical(colnames(direct$data$x), "faclo")
  z <- unname(stats::coef(direct)[3] / sqrt(stats::vcov(direct)[3, 3]))
  step <- .hzr_stepwise_forward_step(obj$fit, scope = ~ fac, data = obj$data,
                                     criterion = "wald", slentry = 0.05)
  # The refit is warm-started and `direct` is not, so they stop at slightly
  # different points on this near-null effect (about 2.5e-4 apart in z).
  expect_equal(step$all_scores$stat[step$all_scores$variable == "fac"], z,
               tolerance = 1e-3)

  # Score must not silently return NA and drop the candidate on the floor.
  expect_error(
    hzr_stepwise(obj$fit, scope = ~ fac, data = obj$data,
                 criterion = "score", direction = "forward", trace = FALSE),
    "not numeric"
  )
})
