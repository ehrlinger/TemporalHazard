# Under criterion = "aic" an entry adds one term to the current model, so the
# candidate model contains the current one and its log-likelihood cannot be
# lower at the optimum. A refit that ends below it did not converge, however
# `converged` reads. Its dAIC was scored all the same, came out positive, and
# the candidate was rejected as if it had been tested, with every counter at 0
# and no warning (#490). Such a candidate is now refused, recorded as
# `loglik_below_base`, and warned about.

.aic_below_fixture <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  ctl <- list(n_starts = 1L)
  set.seed(1)
  fit <- suppressWarnings(hazard(
    Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = ctl
  ))
  list(fit = fit, data = avc, control = ctl)
}

test_that("an AIC entry refuses a refit that ends below its base (#490)", {
  fx <- .aic_below_fixture()
  cand <- suppressWarnings(.hzr_refit_with_scope(
    fx$fit, action = "add", var = "age", phase = "constant",
    data = fx$data, control = fx$control
  ))
  nm <- .hzr_candidate_coef_name(cand, "age", "constant",
                                 current = fx$fit)
  # Known positive: a refit at or above its base is scored as before.
  expect_gte(cand$fit$objective, fx$fit$fit$objective)
  ok <- .hzr_candidate_score("aic", "entry", current = fx$fit,
                             candidate = cand, names = nm)
  expect_true(is.finite(ok$score))
  expect_true(is.na(ok$reason))

  # The same refit, ended one log-likelihood unit below its base.
  low <- cand
  low$fit$objective <- fx$fit$fit$objective - 1
  s <- .hzr_candidate_score("aic", "entry", current = fx$fit,
                            candidate = low, names = nm)
  expect_true(is.na(s$score))
  expect_true(is.na(s$delta_aic))
  expect_identical(s$reason, "loglik_below_base")

  # Optimizer noise is not a failed refit.
  noise <- cand
  noise$fit$objective <- fx$fit$fit$objective - 1e-12
  s <- .hzr_candidate_score("aic", "entry", current = fx$fit,
                            candidate = noise, names = nm)
  expect_true(is.finite(s$score))
  expect_true(is.na(s$reason))
})

.aic_below_screen <- function(fx, ...) {
  msgs <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fx$fit, data = fx$data, direction = "forward",
                 criterion = "aic", trace = FALSE, control = fx$control, ...),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(sw = sw, msgs = msgs)
}

test_that("hzr_stepwise() under AIC counts and warns about it (#490)", {
  skip_on_cran()
  fx <- .aic_below_fixture()
  # The fixture has the property the test needs: this refit reports
  # convergence and ends well below the model it contains.
  cand <- suppressWarnings(.hzr_refit_with_scope(
    fx$fit, action = "add", var = "opmos", phase = "constant",
    data = fx$data, control = fx$control
  ))
  expect_true(cand$fit$converged)
  expect_lt(cand$fit$objective, fx$fit$fit$objective - 1)

  out <- .aic_below_screen(fx, scope = list(constant = ~ opmos + age))
  sw <- out$sw
  # age is scored and misses (dAIC > 0), so nothing enters: the screen
  # completed, and only opmos went untested.
  expect_equal(nrow(sw$steps), 0L)
  expect_false(sw$criteria$stopped_uncomputable)
  expect_equal(sw$criteria$n_uncomputable_scores, 1L)
  expect_equal(sw$criteria$uncomputable_reasons,
               c(loglik_below_base = 1L))
  declined <- grep("declined 1 candidate entry without testing", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(declined, 1L)
  expect_match(declined, "below the current model's", fixed = TRUE)
})

test_that("an AIC screen of only such a candidate stops untested (#490)", {
  skip_on_cran()
  fx <- .aic_below_fixture()
  out <- .aic_below_screen(fx, scope = list(constant = ~ opmos))
  sw <- out$sw
  expect_true(sw$criteria$stopped_uncomputable)
  expect_equal(sw$criteria$uncomputable_reasons,
               c(loglik_below_base = 1L))
  stop_msg <- grep("stopped because no remaining candidate", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(stop_msg, 1L)
  expect_match(stop_msg, "below the current model's", fixed = TRUE)
})
