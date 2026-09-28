# Under criterion = "aic" an entry compares the candidate's AIC with the base
# model's. A multiphase refit drops every row where the candidate is missing,
# so a candidate with NAs was scored on fewer rows than the base, its
# log-likelihood was smaller in magnitude, and a noise variable entered with a
# large negative dAIC (#488). Such a candidate is now refused, recorded as
# `rows_differ`, and not entered.

.aic_rows_fixture <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  set.seed(7)
  avc$z <- stats::rnorm(nrow(avc))
  avc$z[sample(nrow(avc), 60)] <- NA
  avc$w <- stats::rnorm(nrow(avc))
  ctl <- list(n_starts = 2L, maxit = 500L)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = ctl
  )
  list(fit = fit, data = avc, control = ctl)
}

test_that("an AIC entry refuses a candidate fitted on fewer rows (#488)", {
  skip_on_cran()
  fx <- .aic_rows_fixture()
  # The fixture has the property the test needs: z is missing on 60 rows, and
  # the base fit used every row.
  expect_equal(sum(is.na(fx$data$z)), 60L)
  expect_equal(.hzr_fit_rows_used(fx$fit), nrow(fx$data))

  fs <- suppressWarnings(.hzr_stepwise_forward_step(
    fx$fit, scope = list(early = ~ z + w), data = fx$data,
    criterion = "aic", control = fx$control
  ))
  sc <- fs$all_scores
  # Known positive: the complete-data candidate is scored as before.
  expect_true(is.finite(sc$delta_aic[sc$variable == "w"]))
  # z is not scored, so it cannot win.
  expect_true(is.na(sc$score[sc$variable == "z"]))
  expect_true(is.na(sc$delta_aic[sc$variable == "z"]))
  expect_true(fs$accepted)
  expect_equal(fs$variable, "w")
  expect_equal(fs$n_uncomputable, 1L)
  expect_equal(fs$uncomputable_reasons, c(rows_differ = 1L))
})

test_that("an AIC entry screen reports every candidate refused for rows", {
  skip_on_cran()
  fx <- .aic_rows_fixture()
  fs <- suppressWarnings(.hzr_stepwise_forward_step(
    fx$fit, scope = list(early = ~ z), data = fx$data,
    criterion = "aic", control = fx$control
  ))
  # Nothing could be scored, which must read differently from nothing
  # scoring well enough.
  expect_false(fs$accepted)
  expect_identical(fs$stop_reason, "scores_uncomputable")
  expect_equal(fs$uncomputable_reasons, c(rows_differ = 1L))
})

.aic_rows_screen <- function(fx, ...) {
  msgs <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fx$fit, scope = list(early = ~ z + w), data = fx$data,
                 direction = "forward", criterion = "aic", trace = FALSE,
                 control = fx$control, ...),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(sw = sw, msgs = msgs)
}

test_that("hzr_stepwise() under AIC does not enter z and says why (#488)", {
  skip_on_cran()
  fx <- .aic_rows_fixture()
  out <- .aic_rows_screen(fx)
  sw <- out$sw
  expect_equal(sw$steps$variable, "w")
  # z is refused at both iterations, and the second has nothing else to
  # score, so the run stops untested rather than reading as a clean finish.
  expect_equal(sw$criteria$uncomputable_reasons, c(rows_differ = 2L))
  expect_true(sw$criteria$stopped_uncomputable)
  stop_msg <- grep("stopped because no remaining candidate", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(stop_msg, 1L)
  expect_match(stop_msg, "different rows", fixed = TRUE)
})

# Not skipped on CRAN: the fastest block holding both the known positive (w
# enters) and the refusal, and the only one that kills all four mutants.
test_that("a completed AIC screen warns about a refused entry (#488)", {
  fx <- .aic_rows_fixture()
  out <- .aic_rows_screen(fx, max_steps = 1L)
  expect_equal(out$sw$steps$variable, "w")
  expect_false(out$sw$criteria$stopped_uncomputable)
  expect_equal(out$sw$criteria$uncomputable_reasons, c(rows_differ = 1L))
  declined <- grep("declined 1 candidate entry without testing", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(declined, 1L)
  expect_match(declined, "different rows", fixed = TRUE)
})
