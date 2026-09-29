# #490 refused an AIC entry whose refit ended below the current model's
# log-likelihood: the candidate model contains the current one, so such a
# refit did not converge. Under criterion = "wald", and in the score
# criterion's Wald fallback, the same refit still got a Wald p-value. At a
# strict `slentry` it was rejected as if tested, with every counter at 0 and
# no warning (#538). It is now refused with the reason `loglik_below_base`.

.wald_below_fixture <- function() {
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

.wald_below_screen <- function(fx, ...) {
  msgs <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fx$fit, data = fx$data, direction = "forward",
                 trace = FALSE, control = fx$control, ...),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(sw = sw, msgs = msgs)
}

test_that("a Wald entry refuses a refit that ends below its base (#538)", {
  fx <- .wald_below_fixture()
  cand <- suppressWarnings(.hzr_refit_with_scope(
    fx$fit, action = "add", var = "age", phase = "constant",
    data = fx$data, control = fx$control
  ))
  nm <- .hzr_candidate_coef_name(cand, "age", "constant", current = fx$fit)
  # Known positive: a refit at or above its base is tested as before.
  expect_gte(cand$fit$objective, fx$fit$fit$objective)
  ok <- .hzr_candidate_score("wald", "entry", current = fx$fit,
                             candidate = cand, names = nm)
  expect_true(is.finite(ok$score))
  expect_true(is.finite(ok$p_value))
  expect_true(is.na(ok$reason))

  low <- cand
  low$fit$objective <- fx$fit$fit$objective - 1
  s <- .hzr_candidate_score("wald", "entry", current = fx$fit,
                            candidate = low, names = nm)
  expect_true(is.na(s$score))
  expect_true(is.na(s$p_value))
  expect_true(is.na(s$stat))
  expect_identical(s$reason, "loglik_below_base")

  # Optimizer noise is not a failed refit.
  noise <- cand
  noise$fit$objective <- fx$fit$fit$objective - 1e-12
  s <- .hzr_candidate_score("wald", "entry", current = fx$fit,
                            candidate = noise, names = nm)
  expect_true(is.finite(s$score))
  expect_true(is.na(s$reason))

  # A drop is tested on the current model, which has no base to fall below.
  d <- .hzr_candidate_score("wald", "drop", current = cand, names = nm)
  expect_true(is.finite(d$score))
})

test_that("hzr_stepwise() under Wald counts and warns about it (#538)", {
  skip_on_cran()
  fx <- .wald_below_fixture()
  cand <- suppressWarnings(.hzr_refit_with_scope(
    fx$fit, action = "add", var = "opmos", phase = "constant",
    data = fx$data, control = fx$control
  ))
  # The fixture property: the refit reports convergence, ends well below its
  # base, and still has a Wald p-value that misses slentry = 0.05.
  expect_true(cand$fit$converged)
  expect_lt(cand$fit$objective, fx$fit$fit$objective - 1)

  out <- .wald_below_screen(fx, scope = list(constant = ~ opmos + age),
                            criterion = "wald", slentry = 0.05)
  sw <- out$sw
  expect_equal(nrow(sw$steps), 0L)
  expect_false(sw$criteria$stopped_uncomputable)
  expect_equal(sw$criteria$uncomputable_reasons, c(loglik_below_base = 1L))
  # It was not declined for want of a variance, and is not reported so.
  expect_length(sw$criteria$wald_untested_entries, 0L)
  expect_false(any(grepl("without a Wald test", out$msgs, fixed = TRUE)))
  declined <- grep("declined 1 candidate entry without testing", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(declined, 1L)
  expect_match(declined, "below the current model's", fixed = TRUE)
})

test_that("a Wald screen of only such a candidate stops untested (#538)", {
  skip_on_cran()
  fx <- .wald_below_fixture()
  out <- .wald_below_screen(fx, scope = list(constant = ~ opmos),
                            criterion = "wald", slentry = 0.05)
  expect_true(out$sw$criteria$stopped_uncomputable)
  expect_equal(out$sw$criteria$uncomputable_reasons,
               c(loglik_below_base = 1L))
  stop_msg <- grep("stopped because no remaining candidate", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(stop_msg, 1L)
  expect_match(stop_msg, "below the current model's", fixed = TRUE)
})

# The fallback is entered for a candidate whose score information is
# indefinite. That trigger is forced here; the refit and the Wald test that
# follow are the real ones.
test_that("the score criterion's Wald fallback refuses it too (#538)", {
  skip_on_cran()
  fx <- .wald_below_fixture()
  orig_q <- .hzr_score_q
  local_mocked_bindings(.hzr_score_q = function(...) {
    r <- orig_q(...)
    r$stat <- NA_real_
    r$p_value <- NA_real_
    r$reason <- "information_indefinite"
    r
  })
  # Known positive: a candidate above its base is rescued and tested.
  ok <- .wald_below_screen(fx, scope = list(constant = ~ age),
                           criterion = "score", slentry = 0.05)
  expect_equal(ok$sw$criteria$n_wald_fallbacks, 1L)
  expect_length(ok$sw$criteria$uncomputable_reasons, 0L)

  out <- .wald_below_screen(fx, scope = list(constant = ~ opmos + age),
                            criterion = "score", slentry = 0.05)
  cr <- out$sw$criteria
  expect_equal(cr$n_wald_fallbacks, 1L)
  expect_equal(cr$uncomputable_reasons, c(loglik_below_base = 1L))
  declined <- grep("declined 1 candidate entry without testing", out$msgs,
                   fixed = TRUE, value = TRUE)
  expect_length(declined, 1L)
})

# The latest step decides whether an entry is listed as untested for want of
# a variance. A variable with no variance at one step and refused below its
# base at the next was left on that list, and warned about with the wrong
# cause. The two steps' outcomes are marked here so the sequence is exact.
test_that("a later below-base refusal clears an earlier variance entry", {
  skip_on_cran()
  fx <- .wald_below_fixture()
  orig_fs <- .hzr_stepwise_forward_step
  n_call <- 0L
  local_mocked_bindings(.hzr_stepwise_forward_step = function(...) {
    r <- orig_fs(...)
    n_call <<- n_call + 1L
    age <- r$all_scores$variable == "age"
    if (n_call == 1L) {
      # Step 1: age's Wald test had no variance; mal enters.
      r$all_scores$score[age] <- NA_real_
    } else {
      # Step 2: age's refit ended below its base, and nothing enters.
      r$all_scores$score[age] <- NA_real_
      r$all_scores$reason[age] <- "loglik_below_base"
      r$accepted <- FALSE
      r$fit <- list(...)$current
      r$n_uncomputable <- 1L
      r$uncomputable_reasons <- c(loglik_below_base = 1L)
    }
    r
  })
  out <- .wald_below_screen(fx, scope = list(early = ~ age + mal),
                            criterion = "wald")
  # The sequence ran as intended: mal entered at step 1, and step 2 ran.
  expect_equal(n_call, 2L)
  expect_equal(out$sw$steps$variable, "mal")
  expect_length(out$sw$criteria$wald_untested_entries, 0L)
  expect_false(any(grepl("without a Wald test", out$msgs, fixed = TRUE)))
  expect_true(any(grepl("declined 1 candidate entry without testing",
                        out$msgs, fixed = TRUE)))
})
