# tests/testthat/test-stepwise-refit-failure-reason.R
#
# A stepwise candidate whose refit errors used to be reported as "refit failed
# for <var>" and nothing more: the condition message was caught and dropped.
# The concrete case is the duplicate-design-column refusal. A numeric `gb`
# added to a model holding factor `g` (level `b`) makes hazard() stop with a
# message naming the collision, and stepwise showed only the variable name.
#
# The score path never refits per candidate, so it built the colliding design
# itself: the single-distribution builder scored it as though it were fine,
# and the multiphase builder declined it as "not_expandable".

.refit_reason_data <- function() {
  set.seed(7)
  n <- 300
  data.frame(time = rexp(n, 0.2), status = rbinom(n, 1, 0.7),
             g = factor(sample(c("a", "b"), n, TRUE)), gb = runif(n),
             z = rnorm(n))
}

.refit_reason_single <- function(d) {
  hazard(survival::Surv(time, status) ~ g, data = d, dist = "weibull",
         theta = c(mu = 0.2, nu = 1, gb = 0), fit = TRUE)
}

.refit_reason_multi <- function(d) {
  # Exponential times leave phase 'early' unstarted, which hazard() warns
  # about; the fixture only needs a fit holding factor `g`.
  suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
         phases = list(
           early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                fixed = "shapes", formula = ~ g),
           constant = hzr_phase("constant", formula = ~ g)
         ),
         fit = TRUE, control = list(n_starts = 1L, maxit = 200L)))
}

test_that("a candidate refit failure warns with the refit's own error", {
  d <- .refit_reason_data()
  fit <- .refit_reason_single(d)
  expect_warning(
    step <- .hzr_stepwise_forward_step(fit, scope = "gb", data = d,
                                       criterion = "wald", slentry = 0.30),
    "candidate refit failed for gb: .*duplicated column name.*'gb'"
  )
  expect_identical(step$refit_failures, "gb")
  expect_identical(names(step$refit_failure_reasons), "gb")
  expect_match(step$refit_failure_reasons[["gb"]],
               "duplicated column name.*'gb'")
})

test_that("hzr_stepwise() keeps each refit failure's reason on $criteria", {
  d <- .refit_reason_data()
  fit <- .refit_reason_single(d)
  sw <- suppressWarnings(
    hzr_stepwise(fit, scope = "gb", data = d, criterion = "wald",
                 direction = "forward", trace = FALSE)
  )
  expect_identical(sw$criteria$refit_failures, "gb")
  expect_identical(names(sw$criteria$refit_failure_reasons), "gb")
  expect_match(sw$criteria$refit_failure_reasons[["gb"]],
               "duplicated column name.*'gb'")
})

test_that("a refit that ran but did not converge is reported as such", {
  expect_identical(
    .hzr_refit_failure_reason(list(fit = list(converged = FALSE))),
    "the refit did not converge"
  )
  expect_identical(.hzr_refit_failure_reason(simpleError("boom")), "boom")
  expect_identical(.hzr_refit_failure_reason(NULL), "the refit returned no fit")
})

test_that("the single-distribution score path declines a colliding column", {
  # It used to score `gb` against a design with two `gb` columns; `gb` then
  # won the step and failed only at the post-entry refit, which stopped the
  # screen with every other candidate untested.
  d <- .refit_reason_data()
  fit <- .refit_reason_single(d)
  step <- .hzr_stepwise_forward_step(fit, scope = c("gb", "z"), data = d,
                                     criterion = "score", slentry = 1)
  reasons <- stats::setNames(step$all_scores$reason, step$all_scores$variable)
  expect_identical(reasons[["gb"]], "duplicate_column")
  expect_true(is.na(step$all_scores$score[step$all_scores$variable == "gb"]))
  expect_false(is.na(step$all_scores$score[step$all_scores$variable == "z"]))
  expect_identical(step$refit_failures, character())
  expect_true(step$accepted)
  expect_identical(step$variable, "z")
})

test_that("the multiphase score path names the collision, not 'not_expandable'", {
  d <- .refit_reason_data()
  fit <- .refit_reason_multi(d)
  step <- .hzr_stepwise_forward_step(fit, scope = list(early = ~ gb), data = d,
                                     criterion = "score", slentry = 1)
  expect_identical(step$all_scores$reason, "duplicate_column")
  expect_match(TemporalHazard:::.hzr_score_reason_text("duplicate_column"),
               "<factor><level>", fixed = TRUE)
})

test_that("a logical candidate is checked under the name the refit gives it", {
  # model.matrix() names a logical `flag` column `flagTRUE`, so it does not
  # collide with factor `fla`'s dummy `flag`, and hazard() fits it. A check
  # on the bare name declined this strong candidate (p ~ 1e-25) silently.
  set.seed(11)
  n <- 400
  flag <- runif(n) > 0.5
  d <- data.frame(time = rexp(n, 0.2 * exp(1.2 * flag)),
                  status = rbinom(n, 1, 0.8),
                  fla = factor(sample(c("f", "g"), n, TRUE)), flag = flag,
                  w = rnorm(n))
  fit <- hazard(survival::Surv(time, status) ~ fla, data = d, dist = "weibull",
                theta = c(mu = 0.2, nu = 1, b = 0), fit = TRUE)
  step <- .hzr_stepwise_forward_step(fit, scope = c("flag", "w"), data = d,
                                     criterion = "score", slentry = 0.05)
  row <- step$all_scores[step$all_scores$variable == "flag", ]
  expect_true(is.na(row$reason))
  expect_lt(row$p_value, 1e-10)
  expect_true(step$accepted)
  expect_identical(step$variable, "flag")
})

test_that("a score run that completes still warns about a declined collision", {
  d <- .refit_reason_data()
  fit <- .refit_reason_single(d)
  msgs <- character()
  sw <- withCallingHandlers(
    # An impossible slentry ends the run on step 1 with `z` scored and
    # declined on merit, so the screen completes rather than stopping on
    # uncomputable scores, which has its own warning.
    hzr_stepwise(fit, scope = c("gb", "z"), data = d, criterion = "score",
                 direction = "forward", slentry = 1e-12, trace = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(sw$criteria$stopped_uncomputable)
  expect_identical(sw$criteria$uncomputable_reasons[["duplicate_column"]], 1L)
  expect_true(any(grepl("declined 1 candidate score\\(s\\).*<factor><level>",
                        msgs)))
})

test_that("the post-entry refit warning carries the refit's error", {
  d <- .refit_reason_data()
  fit <- .refit_reason_single(d)
  local_mocked_bindings(.hzr_refit_with_scope = function(...) stop("boom"))
  expect_warning(
    step <- .hzr_stepwise_forward_step(fit, scope = "z", data = d,
                                       criterion = "score", slentry = 1),
    "post-entry refit failed for z: boom"
  )
  expect_identical(step$refit_failure_reasons, c(z = "boom"))
})

test_that("the Wald-fallback refit warning carries the refit's error", {
  d <- .refit_reason_data()
  fit <- .refit_reason_single(d)
  local_mocked_bindings(
    .hzr_score_q = function(...) {
      list(stat = NA_real_, df = 1L, p_value = NA_real_,
           reason = "information_indefinite")
    },
    .hzr_refit_with_scope = function(...) stop("boom")
  )
  expect_warning(
    step <- .hzr_stepwise_forward_step(fit, scope = "z", data = d,
                                       criterion = "score", slentry = 1),
    "Wald-fallback refit failed for z: boom"
  )
  expect_identical(step$refit_failures, "z")
  expect_identical(step$refit_failure_reasons, c(z = "boom"))
})

test_that("the post-drop refit warning carries the refit's error", {
  d <- .refit_reason_data()
  fit <- hazard(survival::Surv(time, status) ~ gb + z, data = d,
                dist = "weibull", theta = c(0.2, 1, 0, 0),
                fit = TRUE)
  local_mocked_bindings(.hzr_refit_with_scope = function(...) stop("boom"))
  # slstay = 0: every p-value exceeds it, so the weakest term is dropped.
  expect_warning(
    step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                        slstay = 0),
    "post-drop refit failed for [a-z]+: boom"
  )
  expect_length(step$refit_failures, 1L)
  expect_identical(names(step$refit_failure_reasons), step$refit_failures)
  expect_identical(unname(step$refit_failure_reasons), "boom")
})

test_that("every step return carries refit_failures and its reasons", {
  # hzr_stepwise() reads both with `%||%`, which hid the backward step's
  # no-op and success returns omitting `refit_failure_reasons` (Copilot on
  # #302). A direct caller got a different shape from each step function.
  d <- .refit_reason_data()
  one <- hazard(survival::Surv(time, status) ~ gb, data = d, dist = "weibull",
                theta = c(0.2, 1, 0), fit = TRUE)
  two <- hazard(survival::Surv(time, status) ~ gb + z, data = d,
                dist = "weibull", theta = c(0.2, 1, 0, 0), fit = TRUE)
  wald_none <- .hzr_stepwise_forward_step(one, scope = character(), data = d,
                                          criterion = "wald")
  wald_add  <- .hzr_stepwise_forward_step(one, scope = "z", data = d,
                                          criterion = "wald", slentry = 1)
  score_add <- .hzr_stepwise_forward_step(one, scope = "z", data = d,
                                          criterion = "score", slentry = 1)
  bwd_none  <- .hzr_stepwise_backward_step(two, data = d, criterion = "wald",
                                           slstay = 1)
  bwd_drop  <- .hzr_stepwise_backward_step(two, data = d, criterion = "wald",
                                           slstay = 0)
  expect_true(wald_add$accepted)
  expect_true(bwd_drop$accepted)
  expect_false(bwd_none$accepted)

  # The backward step and the refit-based forward step share one shape; the
  # score path adds its own diagnostics on top of it. The backward step also
  # counts the removals it could not test, under the score path's names
  # (#389).
  uncomputable <- c("n_uncomputable", "uncomputable_reasons")
  expect_identical(setdiff(names(bwd_none), uncomputable), names(wald_none))
  expect_identical(setdiff(names(bwd_drop), uncomputable), names(wald_add))
  expect_true(all(uncomputable %in% names(bwd_none)))
  expect_true(all(uncomputable %in% names(bwd_drop)))
  expect_true(all(uncomputable %in% names(score_add)))
  expect_identical(bwd_drop$n_uncomputable, 0L)
  expect_true(all(names(wald_add) %in% names(score_add)))
  for (s in list(wald_none, wald_add, score_add, bwd_none, bwd_drop)) {
    expect_identical(s$refit_failure_reasons, character())
    expect_identical(s$refit_failures, character())
  }
})
