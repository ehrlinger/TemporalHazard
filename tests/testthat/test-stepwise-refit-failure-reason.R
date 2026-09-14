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
