# test-predict-stored-theta-length.R -- predict() validates the STORED theta
# against the stored design, for every prediction type, before any prediction
# arithmetic (#375, Codex review of #422).
#
# A legacy or hand-edited object whose theta is longer than its design allows
# reached the arithmetic unchecked. `hazard` and `linear_predictor` refused
# downstream, where the design is multiplied as a matrix, but `survival` and
# `cumulative_hazard` recycled the surplus coefficients into an OUTER PRODUCT
# and returned 2n values for n rows, with no error: six predictions for three
# training rows. `hzr_evaluate()` refused the same object. Measured identical
# on main and on the branch before the fix, so this is a gap #422 left open
# rather than one it introduced.

# 40 rows, not the 3 the review's reproduction used. At 3 the fit is
# ill-conditioned (rcond about 4e-14) and predicts Inf and NaN for perfectly
# legitimate reasons, so a known positive asserting finite predictions there
# would be testing the fixture rather than the code. The defect reproduces at
# any n: the surplus coefficients give 2n values for n rows.
ptl_fit <- function() {
  set.seed(1)
  n <- 40
  d <- data.frame(t = stats::rexp(n) + 0.1,
                  s = stats::rbinom(n, 1, 0.7), x = stats::rnorm(n))
  suppressWarnings(hazard(survival::Surv(t, s) ~ x, data = d, dist = "weibull",
                          theta = c(1, 1, 0), fit = TRUE))
}

test_that("predict() refuses a stored theta longer than the design, every type", {
  fit <- ptl_fit()
  fit$fit$theta <- c(1, 1, 0, 5)   # 4 entries for a 3-parameter model
  for (ty in c("hazard", "linear_predictor", "survival", "cumulative_hazard")) {
    expect_error(
      predict(fit, type = ty),
      "'theta' has 4 entries, but this weibull model takes 3",
      fixed = TRUE, info = ty
    )
  }
  # And with newdata, where the surplus produced 4 values for 2 rows.
  nd <- data.frame(time = c(1, 2), x = c(0, 1))
  for (ty in c("survival", "cumulative_hazard")) {
    expect_error(predict(fit, newdata = nd, type = ty),
                 "but this weibull model takes 3", fixed = TRUE, info = ty)
  }
})

test_that("predict() still refuses a stored theta shorter than the design", {
  # The control: the short case already refused before the fix, downstream.
  # It must still refuse, and now name the counts rather than the columns.
  fit <- ptl_fit()
  fit$fit$theta <- c(1, 1)
  expect_error(predict(fit, type = "survival"),
               "'theta' has 2 entries, but this weibull model takes 3",
               fixed = TRUE)
})

test_that("an untampered fit predicts normally, one value per row", {
  # The known positive. Without this the tests above pass for a predict()
  # that refuses everything.
  fit <- ptl_fit()
  for (ty in c("hazard", "survival", "cumulative_hazard")) {
    p <- predict(fit, type = ty)
    expect_length(p, 40L)
    expect_true(all(is.finite(p)), info = ty)
  }
  nd <- data.frame(time = c(1, 2), x = c(0, 1))
  expect_length(predict(fit, newdata = nd, type = "survival"), 2L)
})
