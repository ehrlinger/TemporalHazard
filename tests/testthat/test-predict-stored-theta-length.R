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

test_that("the check does not fire on a time-windowed model, or one with no stored design", {
  # The two cases where the new count is not simply ncol(x), and where getting
  # it wrong refuses a legitimate object rather than letting a bad one through.
  #
  # (a) A time-windowed fit: theta carries one coefficient per covariate PER
  # WINDOW, so the count must be the expanded one. Dropping the expansion
  # would leave the suite green while predict() refused every windowed fit,
  # which is why this asserts the non-firing directly.
  set.seed(4)
  n <- 60
  d <- data.frame(t = stats::rexp(n) + 0.1,
                  s = stats::rbinom(n, 1, 0.7), x = stats::rnorm(n))
  win <- suppressWarnings(hazard(
    survival::Surv(t, s) ~ x, data = d, dist = "weibull",
    time_windows = c(0.5), theta = c(1, 1, 0, 0), fit = TRUE
  ))
  expect_gt(length(win$fit$theta), 3L)   # premise: the design really expanded
  p <- predict(win, type = "survival")
  expect_length(p, n)

  # (b) An object that stored no design but carries a coefficient: position is
  # the only mapping, and it is a documented path, so the check must stay out
  # of its way.
  nodes <- hazard(time = d$t, status = d$s, dist = "exponential", theta = -4)
  nodes$fit$theta <- c(-4, 0.01)
  # suppressWarnings: the object is unfitted, so predict() says its numbers
  # come from starting values. True, and orthogonal to the length check this
  # test is about, so it is kept out of the suite's warning multiset.
  expect_length(
    suppressWarnings(
      predict(nodes, newdata = data.frame(time = c(1, 2), age = 70),
              type = "linear_predictor")
    ),
    2L
  )
})

test_that("predict() warns when it maps newdata columns by position (#422)", {
  # John's decision (2026-09-23): keep the behaviour, say it out loud. An
  # object with no stored design matches newdata's columns to its coefficients
  # in the order supplied, so reordering or renaming them changes the
  # predictions with nothing else changing. The warning names the count, shows
  # the columns it used, and says what to do.
  set.seed(5)
  n <- 30
  d <- data.frame(t = stats::rexp(n) + 0.1, s = stats::rbinom(n, 1, 0.7))
  obj <- hazard(time = d$t, status = d$s, dist = "exponential", theta = -4)
  obj$fit$theta <- c(-4, 0.01)
  nd <- data.frame(time = c(1, 2), age = 70)
  expect_warning(
    suppressWarnings(  # the separate fit = FALSE notice is not what this pins
      withCallingHandlers(
        predict(obj, newdata = nd, type = "linear_predictor"),
        warning = function(w) {
          if (!grepl("BY POSITION", conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
          }
        }
      )
    ),
    NA
  )
  # Pin the text directly, which is what a user reads.
  msgs <- character(0)
  withCallingHandlers(
    predict(obj, newdata = nd, type = "linear_predictor"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  pos <- msgs[grepl("BY POSITION", msgs, fixed = TRUE)]
  expect_length(pos, 1L)
  expect_match(pos, "stored no design matrix", fixed = TRUE)
  expect_match(pos, "1 covariate coefficient BY POSITION", fixed = TRUE)
  expect_match(pos, "age", fixed = TRUE)
  expect_match(pos, "Refit the model", fixed = TRUE)

  # Control: a fit WITH a stored design maps by name and must not warn.
  set.seed(6)
  d2 <- data.frame(t = stats::rexp(n) + 0.1, s = stats::rbinom(n, 1, 0.7),
                   x = stats::rnorm(n))
  ok <- suppressWarnings(hazard(survival::Surv(t, s) ~ x, data = d2,
                                dist = "weibull", theta = c(1, 1, 0),
                                fit = TRUE))
  msgs2 <- character(0)
  withCallingHandlers(
    predict(ok, newdata = data.frame(time = c(1, 2), x = c(0, 1)),
            type = "linear_predictor"),
    warning = function(w) {
      msgs2 <<- c(msgs2, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(msgs2[grepl("BY POSITION", msgs2, fixed = TRUE)], 0L)
})
