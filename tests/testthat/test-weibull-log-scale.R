# The Weibull likelihood on the log scale (#566, PR 2).
#
# .hzr_logl_weibull() formed mu^nu, (mu * t)^nu and t^(nu - 1) on the natural
# scale. Each can underflow or overflow where the log-likelihood is an ordinary
# number. The oracle in each case is the log-scale value, written out here.

.ls_eval <- function(obj, theta) {
  msgs <- character(0)
  ll <- withCallingHandlers(
    hzr_evaluate(obj, theta)$logLik,
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(ll = as.numeric(ll), msgs = msgs)
}

test_that("a left-censored row whose H is subnormal is read on the log scale (#566 acceptance 1)", {
  obj <- hazard(time = 2, status = -1, time_upper = 2, dist = "weibull",
                theta = c(1e-81, 4), fit = FALSE)
  # H = (2e-81)^4 = 1.6e-323 is subnormal, so log(1 - exp(-H)) = log(H) to
  # working precision. An oracle that forms H by exp() first loses the same
  # digits the defect does, so compare against log(H) directly.
  log_h <- 4 * (log(1e-81) + log(2))
  expect_lt(exp(log_h), .Machine$double.xmin)   # premise: H is subnormal
  got <- .ls_eval(obj, c(1e-81, 4))
  # main e7b621c5 returned -743.341, with no warning.
  expect_length(got$msgs, 0L)
  expect_equal(got$ll / log_h, 1, tolerance = 1e-12)
})

test_that("an exact event whose t^(nu - 1) underflows is read on the log scale (#566 acceptance 2)", {
  obj <- hazard(time = 1e-180, status = 1, dist = "weibull",
                theta = c(1e100, 3), fit = FALSE)
  truth <- log(3) + 3 * log(1e100) + 2 * log(1e-180) -
    exp(3 * (log(1e100) + log(1e-180)))
  expect_equal(1e-180^2, 0)   # premise: the natural-scale factor underflows
  got <- .ls_eval(obj, c(1e100, 3))
  # main e7b621c5 returned -Inf with an hzr_evaluate_not_finite warning.
  expect_length(got$msgs, 0L)
  expect_equal(got$ll / truth, 1, tolerance = 1e-12)
})

test_that("a zero-weight event row whose mu^nu underflows leaves the likelihood unchanged (#566 acceptance 3)", {
  both <- hazard(time = c(1, 1), status = c(1, 0), weights = c(0, 1),
                 dist = "weibull", theta = c(1e-200, 2), fit = FALSE)
  alone <- hazard(time = 1, status = 0, dist = "weibull",
                  theta = c(1e-200, 2), fit = FALSE)
  expect_equal(1e-200^2, 0)   # premise: mu^nu underflows
  got_both <- .ls_eval(both, c(1e-200, 2))
  got_alone <- .ls_eval(alone, c(1e-200, 2))
  # main e7b621c5 returned -Inf for `both` (0 * log(0) = NaN on the
  # zero-weight row), with the not-finite warning.
  expect_length(got_both$msgs, 0L)
  expect_identical(got_both$ll, got_alone$ll)
  # Known positive: with weight 1 the event row does move the likelihood.
  # (At mu = 1e-200 hzr_evaluate() refuses a counted event whose mu^nu
  # cannot be represented (#573), so this is read at mu = 1e-100.)
  one <- hazard(time = c(1, 1), status = c(1, 0), weights = c(1, 1),
                dist = "weibull", theta = c(1e-100, 2), fit = FALSE)
  zero <- hazard(time = c(1, 1), status = c(1, 0), weights = c(0, 1),
                 dist = "weibull", theta = c(1e-100, 2), fit = FALSE)
  got_one <- .ls_eval(one, c(1e-100, 2))
  got_zero <- .ls_eval(zero, c(1e-100, 2))
  expect_true(is.finite(got_one$ll))
  expect_equal(got_one$ll - got_zero$ll, log(2) + 2 * log(1e-100),
               tolerance = 1e-12)
})
