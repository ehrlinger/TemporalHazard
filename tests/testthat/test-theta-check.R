# A user-supplied theta is checked before the likelihood sees it, by one
# helper, at every entry point that takes one: hazard(), hzr_evaluate() and
# predict(). The likelihood's own sentinels are internal: a Weibull scale or
# shape <= 0 made .hzr_logl_weibull() return Inf, which reached the user as
# `logLik = Inf` from hzr_evaluate() (read as the best possible fit), and as
# "non-finite value supplied by optim" from hazard() (#383).

tc_data <- function() {
  set.seed(1)
  data.frame(t = rexp(40), s = rbinom(40, 1, 0.7), x = rnorm(40))
}

tc_fit <- function(...) {
  hazard(survival::Surv(t, s) ~ x, data = tc_data(), dist = "weibull",
         theta = c(1, 1, 0), fit = TRUE, ...)
}

test_that("hzr_evaluate() refuses a non-positive Weibull scale or shape, naming it", {
  f <- tc_fit()
  expect_error(hzr_evaluate(f, c(0, 1, 0)), "Weibull scale mu = 0",
               fixed = TRUE)
  expect_error(hzr_evaluate(f, c(1, 0, 0)), "Weibull shape nu = 0",
               fixed = TRUE)
  expect_error(hzr_evaluate(f, c(-1, -2, 0)),
               "Weibull scale mu = -1 and shape nu = -2", fixed = TRUE)
  # Control: the fit's own estimates evaluate to the fit's objective.
  expect_equal(hzr_evaluate(f, unname(f$fit$theta))$logLik,
               f$fit$objective, tolerance = 1e-8)
})

test_that("hazard() refuses a non-positive Weibull start, naming it (#383)", {
  expect_error(
    hazard(survival::Surv(t, s) ~ x, data = tc_data(), dist = "weibull",
           theta = c(0, 1, 0), fit = TRUE),
    "Weibull scale mu = 0", fixed = TRUE
  )
  expect_error(
    hazard(survival::Surv(t, s) ~ x, data = tc_data(), dist = "weibull",
           theta = c(1, -1, 0), fit = TRUE),
    "Weibull shape nu = -1", fixed = TRUE
  )
})

test_that("predict() gives the same refusal for a non-positive Weibull theta", {
  f <- tc_fit()
  f$fit$theta[2] <- 0
  expect_error(predict(f, newdata = data.frame(time = 1, x = 0),
                       type = "survival"),
               "Weibull shape nu = 0", fixed = TRUE)
})

test_that("hzr_evaluate() counts parameters as the likelihood does, not as control says", {
  # control$shape_param_count is read by the score test and the stepwise
  # refit, and deliberately ignored by the likelihood (as by wald.R).
  # hzr_evaluate() evaluates the likelihood, so it must count as the
  # likelihood does. It refused the fit's own theta and accepted a longer
  # one, returning another model's log-likelihood.
  f <- tc_fit(control = list(shape_param_count = 3))
  expect_equal(hzr_evaluate(f, unname(f$fit$theta))$logLik,
               f$fit$objective, tolerance = 1e-8)
  expect_error(hzr_evaluate(f, c(unname(f$fit$theta), 0)),
               "'theta' has 4 entries, but this weibull model takes 3",
               fixed = TRUE)
})

test_that("an unfitted Weibull model refuses a non-positive scale or shape too (#375)", {
  # An unfitted object with such a theta can never be predicted from.
  expect_error(
    hazard(survival::Surv(t, s) ~ x, data = tc_data(), dist = "weibull",
           theta = c(0, 1, 0), fit = FALSE),
    "Weibull scale mu = 0", fixed = TRUE
  )
})

test_that("an unsupported distribution is not told it takes 0 parameters (#375)", {
  # The length check covers the four distributions whose count is known;
  # anything else reaches the refusal that already exists for it.
  d <- tc_data()
  msg <- tryCatch(hazard(time = d$t, status = d$s, dist = "gompertz",
                         theta = c(1, 1), fit = TRUE),
                  error = conditionMessage)
  expect_match(msg, "not yet supported for fitting", fixed = TRUE)
  expect_false(grepl("takes 0", msg, fixed = TRUE))
  expect_s3_class(hazard(time = d$t, status = d$s, dist = "gompertz",
                         theta = c(1, 1), fit = FALSE), "hazard")
})
