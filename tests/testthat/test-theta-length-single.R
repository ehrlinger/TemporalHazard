# A single-distribution theta must have exactly one entry per parameter: the
# distribution's shape parameters, then one coefficient per column of the
# design (#375). The only check was a lower bound against the design's column
# count, so with fit = TRUE a short theta fitted the model with a covariate
# silently dropped, and a long one returned its starting values unfitted.
# The count is the likelihood's own, as the Wald test uses: the likelihood
# ignores control$shape_param_count.

tl_data <- function() {
  set.seed(1)
  data.frame(t = rexp(40), s = rbinom(40, 1, 0.7), x = rnorm(40))
}

test_that("a short or long theta is refused on the formula interface, fitted or not", {
  d <- tl_data()
  for (fit in c(TRUE, FALSE)) {
    for (theta in list(c(1, 1), c(1, 1, 0, 5))) {
      expect_error(
        hazard(survival::Surv(t, s) ~ x, data = d, dist = "weibull",
               theta = theta, fit = fit),
        paste0("'theta' has ", length(theta), " entries, but this weibull ",
               "model takes 3"),
        fixed = TRUE, info = paste("fit =", fit, "length", length(theta))
      )
    }
  }
})

test_that("a short or long theta is refused on the vector interface, fitted or not", {
  d <- tl_data()
  x <- as.matrix(d["x"])
  for (fit in c(TRUE, FALSE)) {
    for (theta in list(c(1, 1), c(1, 1, 0, 5))) {
      expect_error(
        hazard(time = d$t, status = d$s, x = x, dist = "weibull",
               theta = theta, fit = fit),
        paste0("'theta' has ", length(theta), " entries, but this weibull ",
               "model takes 3"),
        fixed = TRUE, info = paste("fit =", fit, "length", length(theta))
      )
    }
  }
})

test_that("the refusal says what the count is made of", {
  d <- tl_data()
  msg <- tryCatch(
    hazard(survival::Surv(t, s) ~ x, data = d, dist = "weibull",
           theta = c(1, 1), fit = TRUE),
    error = conditionMessage
  )
  expect_match(msg, "2 shape parameters", fixed = TRUE)
  expect_match(msg, "1 column", fixed = TRUE)
})

test_that("an exponential intercept-only model refuses a second entry (#375)", {
  d <- tl_data()
  expect_error(
    hazard(survival::Surv(t, s) ~ 1, data = d, dist = "exponential",
           theta = c(1, 2), fit = FALSE),
    "'theta' has 2 entries, but this exponential model takes 1",
    fixed = TRUE
  )
})

test_that("a theta of the right length still fits, on both interfaces", {
  d <- tl_data()
  a <- hazard(survival::Surv(t, s) ~ x, data = d, dist = "weibull",
              theta = c(1, 1, 0), fit = TRUE)
  b <- hazard(time = d$t, status = d$s, x = as.matrix(d["x"]),
              dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_length(a$fit$theta, 3L)
  expect_identical(a$fit$objective, b$fit$objective)
  # The optimizer moved: a fit returned at its start would not.
  expect_false(isTRUE(all.equal(unname(a$fit$theta), c(1, 1, 0))))
})

test_that("time windows count one coefficient per covariate per window", {
  d <- tl_data()
  x <- as.matrix(d["x"])
  # One cut point, two windows: 2 shape parameters + 2 coefficients.
  expect_error(
    hazard(time = d$t, status = d$s, x = x, time_windows = 1,
           dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
    "'theta' has 3 entries, but this weibull model takes 4",
    fixed = TRUE
  )
  fit <- hazard(time = d$t, status = d$s, x = x, time_windows = 1,
                dist = "weibull", theta = c(1, 1, 0, 0), fit = TRUE)
  expect_length(fit$fit$theta, 4L)
})

test_that("control$shape_param_count does not change the count the likelihood uses", {
  # The likelihood always takes the distribution's own shape parameters; the
  # control element is read only by the stepwise refit and the score test.
  d <- tl_data()
  expect_error(
    hazard(survival::Surv(t, s) ~ x, data = d, dist = "weibull",
           theta = c(1, 1, 0, 0), fit = TRUE,
           control = list(shape_param_count = 3)),
    "'theta' has 4 entries, but this weibull model takes 3",
    fixed = TRUE
  )
})
