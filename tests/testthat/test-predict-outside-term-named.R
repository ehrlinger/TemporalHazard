# predict(newdata =) refuses a model term that takes row-level values from
# outside `data`. When newdata's row count differed from the fit's, the
# refusal came from model.frame() -- "variable lengths differ (found for
# 'mal')", naming another variable -- or from the row-count backstop, which
# named no term. It now names the term before model.frame() runs, in all
# three shapes (#409).

ot_setup <- function() {
  set.seed(11)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), mal = rbinom(60, 1, 0.4))
  list(d = d, zz = rnorm(60))
}

test_that("an outside term beside a data column is named, whichever columns newdata has", {
  o <- ot_setup()
  zz <- o$zz
  fit <- hazard(survival::Surv(t, s) ~ zz + mal, data = o$d, dist = "weibull",
                theta = c(mu = 0.5, nu = 1, 0, 0), fit = TRUE)
  want <- paste0("term 'zz' of the model uses row-level values taken from ",
                 "outside `data`")
  for (nd in list(data.frame(zz = c(-1, 1), mal = c(0, 1)),
                  data.frame(mal = c(0, 1)))) {
    msg <- tryCatch(predict(fit, newdata = nd, type = "linear_predictor"),
                    error = conditionMessage)
    expect_match(msg, want, fixed = TRUE)
    expect_false(grepl("variable lengths differ", msg, fixed = TRUE))
  }
})

test_that("an outside term alone is named, not just counted", {
  o <- ot_setup()
  zz <- o$zz
  fit <- hazard(survival::Surv(t, s) ~ zz, data = o$d, dist = "weibull",
                theta = c(mu = 0.5, nu = 1, 0), fit = TRUE)
  msg <- tryCatch(
    predict(fit, newdata = data.frame(zz = c(-1, 1)), type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "term 'zz' of the model uses row-level values", fixed = TRUE)
  expect_match(msg, "Move them into `data` as columns and refit.", fixed = TRUE)
})

test_that("an outside part inside a larger term names that term", {
  o <- ot_setup()
  zz <- o$zz
  fit <- hazard(survival::Surv(t, s) ~ I(mal + zz), data = o$d,
                dist = "weibull", theta = c(mu = 0.5, nu = 1, 0), fit = TRUE)
  expect_error(
    predict(fit, newdata = data.frame(mal = c(0, 1)), type = "linear_predictor"),
    "term 'I(mal + zz)' of the model uses row-level values", fixed = TRUE
  )
})

test_that("constants and knots from outside data still predict at new rows", {
  # Only a variable whose rows differ from newdata's is refused: a scalar
  # cutoff, or a knot vector used as an argument, is not row-level.
  o <- ot_setup()
  d <- o$d
  d$age <- seq(20, 80, length.out = 60)
  cutoff <- 50
  knots <- c(40, 60)
  fit <- hazard(survival::Surv(t, s) ~ I(age > cutoff) + splines::ns(age, knots = knots),
                data = d, dist = "weibull",
                theta = c(mu = 0.5, nu = 1, rep(0, 4)), fit = TRUE)
  nd <- data.frame(age = c(30, 55))
  p <- predict(fit, newdata = nd, type = "linear_predictor")
  expect_length(p, 2L)
  expect_true(all(is.finite(p)))
})

test_that("the row-count backstop still refuses when the pre-check is bypassed", {
  # The pre-check now catches every wrong-rows shape the suite has, so this
  # is the only test that reaches .hzr_check_design_rows(). It proves the
  # backstop is live defence, not code that merely goes unreached.
  local_mocked_bindings(.hzr_refuse_outside_rows = function(...) invisible(NULL))
  o <- ot_setup()
  zz <- o$zz
  fit <- hazard(survival::Surv(t, s) ~ zz, data = o$d, dist = "weibull",
                theta = c(mu = 0.5, nu = 1, 0), fit = TRUE)
  expect_error(
    predict(fit, newdata = data.frame(zz = c(-1, 1)), type = "linear_predictor"),
    "The design rebuilt for the model has 60 rows for 2 row(s) of 'newdata'",
    fixed = TRUE
  )
})

test_that("a vector fit made through a wrapper's formula argument takes newdata like any vector fit (#406)", {
  # A wrapper passing `formula = fml` stores the symbol `fml` in the call
  # even when fml is NULL. The fit is a vector-interface fit (its call has
  # `time =`), but predict() read the non-NULL call$formula as a formula
  # fit saved before its design, and refused a newdata with an extra
  # column. Part of #406; the other sites are in stream C's PR.
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age")])
  wrap <- function(dat) {
    fml <- NULL
    hazard(formula = fml, time = dat$int_dead, status = dat$dead,
           x = cbind(age = dat$age), dist = "weibull",
           theta = c(0.1, 1, 0), fit = TRUE)
  }
  plain <- hazard(time = d$int_dead, status = d$dead, x = cbind(age = d$age),
                  dist = "weibull", theta = c(0.1, 1, 0), fit = TRUE)
  nd <- data.frame(time = c(1, 5), age = c(50, 70), extra = 1)
  expect_identical(predict(wrap(d), newdata = nd, type = "survival"),
                   predict(plain, newdata = nd, type = "survival"))
})

test_that("a formula fit saved without its design still refuses an extra newdata column", {
  # The control for the #406 edit: a genuine pre-design formula fit (no
  # x_design, no kept data frame, no call environment, so the design cannot
  # be rebuilt) still refuses a newdata column it cannot tell from a formula
  # variable. Its call has no `time =`.
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age")])
  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = d,
                dist = "weibull", theta = c(0.1, 1, 0), fit = TRUE)
  fit$data$x_design <- NULL
  fit$fit$x_design <- NULL
  fit$data$frame <- NULL
  fit$call_env <- NULL
  expect_error(
    predict(fit, newdata = data.frame(time = c(1, 5), age = c(50, 70),
                                      extra = 1), type = "survival"),
    "saved by an earlier version of TemporalHazard, without a stored formula design"
  )
})

test_that("a formula fit whose call also named time = is still a formula fit (#406)", {
  # hazard() fits the formula when one is given and ignores `time =` beside
  # it, so `time` in the call does not make a fit vector-interface. Saved
  # without its design, such a fit must still refuse an extra column.
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age")])
  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = d,
                time = d$int_dead * 3, dist = "weibull",
                theta = c(0.1, 1, 0), fit = TRUE)
  fit$data$x_design <- NULL
  fit$fit$x_design <- NULL
  fit$data$frame <- NULL
  expect_error(
    predict(fit, newdata = data.frame(time = c(1, 5), age = c(50, 70),
                                      extra = 1), type = "survival"),
    "saved by an earlier version of TemporalHazard, without a stored formula design"
  )
})
