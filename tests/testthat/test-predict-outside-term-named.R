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
  want <- paste0("term 'zz' of the model does not give one value per row ",
                 "of 'newdata'")
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
  expect_match(msg, "term 'zz' of the model does not give one value per row",
               fixed = TRUE)
  expect_match(msg, "move those into `data` as columns and refit", fixed = TRUE)
})

test_that("an outside part inside a larger term names that term", {
  o <- ot_setup()
  zz <- o$zz
  fit <- hazard(survival::Surv(t, s) ~ I(mal + zz), data = o$d,
                dist = "weibull", theta = c(mu = 0.5, nu = 1, 0), fit = TRUE)
  expect_error(
    predict(fit, newdata = data.frame(mal = c(0, 1)), type = "linear_predictor"),
    "term 'I(mal + zz)' of the model does not give one value per row",
    fixed = TRUE
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

test_that("a formula fit saved without its design refuses an extra newdata column", {
  # Not #409's own change: this pins a refusal that had NO test at all
  # (grep for its text in tests/ found 0 hits, 1 in R/). A genuine
  # pre-design formula fit -- no x_design, no kept frame, no call
  # environment, so .hzr_recover_x_design() cannot rebuild -- still refuses
  # a newdata column it cannot tell from a formula variable.
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

# Three shapes the first version of the pre-check got wrong, found by an
# adversarial review of this PR before it was merged.

ot_long <- paste0("I(log(baseline_creatinine_mg_dl) - ",
                  "log(followup_creatinine_mg_dl) + zz)")

ot_wide_setup <- function() {
  set.seed(12)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7),
                  baseline_creatinine_mg_dl = runif(60, 0.6, 2),
                  followup_creatinine_mg_dl = runif(60, 0.6, 2))
  list(d = d, zz = rnorm(60),
       nd = d[1:2, c("baseline_creatinine_mg_dl",
                     "followup_creatinine_mg_dl")])
}

ot_fit <- function(rhs, d, theta) {
  f <- stats::as.formula(paste("survival::Surv(t, s) ~", rhs))
  environment(f) <- parent.frame()
  hazard(f, data = d, dist = "weibull", theta = theta, fit = TRUE)
}

test_that("a term too long for deparse()'s default cutoff is still named", {
  # The pre-check matched the terms matrix by a deparsed string. deparse()
  # breaks at 60 characters and indents the continuation, so the key never
  # equalled the rowname, no term was found, and the refusal fell back to
  # the unnamed backstop -- silently, for exactly the terms #409 is about.
  o <- ot_wide_setup()
  zz <- o$zz
  # Assert the property the test depends on: a term this test can only
  # exercise if deparse() really does break it over more than one line.
  # A 61-character term does not, and an earlier draft of this test
  # passed on the unfixed code for exactly that reason.
  expect_gt(length(deparse(str2lang(ot_long))), 1L)
  msg <- tryCatch(
    predict(ot_fit(ot_long, o$d, c(mu = 0.5, nu = 1, 0)), newdata = o$nd,
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, paste0("term '", ot_long, "'"), fixed = TRUE)
  expect_false(grepl("The design rebuilt for", msg, fixed = TRUE))
})

test_that("every outside term is named, not only the ones deparse keeps short", {
  # Naming one of two sends the user round the loop twice: they move `yy`
  # into `data`, refit, and are refused again for the term not named.
  o <- ot_wide_setup()
  zz <- o$zz
  yy <- rnorm(60)
  msg <- tryCatch(
    predict(ot_fit(paste("yy +", ot_long), o$d, c(mu = 0.5, nu = 1, 0, 0)),
            newdata = o$nd, type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "'yy'", fixed = TRUE)
  expect_match(msg, paste0("'", ot_long, "'"), fixed = TRUE)
})

test_that("a length-changing term of a data column is not blamed on outside data", {
  # A row count other than newdata's does NOT establish that the values
  # came from outside `data`: unique() and na.omit() change the length of a
  # `data` column. The old message asserted the outside-data cause as fact
  # and told the user to move a column that is already in `data`.
  set.seed(13)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = seq(40, 99))
  stopifnot(length(unique(d$age)) == nrow(d))
  fit <- ot_fit("I(unique(age))", d, c(mu = 0.5, nu = 1, 0))
  msg <- tryCatch(
    predict(fit, newdata = data.frame(age = c(30, 30, 50)),
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "term 'I(unique(age))'", fixed = TRUE)
  expect_match(msg, "does not give one value per row of 'newdata'", fixed = TRUE)
  expect_match(msg, "length-changing", fixed = TRUE)
  # The unconditional remedy is what made it unfollowable here.
  expect_false(grepl("Move them into `data` as columns and refit.", msg,
                     fixed = TRUE))
})
