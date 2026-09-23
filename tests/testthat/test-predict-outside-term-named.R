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

test_that("the backstop's own message propagates when nothing can be named", {
  # The backstop is what REFUSES a short design; naming the term is a
  # diagnosis run afterwards. When the diagnosis names nothing -- here by
  # mocking it to return no term -- the refusal must still reach the user,
  # as the backstop wrote it, rather than being swallowed.
  local_mocked_bindings(.hzr_outside_rows_term = function(...) character(0))
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

# Codex, reviewing this PR at c5144c35, found that the pre-check did not
# merely name a term: by evaluating each variable on its own, before
# model.frame(), it changed which fits predict at all and which errors
# reach the caller. Both shapes below come from that review. The pre-check
# semantics arrived with the original #409 commit, not with the fix to it.

test_that("a term that assigns into model.frame's shared mask still predicts", {
  # model.frame() evaluates every term in ONE mask, so `I(zz <- age)` is
  # what `I(zz^2)` reads. Evaluated separately, the second term finds the
  # unrelated `zz` of the formula's environment instead, and a prediction
  # main makes was refused.
  set.seed(21)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60))
  zz <- rnorm(60)
  f <- survival::Surv(t, s) ~ I(zz <- age) + I(zz^2)
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0, 0),
           fit = TRUE)
  )
  p <- predict(fit, newdata = d[1:2, "age", drop = FALSE],
               type = "linear_predictor")
  # The newdata path must agree with the fitted path on the same two rows.
  expect_equal(unname(p),
               unname(predict(fit, type = "linear_predictor")[1:2]))
})

test_that("an error raised while building the design reaches the caller", {
  # The pre-check evaluated each variable in a tryCatch and discarded the
  # error, then model.frame() evaluated the term again. A function that
  # fails once and then succeeds therefore returned numbers where the
  # error should have propagated.
  set.seed(21)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60))
  boom <- FALSE
  ff <- function(x) {
    if (boom) {
      boom <<- FALSE
      stop("formula boom")
    }
    x
  }
  f <- survival::Surv(t, s) ~ I(ff(age))
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0),
           fit = TRUE)
  )
  boom <- TRUE
  msg <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    error = conditionMessage
  )
  # identical(), not a match: the caller's own condition, unchanged.
  expect_identical(msg, "formula boom")
})

test_that("the design is built once per predict(newdata = ), as before #409", {
  # The pre-check evaluated every model-frame variable an extra time. A
  # term with a side effect, or an expensive one, paid for the diagnosis on
  # every call, including the calls that succeed.
  set.seed(21)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60))
  n_eval <- 0L
  sideg <- function(x) {
    n_eval <<- n_eval + 1L
    x
  }
  f <- survival::Surv(t, s) ~ sideg(age)
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0),
           fit = TRUE)
  )
  n_eval <- 0L
  invisible(predict(fit, newdata = data.frame(age = c(0.4, 0.5)),
                    type = "linear_predictor"))
  # Two: the design, and the row-shifted copy the equivariance check builds.
  expect_identical(n_eval, 2L)
})
