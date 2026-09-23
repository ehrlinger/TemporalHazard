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
  local_mocked_bindings(
    .hzr_outside_rows_wrong = function(...) {
      list(idx = integer(0), vals = list(), term = character(0))
    }
  )
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
  # The fit is suppressed because a contrived collinear formula warns about
  # conditioning; the predict() under test is NOT, because a warning from
  # the rebuild path is exactly what this test has to be able to see.
  expect_no_warning(
    p <- predict(fit, newdata = d[1:2, "age", drop = FALSE],
                 type = "linear_predictor")
  )
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

# Codex, reviewing this PR at c0b2c9e7, found that the diagnosis could replace
# a condition it had no business touching: when an EARLIER term raises and a
# LATER term happens to be row-mismatched, the caller lost their own error and
# got the row-count refusal, blaming a term that had nothing to do with the
# failure. The diagnosis is now substituted only for failures that ARE about
# row counts, identified structurally rather than by message text.

test_that("an unrelated error survives a row-mismatched term beside it", {
  # The shape the raising-once test above cannot reach: it has no second,
  # mismatched term, so the diagnosis finds nothing and re-raises anyway.
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  zz <- rnorm(60)
  armed <- FALSE
  ff <- function(x) {
    if (armed) {
      stop(structure(class = c("ff_boom", "error", "condition"),
                     list(message = "ff exploded", call = NULL)))
    }
    x
  }
  f <- survival::Surv(t, s) ~ I(ff(age)) + zz
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0, 0),
           fit = TRUE)
  )
  armed <- TRUE
  e <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    condition = function(e) e
  )
  # The caller's class, not ours: a tryCatch(ff_boom = ) must still fire.
  expect_identical(class(e), c("ff_boom", "error", "condition"))
  expect_identical(conditionMessage(e), "ff exploded")
})

test_that("a user error with no call survives too, beside a mismatched term", {
  # conditionCall() is what separates a design-build failure from a user's,
  # and this condition has none -- so the only reason the row-count refusal
  # does not swallow it is that OUR backstop carries a class instead.
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  zz <- rnorm(60)
  armed <- FALSE
  ff <- function(x) {
    if (armed) stop("plain user error", call. = FALSE)
    x
  }
  f <- survival::Surv(t, s) ~ I(ff(age)) + zz
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0, 0),
           fit = TRUE)
  )
  armed <- TRUE
  msg <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_identical(msg, "plain user error")
})

test_that("every row-count failure family still names its term", {
  # The three underlying errors the named refusal rests on, measured by
  # neutralising the diagnosis: model.frame.default's "variable lengths
  # differ", model.matrix.default's "length of 'dimnames'", and our own
  # backstop. If the predicate stops recognising any one of them, #409's
  # whole point is lost for that family -- silently, since the user still
  # gets an error.
  set.seed(11)
  n <- 60
  d <- data.frame(t = rexp(n), s = rbinom(n, 1, 0.7), mal = rbinom(n, 1, 0.4),
                  age = seq(40, 40 + n - 1))
  zz <- rnorm(n)
  mk <- function(rhs, th) {
    f <- stats::as.formula(paste("survival::Surv(t, s) ~", rhs))
    environment(f) <- parent.frame()
    suppressWarnings(hazard(f, data = d, dist = "weibull", theta = th,
                            fit = TRUE))
  }
  named <- "does not give one value per row of 'newdata'"
  # model.frame.default: variable lengths differ
  expect_error(
    predict(mk("zz + mal", c(0.5, 1, 0, 0)),
            newdata = d[1:2, "mal", drop = FALSE], type = "linear_predictor"),
    named, fixed = TRUE
  )
  # model.matrix.default: length of 'dimnames'
  expect_error(
    predict(mk("I(unique(age))", c(0.5, 1, 0)),
            newdata = data.frame(age = c(30, 30, 50)),
            type = "linear_predictor"),
    named, fixed = TRUE
  )
  # our own backstop, which carries a class rather than a message to match
  expect_error(
    predict(mk("zz", c(0.5, 1, 0)), newdata = data.frame(zz = c(-1, 1)),
            type = "linear_predictor"),
    named, fixed = TRUE
  )
})

test_that("a term that reads an earlier term's assignment is not named", {
  # model.frame() evaluates every term in ONE mask, in order, so
  # `I(zz <- age)` is what `I(zz^2)` reads. Evaluated in separate masks the
  # diagnosis saw the formula environment's unrelated 60-length `zz` and
  # named `I(zz^2)` too -- sending the user to repair a term that is fine.
  # Base R's own message names exactly `yy` here, so this shape was one
  # where #409's naming was WORSE than what it replaced (#430 review).
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  zz <- rnorm(60)
  yy <- rnorm(60)
  f <- survival::Surv(t, s) ~ I(zz <- age) + I(zz^2) + yy
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull",
           theta = c(mu = 0.5, nu = 1, 0, 0, 0), fit = TRUE)
  )
  msg <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "'yy'", fixed = TRUE)
  expect_false(grepl("I(zz^2)", msg, fixed = TRUE))
})

test_that("a nested model.frame() failure keeps the caller's own error", {
  # Was pinned as a LIMITATION by #430: `conditionCall()` reports
  # `model.frame.default` for both our own frame assembly and one a user's
  # term performed itself, so the diagnosis replaced the real failure with
  # blame for whichever term was row-mismatched. #446 stops classifying the
  # failure and verifies CAUSATION instead, so this now propagates.
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  zz <- rnorm(60)
  armed <- FALSE
  ff <- function(x) {
    if (armed) {
      # Raised by a nested model.frame(), UNCAUGHT, so its conditionCall is
      # `model.frame.default` -- the same callee our own frame assembly
      # reports, which is exactly what main misattributes.
      return(stats::model.frame(~ x + I(1:3)))
    }
    x
  }
  f <- survival::Surv(t, s) ~ I(ff(age)) + zz
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0, 0),
           fit = TRUE)
  )
  armed <- TRUE
  cond <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    condition = function(e) e
  )
  # Assert the variable the NESTED call names. The build's own generic
  # failure also reads "variable lengths differ", but it names `zz`, so a
  # bare "variable lengths differ" cannot tell the two apart -- and a
  # CLASSED condition cannot either: giving it `call = NULL` makes main's
  # old classifier fall through and re-raise, so the test would pass on
  # main and pin nothing. Measured both ways before settling here.
  expect_match(conditionMessage(cond), "found for 'I(1:3)'", fixed = TRUE)
  expect_false(grepl("does not give one value per row",
                     conditionMessage(cond), fixed = TRUE))
})

test_that("a term the model frame rejects keeps its own error, not row blame", {
  # The second shape #446 covers, and it needs no nested model.frame() at
  # all: `ff()` returns a LIST, our own outer frame rejects it, and the
  # diagnosis used to blame `zz` for being row-mismatched. The remedy it
  # prescribed could not even be followed -- `zz` is not a `data` column, so
  # `newdata` cannot supply it.
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  zz <- rnorm(60)
  ff <- function(x) if (length(x) == 2L) as.list(x) else x
  f <- survival::Surv(t, s) ~ I(ff(age)) + zz
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0, 0),
           fit = TRUE)
  )
  msg <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "invalid type (list) for variable", fixed = TRUE)
  expect_false(grepl("does not give one value per row", msg, fixed = TRUE))
})


# The probe matrix, as tests. Term ORDER is permuted deliberately: the
# `drop.terms()` regression was invisible to a four-row table whose every
# shape put the outside variable last, because `drop.terms()` indexed
# `predvars` positionally while `terms()` orders `variables` by first
# appearance and `term.labels` by interaction order.
#
# Which of the three blocks below actually DISCRIMINATE, measured by running
# this file against both implementations rather than assumed:
#   "used only inside an interaction"  -- RED on the drop.terms head
#   "a length-changing function"       -- RED on the drop.terms head
#   "whatever ORDER it is written in"  -- passes on BOTH, and on main
# The third catches neither defect: a plain additive term cannot trigger the
# positional misalignment, which needs a variable that appears only inside an
# interaction. It is kept as a cheap pin on the ordinary case, and labelled
# here so nobody reads the trio as three regression guards (#450 review).

vc_data <- function() {
  set.seed(31)
  data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5),
             mal = rbinom(60, 1, 0.4), w = rnorm(60))
}

vc_fit <- function(rhs, theta, d) {
  zz <- rnorm(nrow(d))
  f <- stats::as.formula(paste("survival::Surv(t, s) ~", rhs))
  environment(f) <- environment()
  suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = theta, fit = TRUE)
  )
}

vc_named <- "does not give one value per row of 'newdata'"

test_that("an outside term is named whatever ORDER it is written in", {
  d <- vc_data()
  nd <- d[1:2, c("age", "mal", "w")]
  shapes <- list(
    list("zz + mal", c(0.5, 1, 0, 0)),
    list("mal + zz", c(0.5, 1, 0, 0)),
    list("mal + zz + w", c(0.5, 1, 0, 0, 0))
  )
  for (sh in shapes) {
    expect_error(
      predict(vc_fit(sh[[1]], sh[[2]], d), newdata = nd,
              type = "linear_predictor"),
      vc_named, fixed = TRUE, info = sh[[1]]
    )
  }
})

test_that("a variable used only inside an interaction is still named", {
  d <- vc_data()
  nd <- d[1:2, c("age", "mal", "w")]
  shapes <- list(
    list("mal:age + zz", c(0.5, 1, 0, 0)),
    list("zz + mal:age", c(0.5, 1, 0, 0)),
    list("mal*age + zz", c(0.5, 1, 0, 0, 0, 0))
  )
  for (sh in shapes) {
    expect_error(
      predict(vc_fit(sh[[1]], sh[[2]], d), newdata = nd,
              type = "linear_predictor"),
      vc_named, fixed = TRUE, info = sh[[1]]
    )
  }
})

test_that("a length-changing function of a data column is still named", {
  # These have NO outside symbol, so an experiment that corrected only
  # outside variables would leave them failing and hand back base R's raw
  # error -- losing #409's message for a case #430 tests. Their own value is
  # corrected like any other, so no special case is needed.
  set.seed(11)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = seq(40, 99),
                  mal = rbinom(60, 1, 0.4))
  dup <- data.frame(age = c(30, 30, 50), mal = c(0, 1, 0))
  shapes <- list(
    list("I(unique(age))", c(0.5, 1, 0)),
    list("I(unique(age)) + zz", c(0.5, 1, 0, 0)),
    list("zz + I(unique(age))", c(0.5, 1, 0, 0)),
    list("mal:age + I(unique(age))", c(0.5, 1, 0, 0))
  )
  for (sh in shapes) {
    expect_error(
      predict(vc_fit(sh[[1]], sh[[2]], d), newdata = dup,
              type = "linear_predictor"),
      vc_named, fixed = TRUE, info = sh[[1]]
    )
  }
})

test_that("a newdata column named like the substitution symbol cannot shadow it", {
  # `model.frame()` resolves `predvars` against `data` BEFORE the terms'
  # environment, so a column named like the symbol the causation experiment
  # binds would shadow it: the experiment would run, report success and have
  # tested nothing -- this package's signature shape (#450 review). The
  # symbol base is lengthened until no `newdata` column starts with it.
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  gg <- function(x) if (length(x) == 2L) as.list(rnorm(60)) else x
  for (nm in c("ordinary", "..hzr_v1")) {
    dd <- d
    dd[[nm]] <- rnorm(60)
    f <- stats::as.formula(
      paste0("survival::Surv(t, s) ~ I(gg(age)) + `", nm, "`")
    )
    environment(f) <- environment()
    fit <- suppressWarnings(
      hazard(f, data = dd, dist = "weibull", theta = c(0.5, 1, 0, 0),
             fit = TRUE)
    )
    msg <- tryCatch(
      predict(fit, newdata = dd[1:2, c("age", nm)], type = "linear_predictor"),
      error = conditionMessage
    )
    # The list-valued term is the real cause either way, so the shadowing
    # column must not change the answer.
    expect_match(msg, "invalid type (list) for variable", fixed = TRUE,
                 info = nm)
  }
})

test_that("a failure that does not reproduce may be misattributed (documented limit)", {
  # PINS A KNOWN LIMITATION, deliberately, so the code and the note in #446
  # cannot drift apart. The diagnosis must evaluate every variable to learn
  # its row count, and that evaluation is itself a second run. A term that
  # fails ONCE has therefore already succeeded by the time any causation
  # test could run, and the evidence is gone: no experiment can tell "failed
  # because of length" from "failed for a reason that went away".
  #
  # Measured step by step: build calls ff (fails), diagnosis calls ff again
  # (succeeds, one-shot), so `I(ff(age))` is no longer row-mismatched and
  # only `zz` is -- which is then named. Loud and rare; recorded in #446
  # rather than chased.
  set.seed(31)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7), age = rnorm(60, 60, 5))
  zz <- rnorm(60)
  boom <- FALSE
  ff <- function(x) {
    if (boom) {
      boom <<- FALSE
      stop("formula boom")
    }
    x
  }
  f <- survival::Surv(t, s) ~ I(ff(age)) + zz
  environment(f) <- environment()
  fit <- suppressWarnings(
    hazard(f, data = d, dist = "weibull", theta = c(mu = 0.5, nu = 1, 0, 0),
           fit = TRUE)
  )
  boom <- TRUE
  msg <- tryCatch(
    predict(fit, newdata = d[1:2, "age", drop = FALSE],
            type = "linear_predictor"),
    error = conditionMessage
  )
  # The limitation, stated as an expectation. If a future change makes the
  # caller's "formula boom" survive here, THIS TEST FAILS and the note in
  # #446 and the roxygen must be updated with it.
  expect_match(msg, "term 'zz' of the model does not give one value per row",
               fixed = TRUE)
})
