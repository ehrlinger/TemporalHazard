# `predict(newdata = )` replaced a failed design build with an error naming
# whichever term did not give one value per row of `newdata`. The naming was
# never verified to explain the failure, so it blamed innocent terms and
# destroyed conditions the caller raised inside their own terms (#446).
#
# It now decides by WHO raised the failure. A condition carrying any class of
# the caller's own passes through UNCHANGED -- object, class, call and every
# field -- and gains no note. A plain base error from the frame build gets our
# note FIRST, with the base text after it. And a term whose value is not a
# legal model-frame column is not diagnosed as a row-count problem at all.

an_setup <- function() {
  set.seed(11)
  d <- data.frame(t = rexp(60), s = rbinom(60, 1, 0.7),
                  mal = rbinom(60, 1, 0.4), age = rnorm(60))
  list(d = d, zz = rnorm(60))
}

an_fit <- function(f, d, k) {
  hazard(f, data = d, dist = "weibull",
         theta = c(mu = 0.5, nu = 1, rep(0, k)), fit = TRUE)
}

an_msg <- function(fit, nd) {
  tryCatch(predict(fit, newdata = nd, type = "linear_predictor"),
           error = conditionMessage)
}

an_first <- function(msg) strsplit(msg, "\n", fixed = TRUE)[[1L]][[1L]]

test_that("the first line names the term, not a column the user supplied", {
  # THE COMMONEST SHAPE. `mal` is supplied correctly; `zz` is the outside
  # vector. Left to model.frame() the failure reads "variable lengths differ
  # (found for 'mal')", which sends the user to repair `mal`. Our statement
  # leads and the base text follows it.
  o <- an_setup()
  zz <- o$zz
  fit <- an_fit(survival::Surv(t, s) ~ zz + mal, o$d, 2)
  msg <- an_msg(fit, data.frame(mal = c(0, 1)))
  expect_match(an_first(msg), "term 'zz' of the model does not give one value",
               fixed = TRUE)
  expect_false(grepl("variable lengths differ", an_first(msg), fixed = TRUE))
  expect_match(msg, "The design build reported: variable lengths differ",
               fixed = TRUE)
})

test_that("a term whose value is not a legal model-frame column keeps its own error", {
  # CASE D. A list-valued term fails on its TYPE. `zz` is innocent, and
  # "move zz into data" could not be followed -- zz is not a column of `data`
  # for `newdata` to supply. The row-count diagnosis declines, structurally,
  # because a variable evaluated to a type no model frame can hold.
  o <- an_setup()
  zz <- o$zz
  ff <- function(x) if (length(x) < 5) as.list(x) else x
  fit <- an_fit(survival::Surv(t, s) ~ I(ff(age)) + zz, o$d, 2)
  msg <- an_msg(fit, data.frame(age = c(0, 1), mal = c(0, 1)))
  expect_equal(msg, "invalid type (list) for variable 'I(ff(age))'")
  # DOCUMENTED CONSEQUENCE: no note is added. The user is told the real
  # failure and nothing else, because nothing else is known to be relevant.
  expect_false(grepl("does not give one value per row", msg, fixed = TRUE))
})

test_that("a caller's condition passes through unchanged, with class, call and fields", {
  # CASE J. Every field survives, so a tryCatch() handler reading $culprit
  # gets the value rather than NULL -- firing on missing data is worse than
  # not firing.
  o <- an_setup()
  zz <- o$zz
  boom <- function(x) {
    if (length(x) >= 5) return(x)
    stop(structure(
      class = c("my_err", "error", "condition"),
      list(message = "payload error", call = quote(model.frame.default(x)),
           culprit = "age", data = 42)
    ))
  }
  fit <- an_fit(survival::Surv(t, s) ~ I(boom(age)) + zz, o$d, 2)
  caught <- tryCatch(
    predict(fit, newdata = data.frame(age = c(0, 1), mal = c(0, 1)),
            type = "linear_predictor"),
    my_err = function(e) e, error = function(e) e
  )
  expect_equal(class(caught), c("my_err", "error", "condition"))
  expect_equal(conditionMessage(caught), "payload error")
  expect_equal(conditionCall(caught), quote(model.frame.default(x)))
  expect_equal(caught$culprit, "age")
  expect_equal(caught$data, 42)
  # DOCUMENTED CONSEQUENCE: no note is appended to a caller's condition.
  expect_false(grepl("does not give one value per row",
                     conditionMessage(caught), fixed = TRUE))
})

test_that("a caller's own conditionMessage method is not defeated", {
  # CASE I. `conditionMessage()` is a GENERIC. Writing our text into
  # `$message` does not mean the user reads it: a class with its own method
  # never looks there. Appending destroyed BOTH statements -- the note was
  # unreachable and the caller's own detail was dropped.
  o <- an_setup()
  zz <- o$zz
  conditionMessage.detail_err <- function(c) paste0("detail: ", c$detail)
  registerS3method("conditionMessage", "detail_err",
                   conditionMessage.detail_err)
  ff <- function(x) {
    if (length(x) >= 5) return(x)
    stop(structure(
      class = c("detail_err", "error", "condition"),
      list(message = "fallback", call = quote(model.frame.default(x)),
           detail = "column 'age' had a list")
    ))
  }
  fit <- an_fit(survival::Surv(t, s) ~ I(ff(age)) + zz, o$d, 2)
  caught <- tryCatch(
    predict(fit, newdata = data.frame(age = c(0, 1), mal = c(0, 1)),
            type = "linear_predictor"),
    detail_err = function(e) e, error = function(e) e
  )
  expect_s3_class(caught, "detail_err")
  expect_equal(conditionMessage(caught), "detail: column 'age' had a list")
})

test_that("a multi-element condition message is not multiplied", {
  # CASE K. `paste0()` vectorises, so appending to a length-2 message emitted
  # a full copy of the note against each element. Passing the object through
  # cannot do that.
  o <- an_setup()
  zz <- o$zz
  ff <- function(x) {
    if (length(x) >= 5) return(x)
    stop(structure(
      class = c("m2", "error", "condition"),
      list(message = c("first half", "second half"),
           call = quote(model.frame.default(x)))
    ))
  }
  fit <- an_fit(survival::Surv(t, s) ~ I(ff(age)) + zz, o$d, 2)
  caught <- tryCatch(
    predict(fit, newdata = data.frame(age = c(0, 1), mal = c(0, 1)),
            type = "linear_predictor"),
    m2 = function(e) e, error = function(e) e
  )
  expect_equal(conditionMessage(caught), c("first half", "second half"))
})

test_that("our own errors carry no call, so no internal names reach the user", {
  # UNGUARDED until this test existed: flipping `call. = FALSE` to TRUE killed
  # nothing, and the mutant emits
  #   Error in .hzr_stop_unmatched_rows(term, where, nrow(newdata), ...)
  # at a CRAN user. Both routes that raise an error of OUR OWN are pinned.
  o <- an_setup()
  zz <- o$zz
  base_route <- an_fit(survival::Surv(t, s) ~ zz + mal, o$d, 2)
  backstop_route <- an_fit(survival::Surv(t, s) ~ zz, o$d, 1)
  for (case in list(list(base_route, data.frame(mal = c(0, 1))),
                    list(backstop_route, data.frame(mal = c(0, 1))))) {
    e <- tryCatch(
      predict(case[[1L]], newdata = case[[2L]], type = "linear_predictor"),
      error = function(e) e
    )
    expect_null(conditionCall(e))
  }
})

test_that("our own row-count backstop is not duplicated", {
  # `~ zz` ALONE reaches the classed backstop, measured rather than assumed:
  # with a `data` column in the formula too, model.frame() raises "variable
  # lengths differ" first and the backstop is never reached.
  #
  # All three assertions pass on `main` as well, so this block pins nothing
  # against main -- it kills the mutant that drops the backstop exception,
  # which was verified. Do not delete it as dead weight.
  o <- an_setup()
  zz <- o$zz
  fit <- an_fit(survival::Surv(t, s) ~ zz, o$d, 1)
  msg <- an_msg(fit, data.frame(mal = c(0, 1)))
  expect_match(msg, "term 'zz' of the model does not give one value per row",
               fixed = TRUE)
  expect_false(grepl("The design build reported", msg, fixed = TRUE))
  expect_false(grepl("The design rebuilt for", msg, fixed = TRUE))
})

test_that("an xlev failure keeps the name and leads with it", {
  # The name is NOT lost when the original failure is about factor levels:
  # `gg` is genuinely outside `data` and row-mismatched, and saying so is true
  # whichever error model.frame() happened to raise first.
  o <- an_setup()
  d <- o$d
  gg <- factor(sample(c("a", "b"), 60, TRUE))
  fit <- an_fit(survival::Surv(t, s) ~ gg, d, 2)
  gg <- factor(sample(c("Y", "Z"), 60, TRUE))
  msg <- an_msg(fit, data.frame(mal = c(0, 1)))
  expect_match(an_first(msg), "term 'gg' of the model does not give one value",
               fixed = TRUE)
  expect_match(msg, "The design build reported: factor gg has new levels",
               fixed = TRUE)
})

test_that("a user error raised outside model.frame() is still untouched", {
  o <- an_setup()
  zz <- o$zz
  n <- 99
  ff <- function(x) {
    n <<- n + 1
    if (n == 1) stop("formula boom")
    x
  }
  fit <- an_fit(survival::Surv(t, s) ~ I(ff(age)) + zz, o$d, 2)
  n <- 0
  msg <- an_msg(fit, data.frame(age = c(0, 1), mal = c(0, 1)))
  expect_equal(msg, "formula boom")
})
