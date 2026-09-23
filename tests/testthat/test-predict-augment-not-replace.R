# predict(newdata = ) names a term that does not give one value per row, but
# it never verified that the named term EXPLAINS the failure, so it replaced
# the real error with row blame: a list-returning term failed, and the message
# sent the user to repair an innocent `zz` (#446). The naming is now APPENDED
# to the original condition rather than substituted for it, so both statements
# are independently true and no causal claim is made.

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

test_that("a term that fails on its type keeps its own error, and zz is not blamed for it", {
  o <- an_setup()
  zz <- o$zz
  # Returns a list only at predict time, so the fit itself is clean.
  ff <- function(x) if (length(x) < 5) as.list(x) else x
  fit <- an_fit(survival::Surv(t, s) ~ I(ff(age)) + zz, o$d, 2)
  msg <- tryCatch(
    predict(fit, newdata = data.frame(age = c(0, 1), mal = c(0, 1)),
            type = "linear_predictor"),
    error = conditionMessage
  )
  # The real failure, verbatim and FIRST: this is what the user must fix.
  expect_equal(strsplit(msg, "\n", fixed = TRUE)[[1L]][[1L]],
               "invalid type (list) for variable 'I(ff(age))'")
  # The naming survives too, but as a separate statement that disclaims cause.
  expect_match(msg, "Separately: term 'zz' of the model does not give one value",
               fixed = TRUE)
  expect_match(msg, "whether or not it caused the error above", fixed = TRUE)
})

test_that("the caller's condition class and call survive the augmentation", {
  o <- an_setup()
  zz <- o$zz
  # A classed condition carrying a model.frame call: indistinguishable from a
  # design failure by conditionCall(), which is exactly the case that must not
  # lose its class (#446, second shape).
  boom <- function(x) {
    if (length(x) >= 5) return(x)
    stop(structure(
      class = c("ff_boom", "error", "condition"),
      list(message = "ff_boom fired", call = quote(model.frame.default(x)))
    ))
  }
  fit <- an_fit(survival::Surv(t, s) ~ I(boom(age)) + zz, o$d, 2)
  caught <- tryCatch(
    predict(fit, newdata = data.frame(age = c(0, 1), mal = c(0, 1)),
            type = "linear_predictor"),
    ff_boom = function(e) e, error = function(e) e
  )
  expect_s3_class(caught, "ff_boom")
  expect_equal(class(caught), c("ff_boom", "error", "condition"))
  expect_equal(conditionCall(caught), quote(model.frame.default(x)))
  expect_equal(strsplit(conditionMessage(caught), "\n", fixed = TRUE)[[1L]][[1L]],
               "ff_boom fired")
})

test_that("our own row-count backstop is not duplicated", {
  # The backstop already says what the note says, so it is replaced, not
  # appended to, and the message carries exactly one statement.
  #
  # `~ zz` ALONE is the shape that reaches it, measured rather than assumed:
  # with a `data` column in the formula too, `model.frame()` raises
  # "variable lengths differ" first and the backstop is never reached.
  o <- an_setup()
  zz <- o$zz
  fit <- an_fit(survival::Surv(t, s) ~ zz, o$d, 1)
  msg <- tryCatch(
    predict(fit, newdata = data.frame(mal = c(0, 1)),
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "term 'zz' of the model does not give one value per row",
               fixed = TRUE)
  expect_false(grepl("Separately:", msg, fixed = TRUE))
  expect_false(grepl("The design rebuilt for", msg, fixed = TRUE))
})

test_that("an xlev failure keeps its own message and still names the term", {
  # The name is NOT lost when the original failure is about factor levels:
  # gg is genuinely outside `data` and row-mismatched, and saying so is true
  # whichever error model.frame() happened to raise first.
  o <- an_setup()
  d <- o$d
  gg <- factor(sample(c("a", "b"), 60, TRUE))
  fit <- an_fit(survival::Surv(t, s) ~ gg, d, 2)
  gg <- factor(sample(c("Y", "Z"), 60, TRUE))
  msg <- tryCatch(
    predict(fit, newdata = data.frame(mal = c(0, 1)), type = "linear_predictor"),
    error = conditionMessage
  )
  expect_match(msg, "new levels", fixed = TRUE)
  expect_match(msg, "Separately: term 'gg' of the model does not give one value",
               fixed = TRUE)
})

test_that("a user error raised outside model.frame() is still untouched", {
  # Unchanged behaviour, pinned so the augmentation cannot start reaching a
  # condition the caller raised in their own call.
  o <- an_setup()
  zz <- o$zz
  n <- 0
  ff <- function(x) {
    n <<- n + 1
    if (n == 1) stop("formula boom")
    x
  }
  n <- 99
  fit <- an_fit(survival::Surv(t, s) ~ I(ff(age)) + zz, o$d, 2)
  n <- 0
  msg <- tryCatch(
    predict(fit, newdata = data.frame(age = c(0, 1), mal = c(0, 1)),
            type = "linear_predictor"),
    error = conditionMessage
  )
  expect_equal(msg, "formula boom")
})
