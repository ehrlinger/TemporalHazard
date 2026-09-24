# A column of `data` named "" stopped every formula-interface fit with
# "attempt to use zero-length variable name", even when the formula never
# named it (#470). No formula can name such a column, so it can play no part
# in the model: the fit must be the fit without it.
#
# Two base-R calls in the formula path fail on a "" name: `list2env()` over
# the data, and `terms(formula, data = data)` for any two-sided formula. The
# `time =`/`status =` interface never failed, because it uses neither.
#
# The column is made with `names(d)[j] <- ""`, never `d[[""]] <- x`, which
# silently names it "V5" and tests nothing. Each test asserts the name first.

ecn_data <- function() {
  set.seed(9)
  d <- data.frame(t = rexp(40) + 0.1, s = rbinom(40, 1, 0.8),
                  age = rnorm(40), mal = rbinom(40, 1, 0.4), x = rnorm(40))
  names(d)[names(d) == "x"] <- ""
  d
}

ecn_fit <- function(formula, d, theta) {
  hazard(formula, data = d, dist = "weibull", theta = theta, fit = TRUE)
}

ecn_warnings <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(msgs = msgs, value = val)
}

test_that("an unused column named \"\" does not change the fit", {
  d <- ecn_data()
  expect_true(any(names(d) == ""))
  clean <- d[, names(d) != ""]
  th <- c(mu = 0.5, nu = 1, 0, 0)
  r <- ecn_warnings(ecn_fit(survival::Surv(t, s) ~ age + mal, d, th))
  # Named columns only: the "" column is irrelevant here, so nothing to say.
  expect_length(r$msgs, 0L)
  a <- r$value
  b <- ecn_fit(survival::Surv(t, s) ~ age + mal, clean, th)
  expect_identical(a$fit$theta, b$fit$theta)
  expect_identical(a$fit$objective, b$fit$objective)
})

test_that("stepwise refits are not broken by a column named \"\"", {
  # The #470 reproducer: the base model is fitted on clean data, and only
  # the data handed to hzr_stepwise() carries the column. Before the fix
  # every candidate refit failed and selection stopped after 0 steps.
  d <- ecn_data()
  expect_true(any(names(d) == ""))
  clean <- d[, names(d) != ""]
  base <- ecn_fit(survival::Surv(t, s) ~ age, clean, c(mu = 0.5, nu = 1, 0))
  run <- function(data) {
    msgs <- character(0)
    sw <- withCallingHandlers(
      hzr_stepwise(base, scope = ~ age + mal, data = data,
                   direction = "forward", criterion = "wald",
                   slentry = 0.9, trace = FALSE),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(msgs = msgs, sw = sw)
  }
  with_col <- run(d)
  without  <- run(clean)
  # The control must take a step, or "same model" would compare two empty
  # screens and prove nothing. Coefficient names in `theta` are blank, so
  # count parameters: the base `~ age` has 3, and a step adds one.
  expect_gt(length(without$sw$fit$theta), length(base$fit$theta))
  expect_identical(with_col$msgs, without$msgs)
  expect_identical(with_col$sw$fit$theta, without$sw$fit$theta)
  expect_identical(with_col$sw$fit$objective, without$sw$fit$objective)
})

test_that("`.` leaves a column named \"\" out, and says so", {
  d <- ecn_data()
  expect_true(any(names(d) == ""))
  clean <- d[, names(d) != ""]
  th <- c(mu = 0.5, nu = 1, 0, 0)
  r <- ecn_warnings(ecn_fit(survival::Surv(t, s) ~ ., d, th))
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "1 column of `data` named \"\" is left out of `.`",
               fixed = TRUE)
  expect_match(r$msgs[[1L]], "Rename it to include it.", fixed = TRUE)
  b <- ecn_fit(survival::Surv(t, s) ~ ., clean, th)
  expect_identical(r$value$fit$theta, b$fit$theta)
  expect_identical(r$value$fit$objective, b$fit$objective)

  # Two such columns: one warning, plural. Built by renaming, because
  # cbind() and `$<-` rebuild the frame and quietly repair an existing ""
  # name, leaving only one.
  d2 <- data.frame(clean, x = 1, y = 2)
  names(d2)[names(d2) %in% c("x", "y")] <- ""
  expect_identical(sum(names(d2) == ""), 2L)
  r2 <- ecn_warnings(ecn_fit(survival::Surv(t, s) ~ ., d2, th))
  expect_length(r2$msgs, 1L)
  expect_match(r2$msgs[[1L]], "2 columns of `data` named \"\" are left out",
               fixed = TRUE)
})
