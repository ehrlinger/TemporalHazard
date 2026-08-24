# hazard()'s `data` argument was accepted and silently ignored on the vector
# path: only the formula path consulted it. That is why the SAS translator's
# emitted hazard(data = AVCS, time = INT_DEAD, ...) failed with
# "object 'INT_DEAD' not found" (#151).

make_df <- function(n = 100) {
  set.seed(4)
  data.frame(tt = stats::rexp(n, 0.3),
             ev = rep(c(1, 0), length.out = n),
             entry = rep(0, n))
}

test_that("bare column names resolve against `data` on the vector path", {
  df <- make_df()
  fit <- hazard(data = df, time = tt, status = ev, dist = "weibull",
                theta = c(0.5, 1), fit = TRUE)
  expect_s3_class(fit, "hazard")
  expect_equal(fit$data$time, df$tt)
  expect_equal(fit$data$status, df$ev)
})

test_that("df$col, a local vector and a bare name all still work", {
  df <- make_df()
  local_time <- df$tt
  a <- hazard(data = df, time = df$tt, status = df$ev, dist = "weibull",
              theta = c(0.5, 1), fit = TRUE)
  b <- hazard(data = df, time = local_time, status = ev, dist = "weibull",
              theta = c(0.5, 1), fit = TRUE)
  d <- hazard(time = df$tt, status = df$ev, dist = "weibull",
              theta = c(0.5, 1), fit = TRUE)
  expect_equal(a$fit$theta, d$fit$theta)
  expect_equal(b$fit$theta, d$fit$theta)
})

test_that("a column masks a same-named caller variable", {
  df <- make_df()
  tt <- rep(999, nrow(df))
  # `tt` is both a column and bound here, so the call is now audibly ambiguous.
  expect_warning(
    fit <- hazard(data = df, time = tt, status = ev, dist = "weibull",
                  theta = c(0.5, 1), fit = TRUE),
    "tt"
  )
  expect_equal(fit$data$time, df$tt)
  expect_false(any(fit$data$time == 999))
})

test_that("time_lower and weights are masked too", {
  df <- make_df()
  df$w <- rep(c(1, 2), length.out = nrow(df))
  fit <- hazard(data = df, time = tt, status = ev, time_lower = entry,
                weights = w, dist = "weibull", theta = c(0.5, 1), fit = TRUE)
  expect_equal(fit$data$weights, df$w)
})

test_that("the formula path is unaffected", {
  df <- make_df()
  fit <- hazard(survival::Surv(tt, ev) ~ 1, data = df, dist = "weibull",
                theta = c(0.5, 1), fit = TRUE)
  expect_s3_class(fit, "hazard")
  expect_equal(fit$data$time, df$tt)
})

# Masking made a wrapper that forwards its own argument by name silently pick
# up the column instead of the vector the caller passed -- a fit over the
# wrong rows, with no error and no warning. The column still wins; the
# ambiguity is now audible.

test_that("a wrapper forwarding its own argument warns about the shadowed name", {
  df <- make_df()
  sub <- df[1:20, ]
  wrap <- function(tt, ev, d) {
    hazard(data = d, time = tt, status = ev, dist = "weibull",
           theta = c(0.5, 1), fit = TRUE)
  }
  expect_warning(fit <- wrap(sub$tt, sub$ev, df), "tt")
  expect_equal(fit$data$time, df$tt)
})

test_that("an unambiguous masked call does not warn", {
  df <- make_df()
  expect_no_warning(
    hazard(data = df, time = tt, status = ev, dist = "weibull",
           theta = c(0.5, 1), fit = TRUE)
  )
})

# `exists(inherits = FALSE)` saw only the immediate frame, so the wrapper case
# above warned only when the shadowed name was an argument of that one
# function. The ordinary top-level shape -- tt <- ...; g <- function(d)
# hazard(data = d, time = tt, ...) -- resolves `tt` from the wrapper's
# enclosure, one link further out, and used to use the column in silence.
# The lexical chain is what eval() itself searches, so it is what the check
# must walk.

test_that("a name in the wrapper's enclosure, not its frame, is ambiguous", {
  df <- make_df()
  sub <- df[1:20, ]
  wrap <- local({
    tt <- sub$tt
    ev <- sub$ev
    function(d) {
      hazard(data = d, time = tt, status = ev, dist = "weibull",
             theta = c(0.5, 1), fit = TRUE)
    }
  })
  expect_warning(fit <- wrap(df), "tt")
  # The column won, as subset()'s rule requires -- the point is that it is
  # now audible.
  expect_equal(fit$data$time, df$tt)
})

# all.vars(quote(df$tt)) is c("df", "tt"), so the `$` form -- the remedy the
# warning itself prescribes -- used to trigger the warning, and named a column
# the fit had never consulted.

test_that("data$col does not warn, and neither does another frame's $ form", {
  df <- make_df()
  other <- data.frame(tt = df$tt * 1000, ev = df$ev)
  expect_no_warning(
    hazard(data = df, time = df$tt, status = df$ev, dist = "weibull",
           theta = c(0.5, 1), fit = TRUE)
  )
  # The column `tt` is not what this call uses; the warning must not claim it.
  expect_no_warning(
    fit <- hazard(data = df, time = other$tt, status = other$ev,
                  dist = "weibull", theta = c(0.5, 1), fit = TRUE)
  )
  expect_equal(fit$data$time, other$tt)
})

# A column named after a base object -- `c`, `t`, `df` -- is common in SAS
# exports. Walking the whole search path would warn on every such call.

test_that("a column named after a base function does not warn", {
  d <- data.frame(c = stats::rexp(40, 0.3), t = rep(c(1, 0), 20))
  f <- function(dd) {
    hazard(data = dd, time = c, status = t, dist = "weibull",
           theta = c(0.5, 1), fit = TRUE)
  }
  expect_no_warning(f(d))
})

test_that("a matrix passed as `data` errors rather than being ignored", {
  m <- as.matrix(make_df())
  expect_error(
    hazard(data = m, time = tt, status = ev, dist = "weibull", fit = FALSE),
    "must be a data frame or a list"
  )
})
