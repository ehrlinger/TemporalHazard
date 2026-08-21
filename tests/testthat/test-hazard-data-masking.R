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
  fit <- hazard(data = df, time = tt, status = ev, dist = "weibull",
                theta = c(0.5, 1), fit = TRUE)
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
