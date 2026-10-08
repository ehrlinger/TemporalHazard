# Bootstrap replicates that passed as successes over work that was not done.
#
# (B) A Weibull fit with a covariate far from zero keeps log(mu) in
# $fit$log_scale because mu itself is 0, Inf or subnormal (#566). The
# single-distribution stepwise warm start copied natural mu, which
# .hzr_check_theta() refuses, so every candidate refit from such a fit failed.
# hzr_stepwise() said so; a select-mode hzr_bootstrap() read none of it and
# pooled each replicate as a success with its candidates "not selected".
#
# (C) hazard() warns, unclassed, when a fit stops on nlm() code 4 or 5 and
# fails the relative-gradient test. Replicates run under suppressWarnings(),
# and the #531 count leaves codes 4 and 5 out, so such replicates were pooled
# with no count and no warning.
#
# The oracle for (B) is the same model with the covariate centred: mu is
# representable there, and the screen is the same screen.

.bs_data <- function() {
  set.seed(1)
  n <- 400
  x <- 1000 + stats::rnorm(n, sd = 5)
  z <- stats::rnorm(n)
  noise <- stats::rnorm(n)
  t <- (stats::rexp(n) / exp(-1.4 + 0.14 * (x - 1000) + 0.8 * z))^(1 / 0.2)
  cens <- stats::rexp(n)^5 * 3
  data.frame(time = pmin(t, cens), dead = as.integer(t <= cens),
             x = x, xc = x - 1000, z = z, noise = noise)
}

# The centred fit, and the uncentred fit started at its maximum.
.bs_fits <- function(d) {
  fc <- suppressWarnings(hazard(survival::Surv(time, dead) ~ xc, data = d,
                                dist = "weibull", theta = c(1, 0.2, 0.1),
                                fit = TRUE))
  th <- unname(coef(fc))
  log_mu0 <- log(th[1]) - 1000 * th[3] / th[2]
  f <- suppressWarnings(hazard(survival::Surv(time, dead) ~ x, data = d,
                               dist = "weibull",
                               theta = c(exp(log_mu0), th[2], th[3]),
                               fit = TRUE))
  list(centred = fc, raw = f)
}

.bs_collect <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, msgs = msgs)
}

test_that("the refit's Weibull start reproduces the base fit's linear predictor (#566)", {
  d <- .bs_data()
  # A base whose mu is exactly 0: log(mu) is far below the double range.
  sw <- suppressWarnings(hzr_stepwise(
    .bs_fits(d)$raw, data = d, scope = ~ z, criterion = "wald",
    direction = "forward", trace = FALSE))
  expect_identical(unname(coef(sw)[[1L]]), 0)          # the premise
  expect_true(isTRUE(sw$fit$log_scale$needed))
  log_mu <- sw$fit$log_scale$theta[[1L]]
  expect_lt(log_mu, log(.Machine$double.xmin) - 10)

  x_new <- as.matrix(d[, c("x", "z", "noise")])
  start <- .hzr_refit_weibull_start(c(unname(sw$fit$theta), 0), sw, x_new,
                                    n_win = 1L)
  expect_true(start[[1L]] >= .Machine$double.xmin)    # representable
  expect_silent(.hzr_check_theta(start, "weibull", n_coef = 3L))
  # The log cumulative hazard's offset, nu * log(mu) + x'beta, per row: the
  # start must carry the base fit's, not merely be accepted.
  nu <- sw$fit$theta[[2L]]
  lp_base <- nu * log_mu + drop(x_new[, 1:2] %*% sw$fit$theta[3:4])
  lp_start <- nu * log(start[[1L]]) + drop(x_new %*% start[3:5])
  # Moving log(mu) alone would leave every row off by the whole shift; the
  # coefficients must carry nearly all of it.
  shift <- nu * (log(start[[1L]]) - log_mu)
  expect_gt(shift, 1)
  expect_lt(max(abs(lp_start - lp_base)), shift / 20)

  # A representable mu is left alone.
  fc <- .bs_fits(d)$centred
  expect_identical(.hzr_refit_weibull_start(c(fc$fit$theta, 0), fc,
                                            as.matrix(d[, c("xc", "z")]), 1L),
                   c(fc$fit$theta, 0))
})

test_that("hzr_stepwise() on a fit whose mu reads 0 tests its candidates (#566)", {
  d <- .bs_data()
  fits <- .bs_fits(d)
  expect_true(isTRUE(fits$raw$fit$log_scale$needed))   # the premise
  sw <- suppressWarnings(hzr_stepwise(
    fits$raw, data = d, scope = ~ z + noise, criterion = "wald",
    direction = "forward", trace = FALSE))
  swc <- suppressWarnings(hzr_stepwise(
    fits$centred, data = d, scope = ~ z + noise, criterion = "wald",
    direction = "forward", trace = FALSE))
  # After z enters, mu is exactly 0, so noise's refit starts from such a fit.
  expect_identical(unname(coef(sw)[[1L]]), 0)
  expect_identical(sw$criteria$refit_failures, character(0))
  expect_identical(sw$steps$variable, swc$steps$variable)
  expect_equal(sw$steps$p_value, swc$steps$p_value, tolerance = 1e-4)
  expect_equal(sw$fit$objective, swc$fit$objective, tolerance = 1e-6)
})

test_that("a select-mode bootstrap on such a fit selects as the centred fit does (#566)", {
  skip_on_cran()
  d <- .bs_data()
  fits <- .bs_fits(d)
  b <- .bs_collect(hzr_bootstrap(
    fits$raw, n_boot = 20, seed = 2, scope = ~ x + z + noise,
    criterion = "wald", direction = "forward"))
  bc <- .bs_collect(hzr_bootstrap(
    fits$centred, n_boot = 20, seed = 2, scope = ~ xc + z + noise,
    criterion = "wald", direction = "forward"))
  # The premise: replicates did report a mu that cannot be represented.
  expect_length(grep("reported a Weibull scale mu", b$msgs, fixed = TRUE), 1L)
  expect_identical(b$value$n_success, 20L)
  expect_identical(bc$value$n_success, 20L)
  # Selections, replicate for replicate.
  sel <- function(r) {
    vapply(split(r$parameter, r$replicate), function(p) {
      paste(sort(intersect(p, c("z", "noise"))), collapse = "+")
    }, character(1))
  }
  s_raw <- sel(b$value$replicates)
  s_ctr <- sel(bc$value$replicates)
  expect_length(s_raw, 20L)
  expect_identical(s_raw, s_ctr)
  # The selections varied, so agreement is not two constant columns.
  expect_gt(length(unique(s_raw)), 1L)
  z_est <- b$value$replicates$estimate[b$value$replicates$parameter == "z"]
  expect_gt(stats::sd(z_est), 0)
  expect_false(any(grepl("candidate refit fail", b$msgs, fixed = TRUE)))
})

test_that("a select-mode bootstrap counts replicates whose candidate refits failed", {
  set.seed(11)
  n <- 150
  z <- stats::rnorm(n)
  t <- stats::rexp(n, 0.3 * exp(0.9 * z))
  df <- data.frame(time = pmin(t, 6), status = as.integer(t <= 6), z = z)
  base <- hazard(survival::Surv(time, status) ~ 1, data = df,
                 dist = "exponential", theta = c(log_rate = 0), fit = TRUE)
  calls <- new.env()
  calls$n <- 0L
  testthat::local_mocked_bindings(.hzr_refit_with_scope = function(...) {
    calls$n <- calls$n + 1L
    stop("planted refit failure")
  })
  b <- .bs_collect(hzr_bootstrap(base, n_boot = 3L, seed = 2L, scope = ~ z,
                                 direction = "forward", criterion = "wald"))
  expect_gt(calls$n, 0L)                               # the premise
  n_ok <- b$value$n_success
  expect_identical(n_ok, 3L)
  hit <- grep("had a candidate refit fail", b$msgs, fixed = TRUE,
              value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, paste0("^", n_ok, " of ", n_ok, " successful replicates"))
  expect_match(hit, paste0("(", n_ok, " stopped on it)"), fixed = TRUE)
})

.bs_exp_data <- function() {
  set.seed(7)
  data.frame(time = stats::rexp(80, 0.4),
             status = rep(c(1, 1, 0), length.out = 80),
             z = stats::rnorm(80))
}

.bs_exp_fit <- function(df) {
  hazard(survival::Surv(time, status) ~ z, data = df, dist = "exponential",
         theta = c(log_rate = 0, z = 0), fit = TRUE)
}

.bs_stuck_code <- function(code, only_intercept = FALSE) {
  real <- .hzr_optim_exponential
  function(...) {
    r <- real(...)
    if (!only_intercept || length(r$par) == 1L) {
      r$polish_code <- code
      r$rel_gradient <- 1e-2
      r$convergence <- 0L
    }
    r
  }
}

test_that("hzr_bootstrap() counts replicates stopped on nlm() code 4 with a failed gradient test", {
  df <- .bs_exp_data()
  base <- .bs_exp_fit(df)
  testthat::local_mocked_bindings(.hzr_optim_exponential = .bs_stuck_code(4L))
  # The premise: a direct fit under the mock warns, unclassed, and stays
  # converged, so a replicate would pool it.
  direct <- .bs_collect(.bs_exp_fit(df))
  expect_true(direct$value$fit$converged)
  expect_identical(direct$value$fit$polish_code, 4L)
  expect_length(grep("fail the relative-gradient test", direct$msgs,
                     fixed = TRUE), 1L)
  b <- .bs_collect(hzr_bootstrap(base, n_boot = 4L, seed = 3L))
  n_ok <- b$value$n_success
  expect_identical(n_ok, 4L)
  expect_identical(b$value$n_failed, 0L)
  hit <- grep("stopped on nlm() code 4 or 5", b$msgs, fixed = TRUE,
              value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, paste0("^", n_ok, " of ", n_ok, " successful replicates"))
  # #531's count leaves codes 4 and 5 out, so it must stay silent here.
  expect_false(any(grepl("may not be a maximum, in the base fit",
                         b$msgs, fixed = TRUE)))
})

test_that("a select-mode bootstrap counts a code-5 stop in the BASE refit", {
  set.seed(11)
  n <- 150
  z <- stats::rnorm(n)
  t <- stats::rexp(n, 0.3 * exp(0.9 * z))
  df <- data.frame(time = pmin(t, 6), status = as.integer(t <= 6), z = z)
  base <- hazard(survival::Surv(time, status) ~ 1, data = df,
                 dist = "exponential", theta = c(log_rate = 0), fit = TRUE)
  testthat::local_mocked_bindings(
    .hzr_optim_exponential = .bs_stuck_code(5L, only_intercept = TRUE))
  b <- .bs_collect(hzr_bootstrap(base, n_boot = 3L, seed = 2L, scope = ~ z,
                                 direction = "forward", criterion = "wald"))
  n_ok <- b$value$n_success
  expect_identical(n_ok, 3L)
  # The premise: the final fits carry z and are clean, so only the base
  # refit can be what is counted.
  expect_true(all(b$value$summary$parameter %in% c("log_rate", "z")))
  expect_equal(b$value$summary$pct[b$value$summary$parameter == "z"], 100)
  hit <- grep("stopped on nlm() code 4 or 5", b$msgs, fixed = TRUE,
              value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, paste0("^", n_ok, " of ", n_ok, " successful replicates"))
})
