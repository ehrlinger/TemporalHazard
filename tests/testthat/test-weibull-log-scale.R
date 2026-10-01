# The Weibull likelihood on the log scale (#566, PR 2).
#
# .hzr_logl_weibull() formed mu^nu, (mu * t)^nu and t^(nu - 1) on the natural
# scale. Each can underflow or overflow where the log-likelihood is an ordinary
# number. The oracle in each case is the log-scale value, written out here.

.ls_eval <- function(obj, theta) {
  msgs <- character(0)
  ll <- withCallingHandlers(
    hzr_evaluate(obj, theta)$logLik,
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(ll = as.numeric(ll), msgs = msgs)
}

.ls_data <- function(alpha, b) {
  set.seed(1)
  n <- 400
  x <- 1000 + stats::rnorm(n, sd = 5)
  t <- (stats::rexp(n) / exp(alpha + b * x))^(1 / 0.2)
  cens <- stats::rexp(n)^5 * 3
  data.frame(time = pmin(t, cens), dead = as.integer(t <= cens),
             x = x, xc = x - 1000)
}

.ls_fit <- function(formula, data, theta) {
  suppressWarnings(hazard(formula, data = data, dist = "weibull",
                          theta = theta, fit = TRUE))
}

test_that("summary() shows the log(mu) a fit kept, with its standard error (#566)", {
  d <- .ls_data(200, -0.2)
  raw <- .ls_fit(survival::Surv(time, dead) ~ x, d, c(exp(50), 0.2, -0.2))
  ctr <- .ls_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.2))
  expect_identical(unname(coef(raw)[[1]]), Inf)   # premise
  tab <- summary(raw)$coefficients
  expect_identical(rownames(tab)[1:2], c(rownames(tab)[1], "log(mu)"))
  expect_equal(nrow(tab), length(coef(raw)) + 1L)
  # mu itself has no standard error; log(mu) does.
  expect_true(is.na(tab$std_error[1]))
  # Oracle: the centred fit. H = exp(nu (log mu + log t) + b x) with
  # x = xc + 1000, so log(mu_raw) = log(mu_ctr) - 1000 b / nu, and its
  # variance is the delta method on the centred fit's (log mu, nu, b).
  th <- unname(coef(ctr))
  lmu_c <- log(th[1])
  want <- lmu_c - 1000 * th[3] / th[2]
  expect_equal(tab$estimate[2] / want, 1, tolerance = 1e-4)
  v <- unname(vcov(ctr))
  dd <- c(1 / th[1], 1, 1)
  v_log <- v * outer(dd, dd)                       # (log mu, nu, b)
  grad <- c(1, 1000 * th[3] / th[2]^2, -1000 / th[2])
  se_want <- sqrt(as.numeric(grad %*% v_log %*% grad))
  expect_equal(tab$std_error[2] / se_want, 1, tolerance = 1e-2)
  expect_true(is.na(tab$z_stat[2]) && is.na(tab$p_value[2]))
  # print() shows the row.
  out <- utils::capture.output(print(summary(raw)))
  expect_true(any(grepl("^log\\(mu\\)", out)))

  # Known negative: a fit whose mu is representable gets no extra row.
  tab_ctr <- summary(ctr)$coefficients
  expect_equal(nrow(tab_ctr), length(coef(ctr)))
  expect_false("log(mu)" %in% rownames(tab_ctr))
})

test_that("a fit whose mu is representable predicts as before (#566)", {
  # Known negative for the log-scale readers: on an ordinary fit the
  # cumulative hazard and its SE are what the natural scale gives.
  d <- .ls_data(-2, 0.001)
  f <- .ls_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0))
  expect_null(f$fit$log_scale$needed)
  nd <- data.frame(time = c(0.5, 2), xc = c(-3, 4))
  th <- unname(coef(f))
  h_nat <- (th[1] * nd$time)^th[2] * exp(th[3] * nd$xc)
  got <- predict(f, newdata = nd, type = "cumulative_hazard", se.fit = TRUE)
  expect_equal(got$fit / h_nat, rep(1, 2), tolerance = 1e-12)
  # The SE, by the natural-scale Jacobian and vcov.
  j_nat <- cbind((th[2] / th[1]) * h_nat,
                 log(th[1] * nd$time) * h_nat,
                 nd$xc * h_nat)
  se_nat <- sqrt(rowSums((j_nat %*% unname(vcov(f))) * j_nat))
  expect_equal(got$se.fit / se_nat, rep(1, 2), tolerance = 1e-10)
})

test_that("a mu just below 1e-154 with a normal variance keeps its SEs (#566, r-reviewer)", {
  # mu^2 underflows here while var(mu) is still a normal double, so the fit
  # does not need its stored log scale, and the derived one must not form
  # 1 / mu^2 (which overflows). 1e8377fe returned NA SEs with a warning; main
  # e7b621c5 returned the centred fit's.
  set.seed(3)
  n <- 300
  x <- 1000 + stats::rnorm(n)
  lmu <- log(1e-155)
  b <- -lmu / 1000
  t <- stats::rexp(n) / exp(lmu + b * x)
  cens <- stats::rexp(n) * 2
  d <- data.frame(time = pmin(t, cens), dead = as.integer(t <= cens),
                  xc = x - 1000)
  fc <- .ls_fit(survival::Surv(time, dead) ~ xc, d, c(1, 1, 0))
  th <- unname(coef(fc))
  shift <- (log(th[1]) - lmu) / (th[3] / th[2])
  d$x2 <- d$xc + shift
  f <- .ls_fit(survival::Surv(time, dead) ~ x2, d, c(1e-155, th[2], th[3]))
  # Premise: mu in the band, its variance normal, no stored scale needed.
  mu <- unname(coef(f))[1]
  expect_true(mu < 1e-154 && mu > .Machine$double.xmin)
  expect_true(vcov(f)[1, 1] >= .Machine$double.xmin)
  expect_null(f$fit$log_scale$needed)
  expect_true(is.infinite(1 / mu^2))
  got <- expect_no_warning(predict(f, newdata = data.frame(
    time = c(0.5, 1), x2 = shift + 0:1), type = "survival", se.fit = TRUE))
  want <- predict(fc, newdata = data.frame(time = c(0.5, 1), xc = 0:1),
                  type = "survival", se.fit = TRUE)
  expect_equal(got$se.fit / want$se.fit, rep(1, 2), tolerance = 1e-3)
})

test_that("a left-censored row whose H is subnormal is read on the log scale (#566 acceptance 1)", {
  obj <- hazard(time = 2, status = -1, time_upper = 2, dist = "weibull",
                theta = c(1e-81, 4), fit = FALSE)
  # H = (2e-81)^4 = 1.6e-323 is subnormal, so log(1 - exp(-H)) = log(H) to
  # working precision. An oracle that forms H by exp() first loses the same
  # digits the defect does, so compare against log(H) directly.
  log_h <- 4 * (log(1e-81) + log(2))
  expect_lt(exp(log_h), .Machine$double.xmin)   # premise: H is subnormal
  got <- .ls_eval(obj, c(1e-81, 4))
  # main e7b621c5 returned -743.341, with no warning.
  expect_length(got$msgs, 0L)
  expect_equal(got$ll / log_h, 1, tolerance = 1e-12)
})

test_that("an interval whose cumulative hazards are subnormal is read on the log scale (#566)", {
  obj <- hazard(time = 1, status = 2, time_lower = 1, time_upper = 2,
                dist = "weibull", theta = c(1e-81, 4), fit = FALSE)
  # log(H(u) - H(l)) = log H(u) + log(1 - (l/u)^4), and -H(l) is negligible.
  truth <- 4 * (log(1e-81) + log(2)) + log(1 - (1 / 2)^4)
  expect_lt(exp(4 * (log(1e-81) + log(2))), .Machine$double.xmin)  # premise
  got <- .ls_eval(obj, c(1e-81, 4))
  expect_length(got$msgs, 0L)
  expect_equal(got$ll / truth, 1, tolerance = 1e-12)
})

test_that("an interval-censored fit whose mu overflows is optimized through log(mu) (#566)", {
  # The optimizer hands left- and interval-censored rows to
  # .hzr_logl_weibull(), which must take its alpha / nu as log(mu): mu itself
  # is Inf here. Oracle: the centred fit, the same model.
  d <- .ls_data(200, -0.2)
  status <- d$dead
  lower <- d$time
  upper <- d$time
  ic <- which(status == 0)[1:5]
  status[ic] <- 2L
  lower[ic] <- d$time[ic] / 2
  fit_on <- function(x, theta) {
    suppressWarnings(hazard(time = d$time, status = status,
                            time_lower = lower, time_upper = upper,
                            x = matrix(x, ncol = 1), dist = "weibull",
                            theta = theta, fit = TRUE))
  }
  s_raw <- fit_on(d$x, c(exp(50), 0.2, -0.2))
  s_ctr <- fit_on(d$xc, c(1, 0.2, -0.2))
  expect_true(any(s_raw$data$status == 2))        # premise: interval rows
  expect_identical(unname(coef(s_raw))[1], Inf)    # premise: mu overflows
  expect_true(is.finite(s_raw$fit$objective))
  expect_equal(s_raw$fit$objective, s_ctr$fit$objective, tolerance = 1e-4)
})

test_that("an exact event whose t^(nu - 1) underflows is read on the log scale (#566 acceptance 2)", {
  obj <- hazard(time = 1e-180, status = 1, dist = "weibull",
                theta = c(1e100, 3), fit = FALSE)
  truth <- log(3) + 3 * log(1e100) + 2 * log(1e-180) -
    exp(3 * (log(1e100) + log(1e-180)))
  expect_equal(1e-180^2, 0)   # premise: the natural-scale factor underflows
  got <- .ls_eval(obj, c(1e100, 3))
  # main e7b621c5 returned -Inf with an hzr_evaluate_not_finite warning.
  expect_length(got$msgs, 0L)
  expect_equal(got$ll / truth, 1, tolerance = 1e-12)
})

test_that("a zero-weight event row whose mu^nu underflows leaves the likelihood unchanged (#566 acceptance 3)", {
  both <- hazard(time = c(1, 1), status = c(1, 0), weights = c(0, 1),
                 dist = "weibull", theta = c(1e-200, 2), fit = FALSE)
  alone <- hazard(time = 1, status = 0, dist = "weibull",
                  theta = c(1e-200, 2), fit = FALSE)
  expect_equal(1e-200^2, 0)   # premise: mu^nu underflows
  got_both <- .ls_eval(both, c(1e-200, 2))
  got_alone <- .ls_eval(alone, c(1e-200, 2))
  # main e7b621c5 returned -Inf for `both` (0 * log(0) = NaN on the
  # zero-weight row), with the not-finite warning.
  expect_length(got_both$msgs, 0L)
  expect_identical(got_both$ll, got_alone$ll)
  # Known positive: with weight 1 the event row does move the likelihood.
  # (At mu = 1e-200 hzr_evaluate() refuses a counted event whose mu^nu
  # cannot be represented (#573), so this is read at mu = 1e-100.)
  one <- hazard(time = c(1, 1), status = c(1, 0), weights = c(1, 1),
                dist = "weibull", theta = c(1e-100, 2), fit = FALSE)
  zero <- hazard(time = c(1, 1), status = c(1, 0), weights = c(0, 1),
                 dist = "weibull", theta = c(1e-100, 2), fit = FALSE)
  got_one <- .ls_eval(one, c(1e-100, 2))
  got_zero <- .ls_eval(zero, c(1e-100, 2))
  expect_true(is.finite(got_one$ll))
  expect_equal(got_one$ll - got_zero$ll, log(2) + 2 * log(1e-100),
               tolerance = 1e-12)
})
