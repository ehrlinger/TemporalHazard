# A Weibull scale mu the double range cannot hold (#566).
#
# The model depends on mu only through nu * log(mu). With a covariate far from
# zero the intercept absorbs -beta * mean(x), so log(mu) can leave exp()'s
# range (about +/-709) at a genuine maximum: mu is then reported as Inf or 0.
# On main fcc6501d the Inf case was silent: predict() returned survival 0 and
# cumulative hazard Inf, and hzr_gof() reported E = Inf, with no warning.
#
# The oracle is the same model with the covariate centred, whose mu is
# representable and whose log-likelihood is identical.

.w566_data <- function(alpha, b) {
  set.seed(1)
  n <- 400
  x <- 1000 + stats::rnorm(n, sd = 5)
  t <- (stats::rexp(n) / exp(alpha + b * x))^(1 / 0.2)
  cens <- stats::rexp(n)^5 * 3
  data.frame(time = pmin(t, cens), dead = as.integer(t <= cens),
             x = x, xc = x - 1000)
}

.w566_fit <- function(formula, data, theta) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    hazard(formula, data = data, dist = "weibull", theta = theta, fit = TRUE),
    warning = function(w) {
      if (inherits(w, "hzr_unrepresentable_scale")) {
        msgs <<- c(msgs, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  # The hzr_unrepresentable_scale warnings this fit raised, and their text.
  list(fit = fit, n_scale = length(msgs), msgs = msgs)
}

test_that("a fit whose mu overflows warns, and its readers refuse (#566)", {
  d <- .w566_data(200, -0.2)
  res <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(50), 0.2, -0.2))
  f <- res$fit
  # Premise: a genuine maximum, the centred fit's, with mu = Inf. The two
  # optimizer runs need not stop at the same digit, so the tolerance is loose.
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.2))
  expect_identical(ctr$n_scale, 0L)
  expect_equal(f$fit$objective, ctr$fit$fit$objective, tolerance = 1e-4)
  expect_identical(unname(coef(f)[[1]]), Inf)
  # It stays converged: the fit is sound, only its scale is not reportable.
  expect_identical(f$fit$converged, TRUE)
  # One warning, and it is the one about mu: the variance of an Inf mu is
  # NaN, which would raise the variance warning if this one did not.
  expect_identical(res$n_scale, 1L)
  expect_match(res$msgs, "scale mu is reported as Inf", fixed = TRUE)

  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  # The centred fit predicts these rows; on main the raw fit returned 0 for
  # both, as if that were the survival.
  s_ctr <- predict(ctr$fit, newdata = nd, type = "survival")
  expect_true(all(s_ctr > 0.05 & s_ctr < 0.95))
  # Every type is refused, including the two that never read mu.
  for (type in c("survival", "cumulative_hazard", "hazard",
                 "linear_predictor")) {
    expect_error(predict(f, newdata = nd, type = type),
                 "scale mu = Inf, which cannot be represented", fixed = TRUE)
  }
  expect_error(predict(f, type = "cumulative_hazard"),
               "scale mu = Inf", fixed = TRUE)
  expect_error(hzr_gof(f), "scale mu = Inf", fixed = TRUE)
  expect_error(hzr_deciles(f, time = 1), "scale mu = Inf", fixed = TRUE)
})

test_that("a fit whose mu underflows to 0 warns too (#566)", {
  d <- .w566_data(-200, 0.2)
  res <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(-700), 0.2, 0.2))
  f <- res$fit
  expect_identical(unname(coef(f)[[1]]), 0)
  expect_identical(f$fit$converged, TRUE)
  expect_identical(res$n_scale, 1L)
  expect_match(res$msgs, "scale mu is reported as 0", fixed = TRUE)
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  expect_error(predict(f, newdata = nd, type = "survival"),
               "scale mu = 0, which must be positive", fixed = TRUE)
})

test_that("an intercept-only object is checked without a stored design (#566)", {
  # predict() skipped the theta check when the object stored no design, so a
  # fitted scale of 0 there returned NA. A fit cannot be built with such a
  # theta, so one is placed on a sound object.
  data(avc, package = "TemporalHazard", envir = environment())
  a <- stats::na.omit(avc)
  f <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1, data = a,
                               dist = "weibull", theta = c(0.01, 0.5),
                               fit = TRUE))
  expect_null(f$data$x)
  nd <- data.frame(time = c(1, 10))
  expect_true(all(is.finite(predict(f, newdata = nd, type = "survival"))))
  for (bad in list(c(Inf, 0.5), c(0, 0.5), c(0.01, Inf))) {
    g <- f
    g$fit$theta[] <- bad
    expect_error(predict(g, newdata = nd, type = "survival"),
                 "'theta' gives Weibull", fixed = TRUE)
  }
})

test_that("mu * time overflowing does not make the cumulative hazard Inf (#566)", {
  # mu is finite here, but mu * time is not: (mu * t)^nu was Inf where the
  # cumulative hazard is finite. On main the second row below gave Inf.
  data(avc, package = "TemporalHazard", envir = environment())
  a <- stats::na.omit(avc)
  obj <- hazard(survival::Surv(int_dead, dead) ~ age, data = a,
                dist = "weibull", theta = c(1e300, 0.2, -0.5), fit = FALSE)
  nd <- data.frame(time = c(1, 1e10), age = c(300, 300))
  expect_false(is.finite(1e300 * 1e10))
  truth <- exp(0.2 * (log(1e300) + log(nd$time)) - 0.5 * nd$age)
  got <- suppressWarnings(predict(obj, newdata = nd, type = "cumulative_hazard"))
  expect_equal(got, truth, tolerance = 1e-10)
  ci <- suppressWarnings(predict(obj, newdata = nd, type = "cumulative_hazard",
                                 se.fit = TRUE))
  expect_equal(ci$fit, truth, tolerance = 1e-10)

  # With a standard error: the centred fit, rewritten exactly in raw-x
  # coordinates, is the same model, so its predictions and their standard
  # errors are the oracle. beta * xc = beta * x - 1000 * beta moves the
  # intercept: log(mu_raw) = log(mu_c) - 1000 * beta / nu, about 300 here, so
  # mu_raw and its variance are finite while mu_raw * 1e200 is not.
  d <- .w566_data(60, -0.06)
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.06))$fit
  raw <- .w566_fit(survival::Surv(time, dead) ~ x, d,
                   c(exp(300), 0.2, -0.06))$fit
  th <- unname(coef(ctr))
  mu_raw <- exp(log(th[1]) - 1000 * th[3] / th[2])
  jac <- rbind(c(mu_raw / th[1], mu_raw * 1000 * th[3] / th[2]^2,
                 -mu_raw * 1000 / th[2]),
               c(0, 1, 0), c(0, 0, 1))
  raw$fit$theta[] <- c(mu_raw, th[2], th[3])
  raw$fit$vcov[] <- jac %*% vcov(ctr) %*% t(jac)
  expect_true(is.finite(mu_raw) && all(is.finite(raw$fit$vcov)))
  expect_false(is.finite(mu_raw * 1e200))
  ndf <- d[1:3, ]
  ndf$time <- c(0.5, 2, 1e200)
  ndc <- ndf
  for (type in c("cumulative_hazard", "survival")) {
    p_raw <- predict(raw, newdata = ndf, type = type, se.fit = TRUE)
    p_ctr <- predict(ctr, newdata = ndc, type = type, se.fit = TRUE)
    expect_true(all(is.finite(p_ctr$fit)) && all(is.finite(p_ctr$se.fit)))
    expect_equal(p_raw$fit, p_ctr$fit, tolerance = 1e-8)
    expect_equal(p_raw$se.fit, p_ctr$se.fit, tolerance = 1e-6)
  }

  # Ordinary values are unchanged: the log-scale form against the power form.
  ord <- hazard(survival::Surv(int_dead, dead) ~ age, data = a,
                dist = "weibull", theta = c(0.0123, 1.7, 0.01), fit = FALSE)
  nd2 <- data.frame(time = c(1e-3, 0.5, 3, 170, 0), age = c(1, 20, 40, 60, 5))
  expect_equal(
    suppressWarnings(predict(ord, newdata = nd2, type = "cumulative_hazard")),
    (0.0123 * nd2$time)^1.7 * exp(0.01 * nd2$age), tolerance = 1e-12
  )
})


test_that("an overflowed variance gives NA standard errors, not wrong ones (#566)", {
  # log(mu) is about 600: mu is finite, but its variance carries mu^2 and is
  # Inf. predict() read that as a fixed parameter and returned a standard
  # error 187 times the centred fit's (16.85 against 0.0901) on main
  # fcc6501d, and summary() showed the SE of mu with no warning.
  d <- .w566_data(115, -0.115)
  raw <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(500), 0.2, -0.115))
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.115))
  # The fit says so: it is the variance here, not mu, that cannot be held.
  expect_identical(raw$n_scale, 1L)
  expect_match(raw$msgs, "variance of the Weibull scale mu", fixed = TRUE)
  expect_true(is.finite(coef(raw$fit)[[1]]))
  expect_identical(unname(vcov(raw$fit)[1, 1]), Inf)
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  p_ctr <- predict(ctr$fit, newdata = nd, type = "cumulative_hazard",
                   se.fit = TRUE)
  expect_warning(
    p_raw <- predict(raw$fit, newdata = nd, type = "cumulative_hazard",
                     se.fit = TRUE),
    "variance that cannot be represented", fixed = TRUE
  )
  # The value is right; its standard error and limits are withheld. The two
  # fits are separate optimizer runs, hence the tolerance.
  expect_equal(p_raw$fit, p_ctr$fit, tolerance = 1e-2)
  expect_true(all(is.na(p_raw$se.fit)))
  expect_true(all(is.na(p_raw$lower)) && all(is.na(p_raw$upper)))
  expect_true(all(is.finite(p_ctr$se.fit)))

  # The linear predictor and the relative hazard do not depend on mu, so
  # their standard errors are |x| * se(beta), and exp(eta) times that,
  # whatever mu's variance is; they are still returned.
  se_eta <- abs(nd$x) * sqrt(vcov(raw$fit)[3, 3])
  lp <- expect_no_warning(predict(raw$fit, newdata = nd,
                                  type = "linear_predictor", se.fit = TRUE))
  expect_equal(lp$se.fit, se_eta, tolerance = 1e-10)
  hz <- expect_no_warning(predict(raw$fit, newdata = nd, type = "hazard",
                                  se.fit = TRUE))
  expect_equal(hz$se.fit, hz$fit * se_eta, tolerance = 1e-10)
})

test_that("an underflowed or NaN variance is withheld too (#566)", {
  # log(mu) is about -400: mu is an ordinary double, but its variance is
  # subnormal while its covariances are not. The sandwich then gave a
  # standard error of 0.1018 where the centred fit gives 0.0583.
  d <- .w566_data(-80, 0.08)
  raw <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(-400), 0.2, 0.08))
  expect_identical(raw$n_scale, 1L)
  expect_match(raw$msgs, "variance of the Weibull scale mu", fixed = TRUE)
  v11 <- unname(vcov(raw$fit)[1, 1])
  expect_true(!is.na(v11) && v11 < .Machine$double.xmin)
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  expect_warning(
    p <- predict(raw$fit, newdata = nd, type = "cumulative_hazard",
                 se.fit = TRUE),
    "variance that cannot be represented", fixed = TRUE
  )
  expect_true(all(is.finite(p$fit)) && all(is.na(p$se.fit)))

  # A NaN variance is not a fixed parameter either: it is Inf - Inf in the
  # delta method. A fixed one is NA, and is still dropped without a warning.
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0.08))$fit
  nan <- ctr
  nan$fit$vcov[1, ] <- NaN
  nan$fit$vcov[, 1] <- NaN
  expect_warning(
    p_nan <- predict(nan, newdata = nd, type = "cumulative_hazard",
                     se.fit = TRUE),
    "variance that cannot be represented", fixed = TRUE
  )
  expect_true(all(is.na(p_nan$se.fit)))
  fixed <- ctr
  fixed$fit$vcov[1, ] <- NA_real_
  fixed$fit$vcov[, 1] <- NA_real_
  p_fixed <- expect_no_warning(predict(fixed, newdata = nd,
                                       type = "cumulative_hazard",
                                       se.fit = TRUE))
  expect_true(all(is.finite(p_fixed$se.fit)))
})

test_that("a subnormal mu is refused like one that reached 0 (#566)", {
  # A subnormal mu is positive, so it passed the positivity rule, but it has
  # lost digits: at 4.94e-324 predict() gave a cumulative hazard of 0.4888
  # where the centred fit gives 0.5243.
  #
  # The value is put in place, not fitted for. The subnormals span only
  # exp(-745) to exp(-708), and a raw-covariate fit is ill-conditioned, so
  # where it stops differs by platform: the same fit gave mu = 5.7e-313 on
  # macOS and 2.4e-306, a normal double, on Linux.
  sub <- 5e-313
  expect_true(sub > 0 && sub < .Machine$double.xmin)
  d <- .w566_data(-150, 0.15)
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0.15))$fit
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  g <- ctr
  g$fit$theta[[1]] <- sub
  expect_error(predict(g, newdata = nd, type = "cumulative_hazard"),
               "scale mu = 5e-313, which cannot be represented", fixed = TRUE)
  expect_error(hzr_evaluate(ctr, c(sub, unname(coef(ctr))[2:3])),
               "scale mu = 5e-313, which cannot be represented", fixed = TRUE)

  # The fit-time warning, with the optimizer's own result carrying the
  # subnormal scale: hazard() reads mu from what its optimizer returns.
  real_optim <- .hzr_optim_weibull
  testthat::local_mocked_bindings(
    .hzr_optim_weibull = function(...) {
      res <- real_optim(...)
      res$par[[1]] <- sub
      res
    }
  )
  res <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0.15))
  expect_identical(unname(coef(res$fit)[[1]]), sub)
  expect_identical(res$n_scale, 1L)
  expect_match(res$msgs, "scale mu is reported as 5e-313", fixed = TRUE)
})

test_that("hzr_bootstrap() counts replicates whose mu cannot be represented (#566)", {
  collect <- function(expr) {
    msgs <- character(0)
    value <- withCallingHandlers(expr, warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, msgs = msgs)
  }
  pattern <- "successful replicates reported a Weibull scale mu"
  d <- .w566_data(200, -0.2)
  f <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(50), 0.2, -0.2))$fit
  bs <- collect(hzr_bootstrap(f, n_boot = 3, seed = 3))
  hit <- grep(pattern, bs$msgs, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  # The count is of replicates that ran, and here every one overflows.
  n_ok <- bs$value$n_success
  expect_gt(n_ok, 0L)
  expect_match(hit, paste0("^", n_ok, " of ", n_ok, " successful"))
  # Known negative: the centred fit's replicates are all representable.
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.2))$fit
  bs_ctr <- collect(hzr_bootstrap(ctr, n_boot = 3, seed = 3))
  expect_false(any(grepl(pattern, bs_ctr$msgs, fixed = TRUE)))
  # A mu that reached exactly 0 is counted too.
  d0 <- .w566_data(-200, 0.2)
  f0 <- .w566_fit(survival::Surv(time, dead) ~ x, d0,
                  c(exp(-700), 0.2, 0.2))$fit
  expect_identical(unname(coef(f0)[[1]]), 0)
  bs0 <- collect(hzr_bootstrap(f0, n_boot = 3, seed = 3))
  hit0 <- grep(pattern, bs0$msgs, fixed = TRUE, value = TRUE)
  expect_length(hit0, 1L)
  n0 <- bs0$value$n_success
  expect_match(hit0, paste0("^", n0, " of ", n0, " successful"))
})

test_that("the rule's edges: the smallest normal double, and a subnormal nu (#566)", {
  xmin <- .Machine$double.xmin
  expect_false(.hzr_unrepresentable(xmin))
  expect_true(.hzr_unrepresentable(xmin / 2))
  expect_false(.hzr_unrepresentable(.Machine$double.xmax))
  expect_true(.hzr_unrepresentable(Inf))
  expect_true(.hzr_unrepresentable(NaN))
  # Zero and negative values are the positivity rule's, which names them.
  expect_false(.hzr_unrepresentable(0))
  expect_false(.hzr_unrepresentable(-1))
  expect_error(.hzr_check_theta(c(1, xmin / 2), "weibull"),
               "shape nu = .* cannot be represented")
  expect_error(.hzr_check_theta(c(0, 1), "weibull"),
               "scale mu = 0, which must be positive", fixed = TRUE)
  expect_null(.hzr_check_theta(c(xmin, 1), "weibull"))
  # The other families report on the optimizer's own scale: no rule applies.
  expect_null(.hzr_check_theta(c(Inf, 1), "lognormal"))
})

test_that("hzr_evaluate() refuses a product mu * time it cannot hold (#566)", {
  # mu is an ordinary double here, so the rule on mu alone passes it, but the
  # likelihood forms (mu * t)^nu and mu * t underflows: it is exactly 0 on
  # some rows and subnormal on others. On a fit of this data whose mu was
  # about 1e-291, hzr_evaluate() returned 23084.15 at the fit's own
  # estimates, whose log-likelihood is 22972.56, with no warning.
  #
  # The parameters are supplied, not fitted: the data are seeded, so which
  # products underflow is the same on every platform.
  set.seed(4)
  n <- 500
  age <- stats::rnorm(n, 60, 10)
  nu <- 0.05
  b <- 0.57
  t <- (stats::rexp(n) / exp(b * (age - 60)))^(1 / nu)
  cc <- (stats::rexp(n) * 2)^(1 / nu)
  d <- data.frame(time = pmin(t, cc), dead = as.integer(t <= cc), age = age,
                  agec = age - 60)
  mu <- 1e-290
  raw <- hazard(survival::Surv(time, dead) ~ age, data = d, dist = "weibull",
                theta = c(mu, nu, b), fit = FALSE)
  # Premise: mu itself is a normal double, and the product is not, on some
  # rows but not all.
  expect_true(mu >= .Machine$double.xmin)
  n_lost <- sum(mu * d$time < .Machine$double.xmin)
  expect_true(n_lost > 0L && n_lost < n)
  expect_error(hzr_evaluate(raw, c(mu, nu, b)),
               paste0("mu * time cannot be represented for ", n_lost, " of ",
                      n, " times"), fixed = TRUE)
  # Known negative: the same model in centred coordinates, where mu is near
  # 1 and no product is lost, evaluates to a finite log-likelihood.
  ctr <- hazard(survival::Surv(time, dead) ~ agec, data = d,
                dist = "weibull", theta = c(1, nu, b), fit = FALSE)
  expect_true(is.finite(as.numeric(hzr_evaluate(ctr, c(1, nu, b))$logLik)))
})

test_that("the product guard reads only the times each row's status uses (#566)", {
  # The likelihood reads `time` on an event or right-censored row, plus its
  # entry time when there is one; the upper bound on a left-censored row; and
  # both bounds on an interval-censored row. A bound it does not read must
  # not refuse the evaluation (Copilot, #573).
  mu <- 1e-290
  nu <- 0.5
  tm <- c(1, 2, 3, 4, 5, 6)
  ev <- function(...) {
    obj <- hazard(time = tm, ..., dist = "weibull", theta = c(mu, nu),
                  fit = FALSE)
    as.numeric(suppressWarnings(hzr_evaluate(obj, c(mu, nu)))$logLik)
  }
  tiny <- 1e-100
  expect_false(mu * tiny >= .Machine$double.xmin)
  right <- c(1, 0, 1, 0, 1, 0)
  reference <- ev(status = right)
  expect_true(is.finite(reference))
  # Unread: an upper bound on event and right-censored rows.
  expect_identical(ev(status = right, time_upper = rep(tiny, 6)), reference)
  # Unread: a lower bound on a left-censored row.
  left <- c(1, 0, 1, 0, 1, -1)
  with_lower <- ev(status = left, time_lower = c(rep(0, 5), tiny),
                   time_upper = tm)
  expect_identical(with_lower, ev(status = left, time_upper = tm))
  # Read: a left-censored row's upper bound, an entry time on an event row,
  # and either bound of an interval-censored row. One time is lost in each.
  lost <- "mu * time cannot be represented for 1 of "
  expect_error(ev(status = left, time_upper = c(tm[1:5], tiny)), lost,
               fixed = TRUE)
  expect_error(ev(status = right, time_lower = c(tiny, rep(0, 5))), lost,
               fixed = TRUE)
  interval <- c(1, 0, 1, 0, 1, 2)
  expect_error(ev(status = interval, time_lower = c(rep(0, 5), tiny),
                  time_upper = tm), lost, fixed = TRUE)

  # A row with weight 0 contributes nothing, whatever its time, so its time
  # is not read either. With weight 1 the same row is refused.
  ev_w <- function(time, weights) {
    obj <- hazard(time = time, status = right, weights = weights,
                  dist = "weibull", theta = c(mu, nu), fit = FALSE)
    as.numeric(suppressWarnings(hzr_evaluate(obj, c(mu, nu)))$logLik)
  }
  w0 <- c(1, 1, 1, 1, 1, 0)
  expect_identical(ev_w(c(tm[1:5], tiny), w0), ev_w(tm, w0))
  expect_error(ev_w(c(tm[1:5], tiny), rep(1, 6)), lost, fixed = TRUE)
})
