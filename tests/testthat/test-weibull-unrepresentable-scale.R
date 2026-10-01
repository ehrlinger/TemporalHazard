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

test_that("a fit whose mu overflows warns, and its readers use log(mu) (#566)", {
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
  # The centred fit predicts these rows; on main fcc6501d the raw fit
  # returned 0 for both, as if that were the survival, and after #573 every
  # reader refused it. Now every reader goes through the log(mu) the fit
  # kept, and agrees with the centred fit, the same model.
  nd_ctr <- transform(nd, xc = x - 1000)
  s_ctr <- predict(ctr$fit, newdata = nd_ctr, type = "survival")
  expect_true(all(s_ctr > 0.05 & s_ctr < 0.95))
  for (type in c("survival", "cumulative_hazard", "hazard",
                 "linear_predictor")) {
    got <- predict(f, newdata = nd, type = type)
    want <- predict(ctr$fit, newdata = nd_ctr, type = type)
    if (type == "linear_predictor") {
      # The linear predictors differ by the centring constant, b * 1000.
      want <- want + unname(coef(ctr$fit)[[3]]) * 1000
    }
    if (type == "hazard") {
      want <- want * exp(unname(coef(ctr$fit)[[3]]) * 1000)
    }
    expect_equal(got / want, rep(1, 2), tolerance = 1e-3, info = type)
  }
  got_ci <- predict(f, newdata = nd, type = "survival", se.fit = TRUE)
  want_ci <- predict(ctr$fit, newdata = nd_ctr, type = "survival",
                     se.fit = TRUE)
  expect_equal(got_ci$se.fit / want_ci$se.fit, rep(1, 2), tolerance = 1e-2)
  expect_equal(got_ci$lower / want_ci$lower, rep(1, 2), tolerance = 1e-2)
  g_raw <- hzr_gof(f)
  g_ctr <- hzr_gof(ctr$fit)
  expect_true(all(is.finite(g_raw$cum_expected)))
  last <- nrow(g_ctr)
  expect_gt(g_ctr$cum_expected[last], 0)
  expect_equal(g_raw$cum_expected[last] / g_ctr$cum_expected[last], 1,
               tolerance = 1e-3)
  dc_raw <- hzr_deciles(f, time = 1)
  dc_ctr <- hzr_deciles(ctr$fit, time = 1)
  expect_equal(dc_raw$expected / dc_ctr$expected, rep(1, nrow(dc_ctr)),
               tolerance = 1e-3)
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
  # Read through log(mu), against the centred fit of the same model. After
  # #573 this refused ("scale mu = 0, which must be positive").
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0.2))
  expect_equal(f$fit$objective, ctr$fit$fit$objective, tolerance = 1e-4)
  got <- predict(f, newdata = nd, type = "survival", se.fit = TRUE)
  want <- predict(ctr$fit, newdata = transform(nd, xc = x - 1000),
                  type = "survival", se.fit = TRUE)
  expect_equal(got$fit / want$fit, rep(1, 2), tolerance = 1e-3)
  expect_equal(got$se.fit / want$se.fit, rep(1, 2), tolerance = 1e-2)
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


test_that("an overflowed variance of mu is read on the log scale, not withheld (#566)", {
  # log(mu) is about 600: mu is finite, but its variance carries mu^2 and is
  # Inf. predict() read that as a fixed parameter and returned a standard
  # error 187 times the centred fit's (16.85 against 0.0901) on main
  # fcc6501d; after #573 it withheld it. The variance of log(mu) is an
  # ordinary number, and the fit keeps it.
  d <- .w566_data(115, -0.115)
  raw <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(500), 0.2, -0.115))
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.115))
  # The fit says so: it is the variance here, not mu, that cannot be held,
  # and the warning says what each reader now does.
  expect_identical(raw$n_scale, 1L)
  expect_match(raw$msgs, "variance of the Weibull scale mu", fixed = TRUE)
  expect_match(raw$msgs, "vcov() gives NA for mu", fixed = TRUE)
  expect_match(raw$msgs, "summary() and predict() use the variance of log(mu)",
               fixed = TRUE)
  expect_true(is.finite(coef(raw$fit)[[1]]))
  # vcov(): mu's row and column are NA, with the reason recorded -- no Inf,
  # no NaN, no false 0. After #573 vcov(raw)[1, 1] was Inf.
  v <- vcov(raw$fit)
  expect_true(all(is.na(v[1, ])) && all(is.na(v[, 1])))
  expect_false(any(is.nan(v)))
  expect_true(all(is.finite(v[-1, -1])))
  expect_match(raw$fit$degraded_causes[["standard_errors"]],
               "log(mu) and its standard error are kept instead", fixed = TRUE)
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  p_ctr <- predict(ctr$fit, newdata = nd, type = "cumulative_hazard",
                   se.fit = TRUE)
  p_raw <- expect_no_warning(predict(raw$fit, newdata = nd,
                                     type = "cumulative_hazard",
                                     se.fit = TRUE))
  # Against the centred fit, the same model: separate optimizer runs, hence
  # the tolerance.
  expect_equal(p_raw$fit / p_ctr$fit, rep(1, 2), tolerance = 1e-2)
  expect_equal(p_raw$se.fit / p_ctr$se.fit, rep(1, 2), tolerance = 1e-2)
  expect_equal(p_raw$lower / p_ctr$lower, rep(1, 2), tolerance = 1e-2)
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

test_that("an underflowed variance of mu is read on the log scale; a NaN one is withheld (#566)", {
  # log(mu) is about -400: mu is an ordinary double, but its variance is
  # subnormal while its covariances are not. The sandwich then gave a
  # standard error of 0.1018 where the centred fit gives 0.0583; after #573
  # it withheld it.
  d <- .w566_data(-80, 0.08)
  raw <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(-400), 0.2, 0.08))
  expect_identical(raw$n_scale, 1L)
  expect_match(raw$msgs, "variance of the Weibull scale mu", fixed = TRUE)
  expect_true(all(is.na(vcov(raw$fit)[1, ])))
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0.08))$fit
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  p <- expect_no_warning(predict(raw$fit, newdata = nd,
                                 type = "cumulative_hazard", se.fit = TRUE))
  p_ctr <- predict(ctr, newdata = nd, type = "cumulative_hazard",
                   se.fit = TRUE)
  expect_equal(p$se.fit / p_ctr$se.fit, rep(1, 2), tolerance = 1e-2)

  # A NaN variance is not a fixed parameter either: it is Inf - Inf in the
  # delta method. A fixed one is NA, and is still dropped without a warning.
  nan <- ctr
  nan$fit$vcov[1, ] <- NaN
  nan$fit$vcov[, 1] <- NaN
  expect_warning(
    p_nan <- predict(nan, newdata = nd, type = "cumulative_hazard",
                     se.fit = TRUE),
    "variance that cannot be represented", fixed = TRUE
  )
  expect_true(all(is.na(p_nan$se.fit)))
  # Marked fixed in fixed_mask: an NA on a parameter that is not fixed is a
  # mask, and withholds the SE (#586).
  fixed <- ctr
  fixed$fit$vcov[1, ] <- NA_real_
  fixed$fit$vcov[, 1] <- NA_real_
  fixed$fit$fixed_mask <- c(TRUE, FALSE, FALSE)
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
  # Warn only (#566, decision B): the result has the same shape as the
  # centred fit's, with no column or element added for log(mu).
  expect_identical(names(bs$value), names(bs_ctr$value))
  for (el in names(bs$value)) {
    if (is.data.frame(bs$value[[el]])) {
      expect_identical(names(bs$value[[el]]), names(bs_ctr$value[[el]]),
                       info = el)
    }
  }
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

test_that("the event hazard's mu^nu is guarded too (#566, Copilot on #573)", {
  # The likelihood forms mu in two ways: mu * time, guarded above, and mu^nu
  # in an exact event's hazard. With mu = 4e-162 and nu = 2, mu * t is 1 but
  # mu^nu is subnormal, and hzr_evaluate() returned -372.0158 where the
  # closed form gives -371.9393, with no warning.
  one_event <- function(mu, nu, t, status = 1) {
    obj <- hazard(time = t, status = status, dist = "weibull",
                  theta = c(mu, nu), fit = FALSE)
    as.numeric(suppressWarnings(hzr_evaluate(obj, c(mu, nu)))$logLik)
  }
  closed <- function(mu, nu, t) {
    log(nu) + nu * log(mu) + (nu - 1) * log(t) - exp(nu * (log(mu) + log(t)))
  }
  expect_true(4e-162^2 < .Machine$double.xmin)
  expect_error(one_event(4e-162, 2, 2.5e161),
               "cannot be represented, so the event hazard", fixed = TRUE)
  # Known negative: mu^nu = 1e-200, a normal double, and the evaluation is
  # the closed form.
  expect_equal(one_event(1e-100, 2, 1e100), closed(1e-100, 2, 1e100),
               tolerance = 1e-12)
  # A right-censored row never forms the hazard, so mu^nu is not read there.
  expect_true(is.finite(one_event(4e-162, 2, 2.5e161, status = 0)))
  # An overflowing mu^nu is not refused: hzr_evaluate() reports a likelihood
  # that is not finite as -Inf, with a warning (test-evaluate-at-parameters.R).
  obj <- hazard(time = 3, status = 1, dist = "weibull", theta = c(2, 1e5),
                fit = FALSE)
  expect_warning(ll <- hzr_evaluate(obj, c(2, 1e5))$logLik,
                 class = "hzr_evaluate_not_finite")
  expect_identical(as.numeric(ll), -Inf)
})

test_that("prediction SEs: a negative variance is named, and an unused one is dropped (#566, Copilot on #573)", {
  data(avc, package = "TemporalHazard", envir = environment())
  a <- stats::na.omit(avc)
  f <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age + mal,
                               data = a, dist = "weibull",
                               theta = c(0.01, 0.5, 0, 0), fit = TRUE))
  nd <- data.frame(time = c(1, 5), age = c(60, 70), mal = c(0, 0))
  se_of <- function(obj, newdata) {
    predict(obj, newdata = newdata, type = "cumulative_hazard",
            se.fit = TRUE)$se.fit
  }
  base_se <- se_of(f, nd)
  expect_true(all(is.finite(base_se) & base_se > 0))

  # A negative variance is not an overflow: the covariance is not positive
  # definite. The standard error is withheld and the warning says why (on
  # main fcc6501d it was reported as 0, with no warning).
  neg <- f
  neg$fit$vcov[1, 1] <- -abs(neg$fit$vcov[1, 1])
  expect_warning(se_neg <- se_of(neg, nd), "negative variance", fixed = TRUE)
  expect_true(all(is.na(se_neg)))

  # mal is 0 in every requested row, so the prediction does not depend on its
  # coefficient: an unrepresentable variance there is dropped exactly, and
  # the SE equals the one with mal treated as fixed.
  inf_mal <- f
  inf_mal$fit$vcov[4, 4] <- Inf
  fixed_mal <- f
  fixed_mal$fit$vcov[4, ] <- NA_real_
  fixed_mal$fit$vcov[, 4] <- NA_real_
  se_inf <- expect_no_warning(se_of(inf_mal, nd))
  expect_equal(se_inf, se_of(fixed_mal, nd), tolerance = 1e-12)
  expect_equal(se_inf, base_se, tolerance = 1e-12)
  # Known positive: with mal = 1 the prediction depends on it, and the SE is
  # withheld.
  expect_warning(se_used <- se_of(inf_mal, transform(nd, mal = 1)),
                 "variance that cannot be represented", fixed = TRUE)
  expect_true(all(is.na(se_used)))
})

test_that("a numerically zero Jacobian column does not exempt a parameter (#566, Copilot on #573)", {
  data(avc, package = "TemporalHazard", envir = environment())
  a <- stats::na.omit(avc)
  f <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age + mal,
                               data = a, dist = "weibull",
                               theta = c(0.01, 0.5, 0, 0), fit = TRUE))
  b_age <- unname(coef(f))[3]
  # eta = -1000 in every row: the relative hazard exp(eta) underflows to 0,
  # so age's Jacobian column is exactly 0, yet the prediction depends on it.
  nd <- data.frame(time = c(1, 5), age = -1000 / b_age, mal = 0)
  hz_of <- function(obj) {
    predict(obj, newdata = nd, type = "hazard", se.fit = TRUE)
  }
  base <- expect_no_warning(hz_of(f))
  expect_identical(base$fit, c(0, 0))
  inf_age <- f
  inf_age$fit$vcov[3, 3] <- Inf
  # On 0ccf028d this dropped age as unused and returned se.fit 0, silently.
  expect_warning(got <- hz_of(inf_age), "variance that cannot be represented",
                 fixed = TRUE)
  expect_true(all(is.na(got$se.fit)))
  # Known negative: mal is 0 in every row, so its Inf variance is still
  # dropped exactly.
  inf_mal <- f
  inf_mal$fit$vcov[4, 4] <- Inf
  expect_identical(expect_no_warning(hz_of(inf_mal))$se.fit, base$se.fit)
})

test_that("a multiphase coefficient whose covariate is 0 in every row is dropped (#566)", {
  set.seed(3)
  n <- 120
  df <- data.frame(time = stats::rexp(n, 0.3),
                   status = stats::rbinom(n, 1, 0.7),
                   x = stats::rnorm(n))
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = df, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                    fixed = "shapes", formula = ~ x),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1)))
  k <- grep("^early.*x$", names(fit$fit$theta))
  expect_length(k, 1L)
  ch <- function(o, x) {
    predict(o, newdata = data.frame(time = c(0.5, 2), x = x),
            type = "cumulative_hazard", se.fit = TRUE)$se.fit
  }
  base <- ch(fit, 0)
  expect_true(all(is.finite(base) & base > 0))
  bad <- fit
  bad$fit$vcov[k, k] <- Inf
  # x is 0 in every row, so the prediction does not depend on its
  # coefficient: the Inf variance drops out exactly.
  expect_equal(expect_no_warning(ch(bad, 0)), base, tolerance = 1e-12)
  # Known positive: with x = 1 it does, and the SE is withheld.
  expect_warning(se1 <- ch(bad, 1), "variance that cannot be represented",
                 fixed = TRUE)
  expect_true(all(is.na(se1)))
})

test_that("a negative variance of mu is not called an overflow at fit time (#566, Copilot on #573)", {
  expect_false(.hzr_variance_unrepresentable(-1))
  expect_false(.hzr_variance_unrepresentable(-1e-320))
  # Known positive: zero and subnormal are underflow.
  expect_true(.hzr_variance_unrepresentable(0))
  expect_true(.hzr_variance_unrepresentable(1e-320))

  d <- .w566_data(-2, 0.001)
  real_optim <- .hzr_optim_weibull
  testthat::local_mocked_bindings(
    .hzr_optim_weibull = function(...) {
      res <- real_optim(...)
      res$vcov[1, 1] <- -abs(res$vcov[1, 1])
      res
    }
  )
  res <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0))
  expect_true(res$fit$fit$vcov[1, 1] < 0)
  expect_identical(res$n_scale, 0L)
  expect_true(is.na(res$fit$fit$se[[1]]))
})

test_that("an unusable vcov withholds SEs before any Jacobian is built (#566, Copilot on #573)", {
  d <- .w566_data(-2, 0.001)
  f <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, 0))$fit
  nd <- d[1:2, ]
  n_jac <- 0L
  real_jac <- .hzr_predict_jacobian_weibull
  testthat::local_mocked_bindings(
    .hzr_predict_jacobian_weibull = function(...) {
      n_jac <<- n_jac + 1L
      real_jac(...)
    }
  )
  # Known positive: a usable vcov builds the Jacobian.
  se_ok <- predict(f, newdata = nd, type = "cumulative_hazard",
                   se.fit = TRUE)$se.fit
  expect_true(all(is.finite(se_ok)))
  expect_identical(n_jac, 1L)
  g <- f
  g$fit$vcov <- g$fit$vcov[1:2, 1:2]
  expect_warning(se_bad <- predict(g, newdata = nd, type = "cumulative_hazard",
                                   se.fit = TRUE)$se.fit,
                 "Variance-covariance matrix is unavailable", fixed = TRUE)
  expect_true(all(is.na(se_bad)))
  expect_identical(n_jac, 1L)
})

test_that("decomposed SEs are screened per phase (#566, Copilot on #573)", {
  set.seed(17)
  df <- data.frame(time = stats::rexp(60, 0.3),
                   status = stats::rbinom(60, 1, 0.6))
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                   fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- suppressWarnings(hazard(survival::Surv(time, status) ~ 1, data = df,
                                 dist = "multiphase", phases = phases,
                                 fit = TRUE, control = list(n_starts = 2)))
  k <- grep("^constant", names(fit$fit$theta))
  expect_length(k, 1L)
  se_by <- function(o) {
    msg <- character(0)
    r <- withCallingHandlers(
      predict(o, newdata = data.frame(time = c(0.5, 1, 2)),
              type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE),
      warning = function(w) {
        msg <<- c(msg, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(se = split(r$se.fit, r$component), msg = msg)
  }
  base <- se_by(fit)
  expect_length(base$msg, 0L)
  expect_true(all(is.finite(unlist(base$se)) & unlist(base$se) > 0))

  bad <- fit
  bad$fit$vcov[k, k] <- Inf
  fx <- fit
  fx$fit$vcov[k, ] <- NA_real_
  fx$fit$vcov[, k] <- NA_real_
  got <- se_by(bad)
  # The early phase does not depend on the constant phase's parameter, so its
  # SE is the one with that parameter fixed, and the one from the fit itself.
  expect_equal(got$se$early, se_by(fx)$se$early, tolerance = 1e-12)
  expect_equal(got$se$early, base$se$early, tolerance = 1e-12)
  # The phase that depends on it, and the total, are withheld, under one
  # warning. On main 263c7057 they were 0 and the early phase's SE, silently.
  expect_true(all(is.na(got$se$constant)))
  expect_true(all(is.na(got$se$total)))
  expect_length(got$msg, 1L)
  expect_match(got$msg, "variance that cannot be represented", fixed = TRUE)

  # An unusable vcov withholds every SE before the Jacobians are built.
  n_jac <- 0L
  real_jac <- .hzr_predict_jacobian_multiphase
  testthat::local_mocked_bindings(
    .hzr_predict_jacobian_multiphase = function(...) {
      n_jac <<- n_jac + 1L
      real_jac(...)
    }
  )
  se_by(fit)
  expect_identical(n_jac, 1L)  # known positive
  short <- fit
  short$fit$vcov <- short$fit$vcov[-k, -k]
  got <- se_by(short)
  expect_true(all(is.na(unlist(got$se))))
  expect_length(got$msg, 1L)
  expect_match(got$msg, "Variance-covariance matrix is unavailable",
               fixed = TRUE)
  expect_identical(n_jac, 1L)
})
