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
  classes <- list()
  fit <- withCallingHandlers(
    hazard(formula, data = data, dist = "weibull", theta = theta, fit = TRUE),
    warning = function(w) {
      classes[[length(classes) + 1L]] <<- class(w)
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit,
       n_scale = sum(vapply(classes, function(k) {
         "hzr_unrepresentable_scale" %in% k
       }, logical(1))))
}

test_that("a fit whose mu overflows warns, and its readers refuse (#566)", {
  d <- .w566_data(200, -0.2)
  res <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(50), 0.2, -0.2))
  f <- res$fit
  # Premise: a genuine maximum, the centred fit's, with mu = Inf.
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.2))
  expect_identical(ctr$n_scale, 0L)
  expect_equal(f$fit$objective, ctr$fit$fit$objective, tolerance = 1e-8)
  expect_lte(f$fit$rel_gradient, .Machine$double.eps^(1 / 3))
  expect_identical(unname(coef(f)[[1]]), Inf)
  # It stays converged: the fit is sound, only its scale is not reportable.
  expect_identical(f$fit$converged, TRUE)
  expect_identical(res$n_scale, 1L)

  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  # The centred fit gives the true predictions for these rows.
  expect_equal(predict(ctr$fit, newdata = nd, type = "survival"),
               c(0.21701522, 0.41077904), tolerance = 1e-6)
  for (type in c("survival", "cumulative_hazard")) {
    expect_error(predict(f, newdata = nd, type = type),
                 "scale mu = Inf, which must be finite", fixed = TRUE)
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
  # log(mu) is about 600: mu is finite, so the fit raises no warning, but its
  # variance carries mu^2 and is Inf. predict() read that as a fixed
  # parameter and returned a standard error 187 times the centred fit's
  # (16.85 against 0.0901) on main fcc6501d.
  d <- .w566_data(115, -0.115)
  raw <- .w566_fit(survival::Surv(time, dead) ~ x, d, c(exp(500), 0.2, -0.115))
  ctr <- .w566_fit(survival::Surv(time, dead) ~ xc, d, c(1, 0.2, -0.115))
  expect_identical(raw$n_scale, 0L)
  expect_true(is.finite(coef(raw$fit)[[1]]))
  expect_identical(unname(vcov(raw$fit)[1, 1]), Inf)
  nd <- d[1:2, ]
  nd$time <- c(0.5, 2)
  p_ctr <- predict(ctr$fit, newdata = nd, type = "cumulative_hazard",
                   se.fit = TRUE)
  expect_warning(
    p_raw <- predict(raw$fit, newdata = nd, type = "cumulative_hazard",
                     se.fit = TRUE),
    "infinite variance", fixed = TRUE
  )
  # The value is right; its standard error and limits are withheld.
  expect_equal(p_raw$fit, p_ctr$fit, tolerance = 1e-3)
  expect_true(all(is.na(p_raw$se.fit)))
  expect_true(all(is.na(p_raw$lower)) && all(is.na(p_raw$upper)))
  expect_true(all(is.finite(p_ctr$se.fit)))
})
