# A row at time 0 has a well-defined likelihood contribution (#341). Right
# censored at 0, it contributes log S(0) = 0. An interval (0, u] contributes
# log(1 - S(u)), which is the left-censored row at u. SAS/C reads an interval
# whose earlier time is 0 the same way (setcoe_obs_loop.c). The single-
# distribution likelihoods refused such a row for the whole data set, so one
# of them turned a fit into a clamped objective at its starting values,
# reported as converged.

zero_rows_data <- function() {
  withr::local_seed(20260917)
  n <- 60
  t <- stats::rweibull(n, shape = 1.3, scale = 4)
  visit <- stats::runif(n, 0.5, 3)
  # Interval-censored between visits; some rows first seen after the event,
  # so their interval opens at 0.
  lower <- floor(t / visit) * visit
  upper <- lower + visit
  data.frame(l = lower, u = upper)
}

zero_rows_starts <- list(exponential = -1, weibull = c(0.3, 1),
                         lognormal = c(1, 0), loglogistic = c(-1, 0))

zero_rows_fit <- function(d, dist) {
  starts <- zero_rows_starts
  hazard(survival::Surv(l, u, type = "interval2") ~ 1, data = d,
         dist = dist, theta = starts[[dist]], fit = TRUE)
}

test_that("a right-censored row at time 0 contributes nothing (#341)", {
  d <- zero_rows_data()
  expect_gt(sum(d$l == 0), 0)
  with_zero <- rbind(d, data.frame(l = 0, u = NA))
  for (dist in c("exponential", "weibull", "lognormal", "loglogistic")) {
    without <- zero_rows_fit(d, dist)
    with <- zero_rows_fit(with_zero, dist)
    expect_true(isTRUE(without$fit$converged), info = dist)
    expect_gt(without$fit$objective, -1e9)
    # A clamped objective stays at the start: the fit must have moved.
    expect_gt(max(abs(unname(with$fit$theta) - zero_rows_starts[[dist]])),
              1e-3)
    expect_equal(with$fit$objective, without$fit$objective,
                 tolerance = 1e-8, info = dist)
    expect_equal(with$fit$theta, without$fit$theta, tolerance = 1e-6,
                 info = dist)
  }
})

test_that("an interval opening at 0 is left censoring at its upper bound (#341)", {
  d <- zero_rows_data()
  left <- d
  left$l[left$l == 0] <- NA
  for (dist in c("exponential", "weibull", "lognormal", "loglogistic")) {
    open <- zero_rows_fit(d, dist)
    ref <- zero_rows_fit(left, dist)
    expect_true(isTRUE(ref$fit$converged), info = dist)
    expect_equal(open$fit$objective, ref$fit$objective, tolerance = 1e-8,
                 info = dist)
    expect_equal(open$fit$theta, ref$fit$theta, tolerance = 1e-6,
                 info = dist)
  }
})

test_that("the exponential matches a one-phase multiphase fit at time 0 (#341)", {
  d <- rbind(zero_rows_data(), data.frame(l = 0, u = NA))
  expo <- zero_rows_fit(d, "exponential")
  mp <- hazard(survival::Surv(l, u, type = "interval2") ~ 1, data = d,
               dist = "multiphase", phases = list(c = hzr_phase("constant")),
               fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  expect_equal(expo$fit$objective, mp$fit$objective, tolerance = 1e-8)
  expect_equal(unname(expo$fit$theta), unname(mp$fit$theta), tolerance = 1e-5)
})

test_that("a lognormal right-censored row at 0 needs no interval rows (#341)", {
  withr::local_seed(20260917)
  t <- stats::rlnorm(80, 1, 0.8)
  cens <- stats::runif(80, 1, 12)
  d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens))
  d0 <- rbind(d, data.frame(time = 0, status = 0))
  fit <- function(dd) {
    hazard(survival::Surv(time, status) ~ 1, data = dd, dist = "lognormal",
           theta = c(1, 0), fit = TRUE)
  }
  without <- fit(d)
  with <- fit(d0)
  expect_equal(with$fit$objective, without$fit$objective, tolerance = 1e-8)
  expect_equal(with$fit$theta, without$fit$theta, tolerance = 1e-6)
  # The analytic score and Hessian see the row too: no NaN, the same vcov.
  expect_true(all(is.finite(with$fit$vcov)))
  expect_equal(with$fit$vcov, without$fit$vcov, tolerance = 1e-6)
})

test_that("one time-0 row leaves every family's fit intact (#341)", {
  skip_on_cran()
  # #336's own event-at-0 control is exponential only, and the issue #341
  # names lognormal. The families do NOT behave alike here, so one family
  # cannot stand for the rest: before this fix, a single right-censored row
  # at time 0, or one (0, u] interval, drove lognormal to the optimizer's
  # -1e10 clamp -- a fit reported `converged` over no likelihood -- while
  # exponential, weibull and loglogistic were unaffected. Testing only the
  # families that already worked would have shown a clean sweep.
  t_ok <- c(0.5, 1.2, 2.3, 3.1, 4.7, 5.5)
  s_ok <- c(1, 1, 0, 1, 0, 1)
  thetas <- list(exponential = 0.1, weibull = c(0.1, 1),
                 lognormal = c(0.1, 1), loglogistic = c(0.1, 1))
  for (d in names(thetas)) {
    fit <- function(...) {
      suppressWarnings(hazard(dist = d, theta = thetas[[d]], fit = TRUE, ...))
    }
    base <- fit(time = t_ok, status = s_ok)
    cens0 <- fit(time = c(0, t_ok), status = c(0, s_ok))
    ivl0 <- fit(time = c(1, t_ok), status = c(2, s_ok),
                time_lower = c(0, t_ok), time_upper = c(1, t_ok))

    # Never the clamp. abs(obj) >= 1e10 is .hzr_optim_generic()'s sentinel
    # for a non-finite objective, and it arrives with converged = TRUE, so
    # a finiteness check alone would pass on exactly the broken case.
    for (f in list(base, cens0, ivl0)) {
      expect_true(is.finite(f$fit$objective), info = d)
      expect_lt(abs(f$fit$objective), 1e10)
    }

    # A right-censored row at time 0 contributes H(0) = 0: the objective is
    # not merely finite, it is UNCHANGED. This is the assertion that says
    # the row was handled rather than silently dropping the other six.
    expect_equal(cens0$fit$objective, base$fit$objective, info = d)

    # The (0, u] interval is real information, so it must MOVE the
    # objective. Asserting it is finite would also pass if the row were
    # dropped; asserting it differs is what proves it was used.
    expect_false(isTRUE(all.equal(ivl0$fit$objective, base$fit$objective)))
  }
})
