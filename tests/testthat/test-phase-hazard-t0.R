# The phase hazard at exactly t = 0 (#672).
#
# Start/stop scoring evaluates the hazard at every segment's start, and each
# patient's first segment starts at 0. For m = 0 and nu > 0 (Case 1L) g was
# formed as G * bt^num1 / rho, which at t = 0 is 0 * Inf: NaN, and the NaN
# took the whole multiphase hazard with it in predict(), with no message.
# The limit is 0 for these shapes, which t = 1e-12 shows.

test_that("Case 1L g and the phase hazard at t = 0 equal their limit", {
  for (nu in c(0.5, 1, 2, 5)) {
    d0 <- hzr_decompos(0, t_half = 1, nu = nu, m = 0)
    d1 <- hzr_decompos(1e-12, t_half = 1, nu = nu, m = 0)
    expect_identical(d0$g, 0, info = paste("nu =", nu))
    expect_lt(d1$g, 1e-10)
    h0 <- hzr_phase_hazard(0, t_half = 1, nu = nu, m = 0)
    expect_identical(h0, 0, info = paste("nu =", nu))
  }
  # Away from 0 the values are unchanged: the log form is the same quantity.
  tt <- c(1e-6, 0.1, 1, 10)
  d <- hzr_decompos(tt, t_half = 2, nu = 1.5, m = 0)
  rho <- 1.5 * 2 * log(2)^1.5
  bt <- 1.5 * tt / rho
  expect_equal(d$g, d$G * bt^(-1 / 1.5 - 1) / rho, tolerance = 1e-13)
})

test_that("predict(type = 'hazard') at time 0 is finite on an early nu > 0, m = 0 phase", {
  set.seed(1)
  n <- 200
  d <- data.frame(time = stats::rexp(n, 0.3) + 0.01)
  d$status <- stats::rbinom(n, 1, 0.7)
  ph <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0,
                               fixed = c("nu", "m")),
             late  = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = d,
                dist = "multiphase", phases = ph,
                theta = c(log(0.3), log(0.5), 1, 0, log(0.1)), fit = TRUE,
                control = list(n_starts = 1))
  nd <- data.frame(time = c(0, 1e-8, 0.1))
  h <- predict(fit, newdata = nd, type = "hazard")
  expect_true(all(is.finite(h)))
  # The early phase contributes 0 at t = 0, so the hazard is the constant
  # phase's alone, as at t = 1e-8.
  expect_equal(h[1], h[2], tolerance = 1e-12)
})

test_that("Weibull predict(se.fit = TRUE) at time 0 is finite (#619 item 22)", {
  set.seed(2)
  t <- stats::rweibull(150, shape = 1.4, scale = 3)
  fit <- hazard(time = t, status = rep(1L, 150), dist = "weibull",
                theta = c(mu = 0.3, nu = 1), fit = TRUE)
  p <- predict(fit, newdata = data.frame(time = c(0, 1)),
               type = "survival", se.fit = TRUE)
  # At t = 0 survival is 1 with no uncertainty. The limits stay NA there by
  # design: they are built on log(H), and H = 0.
  expect_identical(p$fit[1], 1)
  expect_identical(p$se.fit[1], 0)
  expect_true(all(is.na(c(p$lower[1], p$upper[1]))))
  expect_true(all(is.finite(unlist(p[2, ]))))
})
