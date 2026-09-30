# log(1 - G) for a phase far past saturation (#578).
#
# A "hazard" phase's cumulative hazard is -log(1 - G(t)). It was formed from
# 1 - G, and once G rounded to 1 that was clamped at double.xmin: the
# cumulative hazard came back as -log(double.xmin) = 708.396 and the hazard
# g / (1 - G) near 1e290, where both are of order 1 to 100. hzr_decompos()
# now carries log(1 - G) as `log_surv`, from each case's own log-scale terms.
#
# The oracles are closed forms that do not go through hzr_decompos()'s
# branches: for nu = 1, m = 1 (Case 1), 1 - G = e^a / (1 + e^a) with
# a = log(t_half / t), so log(1 - G) = a - log1pexp(a); and for Case 3L,
# 1 - G = exp(-(t / (t_half log(2)^nu))^(-1/nu)).

ls_times <- c(0.05, 0.5, 1, 2, 6)

test_that("log_surv is log(1 - G) wherever 1 - G can be formed (#578)", {
  # All six cases at ordinary scales, where the old form is accurate: the
  # known positive.
  shapes <- list(c(1, 1), c(0.3, 40), c(1.5, 0), c(1, -0.5), c(2, -5),
                 c(0, -1), c(-1, 1), c(-1.5, 0))
  for (sh in shapes) {
    for (th in c(0.2, 1, 5)) {
      d <- hzr_decompos(ls_times, t_half = th, nu = sh[1], m = sh[2])
      ok <- d$G < 1 - 1e-6
      expect_true(any(ok))
      expect_equal(d$log_surv[ok], log(1 - d$G[ok]), tolerance = 1e-9,
                   label = paste0("nu = ", sh[1], ", m = ", sh[2],
                                  ", t_half = ", th))
    }
  }
})

test_that("a saturated Case 1 phase has its true cumulative hazard (#578)", {
  for (lth in c(-5, -30, -40, -100, -300, -700)) {
    d <- hzr_decompos(ls_times, t_half = exp(lth), nu = 1, m = 1)
    a <- lth - log(ls_times)
    expect_equal(d$log_surv, a - hzr_log1pexp(a), tolerance = 1e-12,
                 label = paste("log_t_half =", lth))
    # The phase's cumulative hazard and hazard, as the model reads them. The
    # hazard of this shape is 1 / t once saturated.
    Phi <- hzr_phase_cumhaz(ls_times, t_half = exp(lth), nu = 1, m = 1,
                            type = "hazard")
    expect_equal(Phi, hzr_log1pexp(a) - a, tolerance = 1e-12)
    h <- hzr_phase_hazard(ls_times, t_half = exp(lth), nu = 1, m = 1,
                          type = "hazard")
    if (lth <= -30) expect_equal(h * ls_times, rep(1, 5), tolerance = 1e-8)
  }
  # Not the old clamp: -log(double.xmin) is 708.396.
  Phi <- hzr_phase_cumhaz(1, t_half = exp(-40), nu = 1, m = 1, type = "hazard")
  expect_equal(Phi, 40, tolerance = 1e-12)
})

test_that("a saturated Case 3L phase matches its closed form (#578)", {
  nu <- -1.5
  for (lth in c(-3, -20, -40)) {
    th <- exp(lth)
    d <- hzr_decompos(ls_times, t_half = th, nu = nu, m = 0)
    btnu <- (ls_times / (th * log(2)^nu))^(-1 / nu)
    expect_equal(d$log_surv, -btnu, tolerance = 1e-12,
                 label = paste("log_t_half =", lth))
  }
})

test_that("every case's cumulative hazard and hazard agree as it saturates (#578)", {
  # h is the time derivative of -log(1 - G). Differencing -log_surv in log(t)
  # and comparing with t * h checks the two against each other, over t_half
  # far below the times, in every case, including the two shapes where the
  # first version of this form returned NA.
  skip_if_not_installed("numDeriv")
  shapes <- list(c(1, 1), c(0.2931, 120), c(2, 0.5), c(1.5, 0), c(1, -0.5),
                 c(2, -5), c(0, -1), c(-1, 1), c(-1.5, 0))
  for (sh in shapes) {
    for (lth in c(-3, -30, -60, -300, -700)) {
      th <- exp(lth)
      d <- hzr_decompos(ls_times, t_half = th, nu = sh[1], m = sh[2])
      label <- paste0("nu = ", sh[1], ", m = ", sh[2], ", log_t_half = ", lth)
      expect_true(all(is.finite(d$log_surv)), label = paste(label, "finite"))
      expect_true(all(d$log_surv < 0), label = paste(label, "negative"))
      Phi_t <- vapply(ls_times, function(t) {
        numDeriv::grad(function(lt) {
          -hzr_decompos(exp(lt), t_half = th, nu = sh[1], m = sh[2])$log_surv
        }, log(t), method = "Richardson")
      }, numeric(1))
      expect_equal(ls_times * d$h, Phi_t, tolerance = 1e-6, label = label)
    }
  }
})

test_that("the log_t_half derivative of a saturated hazard phase is right (#578)", {
  # #574's residual: the derivative went noisy from t_half = exp(-30) and to
  # 0 from exp(-34), at nu = 1, m = 1, where it is -1.
  for (lth in c(-25, -30, -34, -40, -100)) {
    pd <- .hzr_phase_derivatives(c(0.5, 1, 2), t_half = exp(lth), nu = 1,
                                 m = 1, type = "hazard")
    expect_equal(pd$dPhi_dlog_thalf, rep(-1, 3), tolerance = 1e-6,
                 label = paste("log_t_half =", lth))
  }
})
