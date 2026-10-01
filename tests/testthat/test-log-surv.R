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
  # Case 2L at m = -3 and m = -0.4 as well as -1: at m = -1 the coefficient
  # mm1 of its hazard's first term is 0, so that term goes untested (Copilot).
  shapes <- list(c(1, 1), c(0.2931, 120), c(2, 0.5), c(1.5, 0), c(1, -0.5),
                 c(2, -5), c(0, -1), c(0, -3), c(0, -0.4), c(-1, 1),
                 c(-1.5, 0))
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

test_that("near t = 0, and where bt underflows, log_surv is 0, not NA (#578)", {
  # Found by review. Case 1L (m = 0) carried log(-log G) = -log(bt) / nu, and
  # where that overflowed -- nu < 1 at the time-0 clamp, or bt underflowing
  # to 0 for a large t_half -- log(1 - G) came back NA where it is 0: a phase
  # that has not started. An entry time of 0 is evaluated at every row with
  # none, so the likelihood read -Inf for every nu below about 1.
  tt <- c(0, 1e-300, 1e-12, 1e-6)
  shapes <- list(c(1, 1), c(0.3, 40), c(1.5, 0), c(0.5, 0), c(0.05, 0),
                 c(1, -0.5), c(2, -5), c(0, -1), c(-1, 1), c(-1.5, 0))
  for (sh in shapes) {
    d <- hzr_decompos(tt, t_half = 1, nu = sh[1], m = sh[2])
    label <- paste0("nu = ", sh[1], ", m = ", sh[2])
    expect_true(all(is.finite(d$log_surv)), label = paste(label, "finite"))
    expect_true(all(d$log_surv <= 0), label = paste(label, "not positive"))
    expect_true(all(is.finite(d$h)), label = paste(label, "hazard finite"))
  }
  # Case 1L with bt underflowing: the phase has not started at these times.
  d <- hzr_decompos(c(1e-3, 0.01, 0.1), t_half = 12, nu = 0.005, m = 0)
  expect_identical(d$log_surv, c(0, 0, 0))
  expect_true(all(is.finite(d$h)))
  expect_identical(hzr_decompos(c(0, 1e-300), t_half = 1, nu = 0.5,
                                m = 0)$log_surv, c(0, 0))
  # Case 2L where t / rho underflows to 0: a large t_half, or a very
  # negative m, at the time-0 entry (second review).
  expect_identical(hzr_decompos(0, t_half = 1e16, nu = 0, m = -1)$log_surv, 0)
  # At m = -50, t / rho is about 2e-325 at the time-0 clamp and underflows,
  # but G = (1 - e^(-t / rho))^(1/50) is still about 3e-7: log(1 - G) is
  # about -3e-7, not 0. Closed form from log(t / rho).
  rho <- 100 / -log1p(-2^-50)
  log_g_cdf <- (log(.Machine$double.xmin) - log(rho)) / 50
  want <- log1p(-exp(log_g_cdf))
  got <- hzr_decompos(0, t_half = 100, nu = 0, m = -50)$log_surv
  expect_lt(want, -1e-8)
  expect_equal(got / want, 1, tolerance = 1e-6)
  # And its hazard there (Copilot, #583). At m = -1 the hazard is exactly
  # 1 / rho, rho = t_half / log(2), at every time; t / rho underflows to 0
  # at these times, which made it NA.
  for (th in c(1e16, 1e300)) {
    h <- hzr_decompos(c(0, 1e-300, 1), t_half = th, nu = 0, m = -1)$h
    expect_equal(h / (log(2) / th), c(1, 1, 1), tolerance = 1e-10,
                 label = paste("t_half =", th))
  }
  # Case 2L at m = -1 has 1 - G = exp(-t / rho), rho = t_half / log(2): the
  # value must keep its relative accuracy as t goes to 0, not cancel away.
  # As a ratio: these values are far below any tolerance, and an absolute
  # comparison would pass anything.
  small <- c(1e-12, 1e-100, 1e-300)
  ratio <- hzr_decompos(small, t_half = 1, nu = 0, m = -1)$log_surv /
    (-small * log(2))
  expect_equal(ratio, c(1, 1, 1), tolerance = 1e-10)
})

test_that("the hazard keeps its value for a tiny nu (#578)", {
  skip_if_not_installed("numDeriv")
  # Found by Copilot. For a tiny |nu|, log(g) and log(1 - G) share a term of
  # size 1 / nu, so their difference lost the hazard: Case 3 at nu = -1e-18,
  # t = 2, t_half = 1, m = 1 gave h = 1 where log(h) is 40.753. Fitted |nu|
  # near 1e-16 occur (#448). Oracle: t * h is the derivative of -log_surv in
  # log(t), and log_surv keeps its accuracy there.
  expect_equal(log(hzr_decompos(2, 1, -1e-18, 1)$h), 40.75338, tolerance = 1e-6)
  # Case 2 saturated as well as tiny nu (Copilot): the exact hazard is
  # dm / (nu * (t_half + t * dm)), dm = (1 - 2^m)^(-nu) - 1.
  for (sh in list(c(1e-18, -1), c(1e-16, -0.5))) {
    th <- 1e-40
    dm <- expm1(-sh[1] * log1p(-2^sh[2]))
    want <- dm / (sh[1] * (th + 1 * dm))
    got <- hzr_decompos(1, t_half = th, nu = sh[1], m = sh[2])$h
    expect_equal(got / want, 1, tolerance = 1e-8,
                 label = paste("Case 2, nu =", sh[1], "m =", sh[2]))
  }
  tt <- c(0.5, 2, 5)
  shapes <- list(c(1e-16, 1), c(-1e-16, 1), c(1e-18, 1), c(-1e-18, 1),
                 c(1e-16, 0), c(1e-14, 2), c(-1e-14, 2), c(1e-12, -0.5))
  for (sh in shapes) {
    d <- hzr_decompos(tt, t_half = 1, nu = sh[1], m = sh[2])
    oracle <- vapply(tt, function(t) {
      numDeriv::grad(function(lt) {
        -hzr_decompos(exp(lt), t_half = 1, nu = sh[1], m = sh[2])$log_surv
      }, log(t))
    }, numeric(1))
    live <- oracle > 1
    expect_true(any(live), label = paste("nu =", sh[1], "m =", sh[2]))
    expect_equal(tt[live] * d$h[live] / oracle[live], rep(1, sum(live)),
                 tolerance = 1e-6,
                 label = paste("nu =", sh[1], "m =", sh[2]))
  }
})

test_that("a hazard phase with entry times and nu below 1 fits (#578)", {
  skip_on_cran() # a multiphase fit
  # Data from the phase's own model: m = 0, nu = 0.5, t_half = 2, with a
  # fifth of the rows entering late. With log_surv NA at the entry time 0,
  # every nu below about 1 read -Inf, and a fit stopped against that edge
  # at nu = 0.9998 reporting convergence (found by review). Before #578 the
  # same call reported a log-likelihood of +23295 at nu = 0.037: the clamp.
  withr::local_seed(2)
  n <- 400
  nu <- 0.5
  rho <- nu * 2 * log(2)^nu
  tt <- (-log(1 - runif(n)))^(-nu) * rho / nu
  cens <- runif(n, 0, 20)
  st <- as.integer(tt <= cens)
  tt <- pmin(tt, cens) + 1e-6
  entry <- ifelse(runif(n) < 0.2, runif(n) * 0.5 * tt, 0)
  d <- data.frame(start = entry, stop = tt, event = st)
  expect_true(any(d$start > 0) && any(d$start == 0))
  fit <- suppressWarnings(hazard(
    survival::Surv(start, stop, event) ~ 1, data = d, dist = "multiphase",
    phases = list(late = hzr_phase("hazard", t_half = 1, nu = 2, m = 0,
                                   fixed = "m")),
    fit = TRUE, control = list(n_starts = 1)
  ))
  expect_lt(fit$fit$objective, 0)
  expect_lte(fit$fit$rel_gradient, .Machine$double.eps^(1 / 3))
  nu_hat <- unname(fit$fit$theta[["late.nu"]])
  expect_gt(nu_hat, 0.3)
  expect_lt(nu_hat, 0.7)
  # The likelihood is defined on the far side of nu = 1 as well.
  p <- coef(fit)
  p[["late.nu"]] <- 0.4
  expect_true(is.finite(as.numeric(hzr_evaluate(fit, p)$logLik)))
})
