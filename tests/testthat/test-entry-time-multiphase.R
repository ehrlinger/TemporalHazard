# Multiphase entry time: `time_lower` on status 0/1 rows (issue #253) -------
#
# On a status 0 or 1 row, `time_lower` is the counting-process entry time only
# when 0 < time_lower < time. `time_lower == time` means no entry: that is the
# mixed-interval layout, where exact and right-censored rows carry
# time_lower = time and only status-2 rows carry a real interval lower bound.
# This is the rule dist = "weibull" has always applied. The multiphase path
# used to subtract H(time_lower) whenever `time_lower` was supplied, so
# time_lower = time made every row enter as it left: on avc the fit reported
# converged with a log-likelihood of +47915.76 against -211.4677 without it.

.mp_case <- function() {
  d <- stats::na.omit(avc)
  ph <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0),
             constant = hzr_phase("constant"))
  list(
    time = d$int_dead, status = d$dead, phases = ph,
    phases_v = .hzr_validate_phases(ph),
    theta = c(log(0.35), log(0.19), 1.44, 1, log(0.03)),
    counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL)
  )
}

.mp_logl <- function(k, theta, time_lower, status = k$status,
                     time_upper = NULL) {
  .hzr_logl_multiphase(theta, time = k$time, status = status,
                       time_lower = time_lower, time_upper = time_upper,
                       phases = k$phases_v, covariate_counts = k$counts,
                       x_list = k$x_list)
}

# Some rows genuine entries, some time_lower = time, the rest 0.
.mp_mixed <- function(k) {
  n <- length(k$time)
  tl <- rep(0, n)
  genuine <- seq(1, n, by = 3)
  at_exit <- seq(2, n, by = 3)
  tl[genuine] <- k$time[genuine] / 2
  tl[at_exit] <- k$time[at_exit]
  ref <- tl
  ref[at_exit] <- 0
  # The layout must actually contain all three kinds, or the test compares
  # nothing.
  stopifnot(any(tl > 0 & tl < k$time), any(tl == k$time & tl > 0),
            any(tl == 0))
  list(tl = tl, ref = ref)
}

test_that("(a) a fit with time_lower = time equals the fit with no time_lower", {
  k <- .mp_case()
  f0 <- .hzr_optim_multiphase(k$time, k$status, phases = k$phases)
  f1 <- .hzr_optim_multiphase(k$time, k$status, time_lower = k$time,
                              phases = k$phases)
  expect_lt(f1$value, 0)
  expect_equal(f1$value, f0$value, tolerance = 1e-6)
  expect_equal(unname(f1$par), unname(f0$par), tolerance = 1e-4)
  expect_equal(f0$value, -211.4677, tolerance = 1e-6)
})

test_that("(b) at fixed theta, time_lower = time gives the NULL logl exactly", {
  k <- .mp_case()
  ll_null <- .mp_logl(k, k$theta, NULL)
  expect_true(is.finite(ll_null))
  expect_identical(.mp_logl(k, k$theta, k$time), ll_null)
})

test_that("(c) mixed layout: time_lower = time rows read as entry 0", {
  k <- .mp_case()
  m <- .mp_mixed(k)
  ll_ref <- .mp_logl(k, k$theta, m$ref)
  expect_true(is.finite(ll_ref))
  expect_identical(.mp_logl(k, k$theta, m$tl), ll_ref)
  # And the genuine entries are still entries: the reference is not the
  # no-entry logl.
  expect_false(isTRUE(all.equal(ll_ref, .mp_logl(k, k$theta, NULL))))

  # The CoE adjustment conserves on the same entry-time scale. The constant
  # phase is the one that can absorb the discrepancy at this theta; with the
  # early phase CoE declines and returns theta unchanged, which would compare
  # nothing.
  fix_pos <- 5L
  total <- sum(k$status == 1)
  coe <- function(tl) {
    .hzr_conserve_events(k$theta, "constant", fix_pos, k$time, k$status,
                         k$phases_v, k$counts, k$x_list, total,
                         time_lower = tl)
  }
  adjusted <- coe(m$ref)
  expect_false(isTRUE(all.equal(adjusted[fix_pos], k$theta[fix_pos])))
  expect_identical(coe(m$tl), adjusted)
})

test_that("(d) a status-2 row carries no entry subtraction", {
  k <- .mp_case()
  th <- k$theta
  status <- c(1, 0, 2)
  time <- c(1.0, 2.0, 1.5)
  lower <- c(0.4, 0.5, 0.8)
  upper <- c(1.0, 2.0, 1.5)
  kk <- list(time = time, phases_v = k$phases_v, counts = k$counts,
             x_list = k$x_list)
  ll <- .mp_logl(kk, th, lower, status = status, time_upper = upper)
  H <- function(t) {
    .hzr_multiphase_cumhaz(t, th, k$phases_v, k$counts, k$x_list)
  }
  h1 <- .hzr_multiphase_hazard(1.0, th, k$phases_v, k$counts, k$x_list)
  hand <- (log(h1) - (H(1.0) - H(0.4))) - (H(2.0) - H(0.5)) +
    log(exp(-H(0.8)) - exp(-H(1.5)))
  expect_equal(ll, hand, tolerance = 1e-12)
})

test_that("(e) gradient and Hessian match numDeriv on the mixed layout", {
  skip_if_not_installed("numDeriv")
  k <- .mp_case()
  m <- .mp_mixed(k)
  for (nm in c("mixed", "at exit")) {
    tl <- if (nm == "mixed") m$tl else k$time
    f <- function(p) .mp_logl(k, p, tl)
    g <- .hzr_gradient_multiphase(k$theta, time = k$time, status = k$status,
                                  time_lower = tl, phases = k$phases_v,
                                  covariate_counts = k$counts,
                                  x_list = k$x_list)
    expect_equal(g, numDeriv::grad(f, k$theta), tolerance = 1e-5, info = nm)
    hs <- .hzr_hessian_multiphase(k$theta, time = k$time, status = k$status,
                                  time_lower = tl, phases = k$phases_v,
                                  covariate_counts = k$counts,
                                  x_list = k$x_list)
    expect_false(is.null(hs), info = nm)
    expect_equal(unname(hs), numDeriv::hessian(function(p) -f(p), k$theta),
                 tolerance = 1e-4, info = nm)
    # The analytic derivatives also agree with the reference layout's.
    if (nm == "mixed") {
      expect_equal(g, .hzr_gradient_multiphase(
        k$theta, time = k$time, status = k$status, time_lower = m$ref,
        phases = k$phases_v, covariate_counts = k$counts, x_list = k$x_list),
        tolerance = 1e-12)
    }
  }
})

test_that("(f) genuine left truncation is unchanged", {
  # Values computed on the code before the #253 fix, where every entry was
  # already 0 < time_lower < time and so read as an entry under both rules.
  set.seed(1)
  n <- 60
  time <- stats::rexp(n, 0.5) + 0.2
  status <- stats::rbinom(n, 1, 0.7)
  entry <- time * stats::runif(n, 0.1, 0.8)
  stopifnot(all(entry > 0 & entry < time))
  ph <- list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0),
             constant = hzr_phase("constant"))
  k <- list(time = time, phases_v = .hzr_validate_phases(ph),
            counts = c(early = 0L, constant = 0L),
            x_list = list(early = NULL, constant = NULL))
  th <- c(log(0.3), log(1), 1, 0, log(0.1))

  expect_equal(.mp_logl(k, th, entry, status = status),
               -95.187915824809963, tolerance = 1e-12)
  expect_equal(
    .hzr_gradient_multiphase(th, time, status, time_lower = entry,
                             phases = k$phases_v, covariate_counts = k$counts,
                             x_list = k$x_list),
    c(12.729520416962995, 2.8494490333870264, -9.2012662706668422,
      -7.3753385293653607e-10, 23.589652067548901),
    tolerance = 1e-10)
  hs <- .hzr_hessian_multiphase(th, time, status, time_lower = entry,
                                phases = k$phases_v,
                                covariate_counts = k$counts,
                                x_list = k$x_list)
  expect_equal(unname(hs[c(1, 2, 3, 5), c(1, 2, 3, 5)]), matrix(c(
    -4.6719823126998747, -2.7924372244408859, 3.2186923025249516,
    8.817994044279871,
    -2.7924372244408859, 9.0112847085199377, 0.31747818316771581,
    2.9159313448562774,
    3.2186923025249516, 0.3174781831677157, -3.5861506233469047,
    -5.6510355499695004,
    8.817994044279871, 2.9159313448562765, -5.6510355499694995,
    -2.2831782603717627), 4, 4), tolerance = 1e-10)
  expect_equal(hs[4, 4], 64807.040011232821, tolerance = 1e-10)

  ft <- suppressWarnings(
    .hzr_optim_multiphase(time, status, time_lower = entry, phases = ph))
  # Optimizer end points are pinned loosely: BFGS can stop at slightly
  # different points on the five CI platforms. The fixed-theta pins above
  # (logl, gradient, Hessian) carry the tight checks.
  expect_equal(ft$value, -61.048341005835724, tolerance = 1e-6)
  # log_mu, log_t_half, nu and log_mu: m sits at ~6e-10, where a relative
  # comparison is meaningless, so it is checked on its own absolute scale.
  expect_equal(unname(ft$par[c(1, 2, 3, 5)]),
               c(-0.87430560608249352, -0.42697704500116312,
                 0.40509325865065615, -0.50755433855562759),
               tolerance = 1e-6)
  expect_lt(abs(ft$par[4]), 1e-6)
})
