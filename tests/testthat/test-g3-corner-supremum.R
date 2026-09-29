# A g3 phase under constraint = "eta_gamma" (FIXGE2) tends, as gamma grows,
# to a corner law: (t/tau)^2 below tau and (t/tau)^(2/alpha) above. The
# likelihood has a finite limit there, and when the data prefer a sharp bend
# the supremum is at gamma = Inf. A fit can stop at a finite gamma, report
# converged = TRUE, a plausible gamma-hat and SE, and warn about nothing
# (#418). The check searches the corner law and fires only on a CERTIFICATE:
# the package's own likelihood, at a finite large gamma, above the fit's.

# Survival times drawn from the exact gamma = Inf corner law, inverted in
# closed form: H = mu (t/tau)^2 below tau, mu (t/tau)^(2/alpha) above.
corner_data <- function(n, seed, mu = 0.7, tau = 1, alpha = 0.5, fu = 20) {
  withr::local_seed(seed)
  e <- stats::rexp(n)
  t <- ifelse(e <= mu, tau * sqrt(e / mu), tau * (e / mu)^(alpha / 2))
  list(time = pmin(t, fu * tau), status = as.numeric(t <= fu * tau))
}

fit_fixge2 <- function(d, tau = 1.5, gamma = 4) {
  hazard(time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
         phases = list(late = hzr_phase("g3", tau = tau, gamma = gamma,
                                        alpha = 0.5, fixed = "alpha",
                                        constraint = "eta_gamma")),
         control = list(n_starts = 1L))
}

corner_records <- function(fit) {
  b <- fit$fit$boundary
  if (!is.list(b)) return(list())
  Filter(function(r) identical(r$mechanism, "g3_corner_supremum"), b)
}

test_that("a gamma-hat short of the corner supremum is recorded and warned (#418)", {
  skip_on_cran()
  d <- corner_data(300, seed = 16)
  w <- list()
  fit <- withCallingHandlers(
    fit_fixge2(d),
    warning = function(cnd) {
      w[[length(w) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }
  )
  # The premise: an ordinary-looking interior fit, which is what makes it
  # dangerous.
  expect_true(fit$fit$converged)
  expect_lt(unname(coef(fit)["late.gamma"]), 1e3)

  rec <- corner_records(fit)
  expect_length(rec, 1L)
  rec <- rec[[1L]]
  expect_identical(rec$phase, "late")
  expect_identical(rec$parameter, "gamma")
  expect_gt(rec$gain, 0.01)
  # The certificate is the package's own likelihood, recomputed here from the
  # recorded point: the check cannot fire on its own search's arithmetic.
  again <- as.numeric(hzr_evaluate(fit, rec$certificate_theta)$logLik)
  expect_equal(again, rec$certificate_loglik, tolerance = 1e-10)
  expect_gt(again, fit$fit$objective + 0.01)
  expect_gt(rec$certificate_theta[["late.gamma"]], 1e3)

  # One warning of the family carries it, classed by mechanism.
  cls <- vapply(w, function(cnd) inherits(cnd, "hzr_g3_corner_supremum"),
                logical(1))
  expect_identical(sum(cls), 1L)
  expect_true(inherits(w[[which(cls)]], "hzr_boundary"))
  expect_match(conditionMessage(w[[which(cls)]]), "gamma", fixed = TRUE)

  # Report only: the estimates are the optimizer's, and the objective is the
  # likelihood at them.
  expect_equal(as.numeric(hzr_evaluate(fit, coef(fit))$logLik),
               fit$fit$objective, tolerance = 1e-8)
})

test_that("a fit at its maximum carries no corner record (#418)", {
  skip_on_cran()
  # The phase-constraint-oracle fixture's eta_gamma data: gamma* = 2, n = 800.
  withr::local_seed(325)
  u <- (stats::rexp(800) / 0.5)^(1 / 1)
  event <- 8 * expm1(0.5 * log1p(u))^(1 / 2)
  censor <- stats::runif(800, 1, 20)
  d <- list(time = pmin(event, censor), status = as.numeric(event <= censor))
  fit <- expect_no_warning(fit_fixge2(d, tau = 6, gamma = 3),
                           class = "hzr_g3_corner_supremum")
  expect_length(corner_records(fit), 0L)
  expect_null(fit$fit$boundary)
})

test_that("a gamma-hat already past 1e3 is left to the loud checks (#418)", {
  skip_on_cran()
  fit <- suppressWarnings(fit_fixge2(corner_data(300, seed = 3)))
  expect_gt(unname(coef(fit)["late.gamma"]), 1e3)
  expect_length(corner_records(fit), 0L)
})

test_that("the certificate is taken off the bend as well as on it (#418)", {
  skip_on_cran()
  # Here tau lands on an observed time, where the finite-gamma hazard blends
  # the two corner slopes; at tau itself the certificate falls short by
  # log(1.5), and only the +-10 / gamma nudge finds the higher likelihood.
  fit <- suppressWarnings(fit_fixge2(corner_data(300, seed = 7)))
  expect_lt(unname(coef(fit)["late.gamma"]), 1e3)
  expect_length(corner_records(fit), 1L)
})

test_that("only the certificate can fire, and only past its margin (#418)", {
  skip_on_cran()
  d <- corner_data(300, seed = 16)
  fit <- suppressWarnings(fit_fixge2(d))
  theta <- coef(fit)
  value <- fit$fit$objective
  record <- function(objective_fn, th = theta) {
    .hzr_g3_corner_record(
      th, value, objective_fn, 1L, .hzr_validate_phases(fit$spec$phases),
      c(late = 0L), list(late = NULL), d$time, d$status, NULL,
      rep(1, length(d$time)))
  }
  # Known positive: the real objective certifies.
  real <- function(p) as.numeric(hzr_evaluate(fit, p)$logLik)
  expect_false(is.null(record(real)))
  # The corner search finds the same higher region either way; what decides
  # is the objective handed in. Below the fit, or within 0.01 of it: nothing.
  expect_null(record(function(p) value - 1))
  expect_null(record(function(p) value + 0.005))
  expect_false(is.null(record(function(p) value + 0.02)))
  # A gamma-hat at 1e3 or beyond is not examined, whatever the objective says.
  big <- theta
  big[["late.gamma"]] <- 1e3
  expect_null(record(function(p) value + 1, th = big))
  # Nor is a fit the optimizer did not report as converged: its gamma-hat is
  # not a claimed optimum, so there is no claim to contradict.
  sweep <- function(converged) {
    .hzr_g3_corner_supremum(
      theta, value, converged, function(p) value + 1,
      .hzr_validate_phases(fit$spec$phases), c(late = 0L), list(late = NULL),
      d$time, d$status, NULL, rep(1, length(d$time)))
  }
  expect_length(sweep(TRUE), 1L)
  expect_length(sweep(FALSE), 0L)
})

test_that("the certificate moves only what the fit estimated (#418)", {
  skip_on_cran()
  # The same data as the known positive, which records a finding with gamma,
  # tau and alpha as in fit_fixge2().
  d <- corner_data(300, seed = 16)
  fit_with <- function(fixed) {
    suppressWarnings(hazard(
      time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
      phases = list(late = hzr_phase("g3", tau = 1.5, gamma = 4, alpha = 0.5,
                                     fixed = fixed,
                                     constraint = "eta_gamma")),
      control = list(n_starts = 1L)))
  }
  # A fixed gamma was never estimated: there is no gamma-hat to contradict,
  # and a certificate at gamma = 1e4 would score another model.
  expect_length(corner_records(fit_with(c("alpha", "gamma"))), 0L)
  # A fixed tau stays where it was fixed, in the search and the certificate.
  # Take the known positive's estimates with tau moved to where its own
  # certificate put the bend, and declare tau fixed there: the corner law
  # still beats that point, and the certificate must keep the tau.
  fit <- fit_with("alpha")
  free_rec <- corner_records(fit)[[1L]]
  th <- coef(fit)
  # Just off the bend, so a search free to move tau would move it.
  th[["late.log_tau"]] <- free_rec$certificate_theta[["late.log_tau"]] + 0.003
  ll <- function(p) as.numeric(hzr_evaluate(fit, p)$logLik)
  held_tau <- .hzr_validate_phases(list(
    late = hzr_phase("g3", tau = 1.5, gamma = 4, alpha = 0.5,
                     fixed = c("alpha", "tau"), constraint = "eta_gamma")))
  rec <- .hzr_g3_corner_record(
    th, ll(th), ll, 1L, held_tau,
    c(late = 0L), list(late = NULL), d$time, d$status, NULL,
    rep(1, length(d$time)))
  expect_false(is.null(rec))
  expect_identical(rec$certificate_theta[["late.log_tau"]],
                   th[["late.log_tau"]])
  expect_identical(rec$certificate_theta[["late.alpha"]], 0.5)
})

test_that("only eta_gamma phases are examined (#418)", {
  skip_on_cran()
  # Paired with the known positive on the same data, which does record one.
  d <- corner_data(300, seed = 16)
  fit <- suppressWarnings(hazard(
    time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
    phases = list(late = hzr_phase("g3", tau = 1.5, gamma = 4, alpha = 0.5,
                                   eta = 0.5, fixed = c("alpha", "eta"))),
    control = list(n_starts = 1L)))
  expect_length(corner_records(fit), 0L)
})
