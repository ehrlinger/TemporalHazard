# test-phase-constraint-oracle.R
# An answer for hzr_phase(constraint = ) that owes nothing to the package
# (#325). The late-phase G3 cumulative hazard, its derivative, the simulation
# and the constrained optimizer are all written out below from the closed
# form, so agreement is between two implementations rather than one
# implementation and a restatement of it.
#
#   H(t) = mu * (((t / tau)^gamma + 1)^(1 / alpha) - 1)^eta

oracle_cumhaz <- function(t, mu, tau, gamma, alpha, eta) {
  s <- (t / tau)^gamma
  mu * expm1(log1p(s) / alpha)^eta
}

oracle_hazard <- function(t, mu, tau, gamma, alpha, eta) {
  s <- (t / tau)^gamma
  u <- expm1(log1p(s) / alpha)
  mu * eta * u^(eta - 1) * exp((1 / alpha - 1) * log1p(s)) / alpha *
    gamma * s / t
}

# Exact inversion of H: u = (E / mu)^(1 / eta), s = (1 + u)^alpha - 1,
# t = tau * s^(1 / gamma), for E ~ Exp(1).
oracle_data <- function(mu, tau, gamma, alpha, eta, n) {
  withr::local_seed(325)
  u <- (stats::rexp(n) / mu)^(1 / eta)
  event <- tau * expm1(alpha * log1p(u))^(1 / gamma)
  censor <- stats::runif(n, 1, 20)
  list(time = pmin(event, censor), status = as.numeric(event <= censor))
}

oracle_loglik <- function(d, mu, tau, gamma, alpha, eta) {
  ev <- d$status == 1
  sum(log(oracle_hazard(d$time[ev], mu, tau, gamma, alpha, eta))) -
    sum(oracle_cumhaz(d$time, mu, tau, gamma, alpha, eta))
}

oracle_cases <- list(
  # FIXGAE2: search log mu, log tau, log gamma, log eta; alpha = gamma*eta/2.
  alpha_gamma_eta = list(
    truth = list(mu = 0.3, tau = 10, gamma = 6, eta = 0.8),
    start = list(mu = 0.2, tau = 8, gamma = 5, eta = 1),
    shapes = function(p) {
      gamma <- exp(p[3])
      eta <- exp(p[4])
      list(mu = exp(p[1]), tau = exp(p[2]), gamma = gamma,
           alpha = gamma * eta / 2, eta = eta)
    },
    phase = function(s) {
      hzr_phase("g3", tau = s$tau, gamma = s$gamma, eta = s$eta,
                constraint = "alpha_gamma_eta")
    }
  ),
  # FIXGE2: search log mu, log tau, log gamma; eta = 2/gamma. alpha is held,
  # because with alpha free as well a sample this size runs gamma to a ridge.
  # Held below 1: PROC HAZARD refuses a fixed alpha above 1 under FIXGE2
  # (SETG31040, #418), and so does hazard().
  eta_gamma = list(
    truth = list(mu = 0.5, tau = 8, gamma = 2, eta = 1),
    start = list(mu = 0.3, tau = 6, gamma = 3, eta = 2 / 3),
    alpha = 0.5,
    shapes = function(p) {
      gamma <- exp(p[3])
      list(mu = exp(p[1]), tau = exp(p[2]), gamma = gamma, alpha = 0.5,
           eta = 2 / gamma)
    },
    phase = function(s) {
      hzr_phase("g3", tau = s$tau, gamma = s$gamma, alpha = 0.5,
                fixed = "alpha", constraint = "eta_gamma")
    }
  )
)

for (constraint in names(oracle_cases)) {
  test_that(paste("hazard() under", constraint, "matches an independent fit"), {
    case <- oracle_cases[[constraint]]
    truth <- case$truth
    alpha_true <- if (constraint == "alpha_gamma_eta") {
      truth$gamma * truth$eta / 2
    } else {
      case$alpha
    }
    d <- oracle_data(truth$mu, truth$tau, truth$gamma, alpha_true, truth$eta,
                     n = 800)

    s0 <- case$start
    p0 <- log(c(s0$mu, s0$tau, s0$gamma,
                if (constraint == "alpha_gamma_eta") s0$eta))
    nll <- function(p) {
      v <- -do.call(oracle_loglik, c(list(d), case$shapes(p)))
      if (is.finite(v)) v else 1e10
    }
    ref <- stats::nlm(nll, p0, gradtol = 1e-10, steptol = 1e-12,
                      iterlim = 2000, hessian = TRUE)
    # Code 3 is nlm stopping on step size along a flat direction; the
    # log-likelihood comparison below is what holds it to the optimum.
    expect_true(ref$code %in% 1:3)
    ref_shapes <- case$shapes(ref$estimate)

    fit <- suppressWarnings(hazard(
      time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
      phases = list(late = case$phase(s0)),
      theta = c(log(s0$mu), log(s0$tau), s0$gamma,
                if (constraint == "alpha_gamma_eta") s0$gamma * s0$eta / 2
                else case$alpha,
                s0$eta),
      control = list(conserve = FALSE, n_starts = 1L)
    ))
    theta <- unname(coef(fit))
    got <- list(mu = exp(theta[1]), tau = exp(theta[2]), gamma = theta[3],
                alpha = theta[4], eta = theta[5])

    expect_true(fit$fit$converged)
    expect_equal(fit$fit$objective, -ref$minimum, tolerance = 1e-8)
    # Estimates are held exactly as tightly as the likelihood above, no
    # tighter: that assertion admits a shortfall of up to
    # dll = 1e-8 * |logLik| (relative tolerance), and near the optimum
    # the largest change in log(element i) that costs no more than dll is
    # sqrt(2 * dll * (J H^-1 J')_ii), with H the oracle's Hessian and J the
    # map from search to log-element scale (linear here, so the unit-step
    # difference is exact). Along the flat direction that is a few 1e-3;
    # across it, far less. A fixed element gets a bound of 0.
    dll <- 1e-8 * abs(ref$minimum)
    log_shapes <- function(p) log(unlist(case$shapes(p)))
    jac <- sapply(seq_along(ref$estimate), function(k) {
      log_shapes(ref$estimate + replace(0 * ref$estimate, k, 1)) -
        log_shapes(ref$estimate)
    })
    bound <- sqrt(2 * dll * diag(jac %*% solve(ref$hessian) %*% t(jac)))
    ratio <- unlist(got) / unlist(ref_shapes)
    # Report the worst excess, not just FALSE, when this fails.
    expect_lte(max(abs(log(ratio)) - bound), 1e-12)
    # The package's objective is the oracle's likelihood at the package's own
    # estimates, not only at the oracle's.
    expect_equal(fit$fit$objective, do.call(oracle_loglik, c(list(d), got)),
                 tolerance = 1e-10)

    # Cumulative hazard, not survival: survival sits near 1 at early times,
    # where a difference in H barely moves it. Worst element, as a ratio;
    # measured up to 2.2e-4, with H as small as 2e-6 at t = 1.
    grid <- c(1, 3, 5, 10, 15)
    cumhaz <- predict(fit, newdata = data.frame(time = grid),
                      type = "cumulative_hazard")
    ref_cumhaz <- do.call(oracle_cumhaz, c(list(grid), ref_shapes))
    expect_true(all(ref_cumhaz > 0))
    expect_lt(max(abs(as.numeric(cumhaz) / ref_cumhaz - 1)), 1e-3)
  })
}
