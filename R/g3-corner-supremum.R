# A g3 phase heading for gamma = Inf under FIXGE2 (#418)
#
# Under constraint = "eta_gamma" (eta = 2 / gamma), G3 tends as gamma grows
# to a corner law: (t / tau)^2 below tau and (t / tau)^(2 / alpha) above.
# The log-likelihood has a finite limit there, and when the data prefer a
# sharp bend the supremum is at gamma = Inf and no finite maximum exists. A
# fit can stop short of it with converged = TRUE, an ordinary gamma-hat and
# SE, and a well-conditioned Hessian, so nothing warned (#418).
#
# The check searches the corner law with every other phase held at the fit,
# and fires only on a CERTIFICATE: the package's own objective, at a finite
# large gamma, above the fit's. A search that under-optimises (the first
# comparator for this did, in 180 of 250 rows) can only miss a case; it
# cannot report one that is not there. PROC HAZARD has no such check: UMSTOP
# accepts on the relative gradient or the step size alone (umstop.c:173-203).

#' Corner-law log shape and log hazard shape of a FIXGE2 g3 phase
#' @param t Times, positive.
#' @param tau,alpha The phase's tau and alpha.
#' @return List of `log_G` and `log_g` (log dG/dt).
#' @keywords internal
#' @noRd
.hzr_g3_corner_law <- function(t, tau, alpha) {
  below <- t <= tau
  r <- log(t / tau)
  power <- ifelse(below, 2, 2 / alpha)
  log_G <- power * r
  list(log_G = log_G, log_g = log_G + log(power) - log(t))
}

#' Search the gamma = Inf corner law of one FIXGE2 g3 phase and certify
#'
#' @param theta The fitted, named full theta.
#' @param value The fit's objective at `theta`.
#' @param objective_fn Function of a full theta returning the same objective
#'   the fit maximised.
#' @param k Index of the phase in `phases`.
#' @param phases,covariate_counts,x_list The fit's phase specification.
#' @param time,status,time_lower,weights The rows the likelihood read.
#' @param delta The gain the certificate must exceed.
#' @return A boundary record, or `NULL`.
#' @keywords internal
#' @noRd
.hzr_g3_corner_record <- function(theta, value, objective_fn, k, phases,
                                  covariate_counts, x_list, time, status,
                                  time_lower, weights, delta = 0.01) {
  nm <- names(phases)[[k]]
  key <- function(p) paste0(nm, ".", p)
  if (!all(key(c("log_mu", "log_tau", "gamma", "alpha", "eta")) %in%
             names(theta))) {
    return(NULL)
  }
  gamma_hat <- unname(theta[[key("gamma")]])
  alpha <- unname(theta[[key("alpha")]])
  if (!is.finite(gamma_hat) || gamma_hat >= 1e3 || !is.finite(alpha) ||
        alpha <= 0) {
    return(NULL)
  }
  # The certificate may move only what the fit estimated. A fixed gamma was
  # never estimated, so there is no gamma-hat to contradict; a fixed tau or
  # alpha stays at its value, or the certificate would score another model.
  fixed <- phases[[k]]$fixed
  if (any(c("gamma", "shapes") %in% fixed)) return(NULL)
  alpha_free <- !any(c("alpha", "shapes") %in% fixed)
  tau_free <- !("tau" %in% fixed)
  w <- if (is.null(weights)) rep(1, length(time)) else weights
  entry <- .hzr_multiphase_entry(time, status, time_lower)
  if (is.null(entry)) entry <- rep(0, length(time))
  ev <- status == 1 & w > 0

  # Every other phase at the fit: this phase's intercept sent to exp(-745).
  split <- .hzr_split_theta(theta, phases, covariate_counts)
  pars <- .hzr_unpack_phase_theta(split[[nm]], phases[[nm]])
  offset <- if (length(pars$beta) && !is.null(x_list[[nm]])) {
    as.numeric(x_list[[nm]] %*% pars$beta)
  } else {
    rep(0, length(time))
  }
  off_theta <- theta
  off_theta[[key("log_mu")]] <- -745
  h_other <- .hzr_multiphase_hazard(time, off_theta, phases,
                                    covariate_counts, x_list)
  cum_other <- .hzr_multiphase_cumhaz(time, off_theta, phases,
                                      covariate_counts, x_list)
  has_entry <- entry > 0
  if (any(has_entry)) {
    cum_other[has_entry] <- cum_other[has_entry] - .hzr_multiphase_cumhaz(
      entry[has_entry], off_theta, phases, covariate_counts,
      lapply(x_list, function(m) {
        if (is.null(m)) NULL else m[has_entry, , drop = FALSE]
      }))
  }
  if (!all(is.finite(h_other[ev])) || !all(is.finite(cum_other))) {
    return(NULL)
  }

  corner_ll <- function(log_mu, tau, a) {
    shape <- .hzr_g3_corner_law(time, tau, a)
    cum <- exp(log_mu + offset + shape$log_G)
    if (any(has_entry)) {
      cum[has_entry] <- cum[has_entry] -
        exp(log_mu + offset[has_entry] +
              .hzr_g3_corner_law(entry[has_entry], tau, a)$log_G)
    }
    haz <- h_other[ev] + exp(log_mu + offset[ev] + shape$log_g[ev])
    v <- sum(w[ev] * log(haz)) - sum(w * (cum_other + cum))
    if (is.finite(v)) v else -Inf
  }
  best_log_mu <- function(tau, a) {
    o <- stats::optimize(function(m) -corner_ll(m, tau, a),
                         pars$log_mu + c(-30, 30))
    c(log_mu = o$minimum, ll = -o$objective)
  }

  # tau over the observed times near tau-hat and their midpoints: under the
  # corner law the tau profile is smooth between observed times, which is
  # what a finite-gamma search lacks.
  tau_hat <- exp(unname(theta[[key("log_tau")]]))
  obs <- sort(unique(time[w > 0]))
  cand <- if (tau_free) {
    obs[obs > tau_hat * exp(-4) & obs < tau_hat * exp(4)]
  } else {
    numeric(0)
  }
  if (length(cand) > 100L) {
    cand <- stats::quantile(cand, seq(0, 1, length.out = 100L), names = FALSE)
  }
  cand <- sort(unique(c(cand, tau_hat)))
  if (length(cand) > 1L) {
    cand <- sort(c(cand, (utils::head(cand, -1L) + utils::tail(cand, -1L)) / 2))
  }
  grid <- vapply(cand, function(tt) best_log_mu(tt, alpha), numeric(2))
  top <- utils::head(order(-grid["ll", ]), 5L)
  polished <- lapply(top, function(j) {
    # Search vector: log_mu, then log tau and logit alpha where free.
    unpack <- function(p) {
      list(log_mu = p[1],
           log_tau = if (tau_free) p[2] else unname(theta[[key("log_tau")]]),
           alpha = if (alpha_free) stats::plogis(p[length(p)]) else alpha)
    }
    start <- c(grid["log_mu", j], if (tau_free) log(cand[[j]]),
               if (alpha_free) stats::qlogis(min(max(alpha, 1e-6), 1 - 1e-6)))
    f <- function(p) {
      u <- unpack(p)
      -corner_ll(u$log_mu, exp(u$log_tau), u$alpha)
    }
    o <- if (length(start) > 1L) {
      stats::optim(start, f, method = "Nelder-Mead",
                   control = list(reltol = 1e-12, maxit = 2000))
    } else {
      list(par = start, value = f(start))
    }
    c(unpack(o$par), ll = -o$value)
  })
  best <- polished[[which.max(vapply(polished, `[[`, numeric(1), "ll"))]]

  # The certificate: the fit's own objective at finite gamma, with tau also
  # nudged by +-10 / gamma. At an observed time lying on tau the finite-gamma
  # hazard is a blend of the two corner slopes (8/3 against 4 at alpha = 0.5),
  # which a certificate taken at tau itself misses by log(1.5).
  certs <- list()
  for (g in c(1e4, 1e6, 1e8)) {
    for (shift in if (tau_free) c(-10, 0, 10) else 0) {
      p <- theta
      p[[key("log_mu")]] <- best$log_mu
      p[[key("log_tau")]] <- best$log_tau + shift / g
      p[[key("gamma")]] <- g
      p[[key("eta")]] <- 2 / g
      p[[key("alpha")]] <- best$alpha
      v <- tryCatch(objective_fn(p), error = function(e) NA_real_)
      if (length(v) == 1L && is.finite(v)) {
        certs[[length(certs) + 1L]] <- list(theta = p, ll = v)
      }
    }
  }
  if (!length(certs)) return(NULL)
  cert <- certs[[which.max(vapply(certs, `[[`, numeric(1), "ll"))]]
  gain <- cert$ll - value
  if (!(gain > delta)) return(NULL)

  list(
    mechanism = "g3_corner_supremum",
    phase = nm,
    parameter = "gamma",
    gain = gain,
    gamma_hat = gamma_hat,
    certificate_theta = cert$theta,
    certificate_loglik = cert$ll,
    detail = paste0(
      "phase '", nm, "' (constraint = \"eta_gamma\") stopped at gamma = ",
      format(gamma_hat, digits = 4), ", but the log-likelihood is ",
      format(gain, digits = 3), " higher at gamma = ",
      format(cert$theta[[key("gamma")]], digits = 3), " (tau = ",
      format(exp(cert$theta[[key("log_tau")]]), digits = 4),
      "), so gamma-hat is not the maximum-likelihood estimate and its ",
      "standard error and Wald interval do not describe one. As gamma grows ",
      "this phase tends to a corner law with a finite likelihood, and the ",
      "likelihood rises in that direction here, so the supremum may lie at ",
      "gamma = Inf. The fitted hazard away from tau is barely affected."
    )
  )
}

#' The corner-supremum check over every FIXGE2 g3 phase of a fit
#' @inheritParams .hzr_g3_corner_record
#' @param converged Did the optimizer report convergence?
#' @return A list of boundary records, possibly empty.
#' @keywords internal
#' @noRd
.hzr_g3_corner_supremum <- function(theta, value, converged, objective_fn,
                                    phases, covariate_counts, x_list, time,
                                    status, time_lower, weights) {
  # Right-censored and exact rows only: the corner search is written for
  # them. Anything else is not examined, which costs power, never a false
  # record, since only the certificate can fire.
  if (!isTRUE(converged) || !is.finite(value) ||
        !all(status %in% c(0, 1))) {
    return(list())
  }
  out <- list()
  for (k in seq_along(phases)) {
    ph <- phases[[k]]
    if (!identical(ph$type, "g3") || !identical(ph$constraint, "eta_gamma")) {
      next
    }
    rec <- tryCatch(
      .hzr_g3_corner_record(theta, value, objective_fn, k, phases,
                            covariate_counts, x_list, time, status,
                            time_lower, weights),
      error = function(e) NULL
    )
    if (!is.null(rec)) out[[length(out) + 1L]] <- rec
  }
  out
}
