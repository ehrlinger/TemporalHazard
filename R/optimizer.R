#' @importFrom stats optim
#' @keywords internal
NULL

# optimizer.R -- Generic parametric hazard optimizer
#
# Consolidates the common optimization pattern shared across all distributions.
# Each distribution provides its own log-likelihood and gradient functions;
# this module handles the objective wrapping, optim() call, and post-fit
# Hessian/vcov computation.

# numDeriv is a Suggests. Wrapped so tests can simulate an install without it:
# requireNamespace() is base R and cannot be mocked directly.
.hzr_numderiv_available <- function() {
  requireNamespace("numDeriv", quietly = TRUE)
}

# Central differences of an objective, for when the analytic gradient is not
# the objective's own (Conservation of Events). A point on the 1e10 clamp makes
# that component NA rather than a number that looks measured. Components in
# `sign_bounded` -- the multiphase shape m -- never straddle 0, where the cdf
# and hazard families meet in a cusp (#251), and follow the rule of
# .hzr_phase_derivatives(): one-sided second order forward from m >= 0; from
# m < 0 a step capped at 1% of |m| (floored at 1e-10), one-sided backward if
# it would still reach 0.
.hzr_fd_gradient <- function(objective, theta, sign_bounded = integer(0)) {
  f0 <- NULL
  vapply(seq_along(theta), function(i) {
    x <- theta[i]
    h <- .Machine$double.eps^(1 / 3) * max(abs(x), 1)
    side <- 0
    if (i %in% sign_bounded) {
      if (x < 0) h <- min(h, max(0.01 * abs(x), 1e-10))
      side <- if (x >= 0 && x - h < 0) {
        1
      } else if (x < 0 && x + h >= 0) {
        -1
      } else {
        0
      }
    }
    at <- function(step) {
      z <- theta
      z[i] <- x + step
      objective(z)
    }
    if (side == 0) {
      up <- at(h)
      down <- at(-h)
      if (up >= 1e10 || down >= 1e10) return(NA_real_)
      return((up - down) / (2 * h))
    }
    if (is.null(f0)) f0 <<- objective(theta)
    near <- at(side * h)
    far <- at(2 * side * h)
    if (f0 >= 1e10 || near >= 1e10 || far >= 1e10) return(NA_real_)
    side * (-3 * f0 + 4 * near - far) / (2 * h)
  }, numeric(1))
}

#' Generic optimizer for parametric hazard likelihoods
#'
#' Maximises a log-likelihood by minimising its negation via \code{stats::optim}.
#' Handles constrained (L-BFGS-B) and unconstrained (BFGS) optimization,
#' post-fit Hessian computation, and vcov extraction.
#'
#' @param logl_fn Function(theta, time, status, time_lower, time_upper, x, ...)
#'   returning the scalar log-likelihood.  Must accept \code{return_gradient}
#'   argument (not used here directly, but the function signature must allow it).
#' @param gradient_fn Function(theta, time, status, time_lower, time_upper, x, ...)
#'   returning the score vector (gradient of the log-likelihood, NOT negated).
#' @param time Numeric vector of follow-up times.
#' @param status Numeric vector of event/censoring indicators.
#' @param time_lower Optional lower bounds for interval censoring.
#' @param time_upper Optional upper bounds for left/interval censoring.
#' @param x Optional design matrix.
#' @param theta_start Starting parameter vector.
#' @param control List of control options (maxit, reltol, abstol).
#' @param use_bounds Logical; if TRUE, use L-BFGS-B with lower bounds on
#'   constrained parameters.  Used for Weibull where mu/nu are on the natural
#'   scale.
#' @param lower_bounds Numeric vector of lower bounds for L-BFGS-B, the same
#'   length as \code{theta_start}.  Only used when \code{use_bounds = TRUE}.
#'   Defaults to \code{rep(1e-6, length(theta_start))} when NULL, but callers
#'   should supply distribution-specific bounds (e.g. \code{c(1e-6, 1e-6,
#'   rep(-Inf, p_cov))} for Weibull so that covariate betas are unconstrained).
#' @param hessian_fn Optional function(theta) returning the Hessian of the
#'   negative log-likelihood at \code{theta} (same scale as
#'   \code{numDeriv::hessian(objective)}).  When it returns a conformant square
#'   matrix (dimension equal to the number of parameters), that matrix is used
#'   for standard errors.  The numerical Hessian is used instead when
#'   \code{hessian_fn} is \code{NULL} (the default), when the function returns
#'   \code{NULL} (e.g. a censoring branch it does not cover analytically), or
#'   when it errors.  A non-NULL, non-conformant return raises a warning and
#'   also falls back to the numerical Hessian.
#' @param gradient_exact Logical; `TRUE` (the default) when `gradient_fn` is
#'   the gradient of the objective being maximised. Conservation of Events
#'   passes `FALSE`: its `gradient_fn` is the partial score at the conserved
#'   theta, which leaves out how the conserved `log_mu` moves with the free
#'   parameters. SAS/C's acceptance test is then computed from finite
#'   differences of the objective itself, with the scale re-solved at every
#'   step as SAS/C does (`setobj.c`). The `nlm()` continuation keeps the
#'   analytic score: with finite differences it walks onto the CoE solve's
#'   discontinuity where no events are left to conserve.
#' @param sign_bounded Integer positions in `theta_start` whose
#'   finite-difference stencil must not cross 0 (the multiphase shape `m`,
#'   where the phase families meet in a cusp). Used only when
#'   `gradient_exact = FALSE`.
#'
#' @return List with par, value (log-likelihood), convergence, counts, message,
#'   hessian, vcov. Includes \code{se_unavailable_reason}.
#' @noRd
.hzr_optim_generic <- function(
    logl_fn,
    gradient_fn,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    theta_start,
    weights = NULL,
    control = list(),
    use_bounds = FALSE,
    lower_bounds = NULL,
    hessian_fn = NULL,
    gradient_exact = TRUE,
    sign_bounded = integer(0)) {

  control <- utils::modifyList(
    list(maxit = 1000, reltol = 1e-5, abstol = 1e-6),
    control
  )

  # Resolve lower bounds: default to 1e-6 on all params only when the caller
  # has not supplied explicit bounds.  Callers like .hzr_optim_weibull pass
  # c(1e-6, 1e-6, -Inf, ...) so that shape params are constrained but betas
  # are free.
  if (use_bounds && is.null(lower_bounds)) {
    lower_bounds <- rep(1e-6, length(theta_start))
  }

  # Negative log-likelihood for minimization
  objective <- function(theta) {
    if (!all(is.finite(theta))) return(1e10)
    # Only penalise infeasible shape params (those with finite lower bounds > -Inf).
    if (use_bounds && any(theta[is.finite(lower_bounds)] <= lower_bounds[is.finite(lower_bounds)])) return(1e10)

    ll <- -logl_fn(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights, return_gradient = FALSE
    )
    if (!is.finite(ll)) return(1e10)
    ll
  }

  # Negative score vector for minimization
  gradient <- function(theta) {
    if (!all(is.finite(theta))) return(rep(0, length(theta)))
    if (use_bounds && any(
      theta[is.finite(lower_bounds)] <= lower_bounds[is.finite(lower_bounds)]
    )) {
      return(rep(0, length(theta)))
    }

    grad <- tryCatch(
      gradient_fn(
        theta = theta, time = time, status = status,
        time_lower = time_lower, time_upper = time_upper,
        x = x, weights = weights
      ),
      error = function(e) rep(0, length(theta))
    )

    grad[!is.finite(grad)] <- 0
    -grad
  }

  # Dispatch to constrained or unconstrained optimizer
  if (use_bounds) {
    result <- optim(
      par = theta_start, fn = objective, gr = gradient,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = rep(Inf, length(theta_start)),
      control = list(
        maxit = control$maxit,
        factr = 1 / control$reltol,
        pgtol = control$abstol
      ),
      hessian = FALSE
    )
  } else {
    result <- optim(
      par = theta_start, fn = objective, gr = gradient,
      method = "BFGS",
      control = list(maxit = control$maxit, reltol = control$reltol),
      hessian = FALSE
    )
  }

  # SAS/C's acceptance test (src/optim/umstop.c): the optimum is accepted
  # only when the relative gradient max_i |g_i| * max(|x_i|, 1) / max(|f|, 1)
  # is at most gradtl = eps^(1/3).  optim()'s BFGS stops on the relative
  # change in the objective instead, and at the default reltol = 1e-5 that
  # accepts a flat ridge well short of the optimum while reporting
  # convergence 0: on the test suite, 70% of converged BFGS stops failed
  # SAS's test.  When BFGS reports convergence and the test fails, polish
  # with stats::nlm() -- R's Dennis-Schnabel UNCMIN, the algorithm SAS/C's
  # optimizer was ported from -- at SAS's tolerances; its typsize = 1 and
  # fscale = 1 defaults are SAS's typx and typf.  The polished point is kept
  # only when it improves the objective.  L-BFGS-B stops on a projected
  # gradient already, so the bounded path is left alone.
  gradtl <- .Machine$double.eps^(1 / 3)
  # NA, never 0, wherever the gradient cannot be trusted. The wrapped
  # gradient() above returns zeros at a clamped or failing point, and a zero
  # there would read as a pass; so this calls gradient_fn itself (or, when
  # gradient_exact is FALSE, differences the objective) and refuses the 1e10
  # sentinel, a non-finite point, and any non-finite component.
  # The acceptance test needs the unsanitised score. The multiphase gradient
  # zeroes components it cannot evaluate, which the optimizer needs and this
  # test must not see: a zero there reads as a pass over a point that was not
  # scored. Asked only of a gradient_fn that can take the request (it names
  # `sanitize` or `...`); the single-distribution gradients take neither and
  # do not sanitise.
  raw_score_arg <- if (any(c("sanitize", "...") %in% names(formals(gradient_fn)))) {
    list(sanitize = FALSE)
  } else {
    list()
  }
  rel_gradient <- function(theta, value) {
    if (!all(is.finite(theta)) || !is.finite(value) || value >= 1e10) {
      return(NA_real_)
    }
    g <- if (gradient_exact) {
      tryCatch(
        do.call(gradient_fn, c(
          list(theta = theta, time = time, status = status,
               time_lower = time_lower, time_upper = time_upper,
               x = x, weights = weights),
          raw_score_arg
        )),
        error = function(e) NULL
      )
    } else {
      .hzr_fd_gradient(objective, theta, sign_bounded)
    }
    if (is.null(g) || length(g) != length(theta) || !all(is.finite(g))) {
      return(NA_real_)
    }
    max(abs(g) * pmax(abs(theta), 1)) / max(abs(value), 1)
  }
  rel_grad <- NA_real_
  polish_code <- NA_integer_
  if (!use_bounds && result$convergence == 0L) {
    rel_grad <- rel_gradient(result$par, result$value)
    if (is.finite(rel_grad) && rel_grad > gradtl) {
      f_nlm <- function(theta) {
        v <- objective(theta)
        attr(v, "gradient") <- gradient(theta)
        v
      }
      polish <- tryCatch(
        suppressWarnings(stats::nlm(
          f_nlm, result$par, gradtol = gradtl,
          steptol = .Machine$double.eps^(2 / 3), iterlim = control$maxit,
          check.analyticals = FALSE
        )),
        error = function(e) NULL
      )
      if (!is.null(polish) && is.finite(polish$minimum) &&
          polish$minimum < result$value) {
        # nlm() drops names; optim() keeps them, and the single-distribution
        # fits hand par straight back to hazard().
        result$par   <- stats::setNames(polish$estimate, names(result$par))
        result$value <- polish$minimum
        polish_code  <- as.integer(polish$code)
        rel_grad     <- rel_gradient(result$par, result$value)
        # counts still describe the BFGS run alone, so say the fit went on.
        result$message <- paste0(
          if (length(result$message)) paste0(result$message, "; ") else "",
          "continued with nlm() for ", polish$iterations,
          " iterations (code ", polish$code, ")"
        )
      }
    }
  }

  # Post-fit Hessian for standard errors.  Prefer the caller's analytic Hessian
  # (on the objective / negative-log-likelihood scale) when supplied; otherwise,
  # or when it declines by returning NULL, fall back to a numerical Hessian.
  hess_result <- NULL
  if (!is.null(hessian_fn)) {
    # A hook that ERRORS is a different thing from a hook that declines by
    # returning NULL, and swallowing the error made the two indistinguishable
    # downstream. Say which happened.
    hess_result <- tryCatch(
      hessian_fn(result$par),
      error = function(e) {
        warning("hessian_fn() errored, so the analytic Hessian was not used: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    # A non-NULL hook result must be a square matrix matching the parameter
    # dimension; a misbehaving hook should not silently produce a wrong vcov.
    p_dim <- length(result$par)
    if (!is.null(hess_result) &&
        !(is.matrix(hess_result) && all(dim(hess_result) == p_dim))) {
      warning("hessian_fn returned a non-conformant result; using numerical Hessian")
      hess_result <- NULL
    }
  }
  # Why standard errors are unavailable, when they are, for the record that
  # print() and summary() show (#242). Set beside each warning below, so the
  # warning and the record name the same cause.
  se_reason <- NA_character_
  if (is.null(hess_result)) {
    # Reached when no analytic Hessian was supplied, or the hook declined by
    # returning NULL, or it returned something non-conformant, or it errored.
    # Each of those is reported at its own site above; this branch only knows
    # that it has no analytic Hessian, so the message says exactly that rather
    # than asserting a reason it cannot distinguish.
    #
    # numDeriv is the documented fallback, but it is a *Suggests*, so it is
    # legitimately absent on a machine installed without Suggests -- and then
    # there is no third option and no standard errors.
    #
    # That combination used to fail silently: rcond and pd came back NA and
    # vcov() returned a bare logical, with nothing anywhere naming numDeriv.
    # The user-visible symptom was diag(vcov(fit)) complaining about an
    # invalid 'nrow', which is unrecognisable from the cause.
    if (!.hzr_numderiv_available()) {
      warning("No analytic Hessian was obtained for this fit, and the ",
              "'numDeriv' fallback is not installed, so standard errors ",
              "cannot be computed. Install it with ",
              "install.packages(\"numDeriv\"); note numDeriv is a Suggests ",
              "dependency, so install.packages() and install_github() do ",
              "not pull it by default.", call. = FALSE)
      se_reason <- "numDeriv not installed and no analytic Hessian"
    } else {
      hess_result <- tryCatch(
        numDeriv::hessian(objective, result$par),
        error = function(e) {
          warning("numDeriv::hessian() failed, so standard errors are ",
                  "unavailable: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (is.null(hess_result)) se_reason <- "numDeriv::hessian() failed"
    }
  }

  # Hardened inversion + conditioning diagnostics (Layer 1).
  inv <- if (is.matrix(hess_result)) {
    .hzr_safe_solve(hess_result)
  } else {
    # Reached only when no Hessian could be produced at all.
    # .hzr_safe_solve() warns on every degenerate path it handles, but this
    # branch bypasses it entirely, so it needs its own voice or the fit
    # returns NA diagnostics mutely.
    warning("No Hessian could be computed for this fit; standard errors, ",
            "rcond and pd are all NA.", call. = FALSE)
    list(vcov = NA, rcond = NA_real_, pd = NA, reason = se_reason)
  }

  list(
    par = result$par,
    value = -result$value,
    convergence = result$convergence,
    counts = result$counts,
    message = result$message,
    hessian = hess_result,
    vcov = inv$vcov,
    rcond = inv$rcond,
    pd = inv$pd,
    se_unavailable_reason = if (is.matrix(inv$vcov)) NA_character_ else inv$reason,
    # SAS/C's relative gradient at the returned point, after any polish; NA
    # when not evaluated (the bounded path, or BFGS did not converge).
    rel_gradient = rel_grad,
    # stats::nlm()'s termination code when the polish ran and was kept.
    polish_code = polish_code
  )
}
