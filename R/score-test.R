# score-test.R -- Score (Q) statistic for stepwise entry candidates.
#
# Replaces the per-candidate model refit in the forward step. At the current
# model's MLE the reduced-model score is zero, so adding a candidate with its
# coefficient pinned at 0 leaves only the candidate's own score component:
#
#   U_beta = dlogL/dbeta at (theta_hat, beta = 0)
#   V_beta = I_bb - I_bt %*% solve(I_tt) %*% I_tb
#   Q      = U_beta^2 / V_beta      ~ chi^2(1)
#
# I_tt is the information of the CURRENT model and does not depend on the
# candidate, so it is inverted once per step (.hzr_score_nuisance()) and reused
# for every candidate. That reuse is what removes the optimizer from the loop.
#
# SAS's Q is deliberately approximate -- its listing states "variances are
# approximate because shaping parameter covariances are ignored" -- so the
# nuisance partition EXCLUDES shape parameters. Matching that is required for
# parity; the efficient score would not match. See
# inst/dev/SCORE-CRITERION-DESIGN.md.
#
# The multiphase path is gated against SAS's own Q values. No SAS reference
# exists for the single-distribution families -- every proc hazard source in
# this repo is multiphase -- so those are held to a numeric oracle instead:
# numDeriv for the score, and agreement with the refit path's Wald chi-square
# where theory requires it.

#' Phases as the fitted model actually used them
#'
#' `.hzr_optim_multiphase()` stores the validated phase list on the fit; fall
#' back to the spec when it is absent.
#'
#' @noRd
.hzr_score_phases <- function(current) {
  if (!is.null(current$fit$phases)) current$fit$phases else current$spec$phases
}

#' Positions of the non-shape (mu / covariate) parameters in theta
#'
#' SAS ignores shaping-parameter covariances during selection, so only these
#' positions enter the nuisance adjustment.
#'
#' The partition is derived from the phase layout, never from parameter names.
#' Each phase's block is `[log_mu, shapes..., betas...]` (see
#' `.hzr_unpack_phase_theta()`), so the shape slots are known by position:
#' `.hzr_phase_n_shape()` of them, immediately after `log_mu`. A name-based
#' rule cannot do this -- a covariate called `m`, `nu`, `gamma`, `alpha` or
#' `eta` produces a theta name like `constant.m` that is indistinguishable
#' from a shape by name alone (a `constant` phase has no shapes at all), and
#' dropping it silently under-adjusts `V_beta` for every other candidate.
#'
#' @noRd
.hzr_score_free_idx <- function(current) {
  if (current$spec$dist != "multiphase") {
    return(.hzr_score_single_free_idx(current))
  }
  # Multiphase: keep each phase's log_mu and its covariate betas; drop exactly
  # the phase's shape slots.
  phases <- .hzr_score_phases(current)
  counts <- current$fit$covariate_counts
  if (is.null(counts)) {
    # A fitted multiphase object always carries counts. They set every phase's
    # start position via .hzr_log_mu_positions(), so substituting zeros would
    # not degrade gracefully -- it would silently return the wrong index set
    # and a wrong V_beta for every candidate in the step.
    stop(
      ".hzr_score_free_idx(): the fitted model has no `covariate_counts`. ",
      "They are required to locate each phase's parameters in theta.",
      call. = FALSE
    )
  }
  starts <- .hzr_log_mu_positions(phases, counts)
  idx <- integer(0)
  for (nm in names(phases)) {
    start <- starts[[nm]]
    n_shape <- .hzr_phase_n_shape(phases[[nm]])
    n_cov <- as.integer(counts[[nm]])
    idx <- c(idx, start)
    if (n_cov > 0L) {
      beta_start <- start + 1L + n_shape
      idx <- c(idx, seq.int(beta_start, beta_start + n_cov - 1L))
    }
  }
  idx
}

#' Size of a single distribution's leading baseline block in theta
#'
#' `.hzr_shape_parameter_count()` returns the size of the WHOLE leading block --
#' intercept and shapes together -- so within it the intercept is always slot 1
#' and the genuine shape slots are `2:n_base`. The documented layouts are
#' `[log_lambda, betas...]` (exponential, so no shape at all), `[mu, nu,
#' betas...]` (weibull), `[log_alpha, log_beta, betas...]` (log-logistic) and
#' `[mu, log_sigma, betas...]` (log-normal).
#'
#' @noRd
.hzr_score_n_base <- function(current) {
  n_base <- .hzr_shape_parameter_count(
    current$spec$dist, control = current$spec$control
  )
  if (!is.finite(n_base) || n_base < 1L) {
    stop(
      ".hzr_score_free_idx(): no known theta layout for dist ",
      sQuote(current$spec$dist), "; refusing to guess which slots are shapes.",
      call. = FALSE
    )
  }
  as.integer(n_base)
}

#' Non-shape positions in a single distribution's theta
#'
#' Keeps the intercept and every covariate beta; drops ONLY the shape slots.
#' The intercept is a nuisance parameter like any other and must stay in the
#' block -- dropping it under-adjusts `V_beta`, which makes `Q` too small and
#' silently keeps real candidates out of the model. Exponential is the clearest
#' case: it has no shape parameter, so nothing is dropped.
#'
#' @noRd
.hzr_score_single_free_idx <- function(current) {
  theta <- current$fit$theta
  n_base <- .hzr_score_n_base(current)
  p_cov <- if (is.null(current$data$x)) 0L else ncol(current$data$x)
  # An empty theta would make a positional rule produce NEGATIVE subscripts,
  # which R happily accepts as a drop-these-rows instruction and turns into a
  # plausible-looking wrong matrix. Check the layout instead of trusting it.
  if (length(theta) != n_base + p_cov) {
    stop(
      ".hzr_score_free_idx(): theta has ", length(theta), " element(s) but the ",
      current$spec$dist, " layout needs ", n_base + p_cov, " (", n_base,
      " baseline + ", p_cov, " covariate).",
      call. = FALSE
    )
  }
  idx <- 1L
  if (p_cov > 0L) idx <- c(idx, seq.int(n_base + 1L, n_base + p_cov))
  idx
}

#' Log-likelihood / gradient entry points for a single distribution
#'
#' SIGN CONVENTION -- both return the POSITIVE log-likelihood scale, despite
#' what some of their roxygen blocks say. `.hzr_optim_generic()` is what negates
#' them for minimisation, and the analytic `hessian_fn` hook it takes is on the
#' negated (objective) scale. So the observed information is a Hessian of
#' `-logl_fn`, not of `logl_fn`. Getting this backwards yields a negative
#' `v_beta`, which the collinearity floor then turns into a silent `NA` for
#' every candidate.
#'
#' @noRd
.hzr_score_logl_fn <- function(dist) {
  switch(
    dist,
    exponential = .hzr_logl_exponential,
    weibull     = .hzr_logl_weibull,
    loglogistic = .hzr_logl_loglogistic,
    lognormal   = .hzr_logl_lognormal,
    stop(".hzr_score_logl_fn(): unsupported dist ", sQuote(dist), ".",
         call. = FALSE)
  )
}

#' @noRd
.hzr_score_gradient_fn <- function(dist) {
  switch(
    dist,
    exponential = .hzr_gradient_exponential,
    weibull     = .hzr_gradient_weibull,
    loglogistic = .hzr_gradient_loglogistic,
    lognormal   = .hzr_gradient_lognormal,
    stop(".hzr_score_gradient_fn(): unsupported dist ", sQuote(dist), ".",
         call. = FALSE)
  )
}

#' Negative log-likelihood of a single-distribution model at `theta`
#'
#' `x` is the design matrix to evaluate against -- the fit's own for the current
#' model, or the expanded one when a candidate is pinned at zero.
#'
#' @noRd
.hzr_score_single_nll <- function(current, x, theta) {
  d <- current$data
  fn <- .hzr_score_logl_fn(current$spec$dist)
  -fn(theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = x, weights = d$weights)
}

#' Numeric observed information for a single distribution
#'
#' Weibull's analytic Hessian is on an internal reparameterisation, not the
#' `(mu, nu, beta)` scale theta is stored on, and there is no analytic Hessian
#' on the expanded design for the others either. A numeric Hessian of the
#' negative log-likelihood is the one form that is uniform across all four
#' families -- and it is the same oracle the analytic Hessians are themselves
#' tested against.
#'
#' @noRd
.hzr_score_single_hessian <- function(current, x, theta) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    # Returning NULL here would NA every candidate and make stepwise report
    # nothing significant -- a plausible-looking wrong answer.
    stop(
      "The score criterion needs the 'numDeriv' package for dist ",
      sQuote(current$spec$dist), ": its observed information is computed ",
      "numerically. Install 'numDeriv', or select with the Wald criterion.",
      call. = FALSE
    )
  }
  h <- tryCatch(
    numDeriv::hessian(
      function(par) .hzr_score_single_nll(current, x, par), theta
    ),
    error = function(e) NULL
  )
  if (is.null(h) || !is.matrix(h) || nrow(h) != length(theta) ||
        !all(is.finite(h))) {
    return(NULL)
  }
  h
}

#' Numeric observed information for a multiphase fit
#'
#' `.hzr_hessian_multiphase()` declines by design when any row is left- or
#' interval-censored (`status` in `{-1, 2}`): the analytic second derivative
#' is not defined for those contributions. Before this fallback existed the
#' `NULL` propagated to `nuisance$ok = FALSE` and every candidate scored `NA`,
#' so the screen stopped having tested nothing -- and said so in the language
#' of a degenerate column, which is a different fault entirely.
#'
#' This costs a numeric Hessian per candidate, which is the per-candidate work
#' the score criterion exists to avoid. It fires only where the analytic form
#' is unavailable, and slower is the correct trade against selecting nothing.
#'
#' @noRd
.hzr_score_multiphase_hessian <- function(current, theta, phases,
                                           covariate_counts, x_list) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    # Returning NULL here would NA every candidate and make stepwise report
    # nothing significant -- a plausible-looking wrong answer.
    stop(
      "The score criterion needs the 'numDeriv' package for a multiphase fit ",
      "with left- or interval-censored rows: the analytic observed ",
      "information is not defined there, so it must be computed numerically. ",
      "Install 'numDeriv', or select with the Wald criterion.",
      call. = FALSE
    )
  }
  d <- current$data
  nll <- function(par) {
    -.hzr_logl_multiphase(
      par, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights, phases = phases,
      covariate_counts = covariate_counts, x_list = x_list
    )
  }
  h <- tryCatch(numDeriv::hessian(nll, theta), error = function(e) NULL)
  if (is.null(h) || !is.matrix(h) || nrow(h) != length(theta) ||
        !all(is.finite(h))) {
    return(NULL)
  }
  # Mirror .hzr_hessian_multiphase(), which names from `theta`, so the two are
  # interchangeable to a caller. `theta` is unnamed on an expanded model, and
  # this then yields NULL dimnames there -- as the analytic path also does.
  dimnames(h) <- list(names(theta), names(theta))
  h
}

#' Observed information (Hessian of the negative log-likelihood)
#'
#' @noRd
.hzr_score_information <- function(current, theta) {
  d <- current$data
  if (current$spec$dist != "multiphase") {
    return(.hzr_score_single_hessian(current, d$x, theta))
  }
  phases <- .hzr_score_phases(current)
  h <- tryCatch(
    .hzr_hessian_multiphase(
      theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = phases,
      covariate_counts = current$fit$covariate_counts,
      x_list = current$fit$x_list
    ),
    error = function(e) NULL
  )
  if (!is.null(h)) {
    return(h)
  }
  .hzr_score_multiphase_hessian(
    current, theta, phases, current$fit$covariate_counts, current$fit$x_list
  )
}

#' Per-step reusable nuisance block
#'
#' @param current Fitted `hazard` object at the step's current model.
#' @return `list(inv, idx, ok)`. `ok = FALSE` means a nuisance block exists but
#'   could not be formed or inverted, and every candidate must return `NA`
#'   rather than a wrongly unadjusted Q. `ok = TRUE` with `inv = NULL` is the
#'   only legitimately unadjusted case: there are no nuisance parameters at all.
#' @noRd
.hzr_score_nuisance <- function(current) {
  idx <- .hzr_score_free_idx(current)
  if (length(idx) == 0L) {
    # No nuisance parameters exist, so there is nothing to adjust for.
    return(list(inv = NULL, idx = idx, ok = TRUE))
  }
  info <- .hzr_score_information(current, theta = current$fit$theta)
  if (is.null(info)) {
    return(list(inv = NULL, idx = idx, ok = FALSE))
  }
  blk <- info[idx, idx, drop = FALSE]
  inv <- tryCatch(solve(blk), error = function(e) NULL)
  list(inv = inv, idx = idx, ok = !is.null(inv))
}

#' Expand the current model's design and theta with one pinned candidate
#'
#' Mirrors `.hzr_optim_multiphase()`'s construction of `x_list` /
#' `covariate_counts` so the expanded model is what a refit via
#' `.hzr_refit_with_scope(action = "add")` would build, with the new
#' coefficient pinned at zero.
#'
#' @return For multiphase, `list(theta, beta_idx, theta_idx, phases, x_list,
#'   covariate_counts)`; for a single distribution, `list(theta, beta_idx,
#'   theta_idx, x)`. `NULL` when the candidate cannot be expanded (already in
#'   scope, unknown phase, row misalignment).
#' @noRd
.hzr_score_expand <- function(current, var, phase, data) {
  if (current$spec$dist != "multiphase") {
    return(.hzr_score_expand_single(current, var, phase, data))
  }
  phases <- .hzr_score_phases(current)
  if (is.null(phase) || !is.character(phase) || length(phase) != 1L ||
        !phase %in% names(phases)) {
    return(NULL)
  }
  if (var %in% .hzr_scope_current_vars(current, phase)) {
    return(NULL)
  }

  new_phases <- phases
  new_phases[[phase]] <- .hzr_phase_update_formula(
    new_phases[[phase]], action = "add", var = var
  )

  d <- current$data
  n_time <- length(d$time)
  nms <- names(new_phases)

  x_list <- vector("list", length(new_phases))
  names(x_list) <- nms
  cov_counts <- stats::setNames(integer(length(new_phases)), nms)

  for (nm in nms) {
    ph <- new_phases[[nm]]
    if (!is.null(ph$formula) && !is.null(data)) {
      mf_j <- tryCatch(
        stats::model.frame(ph$formula, data = data,
                           na.action = stats::na.pass),
        error = function(e) NULL
      )
      if (is.null(mf_j)) return(NULL)
      x_j <- stats::model.matrix(ph$formula, data = mf_j)[, -1L, drop = FALSE]
      x_list[[nm]] <- x_j
      cov_counts[[nm]] <- ncol(x_j)
    } else if (!is.null(d$x)) {
      x_list[[nm]] <- d$x
      cov_counts[[nm]] <- ncol(d$x)
    } else {
      x_list[[nm]] <- NULL
      cov_counts[[nm]] <- 0L
    }
  }

  old_counts <- current$fit$covariate_counts
  if (is.null(old_counts) ||
        cov_counts[[phase]] != old_counts[[phase]] + 1L) {
    # A term that expands to more than one column (factor, spline) breaks the
    # one-zero-per-add layout; the refit path guards this upstream too.
    return(NULL)
  }
  for (nm in nms) {
    xm <- x_list[[nm]]
    if (is.null(xm) || ncol(xm) == 0L) next
    if (nrow(xm) != n_time || anyNA(xm)) return(NULL)
  }

  parts <- .hzr_split_theta(current$fit$theta, phases, old_counts)

  theta_new <- numeric(0)
  theta_idx <- integer(0)
  beta_idx <- NA_integer_
  pos <- 0L
  for (nm in nms) {
    part <- parts[[nm]]
    theta_new <- c(theta_new, part)
    theta_idx <- c(theta_idx, pos + seq_along(part))
    pos <- pos + length(part)
    if (identical(nm, phase)) {
      # The new covariate is the phase's last term, so its coefficient is the
      # last slot of that phase's sub-vector -- matching model.matrix ordering.
      theta_new <- c(theta_new, 0)
      pos <- pos + 1L
      beta_idx <- pos
    }
  }

  names(theta_new) <- unlist(lapply(nms, function(nm) {
    cov_names <- if (cov_counts[[nm]] > 0L && !is.null(x_list[[nm]])) {
      colnames(x_list[[nm]])
    } else {
      character(0)
    }
    .hzr_phase_theta_names(new_phases[[nm]], nm, cov_names)
  }))

  list(
    theta = theta_new,
    beta_idx = beta_idx,
    theta_idx = theta_idx,
    phases = new_phases,
    x_list = x_list,
    covariate_counts = cov_counts
  )
}

#' Expand a single distribution's design and theta with one pinned candidate
#'
#' Mirrors `.hzr_refit_with_scope()`'s single-distribution path: the candidate
#' becomes the formula's last term, so `model.matrix()` puts its column last and
#' its coefficient takes the last slot of theta -- exactly where the refit's
#' `c(theta_old, 0)` warm start puts it. Here it stays pinned at zero.
#'
#' @noRd
.hzr_score_expand_single <- function(current, var, phase, data) {
  if (!is.null(phase)) {
    return(NULL)
  }
  if (var %in% .hzr_scope_current_vars(current)) {
    return(NULL)
  }

  d <- current$data
  xcand <- .hzr_candidate_numeric(data[[var]])
  if (is.null(xcand) || length(xcand) != length(d$time) || anyNA(xcand)) {
    return(NULL)
  }

  new_col <- matrix(as.numeric(xcand), ncol = 1L,
                    dimnames = list(NULL, var))
  x_new <- if (is.null(d$x)) new_col else cbind(d$x, new_col)

  # theta is UNNAMED on these fits (see R/wald.R); index positionally and do
  # not attach names that the rest of the package would not have produced.
  theta_old <- unname(current$fit$theta)
  theta_new <- c(theta_old, 0)

  list(
    theta = theta_new,
    beta_idx = length(theta_new),
    theta_idx = seq_along(theta_old),
    x = x_new
  )
}

#' Score vector of the expanded model at (theta_hat, beta = 0)
#'
#' @noRd
.hzr_score_gradient <- function(current, exp_) {
  d <- current$data
  if (current$spec$dist != "multiphase") {
    fn <- .hzr_score_gradient_fn(current$spec$dist)
    g <- tryCatch(
      fn(exp_$theta, time = d$time, status = d$status,
         time_lower = d$time_lower, time_upper = d$time_upper,
         x = exp_$x, weights = d$weights),
      error = function(e) NULL
    )
    if (is.null(g) || length(g) != length(exp_$theta) || !all(is.finite(g))) {
      return(NULL)
    }
    return(as.numeric(g))
  }
  g <- tryCatch(
    .hzr_gradient_multiphase(
      exp_$theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = exp_$phases,
      covariate_counts = exp_$covariate_counts,
      x_list = exp_$x_list
    ),
    error = function(e) NULL
  )
  if (is.null(g) || length(g) != length(exp_$theta)) return(NULL)
  g
}

#' Observed information of the expanded model at (theta_hat, beta = 0)
#'
#' @noRd
.hzr_score_information_expanded <- function(current, exp_) {
  d <- current$data
  if (current$spec$dist != "multiphase") {
    return(.hzr_score_single_hessian(current, exp_$x, exp_$theta))
  }
  h <- tryCatch(
    .hzr_hessian_multiphase(
      exp_$theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = exp_$phases,
      covariate_counts = exp_$covariate_counts,
      x_list = exp_$x_list
    ),
    error = function(e) NULL
  )
  if (!is.null(h) && is.matrix(h) && nrow(h) == length(exp_$theta)) {
    return(h)
  }
  .hzr_score_multiphase_hessian(
    current, exp_$theta, exp_$phases, exp_$covariate_counts, exp_$x_list
  )
}

#' Score statistic for one entry candidate
#'
#' @param current Fitted `hazard` object (the step's current model).
#' @param var Character scalar; candidate column name in `data`.
#' @param phase Character scalar naming the phase, or `NULL` for
#'   single-distribution fits.
#' @param data Data frame the model was fitted on.
#' @param nuisance Optional result of `.hzr_score_nuisance(current)`; recomputed
#'   when `NULL`. Pass it to reuse across candidates within a step.
#' @return `list(stat, df, p_value)`. `stat`/`p_value` are `NA_real_` for a
#'   degenerate candidate, a collinear candidate, or an unusable nuisance block.
#' @noRd
.hzr_score_q <- function(current, var, phase = NULL, data,
                         nuisance = NULL) {
  # Every NA return carries WHY. The reasons are not interchangeable: a
  # collinear column should be dropped, while an indefinite information matrix
  # usually means the candidate is among the strongest on offer. Reporting the
  # second as the first tells a user to discard their best variable.
  na_result <- function(reason) {
    list(stat = NA_real_, df = 1L, p_value = NA_real_, reason = reason)
  }

  # `data` must be row-aligned with the fit: the candidate column is read from
  # `data` while the score is evaluated on the fit's own stored rows. A caller
  # passing pre-`na.omit()` data would fail the row check inside
  # .hzr_score_expand() for EVERY candidate, and stepwise would report nothing
  # significant -- a plausible-looking wrong answer. Fail loudly instead.
  n_obs <- length(current$data$time)
  if (nrow(data) != n_obs) {
    stop(
      "`data` has ", nrow(data), " rows but the fitted model used ", n_obs,
      ". The score test needs `data` row-aligned with the fit; pass the same ",
      "data frame the model was fitted on (after any NA removal).",
      call. = FALSE
    )
  }

  xcand <- .hzr_candidate_numeric(data[[var]])
  if (is.null(xcand) || anyNA(xcand)) {
    return(na_result("non_numeric"))
  }
  s <- stats::sd(xcand)
  if (!is.finite(s) || s == 0) {
    return(na_result("constant"))
  }
  if (is.null(nuisance)) nuisance <- .hzr_score_nuisance(current)
  # A nuisance block that exists but could not be inverted must not fall
  # through to an unadjusted (too large) v_beta.
  if (!isTRUE(nuisance$ok)) return(na_result("nuisance_singular"))

  exp_ <- .hzr_score_expand(current, var, phase, data)
  if (is.null(exp_)) return(na_result("not_expandable"))

  grad <- .hzr_score_gradient(current, exp_)
  info <- .hzr_score_information_expanded(current, exp_)
  if (is.null(grad) || is.null(info)) return(na_result("no_information"))

  b <- exp_$beta_idx
  u_beta <- grad[b]
  i_bb <- info[b, b]

  if (!is.finite(u_beta) || !is.finite(i_bb)) {
    return(na_result("nonfinite"))
  }

  # The candidate's OWN curvature, before any adjustment for the model. SAS
  # tests this separately and before the tolerance test below -- `q1.c`:
  #
  #     diag = b1[n];
  #     if (diag <= ZERO) { *err = 2; return ZERO; }   /* flag 2 */
  #     ...
  #     if (*qtol <= ZERO) { *err = 3; return ZERO; }  /* flag 3 */
  #
  # and the two are different diagnoses. This one says the candidate's own
  # observed information is not positive; the tolerance test says the
  # candidate is unusable *given what is already in the model*. Folding them
  # together reports a candidate whose own curvature is wrong as though the
  # current model were responsible for it.
  #
  # Reachable, and where this package's other censoring faults live: a
  # multiphase fit with a large share of interval-censored rows drives i_bb
  # slightly negative, because those rows contribute a difference of survival
  # terms rather than a log density.
  if (i_bb <= 0) {
    return(na_result("information_nonpositive"))
  }

  v_beta <- i_bb
  if (!is.null(nuisance$inv) && length(nuisance$idx) > 0L) {
    t_idx <- exp_$theta_idx[nuisance$idx]
    i_bt <- info[b, t_idx, drop = FALSE]
    v_beta <- i_bb - as.numeric(i_bt %*% nuisance$inv %*% t(i_bt))
  }

  if (!is.finite(v_beta)) {
    return(na_result("nonfinite"))
  }

  # v_beta is a Schur complement of the OBSERVED information at beta = 0, and
  # observed information is not positive definite away from a maximum. Two
  # unrelated failures both land on "v_beta is unusable", and they carry
  # OPPOSITE meanings:
  #
  #   |v_beta| <= tol   collinearity. I_bt %*% solve(I_tt) %*% I_tb -> I_bb,
  #                     so v_beta is a difference of two nearly equal numbers
  #                     and its true value is 0. Left unguarded it lands at
  #                     ~1e-14 and Q ~ 1e15 wins the step.
  #   v_beta < -tol     the log-likelihood curves UPWARD in beta at 0, so the
  #                     quadratic approximation the score test rests on does
  #                     not hold there. Reached when the candidate's effect is
  #                     far from zero -- a STRONG candidate -- and also when
  #                     `current` has not truly converged.
  #
  # The magnitude test comes FIRST and the sign test is against -tol, not 0.
  # The `-tol` is redundant given the early return above, and is written out
  # anyway so that the rule holds on its own: reordering these two guards must
  # not be able to reintroduce the defect below. For an exactly collinear candidate the true v_beta is
  # exactly 0, so its computed sign is decided by rounding: an `if (v_beta <
  # 0)` ahead of the floor reported a perfect duplicate as
  # "information_indefinite" -- "this is a strong candidate, keep it" -- on
  # Linux while reporting "collinear" on macOS, from the same code and data.
  # Caught by CI, invisible on one platform. Inside the band the two are not
  # distinguishable, and zero means collinear.
  #
  # i_bb is known positive by the guard above, so tol needs no abs(); an
  # earlier signed floor let a slightly negative v_beta through and returned
  # a negative Q.
  tol <- i_bb * sqrt(.Machine$double.eps)
  if (abs(v_beta) <= tol) {
    return(na_result("collinear"))
  }
  if (v_beta < -tol) {
    return(na_result("information_indefinite"))
  }
  stat <- (u_beta^2) / v_beta

  # The coefficient this Q implies: one Newton step from beta = 0, with
  # standard error 1 / sqrt(v_beta). SAS forms the same quantity and rejects
  # the candidate when it is absurd -- `dqstat.c`:
  #
  #     *qz    = sqrt(q);
  #     *qse   = -(*qz)/d1llad;
  #     *qbeta = (*qz)*(*qse);          /* = u_beta / v_beta */
  #     if (fabs(*qbeta) > 50.0e0) { *qflag = 4; return; }
  #
  # calling it "the model is going to infinity ... legitimate for a variable
  # that is either positive or negative with respect to all the remaining
  # events in its phase. However the model cannot manage this."
  #
  # This is the guard against a Q that is finite, enormous and meaningless.
  # The other two guards do not reach it: v_beta stays positive and well
  # above the collinearity floor. Measured on `avc` with a near-collinear
  # candidate, scored against a model moved 0.25 off its optimum -- which is
  # what a failed refit leaves behind -- Q was 6.5e7 with v_beta = 0.0076 and
  # no reason reported at all. A legitimate candidate at the optimum reached
  # |qbeta| of about 14 in the same setting, so 50 leaves real headroom.
  q_beta <- u_beta / v_beta
  if (!is.finite(q_beta) || abs(q_beta) > 50) {
    return(na_result("coefficient_diverging"))
  }

  list(stat = stat, df = 1L,
       p_value = stats::pchisq(stat, df = 1L, lower.tail = FALSE),
       reason = NA_character_)
}


#' One-line explanation of an unscorable candidate
#'
#' Maps `.hzr_score_q()`'s reason codes to prose for a warning. An unknown code
#' passes through as itself, so a code added later still reports rather than
#' silently becoming an empty string.
#'
#' @noRd
.hzr_score_reason_text <- function(reason) {
  txt <- c(
    information_indefinite = paste(
      "the observed information at beta = 0 was indefinite. That usually means",
      "the candidate's effect is too large for the score test's approximation",
      "at zero, so these are typically STRONG candidates rather than",
      "degenerate ones; `criterion = \"wald\"` tests them. It is also reached",
      "when the current fit has not truly converged"
    ),
    information_nonpositive = paste(
      "the candidate's own observed information was not positive, before any",
      "adjustment for the current model. This is a different fault from",
      "collinearity: the candidate is a poor one in itself, or the fit it",
      "would be added to is not at a maximum"
    ),
    coefficient_diverging = paste(
      "the coefficient the score implies exceeds +/-50, so the fit for this",
      "candidate is running to infinity. That is legitimate for a variable",
      "that separates the remaining events in its phase, and the model cannot",
      "represent it; it is also what a candidate scored against a model that",
      "is NOT at its optimum looks like, so check whether an earlier step's",
      "refit failed"
    ),
    collinear = "the candidate was collinear with the current model",
    constant  = "the candidate column was constant",
    non_numeric = "the candidate column was not numeric, or held NA",
    nuisance_singular = paste(
      "the current model's information matrix could not be inverted, so no",
      "candidate could be scored at that step"
    ),
    no_information = "no observed information was available for the candidate",
    not_expandable = "the candidate could not be added to the model",
    nonfinite = "the score or its variance was not finite"
  )
  out <- unname(txt[reason])
  out[is.na(out)] <- reason[is.na(out)]
  out
}

#' Tally reason codes into a named integer vector, commonest first
#'
#' @noRd
.hzr_tally_reasons <- function(reasons) {
  reasons <- reasons[!is.na(reasons)]
  if (length(reasons) == 0L) {
    return(stats::setNames(integer(0), character(0)))
  }
  tab <- sort(table(reasons), decreasing = TRUE)
  stats::setNames(as.integer(tab), names(tab))
}

#' Add two reason tallies together
#'
#' @noRd
.hzr_merge_reasons <- function(a, b) {
  if (length(b) == 0L) return(a)
  if (length(a) == 0L) return(b)
  keys <- union(names(a), names(b))
  out <- stats::setNames(integer(length(keys)), keys)
  out[names(a)] <- out[names(a)] + as.integer(a)
  out[names(b)] <- out[names(b)] + as.integer(b)
  sort(out, decreasing = TRUE)
}

#' Render a reason tally as "n x <prose>" lines
#'
#' @noRd
.hzr_format_reasons <- function(tally) {
  if (length(tally) == 0L) return("")
  paste0(
    " Causes: ",
    paste0(as.integer(tally), " x ", .hzr_score_reason_text(names(tally)),
           collapse = "; "),
    "."
  )
}
