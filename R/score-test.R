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
# Scope: multiphase only. The single-distribution families land in a follow-up
# task; their branches below return NULL rather than a wrong number.

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
  theta <- current$fit$theta
  if (current$spec$dist != "multiphase") {
    # Single-distribution: theta is [shape..., beta...]; covariate coefficients
    # occupy the trailing ncol(x) slots. Shapes are excluded per SAS.
    p_cov <- if (is.null(current$data$x)) 0L else ncol(current$data$x)
    if (p_cov == 0L) return(integer(0))
    return(seq.int(length(theta) - p_cov + 1L, length(theta)))
  }
  # Multiphase: keep each phase's log_mu and its covariate betas; drop exactly
  # the phase's shape slots.
  phases <- .hzr_score_phases(current)
  counts <- current$fit$covariate_counts
  if (is.null(counts)) {
    counts <- stats::setNames(integer(length(phases)), names(phases))
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

#' Observed information (Hessian of the negative log-likelihood)
#'
#' @noRd
.hzr_score_information <- function(current, theta) {
  d <- current$data
  if (current$spec$dist != "multiphase") {
    # Single-distribution families are out of scope for this task; returning
    # NULL leaves the candidate variance unadjusted rather than wrong.
    return(NULL)
  }
  tryCatch(
    .hzr_hessian_multiphase(
      theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = .hzr_score_phases(current),
      covariate_counts = current$fit$covariate_counts,
      x_list = current$fit$x_list
    ),
    error = function(e) NULL
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
#' @return `list(theta, beta_idx, theta_idx, phases, x_list,
#'   covariate_counts)`, or `NULL` when the candidate cannot be expanded
#'   (already in scope, unknown phase, row misalignment, non-multiphase).
#' @noRd
.hzr_score_expand <- function(current, var, phase, data) {
  if (current$spec$dist != "multiphase") {
    return(NULL)
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

#' Score vector of the expanded model at (theta_hat, beta = 0)
#'
#' @noRd
.hzr_score_gradient <- function(current, exp_) {
  d <- current$data
  if (current$spec$dist != "multiphase") {
    return(NULL)
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
    return(NULL)
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
  if (is.null(h) || !is.matrix(h) || nrow(h) != length(exp_$theta)) {
    return(NULL)
  }
  h
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
  na_result <- list(stat = NA_real_, df = 1L, p_value = NA_real_)

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

  xcand <- data[[var]]
  if (is.null(xcand) || !is.numeric(xcand) ||
        anyNA(xcand) || stats::sd(xcand) == 0) {
    return(na_result)
  }
  if (is.null(nuisance)) nuisance <- .hzr_score_nuisance(current)
  # A nuisance block that exists but could not be inverted must not fall
  # through to an unadjusted (too large) v_beta.
  if (!isTRUE(nuisance$ok)) return(na_result)

  exp_ <- .hzr_score_expand(current, var, phase, data)
  if (is.null(exp_)) return(na_result)

  grad <- .hzr_score_gradient(current, exp_)
  info <- .hzr_score_information_expanded(current, exp_)
  if (is.null(grad) || is.null(info)) return(na_result)

  b <- exp_$beta_idx
  u_beta <- grad[b]
  i_bb <- info[b, b]

  v_beta <- i_bb
  if (!is.null(nuisance$inv) && length(nuisance$idx) > 0L) {
    t_idx <- exp_$theta_idx[nuisance$idx]
    i_bt <- info[b, t_idx, drop = FALSE]
    v_beta <- i_bb - as.numeric(i_bt %*% nuisance$inv %*% t(i_bt))
  }

  # A collinear candidate drives I_bt %*% solve(I_tt) %*% I_tb -> I_bb, so
  # v_beta collapses towards zero from ABOVE: an absolute `v_beta <= 0` test
  # never fires, leaves v_beta ~ 1e-14, and Q ~ 1e15 wins the step -- exactly
  # the failure the NA contract exists to prevent. Floor v_beta relative to the
  # unadjusted i_bb it is a difference of, at the scale where cancellation has
  # eaten the significant digits.
  if (!is.finite(u_beta) || !is.finite(v_beta) || !is.finite(i_bb) ||
        v_beta <= i_bb * sqrt(.Machine$double.eps)) {
    return(na_result)
  }
  stat <- (u_beta^2) / v_beta
  list(stat = stat, df = 1L,
       p_value = stats::pchisq(stat, df = 1L, lower.tail = FALSE))
}
