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
#' rule cannot do this: a covariate called `m`, `nu`, `gamma`, `alpha` or
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
#' `.hzr_shape_parameter_count()` returns the size of the WHOLE leading block
#' (intercept and shapes together), so within it the intercept is always slot 1
#' and the genuine shape slots are `2:n_base`. The documented layouts are
#' `[log_lambda, betas...]` (exponential, so no shape at all), `[mu, nu,
#' betas...]` (weibull), `[log_alpha, log_beta, betas...]` (log-logistic) and
#' `[mu, log_sigma, betas...]` (log-normal).
#'
#' @noRd
.hzr_score_n_base <- function(current) {
  # The likelihood's own count, as hazard() checks theta with: it ignores
  # control$shape_param_count, so a fit made with a misleading one is valid
  # and must not be read with a different layout here (#489).
  n_base <- .hzr_shape_parameter_count(current$spec$dist)
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
#' block; dropping it under-adjusts `V_beta`, which makes `Q` too small and
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
#' SIGN CONVENTION: both return the POSITIVE log-likelihood scale, as their
#' roxygen blocks say. `.hzr_optim_generic()` is what negates
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
#' `x` is the design matrix to evaluate against: the fit's own for the current
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

#' Absorb a NUMERICAL failure, but let a data defect through (#407)
#'
#' The score path deliberately swallows a Hessian it cannot build or invert
#' and reports it as a numerical failure. It used to swallow EVERY error, so
#' a data defect raised inside the likelihood came back as "the information
#' matrix could not be inverted".
#'
#' The re-raise has to happen OUTSIDE the `tryCatch`, not from a handler.
#' Measured: in `tryCatch(expr, hzr_data_error = function(e) stop(e), error
#' = function(e) NULL)` the `stop(e)` is caught by the SIBLING `error`
#' handler of that same call, so the classed error is swallowed anyway and
#' the narrowing is inert while reading as correct.
#'
#' @param expr Expression to evaluate.
#' @return The value, or `NULL` if it failed numerically. An
#'   `hzr_data_error` propagates.
#' @keywords internal
#' @noRd
.hzr_score_try <- function(expr) {
  out <- tryCatch(expr,
                  hzr_data_error = function(e) e,
                  error = function(e) NULL)
  if (inherits(out, "hzr_data_error")) stop(out)
  out
}

#' Numeric observed information for a single distribution
#'
#' Weibull's analytic Hessian is on an internal reparameterisation, not the
#' `(mu, nu, beta)` scale theta is stored on, and there is no analytic Hessian
#' on the expanded design for the others either. A numeric Hessian of the
#' negative log-likelihood is the one form that is uniform across all four
#' families, and it is the same oracle the analytic Hessians are themselves
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
  h <- .hzr_score_try(
    numDeriv::hessian(
      function(par) .hzr_score_single_nll(current, x, par), theta
    )
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
#' so the screen stopped having tested nothing, and said so in the language
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
      covariate_counts = covariate_counts, x_list = x_list,
      objective = .hzr_fit_objective(current)
    )
  }
  h <- .hzr_score_try(numDeriv::hessian(nll, theta))
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
  h <- .hzr_score_try(
    .hzr_hessian_multiphase(
      theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = phases,
      covariate_counts = current$fit$covariate_counts,
      x_list = current$fit$x_list
    )
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

#' Would hazard() refuse this design for a repeated column name?
#'
#' The score path builds each candidate's design itself rather than refitting,
#' so it never meets hazard()'s refusal. Ask the same helper hazard() uses, so
#' the two agree on what counts as a duplicate.
#'
#' @noRd
.hzr_score_duplicate_columns <- function(x) {
  inherits(tryCatch(.hzr_refuse_duplicate_columns(x), error = function(e) e),
           "error")
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
#'   scope, unknown phase, row misalignment), and `list(reason =
#'   "duplicate_column")` when the expanded design repeats a column name,
#'   which hazard() would refuse.
#' @noRd
.hzr_score_expand <- function(current, var, phase, data, term = var) {
  if (current$spec$dist != "multiphase") {
    return(.hzr_score_expand_single(current, var, phase, data, term = term))
  }
  # The phase formula is rebuilt from TEXT, so it is given the candidate's
  # term label, not the column name `var`: pasted, a column `x2 ` read back
  # as `x2`, and the score was computed on the wrong column (#449).
  phases <- .hzr_score_phases(current)
  if (is.null(phase) || !is.character(phase) || length(phase) != 1L ||
        !phase %in% names(phases)) {
    return(NULL)
  }
  if (term %in% .hzr_scope_current_vars(current, phase)) {
    return(NULL)
  }

  new_phases <- phases
  # As in .hzr_refit_with_scope(): a phase with no formula inherits the
  # global terms, and the candidate is added to them (#284).
  inherited <- if (is.null(new_phases[[phase]]$formula)) {
    .hzr_inherited_rhs(current)
  }
  new_phases[[phase]] <- .hzr_phase_update_formula(
    new_phases[[phase]], action = "add", var = term, inherited = inherited
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
      # The fit's own construction, so an intercept-free formula builds as
      # it did at fit time (#303).
      x_j <- .hzr_formula_design(ph$formula, data)$x
      x_list[[nm]] <- x_j
      cov_counts[[nm]] <- ncol(x_j)
    } else if (!is.null(current$fit$x_list[[nm]]) || !is.null(d$x)) {
      # A phase with no formula uses the design the fit used, which under
      # `time_windows` is the window-expanded one; `d$x` is the plain global
      # matrix and would change the phase's columns (#284).
      xm <- current$fit$x_list[[nm]] %||% d$x
      x_list[[nm]] <- xm
      cov_counts[[nm]] <- ncol(xm)
    } else {
      x_list[[nm]] <- NULL
      cov_counts[[nm]] <- 0L
    }
  }

  if (!is.null(x_list[[phase]]) &&
        .hzr_score_duplicate_columns(x_list[[phase]])) {
    return(list(reason = "duplicate_column"))
  }

  old_counts <- current$fit$covariate_counts
  if (is.null(old_counts) ||
        cov_counts[[phase]] != old_counts[[phase]] + 1L) {
    # A term that expands to more than one column (factor, spline) breaks the
    # one-zero-per-add layout; the refit path guards this upstream too.
    return(NULL)
  }
  unchanged <- setdiff(nms, phase)
  if (any(cov_counts[unchanged] != old_counts[unchanged])) {
    # Every phase the step does not touch must keep the design the fit used;
    # a different column count means the expansion rebuilt it wrongly.
    return(NULL)
  }
  for (nm in unchanged) {
    # The same count can still be a different design: a fit saved before
    # #303 holds `~ 0 + o`, for an ordered `o`, as dummies `om`, `oh`, which
    # now rebuild as `o.L`, `o.Q`.
    if (cov_counts[[nm]] > 0L &&
          !identical(colnames(x_list[[nm]]),
                     colnames(current$fit$x_list[[nm]]))) {
      return(NULL)
    }
  }
  for (nm in nms) {
    xm <- x_list[[nm]]
    if (is.null(xm) || ncol(xm) == 0L) next
    if (nrow(xm) != n_time || anyNA(xm)) return(NULL)
  }

  # The candidate's column is not always the phase's last: model.matrix()
  # puts main effects before interactions, so a phase with `age * mal`
  # gains `opmos` before `age:mal`. Find it by name among the phase's
  # columns. A phase with no columns yet takes the only slot; one whose
  # columns cannot be matched by name declines (a score of NA) rather than
  # pin the zero on whichever coefficient happens to be last.
  old_cols <- colnames(current$fit$x_list[[phase]])
  new_cols <- colnames(x_list[[phase]])
  new_col_pos <- cov_counts[[phase]]
  if (old_counts[[phase]] > 0L) {
    hit <- if (!is.null(new_cols) &&
                 length(old_cols) == old_counts[[phase]]) {
      which(!new_cols %in% old_cols)
    } else {
      integer()
    }
    if (length(hit) != 1L) return(NULL)
    new_col_pos <- hit
  }

  parts <- .hzr_split_theta(current$fit$theta, phases, old_counts)

  theta_new <- numeric(0)
  theta_idx <- integer(0)
  beta_idx <- NA_integer_
  pos <- 0L
  for (nm in nms) {
    part <- parts[[nm]]
    if (identical(nm, phase)) {
      # Old entries before the new slot: log_mu, the shapes, and the betas
      # whose columns precede the candidate's. The pinned zero goes there,
      # so every other coefficient keeps its own column.
      at <- length(part) - old_counts[[phase]] + new_col_pos - 1L
      after <- length(part) - at
      theta_new <- c(theta_new, part[seq_len(at)], 0,
                     part[at + seq_len(after)])
      theta_idx <- c(theta_idx, pos + seq_len(at),
                     pos + at + 1L + seq_len(after))
      beta_idx <- pos + at + 1L
      pos <- pos + length(part) + 1L
    } else {
      theta_new <- c(theta_new, part)
      theta_idx <- c(theta_idx, pos + seq_along(part))
      pos <- pos + length(part)
    }
  }

  # Shared with the optimizer's own naming (and with hzr_theta_names()). This
  # path previously fell through to character(0) when a covariate matrix had
  # no colnames, producing FEWER names than parameters -- and `names(x) <-`
  # pads with NA rather than erroring, so the misalignment was silent.
  names(theta_new) <- .hzr_theta_names_list(
    new_phases, .hzr_phase_cov_names(new_phases, cov_counts, x_list))

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
#' its coefficient takes the last slot of theta, exactly where the refit's
#' `c(theta_old, 0)` warm start puts it. Here it stays pinned at zero.
#'
#' @noRd
.hzr_score_expand_single <- function(current, var, phase, data,
                                     term = var) {
  if (!is.null(phase)) {
    return(NULL)
  }
  # By the candidate's term label, never its column name: the model's terms
  # are labels, and the column `age:mal` is spelled like the interaction
  # `age:mal`, which this compared it with and declined (#449). `var` is
  # still the column the values are read from.
  if (term %in% .hzr_scope_current_vars(current)) {
    return(NULL)
  }

  d <- current$data
  xcand <- .hzr_candidate_numeric(data[[var]])
  if (is.null(xcand) || length(xcand) != length(d$time) || anyNA(xcand)) {
    return(NULL)
  }

  new_col <- matrix(as.numeric(xcand), ncol = 1L,
                    dimnames = list(NULL, term))
  x_new <- if (is.null(d$x)) new_col else cbind(d$x, new_col)
  # The refit builds its design with model.matrix(), which names a logical
  # column <var>TRUE. Check the name the refit would create, not `var`: a
  # logical `flag` beside factor `fla`'s dummy `flag` fits fine.
  # model.matrix() names the column by the term label, so a non-syntactic
  # `age:mal` is `` `age:mal` ``, not the interaction's `age:mal` (#449).
  refit_name <- if (is.logical(data[[var]])) paste0(term, "TRUE") else term
  if (refit_name %in% colnames(d$x)) {
    return(list(reason = "duplicate_column"))
  }

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
    g <- .hzr_score_try(
      fn(exp_$theta, time = d$time, status = d$status,
         time_lower = d$time_lower, time_upper = d$time_upper,
         x = exp_$x, weights = d$weights)
    )
    if (is.null(g) || length(g) != length(exp_$theta) || !all(is.finite(g))) {
      return(NULL)
    }
    return(as.numeric(g))
  }
  g <- .hzr_score_try(
    .hzr_gradient_multiphase(
      exp_$theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = exp_$phases,
      covariate_counts = exp_$covariate_counts,
      x_list = exp_$x_list,
      objective = .hzr_fit_objective(current)
    )
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
  h <- .hzr_score_try(
    .hzr_hessian_multiphase(
      exp_$theta, time = d$time, status = d$status,
      time_lower = d$time_lower, time_upper = d$time_upper,
      x = d$x, weights = d$weights,
      phases = exp_$phases,
      covariate_counts = exp_$covariate_counts,
      x_list = exp_$x_list
    )
  )
  if (!is.null(h) && is.matrix(h) && nrow(h) == length(exp_$theta)) {
    return(h)
  }
  .hzr_score_multiphase_hessian(
    current, exp_$theta, exp_$phases, exp_$covariate_counts, exp_$x_list
  )
}

#' Number of rows a fit was estimated on, read off the fit
#'
#' A multiphase fit aligns its phase designs by dropping every row with a
#' missing covariate, and stores the aligned designs in `$fit$x_list`, while
#' `$data$time` keeps the caller's full length. The designs' row count is
#' therefore the rows the fit used. A fit with no phase design columns has
#' nothing to drop, so its `time` length is the count.
#'
#' @param fit A fitted `hazard` object.
#' @return A single integer.
#' @keywords internal
#' @noRd
.hzr_fit_rows_used <- function(fit) {
  xl <- Filter(function(x) !is.null(x) && NCOL(x) > 0L, fit$fit$x_list)
  # The fitter refuses phase designs of different lengths before it stores
  # them, so the stored ones agree; none means nothing was dropped.
  if (length(xl)) NROW(xl[[1L]]) else length(fit$data$time)
}

#' Columns of `data` that differ from the fit's own data frame
#'
#' `$data$frame` is the frame hazard() was given, after its time-0 rows were
#' dropped (#374), so it holds the fit's rows in the fit's order. Every column
#' the two frames share must agree value for value; a column only `data` has
#' (a candidate derived after the fit) cannot be compared and is not. A fit
#' with no stored frame (the vector interface without `data =`) has nothing to
#' compare against and returns no columns.
#'
#' @return Character vector of the shared column names that differ.
#' @noRd
.hzr_score_rows_moved <- function(current, data) {
  frame <- current$data$frame
  if (!is.data.frame(frame)) return(character())
  if (nrow(frame) != nrow(data)) {
    # The caller has already matched nrow(data) to the fit's rows, so a frame
    # of another length is not those rows: on the vector interface `data` may
    # serve only to look names up, and is stored at its own length. It says
    # nothing about row order.
    return(character())
  }
  common <- intersect(names(data), names(frame))
  # Exact: within a tolerance, rows whose values differ only by rounding
  # could be swapped unseen (#515).
  same <- vapply(common, function(nm) {
    isTRUE(all.equal(data[[nm]], frame[[nm]], tolerance = 0,
                     check.attributes = FALSE))
  }, logical(1))
  common[!same]
}

#' The refusal for a `data` whose shared columns differ from the fit's frame
#'
#' @param moved Character vector from `.hzr_score_rows_moved()`.
#' @noRd
.hzr_rows_moved_message <- function(moved) {
  paste0(
    "`data` does not hold the rows the model was fitted on, in the same ",
    "order: column", if (length(moved) > 1L) "s", " ",
    paste0("`", utils::head(moved, 5L), "`", collapse = ", "),
    if (length(moved) > 5L) paste0(" and ", length(moved) - 5L, " more"),
    " differ", if (length(moved) == 1L) "s", " from the data frame given ",
    "to hazard(). Candidates are read from `data` row by row, so a sorted ",
    "or reordered frame scores every candidate against the wrong ",
    "observations. Pass the data frame the model was fitted on, in its ",
    "original row order."
  )
}

#' The per-row inputs a fit stores, other than its event times
#'
#' Everything the likelihood reads row by row: the status, the interval
#' bounds, the weights, and the covariate design columns (the global `x`, and
#' each multiphase phase's design). An input with one value on every row is
#' left out, since rows cannot be mismatched on it.
#'
#' @return A list named for messages, one element per input:
#'   `list(value, column)`, where `value` is the numeric vector and `column`
#'   the name of the column it came from: the design column's name for a
#'   covariate, and for the status, bounds and weights the column the stored
#'   call names (`NA` when the call computes them). An input stored at
#'   another length than `time` (a phase design after rows with a missing
#'   covariate were dropped) keeps that length, so the caller counts it as
#'   unmatched.
#' @noRd
.hzr_fit_row_inputs <- function(current) {
  d <- current$data
  call <- current$call
  input <- function(v, column = NA_character_) {
    if (is.null(v)) NULL else list(value = as.numeric(v), column = column)
  }
  cols <- function(m, label) {
    if (is.null(m) || NCOL(m) == 0L) return(list())
    m <- as.matrix(m)
    nm <- colnames(m) %||% paste0("V", seq_len(ncol(m)))
    stats::setNames(lapply(seq_len(ncol(m)), function(j) input(m[, j], nm[j])),
                    paste0(label, " `", nm, "`"))
  }
  xl <- current$fit$x_list
  inputs <- c(
    list(status = input(d$status, .hzr_status_column(call)),
         `time_lower` = input(d$time_lower, .hzr_call_column(call$time_lower)),
         `time_upper` = input(d$time_upper, .hzr_call_column(call$time_upper)),
         weights = input(d$weights, .hzr_call_column(call$weights))),
    cols(d$x, "covariate"),
    unlist(lapply(names(xl), function(ph) {
      cols(xl[[ph]], paste0("phase `", ph, "` covariate"))
    }), recursive = FALSE)
  )
  inputs <- Filter(Negate(is.null), inputs)
  Filter(function(inp) length(unique(inp$value)) > 1L, inputs)
}

#' The column of `data` a stored call argument names
#'
#' `status = dead` names `dead`, as do `status = d$dead` and
#' `status = d[["dead"]]`. Anything computed names no column, and nor does
#' `..2`, which is what a call forwarded through a wrapper's `...` records.
#'
#' @param expr An argument expression from the stored call, or `NULL`.
#' @return A single string, `NA` when no column is named.
#' @noRd
.hzr_call_column <- function(expr) {
  # `..2` is how a call made through a wrapper's `...` records an argument:
  # it names the wrapper's argument, not a column.
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    return(if (grepl("^\\.\\.(\\.|[0-9]+)$", nm)) NA_character_ else nm)
  }
  if (is.call(expr) && length(expr) == 3L) {
    fn <- expr[[1L]]
    key <- expr[[3L]]
    if (identical(fn, as.name("$")) && (is.symbol(key) || is.character(key))) {
      return(as.character(key))
    }
    if (identical(fn, as.name("[[")) && is.character(key) && length(key) == 1L) {
      return(key)
    }
  }
  NA_character_
}

#' The column of `data` a fit's status came from
#'
#' The vector interface's `status =`, or the event of a two-argument
#' `Surv(time, event)` on the formula interface's left-hand side.
#'
#' @param call The fit's stored call.
#' @return A single string, `NA` when no column is named.
#' @noRd
.hzr_status_column <- function(call) {
  if (!is.null(call$status)) return(.hzr_call_column(call$status))
  f <- call$formula
  if (!is.call(f) || !identical(f[[1L]], as.name("~")) || length(f) != 3L) {
    return(NA_character_)
  }
  lhs <- f[[2L]]
  surv <- is.call(lhs) &&
    deparse(lhs[[1L]]) %in% c("Surv", "survival::Surv") &&
    length(lhs) == 3L && is.null(names(lhs))
  if (surv) .hzr_call_column(lhs[[3L]]) else NA_character_
}

#' Check row order within tied event times against the fit's other inputs
#'
#' A column of `data` equal to the fit's `time` fixes the order of rows except
#' within a tie. Each other per-row input the fit stores is then looked for in
#' `data` under the column it came from, and only there:
#' - in order: it settles the rows it tells apart;
#' - the same (time, value) pairs in another order: the rows were reordered
#'   within a tie, and the screen is refused;
#' - otherwise, or with no such column: it cannot be checked.
#' A column found only by its values is not accepted, as any column that
#' happens to hold them would then vouch for the order.
#'
#' The order is proven when the inputs found tell every row apart, or when
#' every input is found. In the second case rows can differ in position only
#' where they are identical in everything the likelihood reads, so the
#' screen's answer is unchanged (#515).
#'
#' @param time The fit's stored event times.
#' @return `list(proven, moved, unmatched)`: `moved` is `NULL` or
#'   `c(column, input)`; `unmatched` names the inputs not found.
#' @noRd
.hzr_check_rows_within_ties <- function(current, data, time) {
  n <- length(time)
  found <- list(time = time)
  unmatched <- character()
  inputs <- .hzr_fit_row_inputs(current)
  for (nm in names(inputs)) {
    v <- inputs[[nm]]$value
    own <- inputs[[nm]]$column
    label <- if (is.na(own)) nm else paste0(nm, " (column `", own, "`)")
    col <- if (!is.na(own) && own %in% names(data)) data[[own]]
    if (length(v) != n || !(is.numeric(col) || is.logical(col)) ||
          length(col) != n) {
      unmatched <- c(unmatched, label)
      next
    }
    col <- as.numeric(col)
    # Exact: rows that may be swapped must be identical in the input, not
    # merely within a tolerance relative to the whole column.
    if (identical(col, v)) {
      found[[nm]] <- v
    } else if (identical(col[order(time, col)], v[order(time, v)])) {
      return(list(proven = FALSE, moved = c(own, nm), unmatched = unmatched))
    } else {
      unmatched <- c(unmatched, label)
    }
  }
  proven <- !length(unmatched) || anyDuplicated(as.data.frame(found)) == 0L
  list(proven = proven, moved = NULL, unmatched = unmatched)
}

#' Check `data`'s row order once per screen
#'
#' When the refits pair per-row vectors stored on the fit with `data` (the
#' vector interface's response, or `weights` on either interface), a
#' comparable stored frame is compared here for every criterion. An
#' unweighted formula fit rebuilds every per-row input from `data`, so there
#' only the score test needs it, and `.hzr_score_q()` does it.
#'
#' `.hzr_score_rows_moved()` compares against `$data$frame`, and skips when
#' there is none (the vector interface without `data =`) or when it has
#' another row count (a lookup-only `data`). Those skips let a reordered
#' `data` through silently (#487). What the fit always stores is its response,
#' so on that route a column of `data` holding the fit's `time` values is
#' compared with it: in the fit's order, the rows are taken as aligned; as the
#' same values in another order, the screen is refused. With ties in the
#' times, `.hzr_check_rows_within_ties()` checks the fit's other per-row
#' inputs as well (#515). With no such column, or ties it cannot resolve, the
#' order cannot be checked, and a classed `hzr_score_rows_unverified`
#' warning says so. Called once per screen, on the base fit only: every later
#' step's fit is a refit on `data` itself.
#'
#' @return `NULL`, invisibly; called for its error or warning.
#' @noRd
.hzr_check_data_row_order <- function(current, data, score = TRUE) {
  # Every comparison below reads `data` by column name, and of duplicated
  # names `data[[name]]` reads only the first: a second column of that name,
  # holding the rows in another order, would never be seen (#515).
  dup <- unique(names(data)[duplicated(names(data))])
  if (length(dup)) {
    stop(
      "`data` has more than one column named ",
      paste0("`", utils::head(dup, 5L), "`", collapse = ", "),
      if (length(dup) > 5L) paste0(" and ", length(dup) - 5L, " more"),
      ". hzr_stepwise() reads columns by name, so it cannot tell which one ",
      "is meant, nor check that `data` holds the rows the model was fitted ",
      "on in the same order. Give every column of `data` a unique name.",
      call. = FALSE
    )
  }
  frame <- current$data$frame
  # A `data` that shares no column with the frame (derived candidates only)
  # leaves nothing to compare: that is no proof of order, and falls through
  # to the checks below, which warn when they cannot prove it either.
  shares <- is.data.frame(frame) &&
    length(intersect(names(data), names(frame))) > 0L
  why_frame <- NULL
  if (shares && nrow(frame) == nrow(data)) {
    # .hzr_score_q() compares `data` with the frame, but only the score test
    # calls it. Every other criterion refits through .hzr_refit_with_scope(),
    # and that pairs `data`, read by position, with per-row vectors stored
    # on the fit: `time`, `status`, `time_lower` and `time_upper` on the
    # vector interface (which it keys on the same `call$formula`), and
    # `weights` on both. Its other stored inputs, `time_windows`, the phase
    # specs and the objective, are not per row. A formula fit with no stored
    # weights rebuilds every per-row input from `data`, so it alone is
    # consistent with any row order.
    stored_rows <- is.null(current$call$formula) ||
      !is.null(current$data$weights)
    if (stored_rows) {
      moved <- .hzr_score_rows_moved(current, data)
      if (length(moved)) stop(.hzr_rows_moved_message(moved), call. = FALSE)
    }
    # Equal shared columns prove the order only if they tell every row
    # apart: rows reordered within a group of duplicates leave them
    # identical while a derived candidate moves. Then the match is no proof,
    # and the time check below is tried instead. (The score test's own
    # comparison in .hzr_score_q() has the same limit, and relies on this.)
    if (!stored_rows && !score) return(invisible(NULL))
    shared <- frame[match(intersect(names(data), names(frame)), names(frame))]
    n_dup <- sum(duplicated(shared))
    if (n_dup == 0L) return(invisible(NULL))
    why_frame <- paste0(
      "the columns `data` shares with the data frame stored with the fit ",
      "have duplicate rows (", n_dup, " of ", nrow(shared), " repeat an ",
      "earlier row), so rows reordered among duplicates leave them unchanged"
    )
  }
  time <- current$data$time
  why <- if (!is.null(why_frame)) {
    why_frame
  } else if (is.data.frame(frame) && !shares) {
    "`data` shares no column with the data frame stored with the fit"
  } else if (is.data.frame(frame)) {
    paste0("the data frame stored with the fit has ", nrow(frame),
           " rows, not the fit's ", length(time))
  } else {
    "the fit was made without `data =`, so it stores no data frame"
  }
  # With ties in the times, a column holding them in order proves nothing
  # about the rows WITHIN a tie: reordering those leaves it identical. The
  # fit's other per-row inputs are then checked as well (#515).
  tied <- anyDuplicated(time) > 0L
  permuted <- character()
  in_order_tied <- character()
  for (nm in names(data)) {
    col <- data[[nm]]
    if (!is.numeric(col) || length(col) != length(time)) next
    # Exact, as ties are: within all.equal()'s tolerance, rows whose times
    # differ only by rounding would be swapped unseen, and taken as ordered.
    if (identical(as.numeric(col), as.numeric(time))) {
      if (!tied) return(invisible(NULL))
      in_order_tied <- c(in_order_tied, nm)
      next
    }
    # Exact as well: a column merely close to the times is not them, and
    # must not be refused as the times in another order.
    if (identical(sort(as.numeric(col)), sort(as.numeric(time)))) {
      permuted <- c(permuted, nm)
    }
  }
  if (length(permuted)) {
    stop(
      "`data` does not hold the rows the model was fitted on, in the same ",
      "order: column `", permuted[1L], "` holds the fit's event times in ",
      "another order. The score test reads each candidate from `data` row ",
      "by row, so a sorted or reordered frame scores every candidate against ",
      "the wrong observations. Pass the rows in the order the model was ",
      "fitted on.",
      call. = FALSE
    )
  }
  if (length(in_order_tied)) {
    within <- .hzr_check_rows_within_ties(current, data, time)
    if (!is.null(within$moved)) {
      stop(
        "`data` does not hold the rows the model was fitted on, in the same ",
        "order: column `", within$moved[1L], "` holds the fit's ",
        within$moved[2L], " in another order within tied event times. The ",
        "screen reads each candidate from `data` row by row, so rows ",
        "reordered within a tie are scored against the wrong observations. ",
        "Pass the rows in the order the model was fitted on.",
        call. = FALSE
      )
    }
    if (within$proven) return(invisible(NULL))
  }
  held <- if (length(in_order_tied)) {
    paste0(
      "; column `", in_order_tied[1L], "` holds the fit's event times in ",
      "order, but those times have ties (", length(unique(time)), " distinct ",
      "values in ", length(time), " rows), and the tied rows cannot be told ",
      "apart without the fit's ", paste(within$unmatched, collapse = ", "),
      ", which `data` does not hold in order under the column it came from"
    )
  } else {
    ", and no column of `data` holds the fit's event times"
  }
  warning(warningCondition(
    paste0(
      "The row order of `data` could not be checked against the fit: ", why,
      held, ". The score ",
      "test reads each candidate from `data` row by row, so if its rows are ",
      "not in the order the model was fitted on, every candidate is scored ",
      "against the wrong observations. ",
      if (length(in_order_tied)) {
        paste0("Add those to `data`, in the order the model was fitted on ",
               "(an input computed in the call, such as `status = x > 0`, ",
               "has no column to look under), to have it checked.")
      } else if (is.data.frame(frame)) {
        paste0("Include in `data` enough of the columns given to hazard() ",
               "to tell every row apart to have it checked.")
      } else {
        "Refit with `data =` to have it checked."
      }
    ),
    class = "hzr_score_rows_unverified"
  ))
  invisible(NULL)
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
#' @param term The candidate's term label, written into a multiphase phase
#'   formula. Defaults to `var`, which is right only for a syntactic name;
#'   the stepwise step passes the resolved label (#449). `var` may be `NA`
#'   for a candidate that is no column, which is declined as
#'   `not_single_column`.
#' @return `list(stat, df, p_value)`. `stat`/`p_value` are `NA_real_` for a
#'   degenerate candidate, a collinear candidate, or an unusable nuisance block.
#' @noRd
.hzr_score_q <- function(current, var, phase = NULL, data,
                         nuisance = NULL, term = var) {
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
  # Count the rows the FIT used, read off the fit itself: a multiphase fit
  # drops every row with a missing phase covariate but keeps the caller's full
  # `time` in `$data`, so `length(time)` overstated them, this check passed,
  # and every candidate was then labelled `not_expandable` (#372).
  n_time <- length(current$data$time)
  n_obs <- .hzr_fit_rows_used(current)
  if (n_obs != n_time) {
    stop(
      "The base fit dropped ", n_time - n_obs, " rows whose covariate ",
      "values were missing (NA or NaN), so its stored response (", n_time,
      " rows) no longer lines up with the rows it was fitted on (", n_obs,
      "), and no candidate can be scored against it. The values can be ",
      "missing in the data, or made missing by a transform in a model ",
      "formula, such as sqrt() or log() of a negative value, which ",
      "`na.omit()` on the data does not catch. Refit the base model on only ",
      "the rows it used, and pass that same data frame.",
      call. = FALSE
    )
  }
  if (nrow(data) != n_obs) {
    stop(
      "`data` has ", nrow(data), " rows but the fitted model used ", n_obs,
      ". The score test needs `data` row-aligned with the fit; pass the same ",
      "data frame the model was fitted on (after any NA removal, and without ",
      "rows at time 0, which hazard() drops).",
      call. = FALSE
    )
  }
  # The row COUNT matching is not the rows matching. Candidate values are read
  # from `data` by position and scored against the fit's stored rows, so the
  # same rows in another order scored every candidate against the wrong
  # patients and entered a different variable, with no warning (#487).
  moved <- .hzr_score_rows_moved(current, data)
  if (length(moved)) stop(.hzr_rows_moved_message(moved), call. = FALSE)

  # A term that is no column (an interaction, a transform) is not a candidate
  # the score can test; saying `non_numeric` described a column (#449).
  if (is.na(var)) return(na_result("not_single_column"))
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

  exp_ <- .hzr_score_expand(current, var, phase, data, term = term)
  if (is.null(exp_)) return(na_result("not_expandable"))
  if (!is.null(exp_$reason)) return(na_result(exp_$reason))

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


# Reasons a candidate can be rescued by refitting it and testing by Wald.
#
# The first two mean "the quadratic approximation at beta = 0 broke down",
# which is what a LARGE true effect looks like -- so declining them is exactly
# backwards and a refit gives the right answer. SAS's own q1.c says as much
# ("IT IS POSSIBLE THAT THE PROGRAM WILL RETURN A NEGATIVE Q VALUE ... THE USER
# SHOULD USE THE MORE EXPENSIVE Q2 AS AN ALTERNATIVE"); Q2 is named once in the
# C tree and never implemented, and dqstat.c instead declines the candidate
# with p = 1. This is that unbuilt alternative.
#
# The other two are faults of the score test's inputs, not of the candidate
# (#570):
#   information_nonpositive  the candidate's own observed information at
#                            beta = 0 is not positive, the same breakdown one
#                            step earlier in q1.c (its flag 2 against flag 3);
#   nuisance_singular        the CURRENT model's information block could not
#                            be inverted, so no candidate can be adjusted for
#                            it. That is a base fit on a ridge, and it takes
#                            every candidate at the step with it.
# A refit of the extended model has its own Hessian and tests the candidate.
# An earlier version of this comment listed nuisance_singular among the
# reasons "no refit can make testable"; measured for #565 on a base fit up a
# ridge, the refit gave beta 0.847 (se 0.057) against 0.841 (0.057) from an
# independent fit. When the refit's own Hessian fails too, the row is recorded as
# `fallback_no_variance` or as a refit failure, never as tested.
#
# Still kept narrow. The reasons that describe the candidate's COLUMN, or a
# design hazard() would refuse -- collinear, constant, non_numeric,
# not_single_column, duplicate_column, not_expandable -- are NOT here: no
# refit can make those candidates testable, and paying one per degenerate
# candidate would give back the whole speed advantage the score criterion
# exists for. no_information and nonfinite are not here either. They report
# a gradient or information that could not be computed for the expanded
# model, and whether a refit rescues them has not been measured (#577).
#
# The cost: nuisance_singular refits every candidate at its step that reaches
# the nuisance check, which is what `criterion = "wald"` pays at every step.
# .hzr_score_q() returns that reason BEFORE it expands the candidate, so at
# such a step a collinear, duplicate_column or not_expandable candidate is
# refitted too, and its refit fails or yields no variance; only constant,
# non_numeric and not_single_column are screened out ahead of it.
.hzr_score_fallback_reasons <- c("information_indefinite",
                                 "coefficient_diverging",
                                 "nuisance_singular",
                                 "information_nonpositive")

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
      "degenerate ones. `criterion = \"score\"` refits and Wald-tests them",
      "itself, so reaching this reason means that refit ERRORED or did not",
      "converge -- see `refit_failures`. A refit that converged but could not",
      "be Wald-tested reports `fallback_no_variance` instead. It is also",
      "reached when the current fit has not truly converged"
    ),
    fallback_no_variance = paste(
      "the score could not test the candidate, and the Wald refit that would",
      "have rescued it converged but produced no usable variance -- its",
      "Hessian was not invertible, so the coefficient has an estimate but no",
      "standard error. The candidate was therefore tested by NEITHER",
      "criterion. It is typically a STRONG one, since that is what drives the",
      "score's information indefinite in the first place, so a screen",
      "reporting this has understated that variable rather than declined it"
    ),
    information_nonpositive = paste(
      "the candidate's own observed information was not positive, before any",
      "adjustment for the current model. This is a different fault from",
      "collinearity: the candidate is a poor one in itself, or the fit it",
      "would be added to is not at a maximum. `criterion = \"score\"` refits",
      "and Wald-tests such a candidate itself, so reaching this reason means",
      "that refit errored or did not converge -- see `refit_failures`"
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
    not_single_column = paste(
      "the candidate is a term, such as an interaction, and not a single",
      "column of `data`, which is all the score criterion can test.",
      "`criterion = \"wald\"` refits it instead"
    ),
    nuisance_singular = paste(
      "the current model's information matrix could not be inverted, so no",
      "candidate could be scored at that step. `criterion = \"score\"` refits",
      "and Wald-tests each of them itself, so reaching this reason means that",
      "refit errored or did not converge -- see `refit_failures`"
    ),
    no_information = "no observed information was available for the candidate",
    not_expandable = "the candidate could not be added to the model",
    duplicate_column = paste(
      "the candidate's design column has the same name as one already in the",
      "model, and hazard() refuses a design whose column names repeat. A",
      "factor's dummy columns are named <factor><level>, so a numeric `gb`",
      "collides with factor `g`'s level `b`: rename the column, or rename or",
      "relevel the factor"
    ),
    nonfinite = "the score or its variance was not finite",
    rows_differ = paste(
      "under `criterion = \"aic\"`, the candidate's refit was fitted on",
      "different rows from the current model: a multiphase fit drops every",
      "row where a covariate is missing, so the candidate's log-likelihood",
      "summed fewer rows and its AIC could not be compared. Remove or impute",
      "the missing values before the screen, so every model uses the same rows"
    ),
    loglik_below_base = paste(
      "the candidate's refit ended with a log-likelihood below the current",
      "model's, although the candidate model contains the current one. That",
      "cannot happen at the optimum, so the refit did not converge, and",
      "neither its AIC nor its Wald test describes a fitted model. More",
      "starting points (`control$n_starts`) or iterations (`control$maxit`)",
      "may let it converge"
    ),
    wald_no_variance = paste(
      "the model had no usable variance for the coefficient, so its Wald",
      "test could not be computed: there was no variance matrix, the",
      "coefficient's variance was not positive, or a multi-column term's",
      "variance block was singular. A fit with interval- or left-censored",
      "rows takes its variance from numDeriv, so a screen run without",
      "numDeriv installed reports this for every variable"
    )
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
