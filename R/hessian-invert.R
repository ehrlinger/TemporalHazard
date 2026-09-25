#' @keywords internal
NULL

# hessian-invert.R -- Numerically stable Hessian inversion with diagnostics.

# Reciprocal-condition warning threshold (~1.5e-8). Shared by .hzr_safe_solve()
# and the summary() diagnostic note so they never drift apart.
.hzr_rcond_tol <- .Machine$double.eps^0.5

#' Stable Hessian inversion with conditioning diagnostics
#'
#' Inverts a negative-log-likelihood Hessian into a variance-covariance
#' matrix, hardening against the ill-conditioning that arises at high
#' parameter counts.  Symmetrizes the input, checks the reciprocal condition
#' number, inverts via Cholesky (with a \code{solve()} fallback for non-PD
#' Hessians), and guards non-positive variances.  Emits a named warning for
#' every degenerate path.
#'
#' @param H Square numeric Hessian of the negative log-likelihood.
#' @param tol Reciprocal-condition warning threshold.
#' @return A list with:
#'   \code{vcov} (the variance-covariance matrix, or \code{NA} on failure;
#'   diagonals with non-positive variance, and their rows/cols, set to
#'   \code{NA}); \code{rcond} (reciprocal condition number of the symmetrized
#'   Hessian, \code{NA} if unavailable); \code{pd} (\code{TRUE} if the Hessian
#'   was positive-definite, \code{FALSE} if inverted via fallback, \code{NA}
#'   if not invertible); \code{reason} (why \code{vcov} is \code{NA}:
#'   \code{"Hessian has non-finite entries"} or \code{"Hessian not
#'   invertible"}; \code{NA_character_} when a matrix was returned).
#' @noRd
.hzr_safe_solve <- function(H, tol = .hzr_rcond_tol) {
  # (1) Non-finite / non-matrix guard
  if (is.null(H) || !is.matrix(H) || anyNA(H) || any(!is.finite(H))) {
    warning("Hessian contains non-finite entries; standard errors unavailable")
    return(list(vcov = NA, rcond = NA_real_, pd = NA,
                reason = "Hessian has non-finite entries"))
  }

  # (2) Symmetrize (numDeriv Hessians are only symmetric to Richardson tol)
  H <- (H + t(H)) / 2

  # (3) Conditioning check
  rc <- tryCatch(rcond(H), error = function(e) NA_real_)
  if (is.na(rc)) {
    warning("Hessian conditioning could not be assessed; standard errors may be unreliable")
  } else if (rc < tol) {
    warning(sprintf(
      "Hessian is ill-conditioned (rcond = %.3g); standard errors may be unreliable",
      rc
    ))
  }

  # (4) Stable inversion: Cholesky (PD) with solve() fallback (non-PD)
  ch <- tryCatch(chol(H), error = function(e) NULL)
  if (!is.null(ch)) {
    pd <- TRUE
    vcov <- chol2inv(ch)
  } else {
    pd <- FALSE
    vcov <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(vcov)) {
      warning("Hessian not invertible; standard errors unavailable")
      return(list(vcov = NA, rcond = rc, pd = NA,
                  reason = "Hessian not invertible"))
    }
    warning("Hessian is not positive-definite at the optimum; standard errors may be unreliable")
  }
  dimnames(vcov) <- dimnames(H)  # chol2inv() drops names

  # (5) Non-positive-variance guard
  d <- diag(vcov)
  bad <- !is.finite(d) | d <= 0
  if (any(bad)) {
    warning("Non-positive variance estimates; the optimum may not be a proper maximum")
    vcov[bad, ] <- NA_real_
    vcov[, bad] <- NA_real_
  }

  list(vcov = vcov, rcond = rc, pd = pd, reason = NA_character_)
}


# Minimum |correlation| between two estimates before their trade-off counts as
# a ridge rather than ordinary imprecision. Kept next to .hzr_rcond_tol so the
# two thresholds that gate the weak-identification note stay together.
.hzr_ridge_cor_tol <- 0.99

#' Identify a ridge (weak-identification) direction in a fitted model
#'
#' An ill-conditioned Hessian makes standard errors unreliable.  When the
#' ill-conditioning is a \emph{ridge} the point estimates are unreliable too:
#' the likelihood is near-flat along some combination of parameters, so the
#' data pin down only that combination and not the individual values.  This
#' locates the flat direction and names the parameters spanning it.
#'
#' The direction is taken from the \emph{correlation} of the estimates, not
#' from the raw covariance.  Parameters in this package sit on very different
#' scales (an \code{m} of 27 against a \code{nu} of 0.027), and in raw units a
#' direction that moves both equally in statistical terms loads almost
#' entirely on the larger one, so a two-parameter ridge would be reported as a
#' single unidentified parameter.  Standardising first makes the loadings comparable.
#'
#' The pairwise reading cannot see a single parameter the data do not
#' determine, since a correlation matrix normalises its variance away. When
#' it finds nothing and \code{theta} is supplied, a second reading names one
#' parameter that carries the flattest direction on its own (#415); see
#' \code{.hzr_weak_single_parameter()}. Behind the rcond gate a parameter
#' that dominates that direction is named whatever its own standard error.
#'
#' Where two eigenvalues are exactly tied, LAPACK may return any basis of the
#' degenerate eigenspace, so the loadings can spread across both ridges and
#' the reported parameter set be their union.  \code{n_directions} is the
#' honest signal in that case: treat a count above one as "look at the whole
#' coefficient table", not as a precise inventory of two separate ridges.
#'
#' @param vcov Variance-covariance matrix of the estimates.  Rows/columns for
#'   parameters that were not estimated (\code{NA} diagonals) are dropped.
#' @param rcond Reciprocal condition number of the Hessian, as returned by
#'   \code{.hzr_safe_solve()}.
#' @param param_names Optional character vector of parameter names, the same
#'   length as \code{nrow(vcov)}.  Positional labels are used when absent.
#' @param tol Gate on \code{rcond}.  Defaults to the package-wide
#'   ill-conditioning threshold so this never fires where the rcond warning
#'   stays silent.
#' @param cor_tol Minimum absolute correlation for a trade-off to count.
#' @param share Cumulative squared-loading share used to decide how many
#'   parameters span the flat direction.
#' @param theta The fitted parameter vector, aligned with \code{vcov}. Needed
#'   only by the single-parameter reading, which is skipped without it.
#' @param shape_names Names of the g3 shape parameters (\code{gamma},
#'   \code{alpha}, \code{eta}) the single-parameter reading may name; it is
#'   skipped when empty.
#' @return One of three values, which callers must keep distinct:
#'   \code{NULL} when the check ran and found no ridge; \code{NA} when the
#'   check \emph{could not} run, because no usable Hessian was available; and
#'   otherwise a list with \code{params} (names spanning the flat direction),
#'   \code{weights} (their squared loadings), \code{correlation} (the
#'   strongest pairwise correlation among them), \code{rcond}, and
#'   \code{n_directions} (how many independent near-flat directions cleared
#'   the gate, of which this is the flattest).
#'
#'   Collapsing \code{NA} into \code{NULL} is the defect this separation
#'   exists to prevent: \code{numDeriv} is a \code{Suggests} and the
#'   analytic Hessian declines for left- and interval-censored rows by
#'   design, so "no Hessian" is reachable on a real install, and reporting
#'   it as "well identified" is a result-shaped answer over a computation
#'   that never happened.
#'
#'   Directions are examined from flattest to stiffest and the flattest
#'   genuine trade-off is reported, so a ridge is still found when some larger
#'   block of moderately correlated parameters carries more variance than it
#'   does.
#'
#'   The \code{_impl} returns \code{list(weak, reason)}, where \code{reason}
#'   names why \code{weak} is \code{NA}; \code{.hzr_weak_direction()} returns
#'   \code{weak} alone, unchanged, for every existing caller.
#' @noRd
.hzr_weak_direction_impl <- function(vcov, rcond, param_names = NULL,
                                     tol = .hzr_rcond_tol,
                                     cor_tol = .hzr_ridge_cor_tol,
                                     share = 0.9, theta = NULL,
                                     shape_names = NULL) {
  na_because <- function(reason) list(weak = NA, reason = reason)
  looked <- function(weak) list(weak = weak, reason = NA_character_)

  # (1) Gate on the existing ill-conditioning threshold. Each exit below is
  #     either NA ("could not look") or NULL ("looked, nothing there"); see
  #     the @return note on why the two must not be merged.
  #
  #     An absent or unassessable rcond means no Hessian was obtained, so
  #     nothing was examined. A well-conditioned one is a real answer: the
  #     likelihood cannot be near-flat in any direction when it is.
  if (length(rcond) != 1L || is.na(rcond)) {
    return(na_because(if (is.matrix(vcov)) "Hessian condition number unavailable"
                      else "standard errors unavailable"))
  }
  if (rcond >= tol) return(looked(NULL))
  if (is.null(vcov) || !is.matrix(vcov)) return(na_because("standard errors unavailable"))
  # One parameter cannot trade off against another, so there is no ridge to
  # find. That is a conclusion, not a gap.
  if (nrow(vcov) < 2L) return(looked(NULL))

  # (2) Drop parameters that were not estimated (fixed params carry NA rows).
  d <- diag(vcov)
  keep <- which(is.finite(d) & d > 0)
  if (length(keep) < 2L) return(looked(NULL))
  V <- vcov[keep, keep, drop = FALSE]
  # A non-finite entry among parameters that *were* estimated is a gap: the
  # covariance exists but cannot be decomposed.
  if (anyNA(V) || any(!is.finite(V))) return(na_because("covariance has non-finite entries"))

  nms <- if (length(param_names) == nrow(vcov)) {
    as.character(param_names)[keep]
  } else {
    paste0("par", keep)
  }

  # (3) Standardise to a correlation matrix (see note above on scaling).
  s <- sqrt(diag(V))
  R <- V / outer(s, s)
  if (anyNA(R) || any(!is.finite(R))) return(na_because("covariance has non-finite entries"))
  # A negative eigenvalue means this is not a covariance: the Hessian was
  # taken where it is not negative definite, typically short of the optimum.
  # Reading it as a correlation matrix named a ridge from "correlation 1.09"
  # on an identified model (#416), so it is a gap, not a finding. A
  # "correlation" above 1 always leaves a negative eigenvalue, so this one
  # test covers it.
  # The tolerance is the numerical error scale of the eigenvalue computation, n * eps * max|lambda|,
  # not a fixed -sqrt(eps). A fixed cutoff of about -1.49e-08 admitted matrices
  # that are indefinite far beyond rounding error: c(1, 1 + 1e-8, 1 + 1e-8, 1) has a
  # minimum eigenvalue near -1e-08 and an off-diagonal of 1.00000001, and was
  # reported as a ridge rather than declined.
  #
  # The scale is taken from the CORRELATION matrix, which is why it is safe to
  # take it at all. R is already scale-free, so its spectrum is bounded by the
  # dimension and the tolerance cannot depend on the parameters' units. Taking
  # the scale from the covariance instead would make the same fit pass or fail
  # according to whether a time was recorded in days or years, which is exactly
  # what this check must not do.
  #
  # A real correlation matrix is positive semi-definite, so a genuine ridge
  # sits at or above zero and is unaffected: r = 0.99 through exactly 1, and
  # equicorrelated blocks at k = 3, 4 and 5, all have a minimum eigenvalue >= 0.
  e_all <- tryCatch(eigen(R, symmetric = TRUE, only.values = TRUE)$values,
                    error = function(e) NA_real_)
  e_min <- if (anyNA(e_all)) NA_real_ else min(e_all)
  e_tol <- if (anyNA(e_all)) NA_real_ else
    nrow(R) * .Machine$double.eps * max(abs(e_all))
  if (is.na(e_min) || e_min < -e_tol) {
    return(na_because("covariance is not positive definite"))
  }

  # (4) Scan the standardised directions from flattest to stiffest, rather
  #     than gating only the leading one. Taking just the top eigenvector
  #     misses a real ridge whenever a larger *block* of moderately
  #     correlated parameters outranks it: an equicorrelated block of k
  #     parameters has eigenvalue 1 + (k - 1) * r, which passes a perfect
  #     two-parameter ridge's ceiling of 2 as soon as r > 1 / (k - 1) --
  #     0.50 at k = 3, 0.33 at k = 4. The block is the higher-variance
  #     direction and is not a trade-off, so the gate below rejects it, and
  #     the ridge underneath it was never looked at. Returning NULL there
  #     reports "well identified" for a fit that is not, which is the exact
  #     failure this function exists to prevent.
  e <- tryCatch(eigen(R, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) return(na_because("eigendecomposition failed"))

  # eigen() returns values in decreasing order, so this walks flattest first
  # and the first direction that is a genuine trade-off wins. No separate
  # floor on the eigenvalue is needed to keep a stiff direction from being
  # reported: clearing cor_tol requires a strongly correlated pair among the
  # selected parameters, and such a pair always puts a high-variance
  # direction earlier in this same scan, so the flat one is returned first.
  #     The scan runs to the end rather than returning on the first hit. A
  #     second flat direction means a second set of parameters is unidentified
  #     too, and reporting only the first invites reading every parameter it
  #     does not name as identified. Counting them costs one pass over an
  #     already-computed eigendecomposition; the flattest is still what gets
  #     named, since it is the one the data constrain least.
  #
  #     Distinctness is measured on the parameter SET, not on the eigenvector.
  #     A ridge between a and b clears the gate twice -- once on the flat
  #     direction (a - b) and again on its stiff partner (a + b), since the
  #     gate tests the correlation among the selected parameters and both
  #     directions select the same pair. Counting eigenvectors would report a
  #     single ridge as two. Only the first direction over each set of
  #     parameters is counted.
  found <- NULL
  seen <- character(0)

  for (j in seq_along(e$values)) {
    w <- e$vectors[, j]^2

    ord <- order(w, decreasing = TRUE)
    k <- which(cumsum(w[ord]) >= share)[1]
    if (is.na(k)) k <- length(ord)
    idx <- ord[seq_len(k)]
    if (length(idx) < 2L) next

    # (5) Require a genuine trade-off. Without this an uncorrelated but
    #     imprecise parameter set (correlation matrix near identity, where
    #     the eigenvectors are arbitrary) would be reported as a ridge.
    sub <- R[idx, idx, drop = FALSE]
    off <- sub[upper.tri(sub)]
    if (!length(off)) next
    r_max <- off[which.max(abs(off))]
    if (abs(r_max) < cor_tol) next

    key <- paste(sort(idx), collapse = ",")
    if (key %in% seen) next
    seen <- c(seen, key)

    if (is.null(found)) {
      found <- list(params = nms[idx], weights = w[idx],
                    correlation = r_max, rcond = rcond)
    }
  }

  if (is.null(found)) {
    return(looked(.hzr_weak_single_parameter(V, theta, keep, nms, rcond,
                                             cor_tol, shape_names)))
  }
  found$n_directions <- length(seen)
  looked(found)
}

#' A flat direction along ONE parameter (#415)
#'
#' The pairwise scan above cannot see this by construction: a correlation
#' matrix normalises every variance to 1, so a parameter the data do not
#' determine, with no partner to trade off against, looks exactly like a
#' well-determined one. It is read instead from the covariance in a metric
#' where the positive g3 shapes (`gamma`, `alpha`, `eta`) are on the log
#' scale, so their variance is relative, and everything else is as it is
#' (already log-scale, or a coefficient). Only a g3 shape can be NAMED; see
#' the comment at the naming step for why. Scaling every parameter by its own value instead -- the
#' candidate measured in #415 -- named `log_mu` on fits whose `gamma` had run
#' to 1e7, because a log-scale estimate near 0 inflates its relative
#' variance without bound.
#'
#' Fires only when the pairwise scan found nothing (so every pair it names
#' is unchanged), behind the same rcond gate as the caller, when the flattest
#' direction loads on one parameter (`|loading| >= cor_tol`). The gate does
#' the discriminating, as #415 found: there is no magnitude bar, because a
#' likelihood that is flat all the way to a boundary (gamma -> Inf) still
#' reports a finite local curvature there, and a relative standard error of
#' 0.1 to 0.4 at gamma-hat = 1e7 is numerical, not information. Measured on
#' 34 fits (#415's design, FIXGE2, n = 300; identified fits at n = 1500; the
#' ill-scaled-covariate control): it fires on 12 of the 14 gated degenerate
#' fits it can read, naming gamma on 11 and alpha on 1, and on none of the
#' 10 identified fits, whose largest loading is 0.972 even ungated.
#'
#' @return `NULL`, or the `$weak` list with `single = TRUE`, one `params`
#'   entry, its `estimate`, and `se_metric`, its standard error in the
#'   metric above.
#' @noRd
.hzr_weak_single_parameter <- function(V, theta, keep, nms, rcond, cor_tol,
                                       shape_names) {
  if (is.null(theta) || length(theta) < max(keep)) return(NULL)
  if (!length(shape_names)) return(NULL)
  th <- unname(theta)[keep]
  if (anyNA(th) || any(!is.finite(th))) return(NULL)
  # The g3 shapes are named by the CALLER from the phase specs, not guessed
  # from a name suffix: a covariate called `gamma` in a cdf phase is
  # `early.gamma` too.
  positive_shape <- nms %in% shape_names & th > 0
  sc <- ifelse(positive_shape, th, 1)
  S <- V / outer(sc, sc)
  e <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e) || !is.finite(e$values[1])) return(NULL)
  v1 <- e$vectors[, 1]
  j <- which.max(abs(v1))
  # Only a positive g3 shape is named. A coefficient or a log-scale
  # parameter can dominate this direction through its UNITS alone: a
  # covariate recorded in units of 1e-4 opened the rcond gate and named its
  # well-identified coefficient (z = 4.7) as undetermined. Its flatness is
  # not distinguishable from scaling without a natural unit, and the g3
  # shapes are the parameters with one. (`nu` is excluded: it is signed, and
  # 0 is a legitimate limiting value, so relative variance near it is noise.)
  if (abs(v1[j]) < cor_tol || !positive_shape[j]) return(NULL)
  list(params = nms[j], weights = v1[j]^2, correlation = NA_real_,
       rcond = rcond, n_directions = 1L, single = TRUE,
       estimate = th[j], se_metric = sqrt(e$values[1]))
}

# The shape every existing caller relies on: list / NULL / NA, unchanged.
.hzr_weak_direction <- function(vcov, rcond, param_names = NULL,
                                tol = .hzr_rcond_tol,
                                cor_tol = .hzr_ridge_cor_tol,
                                share = 0.9, theta = NULL,
                                shape_names = NULL) {
  .hzr_weak_direction_impl(vcov, rcond, param_names = param_names,
                           tol = tol, cor_tol = cor_tol, share = share,
                           theta = theta, shape_names = shape_names)$weak
}

#' Warning text for a detected ridge direction
#'
#' Shared by the fit-time warning and the \code{summary()} note so the two
#' never drift apart.
#'
#' @param weak A list result from \code{.hzr_weak_direction()}.
#' @return A single string.
#' @noRd
.hzr_weak_direction_message <- function(weak) {
  if (isTRUE(weak$single)) {
    return(paste0(
      "weakly identified fit: '", weak$params, "' = ",
      format(weak$estimate, digits = 4), " carries a near-flat direction ",
      "on its own (loading ", format(sqrt(weak$weights), digits = 3),
      ", Hessian rcond = ", format(weak$rcond, digits = 3), "). Its point ",
      "estimate -- not just its standard error -- is not pinned down by the ",
      "data, and it may be running to a boundary of its range."
    ))
  }
  # Naming one direction and stopping invites the reader to treat every
  # parameter it does not mention as identified. When the scan cleared the
  # gate more than once, say so rather than letting the omission speak.
  more <- if (isTRUE(weak$n_directions > 1L)) {
    paste0(
      " ", weak$n_directions, " near-flat directions were found; this is ",
      "the flattest, and parameters it does not name may be unidentified ",
      "too."
    )
  } else {
    ""
  }
  paste0(
    "weakly identified fit: ",
    paste0("'", weak$params, "'", collapse = " and "),
    " are determined only in combination (correlation ",
    format(weak$correlation, digits = 3),
    ", Hessian rcond = ", format(weak$rcond, digits = 3),
    "). The likelihood is near-flat along that direction, so their ",
    "individual point estimates -- not just their standard errors -- are ",
    "not pinned down by the data.",
    more
  )
}
