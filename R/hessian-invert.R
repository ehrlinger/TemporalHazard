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
#'   if not invertible).
#' @noRd
.hzr_safe_solve <- function(H, tol = .hzr_rcond_tol) {
  # (1) Non-finite / non-matrix guard
  if (is.null(H) || !is.matrix(H) || anyNA(H) || any(!is.finite(H))) {
    warning("Hessian contains non-finite entries; standard errors unavailable")
    return(list(vcov = NA, rcond = NA_real_, pd = NA))
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
      return(list(vcov = NA, rcond = rc, pd = NA))
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

  list(vcov = vcov, rcond = rc, pd = pd)
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
#' entirely on the larger one -- reporting a two-parameter ridge as a single
#' unidentified parameter.  Standardising first makes the loadings comparable.
#'
#' A single imprecise but uncorrelated parameter is deliberately \emph{not}
#' reported: that is ordinary low precision, already covered by the
#' ill-conditioning warning and by the parameter's own standard error.
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
#' @noRd
.hzr_weak_direction <- function(vcov, rcond, param_names = NULL,
                                tol = .hzr_rcond_tol,
                                cor_tol = .hzr_ridge_cor_tol,
                                share = 0.9) {
  # (1) Gate on the existing ill-conditioning threshold. Each exit below is
  #     either NA ("could not look") or NULL ("looked, nothing there"); see
  #     the @return note on why the two must not be merged.
  #
  #     An absent or unassessable rcond means no Hessian was obtained, so
  #     nothing was examined. A well-conditioned one is a real answer: the
  #     likelihood cannot be near-flat in any direction when it is.
  if (length(rcond) != 1L || is.na(rcond)) return(NA)
  if (rcond >= tol) return(NULL)
  if (is.null(vcov) || !is.matrix(vcov)) return(NA)
  # One parameter cannot trade off against another, so there is no ridge to
  # find. That is a conclusion, not a gap.
  if (nrow(vcov) < 2L) return(NULL)

  # (2) Drop parameters that were not estimated (fixed params carry NA rows).
  d <- diag(vcov)
  keep <- which(is.finite(d) & d > 0)
  if (length(keep) < 2L) return(NULL)
  V <- vcov[keep, keep, drop = FALSE]
  # A non-finite entry among parameters that *were* estimated is a gap: the
  # covariance exists but cannot be decomposed.
  if (anyNA(V) || any(!is.finite(V))) return(NA)

  nms <- if (length(param_names) == nrow(vcov)) {
    as.character(param_names)[keep]
  } else {
    paste0("par", keep)
  }

  # (3) Standardise to a correlation matrix (see note above on scaling).
  s <- sqrt(diag(V))
  R <- V / outer(s, s)
  if (anyNA(R) || any(!is.finite(R))) return(NA)

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
  if (is.null(e)) return(NA)

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

  if (is.null(found)) return(NULL)
  found$n_directions <- length(seen)
  found
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
