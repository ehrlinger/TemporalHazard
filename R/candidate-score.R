# candidate-score.R -- Unified criterion wrapper for stepwise selection.
#
# Step 8.2 of STEPWISE-DESIGN.md.  Wraps the Wald helper (Phase 4b sec.8.1)
# and an AIC path behind a single signature so the forward and backward
# drivers can share structure.
#
# Score semantics
# ---------------
# Lower `score` is always "better" in the candidate's favour:
#
#   criterion="score", mode="entry" ->  p_value(Q for the pinned new coef)
#   criterion="wald", mode="entry"  ->  p_value(candidate for new coef)
#   criterion="wald", mode="drop"   ->  1 - p_value(current)  (high p -> easy drop)
#   criterion="aic",  mode="entry"  ->  AIC(candidate) - AIC(current)
#   criterion="aic",  mode="drop"   ->  DeltaAIC_drop ~ W - 2*df  (Wald->LR approx)
#
# `criterion = "score"` is an ENTRY-ONLY criterion: it exists to remove the
# per-candidate refit from the forward step, and the drop path never had one to
# remove.  Following SAS, removal is tested on the current model's Wald
# p-value, so hzr_stepwise() maps score -> wald for the backward step.
#
# For forward selection: argmin(score) picks the best candidate and the
# threshold test is `score < slentry` (Wald) or `score < 0` (AIC).
# For backward selection: argmin(score) picks the variable most worth
# dropping and the threshold test is `score < (1 - slstay)` (Wald) or
# `score < 0` (AIC).
#
# The list always carries both the raw `p_value` and `delta_aic` fields
# so the selection trace can show whichever is most informative
# regardless of the active criterion.

#' AIC for a fitted hazard model
#'
#' @param fit A fitted `hazard` object (`fit = TRUE`).
#' @return Numeric scalar AIC, or `NA_real_` if the log-likelihood is
#'   not available.
#'
#' @keywords internal
#' @noRd
.hzr_aic <- function(fit) {
  if (!inherits(fit, "hazard")) {
    stop("`fit` must be a `hazard` object.", call. = FALSE)
  }
  loglik <- fit$fit$objective
  if (is.null(loglik) || !is.finite(loglik)) {
    return(NA_real_)
  }

  # Count free parameters.  `fixed_mask` is TRUE for fixed entries
  # (see likelihood-multiphase.R: `fixed_mask <- !free_mask`), so free
  # parameters are the `!fixed_mask` count.  Fall back to the full
  # theta length when the mask is absent or length-mismatched.
  mask <- fit$fit$fixed_mask
  k <- if (!is.null(mask) && length(mask) == length(fit$fit$theta)) {
    sum(!as.logical(mask))
  } else {
    length(fit$fit$theta)
  }

  -2 * loglik + 2 * k
}


#' Score a stepwise candidate under Wald or AIC criterion
#'
#' Produces a unified "smaller is better" score for either adding new
#' coefficients to a model (`mode = "entry"`) or dropping coefficients
#' from the current model (`mode = "drop"`), under a Wald or AIC
#' criterion.
#'
#' @param criterion One of `"score"`, `"wald"`, or `"aic"`.  `"score"` is
#'   entry-only (see the file header).
#' @param mode Either `"entry"` or `"drop"`.
#' @param current Fitted `hazard` object that is the starting point for
#'   the step.  Always required (AIC baseline; Wald target for drops).
#' @param candidate Fitted `hazard` object containing the proposed new
#'   coefficient(s).  Required when `mode = "entry"` under the `"wald"` and
#'   `"aic"` criteria; ignored otherwise.  The `"score"` criterion has no
#'   candidate fit -- that is the point of it -- and takes `score` instead.
#' @param score Precomputed `.hzr_score_q()` result (`list(stat, df,
#'   p_value)`).  Required when `criterion = "score"`; ignored otherwise.
#' @param names Character vector of coefficient names being tested.
#'   Must exist in the target fit (`candidate` for entries, `current`
#'   for drops).  Purely descriptive under `criterion = "score"`, which
#'   never looks a coefficient up by name.
#'
#' @return A list with:
#' \describe{
#'   \item{score}{Numeric, smaller is better (see file header).}
#'   \item{criterion}{Echoed input.}
#'   \item{mode}{Echoed input.}
#'   \item{stat}{Wald z (scalar) or chi^2 (joint).}
#'   \item{df}{Integer degrees of freedom.}
#'   \item{p_value}{Always populated when computable.}
#'   \item{delta_aic}{Always populated when computable.}
#'   \item{names}{Echoed input.}
#' }
#'
#' @keywords internal
#' @noRd
.hzr_candidate_score <- function(criterion = c("score", "wald", "aic"),
                                  mode      = c("entry", "drop"),
                                  current,
                                  candidate = NULL,
                                  score     = NULL,
                                  names) {
  criterion <- match.arg(criterion)
  mode      <- match.arg(mode)

  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }

  if (criterion == "score") {
    # Entry-only.  hzr_stepwise() maps score -> wald before the backward step,
    # so reaching here with mode = "drop" is a caller bug, not a user error.
    if (mode != "entry") {
      stop("The score criterion is entry-only; `mode` must be 'entry'.",
           call. = FALSE)
    }
    # `names` is a formal of this function, so spell base::names() out rather
    # than relying on R skipping the non-function binding.
    if (!is.list(score) ||
          !all(c("stat", "df", "p_value") %in% base::names(score))) {
      stop("`score` must be a `.hzr_score_q()` result when criterion = 'score'.",
           call. = FALSE)
    }
    return(list(
      score     = score$p_value,   # smaller p = stronger entry case
      criterion = criterion,
      mode      = mode,
      stat      = score$stat,
      df        = score$df,
      p_value   = score$p_value,
      delta_aic = NA_real_,        # no candidate refit to take an AIC from
      names     = names
    ))
  }

  if (mode == "entry") {
    if (is.null(candidate) || !inherits(candidate, "hazard")) {
      stop("`candidate` must be a fitted `hazard` object when mode = 'entry'.",
           call. = FALSE)
    }
    test_fit <- candidate
  } else {
    test_fit <- current
  }

  wald <- .hzr_wald_p(test_fit, names)

  # Wald -> LR approximation:  W = z^2 when df = 1;
  # otherwise W is already the chi^2 statistic.
  W <- if (wald$df == 1L) {
    if (is.finite(wald$stat)) wald$stat^2 else NA_real_
  } else {
    if (is.finite(wald$stat)) wald$stat else NA_real_
  }

  # AIC components
  if (criterion == "aic" && mode == "entry") {
    aic_cur <- .hzr_aic(current)
    aic_can <- .hzr_aic(candidate)
    delta   <- if (is.finite(aic_cur) && is.finite(aic_can)) {
      aic_can - aic_cur
    } else {
      NA_real_
    }
  } else if (criterion == "aic" && mode == "drop") {
    delta <- if (is.finite(W)) W - 2 * wald$df else NA_real_
  } else if (mode == "drop") {
    # Wald + drop: compute DeltaAIC approx anyway so the trace can show it.
    delta <- if (is.finite(W)) W - 2 * wald$df else NA_real_
  } else {
    # Wald + entry: no DeltaAIC available without refitting;  leave as NA.
    delta <- NA_real_
  }

  score <- if (criterion == "wald") {
    if (mode == "entry") {
      wald$p_value           # smaller p = stronger entry case
    } else {
      if (is.na(wald$p_value)) NA_real_ else 1 - wald$p_value
    }                        # larger p = easier drop, so 1-p is small
  } else {
    delta                    # AIC: smaller DeltaAIC = better outcome
  }

  list(
    score     = score,
    criterion = criterion,
    mode      = mode,
    stat      = wald$stat,
    df        = wald$df,
    p_value   = wald$p_value,
    delta_aic = delta,
    names     = names
  )
}
