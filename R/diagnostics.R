# diagnostics.R -- Model validation and goodness-of-fit utilities
#
# Implements SAS HAZARD macro equivalents for calibration and validation.

#' @importFrom stats predict
NULL

#' Decile-of-risk calibration
#'
#' Partition observations into groups (default 10) by predicted risk and
#' compare observed vs. expected event counts in each group.  Good
#' calibration means the two track each other across the risk spectrum.
#'
#' This reproduces the SAS `deciles.hazard.sas` macro. **All** subjects are
#' ranked by predicted survival at the horizon `time` and split into equal-sized
#' risk groups. Within each group the **expected** event count is the sum of
#' each subject's predicted cumulative hazard at its *own* follow-up time, and
#' the **observed** count is its number of events; under conservation of events
#' the group totals sum to the total observed events. The horizon therefore only
#' stratifies subjects into risk groups -- it does not restrict or exclude any
#' subject, and the expected/observed totals are independent of it.
#'
#' @param object A fitted `hazard` object (with `fit = TRUE`).
#' @param time Numeric scalar: the horizon at which predicted survival is used
#'   to **rank subjects into risk groups** (e.g. `time = 12` ranks by 12-month
#'   predicted survival). It does not restrict the event/expected counts, which
#'   are accumulated over each subject's full follow-up.
#' @param groups Integer: number of risk groups (default 10 for deciles).
#' @param status Optional numeric vector of event indicators (1 = event,
#'   0 = censored).
#'   If `NULL` (default), extracted from the fitted object's stored data.
#' @param event_time Optional numeric vector of observed event/censoring
#'   times.  If `NULL` (default), extracted from the fitted object.
#'
#' @return A data frame with one row per risk group and columns:
#' \describe{
#'   \item{group}{Integer group label (1 = lowest risk, ranked by predicted
#'     survival at \code{time}).}
#'   \item{n}{Number of observations in the group.}
#'   \item{events}{Observed event count in the group (all events over
#'     follow-up).}
#'   \item{expected}{Expected event count: the sum of each subject's predicted
#'     cumulative hazard at its own follow-up time.}
#'   \item{observed_rate}{Observed event rate (events / n).}
#'   \item{expected_rate}{Expected event rate (expected / n).}
#'   \item{chi_sq}{Chi-square contribution: (events - expected)^2 /
#'     expected.}
#'   \item{p_value}{Upper-tail p-value from the chi-square test for
#'     this group (1 df).}
#'   \item{mean_survival}{Mean predicted survival probability at the horizon
#'     in the group.}
#'   \item{mean_cumhaz}{Mean predicted cumulative hazard at follow-up in the
#'     group.}
#' }
#'
#' An attribute `"overall"` is attached with the overall chi-square
#' statistic, degrees of freedom, and p-value.
#'
#' @examples
#' \donttest{
#' data(avc)
#' avc <- na.omit(avc)
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ age + mal,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
#'   fit   = TRUE
#' )
#' cal <- hzr_deciles(fit, time = 120)
#' print(cal)
#' }
#'
#' @seealso [predict.hazard()] for the prediction types used internally.
#' @export
hzr_deciles <- function(object, time, groups = 10L,
                        status = NULL, event_time = NULL) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a fitted hazard object.", call. = FALSE)
  }

  if (is.null(object$fit$theta) ||
      (is.logical(object$fit$converged) && is.na(object$fit$converged))) {
    stop("'object' has no fitted parameters. Refit with fit = TRUE.",
         call. = FALSE)
  }
  if (!is.numeric(time) || length(time) != 1L || time <= 0) {
    stop("'time' must be a single positive number.", call. = FALSE)
  }
  groups <- as.integer(groups)
  if (groups < 2L) {
    stop("'groups' must be at least 2.", call. = FALSE)
  }


  # --- Extract observed data ------------------------------------------------
  if (is.null(status)) {
    status <- object$data$status
  }
  if (is.null(event_time)) {
    event_time <- object$data$time
  }
  n_obs <- length(status)
  if (length(event_time) != n_obs) {
    stop("'status' and 'event_time' must have the same length.",
         call. = FALSE)
  }
  if (any(!is.finite(event_time)) || any(event_time < 0)) {
    stop("'event_time' must be finite and non-negative.", call. = FALSE)
  }
  if (any(!is.finite(status)) || any(!status %in% c(0, 1))) {
    stop("'status' must be coded as 0/1 with finite values.", call. = FALSE)
  }

  n_model <- length(object$data$time)
  if (n_obs != n_model) {
    stop("'status'/'event_time' lengths must match the fitted data length (",
         n_model, ").", call. = FALSE)
  }

  # --- Per-observation predicted cumulative hazard --------------------------
  # SAS %DECILES method (Blackstone/Naftel HAZARD): the "expected events" in a
  # group are the sum of predicted cumulative hazard at each subject's OWN
  # follow-up time, so the total expected equals the observed event count under
  # conservation of events.  Risk grouping is by predicted survival at the
  # calibration horizon `time`; the horizon only stratifies subjects into risk
  # groups -- it does not restrict the event or expected counts, and all
  # subjects are included.
  cumhaz_at <- function(times) {
    if (identical(object$spec$dist, "multiphase")) {
      phases <- object$fit$phases
      if (is.null(phases)) phases <- object$spec$phases
      .hzr_multiphase_cumhaz(times, object$fit$theta, phases,
                             object$fit$covariate_counts, object$fit$x_list,
                             per_phase = FALSE)
    } else if (!is.null(object$data$x) && ncol(object$data$x) > 0) {
      nd <- as.data.frame(object$data$x)
      nd$time <- times
      predict(object, newdata = nd, type = "cumulative_hazard")
    } else {
      predict(object, newdata = data.frame(time = times),
              type = "cumulative_hazard")
    }
  }

  cumhaz_fu    <- cumhaz_at(event_time)        # expected-event contribution
  cumhaz_hor   <- cumhaz_at(rep(time, n_obs))  # risk grouping at the horizon
  survival_hor <- exp(-cumhaz_hor)
  observed     <- as.integer(status == 1)
  n_included   <- n_obs
  n_excluded   <- 0L

  if (groups > n_obs) {
    stop("'groups' must be <= number of observations (", n_obs, ").",
         call. = FALSE)
  }

  # --- Rank into equal-sized groups by predicted risk at the horizon --------
  # Rank by cumulative hazard at the horizon, ascending, so group 1 = lowest
  # risk (highest predicted survival).  This is the same partition SAS PROC
  # RANK GROUPS= produces on predicted survival, with the opposite label order
  # (SAS _DECILE_ 0 = highest risk).  Equal-count bins via PROC RANK's floor
  # rule.  Ties are broken by order of appearance ("first") so that an
  # intercept-only model (all predictions identical) still yields equal-sized
  # groups rather than collapsing into one; for models with distinct
  # predictions the tie rule is irrelevant.
  ranks <- rank(cumhaz_hor, ties.method = "first")
  group <- as.integer(floor(groups * (ranks - 1) / n_obs)) + 1L
  group[group > groups] <- groups

  # --- Aggregate by group ---------------------------------------------------
  result <- data.frame(
    group = seq_len(groups),
    n = integer(groups),
    events = numeric(groups),
    expected = numeric(groups),
    observed_rate = numeric(groups),
    expected_rate = numeric(groups),
    chi_sq = numeric(groups),
    p_value = numeric(groups),
    mean_survival = numeric(groups),
    mean_cumhaz = numeric(groups)
  )

  for (g in seq_len(groups)) {
    idx <- which(group == g)
    ng <- length(idx)
    obs_events <- sum(observed[idx] == 1)
    exp_events <- sum(cumhaz_fu[idx])

    result$n[g] <- ng
    result$events[g] <- obs_events
    result$expected[g] <- exp_events
    result$observed_rate[g] <- if (ng > 0) obs_events / ng else NA_real_
    result$expected_rate[g] <- if (ng > 0) exp_events / ng else NA_real_
    result$mean_survival[g] <- if (ng > 0) mean(survival_hor[idx]) else NA_real_
    result$mean_cumhaz[g] <- if (ng > 0) mean(cumhaz_fu[idx]) else NA_real_

    # Per-group chi-square: (O - E)^2 / E
    if (exp_events > 0) {
      result$chi_sq[g] <- (obs_events - exp_events)^2 / exp_events
      # Upper-tail p-value from chi-square with 1 df
      result$p_value[g] <- stats::pchisq(result$chi_sq[g], df = 1,
                                          lower.tail = FALSE)
    } else {
      result$chi_sq[g] <- NA_real_
      result$p_value[g] <- NA_real_
    }
  }

  # --- Overall chi-square ---------------------------------------------------
  # Only groups with a positive expected count contribute (a zero-expected
  # group has an undefined (O-E)^2/E term and NA chi_sq). Clamp df at 0 so
  # degenerate inputs (<=1 contributing group) never yield a negative df.
  valid <- result$expected > 0 & !is.na(result$chi_sq)
  overall_chi_sq <- sum(result$chi_sq[valid])
  overall_df <- max(0L, sum(valid) - 1L)
  overall_p <- if (overall_df > 0) {
    stats::pchisq(overall_chi_sq, df = overall_df, lower.tail = FALSE)
  } else {
    NA_real_
  }

  attr(result, "overall") <- list(
    chi_sq = overall_chi_sq,
    df = overall_df,
    p_value = overall_p,
    time = time,
    groups = groups,
    total_events = sum(observed == 1),
    total_expected = sum(cumhaz_fu),
    n_included = n_included,
    n_excluded = n_excluded
  )

  class(result) <- c("hzr_deciles", "data.frame")
  result
}

#' Print method for hzr_deciles
#'
#' @param x An `hzr_deciles` object.
#' @param digits Number of decimal places for formatting.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x} of class \code{c("hzr_deciles", "data.frame")},
#'   invisibly. The data frame has one row per risk group and columns:
#'   \code{group} (integer group index, 1 = lowest risk),
#'   \code{n} (group size),
#'   \code{events} (observed event count),
#'   \code{expected} (expected event count from model predictions),
#'   \code{observed_rate}, \code{expected_rate} (events / n),
#'   \code{chi_sq} (per-group (O-E)^2/E contribution),
#'   \code{p_value} (1-df chi-square upper-tail p),
#'   \code{mean_survival}, \code{mean_cumhaz} (mean predicted values in group).
#'   An \code{"overall"} attribute contains the omnibus chi-square test
#'   (fields: \code{chi_sq}, \code{df}, \code{p_value}, \code{time},
#'   \code{groups}, \code{total_events}, \code{total_expected},
#'   \code{n_included}, \code{n_excluded}).
#' @export
print.hzr_deciles <- function(x, digits = 3, ...) {
  ov <- attr(x, "overall")
  cat("Decile-of-risk calibration (risk grouped at time =", ov$time, ")\n")
  if (!is.null(ov$n_included)) {
    cat(ov$n_included, "subjects, all included.\n")
  }
  cat(ov$groups, "groups,", ov$total_events, "observed events,",
      signif(ov$total_expected, digits), "expected\n\n")

  # Format for display
  display <- x
  display$expected <- signif(display$expected, digits)
  display$observed_rate <- signif(display$observed_rate, digits)
  display$expected_rate <- signif(display$expected_rate, digits)
  display$chi_sq <- signif(display$chi_sq, digits)
  display$p_value <- signif(display$p_value, digits)
  display$mean_survival <- signif(display$mean_survival, digits)
  display$mean_cumhaz <- signif(display$mean_cumhaz, digits)
  class(display) <- "data.frame"
  print(display, row.names = FALSE)

  cat("\nOverall: chi-sq =", signif(ov$chi_sq, digits),
      "on", ov$df, "df, p =",
      format.pval(ov$p_value, digits = digits), "\n")

  invisible(x)
}


# =========================================================================
# hzr_gof -- Observed vs. expected goodness-of-fit
# =========================================================================

#' Goodness-of-fit: observed vs. predicted events
#'
#' Compare a fitted hazard model against the nonparametric Kaplan-Meier
#' estimate by computing observed and expected (parametric) event counts
#' at each distinct event time.  This is the R equivalent of the SAS
#' `hazplot.sas` macro and implements the conservation-of-events
#' diagnostic.
#'
#' At each observed event time the function computes:
#' \itemize{
#'   \item The Kaplan-Meier survival and cumulative hazard.
#'   \item The parametric survival and cumulative hazard from the fitted
#'     model (and per-phase components for multiphase models).
#'   \item Cumulative observed events vs. cumulative expected events
#'     (sum of individual cumulative hazards for those exiting the risk
#'     set at each time).
#'   \item The running residual (expected minus observed).
#' }
#'
#' Perfect model fit implies the expected and observed event counts track
#' each other (residual near zero).  This is the conservation-of-events
#' principle.
#'
#' @param object A fitted `hazard` object (with `fit = TRUE`).
#' @param time_grid Optional numeric vector of time points at which to
#'   evaluate the parametric model.
#'   If `NULL` (default), uses the sorted unique event times from the
#'   fitted data.
#'
#' @return A data frame with one row per time point and columns:
#' \describe{
#'   \item{time}{Evaluation time.}
#'   \item{n_risk}{Number at risk (Kaplan-Meier).}
#'   \item{n_event}{Number of events at this time.}
#'   \item{n_censor}{Number censored at this time.}
#'   \item{km_surv}{Kaplan-Meier survival estimate.}
#'   \item{km_cumhaz}{Kaplan-Meier cumulative hazard
#'     (\eqn{-\log(\text{km\_surv})}).}
#'   \item{par_surv}{Parametric survival from the fitted model.}
#'   \item{par_cumhaz}{Parametric cumulative hazard.}
#'   \item{cum_observed}{Cumulative observed events to this time.}
#'   \item{cum_expected}{Cumulative expected events (sum of individual
#'     cumulative hazards for observations exiting the risk set).}
#'   \item{residual}{Expected minus observed
#'     (\code{cum_expected - cum_observed}).}
#' }
#'
#' For multiphase models, additional columns are appended for each
#' phase: \code{par_cumhaz_<phase>}.
#'
#' An attribute `"summary"` is attached with scalar diagnostics:
#' total observed events, total expected events, and the final residual.
#'
#' @examples
#' \donttest{
#' data(avc)
#' avc <- na.omit(avc)
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ age + mal,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
#'   fit   = TRUE
#' )
#' gof <- hzr_gof(fit)
#' print(gof)
#'
#' # Plot observed vs expected events
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(gof, aes(x = time)) +
#'     geom_line(aes(y = cum_observed), colour = "#D55E00") +
#'     geom_line(aes(y = cum_expected), colour = "#0072B2") +
#'     labs(x = "Time", y = "Cumulative events") +
#'     theme_minimal()
#' }
#' }
#'
#' @seealso [hzr_deciles()] for decile-of-risk calibration,
#'   [predict.hazard()] for prediction types.
#' @export
hzr_gof <- function(object, time_grid = NULL) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a fitted hazard object.", call. = FALSE)
  }
  if (is.null(object$fit$theta) ||
      (is.logical(object$fit$converged) && is.na(object$fit$converged))) {
    stop("'object' has no fitted parameters. Refit with fit = TRUE.",
         call. = FALSE)
  }

  # --- Extract observed data ------------------------------------------------
  obs_time   <- object$data$time
  obs_status <- object$data$status
  n_total    <- length(obs_time)

  if (length(obs_status) != n_total) {
    stop("Stored 'time' and 'status' vectors must have the same length.",
         call. = FALSE)
  }
  if (any(!is.finite(obs_time)) || any(obs_time < 0)) {
    stop("Stored event/censoring times must be finite and non-negative.",
         call. = FALSE)
  }
  if (any(!is.finite(obs_status)) || any(!obs_status %in% c(0, 1))) {
    stop("Stored status must be coded as 0/1 with finite values.",
         call. = FALSE)
  }

  # --- Kaplan-Meier via survival::survfit -----------------------------------
  km_fit <- survival::survfit(survival::Surv(obs_time, obs_status) ~ 1)

  # survfit output: time, n.risk, n.event, n.censor, surv
  km_times   <- km_fit$time
  km_n_risk  <- km_fit$n.risk
  km_n_event <- km_fit$n.event
  km_n_censor <- km_fit$n.censor
  km_surv    <- km_fit$surv

  # --- Decide time grid -----------------------------------------------------
  if (is.null(time_grid)) {
    time_grid <- km_times
  }

  # --- Parametric predictions at each time point ----------------------------
  is_multiphase <- (object$spec$dist == "multiphase")

  if (!is.null(object$data$x) && ncol(object$data$x) > 0) {
    # For covariate models, evaluate at covariate means (baseline patient)
    x_means <- colMeans(object$data$x)
    nd <- as.data.frame(t(x_means))
    nd <- nd[rep(1, length(time_grid)), , drop = FALSE]
    nd$time <- time_grid
  } else {
    nd <- data.frame(time = time_grid)
  }

  par_cumhaz <- predict(object, newdata = nd, type = "cumulative_hazard")

  # Phase decomposition for multiphase models
  phase_cumhaz <- NULL
  if (is_multiphase) {
    decomp <- predict(object, newdata = nd, type = "cumulative_hazard",
                      decompose = TRUE)
    # decomp is a matrix; first column is "total", rest are phase names
    phase_cols <- colnames(decomp)[colnames(decomp) != "total"]
    phase_cumhaz <- as.data.frame(decomp[, phase_cols, drop = FALSE])
  }

  par_surv <- exp(-par_cumhaz)

  # --- Interpolate KM at the time grid --------------------------------------
  # Use stepfun-style interpolation for KM (right-continuous)
  km_surv_at_grid <- stats::approx(
    x = c(0, km_times), y = c(1, km_surv),
    xout = time_grid, method = "constant", f = 0, rule = 2
  )$y
  km_cumhaz_at_grid <- -log(pmax(km_surv_at_grid, .Machine$double.xmin))

  # Interpolate n.risk, n.event, n.censor at grid times
  # For event counts, sum events at matching times; 0 otherwise
  km_n_risk_grid <- stats::approx(
    x = c(0, km_times), y = c(n_total, km_n_risk),
    xout = time_grid, method = "constant", f = 0, rule = 2
  )$y
  km_n_event_grid <- rep(0, length(time_grid))
  km_n_censor_grid <- rep(0, length(time_grid))
  for (i in seq_along(km_times)) {
    match_idx <- which(abs(time_grid - km_times[i]) < .Machine$double.eps * 100)
    if (length(match_idx) > 0) {
      km_n_event_grid[match_idx[1]] <- km_n_event[i]
      km_n_censor_grid[match_idx[1]] <- km_n_censor[i]
    }
  }

  # --- Conservation of Events accounting ------------------------------------
  # At each event time, accumulate:
  #   cum_observed: running sum of observed events
  #   cum_expected: running sum of individual cumulative hazards for
  #                 observations exiting the risk set (events + censored)
  #
  # The expected events for observations leaving at time t is:
  #   (n_event + n_censor) * parametric_cumhaz(t)
  # This is the SAS hazplot approach: total * _CUMHAZ at that interval.

  cum_observed <- cumsum(km_n_event_grid)
  interval_expected <- (km_n_event_grid + km_n_censor_grid) * par_cumhaz
  cum_expected <- cumsum(interval_expected)
  residual <- cum_expected - cum_observed

  # --- Assemble result ------------------------------------------------------
  result <- data.frame(
    time         = time_grid,
    n_risk       = km_n_risk_grid,
    n_event      = km_n_event_grid,
    n_censor     = km_n_censor_grid,
    km_surv      = km_surv_at_grid,
    km_cumhaz    = km_cumhaz_at_grid,
    par_surv     = par_surv,
    par_cumhaz   = par_cumhaz,
    cum_observed = cum_observed,
    cum_expected = cum_expected,
    residual     = residual
  )

  # Add phase columns for multiphase
  if (!is.null(phase_cumhaz)) {
    for (ph in names(phase_cumhaz)) {
      result[[paste0("par_cumhaz_", ph)]] <- phase_cumhaz[[ph]]
    }
  }

  # Summary diagnostics
  attr(result, "summary") <- list(
    total_observed = cum_observed[length(cum_observed)],
    total_expected = cum_expected[length(cum_expected)],
    final_residual = residual[length(residual)],
    dist = object$spec$dist,
    n = n_total
  )

  class(result) <- c("hzr_gof", "data.frame")
  result
}

#' Print method for hzr_gof
#'
#' @param x An `hzr_gof` object.
#' @param digits Number of decimal places for formatting.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x} of class \code{c("hzr_gof", "data.frame")},
#'   invisibly. The data frame has one row per time point and columns:
#'   \code{time}, \code{n_risk}, \code{n_event}, \code{n_censor},
#'   \code{km_surv} (Kaplan-Meier survival), \code{km_cumhaz},
#'   \code{par_surv} (parametric survival), \code{par_cumhaz},
#'   \code{cum_observed} (cumulative observed events),
#'   \code{cum_expected} (cumulative expected events from model),
#'   \code{residual} (cum_expected - cum_observed).
#'   Multiphase models additionally include \code{par_cumhaz_<phase>} columns
#'   for per-phase cumulative hazard contributions.
#'   A \code{"summary"} attribute contains scalar diagnostics:
#'   \code{total_observed}, \code{total_expected}, \code{final_residual},
#'   \code{dist}, \code{n}.
#' @export
print.hzr_gof <- function(x, digits = 3, ...) {
  s <- attr(x, "summary")
  cat("Goodness-of-fit: observed vs. expected events\n")
  cat("Distribution:", s$dist, " | n =", s$n, "\n\n")
  cat("Total observed events:", s$total_observed, "\n")
  cat("Total expected events:", round(s$total_expected, digits), "\n")
  cat("Final residual (E - O):", round(s$final_residual, digits), "\n")
  cat("Conservation ratio (E/O):",
      round(s$total_expected / max(s$total_observed, 1), digits), "\n")
  cat("\nUse plot columns: time, km_surv, par_surv, cum_observed,",
      "cum_expected, residual\n")
  invisible(x)
}


# =========================================================================
# hzr_kaplan -- Kaplan-Meier with exact logit confidence limits
# =========================================================================

#' Kaplan-Meier survival with exact logit confidence limits
#'
#' Compute the product-limit (Kaplan-Meier) survival estimate with
#' logit-transformed confidence limits that respect the \eqn{[0, 1]}
#' boundary.
#' This is the R equivalent of the SAS `kaplan.sas` macro.
#'
#' The standard Greenwood confidence interval can exceed \eqn{[0, 1]} in the
#' tails. The logit-transformed interval avoids this by working on the
#' log-odds scale:
#'
#' \deqn{
#'   \text{CL}_{\text{lower}} = S / \bigl(S + (1-S)\,
#'   \exp(z_\alpha\,\text{SI})\bigr)
#' }
#' \deqn{
#'   \text{CL}_{\text{upper}} = S / \bigl(S + (1-S)\,
#'   \exp(-z_\alpha\,\text{SI})\bigr)
#' }
#'
#' where \eqn{\text{SI} = \sqrt{V_P - 1} / (1 - S)}, \eqn{V_P} is the
#' cumulative Greenwood variance product, and \eqn{z_\alpha} is the
#' normal quantile for the requested confidence level.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator (1 = event, 0 = censored).
#' @param conf_level Confidence level for the interval (default 0.95).
#'   The SAS default of 0.68268948 corresponds to a 1-SD interval.
#' @param event_only Logical; if `TRUE` (default), only return rows at
#'   event times (where `n_event > 0`).
#'   If `FALSE`, return rows at all times reported by
#'   `survival::survfit()` (events and censoring times).
#'
#' @return A data frame with one row per time point and columns:
#' \describe{
#'   \item{time}{Event/censoring time.}
#'   \item{n_risk}{Number at risk at start of interval.}
#'   \item{n_event}{Number of events at this time.}
#'   \item{n_censor}{Number censored at this time.}
#'   \item{survival}{Kaplan-Meier survival estimate.}
#'   \item{std_err}{Standard error of survival (Greenwood).}
#'   \item{cl_lower}{Lower confidence limit (logit-transformed).}
#'   \item{cl_upper}{Upper confidence limit (logit-transformed).}
#'   \item{cumhaz}{Cumulative hazard \eqn{= -\log(S)}.}
#'   \item{hazard}{Interval hazard rate
#'     \eqn{= \log(S_{t-1} / S_t) / \Delta t}.}
#'   \item{density}{Probability density estimate
#'     \eqn{= (S_{t-1} - S_t) / \Delta t}.}
#'   \item{life}{Restricted mean survival time (area under curve to
#'     this time).}
#' }
#'
#' @examples
#' data(cabgkul)
#' km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)
#' head(km)
#'
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(km, aes(time)) +
#'     geom_step(aes(y = survival * 100)) +
#'     geom_ribbon(aes(ymin = cl_lower * 100, ymax = cl_upper * 100),
#'                 stat = "identity", alpha = 0.2) +
#'     labs(x = "Months", y = "Survival (%)") +
#'     theme_minimal()
#' }
#' }
#'
#' @references
#' Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete
#' observations. *J Am Stat Assoc* 53(282):457--481.
#' \doi{10.1080/01621459.1958.10501452}
#'
#' Greenwood M (1926). The natural duration of cancer. *Reports on Public
#' Health and Medical Subjects* 33:1--26.
#' @seealso [hzr_gof()] for parametric vs. nonparametric comparison.
#' @export
hzr_kaplan <- function(time, status, conf_level = 0.95,
                       event_only = TRUE) {
  if (!is.numeric(time) || !is.numeric(status)) {
    stop("'time' and 'status' must be numeric vectors.", call. = FALSE)
  }
  if (length(time) != length(status)) {
    stop("'time' and 'status' must have the same length.", call. = FALSE)
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.", call. = FALSE)
  }

  z_alpha <- stats::qnorm(0.5 + 0.5 * conf_level)

  # Use survival::survfit for the core KM computation
  km_fit <- survival::survfit(survival::Surv(time, status) ~ 1)

  n_times <- length(km_fit$time)

  # --- Build output vectors at each event/censoring time --------------------
  km_time    <- km_fit$time
  km_n_risk  <- km_fit$n.risk
  km_n_event <- km_fit$n.event
  km_n_censor <- km_fit$n.censor
  km_surv    <- km_fit$surv

  # --- Greenwood variance product and exact logit CL -----------------------
  # VAR_PROD = prod_i [1/(n_i - d_i) - 1/n_i + 1] for all event times
  # The product is cumulative: each time contributes a multiplicative factor.
  var_prod <- rep(1.0, n_times)
  for (i in seq_len(n_times)) {
    ni <- km_n_risk[i]
    di <- km_n_event[i]
    if (ni > 0 && (ni - di) > 0 && di > 0) {
      factor_i <- 1.0 / (ni - di) - 1.0 / ni + 1.0
      if (i == 1L) {
        var_prod[i] <- factor_i
      } else {
        var_prod[i] <- var_prod[i - 1L] * factor_i
      }
    } else if (i > 1L) {
      var_prod[i] <- var_prod[i - 1L]
    }
  }

  # Standard error (Greenwood)
  var_exact <- km_surv^2 * pmax(var_prod - 1, 0)
  std_err <- sqrt(var_exact)

  # Logit-transformed confidence limits (SAS SI_EXACT formula)
  si_exact <- rep(0.0, n_times)
  idx_valid <- km_surv < 1 & km_surv > 0
  si_exact[idx_valid] <- sqrt(pmax(var_prod[idx_valid] - 1, 0)) /
    (1 - km_surv[idx_valid])

  cl_lower <- km_surv / (km_surv + (1 - km_surv) *
                            exp(z_alpha * si_exact))
  cl_upper <- km_surv / (km_surv + (1 - km_surv) *
                            exp(-z_alpha * si_exact))

  # Edge cases: if S = 1 or S = 0, CL = S

  cl_lower[km_surv >= 1] <- 1
  cl_upper[km_surv >= 1] <- 1
  cl_lower[km_surv <= 0] <- 0
  cl_upper[km_surv <= 0] <- 0

  # --- Cumulative hazard -----------------------------------------------------
  cumhaz <- -log(pmax(km_surv, .Machine$double.xmin))

  # --- Interval hazard and density ------------------------------------------
  lag_surv <- c(1, km_surv[-n_times])
  lag_time <- c(0, km_time[-n_times])
  delta_t <- km_time - lag_time

  hazard <- rep(NA_real_, n_times)
  density <- rep(NA_real_, n_times)
  idx_haz <- km_n_event > 0 & delta_t > 0 & km_surv > 0
  hazard[idx_haz] <- log(lag_surv[idx_haz] / km_surv[idx_haz]) /
    delta_t[idx_haz]
  density[idx_haz] <- (lag_surv[idx_haz] - km_surv[idx_haz]) /
    delta_t[idx_haz]

  # --- Restricted mean survival (life integral) -----------------------------
  # KM survival is a right-continuous step function.  RMST is accumulated
  # as the sum of rectangle areas: each interval contributes dt * S(t-1).
  life <- rep(0, n_times)
  lag_life <- 0
  for (i in seq_len(n_times)) {
    life[i] <- lag_life + delta_t[i] * lag_surv[i]
    lag_life <- life[i]
  }

  # --- Assemble result ------------------------------------------------------
  result <- data.frame(
    time      = km_time,
    n_risk    = km_n_risk,
    n_event   = km_n_event,
    n_censor  = km_n_censor,
    survival  = km_surv,
    std_err   = std_err,
    cl_lower  = cl_lower,
    cl_upper  = cl_upper,
    cumhaz    = cumhaz,
    hazard    = hazard,
    density   = density,
    life      = life
  )

  if (event_only) {
    result <- result[result$n_event > 0, , drop = FALSE]
    rownames(result) <- NULL
  }

  class(result) <- c("hzr_kaplan", "data.frame")
  result
}

#' Print method for hzr_kaplan
#'
#' @param x An `hzr_kaplan` object.
#' @param digits Number of decimal places for formatting.
#' @param n Maximum rows to print (default 20).
#' @param ... Additional arguments (ignored).
#' @return The object \code{x} of class \code{c("hzr_kaplan", "data.frame")},
#'   invisibly. The data frame has one row per event time (or all times when
#'   \code{event_only = FALSE}) and columns:
#'   \code{time} (follow-up time),
#'   \code{n_risk} (number at risk),
#'   \code{n_event} (events in interval),
#'   \code{n_censor} (censored observations in interval),
#'   \code{survival} (Kaplan-Meier survival estimate),
#'   \code{std_err} (Greenwood standard error on log-hazard scale),
#'   \code{cl_lower}, \code{cl_upper} (logit-transformed confidence limits
#'   on the survival scale),
#'   \code{cumhaz} (Nelson-Aalen cumulative hazard),
#'   \code{hazard} (interval hazard estimate),
#'   \code{density} (estimated event density),
#'   \code{life} (life-table life expectancy contribution).
#' @export
print.hzr_kaplan <- function(x, digits = 4, n = 20, ...) {
  cat("Kaplan-Meier estimate with logit confidence limits\n")
  total_events <- sum(x$n_event)
  cat("Events:", total_events, " | Time points:", nrow(x), "\n")
  if (nrow(x) > 0) {
    cat("Survival range:",
        round(min(x$survival), digits), "to",
        round(max(x$survival), digits), "\n")
    cat("RMST at last event:", round(x$life[nrow(x)], digits), "\n\n")
  }
  display <- x
  display$survival <- round(display$survival, digits)
  display$std_err <- round(display$std_err, digits)
  display$cl_lower <- round(display$cl_lower, digits)
  display$cl_upper <- round(display$cl_upper, digits)
  display$cumhaz <- round(display$cumhaz, digits)
  display$hazard <- round(display$hazard, digits)
  display$density <- round(display$density, digits)
  display$life <- round(display$life, digits)
  class(display) <- "data.frame"
  if (nrow(display) > n) {
    print(utils::head(display, n), row.names = FALSE)
    cat("... (", nrow(display) - n, " more rows)\n")
  } else {
    print(display, row.names = FALSE)
  }
  invisible(x)
}


# =========================================================================
# hzr_calibrate -- Variable calibration (logit / Gompertz / Cox transform)
# =========================================================================

#' Calibrate a continuous variable against an outcome
#'
#' Group a continuous covariate into quantile bins, compute the event
#' probability (or hazard rate) per bin, and apply a link transform
#' (logit, Gompertz, or Cox).
#' This is the R equivalent of the SAS `logit.sas` and `logitgr.sas`
#' macros.
#'
#' Use this function before model entry to assess whether the covariate
#' relationship with the outcome is approximately linear on the link
#' scale. If the transformed probabilities are roughly linear against
#' the group means, the covariate can enter the model untransformed.
#' Curvature suggests a transformation (log, quadratic) may improve fit.
#'
#' @param x Numeric vector: the continuous covariate to calibrate.
#' @param event Numeric vector: event indicator (1 = event, 0 = no event).
#' @param groups Integer: number of quantile bins (default 10).
#' @param by Optional factor or character vector for stratified calibration
#'   (SAS `logitgr.sas` functionality). If provided, calibration is
#'   computed within each stratum. Default `NULL` (no stratification).
#' @param link Character: transform to apply to event probabilities.
#'   One of `"logit"` (default), `"gompertz"` (complementary log-log),
#'   or `"cox"`.
#' @param time Optional numeric vector: follow-up time, required when
#'   `link = "cox"`. The Cox link computes
#'   \eqn{\log(\text{events} / \sum \text{time})} (constant hazard rate).
#'
#' @return A data frame with one row per group (or per group-by-stratum
#'   combination) and columns:
#' \describe{
#'   \item{group}{Integer group label.}
#'   \item{by}{Stratum level (only present when `by` is provided).}
#'   \item{n}{Number of observations in the group.}
#'   \item{events}{Number of events.}
#'   \item{mean}{Mean of `x` within the group.}
#'   \item{min}{Minimum of `x` within the group.}
#'   \item{max}{Maximum of `x` within the group.}
#'   \item{prob}{Event probability (events / n), or for Cox link:
#'     events / sum(time).}
#'   \item{link_value}{Transformed probability on the chosen link scale.}
#' }
#'
#' @examples
#' data(avc)
#' avc <- na.omit(avc)
#'
#' # Logit calibration of age
#' cal <- hzr_calibrate(avc$age, avc$dead, groups = 10)
#' print(cal)
#'
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(cal, aes(mean, link_value)) +
#'     geom_point(size = 3) +
#'     geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
#'                 linetype = "dashed") +
#'     labs(x = "Age at repair (months)", y = "Logit(P(death))") +
#'     theme_minimal()
#' }
#' }
#'
#' @seealso [hzr_deciles()] for model-based calibration after fitting.
#' @export
hzr_calibrate <- function(x, event, groups = 10L, by = NULL,
                          link = c("logit", "gompertz", "cox"),
                          time = NULL) {
  link <- match.arg(link)

  if (!is.numeric(x) || !is.numeric(event)) {
    stop("'x' and 'event' must be numeric vectors.", call. = FALSE)
  }
  if (length(x) != length(event)) {
    stop("'x' and 'event' must have the same length.", call. = FALSE)
  }
  groups <- as.integer(groups)
  if (groups < 2L) {
    stop("'groups' must be at least 2.", call. = FALSE)
  }
  if (link == "cox" && is.null(time)) {
    stop("'time' is required when link = 'cox'.", call. = FALSE)
  }
  if (!is.null(time) && length(time) != length(x)) {
    stop("'time' must have the same length as 'x'.", call. = FALSE)
  }

  # --- Assign quantile groups -----------------------------------------------
  ranks <- rank(x, ties.method = "first")
  n <- length(x)
  group <- as.integer(cut(ranks,
                          breaks = seq(0, n, length.out = groups + 1L),
                          include.lowest = TRUE,
                          labels = seq_len(groups)))

  # --- Build data frame for aggregation ------------------------------------
  df <- data.frame(x = x, event = event, group = group)
  if (!is.null(time)) df$time <- time
  if (!is.null(by)) df$by <- by

  # --- Aggregation function -------------------------------------------------
  .calibrate_group <- function(d) {
    ng <- nrow(d)
    ev <- sum(d$event)
    if (link == "cox") {
      prob <- if (sum(d$time) > 0) ev / sum(d$time) else NA_real_
    } else {
      prob <- if (ng > 0) ev / ng else NA_real_
    }

    # Apply link transform
    link_val <- if (is.na(prob) || prob <= 0 || prob >= 1) {
      NA_real_
    } else if (link == "logit") {
      log(prob / (1 - prob))
    } else if (link == "gompertz") {
      log(-log(1 - prob))
    } else {
      # Cox: prob is already events/time (hazard rate)
      if (prob > 0) log(prob) else NA_real_
    }

    data.frame(
      n      = ng,
      events = ev,
      mean   = mean(d$x),
      min    = min(d$x),
      max    = max(d$x),
      prob   = prob,
      link_value = link_val
    )
  }

  # --- Aggregate by group (and optionally by stratum) -----------------------
  if (is.null(by)) {
    result_list <- lapply(split(df, df$group), .calibrate_group)
    result <- do.call(rbind, result_list)
    result$group <- as.integer(rownames(result))
    rownames(result) <- NULL
    result <- result[order(result$group),
                     c("group", "n", "events", "mean", "min", "max",
                       "prob", "link_value")]
  } else {
    splits <- split(df, list(df$by, df$group))
    result_list <- lapply(splits, function(d) {
      if (nrow(d) == 0) return(NULL)
      r <- .calibrate_group(d)
      r$by <- d$by[1]
      r$group <- d$group[1]
      r
    })
    result <- do.call(rbind, Filter(Negate(is.null), result_list))
    rownames(result) <- NULL
    result <- result[order(result$by, result$group),
                     c("group", "by", "n", "events", "mean", "min", "max",
                       "prob", "link_value")]
  }

  attr(result, "link") <- link
  attr(result, "groups") <- groups
  class(result) <- c("hzr_calibrate", "data.frame")
  result
}

#' Print method for hzr_calibrate
#'
#' @param x An `hzr_calibrate` object.
#' @param digits Number of decimal places for formatting.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x} of class \code{c("hzr_calibrate", "data.frame")},
#'   invisibly. The data frame has one row per quantile group and columns:
#'   \code{group} (group index),
#'   \code{n} (group size),
#'   \code{events} (event count),
#'   \code{mean}, \code{min}, \code{max} (covariate summary within group),
#'   \code{prob} (observed event probability),
#'   \code{link_value} (transformed probability on the link scale).
#'   When stratified via the \code{by} argument, a \code{by} column is also
#'   present. Attributes: \code{"link"} (the transform applied,
#'   e.g. \code{"logit"}) and \code{"groups"} (number of quantile bins).
#' @export
print.hzr_calibrate <- function(x, digits = 3, ...) {
  lnk <- attr(x, "link")
  grp <- attr(x, "groups")
  cat("Variable calibration (", lnk, " link, ", grp, " groups)\n",
      sep = "")
  if ("by" %in% names(x)) {
    cat("Stratified by:", paste(unique(x$by), collapse = ", "), "\n")
  }
  cat("\n")
  display <- x
  display$mean <- round(display$mean, digits)
  display$min <- round(display$min, digits)
  display$max <- round(display$max, digits)
  display$prob <- round(display$prob, digits)
  display$link_value <- round(display$link_value, digits)
  class(display) <- "data.frame"
  print(display, row.names = FALSE)
  invisible(x)
}


# =========================================================================
# hzr_nelson -- Wayne Nelson cumulative hazard estimator
# =========================================================================

#' Wayne Nelson cumulative hazard estimator with lognormal confidence limits
#'
#' Compute the Nelson-Aalen cumulative hazard estimate with lognormal
#' confidence limits. Supports weighted events for severity-adjusted
#' analyses of repeated/recurrent events.
#' This is the R equivalent of the SAS `nelsonl.sas` macro.
#'
#' Unlike `survival::survfit()` which uses the Breslow estimator with
#' Greenwood variance, this function uses the Wayne Nelson estimator
#' with lognormal confidence limits that are always non-negative.
#'
#' @param time Numeric vector of follow-up times.
#' @param event Numeric event indicator (1 = event, 0 = censored).
#' @param weight Optional numeric vector of event weights (default 1).
#'   Weights are applied only to events (censored observations contribute
#'   zero weight). Use for severity-weighted repeated events.
#' @param conf_level Confidence level for the interval (default 0.95).
#'
#' @return A data frame with one row per unique event time and columns:
#' \describe{
#'   \item{time}{Event time.}
#'   \item{n_risk}{Number at risk.}
#'   \item{n_event}{Number of events at this time.}
#'   \item{weight_sum}{Sum of event weights at this time.}
#'   \item{cumhaz}{Nelson cumulative hazard estimate.}
#'   \item{std_err}{Standard error.}
#'   \item{cl_lower}{Lower lognormal confidence limit.}
#'   \item{cl_upper}{Upper lognormal confidence limit.}
#'   \item{hazard}{Interval hazard rate.}
#'   \item{cum_events}{Cumulative (weighted) event count.}
#' }
#'
#' @examples
#' data(cabgkul)
#' nel <- hzr_nelson(cabgkul$int_dead, cabgkul$dead)
#' head(nel)
#'
#' @references
#' Nelson W (1972). Theory and applications of hazard plotting for censored
#' failure data. *Technometrics* 14(4):945--966.
#' \doi{10.1080/00401706.1972.10488991}
#'
#' Aalen O (1978). Nonparametric inference for a family of counting processes.
#' *Ann Statist* 6(4):701--726. \doi{10.1214/aos/1176344247}
#' @seealso [hzr_kaplan()] for survival estimation.
#' @export
hzr_nelson <- function(time, event, weight = NULL, conf_level = 0.95) {
  if (!is.numeric(time) || !is.numeric(event)) {
    stop("'time' and 'event' must be numeric vectors.", call. = FALSE)
  }
  n <- length(time)
  if (length(event) != n) {
    stop("'time' and 'event' must have the same length.", call. = FALSE)
  }
  if (is.null(weight)) weight <- rep(1, n)
  if (length(weight) != n) {
    stop("'weight' must have the same length as 'time'.", call. = FALSE)
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.", call. = FALSE)
  }

  z_alpha <- stats::qnorm(0.5 + 0.5 * conf_level)

  # Effective weight: weight * event (censored get 0)
  e_wght <- weight * event

  # Sort by time

  ord <- order(time)
  time_s <- time[ord]
  event_s <- event[ord]
  e_wght_s <- e_wght[ord]

  # Collapse to unique times
  u_times <- sort(unique(time_s[event_s == 1]))
  if (length(u_times) == 0) {
    out <- data.frame(time = numeric(0), n_risk = numeric(0),
                      n_event = numeric(0), weight_sum = numeric(0),
                      cumhaz = numeric(0), std_err = numeric(0),
                      cl_lower = numeric(0), cl_upper = numeric(0),
                      hazard = numeric(0), cum_events = numeric(0))
    class(out) <- c("hzr_nelson", "data.frame")
    return(out)
  }

  out_n <- length(u_times)
  n_risk <- numeric(out_n)
  n_event_out <- numeric(out_n)
  weight_sum <- numeric(out_n)
  cumhaz <- numeric(out_n)

  # Single pass: compute cumhaz, variance, and lognormal CL together
  std_err <- numeric(out_n)
  cl_lower <- numeric(out_n)
  cl_upper <- numeric(out_n)

  run_i_nrisk <- 0
  run_it <- 0
  cum_dist <- 0

  for (k in seq_along(u_times)) {
    t_k <- u_times[k]
    n_risk[k] <- sum(time_s >= t_k)
    at_time <- time_s == t_k & event_s == 1
    n_event_out[k] <- sum(at_time)
    weight_sum[k] <- sum(e_wght_s[at_time])

    # Nelson estimator: dist = sum(weight) / n_risk
    dist_k <- if (n_risk[k] > 0) weight_sum[k] / n_risk[k] else 0
    cum_dist <- cum_dist + dist_k
    cumhaz[k] <- cum_dist

    # Variance and lognormal CL (running accumulators)
    if (n_risk[k] > 0) run_i_nrisk <- run_i_nrisk + 1 / n_risk[k]
    run_it <- run_it + n_event_out[k]

    if (run_it > 0 && cum_dist > 0) {
      var_cef <- run_i_nrisk * cum_dist / run_it
      std_err[k] <- sqrt(var_cef)

      sigma2 <- log(run_i_nrisk / (run_it * cum_dist) + 1)
      sigma <- sqrt(sigma2)
      mu_ln <- log(cum_dist) - sigma2 / 2

      cl_upper[k] <- exp(mu_ln + z_alpha * sigma)
      cl_lower[k] <- exp(mu_ln - z_alpha * sigma)
    }
  }

  # Interval hazard
  lag_cumhaz <- c(0, cumhaz[-out_n])
  lag_time <- c(0, u_times[-out_n])
  delta_t <- u_times - lag_time
  hazard <- rep(NA_real_, out_n)
  idx <- delta_t > 0
  hazard[idx] <- (cumhaz[idx] - lag_cumhaz[idx]) / delta_t[idx]

  result <- data.frame(
    time       = u_times,
    n_risk     = n_risk,
    n_event    = n_event_out,
    weight_sum = weight_sum,
    cumhaz     = cumhaz,
    std_err    = std_err,
    cl_lower   = cl_lower,
    cl_upper   = cl_upper,
    hazard     = hazard,
    cum_events = cumsum(n_event_out)
  )

  class(result) <- c("hzr_nelson", "data.frame")
  result
}

#' @rdname hzr_nelson
#' @param x An `hzr_nelson` object.
#' @param digits Number of decimal places for formatting.
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_nelson <- function(x, digits = 4, ...) {
  cat("Nelson cumulative hazard estimate with lognormal CL\n")
  cat("Events:", sum(x$n_event), " | Time points:", nrow(x), "\n\n")
  display <- x
  for (col in c("cumhaz", "std_err", "cl_lower", "cl_upper", "hazard")) {
    display[[col]] <- round(display[[col]], digits)
  }
  class(display) <- "data.frame"
  print(utils::head(display, 20), row.names = FALSE)
  if (nrow(display) > 20) cat("... (", nrow(display) - 20, " more rows)\n")
  invisible(x)
}


# =========================================================================
# hzr_bootstrap -- Bootstrap inference for hazard models
# =========================================================================

#' Resolve parameter names for a bootstrap replicate's fitted theta
#'
#' Shape parameters are already named in `theta`; covariate betas often
#' come through with empty names. Covariate coefficients occupy the last
#' `ncol(x)` positions of theta -- fill any blanks within that block from
#' the design matrix column names by relative index, so downstream pivots
#' (e.g. `reshape(wide)`) get a distinct column per covariate, even when
#' some betas are already named and others are not.
#'
#' @param fit_obj A fitted `hazard`-like object (has `$fit$theta` and
#'   `$data$x`).
#' @return Character vector of resolved parameter names, same length as
#'   `fit_obj$fit$theta`.
#' @keywords internal
#' @noRd
.hzr_bootstrap_param_names <- function(fit_obj) {
  theta <- fit_obj$fit$theta
  param_names <- names(theta)
  if (is.null(param_names)) {
    param_names <- character(length(theta))
  }
  if (!is.null(fit_obj$data$x)) {
    x_names <- colnames(fit_obj$data$x)
    p <- ncol(fit_obj$data$x)
    n_theta <- length(param_names)
    if (!is.null(x_names) && p > 0L && n_theta >= p) {
      cov_idx <- seq.int(n_theta - p + 1L, n_theta)
      blank_in_block <- !nzchar(param_names[cov_idx])
      param_names[cov_idx[blank_in_block]] <- x_names[blank_in_block]
    }
  }
  still_blank <- !nzchar(param_names)
  if (any(still_blank)) {
    param_names[still_blank] <- paste0("param_", which(still_blank))
  }
  param_names
}

#' Bootstrap resampling for hazard model coefficients
#'
#' Resample data with replacement, refit the hazard model on each
#' replicate, and accumulate coefficient distributions. Returns a tidy
#' data frame of per-replicate estimates with summary statistics.
#' This is the R equivalent of the SAS `bootstrap.hazard.sas` macro.
#'
#' When `scope` is supplied, each replicate instead runs a fresh
#' [hzr_stepwise()] selection on the resampled data (starting from a
#' fixed-shape refit of `object`) instead of refitting `object`'s exact
#' formula. This is the R equivalent of the SAS `%HAZBOOT` macro: fit the
#' hazard shape with no covariates (fixing it via `hzr_phase(..., fixed =
#' "shapes")`), then bootstrap-screen candidate covariates for how often
#' they enter the model. `summary$pct` then reports the selection
#' frequency across replicates, and `summary$mean`/`sd`/`ci_*` describe the
#' coefficient distribution conditional on selection.
#'
#' @param object A fitted `hazard` object (with `fit = TRUE`).
#' @param n_boot Integer: number of bootstrap replicates (default 200).
#' @param fraction Numeric in (0, 1]: fraction of data to sample per
#'   replicate (default 1.0 for full bootstrap; < 1 for bagging).
#' @param seed Optional integer random seed for reproducibility. When
#'   supplied, `set.seed(seed)` is called at function entry, jumping the
#'   global RNG to the seeded state; it is not restored on exit. Pass
#'   `NULL` (the default) to skip the `set.seed()` call and start from
#'   the caller's current RNG state. Note that the bootstrap consumes
#'   random numbers either way, so the global RNG state will advance
#'   during the call -- `seed = NULL` avoids the *reset* at entry, not
#'   the advance during resampling.
#' @param verbose Logical; if `TRUE`, display a text progress bar over the
#'   `n_boot` replicates (via [utils::txtProgressBar()]).
#' @param scope **Experimental.** Candidate variable scope for embedded
#'   stepwise selection during each bootstrap replicate. The argument and the
#'   shape of the object it returns may change in a future release; see the
#'   "Selection mode is experimental" section below. `NULL` (default)
#'   preserves the
#'   original fixed-formula bootstrap: every replicate refits `object`'s
#'   exact model, and `summary$pct` is always ~100. When supplied (a
#'   one-sided formula, character vector, or -- for multiphase fits -- a
#'   named list of one-sided formulas keyed by phase, matching
#'   [hzr_stepwise()]'s `scope`), each replicate runs a fresh
#'   [hzr_stepwise()] selection instead; see Details.
#' @param criterion Entry / retention rule passed through to
#'   [hzr_stepwise()] on each replicate when `scope` is supplied; ignored
#'   when `scope = NULL`. One of `"score"` (default), `"wald"`, or `"aic"`.
#'   `"score"` reproduces C/SAS HAZARD's `SELECTION` statistic and needs no
#'   per-candidate refit, which is what makes a bootstrap screen over many
#'   candidates tractable. Following SAS, the variance used during
#'   *selection* is approximate (shaping-parameter covariances are ignored);
#'   final-model standard errors are unaffected. For single-distribution
#'   fits, `"score"` computes the observed information numerically via the
#'   suggested \pkg{numDeriv} package and errors if it is not installed; a
#'   multiphase fit uses the analytic Hessian instead and does not need it.
#'   See [hzr_stepwise()].
#' @param direction,slentry,slstay,max_steps,max_move,force_in,force_out
#'   Passed through to [hzr_stepwise()] on each replicate when `scope` is
#'   supplied; ignored when `scope = NULL`. See [hzr_stepwise()] for
#'   definitions and defaults.
#' @param ... Additional arguments forwarded to [hzr_stepwise()] (e.g.
#'   `control = list(n_starts = 1)`) when `scope` is supplied; ignored
#'   otherwise.
#'
#' @section Selection mode is experimental:
#'
#' Everything reached through `scope` -- the selection arguments, and the
#' `summary$pct` selection frequencies they produce -- is new and should be
#' treated as unstable. The fixed-formula bootstrap (`scope = NULL`) is not
#' affected and has been stable since 0.9.3.
#'
#' Two reasons to expect change. The first is that the design is still being
#' read off real runs rather than settled in advance: production screens have
#' already moved the defaults once and turned up several ways a screen could
#' report success while selecting nothing.
#'
#' The second is scale, and it is the one to plan around. A screen over a
#' large candidate pool runs for hours, and this function writes nothing
#' until its final replicate, so a run that dies late loses everything. There
#' is no built-in way to split one screen across processes and combine the
#' parts. If you are running at that scale, drive `hzr_bootstrap()` in chunks
#' from your own script and pool the replicates yourself -- deriving each
#' chunk's seed from its chunk number, offsetting replicate ids so a variable
#' selected in two chunks is not counted once, and recomputing frequencies
#' from the pooled replicates rather than averaging across chunks. Whatever
#' eventually covers that inside the package may well change this function's
#' interface.
#'
#' @return A list with class `"hzr_bootstrap"` containing:
#' \describe{
#'   \item{replicates}{Data frame with columns `replicate`, `parameter`,
#'     and `estimate` -- one row per parameter per successful replicate.}
#'   \item{summary}{Data frame with columns `parameter`, `n`, `pct`,
#'     `mean`, `sd`, `min`, `max`, `ci_lower`, `ci_upper` -- one row per
#'     parameter. In `mode = "select"`, `pct` is the selection frequency
#'     and the other statistics are conditional on selection.}
#'   \item{n_success}{Number of successfully converged replicates.}
#'   \item{n_failed}{Number of replicates that failed to converge.}
#'   \item{n_uncomputable_replicates}{Select mode only: number of otherwise
#'     successful replicates whose screen stopped because no remaining
#'     candidate's score statistic could be computed, rather than because no
#'     candidate met `slentry`. Such replicates contribute no selections, so
#'     a non-zero count means every reported selection frequency is
#'     depressed. Always `0` in refit mode.}
#'   \item{mode}{`"refit"` (fixed-formula bootstrap) or `"select"`
#'     (embedded stepwise selection).}
#'   \item{scope}{Only present when `mode == "select"`: the candidate
#'     scope used.}
#' }
#'
#' @examples
#' \donttest{
#' data(avc)
#' avc <- na.omit(avc)
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ age + mal,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5, 0, 0),
#'   fit   = TRUE
#' )
#' bs <- hzr_bootstrap(fit, n_boot = 50, seed = 123)
#' print(bs)
#'
#' # Embedded stepwise selection: screen candidate covariates for how
#' # often they enter the model across resamples (R equivalent of SAS
#' # %HAZBOOT).
#' base <- hazard(
#'   survival::Surv(int_dead, dead) ~ 1,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5),
#'   fit   = TRUE
#' )
#' bs_sel <- hzr_bootstrap(base, n_boot = 20, seed = 123,
#'                          scope = ~ age + mal,
#'                          slentry = 0.3, slstay = 0.2)
#' print(bs_sel)
#' }
#'
#' @seealso [hazard()] for model fitting, [vcov.hazard()] for
#'   Hessian-based standard errors, [hzr_stepwise()] for the selection
#'   procedure used when `scope` is supplied.
#' @export
hzr_bootstrap <- function(object, n_boot = 200L, fraction = 1.0,
                           seed = NULL, verbose = FALSE,
                           scope = NULL,
                           direction = c("both", "forward", "backward"),
                           criterion = c("score", "wald", "aic"),
                           slentry = 0.30, slstay = 0.20,
                           max_steps = 50L, max_move = 4L,
                           force_in = character(), force_out = character(),
                           ...) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a fitted hazard object.", call. = FALSE)
  }
  if (is.null(object$fit$theta) ||
      (is.logical(object$fit$converged) && is.na(object$fit$converged))) {
    stop("'object' has no fitted parameters. Refit with fit = TRUE.",
         call. = FALSE)
  }

  n_boot <- as.integer(n_boot)
  if (n_boot < 1L) stop("'n_boot' must be at least 1.", call. = FALSE)
  if (fraction <= 0 || fraction > 1) {
    stop("'fraction' must be in (0, 1].", call. = FALSE)
  }

  direction <- match.arg(direction)
  criterion <- match.arg(criterion)
  select_mode <- !is.null(scope)

  # `...` exists only to forward stepwise-control arguments (e.g. `control=`)
  # to hzr_stepwise() in select-mode; fixed-refit mode (scope = NULL) has no
  # use for it. Without this check, a mistyped argument (e.g. `verbsoe =
  # TRUE`) would silently be accepted and ignored instead of erroring, as it
  # did before `...` was added for `scope=`.
  extra_args <- list(...)
  if (!select_mode && length(extra_args) > 0) {
    stop("Unused argument(s) in '...': ",
         paste(names(extra_args), collapse = ", "),
         ". '...' is only forwarded to hzr_stepwise() when 'scope' is set.",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # hzr_stepwise() is always called below with trace = FALSE (per-step
  # stepwise output would be too noisy across n_boot replicates; `verbose`
  # controls bootstrap-level progress instead). Strip `trace` from `...`
  # first so a caller-supplied `trace=` doesn't collide with it ("formal
  # argument matched by multiple actual arguments").
  extra_args$trace <- NULL

  # Reconstruct the call components
  cl <- object$call
  # Prefer the evaluated `data` argument stored on the fitted object
  # (`object$data$frame` -- the data frame passed to hazard(), not a
  # model.frame() result): it is guaranteed available and needs no caller-frame
  # lookup. Re-evaluating `cl$data` in `parent.frame()` is fragile -- the
  # original `data` symbol may no longer be in scope (e.g. the fit was built
  # inside a helper that has returned) -- so fall back to it only for objects
  # fitted before the frame was stored on `object$data`.
  orig_data <- if (!is.null(object$data$frame)) {
    object$data$frame
  } else {
    eval(cl$data, envir = parent.frame())
  }
  n_obs <- nrow(orig_data)
  sample_size <- max(1L, as.integer(n_obs * fraction))

  # Observation weights, if any, must be resampled in lockstep with the data.
  # Prefer the weights already evaluated and stored on the fitted object: they
  # are guaranteed aligned with the original rows and need no caller-frame
  # lookup. Re-evaluating `cl$weights` in `parent.frame()` is fragile -- the
  # original symbol/expression may no longer be in scope, or may now resolve to
  # a different value -- so fall back to it only for objects fitted before
  # weights were stored on `object$data`. Each replicate is then rewired to a
  # locally bound, resampled copy (mirroring `data`).
  orig_weights <- if (!is.null(object$data$weights)) {
    object$data$weights
  } else if (!is.null(cl$weights)) {
    eval(cl$weights, envir = parent.frame())
  } else {
    NULL
  }

  # Select-mode: fail loud on a structurally invalid scope (wrong type, an
  # unnamed list for a multiphase fit, an unknown phase name) up front, by
  # running one stepwise search against the untouched data, instead of
  # surfacing it n_boot replicates later. This does NOT catch a merely
  # nonexistent column name -- hzr_stepwise() itself converts that into a
  # per-candidate warning() rather than an error (see
  # inst/dev/BOOTSTRAP-SELECTION-DESIGN.md, "Error handling").
  if (select_mode) {
    # Muffle only .hzr_safe_solve()'s numerical post-fit warnings
    # (ill-conditioned / non-invertible / non-positive-definite Hessian,
    # non-positive variance) -- the per-fit noise this screen aggregates
    # over, which would otherwise fire on the real-data selected model here.
    # Identify them by their originating call, not message text, so all of
    # them are caught (including messages without the word "Hessian") while
    # unrelated warnings still surface: a mistyped scope column
    # ("candidate refit failed for ..."), the optimizer's
    # non-conformant-Hessian note, and all errors -- so a bad scope is caught
    # once, up front.
    withCallingHandlers(
      do.call(hzr_stepwise, c(
        list(
          object, scope = scope, data = orig_data,
          direction = direction, criterion = criterion,
          slentry = slentry, slstay = slstay,
          max_steps = max_steps, max_move = max_move,
          force_in = force_in, force_out = force_out,
          trace = FALSE
        ),
        extra_args
      )),
      warning = function(w) {
        if (any(grepl("hzr_safe_solve", deparse(conditionCall(w)),
                      fixed = TRUE))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  # VECTOR-INTERFACE FITS.
  #
  # hazard() accepts either a formula plus `data`, or bare `time`/`status`
  # vectors. Resampling `data` alone is enough for the formula interface, but
  # NOT for the vector one: the stored call holds `time = d$col` as an
  # *expression*, so each replicate re-evaluates it against the ORIGINAL data
  # and returns the original fit. That produced n_success = n_boot, n_failed
  # = 0, no warning, and n_boot identical replicates -- a bootstrap summary
  # that looks complete and contains nothing.
  #
  # The evaluated vectors are already stored on the object, so they can be
  # resampled by the same index and rewired the same way `data` and `weights`
  # are. `x` is excluded deliberately: a design matrix supplied that way is
  # rebuilt from `data`/`scope` per replicate.
  vector_interface <- is.null(cl$formula) && !is.null(cl$time)
  vec_args <- c("time", "status", "time_lower", "time_upper")
  vec_orig <- if (vector_interface) {
    stats::setNames(lapply(vec_args, function(a) object$data[[a]]), vec_args)
  } else {
    NULL
  }
  if (vector_interface) {
    # EVERY vector argument the call actually passed must have a stored,
    # correctly-sized copy. Rewiring only some of them is worse than rewiring
    # none: the rewired arguments follow the resample while the rest still
    # evaluate against the original data, so row i's time gets paired with
    # row j's status or interval bound. That is silent corruption producing
    # plausible numbers, not an error.
    passed <- vec_args[vapply(vec_args, function(a) !is.null(cl[[a]]), logical(1))]
    have   <- vapply(vec_orig, function(v) !is.null(v) && length(v) == n_obs,
                     logical(1))
    missing_vecs <- passed[!have[passed]]
    if (length(missing_vecs)) {
      stop("hzr_bootstrap(): this fit was built with the vector interface, ",
           "but the evaluated vector(s) ",
           paste(sQuote(missing_vecs), collapse = ", "),
           " are not stored on the object, so replicates cannot be resampled ",
           "consistently. Refit with the formula interface ",
           "(Surv(...) ~ ., data = ...) and bootstrap that.", call. = FALSE)
    }
    vec_orig <- vec_orig[passed]
  }

  # Parameter names from the fitted model. In fixed-refit mode every
  # replicate shares the same theta layout, so names are resolved once, up
  # front. In select-mode, each replicate can select a different variable
  # set, so names are resolved per replicate inside the loop instead.
  param_names <- if (!select_mode) .hzr_bootstrap_param_names(object) else NULL

  # Replicate calls must resolve two different things: `boot_data`/`boot_weights`
  # (locals here) and every other argument the user passed by symbol (theta,
  # phases, control), which only exist in the environment the call was written
  # in. A child of that environment carrying the resample bindings resolves
  # both. Falls back to parent.frame() for objects fitted before call_env was
  # stored.
  eval_env <- object$call_env %||% parent.frame()

  # Accumulate results
  rep_list <- vector("list", n_boot)
  n_success <- 0L
  n_failed <- 0L
  # Replicates whose stepwise screen stopped because no candidate's score
  # could be computed.  Each replicate runs under suppressWarnings(), so the
  # step-level warning never reaches the user here; the count has to be read
  # off the returned objects and reported in aggregate.
  n_uncomputable_reps <- 0L

  # Progress bar over replicates (verbose only). Closed after the loop.
  pb <- if (verbose) utils::txtProgressBar(min = 0, max = n_boot, style = 3)
  if (verbose) on.exit(close(pb), add = TRUE)

  for (b in seq_len(n_boot)) {
    # Resample with replacement
    idx <- sample.int(n_obs, size = sample_size, replace = TRUE)
    boot_data <- orig_data[idx, , drop = FALSE] # nolint: object_usage_linter.
    # boot_weights is referenced via quote() inside eval -- lintr cannot trace it
    boot_weights <- if (is.null(orig_weights)) NULL else orig_weights[idx] # nolint: object_usage_linter.

    # Fresh child of the fitting environment per replicate, carrying this
    # replicate's resample bindings (referenced via quote() in the call below).
    rep_env <- new.env(parent = eval_env)
    assign("boot_data", boot_data, envir = rep_env)
    if (!is.null(orig_weights)) assign("boot_weights", boot_weights, envir = rep_env)
    # Vector-interface arguments follow the same index as the rows.
    if (vector_interface) {
      for (a in names(vec_orig)) {
        assign(paste0("boot_", a), vec_orig[[a]][idx], envir = rep_env)
      }
    }

    # Per-replicate fits routinely hit ill-conditioned Hessians and other
    # numerical warnings on individual resamples; these are not individually
    # actionable (the bootstrap aggregates over replicates) and would swamp
    # the console over n_boot fits. Suppress them here -- structural problems
    # (e.g. a mistyped `scope` column) still surface once from the up-front
    # validation call above, and hard failures are caught below and counted.
    if (select_mode) {
      # Refit the (shape-fixed) base model on the resampled data first, so
      # the stepwise search's entry/retention tests compare candidates
      # against a base likelihood computed on the SAME resampled data --
      # then run a fresh stepwise selection from that base.
      boot_fit <- tryCatch(
        suppressWarnings({
          cl_base <- cl
          cl_base$data <- quote(boot_data)
          if (!is.null(orig_weights)) cl_base$weights <- quote(boot_weights)
          for (a in names(vec_orig)) {
            cl_base[[a]] <- as.name(paste0("boot_", a))
          }
          cl_base$fit <- TRUE
          base_boot <- eval(cl_base, envir = rep_env)
          if (!is.finite(base_boot$fit$objective)) {
            stop("base refit did not converge")
          }
          do.call(hzr_stepwise, c(
            list(
              base_boot, scope = scope, data = boot_data,
              direction = direction, criterion = criterion,
              slentry = slentry, slstay = slstay,
              max_steps = max_steps, max_move = max_move,
              force_in = force_in, force_out = force_out,
              trace = FALSE
            ),
            extra_args
          ))
        }),
        error = function(e) NULL
      )
    } else {
      # Refit using the same call but with resampled data (and weights, if any)
      # (boot_data/boot_weights are referenced via quote() inside eval)
      boot_fit <- tryCatch(
        suppressWarnings({
          cl_boot <- cl
          cl_boot$data <- quote(boot_data)
          if (!is.null(orig_weights)) cl_boot$weights <- quote(boot_weights)
          for (a in names(vec_orig)) {
            cl_boot[[a]] <- as.name(paste0("boot_", a))
          }
          cl_boot$fit <- TRUE
          eval(cl_boot, envir = rep_env)
        }),
        error = function(e) NULL
      )
    }

    if (!is.null(boot_fit) && is.finite(boot_fit$fit$objective)) {
      n_success <- n_success + 1L
      if (select_mode &&
            isTRUE(boot_fit$criteria$stopped_uncomputable)) {
        n_uncomputable_reps <- n_uncomputable_reps + 1L
      }
      theta_b <- boot_fit$fit$theta
      names_b <- if (select_mode) {
        .hzr_bootstrap_param_names(boot_fit)
      } else {
        param_names[seq_along(theta_b)]
      }
      rep_list[[b]] <- data.frame(
        replicate = b,
        parameter = names_b,
        estimate  = as.numeric(theta_b),
        stringsAsFactors = FALSE
      )
    } else {
      n_failed <- n_failed + 1L
    }

    if (verbose) utils::setTxtProgressBar(pb, b)
  }

  # Combine replicates
  replicates <- do.call(rbind, Filter(Negate(is.null), rep_list))
  if (is.null(replicates)) {
    replicates <- data.frame(replicate = integer(0),
                              parameter = character(0),
                              estimate = numeric(0),
                              stringsAsFactors = FALSE)
  }
  rownames(replicates) <- NULL

  # Summary statistics per parameter
  if (nrow(replicates) > 0) {
    summary_list <- lapply(split(replicates, replicates$parameter), function(d) {
      data.frame(
        parameter = d$parameter[1],
        n         = nrow(d),
        pct       = 100 * nrow(d) / n_boot,
        mean      = mean(d$estimate),
        sd        = stats::sd(d$estimate),
        min       = min(d$estimate),
        max       = max(d$estimate),
        ci_lower  = stats::quantile(d$estimate, 0.025),
        ci_upper  = stats::quantile(d$estimate, 0.975),
        stringsAsFactors = FALSE
      )
    })
    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL
    if (select_mode) {
      # No single canonical variable order across replicates (different
      # replicates select different sets) -- rank by selection frequency,
      # most-selected first, ties broken alphabetically.
      summary_df <- summary_df[order(-summary_df$pct, summary_df$parameter), ]
    } else {
      # Sort by parameter order in the original model
      idx_order <- match(summary_df$parameter, param_names)
      summary_df <- summary_df[order(idx_order), ]
    }
  } else {
    summary_df <- data.frame(parameter = character(0), n = integer(0),
                              pct = numeric(0), mean = numeric(0),
                              sd = numeric(0), min = numeric(0),
                              max = numeric(0), ci_lower = numeric(0),
                              ci_upper = numeric(0),
                              stringsAsFactors = FALSE)
  }

  # A selection frequency is the whole deliverable of a select-mode run, so a
  # screen that never selected anything is the result being empty rather than
  # a quiet edge case.  It cannot be read off the object either: the base
  # model's own parameters appear in every replicate by construction, so they
  # fill the summary at pct = 100 and the table looks like a set of perfectly
  # reliable variables.
  if (select_mode && n_success > 0L) {
    selected <- setdiff(unique(replicates$parameter), names(coef(object)))
    if (length(selected) == 0L) {
      warning("Bootstrap selection selected no covariate in any of the ",
              n_success, " successful replicates. The summary holds only the ",
              "base model's parameters, each at pct = 100. Common causes: ",
              "`slentry` stricter than intended; a `scope` naming columns ",
              "absent from the data; or a base fit whose stored call cannot ",
              "be rewritten for the refit.", call. = FALSE)
    }
  }

  if (n_uncomputable_reps > 0L) {
    warning(n_uncomputable_reps, " of ", n_success, " successful replicates ",
            "stopped because the score statistic could not be computed for ",
            "any remaining candidate, rather than because no candidate met ",
            "`slentry`. Those replicates contribute no selections, so every ",
            "reported selection frequency is depressed by them. Usual causes ",
            "are a degenerate or collinear candidate column and resamples on ",
            "which the information matrix is not invertible.", call. = FALSE)
  }

  result <- list(
    replicates = replicates,
    summary    = summary_df,
    n_success  = n_success,
    n_failed   = n_failed,
    n_uncomputable_replicates = n_uncomputable_reps,
    mode       = if (select_mode) "select" else "refit"
  )
  if (select_mode) result$scope <- scope
  class(result) <- "hzr_bootstrap"
  result
}

#' @rdname hzr_bootstrap
#' @param x An `hzr_bootstrap` object.
#' @param digits Number of decimal places for formatting.
#' @export
print.hzr_bootstrap <- function(x, digits = 4, ...) {
  cat("Bootstrap inference for hazard model\n")
  mode <- x$mode %||% "refit"
  cat("Mode:", if (identical(mode, "select")) {
    "embedded stepwise selection"
  } else {
    "fixed refit"
  }, "\n")
  cat("Replicates:", x$n_success, "successful,", x$n_failed, "failed\n\n")
  if (nrow(x$summary) > 0) {
    display <- x$summary
    for (col in c("pct", "mean", "sd", "min", "max",
                   "ci_lower", "ci_upper")) {
      display[[col]] <- round(display[[col]], digits)
    }
    print(display, row.names = FALSE)
  }
  invisible(x)
}


# =========================================================================
# hzr_competing_risks -- Cumulative incidence with Greenwood variance
# =========================================================================

#' Competing risks cumulative incidence
#'
#' Compute cumulative incidence functions for multiple competing event
#' types using the Aalen-Johansen estimator with Greenwood variance.
#' This is the R equivalent of the SAS `markov.sas` macro.
#'
#' Unlike the naive 1 - KM estimator (which overestimates incidence when
#' competing risks exist), this provides the correct marginal cumulative
#' incidence for each event type.
#'
#' @note This estimator is **unweighted**: every observation contributes a
#'   unit count to the at-risk and event tallies. There is no `weights`
#'   argument, so case or inverse-probability weights are not yet supported
#'   for competing-risks incidence (unlike [hazard()], which accepts
#'   `weights`). Pre-expand weighted rows to individual records if an
#'   approximate weighted estimate is needed.
#'
#' @param time Numeric vector of follow-up times.
#' @param event Integer vector of event type indicators:
#'   0 = censored, 1 = event type 1, 2 = event type 2, etc.
#'
#' @return A data frame with one row per unique event time and columns:
#' \describe{
#'   \item{time}{Event time.}
#'   \item{n_risk}{Number at risk.}
#'   \item{n_event_1, n_event_2, ...}{Events of each type at this time.}
#'   \item{n_censor}{Number censored at this time.}
#'   \item{surv}{Overall event-free survival (freedom from all events).}
#'   \item{incid_1, incid_2, ...}{Cumulative incidence for each event type.}
#'   \item{se_surv}{Standard error of overall survival.}
#'   \item{se_1, se_2, ...}{Standard error of each cumulative incidence.}
#' }
#'
#' @examples
#' data(valves)
#' valves_cc <- na.omit(valves)
#' # Combine death and PVE into a competing risks event variable
#' # 0 = censored, 1 = death, 2 = PVE
#' event_cr <- ifelse(valves_cc$dead == 1, 1L,
#'                    ifelse(valves_cc$pve == 1, 2L, 0L))
#' time_cr <- pmin(valves_cc$int_dead, valves_cc$int_pve)
#' cr <- hzr_competing_risks(time_cr, event_cr)
#' head(cr)
#'
#' @references
#' Aalen O, Johansen S (1978). An empirical transition matrix for
#' non-homogeneous Markov chains based on censored observations.
#' *Scand J Statist* 5(3):141--150.
#'
#' Kalbfleisch JD, Prentice RL (1980). *The Statistical Analysis of Failure
#' Time Data.* Wiley, New York.
#' @seealso [hzr_kaplan()] for single-event survival estimation.
#' @export
hzr_competing_risks <- function(time, event) {
  if (!is.numeric(time) || !is.numeric(event)) {
    stop("'time' and 'event' must be numeric vectors.", call. = FALSE)
  }
  if (length(time) != length(event)) {
    stop("'time' and 'event' must have the same length.", call. = FALSE)
  }

  event_types <- sort(setdiff(unique(event), 0))
  n_types <- length(event_types)
  if (n_types == 0) {
    stop("No events found (all observations are censored).", call. = FALSE)
  }

  # Sort by time
  ord <- order(time)
  time_s <- time[ord]
  event_s <- event[ord]
  n_total <- length(time_s)

  # Unique event times (where at least one event occurred)
  u_times <- sort(unique(time_s[event_s > 0]))
  out_n <- length(u_times)

  # Initialize output
  n_risk <- numeric(out_n)
  n_censor_out <- numeric(out_n)
  n_event_mat <- matrix(0, nrow = out_n, ncol = n_types)
  colnames(n_event_mat) <- paste0("n_event_", event_types)

  # Cumulative incidence vectors
  surv <- numeric(out_n)
  incid <- matrix(0, nrow = out_n, ncol = n_types)
  colnames(incid) <- paste0("incid_", event_types)

  # Variance (diagonal of Greenwood matrix -- simplified)
  var_surv <- numeric(out_n)
  var_incid <- matrix(0, nrow = out_n, ncol = n_types)

  prev_surv <- 1.0
  prev_incid <- rep(0, n_types)
  prev_var_surv <- 0
  prev_var_incid <- rep(0, n_types)
  at_risk <- n_total

  for (k in seq_along(u_times)) {
    t_k <- u_times[k]

    # Censored before this time (between previous event time and t_k)
    prev_t <- if (k == 1) 0 else u_times[k - 1]
    censored_between <- sum(event_s == 0 & time_s > prev_t & time_s < t_k)
    at_risk <- at_risk - censored_between

    n_risk[k] <- at_risk

    # Events and censored AT this time
    at_time <- time_s == t_k
    for (j in seq_along(event_types)) {
      n_event_mat[k, j] <- sum(event_s == event_types[j] & at_time)
    }
    n_censor_out[k] <- sum(event_s == 0 & at_time)
    total_events_k <- sum(n_event_mat[k, ])

    # Aalen-Johansen update
    if (at_risk > 0) {
      # Transition probabilities
      d_total <- total_events_k / at_risk
      surv[k] <- prev_surv * (1 - d_total)

      for (j in seq_along(event_types)) {
        d_j <- n_event_mat[k, j] / at_risk
        incid[k, j] <- prev_incid[j] + prev_surv * d_j
      }

      # Greenwood variance update (simplified diagonal)
      if (at_risk > total_events_k && at_risk > 0) {
        greenwood_term <- total_events_k / (at_risk * (at_risk - total_events_k))
        var_surv[k] <- surv[k]^2 *
          (prev_var_surv / max(prev_surv^2, .Machine$double.xmin) +
             greenwood_term)

        for (j in seq_along(event_types)) {
          d_j <- n_event_mat[k, j] / at_risk
          var_incid[k, j] <- prev_var_incid[j] +
            prev_surv^2 * d_j * (1 - d_j) / at_risk
        }
      } else {
        var_surv[k] <- prev_var_surv
        var_incid[k, ] <- prev_var_incid
      }
    } else {
      surv[k] <- prev_surv
      incid[k, ] <- prev_incid
      var_surv[k] <- prev_var_surv
      var_incid[k, ] <- prev_var_incid
    }

    # Decrease at-risk by events + censored AT this time
    at_risk <- at_risk - total_events_k - n_censor_out[k]

    prev_surv <- surv[k]
    prev_incid <- incid[k, ]
    prev_var_surv <- var_surv[k]
    prev_var_incid <- var_incid[k, ]
  }

  # Assemble result
  result <- data.frame(
    time     = u_times,
    n_risk   = n_risk,
    n_event_mat,
    n_censor = n_censor_out,
    surv     = surv,
    incid,
    se_surv  = sqrt(pmax(var_surv, 0)),
    stringsAsFactors = FALSE
  )

  # Add SE columns for each event type
  for (j in seq_along(event_types)) {
    result[[paste0("se_", event_types[j])]] <-
      sqrt(pmax(var_incid[, j], 0))
  }

  class(result) <- c("hzr_competing_risks", "data.frame")
  result
}

#' @rdname hzr_competing_risks
#' @param x An `hzr_competing_risks` object.
#' @param digits Number of decimal places for formatting.
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_competing_risks <- function(x, digits = 4, ...) {
  incid_cols <- grep("^incid_", names(x), value = TRUE)
  cat("Competing risks cumulative incidence\n")
  cat("Event types:", length(incid_cols), " | Time points:", nrow(x), "\n")
  if (nrow(x) > 0) {
    cat("Final survival:", round(x$surv[nrow(x)], digits), "\n")
    for (col in incid_cols) {
      cat("Final", col, ":", round(x[[col]][nrow(x)], digits), "\n")
    }
  }
  cat("\n")
  display <- x
  for (col in c("surv", incid_cols, "se_surv",
                 grep("^se_\\d", names(x), value = TRUE))) {
    if (col %in% names(display)) {
      display[[col]] <- round(display[[col]], digits)
    }
  }
  class(display) <- "data.frame"
  print(utils::head(display, 15), row.names = FALSE)
  if (nrow(display) > 15) cat("... (", nrow(display) - 15, " more rows)\n")
  invisible(x)
}
