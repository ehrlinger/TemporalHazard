#' @importFrom stats optim pnorm
#' @importFrom survival Surv
#' @keywords internal
NULL

# hazard_api.R -- Primary user-facing API for TemporalHazard
#
# This file contains the main entry points:
#   hazard()         -- construct and optionally fit a parametric hazard model
#   predict.hazard() -- generate predictions from a fitted hazard object
#   print.hazard()   -- compact S3 print method
#   coef.hazard()    -- extract fitted parameter vector
#   vcov.hazard()    -- extract variance-covariance matrix
#
# DISTRIBUTION PARAMETERIZATION CONVENTIONS
# -----------------------------------------
# Each distribution stores a flattened theta vector.  Shape/scale parameters
# always come first; covariate coefficients (beta) follow:
#
#   Weibull:     theta = [mu, nu, beta_1, beta_2, ...]      mu>0, nu>0
#   Exponential: theta = [log(lambda), beta_1, beta_2, ...]      lambda>0 via exp()
#   Loglogistic: theta = [log(alpha), log(beta), beta_1, ...]   alpha>0, beta>0 via exp()
#   Lognormal:   theta = [mu, log(sigma), beta_1, beta_2, ...]   sigma>0 via exp(); AFT model
#
# ADDING A NEW DISTRIBUTION
# -------------------------
# 1. Create R/likelihood-<dist>.R with .hzr_logl_<dist>(),
#    .hzr_gradient_<dist>(), and .hzr_optim_<dist>().
# 2. Add an else-if branch in hazard() below under "Distribution dispatch".
# 3. Add an else-if branch in predict.hazard() under "Dispatch by distribution".
# 4. Update the supported-distributions guard in predict.hazard().
# 5. Update the coefficient-extraction block in predict.hazard() for
#    linear_predictor / hazard prediction types.
# 6. Write tests in tests/testthat/test-<dist>-dist.R.
# 7. Generate a golden fixture via data-raw/golden_fixtures.R
#    (.hzr_create_<dist>_golden_fixture()).

#' Build and optionally fit a hazard model
#'
#' Creates a `hazard` object and optionally fits it via maximum likelihood.
#' This mirrors the argument-oriented workflow of the legacy HAZARD C/SAS
#' implementation: supply starting values in `theta` and the function will
#' optimize to produce fitted estimates.
#'
#' @section The fitted model:
#'
#' For `dist = "multiphase"` the cumulative hazard, instantaneous hazard, and
#' survival are
#'
#' \deqn{H(t \mid \mathbf{x}) = \sum_{j=1}^{J} \mu_j(\mathbf{x}) \, \Phi_j(t),
#'       \qquad
#'       h(t \mid \mathbf{x}) = \sum_{j=1}^{J} \mu_j(\mathbf{x}) \, \varphi_j(t),
#'       \qquad
#'       S(t \mid \mathbf{x}) = \exp\!\bigl(-H(t \mid \mathbf{x})\bigr)}
#'
#' where \eqn{\mu_j(\mathbf{x}) = \exp(\alpha_j + \mathbf{x}_j^\top
#' \beta_j)} and the temporal shapes \eqn{\Phi_j}, \eqn{\varphi_j}
#' are set by each phase's `type` (see [hzr_phase()]).  The
#' proportional-hazards single-phase families (`"weibull"`, `"exponential"`)
#' are the special case \eqn{J = 1}, with covariates acting multiplicatively on
#' one temporal shape.  The `"loglogistic"` (proportional-odds) and
#' `"lognormal"` (accelerated-failure-time) families place covariates
#' differently --- on the odds of failure and the log-time location,
#' respectively --- so they are separate parameterizations, not special cases
#' of this additive form.  Parameters are estimated
#' on an unconstrained internal scale (e.g. \eqn{\log\mu}, \eqn{\log t_{1/2}})
#' and transformed back for reporting; see
#' `vignette("mf-mathematical-foundations")`.
#'
#' @section Baseline distributions:
#'
#' The `dist` argument selects the parametric form of the baseline hazard.  The
#' four single-distribution families differ in the *shape* the
#' hazard traces over follow-up; choose by what the risk is expected to do over
#' time.  `"multiphase"` is the general additive model that lets several such
#' shapes coexist.
#'
#' \describe{
#'   \item{`"weibull"` --- monotone rising or falling hazard (default)}{The
#'     workhorse parametric model: \eqn{H(t \mid \mathbf{x}) = (\mu t)^\nu
#'     \exp(\eta)}, with hazard \eqn{h \propto t^{\nu - 1}}.  The single shape
#'     \eqn{\nu} makes risk increase over time (\eqn{\nu > 1}), decrease
#'     (\eqn{\nu < 1}), or stay flat (\eqn{\nu = 1}).  Use it as the default when
#'     a single monotone trend describes the hazard.}
#'   \item{`"exponential"` --- constant hazard}{The memoryless special case
#'     \eqn{\nu = 1}: a time-invariant baseline rate, \eqn{H(t \mid \mathbf{x}) =
#'     \mu t \exp(\eta)}.  Use it when the event rate does not change with
#'     follow-up time (the constant background risk also appears as the
#'     `"constant"` phase in a multiphase model).}
#'   \item{`"loglogistic"` --- unimodal (rise-then-fall) hazard}{A log-logistic
#'     proportional-odds form (covariates act multiplicatively on the odds of
#'     failure, \eqn{\exp(\eta)}, not as an AFT time shift) whose hazard rises to
#'     a single peak and then declines when the shape exceeds 1 (and is monotone
#'     decreasing otherwise), with heavier tails than the log-normal.  Use it
#'     when risk climbs to an early peak and then eases off.}
#'   \item{`"lognormal"` --- early-peaking, resolving hazard}{An
#'     accelerated-failure-time form in which \eqn{\log} time is Gaussian; the
#'     hazard rises to an early peak and then decays toward zero.  Use it for
#'     risk that is concentrated early and resolves over time.}
#'   \item{`"multiphase"` --- additive N-phase hazard}{Sums several phase shapes
#'     into one model, \eqn{H = \sum_j \mu_j(\mathbf{x}) \Phi_j(t)}, so the
#'     overall hazard can fall, level off, and rise again within one fit.
#'     Requires `phases`; see [hzr_phase()] for the available phase shapes.  This
#'     is the form that reproduces the classic C/SAS HAZARD models.}
#' }
#'
#' @param time Numeric follow-up time vector.
#' @param status Numeric or logical event indicator vector.
#' @param time_lower Optional numeric lower bound vector for censoring intervals.
#'   Used when `status == 2` (interval-censored); defaults to `time` if NULL.
#' @param time_upper Optional numeric upper bound vector for censoring intervals.
#'   Used when `status %in% c(-1, 2)`; defaults to `time` if NULL.
#' @param x Optional design matrix (or data frame coercible to matrix).
#' @param formula Optional formula of the form `Surv(time, status) ~ predictors`.
#'   When provided, overrides direct time/status/x arguments and extracts from data.
#'   Example: `hazard(Surv(time, status) ~ x1 + x2, data = df, dist = "weibull", fit = TRUE)`.
#' @param data Optional data frame containing variables referenced in formula.
#' @param time_windows Optional numeric vector of strictly positive cut points for
#'   piecewise time-varying coefficients. When provided, each predictor column in
#'   `x` is expanded into one column per time window so each window gets its own
#'   coefficient.
#' @param theta Optional numeric coefficient vector (starting values for optimization).
#' @param dist Character baseline distribution label (default `"weibull"`).
#'   One of `"weibull"`, `"exponential"`, `"loglogistic"`, `"lognormal"`, or
#'   `"multiphase"`.  The single-distribution families differ in the *shape* the
#'   hazard traces over time; `"multiphase"` builds an additive N-phase hazard
#'   and requires `phases`.  See the **Baseline distributions** section for what
#'   each means and when to use it.
#' @param phases Optional named list of [hzr_phase()] objects specifying the
#'   phases for a multiphase model (`dist = "multiphase"`).  See Examples.
#' @param fit Logical; if TRUE, fit the model via maximum likelihood (default FALSE).
#' @param weights Optional numeric vector of observation weights (non-negative).
#'   Each observation's log-likelihood contribution is multiplied by its weight.
#'   Use for severity-weighted repeated events. Default `NULL` (unit weights).
#'   Implements the SAS `WEIGHT` statement.
#' @param control Named list of control options (see Details).
#' @param ... Additional named arguments retained for parity with legacy calling
#'   conventions.
#'
#' @details
#' Control parameters:
#' - `maxit`: Maximum iterations (default 1000)
#' - `n_starts`: Number of optimization starts for multiphase fits (default 5).
#'   Each start after the first adds random noise to the initial values, drawn
#'   from the ambient RNG stream; call `set.seed()` before fitting for
#'   reproducible results.
#' - `reltol`: Relative parameter change tolerance (default 1e-5)
#' - `abstol`: Absolute gradient norm tolerance (default 1e-6)
#' - `method`: Optimization method: "bfgs" or "nm" (default "bfgs")
#' - `condition`: Condition number control (default 14)
#' - `nocov`, `nocor`: Suppress covariance/correlation output (legacy; no-op in M2)
#'
#' Censoring status coding:
#' - 1: Exact event at time
#' - 0: Right-censored at time
#' - -1: Left-censored with upper bound at time_upper \(or time\)
#' - 2: Interval-censored in the interval \(time_lower, time_upper\)
#'
#' Time-varying coefficients:
#' - If `time_windows` is supplied, predictors are expanded to piecewise window
#'   interactions so each window has its own coefficient vector.
#' - This is implemented as design-matrix expansion, so the existing likelihood
#'   engines remain unchanged.
#'
#' @examples
#' # -- Univariable Weibull ----------------------------------------------
#' set.seed(1)
#' time   <- rexp(50, rate = 0.3)
#' status <- sample(0:1, 50, replace = TRUE, prob = c(0.3, 0.7))
#' fit <- hazard(time = time, status = status,
#'               theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
#' summary(fit)
#'
#' # -- Formula interface with covariates --------------------------------
#' set.seed(1001)
#' n   <- 180
#' dat <- data.frame(
#'   time   = rexp(n, rate = 0.35) + 0.05,
#'   status = rbinom(n, size = 1, prob = 0.6),
#'   age    = rnorm(n, mean = 62, sd = 11),
#'   nyha   = sample(1:4, n, replace = TRUE),
#'   shock  = rbinom(n, size = 1, prob = 0.18)
#' )
#'
#' fit2 <- hazard(
#'   survival::Surv(time, status) ~ age + nyha + shock,
#'   data    = dat,
#'   theta   = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
#'   dist    = "weibull",
#'   fit     = TRUE,
#'   control = list(maxit = 300)
#' )
#' summary(fit2)
#'
#' \donttest{
#' # -- Parametric survival with Kaplan-Meier overlay -----------------
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   # Parametric curve on a fine grid at median covariate profile
#'   t_grid   <- seq(0.05, max(dat$time), length.out = 80)
#'   curve_df <- data.frame(
#'     time = t_grid, age = median(dat$age), nyha = 2, shock = 0
#'   )
#'   curve_df$survival <- predict(fit2, newdata = curve_df,
#'                                type = "survival") * 100
#'
#'   # Kaplan-Meier empirical overlay
#'   km    <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
#'   km_df <- data.frame(time = km$time, survival = km$surv * 100)
#'
#'   ggplot() +
#'     geom_step(data = km_df, aes(time, survival, colour = "Kaplan-Meier")) +
#'     geom_line(data = curve_df, aes(time, survival,
#'                                    colour = "Parametric (Weibull)")) +
#'     scale_colour_manual(
#'       values = c("Parametric (Weibull)" = "#0072B2",
#'                  "Kaplan-Meier"         = "#D55E00")
#'     ) +
#'     scale_y_continuous(limits = c(0, 100)) +
#'     labs(x = "Months after surgery", y = "Freedom from death (%)",
#'          colour = NULL) +
#'     theme_minimal() +
#'     theme(legend.position = "bottom")
#' }
#' }
#'
#' \donttest{
#' # -- Multiphase model (two phases) ---------------------------------
#' fit_mp <- hazard(
#'   survival::Surv(time, status) ~ 1,
#'   data   = dat,
#'   dist   = "multiphase",
#'   phases = list(
#'     early = hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                        fixed = "shapes"),
#'     late  = hzr_phase("cdf", t_half = 5,   nu = 1, m = 0,
#'                        fixed = "shapes")
#'   ),
#'   fit     = TRUE,
#'   control = list(n_starts = 5, maxit = 1000)
#' )
#' summary(fit_mp)
#'
#' # -- Per-phase decomposed cumulative hazard ------------------------
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   t_grid <- seq(0.01, max(dat$time), length.out = 100)
#'   decomp <- predict(fit_mp, newdata = data.frame(time = t_grid),
#'                     type = "cumulative_hazard", decompose = TRUE)
#'
#'   df_long <- data.frame(
#'     time = rep(decomp$time, 3),
#'     cumhaz = c(decomp$total, decomp$early, decomp$late),
#'     component = rep(c("Total", "Early (cdf)", "Late (cdf)"),
#'                     each = nrow(decomp))
#'   )
#'   df_long$component <- factor(df_long$component,
#'     levels = c("Total", "Early (cdf)", "Late (cdf)"))
#'
#'   ggplot2::ggplot(df_long,
#'     ggplot2::aes(x = time, y = cumhaz, colour = component,
#'                  linewidth = component)) +
#'     ggplot2::geom_line() +
#'     ggplot2::scale_colour_manual(values = c(
#'       "Total" = "black", "Early (cdf)" = "#0072B2",
#'       "Late (cdf)" = "#D55E00"
#'     )) +
#'     ggplot2::scale_linewidth_manual(values = c(
#'       "Total" = 1.2, "Early (cdf)" = 0.6, "Late (cdf)" = 0.6
#'     )) +
#'     ggplot2::labs(
#'       x = "Time", y = "Cumulative hazard H(t)",
#'       colour = NULL, linewidth = NULL,
#'       title = "Multiphase decomposition: early + late"
#'     ) +
#'     ggplot2::theme_minimal() +
#'     ggplot2::theme(legend.position = "bottom")
#' }
#' }
#'
#' @seealso
#' [predict.hazard()] for survival/cumulative-hazard predictions,
#' [summary.hazard()] for model summaries,
#' [hzr_phase()] for specifying multiphase temporal shapes.
#'
#' Vignettes with worked examples:
#' \code{vignette("fitting-hazard-models")} --- single-phase through multiphase fitting,
#' \code{vignette("prediction-visualization")} --- prediction types and decomposed hazard plots,
#' \code{vignette("inference-diagnostics")} --- bootstrap CIs and model diagnostics.
#'
#' @references
#' Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of time-varying
#' hazard into phases, each incorporating a separate stream of concomitant
#' information. *J Am Stat Assoc.* 1986;81(395):615--624.
#' \doi{10.1080/01621459.1986.10478314}
#'
#' Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
#' Probability of atrial fibrillation after ablation: Using a parametric
#' nonlinear temporal decomposition mixed effects model. *Stat Methods Med Res.*
#' 2018;27(1):126--141. \doi{10.1177/0962280215623583}
#'
#' @return An object of class \code{hazard}, a named list with components:
#'   \code{call} (the matched call),
#'   \code{spec} (model specification: \code{dist}, \code{control},
#'   \code{time_windows}, \code{phases}),
#'   \code{data} (input data: \code{time}, \code{status}, \code{x},
#'   \code{weights}, etc.),
#'   \code{fit} (optimisation results: \code{theta}, \code{objective},
#'   \code{converged}, \code{se}, \code{vcov}, \code{counts}, \code{message};
#'   all \code{NULL} when \code{fit = FALSE}),
#'   and \code{engine} (implementation tag, \code{"native-r-m2"}).
#' @export
hazard <- function(formula = NULL,
                   data = NULL,
                   time = NULL,
                   status = NULL,
                   time_lower = NULL,
                   time_upper = NULL,
                   x = NULL,
                   time_windows = NULL,
                   theta = NULL,
                   dist = "weibull",
                   phases = NULL,
                   fit = FALSE,
                   weights = NULL,
                   control = list(),
                   ...) {
  # Formula dispatch: if formula is provided, parse it and extract time/status/x from data
  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("'data' is required when 'formula' is provided.", call. = FALSE)
    }

    # For multiphase models, the formula RHS may contain phase-scoped terms of
    # the form `phase_name(var1 + var2)`.  These are not valid R expressions
    # so model.matrix() would fail.  Strip them here: replace the RHS with `1`
    # (no global predictors) so that .hzr_parse_formula only extracts
    # time/status from the LHS.  Covariate routing is handled per-phase via
    # hzr_phase(formula = ...) and resolved inside .hzr_optim_multiphase().
    formula_for_parse <- formula
    if (!is.null(phases) && length(formula) >= 3L) {
      # Check if RHS contains phase-scoped calls of the form `phase_name(...)`.
      # Use a parse-tree walk (not string regex) to avoid false positives when
      # a phase name coincides with a base-R function (e.g., "log", "exp").
      if (.hzr_formula_has_phase_scope(formula[[3L]], names(phases))) {
        formula_for_parse <- stats::reformulate("1", response = formula[[2L]])
      }
    }

    parsed <- .hzr_parse_formula(formula = formula_for_parse, data = data)
    time <- parsed$time
    status <- parsed$status
    time_lower <- parsed$time_lower
    time_upper <- parsed$time_upper
    x <- parsed$x

  }

  # After formula dispatch, require time and status
  if (is.null(time) || is.null(status)) {
    stop("'time' and 'status' are required (either directly or via 'formula').", call. = FALSE)
  }
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
  }

  n <- length(time)
  if (length(status) != n) {
    stop("'status' must have the same length as 'time'.", call. = FALSE)
  }

  # Convert Surv object status to numeric if needed (after formula parsing)
  if (inherits(status, "Surv")) {
    status <- unclass(status)[, 2L]
  }

  # Optional censoring bounds:
  # - status = 1 (event): uses observed event time in `time`
  # - status = 0 (right-censored): censoring time in `time` (or `time_lower`)
  # - status = -1 (left-censored): upper bound in `time` (or `time_upper`)
  # - status = 2 (interval-censored): [time_lower, time_upper] required
  if (!is.null(time_lower)) {
    if (!is.numeric(time_lower) || length(time_lower) != n || any(!is.finite(time_lower)) || any(time_lower < 0)) {
      stop("'time_lower' must be a numeric vector of finite non-negative values matching length(time).", call. = FALSE)
    }
  }

  if (!is.null(time_upper)) {
    if (!is.numeric(time_upper) || length(time_upper) != n || any(!is.finite(time_upper)) || any(time_upper < 0)) {
      stop("'time_upper' must be a numeric vector of finite non-negative values matching length(time).", call. = FALSE)
    }
  }

  if (!is.null(x)) {
    x <- .hzr_as_design_matrix(x, n = n)
  }

  if (!is.null(time_windows)) {
    # We use cut points to define window-specific coefficient blocks.
    if (!is.numeric(time_windows) || any(!is.finite(time_windows)) || anyDuplicated(time_windows)) {
      stop("'time_windows' must be a numeric vector of unique finite cut points.", call. = FALSE)
    }
    time_windows <- sort(time_windows)
    if (any(time_windows <= 0)) {
      stop("'time_windows' cut points must be strictly positive.", call. = FALSE)
    }
    if (is.null(x)) {
      stop("'time_windows' requires predictor matrix 'x'.", call. = FALSE)
    }
  }

  x_fit <- x
  if (!is.null(time_windows) && !is.null(x)) {
    # Expand X -> [X_w1 | X_w2 | ...] where each row is active only in its window.
    x_fit <- .hzr_expand_time_varying_design(x = x, time = time, time_windows = time_windows)
  }

  if (!is.null(theta)) {
    if (!is.numeric(theta) || any(!is.finite(theta))) {
      stop("'theta' must be a finite numeric vector when provided.", call. = FALSE)
    }

    # If x exists, theta must include coefficients for all variates
    # theta = [shape parms ... | covariate coefficients ...]
    # For now, assume theta length determines whether we expect x
    required_coef <- if (is.null(x_fit)) 0L else ncol(x_fit)
    if (!is.null(x_fit) && length(theta) < required_coef) {
      stop("'theta' length must be >= number of required coefficients (", required_coef, ").", call. = FALSE)
    }
  }

  # --- Validate and normalize weights ----------------------------------------
  n_obs <- length(time)
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != n_obs) {
      stop("'weights' must be a numeric vector of length ", n_obs, ".",
           call. = FALSE)
    }
    if (any(weights < 0) || any(!is.finite(weights))) {
      stop("'weights' must be non-negative and finite.", call. = FALSE)
    }
  }

  if (!is.character(dist) || length(dist) != 1 || !nzchar(dist)) {
    stop("'dist' must be a non-empty character scalar.", call. = FALSE)
  }

  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }

  # Multiphase validation
  if (dist == "multiphase") {
    if (is.null(phases)) {
      stop("'phases' is required when dist = 'multiphase'. ",
           "Supply a list of hzr_phase() specifications.", call. = FALSE)
    }
    phases <- .hzr_validate_phases(phases)
  } else if (!is.null(phases)) {
    warning("'phases' is ignored when dist != 'multiphase'.")
    phases <- NULL
  }

  if (any(status %in% c(-1, 2))) {
    if (is.null(time_upper)) {
      warning("'time_upper' not provided; using 'time' as upper bound for left/interval-censored rows.")
    }
    if (any(status == 2) && is.null(time_lower)) {
      warning("'time_lower' not provided; using 'time' as lower bound for interval-censored rows.")
    }
  }

  # fit_state holds the result of optimization (or just starting values if fit=FALSE).
  # Fields:
  #   theta     -- parameter vector (starting values before fit; MLE estimates after)
  #   converged -- TRUE if optimizer reported convergence (code 0); NA if not fitted
  #   objective -- log-likelihood at theta; NA_real_ if not fitted
  #   se        -- standard errors (sqrt of vcov diagonal); NULL / NA if unavailable
  #   gradient  -- score vector at solution (populated by some optimizers); NULL otherwise
  #   vcov      -- variance-covariance matrix; NA if Hessian not invertible
  #   counts    -- c(fn, gr) evaluation counts from optim()
  #   message   -- convergence message string from optim()
  fit_state <- list(
    theta = theta,
    converged = NA,
    objective = NA_real_,
    se = NULL,
    gradient = NULL
  )

  .hzr_run_fit_safely <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        msg <- conditionMessage(w)
        # Muffle only benign numerical noise from optim().  The hardened-inversion
        # diagnostics (ill-conditioned / not positive-definite / not invertible /
        # non-finite) are deliberately allowed to surface so a fit that cannot
        # produce reliable standard errors is never silent.
        if (grepl("NaNs produced", msg, fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  .hzr_safe_se_from_vcov <- function(vcov_mat) {
    if (is.null(vcov_mat) || !is.matrix(vcov_mat) || anyNA(vcov_mat)) {
      return(NA)
    }
    d <- diag(vcov_mat)
    # Guard against small negative/invalid variances from numerical Hessians.
    if (!all(is.finite(d)) || any(d < 0)) {
      return(rep(NA_real_, length(d)))
    }
    sqrt(d)
  }

  # Distribution dispatch -- select the distribution-specific optimizer and fit.
  if (fit && dist == "multiphase") {
    # Multiphase: theta_start is optional (assembled from phase specs if NULL)
    optim_result <- .hzr_run_fit_safely(.hzr_optim_multiphase(
      time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x_fit, theta_start = theta, weights = weights,
      control = control,
      phases = phases, formula_global = formula, data = data
    ))

    fit_state$theta <- optim_result$par
    fit_state$par   <- optim_result$par   # alias for downstream use
    fit_state$objective <- optim_result$value
    fit_state$converged <- (optim_result$convergence == 0)
    fit_state$se <- .hzr_safe_se_from_vcov(optim_result$vcov)
    fit_state$vcov <- optim_result$vcov
    fit_state$rcond <- optim_result$rcond
    fit_state$pd <- optim_result$pd
    fit_state$counts <- optim_result$counts
    fit_state$message <- optim_result$message
    # Store optimizer metadata for predict/summary
    fit_state$phases <- optim_result$phases
    fit_state$covariate_counts <- optim_result$covariate_counts
    fit_state$x_list <- optim_result$x_list
    fit_state$fixed_mask <- optim_result$fixed_mask

  } else if (fit && !is.null(theta)) {
    optim_fn <- switch(
      dist,
      weibull = .hzr_optim_weibull,
      exponential = .hzr_optim_exponential,
      loglogistic = .hzr_optim_loglogistic,
      lognormal = .hzr_optim_lognormal,
      stop("Distribution '", dist, "' not yet supported for fitting.", call. = FALSE)
    )

    optim_result <- .hzr_run_fit_safely(optim_fn(
      time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x_fit, theta_start = theta, weights = weights,
      control = control
    ))

    fit_state$theta <- optim_result$par
    fit_state$par   <- optim_result$par
    fit_state$objective <- optim_result$value
    fit_state$converged <- (optim_result$convergence == 0)
    fit_state$se <- .hzr_safe_se_from_vcov(optim_result$vcov)
    fit_state$vcov <- optim_result$vcov
    fit_state$rcond <- optim_result$rcond
    fit_state$pd <- optim_result$pd
    fit_state$counts <- optim_result$counts
    fit_state$message <- optim_result$message
  }

  # Assemble the hazard S3 object.
  # $call       -- captured call for reproducibility / print
  # $spec       -- model specification (dist, control)
  # $data       -- raw data stored for default predict() / refit
  # $fit        -- optimisation results (see fit_state fields above)
  # $legacy_args -- pass-through ... args for SAS-migration parity
  # $engine     -- implementation tag ("native-r-m2")
  obj <- list(
    call = match.call(),
    spec = list(dist = dist, control = control, time_windows = time_windows,
                phases = phases),
    data = list(
      time = time,
      time_lower = time_lower,
      time_upper = time_upper,
      status = as.numeric(status),
      x = x,
      weights = weights,
      # The evaluated `data` argument as passed to hazard() (formula path; NULL
      # when called with raw vectors). This is the user's data frame, not a
      # model.frame() result. Stored so refit-based tooling such as
      # hzr_bootstrap() can re-run the original call without a caller-frame
      # lookup of the `data` symbol, which fails when that symbol has gone out
      # of scope.
      frame = data
    ),
    fit = fit_state,
    legacy_args = list(...),
    engine = "native-r-m2"
  )

  class(obj) <- "hazard"
  obj
}

#' Predict from a hazard model object
#'
#' Produces prediction outputs from a `hazard` object. Supports multiple prediction
#' types including linear predictor, hazard, survival probability, and cumulative hazard.
#'
#' @param object A `hazard` object.
#' @param newdata Optional matrix or data frame of predictors. For types requiring
#'   time (e.g., "survival", "cumulative_hazard"), newdata should include a `time`
#'   column, or time will be taken from the fitted object's data.
#' @param type Prediction type:
#'   - `"linear_predictor"`: Linear predictor eta = x*beta (not available for multiphase)
#'   - `"hazard"`: Instantaneous hazard. Single-distribution models return the
#'     hazard scale exp(eta); multiphase models return the additive hazard
#'     h(t|x) = sum_j mu_j(x) phi_j'(t) and so require time values (like
#'     `"survival"`/`"cumulative_hazard"`). `decompose` is not supported for
#'     `"hazard"`.
#'   - `"survival"`: Survival probability S(t|x) = exp(-H(t|x))
#'   - `"cumulative_hazard"`: Cumulative hazard H(t|x) at event times
#' @param decompose Logical; if `TRUE` and the model is multiphase, return a
#'   data frame with per-phase cumulative hazard contributions alongside the
#'   total.  Ignored for single-distribution models.  Default `FALSE`.
#' @param se.fit Logical; if `TRUE`, compute delta-method standard errors and
#'   confidence limits for each prediction.  The return value becomes a data
#'   frame with columns `fit`, `se.fit`, `lower`, `upper`.  Default `FALSE`.
#'   CLs are computed on the log-hazard / log-cumhaz scale and on the
#'   log(-log(survival)) scale so lower/upper stay inside the valid range
#'   of each prediction type; `linear_predictor` uses symmetric natural-scale
#'   CLs.
#'   For multiphase models, `se.fit = TRUE` combines with `decompose = TRUE`
#'   when `type = "cumulative_hazard"`: the result is a long data frame with one
#'   row per prediction time and component (`component` in `"total"` plus each
#'   phase name) and columns `fit`, `se.fit`, `lower`, `upper`. Per-phase CLs
#'   use only that phase's parameters, so they do not sum to the total CL.
#'   The combination is not available for `type = "survival"` (per-phase
#'   survival is not additive).
#' @param level Numeric confidence level in `(0, 1)`; default `0.95`.
#'   Only used when `se.fit = TRUE`.
#' @param conf.type Transform for `type = "survival"` confidence limits when
#'   `se.fit = TRUE`: `"log-log"` (default) builds them on `log(-log S)` (the
#'   `survival::survfit` standard); `"logit"` builds them on `logit(1 - S)`,
#'   reproducing SAS HAZARD's `HAZPRED` survival limits. Other types are
#'   unaffected (hazard/cumulative-hazard use a log scale that already matches
#'   HAZPRED). Only used when `se.fit = TRUE`.
#' @param ... Unused; included for S3 compatibility.
#'
#' @details
#' For Weibull models with survival or cumulative_hazard predictions:
#' - Cumulative hazard: H(t|x) = (mu*t)^nu * exp(eta)
#' - Survival: S(t|x) = exp(-H(t|x))
#'
#' Time values must be positive and finite. If newdata contains a `time` column,
#' it will be used; otherwise, the time vector from the fitted object is used.
#' For models fit with `time_windows`, predictions for `type = "linear_predictor"`
#' or `"hazard"` also require time values (via `newdata$time` or fitted-time fallback)
#' so window-specific coefficients can be selected.
#'
#' @return When `se.fit = FALSE` (default), a numeric vector of predictions.
#'   When `se.fit = TRUE`, a data frame with columns `fit`, `se.fit`, `lower`,
#'   `upper` (delta-method point estimate, standard error, and confidence
#'   limits at `level`). For multiphase `type = "cumulative_hazard"` with
#'   `decompose = TRUE`, a long data frame (`time`, `component`, `fit`,
#'   `se.fit`, `lower`, `upper`); with `decompose = TRUE` and `se.fit = FALSE`,
#'   a wide data frame of per-phase contributions.
#' @examples
#' # -- Basic predictions ------------------------------------------------
#' set.seed(1)
#' fit <- hazard(time = rexp(50, 0.3), status = rep(1L, 50),
#'               theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
#' predict(fit, type = "survival")
#' predict(fit, newdata = data.frame(time = c(1, 2, 5)),
#'         type = "cumulative_hazard")
#'
#' # -- Patient-specific survival curves ---------------------------------
#' set.seed(1001)
#' n   <- 180
#' dat <- data.frame(
#'   time   = rexp(n, rate = 0.35) + 0.05,
#'   status = rbinom(n, size = 1, prob = 0.6),
#'   age    = rnorm(n, mean = 62, sd = 11),
#'   nyha   = sample(1:4, n, replace = TRUE),
#'   shock  = rbinom(n, size = 1, prob = 0.18)
#' )
#' fit2 <- hazard(
#'   survival::Surv(time, status) ~ age + nyha + shock,
#'   data  = dat,
#'   theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
#'   dist  = "weibull", fit = TRUE
#' )
#'
#' new_patients <- data.frame(
#'   time = c(0.5, 1.5, 3.0),
#'   age  = c(50, 65, 75),
#'   nyha = c(1, 3, 4),
#'   shock = c(0, 0, 1)
#' )
#' # Compute predictions from the clean covariate frame before adding columns
#' surv   <- predict(fit2, newdata = new_patients, type = "survival")
#' cumhaz <- predict(fit2, newdata = new_patients, type = "cumulative_hazard")
#' new_patients$survival          <- surv
#' new_patients$cumulative_hazard <- cumhaz
#' new_patients
#'
#' \donttest{
#' # -- Grouped survival curves ---------------------------------------
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'
#'   t_grid <- seq(0.05, max(dat$time), length.out = 80)
#'   profiles <- data.frame(
#'     label = c("Low risk (age 50, NYHA I)",
#'               "High risk (age 75, NYHA IV)"),
#'     age   = c(50, 75),
#'     nyha  = c(1, 4),
#'     shock = c(0, 1)
#'   )
#'
#'   curve_list <- lapply(seq_len(nrow(profiles)), function(i) {
#'     nd <- data.frame(
#'       time  = t_grid,
#'       age   = profiles$age[i],
#'       nyha  = profiles$nyha[i],
#'       shock = profiles$shock[i]
#'     )
#'     nd$survival <- predict(fit2, newdata = nd, type = "survival") * 100
#'     nd$profile  <- profiles$label[i]
#'     nd
#'   })
#'   curve_df <- do.call(rbind, curve_list)
#'
#'   ggplot(curve_df, aes(time, survival, colour = profile)) +
#'     geom_line() +
#'     scale_y_continuous(limits = c(0, 100)) +
#'     labs(x = "Months after surgery",
#'          y = "Freedom from death (%)",
#'          title = "Predicted survival by risk profile",
#'          colour = NULL) +
#'     theme_minimal()
#' }
#' }
#'
#' \donttest{
#' # -- Multiphase predictions with decomposition --------------------
#' set.seed(42)
#' n   <- 200
#' dat <- data.frame(
#'   time   = rexp(n, rate = 0.25) + 0.01,
#'   status = rbinom(n, size = 1, prob = 0.65)
#' )
#' fit_mp <- hazard(
#'   survival::Surv(time, status) ~ 1,
#'   data   = dat,
#'   dist   = "multiphase",
#'   phases = list(
#'     early = hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                        fixed = "shapes"),
#'     late  = hzr_phase("cdf", t_half = 5,   nu = 1, m = 0,
#'                        fixed = "shapes")
#'   ),
#'   fit     = TRUE,
#'   control = list(n_starts = 5, maxit = 1000)
#' )
#'
#' t_grid <- seq(0.01, max(dat$time) * 0.9, length.out = 100)
#' nd     <- data.frame(time = t_grid)
#'
#' # Overall survival
#' predict(fit_mp, newdata = nd, type = "survival")
#'
#' # Per-phase decomposed cumulative hazard
#' decomp <- predict(fit_mp, newdata = nd,
#'                   type = "cumulative_hazard", decompose = TRUE)
#' head(decomp)
#' }
#'
#' @seealso
#' [hazard()] for model fitting,
#' [summary.hazard()] for model summaries,
#' [hzr_phase()] for multiphase temporal shapes.
#'
#' \code{vignette("prediction-visualization")} for detailed prediction
#' workflows including decomposed hazard plots and patient-specific curves.
#' @export
predict.hazard <- function(object, newdata = NULL,
                           type = c("hazard", "linear_predictor",
                                    "survival", "cumulative_hazard"),
                           decompose = FALSE,
                           se.fit = FALSE, level = 0.95,
                           conf.type = c("log-log", "logit"), ...) {
  type <- match.arg(type)
  # `conf.type` is validated lazily inside the se.fit survival path (it only
  # affects survival CLs); validating here would error on an otherwise-ignored
  # value (e.g. se.fit = FALSE, or a non-survival type).
  theta <- object$fit$theta
  time_windows <- object$spec$time_windows

  if (is.null(theta)) {
    stop("No coefficients ('theta') are available in 'object'.", call. = FALSE)
  }

  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be TRUE or FALSE.", call. = FALSE)
  }
  if (se.fit) {
    if (!is.numeric(level) || length(level) != 1L ||
          is.na(level) || level <= 0 || level >= 1) {
      stop("'level' must be a single number in (0, 1).", call. = FALSE)
    }
    if (decompose && object$spec$dist == "multiphase" &&
          type != "cumulative_hazard") {
      stop("'se.fit = TRUE' with 'decompose = TRUE' is only supported for ",
           "type = \"cumulative_hazard\" (per-phase survival is not additive). ",
           "Got type = \"", type, "\".", call. = FALSE)
    }
  }

  # -----------------------------------------------------------------------
  # Predictions that do NOT need time (linear_predictor, hazard)
  # -----------------------------------------------------------------------
  # These predictions work purely through the linear predictor eta = X beta.
  # Shape parameters are stripped from theta; only covariate beta's are used.
  # hazard returns exp(eta), which is the relative-hazard multiplier (not the
  # baseline conditional hazard, which would require time + distribution).
  # Multiphase `hazard` is the instantaneous additive hazard h(t) = sum_j
  # mu_j(x) * phi_j'(t); it is time-based, so it is handled in the time-based
  # branch below alongside survival / cumulative_hazard. Only the eta-based
  # types (single-dist hazard = exp(eta); linear_predictor) are handled here.
  if (type %in% c("hazard", "linear_predictor") &&
      !(type == "hazard" && object$spec$dist == "multiphase")) {
    if (object$spec$dist == "multiphase") {
      stop("Prediction type 'linear_predictor' is not supported for multiphase ",
           "models. Use 'survival', 'cumulative_hazard', or 'hazard' instead.",
           call. = FALSE)
    }
    n_pred <- NULL
    pred_time <- NULL
    if (is.null(newdata)) {
      x <- object$data$x
      pred_time <- object$data$time
    } else {
      newdata <- as.data.frame(newdata)
      n_pred <- nrow(newdata)
      if ("time" %in% names(newdata)) {
        pred_time <- newdata$time
      }
      # Remove time column if present (not needed for hazard/linear_predictor)
      newdata <- newdata[, names(newdata) != "time", drop = FALSE]
      if (ncol(newdata) > 0) {
        x <- as.matrix(newdata)
      } else {
        x <- NULL
      }
    }

    if (!is.null(time_windows)) {
      if (is.null(x)) {
        stop("Time-varying coefficients require predictors for '", type, "' predictions.", call. = FALSE)
      }
      if (is.null(pred_time)) {
        stop(
          "Time-varying coefficients require a 'time' column in newdata for '",
          type,
          "' predictions.",
          call. = FALSE
        )
      }
      # Apply the same piecewise expansion used at fit time.
      x <- .hzr_expand_time_varying_design(x = x, time = pred_time, time_windows = time_windows)
    }

    # Extract covariate coefficients by stripping shape parameters from theta.
    n_shape <- .hzr_shape_parameter_count(object$spec$dist)

    if (is.null(x)) {
      if (length(theta) <= n_shape) {
        # Univariable model: no covariates, eta = 0
        if (is.null(n_pred)) n_pred <- length(object$data$time)
        eta <- rep(0, n_pred)
      } else {
        stop("Predictors are required either in the fitted object or via 'newdata'.", call. = FALSE)
      }
    } else {
      beta <- if (length(theta) > n_shape) theta[(n_shape + 1):length(theta)] else theta
      if (ncol(x) != length(beta)) {
        stop("Number of predictor columns (", ncol(x),
             ") must match number of covariate coefficients (", length(beta), ").",
             call. = FALSE)
      }
      eta <- as.numeric(x %*% beta)
    }

    if (type == "linear_predictor") {
      if (!se.fit) return(eta)
      diff_fn <- function(th) {
        if (is.null(x)) return(rep(0, length(eta)))
        n_shape <- .hzr_shape_parameter_count(object$spec$dist)
        beta_cand <- if (length(th) > n_shape) th[(n_shape + 1):length(th)] else th
        as.numeric(x %*% beta_cand)
      }
      return(.hzr_predict_with_se(object = object, type = "linear_predictor",
                                    time = NULL, x = x, level = level,
                                    diff_fn = diff_fn))
    }
    if (!se.fit) return(exp(eta))
    diff_fn <- function(th) {
      if (is.null(x)) return(rep(1, length(eta)))
      n_shape <- .hzr_shape_parameter_count(object$spec$dist)
      beta_cand <- if (length(th) > n_shape) th[(n_shape + 1):length(th)] else th
      exp(as.numeric(x %*% beta_cand))
    }
    return(.hzr_predict_with_se(object = object, type = "hazard",
                                  time = NULL, x = x, level = level,
                                  diff_fn = diff_fn))
  }

  # -----------------------------------------------------------------------
  # Predictions that DO need time (survival, cumulative_hazard)
  # -----------------------------------------------------------------------
  # Each distribution computes H(t|x) in its own way; all then share the
  # same final return statements:
  #   cumulative_hazard -> return(cumhaz)
  #   survival          -> return(exp(-cumhaz))
  #
  # EXCEPTION: log-normal uses a direct Phi(-z) formula for survival and
  # returns early inside its branch, bypassing the shared return below.
  # This is because the log-normal is an AFT model where H(t) = -log Phi(-z)
  # is numerically better computed directly from the normal CDF rather than
  # via exp(-H).
  # `hazard` is included here only for multiphase (single-dist hazard returned
  # above via exp(eta)); it is the instantaneous additive hazard h(t).
  if (type %in% c("survival", "cumulative_hazard", "hazard")) {
    supported <- c("weibull", "exponential", "loglogistic", "lognormal", "multiphase")
    if (!object$spec$dist %in% supported) {
      stop("Prediction type '", type, "' is only supported for ",
           paste(supported, collapse = ", "), " models.", call. = FALSE)
    }

    # --- Multiphase prediction (early return) ---------------------------------
    if (object$spec$dist == "multiphase") {
      # Extract time from newdata or fitted data
      if (!is.null(newdata)) {
        newdata <- as.data.frame(newdata)
        if (!"time" %in% names(newdata)) {
          stop("'newdata' must contain a 'time' column for '", type, "' predictions.",
               call. = FALSE)
        }
        pred_time <- newdata$time
      } else {
        pred_time <- object$data$time
      }

      if (!is.numeric(pred_time) || any(!is.finite(pred_time)) || any(pred_time < 0)) {
        stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
      }

      # Recover phase metadata from fit
      phases <- object$fit$phases
      if (is.null(phases)) phases <- object$spec$phases
      cov_counts <- object$fit$covariate_counts
      x_list <- object$fit$x_list

      # If newdata has covariates, rebuild per-phase design matrices
      if (!is.null(newdata)) {
        nd_covs <- newdata[, names(newdata) != "time", drop = FALSE]
        for (nm in names(phases)) {
          ph <- phases[[nm]]
          if (!is.null(ph$formula) && ncol(nd_covs) > 0) {
            x_list[[nm]] <- stats::model.matrix(ph$formula, data = newdata)[, -1L, drop = FALSE]
          } else if (cov_counts[[nm]] > 0 && ncol(nd_covs) > 0) {
            x_list[[nm]] <- as.matrix(nd_covs)
          } else {
            x_list[[nm]] <- NULL
          }
        }
      }

      # Instantaneous additive hazard is not phase-decomposable through this
      # path (the internal hazard evaluator returns the total only).
      if (type == "hazard" && decompose) {
        stop("decompose = TRUE is not supported for type = 'hazard'.",
             call. = FALSE)
      }

      if (se.fit) {
        if (decompose) {
          return(.hzr_predict_with_se_decomposed( # nolint: object_usage_linter.
            object = object, time = pred_time, x_list = x_list,
            cov_counts = cov_counts, phases = phases, level = level
          ))
        }
        diff_fn <- if (type == "hazard") {
          function(th) {
            .hzr_multiphase_hazard(pred_time, th, phases, cov_counts, x_list)
          }
        } else {
          function(th) {
            .hzr_multiphase_cumhaz(pred_time, th, phases, cov_counts, x_list)
          }
        }
        return(.hzr_predict_with_se(
          object = object, type = type, time = pred_time,
          x_list = x_list, cov_counts = cov_counts, phases = phases,
          level = level, diff_fn = diff_fn, conf_type = conf.type
        ))
      }

      if (type == "hazard") {
        return(.hzr_multiphase_hazard(pred_time, theta, phases, cov_counts,
                                      x_list))
      }

      result <- .hzr_multiphase_cumhaz(
        pred_time, theta, phases, cov_counts, x_list,
        per_phase = decompose
      )

      if (decompose) {
        # Return a data frame with time, total, and per-phase columns
        out <- data.frame(time = pred_time, total = result$total)
        for (nm in names(phases)) {
          out[[nm]] <- result[[nm]]
        }
        if (type == "survival") {
          out$total <- exp(-out$total)
          for (nm in names(phases)) {
            # Per-phase survival contribution is not additive; provide cumhaz
            # columns as-is and only transform total.
          }
        }
        return(out)
      }

      cumhaz <- result
      if (type == "cumulative_hazard") return(cumhaz)
      return(exp(-cumhaz))
    }

    # --- Standard single-distribution prediction ------------------------------
    # Extract time from newdata or use fitted data
    if (!is.null(newdata)) {
      newdata <- as.data.frame(newdata)
      if ("time" %in% names(newdata)) {
        time <- newdata$time
        x <- as.matrix(newdata[, names(newdata) != "time", drop = FALSE])
      } else {
        stop("'newdata' must contain a 'time' column for '", type, "' predictions.", call. = FALSE)
      }
    } else {
      time <- object$data$time
      x <- object$data$x
    }

    if (!is.null(time_windows) && !is.null(x) && ncol(x) > 0) {
      # Keep prediction design matrix consistent with training-time expansion.
      x <- .hzr_expand_time_varying_design(x = x, time = time, time_windows = time_windows)
    }

    if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
      stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
    }

    # Extract shape parameters and covariate coefficients
    n_shape <- .hzr_shape_parameter_count(object$spec$dist)

    if (!is.null(x) && ncol(x) > 0) {
      if (length(theta) < ncol(x) + n_shape) {
        stop("Number of parameters insufficient for predictor columns.", call. = FALSE)
      }
      beta_coef <- theta[(n_shape + 1):length(theta)]
      eta <- as.numeric(x %*% beta_coef)
    } else {
      eta <- rep(0, length(time))
    }

    # Dispatch cumulative hazard computation by distribution.
    # Log-normal is an AFT model, so H is derived from the standardised
    # residual z = (log t - mu - x beta) / sigma.  For se.fit we define a
    # closure `cumhaz_of` that computes H for any candidate theta -- this
    # is the delta-method target for both "cumulative_hazard" and
    # "survival" predictions.
    dist_lbl <- object$spec$dist
    has_cov <- !is.null(x) && ncol(x) > 0

    # Preserve the pre-0.9.8 stop() behaviour on an ill-conditioned MLE.
    # The closures below return NA on negative shape parameters so numeric
    # jacobian perturbations stay robust, but we want a clean error at the
    # point estimate itself.
    if (dist_lbl == "weibull" && (theta[1] <= 0 || theta[2] <= 0)) {
      stop("Weibull shape parameters (mu, nu) must be positive.", call. = FALSE)
    }

    cumhaz_of <- if (dist_lbl == "weibull") {
      function(th) {
        if (th[1] <= 0 || th[2] <= 0) return(rep(NA_real_, length(time)))
        beta_cand <- if (length(th) > 2) th[3:length(th)] else numeric(0)
        eta_cand <- if (has_cov) as.numeric(x %*% beta_cand) else rep(0, length(time))
        (th[1] * time) ^ th[2] * exp(eta_cand)
      }
    } else if (dist_lbl == "exponential") {
      function(th) {
        beta_cand <- if (length(th) > 1) th[2:length(th)] else numeric(0)
        eta_cand <- if (has_cov) as.numeric(x %*% beta_cand) else rep(0, length(time))
        exp(th[1]) * time * exp(eta_cand)
      }
    } else if (dist_lbl == "loglogistic") {
      function(th) {
        beta_cand <- if (length(th) > 2) th[3:length(th)] else numeric(0)
        eta_cand <- if (has_cov) as.numeric(x %*% beta_cand) else rep(0, length(time))
        log(1 + exp(th[1]) * (time ^ exp(th[2])) * exp(eta_cand))
      }
    } else if (dist_lbl == "lognormal") {
      function(th) {
        beta_cand <- if (length(th) > 2) th[3:length(th)] else numeric(0)
        # AFT: covariates shift the location.
        eta_aft <- if (has_cov) th[1] + as.numeric(x %*% beta_cand) else rep(th[1], length(time))
        z <- (log(time) - eta_aft) / exp(th[2])
        -pnorm(-z, log.p = TRUE)
      }
    } else {
      stop("Unknown distribution '", dist_lbl, "'.", call. = FALSE)
    }

    if (se.fit) {
      return(.hzr_predict_with_se(
        object = object, type = type, time = time, x = x,
        level = level, diff_fn = cumhaz_of, conf_type = conf.type
      ))
    }

    cumhaz <- cumhaz_of(theta)
    if (type == "cumulative_hazard") return(cumhaz)
    return(exp(-cumhaz))
  }

  stop("Unknown prediction type: '", type, "'.", call. = FALSE)
}


#' Print method for fitted hazard models
#'
#' Compact one-block summary of a fitted `hazard` object: sample size,
#' number of predictors, distribution, theta vector, and log-likelihood.
#' S3 dispatch only -- users call `print(fit)` rather than invoking this
#' directly.
#'
#' @param x A `hazard` object returned by [hazard()].
#' @param ... Additional arguments (ignored).
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.hazard <- function(x, ...) {
  n <- length(x$data$time)
  p <- if (is.null(x$data$x)) 0L else ncol(x$data$x)
  cat("hazard object\n")
  cat("  observations:", n, "\n")
  cat("  predictors:  ", p, "\n")
  cat("  dist:        ", x$spec$dist, "\n")

  if (x$spec$dist == "multiphase" && !is.null(x$spec$phases)) {
    ph <- x$spec$phases
    cat("  phases:      ", length(ph),
        " (", paste(names(ph), collapse = ", "), ")\n", sep = "")
  }

  cat("  engine:      ", x$engine, "\n")
  if (!anyNA(x$fit$objective)) {
    cat("  log-lik:     ", format(x$fit$objective, digits = 6), "\n")
    cat("  converged:   ", x$fit$converged, "\n")
  }
  invisible(x)
}

#' Summarize a hazard model
#'
#' Returns a compact summary of a `hazard` object, including model metadata,
#' fit diagnostics, and coefficient-level statistics when available.
#'
#' @param object A `hazard` object.
#' @param ... Unused; for S3 compatibility.
#' @return An object of class `summary.hazard`.
#' @examples
#' # -- Single-phase Weibull summary ------------------------------------
#' fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
#'               theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
#' summary(fit)
#'
#' \donttest{
#' # -- Multiphase model summary ----------------------------------------
#' set.seed(42)
#' n   <- 200
#' dat <- data.frame(
#'   time   = rexp(n, rate = 0.25) + 0.01,
#'   status = rbinom(n, size = 1, prob = 0.65)
#' )
#' fit_mp <- hazard(
#'   survival::Surv(time, status) ~ 1,
#'   data   = dat,
#'   dist   = "multiphase",
#'   phases = list(
#'     early = hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                        fixed = "shapes"),
#'     late  = hzr_phase("cdf", t_half = 5,   nu = 1, m = 0,
#'                        fixed = "shapes")
#'   ),
#'   fit     = TRUE,
#'   control = list(n_starts = 5, maxit = 1000)
#' )
#' summary(fit_mp)
#' }
#'
#' @seealso
#' [hazard()] for model fitting, [predict.hazard()] for predictions.
#'
#' \code{vignette("fitting-hazard-models")} for fitting workflows,
#' \code{vignette("inference-diagnostics")} for bootstrap CIs and diagnostics.
#' @export
summary.hazard <- function(object, ...) {
  n <- length(object$data$time)
  p <- if (is.null(object$data$x)) 0L else ncol(object$data$x)
  theta <- object$fit$theta
  vcov_mat <- object$fit$vcov

  coef_table <- NULL
  if (!is.null(theta)) {
    # For multiphase, use the theta names directly (they're already informative)
    if (object$spec$dist == "multiphase") {
      coef_names <- names(theta)
      if (is.null(coef_names)) coef_names <- paste0("param_", seq_along(theta))
    } else {
      coef_names <- .hzr_parameter_names(theta = theta, dist = object$spec$dist, p = p)
    }

    std_error <- rep(NA_real_, length(theta))
    z_stat <- rep(NA_real_, length(theta))
    p_value <- rep(NA_real_, length(theta))

    if (!is.null(vcov_mat) && is.matrix(vcov_mat)) {
      d <- diag(vcov_mat)
      std_error <- sqrt(d)            # NA for fixed params, finite for free
      valid <- is.finite(std_error) & std_error > 0
      z_stat[valid] <- theta[valid] / std_error[valid]
      p_value[valid] <- 2 * pnorm(-abs(z_stat[valid]))
    }

    coef_table <- data.frame(
      estimate = unname(theta),
      std_error = std_error,
      z_stat = z_stat,
      p_value = p_value,
      row.names = coef_names,
      check.names = FALSE
    )
  }

  out <- list(
    call = object$call,
    n = n,
    p = p,
    dist = object$spec$dist,
    engine = object$engine,
    converged = object$fit$converged,
    log_lik = object$fit$objective,
    counts = object$fit$counts,
    message = object$fit$message,
    coefficients = coef_table,
    has_vcov = !is.null(vcov_mat) && is.matrix(vcov_mat),
    rcond = object$fit$rcond,
    pd = object$fit$pd,
    phases = object$spec$phases
  )

  class(out) <- "summary.hazard"
  out
}

#' Print method for hazard summary objects
#'
#' Formatted console display of [summary.hazard()] output: distribution,
#' phase list (for multiphase), coefficient table with standard errors,
#' and log-likelihood.  When the post-fit Hessian is ill-conditioned or not
#' positive-definite, a note warns that the standard errors may be unreliable;
#' when the Hessian could not be inverted at all, a note reports that standard
#' errors are unavailable.  S3 dispatch only -- users call `print(summary(fit))`
#' rather than invoking this directly.
#'
#' @param x A `summary.hazard` object returned by [summary.hazard()].
#' @param ... Additional arguments (ignored).
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.summary.hazard <- function(x, ...) {
  if (x$dist == "multiphase" && !is.null(x$phases)) {
    cat("Multiphase hazard model (", length(x$phases), " phases)\n", sep = "")
  } else {
    cat("hazard model summary\n")
  }
  cat("  observations:", x$n, "\n")
  cat("  predictors:  ", x$p, "\n")
  cat("  dist:        ", x$dist, "\n")

  if (!is.null(x$phases)) {
    for (i in seq_along(x$phases)) {
      nm <- names(x$phases)[i]
      ph <- x$phases[[i]]
      label <- switch(ph$type,
        cdf      = paste0("cdf (", nm, " risk)"),
        hazard   = "hazard (late risk)",
        constant = "constant (flat rate)",
        g3       = "g3 (late risk)")
      cat("  phase ", i, ":      ", nm, " - ", label, "\n", sep = "")
    }
  }

  cat("  engine:      ", x$engine, "\n")

  if (!is.null(x$converged) && !is.na(x$converged)) {
    cat("  converged:   ", x$converged, "\n")
  }
  if (!is.null(x$log_lik) && !is.na(x$log_lik)) {
    cat("  log-lik:     ", format(x$log_lik, digits = 6), "\n")
  }
  if (!is.null(x$rcond) && !is.na(x$rcond) && x$rcond < .hzr_rcond_tol) {
    cat("  Note: Hessian ill-conditioned (rcond = ",
        format(x$rcond, digits = 3),
        "); standard errors may be unreliable.\n", sep = "")
  }
  if (!is.null(x$pd) && !is.na(x$pd) && !isTRUE(x$pd)) {
    cat("  Note: Hessian not positive-definite at the optimum; ",
        "standard errors may be unreliable.\n", sep = "")
  }
  if (!is.null(x$converged) && !is.na(x$converged) && isFALSE(x$has_vcov)) {
    cat("  Note: standard errors unavailable; the Hessian could not be ",
        "inverted.\n", sep = "")
  }
  if (!is.null(x$counts)) {
    fn_count <- x$counts[["function"]] %||% x$counts[["fn"]] %||% NA_integer_
    gr_count <- x$counts[["gradient"]] %||% x$counts[["gr"]] %||% NA_integer_
    if (!is.na(fn_count) || !is.na(gr_count)) {
      cat("  evaluations: ", "fn=", fn_count, ", gr=", gr_count, "\n", sep = "")
    }
  }
  if (!is.null(x$message) && nzchar(x$message)) {
    cat("  message:     ", x$message, "\n")
  }

  if (!is.null(x$coefficients)) {
    if (!is.null(x$phases)) {
      # Group coefficients by phase for readable output
      cat("\nCoefficients (internal scale):\n")
      for (nm in names(x$phases)) {
        prefix <- paste0("^", nm, "\\.")
        rows <- grep(prefix, rownames(x$coefficients))
        if (length(rows) > 0) {
          ph <- x$phases[[nm]]
          label <- switch(ph$type,
            cdf = "cdf", hazard = "hazard", constant = "constant",
            g3 = "g3")
          cat("\n  Phase: ", nm, " (", label, ")\n", sep = "")
          sub_table <- x$coefficients[rows, , drop = FALSE]
          # Strip phase prefix from row names for cleaner display
          rownames(sub_table) <- sub(prefix, "  ", rownames(sub_table))
          print(sub_table)
        }
      }
    } else {
      cat("\nCoefficients:\n")
      print(x$coefficients)
    }
  }

  invisible(x)
}

#' Extract coefficients from hazard model
#'
#' @param object A `hazard` object.
#' @param ... Unused; for S3 compatibility.
#' @examples
#' fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
#'               theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
#' coef(fit)
#' @return A named numeric vector of fitted parameter estimates, or \code{NULL}
#'   if the model has not been fitted (\code{fit = FALSE}).
#' @export
coef.hazard <- function(object, ...) {
  if (is.null(object$fit$theta)) {
    return(NULL)
  }
  object$fit$theta
}

#' Extract variance-covariance matrix from hazard model
#'
#' Returns the estimated variance-covariance matrix of the fitted coefficients.
#'
#' @param object A `hazard` object.
#' @param ... Unused; for S3 compatibility.
#' @examples
#' fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
#'               theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
#' vcov(fit)
#' @return A numeric matrix containing the estimated variance-covariance matrix
#'   of the fitted coefficients, with rows and columns named by the coefficient
#'   labels (phase-prefixed for multiphase models, e.g. \code{early.x}). Rows
#'   and columns for parameters held fixed (e.g. fixed shape parameters) are
#'   \code{NA} because they carry no variance; the finite free-parameter block
#'   is still usable. For Conservation-of-Events fits the conserved phase
#'   \code{log_mu} \emph{normally} carries a variance: it is removed from the
#'   optimizer search but the vcov is the full-information matrix at the optimum
#'   (the CoE solution is the unconstrained MLE). That recomputation requires
#'   \pkg{numDeriv} and an invertible Hessian; if either is unavailable the fit
#'   emits a warning and the conserved \code{log_mu} stays \code{NA} (the rest
#'   of the matrix is unaffected). Returns a scalar \code{NA} only when the
#'   model has not been fitted or no covariance matrix is available.
#' @export
vcov.hazard <- function(object, ...) {
  v <- object$fit$vcov
  if (is.null(v) || !is.matrix(v)) {
    return(NA)
  }
  # Label rows/cols with the coefficient names so callers can align the matrix
  # by name. This matters for multiphase models where the same covariate can
  # enter more than one phase (e.g. early.x vs constant.x): without names the
  # two coefficients are indistinguishable. NA variance rows are retained
  # rather than collapsing the whole matrix to a scalar NA -- a multiphase fit
  # legitimately has NA rows for parameters held fixed (e.g. early shapes),
  # which carry no Hessian-based variance. (The CoE-conserved log_mu is
  # normally NOT NA -- it gets the full-information variance -- but stays NA if
  # that recomputation was unavailable; the fit warns in that case.)
  # The finite free-parameter block is still usable.
  theta <- object$fit$theta
  nm <- names(theta)
  if (is.null(nm) || !all(nzchar(nm))) {
    p <- if (is.null(object$data$x)) 0L else ncol(object$data$x)
    nm <- .hzr_parameter_names(theta = theta, dist = object$spec$dist, p = p)
  }
  if (!is.null(nm) && length(nm) == nrow(v)) {
    dimnames(v) <- list(nm, nm)
  }
  v
}

.hzr_parameter_names <- function(theta, dist, p) {
  theta_names <- names(theta)
  if (!is.null(theta_names) && all(nzchar(theta_names))) {
    return(theta_names)
  }

  if (!is.null(p) && p > 0L && length(theta) == p) {
    return(paste0("beta", seq_along(theta)))
  }

  base_names <- switch(
    dist,
    weibull = c("mu", "nu"),
    exponential = c("log_lambda"),
    loglogistic = c("log_alpha", "log_beta"),
    lognormal = c("mu", "log_sigma"),
    character()
  )

  n_theta <- length(theta)
  n_base <- min(length(base_names), n_theta)
  out <- character(n_theta)

  if (n_base > 0) {
    out[seq_len(n_base)] <- base_names[seq_len(n_base)]
  }
  if (n_theta > n_base) {
    out[seq.int(n_base + 1L, n_theta)] <- paste0("beta", seq_len(n_theta - n_base))
  }
  out
}

#' Coerce x to a validated numeric design matrix
#'
#' Accepts data.frame or numeric matrix input; validates dimensions and
#' finiteness.  Called by hazard() to normalise the x argument before any
#' downstream use.
#'
#' @param x A numeric matrix or data frame with numeric columns.
#' @param n Expected row count (length of time/status vectors); checked if non-NULL.
#' @return A numeric matrix with column names preserved (or added if absent).
#' @keywords internal
#'
.hzr_as_design_matrix <- function(x, n = NULL) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Predictor input must be a numeric matrix or coercible data frame.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop("Predictor matrix contains non-finite values.", call. = FALSE)
  }

  if (!is.null(n) && nrow(x) != n) {
    stop("Predictor rows must match the length of 'time'.", call. = FALSE)
  }

  x
}

#' Expand predictors for piecewise time-varying coefficients
#'
#' Builds a block design matrix with one predictor block per time window. For each
#' observation, only the block corresponding to that row's time window is active;
#' all other blocks are zero.
#'
#' Example with two predictors and one cut point:
#' - input columns: `x1`, `x2`
#' - windows: `(-Inf, c1]`, `(c1, Inf)`
#' - output columns: `x1_w1`, `x2_w1`, `x1_w2`, `x2_w2`
#'
#' @param x Numeric predictor matrix/data.frame.
#' @param time Numeric time vector used to assign rows to windows.
#' @param time_windows Numeric vector of window cut points.
#' @return Expanded numeric matrix with window-specific column names.
#' @noRd
.hzr_expand_time_varying_design <- function(x, time, time_windows) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Time-varying design requires numeric matrix/data.frame predictors.", call. = FALSE)
  }

  if (!is.numeric(time) || length(time) != nrow(x) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must be finite, non-negative, and match predictor rows for time-varying expansion.", call. = FALSE)
  }

  if (!is.numeric(time_windows) || length(time_windows) < 1) {
    stop("'time_windows' must contain at least one cut point.", call. = FALSE)
  }

  time_windows <- sort(unique(time_windows))
  # Bin index k means observation time falls in window k.
  bins <- cut(time, breaks = c(-Inf, time_windows, Inf), right = TRUE, include.lowest = TRUE, labels = FALSE)
  n_bins <- length(time_windows) + 1L

  base_names <- colnames(x)
  if (is.null(base_names)) {
    base_names <- paste0("x", seq_len(ncol(x)))
  }

  out <- matrix(0, nrow = nrow(x), ncol = ncol(x) * n_bins)
  col_ptr <- 1L
  out_names <- character(ncol(out))

  for (k in seq_len(n_bins)) {
    idx <- bins == k
    # Populate only the active window block for rows in this bin.
    block <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    if (any(idx)) {
      block[idx, ] <- x[idx, , drop = FALSE]
    }
    cols <- col_ptr:(col_ptr + ncol(x) - 1L)
    out[, cols] <- block
    out_names[cols] <- paste0(base_names, "_w", k)
    col_ptr <- col_ptr + ncol(x)
  }

  colnames(out) <- out_names
  out
}
