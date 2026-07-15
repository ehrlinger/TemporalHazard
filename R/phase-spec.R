# phase-spec.R -- Phase specification for multiphase hazard models
#
# PURPOSE
# -------
# Defines the hzr_phase() constructor that specifies one phase in an
# N-phase additive cumulative hazard model:
#
#   H(t | x) = sum_j  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)
#
# Each phase is one term in the sum.  The constructor stores:
#   - type:    "cdf" (early risk), "hazard" (late risk), or "constant"
#   - t_half, nu, m: starting values for the decompos shape parameters
#   - formula: optional one-sided formula for phase-specific covariates
#
# The helpers extract metadata needed during likelihood construction:
#   - .hzr_phase_n_shape():     number of shape parameters (3 or 0)
#   - .hzr_phase_theta_names(): named labels for the parameter sub-vector
#
# SAS/C BRIDGE
# ------------
# SAS users specify Early/Constant/Late phases via DIST keywords.  The
# hzr_phase() constructor maps directly:
#   SAS Early (G1) -> hzr_phase("cdf",      t_half, nu, m)
#   SAS Const (G2) -> hzr_phase("constant")
#   SAS Late  (G3) -> hzr_phase("g3",       tau, gamma, alpha, eta)
#
# The "cdf" and "hazard" types use the G1 decomposition (bounded [0,1]).
# The "g3" type uses the G3 decomposition (unbounded power-law), matching
# the C/SAS HAZARD implementation for late-phase rising hazards.
# The "hazard" type (-log(1-G(t))) is available for alternative models.

# ============================================================================
# Constructor
# ============================================================================

#' Specify a single hazard phase
#'
#' Creates an `hzr_phase` object describing one term in a multiphase additive
#' cumulative hazard model.  Pass a list of these to the `phases` argument of
#' [hazard()] when `dist = "multiphase"`.
#'
#' @section Role in the multiphase model:
#'
#' Each phase is one term \eqn{j} in the additive cumulative hazard
#'
#' \deqn{H(t \mid \mathbf{x}) = \sum_{j=1}^{J} \mu_j(\mathbf{x}) \,
#'       \Phi_j(t)}
#'
#' where \eqn{\mu_j(\mathbf{x}) = \exp(\alpha_j + \mathbf{x}_j^\top
#' \beta_j)} is the phase-specific log-linear scale and
#' \eqn{\Phi_j(t)} is the temporal shape selected by `type` (below).  The
#' `t_half`/`nu`/`m` (or g3 `tau`/`gamma`/`alpha`/`eta`) arguments set the
#' starting values for that shape; `formula` attaches the covariates
#' \eqn{\mathbf{x}_j} that enter \eqn{\mu_j}.
#'
#' @section Phase types:
#'
#' The `type` argument chooses the temporal shape \eqn{\Phi_j(t)} for the phase.
#' Each captures a qualitatively different pattern of risk over time; a typical
#' clinical model combines an *early*, a *constant*, and a *late* phase so that
#' the total hazard can fall, level off, and rise again.
#'
#' \describe{
#'   \item{`"cdf"` --- early, resolving risk}{Named for the **c**umulative
#'     **d**istribution **f**unction: the phase contributes \eqn{\Phi(t) = G(t)},
#'     the bounded CDF of the temporal decomposition (\eqn{0} at \eqn{t = 0},
#'     rising to a ceiling of \eqn{1}).  Because it saturates, the *hazard* it
#'     adds, \eqn{\mu\,g(t)}, peaks early and then decays toward zero --- the
#'     signature of a one-time insult that patients either succumb to or survive
#'     past, e.g. peri-operative mortality.  Shape set by `t_half`, `nu`, `m`.
#'     SAS/C equivalent: the Early (G1) phase.}
#'   \item{`"hazard"` --- accumulating aging risk (G1 family)}{Named because the
#'     phase contributes a **cumulative hazard** built from the same G1 family:
#'     \eqn{\Phi(t) = -\log(1 - G(t))}, which is unbounded and monotone
#'     increasing.  Its hazard \eqn{\mu\,h(t)} rises without leveling off, so it
#'     models risk that grows as subjects age.  This is an alternative late-risk
#'     form derived from G1; for the original SAS/C late phase prefer `"g3"`.
#'     Shape set by `t_half`, `nu`, `m`.}
#'   \item{`"g3"` --- late, rising risk (original C/SAS late phase)}{Named for
#'     the **G3** (third) decomposition family used by the original HAZARD
#'     program for the late phase.  It contributes \eqn{\Phi(t) = G_3(t)} from
#'     [hzr_decompos_g3()], an unbounded intensity with its own four-parameter
#'     shape (`tau`, `gamma`, `alpha`, `eta`) that is more flexible than the
#'     G1-derived `"hazard"` form for capturing accelerating late mortality
#'     (e.g. structural valve deterioration years after surgery).  Use this when
#'     reproducing classic three-phase HAZARD models.  SAS/C equivalent: the
#'     Late (G3) phase.}
#'   \item{`"constant"` --- flat background rate}{A time-invariant hazard:
#'     \eqn{\Phi(t) = t}, so the added hazard \eqn{\mu} is constant (the
#'     exponential model).  It represents the steady, ongoing risk present at all
#'     follow-up times, independent of how long ago the time origin was.  Takes
#'     no shape parameters --- only its scale \eqn{\mu} (and any covariates) is
#'     estimated.  SAS/C equivalent: the Constant (G2) phase.}
#' }
#'
#' The shape derivative \eqn{\varphi_j = d\Phi_j/dt} (which forms the
#' instantaneous hazard contribution \eqn{\mu_j\,\varphi_j(t)}) is \eqn{g(t)} for
#' `"cdf"`, \eqn{h(t)} for `"hazard"`, \eqn{g_3(t)} for `"g3"`, and \eqn{1} for
#' `"constant"`.
#'
#' @param type Character; the phase's temporal shape --- one of `"cdf"`
#'   (early resolving risk), `"hazard"` (accumulating G1 aging risk), `"g3"`
#'   (late rising risk, the original C/SAS late phase), or `"constant"` (flat
#'   background rate).  See the **Phase types** section for what each means and
#'   when to use it.
#' @param t_half Positive scalar; initial half-life (time at which
#'   \eqn{G(t_{1/2}) = 0.5}).  Used for `"cdf"` and `"hazard"` phases.
#'   SAS early: `THALF`/`RHO`.
#' @param nu Numeric scalar; initial time exponent.  Used for `"cdf"` and
#'   `"hazard"` phases.  SAS early: `NU`.
#' @param m Numeric scalar; initial shape exponent.  Used for `"cdf"` and
#'   `"hazard"` phases.  SAS early: `M`.
#' @param tau Positive scalar; scale parameter for `"g3"` phases.
#'   SAS late: `TAU`.
#' @param gamma Positive scalar; time exponent for `"g3"` phases.
#'   SAS late: `GAMMA`.
#' @param alpha Non-negative scalar; shape parameter for `"g3"` phases.
#'   When `alpha > 0`, the generic G3 formula is used; `alpha = 0` gives the
#'   exponential limiting case.  SAS late: `ALPHA`.
#' @param eta Positive scalar; outer exponent for `"g3"` phases.
#'   SAS late: `ETA`.
#' @param formula Optional one-sided formula (e.g. `~ age + nyha`) for
#'   phase-specific covariates.  When `NULL` (default), the phase inherits
#'   the global formula from [hazard()].
#' @param fixed Character vector naming shape parameters to hold fixed during
#'   optimization.  Valid names for `"cdf"`/`"hazard"`: `"t_half"`, `"nu"`,
#'   `"m"`, or `"shapes"` (shorthand for all three).  Valid names for `"g3"`:
#'   `"tau"`, `"gamma"`, `"alpha"`, `"eta"`, or `"shapes"` (shorthand for all
#'   four).  Fixed parameters are held at their starting values; only `mu`
#'   (and covariates) are estimated.  Ignored for `"constant"` phases.
#'   This mirrors the SAS/C HAZARD workflow where shapes are typically fixed
#'   and only scale parameters are estimated.
#'
#' @return An S3 object of class `"hzr_phase"` with elements:
#' \describe{
#'   \item{type}{Phase type string.}
#'   \item{t_half}{Initial half-life (cdf/hazard phases).}
#'   \item{nu}{Initial time exponent (cdf/hazard phases).}
#'   \item{m}{Initial shape exponent (cdf/hazard phases).}
#'   \item{tau}{Scale parameter (g3 phases).}
#'   \item{gamma}{Time exponent (g3 phases).}
#'   \item{alpha}{Shape parameter (g3 phases).}
#'   \item{eta}{Outer exponent (g3 phases).}
#'   \item{formula}{Phase-specific formula or `NULL`.}
#'   \item{fixed}{Character vector of fixed parameter names (may be empty).}
#' }
#'
#' @examples
#' # Classic 3-phase Blackstone pattern
#' early <- hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0)
#' const <- hzr_phase("constant")
#' late  <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1)
#'
#' # Fix all shapes (C/SAS-style: only estimate mu)
#' early_fixed <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                           fixed = "shapes")
#' late_fixed  <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
#'                           fixed = "shapes")
#'
#' # Fix only some parameters
#' early_partial <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                             fixed = c("nu", "m"))
#'
#' # Phase with specific covariates
#' early_cov <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                         formula = ~ age + shock)
#'
#' # Use in hazard():
#' # hazard(Surv(time, status) ~ age, data = dat,
#' #        dist = "multiphase",
#' #        phases = list(early = early, constant = const, late = late))
#'
#' @seealso [hazard()] for fitting multiphase models,
#'   [hzr_decompos()] for the underlying parametric family,
#'   [hzr_phase_cumhaz()] and [hzr_phase_hazard()] for computing
#'   \eqn{\Phi(t)} and \eqn{\phi(t)} from these specifications.
#'
#' \code{vignette("fitting-hazard-models")} for multiphase fitting examples,
#' \code{vignette("mf-mathematical-foundations")} for the mathematical framework.
#'
#' @export
hzr_phase <- function(type = c("cdf", "hazard", "constant", "g3"),
                      t_half = 1, nu = 1, m = 0,
                      tau = 1, gamma = 1, alpha = 1, eta = 1,
                      formula = NULL,
                      fixed = character(0)) {

  type <- match.arg(type)

  # --- Validate shape parameters based on type ------------------------------
  if (type == "g3") {
    # G3 late-phase decomposition: 4 parameters (tau, gamma, alpha, eta)
    stopifnot(
      "tau must be a positive scalar" =
        is.numeric(tau) && length(tau) == 1L && tau > 0,
      "gamma must be a positive scalar" =
        is.numeric(gamma) && length(gamma) == 1L && gamma > 0,
      "alpha must be a non-negative scalar" =
        is.numeric(alpha) && length(alpha) == 1L && alpha >= 0 && is.finite(alpha),
      "eta must be a positive scalar" =
        is.numeric(eta) && length(eta) == 1L && eta > 0
    )
  } else if (type != "constant") {
    # G1 early-phase decomposition: 3 parameters (t_half, nu, m)
    stopifnot(
      "t_half must be a positive scalar" =
        is.numeric(t_half) && length(t_half) == 1L &&
        is.finite(t_half) && t_half > 0,
      "nu must be a finite numeric scalar" =
        is.numeric(nu) && length(nu) == 1L && is.finite(nu),
      "m must be a finite numeric scalar" =
        is.numeric(m) && length(m) == 1L && is.finite(m)
    )

    if (m < 0 && nu < 0) {
      stop("Decomposition undefined when both m < 0 and nu < 0 (m = ",
           m, ", nu = ", nu, ").  See ?hzr_decompos for valid cases.",
           call. = FALSE)
    }
  }

  # --- Validate formula if supplied -----------------------------------------
  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) {
      stop("formula must be a one-sided formula (e.g. ~ age + nyha) or NULL.",
           call. = FALSE)
    }
    # Must be one-sided (no LHS)
    if (length(formula) == 3L) {
      stop("Phase formula must be one-sided (e.g. ~ age + nyha), ",
           "not two-sided (response ~ predictors).", call. = FALSE)
    }
  }

  # --- Validate and normalize fixed parameters ------------------------------
  if (length(fixed) > 0 && type == "constant") {
    warning("fixed is ignored for constant phases.", call. = FALSE)
    fixed <- character(0)
  }
  if (length(fixed) > 0) {
    if (type == "g3") {
      # G3: expand "shapes" to tau, gamma, alpha, eta
      if ("shapes" %in% fixed) {
        fixed <- union(setdiff(fixed, "shapes"), c("tau", "gamma", "alpha", "eta"))
      }
      valid_fixed <- c("tau", "gamma", "alpha", "eta")
    } else {
      # G1: expand "shapes" to t_half, nu, m
      if ("shapes" %in% fixed) {
        fixed <- union(setdiff(fixed, "shapes"), c("t_half", "nu", "m"))
      }
      valid_fixed <- c("t_half", "nu", "m")
    }
    bad <- setdiff(fixed, valid_fixed)
    if (length(bad) > 0) {
      stop("Invalid fixed parameter names: ",
           paste(bad, collapse = ", "),
           ". Valid names: ", paste(valid_fixed, collapse = ", "),
           ", or 'shapes' for all.", call. = FALSE)
    }
    fixed <- unique(fixed)
  }

  if (type == "g3") {
    obj <- list(
      type    = type,
      tau     = tau,
      gamma   = gamma,
      alpha   = alpha,
      eta     = eta,
      formula = formula,
      fixed   = fixed
    )
  } else {
    obj <- list(
      type    = type,
      t_half  = if (type == "constant") NA_real_ else t_half,
      nu      = if (type == "constant") NA_real_ else nu,
      m       = if (type == "constant") NA_real_ else m,
      formula = formula,
      fixed   = fixed
    )
  }

  structure(obj, class = "hzr_phase")
}


# ============================================================================
# S3 methods
# ============================================================================

#' @rdname hzr_phase
#' @param x An `hzr_phase` object (for `print.hzr_phase()`).
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_phase <- function(x, ...) {
  label <- switch(x$type,
    cdf      = "cdf (early risk)",
    hazard   = "hazard (late risk)",
    constant = "constant (flat rate)",
    g3       = "g3 (late phase)"
  )
  cat("<hzr_phase>", label, "\n")

  if (x$type == "g3") {
    cat("  tau =", format(x$tau, digits = 4),
        " gamma =", format(x$gamma, digits = 4),
        " alpha =", format(x$alpha, digits = 4),
        " eta =", format(x$eta, digits = 4), "\n")
  } else if (x$type != "constant") {
    cat("  t_half =", format(x$t_half, digits = 4),
        " nu =", format(x$nu, digits = 4),
        " m =", format(x$m, digits = 4), "\n")
  }

  if (x$type != "constant" && length(x$fixed) > 0) {
    cat("  fixed:", paste(x$fixed, collapse = ", "), "\n")
  }

  if (!is.null(x$formula)) {
    cat("  covariates:", deparse(x$formula), "\n")
  }

  invisible(x)
}


# ============================================================================
# Type-checking helper
# ============================================================================

#' Test if an object is an hzr_phase
#'
#' @param x Object to test.
#' @return Logical scalar.
#'
#' @examples
#' is_hzr_phase(hzr_phase("cdf"))
#' is_hzr_phase("not a phase")
#'
#' @export
is_hzr_phase <- function(x) {
  inherits(x, "hzr_phase")
}


# ============================================================================
# Internal helpers for likelihood construction
# ============================================================================

#' Number of shape parameters for a phase
#'
#' Returns 3 (t_half, nu, m) for `"cdf"` and `"hazard"` phases, 0 for
#' `"constant"`.
#'
#' @param phase An `hzr_phase` object.
#' @return Integer: 3 or 0.
#' @keywords internal
.hzr_phase_n_shape <- function(phase) {
  stopifnot(is_hzr_phase(phase))
  switch(phase$type,
    constant = 0L,
    g3       = 4L,
    3L  # cdf, hazard
  )
}


#' Generate parameter names for a phase's sub-vector
#'
#' Builds the named labels for the parameter block that belongs to a single
#' phase in the full theta vector.  The block layout is:
#'
#' \itemize{
#'   \item For `"cdf"`/`"hazard"`: `[log_mu, log_t_half, nu, m, beta_1, ..., beta_p]`
#'   \item For `"constant"`:       `[log_mu, beta_1, ..., beta_p]`
#' }
#'
#' @param phase An `hzr_phase` object.
#' @param phase_name Character label for the phase (e.g. `"early"`, `"phase_1"`).
#' @param covariate_names Character vector of covariate column names that
#'   this phase uses.  Can be length 0 if no covariates.
#' @return Character vector of parameter names.
#' @keywords internal
.hzr_phase_theta_names <- function(phase, phase_name, covariate_names = character(0)) {
  stopifnot(is_hzr_phase(phase))

  # Intercept (always present)
  names_out <- paste0(phase_name, ".log_mu")

  # Shape parameters
  if (phase$type == "g3") {
    names_out <- c(names_out,
      paste0(phase_name, ".log_tau"),
      paste0(phase_name, ".gamma"),
      paste0(phase_name, ".alpha"),
      paste0(phase_name, ".eta")
    )
  } else if (phase$type != "constant") {
    names_out <- c(names_out,
      paste0(phase_name, ".log_t_half"),
      paste0(phase_name, ".nu"),
      paste0(phase_name, ".m")
    )
  }

  # Covariate coefficients
  if (length(covariate_names) > 0L) {
    names_out <- c(names_out,
      paste0(phase_name, ".", covariate_names)
    )
  }

  names_out
}


#' Total number of parameters for a phase
#'
#' Returns the count of free parameters in the theta sub-vector for one phase:
#' 1 (log_mu) + n_shape + n_covariates.
#'
#' @param phase An `hzr_phase` object.
#' @param n_covariates Integer; number of covariate columns this phase uses.
#' @return Integer.
#' @keywords internal
.hzr_phase_n_params <- function(phase, n_covariates = 0L) {
  stopifnot(is_hzr_phase(phase))
  1L + .hzr_phase_n_shape(phase) + as.integer(n_covariates)
}


#' Extract starting values from a phase specification
#'
#' Returns initial theta sub-vector on the estimation (internal) scale:
#' log(mu), log(t_half), nu, m, followed by zeros for covariate coefficients.
#'
#' @param phase An `hzr_phase` object.
#' @param n_covariates Integer; number of covariate columns.
#' @param mu_start Numeric scalar; initial scale parameter (default 0.1).
#' @return Named numeric vector of starting values.
#' @keywords internal
.hzr_phase_start <- function(phase, n_covariates = 0L, mu_start = 0.1) {
  stopifnot(is_hzr_phase(phase))
  stopifnot(
    "mu_start must be a positive scalar" =
      is.numeric(mu_start) && length(mu_start) == 1L && mu_start > 0
  )

  vals <- log(mu_start)  # log_mu

  if (phase$type == "g3") {
    vals <- c(vals, log(phase$tau), phase$gamma, phase$alpha, phase$eta)
  } else if (phase$type != "constant") {
    vals <- c(vals, log(phase$t_half), phase$nu, phase$m)
  }

  # Covariate betas initialized to zero
  if (n_covariates > 0L) {
    vals <- c(vals, rep(0, n_covariates))
  }

  vals
}


#' Build a logical mask of free (non-fixed) parameters in the full theta vector
#'
#' Returns a logical vector the same length as the full theta where `TRUE`
#' means the parameter is free (to be optimized) and `FALSE` means it is
#' held fixed at its starting value.
#'
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @return Logical vector of length `sum(.hzr_phase_n_params(...))`.
#' @keywords internal
.hzr_phase_free_mask <- function(phases, covariate_counts) {
  mask <- logical(0)

  for (nm in names(phases)) {
    ph <- phases[[nm]]
    n_cov <- covariate_counts[[nm]]

    # log_mu is always free
    mask <- c(mask, TRUE)

    # Shape parameters
    if (ph$type == "g3") {
      fixed <- if (is.null(ph$fixed)) character(0) else ph$fixed
      mask <- c(mask,
        !("tau"   %in% fixed),  # log_tau
        !("gamma" %in% fixed),  # gamma
        !("alpha" %in% fixed),  # alpha
        !("eta"   %in% fixed)   # eta
      )
    } else if (ph$type != "constant") {
      fixed <- if (is.null(ph$fixed)) character(0) else ph$fixed
      mask <- c(mask,
        !("t_half" %in% fixed),  # log_t_half
        !("nu"     %in% fixed),  # nu
        !("m"      %in% fixed)   # m
      )
    }

    # Covariate betas are always free
    if (n_cov > 0L) {
      mask <- c(mask, rep(TRUE, n_cov))
    }
  }

  mask
}


#' Validate a list of phase specifications
#'
#' Checks that `phases` is a non-empty named list of `hzr_phase` objects.
#' Auto-names unnamed phases as `phase_1`, `phase_2`, etc.
#'
#' @param phases A list of `hzr_phase` objects.
#' @return The validated (and possibly auto-named) list, invisibly.
#' @keywords internal
.hzr_validate_phases <- function(phases) {
  # Catch bare hzr_phase passed instead of list(hzr_phase(...))
  if (is_hzr_phase(phases)) {
    stop("phases must be a non-empty list of hzr_phase objects, ",
         "not a single hzr_phase.  Wrap it: list(hzr_phase(...)).", call. = FALSE)
  }
  if (!is.list(phases) || length(phases) == 0L) {
    stop("phases must be a non-empty list of hzr_phase objects.", call. = FALSE)
  }

  # Check each element

  for (i in seq_along(phases)) {
    if (!is_hzr_phase(phases[[i]])) {
      stop("phases[[", i, "]] is not an hzr_phase object. ",
           "Use hzr_phase() to create phase specifications.", call. = FALSE)
    }
  }

  # Auto-name unnamed phases
  nms <- names(phases)
  if (is.null(nms)) {
    nms <- character(length(phases))
  }
  empty <- nms == "" | is.na(nms)
  if (any(empty)) {
    nms[empty] <- paste0("phase_", which(empty))
  }

  # Check for duplicate names
  if (anyDuplicated(nms)) {
    stop("Phase names must be unique. Duplicates found: ",
         paste(nms[duplicated(nms)], collapse = ", "), call. = FALSE)
  }

  names(phases) <- nms
  invisible(phases)
}
