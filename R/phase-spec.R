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
#'   \item{`"cdf"`: early, resolving risk}{Named for the **c**umulative
#'     **d**istribution **f**unction: the phase contributes \eqn{\Phi(t) = G(t)},
#'     the bounded CDF of the temporal decomposition (\eqn{0} at \eqn{t = 0},
#'     rising to a ceiling of \eqn{1}).  Because it saturates, the *hazard* it
#'     adds, \eqn{\mu\,g(t)}, peaks early and then decays toward zero (the
#'     signature of a one-time insult that patients either succumb to or survive
#'     past, e.g. peri-operative mortality).  Shape set by `t_half`, `nu`, `m`.
#'     SAS/C equivalent: the Early (G1) phase.}
#'   \item{`"hazard"`: accumulating aging risk (G1 family)}{Named because the
#'     phase contributes a **cumulative hazard** built from the same G1 family:
#'     \eqn{\Phi(t) = -\log(1 - G(t))}, which is unbounded and monotone
#'     increasing.  Its hazard \eqn{\mu\,h(t)} rises without leveling off, so it
#'     models risk that grows as subjects age.  This is an alternative late-risk
#'     form derived from G1; for the original SAS/C late phase prefer `"g3"`.
#'     Shape set by `t_half`, `nu`, `m`.}
#'   \item{`"g3"`: late, rising risk (original C/SAS late phase)}{Named for
#'     the **G3** (third) decomposition family used by the original HAZARD
#'     program for the late phase.  It contributes \eqn{\Phi(t) = G_3(t)} from
#'     [hzr_decompos_g3()], an unbounded intensity with its own four-parameter
#'     shape (`tau`, `gamma`, `alpha`, `eta`) that is more flexible than the
#'     G1-derived `"hazard"` form for capturing accelerating late mortality
#'     (e.g. structural valve deterioration years after surgery).  Use this when
#'     reproducing classic three-phase HAZARD models.  SAS/C equivalent: the
#'     Late (G3) phase.}
#'   \item{`"constant"`: flat background rate}{A time-invariant hazard:
#'     \eqn{\Phi(t) = t}, so the added hazard \eqn{\mu} is constant (the
#'     exponential model).  It represents the steady, ongoing risk present at all
#'     follow-up times, independent of how long ago the time origin was.  Takes
#'     no shape parameters; only its scale \eqn{\mu} (and any covariates) is
#'     estimated.  SAS/C equivalent: the Constant (G2) phase.}
#' }
#'
#' The shape derivative \eqn{\varphi_j = d\Phi_j/dt} (which forms the
#' instantaneous hazard contribution \eqn{\mu_j\,\varphi_j(t)}) is \eqn{g(t)} for
#' `"cdf"`, \eqn{h(t)} for `"hazard"`, \eqn{g_3(t)} for `"g3"`, and \eqn{1} for
#' `"constant"`.
#'
#' @param type Character; the phase's temporal shape, one of `"cdf"`
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
#' @param tau Positive finite scalar; scale parameter for `"g3"` phases.
#'   SAS late: `TAU`.
#' @param gamma Positive finite scalar; time exponent for `"g3"` phases.
#'   SAS late: `GAMMA`.
#' @param alpha Non-negative scalar; shape parameter for `"g3"` phases.
#'   When `alpha > 0`, the generic G3 formula is used; `alpha = 0` gives the
#'   exponential limiting case.  SAS late: `ALPHA`.
#' @param eta Positive finite scalar; outer exponent for `"g3"` phases.
#'   SAS late: `ETA`.
#' @param formula Optional one-sided formula (e.g. `~ age + nyha`) for
#'   phase-specific covariates.  It is evaluated in the `data` given to
#'   [hazard()], so without `data` [hazard()] refuses it, unless it builds
#'   nothing either way: an intercept-only `~ 1` with no global `x`.
#'   The phase's scale parameter plays the role of an intercept, so the
#'   design never has one: removing it (`~ 0 + age`, `~ age - 1`) is
#'   ignored with a warning, and builds the design of `~ age`.
#'   When `NULL` (default), the phase inherits the global design from
#'   [hazard()]: the global formula's covariates, or `x` on the vector
#'   interface.
#' @param fixed Character vector naming shape parameters to hold fixed during
#'   optimization.  Valid names for `"cdf"`/`"hazard"`: `"t_half"`, `"nu"`,
#'   `"m"`, or `"shapes"` (shorthand for all three).  Valid names for `"g3"`:
#'   `"tau"`, `"gamma"`, `"alpha"`, `"eta"`, or `"shapes"` (shorthand for all
#'   four).  Fixed parameters are held at their starting values; only `mu`
#'   (and covariates) are estimated.  Ignored for `"constant"` phases.
#'   This mirrors the SAS/C HAZARD workflow where shapes are typically fixed
#'   and only scale parameters are estimated.
#' @param constraint For `"g3"` phases, a rule that *derives* one shape from
#'   the others rather than estimating it:
#'   \describe{
#'     \item{`"none"`}{(default) every shape is estimated or fixed.}
#'     \item{`"alpha_gamma_eta"`}{\eqn{\alpha = \gamma\eta/2}, so that
#'       \eqn{\gamma\eta/\alpha = 2}. SAS/C: `FIXGAE2`.}
#'     \item{`"eta_gamma"`}{\eqn{\eta = 2/\gamma}, so that
#'       \eqn{\gamma\eta = 2}. SAS/C: `FIXGE2`.}
#'   }
#'   At `alpha = 1`, `hzr_phase()`'s default, the g3 form is
#'   \eqn{(t/\tau)^{\gamma\eta}}{(t/tau)^(gamma*eta)} (see [hzr_decompos_g3()]),
#'   which depends on
#'   \eqn{\gamma} and \eqn{\eta} only through their product. So
#'   \eqn{\gamma} and \eqn{\eta} are not separately identified there: under
#'   `"eta_gamma"` the product is fixed at 2 and \eqn{\gamma} is not
#'   identified at all, and with both estimated only the product is. A fit
#'   started there can report convergence with an arbitrary \eqn{\gamma}.
#'   Start `alpha` away from 1, or fix \eqn{\gamma}, for such a phase.
#'   The derived parameter follows the others at every step of the
#'   optimization, so it is not a free parameter and cannot be named in
#'   `fixed`; `"shapes"` leaves it out. Its starting value is computed from the
#'   others, and a value you supply for it, here or in the `theta` given to
#'   [hazard()], is replaced, with a warning when it differs. Under
#'   `hazard(fit = FALSE)` that replacement is made only when no phase carries
#'   covariates, since otherwise its slot in `theta` is not known until the
#'   design is built; `hazard()` warns when it could not be made. Its
#'   standard error in [vcov()] is the delta-method one, carried from the
#'   parameters it is derived from, so confidence limits from [predict()]
#'   include its uncertainty. Only one constraint per phase: SAS's
#'   `FIXGE2` and `FIXGAE2` together force \eqn{\alpha = 1} and fix `tau`,
#'   `gamma` and `eta`, which is better written as those fixed values.
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
#'   \item{constraint}{The shape constraint (g3 phases).}
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
#' # Derive alpha from gamma and eta (SAS/C FIXGAE2)
#' late_gae2 <- hzr_phase("g3", tau = 14, gamma = 22, eta = 0.18,
#'                         constraint = "alpha_gamma_eta")
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
                      fixed = character(0),
                      constraint = c("none", "alpha_gamma_eta", "eta_gamma")) {

  type <- match.arg(type)
  constraint <- match.arg(constraint)
  if (constraint != "none" && type != "g3") {
    stop("constraint = \"", constraint, "\" applies only to g3 phases.",
         call. = FALSE)
  }

  # --- Validate shape parameters based on type ------------------------------
  if (type == "g3") {
    # G3 late-phase decomposition: 4 parameters (tau, gamma, alpha, eta)
    # Finite as well as positive: an infinite shape was accepted here and only
    # failed later, inside the optimizer.
    stopifnot(
      "tau must be a positive scalar, and finite" =
        is.numeric(tau) && length(tau) == 1L && tau > 0 && is.finite(tau),
      "gamma must be a positive scalar, and finite" =
        is.numeric(gamma) && length(gamma) == 1L && gamma > 0 &&
          is.finite(gamma),
      "alpha must be a non-negative scalar" =
        is.numeric(alpha) && length(alpha) == 1L && alpha >= 0 && is.finite(alpha),
      "eta must be a positive scalar, and finite" =
        is.numeric(eta) && length(eta) == 1L && eta > 0 && is.finite(eta)
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
    # The design is built as if the intercept were present (#303), so say
    # so here, once, rather than on every refit that builds the design.
    tt <- stats::terms(formula, allowDotAsName = TRUE)
    # `~ 0` alone builds no columns either way, so there is nothing to say.
    if (attr(tt, "intercept") == 0L && length(attr(tt, "term.labels")) > 0L) {
      warning("Phase formula `", paste(deparse(formula), collapse = " "),
              "` removes the intercept, which a phase formula cannot do: ",
              "the phase's scale parameter plays the intercept role. The ",
              "design is built as if the intercept were present, so factors ",
              "are coded as they would be with it.", call. = FALSE)
    }
  }

  # --- Validate and normalize fixed parameters ------------------------------
  if (length(fixed) > 0 && type == "constant") {
    warning("fixed is ignored for constant phases.", call. = FALSE)
    fixed <- character(0)
  }
  if (length(fixed) > 0) {
    if (type == "g3") {
      # G3: expand "shapes" to tau, gamma, alpha, eta -- less a derived one,
      # which is computed rather than held (see `constraint`).
      derived <- .hzr_constraint_derived(constraint)
      if ("shapes" %in% fixed) {
        fixed <- union(setdiff(fixed, "shapes"),
                       setdiff(c("tau", "gamma", "alpha", "eta"), derived))
      }
      if (length(derived) && derived %in% fixed) {
        stop(derived, " cannot be fixed under constraint = \"", constraint,
             "\": it is derived from the other shapes. Fix the parameter(s) ",
             "it is derived from instead.", call. = FALSE)
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
    # The derived starting value is computed, as PROC HAZARD's SETG3 does
    # (it rewrote a specified ALPHA = 2 to 1.98 under FIXGAE2). Say so only
    # when the caller supplied a value that the rule then replaced.
    if (constraint == "alpha_gamma_eta") {
      supplied <- if (missing(alpha)) NULL else alpha
      alpha <- gamma * eta / 2
    } else if (constraint == "eta_gamma") {
      supplied <- if (missing(eta)) NULL else eta
      eta <- 2 / gamma
    }
    # The checks above saw the sources, not the derived value. Finite sources
    # can still overflow (alpha = Inf) or underflow, so check it too. A
    # derived alpha must be > 0 even though a supplied one may be 0: alpha = 0
    # selects the exponential limiting case, which a user can choose by
    # fixing it, but reaching it through gamma * eta / 2 underflowing is a
    # silent switch of model family, not a choice, so it is refused.
    if (constraint != "none") {
      value <- if (constraint == "alpha_gamma_eta") alpha else eta
      if (!(is.finite(value) && value > 0)) {
        sources <- if (constraint == "alpha_gamma_eta") {
          paste0("gamma = ", format(gamma, digits = 6), ", eta = ",
                 format(eta, digits = 6))
        } else {
          paste0("gamma = ", format(gamma, digits = 6))
        }
        stop(.hzr_constraint_rule(constraint), " is not a finite positive ",
             "number (", sources, ").", call. = FALSE)
      }
    }
    if (constraint != "none" && !is.null(supplied)) {
      derived <- .hzr_constraint_derived(constraint)
      value <- if (derived == "alpha") alpha else eta
      if (!isTRUE(all.equal(supplied, value))) {
        warning(derived, " = ", format(supplied, digits = 6),
                " was replaced by ", format(value, digits = 6),
                ", the value constraint = \"", constraint, "\" derives.",
                call. = FALSE)
      }
    }
    obj <- list(
      type       = type,
      tau        = tau,
      gamma      = gamma,
      alpha      = alpha,
      eta        = eta,
      formula    = formula,
      fixed      = fixed,
      constraint = constraint
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

  rule <- .hzr_constraint_rule(.hzr_phase_constraint(x))
  if (length(rule)) cat("  derived:", rule, "\n")

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
#' Returns 3 (t_half, nu, m) for `"cdf"` and `"hazard"` phases, 4 (tau,
#' gamma, alpha, eta) for `"g3"`, and 0 for `"constant"`.
#'
#' @param phase An `hzr_phase` object.
#' @return Integer: 3, 4 or 0.
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
#'   \item For `"g3"`:             `[log_mu, log_tau, gamma, alpha, eta, beta_1, ..., beta_p]`
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


#' Parameter names for a phase specification, in `theta` order
#'
#' Returns the names of the `theta` vector a multiphase [hazard()] fit builds
#' from `phases`, in the order `theta` requires.  Use it to check a
#' hand-written starting vector against the specification it is meant to go
#' with, before any fit runs.
#'
#' @details
#' `theta` is positional and its entries are not on a common scale.  For an
#' `early` + `late` pair the layout is
#'
#' ```
#' early.log_mu, early.log_t_half, early.nu, early.m,
#' late.log_mu,  late.log_tau,     late.gamma, late.alpha, late.eta
#' ```
#'
#' so the late phase logs `mu` and `tau` while carrying `gamma`, `alpha` and
#' `eta` on the natural scale.  **Wrapping the wrong element in `log()`
#' produces a fit, not an error**, which is why a comment describing the order
#' is not enough and this returns the real thing.  The order is a property of
#' the specification, so it changes the moment a phase is added, removed or
#' retyped.
#'
#' Phase names come from `names(phases)`; unnamed phases are labelled
#' `phase_1`, `phase_2`, ... by the same validation [hazard()] applies, so the
#' labels here are the labels a fit will use.
#'
#' @param phases A non-empty list of [hzr_phase()] objects, as passed to
#'   [hazard()].
#' @param covariates Optional named list mapping phase name to that phase's
#'   covariate column names, e.g. `list(early = c("age", "sex"))`.  Phases
#'   absent from the list are treated as having no covariates.  Names are used
#'   verbatim, so they must match the columns the fit will see.
#'
#' @return A character vector whose **order is the required `theta` order**.
#'   Its length is the number of parameters the specification implies, so
#'   `length(hzr_theta_names(phases))` is the length `theta` must have.
#'
#' @examples
#' phases <- list(early = hzr_phase("cdf"), late = hzr_phase("g3"))
#' hzr_theta_names(phases)
#'
#' # Check a hand-written starting vector before fitting.
#' theta0 <- c(log(0.05), log(0.2), 0, -0.4, log(0.03), log(1), 1, 1, 1)
#' stopifnot(length(theta0) == length(hzr_theta_names(phases)))
#' setNames(theta0, hzr_theta_names(phases))
#'
#' # With covariates on one phase.
#' hzr_theta_names(phases, covariates = list(early = c("age", "sex")))
#'
#' @seealso [hzr_phase()] for building a specification, [hazard()] for the fit
#'   whose `theta` this names.
#' @export
#' @details Phase names must be unique, and `"total"` and `"time"` are
#'   reserved; this applies the same validation [hazard()] does, so both reject
#'   them identically.
hzr_theta_names <- function(phases, covariates = NULL) {
  # The same validation hazard() applies, so auto-named phases get the same
  # phase_1/phase_2 labels here as they will in the fit. Re-deriving them
  # would let the two drift.
  phases <- .hzr_validate_phases(phases)

  if (is.null(covariates)) {
    covariates <- list()
  }
  if (!is.list(covariates) ||
        (length(covariates) > 0L && is.null(names(covariates)))) {
    stop("'covariates' must be a named list mapping phase name to that ",
         "phase's covariate column names.", call. = FALSE)
  }
  unknown <- setdiff(names(covariates), names(phases))
  if (length(unknown)) {
    # Silently ignoring these would return a name vector that is short by
    # exactly the covariates the caller thought they had asked for.
    stop("'covariates' names no such phase: ",
         paste0("'", unknown, "'", collapse = ", "),
         ". Phases are: ", paste0("'", names(phases), "'", collapse = ", "),
         ".", call. = FALSE)
  }
  for (nm in names(covariates)) {
    if (!is.character(covariates[[nm]])) {
      stop("covariates[['", nm, "']] must be a character vector of column ",
           "names.", call. = FALSE)
    }
  }

  .hzr_theta_names_list(phases, covariates)
}


#' Resolve each phase's covariate column names, with a positional fallback
#'
#' The full theta vector needs one name per covariate column, and `x_list` may
#' carry a matrix with no `colnames`.  Falling through with `character(0)`
#' there produces FEWER names than the phase has parameters, and
#' `names(theta) <- <short vector>` pads with `NA` rather than erroring, so
#' the misalignment is silent.  Substitute positional labels instead.
#'
#' @param phases Named list of `hzr_phase` objects.
#' @param covariate_counts Named list/vector of covariate counts per phase.
#' @param x_list Named list of per-phase covariate matrices, or `NULL`s.
#' @return Named list of character vectors, one per phase.
#' @noRd
.hzr_phase_cov_names <- function(phases, covariate_counts, x_list) {
  stats::setNames(lapply(names(phases), function(nm) {
    k <- covariate_counts[[nm]] %||% 0L
    if (k <= 0L) return(character(0))
    nms <- if (!is.null(x_list[[nm]])) colnames(x_list[[nm]]) else NULL
    if (length(nms) == k) nms else paste0("x", seq_len(k))
  }), names(phases))
}

#' Parameter names for a whole phase list, in theta order
#'
#' The list-level counterpart of [.hzr_phase_theta_names()].  Single source of
#' truth: `hazard()`'s optimiser, the score test's re-expansion and the
#' exported [hzr_theta_names()] all go through this, so the order the public
#' function documents cannot drift from the order a fit actually uses.
#'
#' @param phases Named list of `hzr_phase` objects, already validated.
#' @param cov_names Named list of character vectors, as returned by
#'   `.hzr_phase_cov_names()`.
#' @return Character vector, concatenated in phase order.
#' @noRd
.hzr_theta_names_list <- function(phases, cov_names = NULL) {
  unlist(lapply(names(phases), function(nm) {
    .hzr_phase_theta_names(phases[[nm]], nm, cov_names[[nm]] %||% character(0))
  }), use.names = FALSE)
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
#' log(mu), then log(t_half), nu, m for `"cdf"`/`"hazard"` phases or
#' log(tau), gamma, alpha, eta for `"g3"` (nothing for `"constant"`),
#' followed by zeros for covariate coefficients.
#'
#' @param phase An `hzr_phase` object.
#' @param n_covariates Integer; number of covariate columns.
#' @param mu_start Numeric scalar; initial scale parameter (default 0.1).
#' @return Unnamed numeric vector of starting values.
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
      # A derived shape is not searched over; it is recomputed from the others
      # (.hzr_apply_constraints()), so it is masked like a fixed one.
      fixed <- c(if (is.null(ph$fixed)) character(0) else ph$fixed,
                 .hzr_constraint_derived(.hzr_phase_constraint(ph)))
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


# ============================================================================
# Shape constraints (SAS/C FIXGAE2, FIXGE2)
# ============================================================================
#
# A constrained g3 phase keeps its full theta slot for the derived shape, so
# the layout, the names and every consumer of a full theta are unchanged. The
# optimizer masks the derived slot like a fixed one and recomputes it from the
# free shapes whenever it expands a reduced vector (.hzr_apply_constraints()).
# What a fixed mask alone would get wrong is the calculus: the derived slot
# moves with its sources, so the score, the Hessian and the covariance all
# carry its derivatives (.hzr_constraint_jacobian(),
# .hzr_constraint_curvature()). SAS/C does the same thing in hzd_late_t2p.c,
# which recomputes ALPHA = GAMMA*ETA/2 (or ETA = 2/GAMMA) from theta each time.

#' The constraint a phase carries; `"none"` for one built before the argument
#' existed, so a fit saved by an earlier version still reads.
#' @noRd
.hzr_phase_constraint <- function(ph) {
  if (is.null(ph$constraint)) "none" else ph$constraint
}

#' The shape a constraint derives (`character(0)` for `"none"`).
#' @noRd
.hzr_constraint_derived <- function(constraint) {
  switch(constraint,
         alpha_gamma_eta = "alpha",
         eta_gamma       = "eta",
         character(0))
}

#' The derivation as text, for print methods.
#' @noRd
.hzr_constraint_rule <- function(constraint) {
  switch(constraint,
         alpha_gamma_eta = "alpha = gamma * eta / 2",
         eta_gamma       = "eta = 2 / gamma",
         character(0))
}

#' Each derived theta entry at `theta`, with its first and second derivatives
#'
#' Positions come from the layout (`[log_mu, log_tau, gamma, alpha, eta,
#' betas...]` for g3, see `.hzr_unpack_phase_theta()`), never from names.
#'
#' @return A list with one element per constrained phase: `pos` (the derived
#'   slot), `value`, `src` (the slots it is derived from), `d1` (its gradient
#'   over `src`) and `d2` (its Hessian over `src`).
#' @noRd
.hzr_constraint_terms <- function(theta, phases, covariate_counts) {
  starts <- .hzr_log_mu_positions(phases, covariate_counts)
  terms <- list()
  for (nm in names(phases)) {
    constraint <- .hzr_phase_constraint(phases[[nm]])
    if (constraint == "none") next
    gamma_pos <- starts[[nm]] + 2L
    gamma_ <- theta[[gamma_pos]]
    if (constraint == "alpha_gamma_eta") {
      eta_pos <- starts[[nm]] + 4L
      eta_ <- theta[[eta_pos]]
      terms[[nm]] <- list(
        pos = starts[[nm]] + 3L, value = gamma_ * eta_ / 2,
        src = c(gamma_pos, eta_pos), d1 = c(eta_ / 2, gamma_ / 2),
        d2 = matrix(c(0, 0.5, 0.5, 0), 2L, 2L)
      )
    } else {
      terms[[nm]] <- list(
        pos = starts[[nm]] + 4L, value = 2 / gamma_,
        src = gamma_pos, d1 = -2 / gamma_^2,
        d2 = matrix(4 / gamma_^3, 1L, 1L)
      )
    }
  }
  terms
}

#' Recompute every derived entry of a full theta from its sources
#' @noRd
.hzr_apply_constraints <- function(theta, phases, covariate_counts) {
  for (term in .hzr_constraint_terms(theta, phases, covariate_counts)) {
    theta[term$pos] <- term$value
  }
  theta
}

#' Apply the constraints to a supplied theta, saying what was replaced
#'
#' `hzr_phase()` warns when a value passed for a derived shape is replaced;
#' a `theta` passed to [hazard()] deserves the same, or its derived entry
#' would change without a word.
#' @noRd
.hzr_constrain_supplied_theta <- function(theta, phases, covariate_counts) {
  out <- .hzr_apply_constraints(theta, phases, covariate_counts)
  moved <- which(!mapply(function(a, b) isTRUE(all.equal(a, b)), theta, out))
  if (length(moved)) {
    labels <- if (is.null(names(theta))) paste0("theta[", moved, "]") else
      names(theta)[moved]
    warning("theta ", paste(sprintf("%s = %s was replaced by %s", labels,
                                    format(theta[moved], digits = 6),
                                    format(out[moved], digits = 6)),
                            collapse = "; "),
            ", the value its phase's constraint derives.", call. = FALSE)
  }
  out
}

#' Jacobian of the constrained full theta with respect to itself
#'
#' The identity, except that a derived row holds its derivatives over its
#' sources and a derived column is zero (nothing depends on the derived slot
#' directly). `crossprod(J, score)` is then the score with each derived
#' entry's contribution folded into its sources, and `J V J'` carries a
#' covariance over the searched parameters onto the derived ones.
#' @noRd
.hzr_constraint_jacobian <- function(theta, phases, covariate_counts) {
  jac <- diag(length(theta))
  for (term in .hzr_constraint_terms(theta, phases, covariate_counts)) {
    jac[term$pos, ] <- 0
    jac[term$pos, term$src] <- term$d1
  }
  jac
}

#' Fold each derived entry's score into its sources
#'
#' `crossprod(.hzr_constraint_jacobian(), grad)` element by element, written
#' out so that an `NA` in an unrelated component (the unsanitised score) stays
#' where it is rather than spreading through the zeros of the matrix product.
#' The derived entry itself is left as it was: callers read only free slots.
#' @noRd
.hzr_constraint_score <- function(theta, grad, phases, covariate_counts) {
  for (term in .hzr_constraint_terms(theta, phases, covariate_counts)) {
    grad[term$src] <- grad[term$src] + grad[[term$pos]] * term$d1
  }
  grad
}

#' Second-order term of the constrained Hessian
#'
#' The Hessian of `f(c(theta))` is `J' H J + sum_k df/dtheta_k * d2(theta_k)`,
#' summed over the derived entries `k`. The second part is zero only where the
#' derived entry's own gradient vanishes, which it does not at a constrained
#' optimum, so leaving it out would misstate the curvature the covariance is
#' built from.
#'
#' @param grad Gradient of the same function `H` is the Hessian of.
#' @noRd
.hzr_constraint_curvature <- function(theta, grad, phases, covariate_counts) {
  out <- matrix(0, length(theta), length(theta))
  for (term in .hzr_constraint_terms(theta, phases, covariate_counts)) {
    out[term$src, term$src] <- out[term$src, term$src] +
      grad[[term$pos]] * term$d2
  }
  out
}


#' Does each phase have shape parameters at all?
#'
#' A `constant` phase is `mu` and nothing else. The saturated identifiability
#' message says `mu` survives while the shape parameters go flat, which is
#' vacuous for a phase that has none (the wording defect in #211).
#'
#' This asks whether the parameters *exist*, not whether they are free. A phase
#' whose shapes are pinned with `hzr_phase(fixed = )` still has them, and the
#' message is true and useful there: it is usually why they were pinned.
#'
#' @param phases A list of validated `hzr_phase` objects.
#' @return Logical vector, one element per phase, in `phases` order.
#' @keywords internal
.hzr_phase_has_shape <- function(phases) {
  vapply(phases, function(ph) !identical(ph$type, "constant"),
         logical(1), USE.NAMES = FALSE)
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

  # 'total' is the key .hzr_multiphase_cumhaz() stores the accumulator under, so
  # a phase of that name would lose its own contribution vector to it.
  if (any(nms == "total")) {
    stop("'total' is a reserved phase name. Rename the phase: it collides with ",
         "the accumulator holding the summed cumulative hazard.", call. = FALSE)
  }

  # 'time' is the first column of predict(decompose = TRUE)'s output, which then
  # stores each phase's cumulative hazard under the phase's name, so a phase of
  # that name overwrote the requested times with no error (#224).
  if (any(nms == "time")) {
    stop("'time' is a reserved phase name. Rename the phase: it collides with ",
         "the column of prediction times in the decomposed output of ",
         "predict().", call. = FALSE)
  }

  # Check for duplicate names
  if (anyDuplicated(nms)) {
    stop("Phase names must be unique. Duplicates found: ",
         paste(nms[duplicated(nms)], collapse = ", "), call. = FALSE)
  }

  names(phases) <- nms
  invisible(phases)
}


#' Would a phase formula build any columns in `data`?
#'
#' `FALSE` only for an intercept-only formula such as `~ 1`: no variable, and
#' no term (a constant term such as `log(2)` still builds a column).
#'
#' @param pf A one-sided phase formula.
#' @return A single logical.
#' @keywords internal
#' @noRd
.hzr_phase_formula_has_terms <- function(pf) {
  length(all.vars(pf)) > 0L ||
    length(attr(stats::terms(pf), "term.labels")) > 0L
}


#' Refuse a phase formula that has no `data` to be evaluated in
#'
#' A phase formula is evaluated in `data` and nowhere else, so without
#' `data` it used to be ignored: the phase took the global `x`, or no columns
#' at all, and the fit was a different model with no warning (#299). That
#' covers `~ 1` beside a global `x` too: the phase took `x`, where with
#' `data` it has no columns. Only an intercept-only formula with no `x`
#' builds the same design either way, and it is left alone.
#'
#' @param phases A list of validated `hzr_phase` objects.
#' @param data The `data` argument of [hazard()], possibly `NULL`.
#' @param x The global design matrix, possibly `NULL`.
#' @return `NULL`, invisibly; called for its error.
#' @keywords internal
#' @noRd
.hzr_check_phase_formula_data <- function(phases, data, x = NULL) {
  if (!is.null(data)) {
    return(invisible(NULL))
  }
  has_x <- !is.null(x) && NCOL(x) > 0L
  for (nm in names(phases)) {
    pf <- phases[[nm]]$formula
    has_terms <- !is.null(pf) && .hzr_phase_formula_has_terms(pf)
    if (is.null(pf) || !(has_x || has_terms)) next
    stop("Phase '", nm, "' has the formula `",
         paste(deparse(pf), collapse = " "), "`, but no `data` was ",
         "supplied, and a phase formula is evaluated only in `data`: ",
         "the formula would be ignored",
         if (has_x) ", and the phase would take the global `x` instead",
         ". Pass `data =` with the phase's variables as columns, or use ",
         "hazard(Surv(...) ~ ..., data = ...)",
         if (has_x && !has_terms) ", or drop `x`",
         ".", call. = FALSE)
  }
  invisible(NULL)
}
