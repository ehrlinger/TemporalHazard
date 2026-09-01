# decomposition.R -- Generalized temporal decomposition
#
# PURPOSE
# -------
# The unified parametric family decompos(t; t_half, nu, m) generates all
# temporal phase shapes used in multiphase hazard models.  It produces three
# quantities: G(t) (CDF), g(t) (density -- "early" pattern), and h(t) =
# g(t)/(1-G(t)) (hazard -- "late" pattern).
#
# Originally introduced by Blackstone, Naftel, and Turner (1986, JASA 81:615)
# and extended to longitudinal mixed-effects settings by Rajeswaran et al.
# (2018, Stat Methods Med Res 27:126).  Ported here from the mixhazard
# package with enhanced numerical stability and integrated phase-level
# helpers for the additive cumulative hazard model:
#
#   H(t | x) = sum_j  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)
#
# PARAMETER MAPPING FROM SAS/C HAZARD
# ------------------------------------
# SAS Early (G1): DELTA, RHO/THALF, NU, M  -->  t_half, nu, m
# SAS Late  (G3): TAU, GAMMA, ALPHA, ETA   -->  t_half, nu, m
# SAS Const (G2): (none)                   -->  type = "constant"
#
# The 4-parameter C/SAS parameterizations collapse onto this 3-parameter
# family.
#
# DELTA is NOT among them, and is NOT absorbed by the shape -- it is
# unimplemented, and delta = 0 is assumed throughout. Earlier comments here
# and in R/argument_mapping.R said "absorbed", which made the omission look
# deliberate and harmless. It is neither.
#
# The C DELTA controls a time transformation B(t) = (exp(delta*t) - 1)/delta
# that enters in three separate places, and R computes the delta = 0 branch of
# each:
#
#   1. rho is built from B(tHalf), not from tHalf   (src/common/hzd_set_rho.c)
#   2. the time argument is B(t), not t             (hzd_ln_G1_and_SG1.c:60)
#   3. dBt = delta*T enters lnSG1 ADDITIVELY, so the density carries a factor
#      of exp(delta*T)                              (hzd_ln_G1_and_SG1.c:84)
#
# So a PROC HAZARD job with DELTA != 0 is reproduced here against a different
# function. The condition is detectable -- OUTHAZ carries the estimate and the
# .lst flags block prints "Delta = 0  Yes" -- and the SAS-facing paths check
# it: R/read-outhaz.R stops, .hzr_parse_parms() records it as untranslated,
# and the .lst natural-estimates parser warns. Implementing DELTA is a
# separate question (see #181).

# ============================================================================
# hzr_decompos -- core engine
# ============================================================================

#' Log-scale terms shared by the m > 0 branches of [hzr_decompos()]
#'
#' Cases 1 (`nu > 0`) and 3 (`nu < 0`) evaluate the same intermediate
#' quantities; only the sign of `rho` and the final assembly of `G` differ.
#' Forming them directly overflows well inside the range an optimizer can
#' reach: `2^m` is `Inf` from `m = 1024`, and `bt^(-1/nu)` overflows sooner
#' still (from `m ~ 750` at `t/t_half = 0.5`, `m*nu = 3`).  `btnu` then
#' becomes `Inf` and `btnu^(-1/m)` collapses silently to 0.
#'
#' The `m` in `rho`'s `(2^m - 1)/m` factor cancels the explicit `m`
#' multiplier, so exactly
#'   \deqn{m \, b(t)^{-1/\nu} = (t_{1/2}/t)^{1/\nu} (2^m - 1),}
#' and `log(btnu)` follows from a softplus of the log of that product.  Every
#' intermediate below stays finite for any `m > 0` a double can hold.
#'
#' @param time Numeric vector of positive times.
#' @param t_half Positive scalar.
#' @param nu Nonzero scalar.
#' @param m Positive scalar.
#' @return List with `log_S` (\eqn{= -\log(\mathrm{btnu})/m}) and `log_g`.
#' @keywords internal
.hzr_decompos_g1_logs <- function(time, t_half, nu, m) {
  # log(2^m - 1) = x + log(1 - exp(-x)) for x = m*log(2).  Both tails are
  # exact: hzr_log1mexp() switches to log(-expm1(-x)) for small x, where the
  # naive log1p(-2^(-m)) loses digits from about m = 1e-3 down and returns
  # -Inf once 2^(-m) rounds to 1; and for large m, exp(-x) underflows to 0,
  # leaving the exact m*log(2).
  x_log2   <- m * log(2)
  log_2m1  <- x_log2 + hzr_log1mexp(x_log2)

  # log(rho) = log|nu| + log(t_half) + nu*log((2^m - 1)/m)
  q       <- nu * (log_2m1 - log(m))
  log_rho <- log(abs(nu)) + log(t_half) + q
  log_bt  <- log(time) - log(t_half) - q

  # a = log(m * bt^(-1/nu)) = (1/nu)*log(t_half/t) + log(2^m - 1)
  a        <- log(t_half / time) / nu + log_2m1
  log_btnu <- hzr_log1pexp(a)

  # Softplus excess e = log(btnu) - a = log(1 + exp(-a)), evaluated as
  # log1pexp(-a) so neither tail overflows.  Subtracting it below removes an
  # O(m) cancellation: the naive
  #   log_g = (-1/m - 1)*log_btnu + (-1/nu - 1)*log_bt - log_rho
  # differences two terms of size ~m*log(2) to reach an O(1) answer, losing
  # ~1e-13 relative accuracy by m = 2000.
  e     <- hzr_log1pexp(-a)
  log_S <- -log_btnu / m

  list(log_S = log_S,
       log_g = log_S - e - log(m) - log_bt - log_rho)
}

#' Generalized temporal decomposition
#'
#' Computes the cumulative distribution \eqn{G(t)}, density \eqn{g(t)}, and
#' hazard \eqn{h(t) = g(t)/(1 - G(t))} for the parametric family defined by
#' half-life, time exponent, and shape.  This single function generates all
#' temporal phase shapes used in multiphase hazard models.
#'
#' @section Parameter mapping from SAS/C HAZARD:
#'
#' The original C code used separate parameterizations for early (DELTA,
#' RHO/THALF, NU, M) and late (TAU, GAMMA, ALPHA, ETA) phases.  Both
#' collapse onto the three parameters here.  See
#' [hzr_argument_mapping()] for the full translation table.
#'
#' @section Valid parameter combinations:
#'
#' Six cases are defined by the signs of `nu` and `m`:
#'
#' \tabular{lll}{
#'   **Case** \tab **Sign** \tab **Behavior** \cr
#'   1     \tab m > 0, nu > 0 \tab Standard sigmoidal \cr
#'   1L    \tab m = 0, nu > 0 \tab Exponential-like (Weibull CDF) \cr
#'   2     \tab m < 0, nu > 0 \tab Heavy-tailed \cr
#'   2L    \tab m < 0, nu = 0 \tab Exponential decay \cr
#'   3     \tab m > 0, nu < 0 \tab Bounded cumulative \cr
#'   3L    \tab m = 0, nu < 0 \tab Bounded exponential \cr
#' }
#'
#' The combination m < 0 **and** nu < 0 is undefined and raises an error.
#' `nu = 0` is supported only with m < 0 (Case 2L, the exponential-decay
#' limit); `nu = 0` with m >= 0 has no usable limiting form and raises an
#' error.
#'
#' @section Mathematical form:
#'
#' The construction fixes a rate \eqn{\rho} so that \eqn{G(t_{1/2}) = 0.5}
#' exactly.  For the base case (\eqn{m > 0,\ \nu > 0}):
#'
#' \deqn{\rho = \nu \, t_{1/2} \left(\frac{2^m - 1}{m}\right)^{\!\nu},
#'       \qquad b(t) = \frac{\nu t}{\rho}}
#'
#' The CDF and density are then
#'
#' \deqn{G(t) = \bigl(1 + m \, b(t)^{-1/\nu}\bigr)^{-1/m}, \qquad
#'       g(t) = \bigl(1 + m \, b(t)^{-1/\nu}\bigr)^{-1/m - 1}
#'              \, b(t)^{-1/\nu - 1} / \rho}
#'
#' and the hazard is \eqn{h(t) = g(t) / (1 - G(t))}.  The remaining five cases
#' in the table above arise as limits (\eqn{m \to 0}, \eqn{\nu \to 0}) or sign
#' reflections of this base form; the implementation dispatches to the
#' appropriate branch after inspecting the signs of `nu` and `m`.  See
#' `vignette("mf-mathematical-foundations")` for the full derivation of every
#' case.
#'
#' @param time Numeric vector of times (must be > 0).
#' @param t_half Half-life: time at which \eqn{G(t_{1/2}) = 0.5}.
#'   Must be > 0.
#' @param nu Time exponent controlling rate dynamics.
#'   SAS early: `NU`.  SAS late: relates to `GAMMA`/`ETA`.
#' @param m Shape exponent controlling the distributional form.
#'   SAS early: `M`.  SAS late: relates to `GAMMA`/`ALPHA`.
#'
#' @return A named list with three numeric vectors, each the same length
#'   as `time`:
#' \describe{
#'   \item{G}{Cumulative distribution \eqn{G(t) \in [0, 1]}.}
#'   \item{g}{Density \eqn{g(t) = dG/dt \ge 0}.  The "early" phase
#'     temporal pattern.}
#'   \item{h}{Hazard \eqn{h(t) = g(t)/(1 - G(t)) \ge 0}.  The "late"
#'     phase temporal pattern.}
#' }
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
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#'
#' # Case 1: standard sigmoidal (m > 0, nu > 0)
#' d1 <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)
#' plot(t_grid, d1$G, type = "l", main = "CDF (m=1, nu=2)")
#'
#' # Case 1L: Weibull-like (m = 0, nu > 0)
#' d1L <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
#'
#' # Case 2: heavy-tailed (m < 0, nu > 0)
#' d2 <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = -1)
#'
#' @seealso [hzr_phase_cumhaz()] for the phase-level cumulative hazard
#'   contribution, [hzr_argument_mapping()] for SAS/C parameter mapping,
#'   [hzr_phase()] for specifying phases in [hazard()] models.
#'
#' \code{vignette("mf-mathematical-foundations")} for the full derivation.
#'
#' @export
hzr_decompos <- function(time, t_half, nu, m) {

  # --- Input validation ------------------------------------------------------
  stopifnot(
    "time must be numeric"   = is.numeric(time),
    "t_half must be a positive scalar" =
      is.numeric(t_half) && length(t_half) == 1L && t_half > 0,
    "nu must be a numeric scalar" =
      is.numeric(nu) && length(nu) == 1L && is.finite(nu),
    "m must be a numeric scalar" =
      is.numeric(m) && length(m) == 1L && is.finite(m)
  )

  if (m < 0 && nu < 0) {
    stop("Decomposition undefined when both m < 0 and nu < 0 (m = ",
         m, ", nu = ", nu, ").", call. = FALSE)
  }

  # Clamp time away from zero to avoid division-by-zero / 0^negative

  time <- pmax(time, .Machine$double.xmin)

  # Pre-compute recurring exponent terms
  if (m != 0) mm1 <- -(1 / m) - 1
  if (nu != 0) num1 <- -(1 / nu) - 1

  # --- Case dispatch ---------------------------------------------------------
  if (m > 0 && nu > 0) {
    # Case 1: standard sigmoidal.  Evaluated on the log scale via
    # .hzr_decompos_g1_logs(); the direct form overflowed for large m.
    lg    <- .hzr_decompos_g1_logs(time, t_half, nu, m)
    G     <- exp(lg$log_S)
    g     <- exp(lg$log_g)

  } else if (m == 0 && nu > 0) {
    # Case 1L: Weibull-like (m -> 0 limit)
    rho   <- nu * t_half * (log(2)^nu)
    bt    <- nu * time / rho
    btnu  <- bt^(-1 / nu)
    G     <- exp(-btnu)
    g     <- G * (bt^num1) / rho

  } else if (m < 0 && nu > 0) {
    # Case 2: heavy-tailed.  `1 - 2^m` loses every significant digit once
    # 2^m falls below eps (m < -52): it rounds to exactly 1, so
    # (1 - 2^m)^(-nu) - 1 is 0, rho is Inf, and G collapses silently to 0.
    # The decay starts well before the collapse -- at m = -20 the direct form
    # is already wrong in the 10th digit.  hzr_log1mexp() evaluates
    # log(1 - 2^m) without the cancellation and expm1() recovers the
    # difference, holding full precision down to where 2^m itself underflows.
    log_1m2m <- hzr_log1mexp(-m * log(2))       # log(1 - 2^m)
    dm       <- expm1(-nu * log_1m2m)           # (1 - 2^m)^(-nu) - 1
    log_bt   <- log1p(time * dm / t_half)       # bt = 1 + nu*time/rho
    log_btnu <- hzr_log1mexp(log_bt / nu)       # log(1 - bt^(-1/nu))
    log_rho  <- log(nu) + log(t_half) - log(dm)
    G     <- exp(-log_btnu / m)
    g     <- exp(mm1 * log_btnu + num1 * log_bt - log(-m) - log_rho)

  } else if (m < 0 && nu == 0) {
    # Case 2L: exponential decay (nu -> 0 limit)
    rho   <- -t_half / log(1 - 2^m)
    bt    <- exp(-time / rho)
    btm   <- 1 - bt
    G     <- btm^(-1 / m)
    g     <- -(btm^mm1) * bt / (m * rho)

  } else if (m > 0 && nu < 0) {
    # Case 3: bounded cumulative (C HAZARD g1flag = 5, "Mixed Generic")
    # rho uses the (2^m - 1)/m form (matching Case 1), NOT a bare (2^m - 1).
    # Without the /m divisor the CDF carried a spurious factor of m on the
    # bt^(-1/nu) term, diverging from the C G1 evaluator by up to ~0.2 and
    # breaking continuity with the m -> 0 limit (Case 3L).  With /m the m
    # factors cancel and G reduces to 1 - (bt^(-1/nu) + 1)^(-1/m), matching
    # C exactly (verified against src/common/hzd_ln_G1_and_SG1.c case 5).
    # Shares every intermediate with Case 1 (the -nu factors in rho and bt
    # cancel); only the assembly of G differs.  expm1 keeps 1 - exp(log_S)
    # accurate when log_S is near 0.
    lg    <- .hzr_decompos_g1_logs(time, t_half, nu, m)
    G     <- -expm1(lg$log_S)
    g     <- exp(lg$log_g)

  } else if (m == 0 && nu < 0) {
    # Case 3L: bounded exponential (m -> 0 limit)
    rho   <- -nu * t_half * (log(2)^nu)
    bt    <- -nu * time / rho
    btnu  <- bt^(-1 / nu)
    G     <- 1 - exp(-btnu)
    g     <- exp(-btnu) * (bt^num1) / rho

  } else {
    # Remaining combination: nu == 0 with m >= 0.  The nu -> 0 limit is only
    # defined for m < 0 (Case 2L, exponential decay); for m >= 0 the limit is
    # degenerate (G collapses to a step function), so there is no usable form.
    # Fail loud rather than leaving G/g unassigned (which would raise the
    # cryptic "object 'G' not found").
    stop("Decomposition undefined for nu = ", nu, " with m = ", m, ". ",
         "The nu -> 0 limit is implemented only for m < 0; ",
         "supply a nonzero nu when m >= 0.", call. = FALSE)
  }

  # --- Hazard from density and CDF ------------------------------------------
  # h(t) = g(t) / (1 - G(t)), with guard against G(t) = 1
  one_minus_G <- pmax(1 - G, .Machine$double.xmin)
  h <- g / one_minus_G

  list(G = G, g = g, h = h)
}


# ============================================================================
# hzr_decompos_g3 -- Late-phase (G3) engine
# ============================================================================

#' Late-phase (G3) temporal decomposition
#'
#' Computes the cumulative intensity \eqn{G_3(t)} and its derivative
#' \eqn{g_3(t) = dG_3/dt} for the late-phase parametric family used in the
#' original Blackstone C/SAS HAZARD code.  Unlike [hzr_decompos()] (which
#' computes the early-phase G1 -- a bounded CDF), this function can produce
#' **unbounded** values, making it suitable for modelling increasing late risk.
#'
#' @section Mathematical form:
#'
#' When \eqn{\alpha > 0}:
#' \deqn{G_3(t) = \bigl(\bigl((t/\tau)^\gamma + 1\bigr)^{1/\alpha}
#'                       - 1\bigr)^\eta}
#'
#' When \eqn{\alpha = 0} (limiting exponential case):
#' \deqn{G_3(t) = \bigl(\exp\bigl((t/\tau)^\gamma\bigr) - 1\bigr)^\eta}
#'
#' @section Parameter mapping from SAS/C HAZARD:
#'
#' | SAS name | R argument | Role |
#' |----------|-----------|------|
#' | TAU      | `tau`     | Scale (time at which \eqn{(t/\tau) = 1}) |
#' | GAMMA    | `gamma`   | Power exponent on \eqn{t/\tau} |
#' | ALPHA    | `alpha`   | Shape (0 = exponential limiting case) |
#' | ETA      | `eta`     | Outer power exponent |
#'
#' @param time Numeric vector of times (must be > 0).
#' @param tau Positive scalar scale parameter.
#' @param gamma Positive scalar time exponent.
#' @param alpha Non-negative scalar shape parameter (0 selects limiting case).
#' @param eta Positive scalar outer exponent.
#'
#' @return A named list with two numeric vectors, each the same length
#'   as `time`:
#' \describe{
#'   \item{G3}{Cumulative intensity \eqn{G_3(t) \ge 0} (may exceed 1).}
#'   \item{g3}{Derivative \eqn{g_3(t) = dG_3/dt \ge 0}.}
#' }
#'
#' @references
#' Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of time-varying
#' hazard into phases, each incorporating a separate stream of concomitant
#' information. *J Am Stat Assoc.* 1986;81(395):615--624.
#' \doi{10.1080/01621459.1986.10478314}
#'
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#'
#' # Weibull-like: alpha = 1 gives G3(t) = (t/tau)^(gamma*eta)
#' d <- hzr_decompos_g3(t_grid, tau = 1, gamma = 3, alpha = 1, eta = 1)
#' plot(t_grid, d$G3, type = "l", main = "G3: power law (gamma=3)")
#'
#' # General case with alpha > 0
#' d2 <- hzr_decompos_g3(t_grid, tau = 2, gamma = 2, alpha = 0.5, eta = 1)
#'
#' @seealso [hzr_decompos()] for the early-phase (G1) decomposition,
#'   [hzr_phase_cumhaz()] for phase-level cumulative hazard helpers.
#'
#' @export
hzr_decompos_g3 <- function(time, tau, gamma, alpha, eta) {

  # --- Input validation ------------------------------------------------------
  # Soft checks: return Inf/0 for infeasible params (optimizer needs graceful

  # failure) rather than hard stops.
  if (!is.numeric(time)) stop("time must be numeric")
  if (!is.numeric(tau) || length(tau) != 1L) stop("tau must be a scalar")
  if (!is.numeric(gamma) || length(gamma) != 1L) stop("gamma must be a scalar")
  if (!is.numeric(alpha) || length(alpha) != 1L) stop("alpha must be a scalar")

  if (!is.numeric(eta) || length(eta) != 1L) stop("eta must be a scalar")

  # Return Inf/NaN for infeasible parameters (optimizer will see -Inf logl)
  if (tau <= 0 || gamma <= 0 || eta <= 0 || alpha < 0 || !is.finite(alpha)) {
    return(list(G3 = rep(Inf, length(time)), g3 = rep(NaN, length(time))))
  }

  # Clamp time away from zero
  time <- pmax(time, .Machine$double.xmin)

  # --- Common terms ----------------------------------------------------------
  # Work on log scale for numerical stability, mirroring the C implementation
  ln_t_tau  <- log(time / tau)                    # ln(T/tau)
  ln_t_tau_g <- gamma * ln_t_tau                   # gamma * ln(T/tau)

  # --- Case dispatch: alpha > 0 vs alpha = 0 --------------------------------
  if (alpha > 0) {
    # G3(t) = (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta
    #
    # On log scale:
    #   tGamma = ln((t/tau)^gamma + 1) = ln(exp(ln_t_tau_g) + 1) = log1p(exp(ln_t_tau_g))
    #   tEta   = ln(exp(tGamma/alpha) - 1) = ln(expm1(tGamma/alpha))
    #   lnG3   = eta * tEta

    tGamma <- .log1pexp(ln_t_tau_g)   # ln((t/tau)^gamma + 1)
    inner  <- tGamma / alpha           # ln((t/tau)^gamma + 1) / alpha
    tEta   <- .log_expm1(inner)        # ln(exp(inner) - 1)
    lnG3   <- eta * tEta

    G3 <- exp(lnG3)

    # g3(t) = dG3/dt
    # lnScale = ln(gamma) + ln(eta) - ln(tau) - ln(alpha)
    # lnSG3 = lnScale + (eta-1)*tEta + ((1-alpha)/alpha)*tGamma + (gamma-1)*ln(t/tau)
    lnScale <- log(gamma) + log(eta) - log(tau) - log(alpha)
    lnSG3   <- lnScale + (eta - 1) * tEta +
               ((1 - alpha) / alpha) * tGamma +
               (gamma - 1) * ln_t_tau

    g3 <- exp(lnSG3)

  } else {
    # alpha = 0 limiting case:
    # G3(t) = (exp((t/tau)^gamma) - 1)^eta
    #
    # tGamma = ln(exp(exp(ln_t_tau_g)) - 1)
    #        = ln(exp((t/tau)^gamma) - 1)
    #        = log(expm1((t/tau)^gamma))
    # lnG3 = eta * tGamma

    t_tau_g <- exp(ln_t_tau_g)         # (t/tau)^gamma
    tGamma  <- .log_expm1(t_tau_g)     # ln(exp((t/tau)^gamma) - 1)
    lnG3    <- eta * tGamma

    G3 <- exp(lnG3)

    # lnScale = ln(gamma) + ln(eta) - ln(tau)
    # lnSG3 = lnScale + (eta-1)*tGamma + (t/tau)^gamma + (gamma-1)*ln(t/tau)
    lnScale <- log(gamma) + log(eta) - log(tau)
    lnSG3   <- lnScale + (eta - 1) * tGamma +
               t_tau_g +
               (gamma - 1) * ln_t_tau

    g3 <- exp(lnSG3)
  }

  list(G3 = G3, g3 = g3)
}


# --- Numerically stable helpers for G3 decomposition -----------------------

#' Compute log(1 + exp(x)) with overflow/underflow protection
#' @param x Numeric vector.
#' @return Numeric vector: log(1 + exp(x)).
#' @keywords internal
.log1pexp <- function(x) {
  # For large x: log(1+exp(x)) ~ x
  # For moderate x: use log1p(exp(x))
  # For very negative x: log(1+exp(x)) ~ exp(x)
  out <- x
  big  <- x > 35
  mid  <- !big & x > -35
  small <- x <= -35

  out[big]   <- x[big]
  out[mid]   <- log1p(exp(x[mid]))
  out[small] <- exp(x[small])

  out
}


#' Compute log(exp(x) - 1) with numerical stability
#' @param x Numeric vector (must be > 0 for exp(x) > 1).
#' @return Numeric vector: log(exp(x) - 1).
#' @keywords internal
.log_expm1 <- function(x) {
  # For large x: log(exp(x) - 1) ~ x
  # For small positive x: use log(expm1(x)) for accuracy
  out <- x
  big <- x > 35
  ok  <- !big & x > 1e-10
  tiny <- x <= 1e-10

  out[big] <- x[big]
  out[ok]  <- log(expm1(x[ok]))
  # For very small x: expm1(x) ~ x, so log(expm1(x)) ~ log(x)
  out[tiny] <- log(pmax(x[tiny], .Machine$double.xmin))

  out
}


# ============================================================================
# Phase-level helpers for the additive cumulative hazard model
# ============================================================================

#' Cumulative hazard contribution from a single phase
#'
#' Computes \eqn{\Phi_j(t)} for one phase in the additive model
#' \eqn{H(t|x) = \sum_j \mu_j(x) \Phi_j(t)}.
#'
#' @param time Numeric vector of times (> 0).
#' @param t_half Half-life parameter (> 0).
#' @param nu Time exponent.
#' @param m Shape parameter.
#' @param type Phase type: `"cdf"` (early -- uses \eqn{G(t)}),
#'   `"hazard"` (late -- uses cumulative hazard from \eqn{h(t)}), or
#'   `"constant"` (flat rate -- \eqn{\Phi = t}).
#'
#' @return Numeric vector of cumulative hazard contributions \eqn{\Phi(t)},
#'   same length as `time`.
#'
#' @details
#' - `"cdf"`: \eqn{\Phi(t) = G(t)}.  Bounded \eqn{[0, 1]}.  Models early
#'   risk that resolves over time.
#' - `"hazard"`: \eqn{\Phi(t) = -\log(1 - G(t))}.  Monotone increasing.
#'   Models late or aging risk.  This is the cumulative hazard derived from
#'   the hazard function \eqn{h(t)}, since
#'   \eqn{\int_0^t h(s)\,ds = -\log(1 - G(t))}.
#' - `"constant"`: \eqn{\Phi(t) = t}.  Ignores `t_half`, `nu`, `m`.
#'   Equivalent to exponential (constant hazard rate).
#'
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#' phi_early <- hzr_phase_cumhaz(t_grid, t_half = 2, nu = 2, m = 0,
#'                                type = "cdf")
#' phi_late  <- hzr_phase_cumhaz(t_grid, t_half = 5, nu = 1, m = 0,
#'                                type = "hazard")
#' phi_const <- hzr_phase_cumhaz(t_grid, type = "constant")
#'
#' @seealso [hzr_decompos()] for the underlying parametric family,
#'   [hzr_phase_hazard()] for the instantaneous hazard contribution.
#'
#' @export
hzr_phase_cumhaz <- function(time, t_half = 1, nu = 1, m = 0,
                              type = c("cdf", "hazard", "constant")) {
  type <- match.arg(type)

  if (type == "constant") {
    return(time)
  }

  d <- hzr_decompos(time, t_half = t_half, nu = nu, m = m)

  switch(type,
    cdf    = d$G,
    hazard = {
      # Cumulative hazard from h(t): integral_0^t h(s)ds = -log(1 - G(t))
      one_minus_G <- pmax(1 - d$G, .Machine$double.xmin)
      -log(one_minus_G)
    }
  )
}


# ============================================================================
# Phase-level derivatives for analytic gradient
# ============================================================================

#' Derivatives of phase cumulative and instantaneous hazard w.r.t. shape params
#'
#' Computes \eqn{\Phi_j(t)}, \eqn{\phi_j(t)}, and their derivatives with
#' respect to `t_half`, `nu`, and `m` using central finite differences on
#' [hzr_decompos()].  The `log_t_half` derivative is obtained via the chain
#' rule: \eqn{d\Phi/d(\log t_{1/2}) = t_{1/2} \cdot d\Phi/dt_{1/2}}.
#'
#' @param time Numeric vector of positive times.
#' @param t_half Positive scalar half-life.
#' @param nu Numeric scalar time exponent.
#' @param m Numeric scalar shape exponent.
#' @param type Phase type: `"cdf"`, `"hazard"`, or `"constant"`.
#'
#' @return Named list:
#' \describe{
#'   \item{Phi}{Cumulative hazard contribution \eqn{\Phi(t)}.}
#'   \item{phi}{Instantaneous hazard contribution \eqn{\phi(t) = d\Phi/dt}.}
#'   \item{dPhi_dlog_thalf}{\eqn{d\Phi / d(\log t_{1/2})}.}
#'   \item{dPhi_dnu}{\eqn{d\Phi / d\nu}.}
#'   \item{dPhi_dm}{\eqn{d\Phi / dm}.}
#'   \item{dphi_dlog_thalf}{\eqn{d\phi / d(\log t_{1/2})}.}
#'   \item{dphi_dnu}{\eqn{d\phi / d\nu}.}
#'   \item{dphi_dm}{\eqn{d\phi / dm}.}
#' }
#' @keywords internal
.hzr_phase_derivatives <- function(time, t_half, nu, m,
                                    type = c("cdf", "hazard", "constant")) {
  type <- match.arg(type)
  n <- length(time)

  # Constant phases: Phi = t, phi = 1, no shape dependence

  if (type == "constant") {
    zeros <- rep(0, n)
    return(list(
      Phi = time,
      phi = rep(1, n),
      dPhi_dlog_thalf = zeros,
      dPhi_dnu        = zeros,
      dPhi_dm         = zeros,
      dphi_dlog_thalf = zeros,
      dphi_dnu        = zeros,
      dphi_dm         = zeros
    ))
  }

  # Base evaluation
  d0 <- hzr_decompos(time, t_half = t_half, nu = nu, m = m)

  # Extract Phi and phi from decomposition for the given phase type
  extract <- function(d, tp) {
    if (tp == "cdf") {
      Phi <- d$G
      phi <- d$g
    } else {
      # "hazard": Phi = -log(1 - G), phi = h
      one_minus_G <- pmax(1 - d$G, .Machine$double.xmin)
      Phi <- -log(one_minus_G)
      phi <- d$h
    }
    list(Phi = Phi, phi = phi)
  }

  base <- extract(d0, type)

  # Central-difference derivatives w.r.t. shape parameters
  eps_rel <- (.Machine$double.eps) ^ (1 / 3)  # optimal for central diff

  # Helper: perturbed decomposition (returns NULL if invalid params)
  perturb_decompos <- function(th, n_val, m_val) {
    if (m_val < 0 && n_val < 0) return(NULL)
    if (th <= 0) return(NULL)
    tryCatch(
      hzr_decompos(time, t_half = th, nu = n_val, m = m_val),
      error = function(e) NULL
    )
  }

  # Derivative w.r.t. t_half (then multiply by t_half for log_t_half chain rule)
  h_th <- eps_rel * max(abs(t_half), 1e-4)
  d_plus  <- perturb_decompos(t_half + h_th, nu, m)
  d_minus <- perturb_decompos(t_half - h_th, nu, m)
  if (!is.null(d_plus) && !is.null(d_minus)) {
    e_plus  <- extract(d_plus, type)
    e_minus <- extract(d_minus, type)
    dPhi_dt_half <- (e_plus$Phi - e_minus$Phi) / (2 * h_th)
    dphi_dt_half <- (e_plus$phi - e_minus$phi) / (2 * h_th)
  } else {
    # One-sided fallback
    if (!is.null(d_plus)) {
      e_plus <- extract(d_plus, type)
      dPhi_dt_half <- (e_plus$Phi - base$Phi) / h_th
      dphi_dt_half <- (e_plus$phi - base$phi) / h_th
    } else {
      dPhi_dt_half <- rep(0, n)
      dphi_dt_half <- rep(0, n)
    }
  }
  # Chain rule: d/d(log_t_half) = t_half * d/d(t_half)
  dPhi_dlog_thalf <- t_half * dPhi_dt_half
  dphi_dlog_thalf <- t_half * dphi_dt_half

  # Derivative w.r.t. nu
  h_nu <- eps_rel * max(abs(nu), 1)
  d_plus  <- perturb_decompos(t_half, nu + h_nu, m)
  d_minus <- perturb_decompos(t_half, nu - h_nu, m)
  if (!is.null(d_plus) && !is.null(d_minus)) {
    e_plus  <- extract(d_plus, type)
    e_minus <- extract(d_minus, type)
    dPhi_dnu <- (e_plus$Phi - e_minus$Phi) / (2 * h_nu)
    dphi_dnu <- (e_plus$phi - e_minus$phi) / (2 * h_nu)
  } else if (!is.null(d_plus)) {
    e_plus <- extract(d_plus, type)
    dPhi_dnu <- (e_plus$Phi - base$Phi) / h_nu
    dphi_dnu <- (e_plus$phi - base$phi) / h_nu
  } else if (!is.null(d_minus)) {
    e_minus <- extract(d_minus, type)
    dPhi_dnu <- (base$Phi - e_minus$Phi) / h_nu
    dphi_dnu <- (base$phi - e_minus$phi) / h_nu
  } else {
    dPhi_dnu <- rep(0, n)
    dphi_dnu <- rep(0, n)
  }

  # Derivative w.r.t. m
  h_m <- eps_rel * max(abs(m), 1)
  d_plus  <- perturb_decompos(t_half, nu, m + h_m)
  d_minus <- perturb_decompos(t_half, nu, m - h_m)
  if (!is.null(d_plus) && !is.null(d_minus)) {
    e_plus  <- extract(d_plus, type)
    e_minus <- extract(d_minus, type)
    dPhi_dm <- (e_plus$Phi - e_minus$Phi) / (2 * h_m)
    dphi_dm <- (e_plus$phi - e_minus$phi) / (2 * h_m)
  } else if (!is.null(d_plus)) {
    e_plus <- extract(d_plus, type)
    dPhi_dm <- (e_plus$Phi - base$Phi) / h_m
    dphi_dm <- (e_plus$phi - base$phi) / h_m
  } else if (!is.null(d_minus)) {
    e_minus <- extract(d_minus, type)
    dPhi_dm <- (base$Phi - e_minus$Phi) / h_m
    dphi_dm <- (base$phi - e_minus$phi) / h_m
  } else {
    dPhi_dm <- rep(0, n)
    dphi_dm <- rep(0, n)
  }

  list(
    Phi             = base$Phi,
    phi             = base$phi,
    dPhi_dlog_thalf = dPhi_dlog_thalf,
    dPhi_dnu        = dPhi_dnu,
    dPhi_dm         = dPhi_dm,
    dphi_dlog_thalf = dphi_dlog_thalf,
    dphi_dnu        = dphi_dnu,
    dphi_dm         = dphi_dm
  )
}


#' Instantaneous hazard contribution from a single phase
#'
#' Computes \eqn{\phi_j(t) = d\Phi_j/dt} for one phase -- the derivative of
#' the cumulative hazard contribution returned by [hzr_phase_cumhaz()].
#'
#' @inheritParams hzr_phase_cumhaz
#'
#' @return Numeric vector of instantaneous hazard contributions \eqn{\phi(t)},
#'   same length as `time`.
#'
#' @details
#' - `"cdf"`: \eqn{\phi(t) = g(t)} (density).
#' - `"hazard"`: \eqn{\phi(t) = h(t) = g(t)/(1-G(t))}.
#' - `"constant"`: \eqn{\phi(t) = 1}.
#'
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#' phi_early <- hzr_phase_hazard(t_grid, t_half = 2, nu = 2, m = 0,
#'                                type = "cdf")
#' phi_late  <- hzr_phase_hazard(t_grid, t_half = 5, nu = 1, m = 0,
#'                                type = "hazard")
#'
#' @seealso [hzr_decompos()] for the underlying parametric family,
#'   [hzr_phase_cumhaz()] for the cumulative version.
#'
#' @export
hzr_phase_hazard <- function(time, t_half = 1, nu = 1, m = 0,
                              type = c("cdf", "hazard", "constant")) {
  type <- match.arg(type)

  if (type == "constant") {
    return(rep(1, length(time)))
  }

  d <- hzr_decompos(time, t_half = t_half, nu = nu, m = m)

  switch(type,
    cdf    = d$g,
    hazard = d$h
  )
}
