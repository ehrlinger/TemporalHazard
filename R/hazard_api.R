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

#' Find phases fitted outside the support their parameterisation can carry
#'
#' A DIAGNOSIS run after the optimizer returns, never a constraint: it does
#' not move an estimate, and a fit that trips it still returns its fit. An
#' unbounded phase type (`.hzr_phase_type_unbounded()`) whose fitted `t_half`
#' lies below the first observed time has left the data, and `-log(1 - G)` is
#' then evaluated where `G` is essentially 1 (#444).
#'
#' The trigger is a plain FACT -- `t_half` is below the observed support --
#' with no tuned threshold. The MAGNITUDE does the discriminating and is
#' carried in `$detail`: the ratio, and `1 - G(t_min)`, the phase's remaining
#' mass at the first observed time, which is the mechanism itself rather than
#' a proxy for it. Measured: a suite fit hugging the edge of its data reads
#' 0.053, and the `cabgkul` fit that reported a log-likelihood of +290082
#' reads 3.8e-12.
#'
#' Keyed on the FITTED value, not the starting one. A fit that starts below
#' the data and converges inside it is not this defect.
#'
#' @param theta Fitted parameter vector, named `<phase>.<parameter>`.
#' @param phases The fitted spec's phase list.
#' @param time The observed times.
#' @param fitted Did a fit actually run?
#' @return List with `boundary` (`NA` not examined, `NULL` nothing found, or a
#'   list of records) and `reason` (why, when `NA`).
#' @keywords internal
#' @noRd
.hzr_boundary_check_impl <- function(theta, phases, time, fitted,
                                     time_lower = NULL, time_upper = NULL) {
  na_because <- function(reason) list(boundary = NA, reason = reason)
  if (!isTRUE(fitted)) {
    return(na_because("model not fitted"))
  }
  # A fit with no phases -- every single-distribution fit -- has nothing an
  # unbounded phase type could trip, so this is "examined, found nothing"
  # (NULL) rather than "could not look" (NA). Returning NA here listed
  # boundary_check as a lost capability on every weibull fit in the package,
  # which is noise, not a finding. Compare conservation_of_events, which is
  # likewise reported only where it is applicable.
  if (!length(phases) || is.null(names(phases))) {
    return(list(boundary = NULL, reason = NULL))
  }
  # EVERY observed time, not just `time`. For an interval-censored row the
  # interval is [time_lower, time_upper], and on a left-truncated fit
  # time_lower is the entry time -- both can lie below min(time). The record
  # calls this "the first observed time", so it has to be one (#444).
  all_t <- c(time, time_lower, time_upper)
  t_ok <- all_t[is.finite(all_t) & all_t > 0]
  if (!length(t_ok)) {
    return(na_because("no positive observed times"))
  }
  t_min <- min(t_ok)

  found <- list()
  # BY INDEX, not by name: `phases[[nm]]` resolves a duplicated name to the
  # FIRST element, so `names(phases) == c("a", "a")` iterated twice and
  # emitted two identical records for one phase.
  for (k in seq_along(phases)) {
    nm <- names(phases)[[k]]
    if (!nzchar(nm)) next
    # #448: a G1-based phase whose nu collapses toward 0 becomes a step at
    # t_half. Checked BEFORE the unbounded-type filter: a "cdf" phase is not
    # an unbounded type, so the `next` below would skip it entirely.
    step_rec <- .hzr_phase_step_record(nm, phases[[k]]$type, theta, time,
                                       time_lower = time_lower,
                                       time_upper = time_upper)
    if (!is.null(step_rec)) found[[length(found) + 1L]] <- step_rec
    if (!.hzr_phase_type_unbounded(phases[[k]]$type)) next
    key <- paste0(nm, ".log_t_half")
    if (!key %in% names(theta)) next
    t_half <- exp(unname(theta[[key]]))
    if (!is.finite(t_half) || t_half >= t_min) next
    g <- tryCatch(
      hzr_decompos(t_min, t_half = t_half,
                   nu = unname(theta[[paste0(nm, ".nu")]]),
                   m = unname(theta[[paste0(nm, ".m")]]))$G,
      error = function(e) NA_real_
    )
    found[[length(found) + 1L]] <- list(
      mechanism = "unbounded_phase",
      phase = nm,
      parameter = "t_half",
      detail = paste0(
        "phase '", nm, "' is of the unbounded type 'hazard' and its fitted ",
        "t_half (", format(t_half, digits = 4), ") is below the first ",
        "observed time (", format(t_min, digits = 4), "), a factor of ",
        format(t_min / t_half, digits = 3), ". Its remaining mass there, ",
        "1 - G(t_min), is ", format(1 - g, digits = 4),
        ". ",
        # `1e-6` IS a tuned number, and unlike the TRIGGER -- which stays a
        # plain fact with no threshold -- it decides which of three sentences
        # a reader sees. It is a wording cut-off, not a detection cut-off: no
        # record is created or suppressed by it, and the measured magnitude
        # is printed either way, so a reader who disagrees with the cut-off
        # still has the number. Named here so it is not mistaken for part of
        # the criterion.
        #
        # The magnitude decides this, so the sentence must not assert it
        # unconditionally: an earlier version said "G is near 1 ... may be a
        # supremum" whatever the number was, which contradicts the design
        # this check rests on -- the trigger states a fact, the magnitude
        # discriminates. Reproduced across reachable fits at 1 - G from
        # 0.0067 to 0.22; the second is a fifth of the mass remaining and no
        # runaway at all.
        if (!is.finite(g)) {
          paste0("The remaining mass could not be computed, so whether ",
                 "-log(1 - G) is diverging here is NOT KNOWN.")
        } else if (1 - g < 1e-6) {
          paste0("G is numerically 1 across the observed range, so ",
                 "-log(1 - G) diverges: the objective is unbounded above and ",
                 "a reported optimum is likely a supremum rather than a fit.")
        } else {
          paste0("The phase still carries mass beyond the first observation, ",
                 "so this is a fit sitting outside its data rather than one ",
                 "riding the divergence of -log(1 - G).")
        }
      )
    )
  }
  if (!length(found)) {
    return(list(boundary = NULL, reason = NULL))
  }
  list(boundary = found, reason = NULL)
}


#' The warning text for boundary findings
#' @param records The `$boundary` list.
#' @return A single string.
#' @keywords internal
#' @noRd
.hzr_boundary_message <- function(records) {
  # Each mechanism keeps its own lead-in: a setup hold (#415) was not
  # "fitted outside the observed support", and nor is a cdf phase that
  # collapsed to a step inside the data; saying so would name the wrong cause.
  lead <- c(unbounded_phase = "fitted outside the observed support: ",
            phase_discontinuity = "phase collapsed to a step: ",
            g3_alpha_one = "shape parameters held at setup: ",
            g3_fixge2_alpha_start = "starting value moved at setup: ")
  paste(vapply(records, function(r) {
    paste0(if (r$mechanism %in% names(lead)) lead[[r$mechanism]] else "",
           r$detail)
  }, character(1)), collapse = " ")
}

#' The warning condition for boundary findings
#'
#' Classed `hzr_<mechanism>` for EVERY mechanism present, not only the first
#' record's, so a handler for one mechanism is not defeated by another
#' record listed ahead of it; all inherit `hzr_boundary`.
#' @param records The `$boundary` list.
#' @return A warning condition.
#' @keywords internal
#' @noRd
.hzr_boundary_condition <- function(records) {
  mechanisms <- vapply(records, function(r) r$mechanism, character(1))
  structure(
    class = c(unique(paste0("hzr_", mechanisms)), "hzr_boundary", "warning",
              "condition"),
    list(message = .hzr_boundary_message(records), call = NULL)
  )
}


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
#' differently (on the odds of failure and the log-time location,
#' respectively), so they are separate parameterizations, not special cases
#' of this additive form.  Parameters are estimated
#' on an unconstrained internal scale (e.g. \eqn{\log\mu}, \eqn{\log t_{1/2}})
#' and transformed back for reporting; see
#' `vignette("mf-mathematical-foundations")`.
#'
#' @section Convergence:
#'
#' The optimizer stops when an iteration improves the log-likelihood by less
#' than `control$reltol` relative to its size, and on a flat ridge that can
#' happen well short of the maximum. SAS/C HAZARD accepts an optimum only on a
#' different test, the relative gradient
#' \eqn{\max_i |g_i| \max(|x_i|, 1) / \max(|\ell|, 1) \le \epsilon^{1/3}}{max_i
#' |g_i| max(|x_i|, 1) / max(|l|, 1) <= eps^(1/3)}, about 6e-6, and `hazard()`
#' applies it too: when the optimizer reports convergence and the test fails,
#' the fit is continued with [stats::nlm()] at SAS's tolerances, and the
#' continued point is kept if it improves the log-likelihood.
#'
#' Every fit records the result in `fit$fit$rel_gradient` and, when the
#' continuation improved the fit, its termination code in
#' `fit$fit$polish_code`. `print()` and `summary()` show both.
#' `rel_gradient` is `NA` when the test was not applied (the optimizer did
#' not report convergence) or the gradient cannot be evaluated at the
#' estimates; `NA` is never reported as a pass. Neither is it always a
#' failure: some routes to it, such as a non-converged stop, do say the
#' estimates are unreliable, while others, such as a finite-difference score
#' that needed a point where the log-likelihood is not usable, say nothing
#' against them. Which route it took is recorded in
#' `fit$fit$rel_gradient_reason`, `NA_character_` when the test did run, and
#' `print()` and `summary()` show it. Under Conservation of Events
#' the analytic score omits how the conserved scale moves, so the test is
#' computed from finite differences of the log-likelihood with that scale
#' re-solved, as SAS/C does; the continuation still uses the analytic score,
#' so a CoE fit can honestly end with the test not met. A warning is raised only for code 4, the
#' iteration limit (raise `control$maxit`), and code 5, where the
#' log-likelihood kept rising along some direction and the model may have no
#' maximum. Codes 2 and 3, where SAS/C prints a caution, are recorded without
#' one. The test is relative to the size of the log-likelihood, so a fit that
#' meets it is within SAS's tolerance of the maximum, not exactly at it.
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
#'   \item{`"weibull"`: monotone rising or falling hazard (default)}{The
#'     workhorse parametric model: \eqn{H(t \mid \mathbf{x}) = (\mu t)^\nu
#'     \exp(\eta)}, with hazard \eqn{h \propto t^{\nu - 1}}.  The single shape
#'     \eqn{\nu} makes risk increase over time (\eqn{\nu > 1}), decrease
#'     (\eqn{\nu < 1}), or stay flat (\eqn{\nu = 1}).  Use it as the default when
#'     a single monotone trend describes the hazard.}
#'   \item{`"exponential"`: constant hazard}{The memoryless special case
#'     \eqn{\nu = 1}: a time-invariant baseline rate, \eqn{H(t \mid \mathbf{x}) =
#'     \mu t \exp(\eta)}.  Use it when the event rate does not change with
#'     follow-up time (the constant background risk also appears as the
#'     `"constant"` phase in a multiphase model).}
#'   \item{`"loglogistic"`: unimodal (rise-then-fall) hazard}{A log-logistic
#'     proportional-odds form (covariates act multiplicatively on the odds of
#'     failure, \eqn{\exp(\eta)}, not as an AFT time shift) whose hazard rises to
#'     a single peak and then declines when the shape exceeds 1 (and is monotone
#'     decreasing otherwise), with heavier tails than the log-normal.  Use it
#'     when risk climbs to an early peak and then eases off.}
#'   \item{`"lognormal"`: early-peaking, resolving hazard}{An
#'     accelerated-failure-time form in which \eqn{\log} time is Gaussian; the
#'     hazard rises to an early peak and then decays toward zero.  Use it for
#'     risk that is concentrated early and resolves over time.}
#'   \item{`"multiphase"`: additive N-phase hazard}{Sums several phase shapes
#'     into one model, \eqn{H = \sum_j \mu_j(\mathbf{x}) \Phi_j(t)}, so the
#'     overall hazard can fall, level off, and rise again within one fit.
#'     Requires `phases`; see [hzr_phase()] for the available phase shapes.  This
#'     is the form that reproduces the classic C/SAS HAZARD models.}
#' }
#'
#' @param time Numeric follow-up time vector. A row whose time is 0 is
#'   dropped before fitting, whatever its status, as PROC HAZARD drops it
#'   (`TIME <= 0` is inadmissible there). The time tested is the row's upper
#'   bound: `time` for an exact or right-censored row, `time_upper` for a
#'   left- or interval-censored one; a lower bound or entry time of 0 is
#'   admissible. `hazard()` warns with the count (class
#'   `"hzr_time_zero_dropped"`), and every row stored on the fit is what
#'   remains. `fit$data$dropped_time_zero` is the count and
#'   `fit$data$dropped_time_zero_rows` their positions among the rows given.
#' @param status Numeric or logical event indicator vector, or a
#'   [survival::Surv()] object. A `Surv` is read by its `type`, exactly as the
#'   formula interface reads it, and a `time`, `time_lower` or `time_upper`
#'   that disagrees with it is an error.
#' @param time_lower Optional numeric vector whose role depends on `status`.
#'   Supplying it explicitly is **not** a no-op.
#'   * `status == 2` (interval-censored): the lower bound of the censoring
#'     interval, defaulting to `time`. Every `dist` reads it this way.
#'   * `status %in% c(0, 1)` (right-censored or event): the counting-process
#'     **entry time** when `0 < time_lower < time`, so the row contributes
#'     `H(time) - H(time_lower)`. Every `dist` reads it this way. A value of
#'     `0`, or equal to `time`, means no entry time, and left `NULL` the
#'     entry time is **`0`**.
#'   * `status == -1` (left-censored): not used; the bound is `time_upper`.
#'
#'   `time_lower = time` on the status 0 and 1 rows is the mixed-interval
#'   layout: exact and right-censored rows carry their own time, and only
#'   the status-2 rows carry a real lower bound. A subject cannot enter the
#'   risk set after the moment it leaves, so `hazard()` stops with an error
#'   when a status 0 or 1 row has `time_lower > time`, as SAS HAZARD rejects
#'   a start time after the exit time. It also stops when rows with
#'   `time_lower == time > 0` sit beside rows with genuine entry times: in
#'   counting-process data those are zero-length epochs, which
#'   [hzr_repeated_events()] can emit and which must be adjusted first.
#' @param time_upper Optional numeric upper bound vector for censoring intervals.
#'   Used when `status %in% c(-1, 2)`; defaults to `time` if NULL.
#' @param x Optional design matrix (or data frame coercible to matrix).
#'   Unlike `time` and the other vector arguments, `x` is **not**
#'   data-masked: `hazard(data = df, time = tt, x = age)` errors with
#'   `object 'age' not found` where `time = tt` resolves, so write
#'   `x = df[["age"]]` or use the formula interface.
#' @param formula Optional formula with a `Surv()` object on the left.
#'   Right-censored (`Surv(time, status)`), left-censored
#'   (`Surv(time, event, type = "left")`), interval-censored
#'   (`type = "interval"` or `"interval2"`) and counting-process
#'   (`Surv(start, stop, event)`) forms are all accepted. `Surv()` codes
#'   censoring status with different integers than this package does; the
#'   formula path translates them, so write `Surv()`'s codes here. A plain
#'   `status` vector takes this package's codes; a `Surv` passed as `status`
#'   is translated the same way as here. Every distribution carries its own
#'   intercept: its baseline parameter, the one the covariates add to,
#'   already plays that role, so the design never has one: removing
#'   it (`~ 0 + age`, `~ age - 1`) is ignored with a warning, and builds the
#'   design of `~ age`.
#'   A `.` on the right-hand side stands for every column of `data` that the
#'   `Surv()` term does not use, as in `survival::coxph()`.
#'   When provided, overrides direct time/status/x arguments and extracts from data.
#'   Example: `hazard(Surv(time, status) ~ x1 + x2, data = df, dist = "weibull", fit = TRUE)`.
#' @param data Optional data frame. On the formula path it supplies the model
#'   frame, and `weights` is evaluated in its scope. On the vector path
#'   `time`, `status`, `time_lower`, `time_upper` and `weights` are evaluated
#'   in its scope, the way [base::subset()] and [base::transform()] do: a bare column name resolves to that column, and
#'   anything that is not a column (`df$col`, a local vector, a literal)
#'   falls through to the calling environment. A column of the same name as
#'   a caller variable wins, and because that silently discards the caller's
#'   vector (the way a wrapper forwarding its own argument by name does),
#'   such a name raises a warning naming the symbol and the argument.
#'   The warning reads the names written in the expression. When the
#'   argument is the name alone, it says the column was used. When the name
#'   is part of a larger expression, it does not say which value was read:
#'   the expression may never evaluate the name (an unused function
#'   argument), may rebind it first (a loop variable, an assignment) or may
#'   evaluate it somewhere else (`with()`), and the warning cannot tell. A name chosen at run time, as in `get(nm)` or
#'   `eval(as.name(nm))`, resolves the same way, column first, as it does in
#'   [stats::lm()], but is not checked; the vector path has behaved so since
#'   1.2.2. No second evaluation is made to compare the two values, because
#'   evaluating an expression such as `runif(n)` twice gives two different
#'   vectors.
#'   Masked arguments are validated like any other, so an `NA` in a
#'   masked column errors: an `NA` count on the SAS `ICENSOR`
#'   path reaches `weights` and stops with `'weights' must be
#'   non-negative and finite`.
#'   A named element of `data` that is a **function** is refused, on both
#'   paths, as [stats::lm()] refuses it. Because the mask sits in front of
#'   the calling frame, such an element would be called in place of the
#'   function an expression names -- `weights = rep(1, n)` calling a `rep`
#'   held in `data` -- and the fit would change with nothing to show for it.
#'   An S4 generic and a reference-class generator are functions for this
#'   purpose. A list-column of functions is a list, not a function, and is
#'   unaffected, as is an element with no name, which no expression can look
#'   up -- but an `NA_character_` name is refused, because R binds such an
#'   element under the symbol `` `NA` `` and a call reaches it. Remove the
#'   element: for a vector argument or `weights`, define the helper in the
#'   calling environment; for one used inside the `Surv()` response, compute
#'   the value into a `data` column first, since the response is evaluated
#'   without the formula's environment.
#' @param time_windows Optional numeric vector of strictly positive cut points for
#'   piecewise time-varying coefficients. When provided, each predictor column in
#'   `x` is expanded into one column per time window so each window gets its own
#'   coefficient.
#' @param theta Numeric starting values for optimization: the shape
#'   parameters, then one coefficient per column of the design matrix.
#'   Required when `fit = TRUE` with a single-distribution `dist`; optional for
#'   `dist = "multiphase"`, which assembles its own starting values from
#'   `phases` when `theta` is `NULL`.
#' @param dist Character baseline distribution label (default `"weibull"`).
#'   One of `"weibull"`, `"exponential"`, `"loglogistic"`, `"lognormal"`, or
#'   `"multiphase"`.  The single-distribution families differ in the *shape* the
#'   hazard traces over time; `"multiphase"` builds an additive N-phase hazard
#'   and requires `phases`.  See the **Baseline distributions** section for what
#'   each means and when to use it.
#' @param phases Optional named list of [hzr_phase()] objects specifying the
#'   phases for a multiphase model (`dist = "multiphase"`).  Names must be
#'   unique, and `"total"` and `"time"` are reserved: they label the summed
#'   column and the column of prediction times in the decomposed output of
#'   [predict.hazard()], so a phase of either name is rejected.  See Examples.
#' @param fit Logical; if TRUE, fit the model via maximum likelihood (default FALSE).
#' @param weights Optional numeric vector of observation weights (non-negative).
#'   Each observation's log-likelihood contribution is multiplied by its weight.
#'   Use for severity-weighted repeated events. Default `NULL` (unit weights).
#'   Implements the SAS `WEIGHT` statement. With `data`, a name is looked up
#'   among its columns first and then in the calling environment, on both
#'   interfaces. [stats::lm()] also looks in `data` first, but then in the
#'   formula's environment rather than the caller's, so a formula built
#'   inside another function does not bring that function's variables with
#'   it here. See `data` for the warning raised when a name is both.
#' @param control Named list of control options (see Details).
#' @param objective Which interval-censored contribution the multiphase
#'   likelihood accumulates. `"likelihood"` (default) uses the interval
#'   probability \eqn{\log(S(l) - S(u))}. `"sas"` reproduces what
#'   `PROC HAZARD` accumulates: the event-density term with the instantaneous
#'   hazard replaced by the interval-mean hazard over \eqn{(l, u]}. Applies
#'   only to `dist = "multiphase"`; exact-event and right-censored rows are
#'   unaffected either way.
#'
#'   `"sas"` requires data it can represent: no left-censored rows (`PROC
#'   HAZARD` has no left-censoring statement) and a positive width on every
#'   interval-censored row (the interval-mean hazard divides by \eqn{u - l}).
#'   Both are properties of the data rather than of the fit, so they are
#'   checked when the argument is supplied, including under `fit = FALSE`,
#'   which therefore stops rather than returning an unusable object.
#' @note `objective = "sas"` exists to reproduce legacy `PROC HAZARD` runs and
#'   **must not be used for new analyses**. It is a density, not a probability:
#'   it is inconsistent for wide intervals, where the two forms differ
#'   materially (22 log-likelihood units on the esophagectomy reference fit).
#'   The default is the statistically correct interval likelihood. See
#'   `inst/dev/SAS-INTERVAL-OBJECTIVE-DESIGN.md` for the derivation and the
#'   four-reference evidence.
#'
#' @param ... Additional named arguments retained for parity with legacy calling
#'   conventions.
#'
#' @details
#' Control parameters:
#' - `maxit`: Maximum iterations of the quasi-Newton (BFGS) optimizer
#'   (default 1000), applied to each start. A fit whose optimizer reports
#'   convergence but fails SAS's gradient test is continued with
#'   [stats::nlm()] under the same limit (see "Convergence"), so raising
#'   `maxit` lets that continuation run further too. The Nelder-Mead warm-up
#'   that a multiphase fit with fixed parameters may run first has its own
#'   limit, which `maxit` does not change.
#' - `n_starts`: Number of optimization starts for multiphase fits (default 5).
#'   Each start after the first offsets the initial values. The offsets are
#'   drawn from an internally seeded stream, so a multiphase fit is
#'   reproducible without `set.seed()` and does not advance the caller's RNG
#'   stream.
#' - `start_seed`: Seed selecting the ensemble of multi-start offsets
#'   (default 3). Any whole number within integer range, negative included.
#'   Fits are reproducible at any value; vary it to probe a different set of
#'   starts when a fit is suspected of sitting in a local optimum, and compare
#'   the resulting `objective` values. A fractional value is rejected rather
#'   than truncated, since `3.9` and `3` would otherwise select the same
#'   ensemble without saying so.
#' - `phase_share_tol`: Threshold for the multiphase identifiability warning
#'   (default 1e-8). A phase is reported as having left the model when it
#'   supplies less than this share of the cumulative hazard at every observed
#'   time (it never started; neither its `mu` nor its shape is identified),
#'   or when its contribution varies by less than this relative amount across
#'   them (it finished before the first observation and acts as a constant
#'   offset; `mu` stays identified, the shape parameters do not). A third
#'   condition is a property of the observed times rather than of any phase:
#'   when their own relative range falls below this threshold and no phase's
#'   contribution varies across them, they separate no phase from any other
#'   and none of the shapes is identified. That verdict is withheld when the
#'   fit supplies counting-process entry times or interval bounds that add
#'   enough evaluation points the measure did not see. "Enough" is counted, not
#'   assumed: with the event times tied, each added point supplies one more
#'   functional of the parameters, and one more than the number of added points
#'   must reach the number of free parameters. A zero entry time, or a bound
#'   equal to an event time, adds nothing and does not count. The count is
#'   deliberately conservative: it does not credit the extra functional that a
#'   zero entry time makes separately observable, so a fit can be warned about
#'   while being identified. For the same reason the per-phase measures can
#'   overstate what an interval-censored or left-truncated fit loses, since a
#'   phase flat across the event times may still be identified through the
#'   bounds; the saturated warning says so where such points exist, rather
#'   than claiming the likelihood is unchanged (#228). The measured shares are
#'   kept on the fit as `fit$fit$phase_share`. Raise it to catch marginal
#'   phases, set it to 0 to silence the check.
#' - `reltol`: Relative convergence tolerance on the objective, the negative
#'   log-likelihood (default 1e-5). BFGS stops when an iteration reduces it by
#'   less than `reltol * (|objective| + reltol)`, so the stopping gap grows
#'   with the size of the log-likelihood: about 0.0024 at a log-likelihood of
#'   -240. On a flat surface plain BFGS can stop that far short of the optimum
#'   and still report convergence, so `hazard()` then applies SAS/C HAZARD's
#'   relative-gradient test and, when the stop fails it, continues with
#'   [stats::nlm()]; see the "Convergence" section.
#' - `conserve`: Apply Conservation of Events (**`dist = "multiphase"` only**;
#'   default `TRUE`). CoE counts exact events, so it is **automatically
#'   disabled** whenever any `status` falls outside \{0, 1\} (which interval
#'   or left censoring guarantees) and whenever the model has fewer than two
#'   phases. On a fitted multiphase object the outcome is recorded next to the
#'   request, both fields living under `fit$spec$control`:
#'   - `fit$spec$control$conserve_applied`: logical, whether CoE was actually
#'     applied;
#'   - `fit$spec$control$conserve_disabled_reason`: one of
#'     `"not_requested"`, `"unsupported_censoring"`, `"single_phase"`,
#'     `"no_events"`, `"setup_failed"`, or `NA` when CoE was applied.
#'
#'   Read `fit$spec$control$conserve_applied`, not
#'   `fit$spec$control$conserve`: the latter says only what you asked for.
#' - `shape_param_count`: The number of baseline parameters at the front of
#'   `theta`: the scale and any shape parameters, so 2 for `"weibull"` and 1
#'   for `"exponential"`. For a single-distribution model only. The fit itself does not
#'   use it; the stepwise refit and the score test read it back from
#'   `fit$spec$control`. A multiphase fit derives its own layout, so nothing
#'   reads it there.
#'
#' The elements above are accepted without a warning: `maxit` and `reltol`
#' for every model, `shape_param_count` for a single-distribution model, and
#' `n_starts`, `start_seed`, `phase_share_tol` and `conserve` for
#' `dist = "multiphase"` (#376). Any other element draws one warning that
#' names it and says why it has no effect, and the fit proceeds unchanged,
#' as [stats::optim()] does for unknown `control` names. The element is
#' dropped before the fit, so a name such as `n_starts_extra` cannot be read
#' as `n_starts`, and `fit$spec$control` keeps none of the ignored
#' elements. That covers:
#' - a name no fit reads, such as the misspelling `n_startz`, and an unnamed
#'   element;
#' - `abstol` (read only by a bounded optimizer no fit uses), `method`,
#'   `condition`, `nocov` and `nocor`, which earlier versions documented as
#'   accepted without reading them;
#' - `fix` and `quasi`, which no fit has ever read. A fit given `fix` was
#'   never constrained, so results obtained with it may be affected; hold a
#'   parameter with `hzr_phase(fixed = )` on a multiphase phase. A
#'   single-distribution model has no mechanism for fixing a parameter.
#' - a multiphase element such as `n_starts` given to a single-distribution
#'   fit, and `shape_param_count` given to a multiphase one.
#'
#' No name in `control` is an error; a bad value for an element the fit
#' reads, such as `maxit = "a"`, still stops a fit (`fit = TRUE`) where it
#' is read. [hzr_stepwise()] and
#' [hzr_bootstrap()] pass `control` to every candidate refit, and an error
#' there would count as a failed candidate, so a screen would report success
#' having tested nothing.
#'
#' SAS `PROC HAZARD` options with no `control` equivalent: `NOCOV` and `NOCOR`
#' only suppress printed output, and `hazard()` prints nothing until asked.
#' `CONDITION=` stops SAS's optimizer on a condition-number test that
#' `hazard()` does not have. `QUASI` is `hazard()`'s optimizer already: the
#' fits it runs use BFGS (a multiphase fit may run a Nelder-Mead warm-up
#' first, and a stop that fails SAS's gradient test continues with
#' [stats::nlm()]). SAS jobs often write `STEEPEST QUASI` together, steepest
#' descent first; **there is no steepest-descent option and no two-stage
#' strategy**. The multiphase likelihood is multimodal, so a different descent
#' path can land on a different optimum: a fit translated from a job using
#' `STEEPEST` may not reproduce SAS's estimates, and `hzr_translate_sas()`
#' records the keyword as untranslated rather than dropping it.
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
#' \code{vignette("fitting-hazard-models")}: single-phase through multiphase fitting,
#' \code{vignette("prediction-visualization")}: prediction types and decomposed hazard plots,
#' \code{vignette("inference-diagnostics")}: bootstrap CIs and model diagnostics.
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
#'   \code{converged}, \code{se}, \code{vcov}, \code{counts}, \code{message},
#'   and \code{rel_gradient}, \code{rel_gradient_reason} and
#'   \code{polish_code}, the SAS/C acceptance
#'   test described under "Convergence";
#'   all \code{NULL} when \code{fit = FALSE}; multiphase fits add
#'   \code{starts}, one row per optimisation start with its \code{status}
#'   (\code{"ok"}, \code{"nonconverged"}, \code{"infeasible"},
#'   \code{"nonfinite"} or \code{"error"}), \code{objective} (\code{NA}
#'   unless the start reached a point where the likelihood is defined, and
#'   recorded as the optimizer returned it, before any Conservation of Events
#'   adjustment to the conserved phase's scale), \code{convergence} (the
#'   \code{\link[stats]{optim}} code, \code{0} for success), whether it was
#'   the \code{best} and so the reported fit, and the \code{message} of any
#'   error. A start that stops at \code{maxit} has a finite \code{objective}
#'   and can win the selection, so \code{status} distinguishes it from one
#'   that converged. Every fit carries \code{weak}, which takes one of three
#'   values: a list, when the likelihood is near-flat along a parameter
#'   combination, giving the \code{params} spanning that direction, their
#'   squared loadings (\code{weights}), the strongest pairwise
#'   \code{correlation} among them, the Hessian \code{rcond} and
#'   \code{n_directions}, the number of near-flat directions found (when a
#'   single parameter carries the direction on its own, the list has one
#'   \code{params} entry, \code{single = TRUE}, its \code{estimate}, and
#'   \code{correlation = NA}, and \code{se_metric}, its standard error on
#'   the log scale; only a g3 shape, \code{gamma}, \code{alpha} or
#'   \code{eta}, is named this way);
#'   \code{NULL} when the fit was examined and is well identified; and
#'   \code{NA} when the check could not run because no usable Hessian was
#'   available, which includes an unfitted object and an install without
#'   the suggested \pkg{numDeriv}. Test with \code{is.list(fit$fit$weak)},
#'   not \code{!is.null()}: the \code{NA} case has not been examined and
#'   must not be read as a clean result),
#'   \code{engine} (implementation tag, \code{"native-r-m2"}), and two
#'   fields recording what the fit did not do: \code{degraded}, a character
#'   vector of the steps not performed, in the fixed order
#'   \code{"fitting"}, \code{"standard_errors"},
#'   \code{"conserved_phase_variance"}, \code{"weak_direction_check"},
#'   \code{"boundary_check"}, \code{"conservation_of_events"}, and empty
#'   when nothing was lost; and
#'   \code{degraded_causes}, a character vector with the same names giving
#'   the reason for each. \code{print()} and \code{summary()} always show
#'   them as a "Not done in this run" block, which reads "none" when nothing
#'   was lost. \code{fit$fit$weak} is \code{NA} exactly when
#'   \code{"weak_direction_check"} is listed.
#'
#'   \code{fit$fit$boundary} is its sibling and takes the same three states,
#'   for phases whose fit sits where their parameterisation breaks down:
#'   \code{NULL} when the check ran and found nothing, a list of records when
#'   it found something, and \code{NA} when it did not run --- and it is
#'   \code{NA} exactly when \code{"boundary_check"} is listed in
#'   \code{degraded}, with the reason in \code{degraded_causes}. Each record
#'   carries \code{mechanism}, \code{phase}, \code{parameter} and a
#'   printable \code{detail}. The mechanisms are \code{"unbounded_phase"}, a
#'   \code{"hazard"} phase whose \code{t_half} is below the first observed
#'   time; \code{"phase_discontinuity"}, a \code{"cdf"} or \code{"hazard"}
#'   phase whose shape has collapsed to a step the observed times cannot
#'   resolve (it can lie inside the data); and two made at setup, before the
#'   fit: \code{"g3_alpha_one"} (a g3 phase with \code{alpha} fixed at 1,
#'   re-expressed as PROC HAZARD does, see [hzr_phase()]) and
#'   \code{"g3_fixge2_alpha_start"} (a free \code{alpha} start moved to 2/3
#'   under \code{constraint = "eta_gamma"}). Only rows the likelihood reads
#'   count as observed times. A fit with any record raises one warning whose
#'   classes are \code{"hzr_"} plus each mechanism present, all inheriting
#'   \code{"hzr_boundary"}, so one handler catches the whole family.
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
                   objective = c("likelihood", "sas"),
                   ...) {

  objective <- match.arg(objective)
  # A named scalar such as c(model = "multiphase") is a valid `dist`, but
  # identical() against a bare string is FALSE for it. Dropping the names
  # here, before any read, keeps every later test of `dist` agreeing (#405).
  dist <- unname(dist)
  # The caller's own `data` expression, captured before `data` is reassigned:
  # the ambiguity warning names it in its advice (#401).
  data_arg <- substitute(data)

  # `objective` is a top-level argument rather than a `control` element on
  # purpose: it changes the estimand, and burying that among convergence
  # tolerances makes it easy to miss in review.
  if (objective == "sas" && dist != "multiphase") {
    stop("objective = \"sas\" applies only to dist = \"multiphase\": it ",
         "reproduces PROC HAZARD's interval-censored contribution, and no ",
         "other distribution here is a PROC HAZARD target. Got dist = \"",
         dist, "\".", call. = FALSE)
  }
  x_design <- NULL
  # Formula dispatch: if formula is provided, parse it and extract time/status/x from data
  # `data` masks the calling frame for the argument expressions AND for the
  # formula's Surv() response, and R's function lookup walks past every
  # binding that is not a function. So a function-valued element is CALLED in
  # place of the function an expression names, changing the fit with nothing
  # to show for it (#420). stats::lm() refuses the same shape, less clearly
  # ("cannot coerce class '\"function\"' to a data.frame").
  #
  # This runs BEFORE .hzr_numeric_frame_values() and before the formula/vector
  # branch, both deliberately. That helper replicates each column to `nrow`
  # and dies on a function while doing it, which looks like a guard and is
  # not one: at nrow == 1 there is nothing to replicate, the function
  # survives, and the formula path read it. Running first also means the
  # named refusal is what the user sees, never "attempt to replicate an
  # object of type 'closure'".
  #
  # Only elements a call can REACH. An unnamed or ""-named element cannot be
  # looked up at all, so refusing it would be a false refusal. An
  # `NA_character_` name is NOT in that class, however it reads: R binds such
  # an element under the symbol `NA`, and `` `NA`(x) `` calls it, so
  # exempting it left the whole defect open through one spelling (#443
  # review). Data frames are not exempted -- `data.frame()`, `$<-` and `[[<-`
  # each refuse a function column, but `structure(list(...), class =
  # "data.frame")` carries one and `is.data.frame()` is TRUE for it. A
  # list-column is a list, which the lookup skips, so it still fits.
  # `is.list()` first, and not merely for speed: this guard iterates `data`,
  # so on anything that is not list-like it would answer BEFORE the shape
  # check below and answer wrongly -- `vapply()` raises a coercion error for
  # an S4 object, and for an environment it reports a "function element" and
  # tells the user to remove it, which does not make an environment
  # acceptable `data`. Both are questions about shape, not about functions.
  # `is.list()` is TRUE for a data frame, tibble and data.table alike.
  if (is.list(data)) {
    nm <- names(data)
    # An unnamed list has `nm` NULL, and `is.na(NULL)` is already logical(0),
    # which zeroes the whole vector -- so NULL names need no test of their own.
    fn <- vapply(data, is.function, logical(1)) & (is.na(nm) | nzchar(nm))
    if (any(fn)) {
      stop("'data' holds a function named ",
           paste0("'", unique(nm[fn]), "'", collapse = ", "),
           ". `data` masks the calling frame while hazard() evaluates its ",
           "arguments and the formula's response, so such an element is ",
           "called in place of the function the expression names, changing ",
           "the fit with nothing to show for it. Remove it from 'data'. ",
           "For a vector argument, or the formula's 'weights', define the ",
           "helper in the calling environment instead. For a helper used ",
           "inside the formula's Surv() response, compute the value into a ",
           "'data' column first: the response is evaluated without the ",
           "formula's environment, so a helper defined there is not visible ",
           "to it.", call. = FALSE)
    }
  }

  # Columns of `data` are read before any argument is: Surv() and
  # model.matrix() take a classed numeric's stored doubles too (#231).
  data <- .hzr_numeric_frame_values(data)

  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("'data' is required when 'formula' is provided.", call. = FALSE)
    }

    # A multiphase formula RHS may name a phase as a function, as in
    # `constant(age)`. Such a term is refused, not routed (#275): it used to
    # be stripped with the rest of the RHS and never reached its phase, so
    # the fit converged without it. A phase's covariates belong in
    # hzr_phase(formula = ...). The parse-tree walk (not a string regex)
    # keeps a base-R call such as `log(age)` from matching.
    if (identical(dist, "multiphase") && !is.null(phases) &&
          length(formula) >= 3L) {
      scoped <- Filter(
        function(nm) .hzr_formula_has_phase_scope(formula[[3L]], nm),
        names(phases)
      )
      if (length(scoped) > 0L) {
        stop("The formula names ",
             if (length(scoped) > 1L) "phases " else "phase ",
             paste0("'", scoped, "'", collapse = ", "),
             " as a function, as in `", scoped[[1L]], "(var)`. hazard() ",
             "does not route such terms to a phase, so they would be ",
             "dropped. Give the phase its covariates with ",
             "hzr_phase(..., formula = ~ var) instead, and leave them out ",
             "of the formula.", call. = FALSE)
      }
    }

    parsed <- .hzr_parse_formula(formula = formula, data = data)
    time <- parsed$time
    status <- parsed$status
    time_lower <- parsed$time_lower
    time_upper <- parsed$time_upper
    x <- parsed$x
    x_design <- parsed$x_design

    # `weights` is looked up in `data` first, then the calling frame, by the
    # rule the vector path below applies (#392). This path used to skip
    # `data`: a column-only name was not found, and a name bound both as a
    # column and in the calling frame silently read the calling frame's
    # vector. stats::lm() also looks in `data` first, but falls back to the
    # formula's environment, not the calling frame.
    .hzr_warn_masked_ambiguity(list(weights = substitute(weights)), data,
                               parent.frame(), interface = "formula",
                               data_arg = data_arg)
    weights <- eval(substitute(weights), data, parent.frame())
  }

  # Data masking on the vector path. `data` used to be consulted only by the
  # formula path, so hazard(data = df, time = tt) failed with "object 'tt'
  # not found" while looking like it should work -- the defect behind the SAS
  # translator's unrenderable documents (#151). This is the base-R idiom
  # subset()/transform()/with() use: evaluate the argument expression with
  # `data` as the environment and the caller's frame as its parent, so a
  # column wins, and anything that is not a column (df$col, a local vector,
  # a literal) falls through to the caller unchanged.
  if (is.null(formula) && !is.null(data)) {
    if (!is.data.frame(data) && !is.list(data)) {
      stop("'data' must be a data frame or a list.", call. = FALSE)
    }
    mask_env <- parent.frame()
    .hzr_warn_masked_ambiguity(
      list(time = substitute(time), status = substitute(status),
           time_lower = substitute(time_lower),
           time_upper = substitute(time_upper),
           weights = substitute(weights)),
      data, mask_env, data_arg = data_arg
    )
    time <- eval(substitute(time), data, mask_env)
    status <- eval(substitute(status), data, mask_env)
    time_lower <- eval(substitute(time_lower), data, mask_env)
    time_upper <- eval(substitute(time_upper), data, mask_env)
    weights <- eval(substitute(weights), data, mask_env)
  }

  # After formula dispatch, require time and status
  if (is.null(time) || is.null(status)) {
    stop("'time' and 'status' are required (either directly or via 'formula').", call. = FALSE)
  }
  # See .hzr_numeric_values(): a classed numeric's stored doubles are not its
  # values, and the likelihoods read them raw (#231).
  time <- .hzr_numeric_values(time)
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
  }

  n <- length(time)
  if (length(status) != n) {
    stop("'status' must have the same length as 'time'.", call. = FALSE)
  }
  # Over zero rows every path returned an object, fitted or not, with nothing
  # behind it; a fit even reported converged = TRUE (#231).
  if (n == 0L) {
    stop("hazard() was given no observations: 'time' has length 0. ",
         "Check that `data` (or the subset passed to it) has rows.",
         call. = FALSE)
  }

  # A Surv object passed as `status` is read exactly as the formula path reads
  # it (#226). Its codes are not this package's, and under "interval" and
  # "counting" its second column is not the status at all, so taking that
  # column unchanged fitted left-censored rows as right-censored. The Surv
  # defines the bounds it carries; `time` and any bound the caller also gave
  # must agree with it rather than be silently overridden.
  if (inherits(status, "Surv")) {
    resp <- .hzr_surv_response(status)
    if (!identical(as.numeric(time), as.numeric(resp$time))) {
      stop("'time' does not match the times in the Surv object passed as ",
           "'status'. Pass time = unclass(status)[, 1] -- the stop column, ",
           "[, 2], for Surv(start, stop, event) -- or use the formula ",
           "interface.", call. = FALSE)
    }
    surv_bound <- function(given, from_surv, arg) {
      if (is.null(from_surv)) {
        return(given)
      }
      if (!is.null(given) &&
            !identical(as.numeric(given), as.numeric(from_surv))) {
        stop("'", arg, "' does not match the Surv object passed as ",
             "'status'. Omit it and the bound is taken from the Surv.",
             call. = FALSE)
      }
      from_surv
    }
    time_lower <- surv_bound(time_lower, resp$time_lower, "time_lower")
    time_upper <- surv_bound(time_upper, resp$time_upper, "time_upper")
    status <- resp$status
  }

  # Optional censoring bounds:
  # - status = 1 (event): uses observed event time in `time`
  # - status = 0 (right-censored): censoring time in `time` (or `time_lower`)
  # - status = -1 (left-censored): upper bound in `time` (or `time_upper`)
  # - status = 2 (interval-censored): [time_lower, time_upper] required
  if (!is.null(time_lower)) {
    time_lower <- .hzr_numeric_values(time_lower)
    if (!is.numeric(time_lower) || length(time_lower) != n || any(!is.finite(time_lower)) || any(time_lower < 0)) {
      stop("'time_lower' must be a numeric vector of finite non-negative values matching length(time).", call. = FALSE)
    }
  }

  if (!is.null(time_upper)) {
    time_upper <- .hzr_numeric_values(time_upper)
    if (!is.numeric(time_upper) || length(time_upper) != n || any(!is.finite(time_upper)) || any(time_upper < 0)) {
      stop("'time_upper' must be a numeric vector of finite non-negative values matching length(time).", call. = FALSE)
    }
  }

  # A row at time 0 is deleted, whatever its status, as PROC HAZARD deletes
  # it at input (#374): hazard/src/hazard/readt.c:12-14 marks TIME <= 0 as
  # inadmissible, readobs.c:132-138 leaves it out of the data, and
  # obsstat.c:62-69 notes the count. The rule is on SAS's TIME alone; a
  # lower bound or entry time of 0 is admissible, as readct.c allows CT = 0.
  # Fitting such a row instead made an EVENT at 0 return the optimizer's
  # clamp or log(double.xmax) with converged = TRUE (lognormal, loglogistic,
  # exponential), a value that is not a likelihood.
  # SAS's TIME is the row's UPPER bound. Here that is `time` for an exact or
  # right-censored row, but `time_upper` for a left- or interval-censored one:
  # an interval row's `time` is its LOWER bound (Surv "interval2" maps
  # [l, u] to time = l, time_upper = u), so testing `time` there would drop
  # a legitimate interval opening at 0, which #341 fits as left censoring.
  sas_time <- time
  if (!is.null(time_upper) && length(time_upper) == n) {
    bounded <- !is.na(status) & status %in% c(-1, 2)
    sas_time[bounded] <- time_upper[bounded]
  }
  at_zero <- sas_time == 0
  n_dropped_time_zero <- sum(at_zero)
  dropped_time_zero_rows <- which(at_zero)
  data_model <- NULL
  dropped_frame <- NULL
  if (n_dropped_time_zero > 0L) {
    keep <- !at_zero
    subset_rows <- function(v) {
      if (is.null(v)) return(v)
      if (is.matrix(v) || is.data.frame(v)) {
        if (nrow(v) == n) v[keep, , drop = FALSE] else v
      } else if (length(v) == n) {
        v[keep]
      } else {
        v
      }
    }
    time <- time[keep]
    status <- status[keep]
    time_lower <- subset_rows(time_lower)
    time_upper <- subset_rows(time_upper)
    weights <- subset_rows(weights)
    x <- subset_rows(x)
    if (is.list(data)) {
      data_full <- data
      # A list of columns is accepted as `data` too; subset each column of
      # the row count, as a data frame's rows are.
      data <- if (is.data.frame(data)) {
        subset_rows(data)
      } else {
        lapply(data, subset_rows)
      }
      # Every expression -- the response, scale(age), a factor's levels, a
      # phase formula's terms -- is computed on the data AS GIVEN, and the
      # rows are removed afterwards, as PROC HAZARD removes them after the
      # DATA step and as stats::model.frame() does for `subset =` and
      # `na.action`. (John, 2026-09-25: the alternative, rebuilding without
      # the dropped rows, is circular for a response such as
      # Surv(time - min(time), status), whose zeros depend on those rows.)
      # The global design above was built that way already. Phase designs
      # are built later from `data`, so the optimizer gets the full rows
      # with the kept-row mask, and subsets each phase design after building
      # it (.hzr_multiphase_designs()).
      data_model <- data_full
      attr(data_model, "hzr_rows_kept") <- keep
      if (is.data.frame(data_full)) {
        dropped_frame <- data_full[!keep, , drop = FALSE]
      }
    }
    n <- length(time)
    warning(structure(
      class = c("hzr_time_zero_dropped", "warning", "condition"),
      list(message = paste0(
        n_dropped_time_zero, " row(s) with time = 0 were dropped before ",
        "fitting, as PROC HAZARD drops them (TIME <= 0 is inadmissible, ",
        "readt.c). The fit uses the other ", n, "; the count is in ",
        "fit$data$dropped_time_zero, and row numbers in later messages ",
        "count the rows that remain."
      ), call = NULL)
    ))
  }
  # For status 0/1 rows `time_lower` is the counting-process ENTRY time when
  # 0 < time_lower < time, and every family forms H(time) - H(time_lower)
  # there (the entry rule lives in each likelihood, e.g.
  # .hzr_multiphase_entry()).  time_lower == time means no entry: that is the
  # mixed-interval layout, where only status-2 rows carry a real lower bound.
  # Multiphase used to read it as an entry at exit, and the objective became
  # unbounded (+47915.76, converged = TRUE, issue #136/#253).
  # An entry AFTER exit is a data error.  SAS HAZARD refuses STIME >= TIME
  # (SETCOE960 in setcoe_obs_loop.c), but its STIME is only ever an entry
  # time; here time_lower doubles as the interval bound, so only
  # time_lower > time is refused (issue #253).
  if (!is.null(time_lower)) {
    bad_entry <- status %in% c(0, 1) & time_lower > time
    if (any(bad_entry)) {
      stop(sum(bad_entry), " of ", n, " row(s) with status 0 or 1 have ",
           "'time_lower' > 'time'. On these rows 'time_lower' is the ",
           "counting-process entry time, and a subject cannot enter the ",
           "risk set after it leaves (SAS HAZARD rejects this too). For no ",
           "entry time, leave 'time_lower' as NULL or set it to 0 or to ",
           "'time'.",
           call. = FALSE)
    }
    # time_lower == time > 0 reads as "no entry" only in the mixed-interval
    # layout. Beside genuine entries it is a zero-length counting-process
    # epoch (hzr_repeated_events() emits them), and "no entry" would charge
    # the row its full H(0, time], so the mix is refused (r-reviewer, #253).
    genuine  <- status %in% c(0, 1) & time_lower > 0 & time_lower < time
    zero_len <- status %in% c(0, 1) & time_lower > 0 & time_lower == time
    if (any(genuine) && any(zero_len)) {
      stop(sum(zero_len), " of ", n, " row(s) with status 0 or 1 have ",
           "'time_lower' equal to 'time' while other rows carry genuine ",
           "entry times. These are zero-length counting-process epochs, ",
           "entering and leaving the risk set at the same moment; remove ",
           "or adjust them before fitting (see ?hzr_repeated_events). SAS ",
           "HAZARD refuses them too.",
           call. = FALSE)
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
    # Absent names pass the check on `x`, but the expansion names each window
    # column <name>_w<k>, so two of them become "_w1" (or "NA_w1") twice.
    .hzr_refuse_duplicate_columns(x_fit)
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
    weights <- .hzr_numeric_values(weights)
    if (!is.numeric(weights) || length(weights) != n_obs) {
      stop("'weights' must be a numeric vector of length ", n_obs, ".",
           call. = FALSE)
    }
    if (any(weights < 0) || any(!is.finite(weights))) {
      stop("'weights' must be non-negative and finite.", call. = FALSE)
    }
  }
  # Every likelihood branches on these four codes, and a row with any other
  # code falls through all of them and adds nothing: survival's interval
  # code 3, passed as a plain vector, was silently dropped (#231). NA is left
  # to the completeness check, which names the rows.
  # %in% compares as text, so a character or factor status passed that
  # check, and the single-distribution likelihoods then returned their
  # starting values as a converged fit.
  if (!is.numeric(status) && !is.logical(status)) {
    stop("'status' must be numeric (or logical), not ", class(status)[1L],
         ". Convert it, for example with as.numeric(as.character(status)) ",
         "for a factor.", call. = FALSE)
  }
  # After any Surv translation above, so a Surv is never flattened here.
  status <- .hzr_numeric_values(status)
  bad_status <- !is.na(status) & !(status %in% c(-1, 0, 1, 2))
  if (any(bad_status)) {
    stop("'status' must be coded -1 (left-censored), 0 (right-censored), ",
         "1 (event) or 2 (interval-censored); ", sum(bad_status), " of ", n,
         " row(s) are not, at index/indices ",
         paste(utils::head(which(bad_status), 10L), collapse = ", "),
         if (sum(bad_status) > 10L) ", ..." else "", ". A Surv object's ",
         "codes differ from these: pass it as the response, or as 'status', ",
         "and it is translated.", call. = FALSE)
  }
  # A row adds nothing to the likelihood when its weight is 0, or when it is
  # right-censored at time 0 (H(0) = 0). With no other row the fit returned
  # its starting values, objective 0 and converged = TRUE, as zero rows did.
  contributes <- !(status == 0 & time == 0)
  if (!is.null(weights)) contributes <- contributes & weights > 0
  if (!anyNA(status) && !any(contributes)) {
    stop("hazard() was given no observations that contribute to the ",
         "likelihood: every row has weight 0",
         if (n_dropped_time_zero > 0L) {
           paste0(" or time = 0 (", n_dropped_time_zero, " such row(s) ",
                  "dropped, as PROC HAZARD drops them)")
         },
         ".", call. = FALSE)
  }

  if (!is.character(dist) || length(dist) != 1 || !nzchar(dist)) {
    stop("'dist' must be a non-empty character scalar.", call. = FALSE)
  }

  # Only the multiphase optimizer assembles its own starting values. The
  # single-distribution branch below runs only when theta is supplied, so
  # without this a fit = TRUE call returned an unfitted object -- NULL
  # coefficients, NA objective -- with no error and no warning. `fit` is
  # tested exactly as that dispatch tests it, so no value reaches one and
  # not the other. A zero-length theta carries no starting values either; it
  # used to reach optim() and fail there with a message that did not say so.
  if (fit && dist != "multiphase" && length(theta) == 0L) {
    stop("fit = TRUE needs starting values for dist = \"", dist, "\": ",
         "supply 'theta' (the shape parameters, then one coefficient per ",
         "column of the design matrix), or use fit = FALSE for an unfitted ",
         "model. Only dist = \"multiphase\" assembles its own starting ",
         "values.", call. = FALSE)
  }

  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }
  # Only the names this fit reads go on: consumers read control with `$`,
  # which would partial-match a warned name such as n_starts_extra (#405).
  control <- .hzr_validate_control(control, dist)

  # A single-distribution theta must have one entry per parameter and, for
  # Weibull, a positive scale and shape, fitted or not: an unfitted object
  # of the wrong length can never be predicted from. Here, because x_fit is
  # final only after time-window expansion. Multiphase is checked below
  # (#408), where fit = FALSE may legitimately carry fewer entries.
  if (!is.null(theta) && dist != "multiphase") {
    .hzr_check_theta(theta, dist,
                     n_coef = if (is.null(x_fit)) 0L else ncol(x_fit),
                     windowed = !is.null(time_windows))
  }

  # Multiphase validation
  if (dist == "multiphase") {
    if (is.null(phases)) {
      stop("'phases' is required when dist = 'multiphase'. ",
           "Supply a list of hzr_phase() specifications.", call. = FALSE)
    }
    phases <- .hzr_validate_phases(phases)
    # `.` in a phase formula is written out here, once, before anything reads
    # it (#277). The likelihood, the score test, predict() and the stored
    # spec all build the phase design with model.frame(ph$formula, data),
    # which would expand `.` to every column, response included. It is
    # expanded as the global formula's is, against `data` without the Surv()
    # variables. The vector interface has no Surv() term to name those
    # columns, so there `.` is refused.
    for (nm in names(phases)) {
      pf <- phases[[nm]]$formula
      if (!is.null(pf)) {
        .hzr_refuse_offset(pf, paste0("the formula of phase '", nm, "'"))
      }
      if (is.null(pf) || !"." %in% all.vars(pf)) next
      if (is.null(formula)) {
        stop("Phase '", nm, "' uses `.` in its formula, which needs the ",
             "formula interface: with `time =` and `status =`, hazard() ",
             "cannot tell which columns of `data` hold the response. Write ",
             "the phase's terms out, or use hazard(Surv(...) ~ ..., ",
             "data = ...).", call. = FALSE)
      }
      two_sided <- stats::as.formula(
        call("~", formula[[2L]], pf[[length(pf)]]), env = environment(pf)
      )
      phases[[nm]]$formula <- .hzr_expand_rhs(
        two_sided,
        .hzr_drop_empty_names(data, pf, paste0(" in phase '", nm, "'"))
      )
    }
    .hzr_check_phase_formula_data(phases, data, x)
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

  # Data preconditions of objective = "sas", checked once here rather than
  # inside the objective: both are pure functions of the data, so reporting
  # them through the optimizer's per-start tryCatch framed a data defect as a
  # convergence problem. The guards inside the objective and gradient stay --
  # the gradient is reachable without hazard(). See .hzr_check_sas_data().
  .hzr_check_sas_data(status, time, time_lower, time_upper, objective)

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
  # Under fit = TRUE the optimizer derives a constrained shape from the rest
  # of theta. Unfitted, nothing would, and predict() would evaluate a model
  # off its own constraint. The slots are located the way the optimizer
  # locates them (.hzr_optim_multiphase()): a phase formula against `data`,
  # else the global design, else no covariates (#328).
  # A fit needs one theta entry per parameter. The check near the top only
  # compares the length with the global design, and only as a lower bound,
  # so a longer multiphase theta fitted with the extra entries carried along
  # and a shorter one failed inside the fit on a names() mismatch (#408).
  # Unfitted, theta is returned as supplied, so the constraint block below
  # warns instead.
  if (fit && dist == "multiphase" && !is.null(theta)) {
    per_phase <- .hzr_phase_theta_counts(phases, data, x_fit)
    if (length(theta) != sum(per_phase)) {
      stop(.hzr_theta_length_message(length(theta), per_phase),
           call. = FALSE)
    }
  }
  if (!fit && dist == "multiphase" && !is.null(theta) &&
      any(vapply(phases, function(ph) .hzr_phase_constraint(ph) != "none",
                 logical(1)))) {
    counts <- .hzr_phase_covariate_counts(phases, data, x_fit)
    n_theta <- sum(.hzr_phase_theta_counts(phases, data, x_fit))
    if (length(theta) == n_theta) {
      theta <- .hzr_constrain_supplied_theta(theta, phases, counts)
    } else {
      warning("theta has ", length(theta), " entries but these phases take ",
              n_theta, ", so hzr_phase(constraint = ) could not be applied ",
              "to it; the unfitted object carries theta as supplied.",
              call. = FALSE)
    }
  }

  fit_state <- list(
    theta = theta,
    converged = NA,
    objective = NA_real_,
    se = NULL,
    gradient = NULL,
    # NA rather than NULL, but UNOBSERVABLE TODAY and not a guard: the
    # post-fit block below always overwrites this, on every path including
    # fit = FALSE, so a mutation to NULL changes nothing and no test can
    # catch it. Kept so that an early return added later yields "never
    # examined" instead of an absent field -- stated plainly rather than
    # dressed up as a protection the code does not have (#444).
    boundary = NA
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
    # Scalar NA means there is no variance matrix at all -- the same contract
    # vcov.hazard() keeps. Whenever a matrix *is* present, size the result to it
    # so `$se` stays conformable with `$par`: a multiphase vcov legitimately
    # carries NA rows for parameters held fixed, and collapsing that to a
    # length-1 NA left callers unable to name the SEs against the parameters.
    if (is.null(vcov_mat) || !is.matrix(vcov_mat)) {
      return(NA)
    }
    d <- diag(vcov_mat)
    # Element by element, never all-or-nothing. A parameter held fixed carries
    # an NA variance row and earns an NA standard error; every parameter whose
    # variance *was* computed keeps its own. Collapsing the whole vector on the
    # first NA -- or on one small negative variance from a numerical Hessian --
    # discarded standard errors that summary() reports from the same matrix,
    # leaving `$se` the right length and empty of everything it should carry.
    out <- rep(NA_real_, length(d))
    # rep() starts unnamed; carry over whatever diag() gave us so a named vcov
    # still yields a named `$se`, as the old sqrt(d) path did.
    names(out) <- names(d)
    ok <- is.finite(d) & d >= 0
    out[ok] <- sqrt(d[ok])
    out
  }

  # Filled by whichever optimizer branch runs, and read by the record of what
  # this fit did not do (#242). fit_ran is FALSE exactly when fit = FALSE;
  # a single-distribution fit = TRUE call with no theta now stops (#243)
  # rather than reaching this point unfitted.
  fit_ran <- FALSE
  degraded_reasons <- list()

  # Distribution dispatch -- select the distribution-specific optimizer and fit.
  if (fit && dist == "multiphase") {
    # Multiphase: theta_start is optional (assembled from phase specs if NULL)
    optim_result <- .hzr_run_fit_safely(.hzr_optim_multiphase(
      time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x_fit, theta_start = theta, weights = weights,
      control = control,
      phases = phases, objective = objective,
      # The full rows with their kept-row mask when time-0 rows were
      # dropped, so phase designs are computed on the data as given.
      formula_global = formula, data = data_model %||% data
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
    fit_state$phase_share <- optim_result$phase_share
    fit_state$phases <- optim_result$phases
    fit_state$covariate_counts <- optim_result$covariate_counts
    fit_state$x_list <- optim_result$x_list
    fit_state$x_design <- optim_result$x_design
    fit_state$rows_used <- optim_result$rows_used
    fit_state$fixed_mask <- optim_result$fixed_mask
    fit_state$held <- optim_result$held
    fit_state$starts <- optim_result$starts
    # Applied CoE state, recorded next to the requested one in spec$control
    # below. Kept here first so the assembly reads the optimizer's answer
    # rather than re-deriving it from `status`.
    control$conserve_applied <- optim_result$conserve_applied
    control$conserve_disabled_reason <- optim_result$conserve_disabled_reason
    fit_ran <- TRUE
    degraded_reasons$se <- optim_result$se_unavailable_reason
    degraded_reasons$conserved_variance <- optim_result$conserved_variance_reason

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
    fit_ran <- TRUE
    degraded_reasons$se <- optim_result$se_unavailable_reason
  }

  # SAS/C's acceptance test, applied by .hzr_optim_generic() and polished
  # towards when BFGS stopped short of it. Both results are recorded on
  # every fit and shown by print() and summary(). Only the polish's two hard
  # failures warn -- the ones SAS/C reports as "reached no convergence"
  # (nlm code 4) and "unbounded ... or has a finite asymptote" (code 5). A
  # code 2 or 3 stop (step too small, or no lower point found) is where SAS
  # prints a caution and retries; on the test suite about a third of stops
  # end there, mostly on deliberately awkward fixtures, and warning on each
  # would bury the two that matter.
  if (fit_ran) {
    fit_state$rel_gradient <- optim_result$rel_gradient
    fit_state$rel_gradient_reason <- optim_result$rel_gradient_reason
    fit_state$polish_code  <- optim_result$polish_code
    # Codes 4 and 5 imply a failed test when nlm() and the statistic use the
    # same gradient; under CoE they need not, so the statistic is checked too.
    if (isTRUE(fit_state$converged) &&
        isTRUE(fit_state$polish_code %in% c(4L, 5L)) &&
        !isTRUE(fit_state$rel_gradient <= .Machine$double.eps^(1 / 3))) {
      warning(
        "The optimizer reported convergence, but the estimates fail the ",
        "relative-gradient test SAS/C HAZARD requires (at most ",
        signif(.Machine$double.eps^(1 / 3), 3), "; ",
        if (is.finite(fit_state$rel_gradient)) {
          paste0("here ", signif(fit_state$rel_gradient, 3))
        } else {
          "here the gradient could not be evaluated"
        },
        "). ",
        if (identical(fit_state$polish_code, 5L)) {
          paste0("The likelihood kept rising along a direction in which no ",
                 "maximum was found; the model may not have one.")
        } else {
          paste0("Further optimization stopped at its iteration limit; ",
                 "the estimates may not be at the maximum. Raise ",
                 "control$maxit to continue.")
        },
        call. = FALSE
      )
    }
  }

  # An ill-conditioned Hessian already warns that standard errors are
  # unreliable. That understates a ridge: where the likelihood is near-flat
  # along a parameter combination, the point estimates along it are not
  # determined either, and a fit that reports converged with an ordinary-
  # looking coefficient table gives no sign of it. Name the parameters
  # involved instead of leaving them to be read as estimated quantities.
  # Runs for every distribution -- nothing here knows about phase shapes.
  # optim() hands back an unnamed `par` for the single-distribution fits, so
  # taking names() alone would label the warning "par1"/"par2" while
  # summary() prints the real parameter names against the same numbers, and
  # the reader cannot map one onto the other. Fall back to the same source
  # summary() uses, so the two cannot disagree.
  weak_names <- names(fit_state$par)
  if (is.null(weak_names) && !is.null(fit_state$par)) {
    weak_names <- .hzr_parameter_names(
      theta = fit_state$par, dist = dist,
      p = if (is.null(x_fit)) 0L else ncol(x_fit)
    )
  }
  # is.list(), not !is.null(): the detector now returns NA when it could not
  # look at all (no Hessian, or a non-finite covariance among the estimated
  # parameters), and NULL only when it looked and found nothing. Testing for
  # non-NULL would warn on the NA and print a message built from empty fields.
  # A shape derived by hzr_phase(constraint = ) has a delta-method variance
  # but is an exact function of its sources, so it would read as a ridge the
  # constraint itself created. Only estimated parameters can trade off; every
  # other masked row is NA already.
  weak_vcov <- fit_state$vcov
  masked <- fit_state$fixed_mask
  if (is.matrix(weak_vcov) && length(masked) == nrow(weak_vcov)) {
    weak_vcov[which(masked), ] <- NA_real_
    weak_vcov[, which(masked)] <- NA_real_
  }
  # The g3 shapes the single-parameter reading may name (#415), from the
  # phase specs rather than from a name suffix.
  weak_shapes <- if (dist == "multiphase" && length(phases)) {
    unlist(lapply(seq_along(phases), function(k) {
      if (identical(phases[[k]]$type, "g3")) {
        paste0(names(phases)[[k]], ".", c("gamma", "alpha", "eta"))
      }
    }))
  }
  weak_check <- .hzr_weak_direction_impl(weak_vcov, fit_state$rcond,
                                         weak_names, theta = fit_state$par,
                                         shape_names = weak_shapes)
  fit_state$weak <- weak_check$weak
  degraded_reasons$weak <- weak_check$reason
  if (is.list(fit_state$weak)) {
    warning(.hzr_weak_direction_message(fit_state$weak), call. = FALSE)
  }

  # A phase fitted outside the support its parameterisation can carry (#444).
  # A sibling of $weak, with the same tri-state: NA not examined, NULL
  # examined and nothing found, a list of records otherwise.
  # Only rows the likelihood reads: a weight-0 row is excluded from the fit,
  # and so is a row the multiphase designs drop for an NA covariate, so
  # neither may supply the first observed time (#444) or an endpoint of a
  # step (#448).
  in_fit <- if (is.null(weights)) rep(TRUE, length(time)) else weights > 0
  if (length(fit_state$rows_used) == length(in_fit)) {
    in_fit <- in_fit & fit_state$rows_used
  }
  # And only the bounds the likelihood evaluates on each row: `time_lower` is
  # an entry time for status 0/1 and an interval's lower bound for status 2,
  # and is ignored on a left-censored row; `time_upper` is read only for
  # status -1/2. A supplied bound the likelihood never reads must not
  # stretch or split the span either check looks at.
  rows_with <- function(v, codes) {
    if (length(v) != length(in_fit)) return(v)
    v[in_fit & !is.na(status) & status %in% codes]
  }
  boundary_check <- .hzr_boundary_check_impl(
    theta = fit_state$theta, phases = phases,
    # `time` itself is read only where no explicit bound replaces it: always
    # for status 0/1; for a left-censored row only without `time_upper`; for
    # an interval row as either bound that was not supplied.
    time = rows_with(time, c(0, 1,
                             if (is.null(time_upper)) -1,
                             if (is.null(time_lower) || is.null(time_upper)) 2)),
    fitted = fit_ran,
    time_lower = rows_with(time_lower, c(0, 1, 2)),
    time_upper = rows_with(time_upper, c(-1, 2))
  )
  fit_state$boundary <- boundary_check$boundary
  degraded_reasons$boundary <- boundary_check$reason
  # Holds made at setup (#415) are records of the same family, prepended.
  # When the post-fit check could not run ($boundary NA, e.g. no positive
  # observed times) the field stays NA, as the degraded record requires, and
  # the holds are still announced below rather than dropped.
  boundary_records <- fit_state$boundary
  if (fit_ran && length(optim_held <- fit_state$held)) {
    if (.hzr_is_na_scalar(fit_state$boundary)) {
      boundary_records <- optim_held
    } else {
      fit_state$boundary <- c(optim_held, fit_state$boundary)
      boundary_records <- fit_state$boundary
    }
  }
  if (is.list(boundary_records) && length(boundary_records)) {
    warning(.hzr_boundary_condition(boundary_records))
  }

  # Refit-based tooling (hzr_bootstrap()) re-evaluates $call, so it needs the
  # bindings that call refers to. Capturing parent.frame() wholesale would pin
  # the caller's entire frame to every fitted object -- measured at a 1400x
  # saveRDS bloat, and it would drag unrelated data (including other cohorts)
  # into any saved model. Copy only the symbols the call actually references.
  # Computed here, not inside list(), so parent.frame() unambiguously resolves
  # to hazard()'s caller.
  captured_call <- match.call()
  captured_env <- .hzr_capture_call_env(captured_call, parent.frame())

  # Assemble the hazard S3 object.
  # $call       -- captured call for reproducibility / print
  # $call_env   -- the bindings $call references, copied out of hazard()'s
  #                caller. Refit-based tooling such as hzr_bootstrap()
  #                re-evaluates the stored call; it must do so here, otherwise
  #                arguments passed by symbol (theta, phases, control) cannot
  #                be resolved.
  # $spec       -- model specification (dist, control)
  # $data       -- raw data stored for default predict() / refit
  # $fit        -- optimisation results (see fit_state fields above)
  # $legacy_args -- pass-through ... args for SAS-migration parity
  # $engine     -- implementation tag ("native-r-m2")
  obj <- list(
    call = captured_call,
    call_env = captured_env,
    spec = list(dist = dist, control = control, time_windows = time_windows,
                phases = phases, objective = objective),
    data = list(
      time = time,
      time_lower = time_lower,
      time_upper = time_upper,
      status = as.numeric(status),
      x = x,
      # Terms, factor levels and contrasts of the formula RHS that built `x`
      # (formula path; NULL otherwise), so predict(newdata = ) can rebuild it.
      x_design = x_design,
      weights = weights,
      # Rows dropped for time = 0 before fitting (#374): every stored vector
      # above, and `frame`, is what remains.
      dropped_time_zero = n_dropped_time_zero,
      # Their positions in the rows hazard() was given, so a caller that
      # passes the original data frame on (hzr_stepwise(data = )) can be
      # aligned with the fit.
      dropped_time_zero_rows = dropped_time_zero_rows,
      # And the dropped rows themselves, so hzr_stepwise() can confirm that a
      # frame it is given is the one hazard() was given before trimming it.
      dropped_time_zero_frame = dropped_frame,
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

  # What this fit did not do, and why (#242). The object's state decides which
  # entries appear; the reasons carried up from the optimizer say why. The
  # validator stops if the two disagree, because that is a package bug.
  record <- .hzr_degraded_record(
    vcov = fit_state$vcov, weak = fit_state$weak,
    boundary = fit_state$boundary, control = control,
    dist = dist, fitted = fit_ran,
    fixed_mask = fit_state$fixed_mask, param_names = weak_names,
    reasons = degraded_reasons
  )
  obj$degraded <- record$degraded
  obj$degraded_causes <- record$degraded_causes
  .hzr_validate_degraded(obj, fitted = fit_ran)
  obj
}

#' Predict from a hazard model object
#'
#' Produces prediction outputs from a `hazard` object. Supports multiple prediction
#' types including linear predictor, hazard, survival probability, and cumulative hazard.
#'
#' @param object A `hazard` object. One built with `fit = FALSE` holds the
#'   starting values it was given rather than estimates, so predicting from
#'   it warns (condition class `hzr_unfitted_prediction`); under
#'   `dist = "multiphase"` it is an error instead, because the per-phase
#'   designs are resolved only when the model is fitted. To evaluate a model
#'   at parameters you supply, use [hzr_evaluate()]. The warning is governed
#'   by `options(TemporalHazard.warn_unfitted_prediction = )`, which this
#'   package's own tests set to `FALSE` where they exercise that capability
#'   deliberately; leaving it on is what tells a reader that a number came
#'   from a starting value.
#' @param newdata Optional matrix or data frame of predictors. For types requiring
#'   time (e.g., "survival", "cumulative_hazard"), newdata should include a `time`
#'   column, or time will be taken from the fitted object's data.
#'   Covariates are matched to the model by column name, so their order does
#'   not matter. A formula fit rebuilds its designs from the formulas,
#'   global and per phase, with the factor levels and contrasts the fit
#'   saw, so a factor can be given as a level label. `newdata` may
#'   instead carry the
#'   fit's design-matrix columns by name (`grpyoung` for a factor `grp`);
#'   these are used only when no formula variable is given (a numeric
#'   variable that is itself a column counts only if another column, such
#'   as `I(age^2)`, is built from it). With all the
#'   variables given, the design is rebuilt from them, so a design column
#'   that contradicts one is ignored; some variables beside the design
#'   columns, with others missing, is an error. That error is conservative
#'   in two cases where nothing contradicts: `poly()` design columns given
#'   with the variable they are built from but without another variable,
#'   and a non-syntactic name such as `my age` given both as the design
#'   column `` `my age` `` and as the variable. Give all of the formula's
#'   variables, or only the design columns. A
#'   column the model does not use is
#'   ignored, and a covariate the model needs but `newdata` lacks is an error.
#'   Only the columns of the model's `data` are taken from `newdata`: a
#'   formula constant (`cutoff` in `I(age > cutoff)`, spline knots) comes
#'   from the formula's environment, and a term that uses row-level values
#'   kept outside `data` (`~ zz`, with `zz` a vector in the workspace) is an
#'   error, even when `newdata` has a `zz` column; move it into `data` and
#'   refit.
#'   A term that computes a statistic over `newdata`'s rows warns only for
#'   the functions the check knows (see "How `newdata` is evaluated"). Any
#'   other function, including one you write, is still recomputed from
#'   `newdata`'s rows, silently: no warning means undetected, not safe.
#'   A fit made with an unnamed `x` matrix matches by position. For the types
#'   requiring time, a `newdata` with only a `time` column evaluates the
#'   baseline, with every covariate at 0. Because `time` is then the
#'   prediction time, a model whose formula uses a variable named `time`
#'   (a covariate, or a constant such as `I(age > time)`) cannot be given
#'   those types at `newdata`; rename it and refit.
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
#'
#'   **SAS draws narrower bands than this by default.** `PROC HAZPRED` takes
#'   its width from `CLEVEL`, whose default is `0.68268948`, documented in
#'   the macro source as "(1 sd)", so its `T_ALPHA` multiplier is `1` to
#'   seven decimals (the literal is truncated) and the band is one standard
#'   error, 68.3%, not 95%. Reproducing a SAS figure at this
#'   function's default therefore yields a band about 1.96 times wider than the
#'   one being checked against, with no error and no warning on either side.
#'   Pass the SAS level explicitly to match:
#'
#'   ```r
#'   predict(fit, newdata, type = "survival", se.fit = TRUE,
#'           level = 2 * stats::pnorm(1) - 1, conf.type = "logit")
#'   ```
#'
#'   The default is left at `0.95` deliberately: it is the right R-side
#'   default, and silently adopting SAS's would make this method disagree with
#'   every other R modelling function.
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
#' See the section "How `newdata` is evaluated" for what is recomputed from
#' `newdata` and when `predict()` warns.
#'
#' @section How `newdata` is evaluated:
#' `predict()` evaluates the model's formulas on `newdata` as given, as
#' [stats::predict.lm()] does. It uses the fit's factor levels and
#' contrasts, and the centering, basis and knots that a top-level `scale()`,
#' `poly()`, `ns()` or `bs()` term recorded. Everything else is recomputed
#' from `newdata`, so a prediction can differ from the fit without any
#' error. `predict()` warns, naming the cause, in three such cases. The
#' predicted values are the same with or without the warning.
#'
#' - **A term that computes a statistic over the rows.** In
#'   `I(age - mean(age))`, `I(scale(age)^2)` or
#'   `I(as.integer(factor(grp)))`, the mean, the scaling or the factor
#'   coding comes from `newdata`'s rows, so a row's prediction depends on
#'   which other rows are given. Compute such a variable in the data before
#'   fitting, and supply it in `newdata`. The check knows a fixed list of
#'   R's functions, among them `mean()`, `median()`, `min()`, `max()`,
#'   `quantile()`, `sd()`, `IQR()`, `ave()`, `rank()`, `length()`,
#'   `scale()`, `factor()` and `cut()` with a count of breaks. A function
#'   not on it, including one you write, is recomputed from `newdata`'s
#'   rows just the same, with no warning: the list is a floor, not a
#'   boundary.
#' - **A column of another type than the fit saw.** A numeric column given
#'   as character compares as text (`"154.6" > 50` is `FALSE`), and a
#'   `difftime` in other units is used in those units. The check compares
#'   against the fitting data the fit kept, which fits saved before
#'   TemporalHazard 1.1.0 do not have, so those fits are not checked. A
#'   factor given as its level labels, or an integer for a double, is not
#'   a mismatch.
#' - **A design rebuilt under this session's contrasts.** A formula fit
#'   saved by version 1.2.10 or earlier kept no record of its contrasts,
#'   and its design is rebuilt under `options(contrasts =)`. `predict()`
#'   warns when that option names a function other than `contr.treatment`
#'   or `contr.poly` and the rebuilt design is used. A redefined
#'   `contr.treatment` is not detected.
#'
#' Two cases are not detected:
#'
#' - A constant the formula reads from its environment, such as `cutoff`
#'   in `I(age > cutoff)`, is read when you predict, so a value changed
#'   since the fit is used.
#' - A comparison of strings, such as `I(grp > "b")`, follows the session's
#'   collation (`LC_COLLATE`), which can order strings differently from the
#'   session that fitted the model.
#'
#' @return When `se.fit = FALSE` (default), a numeric vector of predictions.
#'   When `se.fit = TRUE`, a data frame with columns `fit`, `se.fit`, `lower`,
#'   `upper` (delta-method point estimate, standard error, and confidence
#'   limits at `level`). `se.fit` is always the standard error of `fit` on
#'   its own scale: for `type = "survival"` that is `S * se(H)`, the delta
#'   method applied to `S = exp(-H)`. The survival limits are built from
#'   `se(H)` on the `conf.type` scale, so they are not `fit +/- z * se.fit`.
#'   `PROC HAZPRED` prints no standard error, only the limits.
#'   For multiphase `type = "cumulative_hazard"` with
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
  # A multiphase model built with fit = FALSE has no per-phase designs: they
  # are resolved at fit time, and without them .hzr_split_theta() looked up a
  # position that is not there and died with "argument of length 0" (#144).
  # After the argument checks, so a bad `se.fit` still reports itself, and
  # not for linear_predictor, which multiphase refuses for a fitted model too
  # and whose remedy is not "fit it". Other distributions predict from
  # supplied parameters perfectly well and are left alone.
  if (identical(object$spec$dist, "multiphase") &&
        !identical(type, "linear_predictor") &&
        is.null(object$fit$covariate_counts)) {
    stop("This multiphase model was built with fit = FALSE, so it has no ",
         "per-phase design matrices: they are resolved when the model is ",
         "fitted, and predict() cannot rebuild the phases without them. ",
         "Refit with fit = TRUE, or use hzr_evaluate() to evaluate the ",
         "model at parameters you supply.", call. = FALSE)
  }
  # The stored theta is checked against the stored design BEFORE any
  # prediction arithmetic, and for every type, because the downstream checks
  # are not equivalent. `hazard` and `linear_predictor` refuse a wrong length
  # where the design is multiplied as a matrix, but `survival` and
  # `cumulative_hazard` recycled a too-long theta into an outer product and
  # returned 2n values for n rows with no error (Codex review of #422). A
  # per-branch check would have to be repeated four times and kept in step;
  # one check ahead of the dispatch cannot fall out of step.
  #
  # The count is the one the theta was validated against at fit time: the
  # stored design, expanded by the time windows when there are any, which is
  # what `hazard()` and `hzr_evaluate()` both count. Multiphase is excluded
  # here as it is there, since fit = FALSE may legitimately carry fewer
  # entries (#408).
  # Only where there IS a stored design to check against. An object that
  # stored no `x` but carries covariate coefficients is a documented,
  # supported shape: `.hzr_newdata_design()` maps newdata's columns onto
  # those coefficients BY POSITION, because position is the only mapping
  # left. Refusing it here would kill that path (and did: it took
  # test-loglogistic-dist.R's supported case with it). Whether that
  # capability should survive at all is a separate decision, not one to make
  # as a side effect of a length check.
  if (!identical(object$spec$dist, "multiphase") && !is.null(object$data$x)) {
    x_stored <- object$data$x
    if (!is.null(time_windows)) {
      x_stored <- .hzr_expand_time_varying_design(
        x = x_stored, time = object$data$time, time_windows = time_windows
      )
    }
    .hzr_check_theta(theta, object$spec$dist,
                     n_coef = if (is.null(x_stored)) 0L else ncol(x_stored),
                     windowed = !is.null(time_windows))
  }

  # The other families predict from an unfitted object perfectly well, and
  # that is an intended, tested capability -- but the numbers come from the
  # starting values, not from estimates, and saying nothing is the
  # fit = FALSE hollow-chunk defect (#144). Classed, so a caller that meant
  # to supply parameters can muffle exactly this.
  if ((is.null(object$fit$converged) || is.na(object$fit$converged)) &&
        isTRUE(getOption("TemporalHazard.warn_unfitted_prediction", TRUE))) {
    warning(warningCondition(paste0(
      "This model was built with fit = FALSE: these predictions come from ",
      "the starting values it was given, not from estimates. Refit with ",
      "fit = TRUE for a fitted model's predictions, or use hzr_evaluate() ",
      "to evaluate a model at parameters you supply."
    ), class = "hzr_unfitted_prediction"))
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

  # newdata's `time` is the prediction time for the time-based types; a
  # model variable of that name would be misread there (#270). The
  # eta-based types have no prediction time unless there are time windows.
  time_based <- type %in% c("survival", "cumulative_hazard") ||
    identical(object$spec$dist, "multiphase") || !is.null(time_windows)
  if (!is.null(newdata)) {
    # A classed numeric column, such as bit64's integer64, stores doubles
    # that are not its values; read the values, as hazard() reads `data`
    # (#347). One rule for fitting and prediction.
    newdata <- .hzr_numeric_frame_values(newdata)
    # A formula fit saved before its design was stored gets it rebuilt, when
    # the rebuild is exact, so the by-name rules below apply to it (#301).
    # Design-level newdata (hzr_gof(), hzr_deciles()) never uses it, so it
    # skips the rebuild's cost.
    if (!isTRUE(attr(newdata, "hzr_design_columns"))) {
      recovered <- is.null(object$data$x_design)
      object <- .hzr_recover_x_design(object)
      recovered <- recovered && !is.null(object$data$x_design)
      # A legacy design rebuilt under this session's contrasts warns, but
      # only when the rebuild is used: newdata giving the design columns
      # uses them as they are (#335). A newdata the design route refuses
      # stops later, with its own message.
      nd_frame <- as.data.frame(newdata)
      if (recovered && !isTRUE(tryCatch(
        .hzr_uses_design_columns(object, nd_frame),
        error = function(e) TRUE
      ))) {
        .hzr_warn_rebuilt_contrasts(object$data$x_design$contrasts)
      }
      # A column of another type than the fit saw is evaluated as given;
      # warn once per call, naming it (#334).
      .hzr_warn_newdata_types(object, nd_frame)
    }
    .hzr_check_time_covariate(object, as.data.frame(newdata), time_based)
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
      # Covariates by name, not position (#267); NULL when there are none.
      # Without time windows `time` is no prediction time here, so it may
      # be a covariate (#270).
      x <- .hzr_newdata_design(object, newdata,
                               drop_time = !is.null(time_windows))
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
        # Route each phase the way the fit built it: the fit uses a phase's
        # own formula only on the formula interface, and a vector-interface
        # fit ignores hzr_phase(formula = ). .hzr_phase_inherits_global()
        # decides, from the fit's record or, for an older fit, its columns;
        # hzr_gof() calls the same helper, so the two cannot disagree.
        # Built once for every phase that inherits it, so its warnings are
        # given once per call.
        x_global <- NULL
        # A design computing over the rows warns once per call, naming
        # every phase, not once per phase.
        row_dependent <- character(0)
        withCallingHandlers({
          for (nm in names(phases)) {
            ph <- phases[[nm]]
            uses_formula <- !.hzr_phase_inherits_global(object, nm)
            if (uses_formula && ncol(nd_covs) > 0) {
              # The fit's levels, contrasts and columns, not newdata's.
              x_list[[nm]] <- .hzr_phase_newdata_design(object, nm, ph,
                                                        newdata)
            } else if (cov_counts[[nm]] > 0 && ncol(nd_covs) > 0) {
              # A formula-less phase inherits the global design: rebuild that,
              # not every non-time column of newdata (which also carries the
              # phase formulas' variables).
              if (is.null(x_global)) {
                x_global <- .hzr_global_design(object, newdata)
              }
              x_g <- x_global
              # With time windows the fit expanded the inherited design per
              # window (age_w1, age_w2); expand it the same way here, or the
              # rows meet the per-window coefficients unexpanded.
              if (!is.null(time_windows)) {
                x_g <- .hzr_expand_time_varying_design(
                  x = x_g, time = pred_time, time_windows = time_windows
                )
              }
              x_list[[nm]] <- x_g
            } else {
              x_list[[nm]] <- NULL
            }
          }
        }, hzr_row_dependent = function(w) {
          row_dependent <<- c(row_dependent, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
        if (length(row_dependent) > 0L) {
          warning(paste(row_dependent, collapse = "\n"), call. = FALSE)
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
        # By name, before any time-varying expansion below (#267).
        x <- .hzr_newdata_design(object, newdata)
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
    # unname() each result: theta's elements are named (mu, the leading one,
    # is a log rate, a scale or a location, not a shape), and R carries such
    # a name onto the product, through rep() for every n and through any
    # length-1 operand when n == 1, and so into predict() (#309; the
    # multiphase path does the same, #289).
    dist_lbl <- object$spec$dist
    has_cov <- !is.null(x) && ncol(x) > 0

    # The theta check that used to sit here has moved ahead of the type
    # dispatch, so it covers every prediction type rather than the two that
    # reach this line. It raised on the same theta through the same helper, so
    # nothing here can now fire that did not fire earlier.

    cumhaz_of <- if (dist_lbl == "weibull") {
      function(th) {
        if (th[1] <= 0 || th[2] <= 0) return(rep(NA_real_, length(time)))
        beta_cand <- if (length(th) > 2) th[3:length(th)] else numeric(0)
        eta_cand <- if (has_cov) as.numeric(x %*% beta_cand) else rep(0, length(time))
        unname((th[1] * time) ^ th[2] * exp(eta_cand))
      }
    } else if (dist_lbl == "exponential") {
      function(th) {
        beta_cand <- if (length(th) > 1) th[2:length(th)] else numeric(0)
        eta_cand <- if (has_cov) as.numeric(x %*% beta_cand) else rep(0, length(time))
        unname(exp(th[1]) * time * exp(eta_cand))
      }
    } else if (dist_lbl == "loglogistic") {
      function(th) {
        beta_cand <- if (length(th) > 2) th[3:length(th)] else numeric(0)
        eta_cand <- if (has_cov) as.numeric(x %*% beta_cand) else rep(0, length(time))
        unname(log(1 + exp(th[1]) * (time ^ exp(th[2])) * exp(eta_cand)))
      }
    } else if (dist_lbl == "lognormal") {
      function(th) {
        beta_cand <- if (length(th) > 2) th[3:length(th)] else numeric(0)
        # AFT: covariates shift the location.
        eta_aft <- if (has_cov) th[1] + as.numeric(x %*% beta_cand) else rep(th[1], length(time))
        z <- (log(time) - eta_aft) / exp(th[2])
        unname(-pnorm(-z, log.p = TRUE))
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


# The SAS/C acceptance test's result, for print() and summary(). NULL when the
# test was not applied: a fit that did not report convergence, or an object
# with no record of it (imported from SAS, or saved by an earlier version).
# A converged fit whose gradient could not be evaluated says so, because
# printing nothing would read as a test that never ran; the nlm() code is
# shown whenever there is one. It also says WHY, when the fit recorded a
# reason: "not evaluated" alone reads as a failure the fit is hiding, and
# under Conservation of Events it is the ordinary outcome (#351). Objects
# fitted before the reason was recorded carry none, and print as before.
.hzr_format_gradient_test <- function(rel_gradient, polish_code,
                                      converged = TRUE,
                                      reason = NA_character_) {
  if (!isTRUE(converged) || length(rel_gradient) != 1L) return(NULL)
  has_code <- length(polish_code) == 1L && !is.na(polish_code)
  if (is.na(rel_gradient)) {
    has_reason <- length(reason) == 1L && !is.na(reason) && nzchar(reason)
    return(paste0("  gradient:     not evaluated at the estimates",
                  if (has_reason) paste0(": ", reason),
                  if (has_code) paste0(" (nlm code ", polish_code, ")")))
  }
  gradtl <- .Machine$double.eps^(1 / 3)
  verdict <- if (rel_gradient <= gradtl) "met" else "not met"
  paste0("  gradient:     relative ", signif(rel_gradient, 3),
         " (SAS/C requires <= ", signif(gradtl, 3), "; ", verdict,
         if (has_code) paste0(", nlm code ", polish_code), ")")
}

#' Print method for fitted hazard models
#'
#' Compact one-block summary of a fitted `hazard` object: sample size,
#' number of predictors, distribution, theta vector, and log-likelihood,
#' followed by the "Not done in this run" block described in [hazard()].
#' S3 dispatch only: users call `print(fit)` rather than invoking this
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
    cat(.hzr_format_gradient_test(x$fit$rel_gradient, x$fit$polish_code,
                                  converged = x$fit$converged,
                                  reason = x$fit$rel_gradient_reason),
        sep = "\n")
  }
  # Always printed, "none" included (#242).
  cat(.hzr_format_not_done(x$degraded, x$degraded_causes), sep = "\n")
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
      # A shape derived by hzr_phase(constraint = ) keeps its delta-method
      # standard error but was not estimated, so it is not tested against 0.
      masked <- object$fit$fixed_mask
      if (length(masked) == length(valid)) valid <- valid & !(masked %in% TRUE)
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
    rel_gradient = object$fit$rel_gradient,
    rel_gradient_reason = object$fit$rel_gradient_reason,
    polish_code = object$fit$polish_code,
    log_lik = object$fit$objective,
    counts = object$fit$counts,
    message = object$fit$message,
    coefficients = coef_table,
    has_vcov = !is.null(vcov_mat) && is.matrix(vcov_mat),
    rcond = object$fit$rcond,
    pd = object$fit$pd,
    weak = object$fit$weak,
    boundary = object$fit$boundary,
    degraded = object$degraded,
    degraded_causes = object$degraded_causes,
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
#' positive-definite, a note warns that the standard errors may be unreliable,
#' and a further note names the parameters spanning a weakly identified
#' direction when one was found.  A "Not done in this run" block is always
#' printed: it lists each step this fit did not perform, with the reason, and
#' reads "none" when nothing was lost.  S3 dispatch only: users
#' call `print(summary(fit))` rather than invoking this directly.
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
    cat(.hzr_format_gradient_test(x$rel_gradient, x$polish_code,
                                  converged = x$converged,
                                  reason = x$rel_gradient_reason), sep = "\n")
  }
  if (!is.null(x$log_lik) && !is.na(x$log_lik)) {
    cat("  log-lik:     ", format(x$log_lik, digits = 6), "\n")
  }
  if (!is.null(x$rcond) && !is.na(x$rcond) && x$rcond < .hzr_rcond_tol) {
    cat("  Note: Hessian ill-conditioned (rcond = ",
        format(x$rcond, digits = 3),
        "); standard errors may be unreliable.\n", sep = "")
  }
  if (is.list(x$weak)) {
    # Wrapped rather than cat()'d flat: the message names parameters and two
    # diagnostics, and an unwrapped line buries them off the right edge.
    cat(strwrap(paste0("Note: ", .hzr_weak_direction_message(x$weak)),
                width = 76, indent = 2, exdent = 8),
        sep = "\n")
    cat("\n")
  }
  if (is.list(x$boundary)) {
    # Reported for the same reason as $weak: a fit that looks converged and
    # is a supremum is this package's signature defect (#444).
    cat(strwrap(paste0("Note: ", .hzr_boundary_message(x$boundary)),
                width = 76, indent = 2, exdent = 8),
        sep = "\n")
    cat("\n")
  }
  if (!is.null(x$pd) && !is.na(x$pd) && !isTRUE(x$pd)) {
    cat("  Note: Hessian not positive-definite at the optimum; ",
        "standard errors may be unreliable.\n", sep = "")
  }
  # Always printed, "none" included (#242): a line that appears only on bad
  # news cannot be told from a line its author forgot to write. It replaces
  # the notes that said "not examined for a weakly identified direction" and
  # "standard errors unavailable; the Hessian could not be inverted", both of
  # which could name the wrong cause.
  cat(.hzr_format_not_done(x$degraded, x$degraded_causes), sep = "\n")
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

#' Capture the bindings a stored call references
#'
#' Copies out of `envir` only the symbols `cl` actually refers to. Names that
#' resolve from the model's data frame rather than the calling scope (formula
#' column names such as `int_dead`/`dead`) do not exist in `envir` and
#' are skipped.
#'
#' @param cl Matched call, as returned by `match.call()`.
#' @param envir Caller environment to copy bindings out of.
#' @return A new environment holding just the referenced bindings.
#' @noRd
.hzr_capture_call_env <- function(cl, envir) {
  # all.names(), deliberately NOT all.vars(): the call may invoke the user's own
  # helper functions (e.g. one that builds `phases` or `theta`), and those live
  # in the caller's scope, not on globalenv()'s search path. all.vars() returns
  # only variables and omits function names, leaving such helpers unresolvable
  # once the call is re-evaluated.
  #
  # The scope-chain walk below copies only bindings the caller itself owns: it
  # stops at the first environment that is package territory rather than user
  # scope. Package and base functions (`hazard`, `list`, `Surv`, ...) are
  # deliberately left out: R serialises a closure's environment as a namespace
  # reference but its BODY by value, so capturing them would freeze a copy of
  # each function's body into every saved fit. A fit saved under one version and
  # bootstrapped after an upgrade would then run the stale body against the new
  # namespace, throw inside the replicate, and be swallowed by hzr_bootstrap()'s
  # tryCatch into a silent n_success = 0. Left uncaptured, they resolve through
  # `out`'s parent instead.
  #
  # The stop set is globalenv() plus any namespace/base env, NOT globalenv()
  # alone. Called from a script the chain is simply frame -> globalenv(), but
  # under testthat and R CMD check the package namespace sits in the chain
  # BEFORE globalenv(), so a globalenv()-only stop would still capture
  # `hazard`. When `envir` IS globalenv() the loop body never runs and nothing
  # is captured -- correct, since every symbol resolves through the parent.
  is_pkg_env <- function(e) {
    identical(e, globalenv()) || identical(e, emptyenv()) ||
      identical(e, baseenv()) || isNamespace(e)
  }
  syms <- unique(all.names(cl))
  found <- character()
  for (nm in syms) {
    e <- envir
    while (!is_pkg_env(e)) {
      if (exists(nm, envir = e, inherits = FALSE)) {
        found <- c(found, nm)
        break
      }
      e <- parent.env(e)
    }
  }
  # Parented to globalenv(), deliberately NOT to `envir`: parenting to the
  # caller would make its whole frame reachable again through the parent chain
  # and defeat the point of copying only the referenced symbols.
  out <- new.env(parent = globalenv())
  if (length(found)) {
    list2env(mget(found, envir = envir, inherits = TRUE), envir = out)
  }
  out
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

  .hzr_refuse_duplicate_columns(x)
  x
}

#' Refuse a design matrix whose column names repeat
#'
#' A factor's dummy columns are named `<factor><level>`, so a factor `g` with
#' level `b` and a numeric column `gb` both produce a column `gb`. Coefficient
#' names, predict() by name, the stepwise scope and the bootstrap all assume
#' the names are unique, and a duplicate used to fit without a word.
#'
#' @param x Numeric design matrix.
#' @param phase Phase name, for a multiphase phase-specific design; NULL for
#'   the global design.
#' @return `x`, invisibly, when its column names are unique.
#' @noRd
.hzr_refuse_duplicate_columns <- function(x, phase = NULL) {
  # Empty and NA names are absent names, not a repeated one:
  # cbind(a = v1, v2, v3) names its columns c("a", "", "").
  nms <- colnames(x)
  dup <- unique(nms[duplicated(nms) & !is.na(nms) & nzchar(nms)])
  if (length(dup) == 0L) {
    return(invisible(x))
  }
  where <- if (is.null(phase)) "" else paste0(" for phase '", phase, "'")
  stop("In hazard(), the design matrix", where, " has a duplicated column ",
       "name: ", paste0("'", dup, "'", collapse = ", "), ". Each covariate ",
       "needs a unique name, because coefficients and predict(newdata =) ",
       "match covariates by name. A factor's dummy columns are named ",
       "<factor><level>, so factor `g` with level `b` collides with a numeric ",
       "column `gb`: rename the numeric column, or relevel or rename the ",
       "factor. For an 'x' matrix, give its columns unique names.",
       call. = FALSE)
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


#' The values of a classed numeric
#'
#' A classed numeric such as `bit64::integer64` passes `is.numeric()`, but its
#' stored doubles are not its values: `unclass()` of an integer64 1 is
#' 4.94e-324. Arithmetic, `Surv()` and `model.matrix()` read the stored
#' doubles, so a fit over such input returned its starting values as
#' converged (#231). The rule, applied identically wherever the package reads
#' numbers a caller supplied (fitting, and prediction via `newdata`): an
#' object (`is.object()`) that is numeric (`is.numeric()`) is replaced by
#' `as.numeric()`, which dispatches to the class's own method. Everything
#' else is returned unchanged: plain numerics, factors, and
#' `Date`/`POSIXct`/`difftime` (not `is.numeric()`).
#'
#' A `dim` matters only for a column of a data frame, where a matrix column
#' (`I(cbind(p, q))`, a `Surv`) is legitimate and flattening it would change
#' the model: `.hzr_numeric_frame_values()` passes `keep_dim = TRUE`. A
#' single argument such as `time` is one vector whatever its shape, so by
#' default a classed numeric with a `dim` is read as its values too.
#'
#' Under `keep_dim`, the values are still read, and the shape kept
#' (#371): a classed numeric matrix column was left entirely alone, so an
#' integer64 matrix column of `data` fitted on its raw doubles, and the same
#' column in `newdata` predicted on them. Which classes need reading cannot
#' be listed, so it is decided by behaviour: the class's own `as.numeric()`
#' is compared with the stored values as numbers, and the column is replaced
#' only when they differ. A `Surv` column, whose stored doubles ARE its
#' values, is therefore returned unchanged and keeps its class, which it
#' must, since a bare matrix is no longer a response. So is an integer
#' `AsIs` matrix, whose values equal its storage in another mode.
#'
#' @param x Any object.
#' @param keep_dim If `TRUE`, keep the `dim` of an object that has one,
#'   reading its values into a matrix of the same shape.
#' @return `x`, or its values, when the rule applies.
#' @noRd
.hzr_numeric_values <- function(x, keep_dim = FALSE) {
  if (!is.object(x) || !is.numeric(x)) {
    return(x)
  }
  values <- as.numeric(x)
  if (!keep_dim || is.null(dim(x))) {
    return(values)
  }
  # Compared as numbers, not with the storage mode: an integer AsIs matrix
  # has values equal to its storage, in another mode.
  if (identical(values, as.numeric(unclass(x)))) {
    # The class reads as its own storage (a Surv, a classed plain matrix):
    # nothing to read, and replacing it would drop a class that is load
    # bearing.
    return(x)
  }
  array(values, dim(x), dimnames(x))
}


# The `control` elements the fitter reads, by distribution (#376), derived
# from the code rather than the documentation: .hzr_optim_generic() reads
# maxit and reltol; .hzr_optim_multiphase() reads and strips the multiphase
# ones before the optimizer; shape_param_count is not read by the fitter but
# is read back from the stored spec$control by the single-distribution
# stepwise refit and score test (every multiphase path derives its own theta
# layout, so on a multiphase fit nothing reads it; #405). abstol is read only by .hzr_optim_generic()'s bounded
# (L-BFGS-B) branch, which every caller turns off (use_bounds = FALSE), so
# no fit reads it.
.hzr_control_names <- list(
  all = c("maxit", "reltol"),
  single = "shape_param_count",
  multiphase = c("n_starts", "conserve", "phase_share_tol", "start_seed")
)

# Names ?hazard documented as accepted, a translation emitted, or a user
# might reasonably pass, although no fit ever read them, with the reason each
# has no effect (#376). Like every unread name they warn, and the fit
# proceeds: an error inside a stepwise or bootstrap candidate refit would be
# recorded as a failed candidate and empty the screen.
.hzr_control_no_effect <- c(
  abstol = paste0("it is read only by a bounded optimizer that no fit ",
                  "hazard() runs uses; `reltol` is the tolerance that ",
                  "applies"),
  method = paste0("hazard() chooses its optimizer: ",
                  "BFGS, a quasi-Newton method (a multiphase fit may run a ",
                  "Nelder-Mead warm-up first, and a stop that fails SAS's ",
                  "gradient test continues with stats::nlm())"),
  condition = paste0("SAS's CONDITION= has no equivalent: hazard() has no ",
                     "condition-number stop, and reports the Hessian's ",
                     "conditioning after the fit instead"),
  nocov = "it suppresses printed output, and hazard() prints nothing",
  nocor = "it suppresses printed output, and hazard() prints nothing",
  fix = paste0("hazard() has never read it, so a fit given it was the ",
               "unconstrained fit with its \"fixed\" parameters free, and ",
               "results obtained with it may be affected; hold a parameter ",
               "with hzr_phase(fixed = ) on a phase of a ",
               "dist = \"multiphase\" model, since a single-distribution ",
               "model has no mechanism for fixing one"),
  quasi = paste0("hazard() has never read it; its optimizer is ",
                 "BFGS, a quasi-Newton method (a multiphase fit may run a ",
                 "Nelder-Mead warm-up first, and a stop that fails SAS's ",
                 "gradient test continues with stats::nlm())")
)


#' Warn about every `control` element the fit does not read
#'
#' `hazard()` accepted any `control` element, so one it never reads -- a
#' typo such as `n_startz`, or `fix`, which no fitting code has ever read --
#' left the fit as it would have been and said nothing (#376). Every such
#' element now draws one warning that names it and says why it does
#' nothing, and the fit proceeds, as `stats::optim()` does for unknown
#' `control` names. Nothing errors: an error inside a stepwise or bootstrap
#' candidate refit would be recorded as a failed candidate, so the screen
#' would report success having tested nothing.
#'
#' @param control The `control` list, already known to be a list.
#' @param dist The distribution name.
#' @return `control` restricted to the elements this fit reads, so that no
#'   consumer's `$` can partial-match a warned name (`control$n_starts`
#'   would read `n_starts_extra`); warns once about every other element.
#' @keywords internal
#' @noRd
.hzr_validate_control <- function(control, dist) {
  if (length(control) == 0L) {
    return(control)
  }
  nm <- names(control)
  if (is.null(nm)) {
    nm <- rep("", length(control))
  }
  unnamed <- is.na(nm) | !nzchar(nm)
  names_all <- nm
  nm <- nm[!unnamed]
  multiphase <- identical(dist, "multiphase")
  accepted <- c(.hzr_control_names$all,
                if (multiphase) {
                  .hzr_control_names$multiphase
                } else {
                  .hzr_control_names$single
                })
  no_effect <- intersect(nm, names(.hzr_control_no_effect))
  # A name the fitter reads, but only for another model.
  off_path <- intersect(nm, if (multiphase) {
    .hzr_control_names$single
  } else {
    .hzr_control_names$multiphase
  })
  unknown <- setdiff(nm, c(accepted, no_effect, off_path))
  notes <- c(
    if (length(no_effect) > 0L) {
      paste0("control$", no_effect, " (", .hzr_control_no_effect[no_effect],
             ")")
    },
    if (length(off_path) > 0L) {
      paste0("control$", off_path, " (it applies only to ",
             if (multiphase) {
               "single-distribution fits"
             } else {
               "dist = \"multiphase\""
             }, ")")
    },
    if (length(unknown) > 0L) {
      paste0("control$", unknown, " (not an element any fit reads; for dist",
             " = \"", dist, "\" the elements accepted without a warning are ",
             paste(accepted, collapse = ", "), ")")
    },
    if (any(unnamed)) {
      paste0(sum(unnamed), " unnamed element(s) (control is read by name, ",
             "as in list(maxit = 500))")
    }
  )
  if (length(notes) > 0L) {
    warning("'control' element(s) with no effect on this dist = \"", dist,
            "\" fit, ignored: ", paste(notes, collapse = "; "), ".",
            call. = FALSE)
  }
  control[!unnamed & names_all %in% accepted]
}

#' Warn when a masked argument names both a column and a caller variable
#'
#' `hazard()` evaluates an argument expression with `data` as the environment
#' and the caller's frame as its parent, so a column wins (the rule
#' `subset()`, `transform()` and `with()` use, and `stats::lm()` for
#' `weights`). A wrapper that forwards its own argument by name --
#' `f <- function(tt) hazard(data = d, time = tt, ...)` -- reads as "use the
#' caller's vector" and silently gets the column instead: a fit over the wrong
#' rows, no error, no warning. The column still wins, but a name that is
#' BOTH a column and visible from the calling frame is ambiguous enough to
#' say so. Both interfaces call this one helper, so they give one message
#' (#151, #392). The lexical walk stops at the global environment
#' (`.hzr_bound_locally()`): `inherits = FALSE` misses the wrapper case, and
#' `inherits = TRUE` reaches base, where a column named `c`, `t` or `df`
#' would warn on every call.
#'
#' The diagnosis is one sentence on both interfaces; the remedy differs,
#' because only the vector interface can drop `data` to reach the calling
#' frame, while the formula interface requires it.
#'
#' @param exprs Named list of the unevaluated argument expressions.
#' @param data The data frame or list the arguments are masked by.
#' @param env The calling frame.
#' @param interface `"vector"` or `"formula"`, choosing the remedy clause.
#' @param data_arg The caller's unevaluated `data` argument, named in the
#'   advice when it is a plain symbol.
#' @return `NULL`, invisibly; warns naming each ambiguous name.
#' @keywords internal
#' @noRd
.hzr_warn_masked_ambiguity <- function(exprs, data, env,
                                       interface = c("vector", "formula"),
                                       data_arg = NULL) {
  interface <- match.arg(interface)
  ambiguous <- lapply(exprs, function(e) {
    if (is.null(e)) {
      return(character(0))
    }
    # Not all.vars(): it counts the RHS of `$` as a variable, so
    # all.vars(quote(other$tt)) is c("other", "tt") and the warning names
    # `tt` -- a column that was never consulted -- while `data$tt`, the
    # remedy the warning itself prescribes, triggers it.
    # The warning is a diagnostic: an expression too deeply nested to walk
    # (generated code) is fitted unchecked rather than refused (#401 review).
    nms <- tryCatch(.hzr_ambiguity_symbols(e),
                    error = function(err) character(0))
    nms[nms %in% names(data) &
          vapply(nms, .hzr_bound_locally, logical(1), env = env)]
  })
  # Only an argument that is the name alone is known to have read the
  # column. Inside a larger expression the name may never be evaluated (a
  # lazy function argument), may be rebound first (a loop variable, an
  # assignment) or may be evaluated elsewhere (with(), local(envir =)), and
  # no reading of the syntax can tell (#401 review). So that case says
  # nothing about which value was used.
  bare <- vapply(exprs, is.symbol, logical(1))
  .hzr_warn_ambiguous_names(ambiguous[bare], "The column was used. ",
                            interface, data_arg)
  .hzr_warn_ambiguous_names(
    ambiguous[!bare],
    paste0("The name is part of a larger expression, and which value it ",
           "read, if any, is not checked. "),
    interface, data_arg
  )
  invisible(NULL)
}

#' Emit the ambiguity warning for one class of names
#'
#' @param ambiguous Named list (by argument) of ambiguous names.
#' @param outcome The sentence saying what is known about the value used.
#' @inheritParams .hzr_warn_masked_ambiguity
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.hzr_warn_ambiguous_names <- function(ambiguous, outcome, interface,
                                      data_arg) {
  ambiguous <- ambiguous[lengths(ambiguous) > 0L]
  if (length(ambiguous) > 0L) {
    warning(
      "In hazard(), ", paste(sprintf("'%s' (%s)",
                                     unlist(ambiguous, use.names = FALSE),
                                     rep(names(ambiguous),
                                         lengths(ambiguous))),
                             collapse = ", "),
      ": the name is both a column of 'data' and a variable visible from ",
      "the calling frame. ", outcome,
      # Name the caller's own data argument only when it is a plain symbol:
      # `data` itself is usually utils::data() in the caller's frame, and an
      # inline expression or magrittr's `.` cannot be written as a prefix.
      if (is.symbol(data_arg) && !identical(data_arg, quote(.))) {
        paste0("Write ", deparse(data_arg), "$<name> for the column, or ")
      } else {
        paste0("Refer to the column through the data frame passed as ",
               "'data', or ")
      },
      if (interface == "vector") {
        "omit 'data' to use the calling frame's value."
      } else {
        "give the calling frame's value a name that is not a column of 'data'."
      },
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' The names a masked argument looks up, for the ambiguity warning
#'
#' As `.hzr_mask_symbols()`, which skips the name after `$` and `@`, but
#' also skipping both operands of `::` and `:::`: `stats::runif(n)` looks up
#' `n`, never a `stats` or `runif` column, so neither can be ambiguous (#401
#' review). A namespace-qualified call's own arguments are still collected.
#' `.hzr_mask_symbols()` is left as it is: it also feeds the formula's
#' `data_vars` and the `time` check.
#'
#' @param e A language object, symbol or constant.
#' @return Character vector of symbol names, possibly empty.
#' @keywords internal
#' @noRd
.hzr_ambiguity_symbols <- function(e) {
  if (is.symbol(e)) {
    return(as.character(e))
  }
  if (!is.call(e)) {
    return(character(0))
  }
  head <- e[[1L]]
  if (is.symbol(head) && as.character(head) %in% c("::", ":::")) {
    return(character(0))
  }
  if (is.symbol(head) && as.character(head) %in% c("$", "@") &&
        length(e) >= 3L) {
    return(.hzr_ambiguity_symbols(e[[2L]]))
  }
  parts <- as.list(e)[-1L]
  if (!is.symbol(head)) {
    parts <- c(list(head), parts)
  }
  # A function literal's defaults are expressions too, but they sit in a
  # pairlist, which is not a call, so the walk would skip them.
  if (is.symbol(head) && identical(as.character(head), "function")) {
    parts <- c(as.list(e[[2L]]), list(e[[3L]]))
  }
  nms <- unlist(lapply(parts, .hzr_ambiguity_symbols), use.names = FALSE)
  # A missing argument, as in `x[, j]`, is the empty symbol; it names nothing.
  unique(nms[nzchar(nms)])
}

#' Apply `.hzr_numeric_values()` to every column of a data frame or list
#'
#' A column with a `dim` keeps it (`keep_dim = TRUE`), and its values are
#' still read (#371). Columns are
#' replaced in a local copy (`data[] <-`), so a caller's `data.table` is not
#' modified by reference. Anything that is not a list is returned unchanged.
#'
#' @param data A data frame, list, or `NULL`.
#' @return `data` with each column passed through `.hzr_numeric_values()`.
#' @noRd
.hzr_numeric_frame_values <- function(data) {
  if (is.list(data)) {
    data[] <- lapply(data, .hzr_numeric_values, keep_dim = TRUE)
  }
  data
}
