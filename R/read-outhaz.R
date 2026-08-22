#' Read a SAS `outhaz` estimate dataset
#'
#' `PROC HAZARD`'s `outhaz=` dataset stores the converged estimates and the
#' asymptotic variance-covariance matrix at full double precision, where the
#' printed `.lst` carries about seven significant figures. For any quantity the
#' dataset holds, it is the better parity reference: print precision stops
#' being the binding constraint and optimizer convergence takes over.
#'
#' The log-likelihood is *not* stored here; take it from the `.lst`.
#'
#' @param path Path to a `.sas7bdat` written by `outhaz=`, or to an `.rds`
#'   holding the same data frame.
#' @return An `hzr_outhaz` object: a list with `estimates` (named numeric),
#'   `status` (named integer, 1 free / 0 fixed), `vcov` (matrix over free
#'   parameters, dimnames set) and `flags` (named numeric of model-structure
#'   flags). When no parameter is free, `vcov` is a 0x0 matrix rather than
#'   `NULL`, so check its dimensions rather than `is.null()`.
#'
#' @section Experimental:
#' This function is experimental and its return shape is expected to change.
#' The result carries the fitted model, so [predict.hzr_outhaz()] predicts
#' from it -- that is what the `hzr_translate_sas(librefs = )` path emits --
#' but it is not a `hazard` object and none of the other `hazard` methods
#' apply to it. The `_STATUS_` coding is asserted against a synthetic fixture,
#' so a real `OUTHAZ=` file using a different convention would yield an empty
#' `vcov` alongside a fully populated `estimates`.
#' @export
#' @examples
#' f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
#' if (nzchar(f)) str(hzr_read_outhaz(f))
hzr_read_outhaz <- function(path) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    d <- readRDS(path)
  } else {
    if (!requireNamespace("haven", quietly = TRUE)) {
      stop("hzr_read_outhaz() needs the 'haven' package to read a .sas7bdat.",
           call. = FALSE)
    }
    d <- as.data.frame(haven::read_sas(path))
  }

  if (!all(c("_NAME_", "_EST_", "_STATUS_") %in% names(d))) {
    stop("hzr_read_outhaz(): '", path, "' is not an outhaz dataset -- it must ",
         "carry columns _NAME_, _EST_ and _STATUS_. Found: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  }

  nm <- as.character(d[["_NAME_"]])
  est <- stats::setNames(as.numeric(d[["_EST_"]]), nm)
  st_raw <- suppressWarnings(as.integer(d[["_STATUS_"]]))

  # Rows with NA status are model-structure flags (G1FLAG, FIXDEL0, ...),
  # not parameters. They carry no estimate row in the vcov block.
  is_param <- !is.na(st_raw)
  status <- stats::setNames(st_raw[is_param], nm[is_param])
  flags <- est[!is_param]

  free <- names(status)[status == 1L]
  vcov <- as.matrix(d[match(free, nm), free, drop = FALSE])
  dimnames(vcov) <- list(free, free)
  storage.mode(vcov) <- "double"

  structure(
    list(estimates = est[is_param], status = status, vcov = vcov,
         flags = flags),
    class = "hzr_outhaz"
  )
}

# ============================================================================
# Reconstructing a fitted model from an OUTHAZ dataset
# ============================================================================
#
# The layout, and every rule below, is read off the reference implementation
# (`hazard/src/hazard/writeOutputDatafile.c` writes the dataset,
# `hazard/src/hazpred/{hzpm,gethazr,hzpe}.c` reads it back):
#
#   * six flag rows, then `Ntheta = 3p + 11` parameter rows, `p` covariates;
#   * rows 1-8 are the shaping parameters, in a fixed order, stored on the
#     MODEL (natural) scale;
#   * then `E0` and this job's covariate names, `C0` and its covariate
#     aliases, `L0` and its aliases. The intercepts are stored on the
#     ESTIMATION scale, i.e. `E0 = log(mu_early)` (`DTRSMU`);
#   * a phase is in the model when its intercept's `_STATUS_` is 1 -- that is
#     `hzpe.c`'s own test, and it is not what `G1FLAG`/`G3FLAG` mean;
#   * `G1FLAG` (1-6) selects which special case of the G1 shaping function the
#     signs of `M` and `NU` put the fit in, so it is a cross-check on those
#     two estimates, not a phase indicator.
.hzr_outhaz_shape_names <- c("DELTA", "THALF", "NU", "M",
                             "TAU", "GAMMA", "ALPHA", "ETA")

.hzr_outhaz_flag_names <- c("G1FLAG", "FIXDEL0", "FIXMNU1",
                            "G3FLAG", "FIXGE2", "FIXGAE2")

#' G1 special case implied by the early-phase shape estimates
#'
#' The reference program stores this as `G1FLAG`; recomputing it from `M` and
#' `NU` and comparing is the one cheap check that the file's conventions are
#' the ones assumed here. Returns `NA` for the combination the decomposition
#' is undefined on (`m < 0` and `nu < 0`).
#' @noRd
.hzr_outhaz_g1flag <- function(m, nu) {
  if (m > 0 && nu > 0) return(1)
  if (m == 0 && nu > 0) return(2)
  if (m < 0 && nu > 0) return(3)
  if (m < 0 && nu == 0) return(4)
  if (m > 0 && nu < 0) return(5)
  if (m == 0 && nu < 0) return(6)
  NA_real_
}

#' Rebuild a prediction-ready `hazard` object from an `hzr_outhaz` object
#'
#' Returns a `hazard`-classed list carrying exactly what [predict.hazard()]
#' reads for a multiphase model: `$spec$dist`, `$spec$phases`,
#' `$spec$time_windows`, `$fit$theta`, `$fit$phases`,
#' `$fit$covariate_counts`, `$fit$x_list` and `$fit$vcov`. It is a prediction
#' object, not a refittable fit: there is no data in an `OUTHAZ=` dataset, so
#' `$call` and `$data$time` are `NULL` and `predict()` must be given
#' `newdata`.
#'
#' `need_vcov` controls the variance block only. The OUTHAZ covariance matrix
#' is on SAS's ESTIMATION scale (`DCOVAR` differentiates the optimizer's own
#' variables), which is not this package's internal scale: SAS estimates
#' `log|NU|` and `log|M|` where the multiphase likelihood carries `nu` and `m`
#' untransformed. The mapping is the diagonal Jacobian `dtheta_R/dtheta_SAS`
#' applied below. It is only valid where SAS's estimation variable is the
#' plain `log()` of the parameter, and the late phase is often not: see
#' `.hzr_outhaz_late_composite()`. Where the Jacobian is not determinable --
#' a composite late-phase estimation scale, a late parameter derived from an
#' estimated one under `FIXGE2`/`FIXGAE2`, or the `FIXMNU1` constraint that
#' ties `m` to `nu` -- this refuses rather than return standard errors built
#' on the wrong scale.
#'
#' @param object An `hzr_outhaz` object.
#' @param need_vcov Logical; build the variance block (required by
#'   `se.fit = TRUE`).
#' @return A `hazard` object.
#' @noRd
.hzr_outhaz_to_spec <- function(object, need_vcov = FALSE) {
  est <- object$estimates
  status <- object$status
  flags <- object$flags
  nm <- names(est)

  missing_flags <- setdiff(.hzr_outhaz_flag_names, names(flags))
  if (length(missing_flags)) {
    stop("This OUTHAZ dataset is missing the model-structure row(s) ",
         paste(missing_flags, collapse = ", "),
         ", so its phase structure cannot be recovered.", call. = FALSE)
  }

  n_par <- length(est)
  p <- (n_par - 11L) / 3L
  if (n_par < 11L || p != round(p)) {
    stop("This OUTHAZ dataset has ", n_par, " parameter rows; PROC HAZARD ",
         "writes 3p + 11 for p covariates, so the layout is not the one ",
         "this reader knows.", call. = FALSE)
  }
  p <- as.integer(p)
  if (p > 0L) {
    stop("This OUTHAZ dataset carries ", p, " covariate(s) (",
         paste(nm[10:(9 + p)], collapse = ", "), "). Reconstructing a ",
         "covariate model from an OUTHAZ dataset is not supported: predict ",
         "from a hazard() fit instead.", call. = FALSE)
  }
  expected <- c(.hzr_outhaz_shape_names, "E0", "C0", "L0")
  if (!identical(nm, expected)) {
    stop("This OUTHAZ dataset's parameter rows are ",
         paste(nm, collapse = ", "), "; PROC HAZARD writes ",
         paste(expected, collapse = ", "), ".", call. = FALSE)
  }

  # A phase is in the model when its intercept is estimated (hzpe.c).
  present <- c(early = status[["E0"]] == 1L,
               constant = status[["C0"]] == 1L,
               late = status[["L0"]] == 1L)
  if (!any(present)) {
    stop("This OUTHAZ dataset has no phase in the model: none of E0, C0, L0 ",
         "carries _STATUS_ = 1.", call. = FALSE)
  }

  phases <- list()
  theta <- numeric(0)
  theta_sas <- character(0)   # OUTHAZ row each theta slot came from
  theta_jac <- numeric(0)     # d(theta_R) / d(theta_SAS) at the estimates

  if (present[["early"]]) {
    if (est[["DELTA"]] != 0) {
      stop("This OUTHAZ dataset has DELTA = ", format(est[["DELTA"]]),
           ". The early-phase time transformation B(t) = (exp(DELTA t) - 1) ",
           "/ DELTA is only absorbed by hzr_decompos() at DELTA = 0, so this ",
           "fit cannot be reproduced here.", call. = FALSE)
    }
    m_hat <- est[["M"]]
    nu_hat <- est[["NU"]]
    want_flag <- .hzr_outhaz_g1flag(m_hat, nu_hat)
    if (is.na(want_flag) || flags[["G1FLAG"]] != want_flag) {
      stop("This OUTHAZ dataset reports G1FLAG = ", flags[["G1FLAG"]],
           " but M = ", format(m_hat), " and NU = ", format(nu_hat),
           " imply ", if (is.na(want_flag)) "no valid case" else want_flag,
           ". The file's conventions are not the ones this reader assumes.",
           call. = FALSE)
    }
    fixed_early <- c("t_half", "nu", "m")[
      c(status[["THALF"]], status[["NU"]], status[["M"]]) == 0L]
    phases$early <- hzr_phase("cdf", t_half = est[["THALF"]], nu = nu_hat,
                              m = m_hat, fixed = fixed_early)
    theta <- c(theta, est[["E0"]], log(est[["THALF"]]), nu_hat, m_hat)
    theta_sas <- c(theta_sas, "E0", "THALF", "NU", "M")
    # log_mu and log_t_half are the SAS estimation-scale variables already;
    # nu and m are exp() of theirs, so the derivative is the estimate itself.
    theta_jac <- c(theta_jac, 1, 1, nu_hat, m_hat)
  }

  if (present[["constant"]]) {
    phases$constant <- hzr_phase("constant")
    theta <- c(theta, est[["C0"]])
    theta_sas <- c(theta_sas, "C0")
    theta_jac <- c(theta_jac, 1)
  }

  if (present[["late"]]) {
    fixed_late <- c("tau", "gamma", "alpha", "eta")[
      c(status[["TAU"]], status[["GAMMA"]],
        status[["ALPHA"]], status[["ETA"]]) == 0L]
    phases$late <- hzr_phase("g3", tau = est[["TAU"]], gamma = est[["GAMMA"]],
                             alpha = est[["ALPHA"]], eta = est[["ETA"]],
                             fixed = fixed_late)
    theta <- c(theta, est[["L0"]], log(est[["TAU"]]), est[["GAMMA"]],
               est[["ALPHA"]], est[["ETA"]])
    theta_sas <- c(theta_sas, "L0", "TAU", "GAMMA", "ALPHA", "ETA")
    theta_jac <- c(theta_jac, 1, 1, est[["GAMMA"]], est[["ALPHA"]],
                   est[["ETA"]])
  }

  cov_counts <- stats::setNames(rep(0L, length(phases)), names(phases))
  x_list <- stats::setNames(vector("list", length(phases)), names(phases))
  free_slot <- status[theta_sas] == 1L

  fit <- list(
    theta = theta, par = theta,
    converged = TRUE, objective = NA_real_,
    se = NULL, vcov = NULL,
    phases = phases, covariate_counts = cov_counts, x_list = x_list,
    fixed_mask = !free_slot
  )

  if (need_vcov) {
    fit$vcov <- .hzr_outhaz_vcov(object, theta_sas, theta_jac, free_slot,
                                 present)
    fit$se <- sqrt(diag(fit$vcov))
  }

  structure(
    list(
      call = NULL, call_env = NULL,
      spec = list(dist = "multiphase", control = list(),
                  time_windows = NULL, phases = phases),
      data = list(time = NULL, time_lower = NULL, time_upper = NULL,
                  status = NULL, x = NULL, weights = NULL, frame = NULL),
      fit = fit,
      legacy_args = list(),
      engine = "sas-outhaz"
    ),
    class = "hazard"
  )
}

#' What scale PROC HAZARD estimated each late shape on
#'
#' A one-for-one transcription of the three conditionals in
#' `hazard/src/common/hzd_late_p2t.c` (`Common.status` 5, 7 and 6 are GAMMA,
#' ETA and ALPHA):
#'
#' | Parameter | plain `log()` when | otherwise |
#' |---|---|---|
#' | GAMMA | `g_two \|\| g3flag >= 3 \|\| status[ETA] == 1` | `log(gamma*eta - 2)` |
#' | ETA | `g_two \|\| g3flag >= 3` | `log(gamma*eta - 2)` |
#' | ALPHA | `(ga_two && !g_two) \|\| g3flag >= 3` | `log(gamma*eta/alpha - 2)` |
#'
#' Only the plain-`log()` branch has the diagonal derivative
#' `d(parameter)/d(theta_SAS) = parameter` that `.hzr_outhaz_to_spec()`
#' builds. Note which way round the generic case falls: an unconstrained late
#' phase is `G3FLAG = 1` with `FIXGE2 = FIXGAE2 = 0` (`setg3.c`; `G3FLAG = 0`
#' is that routine's initialisation sentinel, not a value it writes), and then
#' ETA and ALPHA are always on a composite scale and GAMMA is whenever ETA is
#' fixed. The composite scales are the common case, not the exotic one.
#'
#' @return Character vector named GAMMA, ETA, ALPHA: `NA` where SAS estimated
#'   the plain `log()`, otherwise the composite expression it estimated.
#' @noRd
.hzr_outhaz_late_composite <- function(status, flags) {
  g_two <- flags[["FIXGE2"]] != 0
  ga_two <- flags[["FIXGAE2"]] != 0
  g3flag <- flags[["G3FLAG"]]
  ge_2 <- "log(GAMMA*ETA - 2)"
  c(
    GAMMA = if (g_two || g3flag >= 3 || status[["ETA"]] == 1L) {
      NA_character_
    } else {
      ge_2
    },
    ETA = if (g_two || g3flag >= 3) NA_character_ else ge_2,
    ALPHA = if ((ga_two && !g_two) || g3flag >= 3) {
      NA_character_
    } else {
      "log(GAMMA*ETA/ALPHA - 2)"
    }
  )
}

#' Map an OUTHAZ covariance block onto this package's parameter scale
#'
#' See `.hzr_outhaz_to_spec()` for why the scales differ. Returns a full
#' `length(theta)` square matrix with `NA` rows and columns for the fixed
#' parameters, which is the convention `.hzr_free_vcov()` expects.
#' @noRd
.hzr_outhaz_vcov <- function(object, theta_sas, theta_jac, free_slot,
                             present) {
  status <- object$status
  flags <- object$flags

  if (flags[["FIXMNU1"]] != 0 && present[["early"]]) {
    stop("This OUTHAZ dataset was fitted under FIXMNU1 (M constrained to ",
         "1/NU). The early-phase covariance is then on a parameterisation ",
         "this reader cannot map onto the multiphase likelihood's, so ",
         "se.fit = TRUE would report standard errors built on the wrong ",
         "scale. Use se.fit = FALSE.", call. = FALSE)
  }
  if (present[["late"]]) {
    g_two <- flags[["FIXGE2"]] != 0
    ga_two <- flags[["FIXGAE2"]] != 0
    g3flag <- flags[["G3FLAG"]]

    composite <- .hzr_outhaz_late_composite(status, flags)
    hit <- names(composite)[!is.na(composite) &
                              status[names(composite)] == 1L]
    if (length(hit)) {
      stop("This OUTHAZ dataset estimates ",
           paste0(hit, " as ", composite[hit], collapse = ", "),
           ", not as log(", paste(hit, collapse = "), log("),
           "). That is what PROC HAZARD's late-phase transform does for ",
           "these flags (G3FLAG = ", g3flag, ", FIXGE2 = ", as.integer(g_two),
           ", FIXGAE2 = ", as.integer(ga_two),
           "); with no late-phase constraint at all it is the generic case, ",
           "not a special one. The stored covariance is on those composite ",
           "scales, which the diagonal Jacobian this reader applies cannot ",
           "map onto the multiphase likelihood's, so se.fit = TRUE would ",
           "report standard errors built on the wrong scale. Use ",
           "se.fit = FALSE.", call. = FALSE)
    }

    # A second way the mapping fails, and the dangerous one because nothing
    # about it is visible in the covariance block: `hzd_late_t2p.c` lines
    # 37-39, 57-59 and 90-94 DERIVE one late parameter from another under the
    # constraints -- eta = 2/gamma or gamma = 2/eta under FIXGE2, and
    # alpha = gamma*eta/2 for a fixed ALPHA under FIXGAE2. The derived
    # parameter carries _STATUS_ = 0, so the reconstruction above marks it
    # fixed and the delta method drops its contribution entirely (for
    # eta = 2/gamma that is d(eta)/d(gamma) = -2/gamma^2). The result is a
    # populated se.fit column that is quietly wrong, so refuse.
    derived <- c(
      if (g_two && status[["GAMMA"]] == 1L) "ETA = 2 / GAMMA",
      if (g_two && status[["ETA"]] == 1L) "GAMMA = 2 / ETA",
      if (ga_two && status[["ALPHA"]] == 0L &&
            (status[["GAMMA"]] == 1L || status[["ETA"]] == 1L))
        "ALPHA = GAMMA * ETA / 2"
    )
    if (length(derived)) {
      stop("This OUTHAZ dataset was fitted with ",
           paste(c("FIXGE2", "FIXGAE2")[c(g_two, ga_two)],
                 collapse = " and "),
           ", which holds ", paste(derived, collapse = " and "),
           ". A derived parameter carries _STATUS_ = 0 and so has no row in ",
           "the stored covariance, but its value moves with an estimated ",
           "one, so se.fit = TRUE would report standard errors that silently ",
           "omit that term. Use se.fit = FALSE.", call. = FALSE)
    }
  }

  n <- length(theta_sas)
  out <- matrix(NA_real_, n, n)
  idx <- which(free_slot)
  if (!length(idx)) return(out)
  sas <- theta_sas[idx]
  miss <- setdiff(sas, rownames(object$vcov))
  if (length(miss)) {
    stop("This OUTHAZ dataset has no covariance row for ",
         paste(miss, collapse = ", "),
         "; se.fit = TRUE needs one (was the fit written with NOCOV?).",
         call. = FALSE)
  }
  j <- theta_jac[idx]
  out[idx, idx] <- object$vcov[sas, sas, drop = FALSE] * outer(j, j)
  out
}

#' Predictions from a fit loaded out of a SAS `OUTHAZ=` dataset
#'
#' Rebuilds the multiphase model the `OUTHAZ=` dataset describes -- which
#' phases are in it, their shapes, and the fitted parameter vector -- and then
#' predicts exactly as [predict.hazard()] does.
#'
#' `newdata` is required. An `OUTHAZ=` dataset holds a converged model and no
#' data, so there is no fitted time vector to fall back on; without `newdata`
#' the prediction would be over nothing.
#'
#' @section What is reconstructed, and what is refused:
#'
#' A phase is in the model when its intercept row (`E0`, `C0`, `L0`) carries
#' `_STATUS_ = 1`, which is the test `PROC HAZPRED` itself applies. The early
#' phase becomes `hzr_phase("cdf")` and the late phase `hzr_phase("g3")`, with
#' the shape estimates read off the `DELTA`/`THALF`/`NU`/`M` and
#' `TAU`/`GAMMA`/`ALPHA`/`ETA` rows.
#'
#' These cases error rather than return a number that looks like a
#' prediction: a dataset carrying covariates (their coefficients cannot be
#' matched to `newdata` columns from the file alone); a non-zero `DELTA` (the
#' early-phase time transformation is not implemented); and a `G1FLAG` that
#' disagrees with the signs of the `M` and `NU` estimates.
#'
#' With `se.fit = TRUE` there are three more, because the stored covariance is
#' on SAS's estimation scale and has to be mapped onto this package's:
#'
#' * a fit constrained by `FIXMNU1`, which ties `M` to `1/NU`;
#' * a fit estimating `GAMMA`, `ALPHA` or `ETA` on one of PROC HAZARD's
#'   *composite* late-phase scales -- `log(GAMMA*ETA - 2)` or
#'   `log(GAMMA*ETA/ALPHA - 2)` rather than `log()` of the parameter. This is
#'   the ordinary unconstrained late phase, not an exotic case: with
#'   `G3FLAG = 1` and no `FIXGE2`/`FIXGAE2`, `ETA` and `ALPHA` are always on a
#'   composite scale and `GAMMA` is whenever `ETA` is fixed. So an
#'   early + constant + late fit with any free late shape gets point
#'   predictions but no standard errors;
#' * a fit under `FIXGE2` or `FIXGAE2` where a late parameter is *derived*
#'   from an estimated one (`ETA = 2/GAMMA`, `GAMMA = 2/ETA`,
#'   `ALPHA = GAMMA*ETA/2`). The derived parameter has no covariance row, so
#'   its contribution to the variance would be dropped without trace.
#'
#' @param object An `hzr_outhaz` object from [hzr_read_outhaz()].
#' @param newdata Data frame with a `time` column. Required.
#' @param type One of `"survival"`, `"hazard"` or `"cumulative_hazard"`.
#' @param se.fit Logical; add delta-method standard errors and confidence
#'   limits, as [predict.hazard()] does.
#' @param ... Passed to [predict.hazard()] (`level`, `conf.type`,
#'   `decompose`).
#' @return What [predict.hazard()] returns: a numeric vector, or a data frame
#'   with `fit`, `se.fit`, `lower` and `upper` when `se.fit = TRUE`.
#' @seealso [hzr_read_outhaz()], [predict.hazard()].
#' @export
#' @examples
#' f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
#' if (nzchar(f)) {
#'   fit <- hzr_read_outhaz(f)
#'   predict(fit, newdata = data.frame(time = c(1, 6, 12)), type = "survival")
#' }
predict.hzr_outhaz <- function(object, newdata,
                               type = c("survival", "hazard",
                                        "cumulative_hazard"),
                               se.fit = FALSE, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    stop("'newdata' is required: an OUTHAZ= dataset carries a fitted model ",
         "and no data, so there are no times to predict at.", call. = FALSE)
  }
  spec <- .hzr_outhaz_to_spec(object, need_vcov = isTRUE(se.fit))
  stats::predict(spec, newdata = newdata, type = type, se.fit = se.fit, ...)
}
