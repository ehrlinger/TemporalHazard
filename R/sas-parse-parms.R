# sas-parse-parms.R -- map a SAS PARMS statement onto hzr_phase() calls and a
# starting-theta call.
#
# PARMS carries the starting values for the multiphase optimizer: MUE/MUC/MUL
# scale each phase, THALF/NU/M shape the early (G1/"cdf") phase, TAU/GAMMA/
# ALPHA/ETA shape the late (G3) phase, bare FIX<param> tokens freeze a
# parameter at its starting value, and bare WEIBULL is setopt(6) /
# SETG3_weibull() in the reference C: it does not add a phase, it constrains
# the existing late phase to alpha = eta = 1, which collapses the general G3
# form to a Weibull cumulative hazard (spec S7.1).
#
# `theta` is the full interleaved starting vector the multiphase engine
# expects -- one block per phase, in the same early -> constant -> late
# order the phases are emitted, each block laid out exactly as
# .hzr_phase_theta_names() names it: log(mu), then that phase's shape
# starts (THALF/TAU logged, NU/M/GAMMA/ALPHA/ETA untransformed), then one
# entry per covariate. See .hzr_parms_theta_block().
#
# Every operand keyword is resolved through .hzr_sas_token() before it is
# routed anywhere; an unresolved keyword is recorded in `untranslated`, never
# guessed at and never silently dropped (see AGENTS.md, "The one thing that
# destroys work").

# Canonical hzr_phase() argument order for each shape family. Both the shape
# arguments and the `fixed=` entries are emitted in this order regardless of
# the order the PARMS operands appeared in, so the result is deterministic.
.hzr_parms_early_arg <- c(THALF = "t_half", NU = "nu", M = "m")
.hzr_parms_late_arg  <- c(TAU = "tau", GAMMA = "gamma", ALPHA = "alpha", ETA = "eta")
.hzr_parms_mu_order  <- c("MUE", "MUC", "MUL")

# hzr_phase() default shape values (mirrored here so the theta starting
# vector agrees with what hzr_phase() itself defaults to when PARMS did not
# supply a given shape parameter). mu has no hzr_phase() argument -- its
# only documented default is .hzr_phase_start()'s mu_start = 0.1, used here
# for a built phase whose PARMS never supplied its MU value.
.hzr_parms_early_default <- c(t_half = 1, nu = 1, m = 0)
.hzr_parms_late_default  <- c(tau = 1, gamma = 1, alpha = 1, eta = 1)
.hzr_parms_default_mu_start <- 0.1

# FIX<param> keywords are their own grammar tokens (FIXTHALF, not FIX+THALF),
# and the R parameter name is not always the lowercased SAS spelling --
# THALF's R name is "t_half". Map explicitly rather than guess. FIXDELTA has
# no owning hzr_phase() parameter and is handled separately below: DELTA is
# unimplemented rather than absorbed, so what matters is its VALUE, not
# whether it was pinned.
.hzr_parms_fix_map <- c(
  FIXTHALF = "t_half", FIXNU = "nu", FIXM = "m",
  FIXTAU = "tau", FIXGAMMA = "gamma", FIXALPHA = "alpha", FIXETA = "eta"
)

#' Order of, and select from, a named list/vector by a canonical key order.
#' @noRd
.hzr_parms_ordered <- function(x, key_order) {
  x[intersect(key_order, names(x))]
}

#' Parse one phase's `EARLY`/`CONSTANT`/`LATE` operand text into covariate
#' names, discarding a `/ options` tail and non-numeric `VAR=VALUE` pairs.
#'
#' The real grammar (`phasevaropt : phasevar phaseval phaseoptspec`, with
#' `phaseoptspec : /*nothing*/ | '/' phaseopts`) is a comma-separated list of
#' `VAR=startvalue` pairs or bare `VAR`s, optionally followed by `/ options`
#' (the `PHOP` family -- `EXCLUDE`/`INCLUDE`/`MOVE`/`ORDER`/`START` --
#' deferred in v1 scope). `x` may be a single raw operand string (from the
#' job parser) or an already-split character vector of bare names (the
#' `.hzr_parse_parms()` `covars=` back-compat interface); both are handled by
#' splitting every element on `/` then `,`, which is a no-op on a plain bare
#' name.
#' @noRd
.hzr_parse_phase_covars <- function(x) {
  names_out <- character(0)
  values_out <- numeric(0)
  bad_construct <- character(0)
  bad_reason <- character(0)
  opt_tail <- character(0)

  for (piece in x) {
    slash <- .idx(piece, "/")
    if (slash > 0L) {
      tail <- trimws(substr(piece, slash + 1L, nchar(piece)))
      if (nzchar(tail)) opt_tail <- c(opt_tail, tail)
      piece <- substr(piece, 1L, slash - 1L)
    }
    parts <- strsplit(piece, ",", fixed = TRUE)[[1L]]
    for (p in parts) {
      p <- trimws(p)
      if (!nzchar(p)) next
      eq <- .idx(p, "=")
      if (eq == 0L) {
        names_out <- c(names_out, p)
        values_out <- c(values_out, NA_real_)
        next
      }
      var <- trimws(substr(p, 1L, eq - 1L))
      val_chr <- trimws(substr(p, eq + 1L, nchar(p)))
      val <- suppressWarnings(as.numeric(val_chr))
      if (is.na(val)) {
        bad_construct <- c(bad_construct, p)
        bad_reason <- c(
          bad_reason,
          sprintf("non-numeric value for phase-statement covariate %s", var)
        )
      } else {
        names_out <- c(names_out, var)
        values_out <- c(values_out, val)
      }
    }
  }

  list(names = names_out, values = values_out, options_tail = opt_tail,
       untranslated_construct = bad_construct, untranslated_reason = bad_reason)
}

#' `fixed=` value: a bare string for one entry, a `c(...)` call for several.
#' @noRd
.hzr_parms_fixed_call <- function(fixed) {
  if (length(fixed) == 0L) return(NULL)
  if (length(fixed) == 1L) return(fixed)
  as.call(c(quote(c), as.list(fixed)))
}

#' Build one phase's block of the full interleaved `theta` starting vector.
#'
#' Mirrors `.hzr_phase_theta_names()`'s layout exactly: `log_mu`, then (for
#' `"early"`/`"late"`) the shape parameters, then one entry per covariate.
#' `THALF`/`TAU` enter as logs (`log_t_half`/`log_tau`); `NU`/`M`/`GAMMA`/
#' `ALPHA`/`ETA` enter untransformed; `MUE`/`MUC`/`MUL` are always logged.
#' @param family One of `"early"`, `"constant"`, `"late"`.
#' @param mu_val Numeric SAS-scale mu starting value for this phase.
#' @param shape Named list of parsed shape starting values (may omit entries
#'   PARMS did not supply -- those fall back to the `hzr_phase()` default).
#' @param covar_vals Numeric vector of covariate starting values, `NA` where
#'   PARMS gave a bare name with no explicit start (defaults to 0, matching
#'   `.hzr_phase_start()`).
#' @return A list of language/numeric elements, in theta order.
#' @noRd
.hzr_parms_theta_block <- function(family, mu_val, shape, covar_vals) {
  block <- list(bquote(log(.(mu_val))))

  if (family == "early") {
    d <- .hzr_parms_early_default
    t_half <- if (!is.null(shape$t_half)) shape$t_half else d[["t_half"]]
    nu     <- if (!is.null(shape$nu))     shape$nu     else d[["nu"]]
    m      <- if (!is.null(shape$m))      shape$m      else d[["m"]]
    block <- c(block, list(bquote(log(.(t_half)))), list(nu), list(m))
  } else if (family == "late") {
    d <- .hzr_parms_late_default
    tau   <- if (!is.null(shape$tau))   shape$tau   else d[["tau"]]
    gamma <- if (!is.null(shape$gamma)) shape$gamma else d[["gamma"]]
    alpha <- if (!is.null(shape$alpha)) shape$alpha else d[["alpha"]]
    eta   <- if (!is.null(shape$eta))   shape$eta   else d[["eta"]]
    block <- c(block,
      list(bquote(log(.(tau)))), list(gamma), list(alpha), list(eta)
    )
  }

  if (length(covar_vals)) {
    covar_vals[is.na(covar_vals)] <- 0
    block <- c(block, as.list(unname(covar_vals)))
  }

  block
}

#' Build one `hzr_phase(...)` call.
#' @noRd
.hzr_parms_phase_call <- function(type, shape, covars, fixed) {
  call_args <- c(list(quote(hzr_phase), type), shape)
  if (length(covars)) {
    call_args <- c(call_args,
                    list(formula = str2lang(paste("~", paste(covars, collapse = " + ")))))
  }
  fixed_call <- .hzr_parms_fixed_call(fixed)
  if (!is.null(fixed_call)) call_args <- c(call_args, list(fixed = fixed_call))
  as.call(call_args)
}

#' Map a SAS `PARMS` statement's operands to phases and a starting theta.
#'
#' @param operands Character vector of `PARMS` tokens, e.g.
#'   `c("MUE=0.2", "THALF=0.15", "NU=1.4", "M=1", "FIXM", "MUC=0.0005")`.
#' @param covars Optional named list of phase covariates, e.g.
#'   `list(early = c("X1", "X2"), constant = , late = )`, from the operands of
#'   the `EARLY` / `CONSTANT` / `LATE` statements.
#' @return `list(phases = <call>, theta = <call>, has_phases = <logical>,
#'   untranslated = <data.frame>)`. `has_phases` is `TRUE` only when at least
#'   one phase was actually built from the operands (i.e. `phases` is not the
#'   empty `list()` call) -- callers use it to decide whether the job
#'   qualifies as multiphase at all.
#' @noRd
.hzr_parse_parms <- function(operands, covars = list()) {
  mu <- list()
  early <- list()
  late <- list()
  fixed_early <- character(0)
  fixed_late <- character(0)
  bad_construct <- character(0)
  bad_reason <- character(0)

  flag_bad <- function(construct, reason) {
    bad_construct <<- c(bad_construct, construct)
    bad_reason <<- c(bad_reason, reason)
  }

  for (op in operands) {
    eq <- .idx(op, "=")

    if (eq > 0L) {
      key <- substr(op, 1L, eq - 1L)
      val <- suppressWarnings(as.numeric(substr(op, eq + 1L, nchar(op))))
      token <- .hzr_sas_token(key, "HAZARD", "PARM")
      if (is.na(token)) {
        flag_bad(op, "unresolved PARMS keyword")
      } else if (is.na(val)) {
        flag_bad(op, sprintf("PARMS value for %s is not numeric", key))
      } else if (token %in% .hzr_parms_mu_order) {
        mu[[token]] <- val
      } else if (token %in% names(.hzr_parms_early_arg)) {
        early[[.hzr_parms_early_arg[[token]]]] <- val
      } else if (token %in% names(.hzr_parms_late_arg)) {
        late[[.hzr_parms_late_arg[[token]]]] <- val
      } else if (token == "DELTA") {
        # DELTA is not absorbed by the shape -- it is unimplemented, and R
        # assumes delta = 0 (see the header of R/decomposition.R). DELTA = 0 is
        # therefore a faithful translation and not a gap, exactly as the
        # print-only PROC options are. A non-zero one is not a missing feature
        # but a WRONG ANSWER: the emitted call fits a different function, with
        # no error. Say which of the two this is; the generic "no phase target"
        # reason fired identically on both and so distinguished nothing.
        if (!identical(val, 0)) {
          flag_bad(op, paste0(
            "DELTA = ", format(val), " is not implemented -- R assumes ",
            "delta = 0, so the emitted call fits a DIFFERENT model than this ",
            "job (rho, the time argument and the density Jacobian all differ)"
          ))
        }
      } else {
        flag_bad(op, "PARMS keyword has no phase target")
      }
      next
    }

    token <- .hzr_sas_token(op, "HAZARD", "PARM")
    if (is.na(token)) {
      flag_bad(op, "unresolved PARMS keyword")
    } else if (token == "WEIBULL") {
      late[["alpha"]] <- 1
      late[["eta"]] <- 1
      fixed_late <- union(fixed_late, c("alpha", "eta"))
    } else if (token %in% names(.hzr_parms_fix_map)) {
      param <- .hzr_parms_fix_map[[token]]
      if (param %in% .hzr_parms_early_arg) {
        fixed_early <- union(fixed_early, param)
      } else {
        fixed_late <- union(fixed_late, param)
      }
    } else if (token == "FIXDELTA") {
      # Pinning DELTA at whatever PARMS set it to. That value is what decides
      # whether this job is reproducible, and a DELTA= operand is flagged
      # above; FIXDELTA on its own leaves it at the SAS default of 0, which is
      # the branch R implements. Mapped, not a gap.
      NULL
    } else {
      flag_bad(op, "PARMS token has no phase target")
    }
  }

  early <- .hzr_parms_ordered(early, unname(.hzr_parms_early_arg))
  late <- .hzr_parms_ordered(late, unname(.hzr_parms_late_arg))
  fixed_early <- intersect(unname(.hzr_parms_early_arg), fixed_early)
  fixed_late <- intersect(unname(.hzr_parms_late_arg), fixed_late)
  mu <- .hzr_parms_ordered(mu, .hzr_parms_mu_order)
  has_muc <- "MUC" %in% names(mu)

  # EARLY/CONSTANT/LATE operand text: comma-separated VAR=VALUE pairs (or
  # bare VARs), optionally followed by a "/ options" tail. Non-numeric values
  # and the options tail are recorded to untranslated, never guessed at; see
  # .hzr_parse_phase_covars(). VAR=VALUE starting values are now mapped into
  # theta (one entry per covariate, appended after that phase's shape block,
  # per .hzr_phase_theta_names()); a bare VAR with no value defaults to 0,
  # matching .hzr_phase_start().
  phase_covars <- list()
  phase_covar_vals <- list()
  for (ph in c("early", "constant", "late")) {
    raw <- covars[[ph]]
    if (is.null(raw)) {
      phase_covars[[ph]] <- character(0)
      phase_covar_vals[[ph]] <- numeric(0)
      next
    }
    parsed <- .hzr_parse_phase_covars(raw)
    phase_covars[[ph]] <- parsed$names
    phase_covar_vals[[ph]] <- parsed$values
    for (i in seq_along(parsed$untranslated_construct)) {
      flag_bad(parsed$untranslated_construct[[i]], parsed$untranslated_reason[[i]])
    }
    if (length(parsed$options_tail)) {
      flag_bad(
        paste("/", paste(parsed$options_tail, collapse = " ")),
        sprintf(
          "%s phase options (EXCLUDE/INCLUDE/MOVE/ORDER/START) are deferred (v1 scope)",
          ph
        )
      )
    }
  }

  # Phases are built, and their theta blocks appended, in the same
  # early -> constant -> late order -- .hzr_phase_theta_names() assigns
  # labels by *position* (phases are auto-named "phase_1", "phase_2", ... by
  # .hzr_validate_phases()), so only this order matters, not the names.
  phase_calls <- list()
  theta_blocks <- list()
  if (length(early)) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "cdf", early, phase_covars$early, fixed_early
    )
    mu_val <- if (!is.null(mu[["MUE"]])) mu[["MUE"]] else .hzr_parms_default_mu_start
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("early", mu_val, early, phase_covar_vals$early)
    )
  }
  if (has_muc) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "constant", list(), phase_covars$constant, character(0)
    )
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("constant", mu[["MUC"]], list(), phase_covar_vals$constant)
    )
  }
  if (length(late)) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "g3", late, phase_covars$late, fixed_late
    )
    mu_val <- if (!is.null(mu[["MUL"]])) mu[["MUL"]] else .hzr_parms_default_mu_start
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("late", mu_val, late, phase_covar_vals$late)
    )
  }

  list(
    phases = as.call(c(quote(list), phase_calls)),
    theta = as.call(c(quote(c), theta_blocks)),
    has_phases = length(phase_calls) > 0L,
    untranslated = .hzr_untranslated_frame(
      line = rep(NA_integer_, length(bad_construct)),
      construct = bad_construct,
      reason = bad_reason
    )
  )
}
