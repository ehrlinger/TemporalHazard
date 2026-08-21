# sas-parse-parms.R -- map a SAS PARMS statement onto hzr_phase() calls and a
# starting-theta call.
#
# PARMS carries the starting values for the multiphase optimizer: MUE/MUC/MUL
# scale each phase (theta is the log of those SAS scales), THALF/NU/M shape
# the early (G1/"cdf") phase, TAU/GAMMA/ALPHA/ETA shape the late (G3) phase,
# bare FIX<param> tokens freeze a parameter at its starting value, and bare
# WEIBULL is setopt(6) / SETG3_weibull() in the reference C: it does not add a
# phase, it constrains the existing late phase to alpha = eta = 1, which
# collapses the general G3 form to a Weibull cumulative hazard (spec S7.1).
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

# FIX<param> keywords are their own grammar tokens (FIXTHALF, not FIX+THALF),
# and the R parameter name is not always the lowercased SAS spelling --
# THALF's R name is "t_half". Map explicitly rather than guess; FIXDELTA has
# no owning hzr_phase() parameter and is left untranslated.
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
  bad_construct <- character(0)
  bad_reason <- character(0)
  opt_tail <- character(0)
  had_value <- FALSE

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
        had_value <- TRUE
      }
    }
  }

  list(names = names_out, had_value = had_value, options_tail = opt_tail,
       untranslated_construct = bad_construct, untranslated_reason = bad_reason)
}

#' `fixed=` value: a bare string for one entry, a `c(...)` call for several.
#' @noRd
.hzr_parms_fixed_call <- function(fixed) {
  if (length(fixed) == 0L) return(NULL)
  if (length(fixed) == 1L) return(fixed)
  as.call(c(quote(c), as.list(fixed)))
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
  # .hzr_parse_phase_covars(). Starting values are noted once per phase (not
  # once per covariate) because they are not yet mapped to theta.
  phase_covars <- list()
  for (ph in c("early", "constant", "late")) {
    raw <- covars[[ph]]
    if (is.null(raw)) {
      phase_covars[[ph]] <- character(0)
      next
    }
    parsed <- .hzr_parse_phase_covars(raw)
    phase_covars[[ph]] <- parsed$names
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
    if (parsed$had_value) {
      flag_bad(ph, "phase covariate starting values are not yet mapped to theta")
    }
  }

  phase_calls <- list()
  if (length(early)) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "cdf", early, phase_covars$early, fixed_early
    )
  }
  if (has_muc) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "constant", list(), phase_covars$constant, character(0)
    )
  }
  if (length(late)) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "g3", late, phase_covars$late, fixed_late
    )
  }

  list(
    phases = as.call(c(quote(list), phase_calls)),
    theta = as.call(c(quote(c), lapply(unname(mu), function(v) bquote(log(.(v)))))),
    has_phases = length(phase_calls) > 0L,
    untranslated = .hzr_untranslated_frame(
      line = rep(NA_integer_, length(bad_construct)),
      construct = bad_construct,
      reason = bad_reason
    )
  )
}
