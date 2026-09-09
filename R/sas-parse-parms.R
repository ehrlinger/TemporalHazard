# sas-parse-parms.R -- map a SAS PARMS statement onto hzr_phase() calls and a
# starting-theta call.
#
# PARMS carries the starting values for the multiphase optimizer: MUE/MUC/MUL
# scale each phase, THALF/NU/M shape the early (G1/"cdf") phase, TAU/GAMMA/
# ALPHA/ETA shape the late (G3) phase, bare FIX<param> tokens freeze a
# parameter at its starting value, and bare WEIBULL is setopt(6) /
# SETG3_weibull() in the reference C: the GENERALIZED Weibull, which admits
# all positive parameter values. It neither adds a phase nor constrains one --
# see the WEIBULL branch below.
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

# PROC HAZARD's OWN shape defaults (src/hazard/stmtprc.c:34-37), used for any shape
# operand a PARMS statement did not name. They are not hzr_phase()'s defaults
# -- nu, m and eta all differ (SAS 2/1/2 against R 1/0/1) -- and the starting
# vector is what the emitted call hands the optimizer, so mirroring R's
# defaults here started a partially-specified job somewhere PROC HAZARD would
# not have started it. On a multimodal likelihood a different start is a
# different answer, not a cosmetic difference (AGENTS.md, "n_starts explores a
# neighbourhood").
#
# TAU is the one entry that is NOT SAS's default. stmtprc.c:33 sets it to
# ZERO, and SETG3() (setg3.c:317) then replaces any non-positive TAU with
# 2*Tmax/3 -- data-dependent, so unreproducible at parse time. The 1 below is
# SETG3_ignore_tau()'s pinned value (setg3.c:378), which is right exactly when
# that branch fires; every other defaulted-TAU late phase gets an untranslated
# row instead. See the TAU guard in .hzr_parse_parms().
#
# mu has no hzr_phase() argument -- its only documented default is
# .hzr_phase_start()'s mu_start = 0.1, used here for a built phase whose PARMS
# never supplied its MU value. That is a separate divergence from SAS (SETG1
# floors a non-positive muE at 1, setg1.c) and, more to the point, a phase
# whose MU is absent is one PROC HAZARD would not build at all -- the
# phase-activation asymmetry tracked in the MU comment further down. Left
# alone here deliberately.
.hzr_parms_early_sas_default <- c(t_half = 1, nu = 2, m = 1)
.hzr_parms_late_sas_default  <- c(tau = 1, gamma = 1, alpha = 1, eta = 2)
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

#' Fill a parsed shape list out to every shape parameter of its family.
#'
#' The one place a default is applied. Both the emitted `hzr_phase()` call and
#' the `theta` starting vector are built from this function's result, so they
#' cannot disagree about what an unspecified operand started at -- filling only
#' one of the two would be worse than defaulting neither, because the printed
#' call would then describe a fit that did not happen.
#'
#' @param shape Named list of the shape values `PARMS` actually supplied.
#' @param defaults Named numeric vector of defaults, in canonical argument
#'   order.
#' @return A named list carrying every name of `defaults`, in that order.
#' @noRd
.hzr_parms_fill_shape <- function(shape, defaults) {
  out <- as.list(defaults)
  out[names(shape)] <- shape
  out[names(defaults)]
}

#' Build one phase's block of the full interleaved `theta` starting vector.
#'
#' Mirrors `.hzr_phase_theta_names()`'s layout exactly: `log_mu`, then (for
#' `"early"`/`"late"`) the shape parameters, then one entry per covariate.
#' `THALF`/`TAU` enter as logs (`log_t_half`/`log_tau`); `NU`/`M`/`GAMMA`/
#' `ALPHA`/`ETA` enter untransformed; `MUE`/`MUC`/`MUL` are always logged.
#' @param family One of `"early"`, `"constant"`, `"late"`.
#' @param mu_val Numeric SAS-scale mu starting value for this phase.
#' @param shape Named list of shape starting values, already filled out to
#'   every parameter of the family by `.hzr_parms_fill_shape()`. Completeness
#'   is asserted rather than defaulted: a short list here would emit a theta
#'   block of the wrong length, which `.hzr_optim_multiphase()` would either
#'   reject or, worse, silently mis-align against the phase's parameter names.
#' @param covar_vals Numeric vector of covariate starting values, `NA` where
#'   PARMS gave a bare name with no explicit start (defaults to 0, matching
#'   `.hzr_phase_start()`).
#' @return A list of language/numeric elements, in theta order.
#' @noRd
.hzr_parms_theta_block <- function(family, mu_val, shape, covar_vals) {
  block <- list(bquote(log(.(mu_val))))

  if (family == "early") {
    stopifnot(all(names(.hzr_parms_early_sas_default) %in% names(shape)))
    block <- c(block, list(bquote(log(.(shape$t_half)))),
               list(shape$nu), list(shape$m))
  } else if (family == "late") {
    stopifnot(all(names(.hzr_parms_late_sas_default) %in% names(shape)))
    block <- c(block, list(bquote(log(.(shape$tau)))),
               list(shape$gamma), list(shape$alpha), list(shape$eta))
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

# ---------------------------------------------------------------------------
# SETG3() trace
# ---------------------------------------------------------------------------
# PROC HAZARD does not optimize from the operands PARMS supplies: SETG3()
# (src/model/setg3.c) rewrites them first, and refuses some jobs outright.
# .hzr_setg3_notes() walks that function and reports what it would do, in the
# C's own order. It changes nothing about the emitted call.
#
# Which of the two happens is a deliberate split, settled once and applied
# throughout:
#
#   * A rewrite that resolves an EXACT non-identifiability is MIRRORED into
#     the emitted call, because the degeneracy is algebra and is just as real
#     in R -- SETG3_ignore_tau()'s TAU pin (setg3.c:378-379) and ETA fix
#     (:405) are the two, both handled at the build site above.
#   * Every other rewrite here keeps PROC HAZARD inside a numerical branch it
#     can evaluate -- the role g3flag plays, which hzr_decompos_g3() does not
#     need because it carries the general four-parameter G3 form. Copying
#     those would import a SAS limitation into R, and reaching a late shape
#     the reference cannot is a purpose of this package. They are RECORDED.
#
# Refusals are recorded for the reason SETG3980 always was: emitting a
# runnable fit for a job PROC HAZARD will not run is an answer that looks like
# a result and is not.
#
# g_two/ga_two (the GAMMA*ETA = 2 and GAMMA*ETA/ALPHA = 2 constraint flags)
# are driven by FIXGE2/FIXGAE2, which .hzr_sas_token() leaves unresolved and
# records as untranslated -- so a job carrying either is never reported clean,
# and this trace takes both as FALSE rather than guessing. Every branch below
# is therefore the !g_two && !ga_two column of the C's own tables.

# What each SETG3 refusal code objects to. The code alone is greppable but
# opaque; a caller reading $untranslated needs to know which operand to change.
.hzr_setg3_refusal <- c(
  "(SETG3900)" = "TAU is fixed at a non-positive value",
  "(SETG3910)" = "GAMMA is fixed at a non-positive value",
  "(SETG3920)" = "ALPHA is fixed at a negative value",
  "(SETG3930)" = "ETA is fixed at a non-positive value",
  "(SETG3960)" = "WEIBULL requires GAMMA > 0",
  "(SETG3970)" = "WEIBULL requires ETA > 0",
  "(SETG3980)" = "WEIBULL rejects a negative ALPHA, and rejects ALPHA = 0 unless ALPHA is fixed",
  "(SETG31020)" = "GAMMA*ETA <= 2 with both GAMMA and ETA fixed, so neither can be adjusted",
  "(SETG31040)" = "GAMMA*ETA/ALPHA <= 2 with ALPHA fixed, so ALPHA cannot be adjusted",
  "(SETG31070)" = "ALPHA is fixed where SETG3 must derive it from GAMMA*ETA",
  "(SETG31090)" = "GAMMA is fixed but non-positive, so SETG3 cannot derive one from ALPHA and ETA",
  "(SETG32020)" = "ETA is fixed but non-positive, so SETG3 cannot derive one from ALPHA and GAMMA",
  "(SETG32050)" = "GAMMA is fixed but non-positive, so SETG3 cannot derive one from ETA",
  "(SETG32080)" = "ETA is fixed but non-positive, so SETG3 cannot derive one from GAMMA",
  "(SETG33010)" = "ETA is fixed but non-positive, with no GAMMA to derive one from",
  "(SETG33020)" = "GAMMA is fixed but non-positive, with no ETA to derive one from"
)
# Deliberately absent, so their omission is not read as an oversight:
#   SETG3950  (setg3.c:319-321) fires when 2*Tmax/3 itself is <= 0. That is a
#             property of the data, not of the PARMS block, so it cannot be
#             decided here -- the TAU row already says the start is
#             data-dependent.
#   SETG3940, SETG3990, SETG31000, SETG31010
#             are reachable only through g_two/ga_two, which FIXGE2/FIXGAE2
#             drive and .hzr_sas_token() records as unresolved.
# SEVEN of the sixteen are unreachable in both languages, because each guards a
# condition the ENTRY checks at setg3.c:269-284 have already refused:
#   SETG31090, SETG32050, SETG33020  want a non-positive GAMMA that is fixed
#                                    -- SETG3910 refuses first.
#   SETG32020, SETG32080, SETG33010  want a non-positive ETA that is fixed
#                                    -- SETG3930 refuses first.
#   SETG31070                        wants a fixed ALPHA on a branch where
#                                    alpha <= 0: negative refuses at SETG3920,
#                                    and zero sets g3flag = 2 (:332-335) so
#                                    SETG3_alpha_gener() is never called.
# They are kept because the C keeps them and because the entry checks are what
# makes them dead -- change those and these wake up. Nine codes can actually
# fire, which an exhaustive search in the tests pins rather than asserts.

#' Plain-language gloss for a SETG3 refusal code.
#' @noRd
.hzr_setg3_refusal_reason <- function(code) {
  # `[[` on an unmatched name in a CHARACTER vector throws "subscript out of
  # bounds" rather than returning NULL, so an is.null() fallback here would be
  # hollow -- right shape, never reachable, and the parser would error instead
  # of degrading if a code were ever added to the trace but not the table.
  hit <- .hzr_setg3_refusal[code]
  if (is.na(hit)) "the operand combination is rejected" else unname(hit)
}

#' `alpha` handling shared by several SETG3 branches (setg3.c:815-852).
#' @return `NULL`, a refusal string, or the substituted alpha.
#' @noRd
.hzr_setg3_alpha_fixup <- function(gamma, eta, alpha, alpha_fixed, weibull) {
  # ga_two is FALSE, so the C falls to the `gteva > TWO || weibul` early
  # return; under WEIBULL this function is a no-op entirely.
  if (weibull) return(NULL)
  if (isTRUE(gamma * eta / alpha > 2)) return(NULL)
  if (alpha_fixed) return("(SETG31040)")
  gamma * eta / 3
}

#' `alpha` handling for the branches that reach SETG3_alpha_gener()
#' (setg3.c:854-871).
#' @noRd
.hzr_setg3_alpha_gener <- function(gamma, eta, alpha_fixed) {
  if (alpha_fixed) return("(SETG31070)")
  gamma * eta / 3
}

#' Walk `SETG3()` and report what it would do to one late phase.
#'
#' @param tau_raw,gamma,alpha,eta The operand values `PARMS` supplied, with
#'   PROC HAZARD's own initializers for any it did not
#'   (`src/hazard/stmtprc.c:34-37`). `tau_raw` is `NA` when `TAU` was absent:
#'   the `stmtprc.c` initializer of 0 does NOT survive to `SETG3()`, because
#'   `src/hazard/readobs.c:153-154` replaces an unspecified `TAU` with
#'   `0.75 * Tmax` on an active late phase, and `readobs()` runs at
#'   `hazard.c:276`, before `hzrg()` reaches `SETG3()` at `:292`. So an absent
#'   `TAU` arrives POSITIVE and cannot raise `SETG3900`.
#' @param fixed Character vector of parameters the job's `FIX*` tokens pinned,
#'   as the user wrote them: these checks run before SETG3 changes any flag.
#' @param weibull Whether the job carries the bare `WEIBULL` keyword.
#' @return `list(refusal = <chr or NULL>, shape = <named numeric or NULL>)`.
#'   `refusal` names the SAS message code for a job PROC HAZARD will not run;
#'   otherwise `shape` gives the values SETG3 would optimize from.
#' @noRd
.hzr_setg3_notes <- function(tau_raw, gamma, alpha, eta, fixed, weibull) {
  fx <- function(p) p %in% fixed
  # `entry` marks the four checks at setg3.c:269-284, the only ones raised
  # BEFORE the TAU rules at :309-323. The other twelve codes fire after those
  # rules have already run, which decides whether a TAU row still applies.
  refuse <- function(code, entry = FALSE) {
    list(refusal = code, entry = entry, shape = NULL)
  }

  # setg3.c:269-284. A FIX* on an operand SAS reads as unspecified is fatal,
  # and it is checked before anything else -- including SETG3_ignore_tau().
  # An absent TAU (NA) is 0.75*Tmax by the time SETG3 sees it -- positive, so
  # only an explicitly non-positive TAU= can refuse here.
  if (isTRUE(tau_raw <= 0) && fx("tau")) return(refuse("(SETG3900)", TRUE))
  if (isTRUE(gamma <= 0) && fx("gamma")) return(refuse("(SETG3910)", TRUE))
  if (isTRUE(alpha < 0) && fx("alpha")) return(refuse("(SETG3920)", TRUE))
  if (isTRUE(eta <= 0) && fx("eta")) return(refuse("(SETG3930)", TRUE))

  # setg3.c:313-315 and 403-421. Reproduced here for the trace only: the
  # emitted call deliberately keeps the user's GAMMA and ETA, because the
  # swap is likelihood-equivalent. The trace needs SAS's values because the
  # dispatch below keys on their signs.
  if (isTRUE(alpha == 1) && fx("alpha")) {
    if (!fx("gamma") && !fx("eta")) fixed <- union(fixed, "eta")
    product <- gamma * eta
    if (fx("eta")) {
      gamma <- if (isTRUE(product > 0)) product else eta
      eta <- 1
    } else {
      eta <- if (isTRUE(product > 0)) product else gamma
      gamma <- 1
    }
    fixed <- union(fixed, c("tau", "alpha"))
  }

  # setg3.c:332-335, with ga_two FALSE.
  g3flag <- if (isTRUE(alpha == 0) && fx("alpha")) 2L else 1L

  if (weibull) {
    # setg3.c:427-482, then the return at :347. The g_two block is skipped
    # (see the header), and SETG3_alpha_fixup() is a no-op under WEIBULL.
    if (isTRUE(gamma <= 0)) return(refuse("(SETG3960)"))
    if (isTRUE(eta <= 0)) return(refuse("(SETG3970)"))
    if (isTRUE(alpha < 0) || (isTRUE(alpha == 0) && g3flag + 2L == 3L)) {
      return(refuse("(SETG3980)"))
    }
    return(list(refusal = NULL, entry = FALSE,
                shape = c(gamma = gamma, alpha = alpha, eta = eta)))
  }

  # setg3.c:873-924 with g_two FALSE: push GAMMA*ETA clear of 2.
  verify_ge_2 <- function() {
    if (!isTRUE(gamma * eta <= 2)) return(NULL)
    if (fx("gamma") && fx("eta")) return("(SETG31020)")
    if (fx("gamma")) eta <<- 3 / gamma else gamma <<- 3 / eta
    NULL
  }

  # setg3.c:359-374. Each branch is the !g_two && !ga_two row of its table.
  bad <- NULL
  if (isTRUE(alpha > 0) && isTRUE(gamma > 0) && isTRUE(eta > 0)) {
    bad <- verify_ge_2()
    if (is.null(bad)) {
      a <- .hzr_setg3_alpha_fixup(gamma, eta, alpha, fx("alpha"), weibull)
      if (is.character(a)) bad <- a else if (!is.null(a)) alpha <- a
    }
  } else if (isTRUE(gamma > 0) && isTRUE(eta > 0)) {
    bad <- verify_ge_2()
    if (is.null(bad) && g3flag == 1L) {
      a <- .hzr_setg3_alpha_gener(gamma, eta, fx("alpha"))
      if (is.character(a)) bad <- a else alpha <- a
    }
  } else if (isTRUE(alpha > 0) && isTRUE(eta > 0)) {
    # gamma <= 0 (setg3.c:534-585)
    if (fx("gamma")) return(refuse("(SETG31090)"))
    gamma <- 3 * alpha / eta
    if (isTRUE(gamma * eta <= 2)) gamma <- 3 / eta
    a <- .hzr_setg3_alpha_fixup(gamma, eta, alpha, fx("alpha"), weibull)
    if (is.character(a)) bad <- a else if (!is.null(a)) alpha <- a
  } else if (isTRUE(alpha > 0) && isTRUE(gamma > 0)) {
    # eta <= 0 (setg3.c:588-636)
    if (fx("eta")) return(refuse("(SETG32020)"))
    eta <- 3 * alpha / gamma
    if (isTRUE(gamma * eta <= 2)) eta <- 3 / gamma
    a <- .hzr_setg3_alpha_fixup(gamma, eta, alpha, fx("alpha"), weibull)
    if (is.character(a)) bad <- a else if (!is.null(a)) alpha <- a
  } else if (isTRUE(eta > 0)) {
    # alpha <= 0, gamma <= 0 (setg3.c:638-675)
    if (fx("gamma")) return(refuse("(SETG32050)"))
    gamma <- 3 / eta
    if (g3flag == 1L) {
      a <- .hzr_setg3_alpha_gener(gamma, eta, fx("alpha"))
      if (is.character(a)) bad <- a else alpha <- a
    }
  } else if (isTRUE(gamma > 0)) {
    # alpha <= 0, eta <= 0 (setg3.c:677-714)
    if (fx("eta")) return(refuse("(SETG32080)"))
    eta <- 3 / gamma
    if (g3flag == 1L) {
      a <- .hzr_setg3_alpha_gener(gamma, eta, fx("alpha"))
      if (is.character(a)) bad <- a else alpha <- a
    }
  } else if (isTRUE(alpha > 0)) {
    # gamma <= 0, eta <= 0 (setg3.c:717-770)
    if (fx("gamma")) return(refuse("(SETG33020)"))
    if (fx("eta")) return(refuse("(SETG33010)"))
    eta <- 2
    gamma <- 1.5 * alpha
    if (isTRUE(gamma * eta <= 2)) gamma <- 3 / eta
    a <- .hzr_setg3_alpha_fixup(gamma, eta, alpha, fx("alpha"), weibull)
    if (is.character(a)) bad <- a else if (!is.null(a)) alpha <- a
  } else {
    # all <= 0 (setg3.c:772-812)
    if (fx("gamma")) return(refuse("(SETG33020)"))
    if (fx("eta")) return(refuse("(SETG33010)"))
    gamma <- 1
    eta <- 3
    if (g3flag == 1L) {
      a <- .hzr_setg3_alpha_gener(gamma, eta, fx("alpha"))
      if (is.character(a)) bad <- a else alpha <- a
    }
  }

  if (!is.null(bad)) return(refuse(bad))
  list(refusal = NULL, entry = FALSE,
       shape = c(gamma = gamma, alpha = alpha, eta = eta))
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
  saw_weibull <- FALSE
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
          # sprintf("%g"), not format(): this string is DATA, not just a
          # message -- it lands in the untranslated frame and is grepped by
          # callers and tests. format() honours getOption("OutDec"), so a
          # session with OutDec = "," would write "DELTA = 0,5" and break both.
          flag_bad(op, paste0(
            "DELTA = ", sprintf("%g", val), " is not implemented -- R assumes ",
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
      # setopt(6) -> SETG3_weibull() (setg3.c:427) is the GENERALIZED Weibull:
      # "NOW HANDLE THE SPECIAL SITUATION OF THE GENERALIZED WEIBULL, WHERE WE
      # ADMIT ALL POSITIVE VALUES OF THE PARAMETERS." It bumps g3flag and
      # validates gamma > 0, eta > 0, alpha >= 0. It assigns nothing and fixes
      # nothing, so neither does this branch -- ALPHA and ETA keep whatever
      # PARMS specified and stay free unless an explicit FIXALPHA/FIXETA pins
      # them. Overwriting them with 1 discarded the user's starting values and
      # silently fitted a smaller model; the listing for
      # hz.ce_cardioversion_repeated.ehb.sas prints both as "Estimated? Yes".
      #
      # g3flag itself has no R counterpart: it selects a numerical branch, and
      # hzr_decompos_g3() handles the general form directly. The GAMMA*ETA = 2
      # and GAMMA*ETA/ALPHA = 2 constraint flags SETG3_weibull() also honours
      # are driven by separate PARMS keywords that this parser does not yet
      # resolve -- they are recorded as untranslated, not assumed absent.
      saw_weibull <- TRUE
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

  # Every shape operand PARMS did not name is filled from PROC HAZARD's own
  # defaults, once, here -- and the same filled list then builds both the
  # hzr_phase() call and that phase's theta block, so the two cannot describe
  # different starting values. A partially specified PARMS block therefore
  # emits a call naming every shape argument explicitly, which is also what
  # makes the divergence readable: the defaulted values are on the page rather
  # than hidden in hzr_phase()'s formals.
  #
  # TAU is the exception, because SETG3() derives it from the data. Its
  # predicate is the C one -- tau <= 0, which catches an explicit TAU=0 as well
  # as an absent TAU (setg3.c:316) -- and it is a divergence only when
  # SETG3_ignore_tau() does not take the branch above it, since that branch
  # pins TAU at 1 (setg3.c:378), exactly the value emitted here.
  # `late` itself is left untouched: `length(late)` is the "did PARMS name a
  # late shape operand" gate below, and a TAU=0 operand still counts as one.
  #
  # Absent and explicitly-non-positive are DIFFERENT cases, and conflating
  # them put a false refusal in $untranslated. PROC HAZARD never sees the
  # stmtprc.c initializer of 0 for an absent TAU: readobs.c:153-154 replaces
  # it with 0.75*Tmax on an active late phase, before hzrg() reaches SETG3()
  # (hazard.c:276 against :292). Only a TAU= the job actually wrote can still
  # be <= 0 at setg3.c:269, and only that one reaches the 2*Tmax/3 assignment
  # at :317. Both defaults are data-dependent; they are different data.
  tau_absent <- is.null(late[["tau"]])
  tau_nonpositive <- !tau_absent && !isTRUE(late[["tau"]] > 0)
  tau_defaulted <- tau_absent || tau_nonpositive
  early_full <- .hzr_parms_fill_shape(early, .hzr_parms_early_sas_default)
  late_full <- .hzr_parms_fill_shape(late, .hzr_parms_late_sas_default)

  # A phase is active in PROC HAZARD iff its MU was specified and positive
  # (parmprc.c:19 -> setparmno.c:14); with phase 3 off, shape.c never calls
  # SETG3() and none of what follows happens. This is the gate every SETG3
  # guard below shares.
  has_late <- !is.null(mu[["MUL"]]) && isTRUE(mu[["MUL"]] > 0)
  # setg3.c:312-314's own predicate for taking the SETG3_ignore_tau() branch.
  # The `g_two && ga_two` disjunct alongside it is driven by PARMS keywords
  # this parser does not resolve; those are recorded as untranslated, so this
  # is the half that can be evaluated here.
  ignore_tau <- isTRUE(late_full[["alpha"]] == 1) && "alpha" %in% fixed_late

  # SETG3's entry checks (setg3.c:269-284) read the job's OWN FIX* flags,
  # before SETG3 sets any of its own, so the trace below needs this snapshot
  # rather than the post-pin set.
  fixed_late_user <- fixed_late

  if (has_late && ignore_tau) {
    # setg3.c:378-379 does TWO things, and mirroring only the first is what
    # made this branch look harmless: it sets TAU to 1 AND fixes it. Leaving
    # TAU free is not a spare degree of freedom here but an unidentified one --
    # at alpha = 1 the G3 form collapses to (t/tau)^(gamma*eta), so mu enters
    # the likelihood only through log_mu - gamma*eta*log_tau and the two are
    # exactly aliased. The emitted fit converged onto that ridge and returned
    # no standard errors for either parameter.
    #
    # So the value AND the fix are both mirrored, which is what the C does and
    # what makes the emitted call identifiable. Where PARMS named a different
    # TAU this overrides the user's own text, so that case is recorded below --
    # PROC HAZARD overrides it too, but a silent rewrite of what the job said
    # is exactly what $untranslated exists to surface.
    late_full[["tau"]] <- 1
    pins <- "tau"
    # setg3.c:406-407: with GAMMA and ETA both free, SETG3_ignore_tau() fixes ETA.
    # Mirrored for the same reason TAU is, and NOT for the reason
    # SETG3_verify_ge_2() is declined below: this one is exact algebra, not a
    # numerical-branch restriction. At alpha = 1 only the product gamma*eta
    # enters the likelihood, so the pair is exactly singular in R too -- the
    # general G3 shape hzr_phase() carries is genuinely degenerate here, and
    # leaving both free relocates the flat ridge the TAU pin removed rather
    # than exercising a shape SAS cannot reach.
    #
    # The VALUE rewrite that follows in the C (setg3.c:411-415: gamma <-
    # gamma*eta, eta <- 1) is still not mirrored, per the decision recorded
    # above: it is likelihood-equivalent -- both parameterisations span the
    # same exponent -- so it changes what is reported, not what is fitted, and
    # rewriting it would put the emitted call at odds with the user's PARMS.
    if (!("gamma" %in% fixed_late) && !("eta" %in% fixed_late)) {
      pins <- c(pins, "eta")
    }
    fixed_late <- intersect(unname(.hzr_parms_late_arg),
                            union(fixed_late, pins))
  } else if (tau_defaulted) {
    late_full[["tau"]] <- .hzr_parms_late_sas_default[["tau"]]
  }

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
      "cdf", early_full, phase_covars$early, fixed_early
    )
    mu_val <- if (!is.null(mu[["MUE"]])) mu[["MUE"]] else .hzr_parms_default_mu_start
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("early", mu_val, early_full, phase_covar_vals$early)
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
      "g3", late_full, phase_covars$late, fixed_late
    )
    mu_val <- if (!is.null(mu[["MUL"]])) mu[["MUL"]] else .hzr_parms_default_mu_start
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("late", mu_val, late_full, phase_covar_vals$late)
    )
  }

  # A MU names its phase. Building a phase only when a *shape* operand
  # appeared let an orphaned MUE/MUL disappear together with the phase it
  # scaled -- a one-phase R model against SAS's two, and no untranslated row
  # to say so. The shape defaults are now PROC HAZARD's own (src/hazard/stmtprc.c:34-37,
  # see .hzr_parms_early_sas_default), but a phase reconstructed from nothing
  # but a MU would still start from a scale SAS does not use, and its tau would
  # still be data-dependent, so the MU stays recorded rather than guessed at.
  # SETG3_ignore_tau() (setg3.c:377-424) fires when ALPHA is fixed at 1 --
  # setg3.c:312-314, before the WEIBULL dispatch. It pins TAU at 1 and rewrites
  # the GAMMA/ETA split, which at alpha = 1, tau = 1 is a reparameterisation of
  # a single exponent: the G3 form collapses to t^(gamma*eta), so only the
  # product is identified. Two consequences, and only one of them is a defect.
  #
  # The reported values diverge: SAS's listing prints the rewritten split, this
  # parser emits what PARMS said. The fit is identical, so that is documented
  # in hzr_translate_sas() rather than mirrored here -- rewriting the emitted
  # call would make it disagree with the user's own PARMS text, and across the
  # public corpus it would change 16 blocks to alter 2 (measured 2026-09-08).
  #
  # With GAMMA and ETA *both* free, though, setg3.c:405-406 fixes ETA and
  # estimates the product as GAMMA alone. That is one fewer estimated parameter
  # than hazard() would use, on a pair that is not separately identifiable --
  # a different model, not a different label for the same one. Recorded.
  # Both guards describe things SETG3() does, so both require a late phase that
  # is actually active -- and in PROC HAZARD a phase is active iff its MU was
  # specified and positive. parmprc.c:19 registers MUL through
  # setparmno(22, 7, 3, ...) and setparmno.c:14 sets C->phase[3] = 1 only when
  # stmtfld() > 0; the four shape operands go through setprmf
  # (parmprc.c:20-23), which never touches C->phase[]. With phase 3 off,
  # stmtprc.c:113-122 zeroes TAU/GAMMA/ALPHA/ETA and shape.c:31 never calls
  # SETG3() at all.
  #
  # So the gate is MUL, not the presence of shape operands: TAU/GAMMA/ETA with
  # no MUL is not a late phase, and warning about SETG3_ignore_tau() or
  # SETG3980 there describes code PROC HAZARD never reaches. A false positive
  # on $untranslated is not harmless -- that frame is how a caller decides
  # whether a translation can be trusted.
  #
  # Note the asymmetry this leaves inside this function, deliberately and
  # tracked separately: the phase-building block below still keys on shape
  # operands (`if (length(late))`), so PARMS carrying TAU/GAMMA with no MUL
  # emits a late phase that PROC HAZARD would not build at all. No corpus job
  # does that, so it is latent, and changing how phases are built is a wider
  # change than gating these two warnings.
  # The GAMMA*ETA divergence this block used to record is gone: setg3.c:405's
  # ETA fix is now mirrored above, so hazard() estimates the same number of
  # free parameters as PROC HAZARD on this branch. What remains unmirrored is
  # the value rewrite (gamma <- gamma*eta, eta <- 1), which is
  # likelihood-equivalent and is documented in hzr_translate_sas() rather than
  # reported per job.

  # Everything SETG3() would do to this phase, in one place. See the header
  # of .hzr_setg3_notes() for why a refusal and a rewrite are both recorded
  # while only the two identifiability fixes are mirrored.
  setg3_refused <- FALSE
  if (length(late) && has_late) {
    setg3 <- .hzr_setg3_notes(
      tau_raw = if (tau_absent) NA_real_ else late[["tau"]],
      gamma = late_full[["gamma"]],
      alpha = late_full[["alpha"]],
      eta = late_full[["eta"]],
      fixed = fixed_late_user,
      weibull = saw_weibull
    )
    if (!is.null(setg3$refusal)) {
      setg3_refused <- isTRUE(setg3$entry)
      flag_bad(
        sprintf("GAMMA=%g ALPHA=%g ETA=%g%s", late_full[["gamma"]],
                late_full[["alpha"]], late_full[["eta"]],
                if (length(fixed_late_user)) {
                  paste0(" fixed:", paste(fixed_late_user, collapse = ","))
                } else {
                  ""
                }),
        paste0("PROC HAZARD refuses this job: SETG3 raises ",
               setg3$refusal, " -- ",
               .hzr_setg3_refusal_reason(setg3$refusal),
               ". hzr_phase() would accept it, so without this the ",
               "translation would emit a runnable fit for a job that does ",
               "not run")
      )
    } else {
      # At alpha = 1 only the product gamma*eta is identified, and the
      # emitted call deliberately keeps the user's split rather than
      # SETG3_ignore_tau()'s (setg3.c:411-415) -- likelihood-equivalent, so
      # compare the product there and each operand otherwise. Without this
      # every ignore_tau job would report a rewrite that changes no fit.
      #
      # The swap is product-preserving only while that product is POSITIVE.
      # setg3.c:411-412 and :416-417 carry a fallback -- when gamma*eta <= 0
      # the C keeps the OTHER operand instead, so gamma = -2, eta = 3 leaves
      # SAS at gamma = 3, eta = 1, a product of 3 rather than -6. Comparing
      # products still catches that (the two differ), which is the point: the
      # comparison must not be read as an assertion that the swap is always
      # harmless. Such a phase is emitted with a non-positive gamma anyway,
      # which hzr_phase() rejects outright.
      emitted <- c(gamma = late_full[["gamma"]], alpha = late_full[["alpha"]],
                   eta = late_full[["eta"]])
      got <- setg3$shape
      if (ignore_tau) {
        emitted <- c(`gamma*eta` = unname(emitted[["gamma"]] * emitted[["eta"]]),
                     alpha = unname(emitted[["alpha"]]))
        got <- c(`gamma*eta` = unname(got[["gamma"]] * got[["eta"]]),
                 alpha = unname(got[["alpha"]]))
      }
      moved <- names(emitted)[!mapply(
        function(a, b) isTRUE(all.equal(a, b)), emitted, got
      )]
      if (length(moved)) {
        # Where the operand is outside hzr_phase()'s own domain the "general
        # shape kept deliberately" framing is simply false: hzr_phase()
        # stopifnot()s gamma > 0, eta > 0 and alpha >= 0 (R/phase-spec.R), so
        # the emitted call cannot even be built. This is the one direction in
        # which PROC HAZARD is the MORE permissive of the two -- it reads a
        # non-positive operand as "derive one for me" -- and saying "a parity
        # run starts elsewhere" would tell the reader the difference is a
        # starting value when the difference is that R will not run this.
        unbuildable <- !isTRUE(late_full[["gamma"]] > 0) ||
          !isTRUE(late_full[["eta"]] > 0) ||
          !isTRUE(late_full[["alpha"]] >= 0)
        flag_bad(
          paste(sprintf("%s=%g", names(emitted), emitted), collapse = " "),
          paste0("SETG3() optimizes from ",
                 paste(sprintf("%s = %g", moved, got[moved]), collapse = ", "),
                 ", not the value(s) emitted here: it ",
                 if (unbuildable) {
                   paste0("reads a non-positive GAMMA/ETA (or negative ALPHA) ",
                          "as a request to derive one. hzr_phase() requires ",
                          "gamma > 0, eta > 0 and alpha >= 0, so the emitted ",
                          "call cannot be built at all -- supply the operand ",
                          "rather than expecting a starting-value difference")
                 } else {
                   paste0("constrains the late shape to keep PROC HAZARD ",
                          "inside a numerical branch it can evaluate. ",
                          "hzr_decompos_g3() carries the general G3 form and ",
                          "needs no such constraint, so the emitted call ",
                          "keeps it deliberately -- but a SAS parity run ",
                          "starts elsewhere")
                 })
        )
      }
    }
  }

  # The emitted call now pins TAU at 1 exactly as SETG3_ignore_tau() does, so
  # an unspecified (or already-1) TAU translates faithfully and needs no row.
  # A job that WROTE a different TAU is a different matter: PROC HAZARD
  # discards that value, and so now does this translator, but a starting value
  # the user typed and neither program uses is worth saying out loud.
  # Only the four ENTRY refusals (setg3.c:269-284) are raised before the TAU
  # rules at :309-323; the other twelve fire from :429 onwards, by which point
  # PROC HAZARD has already applied whichever TAU rule was going to apply. So
  # the suppression below is keyed on `setg3$entry`, not on "refused at all" --
  # reporting a TAU rule alongside an entry refusal would describe code PROC
  # HAZARD never reaches, the same false positive the MUL gate exists to avoid,
  # while suppressing it for a late refusal would hide one that did run.
  if (length(late) && has_late && !setg3_refused && ignore_tau &&
      !is.null(late[["tau"]]) && !isTRUE(late[["tau"]] == 1)) {
    flag_bad(
      paste0("TAU=", sprintf("%g", late[["tau"]])),
      paste0("SETG3_ignore_tau() discards this TAU and runs at tau = 1, ",
             "fixed (setg3.c:378-379), because ALPHA is fixed at 1; the ",
             "emitted phase mirrors that, so the value written here is used ",
             "by neither PROC HAZARD nor the translation")
    )
  } else if (length(late) && has_late && !setg3_refused && tau_defaulted &&
             !ignore_tau) {
    # Which data-dependent default applies depends on whether the job wrote a
    # TAU at all -- see the tau_absent/tau_nonpositive split above.
    # The other SETG3 branch (setg3.c:316-318), reached only when the
    # ignore_tau branch above was NOT taken -- setg3.c:313-316 is an if/else,
    # so a job that fixes ALPHA at 1 never gets this assignment at all.
    # A non-positive TAU -- absent from PARMS, or written as TAU=0 -- is
    # replaced by 2*Tmax/3. It is the one
    # shape default that cannot be reproduced at parse time, because it depends
    # on the data, so the emitted call starts at tau = 1 and this says so
    # rather than letting the difference pass as a translation. Emitting
    # `2 * max(<timevar>) / 3` instead was considered and rejected: SAS's Tmax
    # is taken over the analysis set after exclusions, not over the raw column,
    # so the expression would look exact while being a guess.
    #
    # `length(late)` because a late phase that was never built is already
    # reported by the MUL guard below, and two rows for one absence is noise.
    flag_bad(
      if (tau_absent) "TAU (unspecified)" else
        paste0("TAU=", sprintf("%g", late[["tau"]])),
      paste0("PROC HAZARD starts TAU at ",
             if (tau_absent) {
               "0.75*Tmax (readobs.c:153-154, applied to an unspecified TAU "
             } else {
               "2*Tmax/3 (setg3.c:317, applied to a non-positive TAU "
             },
             "before SETG3 runs), which depends on the data and cannot be ",
             "reproduced at parse time; the emitted phase starts at tau = 1, ",
             "so this fit begins somewhere PROC HAZARD would not and the ",
             "multiphase likelihood is multimodal")
    )
  }

  # sprintf("%g"), not format(): format() honours getOption("OutDec"), so a
  # session with OutDec = "," would record "MUE=0,2" and break every grep --
  # the same trap the DELTA reason string above avoids.
  if (!is.null(mu[["MUE"]]) && !length(early)) {
    flag_bad(paste0("MUE=", sprintf("%g", mu[["MUE"]])),
             "MUE with no early phase shape operand (THALF/NU/M)")
  }
  if (!is.null(mu[["MUL"]]) && !length(late)) {
    flag_bad(paste0("MUL=", sprintf("%g", mu[["MUL"]])),
             "MUL with no late phase shape operand (TAU/GAMMA/ALPHA/ETA)")
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
