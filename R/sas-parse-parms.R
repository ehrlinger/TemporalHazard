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

# The grammar table (.hzr_sas_grammar) is generated from HAZARD's own lexer
# (data-raw/hazard-grammar.R), so a PARMS keyword it does not know is one
# PROC HAZARD's lexer rejects: the job does not run. FIXG1 and FIXG3, for
# instance, are internal flags shape.c:36-41 sets, not options. The prefix is
# kept for callers that grep it.
.hzr_parms_unresolved_reason <- paste0(
  "unresolved PARMS keyword: not in PROC HAZARD's grammar (hazard_l.l), so ",
  "PROC HAZARD rejects this job with a syntax error and it does not run; ",
  "whatever is emitted here translates a job that does not run"
)
# A macro reference (`&EXTRA`) is resolved by SAS before PROC HAZARD reads
# the statement, so the grammar table cannot say what it becomes, and in
# particular cannot say the job fails. Same prefix, neutral consequence.
.hzr_parms_unresolved_macro_reason <- paste0(
  "unresolved PARMS keyword: a SAS macro reference, which SAS resolves ",
  "before PROC HAZARD reads the statement, so this translation cannot tell ",
  "what it becomes"
)
# A SAS macro reference (`&X`) or call (`%CALL`), which SAS expands before
# PROC HAZARD reads the statement: never judged as a syntax error.
.hzr_sas_is_macro <- function(x) grepl("[&%][A-Za-z_]", x)

.hzr_parms_unresolved_why <- function(op) {
  if (.hzr_sas_is_macro(op)) .hzr_parms_unresolved_macro_reason else
    .hzr_parms_unresolved_reason
}

# The lexer's NUMBER (hazard_l.l:34-38). as.numeric() also reads 1E-3 (no
# decimal point before the exponent), 2. and +0.2, none of which PROC HAZARD
# lexes as a number: the job stops with a syntax error.
.hzr_sas_lexer_number <- function(s) {
  grepl("^-?([0-9]+|[0-9]*[.][0-9]+(E[+-]?[0-9]+)?)$", toupper(s))
}

# An operand written with spaces around `=` (`THALF = 0.3`, `THALF =0.3`,
# `THALF= 0.3`) reaches this parser split into pieces. SAS's lexer skips the
# whitespace (hazard_l.l:50) and PROC HAZARD runs the job, so no piece is a
# syntax error; but the operand's value is not read here. Marked from context,
# since a bare number is a piece only after an `=`, and a stray one is not.
.hzr_parms_rejected_piece_reason <- paste0(
  "unresolved PARMS keyword: a piece of an operand written with spaces ",
  "around `=` that is not a value keyword `= NUMBER` even joined ",
  "(hazard_y.y:137-147, hazard_l.l:34-38), so PROC HAZARD rejects this job ",
  "with a syntax error and it does not run"
)
.hzr_parms_macro_piece_reason <- paste0(
  "unresolved PARMS keyword: a piece of an operand written with spaces ",
  "around `=` whose key or value is a SAS macro reference, which SAS ",
  "resolves before PROC HAZARD reads the statement, so this translation ",
  "cannot tell what the operand becomes"
)
.hzr_parms_unresolved_piece_reason <- paste0(
  "unresolved PARMS keyword: a piece of an operand written with spaces ",
  "around `=`, which PROC HAZARD accepts but this translator splits apart, ",
  "so the operand's value was not read"
)
# SAS's lexer skips whitespace (hazard_l.l:32), so `THALF = 0.3` is the same
# operand as `THALF=0.3`: PROC HAZARD reads the value and runs the job. This
# parser splits the statement on whitespace, so it joins the pieces back
# before reading them. Anything that does not join into `KEY=VALUE` is left
# alone for .hzr_parms_spaced_pieces() to judge (#421).
#' Rejoin a macro call whose arguments contain spaces.
#'
#' Operands are split on whitespace, which cuts `%FLAGS(A, B)` into
#' `%FLAGS(A,` and `B)`. Only the first fragment then looks like a macro, and
#' the remainder was classified on its own -- so a job SAS runs collected a
#' false `$untranslated` row and, since 2026-09-22, a false warning (#433
#' review).
#'
#' SAS expands the whole call before `PROC HAZARD` sees any operand, so the
#' call must travel as one token and stay indeterminate.
#'
#' Absorption is bounded by the operands present: an unclosed `%FLAGS(A`
#' takes the rest and stops, rather than looping. That is the right reading
#' anyway, since everything after it is inside the unterminated call.
#' @noRd
.hzr_sas_join_macro_calls <- function(ops) {
  n <- length(ops)
  if (n < 2L) return(ops)
  opens <- function(x) lengths(regmatches(x, gregexpr("(", x, fixed = TRUE)))
  closes <- function(x) lengths(regmatches(x, gregexpr(")", x, fixed = TRUE)))
  out <- character(0)
  i <- 1L
  while (i <= n) {
    op <- ops[[i]]
    if (.hzr_sas_is_macro(op) && opens(op) > closes(op)) {
      j <- i
      acc <- op
      while (j < n && opens(acc) > closes(acc)) {
        j <- j + 1L
        acc <- paste(acc, ops[[j]])
      }
      out <- c(out, acc)
      i <- j + 1L
    } else {
      out <- c(out, op)
      i <- i + 1L
    }
  }
  out
}

.hzr_sas_join_spaced <- function(ops) {
  n <- length(ops)
  if (n < 2L) return(ops)
  out <- character(0)
  i <- 1L
  while (i <= n) {
    op <- ops[[i]]
    if (identical(op, "=") && length(out) && i < n) {
      out[length(out)] <- paste0(out[length(out)], "=", ops[[i + 1L]])
      i <- i + 2L
    } else if (identical(op, "=") && length(out)) {
      # A TRAILING bare `=`: the dangling half of `MAXITER =` with nothing
      # after it. Attaching it to the previous operand makes one construct
      # `MAXITER=` rather than leaving an operand whose key is the empty
      # string, which was then reported as an unknown option with a blank
      # name (#433 review).
      out[length(out)] <- paste0(out[length(out)], "=")
      i <- i + 1L
    } else if (nchar(op) > 1L && endsWith(op, "=") && i < n) {
      out <- c(out, paste0(op, ops[[i + 1L]]))
      i <- i + 2L
    } else if (nchar(op) > 1L && startsWith(op, "=") && length(out)) {
      out[length(out)] <- paste0(out[length(out)], op)
      i <- i + 1L
    } else {
      out <- c(out, op)
      i <- i + 1L
    }
  }
  out
}

.hzr_parms_spaced_pieces <- function(ops) {
  # 0 = not a piece; 1 = a piece of a spaced operand PROC HAZARD accepts
  # (joined, it is a value keyword `= NUMBER`, hazard_y.y:137-147 and
  # hazard_l.l:34-38); 2 = a piece of one it would still reject; 3 = a piece
  # of one whose key or value is a macro reference, which SAS expands before
  # PROC HAZARD reads it (`&KEY = 0.3` may be THALF = 0.3), so it can be
  # judged neither way here (Codex on #365).
  n <- length(ops)
  code <- integer(n)
  tok <- function(x) .hzr_sas_token(x, "HAZARD", "PARM")
  value_key <- function(x) {
    t <- tok(x)
    !is.na(t) && t %in% c(.hzr_parms_mu_order, names(.hzr_parms_early_arg),
                          names(.hzr_parms_late_arg), "DELTA")
  }
  macro <- .hzr_sas_is_macro
  bare_key <- function(k) {
    k >= 1L && code[k] == 0L && !grepl("=", ops[[k]], fixed = TRUE) &&
      (!is.na(tok(ops[[k]])) || macro(ops[[k]]))
  }
  i <- 1L
  while (i <= n) {
    op <- ops[[i]]
    key_i <- if (bare_key(i - 1L)) i - 1L else NA_integer_
    if (identical(op, "=")) {
      val_i <- if (i < n) i + 1L else NA_integer_
      ok <- !is.na(key_i) && value_key(ops[[key_i]]) && !is.na(val_i) &&
        .hzr_sas_lexer_number(ops[[val_i]])
      members <- c(key_i, i, val_i)
    } else if (startsWith(op, "=") && nchar(op) > 1L) {
      ok <- !is.na(key_i) && value_key(ops[[key_i]]) &&
        .hzr_sas_lexer_number(substring(op, 2L))
      members <- c(key_i, i)
    } else if (nchar(op) > 1L && endsWith(op, "=") && i < n) {
      ok <- value_key(substr(op, 1L, nchar(op) - 1L)) &&
        .hzr_sas_lexer_number(ops[[i + 1L]])
      members <- c(i, i + 1L)
    } else {
      i <- i + 1L
      next
    }
    members <- members[!is.na(members)]
    code[members] <- if (any(macro(ops[members]))) 3L else if (ok) 1L else 2L
    i <- max(members) + 1L
  }
  code
}
.hzr_parms_mu_order  <- c("MUE", "MUC", "MUL")

# PROC HAZARD's OWN shape defaults (src/hazard/stmtprc.c:34-37), used for any
# shape operand a PARMS statement did not name. They are not hzr_phase()'s
# defaults -- nu, m and eta all differ (SAS 2/1/2 against R 1/0/1) -- and the
# starting vector is what the emitted call hands the optimizer, so mirroring
# R's defaults here started a partially-specified job somewhere PROC HAZARD
# would not have started it. On a multimodal likelihood a different start is a
# different answer, not a cosmetic difference (AGENTS.md, "n_starts explores a
# neighbourhood").
#
# TAU is the one entry that is NOT SAS's. src/hazard/stmtprc.c:37 sets it to
# ZERO, but that zero never reaches SETG3(): readobs.c:153-154 replaces an
# unspecified TAU with 0.75*Tmax on an active late phase, and readobs() runs
# first (hazard.c:276 against :292). The 1 below is SETG3_ignore_tau()'s pinned
# value (setg3.c:378), right exactly when that branch fires; every other
# defaulted-TAU late phase gets an untranslated row instead, naming whichever
# data-dependent rule applies. See the TAU guards in .hzr_parse_parms().
#
# There is no matching mu default: a phase is built only when its MU was
# specified and positive, so every built phase has a real MU to log.
.hzr_parms_early_sas_default <- c(t_half = 1, nu = 2, m = 1)
.hzr_parms_late_sas_default  <- c(tau = 1, gamma = 1, alpha = 1, eta = 2)

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
#' names, starting values and per-variable options.
#'
#' The real grammar is a comma-separated list of `phasevaropt : phasevar
#' phaseval phaseoptspec` (`hazard_y.y`): `VAR`, an optional `= startvalue`,
#' and an optional `/ options` that belongs to THAT variable alone. The lexer
#' leaves option state at the next comma (`hazard_l.l`, `<PHOP>\\,`), so
#' `AGE, MAL/I, OPMOS` is three covariates. Cutting the list at the first
#' `/` instead dropped every later covariate from the model (#342).
#'
#' Options (lexer aliases in brackets): `EXCLUDE` (`E`), `INCLUDE` (`I`),
#' `START` (`S`), `MOVE=` (`M`), `ORDER=` (`O`). Without a SELECTION
#' statement `setstat.c` puts a bare, `START` or `INCLUDE` variable in the
#' model and leaves an `EXCLUDE` one out, so an excluded variable is returned
#' in `excluded`, not in `names`. `przconc.c` tests EXCLUDE before INCLUDE
#' before START, and so does this. `MOVE=`, `ORDER=` and anything unrecognised
#' are recorded per variable, never dropped silently.
#'
#' `x` may be a single raw operand string (from the job parser) or a
#' character vector of pieces (the `.hzr_parse_parms()` `covars=` back-compat
#' interface); each element is split on `,`.
#' @return `list(names, values, flags, excluded, untranslated_construct,
#'   untranslated_reason)`. `flags` is parallel to `names`: `""`, `"I"` or
#'   `"S"`.
#' @noRd
.hzr_parse_phase_covars <- function(x) {
  names_out <- character(0)
  values_out <- numeric(0)
  flags_out <- character(0)
  bad_construct <- character(0)
  bad_reason <- character(0)
  bad <- function(construct, reason) {
    bad_construct <<- c(bad_construct, construct)
    bad_reason <<- c(bad_reason, reason)
  }
  # Text PROC HAZARD refuses to run. A job carrying any of it produces no
  # estimates, so the caller emits a stop() rather than a fit (#340).
  rejected <- character(0)
  reject <- function(construct, reason) {
    bad(construct, reason)
    rejected <<- c(rejected, paste0(construct, ": ", reason))
  }
  # Every form below sets yysynerr, and initprz.c:75-77 then exits "SYNTAX"
  # before any data are read. Each names its OWN source: the lexer and the
  # grammar reject for different reasons.
  syntax_error <- function(why) {
    paste0("PROC HAZARD stops the job with a syntax error: ", why,
           ", and initprz.c:75-77 exits \"SYNTAX\"")
  }
  # The word rule [.\-_A-Z0-9]+ is active in every lexer state and matches
  # a run such as EI whole, beating E on flex's longest match.
  bad_text <- syntax_error(paste(
    "the option state has no rule for this text, so hazard_l.l:176 reads",
    "it as unexpected and sets yysynerr"))
  # The option state (PHOP) has no "/" rule; only PHVR does (hazard_l.l:151).
  second_slash <- syntax_error(paste(
    "a second \"/\" has no rule in the option state and falls to the",
    "catch-all at hazard_l.l:178, which sets yysynerr"))
  # hazard_y.y:210 phasevaropt starts with the variable, and 220-225 need at
  # least one option after "/"; yyerror.c:19 sets yysynerr on a parse error.
  no_variable <- syntax_error(paste(
    "an option needs a variable before its \"/\" (hazard_y.y:210), so",
    "the parser fails (yyerror.c:19)"))
  no_option <- syntax_error(paste(
    "a \"/\" needs at least one option after it (hazard_y.y:220-225), so",
    "the parser fails (yyerror.c:19)"))
  # hazard_y.y:228-232: MOVE and ORDER take `= NUMBER`; E, I and S take none.
  bad_form <- syntax_error(paste(
    "MOVE= and ORDER= need a number, and E, I and S take no value",
    "(hazard_y.y:228-232), so the parser fails (yyerror.c:19)"))
  # After a phase variable, "=" needs a NUMBER (hazard_y.y:216-218).
  no_value <- syntax_error(paste(
    "a phase variable's \"=\" needs a number (hazard_y.y:216-218), so the",
    "parser fails (yyerror.c:19)"))
  word_value <- syntax_error(paste(
    "the text after \"=\" is not a number to the lexer, so hazard_l.l:176",
    "reads it as unexpected and sets yysynerr"))
  char_value <- syntax_error(paste(
    "a character after \"=\" has no lexer rule and falls to the catch-all",
    "at hazard_l.l:178, which sets yysynerr"))
  # The lexer's NUMBER (hazard_l.l:34-38). as.numeric() also reads Inf, NaN,
  # 1e5, 5. and 0x1A, none of which PROC HAZARD lexes as a number.
  is_number <- .hzr_sas_lexer_number
  # Which rule reads a value that is not a NUMBER. Per input, not per
  # refusal: a character outside the word rule's set falls to the catch-all;
  # a whole name lexes as NAME after a phase variable (hazard_l.l:174-175)
  # and as an option keyword after "/" (hazard_l.l:154-163), and either way
  # the parser finds it where NUMBER belongs; anything else is a word.
  value_error <- function(s, after_slash) {
    s <- toupper(s)
    if (!nzchar(s)) return(if (after_slash) bad_form else no_value)
    if (grepl("[^-._A-Z0-9]", s)) return(char_value)
    if (after_slash && s %in% c("E", "EXCLUDE", "I", "INCLUDE", "S", "START",
                                "M", "MOVE", "O", "ORDER")) return(bad_form)
    if (!after_slash && grepl("^[_A-Z][_A-Z0-9]*$", s)) return(no_value)
    word_value
  }

  for (piece in x) {
    for (p in strsplit(piece, ",", fixed = TRUE)[[1L]]) {
      p <- trimws(p)
      if (!nzchar(p)) next
      opts <- character(0)
      slash <- .idx(p, "/")
      # hazard_l.l has no "/" rule in option state, and a "/" with no name
      # before it is a syntax error, so PROC HAZARD rejects both. Record the
      # item rather than read it as a covariate or an option it is not.
      if (slash == 1L) {
        reject(p, no_variable)
        next
      }
      if (slash > 0L && .idx(substr(p, slash + 1L, nchar(p)), "/") > 0L) {
        reject(p, second_slash)
        next
      }
      if (slash > 0L) {
        opt_txt <- gsub("\\s*=\\s*", "=", substr(p, slash + 1L, nchar(p)))
        opts <- strsplit(trimws(opt_txt), "[[:space:]/]+")[[1L]]
        opts <- toupper(opts[nzchar(opts)])
        # phaseopts needs at least one option (hazard_y.y), so a bare "/" is
        # a syntax error PROC HAZARD rejects, not an option-free covariate.
        if (!length(opts)) {
          reject(p, no_option)
          next
        }
        p <- trimws(substr(p, 1L, slash - 1L))
      }
      eq <- .idx(p, "=")
      var <- if (eq == 0L) p else trimws(substr(p, 1L, eq - 1L))
      val <- NA_real_
      if (eq > 0L) {
        val_chr <- trimws(substr(p, eq + 1L, nchar(p)))
        if (!is_number(val_chr)) {
          reject(p, value_error(val_chr, after_slash = FALSE))
          next
        }
        val <- as.numeric(val_chr)
      }

      flag <- ""
      order_given <- FALSE
      for (o in opts) {
        key <- sub("=.*$", "", o)
        has_val <- grepl("=", o, fixed = TRUE)
        val_ok <- has_val && is_number(sub("^[^=]*=", "", o))
        # Each token must be one that the option state lexes whole: E/I/S (or
        # their long forms) alone, M/O (or MOVE/ORDER) with `= number`
        # (hazard_y.y:228-232). Anything else is a syntax error.
        is_flag <- key %in% c("E", "EXCLUDE", "I", "INCLUDE", "S", "START")
        is_valued <- key %in% c("M", "MOVE", "O", "ORDER")
        if (!is_flag && !is_valued) {
          reject(paste0(var, "/", o), bad_text)
          next
        }
        if ((is_flag && has_val) || (is_valued && !val_ok)) {
          # Which rule rejects the value depends on the value; a flag with
          # "=" or a MOVE/ORDER with none fails in the grammar.
          reject(paste0(var, "/", o),
                 if (is_valued && has_val) {
                   value_error(sub("^[^=]*=", "", o), after_slash = TRUE)
                 } else {
                   bad_form
                 })
          next
        }
        if (key %in% c("O", "ORDER")) order_given <- TRUE
        if (key %in% c("E", "EXCLUDE")) {
          flag <- "E"
        } else if (key %in% c("I", "INCLUDE")) {
          if (flag != "E") flag <- "I"
        } else if (key %in% c("S", "START")) {
          if (!flag %in% c("E", "I")) flag <- "S"
        } else if (key %in% c("M", "MOVE", "O", "ORDER")) {
          long <- if (key %in% c("M", "MOVE")) "MOVE" else "ORDER"
          bad(paste0(var, "/", long, sub("^[^=]*", "", o)),
              sprintf(paste("per-variable %s= option on phase-statement",
                            "covariate %s has no hazard() equivalent"),
                      long, var))
        }
      }
      # ORDER= with /E, /I or /S: przconc.c:45-53 logs "ORDER= and %s are
      # mutually exclusive" and sets semerr, and hazard.c:249-251 exits
      # ("SEMANTIC") before any data are read.
      if (order_given && nzchar(flag)) {
        reject(paste0(var, "/", flag, " ORDER="), paste(
          "PROC HAZARD refuses the job: ORDER= and /E, /I or /S are mutually",
          "exclusive (przconc.c:45-53 sets semerr; hazard.c:249-251 exits",
          "\"SEMANTIC\" before reading any data)"))
      }
      names_out <- c(names_out, var)
      values_out <- c(values_out, val)
      flags_out <- c(flags_out, flag)
    }
  }

  # A repeated covariate is one parameter: setconc.c maps every occurrence to
  # the same slot, and setstat.c runs for each in turn, so the LAST occurrence
  # sets its start value (0 when omitted) and its options. It keeps its first
  # position. Emitting it twice put two entries in theta for one column.
  last <- !duplicated(names_out, fromLast = TRUE)
  first_pos <- match(names_out[last], names_out)
  ord <- order(first_pos)
  names_out <- names_out[last][ord]
  values_out <- values_out[last][ord]
  flags_out <- flags_out[last][ord]
  excluded <- names_out[flags_out == "E"]
  keep <- flags_out != "E"
  names_out <- names_out[keep]
  values_out <- values_out[keep]
  flags_out <- flags_out[keep]

  list(names = names_out, values = values_out, flags = flags_out,
       excluded = excluded, untranslated_construct = bad_construct,
       untranslated_reason = bad_reason, rejected = rejected)
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
#' cannot disagree about what an unspecified operand started at; filling only
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
.hzr_parms_phase_call <- function(type, shape, covars, fixed,
                                  constraint = "none") {
  call_args <- c(list(quote(hzr_phase), type), shape)
  if (length(covars)) {
    call_args <- c(call_args,
                    list(formula = str2lang(paste("~", paste(covars, collapse = " + ")))))
  }
  fixed_call <- .hzr_parms_fixed_call(fixed)
  if (!is.null(fixed_call)) call_args <- c(call_args, list(fixed = fixed_call))
  if (constraint != "none") {
    call_args <- c(call_args, list(constraint = constraint))
  }
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
# are driven by FIXGE2/FIXGAE2, which .hzr_parse_parms() handles itself: it
# maps them onto hzr_phase(constraint = ) for WEIBULL, mirrors (or refuses)
# SETG3_ignore_tau() when that branch takes them, with or without WEIBULL, and
# records every other case -- and does not call this trace for the
# SETG3_ignore_tau() phases. So this trace takes both as FALSE rather than
# guessing. Every branch below is therefore the !g_two && !ga_two column of
# the C's own tables.

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
#   SETG3940, SETG3990, SETG31000
#             are reachable only through g_two/ga_two. The constraint block
#             below raises them itself now that FIXGE2/FIXGAE2 are mapped
#             (#329, #359), so they are absent HERE but not unreachable.
#   SETG31010 is the non-WEIBULL twin of SETG3990 (setg3.c:884-889), on the
#             SETG3_verify_ge_2() path this translator does not trace: a
#             non-WEIBULL job carrying either flag is recorded, not refused.
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
# makes them dead -- change those and these wake up. Nine codes can fire from
# THIS trace, which an exhaustive search in the tests pins rather than
# asserts; the constraint block below raises three more (SETG3940, SETG3990,
# SETG31000), so twelve are reachable through .hzr_parse_parms().

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

#' Would `hzr_phase()` build a `"g3"` phase from these shapes?
#'
#' The refusal message used to assert that it would, from a hand-maintained
#' idea of which codes were "shape" refusals. That drifted: it was false for
#' five of the seven SETG3 classes, not the three the text claimed
#' (`SETG3910`, `SETG3920`, `SETG3930`, and also `SETG3960` and `SETG3970`),
#' because SAS refuses several of them precisely BECAUSE a shape is out of
#' range, and the same value is out of range for `hzr_phase()`.
#'
#' So the sentence is now derived by CONSTRUCTING the phase. It cannot drift
#' again: if `hzr_phase()` changes what it accepts, this answer changes with
#' it (#433 review).
#' @return `TRUE` when the phase builds, `FALSE` when it refuses.
#' @noRd
.hzr_phase_builds <- function(tau, gamma, alpha, eta) {
  tryCatch({
    hzr_phase("g3", tau = tau, gamma = gamma, alpha = alpha, eta = eta)
    TRUE
  }, error = function(e) FALSE)
}

#' The four SETG3 entry refusals, in the C's own order.
#'
#' `setg3.c:269-284` (in `src/model/`, at pin `dad7978`) checks TAU, then
#' GAMMA, then ALPHA, then ETA, and **each one returns immediately**, before
#' `SETG3_ignore_tau()` at `:309-323` and before the WEIBULL branch. So a job
#' that would also trip a later rule is refused by the FIRST of these that
#' matches, and a message naming the later code names a refusal PROC HAZARD
#' never reaches.
#'
#' Kept as one function because two callers need the same order: the trace in
#' `.hzr_setg3_notes()`, and the constraint block, which records its own
#' `SETG3980` and must not do so ahead of an entry refusal (#433 review).
#'
#' An absent TAU (`NA`) is `0.75*Tmax` by the time SETG3 sees it, so only an
#' explicitly non-positive `TAU=` can refuse here.
#' @return The code, or `NULL` when none of the four applies.
#' @noRd
.hzr_setg3_entry_refusal <- function(tau_raw, gamma, alpha, eta, fixed) {
  fx <- function(p) p %in% fixed
  if (isTRUE(tau_raw <= 0) && fx("tau")) return("(SETG3900)")
  if (isTRUE(gamma <= 0) && fx("gamma")) return("(SETG3910)")
  if (isTRUE(alpha < 0) && fx("alpha")) return("(SETG3920)")
  if (isTRUE(eta <= 0) && fx("eta")) return("(SETG3930)")
  NULL
}

.hzr_setg3_notes <- function(tau_raw, gamma, alpha, eta, fixed, weibull) {
  fx <- function(p) p %in% fixed
  # `entry` marks the four checks at setg3.c:269-284, the only ones raised
  # BEFORE the TAU rules at :309-323. The other twelve codes fire after those
  # rules have already run, which decides whether a TAU row still applies.
  refuse <- function(code, entry = FALSE) {
    list(refusal = code, entry = entry, shape = NULL)
  }

  entry_code <- .hzr_setg3_entry_refusal(tau_raw, gamma, alpha, eta, fixed)
  if (!is.null(entry_code)) return(refuse(entry_code, TRUE))

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
#' @param selection `FALSE` for a job with no `SELECTION` statement,
#'   `"screen"` for a forward or two-way screen (bare variables are
#'   candidates, withheld from the phase formulas), or `"backward"` (bare
#'   variables start in the model, as `setstat.c` puts them there when
#'   `H->sw` is 0).
#' @param covars Optional named list of phase covariates, e.g.
#'   `list(early = c("X1", "X2"), constant = , late = )`, from the operands of
#'   the `EARLY` / `CONSTANT` / `LATE` statements.
#' @return `list(phases = <call>, theta = <call>, has_phases = <logical>,
#'   refused = <logical>, untranslated = <data.frame>)`. `has_phases` is
#'   `TRUE` only when at least one phase was actually built from the operands
#'   (i.e. `phases` is not the empty `list()` call); callers use it to decide
#'   whether the job qualifies as multiphase at all. `refused` is `TRUE` only
#'   when no phase is active AND every operand was understood, meaning
#'   `PROC HAZARD` would raise `ERROR 1001` and run nothing
#'   (`src/hazard/modterm.c`); the caller emits a `stop()` in place of the fit.
#'   It is deliberately narrower than `!has_phases`: a statement this parser
#'   could not read is recorded but not refused, because the refusal is a
#'   claim about the reference and not about this parser.
#' @noRd
.hzr_parse_parms <- function(operands, covars = list(), selection = FALSE) {
  mu <- list()
  early <- list()
  late <- list()
  fixed_early <- character(0)
  fixed_late <- character(0)
  saw_weibull <- FALSE
  saw_ge2 <- FALSE
  saw_gae2 <- FALSE
  saw_mnu1 <- FALSE
  bad_construct <- character(0)
  bad_reason <- character(0)
  # Set when an operand could not be read at all -- an unresolved keyword, a
  # non-numeric value, or a keyword with no phase target. It gates the "no
  # phase selected" refusal below, which is a claim about what PROC HAZARD
  # would do and is only sound when this parser actually understood the whole
  # statement. Distinct from `bad_construct`, which also collects the later
  # SEMANTIC guards: `MUE=0 THALF=1` is fully readable and genuinely selects
  # no phase, so it must still refuse even though it flags THALF=1.
  unreadable <- FALSE
  # PARMS text PROC HAZARD rejects with a syntax error (initprz.c:75-77): it
  # joins the phase statements' `rejected`, and .hzr_parse_job() emits a
  # stop() rather than a fit (John's U1 decision, 2026-09-19). A macro
  # reference is not here: SAS expands it first, so it is not known to fail.
  parms_rejected <- character(0)
  # A job PROC HAZARD runs on a model this translation would not emit (U1):
  # .hzr_parse_job() stops on any entry here, naming each.
  not_mirrored <- character(0)
  delta_seen <- NULL

  flag_bad <- function(construct, reason) {
    bad_construct <<- c(bad_construct, construct)
    bad_reason <<- c(bad_reason, reason)
  }
  flag_syntax <- function(construct, reason) {
    flag_bad(construct, reason)
    parms_rejected <<- c(parms_rejected, paste0("PARMS ", construct, ": ", reason))
  }
  flag_unresolved <- function(op) {
    if (.hzr_sas_is_macro(op)) flag_bad(op, .hzr_parms_unresolved_why(op))
    else flag_syntax(op, .hzr_parms_unresolved_why(op))
  }

  # A SETG3 refusal is a job PROC HAZARD stops in shape(), before hzrg() fits
  # anything, so it has to leave this parser as more than prose:
  # .hzr_parse_job() turns `refusal_reason` into the stop() chunk. `refused`
  # is not widened for it -- that field is documented as modterm.c's ERROR
  # 1001 ("no phase selected"), and overloading it would lose that meaning.
  # The row is still recorded, so the document lists what was wrong.
  refusal_reason <- NA_character_
  flag_refusal <- function(construct, reason) {
    flag_bad(construct, reason)
    # The stop() chunk carries the reason and nothing else, and it tells the
    # reader to correct the operands named in it -- so the construct has to
    # travel with the reason or the message names nothing.
    if (is.na(refusal_reason)) {
      refusal_reason <<- paste0(reason, " (PARMS ", construct, ")")
    }
  }

  operands <- .hzr_sas_join_spaced(.hzr_sas_join_macro_calls(operands))
  spaced_piece <- .hzr_parms_spaced_pieces(operands)
  for (i in seq_along(operands)) {
    op <- operands[[i]]
    if (spaced_piece[i] > 0L) {
      unreadable <- TRUE
      if (spaced_piece[i] == 2L) {
        flag_syntax(op, .hzr_parms_rejected_piece_reason)
      } else {
        flag_bad(op, if (spaced_piece[i] == 1L) .hzr_parms_unresolved_piece_reason
                 else .hzr_parms_macro_piece_reason)
      }
      next
    }
    eq <- .idx(op, "=")

    if (eq > 0L) {
      key <- substr(op, 1L, eq - 1L)
      raw <- substr(op, eq + 1L, nchar(op))
      val <- suppressWarnings(as.numeric(raw))
      token <- .hzr_sas_token(key, "HAZARD", "PARM")
      if (is.na(token)) {
        unreadable <- TRUE
        flag_unresolved(op)
      } else if (is.na(val)) {
        unreadable <- TRUE
        if (.hzr_sas_is_macro(raw)) {
          # A macro value SAS expands first: not known to fail.
          flag_bad(op, sprintf(paste0(
            "PARMS value for %s is a SAS macro reference, which SAS resolves ",
            "before PROC HAZARD reads the statement, so this translation ",
            "cannot tell what it becomes"), key))
        } else if (grepl("?", raw, fixed = TRUE)) {
          # A template's placeholder: the lexer has no rule for `?` and its
          # catch-all sets yysynerr (hazard_l.l:178). Filling it from SAS's
          # default and fitting would answer a job that does not run.
          flag_syntax(op, paste0(
            "PARMS value ", raw, " for ", key, " is a template placeholder, ",
            "which PROC HAZARD's lexer rejects (hazard_l.l:178), so the job ",
            "does not run until it is filled in; fill it in and translate ",
            "the job again"))
        } else {
          # A word after `=` is not a NUMBER to the lexer (hazard_l.l:176):
          # PROC HAZARD stops with a syntax error, as for `NU = ABC`.
          flag_syntax(op, paste0(
            "PARMS value ", raw, " for ", key, " is not a number PROC ",
            "HAZARD's lexer reads (hazard_l.l:34-38), so PROC HAZARD rejects ",
            "this job with a syntax error and it does not run"))
        }
      } else if (!.hzr_sas_lexer_number(raw)) {
        unreadable <- TRUE
        flag_syntax(op, paste0(
          "PARMS value ", raw, " for ", key, " is not a number PROC HAZARD's ",
          "lexer reads (hazard_l.l:34-38), so PROC HAZARD rejects this job ",
          "with a syntax error and it does not run"
        ))
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
        # Decided after the loop: DELTA is read only by SETG1()
        # (setg1.c:306), which runs only for an active early phase
        # (shape.c:19-21), so a late-only job ignores it.
        # hazard_y.y:138 is last-wins, so a later DELTA=0 clears an
        # earlier non-zero one, as it does in PROC HAZARD.
        delta_seen <- if (identical(val, 0)) NULL else list(op = op, val = val)
      } else if (token %in% c("WEIBULL", "FIXGE2", "FIXGAE2", "FIXMNU1",
                              "FIXDELTA", names(.hzr_parms_fix_map))) {
        # A flag keyword given a value: the grammar has it as a bare token
        # (hazard_y.y:148-160), so `FIXNU=1` is a syntax error and the job
        # does not run (U1).
        unreadable <- TRUE
        flag_syntax(op, paste0(
          key, " takes no value in PROC HAZARD (hazard_y.y:148-160), so PROC ",
          "HAZARD rejects this job with a syntax error and it does not run"))
      } else {
        unreadable <- TRUE
        flag_bad(op, "PARMS keyword has no phase target")
      }
      next
    }

    token <- .hzr_sas_token(op, "HAZARD", "PARM")
    if (is.na(token)) {
      unreadable <- TRUE
      flag_unresolved(op)
    } else if (token %in% c(.hzr_parms_mu_order, names(.hzr_parms_early_arg),
                            names(.hzr_parms_late_arg), "DELTA")) {
      # A value keyword with no `= NUMBER` after it, and not the first piece
      # of a spaced operand (those are marked above).
      unreadable <- TRUE
      flag_syntax(op, paste0(
        op, " needs a value (", op, "=NUMBER, hazard_y.y:137-147), so PROC ",
        "HAZARD rejects this job with a syntax error and it does not run"
      ))
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
      # are separate PARMS keywords, FIXGE2 and FIXGAE2, handled below.
      saw_weibull <- TRUE
    } else if (token == "FIXGE2") {
      saw_ge2 <- TRUE
    } else if (token == "FIXGAE2") {
      saw_gae2 <- TRUE
    } else if (token == "FIXMNU1") {
      # Recorded below, once whether an early phase is active is known.
      saw_mnu1 <- TRUE
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
  # `late` itself is left untouched, so `tau_absent` still tells an unwritten
  # TAU from a written TAU=0.
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

  # A phase is active iff its MU was specified and *positive*. That is the
  # reference rule rather than an inference from it: parmprc.c:13,18,19
  # registers MUE/MUC/MUL through setparmno(), which sets C->phase[n] = 1 only
  # under `stmtfld(parmno) > ZERO` (setparmno.c:11-14), while the seven shape
  # operands go through setprmf() (parmprc.c:14-17,20-23), which never touches
  # C->phase[] at all. Absence and non-positivity are therefore the same
  # outcome, which is why this tests the value and not merely the name --
  # a presence test lets MUE=0 build a phase PROC HAZARD would not.
  mu_active <- function(key) !is.null(mu[[key]]) && isTRUE(mu[[key]] > 0)
  has_early <- mu_active("MUE")
  has_muc <- mu_active("MUC")
  has_late <- mu_active("MUL")
  # Whether a phase is BUILT, as distinct from whether its MU is active. An
  # active MU with no shape operand builds on PROC HAZARD's own shape defaults
  # (#345), but only when the whole statement was read: if any operand could
  # not be read, a shape operand may have been written in a form this parser
  # split apart (`THALF = 0.3`, which SAS lexes), and building on defaults
  # would fit a different model. Such a phase is not built, and is recorded
  # below. Every build, scope and trace gate reads these, not has_*.
  build_early <- has_early && (length(early) > 0L || !unreadable)
  build_late <- has_late && (length(late) > 0L || !unreadable)
  if (!is.null(delta_seen) && has_early) {
    # sprintf("%g"), not format(): this string is DATA, not just a message --
    # it lands in the untranslated frame and is grepped by callers and tests.
    # format() honours getOption("OutDec"), so a session with OutDec = ","
    # would write "DELTA = 0,5" and break both.
    not_mirrored <- c(not_mirrored, paste0(
      "DELTA = ", sprintf("%g", delta_seen$val), ", which R does not ",
      "implement (it assumes delta = 0); remove DELTA or fit the model by hand"))
    flag_bad(delta_seen$op, paste0(
      "DELTA = ", sprintf("%g", delta_seen$val), " is not implemented -- R ",
      "assumes delta = 0, so the emitted call fits a DIFFERENT model than this ",
      "job (rho, the time argument and the density Jacobian all differ)"))
  }

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

  # FIXGE2 (GAMMA*ETA = 2) and FIXGAE2 (GAMMA*ETA/ALPHA = 2), as
  # SETG3_weibull() applies them (setg3.c:444-481, and SETG3_alpha_fixup() at
  # :815-834), with hzd_late_t2p.c deciding which parameter is derived at each
  # step. Only that branch is traced. Outside WEIBULL SETG3 dispatches to
  # other SETG3_* functions with their own rewrites and refusals; that path is
  # not modelled, and such a job stops saying so (U1). Together, or with ALPHA
  # fixed at 1, they take
  # SETG3_ignore_tau() instead, which is mirrored below.
  late_constraint <- "none"
  ignore_tau_handled <- FALSE
  constraint_flags <- c("FIXGE2", "FIXGAE2")[c(saw_ge2, saw_gae2)]
  if (length(constraint_flags) && !build_late) {
    for (flag in constraint_flags) {
      flag_bad(flag, if (has_late) {
        "the late phase it constrains is not built (see the MUL row)"
      } else {
        "PARMS token has no phase target"
      })
    }
  } else if (length(constraint_flags) &&
               (length(constraint_flags) == 2L || ignore_tau)) {
    # setg3.c:312-314 takes SETG3_ignore_tau() for `g_two && ga_two`, or for
    # ALPHA fixed at 1, and before the WEIBULL dispatch, so with or without
    # WEIBULL. With either flag it then (setg3.c:377-400):
    #   - fixes TAU at 1;
    #   - fixes ALPHA at 1, refusing with SETG3940 if ALPHA is fixed at
    #     anything else;
    #   - fixes GAMMA and ETA at 2 and 1, or 1 and 2 when the job wrote ETA = 2.
    # All four are exact, not a numerical branch: at alpha = 1 the likelihood
    # sees only (t/tau)^(gamma*eta), and gamma*eta = 2 either way. So the
    # values and the fixes are mirrored (see the TAU pin above for why a
    # half-mirror leaves an aliased ridge).
    fx_user <- function(param) param %in% fixed_late_user
    written <- late_full
    tau_raw <- if (tau_absent) NA_real_ else late[["tau"]]
    # The entry refusals (setg3.c:269-284) run before SETG3_ignore_tau(); the
    # SETG3 trace below records those from the unrewritten operands.
    entry_refused <- (isTRUE(tau_raw <= 0) && fx_user("tau")) ||
      (isTRUE(written[["gamma"]] <= 0) && fx_user("gamma")) ||
      (isTRUE(written[["alpha"]] < 0) && fx_user("alpha")) ||
      (isTRUE(written[["eta"]] <= 0) && fx_user("eta"))
    if (!entry_refused && fx_user("alpha") &&
          !isTRUE(written[["alpha"]] == 1)) {
      # PROC HAZARD stops here, so nothing the trace below would describe
      # is ever reached.
      ignore_tau_handled <- TRUE
      flag_refusal(paste(constraint_flags, collapse = " "), paste0(
        "PROC HAZARD refuses this job: SETG3 raises (SETG3940) -- ",
        "SETG3_ignore_tau() must set ALPHA to 1, but ALPHA is fixed at ",
        sprintf("%g", written[["alpha"]]), " (setg3.c:382-385)"))
    } else if (!entry_refused) {
      eta_two <- isTRUE(written[["eta"]] == 2)
      late_full[["tau"]] <- 1
      late_full[["alpha"]] <- 1
      late_full[["gamma"]] <- if (eta_two) 1 else 2
      late_full[["eta"]] <- if (eta_two) 2 else 1
      fixed_late <- unname(.hzr_parms_late_arg)
      ignore_tau_handled <- TRUE
      # A TAU the job wrote is already reported below when ALPHA is fixed at
      # 1; with both flags and ALPHA free nothing else says it.
      moved <- c(
        if (!ignore_tau && !tau_absent && !isTRUE(tau_raw == 1)) {
          sprintf("TAU=%g", tau_raw)
        },
        vapply(c("gamma", "alpha", "eta"), function(param) {
          if (isTRUE(written[[param]] == late_full[[param]])) "" else
            sprintf("%s=%g", toupper(param), written[[param]])
        }, character(1))
      )
      moved <- moved[nzchar(moved)]
      if (length(moved)) {
        flag_bad(paste(moved, collapse = " "), paste0(
          "SETG3_ignore_tau() runs this phase at TAU = 1, ALPHA = 1, GAMMA = ",
          sprintf("%g", late_full[["gamma"]]), ", ETA = ",
          sprintf("%g", late_full[["eta"]]), ", all fixed (setg3.c:377-400), ",
          "because ", if (length(constraint_flags) == 2L) {
            "FIXGE2 and FIXGAE2 are both set"
          } else {
            paste(constraint_flags, "is set with ALPHA fixed at 1")
          },
          "; the emitted phase mirrors that, so the value(s) written here are ",
          "used by neither PROC HAZARD nor the translation"))
      }
    }
  } else if (length(constraint_flags)) {
    not_traced <- if (!saw_weibull) {
      paste0("without WEIBULL, SETG3 applies it through SETG3_verify_ge_2() ",
             "and SETG3_alpha_gener(), which this translator does not trace")
    }
    fx <- function(param) param %in% fixed_late
    gamma_ <- late_full[["gamma"]]
    eta_ <- late_full[["eta"]]
    alpha_ <- late_full[["alpha"]]
    # The operands as this job wrote them (after the shape defaults are
    # filled), kept for the moved-shape record at the end of this block.
    written_late <- late_full
    # SETG3_weibull() refuses on the operands as written (setg3.c:430-440)
    # BEFORE either constraint moves one, so rewriting first would repair a
    # job PROC HAZARD does not run. .hzr_setg3_notes() reports these refusals
    # below, from the unrewritten values, with one exception it cannot see:
    # it takes ga_two as FALSE, so a fixed ALPHA = 0 reads as g3flag 2 + 2 =
    # 4 there, while under FIXGAE2 g3flag is 1 + 2 = 3 and SETG3980 fires.
    alpha_zero_gae2 <- saw_gae2 && isTRUE(alpha_ == 0) && fx("alpha")
    setg3_refuses <- !isTRUE(gamma_ > 0) || !isTRUE(eta_ > 0) ||
      !isTRUE(alpha_ >= 0) || (isTRUE(alpha_ == 0) && !fx("alpha")) ||
      alpha_zero_gae2
    # An entry refusal comes FIRST, because setg3.c:269-284 returns on it
    # before the WEIBULL branch this SETG3980 case lives in. flag_refusal()
    # keeps only the first reason recorded, so recording SETG3980 here
    # unconditionally named a refusal PROC HAZARD never reaches: with
    # GAMMA=0 FIXGAMMA ... ALPHA=0 FIXALPHA FIXGAE2 WEIBULL, SAS raises
    # SETG3910 and never evaluates the alpha rule (#433 review).
    entry_first <- .hzr_setg3_entry_refusal(
      late_full[["tau"]], gamma_, alpha_, eta_, fixed_late)
    entry_construct <- c("(SETG3900)" = "TAU", "(SETG3910)" = "GAMMA",
                         "(SETG3920)" = "ALPHA", "(SETG3930)" = "ETA")
    if (is.null(not_traced) && !is.null(entry_first)) {
      flag_refusal(unname(entry_construct[entry_first]), paste0(
        "PROC HAZARD refuses this job: SETG3 raises ", entry_first, " -- ",
        .hzr_setg3_refusal_reason(entry_first),
        ". It is checked before the WEIBULL rules (setg3.c:269-284), so this ",
        "is the refusal PROC HAZARD reaches first"))
    } else if (is.null(not_traced) && alpha_zero_gae2) {
      flag_refusal("FIXGAE2", paste0(
        "PROC HAZARD refuses this job: SETG3 raises (SETG3980) -- ",
        .hzr_setg3_refusal_reason("(SETG3980)"),
        "; under FIXGAE2 a fixed ALPHA = 0 does not select the exponential ",
        "case (setg3.c:333-335)"))
    }
    if (is.null(not_traced) && setg3_refuses) {
      # Recorded by the SETG3 trace or just above; nothing to translate.
      NULL
    } else if (!is.null(not_traced)) {
      # Without WEIBULL, SETG3 dispatches on the signs of ALPHA, GAMMA and ETA
      # (setg3.c:357-374) to SETG3_all_gt_0(), SETG3_alpha_le_0() and their
      # siblings, each with its own rewrites and refusals. This translation
      # does not model that path. Deriving it by hand went wrong in both
      # directions in two review passes on the U1 branch -- jobs SAS runs came
      # back refused, jobs SAS refuses were fitted -- so it claims neither: the
      # document stops and says it cannot tell (U1). SETG3's entry refusals
      # (setg3.c:269-284), which run before this path, still stop as refusals.
      not_mirrored <- c(not_mirrored, paste0(
        paste(constraint_flags, collapse = " and "), " without WEIBULL: this ",
        "translation does not model PROC HAZARD's constraint path without ",
        "WEIBULL (setg3.c:357-374 and the SETG3_* functions it dispatches to), ",
        "so it cannot tell whether PROC HAZARD refuses this job or which model ",
        "it fits. Add WEIBULL if the job means it, or fit the model by hand"))
      for (flag in constraint_flags) {
        flag_bad(flag, paste0(flag, " is not translated: ", not_traced))
      }
    } else if (saw_ge2) {
      # setg3.c:444-472: make GAMMA*ETA = 2 by moving the operand that is not
      # fixed; if either is fixed both become fixed, and if both are free ETA
      # is marked fixed and hzd_late_t2p.c:37-39 derives it as 2/GAMMA.
      # Exact, as the C is (`gte!=TWO`, setg3.c:450): a product a few ulps
      # from 2 is moved, or refused, there too.
      if (!isTRUE(gamma_ * eta_ == 2)) {
        if (fx("gamma") && fx("eta")) {
          flag_refusal("FIXGE2", paste0(
            "PROC HAZARD refuses this job: SETG3 raises (SETG3990) -- GAMMA ",
            "and ETA are both fixed and GAMMA*ETA = ", sprintf("%g", gamma_ * eta_),
            ", not 2, so neither can be adjusted"))
        } else if (fx("gamma")) {
          late_full[["eta"]] <- 2 / gamma_
        } else {
          late_full[["gamma"]] <- 2 / eta_
        }
      }
      if (fx("gamma") || fx("eta")) {
        fixed_late <- union(fixed_late, c("gamma", "eta"))
      } else {
        late_constraint <- "eta_gamma"
      }
    } else if (fx("alpha") &&
                 !isTRUE(gamma_ * eta_ / late_full[["alpha"]] == 2)) {
      # SETG3_alpha_fixup() (setg3.c:817-826) tests a fixed ALPHA against the
      # constraint before it asks whether GAMMA or ETA is free, so this
      # refusal holds whatever else is fixed.
      flag_refusal("FIXGAE2", paste0(
        "PROC HAZARD refuses this job: SETG3 raises (SETG31000) -- ALPHA is ",
        "fixed at ", sprintf("%g", late_full[["alpha"]]), " where FIXGAE2 ",
        "must move it to GAMMA*ETA/2 = ", sprintf("%g", gamma_ * eta_ / 2)))
    } else if (fx("gamma") && fx("eta")) {
      # SETG3_alpha_fixup(): with no free GAMMA or ETA to derive from, ALPHA is
      # only moved onto the constraint as a starting value and stays free
      # (hzr_parms_ge_1estim() is FALSE, so :832-833 never fixes it). A fixed
      # ALPHA reaching here already sits on the constraint.
      if (!fx("alpha")) late_full[["alpha"]] <- gamma_ * eta_ / 2
    } else {
      # ALPHA is marked fixed and hzd_late_t2p.c:90-94 recomputes it as
      # GAMMA*ETA/2 from theta at every step: derived, not held.
      late_full[["alpha"]] <- gamma_ * eta_ / 2
      fixed_late <- setdiff(fixed_late, "alpha")
      late_constraint <- "alpha_gamma_eta"
    }
    fixed_late <- intersect(unname(.hzr_parms_late_arg), fixed_late)

    # setg3.c:449-467 and :827 move a shape onto the constraint and call
    # hzr_parm_changed(), so PROC HAZARD tells its own reader. The emitted
    # phase is that model, but the job wrote something else, and the SETG3
    # notes trace below covers only its own rewrites -- not this block's,
    # added with the constraint mapping. Record them here so the emitted
    # document says what changed (#359).
    moved_by_flags <- vapply(c("gamma", "alpha", "eta"), function(param) {
      if (isTRUE(written_late[[param]] == late_full[[param]])) "" else
        sprintf("%s=%g -> %g", toupper(param), written_late[[param]],
                late_full[[param]])
    }, character(1))
    moved_by_flags <- moved_by_flags[nzchar(moved_by_flags)]
    # No `refusal_reason` guard here: a refused job never reaches a rewrite,
    # so the two cannot coexist (a mutant allowing both changes nothing).
    if (length(moved_by_flags)) {
      flag_bad(paste(moved_by_flags, collapse = " "), paste0(
        paste(constraint_flags, collapse = " and "),
        " moves the late shape onto the constraint before fitting ",
        "(setg3.c:449-467, :827), as PROC HAZARD does and reports through ",
        "hzr_parm_changed(): ", paste(moved_by_flags, collapse = ", "),
        ". The emitted phase is the model PROC HAZARD fits, not the operands ",
        "written here"))
    }
  }

  # EARLY/CONSTANT/LATE operand text: comma-separated VAR=VALUE pairs (or
  # bare VARs), each with its own optional "/ options". Non-numeric values
  # and options with no hazard() equivalent are recorded to untranslated,
  # never guessed at; see .hzr_parse_phase_covars(). VAR=VALUE starting values are now mapped into
  # theta (one entry per covariate, appended after that phase's shape block,
  # per .hzr_phase_theta_names()); a bare VAR with no value defaults to 0,
  # matching .hzr_phase_start().
  phase_covars <- list()
  # Phase-statement text PROC HAZARD refuses to run (#340).
  rejected <- character(0)
  # Every covariate a phase statement names (not /E, which is excluded and
  # guarded through listwise_only), before SELECTION withholds its
  # candidates from phase_covars. A row about a phase that is not built must
  # name all of them, or a candidate-only phase vanishes without a trace.
  phase_named <- list()
  phase_covar_vals <- list()
  phase_vars <- character(0)
  # Under SELECTION a bare variable starts OUT of the model and is a
  # candidate; /S starts in and may move; /I starts in and never moves
  # (setstat.c with H->sw == 1). Without SELECTION every non-/E variable is
  # simply in the model, which is the `selection = FALSE` path.
  sel_candidates <- list()
  sel_movable <- list()
  # Kept per phase: /I in a phase that is not built pins nothing, because
  # PROC HAZARD skips that phase's variables and their flags (setstat.c:9-12).
  sel_force_in <- list()
  for (ph in c("early", "constant", "late")) {
    raw <- covars[[ph]]
    if (is.null(raw)) {
      phase_covars[[ph]] <- character(0)
      phase_named[[ph]] <- character(0)
      phase_covar_vals[[ph]] <- numeric(0)
      next
    }
    parsed <- .hzr_parse_phase_covars(raw)
    # Candidates are withheld only from a FORWARD or two-way screen. Under
    # BACKWARD, stpwprc.c leaves H->sw at 0, so setstat.c puts a bare
    # variable IN the model and the screen drops from the full set.
    withhold <- identical(selection, "screen")
    keep <- if (withhold) parsed$flags %in% c("I", "S") else
      rep(TRUE, length(parsed$names))
    if (!isFALSE(selection)) {
      sel_candidates[[ph]] <- parsed$names[parsed$flags == ""]
      sel_movable[[ph]] <- parsed$names[parsed$flags %in% c("", "S")]
      sel_force_in[[ph]] <- parsed$names[parsed$flags == "I"]
    }
    phase_covars[[ph]] <- parsed$names[keep]
    phase_named[[ph]] <- parsed$names
    phase_covar_vals[[ph]] <- parsed$values[keep]
    phase_vars <- c(phase_vars, parsed$names, parsed$excluded)
    rejected <- c(rejected, parsed$rejected)
    for (i in seq_along(parsed$untranslated_construct)) {
      flag_bad(parsed$untranslated_construct[[i]], parsed$untranslated_reason[[i]])
    }
  }

  # Phases are built, and their theta blocks appended, in the same
  # early -> constant -> late order -- .hzr_phase_theta_names() assigns
  # labels by *position* (phases are auto-named "phase_1", "phase_2", ... by
  # .hzr_validate_phases()), so only this order matters, not the names. A
  # phase that is not built must therefore contribute no theta block either,
  # or every later block is read against the wrong labels.
  #
  # The MU gate is the whole gate: an active MU whose phase has no shape
  # operand is built on PROC HAZARD's own shape defaults, which
  # .hzr_parms_fill_shape() has already filled in (#345).
  # A shape that is not finite, as written (GAMMA=1e400 reads as Inf) or after
  # a rewrite above (2/ETA, GAMMA*ETA/2), cannot be built: hzr_phase() refuses
  # it. Say so here rather than let the translation read clean over a call
  # that stops.
  for (shape in list(if (build_early) early_full,
                     if (build_late) late_full)) {
    for (param in names(shape)) {
      value <- shape[[param]]
      if (is.numeric(value) && length(value) == 1L && !is.finite(value)) {
        # After SETG3's rewrites, not as written: an operand that is infinite
        # as written but replaced by a rewrite is not flagged, because PROC
        # HAZARD reads it the same way (hazard_l.l:53 scans with sscanf and
        # no range check) and applies the same rewrite, so the emitted model
        # is the one it fits.
        flag_bad(sprintf("%s=%g", toupper(param), value), paste0(
          toupper(param), " is not a finite number after SETG3's rewrites, ",
          "so the emitted hzr_phase() call cannot be built"))
      }
    }
  }

  phase_calls <- list()
  theta_blocks <- list()
  if (build_early) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "cdf", early_full, phase_covars$early, fixed_early
    )
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("early", mu[["MUE"]], early_full,
                             phase_covar_vals$early)
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
  if (build_late) {
    phase_calls[[length(phase_calls) + 1L]] <- .hzr_parms_phase_call(
      "g3", late_full, phase_covars$late, fixed_late, late_constraint
    )
    theta_blocks <- c(theta_blocks,
      .hzr_parms_theta_block("late", mu[["MUL"]], late_full,
                             phase_covar_vals$late)
    )
  }

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
  #
  # Both guards describe things SETG3() does, and shape.c:31 calls SETG3() only
  # under Common.phase[3] == 1, so both require an active late phase. Warning
  # about them for TAU/GAMMA with no MUL would describe code PROC HAZARD never
  # reaches, and a false positive on $untranslated is not harmless -- that
  # frame is how a caller decides whether a translation can be trusted.
  #
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
  # The trace assumes neither constraint flag, so it does not describe a phase
  # SETG3_ignore_tau() ran under one: that phase is fully determined above, or
  # refused there with SETG3940.
  if (build_late && !ignore_tau_handled) {
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
      # The setg3.c trace .hzr_setg3_notes() encodes assumes neither constraint
      # flag is set. With FIXGE2 or FIXGAE2 and no WEIBULL, SAS reaches SETG3
      # down a path the trace does not model, and jobs PROC HAZARD RUNS come
      # back refused here: `MUL=0.2 TAU=1 GAMMA=4 ETA=0.5 FIXGAMMA FIXETA
      # FIXGE2` and `MUL=0.2 TAU=1 GAMMA=4 ETA=1 ALPHA=2 FIXALPHA FIXGAE2` both
      # fit. Such a job is neither refused nor silently passed: the row says
      # what is true, which is that this parser cannot decide it.
      # The entry checks (setg3.c:269-284) run before any constraint or
      # WEIBULL logic, so they are refusals on every path (U1 review).
      not_traced <- length(constraint_flags) && !saw_weibull &&
        !isTRUE(setg3$entry)
      construct <- sprintf("GAMMA=%g ALPHA=%g ETA=%g%s", late_full[["gamma"]],
                           late_full[["alpha"]], late_full[["eta"]],
                           if (length(fixed_late_user)) {
                             paste0(" fixed:",
                                    paste(fixed_late_user, collapse = ","))
                           } else {
                             ""
                           })
      if (not_traced) {
        flag_bad(construct, paste0(
          "this shape reaches SETG3 with ",
          paste(constraint_flags, collapse = " and "),
          " set and no WEIBULL, which this translation does not trace, so ",
          "whether PROC HAZARD refuses the job (it would raise ",
          setg3$refusal, " on the traced path) is not decided here. The ",
          "emitted phase is the one the PARMS statement writes; check the ",
          "SAS log before relying on the fit"
        ))
      } else {
        builds <- .hzr_phase_builds(late_full[["tau"]], late_full[["gamma"]],
                                    late_full[["alpha"]], late_full[["eta"]])
        flag_refusal(construct, paste0(
          "PROC HAZARD refuses this job: SETG3 raises ",
          setg3$refusal, " -- ",
          .hzr_setg3_refusal_reason(setg3$refusal),
          if (builds) {
            paste0(". hzr_phase() accepts this shape, so the translated fit ",
                   "would converge on a job SAS never fits")
          } else {
            paste0(". hzr_phase() will not build this shape either, because ",
                   "the value SAS refuses is also outside the range it ",
                   "accepts, so the document stops at that check rather ",
                   "than fitting")
          }
        ))
      }
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
  if (build_late && !setg3_refused && ignore_tau &&
      !is.null(late[["tau"]]) && !isTRUE(late[["tau"]] == 1)) {
    flag_bad(
      paste0("TAU=", sprintf("%g", late[["tau"]])),
      paste0("SETG3_ignore_tau() discards this TAU and runs at tau = 1, ",
             "fixed (setg3.c:378-379), because ALPHA is fixed at 1; the ",
             "emitted phase mirrors that, so the value written here is used ",
             "by neither PROC HAZARD nor the translation")
    )
  } else if (build_late && !setg3_refused && tau_defaulted &&
             !ignore_tau && !ignore_tau_handled) {
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
    # An active MUL with no shape operand now builds its phase (#345), so this
    # row is the only record of its data-dependent TAU start.
    #
    # Two consequences, and the row names the one that applies. With TAU free,
    # the start differs, and because the multiphase likelihood is multimodal
    # that can change where the fit converges, not only how it gets there.
    # With FIXTAU (reachable only for an unwritten TAU: a written non-positive
    # one is refused as SETG3900), PROC HAZARD holds TAU at that
    # data-dependent value while the emitted phase holds it at 1, which is a
    # different model outright.
    tau_fixed <- "tau" %in% fixed_late
    if (tau_fixed && tau_absent) {
      not_mirrored <- c(not_mirrored, if (unreadable) {
        paste0("FIXTAU with a TAU this translation could not read (see the ",
               "rows above): PROC HAZARD fixes TAU at the value written, or ",
               "at 0.75*Tmax if none was (readobs.c:153-154), while the ",
               "emitted phase would pin it at 1")
      } else {
        paste0("FIXTAU with no TAU written, which PROC HAZARD fixes at ",
               "0.75*Tmax (readobs.c:153-154), a value that depends on the ",
               "data; write TAU= with the value to fix it at")
      })
    }
    flag_bad(
      if (tau_absent) "TAU (unspecified)" else
        paste0("TAU=", sprintf("%g", late[["tau"]])),
      paste0(if (tau_absent && unreadable) {
               paste0("TAU was not read here: this statement has operands ",
                      "the translator could not read (see their rows), so ",
                      "whether and how TAU was written cannot be told. If it ",
                      "was not written, ")
             },
             "PROC HAZARD ", if (tau_fixed) "fixes" else "starts", " TAU at ",
             if (tau_absent && unreadable) {
               "0.75*Tmax (readobs.c:153-154, "
             } else if (tau_absent) {
               "0.75*Tmax (readobs.c:153-154, applied to an unspecified TAU "
             } else {
               "2*Tmax/3 (setg3.c:317, applied to a non-positive TAU "
             },
             "before SETG3 runs), which depends on the data and cannot be ",
             "reproduced at parse time. ",
             if (tau_fixed) {
               paste0("The emitted phase fixes tau = 1 instead, so this is ",
                      "a different model, not only a different start: the ",
                      "estimates will differ from PROC HAZARD's")
             } else {
               paste0("The emitted phase starts at tau = 1. The multiphase ",
                      "likelihood is multimodal, so a different start can ",
                      "converge to a different optimum: the estimates, not ",
                      "only the path to them, may differ from PROC HAZARD's")
             })
    )
  }

  # Everything a MU gate discards has to leave a row behind. A phase that is
  # not built takes its shape operands, its FIX tokens *and* its covariates
  # with it, and silence on any of them hands back a model a phase short of the
  # SAS job with nothing to say so -- the failure this package keeps shipping
  # (see AGENTS.md). PROC HAZARD discards the same material, so the values are
  # right and only the silence would be wrong: stmtprc.c:101-122 zeroes the
  # shape operands and clears their status flags, and setstat.c:9-12 returns
  # early for a phase whose C->phase[] is 0, dropping its covariates.
  #
  # sprintf("%g"), not format(): format() honours getOption("OutDec"), so a
  # session with OutDec = "," would record "MUE=0,2" and break every grep --
  # the same trap the DELTA reason string above avoids.
  dropped <- function(shape, arg_map, fixed, covars) {
    sas <- names(arg_map)[match(names(shape), unname(arg_map))]
    ops <- sprintf("%s=%s", sas, vapply(shape, function(v) sprintf("%g", v), ""))
    fix <- names(.hzr_parms_fix_map)[match(fixed, unname(.hzr_parms_fix_map))]
    paste(c(ops, fix, covars), collapse = " ")
  }

  # (1) A MU that was specified but is not positive. stmtfld(parmno) > ZERO
  # fails, so the phase never activates; the keyword is not an error and is not
  # a phase either.
  for (key in names(mu)) {
    if (!mu_active(key)) {
      flag_bad(paste0(key, "=", sprintf("%g", mu[[key]])),
               paste0(key, " is not positive, so PROC HAZARD leaves that phase ",
                      "inactive (setparmno.c:11): the phase is not built"))
    }
  }

  # (2) Everything belonging to a phase whose MU never activated it.
  if (!has_early) {
    gone <- dropped(early, .hzr_parms_early_arg, fixed_early, phase_named$early)
    if (nzchar(gone)) {
      flag_bad(gone, paste0("early phase material with no active MUE: PROC ",
                            "HAZARD zeroes the shape operands (stmtprc.c:",
                            "101-112) and skips the covariates (setstat.c:9-12)"))
    }
  }
  if (!has_muc && length(phase_named$constant)) {
    flag_bad(paste(phase_named$constant, collapse = " "),
             paste0("constant phase covariates with no active MUC: PROC ",
                    "HAZARD skips them (setstat.c:9-12)"))
  }
  if (!has_late) {
    gone <- dropped(late, .hzr_parms_late_arg, fixed_late, phase_named$late)
    if (nzchar(gone)) {
      flag_bad(gone, paste0("late phase material with no active MUL: PROC ",
                            "HAZARD zeroes the shape operands (stmtprc.c:",
                            "113-122) and skips the covariates (setstat.c:9-12)"))
    }
  }

  # (3) No phase activated at all -- including the job carrying no PARMS
  # statement whatsoever, which arrives here with `operands` empty. PROC
  # HAZARD does not distinguish the two: stmtprc.c:87 zeroes all three phases
  # at init, and the only write that turns one back on is setparmno.c:14,
  # reached from parmprc.c:13,18,19 for MUE/MUC/MUL alone and only under
  # `stmtfld(parmno) > ZERO`. The seven shape operands go through setprmf()
  # (parmprc.c:14-17,20-23), which never touches C->phase[] at all. So "no
  # PARMS" and "PARMS naming no positive MU" are one state, not two cases.
  #
  # That state is refused, not merely defaulted, and modterm() is reached for
  # every job rather than only for a selected multiphase model: its single
  # call site is outmods.c:91, and outmods() sits in main's straight-line
  # sequence (hazard.c:296) with the no-phase test as one of the three
  # disjuncts at outmods.c:89 that fire the call. hazard.c:299-302 then routes
  # errorno 1001 to hzfxit("SEMANTIC"), which exits BEFORE results() -- so the
  # job prints no estimates and never fits. `refused` therefore propagates to
  # .hzr_parse_job(), which emits a stop() in place of the hazard() call: an
  # $untranslated row alone still renders a populated fit, and a fit standing
  # in for a job that produced no result is this package's signature defect.
  #
  # `unreadable` is the gate that keeps this from being a FALSE POSITIVE, and
  # it is load bearing now that the outcome is a hard stop. Testing only the
  # three has_* flags conflates "PARMS selected no phase" with "this parser
  # could not read PARMS". `PARMS MUE = 0.2 THALF = 1` is the case that bites:
  # .hzr_parse_hazard() splits operands on " ", so it yields "MUE", "=", "0.2"
  # and nothing parses -- while HAZARD's own lexer discards whitespace
  # unconditionally (hazard_l.l:32 `ws [ \t\n\r\)]+`, rule at :50), making
  # that a well-formed parmsopt (hazard_y.y:138) whose job RUNS with an active
  # early phase. Refusing it would stop a job PROC HAZARD accepts. The
  # per-operand rows still record what was not understood; only the refusal,
  # which is a claim about the reference, requires having understood it all.
  refused <- !has_early && !has_muc && !has_late && !unreadable
  if (refused) {
    flag_bad(
      # Keyed on `operands`, not on mu/early/late: the SECOND %HAZARD block
      # of dist/examples/hm.dthar.TGA.sas (line 110; its PARMS at 122) carries
      # literal `MUE=? THALF=? NU=?` for the reader to fill in from the
      # stepwise output, which leaves all three empty while the statement
      # plainly existed. The file's FIRST block (line 70, PARMS at 74) is a
      # valid `MUE=0.2 THALF=0.08 NU=1 M=1 FIXM MUC=0.001` activating phases 1
      # and 2 -- so this is a property of that one block, not of the file.
      # (It is unreadable, so it no longer reaches here -- but the wording
      # must not depend on that gate to be true.) A bare `PARMS;` cannot
      # occur: hazard_y.y:133-135 requires at least one parmsopt, so it is a
      # HAZARD syntax error, and character(0) here means the statement was
      # genuinely absent.
      if (length(operands)) {
        "PARMS with no positive MUE, MUC or MUL"
      } else {
        "no PARMS operands (no MUE, MUC or MUL)"
      },
      paste0("PROC HAZARD selects no phase here and refuses the job ",
             "(modterm.c:18-22 raises ERROR 1001, \"No phase selected\"; ",
             "hazard.c:299-302 then exits before results())")
    )
  }

  # FIXMNU1 is a real PROC HAZARD constraint (hazard_y.y:153; parmprc.c:29-38;
  # hzd_early_t2p.c:65-77 derives M = +/-1/NU or NU = +/-1/M at every step)
  # that this translation does not apply. On an active early phase that makes
  # the emitted phase a different model, and the row says the consequence
  # rather than a parse state. Mirroring it is separate work.
  if (saw_mnu1) {
    flag_bad("FIXMNU1", if (build_early) {
      paste0("FIXMNU1 ties M to NU in PROC HAZARD (|M*NU| = 1; ",
             "setg1.c:381-387, hzd_early_t2p.c:65-77), but that constraint is ",
             "not applied here: the emitted early phase does not tie them, a ",
             "different model from PROC HAZARD's")
    } else if (has_early) {
      "the early phase it constrains is not built (see the MUE row)"
    } else {
      "PARMS token has no phase target"
    })
  }

  # (4) An active MU whose phase carries no shape operand is built above, on
  # PROC HAZARD's shape defaults (stmtprc.c:30-37): early tHalf 1, nu 2, m 1,
  # all data-free (setg1.c:343-349 substitutes 1 only for a non-positive
  # tHalf), and late gamma 1, alpha 1, eta 2. The one data-dependent value is
  # the late TAU start, 0.75*Tmax (readobs.c:153-154, before SETG3), recorded by
  # the TAU row above exactly as for a written late phase with no TAU (#345).
  # This used to record the MU instead, when the parser's defaults were not
  # SAS's; they are now. An orphan whose statement was not fully read is the
  # exception (see build_early): it is recorded, not built.
  unread_why <- paste0(
    "with no ", "%s", " shape operand this translator could read: other PARMS ",
    "operands could not be read (see their rows), so whether a shape was ",
    "written cannot be told, and the phase is not built on PROC HAZARD's ",
    "defaults"
  )
  # An operand this parser could not read may be a shape or a FIX flag of a
  # phase it DID build, and the emitted phase then carries SAS's default
  # where the job wrote something else. It cannot be told apart from an
  # operand that changes nothing, so the document stops rather than fitting
  # a model that may not be PROC HAZARD's (U1).
  if (unreadable && (build_early || build_late || has_muc)) {
    not_mirrored <- c(not_mirrored, paste0(
      "operands of this PARMS statement were not read (see the rows above); ",
      "one of them may set a shape or a FIX flag of a phase this translation ",
      "did build, so it cannot tell whether the emitted phases carry the ",
      "values PROC HAZARD uses"))
  }

  # PROC HAZARD fits the phase whatever this parser could read, so a model
  # without it is short a phase: the document stops (U1).
  for (nm in c("MUE", "MUL")) {
    active <- if (nm == "MUE") has_early else has_late
    built <- if (nm == "MUE") build_early else build_late
    if (active && !built) {
      phase <- if (nm == "MUE") "early" else "late"
      flag_bad(paste0(nm, "=", sprintf("%g", mu[[nm]])),
               paste(nm, sprintf(unread_why, phase)))
      not_mirrored <- c(not_mirrored, paste0(
        "an active ", nm, " whose ", phase, " phase this translation could ",
        "not build, because operands of this PARMS statement were not read ",
        "(see the rows above); PROC HAZARD fits that phase, so the emitted ",
        "model would be short of it"))
    }
  }

  # getrisk.c collects every phase-statement variable, of every phase and
  # whatever its options, and readobs.c deletes a row where any is missing.
  # hazard() drops missing rows only for variables in a formula it fits, so
  # the rest -- /E variables, and covariates of a phase that is not built --
  # are returned for the caller to guard. "Built" is the same gate the phase
  # list uses above: the MU alone, since an orphan MU builds on PROC HAZARD's
  # shape defaults (#345). A stricter gate here keyed scope against a phase
  # list it did not match.
  modelled <- c(
    if (build_early) phase_covars$early,
    if (has_muc) phase_covars$constant,
    if (build_late) phase_covars$late
  )
  # Scope is keyed by the name the BASE FIT will carry. The emitted phases
  # list is unnamed, so hazard() auto-names them phase_1, phase_2, ... in
  # build order: keying on "early"/"constant" fails with "Unknown phase(s)
  # in scope". Only built phases have a key. sprintf(), not paste0():
  # paste0("phase_", integer(0)) is "phase_", so a job that builds no phase
  # crashed setNames() instead of reaching its "selects no phase" refusal.
  built <- c(
    if (build_early) "early",
    if (has_muc) "constant",
    if (build_late) "late"
  )
  selection_spec <- if (isFALSE(selection)) NULL else list(
    scope = stats::setNames(
      lapply(built, function(ph) sel_candidates[[ph]] %||% character(0)),
      sprintf("phase_%d", seq_along(built))),
    movable = stats::setNames(
      lapply(built, function(ph) sel_movable[[ph]] %||% character(0)),
      sprintf("phase_%d", seq_along(built))),
    in_model = stats::setNames(
      lapply(built, function(ph) phase_covars[[ph]] %||% character(0)),
      sprintf("phase_%d", seq_along(built))),
    force_in = unique(unlist(sel_force_in[built]))
  )
  list(
    phases = as.call(c(quote(list), phase_calls)),
    theta = as.call(c(quote(c), theta_blocks)),
    listwise_only = setdiff(unique(phase_vars), modelled),
    selection = selection_spec,
    has_phases = length(phase_calls) > 0L,
    refused = refused,
    rejected = c(parms_rejected, rejected),
    # The same vector split by PROVENANCE, because the two halves now get
    # different treatment and telling them apart by their text would be a
    # shape test on a message. `rejected_phase` is the phase statements',
    # which PROC HAZARD has always refused at parse and which the document
    # has always stopped on (#340). `rejected_parms` is the PARMS operands',
    # which used to fit with a row and now fit with a row AND a warning
    # (John's 2026-09-22 decision). `rejected` stays the union: it is read by
    # existing tests and by nothing that needs the distinction.
    rejected_phase = rejected,
    rejected_parms = parms_rejected,
    refusal_reason = refusal_reason,
    # FIXMNU1 on an active early phase: PROC HAZARD fits |M*NU| = 1, which
    # this translation does not mirror (#358), so the model it would emit is
    # a different one. .hzr_parse_job() stops on it (U1).
    not_mirrored = c(not_mirrored, if (saw_mnu1 && has_early) {
      paste0("FIXMNU1, which PROC HAZARD applies as |M*NU| = 1 on the early ",
             "phase and this translation does not mirror (#358); remove ",
             "FIXMNU1 to fit M and NU freely")
    }),
    untranslated = .hzr_untranslated_frame(
      line = rep(NA_integer_, length(bad_construct)),
      construct = bad_construct,
      reason = bad_reason
    )
  )
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
