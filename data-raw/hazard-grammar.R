# data-raw/hazard-grammar.R -- generate .hzr_sas_grammar from HAZARD's own lexer.
#
# WHY GENERATED, NOT HAND-WRITTEN
# -------------------------------
# The SAS HAZARD input language is defined by a lex/yacc grammar in the
# reference implementation (GPL-2, github.com/ehrlinger/hazard):
#
#   src/hazard/hazard_l.l    98 keywords, 11 start conditions
#   src/hazpred/hazpred_l.l  20 keywords
#
# Inferring the keyword table from a corpus of jobs captures only the spellings
# that corpus happens to use. The grammar has alias pairs that fall through to
# one parser token -- QUASI/QUASINEWTON, BW/BACKWARD, SLE/SLENTRY, NOH/NOHAZ,
# and the single-letter phase options E/I/M/O/S -- so a corpus-inferred table
# would silently fail on every spelling the sample did not contain.
#
# LICENCE
# -------
# `hazard` is GPL-2; TemporalHazard is MIT. This script reads the .l files from
# a local checkout and ships ONLY the extracted keyword table. No GPL source
# enters the package. data-raw/ is .Rbuildignore'd.
#
# USAGE
#   Rscript data-raw/hazard-grammar.R [path/to/hazard/checkout]

hazard_repo <- local({
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a)) a[[1L]] else path.expand("~/Documents/GitHub/hazard")
})

# ---------------------------------------------------------------------------
# Extract keyword rules from a lex file.
#
# A rule is  <CTX[,CTX]*>KEYWORD<tab>ACTION  where ACTION is either a bare "|"
# (alias: fall through to the next rule's action) or contains "return TOKEN;".
# Aliases accumulate until a rule supplies the return, which then applies to
# the whole run.
# ---------------------------------------------------------------------------
parse_lex_keywords <- function(path, proc) {
  if (!file.exists(path)) {
    stop("lex file not found: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)

  rule <- "^<([A-Z,]+)>([A-Z][A-Z0-9_]*)[[:space:]]+(.*)$"
  starts <- which(grepl(rule, lines))
  if (!length(starts)) {
    stop("no keyword rules parsed from ", path, call. = FALSE)
  }

  # A lex action may span lines: `<STMT>EARLY { BEGIN PHVR;` puts its
  # `return EARLY;` two lines further down. Reading only the rule's own line
  # drops such keywords silently -- it dropped EARLY, CONSTANT and LATE, the
  # phase covariate statements, from the first version of this table. Take
  # each rule's action to be everything up to the next rule.
  ends <- c(starts[-1L] - 1L, length(lines))

  out <- list()
  pending <- list()

  for (k in seq_along(starts)) {
    ln <- lines[starts[k]]
    p <- regmatches(ln, regexec(rule, ln))[[1L]]
    ctx <- p[[2L]]
    kw <- p[[3L]]
    action <- trimws(paste(
      c(p[[4L]], lines[seq_len(max(0L, ends[k] - starts[k])) + starts[k]]),
      collapse = " "
    ))

    if (grepl("^\\|", action)) {
      pending[[length(pending) + 1L]] <- list(ctx = ctx, kw = kw)
      next
    }

    ret <- regmatches(
      action,
      regexec("return[[:space:]]+([A-Za-z_][A-Za-z0-9_]*)", action)
    )[[1L]]
    if (length(ret) < 2L) {
      pending <- list()
      next
    }
    token <- ret[[2L]]

    run <- c(pending, list(list(ctx = ctx, kw = kw)))
    for (i in seq_along(run)) {
      out[[length(out) + 1L]] <- data.frame(
        proc     = proc,
        context  = run[[i]]$ctx,
        keyword  = run[[i]]$kw,
        token    = token,
        is_alias = i < length(run),
        stringsAsFactors = FALSE
      )
    }
    pending <- list()
  }

  res <- do.call(rbind, out)

  # FAIL LOUD on incomplete extraction. Asserting the table is internally
  # consistent says nothing about whether it covers the source; only this
  # comparison does.
  if (nrow(res) < length(starts)) {
    stop(sprintf(
      "incomplete extraction from %s: %d keyword rules in the source, %d rows produced. Missing: %s",
      basename(path), length(starts), nrow(res),
      paste(setdiff(sub(rule, "\\2", lines[starts]), res$keyword), collapse = ", ")
    ), call. = FALSE)
  }

  res
}

hazard  <- parse_lex_keywords(file.path(hazard_repo, "src/hazard/hazard_l.l"), "HAZARD")
hazpred <- parse_lex_keywords(file.path(hazard_repo, "src/hazpred/hazpred_l.l"), "HAZPRED")
grammar <- rbind(hazard, hazpred)

# ---------------------------------------------------------------------------
# Curated R-side mapping: where each parser token WILL land in the R API.
#
# Status is deliberately three-valued and none of them is "implemented":
# the translator does not exist yet, so claiming implementation here would be
# exactly the kind of output that looks like a result and is not.
#   mapped          -- R target known, translator not yet built
#   unmapped        -- no R target decided
#   no_r_equivalent -- deliberately out of scope, reason recorded
# ---------------------------------------------------------------------------
r_target <- c(
  DATA        = "data",
  TIME        = "time",
  EVENT       = "status",
  OUTHAZ      = "<saveRDS target>",
  INHAZ       = "hzr_read_outhaz()",
  OUT         = "<result object>",
  CONSERVE    = "control$conserve = TRUE",
  NOCONSERVE  = "control$conserve = FALSE",
  MI          = "control$maxit",
  MAXITER     = "control$maxit",
  CONDITION   = "control$condition",
  STEEPEST    = "control$method",
  QUASINEWTON = "control$method",
  NOPRINT     = "<suppress print>",
  NOCOR       = "<print-only: no R effect>",
  NOCOV       = "<print-only: no R effect>",
  NOLOG       = "<print-only: no R effect>",
  NONOTES     = "<print-only: no R effect>",
  MUE         = "theta (early scale)",
  MUC         = "theta (constant scale)",
  MUL         = "theta (late scale)",
  THALF       = "hzr_phase(t_half=)",
  NU          = "hzr_phase(nu=)",
  M           = "hzr_phase(m=)",
  TAU         = "hzr_phase('g3', tau=)",
  GAMMA       = "hzr_phase('g3', gamma=)",
  ALPHA       = "hzr_phase('g3', alpha=)",
  ETA         = "hzr_phase('g3', eta=)",
  FIXM        = "hzr_phase(fixed=)",
  FIXNU       = "hzr_phase(fixed=)",
  FIXTHALF    = "hzr_phase(fixed=)",
  FIXTAU      = "hzr_phase(fixed=)",
  FIXGAMMA    = "hzr_phase(fixed=)",
  FIXALPHA    = "hzr_phase(fixed=)",
  FIXETA      = "hzr_phase(fixed=)",
  FIXDELTA    = "hzr_phase(fixed=)",
  EARLY       = "hzr_phase('cdf', formula=)",
  CONSTANT    = "hzr_phase('constant', formula=)",
  LATE        = "hzr_phase('g3', formula=)",
  ICENSOR     = "status (2) + time_lower/time_upper",
  LCENSOR     = "status (-1)",
  RCENSOR     = "status (0, already carried by EVENT)",
  SELECTION   = "hzr_stepwise()",
  SLENTRY     = "hzr_stepwise(slentry=)",
  SLSTAY      = "hzr_stepwise(slstay=)",
  FORWARD     = "hzr_stepwise(direction=)",
  BACKWARD    = "hzr_stepwise(direction=)",
  STEPWISE    = "hzr_stepwise(direction=)",
  WEIGHT      = "weights",
  PARAMETERS  = "theta + phases",
  PRINTIT     = "<print control: no R effect>"
)

# Keywords deliberately out of scope, with the reason recorded rather than
# left to be rediscovered.
no_r_equivalent <- c(
  JOBID = "SAS job bookkeeping", JOBIX = "SAS job bookkeeping",
  HAZCOUNT = "SAS listing control", OBSCOUNT = "SAS listing control",
  PROC = "SAS syntax", HAZARD = "SAS syntax", HAZPRED = "SAS syntax"
)

grammar$r_target <- unname(r_target[grammar$token])
grammar$implementation_status <- ifelse(
  grammar$token %in% names(no_r_equivalent), "no_r_equivalent",
  ifelse(is.na(grammar$r_target), "unmapped", "mapped")
)
grammar$note <- unname(no_r_equivalent[grammar$token])
grammar <- grammar[order(grammar$proc, grammar$context, grammar$keyword), ]
rownames(grammar) <- NULL

.hzr_sas_grammar <- grammar
save(.hzr_sas_grammar, file = "R/sysdata.rda", version = 2, compress = "xz")

cat(sprintf("keywords            %d (%d distinct parser tokens)\n",
            nrow(grammar), length(unique(grammar$token))))
cat(sprintf("  aliases           %d\n", sum(grammar$is_alias)))
cat(sprintf("  HAZARD / HAZPRED  %d / %d\n",
            sum(grammar$proc == "HAZARD"), sum(grammar$proc == "HAZPRED")))
cat("\nstatus by distinct token:\n")
per_token <- grammar[!duplicated(grammar$token), ]
print(table(per_token$implementation_status))
cat("\nwritten: R/sysdata.rda\n")
