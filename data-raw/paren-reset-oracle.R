# Regenerate tests/testthat/fixtures/paren-reset-oracle.csv (#461).
#
# PROC HAZARD's lexer clears its syntax-error flag at every `(`
# (hazard_l.l:56: `\(  { BEGIN HZRP; yysynerr = 0; yylnctr = 1; }`), and the
# rule has no start condition, so it is active in every lexer state. A syntax
# error raised BEFORE a later `(` therefore no longer stops the job at
# initprz.c:75-77. This grid runs each refusal class the translator emits
# alone, then with a `(` after the offending text, then with it before, and
# records what PROC HAZARD does. Run from the package root:
#
#   Rscript data-raw/paren-reset-oracle.R [path/to/hazard]
#
# The data are the package's own `avc`, written as the XPORT file the binary
# reads from $TMPDIR/hzr.<JOBID>.<JOBIX>.dta (hazard.sas, opnfils.c).
#
# The verdict is read from the fatal-exit REASON together with the fitted
# markers, never from one message:
#   - "refused": a SYNTAX or SEMANTIC fatal exit;
#   - "fits":    no fatal exit and all four fitted markers;
#   - "no_fit":  a fatal exit for another reason (the parse was accepted, and
#                the procedure stopped before fitting) and none of the markers.
# Any other combination says nothing reliable and stops the script. A clean
# job can print no estimates with no fatal exit (`NU=1` with `EARLY LOG()` on
# avc fails its start, #471), so every class carries a clean job measured to
# fit as its known positive.

args <- commandArgs(TRUE)
bin <- if (length(args)) args[[1L]] else
  path.expand("~/Documents/GitHub/hazard/dist/bin/hazard")
stopifnot(file.exists(bin), file.exists("data/avc.rda"))
work <- tempfile("paren-reset-oracle-")
dir.create(work)

env <- new.env()
load("data/avc.rda", envir = env)
a <- env$avc[stats::complete.cases(env$avc), ]
D <- data.frame(TT = a$int_dead, DEAD = a$dead,
                AGE = as.numeric(scale(a$age)),
                LOG = as.numeric(scale(a$opmos)), SEX = a$mal)
n <- nrow(D)

B <- "EVENT DEAD; TIME TT;"
P <- "PARMS MUE=0.2 THALF=1"
job <- function(proc = "", rest) paste0("PROC HAZARD DATA=D", proc, "; ", rest)
# class, where the `(` sits relative to the offending text, the job.
grid <- list(
  # Known positives: clean jobs that fit on avc.
  c("clean", "none", job(rest = paste(B, P, "; EARLY AGE;"))),
  c("clean", "none", job(rest = paste(B, P, "; EARLY LOG();"))),
  c("clean", "none", job(rest = paste(B, P, "; EARLY AGE, SEX, LOG();"))),
  # PARMS lexer and grammar refusals.
  c("parms", "alone",  job(rest = paste(B, P, "NU=ABC; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, P, "NU=ABC; EARLY LOG();"))),
  c("parms", "before", job(rest = paste(B, "EARLY LOG();", P, "NU=ABC;"))),
  c("parms", "alone",  job(rest = paste(B, P, "NU=1E-3; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, P, "NU=1E-3; EARLY LOG();"))),
  c("parms", "before", job(rest = paste(B, "EARLY LOG();", P, "NU=1E-3;"))),
  c("parms", "alone",  job(rest = paste(B, P, "FIXG1; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, P, "FIXG1; EARLY LOG();"))),
  c("parms", "alone",  job(rest = paste(B, P, "NU; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, P, "NU; EARLY LOG();"))),
  c("parms", "alone",  job(rest = paste(B, P, "FIXNU=1; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, P, "FIXNU=1; EARLY LOG();"))),
  c("parms", "alone",  job(rest = paste(B, P, "NU=?; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, P, "NU=?; EARLY LOG();"))),
  c("parms", "alone",  job(rest = paste(B, "PARMS MUE=0.2 = THALF=1; EARLY AGE;"))),
  c("parms", "after",  job(rest = paste(B, "PARMS MUE=0.2 = THALF=1; EARLY LOG();"))),
  # The phase-statement stop (#340).
  c("phase", "alone",  job(rest = paste(B, P, "; EARLY AGE=ABC;"))),
  c("phase", "after",  job(rest = paste(B, P, "; EARLY AGE=ABC; EARLY LOG();"))),
  c("phase", "same",   job(rest = paste(B, P, "; EARLY AGE=ABC, LOG();"))),
  c("phase", "before", job(rest = paste(B, P, "; EARLY LOG(); EARLY AGE=ABC;"))),
  c("phase", "alone",  job(rest = paste(B, P, "; EARLY AGE/X;"))),
  c("phase", "same",   job(rest = paste(B, P, "; EARLY AGE/X, LOG();"))),
  c("phase", "alone",  job(rest = paste(B, P, "; EARLY /I, AGE;"))),
  c("phase", "same",   job(rest = paste(B, P, "; EARLY /I, AGE, LOG();"))),
  # ORDER= with /E is SEMANTIC (przconc.c:45-53), not a syntax error.
  c("semantic", "alone", job(rest = paste(B, P, "; EARLY SEX/E ORDER=1, AGE;"))),
  c("semantic", "after", job(rest = paste(B, P, "; EARLY SEX/E ORDER=1; EARLY LOG();"))),
  # PROC-line refusals.
  c("proc", "alone", job(" MAXITER=5.", paste(B, P, "; EARLY AGE;"))),
  c("proc", "after", job(" MAXITER=5.", paste(B, P, "; EARLY LOG();"))),
  c("proc", "alone", job(" MAXITER=", paste(B, P, "; EARLY AGE;"))),
  c("proc", "after", job(" MAXITER=", paste(B, P, "; EARLY LOG();"))),
  c("proc", "alone", job(" CONDITION=1E5", paste(B, P, "; EARLY AGE;"))),
  c("proc", "after", job(" CONDITION=1E5", paste(B, P, "; EARLY LOG();"))),
  c("proc", "alone", job(" = MAXITER=50", paste(B, P, "; EARLY AGE;"))),
  c("proc", "after", job(" = MAXITER=50", paste(B, P, "; EARLY LOG();"))),
  c("proc", "alone", job(" NOCOV=1", paste(B, P, "; EARLY AGE;"))),
  c("proc", "after", job(" NOCOV=1", paste(B, P, "; EARLY LOG();"))),
  # A one-name statement with other than one operand (#431).
  c("statement", "alone", job(rest = paste("EVENT DEAD SEX; TIME TT;", P, "; EARLY AGE;"))),
  c("statement", "after", job(rest = paste("EVENT DEAD SEX; TIME TT;", P, "; EARLY LOG();"))),
  # SETG3, SETG1 and the no-phase refusal are raised at fit time.
  c("setg", "alone", job(rest = paste(B, "PARMS MUL=0.2 TAU=0 FIXTAU GAMMA=1 ETA=1; LATE AGE;"))),
  c("setg", "after", job(rest = paste(B, "PARMS MUL=0.2 TAU=0 FIXTAU GAMMA=1 ETA=1; LATE LOG();"))),
  c("setg", "alone", job(rest = paste(B, "PARMS MUE=0.2 THALF=0 FIXTHALF; EARLY AGE;"))),
  c("setg", "after", job(rest = paste(B, "PARMS MUE=0.2 THALF=0 FIXTHALF; EARLY LOG();"))),
  c("setg", "alone", job(rest = paste(B, "PARMS MUE=0 THALF=1; EARLY AGE;"))),
  c("setg", "after", job(rest = paste(B, "PARMS MUE=0 THALF=1; EARLY LOG();"))),
  # A phase variable that is not a NAME (#440).
  c("name", "alone",  job(rest = paste(B, P, "; EARLY AGE*SEX;"))),
  c("name", "same",   job(rest = paste(B, P, "; EARLY AGE*SEX, LOG();"))),
  c("name", "after",  job(rest = paste(B, P, "; EARLY AGE*SEX; EARLY LOG();"))),
  c("name", "before", job(rest = paste(B, P, "; EARLY LOG(), AGE*SEX;"))),
  c("name", "before", job(rest = paste(B, P, "; EARLY LOG(); EARLY AGE*SEX;"))),
  c("name", "after",  job(rest = paste(B, P, "; EARLY LOG(); EARLY AGE*SEX; EARLY LOG();"))),
  c("name", "alone",  job(rest = paste(B, P, "; EARLY 1AGE, SEX;"))),
  c("name", "same",   job(rest = paste(B, P, "; EARLY 1AGE, SEX, LOG();"))),
  # Which `(` characters clear the flag, and what after one sets it again.
  c("name", "same",   job(rest = paste(B, P, "; EARLY AGE*SEX, LOG( );"))),
  c("name", "same",   job(rest = paste(B, P, "; EARLY AGE*SEX, LOG() = 0.1;"))),
  c("name", "same",   job(rest = paste(B, P, "; EARLY AGE*SEX, LOG() /I;"))),
  c("name", "same",   job(rest = paste(B, P, "; EARLY AGE*SEX, (LOG);"))),
  # A statement PROC HAZARD does not know, after the `(`, sets the flag again.
  c("unknown", "after", job(rest = paste(B, P, "; EARLY AGE*SEX, LOG(); FOO BAR;")))
)

marks <- c("Estimates for Model Parameters", "Final Results", "Log likelihood",
           "Optimization terminated")
# The covariates PROC HAZARD lists under "Concomitant Information".
concomitant <- function(out) {
  a <- grep("Concomitant Information", out, fixed = TRUE)[1L]
  b <- grep("Model Specifications", out, fixed = TRUE)[1L]
  if (is.na(a) || is.na(b)) return("")
  seg <- out[a:b]
  k <- grep("^@31,", seg)
  nm <- trimws(sub("^\\+1,", "", seg[k + 1L]))
  paste(nm[grepl("^[A-Z_][A-Z0-9_]*$", nm) & !nm %in% c("Phase", "Variable")],
        collapse = " ")
}

rows <- lapply(seq_along(grid), function(k) {
  g <- grid[[k]]
  ix <- sprintf("X%d", k)
  haven::write_xpt(D, file.path(work, sprintf("hzr.J461.%s.dta", ix)),
                   version = 5, name = "HZRCALL")
  ctl <- file.path(work, sprintf("hzr.J461.%s.sas", ix))
  lst <- file.path(work, sprintf("hzr.J461.%s.lst", ix))
  writeLines(c(paste0("(   ", g[[3L]], "   )"),
               sprintf("; OBSCOUNT %d ;", n), "JOBID J461;",
               sprintf("JOBIX %s;", ix)), ctl)
  system(paste0("TMPDIR=", shQuote(work), " ", shQuote(bin), " < ",
                shQuote(ctl), " > ", shQuote(lst), " 2>&1"))
  out <- readLines(lst, warn = FALSE)
  fatal <- regmatches(out, regexpr("fatal exit: [A-Za-z ]+ at", out))
  reason <- if (length(fatal)) {
    trimws(sub(" at$", "", sub("fatal exit: ", "", fatal[[1L]])))
  } else {
    ""
  }
  n_marks <- sum(vapply(marks, function(m) any(grepl(m, out, fixed = TRUE)),
                        logical(1L)))
  verdict <- if (reason %in% c("SYNTAX", "SEMANTIC")) "refused" else
    if (!nzchar(reason) && n_marks == 4L) "fits" else
      if (nzchar(reason) && n_marks == 0L) "no_fit" else NA_character_
  if (is.na(verdict)) {
    stop("row ", k, " has no verdict (", reason, ", ", n_marks, "/4): ",
         g[[3L]])
  }
  data.frame(class = g[[1L]], paren = g[[2L]], job = g[[3L]],
             reason = reason, markers = n_marks, verdict = verdict,
             early_vars = if (verdict == "fits") concomitant(out) else "",
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
version <- grep("C-Version", readLines(file.path(work, "hzr.J461.X1.lst")),
                value = TRUE)[[1L]]
dest <- "tests/testthat/fixtures/paren-reset-oracle.csv"
writeLines(c(
  "# Generated by data-raw/paren-reset-oracle.R; do not edit by hand.",
  paste0("# Binary: ", sub(path.expand("~"), "~", bin, fixed = TRUE), " (",
         trimws(sub("^\\$Note:", "", version)), ")"),
  paste0("# Binary file date: ", format(file.mtime(bin), "%Y-%m-%d")),
  paste0("# Source checkout when generated: ",
         system2("git", c("-C", dirname(dirname(dirname(bin))), "rev-parse",
                          "HEAD"), stdout = TRUE)),
  "# Data: the package's avc (complete cases), written as XPORT.",
  paste0("# Generated: ", format(Sys.Date()))
), dest)
suppressWarnings(utils::write.table(tab, dest, append = TRUE, sep = ",",
                                    row.names = FALSE, qmethod = "double"))
print(tab[, c("class", "paren", "reason", "markers", "verdict", "early_vars")],
      row.names = FALSE)
