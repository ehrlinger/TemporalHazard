#!/usr/bin/env Rscript
# sas-hazard-profile.R -- structural profile of a SAS HAZARD/HAZPRED job tree.
#
# PURPOSE
#   Measure the shape of a corpus of SAS hazard jobs so that a SAS -> R/Quarto
#   translator can be designed against real variety rather than against the
#   handful of curated examples in the public `hazard` repository.
#
#   The lexer in PART 1 is deliberately written to become R/sas-lex.R. If the
#   profiler stripped comments differently from the eventual parser, the
#   profile would describe a different language than the parser parses, and
#   nothing would tell you.
#
# WHAT IT EMITS
#   Counts, frequencies, histograms, and shape classes. It does NOT emit file
#   paths, dataset names, variable names, titles, or any line of source.
#   Librefs are salted-hashed with a salt generated per run and never printed.
#
# DEPENDENCIES
#   Base R only. No packages.
#
# USAGE
#   Rscript sas-hazard-profile.R /path/to/sas/tree > profile.txt
#   Rscript sas-hazard-profile.R --show-librefs /path/to/tree   # local only
#
#   Review profile.txt before sharing it.
#
# EXIT STATUS
#   0 report produced   1 usage error   2 nothing matched (see FAIL LOUD)

# ============================ PART 1: the lexer =============================
# Prototype of R/sas-lex.R. Keep these three functions in sync with it.

.idx <- function(hay, needle) {
  p <- regexpr(needle, hay, fixed = TRUE)
  if (p == -1L) 0L else as.integer(p)
}

.re <- function(hay, pattern) {
  p <- regexpr(pattern, hay)
  if (p == -1L) 0L else as.integer(p)
}

# First occurrence of a PROC/DATA boundary, skipping the one at position 1.
.nth_boundary <- function(hay, needle) {
  p <- .idx(substring(hay, 2L), needle)
  if (p == 0L) 0L else p + 1L
}

# First `;` that is not inside a quoted string. SAS strings use ' or ".
.first_semi <- function(line) {
  pos <- gregexpr(";", line, fixed = TRUE)[[1L]]
  if (pos[1L] == -1L) return(0L)
  for (p in pos) {
    head <- substr(line, 1L, p - 1L)
    nsq <- nchar(gsub("[^']", "", head))
    ndq <- nchar(gsub('[^"]', "", head))
    if (nsq %% 2L == 0L && ndq %% 2L == 0L) return(as.integer(p))
  }
  0L
}

#' Strip SAS comments from a character vector of source lines.
#'
#' Handles `/* ... */` block comments spanning lines, and `* ... ;` / `%* ... ;`
#' statement comments (line-initial form) which also span lines. Both matter:
#' commented-out PARMS statements are common in these jobs and would otherwise
#' inflate every count, and comment prose is where study and patient detail live.
sas_strip_comments <- function(lines) {
  out <- character(length(lines))
  in_block <- FALSE
  in_star <- FALSE

  for (i in seq_along(lines)) {
    line <- lines[i]

    repeat {
      if (in_block) {
        p <- .idx(line, "*/")
        if (p == 0L) { line <- ""; break }
        line <- substring(line, p + 2L)
        in_block <- FALSE
      } else {
        p <- .idx(line, "/*")
        if (p == 0L) break
        rest <- substring(line, p + 2L)
        q <- .idx(rest, "*/")
        if (q == 0L) { line <- substring(line, 1L, p - 1L); in_block <- TRUE; break }
        line <- paste0(substring(line, 1L, p - 1L), " ", substring(rest, q + 2L))
      }
    }

    if (in_star) {
      p <- .first_semi(line)
      if (p == 0L) line <- "" else { line <- substring(line, p + 1L); in_star <- FALSE }
    }

    if (grepl("^[[:space:]]*%?\\*", line)) {
      p <- .first_semi(line)
      if (p == 0L) { line <- ""; in_star <- TRUE } else line <- substring(line, p + 1L)
    }

    out[i] <- line
  }
  out
}

#' Strip `* ... ;` comments that begin mid-line, immediately after a `;`.
#'
#' The line-based pass in sas_strip_comments() only catches the line-initial
#' form. In these jobs a comment frequently follows a statement on the same
#' line, which would otherwise be tokenised as a statement keyword.
sas_strip_inline_comments <- function(txt) {
  repeat {
    p <- .re(txt, "; *%?[*]")
    if (p == 0L) break
    rest <- substring(txt, p + 1L)
    q <- .first_semi(rest)
    if (q == 0L) { txt <- substring(txt, 1L, p); break }
    txt <- paste0(substring(txt, 1L, p), substring(rest, q + 1L))
  }
  txt
}

#' Normalise stripped source into one uppercase, single-spaced string.
sas_normalise <- function(lines) {
  s <- paste(lines, collapse = " ")
  s <- toupper(s)
  s <- gsub("[[:space:]]+", " ", s)
  sas_strip_inline_comments(s)
}

#' Extract PROC HAZARD / PROC HAZPRED blocks from normalised source.
#'
#' Returns a list of list(type=, text=, terminator=). `terminator` records how
#' the block was closed so that a fallback is counted, never hidden.
sas_extract_haz_blocks <- function(txt) {
  blocks <- list()
  s <- txt
  repeat {
    m <- regexpr("PROC HAZ(ARD|PRED) ", s)
    if (m == -1L) break
    tail <- substring(s, m)
    type <- if (grepl("^PROC HAZARD", tail)) "HAZARD" else "HAZPRED"

    # Bound the block at the EARLIEST of: the macro close `)` `;` (whitespace
    # may separate them after normalisation), RUN;, or the next PROC/DATA
    # boundary. A block must NEVER run to end of file: that is how comment
    # prose reaches the token tables, and comment prose is where study and
    # patient detail live.
    cands <- c(
      paren = .re(tail, "\\) *;"),
      run   = .idx(tail, "RUN;"),
      proc  = .nth_boundary(tail, " PROC "),
      data  = .nth_boundary(tail, " DATA ")
    )
    cands <- cands[cands > 0L]
    if (length(cands)) {
      e <- min(cands); term <- names(cands)[which.min(cands)]
    } else {
      e <- nchar(tail) + 1L; term <- "none"
    }

    blocks[[length(blocks) + 1L]] <- list(
      type = type, text = substring(tail, 1L, e - 1L), terminator = term
    )
    s <- substring(tail, e + 1L)
  }
  blocks
}

# ============================ PART 2: redaction =============================
# POLICY (edit here if your boundary differs):
#   librefs, dataset names, variable names, file paths -> never printed raw.
#   SAS keywords and option names -> printed raw; they are language, not data.
.SALT <- paste0(Sys.getpid(), "-", as.numeric(Sys.time()), "-",
                paste(sample(c(letters, 0:9), 16, TRUE), collapse = ""))

redact <- function(s) {
  x <- utf8ToInt(paste0(.SALT, "|", s))
  h <- 5381
  for (ch in x) h <- (h * 33 + ch) %% 99999989
  sprintf("%05x", h %% 1048576L)
}

# ============================ PART 3: known tokens ==========================
KNOWN_OPT <- list(
  HAZARD = c("DATA", "OUTHAZ", "OUT", "CONSERVE", "P", "NOPRINT", "STEEPEST",
             "QUASI", "NEWTON", "MARQUARDT", "CONDITION", "MAXIT", "MAXITER",
             "CONVERGE", "SLENTRY", "SLSTAY", "SCORE", "WALD", "WEIGHT",
             "SEED", "NPHASE", "ALPHA", "COVOUT"),
  HAZPRED = c("DATA", "INHAZ", "OUT", "NOPRINT", "P", "ALPHA")
)
KNOWN_STMT <- c("TIME", "EVENT", "PARMS", "EARLY", "CONSTANT", "LATE",
                "SELECTION", "VAR", "BY", "STRATA", "ID", "OUTPUT", "MODEL",
                "WHERE", "FREQ", "WEIGHT", "LABEL", "FORMAT")
KNOWN_PARM <- c("MUE", "MUC", "MUL", "MU", "THALF", "NU", "M", "DELTA", "TAU",
                "GAMMA", "ALPHA", "ETA", "FIXM", "FIXNU", "FIXDELTA", "FIXTAU",
                "FIXGAMMA", "FIXALPHA", "FIXETA", "FIXTHALF")

# ============================ PART 4: accumulate ============================
acc <- new.env(parent = emptyenv())
acc$opt <- list(); acc$stmt <- list(); acc$parm <- list()
acc$unknown_opt <- list(); acc$unknown_stmt <- list(); acc$unknown_parm <- list()
acc$covar <- list(); acc$grid_class <- list()
acc$outhaz <- character(0); acc$inhaz <- character(0)
acc$n_blocks <- c(HAZARD = 0L, HAZPRED = 0L)
acc$term <- c(paren = 0L, run = 0L, proc = 0L, data = 0L, none = 0L)
acc$files_parsed <- 0L; acc$files_with_blocks <- 0L
acc$grid_ds <- character(0)

bump <- function(slot, key, by = 1L) {
  cur <- acc[[slot]]
  cur[[key]] <- (if (is.null(cur[[key]])) 0L else cur[[key]]) + by
  acc[[slot]] <- cur
}

parse_block <- function(block, type) {
  st <- strsplit(block, ";", fixed = TRUE)[[1L]]
  if (!length(st)) return(invisible(NULL))

  toks <- strsplit(trimws(st[1L]), " ", fixed = TRUE)[[1L]]
  toks <- toks[nzchar(toks) & !(toks %in% c("PROC", "HAZARD", "HAZPRED"))]
  for (tok in toks) {
    eqp <- .idx(tok, "=")
    key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
    val <- if (eqp > 0L) substring(tok, eqp + 1L) else ""
    bump("opt", paste0(type, "|", key))
    if (!(key %in% KNOWN_OPT[[type]])) bump("unknown_opt", paste0(type, "|", key))
    if (key == "OUTHAZ" && nzchar(val)) acc$outhaz <- c(acc$outhaz, val)
    if (key == "INHAZ" && nzchar(val)) acc$inhaz <- c(acc$inhaz, val)
    if (key == "DATA" && type == "HAZPRED" && nzchar(val))
      acc$grid_ds <- c(acc$grid_ds, val)
  }

  for (i in seq_along(st)[-1L]) {
    w <- strsplit(trimws(st[i]), " ", fixed = TRUE)[[1L]]
    w <- w[nzchar(w)]
    if (!length(w)) next
    kw <- w[1L]
    ops <- w[-1L]
    bump("stmt", kw)
    if (!(kw %in% KNOWN_STMT)) bump("unknown_stmt", kw)
    if (kw == "PARMS") {
      for (tok in ops) {
        eqp <- .idx(tok, "=")
        key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
        bump("parm", key)
        if (!(key %in% KNOWN_PARM)) bump("unknown_parm", key)
      }
    }
    if (kw %in% c("EARLY", "CONSTANT", "LATE"))
      bump("covar", sprintf("%s | %d covariate(s)", kw, length(ops)))
  }
  invisible(NULL)
}

# Classify the DATA step that feeds a HAZPRED DATA= grid. Only the CLASS is
# ever reported -- never the step body, which is where labels and titles live.
classify_grids <- function(txt, grid_names) {
  for (nm in unique(grid_names)) {
    p <- .idx(txt, paste0("DATA ", nm, ";"))
    if (p == 0L) { bump("grid_class", "not_found"); next }
    rest <- substring(txt, p)
    ends <- c(.idx(substring(rest, 6L), "DATA "), .idx(substring(rest, 6L), "PROC "))
    ends <- ends[ends > 0L]
    body <- if (length(ends)) substring(rest, 1L, min(ends) + 5L) else rest
    cls <- if (grepl(" DO ", body) && grepl("LOG(", body, fixed = TRUE)) "log_grid"
           else if (grepl(" DO ", body)) "explicit_do"
           else if (grepl(" SET ", body)) "derived_set"
           else "other"
    bump("grid_class", cls)
  }
}

# ============================ PART 5: driver ================================
main <- function() {
args <- commandArgs(trailingOnly = TRUE)
show_librefs <- "--show-librefs" %in% args
args <- args[args != "--show-librefs"]

if (length(args) != 1L) {
  cat("usage: Rscript sas-hazard-profile.R [--show-librefs] <root-dir>\n",
      file = stderr()); quit(status = 1L)
}
root <- args[1L]
if (!dir.exists(root)) {
  cat("not a directory: ", root, "\n", sep = "", file = stderr()); quit(status = 1L)
}

files <- list.files(root, pattern = "\\.sas$", recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
if (!length(files)) {
  cat("FAIL LOUD: no *.sas files under ", root, "\n",
      "Nothing was profiled. This is a failed scan, not an empty result.\n",
      sep = "", file = stderr())
  quit(status = 2L)
}

for (f in files) {
  lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) next
  acc$files_parsed <- acc$files_parsed + 1L
  txt <- sas_normalise(sas_strip_comments(lines))
  blocks <- sas_extract_haz_blocks(txt)
  if (length(blocks)) acc$files_with_blocks <- acc$files_with_blocks + 1L
  before <- length(acc$grid_ds)
  for (b in blocks) {
    acc$n_blocks[b$type] <- acc$n_blocks[b$type] + 1L
    acc$term[b$terminator] <- acc$term[b$terminator] + 1L
    parse_block(b$text, b$type)
  }
  if (length(acc$grid_ds) > before)
    classify_grids(txt, acc$grid_ds[(before + 1L):length(acc$grid_ds)])
}

# ============================ PART 6: report ================================
# Defence in depth. Even if the lexer leaks, only identifier-shaped tokens are
# ever printed; anything else is bucketed by count alone. This guarantees the
# report cannot emit prose regardless of upstream bugs.
safe_token <- function(x) {
  # Keys may carry a "TYPE|TOKEN" prefix; validate only the token half.
  pre <- ifelse(grepl("|", x, fixed = TRUE), sub("[|].*$", "|", x), "")
  tok <- sub("^.*[|]", "", x)
  ok <- grepl("^[A-Z_&%][A-Z0-9_.]*$", tok)
  paste0(pre, ifelse(ok, tok, "<non-identifier token>"))
}

tbl <- function(slot, prefix = NULL, known = NULL, sanitise = TRUE) {
  cur <- acc[[slot]]
  if (!length(cur)) { cat("  (none)\n"); return(invisible(NULL)) }
  keys <- names(cur)
  if (!is.null(prefix)) {
    keep <- startsWith(keys, prefix)
    cur <- cur[keep]; keys <- substring(keys[keep], nchar(prefix) + 1L)
  }
  if (!length(cur)) { cat("  (none)\n"); return(invisible(NULL)) }
  o <- order(-unlist(cur), keys)
  for (i in o) {
    flag <- if (!is.null(known) && !(keys[i] %in% known)) "   <-- UNKNOWN" else ""
    lab <- if (sanitise) safe_token(keys[i]) else keys[i]
    cat(sprintf("  %-22s %6d%s\n", lab, cur[[i]], flag))
  }
  invisible(NULL)
}

cat("SAS HAZARD/HAZPRED corpus profile\n")
cat("=================================\n\n")
cat(sprintf("files found (*.sas)        %d\n", length(files)))
cat(sprintf("files parsed               %d\n", acc$files_parsed))
cat(sprintf("files with a HAZ* block    %d\n", acc$files_with_blocks))
cat(sprintf("PROC HAZARD  blocks        %d\n", acc$n_blocks[["HAZARD"]]))
cat(sprintf("PROC HAZPRED blocks        %d\n", acc$n_blocks[["HAZPRED"]]))

total_blocks <- sum(acc$n_blocks)
if (total_blocks == 0L) {
  cat("\n*** FAIL LOUD ***\n")
  cat(sprintf("Zero HAZ* blocks matched across %d files.\n", acc$files_parsed))
  cat("Every count below would be a confident zero over nothing.\n")
  cat("Check the root path, and check whether these jobs wrap the PROC\n")
  cat("differently: this extractor keys on the literal 'PROC HAZARD ' /\n")
  cat("'PROC HAZPRED ' after comment stripping and whitespace collapse.\n")
  quit(status = 2L)
}
if (acc$term[["run"]] > 0L)
  cat(sprintf("  note: %d block(s) closed on RUN; rather than the macro paren\n",
              acc$term[["run"]]))
if (acc$term[["proc"]] + acc$term[["data"]] > 0L)
  cat(sprintf("  WARNING: %d block(s) had no macro close; bounded at the next\n           PROC/DATA instead. Their option counts may be short.\n",
              acc$term[["proc"]] + acc$term[["data"]]))
if (acc$term[["none"]] > 0L)
  cat(sprintf("  WARNING: %d block(s) reached end of file unbounded\n",
              acc$term[["none"]]))

cat("\n-- PROC options: PROC HAZARD ---------------------------------\n")
tbl("opt", "HAZARD|", KNOWN_OPT$HAZARD)
cat("\n-- PROC options: PROC HAZPRED --------------------------------\n")
tbl("opt", "HAZPRED|", KNOWN_OPT$HAZPRED)
cat("\n-- statements ------------------------------------------------\n")
tbl("stmt", NULL, KNOWN_STMT)
cat("\n-- PARMS keywords --------------------------------------------\n")
tbl("parm", NULL, KNOWN_PARM)
cat("\n-- covariates per phase statement ----------------------------\n")
tbl("covar", sanitise = FALSE)
cat("\n-- HAZPRED prediction-grid DATA-step shape -------------------\n")
tbl("grid_class", sanitise = FALSE)

cat("\n-- UNKNOWN tokens (the gap this profile exists to measure) ----\n")
n_unknown <- length(acc$unknown_opt) + length(acc$unknown_stmt) +
  length(acc$unknown_parm)
if (n_unknown == 0L) {
  cat("  none.\n")
  cat(sprintf("  This zero is meaningful ONLY because %d blocks parsed above.\n",
              total_blocks))
} else {
  cat("  options:\n");    tbl("unknown_opt")
  cat("  statements:\n"); tbl("unknown_stmt")
  cat("  parms:\n");      tbl("unknown_parm")
}

cat("\n-- INHAZ / OUTHAZ resolution ---------------------------------\n")
cat(sprintf("  OUTHAZ= targets written    %d (%d distinct)\n",
            length(acc$outhaz), length(unique(acc$outhaz))))
cat(sprintf("  INHAZ=  sources read       %d (%d distinct)\n",
            length(acc$inhaz), length(unique(acc$inhaz))))
if (length(acc$inhaz)) {
  ok <- acc$inhaz %in% acc$outhaz
  cat(sprintf("  resolved in-tree           %d (%.1f%%)\n",
              sum(ok), 100 * mean(ok)))
  cat(sprintf("  UNRESOLVED                 %d (%.1f%%)\n",
              sum(!ok), 100 * mean(!ok)))
  bad <- acc$inhaz[!ok]
  if (length(bad)) {
    lib <- ifelse(grepl(".", bad, fixed = TRUE),
                  sub("\\..*$", "", bad), "(WORK)")
    cat("\n  unresolved by libref:\n")
    tb <- sort(table(lib), decreasing = TRUE)
    for (k in names(tb))
      cat(sprintf("    %-14s %6d\n",
                  if (show_librefs) k else redact(k), tb[[k]]))
    if (!show_librefs)
      cat("    (librefs salted-hashed; rerun with --show-librefs to see them locally)\n")
  }
} else {
  cat("  no INHAZ= found -- either no prediction jobs, or they name the fit\n")
  cat("  another way. Given HAZPRED blocks were found, treat this as suspicious.\n")
}

cat("\n-- end -------------------------------------------------------\n")
cat("Review this file before sharing it.\n")
}

if (sys.nframe() == 0L) main()
