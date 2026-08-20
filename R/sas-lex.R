# sas-lex.R -- SAS comment stripping and normalisation for the HAZARD job
# translator.
#
# Ported from tools/sas-hazard-profile.R PART 1. The prototype ran clean over
# 990 files across four corpora, but that run predates the fixes below, so it
# is evidence for the general shape of this logic, not a re-validation of it.
# HAZARD's own comment rule (<STMT>\*[^;]*; in hazard_l.l) is quote-agnostic:
# a `* ... ;` comment ends at the first literal `;`, apostrophes and all. Both
# comment-termination call sites here -- the line-initial pass in
# .hzr_sas_strip_comments() and the mid-line pass in
# .hzr_sas_strip_inline_comments() -- must use a plain first-`;` search
# (.idx()), never the quote-aware .first_semi(), to match that rule. An
# earlier revision used .first_semi() in .hzr_sas_strip_inline_comments() and
# silently dropped any statement following a comment that contained an
# apostrophe; see test-sas-lex.R and task-1-report.md for the case that
# surfaced this and the reasoning. Keep this file and
# tools/sas-hazard-profile.R PART 1 in sync on this point -- the profiler's
# measurements are the evidence base for this parser's design.
#
# A second deviation from the prototype: .hzr_sas_normalise() trims
# leading/trailing whitespace from the assembled string (the prototype's
# sas_normalise() did not), so a source that opens with a comment-only line
# does not leave a stray leading space in the returned string.

#' @noRd
.idx <- function(hay, needle) {
  p <- regexpr(needle, hay, fixed = TRUE)
  if (p == -1L) 0L else as.integer(p)
}

#' @noRd
.re <- function(hay, pattern) {
  p <- regexpr(pattern, hay)
  if (p == -1L) 0L else as.integer(p)
}

#' First `;` that is not inside a quoted string. SAS strings use `'` or `"`.
#'
#' For splitting *statements*, where a `;` embedded in a quoted string literal
#' (e.g. `TITLE 'a; b';`) must not be treated as a statement terminator. This
#' is NOT for finding where a `* ... ;` comment ends -- HAZARD's comment rule
#' is quote-agnostic, so using this for comment termination is a bug (see the
#' file header). Nothing in this package calls this function today; it is
#' kept for statement-splitting logic that may need it later.
#'
#' @noRd
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
#'
#' @noRd
.hzr_sas_strip_comments <- function(lines) {
  out <- character(length(lines))
  in_block <- FALSE
  in_star <- FALSE

  for (i in seq_along(lines)) {
    line <- lines[i]

    repeat {
      if (in_block) {
        p <- .idx(line, "*/")
        if (p == 0L) {
          line <- ""
          break
        }
        line <- substring(line, p + 2L)
        in_block <- FALSE
      } else {
        p <- .idx(line, "/*")
        if (p == 0L) break
        rest <- substring(line, p + 2L)
        q <- .idx(rest, "*/")
        if (q == 0L) {
          line <- substring(line, 1L, p - 1L)
          in_block <- TRUE
          break
        }
        line <- paste0(substring(line, 1L, p - 1L), " ", substring(rest, q + 2L))
      }
    }

    # A `* ... ;` statement comment is terminated by the next literal `;`,
    # with no quote-balance check -- HAZARD's own lexer rule is the
    # quote-agnostic <STMT>\*[^;]*; (hazard_l.l). .hzr_sas_strip_inline_comments()
    # below applies the same plain-`;` rule for the mid-line form of this
    # comment; neither pass may treat an apostrophe in comment prose as an
    # unclosed string.
    if (in_star) {
      p <- .idx(line, ";")
      if (p == 0L) {
        line <- ""
      } else {
        line <- substring(line, p + 1L)
        in_star <- FALSE
      }
    }

    if (grepl("^[[:space:]]*%?\\*", line)) {
      p <- .idx(line, ";")
      if (p == 0L) {
        line <- ""
        in_star <- TRUE
      } else {
        line <- substring(line, p + 1L)
      }
    }

    out[i] <- line
  }
  out
}

#' Strip `* ... ;` comments that begin mid-line, immediately after a `;`.
#'
#' The line-based pass in .hzr_sas_strip_comments() only catches the
#' line-initial form. In these jobs a comment frequently follows a statement
#' on the same line, which would otherwise be tokenised as a statement
#' keyword. The comment's own terminating `;` is found with a plain
#' first-`;` search (.idx()), not the quote-aware .first_semi() -- HAZARD's
#' lexer rule (`<STMT>\*[^;]*;`) has no quote awareness, so an apostrophe in
#' comment prose (e.g. "patient's") must not be read as an unclosed string.
#'
#' @noRd
.hzr_sas_strip_inline_comments <- function(txt) {
  repeat {
    p <- .re(txt, "; *%?[*]")
    if (p == 0L) break
    rest <- substring(txt, p + 1L)
    q <- .idx(rest, ";")
    if (q == 0L) {
      txt <- substring(txt, 1L, p)
      break
    }
    txt <- paste0(substring(txt, 1L, p), substring(rest, q + 1L))
  }
  txt
}

#' Normalise stripped source into one uppercase, single-spaced string.
#'
#' @noRd
.hzr_sas_normalise <- function(lines) {
  lines <- .hzr_sas_strip_comments(lines)
  s <- paste(lines, collapse = " ")
  s <- toupper(s)
  s <- gsub("[[:space:]]+", " ", s)
  trimws(.hzr_sas_strip_inline_comments(s))
}

#' Extract PROC HAZARD / PROC HAZPRED blocks from normalised source.
#'
#' Blocks are delimited by parentheses, not by the macro call's name: HAZARD's
#' own lexer treats `)` as whitespace, so a parenthesised group -- not the
#' `%HAZARD(`/`%HAZPRED(` spelling -- is what owns delimitation. Anchoring on
#' the literal macro call text missed macro-argument fragments: files that are
#' `%INCLUDE`d into a `%HAZARD(...)` call held in a different file, so the
#' fragment itself begins with a bare `(` and the macro name is never present
#' in this text at all. The fix locates each `PROC HAZARD `/`PROC HAZPRED `
#' occurrence directly, then scans *backwards* for the nearest unmatched `(`
#' that opens the group containing it, then forwards from that paren to its
#' balancing `)` (unchanged balanced-paren logic -- the corpus confirms 0
#' unterminated). A PROC with no enclosing paren at all is still returned,
#' never dropped: its text is bounded at the next `PROC `/`DATA `/`RUN;`
#' boundary, whichever comes first. If none of those follow, the block
#' extends to the end of the normalised text -- this is safe because this
#' function only ever runs on output from `.hzr_sas_normalise()`, which has
#' already stripped all comments, so there is no trailing comment prose left
#' to sweep in.
#' @noRd
.hzr_sas_blocks <- function(txt) {
  out <- list()
  procs <- gregexpr("PROC HAZ(ARD|PRED) ", txt)[[1L]]
  if (procs[1L] == -1L) return(out)
  proc_lens <- attr(procs, "match.length")

  for (k in seq_along(procs)) {
    proc_at <- procs[k]
    proc <- if (identical(substring(txt, proc_at, proc_at + 10L), "PROC HAZARD")) {
      "HAZARD"
    } else {
      "HAZPRED"
    }

    # Scan backwards from the PROC keyword for the nearest unmatched `(` --
    # the paren that opens the group containing this PROC statement.
    open_at <- NA_integer_
    balance <- 0L
    if (proc_at > 1L) {
      for (i in seq(proc_at - 1L, 1L)) {
        ch <- substr(txt, i, i)
        if (ch == ")") {
          balance <- balance + 1L
        } else if (ch == "(") {
          if (balance == 0L) {
            open_at <- i
            break
          }
          balance <- balance - 1L
        }
      }
    }

    if (is.na(open_at)) {
      # No enclosing paren anywhere before this PROC. Bound the text at the
      # next PROC / DATA / RUN; boundary, never dropping the block. If none
      # of those follow either, the block genuinely extends to the end of
      # txt -- that is safe here because .hzr_sas_blocks() only ever sees
      # output from .hzr_sas_normalise(), which has already stripped all
      # comments, so the end-of-file-prose hazard this bounding rule exists
      # to prevent cannot arise at this point.
      search_from <- proc_at + proc_lens[k]
      rest <- substring(txt, search_from)
      b <- regexpr("PROC |DATA |RUN;", rest)
      body <- if (b == -1L) {
        substring(txt, proc_at)
      } else {
        substring(txt, proc_at, search_from + b - 2L)
      }
      term <- "none"
    } else {
      depth <- 0L
      close_at <- NA_integer_
      for (i in seq(open_at, nchar(txt))) {
        ch <- substr(txt, i, i)
        if (ch == "(") depth <- depth + 1L
        if (ch == ")") {
          depth <- depth - 1L
          if (depth == 0L) {
            close_at <- i
            break
          }
        }
      }
      term <- if (is.na(close_at)) "none" else "paren"
      body <- substring(txt, open_at + 1L,
                        if (is.na(close_at)) nchar(txt) else close_at - 1L)
    }

    out[[length(out) + 1L]] <- list(proc = proc, text = trimws(body),
                                    terminator = term)
  }
  out
}
