# sas-lex.R -- SAS comment stripping and normalisation for the HAZARD job
# translator.
#
# Ported from tools/sas-hazard-profile.R PART 1, which ran clean over 990
# files across four corpora and matches HAZARD's own comment rule
# (<STMT>\*[^;]*; in hazard_l.l). Do not change the logic here without
# re-validating against that corpus.
#
# One deliberate deviation from the prototype: comment termination in
# .hzr_sas_strip_comments() uses a plain next-`;` search rather than the
# quote-aware .first_semi() the prototype reused for that purpose. The
# reference lexer's own rule (<STMT>\*[^;]*;) has no quote awareness, so a
# comment containing a single apostrophe (e.g. "patient's") does not swallow
# the next statement. See test-sas-lex.R and task-1-report.md for the case
# that surfaced this and the reasoning.
#
# A second deviation: .hzr_sas_normalise() trims leading/trailing whitespace
# from the assembled string (the prototype's sas_normalise() did not), so a
# source that opens with a comment-only line does not leave a stray leading
# space in the returned string.

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

# First `;` that is not inside a quoted string. SAS strings use ' or ".
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
    # quote-agnostic <STMT>\*[^;]*; (hazard_l.l). Unlike .hzr_sas_strip_inline_comments(),
    # which uses the quote-aware .first_semi() to find the semicolon that
    # ends the *statement following* a comment, comment termination itself
    # must not treat an apostrophe in comment prose as an unclosed string.
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
#' keyword.
#'
#' @noRd
.hzr_sas_strip_inline_comments <- function(txt) {
  repeat {
    p <- .re(txt, "; *%?[*]")
    if (p == 0L) break
    rest <- substring(txt, p + 1L)
    q <- .first_semi(rest)
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
  s <- trimws(s)
  trimws(.hzr_sas_strip_inline_comments(s))
}
