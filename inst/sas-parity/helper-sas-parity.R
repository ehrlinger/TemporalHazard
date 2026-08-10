# Shared parser for SAS HAZARD .lst capture files.  testthat auto-sources
# this helper before running test-sas-parity.R.
#
# Reference fixtures live at:
#   ~/Documents/GitHub/hazard/examples/*.lst   (paired with *.sas)
#
# Capture vintage: v4.4 macOS (2026-04-27).  Final values refresh when
# v4.4.6 captures land — the parser itself is version-stable.
#
# Encoding notes (verified across all 18 fixtures):
#   - CRLF line endings  ->  strip \r
#   - SAS form-feed page breaks  ->  convert \f to \n
#   - "Non-ISO extended-ASCII" classification causes grep -c to silently
#     return 0 without -a; R's readLines handles this correctly when the
#     normalized text is split on \n.
#
# Public API:
#   .hzr_parse_sas_lst(path) -> list(fits = list(<one per Final Results:>))
#
#   Each fit element:
#     loglik        scalar, final "Log likelihood = N" before "Final Results:"
#     n_obs         integer
#     n_events      integer
#     n_censored    integer
#     params        data.frame(phase, name, label, estimate, se, z, p)
#                   from "Parameter Estimate Summary"
#     natural       data.frame(phase, name, fixed, estimate)
#                   from "Estimates for Model Parameters"
#     vcov          numeric matrix (NULL if NOCOV)
#     correlation   numeric matrix (NULL if NOCOR)
#     events_conserved  scalar or NA

.hzr_read_lst <- function(path) {
  raw <- readLines(path, warn = FALSE, encoding = "latin1")
  # Convert form-feed page breaks to newlines, then re-split.
  joined <- paste(raw, collapse = "\n")
  joined <- gsub("\f", "\n", joined, fixed = TRUE)
  joined <- gsub("\r", "",  joined, fixed = TRUE)
  strsplit(joined, "\n", fixed = TRUE)[[1]]
}

.hzr_extract_loglik <- function(lines) {
  m <- regmatches(
    lines,
    regexec("Log likelihood\\s*=\\s*(-?[0-9.Ee+-]+)", lines)
  )
  vals <- vapply(
    m,
    function(x) if (length(x) >= 2) as.numeric(x[2]) else NA_real_,
    numeric(1)
  )
  vals[!is.na(vals)]
}

.hzr_extract_obs_counts <- function(lines) {
  # "There are  5880 observations available for analysis with:"
  # "          545 events"
  # "         5335 Right Censored Observations"
  obs_line <- grep("observations available for analysis", lines, value = TRUE)[1]
  ev_line  <- grep("\\bevents\\s*$",                       lines, value = TRUE)[1]
  rc_line  <- grep("Right Censored Observations",          lines, value = TRUE)[1]

  pull <- function(pat, s) {
    if (is.na(s)) return(NA_integer_)
    m <- regmatches(s, regexec(pat, s))[[1]]
    if (length(m) >= 2) as.integer(m[2]) else NA_integer_
  }

  list(
    n_obs      = pull("There are\\s+([0-9]+)\\s+observations",          obs_line),
    n_events   = pull("([0-9]+)\\s+events",                              ev_line),
    n_censored = pull("([0-9]+)\\s+Right Censored Observations",         rc_line)
  )
}

.hzr_extract_events_conserved <- function(lines) {
  s <- grep("Number of events conserved", lines, value = TRUE)
  if (!length(s)) return(NA_real_)
  m <- regmatches(s[1], regexec("=\\s*(-?[0-9.Ee+-]+)", s[1]))[[1]]
  if (length(m) >= 2) as.numeric(m[2]) else NA_real_
}

# Parameter Estimate Summary table.  Layout:
#
#   Phase      Parameter   Label              Estimate    Std error      Z        Prob>|Z|
#   ---------- ...
#   Early:     E0                              -3.77955    0.09381214   -40.288   <.0001
#              ---------------- ...
#   Constant:  C0                              -7.2258     0.09312647   -77.591   <.0001
#
# - Phase label appears only on the first row of each phase; subsequent
#   rows in the same phase have a blank Phase column.
# - Inner "---" separator rows divide phases.
# - Outer "===" / long "---" rule bounds the table.
.hzr_extract_param_summary <- function(lines) {
  start <- grep("Parameter Estimate Summary", lines)
  if (!length(start)) return(NULL)
  start <- start[1]

  # Find the header row (starts with "Phase").
  hdr <- which(grepl("^\\s*Phase\\s+Parameter", lines))
  hdr <- hdr[hdr > start][1]
  if (is.na(hdr)) return(NULL)

  # Body starts after the rule line under the header.
  body_start <- hdr + 2L

  # Body ends at the next blank-rule (long dashes) followed by blank lines
  # before "Estimates for Model Parameters".
  end_marker <- grep("Estimates for Model Parameters", lines)
  end_marker <- end_marker[end_marker > hdr][1]
  if (is.na(end_marker)) return(NULL)

  block <- lines[body_start:(end_marker - 1L)]
  block <- block[nzchar(trimws(block))]
  # Drop separator rows (only dashes / spaces).
  block <- block[!grepl("^[\\s-]+$", block, perl = TRUE)]

  current_phase <- NA_character_
  rows <- list()

  for (ln in block) {
    # Phase label appears at the very start of the line, ending in ":".
    m_phase <- regmatches(ln, regexec("^\\s+(Early|Constant|Late):", ln))[[1]]
    if (length(m_phase) >= 2) {
      current_phase <- m_phase[2]
      ln <- sub("^\\s*(Early|Constant|Late):\\s*", "", ln)
    } else {
      # Continuation row — strip leading whitespace.
      ln <- sub("^\\s+", "", ln)
    }

    # Tail of each row: estimate, se, z, p.  Use a trailing-fields regex
    # so we capture them regardless of label column width.
    # Trailing fields: estimate, se, z, p.  p-value can be "<.0001" or a
    # decimal like "0.0032" / ".0032".
    tail_re <- "(-?[0-9.Ee+-]+)\\s+(-?[0-9.Ee+-]+)\\s+(-?[0-9.Ee+-]+)\\s+(<\\.?[0-9]+|[0-9.]+)\\s*$"
    m_tail <- regmatches(ln, regexec(tail_re, ln))[[1]]
    if (length(m_tail) < 5) next

    # Head: parameter name and optional label.  Strip the matched tail.
    head <- sub(tail_re, "", ln)
    head <- trimws(head)
    # First token is the parameter name; remainder is the label.
    head_split <- regmatches(head, regexec("^(\\S+)\\s*(.*)$", head))[[1]]
    name  <- if (length(head_split) >= 2) head_split[2] else NA_character_
    label <- if (length(head_split) >= 3) trimws(head_split[3]) else ""

    p_str <- m_tail[5]
    p_val <- if (grepl("^<", p_str)) {
      as.numeric(sub("^<", "", p_str))
    } else {
      as.numeric(p_str)
    }

    rows[[length(rows) + 1L]] <- data.frame(
      phase    = current_phase,
      name     = name,
      label    = label,
      estimate = as.numeric(m_tail[2]),
      se       = as.numeric(m_tail[3]),
      z        = as.numeric(m_tail[4]),
      p        = p_val,
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Estimates for Model Parameters table — natural-scale values, all params
# (free + fixed).  Layout:
#
#   Phase       Parameter     Fixed?     Estimate
#   Early:      DELTA         Yes                  0
#               THALF         Yes                0.2
#                (RHO                            0.2)
#               ...
#               MUE                       0.02283304
.hzr_extract_natural <- function(lines) {
  start <- grep("Estimates for Model Parameters", lines)
  if (!length(start)) return(NULL)
  start <- start[1]

  hdr <- which(grepl("^\\s*Phase\\s+Parameter\\s+Fixed\\?", lines))
  hdr <- hdr[hdr > start][1]
  if (is.na(hdr)) return(NULL)

  body_start <- hdr + 2L

  end_marker <- grep("Asymptotic Variance-Covariance Matrix|Asymptotic Correlation Matrix",
                     lines)
  end_marker <- end_marker[end_marker > hdr][1]
  end_idx <- if (is.na(end_marker)) length(lines) else end_marker - 1L

  block <- lines[body_start:end_idx]
  block <- block[nzchar(trimws(block))]
  block <- block[!grepl("^[\\s-]+$", block, perl = TRUE)]

  current_phase <- NA_character_
  rows <- list()
  for (ln in block) {
    # Skip RHO informational rows: " (RHO                            0.2)"
    if (grepl("^\\s*\\(RHO", ln)) next

    # Natural-scale table is deeply indented (~40+ leading spaces); accept
    # arbitrary leading whitespace before the phase label.
    m_phase <- regmatches(ln, regexec("^\\s+(Early|Constant|Late):", ln))[[1]]
    if (length(m_phase) >= 2) {
      current_phase <- m_phase[2]
      ln <- sub("^\\s*(Early|Constant|Late):\\s*", "", ln)
    } else {
      ln <- sub("^\\s+", "", ln)
    }

    # Tail: optional "Yes"/"No", then numeric estimate.
    tail_re <- "^(\\S+)\\s+(?:(Yes|No)\\s+)?(-?[0-9.Ee+-]+)\\s*$"
    m <- regmatches(ln, regexec(tail_re, ln, perl = TRUE))[[1]]
    if (length(m) < 4) next

    rows[[length(rows) + 1L]] <- data.frame(
      phase    = current_phase,
      name     = m[2],
      fixed    = identical(m[3], "Yes"),
      estimate = as.numeric(m[4]),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Symmetric matrix parser shared by vcov and correlation blocks.
#
# Wide matrices (>~8 params) wrap both the header row and each data row
# across multiple lines.  Pages are also broken by form-feed inserting a
# title block, so we collect column names from any line of header-shape
# until a separator rule appears, and we accumulate numerics from
# continuation lines until the next labelled row or section.
.hzr_extract_matrix <- function(lines, marker) {
  start <- grep(marker, lines)
  if (!length(start)) return(NULL)
  start <- start[1]

  # Walk forward to collect the (possibly multi-line) header.  Header lines
  # contain only parameter names (identifiers) separated by 2+ spaces;
  # they terminate at the rule line of dashes.
  i <- start + 1L
  col_names <- character()
  saw_header <- FALSE
  while (i <= length(lines)) {
    ln <- lines[i]
    if (grepl("^\\s*-{10,}", ln)) {
      if (saw_header) break
    } else if (grepl("^\\s+[A-Za-z_][A-Za-z0-9_]*(\\s{2,}[A-Za-z_][A-Za-z0-9_]*)+\\s*$", ln)) {
      col_names <- c(col_names, strsplit(trimws(ln), "\\s{2,}")[[1]])
      saw_header <- TRUE
    } else if (saw_header && !nzchar(trimws(ln))) {
      # blank line inside header — keep going
    } else if (saw_header) {
      break
    }
    i <- i + 1L
  }
  if (!length(col_names)) return(NULL)
  k <- length(col_names)

  body <- lines[i:length(lines)]
  stop_pat <- "Asymptotic Correlation Matrix|Final Results|Initial Summary|Estimates for Model"
  if (!grepl("Correlation", marker)) {
    # vcov block: stop at the correlation block.
  } else {
    stop_pat <- "Final Results|Initial Summary|Estimates for Model"
  }
  stop_idx <- grep(stop_pat, body)
  if (length(stop_idx)) body <- body[seq_len(stop_idx[1] - 1L)]

  M <- matrix(NA_real_, k, k)
  # Use unique row/col names so duplicate parameters (e.g. STATUS appearing
  # in both Early and Constant phases) get distinct slots.  Track phase
  # context in the body so callers can resolve "Early:STATUS" vs
  # "Constant:STATUS" via the rownames suffix.
  rownames(M) <- make.unique(col_names, sep = "__")
  colnames(M) <- rownames(M)

  # The body presents rows in the same order as header columns (verified
  # across all in-scope fixtures).  Track a single row cursor that
  # advances each time we see a new labelled row — this is robust to
  # duplicate parameter names because position, not name, drives the
  # routing.  Phase delimiters in the body act as sanity check only.
  row_cursor <- 0L
  pending <- list()
  current_row  <- NA_integer_
  current_name <- NA_character_
  current_phase <- NA_character_
  current_vals <- numeric()
  flush <- function() {
    if (!is.na(current_row) && length(current_vals)) {
      pending[[length(pending) + 1L]] <<- list(
        row = current_row, vals = current_vals,
        name = current_name, phase = current_phase
      )
    }
    current_row  <<- NA_integer_
    current_name <<- NA_character_
    current_vals <<- numeric()
  }

  for (ln in body) {
    # Skip rule lines.
    if (grepl("^\\s*-{10,}", ln)) next
    # Phase-label rows update phase context; don't reset row tracking.
    m_phase <- regmatches(ln, regexec("^\\s+(Early|Constant|Late):\\s*$", ln))[[1]]
    if (length(m_phase) >= 2) {
      current_phase <- m_phase[2]
      next
    }
    # Skip page-break title lines from form-feed inserts.
    if (!nzchar(trimws(ln))) next

    # New labelled row: leading identifier followed by numerics.
    toks <- regmatches(ln, regexec("^\\s+([A-Za-z_][A-Za-z0-9_]*)\\s+(.+)$", ln))[[1]]
    # Only treat as a real row if the identifier is one of the matrix's
    # parameter names.  Form-feed page breaks inject title-block lines
    # like "Primary Valve Operations ..." or "Multivariable Analysis ..."
    # that start with an identifier but aren't rows; they must not
    # advance the cursor or be parsed as data.
    is_real_row <- length(toks) >= 3 &&
      !grepl("^[+-]?[0-9.]", toks[2]) &&
      toks[2] %in% col_names
    if (is_real_row) {
      flush()
      rn <- toks[2]
      # Cursor-based routing handles duplicate parameter names (e.g.
      # STATUS or AGE_COP appearing in both Early and Constant phases).
      # Sanity-check that col_names[cursor] matches; otherwise pick the
      # next at-or-after-cursor match to keep position-based routing.
      row_cursor <- row_cursor + 1L
      if (row_cursor <= k && identical(col_names[row_cursor], rn)) {
        current_row <- row_cursor
      } else {
        candidates <- which(col_names == rn)
        forward <- candidates[candidates >= row_cursor]
        current_row <- if (length(forward)) forward[1] else
                       if (length(candidates)) candidates[1] else NA_integer_
        if (!is.na(current_row)) row_cursor <- current_row
      }
      current_name <- rn
      if (!is.na(current_row)) {
        rest <- toks[3]
        nums <- suppressWarnings(as.numeric(
          regmatches(rest, gregexpr("-?[0-9.]+(?:[Ee][+-]?[0-9]+)?", rest, perl = TRUE))[[1]]
        ))
        current_vals <- nums[!is.na(nums)]
      }
    } else if (!is.na(current_row) && length(toks) < 3) {
      # Continuation line: pure numerics, no leading identifier.  Title
      # lines (Primary, Multivariable, ...) have identifiers and we
      # explicitly skip them here by requiring the leading-identifier
      # regex to have *not* matched.
      nums <- suppressWarnings(as.numeric(
        regmatches(ln, gregexpr("-?[0-9.]+(?:[Ee][+-]?[0-9]+)?", ln, perl = TRUE))[[1]]
      ))
      current_vals <- c(current_vals, nums[!is.na(nums)])
    }
  }
  flush()

  for (p in pending) {
    n <- min(length(p$vals), k)
    if (n > 0L) M[p$row, seq_len(n)] <- p$vals[seq_len(n)]
  }

  for (i2 in seq_len(k)) for (j2 in seq_len(k)) {
    if (is.na(M[i2, j2]) && !is.na(M[j2, i2])) M[i2, j2] <- M[j2, i2]
  }
  M
}

# Split the file into per-fit blocks on "Initial Summary:" anchors.
.hzr_split_fits <- function(lines) {
  anchors <- grep("Initial Summary:", lines)
  if (!length(anchors)) return(list(lines))
  ends <- c(tail(anchors, -1L) - 1L, length(lines))
  mapply(function(s, e) lines[s:e], anchors, ends, SIMPLIFY = FALSE)
}

.hzr_parse_sas_lst <- function(path) {
  lines <- .hzr_read_lst(path)
  fits  <- .hzr_split_fits(lines)
  parsed <- lapply(fits, function(block) {
    obs <- .hzr_extract_obs_counts(block)
    ll  <- .hzr_extract_loglik(block)
    list(
      loglik            = if (length(ll)) ll[length(ll)] else NA_real_,
      loglik_trace      = ll,
      n_obs             = obs$n_obs,
      n_events          = obs$n_events,
      n_censored        = obs$n_censored,
      events_conserved  = .hzr_extract_events_conserved(block),
      params            = .hzr_extract_param_summary(block),
      natural           = .hzr_extract_natural(block),
      vcov              = .hzr_extract_matrix(block, "Asymptotic Variance-Covariance Matrix"),
      correlation       = .hzr_extract_matrix(block, "Asymptotic Correlation Matrix")
    )
  })
  list(fits = parsed)
}

# Parse the %DECILES goodness-of-fit table printed by the SAS deciles.hazard
# macro. The block is introduced by a spaced-out title "D E C I L E ..." and
# has one row per group: the leading token is the _DECILE_ index ("." for the
# overall row, then 0..groups-1), followed by numeric columns whose first four
# are CASES, EXPECTED (sum of cumulative hazard), CUM_HAZ (mean), ACTUAL
# (observed events). Returns a data frame with columns decile/cases/expected/
# actual (NA decile = overall row), or NULL if no decile block is present.
.hzr_parse_sas_deciles <- function(path) {
  lines <- .hzr_read_lst(path)
  start <- grep("D E C I L E   A N A L Y S I S", lines)
  if (!length(start)) return(NULL)
  body <- lines[start[1]:length(lines)]

  rows <- list()
  started <- FALSE
  num_only <- "^[0-9eE+.[:space:]-]+$"
  for (ln in body) {
    if (!nzchar(trimws(ln))) next
    m <- regmatches(ln, regexec("^[[:space:]]*([.]|[0-9]+)[[:space:]]+(.+)$", ln))[[1]]
    is_row <- length(m) >= 3 && grepl(num_only, m[3])
    if (is_row) {
      nums <- suppressWarnings(as.numeric(strsplit(trimws(m[3]),
                                                   "[[:space:]]+")[[1]]))
      nums <- nums[!is.na(nums)]
      if (length(nums) >= 4) {
        rows[[length(rows) + 1L]] <- data.frame(
          decile   = if (m[2] == ".") NA_integer_ else as.integer(m[2]),
          cases    = nums[1],
          expected = nums[2],
          actual   = nums[4]
        )
        started <- TRUE
        next
      }
    }
    # A non-blank, non-data line after the table has started ends the block.
    if (started) break
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Parse a SAS %KAPLAN / %NELSONT life-table block. These are printed by the
# kaplan.sas / nelsont.sas actuarial macros as a wide table with a header line
#   Obs INT_DEAD NUMBER CENSORED DEAD CUM_SURV SE_EXACT CL_LOWER ... PROPLIFE
# (KAPLAN) or a slightly different column set (NELSONT adds NUMRISK, drops
# SE_EXACT). Header-driven so it tolerates either layout. `which` selects the
# block: "kaplan" (first block with an SE_EXACT column), "nelson" (block with a
# NUMRISK column), or a `_CATG=<x>` stratum marker like "catg0"/"catg1".
# Returns a data frame named by the header columns (leading Obs dropped, "."
# -> NA), or NULL if the block is absent.
.hzr_parse_sas_lifetable <- function(path, which = c("kaplan", "nelson",
                                                     "catg0", "catg1")) {
  which <- match.arg(which)
  lines <- .hzr_read_lst(path)
  headers <- grep("^\\s*Obs\\s+INT_DEAD\\s+NUMBER", lines)
  if (!length(headers)) return(NULL)

  pick <- NULL
  if (which == "nelson") {
    pick <- headers[grepl("NUMRISK", lines[headers])][1]
  } else if (which == "kaplan") {
    # First KAPLAN-style header (has SE_EXACT, no NUMRISK).
    cand <- headers[grepl("SE_EXACT", lines[headers]) &
                      !grepl("NUMRISK", lines[headers])]
    pick <- cand[1]
  } else {
    # Stratum: the header following the matching _CATG=<n> rule line.
    n <- sub("catg", "", which)
    catg <- grep(paste0("_CATG=", n, "\\b"), lines)
    if (length(catg)) {
      after <- headers[headers > catg[1]]
      pick <- after[1]
    }
  }
  if (is.null(pick) || is.na(pick)) return(NULL)

  cols <- strsplit(trimws(lines[pick]), "\\s+")[[1]]
  rows <- list()
  for (ln in lines[(pick + 1L):length(lines)]) {
    if (grepl("^\\s*-{5,}", ln)) break          # summary rule ends the block
    if (!nzchar(trimws(ln))) next
    toks <- strsplit(trimws(ln), "\\s+")[[1]]
    # Data rows lead with the integer Obs counter.
    if (!grepl("^[0-9]+$", toks[1])) {
      if (length(rows)) break else next
    }
    if (length(toks) < length(cols)) next
    vals <- toks[seq_along(cols)]
    vals[vals == "."] <- NA
    rows[[length(rows) + 1L]] <- as.numeric(vals)
  }
  if (!length(rows)) return(NULL)
  df <- as.data.frame(do.call(rbind, rows))
  names(df) <- cols
  df[["Obs"]] <- NULL
  df
}

# Parse the HAZPRED "digital nomogram" prediction table printed by jobs like
# hp.death.AVC.sas: a PROC PRINT with header
#   Obs YEARS MONTHS _SURVIV _CLLSURV _CLUSURV _HAZARD _CLLHAZ _CLUHAZ
# Returns a data frame (Obs dropped) with the leading underscores stripped from
# the column names: YEARS, MONTHS, SURVIV, CLLSURV, CLUSURV, HAZARD, CLLHAZ,
# CLUHAZ. NULL if the table is absent.
.hzr_parse_sas_nomogram <- function(path) {
  lines <- .hzr_read_lst(path)
  h <- grep("YEARS[[:space:]]+MONTHS[[:space:]]+_SURVIV", lines)
  if (!length(h)) return(NULL)
  cols <- c("YEARS", "MONTHS", "SURVIV", "CLLSURV", "CLUSURV",
            "HAZARD", "CLLHAZ", "CLUHAZ")
  rows <- list()
  for (ln in lines[(h[1] + 1L):length(lines)]) {
    if (!nzchar(trimws(ln))) next
    toks <- strsplit(trimws(ln), "[[:space:]]+")[[1]]
    if (!grepl("^[0-9]+$", toks[1])) {
      if (length(rows)) break else next
    }
    if (length(toks) < 9L) next
    rows[[length(rows) + 1L]] <- as.numeric(toks[2:9])
  }
  if (!length(rows)) return(NULL)
  df <- as.data.frame(do.call(rbind, rows))
  names(df) <- cols
  df
}

# Parse the SAS %HAZBOOT bootstrap output printed by bs.death.AVC.sas: one
# per-phase coefficient table (Early/Constant/Late) where each row is one
# bootstrap resample's *selected* model -- a "." marks a covariate that was
# NOT selected in that resample.  Header looks like
#   Obs  E0  AGE STATUS COM_IV OPMOS OP_AGE RESAMPL ORIFICE MAL INC_SURG
# (intercept E0/C0/L0 identifies the phase: early/constant/late).
#
# Returns a named list with one data frame per present phase (covariate
# coefficients, NA where unselected; plus the RESAMPL id), and a
# `selection_freq` attribute: per-phase named integer vector counting, over
# the resamples, how often each covariate was selected.  NULL if no bootstrap
# table is present.  A phase whose covariates are all unselected (e.g. the
# placeholder Late table of a 2-phase fit) is dropped.
.hzr_parse_sas_bootstrap <- function(path) {
  lines <- .hzr_read_lst(path)
  phase_for <- c(E0 = "early", C0 = "constant", L0 = "late")

  parse_phase <- function(intercept) {
    h <- grep(paste0("Obs[[:space:]]+", intercept, "\\b.*RESAMPL"), lines)
    if (!length(h)) return(NULL)
    cols <- strsplit(trimws(lines[h[1]]), "[[:space:]]+")[[1]]
    rows <- list()
    for (ln in lines[(h[1] + 1L):length(lines)]) {
      if (!nzchar(trimws(ln))) next        # tolerate blank/pagination rows
      toks <- strsplit(trimws(ln), "[[:space:]]+")[[1]]
      if (!grepl("^[0-9]+$", toks[1])) {
        if (length(rows)) break else next
      }
      if (length(toks) < length(cols)) next
      vals <- toks[seq_along(cols)]
      vals[vals == "."] <- NA
      rows[[length(rows) + 1L]] <- suppressWarnings(as.numeric(vals))
    }
    if (!length(rows)) return(NULL)
    df <- as.data.frame(do.call(rbind, rows))
    names(df) <- cols
    df[["Obs"]] <- NULL
    df
  }

  out <- list()
  freq <- list()
  for (ic in names(phase_for)) {
    df <- parse_phase(ic)
    if (is.null(df)) next
    covs <- setdiff(names(df), c(ic, "RESAMPL"))
    f <- vapply(df[covs], function(x) sum(!is.na(x)), integer(1))
    if (sum(f) == 0L) next                     # drop empty placeholder phase
    out[[phase_for[[ic]]]] <- df
    freq[[phase_for[[ic]]]] <- f
  }
  if (!length(out)) return(NULL)
  attr(out, "selection_freq") <- freq
  out
}

# Parse the multivariable HAZPRED "digital nomogram" printed by the
# patient-specific prediction jobs (hp.death.AVC.hm1 / hm2).  Unlike the null
# model nomogram above, these PROC PRINT tables carry the full covariate
# profile on every row and are printed BY a grouping variable (MAL or COM_IV),
# so each row is a self-contained (time + covariates -> predictions) record:
#   Obs MONTHS YEARS OPYEAR OPMOS AGE COM_IV MAL INC_SURG ORIFICE STATUS
#       _SURVIV _CLLSURV _CLUSURV _HAZARD _CLLHAZ _CLUHAZ
# Header-driven (reads column names from the printed header), so it tolerates
# layout differences between fixtures.  Parsing is scoped to the data block
# immediately following each matching header line: within a block, data rows
# lead with the integer Obs counter; the block ends at the first non-data line
# after data have started (so unrelated `Obs ...` tables elsewhere in the
# listing are never ingested).  Page breaks / BY-groups repeat the header, and
# all their blocks are concatenated.  Returns a data frame with the leading "_"
# stripped from the prediction columns and Obs dropped; NULL if absent.
.hzr_parse_sas_nomogram_mv <- function(path) {
  lines <- .hzr_read_lst(path)
  h <- grep("MONTHS[[:space:]].*_SURVIV", lines)
  if (!length(h)) return(NULL)
  cols <- strsplit(trimws(lines[h[1]]), "[[:space:]]+")[[1]]
  ncol <- length(cols)
  rows <- list()
  for (hi in h) {
    if (hi >= length(lines)) next
    started <- FALSE
    for (ln in lines[(hi + 1L):length(lines)]) {
      if (!nzchar(trimws(ln))) next                 # tolerate blank spacer rows
      toks <- strsplit(trimws(ln), "[[:space:]]+")[[1]]
      if (!grepl("^[0-9]+$", toks[1])) {            # data rows lead with Obs
        if (started) break else next                # block ends after data
      }
      if (length(toks) < ncol) next
      vals <- toks[seq_len(ncol)]
      vals[vals == "."] <- NA
      rows[[length(rows) + 1L]] <- suppressWarnings(as.numeric(vals))
      started <- TRUE
    }
  }
  if (!length(rows)) return(NULL)
  df <- as.data.frame(do.call(rbind, rows))
  names(df) <- sub("^_", "", cols)
  df[["Obs"]] <- NULL
  df
}

# Fit the saved multivariable both-phase "HMDEATH" model used by the AVC
# patient-specific prediction fixtures (hm.death.AVC final model; consumed by
# hp.death.AVC.hm1 / hm2).  Early phase carries 6 covariates, Constant phase 3
# (STATUS appears in both).  Deterministic warm start from the .sas PARMS /
# coefficient statements (documented starting values, not the converged
# answer).  Returns the fitted hazard object.
.hzr_fit_avc_hmdeath <- function() {
  # Load via an explicit env + `$` so object_usage_linter sees a binding
  # (a bare `avc` from data() reads as an undefined global under lint_package).
  e <- new.env()
  utils::data("avc", package = "TemporalHazard", envir = e)
  d <- e$avc
  # SAS PROC STANDARD ... REPLACE fills missing INC_SURG with the column mean.
  d$inc_surg[is.na(d$inc_surg)] <- mean(d$inc_surg, na.rm = TRUE)
  start_thalf <- 0.1905077
  start_nu    <- 1.437416
  theta0 <- c(
    log(0.3504743), log(start_thalf), start_nu, 1,
    -0.03205774, 1.336675, 0.6872028, -0.01963377, 0.0002086689, 0.5169533,
    log(4.391673e-07), 1.375285, 3.11765, 1.054988
  )
  hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = start_thalf, nu = start_nu, m = 1,
                           fixed = "m",
                           formula = ~ age + com_iv + mal + opmos + op_age + status),
      constant = hzr_phase("constant", formula = ~ inc_surg + orifice + status)
    ),
    theta   = theta0,
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 2000, conserve = TRUE)
  )
}

# Parse the stratified "Predict number of deaths" observed-vs-expected table
# printed by hs.death.AVC.hm1.sas: one row per stratum (PROC SUMMARY BY COM_IV)
# with the printed header
#   Obs STATUS INC_SURG OPMOS AGE MAL COM_IV ORIFICE DEAD INT_DEAD PROB
#       CUM_HAZ DEAD PEXPECT EXPECTED ACTUAL
# (note the duplicated DEAD column, so we read by position).  Returns a data
# frame with one row per stratum: com_iv, pexpect (Sum 1 - S), expected
# (Sum cumulative hazard), actual (observed deaths).  NULL if absent.
.hzr_parse_sas_calibration <- function(path) {
  lines <- .hzr_read_lst(path)
  h <- grep("COM_IV.*PEXPECT.*EXPECTED.*ACTUAL", lines)
  if (!length(h)) return(NULL)
  rows <- list()
  for (ln in lines[(h[1] + 1L):length(lines)]) {
    if (!nzchar(trimws(ln))) next                  # SAS double-spaces the rows
    toks <- strsplit(trimws(ln), "[[:space:]]+")[[1]]
    if (!grepl("^[0-9]+$", toks[1])) {
      if (length(rows)) break else next
    }
    if (length(toks) < 16L) next
    v <- suppressWarnings(as.numeric(toks[1:16]))
    rows[[length(rows) + 1L]] <- data.frame(
      com_iv   = v[7],
      pexpect  = v[14],
      expected = v[15],
      actual   = v[16]
    )
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Parse the per-stratum mean-survival "digital" table printed by
# hs.death.AVC.hm1.sas: blocks delimited by "Interventricular communication=<g>"
# rule lines, each with header  Obs YEARS NSURVIV MSURVIV MCLLSURV MCLUSURV
# (NSURVIV = stratum size; MSURVIV/MCLLSURV/MCLUSURV = mean of per-subject
# survival and its CL across the stratum).  Returns a long data frame with a
# leading `com_iv` column; NULL if absent.
.hzr_parse_sas_strata_survival <- function(path) {
  lines <- .hzr_read_lst(path)
  rules <- grep("Interventricular communication=([0-9])", lines)
  if (!length(rules)) return(NULL)
  out <- list()
  for (r in rules) {
    g <- as.integer(sub(".*communication=([0-9]).*", "\\1", lines[r]))
    hdr <- grep("Obs[[:space:]]+YEARS[[:space:]]+NSURVIV", lines)
    hdr <- hdr[hdr > r][1]
    if (is.na(hdr)) next
    started <- FALSE
    for (ln in lines[(hdr + 1L):length(lines)]) {
      if (!nzchar(trimws(ln))) next                # tolerate blank spacer rows
      toks <- strsplit(trimws(ln), "[[:space:]]+")[[1]]
      if (!grepl("^[0-9]+$", toks[1])) {
        if (started) break else next
      }
      if (length(toks) < 6L) next
      started <- TRUE
      v <- suppressWarnings(as.numeric(toks[2:6]))
      out[[length(out) + 1L]] <- data.frame(
        com_iv = g, YEARS = v[1], NSURVIV = v[2],
        MSURVIV = v[3], MCLLSURV = v[4], MCLUSURV = v[5]
      )
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Default discovery: ~/Documents/GitHub/hazard/examples/ if it exists,
# else NULL.  Skip tests when fixtures are unavailable.
.hzr_sas_fixture_dir <- function() {
  env <- Sys.getenv("HAZARD_EXAMPLES_DIR", "")
  if (nzchar(env) && dir.exists(env)) return(env)
  default <- path.expand("~/Documents/GitHub/hazard/examples")
  if (dir.exists(default)) return(default)
  NA_character_
}

# ---------------------------------------------------------------------------
# OMC (PRIMISOL) dataset derivation
# ---------------------------------------------------------------------------
# The shipped `omc` R dataset has only the columns extracted at the source
# (study, te1, te2, te3, int_dead, dead, opdjul).  The OMC parity fixtures
# derive INT_TE / TE / STARTTME / CENSORED / NOPREVTE / MORBID via a
# multi-observation-per-patient DATA step that also needs te1djul/te2djul/
# te3djul, reop flags, fupdjul, and grade columns.  We reproduce that
# derivation here from the raw flat file.
#
# .hzr_omc_raw_path()     -- locate raw data file; returns NA_character_ if absent
# .hzr_derive_primisol()  -- parse + expand rows; returns list(te, te_mod, tm)

.hzr_omc_raw_path <- function() {
  env <- Sys.getenv("HAZARD_OMC_RAW", "")
  if (nzchar(env) && file.exists(env)) return(env)
  fixture_dir <- .hzr_sas_fixture_dir()
  if (!is.na(fixture_dir)) {
    cand <- file.path(fixture_dir, "data", "omc")
    if (file.exists(cand)) return(cand)
  }
  NA_character_
}

# Parse the SAS fixed-width raw file into a one-row-per-patient data frame.
# Column positions exactly match the SAS INPUT statement in hz.te123.OMC.sas.
.hzr_read_omc_raw <- function(path) {
  lines <- readLines(path, warn = FALSE)

  fld <- function(s, e = s) {
    v <- trimws(substr(lines, s, e))
    v[v == "."] <- NA_character_
    v
  }
  num <- function(s, e = s) suppressWarnings(as.numeric(fld(s, e)))
  chr <- function(s, e = s) fld(s, e)

  data.frame(
    study    = chr(1, 10),
    te1      = num(12),
    te2      = num(14),
    te3      = num(16),
    te1djul  = num(18, 25),
    te2djul  = num(27, 33),
    te3djul  = num(35, 41),
    avp_rp1  = chr(43),
    avp_rp2  = chr(45),
    avp_rp3  = chr(47),
    mvp_rp1  = chr(49),
    mvp_rp2  = chr(51),
    mvp_rp3  = chr(53),
    tvp_rp1  = chr(55),
    tvp_rp2  = chr(57),
    tvp_rp3  = chr(59),
    rp1djul  = num(61, 67),
    rp2djul  = num(69, 75),
    rp3djul  = num(77, 83),
    int_dead = num(111, 120),
    opdjul   = num(122, 128),
    dead     = num(130),
    fupdjul  = num(132, 138),
    te1grade = num(140),
    te2grade = num(142),
    te3grade = num(144),
    stringsAsFactors = FALSE
  )
}

# Reproduce the SAS PRIMISOL DATA step from hz.te123.OMC.sas and
# hz.tm123.OMC.sas.  Returns a list with three data frames:
#
#   te      -- for hz.te123.OMC fit 1 (left-truncated intervals, NOPREVTE)
#   te_mod  -- for hz.te123.OMC fit 2 (INT_TE adjusted to relative time,
#               + NOTE2 = NOPREVTE^2, NOTEE = exp(NOPREVTE))
#   tm      -- for hz.tm123.OMC    (same intervals, MORBID severity weight)
.hzr_derive_primisol <- function(raw_path) {
  d <- .hzr_read_omc_raw(raw_path)

  # IF TE1=1 AND TE1DJUL=. THEN DELETE
  d <- d[!(d$te1 == 1 & is.na(d$te1djul)), ]

  # Replace missing TE indicators with 0
  d$te1[is.na(d$te1)] <- 0L
  d$te2[is.na(d$te2)] <- 0L
  d$te3[is.na(d$te3)] <- 0L

  # Valve replacement: RP3 checked first, then RP2 overwrites, then RP1
  # overwrites (SAS IF statements execute in order 3->2->1; RP1 takes
  # precedence because it runs last and overwrites any earlier assignment).
  is_r <- function(x) !is.na(x) & trimws(x) == "R"
  n    <- nrow(d)
  REPL    <- rep(0L, n)
  REPLDJUL <- d$fupdjul   # default: censor at last follow-up

  for (rp_num in c(3L, 2L, 1L)) {
    avp <- is_r(d[[paste0("avp_rp", rp_num)]])
    mvp <- is_r(d[[paste0("mvp_rp", rp_num)]])
    tvp <- is_r(d[[paste0("tvp_rp", rp_num)]])
    idx <- avp | mvp | tvp
    REPL[idx]     <- 1L
    REPLDJUL[idx] <- d[[paste0("rp", rp_num, "djul")]][idx]
  }
  d$REPL    <- REPL
  d$REPLDJUL <- REPLDJUL

  # Censor TEs that occurred at or after valve replacement
  for (j in 1:3) {
    te_col   <- paste0("te",  j)
    djul_col <- paste0("te", j, "djul")
    mask <- d$REPL == 1L & d[[te_col]] == 1L &
            !is.na(d[[djul_col]]) & d[[djul_col]] >= d$REPLDJUL
    d[[te_col]][mask] <- 0L
  }

  # INT_REPL: time to valve replacement (or censoring if no replacement)
  d$INT_REPL <- d$int_dead
  repl_idx <- d$REPL == 1L
  d$INT_REPL[repl_idx] <-
    (d$REPLDJUL[repl_idx] - d$opdjul[repl_idx]) * 12 / 365.2425

  # -- Multi-row expansion -------------------------------------------------
  rows_te <- vector("list", n * 4L)   # pre-allocate generously
  rows_tm <- vector("list", n * 4L)
  nr_te <- 0L
  nr_tm <- 0L

  append_te <- function(study, te1, te2, te3, int_te, starttme, censored, te, noprevte) {
    nr_te <<- nr_te + 1L
    rows_te[[nr_te]] <<- data.frame(
      study = study, te1 = te1, te2 = te2, te3 = te3,
      int_te = int_te, starttme = starttme, censored = censored,
      te = te, noprevte = noprevte, stringsAsFactors = FALSE
    )
  }
  append_tm <- function(study, te1, te2, te3, int_te, starttme, censored, te, morbid) {
    nr_tm <<- nr_tm + 1L
    rows_tm[[nr_tm]] <<- data.frame(
      study = study, te1 = te1, te2 = te2, te3 = te3,
      int_te = int_te, starttme = starttme, censored = censored,
      te = te, morbid = morbid, stringsAsFactors = FALSE
    )
  }

  for (i in seq_len(nrow(d))) {
    r       <- d[i, ]
    study   <- trimws(r$study)
    te1_val <- r$te1
    te2_val <- r$te2
    te3_val <- r$te3

    # ---- Censored summary row (always output) ----------------------------
    # STARTTME advances through any events that preceded the censored interval.
    starttme <- 0
    noprevte <- 0L

    if (!is.na(te1_val) && te1_val == 1L) {
      starttme <- 12 * (r$te1djul - r$opdjul) / 365.2425
      noprevte <- 1L
    }
    # Special hospital-acquired TEs for two patients
    if (study == "039") starttme <- (3 / 24) * 12 / 365.2425
    if (study == "096") starttme <- (2 / 24) * 12 / 365.2425
    if (!is.na(te2_val) && te2_val == 1L) {
      starttme <- 12 * (r$te2djul - r$opdjul) / 365.2425
      noprevte <- 2L
    }
    if (!is.na(te3_val) && te3_val == 1L) {
      starttme <- 12 * (r$te3djul - r$opdjul) / 365.2425
      noprevte <- 3L
    }

    append_te(study, te1_val, te2_val, te3_val, r$INT_REPL, starttme, 1L, 0L, noprevte)
    append_tm(study, te1_val, te2_val, te3_val, r$INT_REPL, starttme, 1L, 0L, 0)

    # ---- TE1 event row ---------------------------------------------------
    if (!is.na(te1_val) && te1_val == 1L) {
      ite <- 12 * (r$te1djul - r$opdjul) / 365.2425
      if (study == "039") ite <- (3 / 24) * 12 / 365.2425
      if (study == "096") ite <- (2 / 24) * 12 / 365.2425
      append_te(study, te1_val, te2_val, te3_val, ite, 0, 0L, 1L, 0L)
      m1 <- if (is.na(r$te1grade)) 1.70 else r$te1grade
      append_tm(study, te1_val, te2_val, te3_val, ite, 0, 0L, 1L, m1)
    }

    # ---- TE2 event row ---------------------------------------------------
    if (!is.na(te2_val) && te2_val == 1L) {
      s2  <- 12 * (r$te1djul - r$opdjul) / 365.2425
      if (study == "039") s2 <- (3 / 24) * 12 / 365.2425
      if (study == "096") s2 <- (2 / 24) * 12 / 365.2425
      ite <- 12 * (r$te2djul - r$opdjul) / 365.2425
      append_te(study, te1_val, te2_val, te3_val, ite, s2, 0L, 1L, 1L)
      m2 <- if (is.na(r$te2grade)) 2.17 else r$te2grade
      append_tm(study, te1_val, te2_val, te3_val, ite, s2, 0L, 1L, m2)
    }

    # ---- TE3 event row ---------------------------------------------------
    if (!is.na(te3_val) && te3_val == 1L) {
      s2  <- 12 * (r$te2djul - r$opdjul) / 365.2425
      ite <- 12 * (r$te3djul - r$opdjul) / 365.2425
      append_te(study, te1_val, te2_val, te3_val, ite, s2, 0L, 1L, 2L)
      m3 <- if (is.na(r$te3grade)) 4.00 else r$te3grade
      append_tm(study, te1_val, te2_val, te3_val, ite, s2, 0L, 1L, m3)
    }
  }

  te <- do.call(rbind, rows_te[seq_len(nr_te)])
  tm <- do.call(rbind, rows_tm[seq_len(nr_tm)])

  # Drop zero/negative-duration rows (e.g. study 173: died same day as TE,
  # so the trailing censored row has STARTTME == INT_TE).  SAS HAZARD drops
  # these automatically when LCENSOR STARTTME is in effect.
  te <- te[te$int_te > te$starttme, ]
  tm <- tm[tm$int_te > tm$starttme, ]

  # hz.te123.OMC fit 2: modulated renewal — INT_TE relative to prior event
  te_mod <- te
  te_mod$int_te <- te_mod$int_te - te_mod$starttme
  te_mod$note2  <- te_mod$noprevte^2
  te_mod$notee  <- exp(te_mod$noprevte)
  # Drop any zero-duration modulated rows (relative interval would be 0)
  te_mod <- te_mod[te_mod$int_te > 0, ]

  list(te = te, te_mod = te_mod, tm = tm)
}
