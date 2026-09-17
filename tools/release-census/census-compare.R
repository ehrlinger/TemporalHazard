## Behavioural census: compare two runs.
##
## Usage:
##   Rscript tools/release-census/census-compare.R <old.rds> <new.rds> <out.txt>
##
## Exit status:
##   0  the comparison ran AND every known positive was detected
##   1  a known positive compared identical, or the harness produced nothing
##
## A non-zero exit is the gate failing. "All identical" with no demonstrated
## detection is not a pass -- it is an untested comparator, so the known
## positives below are checked first and the script refuses to report
## anything if they did not fire.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: census-compare.R <old.rds> <new.rds> <out.txt>", call. = FALSE)
}
old <- readRDS(args[[1]])
new <- readRDS(args[[2]])
out <- args[[3]]

con <- file(out, open = "wt")
say <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = con)
}

say("# Behavioural census: ", old$env$version, " -> ", new$env$version)
say("")

# ---------------------------------------------------------------------------
# Environment reconciliation -- a difference here invalidates comparisons
# ---------------------------------------------------------------------------

say("## Environment")
env_fields <- c("r_version", "platform", "blas", "lapack", "survival",
                "numDeriv", "qhsstudies_mounted", "hazard_repo",
                "hazard_examples_dir")
env_mismatch <- character()
for (f in env_fields) {
  a <- paste(as.character(old$env[[f]]), collapse = ",")
  b <- paste(as.character(new$env[[f]]), collapse = ",")
  flag <- if (identical(a, b)) "same" else "**DIFFERS**"
  say(sprintf("  %-20s %-40s %s", f, a, flag))
  if (!identical(a, b)) env_mismatch <- c(env_mismatch, f)
}
say(sprintf("  %-20s %s", "RNGkind",
            paste(old$env$rng, collapse = ",")))
if (length(env_mismatch)) {
  say("")
  say("  !! ENVIRONMENT DIFFERS in: ", paste(env_mismatch, collapse = ", "))
  say("  !! Numeric differences below may be environmental, not code.")
}
say("")
say("## SAS fixtures")
say("  /Volumes/qhsstudies mounted: ", old$env$qhsstudies_mounted,
    " (old run) / ", new$env$qhsstudies_mounted, " (new run)")
if (!isTRUE(old$env$qhsstudies_mounted)) {
  say("  NOT MOUNTED: no SAS-fixture case ran. The census covers the R")
  say("  behavioural surface only; SAS parity is NOT covered by this run.")
}
say("")

# ---------------------------------------------------------------------------
# Coverage, asserted before any value is compared
# ---------------------------------------------------------------------------

all_names <- union(names(old$results), names(new$results))
say("## Coverage")
say("  cases declared (old run): ", old$counts$declared,
    "   ok ", old$counts$ok, "  error ", old$counts$error,
    "  absent ", old$counts$absent)
say("  cases declared (new run): ", new$counts$declared,
    "   ok ", new$counts$ok, "  error ", new$counts$error,
    "  absent ", new$counts$absent)
say("  union of case names     : ", length(all_names))

only_old <- setdiff(names(old$results), names(new$results))
only_new <- setdiff(names(new$results), names(old$results))
if (length(only_old)) say("  !! cases present only in the old run: ",
                          paste(only_old, collapse = ", "))
if (length(only_new)) say("  !! cases present only in the new run: ",
                          paste(only_new, collapse = ", "))
say("")

# ---------------------------------------------------------------------------
# Component comparison
# ---------------------------------------------------------------------------

compare_probe <- function(a, b) {
  # Returns a character vector of differing component names, or character(0).
  if (identical(a, b)) return(character())
  ka <- union(names(a), names(b))
  diffs <- character()
  for (k in ka) {
    if (!identical(a[[k]], b[[k]])) diffs <- c(diffs, k)
  }
  if (!length(diffs)) diffs <- "(identical components, differing structure)"
  diffs
}

max_rel <- function(a, b) {
  # Worst relative discrepancy, for reporting magnitude only.
  a <- suppressWarnings(as.numeric(unlist(a)))
  b <- suppressWarnings(as.numeric(unlist(b)))
  if (!length(a) || length(a) != length(b)) return(NA_real_)
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  scale <- pmax(abs(a[ok]), abs(b[ok]), 1e-300)
  max(abs(a[ok] - b[ok]) / scale)
}

rows <- list()
for (nm in all_names) {
  o <- old$results[[nm]]
  n <- new$results[[nm]]
  if (is.null(o) || is.null(n)) {
    rows[[nm]] <- list(name = nm, outcome = "MISSING-CASE",
                       detail = "case not declared in both runs",
                       kp = if (is.null(o)) n$kp else o$kp)
    next
  }
  kp <- if (!is.null(n$kp)) n$kp else o$kp
  os <- o$status
  ns <- n$status

  if (os == "absent" || ns == "absent") {
    rows[[nm]] <- list(name = nm, outcome = "ABSENT",
                       detail = paste0("old=", os, " new=", ns, "; missing: ",
                                       paste(c(o$missing, n$missing),
                                             collapse = ",")),
                       kp = kp)
  } else if (os == "error" && ns == "error") {
    same <- identical(o$message, n$message)
    rows[[nm]] <- list(name = nm,
                       outcome = if (same) "ERROR-BOTH-SAME" else "ERROR-BOTH-DIFFERENT",
                       detail = paste0("old: ", o$message, " || new: ", n$message),
                       kp = kp)
  } else if (os == "error" && ns == "ok") {
    rows[[nm]] <- list(name = nm, outcome = "ERROR-OLD-ONLY",
                       detail = paste0("old errored: ", o$message),
                       kp = kp)
  } else if (os == "ok" && ns == "error") {
    rows[[nm]] <- list(name = nm, outcome = "NOW-ERRORS",
                       detail = paste0("new errors: ", n$message),
                       kp = kp)
  } else {
    d <- compare_probe(o$probe, n$probe)
    if (!length(d)) {
      rows[[nm]] <- list(name = nm, outcome = "IDENTICAL", detail = "",
                         kp = kp)
    } else {
      mr <- max_rel(o$probe[d], n$probe[d])
      rows[[nm]] <- list(name = nm, outcome = "DIFFERS",
                         detail = paste0("components: ", paste(d, collapse = ", "),
                                         "; max rel discrepancy ",
                                         format(mr, digits = 4)),
                         maxrel = mr, kp = kp)
    }
  }
  # Warning text is behaviour a user sees, tracked separately from values.
  ow <- if (is.null(o$warnings)) character() else o$warnings
  nw <- if (is.null(n$warnings)) character() else n$warnings
  rows[[nm]]$warn_changed <- !identical(sort(ow), sort(nw))
  rows[[nm]]$warn_detail <- if (rows[[nm]]$warn_changed) {
    paste0("old warnings: [", paste(ow, collapse = " | "),
           "] new warnings: [", paste(nw, collapse = " | "), "]")
  } else ""
}

# ---------------------------------------------------------------------------
# KNOWN POSITIVES -- checked BEFORE any conclusion is drawn
# ---------------------------------------------------------------------------

say("## Known positives (the census proving it can detect a change)")
kp_rows <- Filter(function(r) !is.null(r$kp), rows)
if (!length(kp_rows)) {
  say("  FAIL: no known-positive case is declared. A census with no planted")
  say("  detection is a hollow result and this script will not certify it.")
  close(con)
  quit(status = 1L, save = "no")
}
kp_fail <- character()
detected <- c("DIFFERS", "NOW-ERRORS", "ERROR-OLD-ONLY", "ERROR-BOTH-DIFFERENT")
for (r in kp_rows) {
  hit <- r$outcome %in% detected
  say(sprintf("  %-28s %-22s %-8s %s", r$name, r$kp,
              if (hit) "DETECTED" else "NOT SEEN", r$outcome))
  if (!hit) kp_fail <- c(kp_fail, r$name)
}
say("")
if (length(kp_fail)) {
  say("  FAIL: these deliberate changes were NOT detected: ",
      paste(kp_fail, collapse = ", "))
  say("  The comparator is not proven to work, so its 'identical' verdicts")
  say("  carry no information. Fix the harness before reading anything below.")
  close(con)
  quit(status = 1L, save = "no")
}
say("  PASS: every planted change was detected. Verdicts below are meaningful.")
say("")

# ---------------------------------------------------------------------------
# Outcome tally and the difference table
# ---------------------------------------------------------------------------

outcomes <- vapply(rows, function(r) r$outcome, character(1))
say("## Outcome tally")
tt <- sort(table(outcomes), decreasing = TRUE)
for (i in seq_along(tt)) say(sprintf("  %-24s %d", names(tt)[i], tt[[i]]))
say(sprintf("  %-24s %d", "TOTAL", length(rows)))
say("")

say("## Every non-identical case (classify each as (a) NEWS, (b) issue, (c) UNEXPLAINED)")
say("")
nd <- Filter(function(r) r$outcome != "IDENTICAL", rows)
if (!length(nd)) {
  say("  none")
} else {
  for (r in nd) {
    say("- ", r$name, "  [", r$outcome, "]",
        if (!is.null(r$kp)) paste0("  (planted: ", r$kp, ")") else "")
    if (nzchar(r$detail)) say("    ", r$detail)
    if (isTRUE(r$warn_changed)) say("    WARNINGS CHANGED: ", r$warn_detail)
    say("    CLASSIFY: ")
  }
}
say("")

say("## Identical cases (no behavioural change across the release)")
idn <- names(Filter(function(r) r$outcome == "IDENTICAL", rows))
say("  ", length(idn), " of ", length(rows), ": ",
    paste(idn, collapse = ", "))
say("")

say("## Cases whose warning text changed but whose values did not")
wc <- Filter(function(r) isTRUE(r$warn_changed) && r$outcome == "IDENTICAL", rows)
if (!length(wc)) {
  say("  none")
} else {
  for (r in wc) say("- ", r$name, ": ", r$warn_detail)
}

close(con)
cat("\nwrote:", out, "\n")
