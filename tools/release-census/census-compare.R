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
if (length(args) < 3L || length(args) > 4L) {
  stop("usage: census-compare.R <old.rds> <new.rds> <out.txt> [cases.R]",
       call. = FALSE)
}
old <- readRDS(args[[1]])
new <- readRDS(args[[2]])
out <- args[[3]]

# Which cases are DECLARED as the gate's known positive.
#
# Read from the cases file, not from the run records. The declaration is a
# property of the case definition, and a comparator that depends on the runner
# having copied the field forward has an extra failure mode -- which is not
# hypothetical: the first gated run reported "no case is declared as the gate's
# known positive" and exited 1 because census-run.R stored `kp` but not
# `gate_kp`. The gate was working; the assertion about the gate could not see
# its own declaration. Reading the source of truth removes that link entirely,
# and makes the comparison a pure function of the saved probes plus the
# declarations.
#
# Sourcing the cases file is safe without the package loaded: `census_cases()`
# only builds a list of closures, and the package is touched inside their
# bodies, which are not evaluated here.
declared_gate_kp <- character()
cases_file <- if (length(args) == 4L) {
  args[[4]]
} else {
  f <- file.path(dirname(sub("^--file=", "", grep("^--file=",
                 commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
                 "census-cases.R")
  if (length(f) == 1L && !is.na(f) && file.exists(f)) f else NA_character_
}
if (!is.na(cases_file) && file.exists(cases_file)) {
  local({
    e <- new.env(parent = globalenv())
    sys.source(cases_file, envir = e)
    cs <- e$census_cases()
    declared_gate_kp <<- names(cs)[vapply(cs, function(c) !is.null(c$gate_kp),
                                          logical(1))]
  })
}

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
                "hazard_examples_dir", "rng")
env_mismatch <- character()
for (f in env_fields) {
  a <- paste(as.character(old$env[[f]]), collapse = ",")
  b <- paste(as.character(new$env[[f]]), collapse = ",")
  flag <- if (identical(a, b)) "same" else "**DIFFERS**"
  say(sprintf("  %-20s %-40s %s", f, a, flag))
  if (!identical(a, b)) env_mismatch <- c(env_mismatch, f)
}
if (length(env_mismatch)) {
  say("")
  say("  !! ENVIRONMENT DIFFERS in: ", paste(env_mismatch, collapse = ", "))
  say("  !! A comparison across environments cannot tell code from")
  say("  !! environment, so it is not run. Re-run both sides alike.")
  quit(status = 2)
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

# `gate` carries optimum-quality evidence, not a result. It is dropped before
# comparing, so its presence cannot change any verdict -- and in particular
# cannot change the verdict of a case that does not use it.
GATE_KEY <- "gate"

compare_probe <- function(a, b) {
  # Returns a character vector of differing component names, or character(0).
  a[[GATE_KEY]] <- NULL
  b[[GATE_KEY]] <- NULL
  if (identical(a, b)) return(character())
  ka <- setdiff(union(names(a), names(b)), GATE_KEY)
  diffs <- character()
  for (k in ka) {
    if (!identical(a[[k]], b[[k]])) diffs <- c(diffs, k)
  }
  if (!length(diffs)) diffs <- "(identical components, differing structure)"
  diffs
}

# SAS/C accepts an optimum only when the relative gradient is at most
# eps^(1/3), about 6.055e-06. That is the criterion the package itself applies
# (R/optimizer.R), so it is derived the same way, not copied as a rounded
# literal that would accept gradients the package rejects.
GRADIENT_TOL <- .Machine$double.eps^(1 / 3)
RCOND_FLOOR <- 1e-8

gate_verdict <- function(g) {
  # Why this is in two parts: `fit$fit$rel_gradient` DOES NOT EXIST before
  # 1.2.11, so a gate that simply required the field would mark every row
  # uninterpretable on any older baseline regardless of how good the fit was.
  #
  #   part 1, version-independent: converged, every standard error finite,
  #           rcond above a floor. An NA standard error means the Hessian
  #           could not be inverted, so the point is not a proper interior
  #           maximum whatever `converged` reports.
  #   part 2, where the field exists: the relative-gradient test itself.
  #
  # In part 1, absent or NA is always "cannot confirm", never "passed"; so
  # is a gradient that was recorded as NA. A side that predates the field
  # is judged on part 1 alone, and the verdict line says so.
  if (is.null(g)) {
    return(list(ok = FALSE, why = "no gate data recorded"))
  }
  why <- character()
  if (!isTRUE(g$converged)) why <- c(why, "converged is not TRUE")
  if (!isTRUE(g$se_finite)) why <- c(why, "a standard error is NA (Hessian not invertible)")
  rc <- g$rcond
  if (!is.numeric(rc) || !is.finite(rc) || rc <= RCOND_FLOOR) {
    why <- c(why, paste0("rcond ", format(rc, digits = 3), " at or below ",
                         format(RCOND_FLOOR)))
  }
  if (isTRUE(g$rel_gradient_present)) {
    rg <- g$rel_gradient
    if (!is.numeric(rg) || !is.finite(rg) || rg > GRADIENT_TOL) {
      why <- c(why, paste0("rel_gradient ", format(rg, digits = 3),
                           " fails the SAS/C test (<= ",
                           format(GRADIENT_TOL), ")"))
    }
  }
  list(ok = !length(why), why = paste(why, collapse = "; "))
}

fmt_gate <- function(g) {
  if (is.null(g)) return("gate: none")
  paste0("converged ", isTRUE(g$converged),
         ", se finite ", isTRUE(g$se_finite),
         ", rcond ", format(g$rcond, digits = 3),
         ", rel_gradient ",
         if (isTRUE(g$rel_gradient_present)) format(g$rel_gradient, digits = 3)
         else "ABSENT (version predates the gradient test)")
}

max_disc <- function(a, b) {
  # Worst ABSOLUTE and worst RELATIVE discrepancy, and whether the relative
  # one is driven by values too small for a ratio to mean anything.
  #
  # A pure relative measure is actively misleading here, and the first run of
  # this census proved it: mp3_cdf reported "max rel discrepancy 1", which
  # looked like a parameter that had completely changed. It was `early.m`
  # moving from 4.53e-09 to 2.64e-11 -- both zero to any practical precision,
  # while every other parameter in that fit agreed to about 1e-9. Scaling by
  # pmax(|a|, |b|, 1e-300) turns two numerically-zero values into a ratio of
  # 0.994. This is the waldo absolute-branch trap the other way round: there,
  # an absolute tolerance hid a real difference; here, a relative one
  # manufactures a fake one.
  #
  # So both are reported, plus a floor flag, and neither is called "the"
  # discrepancy.
  a <- suppressWarnings(as.numeric(unlist(a)))
  b <- suppressWarnings(as.numeric(unlist(b)))
  none <- list(abs = NA_real_, rel = NA_real_, tiny = NA)
  if (!length(a) || length(a) != length(b)) return(none)
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(none)
  av <- a[ok]
  bv <- b[ok]
  ad <- abs(av - bv)
  rel <- ad / pmax(abs(av), abs(bv), 1e-300)
  i <- which.max(rel)
  list(
    abs = max(ad),
    rel = max(rel),
    # TRUE when the worst RELATIVE discrepancy comes from a pair whose
    # absolute difference is negligible -- i.e. the ratio is noise.
    tiny = ad[i] < 1e-8
  )
}

fmt_disc <- function(d) {
  if (is.na(d$abs)) return("discrepancy not measurable on these components")
  out <- sprintf("max abs %.4g, max rel %.4g", d$abs, d$rel)
  if (isTRUE(d$tiny)) {
    out <- paste0(out, " (the max-rel pair differs by ", format(d$abs, digits = 3),
                  " in absolute terms: the RATIO is noise, not a real move)")
  }
  out
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
  kp_gate <- if (!is.null(n$gate_kp)) n$gate_kp else o$gate_kp
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
  } else if (!is.null(o$probe[[GATE_KEY]]) || !is.null(n$probe[[GATE_KEY]])) {
    # A gated row. Both sides must reach a proper optimum or the row is NOT
    # classified: comparing two fits that are not at an optimum compares two
    # arbitrary stopping points, and "differs" is not a fact about the model.
    # This is the structural fix for the withdrawn mp_cov_formula finding.
    go <- gate_verdict(o$probe[[GATE_KEY]])
    gn <- gate_verdict(n$probe[[GATE_KEY]])
    gates <- paste0("old gate [", fmt_gate(o$probe[[GATE_KEY]]), "]",
                    " | new gate [", fmt_gate(n$probe[[GATE_KEY]]), "]")
    if (!go$ok || !gn$ok) {
      rows[[nm]] <- list(
        name = nm, outcome = "UNINTERPRETABLE",
        detail = paste0(
          "NOT CLASSIFIED -- ",
          if (!go$ok) paste0("old side: ", go$why, ". ") else "",
          if (!gn$ok) paste0("new side: ", gn$why, ". ") else "",
          "Comparing fits that are not at an optimum compares stopping ",
          "points, not models. ", gates),
        gate_kp = kp_gate, kp = kp)
    } else {
      d <- compare_probe(o$probe, n$probe)
      md <- max_disc(o$probe[d], n$probe[d])
      oo <- o$probe$objective
      no <- n$probe$objective
      obj <- if (is.numeric(oo) && is.numeric(no) && length(oo) == 1L &&
                 length(no) == 1L && identical(oo, no)) {
        "; objective identical"
      } else if (is.numeric(oo) && is.numeric(no) && length(oo) == 1L &&
                 length(no) == 1L) {
        sprintf("; objective moved by %.3g", no - oo)
      } else {
        ""
      }
      rows[[nm]] <- list(
        name = nm,
        outcome = if (!length(d)) "IDENTICAL-GATED" else "DIFFERS-GATED",
        detail = paste0(
          if (length(d)) paste0("components: ", paste(d, collapse = ", "),
                                "; ", fmt_disc(md), obj, ". ") else "",
          if (isTRUE(o$probe[[GATE_KEY]]$rel_gradient_present) &&
                isTRUE(n$probe[[GATE_KEY]]$rel_gradient_present)) {
            "BOTH SIDES AT A PROPER OPTIMUM, so this verdict is about the model. "
          } else {
            paste0("Both sides pass converged, finite SEs and rcond; the ",
                   "SAS/C gradient test ran only where the version records ",
                   "it (see the gates). ")
          },
          gates),
        gate_kp = kp_gate, kp = kp)
    }
  } else {
    d <- compare_probe(o$probe, n$probe)
    if (!length(d)) {
      rows[[nm]] <- list(name = nm, outcome = "IDENTICAL", detail = "",
                         kp = kp)
    } else {
      md <- max_disc(o$probe[d], n$probe[d])
      # The objective is the log-likelihood, so a fit whose objective ROSE
      # moved closer to the maximum. That distinction decides whether a
      # difference is an improvement or a regression, and it is the first
      # thing to look at before reading any discrepancy size.
      obj <- ""
      oo <- o$probe$objective
      no <- n$probe$objective
      if (is.numeric(oo) && is.numeric(no) && length(oo) == 1L &&
          length(no) == 1L && is.finite(oo) && is.finite(no)) {
        obj <- if (identical(oo, no)) {
          "; objective unchanged"
        } else if (no > oo) {
          sprintf("; objective IMPROVED by %.3g", no - oo)
        } else {
          sprintf("; objective WORSE by %.3g <-- look at this one", oo - no)
        }
      }
      rows[[nm]] <- list(name = nm, outcome = "DIFFERS",
                         detail = paste0("components: ", paste(d, collapse = ", "),
                                         "; ", fmt_disc(md), obj),
                         maxrel = md$rel, kp = kp)
    }
  }
  # Warning text is behaviour a user sees, tracked separately from values.
  ow <- if (is.null(o$warnings)) character() else o$warnings
  nw <- if (is.null(n$warnings)) character() else n$warnings
  rows[[nm]]$warn_changed <- !identical(sort(ow), sort(nw))
  rows[[nm]]$warn_detail <- if (rows[[nm]]$warn_changed) {
    paste0("old warnings: [", paste(ow, collapse = " | "),
           "] new warnings: [", paste(nw, collapse = " | "), "]")
  } else {
    ""
  }
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
# THE GATE'S OWN KNOWN POSITIVE
# ---------------------------------------------------------------------------
# A gate that refuses to classify uninterpretable rows is worth nothing unless
# some row actually trips it. Same discipline as the planted changes above.

gated <- Filter(function(r) {
  r$outcome %in% c("UNINTERPRETABLE", "DIFFERS-GATED", "IDENTICAL-GATED")
}, rows)
say("## Gradient gate on the appended multiphase cases")
if (!length(gated)) {
  say("  no gated cases in this run (an older cases file?)")
} else {
  for (r in gated) {
    say(sprintf("  %-28s %s", r$name, r$outcome))
  }
  say("")
  # Declared in the cases file (authoritative), with the run record as a
  # fallback for an rds written before census-run.R propagated the field.
  gate_kps <- Filter(function(r) {
    !is.null(r$gate_kp) || r$name %in% declared_gate_kp
  }, gated)
  if (!length(gate_kps)) {
    say("  FAIL: no case is declared as the gate's known positive, so the")
    say("  gate is untested. A gate no row trips cannot be trusted to stop")
    say("  the next uninterpretable row.")
    close(con)
    quit(status = 1L, save = "no")
  }
  bad <- Filter(function(r) r$outcome != "UNINTERPRETABLE", gate_kps)
  for (r in gate_kps) {
    say(sprintf("  known positive %-26s %-18s %s", r$name, r$outcome,
                if (r$outcome == "UNINTERPRETABLE") "CORRECTLY REFUSED"
                else "*** GATE FAILED TO REFUSE ***"))
  }
  if (length(bad)) {
    say("")
    say("  FAIL: a case built to be unidentifiable was classified anyway.")
    say("  The gate does not work, so every gated verdict here is suspect.")
    close(con)
    quit(status = 1L, save = "no")
  }
  say("  PASS: the gate refused the case built to be unidentifiable.")
}
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
