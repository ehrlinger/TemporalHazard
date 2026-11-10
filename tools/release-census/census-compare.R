## Behavioural census: compare two runs.
##
## Usage:
##   Rscript tools/release-census/census-compare.R <old.rds> <new.rds> <out.txt> \
##           [cases.R]
##
## Exit status:
##   0  the comparison ran, every planted perturbation was flagged and every
##      known positive was detected
##   1  a known positive was not detected, a planted perturbation was missed,
##      or the harness produced nothing
##   2  the two runs' environments differ
##   3  the two runs do not cover the same cases
##   4  the two runs are not provably different builds (same or missing sha)
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
declared_kp_expect <- list()
cases_file <- if (length(args) == 4L) {
  args[[4]]
} else {
  f <- file.path(dirname(sub("^--file=", "", grep("^--file=",
                 commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
                 "census-cases.R")
  if (length(f) == 1L && !is.na(f) && file.exists(f)) f else NA_character_
}
# The cases file is REQUIRED: coverage and the known positives are checked
# against the declared battery, not against whatever the two saved runs
# happen to hold. Two runs that are stale or incomplete in the same way agree
# with each other perfectly.
declared_names <- NULL
declared_kp <- character()
if (is.na(cases_file) || !file.exists(cases_file)) {
  cat("census-compare: no cases file; the declared battery cannot be checked.\n")
  quit(status = 3L, save = "no")
}
if (!is.na(cases_file) && file.exists(cases_file)) {
  local({
    e <- new.env(parent = globalenv())
    sys.source(cases_file, envir = e)
    cs <- e$census_cases()
    declared_names <<- names(cs)
    declared_kp <<- names(cs)[vapply(cs, function(c) !is.null(c[["kp"]]),
                                     logical(1))]
    declared_gate_kp <<- names(cs)[vapply(cs, function(c) !is.null(c[["gate_kp"]]),
                                          logical(1))]
    declared_kp_expect <<- Filter(Negate(is.null),
                                  lapply(cs, function(c) c[["kp_expect"]]))
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
# Build identity -- the two runs must be different commits
# ---------------------------------------------------------------------------
# Version cannot tell two refs apart: v1.3.0 and a 1.3.1 branch before its
# bump are both 1.3.0, and version-keyed paths once made both runs load the
# same library, comparing the new code with itself. Each install carries the
# commit it was built from (install-ref.sh stamps `CensusSha:`), and a pair
# whose stamps are equal or missing is not compared.
say("## Build identity")
o_sha <- old$env$build_sha
n_sha <- new$env$build_sha
say("  old build sha: ", if (is.null(o_sha)) "MISSING" else o_sha)
say("  new build sha: ", if (is.null(n_sha)) "MISSING" else n_sha)
if (is.null(o_sha) || is.null(n_sha) || is.na(o_sha) || is.na(n_sha) ||
      !nzchar(o_sha) || !nzchar(n_sha) || identical(o_sha, n_sha)) {
  say("  !! The two runs are not provably different builds, so a comparison")
  say("  !! could be the same code against itself. Not compared.")
  close(con)
  quit(status = 4L, save = "no")
}
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
say("  cases declared (cases file): ", length(declared_names),
    "   known positives declared: ", length(declared_kp))
# Each run must cover the DECLARED battery, not merely agree with the other.
miss_decl <- list(old = setdiff(declared_names, names(old$results)),
                  new = setdiff(declared_names, names(new$results)))
extra_decl <- list(old = setdiff(names(old$results), declared_names),
                   new = setdiff(names(new$results), declared_names))
for (side in c("old", "new")) {
  if (length(miss_decl[[side]])) say("  !! declared but not in the ", side,
                                     " run: ", paste(miss_decl[[side]], collapse = ", "))
  if (length(extra_decl[[side]])) say("  !! in the ", side, " run but not declared: ",
                                      paste(extra_decl[[side]], collapse = ", "))
}
# Both runs read the same cases file, so the sets must match. Comparing only
# the intersection would certify an incomplete census.
if (length(only_old) || length(only_new) ||
      length(unlist(miss_decl)) || length(unlist(extra_decl)) ||
      !identical(old$counts$declared, new$counts$declared) ||
      !identical(as.integer(old$counts$declared), length(declared_names))) {
  say("  !! The two runs do not cover the same cases, so the census is ")
  say("  !! incomplete and is not compared. Re-run both from one cases file.")
  quit(status = 3)
}
say("")

# ---------------------------------------------------------------------------
# Component comparison
# ---------------------------------------------------------------------------

# `gate` carries optimum-quality evidence, not a result. It is dropped before
# comparing, so its presence cannot change any verdict -- and in particular
# cannot change the verdict of a case that does not use it.
GATE_KEY <- "gate"
# `limit` carries the #683 surface (fit$fit$limit, the printed Limit check
# line, the bootstrap restart counts). It is reported in its own section and
# dropped here, so a fit that only gained the field still reads IDENTICAL on
# its numbers.
LIMIT_KEY <- "limit"
SIDE_KEYS <- c(GATE_KEY, LIMIT_KEY)

drop_side <- function(p) {
  for (k in SIDE_KEYS) p[[k]] <- NULL
  p
}

compare_probe <- function(a, b) {
  # Returns a character vector of differing component names, or character(0).
  a <- drop_side(a)
  b <- drop_side(b)
  if (identical(a, b)) return(character())
  ka <- union(names(a), names(b))
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
# The package's own conditioning threshold (.hzr_rcond_tol, R/hessian-invert.R),
# derived the same way: the rounded 1e-8 was looser and certified fits the
# package flags as ill-conditioned.
RCOND_FLOOR <- .Machine$double.eps^0.5

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

# A bootstrap with no successful replicate has nothing to compare: its
# replicate frame is empty and its spread NULL on both sides, so identical()
# says TRUE over nothing. Such a pair is NOT-COMPARED, never IDENTICAL.
boot_hollow <- function(p) {
  if (!identical(p$kind, "boot")) return(FALSE)
  !(is.numeric(p$n_success) && length(p$n_success) == 1L &&
      isTRUE(p$n_success > 0) && is.data.frame(p$replicates) &&
      nrow(p$replicates) > 0L)
}

rows <- list()
for (nm in all_names) {
  o <- old$results[[nm]]
  n <- new$results[[nm]]
  if (is.null(o) || is.null(n)) {
    rows[[nm]] <- list(name = nm, outcome = "MISSING-CASE",
                       detail = "case not declared in both runs",
                       kp = if (is.null(o)) n[["kp"]] else o[["kp"]])
    next
  }
  kp <- if (!is.null(n[["kp"]])) n[["kp"]] else o[["kp"]]
  kp_gate <- if (!is.null(n[["gate_kp"]])) n[["gate_kp"]] else o[["gate_kp"]]
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
  } else if (boot_hollow(o$probe) || boot_hollow(n$probe)) {
    # Both sides empty: nothing compared, never IDENTICAL. One side empty: the
    # replicates still cannot be compared, but 0 successes against some is
    # itself a change, and is labelled as one.
    both <- boot_hollow(o$probe) && boot_hollow(n$probe)
    rows[[nm]] <- list(
      name = nm, outcome = if (both) "NOT-COMPARED" else "NOT-COMPARED-ONE-SIDE",
      detail = paste0("bootstrap with no successful replicate on ",
                      if (both) "either side" else "one side",
                      ": n_success old ", format(o$probe$n_success), ", new ",
                      format(n$probe$n_success),
                      ". No replicate was compared, so this is not IDENTICAL."),
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
# PLANTED PERTURBATIONS -- the comparator must flag a change in every probe
# ---------------------------------------------------------------------------
# Independent of any release range. For every compared probe, a copy is made
# with ONE leaf nudged (a double by 4 ulps relative, or its smallest
# normal if zero; an integer by 1; a string by one trailing character; a
# logical flipped) and compare_probe() must flag the copy. A probe with no
# leaf to nudge, or a nudge the comparator misses, means an IDENTICAL verdict
# on that row would carry no information, so the script fails.
plant <- function(x) {
  # Within one component, a numeric leaf is preferred, so a data frame of
  # replicates is nudged in `estimate`, not in its `parameter` labels.
  if (is.list(x) && is_num_leaf(x)) {
    for (i in seq_along(x)) {
      if (!is_num_leaf(x[[i]])) next
      r <- plant(x[[i]])
      if (r$done) {
        x[[i]] <- r$value
        return(list(value = x, done = TRUE))
      }
    }
  }
  if (is.list(x)) {
    for (i in seq_along(x)) {
      r <- plant(x[[i]])
      if (r$done) {
        x[[i]] <- r$value
        return(list(value = x, done = TRUE))
      }
    }
    return(list(value = x, done = FALSE))
  }
  if (is.double(x) && length(x)) {
    j <- which(is.finite(x))[1]
    if (!is.na(j)) {
      x[j] <- if (x[j] == 0) .Machine$double.xmin else x[j] * (1 + 4 * .Machine$double.eps)
      return(list(value = x, done = TRUE))
    }
  }
  if (is.integer(x) && length(x) && !is.na(x[1])) {
    x[1] <- x[1] + 1L
    return(list(value = x, done = TRUE))
  }
  if (is.character(x) && length(x) && !is.na(x[1])) {
    x[1] <- paste0(x[1], "~")
    return(list(value = x, done = TRUE))
  }
  if (is.logical(x) && length(x) && !is.na(x[1])) {
    x[1] <- !x[1]
    return(list(value = x, done = TRUE))
  }
  list(value = x, done = FALSE)
}

say("## Planted perturbations (the comparator proving it sees a 4-ulp change)")
# One perturbation PER COMPARED COMPONENT, and compare_probe() must name that
# component. A single plant on the first leaf found only `kind`, a string
# tag, so the self-test proved nothing about theta, se or the objective.
# Numeric probes must have at least one numeric component perturbed; only a
# probe of kind "chr" may be all text.
is_num_leaf <- function(x) {
  if (is.list(x)) return(any(vapply(x, is_num_leaf, logical(1))))
  (is.double(x) || is.integer(x)) && any(!is.na(x))
}
compared <- c("IDENTICAL", "DIFFERS", "IDENTICAL-GATED", "DIFFERS-GATED")
plant_n <- 0L
plant_fail <- character()
for (r in Filter(function(r) r$outcome %in% compared, rows)) {
  for (side in c("old", "new")) {
    p <- drop_side(if (side == "old") old$results[[r$name]]$probe else new$results[[r$name]]$probe)
    n_num <- 0L
    for (k in setdiff(names(p), "kind")) {
      pl <- plant(p[[k]])
      if (!pl$done) next
      p2 <- p
      p2[[k]] <- pl$value
      plant_n <- plant_n + 1L
      if (is_num_leaf(p[[k]])) n_num <- n_num + 1L
      if (!k %in% compare_probe(p, p2)) {
        plant_fail <- c(plant_fail, paste0(r$name, " (", side, "): change to `",
                                           k, "` missed"))
      }
    }
    if (!identical(p$kind, "chr") && n_num == 0L) {
      plant_fail <- c(plant_fail, paste0(r$name, " (", side,
                                         "): no numeric component to perturb"))
    }
  }
}
say("  components perturbed: ", plant_n, "  missed or unperturbable: ",
    length(plant_fail))
if (!plant_n) {
  say("  FAIL: no compared probe to perturb. An empty comparison is not a pass.")
  close(con)
  quit(status = 1L, save = "no")
}
if (length(plant_fail)) {
  for (f in plant_fail) say("  !! ", f)
  say("  FAIL: the comparator cannot be shown to see a change on these rows.")
  close(con)
  quit(status = 1L, save = "no")
}
say("  PASS: every compared component flagged its planted perturbation.")
say("")

# What a known positive's NEW side must show (`kp_expect`, census-cases.R).
# Returns the unmet expectations; character(0) when all are met.
kp_unmet <- function(name, expect) {
  if (is.null(expect)) return(character())
  n <- new$results[[name]]
  lim <- n$probe[[LIMIT_KEY]]
  trig <- if (identical(n$status, "error")) n$trigger else lim$trigger
  miss <- character()
  if (!is.null(expect$limit_outcome) &&
        !identical(lim$outcome, expect$limit_outcome)) {
    miss <- c(miss, paste0("limit outcome '", format(lim$outcome),
                           "', expected '", expect$limit_outcome, "'"))
  }
  if (!is.null(expect$warn_class) &&
        !expect$warn_class %in% n$warn_classes) {
    miss <- c(miss, paste0("no warning of class ", expect$warn_class))
  }
  if (!is.null(expect$error_class) &&
        !identical(n$error_class, expect$error_class)) {
    miss <- c(miss, paste0("error class '", format(n$error_class),
                           "', expected ", expect$error_class))
  }
  if (!is.null(expect$trigger) &&
        !any(startsWith(as.character(trig), paste0(expect$trigger, " ")))) {
    miss <- c(miss, paste0("no trigger on ", expect$trigger))
  }
  if (!is.null(expect$n_restarted) &&
        !(is.numeric(lim$n_restarted) && isTRUE(lim$n_restarted >= expect$n_restarted))) {
    miss <- c(miss, paste0("n_restarted ", format(lim$n_restarted),
                           ", expected at least ", expect$n_restarted))
  }
  if (!is.null(expect$sas_rel_tol)) {
    so <- lim$sas_objective
    ob <- n$probe$objective
    if (!(is.numeric(so) && is.numeric(ob) &&
            isTRUE(abs(ob / so - 1) <= expect$sas_rel_tol))) {
      miss <- c(miss, paste0("objective ", format(ob, digits = 8),
                             " not within ", expect$sas_rel_tol,
                             " of SAS ", format(so)))
    }
  }
  miss
}

# ---------------------------------------------------------------------------
# KNOWN POSITIVES -- checked BEFORE any conclusion is drawn
# ---------------------------------------------------------------------------

say("## Known positives (the census proving it can detect a change)")
# Known positives come from the DECLARATION. A run record that lost or
# never carried its `kp` cannot drop a positive from this check.
kp_rows <- Filter(function(r) r$name %in% declared_kp, rows)
if (!setequal(vapply(kp_rows, function(r) r$name, ""), declared_kp)) {
  say("  FAIL: declared known positives without a compared row: ",
      paste(setdiff(declared_kp, vapply(kp_rows, function(r) r$name, "")),
            collapse = ", "))
  close(con)
  quit(status = 1L, save = "no")
}
if (!length(kp_rows)) {
  say("  FAIL: no known-positive case is declared. A census with no planted")
  say("  detection is a hollow result and this script will not certify it.")
  close(con)
  quit(status = 1L, save = "no")
}
kp_fail <- character()
detected <- c("DIFFERS", "NOW-ERRORS", "ERROR-OLD-ONLY", "ERROR-BOTH-DIFFERENT",
              "NOT-COMPARED-ONE-SIDE")
for (r in kp_rows) {
  expect <- declared_kp_expect[[r$name]]
  if (is.null(expect)) expect <- new$results[[r$name]][["kp_expect"]]
  unmet <- kp_unmet(r$name, expect)
  hit <- r$outcome %in% detected && !length(unmet)
  if (is.null(expect)) unmet <- "no kp_expect declared: a known positive must say what it expects"
  hit <- r$outcome %in% detected && !length(unmet)
  say(sprintf("  %-24s %-8s %-14s %s", r$name,
              if (hit) "DETECTED" else "NOT SEEN", r$outcome,
              format(r[["kp"]])))
  if (length(unmet)) say("      unmet on the new side: ", paste(unmet, collapse = "; "))
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
    !is.null(r[["gate_kp"]]) || r$name %in% declared_gate_kp
  }, rows)
  if (!length(gate_kps)) {
    say("  FAIL: no case is declared as the gate's known positive, so the")
    say("  gate is untested. A gate no row trips cannot be trusted to stop")
    say("  the next uninterpretable row.")
    close(con)
    quit(status = 1L, save = "no")
  }
  # The gate is judged on every side that returned a probe. A side that
  # errored (#683 refuses this case on the branch) leaves the gate untested
  # there, so at least one side must carry a probe, and every side that does
  # must be refused by gate_verdict().
  gate_refuses <- function(r) {
    sides <- Filter(function(rec) identical(rec$status, "ok"),
                    list(old$results[[r$name]], new$results[[r$name]]))
    length(sides) > 0L && all(vapply(sides, function(rec) {
      !gate_verdict(rec$probe[[GATE_KEY]])$ok
    }, logical(1)))
  }
  bad <- Filter(function(r) !gate_refuses(r), gate_kps)
  for (r in gate_kps) {
    say(sprintf("  known positive %-26s %-18s %s", r$name, r$outcome,
                if (gate_refuses(r)) {
                  paste0("CORRECTLY REFUSED on every side with a probe",
                         if (r$outcome == "NOW-ERRORS") " (new side errored)" else "")
                } else {
                  "*** GATE FAILED TO REFUSE ***"
                }))
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
# The #683 limit surface, reported apart from the bit-identity verdicts
# ---------------------------------------------------------------------------
# A field the old version does not have is reported as ABSENT on that side,
# never as a match: "absent -> NULL" is a change of surface, not agreement.
fmt_side <- function(rec) {
  if (is.null(rec)) return("no record")
  if (identical(rec$status, "error")) {
    return(paste0("ERROR ", format(rec$error_class),
                  if (length(rec$trigger)) paste0(" [", paste(rec$trigger, collapse = "; "), "]") else ""))
  }
  if (!identical(rec$status, "ok")) return(rec$status)
  lim <- rec$probe[[LIMIT_KEY]]
  if (is.null(lim)) return("no limit probe")
  if (identical(rec$probe$kind, "boot")) {
    return(paste0("n_restarted ", format(lim$n_restarted),
                  ", restart_reasons ",
                  if (is.character(lim$restart_reasons) &&
                        identical(lim$restart_reasons, "absent")) "absent"
                  else if (!length(lim$restart_reasons)) "none"
                  else paste(paste0(names(lim$restart_reasons), " x",
                                    lim$restart_reasons), collapse = "; "),
                  ", restarted column ", lim$restarted_column))
  }
  paste0("$limit ", lim$field,
         if (length(lim$trigger)) paste0(" [", paste(lim$trigger, collapse = "; "), "]") else "",
         " | print: ", lim$print_limit)
}
say("## #683 limit surface (reported, not part of the bit-identity verdicts)")
lim_rows <- 0L
for (nm in all_names) {
  o <- old$results[[nm]]
  n <- new$results[[nm]]
  has <- function(rec) {
    !is.null(rec) && (!is.null(rec$probe[[LIMIT_KEY]]) ||
                        !is.null(rec$error_class))
  }
  if (!has(o) && !has(n)) next
  ol <- o$probe[[LIMIT_KEY]]
  nl <- n$probe[[LIMIT_KEY]]
  # Error rows are only interesting when they carry a class the old side
  # did not raise (a new refusal); a plain error on both sides is skipped.
  if (is.null(ol) && is.null(nl) &&
        identical(o$error_class, n$error_class)) next
  lim_rows <- lim_rows + 1L
  o_abs <- is.null(ol) || identical(ol$field, "absent")
  n_abs <- is.null(nl) || identical(nl$field, "absent")
  state <- if (o_abs && n_abs) {
    "ABSENT-BOTH"
  } else if (o_abs) {
    "ABSENT-OLD"
  } else if (identical(fmt_side(o), fmt_side(n))) {
    "SAME"
  } else {
    "CHANGED"
  }
  say("- ", nm, "  [", state, "]")
  say("    old: ", fmt_side(o))
  say("    new: ", fmt_side(n))
  new_cls <- setdiff(n$warn_classes, o$warn_classes)
  if (length(new_cls)) say("    warning classes new on this side: ",
                           paste(unique(new_cls), collapse = ", "))
  if (!is.null(ol$print_other) && !is.null(nl$print_other)) {
    say("    print() apart from the limit line: ",
        if (identical(ol$print_other, nl$print_other)) "identical" else "CHANGED")
  }
  if (is.numeric(nl$sas_objective)) {
    say("    SAS optimum ", format(nl$sas_objective), "; new objective ",
        format(n$probe$objective, digits = 8), "; old objective ",
        format(o$probe$objective, digits = 8))
  }
}
if (!lim_rows) say("  none")
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
        if (!is.null(r[["kp"]])) paste0("  (planted: ", r[["kp"]], ")") else "")
    if (nzchar(r$detail)) say("    ", r$detail)
    if (isTRUE(r$warn_changed)) say("    WARNINGS CHANGED: ", r$warn_detail)
    say("    CLASSIFY: ")
  }
}
say("")

# Identical bootstraps are not clean bootstraps: two versions that both
# return the same repeated fit agree perfectly. Each side is read on its own.
say("## Bootstrap spread (read per side, whatever the verdict)")
spread_bad <- character()
for (nm in all_names) {
  for (side in c("old", "new")) {
    rec <- if (side == "old") old$results[[nm]] else new$results[[nm]]
    p <- rec$probe
    if (!identical(rec$status, "ok") || !identical(p$kind, "boot")) next
    ns <- p$n_success
    why <- if (!is.numeric(ns) || !isTRUE(ns >= 2)) {
      paste0("TOO FEW REPLICATES (n_success ", format(ns), "): spread not measurable")
    } else if (isTRUE(p$spread_is_zero)) {
      paste0("ZERO SPREAD on free parameter(s) ",
             paste(p$zero_spread_free, collapse = ", "),
             ": the replicates do not vary")
    } else if (is.null(p$zero_spread_free)) {
      "spread not recorded by this probe"
    }
    if (!is.null(why)) {
      spread_bad <- c(spread_bad, nm)
      say("  !! ", nm, " (", side, "): ", why)
    }
  }
}
spread_bad <- unique(spread_bad)
if (!length(spread_bad)) say("  every bootstrap varies on every free parameter, on both sides")
say("")

say("## Identical cases (no behavioural change across the release)")
idn_all <- names(Filter(function(r) r$outcome == "IDENTICAL", rows))
idn <- setdiff(idn_all, spread_bad)
if (length(intersect(spread_bad, idn_all))) {
  say("  (IDENTICAL but listed under Bootstrap spread, so not counted here: ",
      paste(intersect(spread_bad, idn_all), collapse = ", "), ")")
}
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
