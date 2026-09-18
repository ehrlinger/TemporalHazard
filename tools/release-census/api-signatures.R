## CRAN-facing upgrade review: signature and default-value diff.
##
## Usage:
##   Rscript tools/release-census/api-signatures.R <old-ref> <new-ref>
##
## Compares the FORMALS of every exported function between two refs, by
## parsing the R sources at each ref rather than installing them. A silently
## changed default is a breaking change that the version number does not
## announce, and it is invisible in an export-name diff -- which is why this
## exists alongside api-inventory.sh.
##
## Nothing is checked out: each ref is read through `git show`.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: api-signatures.R <old-ref> <new-ref>")
old_ref <- args[[1]]
new_ref <- args[[2]]

git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    stop("git ", paste(..., collapse = " "), " failed: ",
         paste(out, collapse = "\n"), call. = FALSE)
  }
  out
}

exports_at <- function(ref) {
  ns <- git("show", paste0(ref, ":NAMESPACE"))
  e <- grep("^export\\(", ns, value = TRUE)
  sort(sub("\\)$", "", sub("^export\\(", "", e)))
}

# Every top-level `name <- function(...)` in R/ at that ref, as a named list
# of formals. Parsed, not evaluated: evaluating the sources would need the
# package's dependencies and could run code.
formals_at <- function(ref) {
  files <- git("ls-tree", "--name-only", paste0(ref, ":R"))
  files <- grep("\\.R$", files, value = TRUE)
  out <- list()
  for (f in files) {
    txt <- git("show", paste0(ref, ":R/", f))
    exprs <- tryCatch(parse(text = txt), error = function(e) NULL)
    if (is.null(exprs)) {
      warning("could not parse R/", f, " at ", ref, call. = FALSE)
      next
    }
    for (ex in exprs) {
      if (!is.call(ex)) next
      op <- as.character(ex[[1]])
      if (length(op) != 1L || !op %in% c("<-", "=")) next
      lhs <- ex[[2]]
      rhs <- ex[[3]]
      if (!is.symbol(lhs)) next
      if (!is.call(rhs) || as.character(rhs[[1]])[1] != "function") next
      out[[as.character(lhs)]] <- rhs[[2]]
    }
  }
  out
}

fmt_formals <- function(fm) {
  if (is.null(fm)) return(NA_character_)
  if (!length(fm)) return("()")
  parts <- vapply(seq_along(fm), function(i) {
    nm <- names(fm)[i]
    d <- fm[[i]]
    if (missing(d) || identical(d, quote(expr = ))) return(nm)
    paste0(nm, " = ", paste(deparse(d), collapse = " "))
  }, character(1))
  paste0("(", paste(parts, collapse = ", "), ")")
}

old_exp <- exports_at(old_ref)
new_exp <- exports_at(new_ref)
old_fm <- formals_at(old_ref)
new_fm <- formals_at(new_ref)

cat("== exports:", old_ref, length(old_exp), "->", new_ref, length(new_exp), "==\n")
cat("== parsed top-level functions:", length(old_fm), "->", length(new_fm), "==\n\n")

# A known positive for the signature differ: a name whose formals we perturb
# must be reported. Without this, "no signature changed" is an untested claim.
probe <- intersect(old_exp, names(old_fm))[1]
if (is.na(probe) || is.null(probe)) stop("no exported function parsed; harness broken")
fake <- old_fm
fake[[probe]] <- as.pairlist(alist(zzz_sentinel = ))
if (identical(fmt_formals(fake[[probe]]), fmt_formals(old_fm[[probe]]))) {
  stop("KNOWN POSITIVE FAILED: the signature differ cannot see a change")
}
cat("known positive: signature differ detects a perturbed formals list on '",
    probe, "' -- PASS\n\n", sep = "")

cat("== CHANGED SIGNATURES on exports present in BOTH refs ==\n")
common <- intersect(old_exp, new_exp)
n_changed <- 0L
n_unparsed <- character()
for (nm in common) {
  a <- old_fm[[nm]]
  b <- new_fm[[nm]]
  if (is.null(a) || is.null(b)) {
    # Not a top-level `name <- function()`: an S3 method, or built another
    # way. Recorded, not silently skipped.
    n_unparsed <- c(n_unparsed, nm)
    next
  }
  sa <- fmt_formals(a)
  sb <- fmt_formals(b)
  if (!identical(sa, sb)) {
    n_changed <- n_changed + 1L
    cat("\n* ", nm, "\n", sep = "")
    cat("    ", old_ref, ": ", sa, "\n", sep = "")
    cat("    ", new_ref, ": ", sb, "\n", sep = "")
    # Name the argument-level delta, which is what a user actually hits.
    na <- names(a)
    nb <- names(b)
    added <- setdiff(nb, na)
    removed <- setdiff(na, nb)
    if (length(added)) cat("    ARGS ADDED  : ", paste(added, collapse = ", "), "\n", sep = "")
    if (length(removed)) cat("    ARGS REMOVED: ", paste(removed, collapse = ", "), "\n", sep = "")
    for (k in intersect(na, nb)) {
      da <- if (identical(a[[k]], quote(expr = ))) "<no default>" else paste(deparse(a[[k]]), collapse = " ")
      db <- if (identical(b[[k]], quote(expr = ))) "<no default>" else paste(deparse(b[[k]]), collapse = " ")
      if (!identical(da, db)) {
        cat("    DEFAULT CHANGED: ", k, ": ", da, " -> ", db, "\n", sep = "")
      }
    }
    # Argument ORDER matters for positional calls.
    if (!identical(intersect(na, nb), intersect(nb, na))) {
      cat("    ARG ORDER CHANGED among the shared arguments\n")
    }
  }
}
if (n_changed == 0L) cat("  (none)\n")
cat("\n  exports compared: ", length(common),
    "; signatures changed: ", n_changed,
    "; not parsed as top-level functions: ", length(n_unparsed), "\n", sep = "")
if (length(n_unparsed)) {
  cat("  not parsed (checked by hand, not dropped): ",
      paste(n_unparsed, collapse = ", "), "\n", sep = "")
}

cat("\n== NEW exports and their signatures ==\n")
for (nm in setdiff(new_exp, old_exp)) {
  cat("+ ", nm, " ", fmt_formals(new_fm[[nm]]), "\n", sep = "")
}
cat("\n== REMOVED exports ==\n")
rm_exp <- setdiff(old_exp, new_exp)
if (!length(rm_exp)) cat("  (none)\n") else for (nm in rm_exp) cat("- ", nm, "\n", sep = "")
