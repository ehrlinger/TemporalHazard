#!/usr/bin/env Rscript
# feature-matrix.R -- the data behind the feature test matrix article (#682).
#
# WHY
#   A hand-kept table of what is tested goes stale, and a stale "tested" is an
#   output that looks like a result and is not. So the counts are read from a
#   saved test run, and the run is checked against the tree before anything
#   is written. Only the status column and the notes are human judgement; they
#   live in tools/feature-matrix-map.csv beside the file patterns.
#
# WHAT IT REFUSES
#   - a results file that is not one full run of this tree: a test file on
#     disk with no rows in it, or a test name in a file that the run does not
#     carry, or the reverse;
#   - a run with any failure or error;
#   - a test file that no feature claims, and a feature pattern that matches
#     no file;
#   - a cited default-start test that is not in the run.
#
# USAGE (from the package root)
#   Rscript tools/feature-matrix.R <results.rds>
#   <results.rds> is as.data.frame(devtools::test()), as tools/release-suite.R
#   saves it. Writes vignettes/articles/feature-matrix.csv.

refuse <- function(...) {
  message("FEATURE MATRIX REFUSED: ", ...)
  quit(save = "no", status = 1)
}

if (!file.exists("DESCRIPTION") ||
      read.dcf("DESCRIPTION", "Package")[1, 1] != "TemporalHazard") {
  refuse("run this from the TemporalHazard package root.")
}
results_path <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(results_path) || !file.exists(results_path)) {
  refuse("give the saved test results (an .rds) as the first argument.")
}

res <- readRDS(results_path)
need <- c("file", "test", "passed", "failed", "skipped", "error")
if (!is.data.frame(res) || !all(need %in% names(res)) || nrow(res) == 0L) {
  refuse("the results file is not as.data.frame(devtools::test()).")
}
if (any(res$failed > 0L) || any(res$error)) {
  refuse(sum(res$failed > 0L | res$error), " test(s) failed or errored ",
         "in that run; a matrix over a red run would read as tested.")
}

# The run must be of this tree: every test file, and every test in it.
test_dir <- file.path("tests", "testthat")
files <- list.files(test_dir, pattern = "^test-.*\\.R$")
# A name built at run time (paste() in a loop) cannot be read here; such
# calls are counted, and the run may carry that many names unseen.
literal_tests <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  out <- character()
  computed <- 0L
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1L]], as.name("test_that"))) {
        if (is.character(e[[2L]])) {
          out <<- c(out, e[[2L]])
        } else {
          computed <<- computed + 1L
        }
      }
      for (a in as.list(e)[-1L]) if (!missing(a)) walk(a)
    }
  }
  for (e in exprs) walk(e)
  structure(out, computed = computed)
}
for (f in files) {
  on_disk <- literal_tests(file.path(test_dir, f))
  in_run <- unique(res$test[res$file == f])
  if (!length(in_run)) refuse(f, " has no rows in the run.")
  gone <- setdiff(on_disk, in_run)
  extra <- setdiff(in_run, on_disk)
  # A computed name in a loop can yield several tests, so any count of
  # unseen names is accepted from a file that computes some -- but only then.
  if (length(gone) || (length(extra) && !attr(on_disk, "computed"))) {
    refuse(f, " differs from the run (", length(gone), " test(s) not run, ",
           length(extra), " not on disk), so the run is of another tree. ",
           "First: ", c(gone, extra)[1L])
  }
}
stale <- setdiff(unique(res$file), files)
if (length(stale)) refuse("the run has files no longer on disk: ",
                          paste(stale, collapse = ", "))

map <- utils::read.csv(file.path("tools", "feature-matrix-map.csv"),
                       stringsAsFactors = FALSE, na.strings = character())
stem <- sub("\\.R$", "", sub("^test-", "", files))
claimed <- logical(length(files))
match_files <- function(pattern, what) {
  if (!nzchar(pattern)) return(character())
  hit <- grepl(pattern, stem, perl = TRUE)
  if (!any(hit)) refuse(what, " pattern '", pattern, "' matches no test file.")
  files[hit]
}

rows <- lapply(seq_len(nrow(map)), function(i) {
  m <- map[i, ]
  ref <- match_files(m$parity, paste(m$feature, "reference"))
  # A feature's reference comparisons are tests of it too.
  own <- union(match_files(m$files, m$feature), ref)
  claimed[files %in% own] <<- TRUE
  r <- res[res$file %in% own, ]
  ref_rows <- res[res$file %in% ref, ]
  ref_ran <- sum(!ref_rows$skipped)
  reference <- if (!length(ref)) {
    "none"
  } else {
    paste0(length(ref), " file(s); ", ref_ran, " of ", nrow(ref_rows),
           " tests ran")
  }
  default <- "none"
  if (nzchar(m$from_default)) {
    parts <- strsplit(m$from_default, "::", fixed = TRUE)[[1L]]
    hit <- res[res$file == parts[1L] & res$test == parts[2L], ]
    if (nrow(hit) != 1L) refuse("cited test not in the run: ", m$from_default)
    default <- if (hit$skipped) "cited test skipped" else "cited test passed"
  }
  data.frame(
    area = m$area, feature = m$feature, status = m$status,
    files = length(own), tests = nrow(r), expectations = sum(r$passed),
    skipped = sum(r$skipped), reference = reference,
    default_start = default, notes = m$notes,
    stringsAsFactors = FALSE
  )
})
if (!all(claimed)) {
  refuse("no feature claims ", paste(files[!claimed], collapse = ", "),
         ". Add it to a feature's pattern in tools/feature-matrix-map.csv.")
}
out <- do.call(rbind, rows)

sha <- system2("git", c("rev-parse", "--short=8", "HEAD"), stdout = TRUE)
dirty <- length(system2("git", c("status", "--porcelain", "--", "R", "tests"),
                        stdout = TRUE)) > 0L
target <- file.path("vignettes", "articles", "feature-matrix.csv")
dir.create(dirname(target), showWarnings = FALSE)
con <- file(target, "w")
writeLines(c(
  paste0("# measured_at: ", sha, if (dirty) " (uncommitted changes)"),
  paste0("# measured_on: ", format(file.mtime(results_path), "%Y-%m-%d")),
  paste0("# version: ", read.dcf("DESCRIPTION", "Version")[1, 1]),
  paste0("# suite: ", nrow(res), " tests, ", sum(res$passed),
         " expectations, ", sum(res$skipped), " skipped")
), con)
utils::write.csv(out, con, row.names = FALSE)
close(con)
message("Wrote ", target, ": ", nrow(out), " features at ", sha,
        if (dirty) " (uncommitted changes)")
