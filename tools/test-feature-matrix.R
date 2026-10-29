#!/usr/bin/env Rscript
# test-feature-matrix.R -- each refusal in tools/feature-matrix.R, planted.
#
# The generator is a tool, not package code, so it is tested here rather than
# in tests/testthat (tools/ is not in the built package). Each case starts
# from a synthetic run that the generator accepts, plants one fault in a
# scratch clone, and requires the refusal it names. The accepted run is the
# known positive: without it a refusal could come from a broken fixture.
#
# USAGE (from the package root)
#   Rscript tools/test-feature-matrix.R
# Prints one line per case and "ALL n CASES PASSED"; exits 1 on any failure.

root <- normalizePath(".")
scratch <- tempfile("feature-matrix-test-")
# system2() pastes its arguments into a shell command unquoted, so quote each.
git <- function(...) {
  out <- suppressWarnings(system2("git", shQuote(c("-C", scratch, ...)),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(out, "status"))) {
    stop("git ", paste(c(...), collapse = " "), ": ", paste(out, collapse = " "))
  }
  invisible(out)
}
commit_all <- function(msg) {
  git("add", "-A")
  git("-c", "user.name=test", "-c", "user.email=test@example.invalid",
      "commit", "-q", "--allow-empty", "-m", msg)
}

system2("git", shQuote(c("clone", "-q", "--no-hardlinks", root, scratch)))
# Bring this tree's tools/ in, so uncommitted generator changes are tested.
invisible(file.copy(list.files(file.path(root, "tools"), full.names = TRUE),
                    file.path(scratch, "tools"), overwrite = TRUE))
commit_all("test base")
setwd(scratch)
source(file.path("tools", "feature-matrix.R"))

# A run of every test on disk, all passing, with this tree's provenance.
synthetic_run <- function() {
  computed <- utils::read.csv(file.path("tools", "feature-matrix-computed.csv"),
                              stringsAsFactors = FALSE)
  files <- list.files(file.path("tests", "testthat"),
                      pattern = "^test-.*\\.R$")
  rows <- lapply(files, function(f) {
    lit <- matrix_literal_tests(file.path("tests", "testthat", f))
    t <- c(lit$names, computed$test[computed$file == f])
    data.frame(file = rep(f, length(t)), test = t, passed = 1L,
               failed = 0L, skipped = FALSE, error = FALSE,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  attr(res, "provenance") <- suite_provenance(Sys.time())
  res
}

outcome <- function(expr) {
  tryCatch({
    force(expr)
    "accepted"
  }, feature_matrix_refused = function(e) conditionMessage(e))
}

failures <- 0L
n_cases <- 0L
check <- function(name, got, pattern) {
  n_cases <<- n_cases + 1L
  ok <- if (is.na(pattern)) identical(got, "accepted") else grepl(pattern, got,
                                                                   fixed = TRUE)
  if (!ok) failures <<- failures + 1L
  cat(if (ok) "ok   " else "FAIL ", name, if (!ok) paste0(": got ", got), "\n")
}

base <- synthetic_run()
check("a complete run of this tree is accepted", outcome(feature_matrix(base)),
      NA)
out <- feature_matrix(base)

# Finding 4: the tests column counts only what ran.
sk <- base
sk$skipped[sk$file == "test-repeating-events.R"][1L] <- TRUE
row <- feature_matrix(sk)
row <- row[row$feature == "Time-varying covariates (start/stop epochs)", ]
was <- out[out$feature == "Time-varying covariates (start/stop epochs)", ]
n_cases <- n_cases + 1L
if (identical(row$tests_run, was$tests_run - 1L) && identical(row$skipped, 1L)) {
  cat("ok    a skipped test is not counted as run\n")
} else {
  failures <- failures + 1L
  cat("FAIL  a skipped test is not counted as run\n")
}

# Finding 1: helper tests passing beside skipped comparisons count zero.
refs <- utils::read.csv(file.path("tools", "feature-matrix-references.csv"),
                        stringsAsFactors = FALSE)
w <- refs[refs$feature == "Observation weights", ]
hp <- base
hp$skipped[paste(hp$file, hp$test) %in% paste(w$file, w$test)] <- TRUE
helpers <- hp$file == "test-weights-sas-parity.R" & !hp$skipped
row <- feature_matrix(hp)
row <- row[row$feature == "Observation weights", ]
n_cases <- n_cases + 1L
if (any(helpers) && startsWith(row$reference, paste0("0 of ", nrow(w)))) {
  cat("ok    validator tests beside skipped comparisons count 0 references\n")
} else {
  failures <- failures + 1L
  cat("FAIL  validator tests beside skipped comparisons count 0 references:",
      row$reference, "\n")
}

red <- base
red$failed[1L] <- 1L
check("a run with a failure is refused", outcome(feature_matrix(red)),
      "failed or errored")

bare <- base
attr(bare, "provenance") <- NULL
check("a run without provenance is refused", outcome(feature_matrix(bare)),
      "carries no provenance")

dirty_run <- base
attr(dirty_run, "provenance")$dirty <- " M R/hazard_api.R"
check("a run made on a dirty tree is refused",
      outcome(feature_matrix(dirty_run)), "uncommitted changes")

dep <- base
attr(dep, "provenance")$dependencies[["survival"]] <- "0.0.1"
check("a run under other dependency versions is refused",
      outcome(feature_matrix(dep)), "dependencies differ")

# Finding 3: one generated test missing, or one invented, is refused.
computed <- utils::read.csv(file.path("tools", "feature-matrix-computed.csv"),
                            stringsAsFactors = FALSE)
gen <- computed[computed$file == "test-phase-constraint-oracle.R", ]
drop <- base[!(base$file == gen$file[1L] & base$test == gen$test[1L]), ]
check("a generated test missing from the run is refused",
      outcome(feature_matrix(drop)), "expected test(s) not run")
# Added, not renamed: a rename also removes an expected name, and the missing
# check would refuse it without the extra check ever being needed.
inv <- base[base$file == gen$file[1L], ][1L, ]
inv$test <- "invented"
inv <- rbind(base, inv)
attr(inv, "provenance") <- attr(base, "provenance")
check("an unexpected generated name is refused",
      outcome(feature_matrix(inv)),
      "0 expected test(s) not run, 1 run but not expected")

stale <- base
stale$test[stale$file == "test-wald.R"][1L] <- "renamed"
check("a renamed test is refused", outcome(feature_matrix(stale)),
      "differs from the run")

# Finding 2: a run of an earlier commit, after R/ changed, is refused.
cat("# a later change\n", file = file.path("R", "zzz-feature-matrix-test.R"))
check("an uncommitted change to R/ now is refused",
      outcome(feature_matrix(base)), "the tree has uncommitted changes")
commit_all("change R/")
check("a run of an earlier commit after R/ changed is refused",
      outcome(feature_matrix(base)), "file(s) a test reads have changed")
unlink(file.path("R", "zzz-feature-matrix-test.R"))
commit_all("revert R/")
cat("note\n", file = "NOTES-feature-matrix-test.md")
commit_all("docs only")
check("a run of an earlier commit after a docs-only change is accepted",
      outcome(feature_matrix(base)), NA)

writeLines('test_that("x", expect_true(TRUE))',
           file.path("tests", "testthat", "test-zzz-unclaimed.R"))
commit_all("new test file")
fresh <- synthetic_run()
check("a test file no feature claims is refused",
      outcome(feature_matrix(fresh)), "no feature claims test-zzz-unclaimed.R")
unlink(file.path("tests", "testthat", "test-zzz-unclaimed.R"))
commit_all("drop test file")
base <- synthetic_run()

map_path <- file.path("tools", "feature-matrix-map.csv")
map_text <- readLines(map_path)
writeLines(sub("^repeating-events$", "^repeating-eventsx$", map_text,
               fixed = TRUE), map_path)
check("a pattern matching no file is refused", outcome(feature_matrix(base)),
      "matches no test file")
writeLines(map_text, map_path)

ref_path <- file.path("tools", "feature-matrix-references.csv")
ref_text <- readLines(ref_path)
writeLines(sub("matches C output", "matches C outputx", ref_text, fixed = TRUE),
           ref_path)
check("a cited reference test not in the run is refused",
      outcome(feature_matrix(base)), "cited reference test not in the run")
writeLines(ref_text, ref_path)

writeLines(sub("avc model (#565)", "avc model (#999)", map_text, fixed = TRUE),
           map_path)
check("a cited default-start test not in the run is refused",
      outcome(feature_matrix(base)), "cited default-start test not in the run")
writeLines(map_text, map_path)

setwd(root)
unlink(scratch, recursive = TRUE)
if (failures) {
  message(failures, " OF ", n_cases, " CASES FAILED")
  quit(save = "no", status = 1)
}
message("ALL ", n_cases, " CASES PASSED")
