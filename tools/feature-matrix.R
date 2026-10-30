#!/usr/bin/env Rscript
# feature-matrix.R -- the data behind the feature test matrix article (#682).
#
# WHY
#   A hand-kept table of what is tested goes stale, and a stale "tested" is an
#   output that looks like a result and is not. So the counts are read from a
#   saved test run, and the run is checked against the tree before anything
#   is written. Only the status column and the notes are human judgement.
#
# INPUTS (all in tools/)
#   feature-matrix-map.csv         one row per feature: test-file pattern, a
#                                  cited default-start test, status, notes
#   feature-matrix-references.csv  every test that compares a feature with an
#                                  external reference, cited by file and name
#   feature-matrix-computed.csv    the names of tests whose names are built
#                                  at run time, which no parser can read
#   feature-matrix-assigned.csv    tests that count towards one feature from a
#                                  file no feature's pattern claims, by name
#
# WHAT IT REFUSES
#   - a run without provenance (tools/save-suite.R or tools/release-suite.R
#     attaches it), a run made on an uncommitted tree, a run whose commit
#     differs from HEAD in anything a test reads (suite-provenance.R names the
#     paths), a tree with uncommitted changes there now, and a run made under
#     other installed dependency versions;
#   - a run with any failure or error;
#   - a run whose tests differ from the tests on disk, computed names included;
#   - a test that no feature counts, and a pattern that matches nothing;
#   - a cited or assigned test, or a default-start test, not in the run.
#
# USAGE (from the package root)
#   Rscript tools/save-suite.R results.rds      # or tools/release-suite.R
#   Rscript tools/feature-matrix.R results.rds
#   Writes vignettes/articles/feature-matrix.csv.

source(file.path("tools", "suite-provenance.R"))

matrix_refuse <- function(...) {
  stop(structure(
    class = c("feature_matrix_refused", "error", "condition"),
    list(message = paste0("FEATURE MATRIX REFUSED: ", ...), call = NULL)
  ))
}

# Literal test_that() names in a file, and how many calls compute theirs.
matrix_literal_tests <- function(path) {
  out <- character()
  computed <- 0L
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    if (identical(e[[1L]], as.name("test_that"))) {
      if (is.character(e[[2L]])) out <<- c(out, e[[2L]])
      else computed <<- computed + 1L
    }
    for (a in as.list(e)[-1L]) if (!missing(a)) walk(a)
  }
  for (e in parse(path, keep.source = FALSE)) walk(e)
  list(names = out, computed = computed)
}

matrix_check_provenance <- function(prov) {
  if (!is.list(prov) || is.null(prov$commit)) {
    matrix_refuse("the run carries no provenance. Make it with ",
                  "tools/save-suite.R or tools/release-suite.R.")
  }
  if (length(prov$dirty)) {
    matrix_refuse("the run was made on a tree with uncommitted changes (",
                  prov$dirty[1L], "), so no commit describes it.")
  }
  now <- suite_dirty()
  if (length(now)) {
    matrix_refuse("the tree has uncommitted changes a test reads (", now[1L],
                  "). Commit them and run the suite again.")
  }
  head <- suite_git("rev-parse", "HEAD")
  if (!identical(prov$commit, head)) {
    changed <- suite_git("diff", "--name-only", prov$commit, head, "--",
                         suite_paths)
    if (length(changed)) {
      matrix_refuse("the run is of ", substr(prov$commit, 1L, 8L), ", and ",
                    length(changed), " file(s) a test reads have changed ",
                    "since (", changed[1L], "). Run the suite again.")
    }
  }
  deps <- suite_dependencies()
  was <- prov$dependencies
  common <- intersect(names(deps), names(was))
  moved <- common[!mapply(identical, deps[common], was[common])]
  if (length(moved) || !setequal(names(deps), names(was))) {
    first <- c(moved, setdiff(union(names(deps), names(was)), common))[1L]
    matrix_refuse("the installed dependencies differ from the run's (",
                  first, ": ", was[first], " then, ", deps[first], " now).")
  }
  invisible(TRUE)
}

matrix_check_tree <- function(res, files, test_dir, computed) {
  for (f in files) {
    lit <- matrix_literal_tests(file.path(test_dir, f))
    expected <- c(lit$names, computed$test[computed$file == f])
    if (lit$computed > 0L && !any(computed$file == f)) {
      matrix_refuse(f, " builds ", lit$computed, " test name(s) at run time ",
                    "and tools/feature-matrix-computed.csv lists none.")
    }
    in_run <- res$test[res$file == f]
    if (!length(in_run)) matrix_refuse(f, " has no rows in the run.")
    gone <- setdiff(expected, in_run)
    extra <- setdiff(in_run, expected)
    if (length(gone) || length(extra)) {
      matrix_refuse(f, " differs from the run (", length(gone),
                    " expected test(s) not run, ", length(extra),
                    " run but not expected), so the run is of another tree. ",
                    "First: ", c(gone, extra)[1L])
    }
  }
  stale <- setdiff(unique(res$file), files)
  if (length(stale)) {
    matrix_refuse("the run has files no longer on disk: ",
                  paste(stale, collapse = ", "))
  }
  invisible(TRUE)
}

# One row per feature. `res` is a saved run, checked here before use.
feature_matrix <- function(res, root = ".") {
  need <- c("file", "test", "passed", "failed", "skipped", "error")
  if (!is.data.frame(res) || !all(need %in% names(res)) || nrow(res) == 0L) {
    matrix_refuse("the results are not as.data.frame(devtools::test()).")
  }
  matrix_check_provenance(attr(res, "provenance"))
  if (any(res$failed > 0L) || any(res$error)) {
    matrix_refuse(sum(res$failed > 0L | res$error), " test(s) failed or ",
                  "errored in that run; a matrix over a red run would read ",
                  "as tested.")
  }

  read <- function(name) {
    utils::read.csv(file.path(root, "tools", name), stringsAsFactors = FALSE,
                    na.strings = character())
  }
  map <- read("feature-matrix-map.csv")
  refs <- read("feature-matrix-references.csv")
  computed <- read("feature-matrix-computed.csv")
  assigned <- read("feature-matrix-assigned.csv")

  test_dir <- file.path(root, "tests", "testthat")
  files <- list.files(test_dir, pattern = "^test-.*\\.R$")
  matrix_check_tree(res, files, test_dir, computed)

  unknown <- setdiff(c(refs$feature, assigned$feature), map$feature)
  if (length(unknown)) matrix_refuse("references or assignments name an ",
                                     "unknown feature: ", unknown[1L])
  cited <- function(file, test, what) {
    hit <- res[res$file == file & res$test == test, ]
    if (nrow(hit) != 1L) matrix_refuse(what, " not in the run: ", file, "::",
                                       test)
    hit
  }

  stem <- sub("\\.R$", "", sub("^test-", "", files))
  key <- function(d) paste(d$file, d$test, sep = "::")
  counted <- logical(nrow(res))
  rows <- lapply(seq_len(nrow(map)), function(i) {
    m <- map[i, ]
    hit <- grepl(m$files, stem, perl = TRUE)
    if (!any(hit)) matrix_refuse(m$feature, " pattern '", m$files,
                                 "' matches no test file.")
    fr <- refs[refs$feature == m$feature, ]
    fa <- assigned[assigned$feature == m$feature, ]
    ref_rows <- do.call(rbind, c(
      list(res[0L, ]),
      Map(cited, fr$file, fr$test, "cited reference test")
    ))
    Map(cited, fa$file, fa$test, "assigned test")
    # A feature's reference comparisons are tests of it too, but only the
    # cited ones: the rest of a cited file tests other features (#689).
    mine <- res$file %in% files[hit] | key(res) %in% c(key(fr), key(fa))
    counted[mine] <<- TRUE
    r <- res[mine, ]
    own <- unique(r$file)
    kinds <- paste(sort(unique(fr$reference)), collapse = ", ")
    # The references file is a selection, so an empty one is "none cited",
    # not a finding that no comparison exists (#689).
    reference <- if (!nrow(fr)) {
      "none cited"
    } else {
      paste0(sum(!ref_rows$skipped), " of ", nrow(ref_rows), " ran (", kinds,
             ")")
    }
    default <- "none"
    if (nzchar(m$from_default)) {
      parts <- strsplit(m$from_default, "::", fixed = TRUE)[[1L]]
      d <- cited(parts[1L], parts[2L], "cited default-start test")
      default <- if (d$skipped) "cited test skipped" else "cited test passed"
    }
    data.frame(
      area = m$area, feature = m$feature, status = m$status,
      files = length(own), tests_run = sum(!r$skipped),
      expectations = sum(r$passed), skipped = sum(r$skipped),
      reference = reference, default_start = default, notes = m$notes,
      stringsAsFactors = FALSE
    )
  })
  if (!all(counted)) {
    lost <- res[!counted, ]
    matrix_refuse("no feature claims ", nrow(lost), " test(s) in ",
                  paste(unique(lost$file), collapse = ", "), " (first: ",
                  lost$test[1L], "). Add the file to a feature's pattern in ",
                  "tools/feature-matrix-map.csv, or the test to ",
                  "tools/feature-matrix-assigned.csv.")
  }
  do.call(rbind, rows)
}

write_feature_matrix <- function(res, target) {
  out <- feature_matrix(res)
  prov <- attr(res, "provenance")
  dir.create(dirname(target), showWarnings = FALSE)
  con <- file(target, "w")
  on.exit(close(con))
  writeLines(c(
    paste0("# measured_at: ", substr(prov$commit, 1L, 8L)),
    paste0("# measured_on: ", format(prov$finished, "%Y-%m-%d")),
    paste0("# version: ", prov$version),
    paste0("# suite: ", sum(!res$skipped), " tests run, ", sum(res$passed),
           " expectations, ", sum(res$skipped), " tests skipped")
  ), con)
  utils::write.csv(out, con, row.names = FALSE)
  invisible(out)
}

if (sys.nframe() == 0L) {
  if (!file.exists("DESCRIPTION") ||
        read.dcf("DESCRIPTION", "Package")[1, 1] != "TemporalHazard") {
    message("FEATURE MATRIX REFUSED: run this from the package root.")
    quit(save = "no", status = 1)
  }
  results_path <- commandArgs(trailingOnly = TRUE)[1]
  if (is.na(results_path) || !file.exists(results_path)) {
    message("FEATURE MATRIX REFUSED: give the saved results (an .rds).")
    quit(save = "no", status = 1)
  }
  target <- file.path("vignettes", "articles", "feature-matrix.csv")
  out <- tryCatch(
    write_feature_matrix(readRDS(results_path), target),
    feature_matrix_refused = function(e) {
      message(conditionMessage(e))
      quit(save = "no", status = 1)
    }
  )
  message("Wrote ", target, ": ", nrow(out), " features")
}
