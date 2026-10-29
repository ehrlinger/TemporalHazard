# suite-provenance.R -- what produced a saved test run (#682).
#
# Sourced by tools/save-suite.R and tools/release-suite.R, which attach the
# record to the results they save, and by tools/feature-matrix.R, which
# refuses a run whose record does not match the tree it is asked to describe.
# Test names cannot tell a run of one commit from a run of the next: a change
# to R/, a test body, a fixture or a dependency leaves them all as they were.

# Paths whose content a test run depends on. A change anywhere else (the
# article, NEWS, tools/) cannot change a test's result.
suite_paths <- c("R", "tests", "inst", "src", "data", "DESCRIPTION",
                 "NAMESPACE")

suite_git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("git ", paste(c(...), collapse = " "), " failed: ",
         paste(out, collapse = " "), call. = FALSE)
  }
  out
}

# Uncommitted changes, tracked or not, under the paths a run depends on.
suite_dirty <- function() {
  suite_git("status", "--porcelain", "--untracked-files=all", "--",
            suite_paths)
}

# Installed versions of R and of every package DESCRIPTION names.
suite_dependencies <- function() {
  d <- read.dcf("DESCRIPTION")
  fields <- intersect(c("Depends", "Imports", "Suggests", "LinkingTo"),
                      colnames(d))
  pkgs <- unlist(lapply(fields, function(f) {
    trimws(sub("\\(.*", "", strsplit(d[1, f], ",")[[1L]]))
  }))
  pkgs <- sort(setdiff(unique(pkgs[nzchar(pkgs)]), "R"))
  vers <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) NA_character_)
  }, character(1))
  c(R = paste(R.version$major, R.version$minor, sep = "."), vers)
}

suite_provenance <- function(started) {
  list(
    commit = suite_git("rev-parse", "HEAD"),
    dirty = suite_dirty(),
    version = read.dcf("DESCRIPTION", "Version")[1, 1],
    started = started,
    finished = Sys.time(),
    not_cran = Sys.getenv("NOT_CRAN"),
    dependencies = suite_dependencies()
  )
}
