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

# A hash of what the uncommitted changes say, not only which paths they touch:
# an already dirty file edited again keeps the same status line (#689).
suite_dirty_content <- function(dirty = suite_dirty()) {
  if (!length(dirty)) return("")
  f <- tempfile()
  on.exit(unlink(f))
  untracked <- suite_git("ls-files", "--others", "--exclude-standard", "--",
                         suite_paths)
  writeLines(c(suite_git("diff", "HEAD", "--binary", "--", suite_paths),
               if (length(untracked)) {
                 paste(untracked, suite_git("hash-object", "--", untracked))
               }), f)
  unname(tools::md5sum(f))
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

# What a run depends on: the commit, uncommitted changes, the version and the
# installed dependencies.
suite_identity <- function() {
  dirty <- suite_dirty()
  list(
    commit = suite_git("rev-parse", "HEAD"),
    dirty = dirty,
    dirty_content = suite_dirty_content(dirty),
    version = read.dcf("DESCRIPTION", "Version")[1, 1],
    dependencies = suite_dependencies()
  )
}

suite_provenance <- function(started, identity = suite_identity()) {
  c(identity, list(started = started, finished = Sys.time(),
                   not_cran = Sys.getenv("NOT_CRAN")))
}

# Runs `run()` and attaches its provenance. The identity is read before the
# run and again after it: read only after, a run that started on an older or
# dirty tree, committed or cleaned while it ran, would carry the later clean
# commit (#689). A run whose identity moved is refused, not saved.
suite_record <- function(run) {
  before <- suite_identity()
  started <- Sys.time()
  res <- run()
  after <- suite_identity()
  moved <- names(before)[!mapply(identical, before, after[names(before)])]
  if (length(moved)) {
    stop(structure(
      class = c("suite_changed", "error", "condition"),
      list(message = paste0(
        "the tree changed while the suite ran (", paste(moved, collapse = ", "),
        "), so no single state of it produced these results. Leave the tree ",
        "alone and run again."
      ), call = NULL)
    ))
  }
  attr(res, "provenance") <- suite_provenance(started, before)
  res
}
