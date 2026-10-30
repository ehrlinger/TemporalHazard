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

# What a run depends on: the commit, uncommitted changes, the version and the
# installed dependencies.
suite_identity <- function() {
  list(
    commit = suite_git("rev-parse", "HEAD"),
    dirty = suite_dirty(),
    version = read.dcf("DESCRIPTION", "Version")[1, 1],
    dependencies = suite_dependencies()
  )
}

suite_provenance <- function(started, identity = suite_identity()) {
  c(identity, list(started = started, finished = Sys.time(),
                   not_cran = Sys.getenv("NOT_CRAN")))
}

suite_refuse <- function(...) {
  stop(structure(
    class = c("suite_refused", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# Runs `run()` in a private export of HEAD and attaches its provenance. A
# comparison of the tree before and after a run cannot see an edit undone
# between the two, so the run does not read the tree at all (#689): nothing
# done to it while the suite runs, a commit included, can reach the results.
# The export holds only what is committed, so uncommitted changes are refused
# rather than left silently untested. Installed packages cannot be exported,
# so they are compared, and a change between the two reads is the limit of it.
suite_record <- function(run) {
  before <- suite_identity()
  if (length(before$dirty)) {
    suite_refuse("the tree has uncommitted changes a test reads (",
                 before$dirty[1L], "). The suite runs on an export of the ",
                 "commit and would not see them. Commit them and run again.")
  }
  export <- tempfile("suite-export-")
  on.exit(unlink(export, recursive = TRUE), add = TRUE)
  dir.create(export)
  archive <- paste0(export, ".tar")
  on.exit(unlink(archive), add = TRUE)
  suite_git("archive", "--format=tar", "-o", shQuote(archive), before$commit)
  utils::untar(archive, exdir = export)

  started <- Sys.time()
  home <- setwd(export)
  res <- tryCatch(run(), finally = setwd(home))
  deps <- suite_dependencies()
  if (!identical(deps, before$dependencies)) {
    suite_refuse("the installed packages changed while the suite ran, so no ",
                 "single set of them produced these results. Run again.")
  }
  attr(res, "provenance") <- suite_provenance(started, before)
  res
}
