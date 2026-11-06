## Behavioural census: run the case battery under ONE installed version.
##
## Usage:
##   Rscript tools/release-census/census-run.R <lib-dir> <expected-version> \
##           <cases-file> <out-rds>
##
## The library guard is the first thing that runs and it is not optional. The
## R system library on this machine still holds TemporalHazard 1.1.0 (June
## 2026), so an R process that loses its intended library silently loads that
## instead and the census then compares 1.1.0 against 1.1.0 and reports
## "all identical" -- a populated result over nothing, which is this package's
## signature defect.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: census-run.R <lib-dir> <expected-version> <cases-file> <out-rds>",
       call. = FALSE)
}
lib <- normalizePath(args[[1]], mustWork = TRUE)
expected <- args[[2]]
cases_file <- normalizePath(args[[3]], mustWork = TRUE)
out_rds <- args[[4]]

# ---------------------------------------------------------------------------
# Library guard
# ---------------------------------------------------------------------------

.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages(library(TemporalHazard))

pkg_dir <- normalizePath(find.package("TemporalHazard"))
pkg_lib <- normalizePath(dirname(pkg_dir))
version <- as.character(utils::packageVersion("TemporalHazard"))

cat("=== LIBRARY GUARD ===\n")
cat("requested lib   :", lib, "\n")
cat("loaded from     :", pkg_lib, "\n")
cat("package dir     :", pkg_dir, "\n")
cat("packageVersion  :", version, "\n")
cat("expected version:", expected, "\n")
cat(".libPaths()     :", paste(.libPaths(), collapse = " | "), "\n")

if (!identical(pkg_lib, lib)) {
  stop("LIBRARY GUARD FAILED: TemporalHazard loaded from '", pkg_lib,
       "', not the requested '", lib, "'. Refusing to census the wrong tree.",
       call. = FALSE)
}
if (!identical(version, expected)) {
  stop("LIBRARY GUARD FAILED: loaded TemporalHazard ", version,
       ", expected ", expected, ".", call. = FALSE)
}
cat("LIBRARY GUARD: PASS\n\n")

# ---------------------------------------------------------------------------
# Environment record: anything that could move a numeric result
# ---------------------------------------------------------------------------

env_record <- list(
  version = version,
  lib = lib,
  r_version = R.version.string,
  platform = R.version$platform,
  blas = tryCatch(sessionInfo()$BLAS, error = function(e) NA_character_),
  lapack = tryCatch(sessionInfo()$LAPACK, error = function(e) NA_character_),
  survival = as.character(utils::packageVersion("survival")),
  numDeriv = tryCatch(as.character(utils::packageVersion("numDeriv")),
                      error = function(e) "ABSENT"),
  rng = RNGkind(),
  # SAS fixtures. Recorded, not assumed: the report must say whether they
  # were available rather than leave a silent skip.
  qhsstudies_mounted = dir.exists("/Volumes/qhsstudies"),
  hazard_repo = Sys.getenv("HAZARD_REPO", unset = ""),
  hazard_examples_dir = Sys.getenv("HAZARD_EXAMPLES_DIR", unset = "")
)
cat("=== ENVIRONMENT ===\n")
str(env_record)
cat("\n")
if (identical(env_record$numDeriv, "ABSENT")) {
  cat("NOTE: numDeriv is ABSENT. The analytic Hessian declines for left- and\n",
      "      interval-censored rows by design and falls back to numDeriv, so\n",
      "      those fits will produce no standard errors here. This is an\n",
      "      environment difference, not a code difference -- but it must be\n",
      "      THE SAME in both runs or the vcov comparison is meaningless.\n",
      sep = "")
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

source(cases_file)

dat <- census_data()
cases <- census_cases()

# Which exported names exist under THIS version. A case needing a name that
# does not exist is ABSENT, which is reported as its own outcome.
exports <- ls(asNamespace("TemporalHazard"))
exported <- getNamespaceExports("TemporalHazard")

results <- list()
n_ok <- 0L
n_error <- 0L
n_absent <- 0L

for (nm in names(cases)) {
  case <- cases[[nm]]
  needs <- if (is.null(case$needs)) character() else case$needs
  missing <- setdiff(needs, exported)

  if (length(missing)) {
    results[[nm]] <- list(status = "absent", missing = missing,
                          kp = case$kp, gate_kp = case$gate_kp)
    n_absent <- n_absent + 1L
    cat(sprintf("[ABSENT] %-28s needs %s\n", nm, paste(missing, collapse = ",")))
    next
  }

  warns <- character()
  val <- withCallingHandlers(
    tryCatch(case$fn(dat), error = function(e) {
      structure(list(message = conditionMessage(e)), class = "census_error")
    }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(val, "census_error")) {
    results[[nm]] <- list(status = "error", message = val$message,
                          warnings = warns, kp = case$kp, gate_kp = case$gate_kp)
    n_error <- n_error + 1L
    cat(sprintf("[ERROR ] %-28s %s\n", nm,
                substr(gsub("\\s+", " ", val$message), 1, 110)))
  } else {
    results[[nm]] <- list(status = "ok", probe = val, warnings = warns,
                          kp = case$kp, gate_kp = case$gate_kp)
    n_ok <- n_ok + 1L
    cat(sprintf("[OK    ] %-28s %s\n", nm,
                if (length(warns)) sprintf("(%d warning(s))", length(warns)) else ""))
  }
}

cat("\n=== COVERAGE (this run) ===\n")
cat("cases declared:", length(cases), "\n")
cat("ok            :", n_ok, "\n")
cat("error         :", n_error, "\n")
cat("absent        :", n_absent, "\n")
stopifnot(n_ok + n_error + n_absent == length(cases))

# A run in which nothing succeeded is a broken harness, not a result.
if (n_ok == 0L) {
  stop("HARNESS FAILURE: zero cases succeeded under ", version,
       ". An empty census is not a pass.", call. = FALSE)
}

saveRDS(list(env = env_record, results = results,
             counts = list(declared = length(cases), ok = n_ok,
                           error = n_error, absent = n_absent)),
        out_rds)
cat("\nwrote:", out_rds, "\n")
