#!/usr/bin/env Rscript
# oracle.R -- bit-identical differential oracle for the phase decomposition.
#
# WHY
#   A performance refactor of code already known to be correct claims to
#   change no answer, so the test that fits it is exact agreement between the
#   old code and the new (inst/dev/CORRECTNESS-STRATEGY.md, section 2). This
#   evaluates one grid on both and compares every value with identical(). It
#   runs both sides on the same machine, so last-bit differences between
#   platforms' maths libraries cancel, which a frozen reference file could not
#   promise. The 1.2.14 speed-up (#631) was gated with this grid.
#
# WHAT IT COVERS
#   hzr_phase_cumhaz(), hzr_phase_hazard(), .hzr_phase_derivatives() and
#   hzr_decompos() over 216 shapes (t_half x nu x m, every case and the error
#   cells) and seven kinds of time vector (tiny to huge values, 0, NA, NaN,
#   Inf, negative, named, matrix, integer, classed), for "cdf" and "hazard"
#   phases; plus hzr_log1pexp() over awkward inputs. Values, attributes,
#   errors and warnings all count. Not covered: g3 phases, whole fits and the
#   bootstrap.
#
# USAGE (from the package root)
#   Rscript tools/oracle.R [base-ref]
#   Compares the working tree, uncommitted changes included, with base-ref
#   (default origin/main), checked out in a temporary worktree. Exits 0 and
#   prints "ORACLE: IDENTICAL" only when every value matches; otherwise it
#   lists the differing shapes and fields and exits 1; it exits 2, "ORACLE NOT
#   RUN", when it cannot make the comparison. An intended numeric
#   change will differ, and the listing is how to check it differs only where
#   it should.

shapes <- expand.grid(
  t_half = c(1e-8, 0.3, 1, 50),
  nu = c(-2, -0.01, 0, 1e-4, 0.5, 3),
  m = c(-1100, -60, -1, -1e-6, 0, 1e-6, 0.7, 5, 2000)
)

evaluate_grid <- function(tree, out) {
  suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE,
                                     helpers = FALSE))
  ns <- asNamespace("TemporalHazard")
  times <- list(
    plain = c(1e-300, 1e-12, 1e-5, 0.01, 0.5, 1, 2, 10, 1e3, 1e8, 1e300),
    edge  = c(0, NA, NaN, Inf, 1, -1),
    named = c(a = 0.1, b = 1, c = 30),
    mat   = matrix(c(0.1, 1, 5, 50), 2),
    int   = 1:5,
    cls   = structure(c(0.2, 2, 20), class = "foo", units = "d"),
    one   = 3
  )
  res <- list()
  warn <- list()
  for (tn in names(times)) {
    for (i in seq_len(nrow(shapes))) {
      for (ty in c("cdf", "hazard")) {
        s <- shapes[i, ]
        key <- paste(tn, i, ty)
        f <- function(expr) {
          withCallingHandlers(
            tryCatch(expr, error = function(e) {
              paste("ERR", conditionMessage(e))
            }),
            warning = function(w) {
              warn[[key]] <<- c(warn[[key]], conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          )
        }
        tm <- times[[tn]]
        res[[key]] <- list(
          ch = f(ns$hzr_phase_cumhaz(tm, s$t_half, s$nu, s$m, type = ty)),
          hz = f(ns$hzr_phase_hazard(tm, s$t_half, s$nu, s$m, type = ty)),
          pd = f(ns$.hzr_phase_derivatives(as.numeric(tm), s$t_half, s$nu,
                                           s$m, type = ty)),
          dc = if (ty == "hazard") f(ns$hzr_decompos(tm, s$t_half, s$nu, s$m))
        )
      }
    }
  }
  lp <- list(
    plain = c(-800, -50, -30.5, -1, 0, 1e-320, 1, 40, 800, NA, NaN, Inf, -Inf),
    int = c(-3L, 0L, NA, 5L), mat = matrix(c(-1, 2, NA, 4), 2),
    named = c(a = -1, b = 2), attr = structure(c(-1, 2), foo = "x"),
    cls = structure(c(-1, 2), class = "bar"), lgl = c(TRUE, NA),
    allna = c(NA_real_, NA_real_), na1 = NA, empty = numeric(0),
    df = data.frame(a = c(-1, 2)), cplx = NULL
  )
  res$lp <- lapply(lp, function(x) {
    tryCatch(ns$hzr_log1pexp(x), error = function(e) conditionMessage(e))
  })
  saveRDS(list(res = res, warn = warn), out)
}

args <- commandArgs(trailingOnly = TRUE)

# Worker mode: one side of the comparison, in a fresh R process, so the two
# trees never share a loaded namespace.
if (length(args) == 3L && args[1] == "--worker") {
  evaluate_grid(args[2], args[3])
  quit(save = "no", status = 0)
}

# A refusal is signalled as a condition so that run_oracle()'s on.exit()
# removes the temporary worktree on every path; on.exit() at the top level of
# an Rscript script never runs.
stop_oracle <- function(...) {
  stop(structure(class = c("oracle_not_run", "error", "condition"),
                 list(message = paste0(...), call = NULL)))
}

run_oracle <- function(args) {
  if (!file.exists("DESCRIPTION") ||
        read.dcf("DESCRIPTION", "Package")[1, 1] != "TemporalHazard") {
    stop_oracle("run this from the TemporalHazard package root.")
  }
  base_ref <- if (length(args)) args[1] else "origin/main"
  base_sha <- suppressWarnings(system2("git", c("rev-parse", "--short", base_ref),
                                       stdout = TRUE, stderr = FALSE))
  if (!length(base_sha)) stop_oracle("cannot resolve '", base_ref, "'.")
  head_sha <- system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)
  dirty <- length(system2("git", c("status", "--porcelain",
                                   "--untracked-files=no"), stdout = TRUE)) > 0L

  work <- tempfile("oracle-")
  dir.create(work)
  base_tree <- file.path(work, "base")
  on.exit({
    system2("git", c("worktree", "remove", "--force", base_tree),
            stdout = FALSE, stderr = FALSE)
    unlink(work, recursive = TRUE)
  }, add = TRUE)
  if (system2("git", c("worktree", "add", "--detach", "--quiet", base_tree,
                       base_ref)) != 0L) {
    stop_oracle("could not check out '", base_ref, "' in a worktree.")
  }

  script <- normalizePath(file.path("tools", "oracle.R"))
  side <- function(tree, out) {
    system2(file.path(R.home("bin"), "Rscript"),
            c(script, "--worker", shQuote(tree), shQuote(out)))
  }
  out_base <- file.path(work, "base.rds")
  out_new  <- file.path(work, "new.rds")
  if (side(base_tree, out_base) != 0L ||
        side(normalizePath("."), out_new) != 0L) {
    stop_oracle("a side of the comparison failed to evaluate.")
  }
  a <- readRDS(out_base)
  b <- readRDS(out_new)

  cat(sprintf("base %s (%s)  vs  working tree at %s%s\n", base_ref, base_sha,
              head_sha, if (dirty) " plus uncommitted changes" else ""))
  keys <- names(a$res)
  same <- mapply(identical, a$res[keys], b$res[keys])
  if (!identical(sort(keys), sort(names(b$res)))) {
    same[] <- FALSE
    cat("The two sides evaluated different key sets.\n")
  }
  warn_same <- identical(a$warn, b$warn)
  cat(sprintf("%d keys compared; %d differ; warnings %s.\n", length(keys),
              sum(!same), if (warn_same) "identical" else "DIFFER"))
  if (all(same) && warn_same) {
    cat("ORACLE: IDENTICAL\n")
    return(0L)
  }

  diff_keys <- keys[!same]
  grid_keys <- setdiff(diff_keys, "lp")
  if (length(grid_keys)) {
    i <- as.integer(vapply(strsplit(grid_keys, " "), `[`, "", 2))
    cat("\nDiffering grid keys by shape (nu by m):\n")
    print(table(nu = shapes$nu[i], m = shapes$m[i]))
    fields <- unlist(lapply(grid_keys, function(k) {
      x <- a$res[[k]]
      names(x)[!mapply(identical, x, b$res[[k]][names(x)])]
    }))
    cat("\nBy field (ch cumhaz, hz hazard, pd derivatives, dc decompos):\n")
    print(table(fields))
  }
  if ("lp" %in% diff_keys) cat("\nhzr_log1pexp() results differ.\n")
  cat(sprintf("\nORACLE: %d of %d keys differ\n", sum(!same), length(keys)))
  1L
}

status <- tryCatch(run_oracle(args), oracle_not_run = function(e) {
  message("ORACLE NOT RUN: ", conditionMessage(e))
  2L
})
quit(save = "no", status = status)
