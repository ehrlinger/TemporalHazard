# SETG1's rewrites of an early phase that starts at M = 0 (#471).
#
# PROC HAZARD does not fit `PARMS M=0 NU=1` with M free: SETG1 fixes M at 0
# (setg1.c:713-731), and for NU < 0 likewise (:651-667), unless NU is fixed
# and M free, when it moves M to 1 instead. With M fixed and NU = 0 it moves
# NU to 1 (:700-707). The oracle is the HAZARD binary's own "Initial Parameter
# Values" table and, where the fit runs, its final estimates
# (data-raw/setg1-mzero-oracle.R).

.mz_oracle <- function() {
  utils::read.csv(test_path("fixtures", "setg1-mzero-oracle.csv"),
                  comment.char = "#", stringsAsFactors = FALSE)
}
.mz_data <- function() {
  utils::read.csv(test_path("fixtures", "setg1-mzero-data.csv"),
                  comment.char = "#")
}
.mz_job <- function(parms) {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste0("%HAZARD( PROC HAZARD DATA=D; EVENT DEAD; TIME TT; ",
                    "PARMS ", parms, "; );"), f)
  w <- character(0)
  job <- withCallingHandlers(hzr_translate_sas(f), warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  job$warnings <- w
  job
}
# The emitted early phase, as hzr_phase() builds it.
.mz_phase <- function(job) eval(job$calls$fit[[3L]]$phases[[2L]])
# Run every emitted chunk on D, as a rendered document would.
.mz_run <- function(job, D) {
  env <- new.env(parent = .render_parent())
  env$D <- D
  utils::capture.output(for (nm in names(job$calls)) {
    suppressWarnings(eval(job$calls[[nm]], env))
  })
  env$fit
}

test_that("the emitted early phase starts and fixes as SETG1 leaves it (#471)", {
  oracle <- .mz_oracle()
  # Coverage before comparison: each rewrite, and the known negative.
  expect_true(all(c("mzero_nupos", "mzero_nuneg", "mzero_nuzero_fixm",
                    "mzero_nupos_fixnu", "mzero_nuneg_fixnu",
                    "control_mone") %in% oracle$id))
  expect_true(any(oracle$m_estimated == "No" & !grepl("FIXM", oracle$parms)))
  expect_true(any(oracle$m_used == 1 & grepl("M=0", oracle$parms)))
  for (k in seq_len(nrow(oracle))) {
    o <- oracle[k, ]
    ph <- .mz_phase(.mz_job(o$parms))
    fixed <- ph$fixed %||% character(0)
    expect_equal(ph$nu, o$nu_used, info = o$id)
    expect_equal(ph$m, o$m_used, info = o$id)
    expect_identical("m" %in% fixed, o$m_estimated == "No", info = o$id)
    expect_identical("nu" %in% fixed, o$nu_estimated == "No", info = o$id)
  }
})

test_that("each rewrite is recorded once, and no refusal is emitted (#471)", {
  oracle <- .mz_oracle()
  for (k in seq_len(nrow(oracle))) {
    o <- oracle[k, ]
    job <- .mz_job(o$parms)
    written_m_free <- !grepl("FIXM", o$parms)
    written_nu <- as.numeric(sub(".*NU=(-?[0-9.]+).*", "\\1", o$parms))
    written_m <- as.numeric(sub(".*M=(-?[0-9.]+).*", "\\1", o$parms))
    rewritten <- (written_m_free && o$m_estimated == "No") ||
      o$m_used != written_m || o$nu_used != written_nu
    expect_identical(NROW(job$untranslated), as.integer(rewritten),
                     info = o$id)
    expect_length(job$warnings, as.integer(rewritten))
    expect_false(any(grepl("^refusal", names(job$calls))), info = o$id)
    if (rewritten) {
      expect_match(job$untranslated$reason, "setg1.c:", fixed = TRUE,
                   info = o$id)
    }
  }
})

test_that("#471's three spellings emit one model (#471)", {
  calls <- lapply(c("MUE=0.2 THALF=1 M=0 NU=1", "MUE=0.2 THALF=1 M=0 NU=0 FIXM",
                    "MUE=0.2 THALF=1 M=0 NU=1 FIXM"),
                  function(p) .mz_job(p)$calls$fit)
  expect_identical(calls[[2L]], calls[[1L]])
  expect_identical(calls[[3L]], calls[[1L]])
  # The binary agrees they are one job: one log-likelihood on synth.
  oracle <- .mz_oracle()
  ll <- oracle$synth_ll[startsWith(oracle$id, "issue_")]
  expect_length(ll, 3L)
  expect_true(all(is.finite(ll)))
  expect_identical(length(unique(ll)), 1L)
})

test_that("the emitted fit reproduces PROC HAZARD's estimates (#471)", {
  skip_on_cran()
  oracle <- .mz_oracle()
  D <- .mz_data()
  expect_identical(nrow(D), 300L)
  runs <- oracle[oracle$mzero == "runs", ]
  # The M-fixed rows, the M-moved rows and the both-fixed row all ran.
  expect_gte(nrow(runs), 8L)
  for (k in seq_len(nrow(runs))) {
    o <- runs[k, ]
    fit <- .mz_run(.mz_job(o$parms), D)
    th <- fit$fit$theta
    # The listing prints the log-likelihood to 3 decimals (or fewer).
    expect_lt(abs(fit$fit$objective - o$mzero_ll), 0.006)
    expect_equal(exp(th[[2L]]), o$mzero_thalf, tolerance = 1e-4, info = o$id)
    expect_equal(th[[3L]], o$mzero_nu, tolerance = 1e-4, info = o$id)
    expect_equal(exp(th[[1L]]), o$mzero_mue, tolerance = 1e-4, info = o$id)
    if (o$m_estimated == "No") expect_identical(th[[4L]], o$mzero_m)
  }
  # KNOWN POSITIVE for the comparison: the model main emitted for
  # `M=0 NU=1`, M free, is a different model, and the check above can tell.
  free <- suppressWarnings(hazard(
    data = D, time = TT, status = DEAD, fit = TRUE, dist = "multiphase",
    phases = list(hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0)),
    theta = c(log(1), log(0.5), 1, 0)))
  ll <- oracle$mzero_ll[oracle$id == "mzero_nupos"]
  expect_gt(abs(free$fit$objective - ll), 0.1)
})
