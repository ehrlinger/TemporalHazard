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
  # Oracle sanity check, not code coverage: the binary agrees they are one
  # job, with one log-likelihood on synth. Only the fixture can fail this.
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
    if (o$m_estimated == "No") {
      expect_identical(th[[4L]], o$mzero_m)
    } else {
      # Where M is estimated (moved to 1, then free), compare it too. The
      # likelihood is flat in M near 0, so the scale is max(1, |M|), not |M|.
      expect_lt(abs(th[[4L]] - o$mzero_m), 1e-4 * max(1, abs(o$mzero_m)))
    }
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

# --- SETG1's other rewrites (#601) -------------------------------------------
# M < 0 and NU < 0 has its signs flipped (setg1.c:581-613); NU = 0 with M
# nonzero and both free has NU fixed at 0, and M's sign flipped when M > 0
# (:631-634, :763-770); M = 0, NU = 0 with NU fixed starts M at 1 (:692-699).
# The oracle is the binary on a 36-job grid: M in {-1, 0, 0.5}, NU in
# {-1, 0, 1}, with no flag, FIXM, FIXNU and both.

.mz_grid <- function() {
  utils::read.csv(test_path("fixtures", "setg1-rewrite-grid.csv"),
                  comment.char = "#", stringsAsFactors = FALSE)
}
# The value of `key=` in each PARMS string. A string it cannot read stops the
# test rather than becoming NA, which would silently drop rows from the
# cases a test selects; an exponent is read whole, not truncated at the `e`.
.mz_operand <- function(parms, key) {
  pat <- paste0("(^|.* )", key, "=([-+]?[0-9.]+([eE][-+]?[0-9]+)?)( .*)?$")
  hit <- grepl(pat, parms)
  if (!all(hit)) {
    stop("`", key, "=` could not be read from: ",
         paste(parms[!hit], collapse = "; "), call. = FALSE)
  }
  as.numeric(sub(pat, "\\2", parms))
}

test_that("every emitted early phase is the one SETG1 leaves, on the grid (#601)", {
  grid <- .mz_grid()
  # Coverage before comparison: the whole grid, each case the change reaches,
  # and the jobs SETG1 refuses.
  expect_identical(nrow(grid), 36L)
  m <- .mz_operand(grid$parms, "M")
  nu <- .mz_operand(grid$parms, "NU")
  refused <- grid$mzero == "refused"
  expect_identical(sum(refused), 3L)
  expect_true(any(m < 0 & nu < 0 & !refused))
  expect_true(any(m > 0 & nu == 0 & grid$m_used < 0, na.rm = TRUE))
  expect_true(any(nu == 0 & m != 0 & grid$nu_estimated == "No" &
                    !grepl("FIXNU", grid$parms), na.rm = TRUE))
  for (k in seq_len(nrow(grid))) {
    o <- grid[k, ]
    job <- .mz_job(o$parms)
    if (refused[k]) {
      expect_true(any(grepl("^refusal", names(job$calls))), info = o$parms)
      next
    }
    # The phase builds: unmirrored, M < 0 with NU < 0 errored here.
    expect_no_error(ph <- .mz_phase(job))
    fixed <- ph$fixed %||% character(0)
    expect_equal(ph$nu, o$nu_used, info = o$parms)
    expect_equal(ph$m, o$m_used, info = o$parms)
    expect_identical("m" %in% fixed, o$m_estimated == "No", info = o$parms)
    expect_identical("nu" %in% fixed, o$nu_estimated == "No", info = o$parms)
  }
})

test_that("the sign-flipped and NU-fixed fits reproduce the binary (#601)", {
  skip_on_cran()
  grid <- .mz_grid()
  D <- .mz_data()
  m <- .mz_operand(grid$parms, "M")
  nu <- .mz_operand(grid$parms, "NU")
  free_m <- !grepl("FIXM", grid$parms)
  free_nu <- !grepl("FIXNU", grid$parms)
  touched <- (m < 0 & nu < 0) | (nu == 0 & m != 0 & free_m) |
    (m == 0 & nu == 0 & !free_nu & free_m)
  runs <- grid[touched & grid$mzero == "runs", ]
  expect_gte(nrow(runs), 5L)
  for (k in seq_len(nrow(runs))) {
    o <- runs[k, ]
    fit <- .mz_run(.mz_job(o$parms), D)
    expect_lt(abs(fit$fit$objective - o$mzero_ll), 0.006)
    # Some of these optima lie on a flat ridge (M = -1 NU = -1 FIXM: THALF
    # differs by 5e-4 relative, the log-likelihood by 1e-6), so the
    # estimates are compared by what they ARE: SAS's printed point,
    # evaluated by R's likelihood, is R's maximum.
    sas <- c(log(o$mzero_mue), log(o$mzero_thalf), o$mzero_nu, o$mzero_m)
    at_sas <- suppressWarnings(hzr_evaluate(fit, theta = sas))$logLik
    expect_lt(abs(at_sas - fit$fit$objective), 1e-4)
  }
})

test_that("a rewrite row is conditional when a macro could change SETG1's case (#601)", {
  row <- function(parms, start) {
    r <- .mz_job(parms)$untranslated$reason
    r[startsWith(r, start)]
  }
  plain <- row("MUE=0.2 THALF=1 M=0 NU=1", "SETG1 fixes")
  macro <- row("MUE=0.2 THALF=1 M=0 NU=1 &FLAGS", "SETG1 fixes")
  expect_length(plain, 1L)
  expect_length(macro, 1L)
  expect_match(plain, "as PROC HAZARD does:", fixed = TRUE)
  expect_no_match(plain, "macro", fixed = TRUE)
  expect_match(macro, "as PROC HAZARD does if the macro reference",
               fixed = TRUE)
  # The same for a moved starting value.
  plain <- row("MUE=0.2 THALF=1 M=-1 NU=-1", "SETG1 replaces")
  macro <- row("MUE=0.2 THALF=1 M=-1 NU=-1 &FLAGS", "SETG1 replaces")
  expect_length(plain, 1L)
  expect_match(plain, "as PROC HAZARD does, and", fixed = TRUE)
  expect_match(macro, "if the macro reference", fixed = TRUE)
})
