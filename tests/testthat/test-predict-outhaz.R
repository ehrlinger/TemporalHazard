# Predicting from a model loaded out of a SAS `OUTHAZ=` dataset.
#
# The fixture is synthetic -- real SAS layout, invented numbers (see
# data-raw/outhaz_fixture.R). It is an early + constant fit with no
# covariates: E0 and C0 carry _STATUS_ = 1, L0 does not.
#
# Every value assertion below recomputes the prediction from the *file's own
# numbers* through hzr_decompos(), not through the reconstruction being
# tested. A reconstruction that mis-orders theta, or reads E0 as mu rather
# than log(mu), predicts smoothly and wrongly; only an independent
# computation can see that.

fixture_path <- function() {
  system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
}

# Write a copy of the fixture with `edits` applied to its data frame.
outhaz_variant <- function(edits) {
  d <- readRDS(fixture_path())
  d <- edits(d)
  path <- withr::local_tempfile(fileext = ".rds", .local_envir = parent.frame())
  saveRDS(d, path)
  path
}

# H(t) for the fixture's early + constant model, computed straight from the
# stored estimates: H(t) = exp(E0) G(t) + exp(C0) t.
fixture_cumhaz <- function(time, est, early = TRUE, constant = TRUE) {
  h <- rep(0, length(time))
  if (early) {
    g <- hzr_decompos(time, t_half = est[["THALF"]], nu = est[["NU"]],
                      m = est[["M"]])$G
    h <- h + exp(est[["E0"]]) * g
  }
  if (constant) h <- h + exp(est[["C0"]]) * time
  h
}

test_that("hzr_read_outhaz() returns a classed object", {
  got <- hzr_read_outhaz(fixture_path())
  expect_s3_class(got, "hzr_outhaz")
})

test_that("an all-fixed parameter set yields a 0x0 vcov, not NULL", {
  obj <- hzr_read_outhaz(fixture_path())
  # Documented hollow-object shape: check dimensions, never is.null().
  expect_true(is.matrix(obj$vcov))
  expect_equal(dim(obj$vcov), c(4L, 4L))
})

test_that("predict() works on a loaded OUTHAZ fit", {
  obj <- hzr_read_outhaz(fixture_path())
  got <- predict(obj, newdata = data.frame(time = c(1, 6, 12)),
                 type = "survival")
  expect_length(got, 3L)
  expect_true(all(got >= 0 & got <= 1))
  expect_true(all(diff(got) < 0))
})

test_that("the reconstructed model is the one the file describes", {
  obj <- hzr_read_outhaz(fixture_path())
  tt <- c(0.5, 1, 6, 12, 60)
  want <- fixture_cumhaz(tt, obj$estimates)

  expect_equal(predict(obj, newdata = data.frame(time = tt),
                       type = "cumulative_hazard"), want)
  expect_equal(predict(obj, newdata = data.frame(time = tt),
                       type = "survival"), exp(-want))

  # The two phases must both be contributing: strip the constant phase's
  # contribution and the numbers have to move. Without this the assertions
  # above pass on an early-only reconstruction.
  early_only <- fixture_cumhaz(tt, obj$estimates, constant = FALSE)
  expect_true(all(abs(want - early_only) > 1e-8))
})

test_that("phase membership follows _STATUS_ on E0/C0/L0, not the flags", {
  # G1FLAG is 2 in the fixture and stays 2; only C0's status changes. If
  # phase membership were read off the flags this would not move.
  path <- outhaz_variant(function(d) {
    d[["_STATUS_"]][d[["_NAME_"]] == "C0"] <- 0
    d
  })
  obj <- hzr_read_outhaz(path)
  full <- hzr_read_outhaz(fixture_path())
  tt <- c(1, 6, 12)

  got <- predict(obj, newdata = data.frame(time = tt),
                 type = "cumulative_hazard")
  expect_equal(got, fixture_cumhaz(tt, full$estimates, constant = FALSE))
  expect_true(all(got < predict(full, newdata = data.frame(time = tt),
                                type = "cumulative_hazard")))
})

test_that("a dataset with no phase in the model errors", {
  path <- outhaz_variant(function(d) {
    d[["_STATUS_"]][d[["_NAME_"]] %in% c("E0", "C0", "L0")] <- 0
    d
  })
  expect_error(
    predict(hzr_read_outhaz(path), newdata = data.frame(time = 1)),
    "no phase in the model"
  )
})

test_that("G1FLAG is checked against the M and NU estimates", {
  path <- outhaz_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "G1FLAG"] <- 1
    d
  })
  expect_error(
    predict(hzr_read_outhaz(path), newdata = data.frame(time = 1)),
    "G1FLAG"
  )
})

test_that("a non-zero DELTA is refused, not absorbed", {
  path <- outhaz_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "DELTA"] <- 0.25
    d
  })
  expect_error(
    predict(hzr_read_outhaz(path), newdata = data.frame(time = 1)),
    "DELTA"
  )
})

test_that("newdata is required -- an OUTHAZ dataset carries no times", {
  obj <- hzr_read_outhaz(fixture_path())
  expect_error(predict(obj), "newdata")
  expect_error(predict(obj, newdata = NULL), "newdata")
})

test_that("a dataset carrying covariates is refused by name", {
  # PROC HAZARD writes 3p + 11 parameter rows; the E-block carries the real
  # covariate names, the C and L blocks positional aliases.
  params <- c("DELTA", "THALF", "NU", "M", "TAU", "GAMMA", "ALPHA", "ETA",
              "E0", "AGE", "C0", "C01", "L0", "L01")
  flags <- c(G1FLAG = 2, FIXDEL0 = 1, FIXMNU1 = 0,
             G3FLAG = 0, FIXGE2 = 0, FIXGAE2 = 0)
  d <- data.frame(
    `_NAME_` = c(names(flags), params),
    `_EST_` = c(unname(flags), rep(0, length(params))),
    `_STATUS_` = c(rep(NA_real_, length(flags)), rep(1, length(params))),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (p in params) {
    d[[p]] <- c(rep(NA_real_, length(flags)), rep(0, length(params)))
  }
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(d, path)

  expect_error(
    predict(hzr_read_outhaz(path), newdata = data.frame(time = 1)),
    "AGE"
  )
})

test_that("se.fit maps the OUTHAZ covariance onto this package's scale", {
  # The stored covariance is on SAS's ESTIMATION scale: its free parameters
  # are (log THALF, log NU, E0, C0), where E0 and C0 are already log(mu).
  # Differentiate H(t) in exactly those coordinates and sandwich by the raw
  # stored matrix -- nothing here touches the reconstruction's own Jacobian.
  obj <- hzr_read_outhaz(fixture_path())
  est <- obj$estimates
  tt <- c(1, 6, 12, 60)

  h_of <- function(s) {
    exp(s[[3L]]) *
      hzr_decompos(tt, t_half = exp(s[[1L]]), nu = exp(s[[2L]]),
                   m = est[["M"]])$G +
      exp(s[[4L]]) * tt
  }
  s0 <- c(log(est[["THALF"]]), log(est[["NU"]]), est[["E0"]], est[["C0"]])
  jac <- vapply(seq_along(s0), function(k) {
    step <- 1e-6 * max(1, abs(s0[[k]]))
    hi <- s0
    lo <- s0
    hi[[k]] <- hi[[k]] + step
    lo[[k]] <- lo[[k]] - step
    (h_of(hi) - h_of(lo)) / (2 * step)
  }, numeric(length(tt)))

  v <- obj$vcov[c("THALF", "NU", "E0", "C0"), c("THALF", "NU", "E0", "C0")]
  want_se <- sqrt(rowSums((jac %*% v) * jac))
  # A zero here would mean the Jacobian collapsed and the comparison is empty.
  expect_true(all(want_se > 0))

  got <- predict(obj, newdata = data.frame(time = tt), type = "survival",
                 se.fit = TRUE, conf.type = "logit")
  expect_s3_class(got, "data.frame")
  expect_equal(got$fit, exp(-fixture_cumhaz(tt, est)))
  expect_equal(got$se.fit, want_se, tolerance = 1e-6)
  expect_true(all(got$lower < got$fit & got$fit < got$upper))
})

test_that("se.fit is refused rather than reported on the wrong scale", {
  # FIXMNU1 ties M to 1/NU, so the early block's covariance is not the
  # diagonal reparameterisation the mapping assumes.
  path <- outhaz_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "FIXMNU1"] <- 1
    d
  })
  obj <- hzr_read_outhaz(path)
  nd <- data.frame(time = c(1, 6))
  expect_error(predict(obj, newdata = nd, se.fit = TRUE), "FIXMNU1")
  # The point predictions are unaffected and must still work.
  expect_length(predict(obj, newdata = nd), 2L)
})

test_that("a librefs-resolved INHAZ job renders", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  file.copy(fixture_path(), file.path(dir, "hzdeath.rds"))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c(
    "DATA PGRID; DO MONTHS = 1 TO 12 BY 1; OUTPUT; END; RUN;",
    "%HAZPRED( PROC HAZPRED DATA=PGRID INHAZ=EX.HZDEATH OUT=P; TIME MONTHS; );"
  ), f)
  job <- suppressWarnings(
    hzr_translate_sas(f, librefs = c(EX = file.path(dir, "hzdeath.rds"))))
  got <- render_sim(job)
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
})
