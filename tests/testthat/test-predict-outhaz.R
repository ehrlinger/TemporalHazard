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

test_that("vcov is a matrix sized by the free parameters, not NULL", {
  # The fixture has four FREE parameters, so this asserts 4x4. It is not the
  # all-fixed case -- an all-fixed OUTHAZ would give a 0x0 matrix, and no
  # fixture here has one.
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

# ---------------------------------------------------------------------------
# The late phase
# ---------------------------------------------------------------------------
# outhaz-fixture.rds carries L0 _STATUS_ = 0, so nothing above ever reaches
# the g3 reconstruction block or the late-phase scale gate. This second
# fixture is an early + constant + late fit (see data-raw/outhaz_fixture.R):
# G3FLAG = 1, FIXGE2 = FIXGAE2 = 0, with TAU free and GAMMA/ALPHA/ETA fixed.

late_path <- function() {
  system.file("extdata", "outhaz-late-fixture.rds",
              package = "TemporalHazard")
}

late_variant <- function(edits) {
  d <- edits(readRDS(late_path()))
  path <- withr::local_tempfile(fileext = ".rds", .local_envir = parent.frame())
  saveRDS(d, path)
  path
}

# H(t) for the late fixture's three-phase model, from the stored estimates.
late_cumhaz <- function(time, est, early = TRUE, constant = TRUE,
                        late = TRUE) {
  h <- rep(0, length(time))
  if (early) {
    h <- h + exp(est[["E0"]]) *
      hzr_decompos(time, t_half = est[["THALF"]], nu = est[["NU"]],
                   m = est[["M"]])$G
  }
  if (constant) h <- h + exp(est[["C0"]]) * time
  if (late) {
    h <- h + exp(est[["L0"]]) *
      hzr_decompos_g3(time, tau = est[["TAU"]], gamma = est[["GAMMA"]],
                      alpha = est[["ALPHA"]], eta = est[["ETA"]])$G
  }
  h
}

test_that("the late fixture puts all three phases in the model", {
  obj <- hzr_read_outhaz(late_path())
  expect_equal(obj$status[c("E0", "C0", "L0")],
               c(E0 = 1L, C0 = 1L, L0 = 1L))
  spec <- .hzr_outhaz_to_spec(obj)
  expect_equal(names(spec$spec$phases), c("early", "constant", "late"))
  expect_equal(spec$spec$phases$late$type, "g3")
  # 1 early intercept + 3 early shapes + 1 constant + 1 late intercept
  # + 4 late shapes.
  expect_length(spec$fit$theta, 10L)
})

test_that("the reconstructed late phase is the one the file describes", {
  obj <- hzr_read_outhaz(late_path())
  tt <- c(1, 6, 12, 60)
  want <- late_cumhaz(tt, obj$estimates)

  expect_equal(predict(obj, newdata = data.frame(time = tt),
                       type = "cumulative_hazard"), want)
  expect_equal(predict(obj, newdata = data.frame(time = tt),
                       type = "survival"), exp(-want))

  # The late phase must actually be contributing: drop it and the numbers
  # have to move. Without this the assertions above pass on an
  # early + constant reconstruction that ignored L0 entirely.
  no_late <- late_cumhaz(tt, obj$estimates, late = FALSE)
  expect_true(max(abs(want - no_late)) > 1e-3)
})

test_that("se.fit maps the late block's covariance onto this scale", {
  # SAS's estimation variables here are (log THALF, log NU, log TAU, E0, C0,
  # L0) -- TAU is the one late parameter PROC HAZARD always estimates as
  # log(TAU). Differentiate H(t) in those coordinates and sandwich the raw
  # stored matrix; nothing here uses the reconstruction's own Jacobian.
  obj <- hzr_read_outhaz(late_path())
  est <- obj$estimates
  tt <- c(1, 6, 12, 60)
  free <- c("THALF", "NU", "TAU", "E0", "C0", "L0")

  h_of <- function(s) {
    e <- est
    e[["THALF"]] <- exp(s[[1L]])
    e[["NU"]] <- exp(s[[2L]])
    e[["TAU"]] <- exp(s[[3L]])
    e[["E0"]] <- s[[4L]]
    e[["C0"]] <- s[[5L]]
    e[["L0"]] <- s[[6L]]
    late_cumhaz(tt, e)
  }
  s0 <- c(log(est[["THALF"]]), log(est[["NU"]]), log(est[["TAU"]]),
          est[["E0"]], est[["C0"]], est[["L0"]])
  jac <- vapply(seq_along(s0), function(k) {
    step <- 1e-6 * max(1, abs(s0[[k]]))
    hi <- s0
    lo <- s0
    hi[[k]] <- hi[[k]] + step
    lo[[k]] <- lo[[k]] - step
    (h_of(hi) - h_of(lo)) / (2 * step)
  }, numeric(length(tt)))

  # The log(TAU) column must not be numerically dead, or this test would
  # pass on a reconstruction that dropped the late block from the variance.
  expect_true(abs(jac[length(tt), 3L]) > 1e-6)

  v <- obj$vcov[free, free]
  want_se <- sqrt(rowSums((jac %*% v) * jac))
  expect_true(all(want_se > 0))

  got <- predict(obj, newdata = data.frame(time = tt), type = "survival",
                 se.fit = TRUE, conf.type = "logit")
  expect_equal(got$fit, exp(-late_cumhaz(tt, est)))
  expect_equal(got$se.fit, want_se, tolerance = 1e-6)
})

test_that(".hzr_outhaz_late_composite() matches hzd_late_p2t.c", {
  # One row per branch of the three conditionals in
  # hazard/src/common/hzd_late_p2t.c. The expected values are read off the C
  # by hand, not recomputed from the R, so an inverted condition shows up
  # here rather than passing silently.
  ge_2 <- "log(GAMMA*ETA - 2)"
  gea_2 <- "log(GAMMA*ETA/ALPHA - 2)"
  cases <- list(
    # g_two, ga_two, g3flag, eta free, GAMMA, ETA, ALPHA
    list(0, 0, 1, 0L, ge_2, ge_2, gea_2),
    list(0, 0, 1, 1L, NA_character_, ge_2, gea_2),
    list(1, 0, 2, 0L, NA_character_, NA_character_, gea_2),
    list(1, 0, 2, 1L, NA_character_, NA_character_, gea_2),
    list(0, 1, 3, 0L, NA_character_, NA_character_, NA_character_),
    list(1, 1, 4, 0L, NA_character_, NA_character_, NA_character_),
    # setg3.c never writes g_two and ga_two together below G3FLAG = 4, but
    # hzd_late_p2t.c line 60 guards ALPHA with `ga_two && !g_two` all the
    # same. This row is the only one that isolates that guard.
    list(1, 1, 2, 0L, NA_character_, NA_character_, gea_2),
    list(0, 0, 5, 0L, NA_character_, NA_character_, NA_character_),
    list(1, 0, 6, 0L, NA_character_, NA_character_, NA_character_)
  )
  for (cs in cases) {
    got <- .hzr_outhaz_late_composite(
      status = c(ETA = cs[[4L]]),
      flags = c(FIXGE2 = cs[[1L]], FIXGAE2 = cs[[2L]], G3FLAG = cs[[3L]])
    )
    expect_equal(
      got,
      c(GAMMA = cs[[5L]], ETA = cs[[6L]], ALPHA = cs[[7L]]),
      info = paste("g_two =", cs[[1L]], "ga_two =", cs[[2L]],
                   "g3flag =", cs[[3L]], "eta free =", cs[[4L]])
    )
  }
})

test_that("an unconstrained free late shape is refused by name, not by flag", {
  # G3FLAG = 1 with FIXGE2 = FIXGAE2 = 0 is the GENERIC late phase -- no
  # constraint at all -- and PROC HAZARD still estimates ETA as
  # log(GAMMA*ETA - 2). The refusal has to say that, not claim a constraint
  # while printing zeros for every constraint flag.
  nd <- data.frame(time = c(1, 6, 12))
  for (par in c("GAMMA", "ALPHA", "ETA")) {
    path <- late_variant(function(d) {
      d[["_STATUS_"]][d[["_NAME_"]] == par] <- 1
      d
    })
    obj <- hzr_read_outhaz(path)
    err <- expect_error(predict(obj, newdata = nd, se.fit = TRUE))
    msg <- conditionMessage(err)
    expect_match(msg, par, fixed = TRUE)
    expect_match(msg, "GAMMA*ETA", fixed = TRUE)
    expect_false(grepl("under a constrained late-phase transformation", msg,
                       fixed = TRUE))
    # Point predictions are unaffected and must still work.
    expect_length(predict(obj, newdata = nd), 3L)
  }
})

test_that("freeing ETA moves GAMMA onto the plain log scale", {
  # hzd_late_p2t.c line 30: GAMMA takes log(GAMMA) when ETA is estimated.
  # So with both free it is ETA, and only ETA, that cannot be mapped. This
  # is the assertion that fails if the condition is inverted or dropped.
  path <- late_variant(function(d) {
    d[["_STATUS_"]][d[["_NAME_"]] %in% c("GAMMA", "ETA")] <- 1
    d
  })
  err <- expect_error(
    predict(hzr_read_outhaz(path), newdata = data.frame(time = 1),
            se.fit = TRUE))
  msg <- conditionMessage(err)
  expect_match(msg, "estimates ETA as log(GAMMA*ETA - 2)", fixed = TRUE)
  expect_false(grepl("GAMMA as ", msg, fixed = TRUE))
})

test_that("a late parameter derived under FIXGE2 is refused, not dropped", {
  # hzd_late_t2p.c lines 37-39: with FIXGE2 and GAMMA estimated,
  # ETA = 2/GAMMA. ETA carries _STATUS_ = 0 and so has no covariance row,
  # and d(ETA)/d(GAMMA) = -2/GAMMA^2 would be silently dropped from the
  # delta method -- a populated se.fit column that is quietly wrong.
  path <- late_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "FIXGE2"] <- 1
    d[["_EST_"]][d[["_NAME_"]] == "G3FLAG"] <- 2
    d[["_EST_"]][d[["_NAME_"]] == "GAMMA"] <- 2 / sqrt(3)
    d[["_STATUS_"]][d[["_NAME_"]] == "GAMMA"] <- 1
    d
  })
  obj <- hzr_read_outhaz(path)
  nd <- data.frame(time = c(1, 6))
  err <- expect_error(predict(obj, newdata = nd, se.fit = TRUE))
  expect_match(conditionMessage(err), "ETA = 2 / GAMMA", fixed = TRUE)
  expect_match(conditionMessage(err), "FIXGE2", fixed = TRUE)
  expect_length(predict(obj, newdata = nd), 2L)
})

test_that("a fixed ALPHA derived under FIXGAE2 is refused, not dropped", {
  # hzd_late_t2p.c lines 90-94: with FIXGAE2 and ALPHA fixed,
  # ALPHA = GAMMA*ETA/2, so ALPHA moves with the estimated GAMMA and ETA.
  path <- late_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "FIXGAE2"] <- 1
    d[["_EST_"]][d[["_NAME_"]] == "G3FLAG"] <- 3
    d[["_EST_"]][d[["_NAME_"]] == "ALPHA"] <- exp(1) * sqrt(3) / 2
    d[["_STATUS_"]][d[["_NAME_"]] %in% c("GAMMA", "ETA")] <- 1
    d
  })
  err <- expect_error(
    predict(hzr_read_outhaz(path), newdata = data.frame(time = 1),
            se.fit = TRUE))
  expect_match(conditionMessage(err), "ALPHA = GAMMA * ETA / 2", fixed = TRUE)
})

test_that("FIXGE2 with every late shape fixed still yields se.fit", {
  # The derived-parameter refusal must be about a derived parameter moving
  # with an ESTIMATED one, not about the flag being set. With no late shape
  # free there is nothing to propagate and the mapping is exact, so this
  # must not error -- otherwise the gate is just refusing on FIXGE2.
  path <- late_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "FIXGE2"] <- 1
    d[["_EST_"]][d[["_NAME_"]] == "G3FLAG"] <- 2
    d[["_EST_"]][d[["_NAME_"]] == "GAMMA"] <- 2 / sqrt(3)
    d
  })
  got <- predict(hzr_read_outhaz(path), newdata = data.frame(time = c(1, 6)),
                 se.fit = TRUE, type = "survival")
  expect_s3_class(got, "data.frame")
  expect_true(all(got$se.fit > 0))
})

test_that("a plain-log late shape IS mapped -- the gate is not blanket refusal", {
  # G3FLAG = 5 is setg3.c's EXPONENTIAL case (ALPHA <= 0, GETWO and GAETWO
  # both false). g3flag >= 3, so hzd_late_p2t.c estimates log(GAMMA) and
  # log(ETA) plainly and the diagonal Jacobian d(par)/d(theta) = par is
  # exact. This is the positive control for the whole late-phase gate: if it
  # errored, the refusals above would prove nothing, and if the Jacobian
  # factors were dropped the se.fit comparison below would move.
  free <- c("THALF", "NU", "TAU", "GAMMA", "ETA", "E0", "C0", "L0")
  lt <- matrix(0, 8, 8, dimnames = list(free, free))
  diag(lt) <- 1 / sqrt(c(2, 3, 5, 7, 11, 13, 17, 19))
  lt[lower.tri(lt)] <- 1 / seq(23, by = 6, length.out = 28)
  vc <- lt %*% t(lt)

  path <- late_variant(function(d) {
    d[["_EST_"]][d[["_NAME_"]] == "G3FLAG"] <- 5
    d[["_EST_"]][d[["_NAME_"]] == "ALPHA"] <- 0
    d[["_STATUS_"]][d[["_NAME_"]] %in% c("GAMMA", "ETA")] <- 1
    pr <- which(!is.na(d[["_STATUS_"]]))
    nmp <- d[["_NAME_"]][pr]
    for (a in free) {
      for (b in free) d[[b]][pr][match(a, nmp)] <- vc[a, b]
    }
    d
  })
  obj <- hzr_read_outhaz(path)
  est <- obj$estimates
  tt <- c(1, 6, 12, 60)

  h_of <- function(s) {
    e <- est
    e[["THALF"]] <- exp(s[[1L]])
    e[["NU"]] <- exp(s[[2L]])
    e[["TAU"]] <- exp(s[[3L]])
    e[["GAMMA"]] <- exp(s[[4L]])
    e[["ETA"]] <- exp(s[[5L]])
    e[["E0"]] <- s[[6L]]
    e[["C0"]] <- s[[7L]]
    e[["L0"]] <- s[[8L]]
    late_cumhaz(tt, e)
  }
  s0 <- c(log(est[["THALF"]]), log(est[["NU"]]), log(est[["TAU"]]),
          log(est[["GAMMA"]]), log(est[["ETA"]]),
          est[["E0"]], est[["C0"]], est[["L0"]])
  jac <- vapply(seq_along(s0), function(k) {
    step <- 1e-6 * max(1, abs(s0[[k]]))
    hi <- s0
    lo <- s0
    hi[[k]] <- hi[[k]] + step
    lo[[k]] <- lo[[k]] - step
    (h_of(hi) - h_of(lo)) / (2 * step)
  }, numeric(length(tt)))

  # The two shape columns must carry real signal, or dropping their Jacobian
  # factor would be invisible here.
  expect_true(all(abs(jac[length(tt), 4:5]) > 1e-3))

  want_se <- sqrt(rowSums((jac %*% vc[free, free]) * jac))
  got <- predict(obj, newdata = data.frame(time = tt), type = "survival",
                 se.fit = TRUE, conf.type = "logit")
  expect_equal(got$fit, exp(-late_cumhaz(tt, est)))
  expect_equal(got$se.fit, want_se, tolerance = 1e-6)
})

# predict.hzr_outhaz() and predict.hazard() are two methods of one generic in
# one package. When their argument lists disagreed, positional slot 4 was
# `decompose` on one and `se.fit` on the other, and `level`/`conf.type` fell
# into `...` here -- so a misspelt `conf.type` silently returned log-log
# limits where the SAS job being reproduced asked for logit ones.

test_that("predict.hzr_outhaz() takes predict.hazard()'s arguments in order", {
  fo <- formals(predict.hzr_outhaz)
  fh <- formals(predict.hazard)
  expect_equal(names(fo),
               c("object", "newdata", "type", "decompose", "se.fit",
                 "level", "conf.type", "..."))
  expect_equal(names(fo), names(fh))
  # Same default type, so an unnamed `type` means the same thing in both.
  expect_equal(eval(fo$type)[[1L]], eval(fh$type)[[1L]])
  expect_equal(fo$level, fh$level)
  expect_equal(eval(fo$conf.type), eval(fh$conf.type))
})

test_that("a positional predict() call means the same for both methods", {
  obj <- hzr_read_outhaz(fixture_path())
  nd <- data.frame(time = c(1, 6, 12))
  # Slot 4 is `decompose`, slot 5 is `se.fit` -- as in predict.hazard().
  got <- predict(obj, nd, "survival", FALSE, TRUE)
  want <- predict(obj, newdata = nd, type = "survival", se.fit = TRUE)
  expect_equal(got, want)
  expect_s3_class(got, "data.frame")
  # Not a positional call that quietly did nothing: it really used se.fit.
  expect_true(all(got$se.fit > 0))
})

test_that("decompose = TRUE is refused, not silently ignored", {
  obj <- hzr_read_outhaz(fixture_path())
  nd <- data.frame(time = c(1, 6, 12))
  err <- expect_error(
    predict(obj, newdata = nd, type = "cumulative_hazard", decompose = TRUE)
  )
  expect_match(conditionMessage(err), "OUTHAZ=", fixed = TRUE)
  expect_match(conditionMessage(err), "decompose", fixed = TRUE)
  # The undecomposed prediction is what would have come back instead.
  expect_length(
    predict(obj, newdata = nd, type = "cumulative_hazard", decompose = FALSE),
    3L
  )
})

test_that("a misspelt conf.type errors instead of changing the answer", {
  obj <- hzr_read_outhaz(fixture_path())
  nd <- data.frame(time = c(1, 6, 12))
  err <- expect_error(
    predict(obj, newdata = nd, type = "survival", se.fit = TRUE,
            conf_type = "logit")
  )
  expect_match(conditionMessage(err), "conf_type", fixed = TRUE)
  # The two transforms really do give different limits, so swallowing the
  # typo would have been a wrong answer, not a harmless one.
  logit <- predict(obj, newdata = nd, type = "survival", se.fit = TRUE,
                   conf.type = "logit")
  loglog <- predict(obj, newdata = nd, type = "survival", se.fit = TRUE,
                    conf.type = "log-log")
  expect_false(isTRUE(all.equal(logit$lower, loglog$lower)))
})

test_that("conf.type is validated where it is used, not before", {
  # predict.hazard() defers this match.arg() to the survival-SE path, so an
  # otherwise-ignored value must not error here either: the two methods of
  # one generic disagreeing about which calls are legal is a surprise the
  # bilingual reader pays for. Rejecting an unknown argument NAME is the
  # separate half, pinned in the test above and again below.
  obj <- hzr_read_outhaz(fixture_path())
  nd <- data.frame(time = c(1, 6, 12))

  ref_haz <- predict(obj, newdata = nd, type = "hazard")
  expect_equal(predict(obj, newdata = nd, type = "hazard",
                       conf.type = "unused"),
               ref_haz)
  ref_pt <- predict(obj, newdata = nd, type = "survival")
  expect_equal(predict(obj, newdata = nd, type = "survival",
                       conf.type = "unused"),
               ref_pt)
  # Not vacuous: the reference values must be real predictions, not a
  # degenerate constant that any two calls would agree on.
  expect_true(all(is.finite(ref_haz)))
  expect_gt(length(unique(ref_pt)), 1L)

  # On the one path that reads it, a bad value still errors.
  expect_error(
    predict(obj, newdata = nd, type = "survival", se.fit = TRUE,
            conf.type = "unused")
  )
  # And an unknown argument NAME errors regardless of the type, rather than
  # vanishing into `...`.
  err <- expect_error(
    predict(obj, newdata = nd, type = "hazard", conf_type = "logit")
  )
  expect_match(conditionMessage(err), "conf_type", fixed = TRUE)
})
