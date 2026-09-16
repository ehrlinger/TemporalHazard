# predict(newdata = ) on a multiphase fit saved before the phase design was
# stored (no fit$x_design), whose phase formula has a data-dependent term
# (#307).
#
# scale(), poly() and ns() take their centering or basis from the fitting data.
# Rebuilding such a phase from `newdata` alone recomputed them from newdata's
# own rows: silently wrong for any newdata but the fitting rows, and a
# zero-length prediction for one row of a scale() phase.
#
# Reference: the fit's own stored design. predict() without newdata evaluates
# fit$x_list at the stored times, so its values at a subset of the fitting
# rows are what predict(newdata = those rows) must return. The subset matters:
# at all the fitting rows, recomputing from newdata reproduces the fit, and a
# test there cannot fail.

legacy_fit <- function(term, keep_frame) {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = stats::as.formula(paste("~", term))),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  )
  # A fit saved by 1.1.0 to 1.2.10 kept its data but no phase design; one
  # saved by 1.0.3 or earlier kept neither.
  fit$fit$x_design <- NULL
  if (!keep_frame) fit$data$frame <- NULL
  list(fit = fit, data = d)
}

at_rows <- function(d, rows) {
  data.frame(time = d$int_dead[rows], age = d$age[rows])
}

test_that("a legacy fit that kept its data rebuilds data-dependent phase terms as fitted", {
  skip_on_cran()  # multiphase fits
  for (term in c("scale(age)", "poly(age, 2)", "splines::ns(age, df = 3)")) {
    lf <- legacy_fit(term, keep_frame = TRUE)
    fitted <- predict(lf$fit, type = "cumulative_hazard")
    for (rows in list(c(1, 50, 200), 50L)) {
      got <- predict(lf$fit, newdata = at_rows(lf$data, rows),
                     type = "cumulative_hazard")
      # Length first: one row of a scale() phase came back empty.
      expect_length(got, length(rows))
      expect_equal(got / fitted[rows], rep(1, length(rows)),
                   tolerance = 1e-8, ignore_attr = TRUE,
                   label = paste(term, "at rows", paste(rows, collapse = ",")))
    }
  }
})

test_that("a legacy fit without its data refuses data-dependent phase terms", {
  skip_on_cran()  # multiphase fits
  for (term in c("scale(age)", "poly(age, 2)", "splines::ns(age, df = 3)")) {
    lf <- legacy_fit(term, keep_frame = FALSE)
    expect_error(
      predict(lf$fit, newdata = at_rows(lf$data, c(1, 50, 200)),
              type = "cumulative_hazard"),
      "refit",
      label = term
    )
  }
})

test_that("a legacy fit whose kept data no longer reproduces it is not trusted", {
  skip_on_cran()  # multiphase fits
  # The kept data is used only when it rebuilds the fitted columns. Altered
  # after the fit (one row corrected, say), it would give scale() the wrong
  # centering, silently; the term is refused as for a fit without its data
  # instead. (A rescaling of `age` would not be seen: scale() of 2 * age is
  # scale() of age, so the rebuilt columns still match.)
  lf <- legacy_fit("scale(age)", keep_frame = TRUE)
  lf$fit$data$frame$age[1] <- lf$fit$data$frame$age[1] + 50
  expect_error(
    predict(lf$fit, newdata = at_rows(lf$data, c(1, 50, 200)),
            type = "cumulative_hazard"),
    "refit"
  )
})

legacy_fit_on <- function(term, d, keep_frame, env = parent.frame()) {
  # `env` is where the formula is made, as a user's formula would be, so a
  # constant such as `cutoff` resolves from the calling test.
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = stats::as.formula(paste("~", term),
                                                    env = env)),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  current <- fit
  fit$fit$x_design <- NULL
  if (!keep_frame) fit$data$frame <- NULL
  list(fit = fit, current = current)
}

rows_of <- function(d, rows) {
  data.frame(time = d$int_dead[rows], age = d$age[rows],
             opmos = d$opmos[rows], inc_surg = d$inc_surg[rows])
}

test_that("a data-dependent term is refused whatever newdata's row order", {
  skip_on_cran()  # multiphase fits
  # A probe row moves one column past its range; for opmos that is beyond
  # cut()'s last break, an NA cell, and scale() must still be caught by the
  # probes on age.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("scale(age) + cut(opmos, c(0, 50, 100, 200))", d,
                      keep_frame = FALSE)
  i <- which(d$opmos > 100 & d$opmos * 2 + 1 > 200)[1]
  for (rows in list(c(i, 50, 200), c(200, 50, i))) {
    expect_error(
      predict(lf$fit, newdata = rows_of(d, rows), type = "cumulative_hazard"),
      "refit"
    )
  }
})

test_that("the range probes and the copy probe each catch what the other cannot", {
  skip_on_cran()  # multiphase fits
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  # One row: a copy of it gives two identical rows, so scale() is NaN either
  # way; the range probes set age one below and one above its single value.
  lf1 <- legacy_fit_on("scale(age)", d, keep_frame = FALSE)
  expect_error(
    predict(lf1$fit, newdata = rows_of(d, 50L), type = "cumulative_hazard"),
    "refit"
  )
  # The same column feeds factor() and scale(): a range probe on inc_surg is
  # a new factor level, so both are skipped, and only the copy of row 1 moves
  # scale()'s centering.
  lf2 <- legacy_fit_on("factor(inc_surg) + scale(inc_surg)", d,
                       keep_frame = FALSE)
  lv <- sort(unique(d$inc_surg), decreasing = TRUE)
  rows <- vapply(lv, function(v) which(d$inc_surg == v)[1], integer(1))
  expect_error(
    predict(lf2$fit, newdata = rows_of(d, rows), type = "cumulative_hazard"),
    "refit"
  )
})

adversarial_nd <- data.frame(time = c(1, 2, 3), age = c(10, 20, 30),
                             opmos = c(30, 60, 90), inc_surg = c(1, 2, 3))

test_that("a data-dependent term is refused however it reads the data", {
  skip_on_cran()  # multiphase fits
  # Each reads the data a different way: its minimum, its maximum, its
  # median, bins over its range, its mean alone, and boundary knots from its
  # range. A probe past only one end, or a copy of a row, misses some.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  for (term in c("I(age - min(age))", "I(age / max(age))",
                 "I(age > median(age))", "cut(age, 3)",
                 "scale(age, scale = FALSE)",
                 "splines::ns(age, knots = numeric(0))",
                 "splines::bs(age, df = 3)")) {
    lf <- legacy_fit_on(term, d, keep_frame = FALSE)
    expect_error(
      predict(lf$fit, newdata = adversarial_nd, type = "cumulative_hazard"),
      "refit",
      label = term
    )
  }
})

test_that("row-wise and fully specified phase terms still predict as fitted", {
  skip_on_cran()  # multiphase fits
  # These take nothing from the data they are built on, so a current fit's
  # prediction is the reference.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  cutoff <- 100
  for (term in c("log(age)", "I(age > cutoff)",
                 "cut(opmos, c(0, 50, 100, 200))",
                 "poly(age, 2, raw = TRUE)",
                 "splines::ns(age, knots = 50, Boundary.knots = c(0, 400))")) {
    lf <- legacy_fit_on(term, d, keep_frame = FALSE)
    want <- predict(lf$current, newdata = adversarial_nd,
                    type = "cumulative_hazard")
    got <- predict(lf$fit, newdata = adversarial_nd,
                   type = "cumulative_hazard")
    expect_length(got, 3L)
    expect_equal(got / want, rep(1, 3), tolerance = 1e-8, ignore_attr = TRUE,
                 label = term)
  }
})

test_that("a column with one distinct value is still probed", {
  skip_on_cran()  # multiphase fits
  # The range probes set it one below and one above its single value.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  nd <- adversarial_nd
  nd$age <- 20
  lf <- legacy_fit_on("scale(age)", d, keep_frame = FALSE)
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
  lf <- legacy_fit_on("log(age)", d, keep_frame = FALSE)
  want <- predict(lf$current, newdata = nd, type = "cumulative_hazard")
  expect_equal(predict(lf$fit, newdata = nd, type = "cumulative_hazard") / want,
               rep(1, 3), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("a term that cannot be built from newdata alone gets the refit advice", {
  skip_on_cran()  # multiphase fits
  # poly() needs more unique points than its degree, so one row cannot be
  # built at all; the refusal says why rather than the transform's own error.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("poly(age, 2)", d, keep_frame = FALSE)
  expect_error(
    predict(lf$fit, newdata = rows_of(d, 50L), type = "cumulative_hazard"),
    "refit"
  )
})

test_that("a legacy phase whose factor columns come from newdata is refused", {
  skip_on_cran()  # multiphase fits
  # cut(age, 3) takes its bins from newdata's range, so the rebuilt columns
  # are not the fitted ones.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("cut(age, 3)", d, keep_frame = FALSE)
  expect_error(
    predict(lf$fit, newdata = rows_of(d, c(1, 50, 200)),
            type = "cumulative_hazard"),
    "refit"
  )
})

test_that("a legacy factor() phase predicts with all its levels and refuses without", {
  skip_on_cran()  # multiphase fits
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("factor(inc_surg)", d, keep_frame = FALSE)
  fitted <- predict(lf$current, type = "cumulative_hazard")
  lv <- sort(unique(d$inc_surg), decreasing = TRUE)
  # The range probes on inc_surg are new levels, so they are skipped.
  rows <- vapply(lv, function(v) which(d$inc_surg == v)[1], integer(1))
  got <- predict(lf$fit, newdata = rows_of(d, rows), type = "cumulative_hazard")
  expect_length(got, length(rows))
  expect_equal(got / fitted[rows], rep(1, length(rows)),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_error(
    predict(lf$fit, newdata = rows_of(d, rows[1:3]),
            type = "cumulative_hazard"),
    "refit"
  )
})

test_that("a missing value in newdata gives an NA row for a legacy fit", {
  skip_on_cran()  # multiphase fits
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("log(age)", d, keep_frame = FALSE)
  fitted <- predict(lf$current, type = "cumulative_hazard")
  nd <- rows_of(d, c(1, 50, 200))
  nd$age[2] <- NA
  got <- predict(lf$fit, newdata = nd, type = "cumulative_hazard")
  expect_length(got, 3L)
  expect_true(is.na(got[2]))
  expect_equal(got[c(1, 3)] / fitted[c(1, 200)], c(1, 1),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("a legacy fit that dropped a missing row is checked against its data", {
  skip_on_cran()  # multiphase fits
  # The fit dropped the NA row; its kept data still has it. Untouched, the
  # rebuilt design matches the fit; altered afterwards, it is not trusted.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  d$age[5] <- NA
  lf <- legacy_fit_on("scale(age)", d, keep_frame = TRUE)
  nd <- rows_of(d, c(1, 50, 200))
  want <- predict(lf$current, newdata = nd, type = "cumulative_hazard")
  expect_equal(predict(lf$fit, newdata = nd, type = "cumulative_hazard") / want,
               rep(1, 3), tolerance = 1e-8, ignore_attr = TRUE)
  lf$fit$data$frame$age[1] <- lf$fit$data$frame$age[1] + 50
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
})

test_that("a legacy fit keeps predicting plain row-wise phase terms", {
  skip_on_cran()  # multiphase fits
  for (keep in c(TRUE, FALSE)) {
    lf <- legacy_fit("log(age)", keep_frame = keep)
    fitted <- predict(lf$fit, type = "cumulative_hazard")
    rows <- c(1, 50, 200)
    got <- predict(lf$fit, newdata = at_rows(lf$data, rows),
                   type = "cumulative_hazard")
    expect_length(got, length(rows))
    expect_equal(got / fitted[rows], rep(1, length(rows)),
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})
