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

test_that("a legacy fit keeps predicting plain row-wise phase terms", {
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
