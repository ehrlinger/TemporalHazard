# predict(newdata = ) on a multiphase fit saved before the phase design was
# stored (no fit$x_design), whose phase formula has a data-dependent term
# (#307).
#
# scale(), poly() and ns() take their centering or basis from the fitting data.
# Rebuilding such a phase from `newdata` alone recomputed them from newdata's
# own rows: silently wrong for any newdata but the fitting rows, and a
# zero-length prediction for one row of a scale() phase.
#
# A fit that kept its fitting data (data$frame, 1.1.0 onward) is rebuilt from
# it. One that kept neither (1.0.3 or earlier) is rebuilt only from a closed
# formula: data columns and known elementwise operations. Anything else is
# refused, since nothing is left to check the rebuild against.
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

legacy_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  d$male <- d$mal == 1
  d$opdate <- as.Date("1970-01-01") + round(d$opmos * 30)
  d$grp <- c("A", "B", "C")[(d$inc_surg %% 3) + 1]
  d$ni <- d$inc_surg
  d
}

legacy_fit_on <- function(term, d, keep_frame, env = parent.frame()) {
  # `env` is where the formula is made, as a user's formula would be, so a
  # constant such as `cutoff` resolves from the calling test.
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = if (inherits(term, "formula")) term else
                          stats::as.formula(paste("~", term), env = env)),
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

# Six rows covering every inc_surg / ni level (0 to 5).
legacy_nd <- function() {
  data.frame(time = 1:6, age = c(10, 20, 30, 40, 50, 60),
             opmos = c(30, 60, 90, 120, 150, 180),
             inc_surg = 0:5, ni = 0:5,
             male = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE),
             opdate = as.Date(c("1975-01-01", "1976-01-01", "1977-01-01",
                                "1978-01-01", "1979-01-01", "1980-01-01")),
             grp = c("B", "C", "C", "B", "C", "B"))
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
  # centering, silently; the fit is treated as having no data instead.
  # (A rescaling of `age` would not be seen: scale() of 2 * age is scale() of
  # age, so the rebuilt columns still match.)
  lf <- legacy_fit("scale(age)", keep_frame = TRUE)
  lf$fit$data$frame$age[1] <- lf$fit$data$frame$age[1] + 50
  expect_error(
    predict(lf$fit, newdata = at_rows(lf$data, c(1, 50, 200)),
            type = "cumulative_hazard"),
    "refit"
  )
})

test_that("a data-dependent term is refused whatever newdata's row order", {
  skip_on_cran()  # multiphase fits
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

test_that("without its data, a phase term that is not closed is refused", {
  skip_on_cran()  # multiphase fits
  # Each reads the data it is built on, or something outside the formula:
  # its minimum, maximum, median or mean; bins over its range; boundary
  # knots from its range; factor codes or a date origin inside I(); a
  # constant from the formula's environment; or a basis that, though fixed
  # here, is not an elementwise operation. With no fitting data, none can
  # be checked, so all are refused.
  d <- legacy_data()
  nd <- legacy_nd()
  cutoff <- 100
  # A user's own sqrt(), which reads the data, is not base R's.
  sqrt <- function(x) x - mean(x)
  for (term in c("I(age - min(age))", "I(age / max(age))",
                 "I(age > median(age))", "I(age > mean(age))",
                 "cut(age, 3)", "scale(age, scale = FALSE)",
                 "splines::ns(age, knots = numeric(0))",
                 "splines::bs(age, df = 3)",
                 "I(male - mean(male))",
                 "I(as.numeric(opdate - min(opdate)))",
                 "I(as.numeric(factor(grp)))",
                 "factor(ni) + I(ni - mean(ni))",
                 "I(age > cutoff)", "factor(age > mean(age))", "sqrt(age)",
                 "cut(opmos, c(0, cutoff, 1000))", "cut(age, 3, labels = FALSE)",
                 "cut(age, c(3), labels = FALSE)",
                 "poly(age, 2, raw = TRUE)",
                 "splines::ns(age, knots = 50, Boundary.knots = c(0, 400))")) {
    lf <- legacy_fit_on(term, d, keep_frame = FALSE)
    expect_error(
      predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
      "refit",
      label = term
    )
  }
})

test_that("without its data, a closed phase formula still predicts as fitted", {
  skip_on_cran()  # multiphase fits
  # Data columns and elementwise operations take nothing from the rows they
  # are built on, so a current fit's prediction is the reference.
  d <- legacy_data()
  nd <- legacy_nd()
  for (term in c("log(age)", "I(2 * age + 1)", "I(age > 50)",
                 "sqrt(age) + exp(-opmos / 100)",
                 "cut(opmos, c(0, 50, 100, 200))")) {
    lf <- legacy_fit_on(term, d, keep_frame = FALSE)
    want <- predict(lf$current, newdata = nd, type = "cumulative_hazard")
    got <- predict(lf$fit, newdata = nd, type = "cumulative_hazard")
    expect_length(got, nrow(nd))
    expect_equal(got / want, rep(1, nrow(nd)), tolerance = 1e-8,
                 ignore_attr = TRUE, label = term)
  }
})

test_that("a formula that does not survive deparse is refused", {
  skip_on_cran()  # multiphase fits
  # A Date inlined into the formula passes the symbol and function checks,
  # but is not a literal the formula's text can carry.
  d <- legacy_data()
  f <- eval(bquote(~ I(opdate > .(as.Date("1975-06-01")))))
  lf <- legacy_fit_on(f, d, keep_frame = FALSE)
  expect_error(
    predict(lf$fit, newdata = legacy_nd(), type = "cumulative_hazard"),
    "refit"
  )
})

test_that("a column with one distinct value is refused or predicted by the formula alone", {
  skip_on_cran()  # multiphase fits
  d <- legacy_data()
  nd <- legacy_nd()
  nd$age <- 20
  lf <- legacy_fit_on("scale(age)", d, keep_frame = FALSE)
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
  lf <- legacy_fit_on("log(age)", d, keep_frame = FALSE)
  want <- predict(lf$current, newdata = nd, type = "cumulative_hazard")
  expect_equal(predict(lf$fit, newdata = nd, type = "cumulative_hazard") / want,
               rep(1, nrow(nd)), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("a closed term that cannot be built from newdata alone gets the refit advice", {
  skip_on_cran()  # multiphase fits
  # factor() of one row has one level, so its contrasts cannot be built; the
  # refusal says why rather than the model matrix's own error.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("factor(inc_surg)", d, keep_frame = FALSE)
  expect_error(
    predict(lf$fit, newdata = rows_of(d, 50L), type = "cumulative_hazard"),
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

test_that("without its data, a column shadowing a constant or a coding the names do not carry is refused", {
  skip_on_cran()  # multiphase fits
  # Each passes the column-name check with different values.
  d <- legacy_data()
  d$og <- ordered(d$grp, levels = c("A", "B", "C"))
  nd <- legacy_nd()
  # T and pi were base R's in the fit; newdata's columns of those names
  # would take their place.
  nd$T <- FALSE
  nd$pi <- 1
  # Every level present, so a missing level is not what refuses.
  nd$grp <- c("A", "B", "C", "A", "B", "C")
  for (term in c("I(age > 50 & T)", "I(age * pi)")) {
    lf <- legacy_fit_on(term, d, keep_frame = FALSE)
    expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
                 "refit", label = term)
  }
  # An ordered factor codes by polynomial contrasts, whose columns (.L, .Q)
  # do not name the levels, so reversed levels would pass.
  nd$og <- factor(nd$grp, levels = c("C", "B", "A"), ordered = TRUE)
  lf <- legacy_fit_on("factor(og)", d, keep_frame = FALSE)
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
  # A fit made under other contrasts: Helmert and sum coding both name the
  # columns 1 and 2, and a legacy fit does not record which it used.
  old <- options(contrasts = c("contr.helmert", "contr.poly"))
  on.exit(options(old), add = TRUE)
  lf <- legacy_fit_on("factor(grp)", d, keep_frame = FALSE)
  options(contrasts = c("contr.sum", "contr.poly"))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
})
