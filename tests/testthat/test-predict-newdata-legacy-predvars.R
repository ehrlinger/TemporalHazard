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
# it, when that can be checked: a formula closed over the kept columns, no
# term coded by contrasts, and a rebuild equal to the fitted columns. Anything
# else, and every fit that kept neither design nor data (1.0.3 or earlier), is
# refused at newdata: such a fit cannot say which of its formula's names were
# data columns and which were constants.
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
  d$grp <- c("A", "B", "C")[(d$inc_surg %% 3) + 1]
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

# The gap between two fitted ages, just above 100.
age_gap <- function(d) {
  ages <- sort(unique(d$age))
  i <- which(ages > 100)[1L]
  c(ages[i - 1L], ages[i])
}

test_that("a legacy fit that kept its data rebuilds data-dependent phase terms as fitted", {
  skip_on_cran()  # multiphase fits
  for (term in c("scale(age)", "poly(age, 2)", "splines::ns(age, df = 3)",
                 "log(age)")) {
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

test_that("a missing value in newdata gives an NA row for a legacy fit", {
  skip_on_cran()  # multiphase fits
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  lf <- legacy_fit_on("log(age)", d, keep_frame = TRUE)
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

test_that("kept data changed after the fit is not trusted", {
  skip_on_cran()  # multiphase fits
  # Altered after the fit (one row corrected, say), the kept data would give
  # scale() another centering, silently. The comparison is exact: an
  # approximate one let a change to one age of 1.8e-4 through. (A rescaling
  # of `age` would not be seen: scale() of 2 * age is scale() of age.)
  for (change in c(50, 1e-4)) {
    lf <- legacy_fit("scale(age)", keep_frame = TRUE)
    lf$fit$data$frame$age[1] <- lf$fit$data$frame$age[1] + change
    expect_error(
      predict(lf$fit, newdata = at_rows(lf$data, c(1, 50, 200)),
              type = "cumulative_hazard"),
      "refit", label = paste("change", change)
    )
  }
})

test_that("kept data is not trusted when the formula reads outside it", {
  skip_on_cran()  # multiphase fits
  # The kept data proves only the fitting rows. A cutoff moved within a gap
  # between fitted ages leaves every fitted row the same, so the rebuild
  # matches the fit and still scores a row in the gap with the new cutoff.
  d <- legacy_data()
  gap <- age_gap(d)
  e <- new.env()
  e$cutoff <- gap[1] + 0.25 * diff(gap)
  lf <- legacy_fit_on(stats::as.formula("~ I(age * (age > cutoff))", env = e),
                      d, keep_frame = TRUE)
  e$cutoff <- gap[1] + 0.75 * diff(gap)
  expect_identical(unname(lf$fit$fit$x_list$early[, 1L]),
                   d$age * (d$age > e$cutoff))
  nd <- data.frame(time = 2, age = mean(gap))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
  # A user's log() whose state changes only between two fitted ages: every
  # fitted row, and so the rebuilt design, is unchanged.
  e2 <- new.env()
  e2$k <- gap[2]
  e2$hi <- gap[2]
  e2$log <- function(x) base::log(x) + (x > k & x < hi)
  environment(e2$log) <- e2
  lf <- legacy_fit_on(stats::as.formula("~ log(age)", env = e2), d,
                      keep_frame = TRUE)
  e2$k <- gap[1]
  expect_identical(unname(lf$fit$fit$x_list$early[, 1L]), e2$log(d$age))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit")
  # A function not on the list, the user's own or R's, qualified or not: its
  # rebuild from the kept data matches the fit, but it may read state, or
  # newdata's rows, at predict time.
  e3 <- new.env()
  e3$k <- gap[1] + 0.25 * diff(gap)
  e3$thr <- function(x) x * (x > k)
  environment(e3$thr) <- e3
  lf <- legacy_fit_on(stats::as.formula("~ I(thr(age))", env = e3), d,
                      keep_frame = TRUE)
  e3$k <- gap[1] + 0.75 * diff(gap)
  expect_identical(unname(lf$fit$fit$x_list$early[, 1L]), e3$thr(d$age))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit", label = "user's thr()")
  for (term in c("I(age - mean(age))", "I(age - stats::median(age))")) {
    lf <- legacy_fit_on(term, d, keep_frame = TRUE)
    expect_error(predict(lf$fit, newdata = at_rows(d, c(1, 50, 200)),
                         type = "cumulative_hazard"),
                 "refit", label = term)
  }
})

test_that("kept data is not trusted for a formula its text does not carry", {
  skip_on_cran()  # multiphase fits
  d <- legacy_data()
  nd <- at_rows(d, c(1, 50, 200))
  # -0 prints as 0; a classed constant prints as a call, not as itself.
  for (f in list(eval(bquote(~ I(age + exp(1 / .(-0))))),
                 eval(bquote(~ I(age * .(structure(2, class = "myunit"))))))) {
    lf <- legacy_fit_on(f, d, keep_frame = TRUE)
    expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
                 "refit", label = deparse(f))
  }
})

test_that("kept data is not trusted for a term coded by contrasts", {
  skip_on_cran()  # multiphase fits
  # The fit did not record its contrasts. A level whose rows a term
  # multiplies by 0 (g in x:g) has its code shown by no fitted value, so
  # another contrasts option, or a contrasts function that names its
  # columns as treatment coding does, would recode it silently; so would a
  # level order moved in rows the fit dropped, under numbered contrasts.
  d <- legacy_data()
  lf <- legacy_fit_on("grp", d, keep_frame = TRUE)
  nd <- data.frame(time = 1:3, grp = c("A", "B", "C"))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit", label = "grp")

  d$x <- d$age / 100
  d$x[d$opmos > 100] <- 0
  recode <- function(n, contrasts = TRUE, sparse = FALSE) {
    m <- stats::contr.treatment(n)
    m[3, ] <- c(1, 1)
    m
  }
  assign("contr.recode307", recode, envir = globalenv())
  on.exit(rm("contr.recode307", envir = globalenv()), add = TRUE)
  old <- options(contrasts = c("contr.recode307", "contr.poly"))
  on.exit(options(old), add = TRUE)
  lf <- legacy_fit_on("x + x:cut(opmos, c(0, 50, 100, 1000))", d,
                      keep_frame = TRUE)
  options(contrasts = c("contr.treatment", "contr.poly"))
  nd <- data.frame(time = 1, x = 2, opmos = c(30, 75, 150))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit", label = "recoding contrasts function")

  set.seed(1)
  d$z <- d$age
  d$g <- sample(c("c", "d", "e"), nrow(d), replace = TRUE)
  d$g[1:4] <- c("a", "a", "b", "b")
  d$z[1:4] <- NA
  d$g <- factor(d$g, levels = c("a", "b", "c", "d", "e"))
  options(contrasts = c("contr.sum", "contr.poly"))
  lf <- legacy_fit_on("g + z", d, keep_frame = TRUE)
  nd <- data.frame(time = 5, g = c("a", "b", "c"), z = 100)
  want <- predict(lf$current, newdata = nd, type = "cumulative_hazard")
  expect_false(isTRUE(all.equal(want[1], want[2])))
  lf$fit$data$frame$g <- factor(as.character(d$g),
                                levels = c("b", "a", "c", "d", "e"))
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit", label = "levels moved in dropped rows")
})

test_that("a legacy fit with duplicated phase column names is refused", {
  skip_on_cran()  # multiphase fits
  # Fits made before duplicated names were refused (#296) can hold a factor
  # dummy `gb` beside a numeric `gb`. Selecting the columns by name takes the
  # first `gb` twice.
  d <- legacy_data()
  d$g <- factor(c("a", "b")[(d$inc_surg %% 2) + 1])
  d$gb <- d$age / 100
  local_mocked_bindings(.hzr_refuse_duplicate_columns =
                          function(x, phase = NULL) invisible(NULL))
  lf <- legacy_fit_on("g + gb", d, keep_frame = TRUE)
  expect_identical(colnames(lf$current$fit$x_list$early), c("gb", "gb"))
  rows <- c(1, 2, 3, 50)
  nd <- data.frame(time = d$int_dead[rows], g = d$g[rows], gb = d$gb[rows])
  none <- lf$fit
  none$data$frame <- NULL
  for (f in list(lf$current, lf$fit, none)) {
    expect_error(predict(f, newdata = nd, type = "cumulative_hazard"),
                 "duplicated")
  }
})

test_that("a legacy fit without its data refuses to rebuild any phase formula", {
  skip_on_cran()  # multiphase fits
  # Such a fit cannot say which of its formula's names were data columns: a
  # constant `k` that is gone at predict time (rm(), a new session) would be
  # taken from a newdata column of that name. Every rebuild is refused.
  d <- legacy_data()
  rows <- c(1, 50, 200)
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows],
                   opmos = d$opmos[rows], male = d$male[rows])
  for (term in c("log(age)", "I(2 * age + 1)", "scale(age)", "poly(age, 2)",
                 "splines::ns(age, df = 3)", "male",
                 "cut(opmos, c(0, 50, 100, 1000))")) {
    lf <- legacy_fit_on(term, d, keep_frame = FALSE)
    expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
                 "refit", label = term)
  }
  e <- new.env()
  e$k <- 0.01
  lf <- legacy_fit_on(stats::as.formula("~ I(age * k)", env = e), d,
                      keep_frame = FALSE)
  rm("k", envir = e)
  nd$k <- 0.5
  expect_error(predict(lf$fit, newdata = nd, type = "cumulative_hazard"),
               "refit", label = "constant gone, newdata column of its name")

  # Still predicts without newdata, and from the fitted design columns.
  lf <- legacy_fit_on("log(age)", d, keep_frame = FALSE)
  fitted <- predict(lf$current, type = "cumulative_hazard")
  expect_equal(predict(lf$fit, type = "cumulative_hazard"), fitted)
  cols <- data.frame(time = d$int_dead[rows], log(d$age[rows]))
  names(cols)[2] <- "log(age)"
  got <- predict(lf$fit, newdata = cols, type = "cumulative_hazard")
  expect_length(got, length(rows))
  expect_equal(got / fitted[rows], rep(1, length(rows)), tolerance = 1e-8,
               ignore_attr = TRUE)
})
