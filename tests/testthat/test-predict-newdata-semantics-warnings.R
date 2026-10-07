# predict(newdata = ) evaluates the formula on newdata as given, as
# predict.lm() does. Three cases where that silently differs from the fit
# now warn, naming the cause, and predict the same values as before
# (#331, #334, #335).

# This file predicts from models built with fit = FALSE on purpose, so the
# warning that those numbers come from starting values is switched off for
# this file only (#398). A file that does not expect the warning sees it as
# an ordinary leaked warning.
withr::local_options(TemporalHazard.warn_unfitted_prediction = FALSE)

.sw_data <- local({
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  d$dt <- as.difftime(d$age, units = "days")
  d
})

.sw_weibull <- function(rhs, beta) {
  f <- stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs),
                         env = globalenv())
  eval(bquote(hazard(.(f), data = .sw_data, dist = "weibull",
                     theta = .(c(mu = 0.01, nu = 0.5, beta)))))
}

.sw_cumhaz <- function(time, eta) (0.01 * time)^0.5 * exp(eta)

.sw_multiphase <- function(rhs) {
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = .sw_data, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = stats::as.formula(paste("~", rhs),
                                                    env = globalenv())),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
}

# ---- #331: a data-dependent function recomputed over newdata's rows -------

test_that("a term computing min() over newdata warns and predicts as before", {
  w <- .sw_weibull("I(age - min(age))", 0.004)
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  expect_warning(
    got <- predict(w, newdata = nd, type = "cumulative_hazard"),
    "I\\(age - min\\(age\\)\\).*min\\(\\)"
  )
  # min() is taken over newdata's two rows: 60, not the fitted minimum.
  expect_equal(got, .sw_cumhaz(c(1, 2), 0.004 * (c(60, 90) - 60)),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(min(.sw_data$age), 60)))
})

test_that("mean(), median() and max() inside I() warn too", {
  nd <- data.frame(time = 2, age = c(60, 90))
  for (fn in c("mean", "median", "max")) {
    w <- .sw_weibull(paste0("I(age / ", fn, "(age))"), 0.1)
    expect_warning(predict(w, newdata = nd, type = "cumulative_hazard"),
                   paste0(fn, "\\(\\)"), label = fn)
  }
})

test_that("scale() nested inside I() warns; scale() as a term does not", {
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  nested <- .sw_weibull("I(scale(age)^2)", 0.1)
  expect_warning(predict(nested, newdata = nd, type = "cumulative_hazard"),
                 "scale\\(\\)")
  top <- .sw_weibull("scale(age) + poly(age, 2)", c(0.1, 0.1, 0.1))
  expect_no_warning(predict(top, newdata = nd, type = "cumulative_hazard"))
})

test_that("factor() nested inside a term warns; a factor term does not", {
  nd <- data.frame(time = 2, grp = c("young", "old"))
  nested <- .sw_weibull("I(as.integer(factor(grp)))", 0.3)
  expect_warning(
    got <- predict(nested, newdata = nd, type = "cumulative_hazard"),
    "factor\\(\\)"
  )
  # Levels are recomputed from newdata's rows: "old" codes 1, "young" 2.
  expect_equal(got, .sw_cumhaz(2, 0.3 * c(2, 1)), tolerance = 1e-12)
  top <- .sw_weibull("factor(grp)", 0.7)
  expect_no_warning(predict(top, newdata = nd, type = "cumulative_hazard"))
})

test_that("a phase formula computing min() over newdata warns", {
  m <- .sw_multiphase("I(age - min(age))")
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  expect_warning(predict(m, newdata = nd, type = "cumulative_hazard"),
                 "phase 'early'.*min\\(\\)")
  ctrl <- .sw_multiphase("scale(age)")
  expect_no_warning(predict(ctrl, newdata = nd, type = "cumulative_hazard"))
})

# ---- #334: a newdata column of another type ------------------------------

test_that("a numeric covariate given as character warns and predicts as before", {
  w <- .sw_weibull("I(age > 50)", 0.4)
  nd <- data.frame(time = 2, age = c("154.6", "60"))
  expect_warning(
    got <- predict(w, newdata = nd, type = "cumulative_hazard"),
    "'age'.*character.*numeric"
  )
  # The comparison runs on strings: "154.6" > 50 is FALSE.
  expect_equal(got, .sw_cumhaz(2, 0.4 * c(FALSE, TRUE)), tolerance = 1e-12)
})

test_that("a difftime in other units warns and is used in those units", {
  w <- .sw_weibull("I(as.numeric(dt) / 7)", 0.01)
  nd <- data.frame(time = c(2, 2))
  nd$dt <- as.difftime(c(2, 10), units = "weeks")
  expect_warning(
    got <- predict(w, newdata = nd, type = "cumulative_hazard"),
    "'dt'.*weeks.*days"
  )
  expect_equal(got, .sw_cumhaz(2, 0.01 * c(2, 10) / 7), tolerance = 1e-12)
})

test_that("ordinary newdata types do not warn", {
  nd <- data.frame(time = c(1, 2), age = c(60L, 90L), grp = c("old", "young"))
  nd$dt <- as.difftime(c(3, 4), units = "days")
  w <- .sw_weibull("age + grp + I(as.numeric(dt))", c(0.004, 0.7, 0.001))
  # An integer for a double, and a label for a factor.
  expect_no_warning(predict(w, newdata = nd, type = "cumulative_hazard"))
  m <- .sw_multiphase("age + grp")
  expect_no_warning(predict(m, newdata = nd, type = "cumulative_hazard"))
})

# ---- #335: a legacy design rebuilt under this session's contrasts ---------

test_that("a legacy global design rebuilt under a custom contrasts function warns", {
  d <- .sw_data
  set.seed(1)
  d$g <- factor(sample(c("A", "B", "C"), nrow(d), TRUE))
  d$x <- ifelse(d$g == "C", 0, d$age / 100)
  w <- eval(bquote(hazard(survival::Surv(int_dead, dead) ~ x + x:g,
                          data = .(d), dist = "weibull",
                          theta = c(mu = 0.01, nu = 0.5, 0.3, 0.4, 0.5))))
  leg <- w
  leg$data$x_design <- NULL
  nd <- data.frame(time = 2, x = 0.5, g = c("A", "B", "C"))

  # Control: the default contrasts rebuild silently.
  expect_no_warning(
    ref <- predict(leg, type = "cumulative_hazard", newdata = nd)
  )
  expect_identical(ref, predict(w, type = "cumulative_hazard", newdata = nd))

  contr.custom <- function(n, contrasts = TRUE, sparse = FALSE) {
    m <- stats::contr.treatment(n, contrasts = contrasts, sparse = sparse)
    m[nrow(m), ] <- 1
    m
  }
  assign("contr.custom", contr.custom, envir = globalenv())
  withr::defer(rm("contr.custom", envir = globalenv()))
  withr::local_options(contrasts = c(unordered = "contr.custom",
                                     ordered = "contr.poly"))
  expect_warning(
    got <- predict(leg, type = "cumulative_hazard", newdata = nd),
    "contr\\.custom"
  )
  # The value is still the rebuild's: C is coded (1, 1), about 22% high.
  expect_equal(unname(got[1:2]), unname(ref[1:2]), tolerance = 1e-12)
  expect_gt(got[3] / ref[3], 1.2)
})

test_that("the #335 reproducer warns and keeps its measured values", {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  set.seed(1)
  d$g <- factor(sample(c("A", "B", "C"), nrow(d), TRUE))
  d$x <- ifelse(d$g == "C", 0, d$age / 100)
  w <- eval(bquote(hazard(survival::Surv(int_dead, dead) ~ x + x:g,
                          data = .(d), dist = "weibull",
                          theta = c(mu = 0.01, nu = 0.5, 0.3, 0.4, 0.5))))
  leg <- w
  leg$data$x_design <- NULL
  nd <- data.frame(time = 2, x = 0.5, g = c("A", "B", "C"))
  truth <- c(0.16430817, 0.20068646, 0.21097587)
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               truth, tolerance = 1e-7)

  contr.custom <- function(n, contrasts = TRUE, sparse = FALSE) {
    m <- stats::contr.treatment(n, contrasts = contrasts, sparse = sparse)
    m[nrow(m), ] <- 1
    m
  }
  assign("contr.custom", contr.custom, envir = globalenv())
  withr::defer(rm("contr.custom", envir = globalenv()))
  withr::local_options(contrasts = c(unordered = "contr.custom",
                                     ordered = "contr.poly"))
  expect_warning(
    got <- predict(leg, type = "cumulative_hazard", newdata = nd),
    "contr\\.custom"
  )
  expect_equal(unname(got), c(truth[1:2], 0.25768651), tolerance = 1e-7)
  # A current fit records its contrasts and is right, with no warning.
  expect_no_warning(
    cur <- predict(w, type = "cumulative_hazard", newdata = nd)
  )
  expect_equal(unname(cur), truth, tolerance = 1e-7)
})

# ---- once per user call ----------------------------------------------------

newdata_warnings <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  msgs[grepl("over the rows of 'newdata'|in the data the model was fitted",
             msgs)]
}

test_that("a global design inherited by two phases warns once per call", {
  inherit_fit <- function(rhs) {
    suppressWarnings(hazard(
      stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs),
                        env = globalenv()),
      data = .sw_data, dist = "multiphase",
      phases = list(
        early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
        constant = hzr_phase("constant")
      ),
      fit = TRUE
    ))
  }
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  expect_length(newdata_warnings(
    predict(inherit_fit("I(age - mean(age))"), newdata = nd,
            type = "cumulative_hazard")
  ), 1L)
  nd$age <- as.character(nd$age)
  expect_length(newdata_warnings(
    predict(inherit_fit("I(age > 50)"), newdata = nd, type = "survival")
  ), 1L)
})

test_that("hzr_gof() and hzr_deciles() on such a fit do not warn", {
  w <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ I(age - mean(age)), data = .sw_data,
    dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004), fit = TRUE
  ))
  expect_length(newdata_warnings(hzr_gof(w)), 0L)
  expect_length(newdata_warnings(hzr_deciles(w, time = 5)), 0L)
  # Known positive: the same fit does warn at user newdata.
  expect_length(newdata_warnings(
    predict(w, newdata = data.frame(time = 2, age = c(60, 90)),
            type = "cumulative_hazard")
  ), 1L)
})

# ---- r-reviewer round 1 ------------------------------------------------------

test_that("two phase formulas computing over the rows warn once per call", {
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = .sw_data, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ I(age - mean(age))),
      constant = hzr_phase("constant", formula = ~ I(age - mean(age)))
    ),
    fit = TRUE
  ))
  nd <- data.frame(time = c(1, 20), age = c(60, 90))
  got <- newdata_warnings(predict(m, newdata = nd, type = "cumulative_hazard"))
  expect_length(got, 1L)
  expect_match(got, "phase 'early'")
  expect_match(got, "phase 'constant'")
})

test_that("raw poly() and an unscaled scale() term do not warn", {
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  raw <- .sw_weibull("poly(age, 2, raw = TRUE)", c(0.01, 0.01))
  expect_no_warning(predict(raw, newdata = nd, type = "cumulative_hazard"))
  unscaled <- .sw_weibull("scale(age, center = FALSE, scale = FALSE)", 0.01)
  expect_no_warning(
    predict(unscaled, newdata = nd, type = "cumulative_hazard")
  )
})

test_that("ave(), IQR(), mad() and cut() with a named break count warn", {
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  for (rhs in c("I(age - ave(age))", "I(age / IQR(age))",
                "I(age / mad(age))")) {
    w <- .sw_weibull(rhs, 0.001)
    expect_warning(predict(w, newdata = nd, type = "cumulative_hazard"),
                   "over the rows", label = rhs)
  }
  kk <- 3
  assign("kk", kk, envir = globalenv())
  withr::defer(rm("kk", envir = globalenv()))
  w <- .sw_weibull("I(as.integer(cut(age, kk)))", 0.01)
  expect_warning(predict(w, newdata = nd, type = "cumulative_hazard"),
                 "cut\\(\\)")
  fixed <- .sw_weibull("I(as.integer(cut(age, c(-Inf, 100, Inf))))", 0.01)
  expect_no_warning(predict(fixed, newdata = nd, type = "cumulative_hazard"))
})

test_that("legacy contrasts do not warn when the design columns are used", {
  d <- .sw_data
  w <- eval(bquote(hazard(survival::Surv(int_dead, dead) ~ age + grp,
                          data = .(d), dist = "weibull",
                          theta = c(mu = 0.01, nu = 0.5, 0.004, 0.7))))
  leg <- w
  leg$data$x_design <- NULL
  assign("contr.same", stats::contr.treatment, envir = globalenv())
  withr::defer(rm("contr.same", envir = globalenv()))
  withr::local_options(contrasts = c(unordered = "contr.same",
                                     ordered = "contr.poly"))
  expect_no_warning(predict(leg, type = "cumulative_hazard",
                            newdata = data.frame(time = 1, age = 60,
                                                 grpyoung = 0)))
  # Known positive: the same session warns when the design is rebuilt.
  expect_warning(predict(leg, type = "cumulative_hazard",
                         newdata = data.frame(time = 1, age = 60,
                                              grp = "old")),
                 "contr\\.same")
})

# ---- predict() must not assign into newdata ---------------------------------

# A guard class whose assignment methods stop. It catches ANY assignment into
# newdata on these routes, not only one that a reference-semantics object such
# as a data.table would carry back to the caller.
`[[<-.hzrtest_noassign` <- function(x, i, value) stop("assigned into newdata")
`$<-.hzrtest_noassign` <- function(x, i, value) stop("assigned into newdata")

sw_guard <- function(df) {
  registerS3method("[[<-", "hzrtest_noassign", `[[<-.hzrtest_noassign`)
  registerS3method("$<-", "hzrtest_noassign", `$<-.hzrtest_noassign`)
  structure(df, class = c("hzrtest_noassign", "data.frame"))
}

test_that("the newdata guard class stops an assignment", {
  # Known positive: without this, the tests below could pass over a guard
  # that never fires.
  nd <- sw_guard(data.frame(x = 1))
  expect_error(nd$x <- 2, "assigned into newdata")
  expect_error(nd[["x"]] <- 2, "assigned into newdata")
})

test_that("no warning route assigns into newdata", {
  plain_chr <- data.frame(time = c(1, 2), age = c("60", "90"),
                          grp = c("old", "young"))
  plain_num <- data.frame(time = c(1, 2), age = c(60, 90),
                          grp = c("old", "young"))

  types <- .sw_weibull("I(age > 50) + grp", c(0.4, 0.7))
  expect_warning(want <- predict(types, newdata = plain_chr,
                                 type = "cumulative_hazard"), "'age'")
  expect_warning(got <- predict(types, newdata = sw_guard(plain_chr),
                                type = "cumulative_hazard"), "'age'")
  expect_equal(got, want, tolerance = 1e-12)

  rows <- .sw_weibull("I(age - mean(age)) + grp", c(0.004, 0.7))
  expect_warning(want <- predict(rows, newdata = plain_num,
                                 type = "cumulative_hazard"), "mean\\(\\)")
  expect_warning(got <- predict(rows, newdata = sw_guard(plain_num),
                                type = "cumulative_hazard"), "mean\\(\\)")
  expect_equal(got, want, tolerance = 1e-12)

  phase <- .sw_multiphase("I(age - mean(age))")
  nd <- data.frame(time = c(1, 2), age = c(60, 90))
  expect_warning(want <- predict(phase, newdata = nd,
                                 type = "cumulative_hazard"), "mean\\(\\)")
  expect_warning(got <- predict(phase, newdata = sw_guard(nd),
                                type = "cumulative_hazard"), "mean\\(\\)")
  expect_equal(got, want, tolerance = 1e-12)

  legacy <- .sw_weibull("age + grp", c(0.004, 0.7))
  legacy$data$x_design <- NULL
  contr.custom <- function(n, contrasts = TRUE, sparse = FALSE) {
    m <- stats::contr.treatment(n, contrasts = contrasts, sparse = sparse)
    m[nrow(m), ] <- 1
    m
  }
  assign("contr.custom", contr.custom, envir = globalenv())
  withr::defer(rm("contr.custom", envir = globalenv()))
  withr::local_options(contrasts = c(unordered = "contr.custom",
                                     ordered = "contr.poly"))
  one <- data.frame(time = 2, age = 60, grp = "old")
  expect_warning(want <- predict(legacy, newdata = one,
                                 type = "cumulative_hazard"), "contr\\.custom")
  expect_warning(got <- predict(legacy, newdata = sw_guard(one),
                                type = "cumulative_hazard"), "contr\\.custom")
  expect_equal(got, want, tolerance = 1e-12)
})

# ---- #347: a classed numeric column is read as its values -------------------

test_that("a classed numeric newdata column predicts as its values", {
  # A classed numeric whose stored doubles are not its values, as with
  # bit64's integer64 (#347). The class is local, so no dependency: its
  # as.double() method returns the values the doubles do not hold.
  registerS3method("as.double", "hzr_test_wrapped",
                   function(x, ...) attr(x, "values"))
  wrap <- function(v) {
    structure(rep(9e-300, length(v)), values = v, class = "hzr_test_wrapped")
  }
  plain <- data.frame(time = c(1, 2, 5), age = c(60, 90, 30))
  wrapped <- plain
  wrapped$age <- wrap(plain$age)
  wrapped_time <- plain
  wrapped_time$time <- wrap(plain$time)

  w <- .sw_weibull("age", 0.004)
  m <- .sw_multiphase("age")
  for (fit in list(w, m)) {
    for (type in c("survival", "cumulative_hazard")) {
      want <- predict(fit, newdata = plain, type = type)
      # The rows differ from one another, so reading the stored doubles
      # instead of the values cannot give this answer by accident.
      expect_false(isTRUE(all.equal(want, rep(want[[1L]], 3L))))
      expect_identical(predict(fit, newdata = wrapped, type = type), want)
      expect_identical(predict(fit, newdata = wrapped_time, type = type),
                       want)
    }
  }
})
