# Which route builds the covariate design at newdata (#272).
#
# newdata may carry the formula's variables (grp = "old") or the fit's
# design columns (grpyoung = 1).  When it carried both and they disagreed,
# the design columns silently won: grp = "old" plus grpyoung = 1 gave the
# "young" value.  User newdata is now rebuilt from the formula's variables
# whenever they are all present; a design-named column is then an unused
# extra.  Design columns are used only when the variables are absent (a fit
# saved before the design was stored, or a caller passing design columns),
# and hzr_deciles() / hzr_gof(), which evaluate at the fitted design, say so.

.dp_avc <- local({
  data(avc, package = "TemporalHazard")
  d <- na.omit(avc)
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  d
})

.dp_th <- c(mu = 0.01, nu = 0.5, b_age = 0.004, b_young = 0.7)

.dp_weibull <- function(eta, time) (0.01 * time)^0.5 * exp(eta)

test_that("a contradicting design column does not override the variable", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
              dist = "weibull", theta = .dp_th)
  old <- .dp_weibull(0.004 * 60, 2)
  young <- .dp_weibull(0.004 * 60 + 0.7, 2)
  # The two answers differ, so the check below can fail.
  expect_gt(young / old, 1.5)

  for (type in c("cumulative_hazard", "survival")) {
    want <- if (type == "survival") exp(-old) else old
    got <- predict(w, type = type,
                   newdata = data.frame(time = 2, age = 60, grp = "old",
                                        grpyoung = 1))
    expect_equal(unname(got), want, tolerance = 1e-12, label = type)
  }
  expect_equal(
    predict(w, type = "linear_predictor",
            newdata = data.frame(age = 60, grp = "old", grpyoung = 1)),
    0.004 * 60, tolerance = 1e-12)
})

test_that("some variables beside all design columns is refused, not guessed", {
  # With `sex` missing the variables cannot be rebuilt, and taking the
  # design columns would silently ignore grp = "old" (it gave 0.362, the
  # "young" value).  Neither answer is safe, so it stops.
  d <- .dp_avc
  d$sex <- factor(ifelse(d$mal == 1, "M", "F"))
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp + sex, data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.7, 0.2))
  expect_error(
    predict(w, type = "cumulative_hazard",
            newdata = data.frame(time = 2, age = 60, grp = "old",
                                 grpyoung = 1, sexM = 0)),
    "gives the formula variable\\(s\\) 'grp'.*lacks 'sex'"
  )
  # A transform: log(age) is a design column, age the missing variable.
  w2 <- hazard(survival::Surv(int_dead, dead) ~ log(age) + grp, data = d,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.2, 0.7))
  nd2 <- data.frame(time = 2, check.names = FALSE, `log(age)` = log(60),
                    grp = "old", grpyoung = 1)
  expect_error(predict(w2, type = "cumulative_hazard", newdata = nd2),
               "gives the formula variable\\(s\\) 'grp'.*lacks 'age'")
})

test_that("a variable that other columns are built from is not left stale", {
  # Copy the design, change `age`: I(age^2) and age:grpyoung still hold the
  # old values.  Taking them as given gave 0.362 where age * grp with
  # grp "young" gives 0.660.  `age` is a formula variable feeding other
  # columns, so with `grp` missing this is a mix, and it stops.
  d <- .dp_avc
  w <- hazard(survival::Surv(int_dead, dead) ~ age * grp, data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.7, 0.01))
  nd <- data.frame(time = 2, age = 60, grpyoung = 1, `age:grpyoung` = 0,
                   check.names = FALSE)
  expect_error(predict(w, type = "cumulative_hazard", newdata = nd),
               "gives the formula variable\\(s\\) 'age'.*lacks 'grp'")
  w2 <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2) + grp, data = d,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0, 0.7))
  nd2 <- data.frame(time = 2, age = 60, `I(age^2)` = 0, grpyoung = 1,
                    check.names = FALSE)
  expect_error(predict(w2, type = "cumulative_hazard", newdata = nd2),
               "gives the formula variable\\(s\\) 'age'.*lacks 'grp'")
  # Given grp as well, the design is rebuilt and nothing is stale.
  nd3 <- data.frame(time = 2, age = 60, grp = "young")
  expect_equal(unname(predict(w, type = "cumulative_hazard", newdata = nd3)),
               .dp_weibull(0.004 * 60 + 0.7 + 0.01 * 60, 2), tolerance = 1e-12)
})

test_that("a fit saved before the design was stored refuses extra columns", {
  # A formula fit from 1.2.10 or earlier has no x_design, so nothing can
  # tell a formula variable from an unused column.  Taking the design
  # columns beside a contradicting grp swapped the answer; main (4b68020)
  # stopped instead ("Number of parameters insufficient ..."), and so does
  # this now.  Values pinned from main at 4b68020 on the same object.
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
              dist = "weibull", theta = .dp_th)
  w$data$x_design <- NULL
  expect_error(
    predict(w, type = "cumulative_hazard",
            newdata = data.frame(time = 2, age = 60, grp = c("old", "young"),
                                 grpyoung = c(1, 0))),
    "earlier version.*also has 'grp'.*pass only the design columns"
  )
  # Design-only newdata stays correct: the truth, which main also gave.
  got <- predict(w, type = "cumulative_hazard",
                 newdata = data.frame(grpyoung = c(0, 1), time = 2, age = 60))
  expect_equal(unname(got), .dp_weibull(0.004 * 60 + 0.7 * c(0, 1), 2),
               tolerance = 1e-12)
})

test_that("hzr_gof() and hzr_deciles() still run on a fit saved before x_design", {
  # Their newdata is design-level and marked, so the legacy refusal must
  # not touch it.  Values pinned from b0aa955 (before the refusal) on the
  # same legacy object; a fitted model, so 1e-4 (relative: all >= 0.02).
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2) + grp,
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0, 0, 0), fit = TRUE)
  w$data$x_design <- NULL
  g <- hzr_gof(w)
  expect_length(g$par_cumhaz, 270L)
  expect_equal(unname(g$par_cumhaz[c(1, 135, 270)]),
               c(0.02289393213, 0.2218081341, 0.3079777839), tolerance = 1e-4)
  expect_equal(sum(hzr_deciles(w, time = 60)$expected), 67.99999988,
               tolerance = 1e-4)
})

test_that("a vector-interface fit still ignores a column it does not use", {
  # No formula, so no variable that could contradict: an extra column is
  # just unused, as for a fit with a stored design.
  d <- .dp_avc
  v <- hazard(time = d$int_dead, status = d$dead,
              x = cbind(age = d$age, mal = d$mal),
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  got <- predict(v, type = "cumulative_hazard",
                 newdata = data.frame(time = 2, age = 60, mal = 1, junk = 9))
  expect_equal(unname(got), .dp_weibull(0.004 * 60 + 0.3, 2), tolerance = 1e-12)
})

test_that("design columns alone are still taken by name", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
              dist = "weibull", theta = .dp_th)
  got <- predict(w, type = "cumulative_hazard",
                 newdata = data.frame(grpyoung = c(0, 1), time = 2, age = 60))
  expect_equal(unname(got),
               .dp_weibull(0.004 * 60 + 0.7 * c(0, 1), 2), tolerance = 1e-12)
  # Pinned from c4678f1, before #272: this route is unchanged.
  expect_equal(unname(got), c(0.179781779, 0.3620360441), tolerance = 1e-9)
})

test_that("hzr_gof() still evaluates at the design-column means", {
  # For I(age^2) the mean of the design column is mean(age^2), not
  # mean(age)^2: the rebuild route would give the latter.
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2), data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  th <- w$fit$theta
  x_bar <- colMeans(w$data$x)
  expect_gt(abs(x_bar[[2]] / x_bar[[1]]^2 - 1), 0.1)   # the routes differ
  gof <- hzr_gof(w)
  want <- (th[[1]] * gof$time)^th[[2]] * exp(sum(x_bar * th[3:4]))
  expect_equal(unname(gof$par_cumhaz), want, tolerance = 1e-10)
  # Pinned from c4678f1, before #272. A fitted, poorly scaled model, so
  # 1e-4 to absorb cross-platform optimizer noise; the rebuild route would
  # move it by about 31%.
  expect_length(gof$par_cumhaz, 270L)
  expect_equal(unname(gof$par_cumhaz[c(1, 135, 270)]),
               c(0.02303918854, 0.2220404487, 0.3080652665),
               tolerance = 1e-4)
})

test_that("hzr_deciles() still evaluates each subject at its own design row", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2), data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  th <- w$fit$theta
  dec <- hzr_deciles(w, time = 60)
  # Expected events are the sum of each subject's H at their own follow-up.
  h <- (th[[1]] * w$data$time)^th[[2]] *
    exp(as.numeric(w$data$x %*% th[3:4]))
  expect_equal(sum(dec$expected), sum(h), tolerance = 1e-8)
  # Pinned from c4678f1, before #272. 1e-4 (relative: the value is ~68)
  # absorbs cross-platform optimizer noise on a poorly scaled fit. This
  # pins the value; the check against the design rows above pins the route,
  # since here the two routes happen to agree.
  expect_equal(sum(dec$expected), 67.99859561, tolerance = 1e-4)
})

test_that("hzr_deciles() uses the fitted rows even if a formula constant changed", {
  # ~ age + I(age > k): its design frame carries `age`, a formula variable,
  # so without the design-level marker it would be rebuilt, and I(age > k)
  # re-evaluated with the current k.  The fitted rows must win.
  d <- .dp_avc
  k <- 50
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age > k), data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  th <- w$fit$theta
  want <- sum((th[[1]] * w$data$time)^th[[2]] *
                exp(as.numeric(w$data$x %*% th[3:4])))
  k <- 500   # now every I(age > k) would be FALSE on a rebuild
  expect_gt(sum(w$data$x[, 2]), 0)            # the fitted column is not all 0
  expect_equal(sum(hzr_deciles(w, time = 60)$expected), want,
               tolerance = 1e-8)
})

test_that("a multiphase fit saved before the design was stored refuses too", {
  skip_on_cran()  # a multiphase fit
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  m$data$x_design <- NULL
  expect_error(
    predict(m, type = "cumulative_hazard",
            newdata = data.frame(time = 2, age = 60, grp = c("old", "young"),
                                 grpyoung = c(1, 0))),
    "earlier version.*also has 'grp'.*pass only the design columns"
  )
  # Pinned from main at 4b68020 on the same object (a fitted model: 1e-4).
  got <- predict(m, type = "cumulative_hazard",
                 newdata = data.frame(time = 2, age = 60, grpyoung = c(0, 1)))
  expect_equal(unname(got), c(0.0140056, 0.320817), tolerance = 1e-4)
})

test_that("multiphase time_windows: the rebuilt global design is expanded", {
  skip_on_cran()  # a multiphase fit
  # The global design a formula-less phase inherits is window-expanded at
  # fit time (age_w1, age_w2).  Rebuilt unexpanded at newdata, two rows met
  # two per-window coefficients and %*% returned four values, silently.
  # Reference: predict() at the fitted data, which uses the stored design.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    time_windows = 12,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  # One subject on each side of the cut, so both windows are exercised.
  rows <- c(which(d$int_dead <= 12)[1], which(d$int_dead > 12)[1])
  want <- predict(m, type = "cumulative_hazard")[rows]
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  got <- predict(m, newdata = nd, type = "cumulative_hazard")
  expect_length(got, 2L)
  expect_equal(unname(got), unname(want), tolerance = 1e-10)
})

test_that("a phase-formula factor with reordered levels codes as the fit did", {
  skip_on_cran()  # a multiphase fit
  # Global ~ age plus a phase formula ~ grp.  #292 makes this phase path
  # reachable at newdata; before #290 kept the fit's levels, reversed
  # levels came back swapped (0.213738/0.21485) with no error.  The answer
  # must be the structural one, whatever order or type grp arrives in.
  set.seed(1)
  pf <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = .dp_avc, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ grp),
      constant = hzr_phase("constant")),
    fit = TRUE))
  th <- pf$fit$theta
  tt <- c(2, 20)
  base <- predict(pf, newdata = data.frame(time = tt), type = "cumulative_hazard",
                  decompose = TRUE)
  # Row 1 is "old" (early: no grp effect), row 2 "young"; the constant
  # phase inherits the global age effect at age 60.
  ref <- base$early * exp(th[["early.grpyoung"]] * c(0, 1)) +
    base$constant * exp(sum(th[grep("^constant\\.age", names(th))]) * 60)
  expect_gt(abs(th[["early.grpyoung"]]), 0.5)          # levels matter here
  for (g in list(factor(c("old", "young"), levels = c("young", "old")),
                 c("old", "young"),
                 factor(c("old", "young"), levels = c("old", "young")))) {
    got <- predict(pf, newdata = data.frame(time = tt, age = 60, grp = g),
                   type = "cumulative_hazard")
    expect_equal(unname(got), ref, tolerance = 1e-10)
  }
  # The structural values, pinned (a fitted model: 1e-4, relative).
  expect_equal(ref, c(0.05608905, 0.3599925), tolerance = 1e-4)
})

test_that("a phase formula the fit did not use is not rebuilt at newdata", {
  skip_on_cran()  # multiphase fits
  # The fit uses a phase's own formula only on the formula interface; a
  # vector-interface fit ignores hzr_phase(formula = ), the phase inherits
  # the global x, and attr(fit$x_list, "from_formula") records FALSE.
  # predict() routed on !is.null(ph$formula) alone and rebuilt the unused
  # formula: log(age) times a coefficient fitted on age, and under
  # time_windows eight values for four rows -- on main 7e50e2b as well.
  # Reference: predict() at the fitted rows (the stored design).
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  for (tw in list(NULL, 12)) {
    set.seed(1)
    f <- suppressWarnings(hazard(
      time = d$int_dead, status = d$dead, x = cbind(age = d$age),
      dist = "multiphase", time_windows = tw,
      phases = list(
        early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes", formula = ~ log(age)),
        constant = hzr_phase("constant")),
      fit = TRUE))
    expect_false(any(attr(f$fit$x_list, "from_formula")))
    want <- predict(f, type = "cumulative_hazard")[rows]
    got <- predict(f, newdata = nd, type = "cumulative_hazard")
    label <- if (is.null(tw)) "no windows" else "time_windows = 12"
    expect_length(got, length(rows))
    expect_equal(unname(got), unname(want), tolerance = 1e-10, label = label)
  }
})

test_that("a fit without the from_formula record routes by formula and data", {
  skip_on_cran()  # multiphase fits
  # An older fit has no from_formula record.  Then a phase counts as
  # using its formula only if the fit had data (the formula interface):
  # the vector fit still takes the global route, a formula fit still takes
  # its phase formula.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  rows <- c(1, 50, 150)
  set.seed(1)
  vf <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ log(age)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  attr(vf$fit$x_list, "from_formula") <- NULL
  expect_equal(
    unname(predict(vf, type = "cumulative_hazard",
                   newdata = data.frame(time = d$int_dead[rows],
                                        age = d$age[rows]))),
    unname(predict(vf, type = "cumulative_hazard")[rows]), tolerance = 1e-10)
  set.seed(1)
  ff <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ grp),
      constant = hzr_phase("constant")),
    fit = TRUE))
  attr(ff$fit$x_list, "from_formula") <- NULL
  expect_equal(
    unname(predict(ff, type = "cumulative_hazard",
                   newdata = data.frame(time = d$int_dead[rows],
                                        age = d$age[rows],
                                        grp = d$grp[rows]))),
    unname(predict(ff, type = "cumulative_hazard")[rows]), tolerance = 1e-10)
})

test_that("time_windows with a factor and several covariates matches the fit", {
  skip_on_cran()  # a multiphase fit
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age + grp, data = d,
    dist = "multiphase", time_windows = 12,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows],
                   grp = as.character(d$grp[rows]))
  got <- predict(m, newdata = nd, type = "cumulative_hazard")
  expect_length(got, length(rows))
  expect_equal(unname(got),
               unname(predict(m, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
  expect_true(length(unique(d$grp[rows])) == 2L ||
                any(grepl("grpyoung", colnames(m$fit$x_list$constant))))
})

test_that("a multiphase global design takes the variable over its column", {
  skip_on_cran()  # a multiphase fit
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ grp, data = .dp_avc,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  as_old <- predict(fit, type = "cumulative_hazard",
                    newdata = data.frame(time = c(1, 4), grp = "old"))
  both <- predict(fit, type = "cumulative_hazard",
                  newdata = data.frame(time = c(1, 4), grp = "old",
                                       grpyoung = 1))
  as_young <- predict(fit, type = "cumulative_hazard",
                      newdata = data.frame(time = c(1, 4), grp = "young"))
  expect_gt(max(abs(as_young / as_old - 1)), 0.05)
  expect_equal(both, as_old, tolerance = 1e-12)
})
