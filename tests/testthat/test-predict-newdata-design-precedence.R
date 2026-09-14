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
    # The fit recorded that it ignored the formula (not merely no record).
    expect_identical(unname(attr(f$fit$x_list, "from_formula")),
                     c(FALSE, FALSE))
    want <- predict(f, type = "cumulative_hazard")[rows]
    got <- predict(f, newdata = nd, type = "cumulative_hazard")
    label <- if (is.null(tw)) "no windows" else "time_windows = 12"
    expect_length(got, length(rows))
    expect_equal(unname(got), unname(want), tolerance = 1e-10, label = label)
  }
})

test_that("a fit without the from_formula record routes by formula and columns", {
  skip_on_cran()  # multiphase fits
  # Every fit saved before the record existed (it is in no tag through
  # v1.2.9) is routed by the fallback hzr_gof() shares: a phase uses its
  # formula iff it has one AND its stored columns are not the inherited
  # ones.  The vector fit still takes the global route, a formula fit
  # still takes its phase formula.
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

test_that("a 1.0.3-era fit (no record, no frame, no design) keeps its phase formula", {
  skip_on_cran()  # a multiphase fit
  # Fitted on its phase formula ~ log(age), but saved before x_design,
  # data$frame (1.1.0) and the from_formula record existed.  A frame-based
  # fallback sent it down the global route (age times a log(age)
  # coefficient); the column test keeps it on its phase formula, as main
  # 7e50e2b did.  Values pinned from main on the same stripped object.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  f <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ log(age)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  f$data$x_design <- NULL
  f$fit$x_design <- NULL
  f$data$frame <- NULL
  attr(f$fit$x_list, "from_formula") <- NULL
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  got <- predict(f, type = "cumulative_hazard",
                 newdata = data.frame(time = d$int_dead[rows], age = d$age[rows]))
  expect_equal(unname(got), unname(predict(f, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
  expect_equal(unname(got),
               c(0.01628361325, 0.08219660741, 0.1292144135, 0.1674520008),
               tolerance = 1e-4)
})

test_that("time_windows: an unused phase formula without the record goes global", {
  skip_on_cran()  # a multiphase fit
  # The phase inherited the window-expanded global x (age_w1, age_w2).
  # The fallback compares with those window names, not with data$x's age,
  # so it keeps the phase on the global route: four values, matching.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  f <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    dist = "multiphase", time_windows = 12,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ log(age)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  attr(f$fit$x_list, "from_formula") <- NULL
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  got <- predict(f, type = "cumulative_hazard",
                 newdata = data.frame(time = d$int_dead[rows], age = d$age[rows]))
  expect_length(got, 4L)
  expect_equal(unname(got), unname(predict(f, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
})

test_that("coinciding global and phase formulas give the same answer either way", {
  skip_on_cran()  # a multiphase fit
  # Global ~ age with phase ~ age: without the record the phase's columns
  # equal the inherited ones, so the fallback takes the global route; with
  # the record it takes the phase route.  Both must be the fitted answer.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  f <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant")),
    fit = TRUE))
  rows <- c(1, 50, 150, 250)
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  want <- unname(predict(f, type = "cumulative_hazard")[rows])
  via_phase <- unname(predict(f, newdata = nd, type = "cumulative_hazard"))
  g <- f
  attr(g$fit$x_list, "from_formula") <- NULL
  via_global <- unname(predict(g, newdata = nd, type = "cumulative_hazard"))
  expect_equal(via_phase, want, tolerance = 1e-10)
  expect_equal(via_global, want, tolerance = 1e-10)
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
  # grp arrives with its levels reversed, so ignoring the fit's levels
  # would code it wrong; both levels are present among the rows.
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows],
                   grp = factor(as.character(d$grp[rows]),
                                levels = c("young", "old")))
  expect_setequal(as.character(nd$grp), c("old", "young"))
  got <- predict(m, newdata = nd, type = "cumulative_hazard")
  expect_length(got, length(rows))
  expect_equal(unname(got),
               unname(predict(m, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
})

test_that("an extra newdata column cannot mask a formula constant", {
  # cutoff is a constant of the formula's environment, not a fitting-data
  # variable, so a same-named newdata column is an unused extra.  Passing
  # all of newdata to model.frame() let it replace the constant: an extra
  # cutoff = 0 turned I(30 > 50) into I(30 > 0), silently (0.425 for 0.191).
  cutoff <- 50
  w <- hazard(survival::Surv(int_dead, dead) ~ I(age > cutoff) + mal,
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, b_gt = 0.8, b_mal = 0.3))
  want_eta <- 0.8 * (30 > 50) + 0.3
  nd <- data.frame(time = 2, age = 30, mal = 1, cutoff = 0)
  expect_equal(unname(predict(w, newdata = nd, type = "cumulative_hazard")),
               .dp_weibull(want_eta, 2), tolerance = 1e-12)
  expect_equal(predict(w, newdata = nd[, -1], type = "linear_predictor"),
               want_eta, tolerance = 1e-12)
  # The extra column is ignored, as documented: same as without it.
  expect_equal(predict(w, newdata = nd, type = "cumulative_hazard"),
               predict(w, newdata = nd[, c("time", "age", "mal")],
                       type = "cumulative_hazard"), tolerance = 1e-12)
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
