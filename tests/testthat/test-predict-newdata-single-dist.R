# predict(newdata = ) for the single-distribution models matched covariates
# to coefficients by column POSITION, so newdata with its columns in another
# order returned a wrong value with no error (#267): a Weibull cumulative
# hazard of 0.23 became 4.3e14, and its survival became 0.  Covariates are
# now matched by name.
#
# Every fixture is unfitted with a supplied theta, so each expected value is
# exact: the prediction must equal the family's closed form at that theta.

.sd_avc <- local({
  data(avc, package = "TemporalHazard")
  d <- na.omit(avc)
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  d
})

.sd_theta <- list(
  weibull     = c(mu = 0.01, nu = 0.5, b_age = 0.004, b_mal = 0.6),
  exponential = c(log_lambda = -4, b_age = 0.004, b_mal = 0.6),
  loglogistic = c(a = -4, b = -0.5, b_age = 0.004, b_mal = 0.6),
  lognormal   = c(mu = 4, log_sigma = 0.8, b_age = -0.004, b_mal = -0.6)
)

# The families' cumulative hazards, written out from their definitions.
.sd_cumhaz <- function(dist, th, time, eta) {
  switch(dist,
    weibull     = (th[[1]] * time)^th[[2]] * exp(eta),
    exponential = exp(th[[1]]) * time * exp(eta),
    loglogistic = log1p(exp(th[[1]]) * time^exp(th[[2]]) * exp(eta)),
    lognormal   = -stats::pnorm(-(log(time) - th[[1]] - eta) / exp(th[[2]]),
                                log.p = TRUE)
  )
}

.sd_obj <- function(dist, formula = survival::Surv(int_dead, dead) ~ age + mal,
                    theta = .sd_theta[[dist]]) {
  hazard(formula, data = .sd_avc, dist = dist, theta = theta)
}

.sd_time <- c(1, 5)
.sd_in_order <- data.frame(time = .sd_time, age = 60, mal = 1)
.sd_reordered <- data.frame(mal = 1, time = .sd_time, age = 60)

# In-order values at .sd_in_order, pinned from the closed form above.
.sd_pinned <- list(
  weibull     = c(0.231637, 0.517955),
  exponential = c(0.0424257, 0.212129),
  loglogistic = c(0.0415504, 0.106708),
  lognormal   = c(0.0810163, 0.278381)
)

test_that("reordered columns give the closed-form answer, all four families", {
  for (dist in names(.sd_theta)) {
    th <- .sd_theta[[dist]]
    obj <- .sd_obj(dist)
    want <- .sd_cumhaz(dist, th, .sd_time,
                       eta = th[["b_age"]] * 60 + th[["b_mal"]] * 1)
    expect_equal(want, .sd_pinned[[dist]], tolerance = 1e-5, label = dist)

    for (nd in list(.sd_in_order, .sd_reordered)) {
      expect_equal(predict(obj, newdata = nd, type = "cumulative_hazard"),
                   want, tolerance = 1e-12, label = dist)
      expect_equal(predict(obj, newdata = nd, type = "survival"),
                   exp(-want), tolerance = 1e-12, label = dist)
    }
    # The two orders are far apart when read positionally, so the check
    # above can fail.
    swapped <- .sd_cumhaz(dist, th, .sd_time,
                          eta = th[["b_age"]] * 1 + th[["b_mal"]] * 60)
    expect_gt(max(abs(swapped / want - 1)), 0.5, label = dist)
  }
})

test_that("reordered columns give the same linear predictor and hazard", {
  obj <- .sd_obj("weibull")
  th <- .sd_theta$weibull
  eta <- th[["b_age"]] * 60 + th[["b_mal"]] * 1
  nd <- .sd_reordered[, c("mal", "age")]
  expect_equal(predict(obj, newdata = nd, type = "linear_predictor"),
               rep(eta, 2), tolerance = 1e-12)
  expect_equal(predict(obj, newdata = nd, type = "hazard"),
               rep(exp(eta), 2), tolerance = 1e-12)
})

test_that("a missing covariate column is an error naming it", {
  obj <- .sd_obj("weibull")
  expect_error(
    predict(obj, newdata = data.frame(time = 1, age = 60),
            type = "cumulative_hazard"),
    "lacks the covariate column\\(s\\) 'mal'"
  )
  expect_error(
    predict(obj, newdata = data.frame(age = 60), type = "linear_predictor"),
    "lacks the covariate column\\(s\\) 'mal'"
  )
})

test_that("a missing column is an error even when a same-named object exists", {
  # The formula's environment holds a `mal`.  A newdata without the column
  # must not pick that value up: it is a covariate of the fit's data, not a
  # constant of the formula.
  mal <- 0
  obj <- hazard(survival::Surv(int_dead, dead) ~ age + mal, data = .sd_avc,
                dist = "weibull", theta = .sd_theta$weibull)
  expect_error(
    predict(obj, newdata = data.frame(time = 1, age = 60),
            type = "cumulative_hazard"),
    "lacks the covariate column\\(s\\) 'mal'"
  )
})

test_that("design-column newdata is taken by name, as the diagnostics pass it", {
  # hzr_deciles() and hzr_gof() build newdata from the fitted design matrix,
  # so a factor arrives as `grpyoung` and a transform as `log(age)`: design
  # columns, not the formula's variables.
  th <- c(mu = 0.01, nu = 0.5, b1 = 0.4, b2 = 0.7)
  for (f in list(survival::Surv(int_dead, dead) ~ log(age) + grp,
                 survival::Surv(int_dead, dead) ~ grp + log(age))) {
    obj <- hazard(f, data = .sd_avc, dist = "weibull", theta = th)
    x <- obj$data$x
    nd <- as.data.frame(x[c(3, 1, 2), , drop = FALSE])
    nd$time <- 2
    nd <- nd[, rev(names(nd))]
    want <- .sd_cumhaz("weibull", th, 2,
                       eta = as.numeric(x[c(3, 1, 2), ] %*% th[3:4]))
    expect_equal(predict(obj, newdata = nd, type = "cumulative_hazard"),
                 want, tolerance = 1e-12)
  }
})

test_that("hzr_deciles() and hzr_gof() run on factor and transformed fits", {
  for (f in list(survival::Surv(int_dead, dead) ~ age + grp,
                 survival::Surv(int_dead, dead) ~ log(age) + mal)) {
    obj <- hazard(f, data = .sd_avc, dist = "weibull",
                  theta = c(mu = 0.01, nu = 0.5, b1 = 0, b2 = 0), fit = TRUE)
    th <- obj$fit$theta
    expect_true(isTRUE(obj$fit$converged))
    expect_no_error(hzr_deciles(obj, time = 60))
    gof <- hzr_gof(obj)
    # hzr_gof() evaluates at the design-matrix means.
    eta <- sum(colMeans(obj$data$x) * th[3:4])
    want <- .sd_cumhaz("weibull", th, gof$time, eta = eta)
    expect_equal(unname(gof$par_cumhaz), want, tolerance = 1e-10)
  }
})

test_that("a fit saved before the design was stored matches design columns", {
  th <- c(mu = 0.01, nu = 0.5, b_age = 0.004, b_young = 0.7)
  obj <- .sd_obj("weibull", survival::Surv(int_dead, dead) ~ age + grp,
                 theta = th)
  obj$data$x_design <- NULL   # as saved by an earlier version
  want <- .sd_cumhaz("weibull", th, 2, eta = th[["b_age"]] * 60 + th[[4]])
  # unname(): at a single time the Weibull branch carries theta[1]'s name.
  expect_equal(unname(predict(obj, newdata = data.frame(grpyoung = 1,
                                                        time = 2, age = 60),
                              type = "cumulative_hazard")),
               want, tolerance = 1e-12)
  expect_error(
    predict(obj, newdata = data.frame(time = 2, age = 60, grp = "young"),
            type = "cumulative_hazard"),
    "lacks the covariate column\\(s\\) 'grpyoung'"
  )
})

test_that("a column the model does not use is ignored", {
  obj <- .sd_obj("exponential")
  base <- predict(obj, newdata = .sd_in_order, type = "cumulative_hazard")
  extra <- cbind(unused = 99, .sd_in_order, label = "q")
  expect_equal(predict(obj, newdata = extra, type = "cumulative_hazard"),
               base, tolerance = 1e-12)
})

test_that("a factor covariate can be given as a label", {
  th <- c(mu = 0.01, nu = 0.5, b_age = 0.004, b_young = 0.7)
  obj <- .sd_obj("weibull", survival::Surv(int_dead, dead) ~ age + grp,
                 theta = th)
  want <- .sd_cumhaz("weibull", th, .sd_time,
                     eta = th[["b_age"]] * 60 + th[["b_young"]])
  nd <- data.frame(grp = "young", time = .sd_time, age = 60)
  expect_equal(predict(obj, newdata = nd, type = "cumulative_hazard"),
               want, tolerance = 1e-12)
  expect_equal(predict(obj, newdata = nd, type = "linear_predictor"),
               rep(th[["b_age"]] * 60 + th[["b_young"]], 2), tolerance = 1e-12)
})

test_that("time-varying coefficients match by name before the expansion", {
  # One cut at 3: time 1 falls in window 1, time 5 in window 2, and the
  # expanded columns are age_w1, mal_w1, age_w2, mal_w2.
  th <- c(mu = 0.01, nu = 0.5, age_w1 = 0.004, mal_w1 = 0.6,
          age_w2 = -0.002, mal_w2 = 0.2)
  obj <- hazard(survival::Surv(int_dead, dead) ~ age + mal, data = .sd_avc,
                dist = "weibull", theta = th, time_windows = 3)
  eta <- c(th[["age_w1"]] * 60 + th[["mal_w1"]],
           th[["age_w2"]] * 60 + th[["mal_w2"]])
  want <- .sd_cumhaz("weibull", th, .sd_time, eta = eta)
  expect_equal(predict(obj, newdata = .sd_reordered,
                       type = "cumulative_hazard"),
               want, tolerance = 1e-12)
  expect_equal(predict(obj, newdata = .sd_reordered[, c("time", "mal", "age")],
                       type = "linear_predictor"),
               eta, tolerance = 1e-12)
})

test_that("the vector interface matches a named x by name", {
  d <- .sd_avc
  th <- .sd_theta$lognormal
  obj <- hazard(time = d$int_dead, status = d$dead,
                x = cbind(age = d$age, mal = d$mal),
                dist = "lognormal", theta = th)
  want <- .sd_cumhaz("lognormal", th, .sd_time,
                     eta = th[["b_age"]] * 60 + th[["b_mal"]] * 1)
  expect_equal(predict(obj, newdata = .sd_reordered,
                       type = "cumulative_hazard"),
               want, tolerance = 1e-12)
})
