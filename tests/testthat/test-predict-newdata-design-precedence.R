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
  # Pinned from c4678f1, before #272 (a fitted model, so 1e-6).
  expect_length(gof$par_cumhaz, 270L)
  expect_equal(unname(gof$par_cumhaz[c(1, 135, 270)]),
               c(0.02303918854, 0.2220404487, 0.3080652665),
               tolerance = 1e-6)
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
  # Pinned from c4678f1, before #272 (a fitted model, so 1e-6).
  expect_equal(sum(dec$expected), 67.99859561, tolerance = 1e-6)
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
