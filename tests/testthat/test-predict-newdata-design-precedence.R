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

test_that("design columns alone are still taken by name", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
              dist = "weibull", theta = .dp_th)
  got <- predict(w, type = "cumulative_hazard",
                 newdata = data.frame(grpyoung = c(0, 1), time = 2, age = 60))
  expect_equal(unname(got),
               .dp_weibull(0.004 * 60 + 0.7 * c(0, 1), 2), tolerance = 1e-12)
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
