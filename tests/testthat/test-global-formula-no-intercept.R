# A global formula without an intercept builds the design of the same
# formula with one (#337). Every distribution carries its own intercept (its
# scale), so a design without one gave a factor a column for every level,
# collinear with the scale: the Weibull lost its standard errors, and a
# multiphase fit inheriting the design stopped 9.5 log-likelihood units short.
# The phase-formula counterpart is #303.

global_ni_data <- function() {
  d <- na.omit(avc[, c("int_dead", "dead", "age", "status")])
  d$grp <- factor(d$status)
  d
}

global_ni_parse <- function(rhs, d) {
  f <- stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs))
  .hzr_parse_formula(f, d)$x
}

test_that("a global formula without an intercept builds the intercept design (#337)", {
  d <- global_ni_data()
  pairs <- list(c("0 + grp", "grp"), c("grp - 1", "grp"),
                c("0 + age + grp", "age + grp"), c("0 + age", "age"))
  for (p in pairs) {
    x <- suppressWarnings(global_ni_parse(p[[1]], d))
    expect_identical(x, global_ni_parse(p[[2]], d), info = p[[1]])
  }
  expect_identical(colnames(suppressWarnings(global_ni_parse("0 + grp", d))),
                   c("grp2", "grp3", "grp4"))
})

test_that("the removal warns once, and only where it is ignored (#337)", {
  d <- global_ni_data()
  for (rhs in c("0 + grp", "grp - 1", "0 + age")) {
    # By CLASS as well as message: stepwise's muffler keys on the class
    # (stepwise-refit.R), so a warning that kept its wording and lost its
    # class would break the muffling with every message assertion passing.
    expect_warning(global_ni_parse(rhs, d), "removes the intercept",
                   class = "hzr_intercept_removed", info = rhs)
  }
  for (rhs in c("grp", "age + grp", "1", "0")) {
    expect_no_warning(global_ni_parse(rhs, d))
  }
})

test_that("a fit without a global intercept is the fit with one (#337)", {
  skip_on_cran()
  d <- global_ni_data()
  control <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = d,
                    dist = "weibull", theta = c(1, 1, 0, 0, 0, 0), fit = TRUE)
  expect_warning(
    fit <- hazard(survival::Surv(int_dead, dead) ~ 0 + age + grp, data = d,
                  dist = "weibull", theta = c(1, 1, 0, 0, 0, 0), fit = TRUE),
    "removes the intercept"
  )
  expect_identical(fit$fit$theta, control$fit$theta)
  expect_identical(fit$fit$vcov, control$fit$vcov)
  expect_true(all(is.finite(diag(fit$fit$vcov))))
  nd_w <- data.frame(time = c(0.5, 2), age = c(40, 70),
                     grp = factor(c("1", "3"), levels = levels(d$grp)))
  expect_identical(predict(fit, newdata = nd_w, type = "cumulative_hazard"),
                   predict(control, newdata = nd_w,
                           type = "cumulative_hazard"))

  phases <- list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                   fixed = "m"),
                 constant = hzr_phase("constant"))
  mp <- function(f) {
    hazard(f, data = d, dist = "multiphase", phases = phases, fit = TRUE,
           control = list(n_starts = 1L, conserve = FALSE))
  }
  mp_control <- suppressWarnings(mp(survival::Surv(int_dead, dead) ~ grp))
  mp_fit <- suppressWarnings(mp(survival::Surv(int_dead, dead) ~ 0 + grp))
  expect_identical(mp_fit$fit$objective, mp_control$fit$objective)
  expect_identical(mp_fit$fit$theta, mp_control$fit$theta)
  nd <- data.frame(time = c(0.5, 2), grp = factor(c("2", "4"),
                                                 levels = levels(d$grp)))
  expect_identical(predict(mp_fit, newdata = nd, type = "cumulative_hazard"),
                   predict(mp_control, newdata = nd,
                           type = "cumulative_hazard"))
})

test_that("a stepwise refit does not repeat the warning (#337)", {
  skip_on_cran()
  d <- na.omit(avc[, c("int_dead", "dead", "age", "mal", "orifice")])
  expect_warning(
    base <- hazard(survival::Surv(int_dead, dead) ~ 0 + age, data = d,
                   dist = "multiphase",
                   phases = list(early = hzr_phase("cdf", t_half = 0.15,
                                                   nu = 1.4, m = 1,
                                                   fixed = "m"),
                                 constant = hzr_phase("constant")),
                   fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)),
    "removes the intercept"
  )
  # A multiphase refit passes the global formula through unchanged, so
  # every step re-parses `~ 0 + age`: the steps taken are the refits.
  # Count by MESSAGE, not by class. Counting the class alone cannot fail
  # the way this test needs to: if the class were renamed, the muffler in
  # stepwise-refit.R would stop muffling and the warning WOULD repeat --
  # but a class counter would never see it, the count would stay 0, and
  # the catch-all handler below would swallow the repeats. The test would
  # pass while asserting the opposite of the truth. Matching the message
  # catches a repeat however it is classed.
  repeated <- 0L
  classed <- 0L
  count_warnings <- function(w) {
    if (grepl("removes the intercept", conditionMessage(w))) {
      repeated <<- repeated + 1L
      if (inherits(w, "hzr_intercept_removed")) classed <<- classed + 1L
    }
    invokeRestart("muffleWarning")
  }
  # Known positive: the same handler, on a call that DOES warn, must count
  # one of each. Without this, `repeated == 0` below is equally consistent
  # with "nothing repeated" and "this handler never fires at all".
  invisible(withCallingHandlers(global_ni_parse("0 + age", d),
                                warning = count_warnings))
  expect_identical(repeated, 1L)
  expect_identical(classed, 1L)

  repeated <- 0L
  classed <- 0L
  sw <- withCallingHandlers(
    suppressMessages(hzr_stepwise(base, scope = list(early = ~ mal,
                                                     constant = ~ orifice),
                                  data = d, trace = FALSE, max_steps = 2L)),
    warning = count_warnings
  )
  expect_gte(nrow(sw$steps), 1L)
  expect_identical(repeated, 0L)
})
