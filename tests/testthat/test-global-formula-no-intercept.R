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
    # Exactly ONCE per call: expect_warning() alone passes on a duplicate.
    ws <- list()
    withCallingHandlers(global_ni_parse(rhs, d), warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    expect_length(ws, 1L)
    expect_s3_class(ws[[1L]], "hzr_intercept_removed")
    expect_match(conditionMessage(ws[[1L]]), "removes the intercept",
                 info = rhs)
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

test_that("a single-distribution no-op refit does not repeat the warning (#337)", {
  skip_on_cran()
  # The single-distribution refit (.hzr_refit_with_scope) muffles the intercept
  # warning too, and until now nothing tested it. An ordinary forward or
  # backward step cannot: .hzr_formula_update() rebuilds the right-hand side
  # from the term labels, which never carry `0` or `- 1`, so every real add or
  # drop restores the intercept and there is nothing to muffle -- a test built
  # on one would pass with the muffler removed. The muffler is live only where
  # the formula comes back UNCHANGED: adding a variable already present, or
  # dropping one that is absent. Those refits re-parse `~ 0 + x1`.
  set.seed(337)
  n <- 200
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  df$time <- stats::rexp(n) * exp(-0.8 * df$x1 - 0.9 * df$x2)
  df$status <- 1L
  seen <- 0L
  count_warnings <- function(w) {
    if (grepl("removes the intercept", conditionMessage(w))) seen <<- seen + 1L
    invokeRestart("muffleWarning")
  }
  # Known positive: the same handler over the fit that DOES warn counts one.
  base <- withCallingHandlers(
    hazard(survival::Surv(time, status) ~ 0 + x1, data = df,
           dist = "weibull", theta = c(0.5, 1, 0), fit = TRUE),
    warning = count_warnings
  )
  expect_identical(seen, 1L)

  for (case in list(c("add", "x1"), c("drop", "x2"))) {
    seen <- 0L
    refit <- withCallingHandlers(
      .hzr_refit_with_scope(base, case[[1]], case[[2]], data = df),
      warning = count_warnings
    )
    # The no-op path was taken: the formula is still intercept-free. Without
    # this, a zero below could come from a path that had restored it.
    expect_identical(attr(stats::terms(refit$call$formula), "intercept"), 0L,
                     info = paste(case, collapse = " "))
    expect_identical(seen, 0L, info = paste(case, collapse = " "))
  }
})
