# predict(se.fit = TRUE) over a covariance it cannot use (#586).
#
# Two silent defects on main e7b621c5, found by the r-reviewer pass over
# #573's prediction-SE path:
# 1. .hzr_safe_solve() stores a non-positive variance as an NA row and column,
#    the same mark a fixed parameter carries. The sandwich dropped such a
#    parameter as if it were known, and returned a finite, too-small SE.
# 2. An indefinite covariance with a positive diagonal passed every screen; a
#    negative quadratic form was clamped to 0 (se.fit 0, a zero-width CI), and
#    a positive one was reported though it is not a variance.

.p586_weibull <- function() {
  set.seed(1)
  n <- 300
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.5)
  tt <- stats::rweibull(n, shape = 1.5, scale = 5 * exp(-0.3 * x1))
  cc <- stats::runif(n, 0, 10)
  dat <- data.frame(time = pmin(tt, cc), status = as.integer(tt <= cc),
                    x1 = x1, x2 = x2)
  hazard(survival::Surv(time, status) ~ x1 + x2, data = dat, dist = "weibull",
         theta = c(mu = 0.2, nu = 1, x1 = 0, x2 = 0), fit = TRUE)
}

.p586_nd <- data.frame(time = c(1, 3, 6), x1 = c(0.5, -1, 1), x2 = c(1, 0, 1))

# Run a predict() call; return its value and the warnings it raised.
.p586_run <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(val = val, msgs = msgs)
}

test_that("a masked Weibull variance withholds the SEs that read it (#586)", {
  fw <- .p586_weibull()
  base <- predict(fw, newdata = .p586_nd, type = "survival", se.fit = TRUE)
  expect_true(all(is.finite(base$se.fit) & base$se.fit > 0))

  # mu's variance masked to NA, as .hzr_safe_solve() leaves a non-positive one.
  # A Weibull fit fixes nothing, so this NA is a mask, not a fixed parameter.
  # On main e7b621c5 the SEs came back as 0.0221 / 0.0221 / 0.0524 against
  # the 0.0175 / 0.0345 / 0.0392 of the fit itself, with no warning.
  m <- fw
  m$fit$vcov[1, ] <- NA_real_
  m$fit$vcov[, 1] <- NA_real_
  got <- .p586_run(predict(m, newdata = .p586_nd, type = "survival",
                           se.fit = TRUE))
  expect_true(all(is.na(got$val$se.fit)))
  expect_true(all(is.na(got$val$lower) & is.na(got$val$upper)))
  expect_length(got$msgs, 1L)
  expect_match(got$msgs, "no variance for an estimated parameter", fixed = TRUE)
  expect_match(got$msgs, "mu", fixed = TRUE)
  # The point estimate is untouched.
  expect_equal(got$val$fit, base$fit)

  # Known negative: the relative hazard never reads mu, so its SEs stand.
  hz_base <- predict(fw, newdata = .p586_nd, type = "hazard", se.fit = TRUE)
  hz <- .p586_run(predict(m, newdata = .p586_nd, type = "hazard",
                          se.fit = TRUE))
  expect_length(hz$msgs, 0L)
  expect_equal(hz$val$se.fit, hz_base$se.fit, tolerance = 1e-12)

  # A masked coefficient, read on a row where its covariate is nonzero: main
  # returned an SE of exactly 0 on row 2 (x1 = -1, x2 = 0), with no warning.
  b <- fw
  b$fit$vcov[3, ] <- NA_real_
  b$fit$vcov[, 3] <- NA_real_
  hz3 <- .p586_run(predict(b, newdata = .p586_nd, type = "hazard",
                           se.fit = TRUE))
  expect_true(all(is.na(hz3$val$se.fit)))
  expect_match(hz3$msgs, "no variance for an estimated parameter", fixed = TRUE)
})

test_that("a real multiphase fit with masked free variances withholds its SEs (#586)", {
  set.seed(494)
  d <- data.frame(time = stats::rexp(300, 0.2),
                  status = rep(c(1, 1, 0), length.out = 300),
                  mal = rep(c(0, 1), each = 150))
  fit <- suppressWarnings(hazard(
    time = d$time, status = d$status, data = d, dist = "multiphase",
    phases = list(hzr_phase("cdf", t_half = 0.15, nu = 1, m = 1,
                            formula = ~ mal),
                  hzr_phase("constant")),
    theta = c(log(0.2), log(0.15), 1, 1, 0, log(5e-04)), fit = TRUE))
  # Premise: nothing is fixed, yet some variances are NA -- masked by the
  # Hessian inversion -- and the fit records that its SEs are incomplete.
  expect_false(any(fit$fit$fixed_mask))
  masked <- which(is.na(diag(fit$fit$vcov)))
  expect_gt(length(masked), 0L)
  expect_true("standard_errors" %in% fit$degraded)

  nd <- data.frame(time = c(1, 6, 12), mal = 1)
  got <- .p586_run(predict(fit, newdata = nd, type = "cumulative_hazard",
                           se.fit = TRUE))
  # On main e7b621c5 these were finite, built from the unmasked parameters
  # alone, with no warning.
  expect_true(all(is.na(got$val$se.fit)))
  expect_true(any(grepl("no variance for an estimated parameter", got$msgs,
                        fixed = TRUE)))
  # Decomposed, the masked parameters are all the early phase's, so the early
  # phase and the total are withheld; the constant phase reads only its own
  # log_mu, whose variance is there, and keeps its SE.
  expect_true(all(masked <= 5L))
  dec <- .p586_run(predict(fit, newdata = nd, type = "cumulative_hazard",
                           se.fit = TRUE, decompose = TRUE))
  se_by <- split(dec$val$se.fit, dec$val$component)
  expect_true(all(is.na(se_by[[1]])))   # total
  expect_true(all(is.na(se_by[[2]])))   # early
  expect_true(all(is.finite(se_by[[3]]) & se_by[[3]] > 0))
  expect_length(dec$msgs[grepl("no variance", dec$msgs, fixed = TRUE)], 1L)
})

test_that("an indefinite covariance withholds the SEs that read it (#586)", {
  fw <- .p586_weibull()
  v <- fw$fit$vcov
  # Positive diagonal, indefinite matrix: the x1-x2 covariance is three times
  # what a correlation of 1 allows.
  ind <- fw
  c34 <- -3 * sqrt(v[3, 3] * v[4, 4])
  ind$fit$vcov[3, 4] <- ind$fit$vcov[4, 3] <- c34
  expect_lt(min(eigen(ind$fit$vcov[3:4, 3:4], only.values = TRUE)$values), 0)

  # On main e7b621c5 rows 1 and 3 came back with se.fit 0 and a zero-width
  # interval, and row 2 with a finite SE, all without a warning.
  lp <- .p586_run(predict(ind, newdata = .p586_nd, type = "linear_predictor",
                          se.fit = TRUE))
  expect_true(all(is.na(lp$val$se.fit)))
  expect_length(lp$msgs, 1L)
  expect_match(lp$msgs, "not positive definite over the parameters",
               fixed = TRUE)

  # Indefinite but with every row's quadratic form positive: still not a
  # variance, so still withheld.
  nd_pos <- data.frame(time = 1, x1 = 1, x2 = 0.01)
  q <- as.numeric(c(1, 0.01) %*% ind$fit$vcov[3:4, 3:4] %*% c(1, 0.01))
  expect_gt(q, 0)   # premise: main would have reported sqrt(q)
  pos <- .p586_run(predict(ind, newdata = nd_pos, type = "linear_predictor",
                           se.fit = TRUE))
  expect_true(is.na(pos$val$se.fit))

  # An infinite covariance is not a covariance either.
  inf <- fw
  inf$fit$vcov[1, 3] <- inf$fit$vcov[3, 1] <- Inf
  si <- .p586_run(predict(inf, newdata = .p586_nd, type = "survival",
                          se.fit = TRUE))
  expect_true(all(is.na(si$val$se.fit)))
  expect_length(si$msgs, 1L)
})

test_that("indefiniteness or an infinite covariance outside the prediction does not withhold it (#586, #587)", {
  fw <- .p586_weibull()
  v <- fw$fit$vcov
  ind <- fw
  ind$fit$vcov[3, 4] <- ind$fit$vcov[4, 3] <- -3 * sqrt(v[3, 3] * v[4, 4])
  # x2 is 0 in every row, so the linear predictor reads only x1's block,
  # which is fine: se = |x1| * se(b1), exactly.
  nd0 <- data.frame(time = c(1, 3), x1 = c(0.5, -1), x2 = 0)
  lp <- .p586_run(predict(ind, newdata = nd0, type = "linear_predictor",
                          se.fit = TRUE))
  expect_length(lp$msgs, 0L)
  expect_equal(lp$val$se.fit, abs(nd0$x1) * sqrt(v[3, 3]), tolerance = 1e-12)

  # #587 item 1: an infinite covariance between mu and nu, neither of which
  # the relative hazard reads, gave se.fit NaN with no warning on main.
  inf <- fw
  inf$fit$vcov[1, 2] <- inf$fit$vcov[2, 1] <- -Inf
  base <- predict(fw, newdata = .p586_nd, type = "hazard", se.fit = TRUE)
  hz <- .p586_run(predict(inf, newdata = .p586_nd, type = "hazard",
                          se.fit = TRUE))
  expect_length(hz$msgs, 0L)
  expect_equal(hz$val$se.fit, base$se.fit, tolerance = 1e-12)
})

test_that("an NA variance on a fixed parameter is still dropped as known (#586)", {
  set.seed(17)
  df <- data.frame(time = stats::rexp(60, 0.3),
                   status = stats::rbinom(60, 1, 0.6))
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = df, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                    fixed = "shapes"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 2)))
  fixed <- which(fit$fit$fixed_mask)
  expect_length(fixed, 3L)
  expect_true(all(is.na(diag(fit$fit$vcov)[fixed])))
  got <- .p586_run(predict(fit, newdata = data.frame(time = c(0.5, 2)),
                           type = "cumulative_hazard", se.fit = TRUE))
  expect_length(got$msgs, 0L)
  expect_true(all(is.finite(got$val$se.fit) & got$val$se.fit > 0))
})

test_that("a failed CoE variance recompute warns once and keeps the SEs (#586)", {
  set.seed(17)
  df <- data.frame(time = stats::rexp(60, 0.3),
                   status = stats::rbinom(60, 1, 0.6))
  # With the shapes fixed, the search estimates one log_mu and conserves the
  # other; the recompute inverts a 2 x 2 Hessian over both. Fail only that
  # inversion.
  real_solve <- .hzr_safe_solve
  testthat::local_mocked_bindings(
    .hzr_safe_solve = function(H, ...) {
      if (is.matrix(H) && nrow(H) == 2L) {
        return(list(vcov = NA, rcond = NA_real_, pd = NA,
                    reason = "mocked recompute failure"))
      }
      real_solve(H, ...)
    }
  )
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = df, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                    fixed = "shapes"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(conserve = TRUE, n_starts = 1)))
  # Premise: CoE ran, its recompute failed, and the reduced vcov stands with
  # the conserved log_mu held fixed beside the three shapes.
  expect_true(fit$spec$control$conserve_applied)
  expect_true("conserved_phase_variance" %in% fit$degraded)
  expect_true(is.matrix(fit$fit$vcov))
  expect_equal(sum(fit$fit$fixed_mask), 4L)

  nd <- data.frame(time = c(0.5, 2, 6))
  got <- .p586_run(predict(fit, newdata = nd, type = "cumulative_hazard",
                           se.fit = TRUE))
  expect_true(all(is.finite(got$val$se.fit) & got$val$se.fit > 0))
  expect_length(got$msgs, 1L)
  expect_match(got$msgs, "conserved phase", fixed = TRUE)
  dec <- .p586_run(predict(fit, newdata = nd, type = "cumulative_hazard",
                           se.fit = TRUE, decompose = TRUE))
  expect_length(dec$msgs, 1L)
  expect_match(dec$msgs, "conserved phase", fixed = TRUE)
  # Known negative: no warning without se.fit.
  expect_length(.p586_run(predict(fit, newdata = nd,
                                  type = "cumulative_hazard"))$msgs, 0L)
})
