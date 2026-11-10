# predict(se.fit = TRUE): the "jacobian" attribute and `average = TRUE`
# (#731). Every fit here is a real optimisation, so the file is skipped on
# CRAN.
skip_on_cran()

avg_toy <- function(seed = 731, n = 80) {
  set.seed(seed)
  data.frame(time = rexp(n, 0.3), status = rbinom(n, 1, 0.7), x = rnorm(n))
}

fit_weibull <- function() {
  hazard(survival::Surv(time, status) ~ x, data = avg_toy(),
         dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
}

# Shapes fixed, and x on the constant phase only: with x on both phases too,
# early.log_mu and early.x are determined only in combination on this data.
fit_multiphase <- function() {
  hazard(survival::Surv(time, status) ~ x, data = avg_toy(n = 200),
         dist = "multiphase",
         phases = list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                         fixed = "shapes", formula = ~1),
                       constant = hzr_phase("constant")),
         fit = TRUE, control = list(n_starts = 2))
}

# The reference: the delta-method SE of `stat(predictions)` from a numerical
# gradient over the STORED parameters and the stored vcov, perturbing
# object$fit$theta and re-predicting. It shares nothing with the analytic
# path but predict() itself: a Weibull fit is perturbed on the mu scale, not
# the log(mu) scale the attribute reports, and a parameter with no variance
# (held fixed) is not perturbed.
numeric_se <- function(fit, newdata, type, stat) {
  theta <- fit$fit$theta
  V <- fit$fit$vcov
  free <- which(!is.na(diag(V)))
  f <- function(th_free) {
    obj <- fit
    obj$fit$theta[free] <- th_free
    stat(predict(obj, newdata = newdata, type = type))
  }
  g <- numDeriv::grad(f, theta[free])
  sqrt(drop(t(g) %*% V[free, free] %*% g))
}

hetero_nd <- function(times) {
  xs <- c(-1.5, -0.3, 0.4, 1.8)
  expand.grid(x = xs, time = times)
}

for (dist in c("weibull", "multiphase")) {
  make_fit <- if (dist == "weibull") fit_weibull else fit_multiphase

  test_that(paste(dist, "- one row: the average is that row"), {
    fit <- make_fit()
    nd <- data.frame(time = 2, x = 0.7)
    for (type in c("survival", "cumulative_hazard")) {
      for (ct in c("log-log", "logit")) {
        row <- predict(fit, nd, type = type, se.fit = TRUE, conf.type = ct)
        avg <- predict(fit, nd, type = type, se.fit = TRUE, conf.type = ct,
                       average = TRUE)
        expect_named(avg, c("time", "n", "fit", "se.fit", "lower", "upper"))
        expect_equal(avg$n, 1L)
        expect_true(is.finite(row$se.fit) && row$se.fit > 0)
        expect_equal(unlist(avg[c("fit", "se.fit", "lower", "upper")]),
                     unlist(row[c("fit", "se.fit", "lower", "upper")]),
                     tolerance = 1e-12)
      }
    }
  })

  test_that(paste(dist, "- n identical rows: the average is one row"), {
    fit <- make_fit()
    one <- data.frame(time = c(1, 3), x = 0.7)
    many <- one[rep(1:2, each = 6), ]
    row <- predict(fit, one, type = "survival", se.fit = TRUE)
    avg <- predict(fit, many, type = "survival", se.fit = TRUE,
                   average = TRUE)
    expect_equal(avg$time, c(1, 3))
    expect_equal(avg$n, c(6L, 6L))
    expect_equal(avg$se.fit, row$se.fit, tolerance = 1e-10)
    expect_equal(avg$lower, row$lower, tolerance = 1e-10)
  })

  test_that(paste(dist, "- heterogeneous rows: bounded and numerically right"),
            {
    fit <- make_fit()
    nd <- hetero_nd(c(0.5, 2, 5))
    for (type in c("survival", "cumulative_hazard")) {
      rows <- predict(fit, nd, type = type, se.fit = TRUE)
      avg <- predict(fit, nd, type = type, se.fit = TRUE, average = TRUE)
      # Coverage: every row, every time, nothing else.
      expect_equal(avg$time, sort(unique(nd$time)))
      expect_equal(sum(avg$n), nrow(nd))
      expect_equal(avg$n, as.integer(table(nd$time)[as.character(avg$time)]),
                   ignore_attr = TRUE)
      for (k in seq_along(avg$time)) {
        at <- nd$time == avg$time[k]
        expect_equal(avg$fit[k], mean(rows$fit[at]), tolerance = 1e-12)
        # Strictly below: the rows differ, so they are not perfectly
        # correlated, and the bound must not be met with equality.
        expect_lt(avg$se.fit[k], mean(rows$se.fit[at]))
        ref <- numeric_se(fit, nd[at, ], type, mean)
        expect_equal(avg$se.fit[k] / ref, 1, tolerance = 1e-5)
      }
      expect_true(all(avg$lower < avg$fit & avg$fit < avg$upper))
    }
  })

  test_that(paste(dist, "- the attribute gives Var(E) and the per-row SE"), {
    fit <- make_fit()
    nd <- data.frame(time = c(0.4, 1.1, 2.5, 4, 6), x = c(-1, 0, 0.5, 1, 2))
    res <- predict(fit, nd, type = "cumulative_hazard", se.fit = TRUE)
    j <- attr(res, "jacobian")
    expect_named(j, c("J", "vcov", "names"))
    expect_equal(dim(j$J), c(nrow(nd), length(j$names)))
    expect_equal(dim(j$vcov), rep(length(j$names), 2L))
    expect_equal(colnames(j$J), j$names)
    # The per-row SE is the attribute's own sandwich.
    expect_equal(unname(sqrt(rowSums((j$J %*% j$vcov) * j$J))), res$se.fit,
                 tolerance = 1e-12)
    g <- colSums(j$J)
    se_e <- sqrt(drop(t(g) %*% j$vcov %*% g))
    ref <- numeric_se(fit, nd, "cumulative_hazard", sum)
    expect_equal(se_e / ref, 1, tolerance = 1e-5)

    # Survival: J is of S itself, not of H.
    res_s <- predict(fit, nd, type = "survival", se.fit = TRUE)
    js <- attr(res_s, "jacobian")
    expect_equal(js$J, -res_s$fit * j$J, tolerance = 1e-12)
  })
}

test_that("Weibull attribute is on the log(mu) scale", {
  fit <- fit_weibull()
  res <- predict(fit, data.frame(time = 2, x = 0.5), type = "survival",
                 se.fit = TRUE)
  # summary()'s row names, with the first on the scale predictions use.
  expect_equal(rownames(summary(fit)$coefficients), c("mu", "nu", "beta1"))
  expect_equal(attr(res, "jacobian")$names, c("log(mu)", "nu", "beta1"))
  # A coefficient whose covariate is 0 in every row is not read, so it is
  # not a column.
  res0 <- predict(fit, data.frame(time = 2, x = 0), type = "survival",
                  se.fit = TRUE)
  expect_equal(attr(res0, "jacobian")$names, c("log(mu)", "nu"))
})

test_that("multiphase attribute leaves out the fixed shapes", {
  fit <- fit_multiphase()
  res <- predict(fit, data.frame(time = 2, x = 0.3),
                 type = "cumulative_hazard", se.fit = TRUE)
  nms <- attr(res, "jacobian")$names
  fixed <- names(fit$fit$theta)[fit$fit$fixed_mask %in% TRUE &
                                  is.na(diag(fit$fit$vcov))]
  expect_gt(length(fixed), 0L)
  expect_length(intersect(nms, fixed), 0L)
  expect_setequal(nms, setdiff(names(fit$fit$theta), fixed))
})

test_that("an averaged result carries one Jacobian row per time", {
  fit <- fit_weibull()
  nd <- hetero_nd(c(1, 4))
  avg <- predict(fit, nd, type = "survival", se.fit = TRUE, average = TRUE)
  j <- attr(avg, "jacobian")
  expect_equal(nrow(j$J), 2L)
  expect_equal(unname(sqrt(rowSums((j$J %*% j$vcov) * j$J))), avg$se.fit,
               tolerance = 1e-12)
})

test_that("no attribute when the vcov is unusable", {
  fit <- fit_weibull()
  fit$fit$vcov <- NULL
  expect_warning(
    res <- predict(fit, data.frame(time = 2, x = 0), type = "survival",
                   se.fit = TRUE),
    "unavailable"
  )
  expect_true(is.na(res$se.fit))
  expect_null(attr(res, "jacobian"))
  expect_warning(
    avg <- predict(fit, hetero_nd(c(1, 2)), type = "survival", se.fit = TRUE,
                   average = TRUE),
    "unavailable"
  )
  expect_equal(avg$n, c(4L, 4L))
  expect_true(all(is.na(avg$se.fit)))
  expect_null(attr(avg, "jacobian"))
})

test_that("average = TRUE refuses what it cannot average", {
  fit <- fit_weibull()
  nd <- data.frame(time = 2, x = 0)
  expect_error(predict(fit, nd, type = "survival", average = TRUE),
               "se.fit = TRUE")
  expect_error(predict(fit, nd, type = "hazard", se.fit = TRUE,
                       average = TRUE), "survival")
  expect_error(predict(fit, nd, type = "survival", se.fit = TRUE,
                       average = NA), "average")
  mp <- fit_multiphase()
  expect_error(predict(mp, data.frame(time = 2, x = 0),
                       type = "cumulative_hazard", se.fit = TRUE,
                       decompose = TRUE, average = TRUE), "decompose")
})

test_that("J rows carry the result's row names, for realignment", {
  fit <- fit_weibull()
  nd <- data.frame(time = c(0.5, 1, 2, 4), x = c(-1, 0.2, 0.9, 1.6))
  res <- predict(fit, nd, type = "cumulative_hazard", se.fit = TRUE)
  J <- attr(res, "jacobian")$J
  expect_equal(rownames(J), rownames(res))
  sub <- res[c(4, 2), ]
  Js <- J[rownames(sub), , drop = FALSE]
  V <- attr(res, "jacobian")$vcov
  expect_equal(unname(sqrt(rowSums((Js %*% V) * Js))), sub$se.fit, tolerance = 1e-12)
})

test_that("average = TRUE needs newdata", {
  fit <- fit_weibull()
  expect_error(predict(fit, type = "survival", se.fit = TRUE, average = TRUE),
               "needs 'newdata'")
})

test_that("an underflowed mean survival has NA limits, not NaN", {
  fit <- fit_weibull()
  nd <- data.frame(time = c(1e4, 1e4), x = c(3, 3))
  avg <- predict(fit, nd, type = "survival", se.fit = TRUE, average = TRUE)
  expect_equal(avg$fit, 0)
  expect_true(is.na(avg$lower) && !is.nan(avg$lower))
  expect_true(is.na(avg$upper) && !is.nan(avg$upper))
})

# The numeric-Jacobian families, checked against the same reference.
for (dist in c("exponential", "loglogistic", "lognormal")) {
  test_that(paste(dist, "- averaged SE matches the numerical gradient"), {
    theta0 <- switch(dist, exponential = c(-1, 0),
                     loglogistic = c(-1, 0, 0), lognormal = c(1, 0, 0))
    fit <- hazard(survival::Surv(time, status) ~ x, data = avg_toy(),
                  dist = dist, theta = theta0, fit = TRUE)
    nd <- hetero_nd(c(0.5, 2, 5))
    for (type in c("survival", "cumulative_hazard")) {
      avg <- predict(fit, nd, type = type, se.fit = TRUE, average = TRUE)
      expect_equal(sum(avg$n), nrow(nd))
      for (k in seq_along(avg$time)) {
        at <- nd$time == avg$time[k]
        ref <- numeric_se(fit, nd[at, ], type, mean)
        expect_gt(ref, 0)
        expect_equal(avg$se.fit[k] / ref, 1, tolerance = 1e-5)
      }
    }
  })
}
