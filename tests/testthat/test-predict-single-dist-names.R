# predict() on a single-distribution fit returns unnamed values (#309).
# The shape parameters are named elements of theta, and R carries such a
# name onto arithmetic through rep() and through any length-1 operand: a
# lognormal fit named every survival and cumulative-hazard value "mu", and
# every family did so for one row of newdata. `expected` was produced by
# predict() on main at e4059a5, before the fix, so these tests also check
# that the fix does not move a number.

names_d <- data.frame(
  time = c(1, 2, 3, 4, 5, 6, 7, 8),
  status = c(1, 0, 1, 1, 0, 1, 0, 1),
  age = c(-1, 0.5, 0, 1.5, -0.3, 0.8, -1.2, 0.2)
)
names_shape <- list(
  weibull = c(mu = 0.2, nu = 1.3),
  exponential = c(mu = log(0.2)),
  loglogistic = c(mu = log(0.1), nu = 0.2),
  lognormal = c(mu = 1.5, sigma = -0.1)
)
names_nd <- data.frame(time = c(0.5, 2, 9), age = c(-0.5, 0, 1.2))
names_types <- c("survival", "cumulative_hazard", "hazard",
                 "linear_predictor")

names_expected <- list(
  weibull.nocov.survival =
    c(0.95111649804834897, 0.73796187384489198, 0.11682130584749517),
  weibull.nocov.cumulative_hazard =
    c(0.050118723362727227, 0.30386311717294956, 2.1471098121453833),
  weibull.nocov.hazard = c(1, 1, 1),
  weibull.nocov.linear_predictor = c(0, 0, 0),
  weibull.cov.survival =
    c(0.9577796049265147, 0.73796187384489198, 0.046073577474172824),
  weibull.cov.cumulative_hazard =
    c(0.043137584966540223, 0.30386311717294956, 3.0775156500391043),
  weibull.cov.hazard = c(0.86070797642505781, 1, 1.4333294145603401),
  weibull.cov.linear_predictor =
    c(-0.14999999999999999, 0, 0.35999999999999999),
  exponential.nocov.survival =
    c(0.90483741803595952, 0.67032004603563933, 0.16529888822158653),
  exponential.nocov.cumulative_hazard =
    c(0.10000000000000001, 0.40000000000000002, 1.8),
  exponential.nocov.hazard = c(1, 1, 1),
  exponential.nocov.linear_predictor = c(0, 0, 0),
  exponential.cov.survival =
    c(0.91752927001138107, 0.67032004603563933, 0.075774538518747589),
  exponential.cov.cumulative_hazard =
    c(0.086070797642505789, 0.40000000000000002, 2.5799929462086122),
  exponential.cov.hazard = c(0.86070797642505781, 1, 1.4333294145603401),
  exponential.cov.linear_predictor =
    c(-0.14999999999999999, 0, 0.35999999999999999),
  loglogistic.nocov.survival =
    c(0.95887706856295307, 0.81091601574655547, 0.40585894941602102),
  loglogistic.nocov.cumulative_hazard =
    c(0.041992399423863604, 0.20959078664514447, 0.90174959497453222),
  loglogistic.nocov.hazard = c(1, 1, 1),
  loglogistic.nocov.linear_predictor = c(0, 0, 0),
  loglogistic.cov.survival =
    c(0.96440125183927583, 0.81091601574655547, 0.32276123199947904),
  loglogistic.cov.cumulative_hazard =
    c(0.036247834626113078, 0.20959078664514447, 1.1308424489330255),
  loglogistic.cov.hazard = c(0.86070797642505781, 1, 1.4333294145603401),
  loglogistic.cov.linear_predictor =
    c(-0.14999999999999999, 0, 0.35999999999999999),
  lognormal.nocov.survival =
    c(0.99232052174280183, 0.81372587662777718, 0.22048616411531335),
  lognormal.nocov.cumulative_hazard =
    c(0.0077091172893927404, 0.2061317305958188, 1.5119203338403171),
  lognormal.nocov.hazard = c(1, 1, 1),
  lognormal.nocov.linear_predictor = c(0, 0, 0),
  lognormal.cov.survival =
    c(0.98802800709118632, 0.81372587662777718, 0.35468929383240411),
  lognormal.cov.cumulative_hazard =
    c(0.012044234377759399, 0.2061317305958188, 1.0365131013840527),
  lognormal.cov.hazard = c(0.86070797642505781, 1, 1.4333294145603401),
  lognormal.cov.linear_predictor =
    c(-0.14999999999999999, 0, 0.35999999999999999)
)

test_that("single-distribution predict() returns unnamed values (#309)", {
  n_checked <- 0L
  for (dist in names(names_shape)) {
    for (cov in c(FALSE, TRUE)) {
      f <- hazard(
        if (cov) survival::Surv(time, status) ~ age else
          survival::Surv(time, status) ~ 1,
        data = names_d, dist = dist,
        theta = c(names_shape[[dist]], if (cov) c(age = 0.3)), fit = FALSE
      )
      nd <- if (cov) names_nd else names_nd["time"]
      for (type in names_types) {
        key <- paste(dist, if (cov) "cov" else "nocov", type, sep = ".")
        for (rows in list(1L, 1:3)) {
          lbl <- paste0(key, ", ", length(rows), " row(s)")
          nd_rows <- nd[rows, , drop = FALSE]

          p <- predict(f, newdata = nd_rows, type = type)
          expect_null(names(p), label = paste("names of", lbl))
          expect_equal(p, names_expected[[key]][rows], tolerance = 1e-12,
                       label = lbl)

          # The se.fit data frame took the name as its row name. An
          # unfitted model has no vcov, so its SEs are NA, with a warning.
          expect_warning(
            s <- predict(f, newdata = nd_rows, type = type, se.fit = TRUE),
            "Variance-covariance matrix is unavailable"
          )
          expect_identical(rownames(s), as.character(seq_along(rows)),
                           label = paste("se.fit row names of", lbl))
          expect_equal(s$fit, names_expected[[key]][rows], tolerance = 1e-12,
                       label = paste("se.fit fit of", lbl))
          n_checked <- n_checked + 1L
        }
      }
    }
  }
  # 4 families x 2 designs x 4 types x 2 row counts: the grid ran in full.
  expect_identical(n_checked, 64L)
})

test_that("a fitted lognormal predicts unnamed values (#309 reproducer)", {
  set.seed(1)
  d <- data.frame(time = rexp(200, 0.1), status = rbinom(200, 1, 0.7))
  f <- hazard(survival::Surv(time, status) ~ 1, data = d, dist = "lognormal",
              theta = c(mu = 2, sigma = 1), fit = TRUE)
  expect_true(isTRUE(f$fit$converged))
  nd <- data.frame(time = c(2, 5))
  for (type in c("survival", "cumulative_hazard")) {
    p <- predict(f, newdata = nd, type = type)
    expect_length(p, 2L)
    expect_null(names(p), label = type)
  }
})
