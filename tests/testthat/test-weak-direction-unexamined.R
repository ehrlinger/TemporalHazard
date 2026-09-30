# .hzr_weak_direction_impl() returns NULL when it looked and found no ridge,
# and NA when it could not look. With fewer than two finite variances it
# returned NULL. That exit is reached only under an ill-conditioned Hessian,
# and .hzr_safe_solve() masks a failed variance to NA, so it means the
# covariance mostly failed: there was nothing to decompose, and the fit was
# reported as examined and clean (#570).

wd_vcov <- function(d) {
  v <- diag(d, nrow = length(d))
  bad <- !is.finite(d)
  v[bad, ] <- NA_real_
  v[, bad] <- NA_real_
  v
}

test_that("fewer than two finite variances is NA, with its reason (#570)", {
  nms <- c("a", "b")
  one_failed <- wd_vcov(c(NA, 4))
  # Known positive: well conditioned, the same matrix is a real "no ridge".
  expect_null(.hzr_weak_direction(one_failed, rcond = 0.3, param_names = nms))

  r <- .hzr_weak_direction_impl(one_failed, rcond = 1e-10, param_names = nms)
  expect_identical(r$weak, NA)
  expect_match(r$reason, "fewer than two parameters have a finite positive",
               fixed = TRUE)
  expect_identical(.hzr_weak_direction(one_failed, rcond = 1e-10,
                                       param_names = nms), NA)
  # All of them failed.
  expect_identical(.hzr_weak_direction(wd_vcov(c(NA, NA)), rcond = 1e-10,
                                       param_names = nms), NA)
  # A non-positive variance that was not masked is a failed one too.
  expect_identical(.hzr_weak_direction(diag(c(-1, 4)), rcond = 1e-10,
                                       param_names = nms), NA)
  # Three parameters, one usable variance.
  expect_identical(
    .hzr_weak_direction(wd_vcov(c(NA, NA, 4)), rcond = 1e-10,
                        param_names = c("a", "b", "c")), NA)

  # Two usable variances are enough to examine: that is a conclusion again.
  expect_null(.hzr_weak_direction(wd_vcov(c(NA, 1, 4)), rcond = 1e-10,
                                  param_names = c("a", "b", "c")))
})

test_that("a fit whose covariance mostly failed says it was not examined (#570)", {
  d <- data.frame(tt = c(1.2, 0.4, 2.2, 3.1, 0.9, 1.7, 2.8, 0.6, 1.1, 3.6),
                  ev = c(1, 1, 0, 1, 1, 0, 1, 1, 1, 0),
                  z = c(0.1, -1.2, 0.7, 0.3, -0.4, 1.5, -0.8, 0.2, 0.9, -0.1))
  fit_with <- function() {
    suppressWarnings(hazard(survival::Surv(tt, ev) ~ z, data = d,
                            dist = "weibull", theta = c(0.5, 1, 0),
                            fit = TRUE))
  }
  # Known positive: the fit as it is was examined, and is clean.
  clean <- fit_with()
  expect_null(clean$fit$weak)
  expect_false("weak_direction_check" %in% names(clean$degraded_causes))

  # The inversion is made to fail on all but the last variance, under an
  # ill-conditioned Hessian: what .hzr_safe_solve() returns for such a fit.
  orig <- .hzr_safe_solve
  local_mocked_bindings(.hzr_safe_solve = function(H, ...) {
    r <- orig(H, ...)
    if (is.matrix(r$vcov)) {
      bad <- seq_len(nrow(r$vcov) - 1L)
      r$vcov[bad, ] <- NA_real_
      r$vcov[, bad] <- NA_real_
      r$rcond <- 1e-12
    }
    r
  })
  fit <- fit_with()
  expect_identical(sum(is.finite(diag(fit$fit$vcov))), 1L)
  expect_identical(fit$fit$weak, NA)
  expect_match(fit$degraded_causes[["weak_direction_check"]],
               "fewer than two", fixed = TRUE)
  # What the reader sees: no ridge is named, and the check is listed as not
  # done.
  out <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
  expect_no_match(out, "weakly identified", fixed = TRUE)
  expect_match(out, "fewer than two", fixed = TRUE)
})
