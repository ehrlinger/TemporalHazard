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

})

test_that("NULL needs every estimated parameter to have been examined (#570)", {
  abc <- c("a", "b", "c")
  three <- wd_vcov(c(NA, 1, 4))
  # `a` was estimated and its variance failed. It is dropped, nothing is
  # found between b and c, and that is not a clean result: `a` was never
  # looked at, and is the likeliest cause of the ill-conditioning.
  r <- .hzr_weak_direction_impl(three, rcond = 1e-10, param_names = abc)
  expect_identical(r$weak, NA)
  expect_match(r$reason, "not every direction was examined", fixed = TRUE)
  r <- .hzr_weak_direction_impl(three, rcond = 1e-10, param_names = abc,
                                fixed_mask = c(FALSE, FALSE, FALSE))
  expect_identical(r$weak, NA)

  # The same matrix with `a` held FIXED: its NA row is by design, b and c
  # are all there is to examine, and nothing found is a conclusion.
  r <- .hzr_weak_direction_impl(three, rcond = 1e-10, param_names = abc,
                                fixed_mask = c(TRUE, FALSE, FALSE))
  expect_null(r$weak)
  expect_true(is.na(r$reason))

  # A ridge found among the examined parameters stands, whatever was dropped.
  ridge <- matrix(NA_real_, 3, 3)
  ridge[2:3, 2:3] <- matrix(c(1, 0.999, 0.999, 1), 2, 2)
  w <- .hzr_weak_direction(ridge, rcond = 1e-10, param_names = abc)
  expect_true(is.list(w))
  expect_setequal(w$params, c("b", "c"))
})

test_that("hazard() tells a fixed parameter's row from a failed one (#570)", {
  skip_on_cran() # a multiphase fit
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  # Shapes held fixed, so their rows are NA by design. The Hessian is
  # reported ill-conditioned, which sends the fit through the examination.
  orig <- .hzr_safe_solve
  local_mocked_bindings(.hzr_safe_solve = function(H, ...) {
    r <- orig(H, ...)
    r$rcond <- 1e-12
    r
  })
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  v <- diag(fit$fit$vcov)
  fixed <- fit$fit$fixed_mask
  # The premise: fixed rows carry NA, and every estimated one has a variance.
  expect_true(any(fixed))
  expect_true(all(is.na(v[fixed])))
  expect_true(all(is.finite(v[!fixed]) & v[!fixed] > 0))
  expect_lt(fit$fit$rcond, 1e-8)
  # So the fit was examined in full, and the fixed rows do not make it NA.
  expect_false(.hzr_is_na_scalar(fit$fit$weak))
  expect_false("weak_direction_check" %in% names(fit$degraded_causes))
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
