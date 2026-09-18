# A backward screen whose removal tests cannot be computed (#389).
#
# Removal is tested on the current model's Wald p-value, which needs a
# variance for the coefficient. A multiphase fit with interval-censored rows
# gets that variance from numDeriv (a Suggests), so without it every Wald p
# is NA. The backward step dropped NA-scored candidates from `eligible` and
# returned the same empty result as "nothing met slstay": the screen stopped
# after 0 steps with no warning and n_uncomputable_scores = 0.

.ic_backward_data <- function() {
  data(avc, envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
  set.seed(3)
  ic <- sample(which(d$dead == 1), 40)
  d$lo <- d$int_dead
  d$hi <- d$int_dead
  d$lo[ic] <- d$int_dead[ic] * 0.7
  d$hi[ic] <- d$int_dead[ic] * 1.3
  d$st <- ifelse(d$dead == 1, 1L, 0L)
  d$st[ic] <- 2L
  d
}

.ic_backward_fit <- function(d) {
  suppressWarnings(hazard(
    time = d$hi, status = d$st, time_lower = d$lo, time_upper = d$hi,
    data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ age + mal + com_iv),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L)
  ))
}

.mask_numderiv <- function(env = parent.frame()) {
  orig <- base::requireNamespace
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "numDeriv")) FALSE else orig(package, ...)
    },
    .package = "base", .env = env
  )
}

test_that("control: with variances available the backward screen drops a variable", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- .ic_backward_data()
  fit <- .ic_backward_fit(d)
  expect_true(any(is.finite(fit$fit$se)))

  sw <- hzr_stepwise(fit, data = d, direction = "backward",
                     criterion = "wald", slstay = 0.05, trace = FALSE)
  expect_gt(nrow(sw$steps), 0L)
  expect_true("age" %in% sw$steps$variable[sw$steps$action == "drop"])
  expect_identical(sw$criteria$n_uncomputable_scores, 0L)
  expect_false(sw$criteria$stopped_uncomputable)
})

test_that("a backward screen with no Wald p-values warns and counts them", {
  skip_on_cran()
  .mask_numderiv()
  d <- .ic_backward_data()
  fit <- .ic_backward_fit(d)
  expect_false(any(is.finite(fit$fit$se)))

  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(fit, data = d, direction = "backward",
                       criterion = "wald", slstay = 0.05, trace = FALSE)
  )
  expect_identical(nrow(sw$steps), 0L)
  # All three in-model candidates were NA, and all three are counted.
  expect_identical(sw$criteria$n_uncomputable_scores, 3L)
  expect_identical(sw$criteria$uncomputable_reasons,
                   c(wald_no_variance = 3L))
  expect_true(sw$criteria$stopped_uncomputable)
  expect_true(any(grepl("could not be computed", w)))
  expect_true(any(grepl("slstay", w)))
})

test_that("a forced-in candidate with no variance is not counted", {
  skip_on_cran()
  .mask_numderiv()
  d <- .ic_backward_data()
  fit <- .ic_backward_fit(d)

  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(fit, data = d, direction = "backward",
                       criterion = "wald", slstay = 0.05, trace = FALSE,
                       force_in = c("age", "mal", "com_iv"))
  )
  # Nothing was a removal candidate, so nothing went untested.
  expect_identical(sw$criteria$n_uncomputable_scores, 0L)
  expect_false(sw$criteria$stopped_uncomputable)
  expect_false(any(grepl("could not be computed", w)))
})

test_that("a variable kept only because its removal test was NA is reported", {
  obj <- .fit_overfitted()
  orig <- .hzr_candidate_score
  local_mocked_bindings(
    .hzr_candidate_score = function(...) {
      a <- list(...)
      s <- orig(...)
      # A single-distribution drop names its coefficient by position
      # (`beta<k>`), so which variable it is depends on the step.
      cols <- colnames(a$current$data$x)
      var <- cols[match(a$names, paste0("beta", seq_along(cols)))]
      if (identical(var, "x3")) {
        s$score <- NA_real_
        s$p_value <- NA_real_
        s$stat <- NA_real_
      }
      s
    }
  )
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(obj$fit, data = obj$data, direction = "backward",
                       criterion = "wald", slstay = 0.20, trace = FALSE)
  )
  # The screen ran (x2 can be tested), so it did not stop uncomputable ...
  expect_false(sw$criteria$stopped_uncomputable)
  expect_false("x3" %in% sw$steps$variable)
  expect_true("x3" %in% colnames(sw$data$x))
  # ... but x3 was kept without a test, once per backward step taken.
  expect_identical(sw$criteria$uncomputable_reasons[["wald_no_variance"]],
                   sw$criteria$n_uncomputable_scores)
  expect_identical(sw$criteria$n_uncomputable_scores,
                   nrow(sw$steps) + 1L)
  expect_true(any(grepl("kept without being tested", w)))
})
