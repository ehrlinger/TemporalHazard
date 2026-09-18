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
  # The stop warning covers them; they are not named a second time.
  expect_false(any(grepl("without a Wald test", w)))
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

# Removal tests for x3 only come back NA. A single-distribution drop names its
# coefficient by position (`beta<k>`), so which variable it is depends on the
# step: resolve it through the current fit's design columns.
.mask_x3_removal <- function(env = parent.frame(), only_while = NULL) {
  orig <- .hzr_candidate_score
  local_mocked_bindings(
    .hzr_candidate_score = function(...) {
      a <- list(...)
      s <- orig(...)
      if (identical(a$mode, "drop")) {
        cols <- colnames(a$current$data$x)
        var <- cols[match(a$names, paste0("beta", seq_along(cols)))]
        if (identical(var, "x3") &&
              (is.null(only_while) || only_while %in% cols)) {
          s$score <- NA_real_
          s$p_value <- NA_real_
          s$stat <- NA_real_
        }
      }
      s
    },
    .env = env
  )
}

test_that("a variable kept only because its removal test was NA is reported", {
  obj <- .fit_overfitted()
  .mask_x3_removal()
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(obj$fit, data = obj$data, direction = "backward",
                       criterion = "wald", slstay = 0.20, trace = FALSE)
  )
  # The screen ran (x2 can be tested), so it did not stop uncomputable ...
  expect_false(sw$criteria$stopped_uncomputable)
  expect_false("x3" %in% sw$steps$variable)
  expect_true("x3" %in% colnames(sw$data$x))
  # ... but x3 was kept without a test. The counter counts attempts, one per
  # backward step; the warning names the variable once.
  expect_identical(sw$criteria$uncomputable_reasons[["wald_no_variance"]],
                   sw$criteria$n_uncomputable_scores)
  expect_identical(sw$criteria$n_uncomputable_scores,
                   nrow(sw$steps) + 1L)
  expect_identical(sum(grepl("without a Wald test", w)), 1L)
  expect_true(any(grepl(
    "decided 1 variable\\(s\\) without a Wald test; kept in the model with its removal untested: x3\\.",
    w
  )))
})

test_that("a removal untested at one step and tested at a later one is not reported", {
  # x3's test is NA only while x2 is still in the model; once x2 is dropped
  # x3 is tested, so the screen decided it on a test.
  obj <- .fit_overfitted()
  .mask_x3_removal(only_while = "x2")
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(obj$fit, data = obj$data, direction = "backward",
                       criterion = "wald", slstay = 0.20, trace = FALSE)
  )
  expect_identical(sw$steps$variable[1], "x2")
  expect_gt(nrow(sw$steps), 0L)
  expect_identical(sw$criteria$n_uncomputable_scores, 1L)
  expect_false(any(grepl("without a Wald test", w)))
})

test_that("a two-way stop names the half that could not be tested", {
  # The entry is tested and rejected (noise x2, unreachable slentry) and the
  # only removal cannot be tested: the screen stopped for want of a removal
  # test, and must not claim its entries went untested.
  obj <- .fit_overfitted()
  base <- hazard(Surv(time, status) ~ x3, data = obj$data,
                 theta = c(0.5, 1.0, 0), dist = "weibull", fit = TRUE)
  .mask_x3_removal()
  msgs <- character()
  out <- utils::capture.output(sw <- withCallingHandlers(
    hzr_stepwise(base, scope = "x2", data = obj$data,
                 direction = "both", criterion = "wald", slentry = 1e-12,
                 slstay = 0.20, trace = TRUE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  msgs <- c(msgs, out)
  expect_identical(nrow(sw$steps), 0L)
  expect_true(sw$criteria$stopped_uncomputable)
  expect_identical(sw$criteria$wald_untested_removals, "x3")
  expect_identical(sw$criteria$wald_untested_entries, character())
  expect_true(any(grepl("could be COMPUTED for removal -- none was tested",
                        msgs)))
  expect_true(any(grepl("could be tested for removal: .*could not be computed",
                        msgs)))
  expect_false(any(grepl("for entry", msgs)))
})

test_that("a two-way screen that recovers is not reported as stopped", {
  # x3's entry test is NA on the first forward step only. The screen goes on
  # (x2 is dropped), tests x3 at the next iteration and stops on merit.
  obj <- .fit_overfitted()
  base <- hazard(Surv(time, status) ~ x1 + x2, data = obj$data,
                 theta = c(0.5, 1.0, 0, 0), dist = "weibull", fit = TRUE)
  orig <- .hzr_candidate_score
  seen <- new.env()
  seen$n <- 0L
  local_mocked_bindings(
    .hzr_candidate_score = function(...) {
      a <- list(...)
      s <- orig(...)
      if (identical(a$mode, "entry") && seen$n == 0L) {
        seen$n <- 1L
        s$score <- NA_real_
        s$p_value <- NA_real_
        s$stat <- NA_real_
      }
      s
    }
  )
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(base, scope = c("x1", "x2", "x3"), data = obj$data,
                       direction = "both", criterion = "wald",
                       slentry = 0.01, slstay = 0.20, trace = FALSE)
  )
  expect_identical(seen$n, 1L)
  expect_identical(sw$steps$variable, "x2")
  expect_identical(sw$criteria$uncomputable_reasons[["wald_no_variance"]], 1L)
  expect_false(sw$criteria$stopped_uncomputable)
  expect_identical(sw$criteria$wald_untested_entries, character())
  expect_false(any(grepl("could not be computed|without a Wald test", w)))
})

test_that("a forward Wald screen with no variances warns and counts them", {
  skip_on_cran()
  .mask_numderiv()
  d <- .ic_backward_data()
  fit <- suppressWarnings(hazard(
    time = d$hi, status = d$st, time_lower = d$lo, time_upper = d$hi,
    data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L)
  ))
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(fit, data = d, scope = list(early = ~ mal + com_iv),
                       direction = "forward", criterion = "wald",
                       slentry = 0.05, trace = FALSE)
  )
  expect_identical(nrow(sw$steps), 0L)
  expect_identical(sw$criteria$n_uncomputable_scores, 2L)
  expect_identical(sw$criteria$uncomputable_reasons,
                   c(wald_no_variance = 2L))
  expect_true(sw$criteria$stopped_uncomputable)
  expect_true(any(grepl("could be tested for entry: .*could not be computed",
                        w)))
  expect_false(any(grepl("without a Wald test", w)))
})

test_that("control: with variances the same forward Wald screen enters com_iv", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- .ic_backward_data()
  fit <- suppressWarnings(hazard(
    time = d$hi, status = d$st, time_lower = d$lo, time_upper = d$hi,
    data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L)
  ))
  sw <- suppressWarnings(
    hzr_stepwise(fit, data = d, scope = list(early = ~ mal + com_iv),
                 direction = "forward", criterion = "wald",
                 slentry = 0.05, trace = FALSE)
  )
  expect_true("com_iv" %in% sw$steps$variable[sw$steps$action == "enter"])
  expect_identical(sw$criteria$n_uncomputable_scores, 0L)
})

test_that("hzr_bootstrap() warns when replicates decide a variable untested", {
  obj <- .fit_overfitted()
  .mask_x3_removal()
  w <- testthat::capture_warnings(
    boot <- hzr_bootstrap(obj$fit, n_boot = 3, seed = 1,
                          scope = c("x1", "x2", "x3"),
                          direction = "backward", criterion = "wald",
                          slstay = 0.20)
  )
  expect_identical(boot$n_uncomputable_replicates, 0L)
  expect_gt(boot$uncomputable_reasons[["wald_no_variance"]], 0L)
  expect_true(any(grepl(
    paste0("^", boot$n_success, " of ", boot$n_success,
           " successful replicates decided a variable without a Wald test"),
    w
  )))
})

test_that("a forward candidate whose refit failed is a refit failure, not untested", {
  obj <- .fit_stepwise_base()
  local_mocked_bindings(.hzr_refit_with_scope = function(...) stop("boom"))
  step <- suppressWarnings(
    .hzr_stepwise_forward_step(obj$fit, scope = c("x1", "x2"), data = obj$data,
                               criterion = "wald", slentry = 0.05)
  )
  expect_identical(step$refit_failures, c("x1", "x2"))
  expect_identical(step$n_uncomputable, 0L)
  expect_null(step$stop_reason)
})

test_that("a forward entry left out only because its Wald test was NA is reported", {
  obj <- .fit_stepwise_base()
  orig <- .hzr_candidate_score
  local_mocked_bindings(
    .hzr_candidate_score = function(...) {
      a <- list(...)
      s <- orig(...)
      if (identical(a$mode, "entry")) {
        cols <- colnames(a$candidate$data$x)
        var <- cols[match(a$names, paste0("beta", seq_along(cols)))]
        if (identical(var, "x3")) {
          s$score <- NA_real_
          s$p_value <- NA_real_
          s$stat <- NA_real_
        }
      }
      s
    }
  )
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(obj$fit, scope = c("x1", "x2", "x3"), data = obj$data,
                       direction = "forward", criterion = "wald",
                       slentry = 0.05, trace = FALSE)
  )
  # x1 is tested and enters, and noise x2 is tested and rejected at the last
  # step, so the run did not stop for want of a test ...
  expect_identical(sw$steps$variable, "x1")
  expect_false(sw$criteria$stopped_uncomputable)
  # ... but x3 was left out untested, and is named once.
  expect_identical(sw$criteria$uncomputable_reasons[["wald_no_variance"]], 2L)
  expect_identical(sum(grepl("without a Wald test", w)), 1L)
  expect_true(any(grepl("left out with its entry untested: x3\\.", w)))
})

test_that("a failed entry refit is not reported as a variable left out untested", {
  # x2's refit fails, which is reported as a refit failure with its own
  # reason; it must not also be named with the no-variance diagnosis.
  obj <- .fit_stepwise_base()
  orig <- .hzr_refit_with_scope
  local_mocked_bindings(
    .hzr_refit_with_scope = function(current, action, var, ...) {
      if (identical(action, "add") && identical(var, "x2")) stop("boom")
      orig(current, action = action, var = var, ...)
    }
  )
  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(obj$fit, scope = c("x1", "x2"), data = obj$data,
                       direction = "forward", criterion = "wald",
                       slentry = 0.05, trace = FALSE)
  )
  expect_true("x2" %in% sw$criteria$refit_failures)
  expect_true(any(grepl("candidate refit failed for x2: boom", w)))
  expect_false(any(grepl("without a Wald test", w)))
})

test_that("hzr_bootstrap() counts only replicates that left a variable untested", {
  # x3's removal test is NA only while x2 is in the model. A replicate that
  # drops x2 then tests x3 decided it on a test; one that keeps x2 did not.
  # The bootstrap's count must match the replicates' own screens.
  obj <- .fit_overfitted()
  .mask_x3_removal(only_while = "x2")
  screens <- list()
  orig_sw <- hzr_stepwise
  local_mocked_bindings(
    hzr_stepwise = function(...) {
      r <- orig_sw(...)
      screens[[length(screens) + 1L]] <<- list(
        listed = length(c(r$criteria$wald_untested_removals,
                          r$criteria$wald_untested_entries)) > 0L,
        na = r$criteria$n_uncomputable_scores > 0L
      )
      r
    }
  )
  w <- testthat::capture_warnings(
    boot <- hzr_bootstrap(obj$fit, n_boot = 3, seed = 1,
                          scope = c("x1", "x2", "x3"),
                          direction = "backward", criterion = "wald",
                          slstay = 0.20)
  )
  expect_identical(boot$n_failed, 0L)
  reps <- utils::tail(screens, boot$n_success)
  n_listed <- sum(vapply(reps, `[[`, logical(1L), "listed"))
  n_na <- sum(vapply(reps, `[[`, logical(1L), "na"))
  # At least one replicate had an NA test and then tested the variable.
  expect_gt(n_na, n_listed)
  hit <- grepl("successful replicates decided a variable without a Wald test",
               w)
  if (n_listed > 0L) {
    expect_true(any(grepl(paste0("^", n_listed, " of ", boot$n_success,
                                 " successful replicates decided"), w)))
  } else {
    expect_false(any(hit))
  }
})
