# The score criterion refits a candidate it cannot score and tests it by Wald,
# but only for the reasons in .hzr_score_fallback_reasons. Two reasons were
# left out as "no refit can make those candidates testable", and that is wrong
# for both (#570):
#
#   nuisance_singular        the CURRENT model's information block could not
#                            be formed or inverted, so no candidate can be
#                            adjusted for it. It says nothing about the
#                            candidate, and it hits every candidate at the
#                            step.
#   information_nonpositive  the candidate's own observed information at
#                            beta = 0 is not positive.
#
# A refit of the extended model has its own Hessian, and tests the candidate.

fb_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)
}

fb_fit <- function(d) {
  suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                          dist = "weibull", theta = c(mu = 0.01, nu = 0.5),
                          fit = TRUE))
}

fb_scope <- ~ age + mal + com_iv

fb_screen <- function(fit, d, criterion, scope = fb_scope, ...) {
  suppressWarnings(hzr_stepwise(fit, scope = scope, data = d,
                                direction = "forward", criterion = criterion,
                                slentry = 0.30, trace = FALSE, ...))
}

test_that("a real fit whose score information cannot be formed is refitted (#570)", {
  skip_if_not_installed("numDeriv")
  # avc, a Weibull base. Once com_iv has entered, mu is so small that the
  # numeric Hessian the score test needs comes back non-finite: there is no
  # nuisance block, and every remaining candidate reads `nuisance_singular`.
  # The fit itself is ordinary.
  d <- fb_data()
  fit <- fb_fit(d)
  wald <- fb_screen(fit, d, "wald")
  # Known positive: the Wald criterion, which refits every candidate, enters
  # all three.
  expect_identical(wald$steps$variable, c("com_iv", "mal", "age"))

  step1 <- fb_screen(fit, d, "score", max_steps = 1L)
  expect_identical(step1$steps$variable, "com_iv")
  expect_identical(step1$steps$stat_type, "score_q")
  # A stepwise result is the fitted model.
  m1 <- step1
  expect_s3_class(m1, "hazard")
  expect_true(m1$fit$converged)
  skip_if(isTRUE(.hzr_score_nuisance(m1)$ok),
          "this platform forms the score information for this fit")
  for (v in c("age", "mal")) {
    expect_identical(.hzr_score_q(m1, v, phase = NULL, data = d)$reason,
                     "nuisance_singular")
  }

  score <- fb_screen(fit, d, "score")
  # The screen stopped after com_iv, having tested nothing further. It now
  # reaches the Wald criterion's selection, with its p-values.
  expect_false(score$criteria$stopped_uncomputable)
  expect_length(score$criteria$uncomputable_reasons, 0L)
  expect_identical(score$steps$variable, wald$steps$variable)
  expect_identical(score$steps$stat_type, c("score_q", "wald_z", "wald_z"))
  expect_equal(log(score$steps$p_value[2:3]), log(wald$steps$p_value[2:3]),
               tolerance = 1e-6)
  # mal and age at step 2, age at step 3.
  expect_identical(score$criteria$n_wald_fallbacks, 3L)
})

test_that("a singular nuisance block sends every candidate to a refit (#570)", {
  # The same reason, forced at the first step on every platform. The refit
  # and the Wald test are real.
  d <- fb_data()
  fit <- fb_fit(d)
  wald <- fb_screen(fit, d, "wald", max_steps = 1L)
  expect_identical(wald$steps$variable, "com_iv")
  expect_gt(wald$steps$p_value, 0)

  orig <- .hzr_score_nuisance
  local_mocked_bindings(.hzr_score_nuisance = function(current) {
    r <- orig(current)
    list(inv = NULL, idx = r$idx, ok = FALSE)
  })
  # The premise: with no usable nuisance block, no candidate is scored.
  for (v in c("age", "mal", "com_iv")) {
    expect_identical(.hzr_score_q(fit, v, phase = NULL, data = d)$reason,
                     "nuisance_singular")
  }

  sw <- fb_screen(fit, d, "score", max_steps = 1L)
  # All three were refitted and tested; the screen did not stop untested.
  expect_identical(sw$criteria$n_wald_fallbacks, 3L)
  expect_false(sw$criteria$stopped_uncomputable)
  expect_length(sw$criteria$uncomputable_reasons, 0L)
  # The rescued test is the Wald criterion's own: same variable, same p.
  expect_identical(sw$steps$variable, "com_iv")
  expect_identical(sw$steps$stat_type, "wald_z")
  expect_equal(log(sw$steps$p_value), log(wald$steps$p_value),
               tolerance = 1e-6)
})

test_that("a candidate with non-positive information is refitted (#570)", {
  skip_if_not_installed("numDeriv")
  d <- fb_data()
  fit <- fb_fit(d)
  wald <- fb_screen(fit, d, "wald", max_steps = 1L)
  expect_identical(wald$steps$variable, "com_iv")

  # com_iv alone is given the reason; the others are scored as they are. A
  # screen that leaves com_iv untested enters the best of the rest.
  orig_q <- .hzr_score_q
  local_mocked_bindings(.hzr_score_q = function(current, var, ...) {
    r <- orig_q(current, var, ...)
    if (identical(var, "com_iv")) {
      r$stat <- NA_real_
      r$p_value <- NA_real_
      r$reason <- "information_nonpositive"
    }
    r
  })
  others <- vapply(c("age", "mal"), function(v) {
    .hzr_score_q(fit, v, phase = NULL, data = d)$p_value
  }, numeric(1))
  expect_lt(min(others), 0.30)
  expect_identical(.hzr_score_q(fit, "com_iv", phase = NULL, data = d)$reason,
                   "information_nonpositive")

  sw <- fb_screen(fit, d, "score", max_steps = 1L)
  expect_identical(sw$criteria$n_wald_fallbacks, 1L)
  expect_identical(sw$steps$variable, "com_iv")
  expect_equal(log(sw$steps$p_value), log(wald$steps$p_value),
               tolerance = 1e-6)
})

test_that("the information fault itself, at its source, is rescued (#570)", {
  skip_if_not_installed("numDeriv")
  d <- fb_data()
  fit <- fb_fit(d)
  # Each candidate's own diagonal of the observed information is flipped, as
  # test-score-uncomputable-reasons.R does, so .hzr_score_q() reaches the
  # reason by its own guard.
  orig <- .hzr_score_information_expanded
  local_mocked_bindings(
    .hzr_score_information_expanded = function(current, exp_) {
      info <- orig(current, exp_)
      if (!is.null(info)) info[exp_$beta_idx, exp_$beta_idx] <- -1e-3
      info
    }
  )
  expect_identical(.hzr_score_q(fit, "mal", phase = NULL, data = d)$reason,
                   "information_nonpositive")
  sw <- fb_screen(fit, d, "score", max_steps = 1L)
  expect_identical(sw$criteria$n_wald_fallbacks, 3L)
  expect_identical(sw$steps$variable, "com_iv")
})

test_that("reasons a refit cannot rescue still cost no refit (#570)", {
  d <- fb_data()
  d$flat <- 1
  fit <- fb_fit(d)
  n_refit <- 0L
  orig <- .hzr_refit_with_scope
  local_mocked_bindings(.hzr_refit_with_scope = function(...) {
    n_refit <<- n_refit + 1L
    orig(...)
  })
  sw <- fb_screen(fit, d, "score", scope = ~ flat)
  expect_identical(n_refit, 0L)
  expect_identical(sw$criteria$n_wald_fallbacks, 0L)
  expect_true(sw$criteria$stopped_uncomputable)
  expect_identical(sw$criteria$uncomputable_reasons, c(constant = 1L))

  # The set is the four reasons that describe the score test's approximation
  # or the current fit, and none that describes the candidate's column.
  expect_setequal(.hzr_score_fallback_reasons,
                  c("information_indefinite", "coefficient_diverging",
                    "nuisance_singular", "information_nonpositive"))
})

test_that("a rescue that fails under a new reason is reported as untested (#570)", {
  # The refit is what makes these candidates testable, so one whose refit
  # fails was tested by neither criterion. The warning that says so was
  # keyed on `information_indefinite` alone, and in hzr_bootstrap(), whose
  # replicates run with warnings suppressed, nothing else reports it.
  d <- fb_data()
  fit <- fb_fit(d)
  orig_n <- .hzr_score_nuisance
  orig_r <- .hzr_refit_with_scope
  local_mocked_bindings(
    .hzr_score_nuisance = function(current) {
      r <- orig_n(current)
      list(inv = NULL, idx = r$idx, ok = FALSE)
    },
    .hzr_refit_with_scope = function(current, action = c("add", "drop"),
                                     var, ...) {
      if (identical(var, "mal")) stop("forced refit failure")
      orig_r(current, action = action, var = var, ...)
    }
  )

  w <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fit, scope = fb_scope, data = d, direction = "forward",
                 criterion = "score", slentry = 0.30, max_steps = 1L,
                 trace = FALSE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  # The run completed: com_iv entered on its rescued test, mal went untested.
  expect_identical(sw$steps$variable, "com_iv")
  expect_false(sw$criteria$stopped_uncomputable)
  expect_identical(sw$criteria$uncomputable_reasons,
                   c(nuisance_singular = 1L))
  expect_true("mal" %in% sw$criteria$refit_failures)
  expect_true(any(grepl("Wald-fallback refit failed for mal", w,
                        fixed = TRUE)))
  neither <- grep("NEITHER criterion", w, fixed = TRUE, value = TRUE)
  expect_length(neither, 1L)
  expect_match(neither, "1 candidate(s)", fixed = TRUE)

  # The bootstrap: each replicate completes its one step with mal untested.
  bw <- character()
  bs <- withCallingHandlers(
    hzr_bootstrap(fit, n_boot = 3, seed = 570, scope = fb_scope,
                  direction = "forward", slentry = 0.30, max_steps = 1L),
    warning = function(x) {
      bw <<- c(bw, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(bs$n_success, 3L)
  expect_identical(bs$n_uncomputable_replicates, 0L)
  expect_identical(bs$uncomputable_reasons[["nuisance_singular"]], 3L)
  # The up-front screen on the data warns for itself; this is the replicates'.
  b_neither <- grep("replicates were tested by NEITHER criterion", bw,
                    fixed = TRUE, value = TRUE)
  expect_length(b_neither, 1L)
  expect_match(b_neither, "3 candidate score(s) across 3 replicates",
               fixed = TRUE)
})

test_that("every refitted reason reports its failed rescue (#570)", {
  # The warning reads the fallback set, so each reason in it is covered, not
  # only the one the test above forces. `mal` is given each reason in turn,
  # and its refit fails.
  d <- fb_data()
  fit <- fb_fit(d)
  expect_gt(length(.hzr_score_fallback_reasons), 2L)
  for (reason in .hzr_score_fallback_reasons) {
    orig_q <- .hzr_score_q
    orig_r <- .hzr_refit_with_scope
    local({
      local_mocked_bindings(
        .hzr_score_q = function(current, var, ...) {
          r <- orig_q(current, var, ...)
          if (identical(var, "mal")) {
            r$stat <- NA_real_
            r$p_value <- NA_real_
            r$reason <- reason
          }
          r
        },
        .hzr_refit_with_scope = function(current, action = c("add", "drop"),
                                         var, ...) {
          if (identical(var, "mal")) stop("forced refit failure")
          orig_r(current, action = action, var = var, ...)
        }
      )
      w <- character()
      sw <- withCallingHandlers(
        hzr_stepwise(fit, scope = fb_scope, data = d, direction = "forward",
                     criterion = "score", slentry = 0.30, max_steps = 1L,
                     trace = FALSE),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        }
      )
      expect_identical(sw$criteria$uncomputable_reasons,
                       stats::setNames(1L, reason), label = reason)
      expect_length(grep("NEITHER criterion", w, fixed = TRUE), 1L)
    })
  }
})

test_that("the bootstrap counts Wald ENTRIES, not candidates Wald-tested (#570)", {
  # With no usable nuisance block every candidate is refitted and tested, so
  # the screen's n_wald_fallbacks is positive even when none enters. The
  # bootstrap's count and its warning are about entries.
  d <- fb_data()
  fit <- fb_fit(d)
  orig_n <- .hzr_score_nuisance
  local_mocked_bindings(.hzr_score_nuisance = function(current) {
    r <- orig_n(current)
    list(inv = NULL, idx = r$idx, ok = FALSE)
  })
  boot <- function(slentry) {
    w <- character()
    bs <- withCallingHandlers(
      hzr_bootstrap(fit, n_boot = 3, seed = 570, scope = fb_scope,
                    direction = "forward", slentry = slentry, max_steps = 1L),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      }
    )
    list(bs = bs, entered = any(grepl("entered at least one variable on a Wald",
                                      w, fixed = TRUE)))
  }
  # Known positive: at an ordinary slentry each replicate enters a variable
  # on its rescued Wald test.
  yes <- boot(0.30)
  expect_identical(yes$bs$n_wald_fallback_replicates, 3L)
  expect_identical(yes$bs$n_wald_fallbacks, 3L)
  expect_true(yes$entered)
  # An slentry nothing can meet: every candidate is still Wald-tested, and
  # none enters.
  no <- boot(1e-300)
  expect_identical(no$bs$n_success, 3L)
  expect_identical(no$bs$n_wald_fallback_replicates, 0L)
  expect_identical(no$bs$n_wald_fallbacks, 0L)
  expect_false(no$entered)

  # Under criterion = "wald" every entry is a Wald z by design: no fallback.
  wald <- suppressWarnings(hzr_bootstrap(
    fit, n_boot = 3, seed = 570, scope = fb_scope, direction = "forward",
    criterion = "wald", slentry = 0.30, max_steps = 1L
  ))
  expect_identical(wald$n_success, 3L)
  expect_gt(sum(wald$summary$n[wald$summary$parameter %in%
                                 c("age", "mal", "com_iv")]), 0)
  expect_identical(wald$n_wald_fallback_replicates, 0L)
})

test_that("the two reasons' texts say the refit was tried (#570)", {
  for (r in c("nuisance_singular", "information_nonpositive")) {
    txt <- .hzr_score_reason_text(r)
    expect_match(txt, "refit", fixed = TRUE)
    expect_match(txt, "refit_failures", fixed = TRUE)
  }
})
