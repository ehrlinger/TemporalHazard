# A drop that removes no column is refused (#320).
#
# Under treatment contrasts, `model.matrix()` codes an interaction whose main
# effect is absent with a full set of dummies: `~ z:f` gives `z:fa, z:fb`,
# which spans what `z, z:fb` spans. So dropping the term `z` from `~ z + z:f`
# removes no column, and the "reduced" model is the model it started from:
# the same column space, the same coefficient count, the same likelihood.
#
# The step used to accept that drop and report a p-value for it. An uncapped
# run then stopped at the next step, because `z:f` had become two columns --
# loud, but by accident. A run that ended right after the drop (`max_steps`)
# returned it as a result, with `$steps` and the final model disagreeing.
#
# The forward step has the mirror of this guard: a candidate must ADD a
# column (`.hzr_entered_coef_name()`, #306). These tests pin the drop side.

nod_data <- function(n = 400L, seed = 5L) {
  set.seed(seed)
  d <- data.frame(z = stats::rnorm(n),
                  f = factor(sample(c("a", "b"), n, replace = TRUE)))
  d$time <- stats::rexp(n) * exp(-0.5 * d$z * (d$f == "b"))
  d$status <- 1L
  d
}

nod_multi <- function(d, formula = ~ z + z:f) {
  hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
         phases = list(constant = hzr_phase("constant", formula = formula)),
         fit = TRUE)
}

test_that("a multiphase drop that removes no column is refused, not accepted", {
  d <- nod_data()
  fit <- nod_multi(d)
  # The fixture's mechanism: dropping `z` reparameterises rather than reduces.
  expect_identical(colnames(fit$fit$x_list$constant), c("z", "z:fb"))

  expect_warning(
    step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                        slstay = 0.2),
    "removes no column"
  )
  # `z` was the chosen drop: the guard fires on the decision, not before it.
  row <- step$all_scores[step$all_scores$variable == "z", ]
  expect_equal(nrow(row), 1L)
  expect_gt(row$p_value, 0.2)

  expect_false(step$accepted)
  expect_true(is.na(step$variable))
  # The returned fit still carries `z`. Under the bug it is the refitted,
  # reparameterised model, whose coefficients are `constant.z:fa`/`:fb`, so
  # these fail there. (`null_result()` returns `fit = current`, so comparing
  # that object with `fit` would compare it with itself and could not fail.)
  # They come before the reasons lookup below, which errors under the bug --
  # an error ends the block, and assertions after it would never run.
  expect_true("constant.z" %in% names(stats::coef(step$fit)))
  expect_identical(.hzr_scope_current_vars(step$fit, "constant"),
                   c("z", "z:f"))
  expect_identical(step$refit_failures, "z@constant")
  expect_match(step$refit_failure_reasons[["z@constant"]], "removes no column")
})

test_that("a capped run records no drop and returns the unchanged model", {
  # This is the silent case: the run ends at the cap right after the bad drop,
  # so nothing downstream errors and `$steps` is all the reader sees.
  d <- nod_data()
  fit <- nod_multi(d)
  sw <- suppressWarnings(hzr_stepwise(fit, data = d, direction = "backward",
                                      criterion = "wald", slstay = 0.2,
                                      max_steps = 1L, trace = FALSE))
  expect_equal(sum(sw$steps$action == "drop"), 0L)
  expect_identical(colnames(sw$fit$x_list$constant), c("z", "z:fb"))
  # The phase formula still holds both terms. Equal objective and equal
  # coefficient count are the SIGNATURE of this bug, so asserting them would
  # hold whether or not the guard fired; the surviving term does not.
  expect_identical(.hzr_scope_current_vars(sw, "constant"), c("z", "z:f"))
})

test_that("an uncapped run stops on the guard, not on the next step's term", {
  # Before the guard, the run reached a second step where `z:f` had become
  # `z:fa, z:fb` and errored with "expands to multiple coefficients".
  d <- nod_data()
  fit <- nod_multi(d)
  w <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fit, data = d, direction = "backward", criterion = "wald",
                 slstay = 0.2, trace = FALSE),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("removes no column", w)))
  expect_equal(sum(sw$steps$action == "drop"), 0L)
  expect_identical(colnames(sw$fit$x_list$constant), c("z", "z:fb"))
})

test_that("a drop that does remove a column still happens", {
  # The guard must not refuse ordinary work: `w` owns its own column.
  d <- nod_data()
  d$w <- stats::rnorm(nrow(d))
  fit <- nod_multi(d, formula = ~ z + w)
  expect_identical(colnames(fit$fit$x_list$constant), c("z", "w"))
  step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                      slstay = 0.2)
  expect_true(step$accepted)
  expect_identical(step$variable, "w")
  expect_identical(colnames(step$fit$fit$x_list$constant), "z")
  expect_length(step$refit_failures, 0L)
})

test_that("the single-distribution path refuses the same drop for the same reason (#323)", {
  # The single-distribution refit warm-starts from `theta_old[-drop_idx]`, one
  # element shorter than a no-op drop's design, so it used to fail to conform
  # first and report "non-conformable arguments": a linear-algebra symptom,
  # not the cause. The reduced design is now decided before the refit.
  d <- nod_data()
  fit <- hazard(survival::Surv(time, status) ~ z + z:f, data = d,
                dist = "weibull", theta = c(0.5, 1, 0, 0), fit = TRUE)
  expect_identical(colnames(fit$data$x), c("z", "z:fb"))
  expect_warning(
    step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                        slstay = 0.2),
    "Stepwise backward: dropping z removes no column", fixed = TRUE
  )
  expect_false(step$accepted)
  expect_identical(step$refit_failures, "z")

  mp_step <- suppressWarnings(
    .hzr_stepwise_backward_step(nod_multi(d), data = d, criterion = "wald",
                                slstay = 0.2)
  )
  expect_identical(unname(step$refit_failure_reasons[["z"]]),
                   unname(mp_step$refit_failure_reasons[["z@constant"]]))
})

test_that("a single-distribution drop that does remove a column still happens", {
  d <- nod_data()
  d$w <- stats::rnorm(nrow(d))
  fit <- hazard(survival::Surv(time, status) ~ z + w, data = d,
                dist = "weibull", theta = c(0.5, 1, 0, 0), fit = TRUE)
  step <- .hzr_stepwise_backward_step(fit, data = d, criterion = "wald",
                                      slstay = 0.2)
  expect_true(step$accepted)
  expect_identical(step$variable, "w")
  expect_identical(colnames(step$fit$data$x), "z")
  expect_length(step$refit_failures, 0L)
})
