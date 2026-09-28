# A multiphase refit drops every row where a variable in the model is missing.
# So entering a variable with NAs shrinks the sample every later step is tested
# on, and dropping it grows the sample back. Nothing said so (#519). Each step
# now records the rows it was fitted on in `$steps$n_rows`, and a step that
# changes them warns with class `hzr_stepwise_sample_changed`.

.s519_fixture <- function(z_effect, v_effect = 0.8) {
  data(avc, package = "TemporalHazard", envir = environment())
  avc <- avc[!is.na(avc$int_dead) & !is.na(avc$dead), ]
  set.seed(7)
  na_rows <- sample(nrow(avc), 60)
  avc$z <- stats::rnorm(nrow(avc)) + z_effect * avc$dead
  avc$z[na_rows] <- NA
  avc$v <- stats::rnorm(nrow(avc)) + v_effect * avc$dead
  avc$u <- stats::rnorm(nrow(avc))
  avc
}

.s519_phases <- function(f = NULL) {
  list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                         fixed = "shapes"),
    constant = hzr_phase("constant", formula = f)
  )
}

.s519_fit <- function(data, f = NULL) {
  set.seed(1)
  hazard(Surv(int_dead, dead) ~ 1, data = data, dist = "multiphase",
         phases = .s519_phases(f), fit = TRUE,
         control = list(n_starts = 1L, maxit = 500L))
}

# Runs the screen and keeps every warning, split into the sample-change ones
# and the rest.
.s519_screen <- function(fit, data, ...) {
  changed <- list()
  other <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fit, data = data, slentry = 0.05, slstay = 0.05,
                 trace = FALSE,
                 control = list(n_starts = 1L, maxit = 500L), ...),
    warning = function(w) {
      if (inherits(w, "hzr_stepwise_sample_changed")) {
        changed[[length(changed) + 1L]] <<- w
      } else {
        other <<- c(other, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  list(sw = sw, changed = changed, other = other)
}

test_that("an entry that shrinks the sample warns with the counts (#519)", {
  # The issue's shape: a complete variable enters, then one missing on 60
  # rows enters and the model moves to 250 rows. On main this ran silently.
  d <- .s519_fixture(z_effect = 3, v_effect = 2)
  expect_equal(sum(is.na(d$z)), 60L)
  base <- .s519_fit(d)
  expect_equal(sum(.hzr_fit_row_mask(base)), 310L)

  out <- .s519_screen(base, d, scope = list(constant = ~ z + v),
                      direction = "forward", criterion = "wald")
  st <- out$sw$steps
  expect_equal(st$variable, c("v", "z"))
  expect_equal(st$n_rows, c(310L, 250L))
  expect_equal(sum(.hzr_fit_row_mask(out$sw)), 250L)
  # One change, one warning, naming the step, the variable and both counts.
  expect_length(out$changed, 1L)
  msg <- conditionMessage(out$changed[[1L]])
  expect_match(msg, "step 2 (enter z)", fixed = TRUE)
  expect_match(msg, "from 310 to 250", fixed = TRUE)
  expect_match(stepwise_trace(out$sw), "rows fitted changed: 310 -> 250",
               fixed = TRUE, all = FALSE)
})

test_that("a drop that grows the sample warns, and the next does not (#519)", {
  d <- .s519_fixture(z_effect = -0.3, v_effect = 0.5)
  base <- .s519_fit(d, ~ z + v)
  expect_equal(sum(.hzr_fit_row_mask(base)), 250L)

  out <- .s519_screen(base, d, direction = "backward", criterion = "wald")
  st <- out$sw$steps
  # z leaves first and the model moves back to 310 rows; v then leaves on
  # those same 310, which is no change and must not warn.
  expect_equal(st$variable, c("z", "v"))
  expect_equal(st$n_rows, c(310L, 310L))
  expect_length(out$changed, 1L)
  msg <- conditionMessage(out$changed[[1L]])
  expect_match(msg, "step 1 (drop z)", fixed = TRUE)
  expect_match(msg, "from 250 to 310", fixed = TRUE)
})

test_that("a complete-data screen records n_rows and does not warn (#519)", {
  d <- .s519_fixture(z_effect = 0)
  base <- .s519_fit(d)
  out <- .s519_screen(base, d, scope = list(constant = ~ v + u),
                      direction = "forward", criterion = "wald")
  st <- out$sw$steps
  # Known positive: a step was taken, so there was a step that could warn.
  expect_equal(st$variable, "v")
  expect_equal(st$n_rows, 310L)
  expect_length(out$changed, 0L)
})

test_that("an empty `$steps` still carries the n_rows column (#519)", {
  d <- .s519_fixture(z_effect = 0)
  base <- .s519_fit(d)
  out <- .s519_screen(base, d, scope = list(constant = ~ u),
                      direction = "forward", criterion = "wald")
  expect_equal(nrow(out$sw$steps), 0L)
  expect_identical(out$sw$steps$n_rows, integer())
})

test_that("a frozen row records the rows of the model it froze in (#519)", {
  # max_move = 0 freezes each variable on its first move, so each entry is
  # followed by a `frozen` row, written by record_freeze(), not record_step().
  # z's frozen row comes after the sample changed, so it must read 250, not
  # the base model's 310.
  d <- .s519_fixture(z_effect = 3, v_effect = 2)
  base <- .s519_fit(d)
  out <- .s519_screen(base, d, scope = list(constant = ~ z + v),
                      direction = "forward", criterion = "wald",
                      max_move = 0L)
  st <- out$sw$steps
  expect_equal(st$action, c("enter", "frozen", "enter", "frozen"))
  expect_equal(st$variable, c("v", "v", "z", "z"))
  expect_equal(st$n_rows, c(310L, 310L, 250L, 250L))
  expect_equal(st$n_rows[4L], sum(.hzr_fit_row_mask(out$sw)))
  # Freezing changes no row, so it adds no warning of its own.
  expect_length(out$changed, 1L)
})
