# In a two-way screen a variable that reached max_move on ENTRY was frozen at
# once, but the iteration's protected sets had been fixed before the forward
# step, so the backward step that followed could still drop it: the trace
# read FROZEN then DROP, and $scope$frozen named a variable the final model
# did not contain (#580). A freeze now takes effect after the backward step,
# so a variable that enters and is dropped in the same iteration is frozen
# OUT. Outside NOSTEPWISE, PROC HAZARD counts only exits
# (src/hazard/hazrd4.c:361-377), so it freezes a variable only as it leaves.

fz_data <- function(env = parent.frame()) {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  withr::local_seed(580, .local_envir = env)
  d$noise <- stats::rnorm(nrow(d))
  d
}

fz_fit <- function(d) {
  suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                          dist = "weibull", theta = c(mu = 0.01, nu = 0.5),
                          fit = TRUE))
}

fz_screen <- function(fit, d, max_move, direction = "both", scope = ~ noise) {
  suppressWarnings(hzr_stepwise(fit, scope = scope, data = d,
                                direction = direction, criterion = "score",
                                slentry = 0.9, slstay = 0.001,
                                max_move = max_move, max_steps = 30L,
                                trace = FALSE))
}

fz_has <- function(sw, var) any(grepl(var, names(sw$fit$theta), fixed = TRUE))

test_that("a variable frozen in a two-way screen is frozen where it ends (#580)", {
  skip_on_cran() # eight short two-way screens
  d <- fz_data()
  fit <- fz_fit(d)
  # The premise: noise enters on its score (p about 0.55 < slentry 0.9) and
  # fails slstay 0.001 once in, so it moves on every iteration.
  q <- .hzr_score_q(fit, "noise", phase = NULL, data = d)
  expect_gt(q$p_value, 0.001)
  expect_lt(q$p_value, 0.9)

  expected <- list(
    `1` = c("enter", "drop", "frozen"),
    `2` = c("enter", "drop", "enter", "drop", "frozen"),
    `3` = c("enter", "drop", "enter", "drop", "frozen"),
    `4` = c("enter", "drop", "enter", "drop", "enter", "drop", "frozen")
  )
  for (mm in names(expected)) {
    sw <- fz_screen(fit, d, as.integer(mm))
    # Freezing on an entry (max_move 2 and 4) read enter, FROZEN, drop.
    expect_identical(sw$steps$action, expected[[mm]], label = paste("max_move", mm))
    # Frozen out: the frozen set agrees with the final model.
    expect_identical(sw$scope$frozen, "noise")
    expect_false(fz_has(sw, "noise"))
    expect_false(sw$criteria$hit_max_steps)
  }
})

test_that("a freeze still holds a variable that stays in (#580)", {
  skip_on_cran() # three short screens
  d <- fz_data()
  fit <- fz_fit(d)
  # Known positive for the detector, and the forward-only case the issue
  # does not touch: noise enters and nothing can drop it.
  fwd <- fz_screen(fit, d, 4L, direction = "forward")
  expect_identical(fwd$steps$action, "enter")
  expect_true(fz_has(fwd, "noise"))

  # A strong variable that no step drops, frozen on its entry at
  # max_move = 0: frozen in, and still in at the end.
  sw <- suppressWarnings(hzr_stepwise(fit, scope = ~ com_iv, data = d,
                                      direction = "both", criterion = "score",
                                      slentry = 0.3, slstay = 0.2,
                                      max_move = 0L, max_steps = 30L,
                                      trace = FALSE))
  expect_identical(sw$steps$action, c("enter", "frozen"))
  expect_identical(sw$scope$frozen, "com_iv")
  expect_true(fz_has(sw, "com_iv"))
})
