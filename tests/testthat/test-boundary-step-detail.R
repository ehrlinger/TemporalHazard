# Unit tests for .hzr_phase_step_detail() (#448).  `g_fn` is injected, so the
# synthetic cases pin the logic; the last two run the REAL decomposition, so a
# shared wrong assumption between test and implementation cannot hide.

test_that("a phase that rises inside one observation gap is reported", {
  tt <- c(1, 2, 3, 4)
  # numerically 0 up to t = 2, numerically 1 from t = 3: the whole rise sits
  # between two adjacent observed times.
  step <- function(x) ifelse(x < 2.5, 0, 1)
  out <- .hzr_phase_step_detail(tt, t_half = 2.5, nu = 1e-16, g_fn = step)
  expect_type(out, "list")
  expect_identical(out$parameter, "nu")
  expect_match(out$detail, "step at this data's resolution", fixed = TRUE)
  expect_match(out$detail, "1e-16", fixed = TRUE)   # the magnitude is carried
})

test_that("a smooth phase is NOT reported", {
  tt <- seq(0.5, 5, by = 0.5)
  smooth <- function(x) pnorm(x, mean = 2.5, sd = 1)
  expect_null(.hzr_phase_step_detail(tt, t_half = 2.5, nu = 1.4, g_fn = smooth))
})

test_that("degenerate inputs decline rather than guess", {
  step <- function(x) ifelse(x < 2.5, 0, 1)
  expect_null(.hzr_phase_step_detail(c(1), 2.5, 1e-16, step))        # one time
  expect_null(.hzr_phase_step_detail(c(1, 1, 1), 2.5, 1e-16, step))  # one distinct
  expect_null(.hzr_phase_step_detail(c(1, 2), NA_real_, 1e-16, step))
  expect_null(.hzr_phase_step_detail(c(1, 2), 2.5, NA_real_, step))
  expect_null(.hzr_phase_step_detail(c(1, 2), 2.5, 1e-16, g_fn = "not a fn"))
  # a g_fn that errors must decline, not propagate
  expect_null(.hzr_phase_step_detail(c(1, 2), 2.5, 1e-16,
                                     function(x) stop("boom")))
})

test_that("the REAL decomposition at #448's fitted nu is reported", {
  # #448's fit: nu = -1.43915e-16, t_half = 0.0657098, m = 1.
  th <- 0.0657098
  nu <- -1.43915e-16
  gfn <- function(x) hzr_decompos(x, t_half = th, nu = nu, m = 1)$G
  tt <- c(th * 0.5, th * 2, th * 4)
  out <- .hzr_phase_step_detail(tt, t_half = th, nu = nu, g_fn = gfn)
  expect_type(out, "list")
  expect_match(out$detail, "unidentified", fixed = TRUE)
})

test_that("the REAL decomposition at a healthy nu is NOT reported", {
  th <- 0.0657098
  gfn <- function(x) hzr_decompos(x, t_half = th, nu = 1.4, m = 1)$G
  tt <- c(th * 0.5, th * 2, th * 4)
  # the census measured |nu| >= 0.8175 on every healthy suite fit
  expect_null(.hzr_phase_step_detail(tt, t_half = th, nu = 1.4, g_fn = gfn))
})

test_that("interval bounds count as observed times", {
  step <- function(x) ifelse(x < 0.15, 0, 1)
  # `time` alone has only ONE positive value, so the old signature declined.
  # The interval bounds bracket the rise, so the phase IS a step at this
  # data's resolution and must be reported.
  expect_null(.hzr_phase_step_detail(c(5), t_half = 0.15, nu = 1e-16,
                                     g_fn = step))
  out <- .hzr_phase_step_detail(c(5), t_half = 0.15, nu = 1e-16, g_fn = step,
                                time_lower = c(0.1), time_upper = c(0.2))
  expect_type(out, "list")
  expect_match(out$detail, "0.1", fixed = TRUE)
})
