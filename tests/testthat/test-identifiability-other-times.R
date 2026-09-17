# The saturated-phase warning measures a phase's variation across the event
# times, and said the likelihood is unchanged whether its shapes are pinned
# or fitted. An interval-censored or left-truncated likelihood also evaluates
# the phase at the interval bounds and the entry times, where a phase flat
# across the event times can be climbing, and there the claim is false: on the
# fit below by 103 log-likelihood units (#228). The verdict is unchanged --
# what those further points are worth is not measured here -- but the warning
# no longer claims more than was measured.

other_times_setup <- function(t_half = 1e-6, nu = 0, m = -0.4) {
  phases <- .hzr_validate_phases(
    list(early = hzr_phase("cdf", t_half = t_half, nu = nu, m = m),
         constant = hzr_phase("constant"))
  )
  list(phases = phases,
       counts = c(early = 0L, constant = 0L),
       x_list = list(early = NULL, constant = NULL),
       theta = c(log(0.045), log(t_half), nu, m, log(0.036)))
}

other_times_data <- function() {
  data("avc", package = "TemporalHazard", envir = environment())
  t <- avc$int_dead
  t[t > 0]
}

test_that("the saturated warning does not claim the likelihood is unchanged when other points exist (#228)", {
  s <- other_times_setup()
  time <- other_times_data()
  lower <- rep(c(2e-7, 1e-6, 3e-6), length.out = length(time))
  loglik <- function(t_half) {
    theta <- c(log(0.045), log(t_half), 0, -0.4, log(0.036))
    .hzr_logl_multiphase(theta, time, rep(2L, length(time)),
                         time_lower = lower, time_upper = time,
                         phases = s$phases, covariate_counts = s$counts,
                         x_list = s$x_list)
  }
  # The shape the old message called unchanged moves the likelihood far.
  expect_gt(abs(loglik(1e-6) - loglik(1e-7)), 10)

  w <- capture_warnings(
    .hzr_check_phase_identifiability(s$theta, time, s$phases, s$counts,
                                     s$x_list,
                                     other_times = c(lower, time))
  )
  expect_length(w, 1L)
  expect_match(w, "constant across the observed times")
  expect_match(w, "'mu' remains identified but the shape parameters do not")
  expect_match(w, "at least not by those times")
  expect_match(w, "further times \\(counting-process entry times or interval")
  expect_false(grepl("likelihood is unchanged", w))
})

test_that("with nothing else evaluated, the warning is unchanged (#228)", {
  s <- other_times_setup()
  time <- other_times_data()
  w <- capture_warnings(
    .hzr_check_phase_identifiability(s$theta, time, s$phases, s$counts,
                                     s$x_list)
  )
  expect_length(w, 1L)
  expect_match(w, "the likelihood is unchanged whether they are pinned or")
  expect_false(grepl("further time", w))

  # A bound of 0, or one equal to an event time, adds no evaluation point:
  # Lambda(0) is 0 for every phase, and the event times are already measured.
  w0 <- capture_warnings(
    .hzr_check_phase_identifiability(s$theta, time, s$phases, s$counts,
                                     s$x_list,
                                     other_times = c(rep(0, length(time)),
                                                     time))
  )
  expect_match(w0, "the likelihood is unchanged whether they are pinned or")
})

test_that("the verdict and the reported shares are unchanged (#228)", {
  # Only the wording moves: which phases are reported, and fit$phase_share,
  # stay as they were.
  s <- other_times_setup()
  time <- other_times_data()
  sh <- suppressWarnings(
    .hzr_check_phase_identifiability(s$theta, time, s$phases, s$counts,
                                     s$x_list,
                                     other_times = c(rep(2e-7, length(time)),
                                                     time))
  )
  expect_identical(colnames(sh), c("phase", "share", "variation"))
  expect_identical(sh$variation[sh$phase == "early"], 0)
  expect_gt(sh$share[sh$phase == "early"], 0.5)
  expect_null(attr(sh, "time_variation"))

  # A healthy fit still raises nothing, with or without other times.
  healthy <- other_times_setup(t_half = 0.003)
  expect_no_warning(
    .hzr_check_phase_identifiability(healthy$theta, time, healthy$phases,
                                     healthy$counts, healthy$x_list,
                                     other_times = c(rep(2e-7, length(time)),
                                                     time))
  )
})
