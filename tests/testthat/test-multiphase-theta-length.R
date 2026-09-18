# A supplied multiphase `theta` must have one entry per parameter: each
# phase's log_mu and free shapes, then one coefficient per column of the
# design the phase uses (its own `formula`, or the global one it inherits).
# The only check compared the length with the global design's column count,
# and only as a lower bound, so a longer `theta` fitted with the extra
# entries carried along, and a shorter one failed inside the fit with
# "'names' attribute [13] must be the same length as the vector [11]" (#408).

theta_len_data <- function() {
  withr::local_seed(11)
  n <- 200
  d <- data.frame(f1 = factor(sample(c("a", "b", "c"), n, TRUE)),
                  age = stats::rnorm(n))
  d$time <- stats::rexp(n, 0.3)
  d$dead <- stats::rbinom(n, 1, 0.7)
  d
}

theta_len_phases <- function(own_formula = NULL) {
  list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "m",
                         formula = own_formula),
       constant = hzr_phase("constant"))
}

theta_len_fit <- function(theta, phases = theta_len_phases(), fit = TRUE) {
  hazard(survival::Surv(time, dead) ~ f1, data = theta_len_data(),
         dist = "multiphase", phases = phases, theta = theta, fit = fit,
         control = list(n_starts = 1L, conserve = FALSE))
}

# Both phases inherit `~ f1` (columns f1b, f1c): early takes log_mu,
# log_t_half, nu, m and two coefficients (6), constant takes log_mu and two
# coefficients (3), so the model takes 9.
theta_ok <- c(0, log(0.5), 1, 1, 0, 0, log(0.3), 0, 0)

test_that("a multiphase theta of the right length fits (#408)", {
  fit <- suppressWarnings(theta_len_fit(theta_ok))
  expect_length(fit$fit$theta, 9L)
  expect_true(is.finite(fit$fit$objective))
})

test_that("a longer multiphase theta stops, naming both lengths (#408)", {
  expect_error(theta_len_fit(c(theta_ok, 0, 0)),
               "'theta' has 11 entries, but this model takes 9")
})

test_that("a shorter multiphase theta stops with the same message (#408)", {
  # Was "'names' attribute [9] must be the same length as the vector [7]",
  # raised from inside the fit.
  expect_error(theta_len_fit(theta_ok[1:7]),
               "'theta' has 7 entries, but this model takes 9")
})

test_that("a phase with its own formula is counted from that formula (#408)", {
  # early: log_mu, log_t_half, nu, m, age (5); constant inherits ~ f1 (3).
  ph <- theta_len_phases(~ age)
  own_ok <- c(0, log(0.5), 1, 1, 0, log(0.3), 0, 0)
  fit <- suppressWarnings(theta_len_fit(own_ok, phases = ph))
  expect_length(fit$fit$theta, 8L)
  # The old design's length for this model (the inherited count, 9) is one
  # too many here, so counting every phase from the global design would
  # pass this and fail the fit above.
  expect_error(theta_len_fit(c(own_ok, 0), phases = ph),
               "'theta' has 9 entries, but this model takes 8")
})

test_that("an unfitted multiphase model still carries theta as supplied (#408)", {
  # fit = FALSE returns theta itself, so the caller sees what was passed;
  # the stop is for a fit, where extra entries would corrupt the result.
  obj <- theta_len_fit(c(theta_ok, 0, 0), fit = FALSE)
  expect_identical(obj$fit$theta, c(theta_ok, 0, 0))
})

test_that("no theta still means the fit builds its own starts (#408)", {
  fit <- suppressWarnings(theta_len_fit(NULL))
  expect_length(fit$fit$theta, 9L)
})
