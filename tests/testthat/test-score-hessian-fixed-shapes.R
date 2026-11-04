# The score test's numeric Hessian differentiates only what it reads (#712).
#
# With a left- or interval-censored row the analytic Hessian declines and the
# score path falls back to a numeric one. That one differentiated every
# parameter, fixed shapes included. A PROC HAZARD early phase fixes nu = 0,
# and one step below it the negative log-likelihood is infinite, so the matrix
# was non-finite, the nuisance block was discarded, every candidate scored NA
# as `nuisance_singular`, and the screen ran on Wald refits: a different
# screen from SAS's, which scores the same candidates from the same block.

valves_fit <- function(n_interval) {
  data(valves, package = "TemporalHazard", envir = environment())
  d <- valves[!is.na(valves$nyha), ]
  st <- d$dead
  ev <- which(st == 1)[seq_len(n_interval)]
  st[ev] <- 2
  d$st <- st
  d$tl <- d$int_dead
  d$tl[ev] <- 0.8 * d$int_dead[ev]
  # The translated base model of the hazard repo's hm.deadp.VALVES job: the
  # early shapes fixed, nu on its boundary at 0.
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = d$st, time_lower = d$tl, time_upper = d$int_dead,
    data = d, fit = TRUE, dist = "multiphase",
    phases = list(hzr_phase("cdf", t_half = 0.6781774, nu = 0, m = -2.15596,
                            fixed = c("t_half", "nu", "m")),
                  hzr_phase("constant")),
    theta = c(log(0.09596854), log(0.6781774), 0, -2.15596, log(0.00349094))
  ))
  list(fit = fit, data = d)
}

test_that("a fixed shape on its boundary does not void the numeric Hessian", {
  m <- valves_fit(20L)
  # Known positive for the premise: the full numeric Hessian is non-finite,
  # so the old path had nothing to use.
  th <- m$fit$fit$theta
  expect_identical(unname(th[["phase_1.nu"]]), 0)
  full <- .hzr_score_multiphase_hessian(m$fit, th, .hzr_score_phases(m$fit),
                                         m$fit$fit$covariate_counts, m$fit$fit$x_list)
  expect_null(full)
  q <- .hzr_score_q(m$fit, "nyha", "phase_1", m$data)
  expect_true(is.finite(q$stat))
  expect_gt(q$stat, 50)
})

test_that("the restricted block equals the analytic Hessian's block", {
  # Right-censored only, so the analytic Hessian exists to compare against.
  m <- valves_fit(0L)
  th <- m$fit$fit$theta
  idx <- .hzr_score_free_idx(m$fit)
  an <- .hzr_hessian_multiphase(
    th, time = m$fit$data$time, status = m$fit$data$status,
    time_lower = m$fit$data$time_lower, time_upper = m$fit$data$time_upper,
    x = m$fit$data$x, weights = m$fit$data$weights,
    phases = .hzr_score_phases(m$fit),
    covariate_counts = m$fit$fit$covariate_counts, x_list = m$fit$fit$x_list
  )
  num <- .hzr_score_multiphase_hessian(m$fit, th, .hzr_score_phases(m$fit),
                                        m$fit$fit$covariate_counts, m$fit$fit$x_list,
                                        positions = idx)
  expect_equal(num[idx, idx] / an[idx, idx], matrix(1, length(idx), length(idx)),
               tolerance = 1e-5, ignore_attr = TRUE)
  # What it did not differentiate it does not claim to know.
  expect_true(all(is.na(num[-idx, ])))
  expect_true(all(is.na(num[, -idx])))
})

test_that("a score screen on interval-censored rows scores rather than refits", {
  skip_on_cran() # a stepwise screen
  m <- valves_fit(20L)
  sw <- suppressWarnings(hzr_stepwise(
    m$fit, scope = list(phase_1 = ~ nyha + i_path + mitral, phase_2 = ~ nyha + i_path),
    data = m$data, direction = "forward", criterion = "score",
    slentry = 0.2, max_steps = 1L, trace = FALSE
  ))
  # Known positive: the screen took its step.
  expect_identical(nrow(sw$steps), 1L)
  expect_identical(sw$criteria$n_wald_fallbacks, 0L)
})
