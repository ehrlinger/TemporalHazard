# A phase's time scale is carried on the log scale (log_t_half, log_tau), and
# an optimizer step can take it past exp()'s range: t_half or tau becomes Inf
# or 0. That is an infeasible point, which the likelihood must penalise so
# the optimizer backs away. It was an error instead, on the log-likelihood,
# the score and the Conservation of Events solve, so a single-start fit
# stopped with "t_half must be a positive scalar" (#262).

scale_data <- function() {
  withr::local_seed(3)
  n <- 300
  tt <- c(stats::rexp(n / 2, 3), stats::rexp(n / 2, 0.15))
  cens <- stats::runif(n, 1, 15)
  list(time = pmin(tt, cens), status = as.integer(tt <= cens))
}

test_that("an out-of-range time scale is infeasible, not an error (#262)", {
  d <- scale_data()
  layouts <- list(
    cdf = list(phases = list(e = hzr_phase("cdf", t_half = 3, nu = 1, m = 0),
                             c = hzr_phase("constant")),
               theta = function(ls) c(log(0.1), ls, 1, 0, log(0.05))),
    g3 = list(phases = list(g = hzr_phase("g3"), c = hzr_phase("constant")),
              theta = function(ls) c(log(0.1), ls, 3, 1, 1, log(0.05)))
  )
  for (nm in names(layouts)) {
    ph <- layouts[[nm]]$phases
    counts <- stats::setNames(c(0L, 0L), names(ph))
    x_list <- stats::setNames(list(NULL, NULL), names(ph))
    conserved <- length(layouts[[nm]]$theta(0))
    # Which cells were what on the previous version: cdf/-800 (t_half = 0)
    # raised "t_half must be a positive scalar" on all three paths, the
    # defect #262 was filed for. cdf/+800 and g3/-800 already returned -Inf
    # from the objective; the score raised at cdf/+800 and returned zeros at
    # g3/-800. g3/+800 returned a finite objective (the late phase switched
    # off) with a score that raised, and is the deliberate new refusal.
    for (ls in c(800, -800)) {
      th <- layouts[[nm]]$theta(ls)
      info <- paste(nm, ls)
      expect_identical(
        .hzr_logl_multiphase(th, d$time, d$status, phases = ph,
                             covariate_counts = counts, x_list = x_list),
        -Inf, info = info
      )
      expect_identical(
        .hzr_gradient_multiphase(th, d$time, d$status, phases = ph,
                                 covariate_counts = counts, x_list = x_list),
        rep(0, length(th)), info = info
      )
      expect_true(all(is.na(
        .hzr_gradient_multiphase(th, d$time, d$status, phases = ph,
                                 covariate_counts = counts, x_list = x_list,
                                 sanitize = FALSE)
      )), info = info)
      expect_identical(
        .hzr_conserve_events(th, "c", conserved, d$time, d$status, ph,
                             counts, x_list, sum(d$status)),
        th, info = info
      )
    }
    # A feasible scale is untouched: finite likelihood, finite score.
    th <- layouts[[nm]]$theta(log(3))
    expect_true(is.finite(.hzr_logl_multiphase(
      th, d$time, d$status, phases = ph, covariate_counts = counts,
      x_list = x_list)))
  }
})

test_that("the single-start fit from #262 fits (#262)", {
  skip_on_cran()
  withr::local_seed(3)
  n <- 300
  tt <- c(stats::rexp(n / 2, 3), stats::rexp(n / 2, 0.15))
  cens <- stats::runif(n, 1, 15)
  d <- data.frame(time = pmin(tt, cens), status = as.integer(tt <= cens),
                  z = stats::rnorm(n))
  fit <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
    phases = list(
      e1 = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 0, formula = ~ z,
                     fixed = "nu"),
      e2 = hzr_phase("cdf", t_half = 3, nu = 1, m = 0, fixed = "m"),
      late = hzr_phase("g3")
    ),
    fit = TRUE, control = list(n_starts = 1)
  ))
  expect_true(is.finite(fit$fit$objective))
  expect_gt(fit$fit$objective, -1e9)
})
