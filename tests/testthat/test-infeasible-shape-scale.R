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
              theta = function(ls) c(log(0.1), ls, 3, 1, 1, log(0.05))),
    # A separate selector in .hzr_phase_scale_feasible(): it shares cdf's
    # log_t_half, so mapping it to log_tau, or dropping it, must fail here.
    hazard = list(phases = list(h = hzr_phase("hazard", t_half = 3, nu = 1,
                                              m = 0),
                                c = hzr_phase("constant")),
                  theta = function(ls) c(log(0.1), ls, 1, 0, log(0.05)))
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
    # A feasible scale is untouched: finite likelihood, finite score, and a
    # Conservation of Events solve that actually runs. Each is asserted so
    # that a regression which quietly treated a FEASIBLE scale as infeasible
    # would fail -- the likelihood alone would not catch a zero-filled score
    # or a skipped solve.
    th <- layouts[[nm]]$theta(log(3))
    expect_true(is.finite(.hzr_logl_multiphase(
      th, d$time, d$status, phases = ph, covariate_counts = counts,
      x_list = x_list)), info = nm)
    # Unsanitised, so an unevaluable component shows as NA rather than 0;
    # and not all zero, since a feasible non-optimal point has a gradient.
    g <- .hzr_gradient_multiphase(th, d$time, d$status, phases = ph,
                                  covariate_counts = counts, x_list = x_list,
                                  sanitize = FALSE)
    expect_true(all(is.finite(g)), info = nm)
    expect_false(isTRUE(all(g == 0)), info = nm)
    # The solve moves the conserved phase's log_mu; the infeasible branch
    # above returns theta unchanged, so equality here would mean it was
    # skipped.
    solved <- .hzr_conserve_events(th, "c", conserved, d$time, d$status, ph,
                                   counts, x_list, sum(d$status))
    expect_true(all(is.finite(solved)), info = nm)
    expect_false(identical(solved, th), info = nm)
  }
})

test_that("the single-start fit from #262 completes instead of erroring (#262)", {
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
  # What #262 claims, and what main fails: the single start is not discarded
  # as errored. On main this fit stops with "no usable fit from 1 start: 1
  # errored"; here the start completes. Its status is "ok" or
  # "nonconverged": the loop writes "ok" only for optim()'s convergence code
  # 0, which is the optimizer's flag, not #262's claim. "error" (main) and
  # "infeasible" (a finite start the guard wrongly penalised) both fail.
  expect_length(fit$fit$starts$status, 1L)
  expect_true(fit$fit$starts$status %in% c("ok", "nonconverged"))
  expect_true(is.finite(fit$fit$objective))
  expect_gt(fit$fit$objective, -1e9)
  # Deliberately NOT asserted: convergence. This fit reports converged = TRUE,
  # but its relative gradient (about 0.006) is roughly 1000 times over the
  # limit SAS/C HAZARD accepts, and its Hessian is not positive definite, so
  # it is not a verified optimum -- on a three-phase model these data
  # probably cannot identify. #262 claims the start stops ERRORING, not that
  # the fit reaches an optimum. Pinning `converged` would assert the
  # optimizer's own flag, which says nothing about #262 here (see #351:
  # `converged` is TRUE when a fit fails the relative-gradient test).
})

test_that("a corrupt interval bound still stops the fit whatever the scale (#262 x #394)", {
  # #394 guards the DATA (interval bounds); #262 guards the PARAMETERS (a
  # time scale stepped out of range). Disjoint inputs -- but #394's comment
  # says order is load-bearing: a corrupt row must STOP the fit, not return
  # -Inf for the optimizer to walk away from, and #262's guard returns -Inf.
  # Inside .hzr_logl_multiphase() #262's return does come first, so a DIRECT
  # call with both defects returns -Inf. That internal order is not pinned,
  # since checking the bounds first would be an improvement. The user-facing
  # guarantee is pinned instead: hazard() checks the data at entry
  # (.hzr_check_sas_data) before any likelihood is evaluated. An NA bound is
  # refused earlier still, by input validation (#232's test), so the corrupt
  # row that reaches this check through hazard() is a finite interval whose
  # upper bound does not exceed its lower one.
  tt <- c(1, 2, 3, 4, 5, 6)
  st <- c(1, 0, 2, 1, 2, 1)
  lo <- c(0, 0, 3, 0, 2, 0)   # row 3: lower == upper, a zero-width interval
  ph <- list(e = hzr_phase("cdf", t_half = 3, nu = 1, m = 0),
             c = hzr_phase("constant"))
  fit_at <- function(log_t_half) {
    hazard(time = tt, status = st, time_lower = lo, time_upper = tt,
           dist = "multiphase", phases = ph,
           theta = c(log(0.1), log_t_half, 1, 0, log(0.05)),
           objective = "sas", fit = TRUE)
  }
  msg <- "requires upper > lower on every interval-censored row"
  # Out-of-range scale AND a corrupt interval: the data defect wins.
  expect_error(fit_at(-800), msg)
  # Control: a feasible scale stops identically, so the stop is the data's,
  # not a side effect of the scale.
  expect_error(fit_at(log(3)), msg)
})
