# Conservation of Events solves the conserved phase's log_mu so it absorbs
# the events the other phases leave. Where none are left there is nothing to
# solve, the objective falls back to that log_mu's starting value, and a fit
# can stop short of a higher likelihood with the conserved phase switched off
# (#261). The check certifies such a point with the fit's own objective.

coe_repro <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  hazard(survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
         phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
                       constant = hzr_phase("constant"),
                       late = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1,
                                        eta = 1)),
         fit = TRUE, control = list(n_starts = 1))
}

coe_records <- function(fit) {
  b <- fit$fit$boundary
  if (!is.list(b)) return(list())
  Filter(function(r) identical(r$mechanism, "coe_no_events_left"), b)
}

test_that("a CoE fit short of its boundary supremum is recorded and warned (#261)", {
  skip_on_cran()
  w <- list()
  fit <- withCallingHandlers(coe_repro(), warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd
    invokeRestart("muffleWarning")
  })
  # The premise, as #261 measured it: CoE applied, the constant phase
  # conserved, an ordinary-looking optimum at -206.7432.
  expect_true(fit$spec$control$conserve_applied)
  expect_equal(fit$fit$objective, -206.7432, tolerance = 1e-4)

  rec <- coe_records(fit)
  expect_length(rec, 1L)
  rec <- rec[[1L]]
  expect_identical(rec$phase, "constant")
  expect_identical(rec$parameter, "log_mu")
  # The certificate is the fit's own likelihood, recomputed from the point,
  # and it is where #261 put the supremum: above -206.7432, at -206.667.
  again <- as.numeric(hzr_evaluate(fit, rec$certificate_theta)$logLik)
  expect_equal(again, rec$certificate_loglik, tolerance = 1e-10)
  expect_gt(again, fit$fit$objective + 0.01)
  expect_lt(abs(again - (-206.669)), 1e-3)
  # The conserved phase is switched off there; only intercepts moved.
  expect_lt(rec$certificate_theta[["constant.log_mu"]], -600)
  shapes <- setdiff(names(coef(fit)), c("early.log_mu", "constant.log_mu",
                                        "late.log_mu"))
  expect_identical(rec$certificate_theta[shapes], coef(fit)[shapes])

  cls <- vapply(w, function(cnd) inherits(cnd, "hzr_coe_no_events_left"),
                logical(1))
  expect_identical(sum(cls), 1L)
  expect_true(inherits(w[[which(cls)]], "hzr_boundary"))
  expect_match(conditionMessage(w[[which(cls)]]), "'constant'", fixed = TRUE)
  expect_match(conditionMessage(w[[which(cls)]]), "starting value",
               fixed = TRUE)
  # Report only.
  expect_equal(as.numeric(hzr_evaluate(fit, coef(fit))$logLik),
               fit$fit$objective, tolerance = 1e-8)
})

test_that("a CoE fit at its maximum carries no boundary record (#261)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard", envir = environment())
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1)))
  expect_true(fit$spec$control$conserve_applied)
  # The check ran: it needs a converged fit and a finite objective.
  expect_true(fit$fit$converged)
  expect_true(is.finite(fit$fit$objective))
  expect_length(coe_records(fit), 0L)
})

test_that("only the certificate fires, past its margin, on a converged value (#261)", {
  skip_on_cran()
  fit <- suppressWarnings(coe_repro())
  theta <- unname(coef(fit))
  value <- fit$fit$objective
  phases <- .hzr_validate_phases(fit$spec$phases)
  counts <- c(early = 0L, constant = 0L, late = 0L)
  d <- fit$data
  pos <- .hzr_log_mu_positions(phases, counts)
  record <- function(objective_fn, converged = TRUE, v = value) {
    .hzr_coe_boundary_record(
      theta, v, converged, objective_fn, "constant", pos[["constant"]], pos,
      phases, counts, list(early = NULL, constant = NULL, late = NULL),
      d$time, d$status, NULL, NULL, sum(d$status == 1))
  }
  real <- function(p) as.numeric(hzr_evaluate(fit, p)$logLik)
  expect_false(is.null(record(real)))
  expect_null(record(function(p) value - 1))
  expect_null(record(function(p) value + 0.005))
  expect_false(is.null(record(function(p) value + 0.02)))
  expect_null(record(real, converged = FALSE))
  expect_null(record(real, v = NA_real_))
  # Both boundary points are scored: switched off with the others held, and
  # with the others rescaled to conserve the events again. Here the second
  # is the higher, so a check with only the first finds less.
  held_only <- function(p) {
    if (identical(p[pos[["early"]]], theta[pos[["early"]]])) real(p) else -Inf
  }
  rescaled_only <- function(p) {
    if (identical(p[pos[["early"]]], theta[pos[["early"]]])) -Inf else real(p)
  }
  expect_gt(record(rescaled_only)$gain, record(real)$gain - 1e-12)
  expect_lt(record(held_only)$gain %||% -Inf, record(real)$gain)
  # "Sent to zero" is literal: the conserved phase contributes exactly
  # nothing at the certificate, even at an extreme time.
  cert <- record(real)$certificate_theta
  xl0 <- list(early = NULL, constant = NULL, late = NULL)
  expect_true(all(.hzr_multiphase_cumhaz(
    c(d$time, 1e300), cert, phases, counts, xl0,
    per_phase = TRUE)[["constant"]] == 0))
  # And where a linear predictor would keep the phase alive at that floor,
  # nothing is certified rather than a live phase called off.
  counts1 <- c(early = 0L, constant = 1L, late = 0L)
  th1 <- append(theta, 1, after = pos[["constant"]])
  pos1 <- .hzr_log_mu_positions(phases, counts1)
  xl1 <- list(early = NULL, constant = matrix(2e4, length(d$time), 1L),
              late = NULL)
  expect_null(.hzr_coe_boundary_record(
    th1, value, TRUE, function(p) value + 1, "constant", pos1[["constant"]],
    pos1, phases, counts1, xl1, d$time, d$status, NULL, NULL,
    sum(d$status == 1)))

  # The held-only point is scored too: with the rescaled one refused, it
  # alone still certifies here.
  expect_false(is.null(record(held_only)))

  # The rescaled point conserves the events on the fit's own footing:
  # weighted, and net of entry-time cumulative hazard.
  w <- rep(c(1, 2), length.out = length(d$time))
  entry <- ifelse(seq_along(d$time) %% 5 == 0, d$time / 2, 0)
  total <- sum(w[d$status == 1])
  grab <- NULL
  .hzr_coe_boundary_record(
    theta, value, TRUE, function(p) {
      if (!identical(p[pos[["early"]]], theta[pos[["early"]]])) grab <<- p
      value - 1
    }, "constant", pos[["constant"]], pos, phases, counts,
    list(early = NULL, constant = NULL, late = NULL), d$time, d$status,
    entry, w, total)
  expect_false(is.null(grab))
  xl <- list(early = NULL, constant = NULL, late = NULL)
  net <- .hzr_multiphase_cumhaz(d$time, grab, phases, counts, xl) -
    .hzr_multiphase_cumhaz(entry, grab, phases, counts, xl)
  expect_equal(sum(w * net), total, tolerance = 1e-8)
})
