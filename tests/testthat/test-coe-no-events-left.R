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

coe_records <- function(fit, mechanism = "coe_no_events_left") {
  b <- fit$fit$boundary
  if (!is.list(b)) return(list())
  Filter(function(r) identical(r$mechanism, mechanism), b)
}

# Where the fit from coe_repro() stopped until #565: an ordinary-looking
# optimum at -206.7432, short of the boundary. The optimizer was then given
# the partial score under CoE; with the gradient of the objective the same
# call runs on towards the boundary, so no real fit stops here any more. The
# certificate is a function of a point, and this is the point #261 measured.
coe_stall_theta <- c(
  -1.38267369320279, -1.63637720507834, 1.82219933970156, 0.631796526964141,
  -10.6990572711889, -14.7979393759282, 1.16243431065041, 1.84633913620587,
  0.0639971149827021, 0.117921229256013
)

# The certificate at a point, scored by the fit's own likelihood.
coe_certificate <- function(fit, theta, value, objective_fn = NULL,
                            converged = TRUE) {
  phases <- .hzr_validate_phases(fit$spec$phases)
  counts <- c(early = 0L, constant = 0L, late = 0L)
  pos <- .hzr_log_mu_positions(phases, counts)
  d <- fit$data
  if (is.null(objective_fn)) {
    objective_fn <- function(p) as.numeric(hzr_evaluate(fit, p)$logLik)
  }
  .hzr_coe_boundary_record(
    theta, value, converged, objective_fn, "constant", pos[["constant"]], pos,
    phases, counts, list(early = NULL, constant = NULL, late = NULL),
    d$time, d$status, NULL, NULL, sum(d$status == 1))
}

test_that("a fit that ends short of the boundary is recorded and warned, end to end (#261)", {
  skip_on_cran()
  # No real fit stops at #261's point any more (#565), so the optimizer's
  # return is pinned there: the run is real, then its estimates are replaced
  # by the stall point. Everything after the optimizer -- the conservation
  # step, the certificate, the record on $boundary and hazard()'s warning --
  # runs as in any fit. The stall point has every parameter free except the
  # conserved constant.log_mu, the fifth.
  orig <- .hzr_optim_generic
  local_mocked_bindings(.hzr_optim_generic = function(...) {
    a <- list(...)
    r <- orig(...)
    r$par[] <- coe_stall_theta[-5L]
    r$value <- -a$logl_fn(r$par, a$time, a$status, a$time_lower,
                          a$time_upper, a$x, weights = a$weights)
    r$convergence <- 0L
    r
  })
  w <- list()
  fit <- withCallingHandlers(coe_repro(), warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd
    invokeRestart("muffleWarning")
  })
  # The premise: the fit is at the stall point.
  expect_equal(fit$fit$objective, -206.7432, tolerance = 1e-6)
  rec <- coe_records(fit)
  expect_length(rec, 1L)
  expect_identical(rec[[1L]]$phase, "constant")
  expect_gt(rec[[1L]]$gain, 0.01)
  cls <- vapply(w, function(cnd) inherits(cnd, "hzr_coe_no_events_left"),
                logical(1))
  expect_identical(sum(cls), 1L)
  expect_true(inherits(w[[which(cls)]], "hzr_boundary"))
})

test_that("a point short of the boundary supremum is certified (#261)", {
  skip_on_cran()
  fit <- suppressWarnings(coe_repro())
  theta <- stats::setNames(coe_stall_theta, names(coef(fit)))
  value <- as.numeric(hzr_evaluate(fit, theta)$logLik)
  # The premise, as #261 measured it: the constant phase conserved, an
  # ordinary-looking optimum at -206.7432.
  expect_true(fit$spec$control$conserve_applied)
  expect_equal(value, -206.7432, tolerance = 1e-6)

  rec <- coe_certificate(fit, theta, value)
  expect_identical(rec$mechanism, "coe_no_events_left")
  expect_identical(rec$phase, "constant")
  expect_identical(rec$parameter, "log_mu")
  # The certificate is the fit's own likelihood, recomputed from the point,
  # and it is where #261 put the supremum: above -206.7432, at -206.667.
  again <- as.numeric(hzr_evaluate(fit, rec$certificate_theta)$logLik)
  expect_equal(again, rec$certificate_loglik, tolerance = 1e-10)
  expect_gt(again, value + 0.01)
  expect_lt(abs(again - (-206.669)), 1e-3)
  # The conserved phase is switched off there; only intercepts moved.
  expect_lt(rec$certificate_theta[["constant.log_mu"]], -600)
  shapes <- setdiff(names(theta), c("early.log_mu", "constant.log_mu",
                                    "late.log_mu"))
  expect_identical(rec$certificate_theta[shapes], theta[shapes])

  # A certificate carries no `warned_by`, so hazard() announces it.
  expect_null(rec$warned_by)
  cnd <- .hzr_boundary_condition(list(rec))
  expect_s3_class(cnd, "hzr_coe_no_events_left")
  expect_s3_class(cnd, "hzr_boundary")
  expect_match(conditionMessage(cnd), "'constant'", fixed = TRUE)
  expect_match(conditionMessage(cnd), "starting value", fixed = TRUE)
})

test_that("a CoE fit that reaches the boundary is recorded, and warned once (#565)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard", envir = environment())
  fit_with <- function(conserve) {
    w <- list()
    fit <- withCallingHandlers(
      hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
             dist = "multiphase",
             phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
                           constant = hzr_phase("constant"),
                           late = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1,
                                            eta = 1)),
             fit = TRUE,
             # The fitted share is about 4e-9; a threshold well above it keeps
             # the fixture off the default's edge.
             control = list(n_starts = 1, conserve = conserve,
                            phase_share_tol = 1e-5)),
      warning = function(cnd) {
        w[[length(w) + 1L]] <<- cnd
        invokeRestart("muffleWarning")
      })
    list(fit = fit, w = w)
  }
  on <- fit_with(TRUE)
  fit <- on$fit
  expect_true(fit$spec$control$conserve_applied)
  # Past the point #261 measured, towards the boundary.
  expect_gt(fit$fit$objective, -206.7432 + 0.1)
  share <- fit$fit$phase_share
  expect_lt(share$share[share$phase == "constant"], 1e-5)

  rec <- coe_records(fit, "coe_phase_vanished")
  expect_length(rec, 1L)
  # One mechanism, one record shape: this is not the certificate's.
  expect_length(coe_records(fit), 0L)
  rec <- rec[[1L]]
  expect_null(rec$certificate_theta)
  expect_identical(rec$phase, "constant")
  expect_identical(rec$parameter, "log_mu")
  expect_identical(rec$share, share$share[share$phase == "constant"])
  expect_identical(rec$tol, 1e-5)
  expect_identical(rec$warned_by, "phase_share")
  expect_match(rec$detail, "'constant'", fixed = TRUE)

  # The identifiability warning says it, once, and nothing says it again.
  msgs <- vapply(on$w, conditionMessage, character(1))
  expect_identical(
    sum(grepl("Phase 'constant' contributes at most", msgs, fixed = TRUE)), 1L
  )
  boundary <- vapply(on$w, function(cnd) inherits(cnd, "hzr_boundary"),
                     logical(1))
  expect_identical(sum(boundary), 0L)
  # Report only.
  expect_equal(as.numeric(hzr_evaluate(fit, coef(fit))$logLik),
               fit$fit$objective, tolerance = 1e-8)

  # Without conservation no phase is the conserved one: nothing to record.
  off <- fit_with(FALSE)
  expect_false(off$fit$spec$control$conserve_applied)
  expect_length(coe_records(off$fit, "coe_phase_vanished"), 0L)
})

test_that("only the conserved phase, under the threshold, is recorded (#565)", {
  shares <- data.frame(phase = c("early", "constant", "late"),
                       share = c(0.9, 4e-9, 2e-12), variation = NA_real_)
  rec <- .hzr_coe_vanished_record(shares, "constant", 1e-8)
  expect_identical(rec$mechanism, "coe_phase_vanished")
  expect_identical(rec$phase, "constant")
  expect_identical(rec$share, 4e-9)
  # Another phase being absent is the identifiability check's business.
  expect_null(.hzr_coe_vanished_record(shares, "early", 1e-8))
  # At or above the threshold, not measurable, or the check switched off.
  expect_null(.hzr_coe_vanished_record(shares, "constant", 4e-9))
  expect_null(.hzr_coe_vanished_record(shares, "constant", 0))
  shares$share[2] <- NA_real_
  expect_null(.hzr_coe_vanished_record(shares, "constant", 1e-8))
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
  # At the point short of the boundary, not at the fit's own estimates,
  # which since #565 lie past it.
  theta <- coe_stall_theta
  value <- as.numeric(hzr_evaluate(fit, theta)$logLik)
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

  # Nor where the phase's cumulative hazard underflows to 0 but its event
  # hazard does not: a linear predictor of 9300 leaves mu = exp(-700), whose
  # cumulative hazard at times near 1e-30 is 0 while its hazard is not.
  # The other phases are switched off too, so nothing masks that hazard.
  xl2 <- list(early = NULL, constant = matrix(9300, length(d$time), 1L),
              late = NULL)
  th2 <- th1
  th2[c(pos1[["early"]], pos1[["late"]])] <- -1e4
  expect_null(.hzr_coe_boundary_record(
    th2, value, TRUE, function(p) value + 1, "constant", pos1[["constant"]],
    pos1, phases, counts1, xl2, d$time * 1e-30, d$status, NULL, NULL,
    sum(d$status == 1)))

  # And the converse: at very late times the same phase's hazard is lost
  # beside the others', but its cumulative hazard is not.
  expect_null(.hzr_coe_boundary_record(
    th1, value, TRUE, function(p) value + 1, "constant", pos1[["constant"]],
    pos1, phases, counts1, xl2, d$time * 1e300, d$status, NULL, NULL,
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
