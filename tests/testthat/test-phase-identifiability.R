# A phase can leave the model without leaving the output. These tests pin the
# two ways it happens, which have DIFFERENT consequences:
#
#   absent    -- the phase never starts; mu and shape are both unidentified.
#   saturated -- the phase finished before the first observation; it is then a
#                constant offset, so mu stays well identified and only the
#                shape parameters go flat.
#
# The second is the one that misled a parity investigation (temporal_hazard#143,
# where the report had it as "mu stops being identified" -- mu is in fact the
# one parameter that survives).

ident_ns <- function(f) get(f, asNamespace("TemporalHazard"))

ident_setup <- function(t_half, nu = 0, m = -0.4) {
  ph <- ident_ns(".hzr_validate_phases")(
    list(early    = hzr_phase("cdf", t_half = t_half, nu = nu, m = m),
         constant = hzr_phase("constant")))
  list(phases = ph,
       cc     = c(early = 0L, constant = 0L),
       xl     = list(early = NULL, constant = NULL),
       theta  = c(log(0.045), log(t_half), nu, m, log(0.036)))
}

ident_time <- function() {
  data("avc", package = "TemporalHazard", envir = environment())
  t <- avc$int_dead
  t[t > 0]
}

ident_check <- function(t_half, nu = 0, m = -0.4, tol = 1e-8) {
  s <- ident_setup(t_half, nu, m)
  ident_ns(".hzr_check_phase_identifiability")(
    s$theta, ident_time(), s$phases, s$cc, s$xl, tol = tol)
}


test_that("a healthy two-phase fit raises nothing", {
  # t_half = 0.003 is temporal_hazard#143's own reproducer. The early phase is
  # still climbing across the observed range, so nothing is unidentified -- the
  # guard must not fire here, or it would confirm the report it was written to
  # correct.
  expect_no_warning(ident_check(0.003))
  sh <- suppressWarnings(ident_check(0.003))
  expect_gt(sh$variation[sh$phase == "early"], 0.5)
})

test_that("a saturated phase warns, and says mu is the survivor", {
  w <- capture_warnings(ident_check(1e-6))
  expect_length(w, 1L)
  expect_match(w, "constant across the observed times")
  expect_match(w, "'mu' remains identified but the shape parameters do not")
  expect_match(w, "'early'")
})

test_that("saturation is detected as exactly zero variation", {
  sh <- suppressWarnings(ident_check(1e-6))
  expect_identical(sh$variation[sh$phase == "early"], 0)
  # ...and mu is NOT flagged absent: the phase still supplies most of Lambda.
  expect_gt(sh$share[sh$phase == "early"], 0.5)
})

test_that("an absent phase warns about mu AND shape", {
  # Half-life far beyond follow-up: G(t) never lifts off, so the phase
  # contributes essentially none of Lambda anywhere.
  w <- capture_warnings(ident_check(1e8))
  expect_true(any(grepl("has not started by the end of follow-up", w)))
  expect_true(any(grepl("neither its 'mu' nor its shape", w)))
})

test_that("an absent phase is reported once, as absent rather than as flat", {
  # A phase that never starts is trivially also flat. Reporting both would be
  # two warnings for one condition, and "absent" is the more informative.
  w <- capture_warnings(ident_check(1e8))
  expect_length(w, 1L)
})

test_that("the tolerance is honoured in both directions", {
  expect_no_warning(ident_check(1e-6, tol = 0))       # nothing can fall below 0
  expect_warning(ident_check(0.003, tol = 0.9),       # 0.843 < 0.9
                 "constant across the observed times")
})

test_that("variation is withheld, not guessed, for a phase carrying covariates", {
  # With covariates mu differs by row, so the spread of the contribution mixes
  # covariate variation with the shape's. Reporting it would invite reading
  # covariate spread as shape identifiability.
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- ident_ns(".hzr_validate_phases")(
    list(early    = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4),
         constant = hzr_phase("constant")))
  xm <- matrix(as.numeric(seq_len(nrow(avc)) %% 2), ncol = 1,
               dimnames = list(NULL, "grp"))
  sh <- ident_ns(".hzr_phase_shares")(
    c(log(0.045), log(1e-6), 0, -0.4, 0.3, log(0.036)),
    avc$int_dead,
    ph, c(early = 1L, constant = 0L), list(early = xm, constant = NULL))
  expect_true(is.na(sh$variation[sh$phase == "early"]))
  expect_false(is.na(sh$share[sh$phase == "early"]))
})

test_that("the fit carries the shares so the warning is checkable", {
  skip_if_not_installed("survival")
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- list(early = hzr_phase("cdf", t_half = 0.003, nu = 0, m = -0.4,
                               fixed = c("nu", "m", "t_half")),
             constant = hzr_phase("constant"))
  f <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
              dist = "multiphase", phases = ph,
              theta = c(log(0.045), log(0.003), 0, -0.4, log(0.036)),
              fit = TRUE, control = list(n_starts = 1, conserve = TRUE))
  expect_s3_class(f$fit$phase_share, "data.frame")
  expect_setequal(f$fit$phase_share$phase, c("early", "constant"))
})


# --- control$phase_share_tol validation -------------------------------------
# NA is the dangerous input, not the obviously-wrong one: `shares < NA` is NA,
# which() drops it, and the guard would pass every phase while looking like it
# ran. A check written against silent failure must not fail silently itself.

ident_fit_tol <- function(tol) {
  data("avc", package = "TemporalHazard", envir = environment())
  ph <- list(early = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4,
                               fixed = c("nu", "m", "t_half")),
             constant = hzr_phase("constant"))
  hazard(survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
         phases = ph, theta = c(log(0.045), log(1e-6), 0, -0.4, log(0.036)),
         fit = TRUE,
         control = list(n_starts = 1, conserve = TRUE, phase_share_tol = tol))
}

test_that("a bad phase_share_tol is rejected, naming the option", {
  skip_if_not_installed("survival")
  for (bad in list(NA_real_, NA, "1e-8", -1, c(1e-8, 1e-9), numeric(0), Inf)) {
    # fixed = TRUE: the option name contains a `$`, and escaping it for a
    # regex is exactly the kind of quiet mismatch this test exists to catch.
    expect_error(ident_fit_tol(bad), "'control$phase_share_tol'",
                 fixed = TRUE, info = deparse(bad))
  }
})

test_that("phase_share_tol = 0 silences the check rather than erroring", {
  skip_if_not_installed("survival")
  expect_no_warning(ident_fit_tol(0))
})

test_that("a valid phase_share_tol still warns", {
  skip_if_not_installed("survival")
  expect_warning(ident_fit_tol(1e-8), "constant across the observed times")
})


# --- return type is stable --------------------------------------------------

test_that("shares are a data.frame even when nothing can be measured", {
  # No observed time carries a positive total, so there is nothing to measure.
  # The frame must still have the same shape -- fit$phase_share is one type for
  # a caller to handle, not two -- and NA must mean "not measured" rather than
  # being confusable with a measured zero.
  s <- ident_setup(0.003)
  sh <- ident_ns(".hzr_phase_shares")(s$theta, numeric(0), s$phases, s$cc, s$xl)
  expect_s3_class(sh, "data.frame")
  expect_identical(sh$phase, c("early", "constant"))
  expect_true(all(is.na(sh$share)))
  expect_true(all(is.na(sh$variation)))
})

test_that("unmeasurable shares warn about nothing", {
  s <- ident_setup(0.003)
  expect_no_warning(
    ident_ns(".hzr_check_phase_identifiability")(
      s$theta, numeric(0), s$phases, s$cc, s$xl))
})


# ---------------------------------------------------------------------------
# Degenerate observed times (#211).
#
# `variation` collapsing does NOT mean the times are degenerate: a saturated
# phase is flat across times that are perfectly well spread. The condition is
# therefore asked of the TIMES -- their own relative range -- and not inferred
# from the phases. Two earlier attempts inferred it, and each reported a cause
# that had not occurred, which is the defect #211 is about.
# ---------------------------------------------------------------------------

ident_warn <- function(times, tol = 1e-8, other_times = NULL,
                       t_half = 0.5) {
  s <- ident_setup(t_half = t_half)
  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      s$theta, times, s$phases, s$cc, s$xl, tol = tol,
      other_times = other_times),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  w
}

# One saturated phase and one absent phase over well-spread times. Both
# diagnostics are correct here and neither is about the times.
ident_mixed <- function(times, other_times = NULL) {
  ph <- ident_ns(".hzr_validate_phases")(
    list(early = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4),
         late  = hzr_phase("cdf", t_half = 1e12, nu = 0, m = -0.4)))
  th <- c(log(0.045), log(1e-6), 0, -0.4, log(1e-12), log(1e12), 0, -0.4)
  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      th, times, ph, c(early = 0L, late = 0L),
      list(early = NULL, late = NULL), other_times = other_times),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  w
}

test_that("tied observation times are reported as their own condition", {
  w <- ident_warn(rep(2, 20))

  expect_length(w, 1L)
  expect_match(w, "span a relative range below")

  # The three claims that were wrong for this input.
  expect_no_match(w, "shape parameters")
  expect_no_match(w, "finished before the first observation")
  expect_no_match(w, "half-life")
})

test_that("a single-row fit reports the same condition", {
  w <- ident_warn(2)
  expect_length(w, 1L)
  expect_match(w, "span a relative range below")
})

test_that("near-tied times reach the same condition as exact ties", {
  # Counting distinct times missed this: `unique()` separates them while the
  # measures are degenerate, and it produced #211's wording verbatim.
  w <- ident_warn(c(100, 100.000001))
  expect_length(w, 1L)
  expect_match(w, "span a relative range below")
  expect_no_match(w, "shape parameters")
})

test_that("a saturated phase over well-spread times is still saturated", {
  # Regression guard. Inferring degeneracy from the phases made this fire the
  # tied-times message on 25 evenly spaced times -- and with a single phase,
  # "every phase is flat" IS "this phase saturated", so single-phase
  # saturation could never be reported correctly.
  ph <- ident_ns(".hzr_validate_phases")(
    list(early = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4)))
  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      c(log(0.045), log(1e-6), 0, -0.4), seq(1, 10, length.out = 25),
      ph, c(early = 0L), list(early = NULL)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "constant across the observed times")
  expect_no_match(w, "span a relative range below")
})

test_that("an absent phase over well-spread times does not mask a saturated one", {
  # `variation` is normalised by the phase's OWN maximum, so an absent phase
  # varies across nearly its full range (0.997 here) -- folding it into a
  # flatness verdict wrongly declared these 25 distinct times degenerate.
  w <- ident_mixed(seq(1, 10, length.out = 25))
  expect_length(w, 2L)
  expect_match(w, "contributes at most", all = FALSE)
  expect_match(w, "constant across the observed times", all = FALSE)
  expect_false(any(grepl("span a relative range below", w)))
})

test_that("an absent phase is reported even when only other times vary", {
  # `share` does not depend on how the times are spread, so the absent test is
  # never withheld -- a phase that has not started by the end of follow-up is
  # unidentified whatever the entry times do.
  w <- ident_mixed(rep(2, 20), other_times = seq(0.01, 1.9, length.out = 20))
  expect_match(w, "contributes at most", all = FALSE)
})

test_that("other evaluation times withhold the flatness verdict", {
  # A left-truncated fit evaluates Lambda(stop) - Lambda(start), so varying
  # entry times identify the shapes even when every `stop` is tied. Measures
  # taken over `time` alone cannot see that, so they are withheld rather than
  # guessed -- the same rule `variation` already follows for covariates.
  w <- ident_warn(rep(2, 20), other_times = seq(0.01, 1.9, length.out = 20))
  expect_length(w, 0L)
})

test_that("the verdict counts functionals against free parameters", {
  # With exits tied, k added evaluation points give k + 1 functionals of theta,
  # and identifying p free parameters needs k + 1 >= p. ident_setup() is the
  # two-phase model, p = 5, so four added points are needed. Two thresholds
  # were guessed before this: one point (silenced a constant entry time) and
  # two (off by a factor of two for exactly this model, in the silent
  # direction).
  expect_match(ident_warn(rep(2, 20), other_times = rep(0.5, 20)),
               "span a relative range below")                      # k = 1
  expect_match(ident_warn(rep(2, 20), other_times = c(0.5, 1.5)),
               "span a relative range below")                      # k = 2
  expect_match(ident_warn(rep(2, 20), other_times = c(0.5, 1.0, 1.5)),
               "span a relative range below")                      # k = 3
  expect_length(ident_warn(rep(2, 20),
                           other_times = c(0.4, 0.8, 1.2, 1.6)), 0L)  # k = 4

  # Points that add nothing do not count toward k. Lambda(0) is 0 for every
  # phase type, and a bound equal to a measured time restates it -- both are
  # synthesised for every row by Surv(type = "interval").
  # Four values of which one is excluded leaves k = 3, which is short; without
  # the exclusion it would be k = 4 and silent, so each case discriminates.
  expect_match(ident_warn(rep(2, 20), other_times = c(0, 0.4, 0.8, 1.2)),
               "span a relative range below")
  expect_match(ident_warn(rep(2, 20), other_times = c(2, 0.4, 0.8, 1.2)),
               "span a relative range below")
})

test_that("pinned shapes lower the bar, because they are not free", {
  # The count is of FREE parameters, not of theta's length. Pinning the early
  # phase's three shapes leaves p = 2 (the two log_mu), so a single added
  # evaluation point supplies the two functionals needed and nothing is
  # warned about -- where the same model with those shapes free would warn.
  ph <- ident_ns(".hzr_validate_phases")(
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 0, m = -0.4,
                           fixed = c("t_half", "nu", "m")),
         constant = hzr_phase("constant")))
  th <- c(log(0.045), log(0.5), 0, -0.4, log(0.036))
  cc <- c(early = 0L, constant = 0L)
  xl <- list(early = NULL, constant = NULL)
  expect_length(th, 5L)                                       # theta is long
  expect_identical(sum(ident_ns(".hzr_phase_free_mask")(ph, cc)), 2L)  # p is not

  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      th, rep(2, 20), ph, cc, xl, other_times = rep(0.5, 20)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 0L)

  # Same data, same one added point, shapes free: p = 5 and it warns.
  expect_match(ident_warn(rep(2, 20), other_times = rep(0.5, 20)),
               "span a relative range below")
})

test_that("a model whose mu is identified at one time point says nothing", {
  # A lone `constant` phase has p = 1, so the single functional the tied times
  # supply determines its mu exactly -- verified against the closed form
  # D/(nT). Warning here would be #211's defect reproduced by its own fix.
  skip_on_cran()  # end-to-end fit; CI runs it, --as-cran need not
  skip_if_not_installed("survival")
  set.seed(3)
  d <- data.frame(t = rep(2, 200), ev = rbinom(200, 1, 0.4))
  w <- character(0)
  f <- withCallingHandlers(
    hazard(survival::Surv(t, ev) ~ 1, data = d, dist = "multiphase",
           phases = list(base = hzr_phase("constant")), fit = TRUE,
           control = list(n_starts = 1)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("span a relative range below", w)))
  expect_equal(unname(exp(coef(f))), sum(d$ev) / (200 * 2), tolerance = 1e-5)
})

test_that("a model with no shape parameters is not told its shapes are flat", {
  # Two `constant` phases DO warrant the warning -- only mu1 + mu2 is
  # identified -- but they have no shapes, so the sentence must not name any.
  # That wording, on a shapeless model, is #211 verbatim.
  ph <- ident_ns(".hzr_validate_phases")(
    list(a = hzr_phase("constant"), b = hzr_phase("constant")))
  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      c(log(0.02), log(0.03)), rep(2, 20), ph, c(a = 0L, b = 0L),
      list(a = NULL, b = NULL)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "cannot be told apart")
  expect_no_match(w, "shape")
})

test_that("the same data written as interval-censored gets the same verdict", {
  skip_on_cran()  # end-to-end fit; CI runs it, --as-cran need not
  # End-to-end form of the above. Surv(type = "interval") synthesises
  # time_lower = 0 and time_upper = time, so a representation change that adds
  # no information silenced #211's own input.
  skip_if_not_installed("survival")
  d <- data.frame(t = rep(2, 40), ev = rep(1, 40))
  ph <- list(early = hzr_phase("cdf"), constant = hzr_phase("constant"))
  # Both through the FORMULA interface: only .hzr_parse_formula() synthesises
  # time_lower = 0 and time_upper = time, and it is that synthesis this test
  # exists for. On the vector path `unclass(status)[, 2L]` reads a Surv matrix
  # positionally, so both spellings would collapse to the same fit and the
  # test would compare a computation with itself.
  expect_warning(
    hazard(survival::Surv(t, ev) ~ 1, data = d, dist = "multiphase",
           phases = ph, fit = TRUE, control = list(n_starts = 1)),
    "span a relative range below")
  expect_warning(
    hazard(survival::Surv(t, t, ev, type = "interval") ~ 1, data = d,
           dist = "multiphase", phases = ph, fit = TRUE,
           control = list(n_starts = 1)),
    "span a relative range below")
})

test_that("a constant entry time does not rescue tied exit times", {
  skip_on_cran()  # end-to-end fit; CI runs it, --as-cran need not
  # Every row shares the same entry and the same exit, so the likelihood
  # depends on theta through two numbers only. Five free parameters, so the
  # shapes are unidentified and the warning is correct -- contrast the
  # 200-distinct-entry fit below, which is withheld.
  skip_if_not_installed("survival")
  set.seed(9)
  d <- data.frame(start = rep(1, 200), stop = rep(2, 200),
                  ev = rbinom(200, 1, 0.4))
  w <- character(0)
  withCallingHandlers(
    hazard(survival::Surv(start, stop, ev) ~ 1, data = d, dist = "multiphase",
           phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 0, m = -0.4),
                         constant = hzr_phase("constant")),
           theta = c(log(0.045), log(0.5), 0, -0.4, log(0.036)), fit = TRUE,
           control = list(n_starts = 1, conserve = FALSE)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("span a relative range below", w)))
})

test_that("degenerate times with an absent phase still warn", {
  # Pins `flat_all` rather than `flat` in the predicate: an absent phase is
  # excluded from `flat`, so using that set would let one dead phase block the
  # verdict. Previously this choice was covered only incidentally.
  w <- ident_mixed(rep(2, 20))
  expect_match(w, "span a relative range below", all = FALSE)
})

test_that("tol = 0 silences the check, as ?hazard documents", {
  expect_length(ident_warn(rep(2, 20), tol = 0), 0L)
  expect_length(ident_warn(c(100, 100.000001), tol = 0), 0L)
  # ...and the same inputs are NOT silent at the default tolerance, so those
  # assertions cannot pass merely because nothing would have warned anyway.
  expect_length(ident_warn(rep(2, 20)), 1L)
  expect_length(ident_warn(c(100, 100.000001)), 1L)
})

test_that("distinguishable times still get the ordinary saturated message", {
  w <- capture_warnings(ident_check(1e-6))
  expect_match(w, "constant across the observed times", all = FALSE)
  expect_false(any(grepl("span a relative range below", w)))
})

test_that("bunched times with a phase that still varies are not degenerate", {
  # A relative range below tol is NOT sufficient. Times spanning 1 unit at
  # t ~ 1e9 are bunched, but a steep enough g3 is fully identified across them:
  # moving gamma from 1e9 to 2e9 moves the log-likelihood by ~11 units. Calling
  # that "nothing separates one phase from another" is the same wrong answer
  # #211 is about, so the verdict requires BOTH a degenerate range and no phase
  # still varying. `const` is flat here and has no shape parameters, so nothing
  # is said at all.
  ph <- ident_ns(".hzr_validate_phases")(
    list(const = hzr_phase("constant"),
         late  = hzr_phase("g3", tau = 1e9 + 0.5, gamma = 1e9, alpha = 1, eta = 1)))
  th <- c(log(1e-12), 0, log(1e9 + 0.5), 1e9, 1, 1)
  tt <- seq(1e9, 1e9 + 1, length.out = 20)
  cc <- c(const = 0L, late = 0L)
  xl <- list(const = NULL, late = NULL)

  sh <- ident_ns(".hzr_phase_shares")(th, tt, ph, cc, xl)
  expect_lt(attr(sh, "time_variation"), 1e-8)      # the range IS degenerate
  expect_gt(sh$variation[sh$phase == "late"], 0.5)  # yet a phase varies

  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(th, tt, ph, cc, xl),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 0L)

  ll <- function(gamma) {
    ident_ns(".hzr_logl_multiphase")(
      theta = c(log(1e-12), 0, log(1e9 + 0.5), gamma, 1, 1), time = tt,
      status = rep(1, 20), weights = rep(1, 20), phases = ph,
      covariate_counts = cc, x_list = xl)
  }
  expect_gt(abs(ll(2e9) - ll(1e9)), 1)
})

test_that("a phase with no shape parameters is never told its shape is flat", {
  # #211's specific complaint. The filter is not dead code: a constant phase's
  # relative range tracks the times' own closely but not exactly in floating
  # point, so it can land on the other side of tol while the times do not.
  has <- ident_ns(".hzr_phase_has_shape")
  expect_identical(has(list(hzr_phase("constant"))), FALSE)
  expect_identical(has(list(hzr_phase("cdf"))), TRUE)
  expect_identical(has(list(hzr_phase("g3"))), TRUE)
  # Pinned shapes still EXIST, and the saturated message is true of them --
  # which is what the pre-existing phase_share_tol test relies on.
  expect_identical(
    has(list(hzr_phase("cdf", t_half = 1, nu = 0, m = 0,
                       fixed = c("t_half", "nu", "m")))),
    TRUE)
})

test_that("exactly tied times warn even when every phase carries covariates", {
  # `variation` is withheld for a covariate-carrying phase, so nothing is
  # measurable and any phase-derived verdict says nothing -- on the maximally
  # degenerate input. The times are the only evidence there is.
  set.seed(211)
  ph <- ident_ns(".hzr_validate_phases")(
    list(a = hzr_phase("cdf"), b = hzr_phase("constant")))
  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      c(log(0.045), log(0.5), 0, -0.4, 0.1, log(0.036), 0.1),
      rep(2, 20), ph, c(a = 1L, b = 1L),
      list(a = matrix(rnorm(20)), b = matrix(rnorm(20)))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "span a relative range below")
})

test_that("entry times withhold the degenerate verdict, not the saturated one", {
  # Withholding is justified only for the verdict measured over `time` alone.
  # Gating the per-phase message on it too meant adding a `start` column
  # silently removed a correct diagnostic.
  ph <- ident_ns(".hzr_validate_phases")(
    list(early = hzr_phase("cdf", t_half = 1e-6, nu = 0, m = -0.4),
         constant = hzr_phase("constant")))
  w <- character(0)
  withCallingHandlers(
    ident_ns(".hzr_check_phase_identifiability")(
      c(log(0.045), log(1e-6), 0, -0.4, log(0.036)),
      seq(1, 10, length.out = 25), ph, c(early = 0L, constant = 0L),
      list(early = NULL, constant = NULL),
      other_times = seq(0.1, 0.9, length.out = 25)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_match(w, "constant across the observed times", all = FALSE)
})

test_that("the time measure is taken over the measured rows only", {
  # Rows with a zero total are dropped by `ok`, and the measure follows them:
  # over all of `time` this spans a full relative range, over the measured
  # rows it is degenerate.
  s <- ident_setup(t_half = 0.5)
  sh <- ident_ns(".hzr_phase_shares")(
    s$theta, c(0, 0, 3, 3), s$phases, s$cc, s$xl)
  expect_identical(attr(sh, "time_variation"), 0)
})

test_that("all-zero times are a measured zero, not an unmeasured NA", {
  # Reachable rather than defensive: a g3 phase returns G3(0) > 0, so `ok` is
  # all TRUE and the maximum is 0. Dividing would give NaN.
  ph <- ident_ns(".hzr_validate_phases")(
    list(late = hzr_phase("g3"), constant = hzr_phase("constant")))
  sh <- ident_ns(".hzr_phase_shares")(
    c(log(0.045), log(1), 1, 1, 1, log(0.036)), rep(0, 10), ph,
    c(late = 0L, constant = 0L), list(late = NULL, constant = NULL))
  expect_identical(attr(sh, "time_variation"), 0)
})

test_that("the degenerate branch still returns the shares, without the measure", {
  # summary() and fit$phase_share read this; the internal time measure must
  # not travel out on it.
  s <- ident_setup(t_half = 0.5)
  sh <- suppressWarnings(ident_ns(".hzr_check_phase_identifiability")(
    s$theta, rep(2, 20), s$phases, s$cc, s$xl))
  expect_s3_class(sh, "data.frame")
  expect_identical(sh$phase, c("early", "constant"))
  expect_null(attr(sh, "time_variation"))
  expect_null(attr(sh, "time_unique"))
})

test_that("the check reaches the degenerate branch through hazard()", {
  skip_on_cran()  # end-to-end fit; CI runs it, --as-cran need not
  # Every other test here calls the internal directly, so a change to the
  # arguments at the call site would go uncaught.
  skip_if_not_installed("survival")
  d <- data.frame(t = rep(2, 40), ev = rep(c(1, 0), 20))
  expect_warning(
    hazard(survival::Surv(t, ev) ~ 1, data = d, dist = "multiphase",
           phases = list(early = hzr_phase("cdf"), constant = hzr_phase("constant")),
           fit = TRUE, control = list(n_starts = 1)),
    "span a relative range below")
})

test_that("interval bounds count as other times at the call site", {
  skip_on_cran()  # end-to-end fit; CI runs it, --as-cran need not
  # The call site gates each bound by the status that makes it live; a version
  # passing only time_lower is undetectable without an interval-censored case.
  #
  # The shapes are PINNED and theta supplied, so the fit is deterministic and
  # this asserts argument plumbing rather than where an unconstrained optimiser
  # happens to land. The earlier version used the bare defaults and failed only
  # on Windows, with "t_half must be a positive scalar" -- a platform-dependent
  # assertion masquerading as a check of the call site.
  skip_if_not_installed("survival")
  n <- 40
  d <- data.frame(lo = rep(2, n), up = seq(2.1, 6, length.out = n),
                  st = rep(2, n))
  ph <- list(
    early = hzr_phase("cdf", t_half = 1, nu = 0, m = -0.4,
                      fixed = c("t_half", "nu", "m")),
    constant = hzr_phase("constant"))
  w <- character(0)
  withCallingHandlers(
    hazard(time = rep(2, n), status = d$st, time_lower = d$lo,
           time_upper = d$up, dist = "multiphase", phases = ph,
           theta = c(log(0.05), log(1), 0, -0.4, log(0.03)),
           fit = TRUE, control = list(n_starts = 1)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("span a relative range below", w)))
})

test_that("a bound the objective never reads does not count", {
  # time_upper is read only for left-censored and interval rows, so on tied
  # right-censored data it changes nothing: .hzr_logl_multiphase() is
  # bit-identical with and without it. Counting it anyway silenced the warning
  # on data that needed it -- the silent failure this guard exists to prevent.
  ph <- ident_ns(".hzr_validate_phases")(
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 0, m = -0.4),
         constant = hzr_phase("constant")))
  th <- c(log(0.045), log(0.5), 0, -0.4, log(0.036))
  cc <- c(early = 0L, constant = 0L)
  xl <- list(early = NULL, constant = NULL)
  tt <- rep(2, 20)
  st <- rep(0, 20)
  up <- seq(0.2, 1.8, length.out = 20)

  ll <- function(u) {
    ident_ns(".hzr_logl_multiphase")(
      theta = th, time = tt, status = st, time_upper = u,
      weights = rep(1, 20), phases = ph, covariate_counts = cc, x_list = xl)
  }
  expect_identical(ll(NULL), ll(up))          # the objective cannot see it

  skip_if_not_installed("survival")
  saw <- function(u) {
    w <- character(0)
    withCallingHandlers(
      hazard(time = tt, status = st, time_upper = u, dist = "multiphase",
             phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 0, m = -0.4),
                           constant = hzr_phase("constant")),
             theta = th, fit = TRUE, control = list(n_starts = 1)),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    any(grepl("span a relative range below", w))
  }
  expect_true(saw(NULL))
  expect_true(saw(up))                         # ...so neither can the guard
})


test_that("a left-truncated fit through hazard() is not called unidentified", {
  skip_on_cran()  # end-to-end fit; CI runs it, --as-cran need not
  skip_if_not_installed("survival")
  set.seed(211)
  d <- data.frame(start = runif(200, 0.01, 1.9), stop = rep(2, 200),
                  ev = rbinom(200, 1, 0.4))
  w <- character(0)
  withCallingHandlers(
    hazard(survival::Surv(start, stop, ev) ~ 1, data = d, dist = "multiphase",
           phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 0, m = -0.4),
                         constant = hzr_phase("constant")),
           theta = c(log(0.045), log(0.5), 0, -0.4, log(0.036)), fit = TRUE,
           control = list(n_starts = 1, conserve = FALSE)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("span a relative range below", w)))
})
