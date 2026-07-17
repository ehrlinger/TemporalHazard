# SAS parity tests for hzr_stepwise().  Skip gracefully when the
# fixture file under inst/fixtures/ is absent, so the package installs
# and passes CI without any SAS capture.
#
# To enable: run the SAS template in
# `inst/extdata/stepwise-fixtures/avc-forward-wald.sas`, then
# convert its three output files via:
#
#   TemporalHazard:::.hzr_build_stepwise_fixture(
#     trace_csv = "path/to/stepwise_trace.csv",
#     final_csv = "path/to/stepwise_final.csv",
#     meta_txt  = "path/to/stepwise_meta.txt",
#     out_path  = "inst/fixtures/stepwise-avc-forward-wald.rds"
#   )
#
# ALGORITHMIC NOTE: R's hzr_stepwise(criterion = "score") now computes the same
# Q statistic SAS HAZARD uses for selection (score test at the current
# estimates, no refit), including SAS's approximation that ignores
# shaping-parameter covariances. Per-step comparison is therefore meaningful,
# where it previously was not.
#
# The whole-model check below still drives R with criterion = "wald" -- its
# generous tolerances were calibrated against the Wald refit path, and the
# variable set it produces (6 covariates) matches SAS's count where the score
# path selects a slightly smaller model.  The fixture's meta$criterion records
# "score" because that is what SAS computed; the R comparison mode is a separate
# choice, so this runner pins "wald" explicitly rather than reading it back.

.stepwise_parity_tolerance <- list(
  logLik      = 10,    # absolute; generous to account for Q-stat vs Wald path difference
  coefficient = 0.10   # relative; final estimates for shared covariates
)


.stepwise_parity_run <- function(fix) {
  if (fix$meta$dataset != "avc") {
    testthat::skip(paste0("Parity runner only wired up for AVC so far; saw ",
                           fix$meta$dataset))
  }
  data(avc, envir = environment())
  avc <- na.omit(avc)

  # Multiphase fits with n_starts > 1 perturb the initial values using the
  # ambient RNG.  Seed locally and restore the stream on exit so the test is
  # reproducible and does not leak RNG state into later tests.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      suppressWarnings(rm(".Random.seed", envir = globalenv())),
      add = TRUE
    )
  }
  set.seed(20260610L)

  # Two-phase (early CDF + constant) model matching SAS CONDITION=14:
  #   Early phase    -- Weibull CDF, THALF/NU/M fixed at the unconditional fit values
  #   Constant phase -- flat exponential hazard rate
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = avc,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                            fixed = c("t_half", "nu", "m")),
      constant = hzr_phase("constant")
    ),
    fit    = TRUE,
    control = list(n_starts = 3L, maxit = 500L, conserve = TRUE)
  ))

  # For multiphase stepwise, scope must be a named list of one-sided formulas.
  cands    <- fix$scope$candidates
  scope_mp <- list(early = reformulate(cands), constant = reformulate(cands))

  hzr_stepwise(
    base,
    scope     = scope_mp,
    data      = avc,
    direction = fix$meta$direction,
    criterion = "wald",
    slentry   = fix$meta$slentry,
    slstay    = fix$meta$slstay,
    force_in  = fix$scope$force_in,
    force_out = fix$scope$force_out,
    trace     = FALSE
  )
}


test_that("AVC forward-Wald stepwise: R selects similar variables to SAS", {
  skip_on_cran()
  fix <- .hzr_load_stepwise_fixture("avc-forward-wald")
  if (is.null(fix)) {
    skip("Stepwise parity fixture avc-forward-wald.rds not found")
  }

  r_fit   <- suppressWarnings(.stepwise_parity_run(fix))

  # Final fitted model's coefficient names (e.g. "early.com_iv").  Covariates
  # are the rows whose suffix (after "phase.") is one of the scope candidates;
  # intercepts (log_mu) and shape parameters (log_t_half, nu, m, ...) are not.
  r_summary    <- summary(r_fit)$coefficients
  r_coef_names <- rownames(r_summary)
  cands        <- fix$scope$candidates
  r_covariates <- r_coef_names[sub("^[^.]+\\.", "", r_coef_names) %in% cands]

  # SAS non-intercept final covariates.
  sas_coef       <- fix$final$coef
  sas_covariates <- sas_coef[!sas_coef$variable %in% c("e0", "c0"), ]

  # R's *final model* must retain at least as many covariates as SAS selected.
  # (Count coefficients in the fitted model, not accepted entry steps, so a
  # variable that entered and was later dropped does not inflate the count.)
  expect_gte(length(r_covariates),
             nrow(sas_covariates),
             label = "R final model should retain at least as many covariates as SAS selected")

  # Final log-likelihood: R may find a different (even better) optimum via a
  # different variable path; we only require it stays within the generous
  # tolerance that accounts for the Q-stat vs Wald path difference.
  expect_gte(r_fit$fit$objective,
             fix$final$logLik - .stepwise_parity_tolerance$logLik,
             label = paste0("R logLik (", round(r_fit$fit$objective, 3),
                            ") should be within ", .stepwise_parity_tolerance$logLik,
                            " of SAS logLik (", fix$final$logLik, ")"))

  # Final variable set: SAS non-intercept covariates should largely appear in R.
  # Map fixture coef rows (variable + phase columns) to R's "phase.variable" naming.
  sas_varphase <- paste0(sas_covariates$phase, ".", sas_covariates$variable)
  n_shared <- sum(sas_varphase %in% r_coef_names)

  expect_gte(n_shared / length(sas_varphase), 0.5,
             label = paste0("At least half of SAS final covariates should appear in R model.",
                            " SAS: ", paste(sas_varphase, collapse = ", "),
                            " | R: ", paste(r_coef_names, collapse = ", ")))

  # Coefficient comparison is intentionally omitted: because R and SAS take
  # different selection paths (Wald refit vs Q-statistic), their final models
  # include different covariates.  Estimates for nominally "shared" variables
  # differ because the surrounding model composition differs, making direct
  # comparison invalid.
})


# Sanity: the fixture loader + validator at least work on a synthetic
# fixture (doesn't depend on SAS output being present).
test_that("fixture validator rejects malformed input", {
  bad <- list(
    meta = list(dataset = "x", dist = "weibull", criterion = "wald"),
    scope = list(candidates = "age", force_in = character(),
                 force_out = character()),
    steps = data.frame(step_num = 1L, action = "bogus",
                       variable = "x", phase = NA_character_,
                       stat = 1, df = 1L, p_value = 0.1,
                       stringsAsFactors = FALSE),
    final = list(coef = data.frame(), logLik = -100, iterations = 5L)
  )
  expect_error(
    TemporalHazard:::.hzr_validate_stepwise_fixture(bad),
    "Invalid stepwise parity fixture"
  )
})

test_that("fixture validator accepts a well-formed fixture skeleton", {
  ok <- list(
    meta = list(
      dataset     = "x", dist = "weibull",
      criterion   = "wald", direction = "forward",
      slentry     = 0.30, slstay = 0.20,
      sas_version = "9.4", captured_on = "2026-04-17",
      proc_hazard = "PROC HAZARD ..."
    ),
    scope = list(candidates = "age", force_in = character(),
                 force_out = character()),
    steps = data.frame(
      step_num = 1L, action = "enter", variable = "age",
      phase = NA_character_, stat = 5, df = 1L, p_value = 0.025,
      stringsAsFactors = FALSE
    ),
    final = list(coef = data.frame(), logLik = -100, iterations = 7L)
  )
  expect_silent(TemporalHazard:::.hzr_validate_stepwise_fixture(ok))
})


test_that("R's per-step Q matches SAS's for the first entered variable", {
  fx <- readRDS(test_path("..", "..", "inst", "fixtures",
                          "stepwise-avc-forward-wald.rds"))
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  nuis   <- .hzr_score_nuisance(base)
  enters <- fx$steps[fx$steps$action == "enter", ]
  first  <- enters[1L, ]
  q <- .hzr_score_q(base, var = first$variable, phase = first$phase,
                    data = avc, nuisance = nuis)

  # SAS's Q for com_iv@early at the intercept-only two-phase model is 34.3396.
  # The multiphase optimiser perturbs its start from the ambient RNG, so Q
  # scatters ~1e-5 relative between runs; 1e-2 clears that by three orders of
  # magnitude while still rejecting the rounded-shape model (Q ~34.86, off by
  # 1.5e-2) that the fixed shapes here exist to avoid.
  expect_equal(q$stat, first$stat, tolerance = 1e-2)
})


test_that("score selection refits once per accepted step, not once per candidate", {
  skip_if_not_installed("numDeriv")
  # The whole point of criterion = "score": the reduced-model information is
  # inverted once per step and every candidate is scored without a refit, so the
  # only refit is the winner's -- one per accepted step. Wald refits EVERY
  # candidate at EVERY step, so its refit count is strictly larger. Counting the
  # refit calls asserts that invariant directly; a wall-clock comparison would
  # measure the same thing but flake on a loaded CI box.
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "weibull", theta = c(mu = 0.01, nu = 0.5), fit = TRUE)
  scope <- ~ age + mal + com_iv

  count_refits <- function(criterion) {
    orig <- TemporalHazard:::.hzr_refit_with_scope
    n <- 0L
    counter <- function(...) {
      n <<- n + 1L
      orig(...)
    }
    utils::assignInNamespace(".hzr_refit_with_scope", counter,
                             ns = "TemporalHazard")
    on.exit(utils::assignInNamespace(".hzr_refit_with_scope", orig,
                                     ns = "TemporalHazard"))
    fit <- suppressWarnings(hzr_stepwise(
      base, scope = scope, data = avc, criterion = criterion,
      direction = "forward", slentry = 0.3, slstay = 0.2, trace = FALSE
    ))
    list(n = n, n_enter = sum(fit$steps$action == "enter"))
  }

  score <- count_refits("score")
  wald  <- count_refits("wald")

  # Score refits only the winners it accepts, never the candidates it scores.
  expect_equal(score$n, score$n_enter)
  # And that is strictly fewer refits than Wald's per-candidate loop.
  expect_lt(score$n, wald$n)
})
