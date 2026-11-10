## Behavioural census: the fixed case battery.
##
## Sourced by census-run.R under EACH installed version. Nothing in here may
## depend on the package version: the data is generated locally from fixed
## seeds rather than taken from the package's own datasets, because a shipped
## dataset that changed between versions would make every downstream
## difference unattributable.
##
## A case returns a PROBE: a small list of plain numerics. Whole fit objects
## are never compared, because `#242` added a "what this fit did not do" field
## to every object and `identical()` on the object would then report a
## difference for every single case -- a result that is technically true and
## says nothing. The probe is the thing a user reads: theta, vcov, the
## objective, and the prediction.
##
## Each case declares `needs`: the exported names it requires. A case whose
## `needs` are absent under a version is recorded as ABSENT, which is a
## different outcome from an error, and is never silently dropped.
##
## ---------------------------------------------------------------------------
## KNOWN LIMITATION OF THIS BATTERY: the multiphase cases run on data that
## cannot identify a multiphase model.
## ---------------------------------------------------------------------------
##
## `census_data()$basic` is plain exponential draws with no early/late/
## background structure. Every multiphase case below therefore FAILS the SAS/C
## relative-gradient test -- `mp_cov_formula` stops at rel_gradient 0.0181
## against a required 6.06e-06, about 3000x outside tolerance, and across
## `start_seed` 1:10, 42 and 2026 not one of 12 fits meets it. Several cases
## also have a phase contributing 1e-9 or less of the cumulative hazard, and
## the fits warn about exactly that.
##
## The consequence for reading this census, which cost one withdrawn finding
## to learn: a multiphase DIFFERS row means **this fit's stopping point
## moved**. It does NOT mean the optimum moved, and it is NOT evidence about
## the likelihood. Comparing two fits that both fail the convergence criterion
## compares two arbitrary stopping points, so neither "better" nor "worse" is
## a property either one has.
##
## The single-distribution cases, the exported-pure-function cases and the
## predict cases are unaffected: those fits converge, and their verdicts carry
## their ordinary meaning.
##
## TO FIX, before the multiphase rows are leaned on: rebuild the multiphase
## cases on data with genuine phase structure so the fits converge, then a
## multiphase difference will mean something. Left undone on purpose -- it
## changes the fixed case battery, and the value of a fixed battery is that it
## does not drift mid-programme. Raise it with the maintainer first.

# ---------------------------------------------------------------------------
# Deterministic data
# ---------------------------------------------------------------------------

census_data <- function() {
  d <- list()

  set.seed(1001)
  n <- 150
  z <- stats::rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  d$basic <- data.frame(
    time = stats::rexp(n, rate = 0.3) + 0.01,
    status = stats::rbinom(n, 1L, 0.7),
    age = stats::runif(n, 40, 80),
    z = z,
    g = g
  )

  # Left truncation: a real entry-time cohort. Some rows enter at the origin,
  # most do not, so the entry term cannot be a no-op (#253's known positive).
  set.seed(253)
  n <- 200
  z2 <- stats::rnorm(n)
  entry <- stats::runif(n, 0, 2)
  entry[1:40] <- 0
  rate <- 0.3 * exp(0.5 * z2)
  t_event <- entry + stats::rexp(n, rate)
  t_cens <- entry + stats::runif(n, 0, 6)
  d$trunc <- data.frame(
    time = pmin(t_event, t_cens),
    status = as.integer(t_event <= t_cens),
    entry = entry,
    z = z2
  )

  # Interval censoring: status 2 rows with a genuine bracket.
  set.seed(1003)
  n <- 120
  lo <- stats::runif(n, 0.1, 3)
  wid <- stats::runif(n, 0.2, 1.5)
  st <- c(rep(1L, 40), rep(0L, 40), rep(2L, 40))
  d$interval <- data.frame(
    time = lo + wid,
    time_lower = lo,
    time_upper = lo + wid,
    status = st,
    z = stats::rnorm(n)
  )

  # Weights: integer-ish case weights, not all 1.
  set.seed(1004)
  d$weighted <- d$basic
  d$weighted$w <- sample(c(1, 1, 2, 3), nrow(d$basic), replace = TRUE)

  # -------------------------------------------------------------------------
  # APPENDED 2026-09-17: data with GENUINE phase structure.
  #
  # The `basic` cohort above is plain exponential draws and cannot identify a
  # multiphase model, which is why every original multiphase case fails the
  # gradient test. These cohorts are simulated from an actual multiphase
  # cumulative hazard, inverted numerically:
  #
  #   H(t) = mu_e * (1 - exp(-t / tau_e))  +  mu_bg * t
  #          \_______ early, saturating _/    \_ background _/
  #
  # so the early hazard is (mu_e/tau_e) * exp(-t/tau_e), decaying on timescale
  # tau_e, against a flat background: two genuinely separated timescales. The
  # parameters were chosen by sweeping n, mu_e, tau_e and t_half and keeping a
  # configuration whose fit MEETS SAS/C's relative-gradient test, rather than
  # assumed to work.
  # -------------------------------------------------------------------------
  sim_phased <- function(n, mu_e, tau_e, mu_bg, fu, seed) {
    set.seed(seed)
    grid <- seq(0, fu * 3, length.out = 200000)
    hh <- mu_e * (1 - exp(-grid / tau_e)) + mu_bg * grid
    u <- stats::runif(n)
    # hh is increasing in grid, so approx() inverts it directly.
    t_event <- stats::approx(hh, grid, xout = -log(u), rule = 2)$y
    t_cens <- stats::runif(n, 0.5 * fu, 1.5 * fu)
    data.frame(
      time = pmin(t_event, t_cens),
      status = as.integer(t_event <= t_cens),
      z = stats::rnorm(n),
      w = sample(c(1, 1, 2, 3), n, replace = TRUE)
    )
  }
  # Meets the gradient test on main (rel_gradient 1.06e-06) and passes the
  # version-independent optimum check on v1.2.9.
  d$phased <- sim_phased(400L, 0.8, 0.3, 0.05, 15, seed = 7001L)
  # Three phases on two timescales: deliberately NOT identifiable, and kept as
  # the GATE'S OWN known positive -- see mpc_3phase_unidentified below.
  d$phased3 <- sim_phased(800L, 1.2, 0.2, 0.04, 20, seed = 7003L)

  # Stepwise candidate pool: a few real signals and a few null ones.
  set.seed(1005)
  n <- 220
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- stats::rnorm(n)
  x4 <- stats::rnorm(n)
  lp <- 0.8 * x1 - 0.5 * x2
  d$step <- data.frame(
    time = stats::rexp(n, rate = 0.2 * exp(lp)) + 0.01,
    status = stats::rbinom(n, 1L, 0.75),
    x1 = x1, x2 = x2, x3 = x3, x4 = x4
  )

  # -------------------------------------------------------------------------
  # APPENDED for the 1.3.1 cut (#683): the limit fixtures.
  #
  # `avc` and `valves` are the package's own datasets, which the header above
  # warns against: a dataset that changed between versions would make every
  # difference on these cases unattributable. They are used because the #448
  # and #729 fixtures ARE those datasets, and the `limit_fixture_data` case
  # compares a checksum of every column read, so a changed dataset shows as
  # a DIFFERS row of its own instead of hiding inside a fit.
  # -------------------------------------------------------------------------
  e <- new.env()
  utils::data(list = c("avc", "valves"), package = "TemporalHazard", envir = e)
  k <- is.finite(e$avc$int_dead) & is.finite(e$avc$dead) & e$avc$int_dead > 0
  d$avc <- data.frame(time = e$avc$int_dead[k], status = e$avc$dead[k])
  v <- e$valves
  v$nyha[is.na(v$nyha)] <- mean(v$nyha, na.rm = TRUE)
  v$mv_nyha <- v$mitral * v$nyha
  d$valves <- v
  # The #444 / #683 test cohort (test-unbounded-phase.R, ub_data()).
  set.seed(1)
  d$ub <- data.frame(t = stats::rweibull(200, 1.4, 5),
                     s = stats::rbinom(200, 1, 0.8))
  # The weak-only ridge cohort (test-limit-restart.R, weak_only_fit()).
  set.seed(42)
  d$ridge <- data.frame(time = stats::rexp(100, rate = 0.5) + 0.01,
                        status = sample(0:1, 100, replace = TRUE,
                                        prob = c(0.2, 0.8)))

  d
}

# ---------------------------------------------------------------------------
# Probe builders -- what actually gets compared
# ---------------------------------------------------------------------------

limit_trigger_text <- function(trig) {
  if (!is.list(trig)) return(character())
  vapply(trig, function(r) {
    paste0(r$phase, ".", r$parameter, " (", r$source, ")")
  }, character(1))
}

limit_surface <- function(fit) {
  # The #683 surface of a fit: `fit$fit$limit` and the printed `Limit check:`
  # line. Kept under the `limit` key, which census-compare.R drops before the
  # bit-identity comparison (as it drops `gate`) and reports in a section of
  # its own. Were it compared, every multiphase fit would read DIFFERS on the
  # new side for gaining a field, and "did the well-identified fits stay
  # bit-identical" could not be answered.
  #
  # `fit$fit$limit` has four states, and they are kept apart: the element is
  # ABSENT (a version before #683), NULL (checked, nothing at a limit; the
  # branch stores it with `["limit"] <- list(NULL)` so the name is present),
  # NA (not examined), or a list whose `outcome` is "restarted" or "kept".
  # "absent" is never read as a match for any of the others.
  present <- "limit" %in% names(fit$fit)
  lim <- fit$fit$limit
  txt <- tryCatch(utils::capture.output(print(fit)),
                  error = function(e) paste("print-error:", conditionMessage(e)))
  is_line <- grepl("Limit check:", txt, fixed = TRUE)
  is_rec <- present && is.list(lim)
  list(
    field = if (!present) {
      "absent"
    } else if (is.null(lim)) {
      "NULL (checked, nothing at a limit)"
    } else if (!is.list(lim)) {
      "NA (not examined)"
    } else {
      paste0("list, outcome ", lim$outcome)
    },
    outcome = if (is_rec) lim$outcome else NA_character_,
    trigger = if (is_rec) limit_trigger_text(lim$trigger) else character(),
    grid_rows = if (is_rec && is.data.frame(lim$grid)) nrow(lim$grid) else NA_integer_,
    winner = if (is_rec) lim$winner else NA_integer_,
    original_objective = if (is_rec) lim$original$objective else NA_real_,
    print_limit = if (any(is_line)) trimws(txt[is_line]) else "absent",
    # Everything print() shows apart from the limit line, so the report can
    # say whether print() changed in any other way.
    print_other = txt[!is_line]
  )
}

probe_fit <- function(fit) {
  # Only the numeric result surface. Not the object, not the call, not the
  # degraded-reasons list: those changed by design this cycle. The #683
  # surface rides along under `limit`, which is reported, not compared.
  v <- tryCatch(stats::vcov(fit), error = function(e) paste("vcov-error:", conditionMessage(e)))
  list(
    kind = "fit",
    theta = fit$fit$theta,
    theta_names = names(fit$fit$theta),
    objective = fit$fit$objective,
    converged = fit$fit$converged,
    se = fit$fit$se,
    vcov = v,
    rcond = fit$fit$rcond,
    limit = limit_surface(fit)
  )
}

probe_fit_gated <- function(fit) {
  # probe_fit() plus a `gate` element carrying optimum-quality evidence.
  #
  # `gate` is deliberately NOT a compared component: census-compare.R drops
  # the key before comparing, so adding it here cannot change the verdict of
  # any case that does not use it. That matters -- putting `rel_gradient`
  # into the compared probe would have turned every existing IDENTICAL row
  # into a DIFFERS row, because the field does not exist before 1.2.11, and
  # would have silently invalidated results already measured.
  p <- probe_fit(fit)
  se <- fit$fit$se
  rc <- fit$fit$rcond
  rg <- fit$fit$rel_gradient
  p$gate <- list(
    converged = isTRUE(fit$fit$converged),
    # An NA standard error means the Hessian could not be inverted, so the
    # point is not a proper interior maximum whatever `converged` says.
    se_finite = !is.null(se) && is.numeric(se) && length(se) > 0L &&
      all(is.finite(se)),
    rcond = if (is.numeric(rc) && length(rc) == 1L) rc else NA_real_,
    # NULL where the version predates the gradient test; NA where the test
    # was not applied. Both are "cannot confirm", never "passed".
    rel_gradient = if (is.numeric(rg) && length(rg) == 1L) rg else NA_real_,
    rel_gradient_present = is.numeric(rg) && length(rg) == 1L
  )
  p
}

probe_num <- function(x) {
  list(kind = "num", value = x)
}

probe_chr <- function(x) {
  list(kind = "chr", value = as.character(x))
}

probe_boot <- function(b, fit) {
  # `replicates` is a long data frame: replicate, parameter, estimate. The
  # REPLICATES are compared, never `$summary` -- a summary statistic of 500
  # identical fits once reported n_success = 500 with no warning.
  #
  # #683 added a `restarted` column to `replicates`. It is moved out of the
  # compared frame into the `limit` key, so a bootstrap whose replicates are
  # bit-identical is not reported as DIFFERS for gaining a column; whether
  # it changed is reported in the limit section.
  #
  # A bootstrap with no successful replicate compares nothing. census-compare.R
  # refuses to call such a pair IDENTICAL (it reports NOT-COMPARED); this probe
  # records what it needs to decide that, and computes no spread over nothing.
  reps <- b$replicates
  restarted_col <- NULL
  if (is.data.frame(reps) && "restarted" %in% names(reps)) {
    restarted_col <- reps$restarted
    reps$restarted <- NULL
  }
  n_rows <- if (is.data.frame(reps)) nrow(reps) else 0L
  spread <- if (n_rows > 0L) {
    tapply(reps$estimate, reps$parameter, stats::sd)
  } else {
    NULL
  }
  # A FIXED parameter has zero spread by design (early.m in the #683 cases),
  # so the zero-spread check reads only the free ones. The fit's own mask
  # says which are fixed; without it every parameter counts as free.
  fm <- fit$fit$fixed_mask
  th <- names(fit$fit$theta)
  fixed <- if (length(fm) == length(th)) th[as.logical(fm)] else character()
  free_spread <- spread[setdiff(names(spread), fixed)]
  list(
    kind = "boot",
    replicates = reps,
    n_success = b$n_success,
    n_failed = b$n_failed,
    failure_reasons = b$failure_reasons,
    mode = b$mode,
    # A free parameter must actually vary across resamples. An exactly-zero
    # spread means the replicates are identical, which is the defect this
    # field exists to make visible rather than to hide.
    spread = spread,
    fixed_params = fixed,
    # Reported by census-compare.R on each side, even when both sides match:
    # identical replicates on both sides are still identical replicates.
    zero_spread_free = names(free_spread)[!is.na(free_spread) & free_spread == 0],
    spread_is_zero = length(free_spread) > 0L &&
      any(!is.na(free_spread) & free_spread == 0),
    limit = list(
      field = if ("n_restarted" %in% names(b)) "present" else "absent",
      n_restarted = if ("n_restarted" %in% names(b)) b$n_restarted else "absent",
      restart_reasons = if ("restart_reasons" %in% names(b)) {
        b$restart_reasons
      } else {
        "absent"
      },
      restarted_column = if (is.null(restarted_col)) {
        "absent"
      } else {
        paste0(sum(restarted_col), " of ", length(restarted_col),
               " rows TRUE")
      }
    )
  )
}

probe_stepwise <- function(sw, trace_txt) {
  # hzr_stepwise() returns the final FIT, so `sw$fit$fit$theta` is NULL and
  # an identical(NULL, NULL) control built on it is hollow. Read `sw$fit`.
  th <- sw$fit$theta
  stopifnot("stepwise probe read a NULL theta" = !is.null(th))
  list(
    kind = "stepwise",
    theta = th,
    coef_names = names(stats::coef(sw)),
    objective = sw$fit$objective,
    trace = as.character(trace_txt)
  )
}

# ---------------------------------------------------------------------------
# The battery
# ---------------------------------------------------------------------------
#
# `kp` marks a case that MUST differ between the two refs compared, and
# `kp_expect` names what the new side must show for it (see GROUP 3). The
# comparison script fails when a `kp` case compares identical or its
# expectation is not met: that is the census proving it can detect a change,
# rather than asserting it found none.
#
# The known positives are tied to a release range, because they are changes
# that range makes. The six from the 1.3.0 cycle (#253 x3, #300, the G3
# underflow and the no-theta contract) were v1.2.9 bugs, fixed before 1.3.0,
# so between 1.3.0 and any later ref they compare identical and every run
# failed at "Known positives" before reporting anything. Those cases stay as
# CONTROLS (a difference on one now is unexplained) and the known positives
# for v1.3.0 -> 1.3.1 are the #683 cases in GROUP 3. Replace them at the next
# cut. census-compare.R also plants a perturbation in every compared probe,
# independent of any range, and fails if the comparator misses one.

census_cases <- function() {
  list(

    ## -- single distributions, no covariates --------------------------------
    exp_basic = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "exponential", theta = c(log(0.3)), fit = TRUE))
    }),
    wei_basic = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "weibull", theta = c(0.3, 1), fit = TRUE))
    }),
    lnorm_basic = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "lognormal", theta = c(0, 1), fit = TRUE))
    }),
    llogis_basic = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "loglogistic", theta = c(0, 1), fit = TRUE))
    }),

    ## -- covariates, both interfaces ---------------------------------------
    wei_cov_formula = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(survival::Surv(time, status) ~ age + z, data = d$basic,
                       dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE))
    }),
    wei_cov_vector = list(needs = "hazard", fn = function(d) {
      xm <- cbind(age = d$basic$age, z = d$basic$z)
      probe_fit(hazard(time = d$basic$time, status = d$basic$status, x = xm,
                       dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE))
    }),
    wei_cov_factor = list(needs = "hazard", fn = function(d) {
      # z + a 3-level factor g -> 3 design columns, so 2 + 3 = 5 values.
      probe_fit(hazard(survival::Surv(time, status) ~ z + g, data = d$basic,
                       dist = "weibull", theta = c(0.3, 1, 0, 0, 0),
                       fit = TRUE))
    }),
    exp_cov_formula = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(survival::Surv(time, status) ~ z, data = d$basic,
                       dist = "exponential", theta = c(log(0.3), 0),
                       fit = TRUE))
    }),

    ## -- left truncation ---------------------------------------------------
    # CONTROLS since 1.3.0 (known positives for v1.2.9 -> 1.3.0, #253):
    # exponential, log-logistic and log-normal fits ignored `time_lower` as an
    # entry time before 1.3.0.
    exp_lefttrunc = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$trunc$time, status = d$trunc$status,
                       time_lower = d$trunc$entry,
                       x = cbind(z = d$trunc$z),
                       dist = "exponential", theta = c(log(0.25), 0.4),
                       fit = TRUE))
    }),
    lnorm_lefttrunc = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$trunc$time, status = d$trunc$status,
                       time_lower = d$trunc$entry,
                       x = cbind(z = d$trunc$z),
                       dist = "lognormal", theta = c(0, 0, 0), fit = TRUE))
    }),
    llogis_lefttrunc = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$trunc$time, status = d$trunc$status,
                       time_lower = d$trunc$entry,
                       x = cbind(z = d$trunc$z),
                       dist = "loglogistic", theta = c(0, 0, 0), fit = TRUE))
    }),
    # Weibull already honoured entry times, so this one is a CONTROL: it is
    # expected to be identical. A difference here is a category (c) finding.
    wei_lefttrunc = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$trunc$time, status = d$trunc$status,
                       time_lower = d$trunc$entry,
                       x = cbind(z = d$trunc$z),
                       dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE))
    }),

    ## -- interval censoring ------------------------------------------------
    wei_interval = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$interval$time, status = d$interval$status,
                       time_lower = d$interval$time_lower,
                       time_upper = d$interval$time_upper,
                       dist = "weibull", theta = c(0.3, 1), fit = TRUE))
    }),
    mp_interval = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$interval$time, status = d$interval$status,
                       time_lower = d$interval$time_lower,
                       time_upper = d$interval$time_upper,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),

    ## -- multiphase --------------------------------------------------------
    mp2_cdf = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),
    mp3_cdf = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 0.5,
                                                       nu = 1, m = 0),
                                     mid = hzr_phase("cdf", t_half = 3,
                                                     nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),
    mp_g3 = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     late = hzr_phase("hazard", t_half = 4,
                                                      nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),
    mp_cov_formula = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(survival::Surv(time, status) ~ age, data = d$basic,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),
    mp_lefttrunc = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$trunc$time, status = d$trunc$status,
                       time_lower = d$trunc$entry,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),
    mp_objective_sas = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$interval$time, status = d$interval$status,
                       time_lower = d$interval$time_lower,
                       time_upper = d$interval$time_upper,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       objective = "sas", fit = TRUE))
    }),

    ## -- weights -----------------------------------------------------------
    wei_weighted = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(survival::Surv(time, status) ~ z, data = d$weighted,
                       weights = d$weighted$w, dist = "weibull",
                       theta = c(0.3, 1, 0), fit = TRUE))
    }),
    mp_weighted = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(survival::Surv(time, status) ~ z, data = d$weighted,
                       weights = d$weighted$w, dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),

    ## -- conservation of events -------------------------------------------
    # CoE auto-disables outside status {0, 1} (coe_supported_data), so this
    # uses the 0/1 cohort deliberately.
    mp_coe = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       control = list(conservation_of_events = TRUE),
                       fit = TRUE))
    }),

    ## -- stepwise ----------------------------------------------------------
    stepwise_score = list(needs = c("hazard", "hzr_stepwise"), fn = function(d) {
      base <- hazard(survival::Surv(time, status) ~ 1, data = d$step,
                     dist = "weibull", theta = c(0.3, 1), fit = TRUE)
      trace_txt <- utils::capture.output(
        sw <- hzr_stepwise(base, scope = ~ x1 + x2 + x3 + x4, data = d$step,
                           direction = "forward", criterion = "score")
      )
      probe_stepwise(sw, trace_txt)
    }),
    stepwise_wald = list(needs = c("hazard", "hzr_stepwise"), fn = function(d) {
      base <- hazard(survival::Surv(time, status) ~ 1, data = d$step,
                     dist = "weibull", theta = c(0.3, 1), fit = TRUE)
      trace_txt <- utils::capture.output(
        sw <- hzr_stepwise(base, scope = ~ x1 + x2 + x3 + x4, data = d$step,
                           direction = "forward", criterion = "wald")
      )
      probe_stepwise(sw, trace_txt)
    }),
    stepwise_backward = list(needs = c("hazard", "hzr_stepwise"), fn = function(d) {
      base <- hazard(survival::Surv(time, status) ~ x1 + x2 + x3 + x4,
                     data = d$step, dist = "weibull",
                     theta = c(0.3, 1, 0, 0, 0, 0), fit = TRUE)
      trace_txt <- utils::capture.output(
        sw <- hzr_stepwise(base, scope = ~ x1 + x2 + x3 + x4, data = d$step,
                           direction = "backward", criterion = "wald")
      )
      probe_stepwise(sw, trace_txt)
    }),

    ## -- bootstrap ---------------------------------------------------------
    # 20 replicates at a fixed seed. Compared as REPLICATES, never as a
    # summary statistic of them: 500 identical fits once reported
    # n_success = 500 with sd exactly 0 and no warning.
    boot_wei = list(needs = c("hazard", "hzr_bootstrap"), fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ z, data = d$basic,
                    dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
      b <- hzr_bootstrap(fit, n_boot = 20L, seed = 20250917L)
      probe_boot(b, fit)
    }),
    boot_vector = list(needs = c("hazard", "hzr_bootstrap"), fn = function(d) {
      # No `x`: a design matrix passed directly cannot be resampled through
      # the stored call, and hzr_bootstrap() refuses it (see boot_vector_x).
      fit <- hazard(time = d$basic$time, status = d$basic$status,
                    dist = "weibull", theta = c(0.3, 1), fit = TRUE)
      probe_boot(hzr_bootstrap(fit, n_boot = 20L, seed = 20250917L), fit)
    }),
    boot_vector_x = list(needs = c("hazard", "hzr_bootstrap"),
                         fn = function(d) {
      fit <- hazard(time = d$basic$time, status = d$basic$status,
                    x = cbind(z = d$basic$z), dist = "weibull",
                    theta = c(0.3, 1, 0), fit = TRUE)
      probe_boot(hzr_bootstrap(fit, n_boot = 20L, seed = 20250917L), fit)
    }),

    ## -- predict(newdata = ), both interfaces ------------------------------
    predict_formula_newdata = list(needs = "hazard", fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ age + z, data = d$basic,
                    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE)
      nd <- data.frame(age = c(50, 65, 80), z = c(-1, 0, 1), time = c(1, 2, 3))
      probe_num(list(
        lp = stats::predict(fit, newdata = nd, type = "linear_predictor"),
        hz = stats::predict(fit, newdata = nd, type = "hazard"),
        sv = stats::predict(fit, newdata = nd, type = "survival"),
        ch = stats::predict(fit, newdata = nd, type = "cumulative_hazard")
      ))
    }),
    predict_vector_newdata = list(needs = "hazard", fn = function(d) {
      fit <- hazard(time = d$basic$time, status = d$basic$status,
                    x = cbind(age = d$basic$age, z = d$basic$z),
                    dist = "weibull", theta = c(0.3, 1, 0, 0), fit = TRUE)
      nd <- data.frame(age = c(50, 65, 80), z = c(-1, 0, 1), time = c(1, 2, 3))
      probe_num(list(
        lp = stats::predict(fit, newdata = nd, type = "linear_predictor"),
        hz = stats::predict(fit, newdata = nd, type = "hazard")
      ))
    }),
    predict_mp_newdata = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ age, data = d$basic,
                    dist = "multiphase",
                    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1,
                                                    m = 0),
                                  bg = hzr_phase("constant")),
                    fit = TRUE)
      nd <- data.frame(age = c(50, 65, 80), time = c(1, 2, 3))
      probe_num(list(
        hz = stats::predict(fit, newdata = nd, type = "hazard"),
        sv = stats::predict(fit, newdata = nd, type = "survival")
      ))
    }),
    # CONTROL since 1.3.0 (known positive for v1.2.9 -> 1.3.0, #300): an
    # exponential fit with no covariates used its
    # log rate as a coefficient for whatever column `newdata` carried, and
    # returned -280 for age = 70 where the answer is 0.
    predict_exp_nocov_newdata = list(needs = "hazard", fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ 1, data = d$basic,
                    dist = "exponential", theta = c(log(0.3)), fit = TRUE)
      nd <- data.frame(age = c(50, 70), time = c(1, 2))
      probe_num(stats::predict(fit, newdata = nd, type = "linear_predictor"))
    }),
    predict_se_survival = list(needs = "hazard", fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ z, data = d$basic,
                    dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
      nd <- data.frame(z = c(-1, 0, 1), time = c(1, 2, 3))
      probe_num(stats::predict(fit, newdata = nd, type = "survival",
                               se.fit = TRUE))
    }),

    ## -- exported pure functions -------------------------------------------
    # CONTROL since 1.3.0 (known positive for v1.2.9 -> 1.3.0): with a large
    # gamma, (t/tau)^gamma underflows and
    # hzr_decompos_g3() used to clamp to the smallest double, freezing G3 and
    # putting log(g3) wrong by more than 100.
    decompos_g3_underflow = list(needs = "hzr_decompos_g3", fn = function(d) {
      tt <- c(1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2)
      probe_num(hzr_decompos_g3(time = tt, tau = 10, gamma = 300,
                                alpha = 1, eta = 1))
    }),
    decompos_g3_ordinary = list(needs = "hzr_decompos_g3", fn = function(d) {
      tt <- seq(0.1, 20, by = 0.1)
      probe_num(hzr_decompos_g3(time = tt, tau = 5, gamma = 1.5,
                                alpha = 1, eta = 1))
    }),
    decompos_cdf = list(needs = "hzr_decompos", fn = function(d) {
      probe_num(hzr_decompos(time = seq(0.05, 20, by = 0.05),
                             t_half = 2, nu = 1.3, m = 0.4))
    }),
    phase_hazard_grid = list(needs = "hzr_phase_hazard", fn = function(d) {
      tt <- seq(0.05, 20, by = 0.05)
      probe_num(list(
        cdf = hzr_phase_hazard(tt, t_half = 2, nu = 1.3, m = 0.4, type = "cdf"),
        haz = hzr_phase_hazard(tt, t_half = 2, nu = 1.3, m = 0.4,
                               type = "hazard"),
        const = hzr_phase_hazard(tt, type = "constant")
      ))
    }),
    phase_cumhaz_grid = list(needs = "hzr_phase_cumhaz", fn = function(d) {
      tt <- seq(0.05, 20, by = 0.05)
      probe_num(list(
        cdf = hzr_phase_cumhaz(tt, t_half = 2, nu = 1.3, m = 0.4, type = "cdf"),
        haz = hzr_phase_cumhaz(tt, t_half = 2, nu = 1.3, m = 0.4,
                               type = "hazard"),
        const = hzr_phase_cumhaz(tt, type = "constant")
      ))
    }),
    math_primitives = list(needs = c("hzr_log1pexp", "hzr_log1mexp",
                                     "hzr_clamp_prob"), fn = function(d) {
      xs <- c(-800, -100, -1e-8, 0, 1e-8, 1, 40, 800)
      probe_num(list(
        l1pe = hzr_log1pexp(xs),
        l1me = hzr_log1mexp(-abs(xs) - 1e-12),
        clamp = hzr_clamp_prob(c(-1, 0, 1e-20, 0.5, 1 - 1e-20, 1, 2))
      ))
    }),
    theta_names = list(needs = c("hzr_theta_names", "hzr_phase"),
                       fn = function(d) {
      probe_chr(hzr_theta_names(
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        covariates = list(early = c("age", "z"), bg = character())))
    }),

    ## -- diagnostics -------------------------------------------------------
    gof_wei = list(needs = c("hazard", "hzr_gof"), fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ z, data = d$basic,
                    dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
      g <- hzr_gof(fit)
      # hzr_gof() returns a data frame (cum_observed, cum_expected, residual,
      # ...) with a `summary` attribute. This probe used to read $observed,
      # $expected, $ratio and $conservation_ratio, none of which exist, so it
      # compared four NULLs and reported IDENTICAL over nothing.
      tab <- as.list(as.data.frame(g))
      stopifnot("gof probe read no rows" = length(tab$cum_expected) > 0L)
      probe_num(list(table = tab, summary = attr(g, "summary")))
    }),
    kaplan = list(needs = c("hazard", "hzr_kaplan"), fn = function(d) {
      k <- hzr_kaplan(time = d$basic$time, status = d$basic$status)
      probe_num(list(time = k$time, surv = k$surv))
    }),
    nelson = list(needs = c("hazard", "hzr_nelson"), fn = function(d) {
      nl <- hzr_nelson(time = d$basic$time, event = d$basic$status)
      probe_num(list(time = nl$time, cumhaz = nl$cumhaz))
    }),
    deciles = list(needs = c("hazard", "hzr_deciles"), fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ z, data = d$basic,
                    dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
      probe_num(unclass(hzr_deciles(fit, time = 2,
                                    status = d$basic$status,
                                    event_time = d$basic$time)))
    }),
    calibrate = list(needs = c("hazard", "hzr_calibrate"), fn = function(d) {
      fit <- hazard(survival::Surv(time, status) ~ z, data = d$basic,
                    dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
      # hzr_calibrate() grades probabilities against events, so feed it the
      # model's own predicted event probability rather than the fit object.
      s <- stats::predict(fit, type = "survival")
      probe_num(unclass(hzr_calibrate(x = 1 - s, event = d$basic$status,
                                      groups = 5L)))
    }),

    ## -- translator --------------------------------------------------------
    # The translator's printed result for a small job: the call names and the
    # coverage it reports. The generated calls are NOT evaluated here (they
    # need the job's data bound); the corpus test in tests/ executes emitted
    # code. The case name is kept so it compares across versions.
    translate_sas_fit = list(needs = "hzr_translate_sas", fn = function(d) {
      src <- paste(
        "%hazard( proc hazard data = one;",
        "  time t;",
        "  event e;",
        "  parms muc = 0.1 mue = 0.2 nue = 1 taue = 1;",
        ");",
        sep = "\n"
      )
      # The translator reads a file, so the fixture is written to one. The
      # tempfile path is scrubbed from the compared output: it changes every
      # run and would report a difference on every case.
      f <- tempfile(fileext = ".sas")
      on.exit(unlink(f), add = TRUE)
      writeLines(src, f)
      scrub <- function(x) {
        gsub(basename(f), "<TMPFILE>", gsub(f, "<TMPFILE>", x, fixed = TRUE),
             fixed = TRUE)
      }
      # The ERROR path has to be scrubbed too, not just the printed output.
      # The translator names the file in its error message, so an un-scrubbed
      # failure records a different message on every run: two runs that behave
      # identically then compare as ERROR-BOTH-DIFFERENT purely because of a
      # random tempfile name, which fabricates a behavioural difference.
      #
      # Found by the Group 1 integrity check, which flagged this case as the
      # single row that moved between two runs of the SAME code -- probe
      # identical, warnings identical, message differing only in
      # "file842353d00f3c.sas" vs "file49174a18b25d.sas".
      # Warnings carry the path too ("14 untranslated construct(s) in
      # fileXXXX.sas"), and the runner captures warnings outside this function,
      # so they are muffled here and re-raised scrubbed.
      tr <- withCallingHandlers(
        tryCatch(hzr_translate_sas(f),
                 error = function(e) {
                   stop(scrub(conditionMessage(e)), call. = FALSE)
                 }),
        warning = function(w) {
          warning(scrub(conditionMessage(w)), call. = FALSE)
          invokeRestart("muffleWarning")
        }
      )
      probe_chr(scrub(utils::capture.output(print(tr))))
    }),
    argument_mapping = list(needs = "hzr_argument_mapping", fn = function(d) {
      probe_chr(utils::capture.output(print(hzr_argument_mapping())))
    }),

    ## -- deliberate contract changes, probed on purpose --------------------
    # CONTROL since 1.3.0 (known positive for v1.2.9 -> 1.3.0): `hazard(fit =
    # TRUE)` with no `theta` used to return an unfitted object with no error.
    # It has been an error since; expected ERROR-BOTH-SAME.
    fit_true_no_theta = list(needs = "hazard", fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "weibull", fit = TRUE))
    }),
    # The multiphase path is explicitly unaffected: it assembles its own start
    # from `phases`. This is the CONTROL for the case above -- a difference
    # here would contradict the NEWS bullet and is a category (c) finding.
    mp_no_theta = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(time = d$basic$time, status = d$basic$status,
                       dist = "multiphase",
                       phases = list(early = hzr_phase("cdf", t_half = 1,
                                                       nu = 1, m = 0),
                                     bg = hzr_phase("constant")),
                       fit = TRUE))
    }),

    ## -- SAS fixtures ------------------------------------------------------
    # Gated on HAZARD_EXAMPLES_DIR, the same variable R-CMD-check.yaml sets
    # for the one runner that has the `ehrlinger/hazard` checkout. Without it
    # this case records ERROR, with a message saying the fixture is
    # unavailable (census-run.R reserves ABSENT for a missing export), rather
    # than disappearing; read an ERROR-BOTH-SAME here as "not run". A
    # mounted /Volumes/qhsstudies is NOT sufficient on its own: the fixtures
    # the translator corpus reads live in the hazard repo, not on the volume.
    sas_outhaz_fixture = list(needs = "hzr_read_outhaz", fn = function(d) {
      dir <- Sys.getenv("HAZARD_EXAMPLES_DIR", unset = "")
      if (!nzchar(dir) || !dir.exists(dir)) {
        stop("HAZARD_EXAMPLES_DIR unset or absent: no SAS fixture to read. ",
             "This is a fixture-availability skip, not a behavioural result.",
             call. = FALSE)
      }
      f <- list.files(dir, pattern = "\\.lst$", recursive = TRUE,
                      full.names = TRUE)
      if (!length(f)) {
        stop("no .lst fixture under HAZARD_EXAMPLES_DIR = ", dir,
             call. = FALSE)
      }
      # Sorted so the same fixture is read under both versions.
      probe_chr(utils::capture.output(str(hzr_read_outhaz(sort(f)[[1]]))))
    }),

    ## =====================================================================
    ## GROUP 2, APPENDED 2026-09-17: gradient-gated multiphase cases.
    ##
    ## Every case ABOVE this line is GROUP 1 and is unchanged, byte for byte,
    ## from the runs already measured -- so their verdicts stay comparable
    ## across those runs and nothing already reported is re-interpreted.
    ##
    ## These cases differ from Group 1's multiphase rows in two ways:
    ##   - they run on `phased`, which has genuine early/background structure;
    ##   - they carry `gate = TRUE`, so census-compare.R REFUSES to classify
    ##     them as DIFFERS or IDENTICAL unless both sides reach a proper
    ##     optimum. A row that cannot be interpreted says so in the output
    ##     instead of appearing as a verdict.
    ##
    ## So a DIFFERS or IDENTICAL verdict on an `mpc_` row means something that
    ## the same verdict on a Group 1 multiphase row does not.
    ## =====================================================================

    mpc_2phase = list(needs = c("hazard", "hzr_phase"), gate = TRUE,
                      group = "appended", fn = function(d) {
      probe_fit_gated(hazard(
        time = d$phased$time, status = d$phased$status, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        fit = TRUE))
    }),
    mpc_2phase_cov = list(needs = c("hazard", "hzr_phase"), gate = TRUE,
                          group = "appended", fn = function(d) {
      probe_fit_gated(hazard(
        time = d$phased$time, status = d$phased$status,
        x = cbind(z = d$phased$z), dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        fit = TRUE))
    }),
    mpc_2phase_weighted = list(needs = c("hazard", "hzr_phase"), gate = TRUE,
                               group = "appended", fn = function(d) {
      probe_fit_gated(hazard(
        time = d$phased$time, status = d$phased$status,
        weights = d$phased$w, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        fit = TRUE))
    }),
    mpc_2phase_formula = list(needs = c("hazard", "hzr_phase"), gate = TRUE,
                              group = "appended", fn = function(d) {
      probe_fit_gated(hazard(
        survival::Surv(time, status) ~ z, data = d$phased, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        fit = TRUE))
    }),

    ## THE GATE'S OWN KNOWN POSITIVE.
    ##
    ## Three phases on data with only two timescales: not identifiable, and
    ## deliberately so. On main it reaches rcond 2.4e-21, an NA standard error
    ## and rel_gradient 2.34 against a required 6.06e-06. It MUST come out as
    ## UNINTERPRETABLE rather than as DIFFERS or IDENTICAL, and
    ## census-compare.R fails if it does not -- the same discipline as the
    ## planted behavioural changes. A gate with no case that trips it is an
    ## untested gate.
    mpc_3phase_unidentified = list(needs = c("hazard", "hzr_phase"),
                                   gate = TRUE, group = "appended",
                                   gate_kp = "3 phases, 2 timescales",
                                   fn = function(d) {
      probe_fit_gated(hazard(
        time = d$phased3$time, status = d$phased3$status, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 0.4, nu = 1, m = 0),
                      mid = hzr_phase("cdf", t_half = 4, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        fit = TRUE))
    }),

    ## =====================================================================
    ## GROUP 3, APPENDED for the 1.3.1 cut: the #683 / #729 limit cases.
    ##
    ## One case for each limit #683 handles, each built from the fixture the
    ## branch's own tests use, plus two controls. The known positives here are
    ## the changes v1.3.0 -> 1.3.1 makes ON PURPOSE; `kp_expect` names what the
    ## NEW side must show, so a known positive that merely differs for some
    ## other reason does not count as detected:
    ##   limit_outcome  new `fit$fit$limit$outcome`
    ##   warn_class     a warning of this class was raised on the new side
    ##   error_class    the new side errored with this class
    ##   trigger        a trigger names this `phase.parameter`
    ##   n_restarted    the new bootstrap reports at least this many restarts
    ##   sas_rel_tol    the new objective is within this of the SAS optimum
    ## =====================================================================

    # Checksums of every fixture column the cases below read. IDENTICAL here
    # is what lets a difference on them be read as code, not data.
    limit_fixture_data = list(needs = "hazard", fn = function(d) {
      probe_num(list(avc = colSums(d$avc), ub = colSums(d$ub),
                     ridge = colSums(d$ridge),
                     valves = colSums(d$valves[, c(
                       "int_dead", "dead", "age_cop", "nyha", "mv_nyha",
                       "double_", "ao_pinc", "black", "i_path")])))
    }),

    # t_half PAST THE DATA. The start puts the early rise a million times
    # beyond the last time; 1.3.0 returns wherever it stops, the branch
    # records phase_beyond_data and restarts from the grid.
    limit_past_data = list(
      needs = c("hazard", "hzr_phase"), kp = "#683 restart: t_half past the data",
      kp_expect = list(limit_outcome = "restarted",
                       warn_class = "hzr_limit_restarted",
                       trigger = "early.t_half"),
      fn = function(d) {
        probe_fit(hazard(
          time = d$ub$t, status = d$ub$s, dist = "multiphase",
          phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1,
                                          fixed = "m"),
                        late = hzr_phase("constant")),
          theta = c(log(0.1), log(max(d$ub$t) * 1e6), 1, 1, log(0.05)),
          fit = TRUE, control = list(n_starts = 1L)))
      }),

    # t_half BELOW THE FIRST OBSERVED TIME (#444): a "hazard" phase started a
    # thousand times below the data. The branch refuses it.
    limit_below_first = list(
      needs = c("hazard", "hzr_phase"),
      kp = "#683 refusal: t_half below the first time (#444)",
      kp_expect = list(error_class = "hzr_limit_refusal",
                       trigger = "early.t_half"),
      fn = function(d) {
        probe_fit(hazard(
          time = d$ub$t, status = d$ub$s, dist = "multiphase",
          phases = list(early = hzr_phase("hazard", t_half = 0.5, nu = 1, m = 0),
                        late = hzr_phase("constant")),
          theta = c(log(0.1), log(min(d$ub$t) / 1000), 1, 0, log(0.05)),
          fit = TRUE, control = list(n_starts = 1L)))
      }),

    # nu TO 0 AS A STEP (#448 fixture): the best objective is at the step, so
    # the branch refuses rather than return a clean start 130 units worse.
    limit_step_448 = list(
      needs = c("hazard", "hzr_phase"),
      kp = "#683 refusal: nu to 0, a step (#448)",
      kp_expect = list(error_class = "hzr_limit_refusal",
                       trigger = "early.nu"),
      fn = function(d) {
        probe_fit(hazard(
          time = d$avc$time, status = d$avc$status, dist = "multiphase",
          phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                          fixed = "m"),
                        late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1,
                                         eta = 1,
                                         constraint = "alpha_gamma_eta")),
          fit = TRUE, control = list(n_starts = 1L)))
      }),

    # #729: hm.death.patient from the blind start t_half = 10, nu = 0.2. SAS
    # reaches -1792.07; 1.3.0 drifts to about -2036.66 with converged = TRUE.
    # The SAS optimum is read from the listing, so the case needs the
    # ehrlinger/hazard fixtures (HAZARD_EXAMPLES_DIR) and records an ERROR on
    # both sides without them, which fails this known positive loudly.
    limit_729_blind = list(
      needs = c("hazard", "hzr_phase"), kp = "#683/#729 restart: blind start",
      kp_expect = list(limit_outcome = "restarted",
                       warn_class = "hzr_limit_restarted",
                       trigger = "early.t_half", sas_rel_tol = 1e-5),
      fn = function(d) {
        dir <- Sys.getenv("HAZARD_EXAMPLES_DIR", unset = "")
        lst <- file.path(dir, "hm.death.patient.lst")
        if (!nzchar(dir) || !file.exists(lst)) {
          stop("HAZARD_EXAMPLES_DIR unset or without hm.death.patient.lst: ",
               "the #729 case cannot be judged against SAS. This is a ",
               "fixture-availability skip, not a behavioural result.",
               call. = FALSE)
        }
        ll <- grep("Log likelihood =", readLines(lst), value = TRUE)
        sas_ll <- as.numeric(sub(".*=", "", ll[[length(ll)]]))
        stopifnot("no SAS log likelihood parsed" = is.finite(sas_ll))
        p <- probe_fit(hazard(
          survival::Surv(int_dead, dead) ~ 1, data = d$valves,
          dist = "multiphase",
          phases = list(
            early = hzr_phase("cdf", t_half = 10, nu = 0.2, m = 1, fixed = "m",
                              formula = ~ age_cop + nyha + mv_nyha + double_ +
                                ao_pinc),
            constant = hzr_phase("constant",
                                 formula = ~ age_cop + black + i_path)),
          fit = TRUE,
          control = list(n_starts = 1L, maxit = 2000, conserve = TRUE)))
        p$limit$sas_objective <- sas_ll
        p
      }),

    # CONTROL: a weak-only ridge. The branch tries the grid, finds nothing
    # better and keeps the original, so the numbers must be bit-identical to
    # 1.3.0 while the limit surface reads "kept".
    limit_weak_kept = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(
        time = d$ridge$time, status = d$ridge$status, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                      const = hzr_phase("constant")),
        fit = TRUE, control = list(n_starts = 5, maxit = 500)))
    }),

    # CONTROL: a clean multiphase fit; the branch records `$limit` NULL.
    limit_clean = list(needs = c("hazard", "hzr_phase"), fn = function(d) {
      probe_fit(hazard(
        time = d$avc$time, status = d$avc$status, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                                        fixed = "m"),
                      constant = hzr_phase("constant")),
        fit = TRUE, control = list(n_starts = 1L)))
    }),

    # The bootstrap of the past-the-data fit: on the branch every replicate
    # meets the limit, so some restart and the rest are refused.
    boot_limit_past_data = list(
      needs = c("hazard", "hzr_phase", "hzr_bootstrap"),
      kp = "#683 bootstrap: restarted replicates counted",
      kp_expect = list(n_restarted = 1L),
      fn = function(d) {
        fit <- hazard(
          time = d$ub$t, status = d$ub$s, dist = "multiphase",
          phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1,
                                          fixed = "m"),
                        late = hzr_phase("constant")),
          theta = c(log(0.1), log(max(d$ub$t) * 1e6), 1, 1, log(0.05)),
          fit = TRUE, control = list(n_starts = 1L))
        probe_boot(hzr_bootstrap(fit, n_boot = 8L, seed = 683L), fit)
      }),

    # CONTROL: the bootstrap of a well-identified multiphase fit. The
    # replicates must be bit-identical; the branch adds n_restarted = 0.
    boot_mp_phased = list(needs = c("hazard", "hzr_phase", "hzr_bootstrap"),
                          fn = function(d) {
      fit <- hazard(
        time = d$phased$time, status = d$phased$status, dist = "multiphase",
        phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
                      bg = hzr_phase("constant")),
        fit = TRUE)
      probe_boot(hzr_bootstrap(fit, n_boot = 10L, seed = 20261110L), fit)
    })
  )
}
