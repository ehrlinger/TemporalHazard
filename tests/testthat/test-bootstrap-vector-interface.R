# hzr_bootstrap() must resample a fit built with the vector interface.
#
# hazard() accepts either a formula plus `data`, or bare time/status vectors.
# Resampling `data` alone is enough for the formula interface; for the vector
# one the stored call holds `time = d$col` as an EXPRESSION, so each replicate
# re-evaluated it against the original data and returned the original fit --
# n_success = n_boot, n_failed = 0, no warning, and n_boot identical
# replicates. A bootstrap summary that looked complete and contained nothing.

fixture <- function() {
  e <- new.env()
  utils::data("cabgkul", package = "TemporalHazard", envir = e)
  e$cabgkul
}

phases_fixed <- function() {
  list(
    early = hzr_phase("cdf", t_half = 0.19, nu = 1.4, m = 1,
                      fixed = c("t_half", "nu", "m")),
    late  = hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1.7,
                      fixed = c("tau", "gamma", "alpha", "eta"))
  )
}

sd_log_mu <- function(fit, n_boot = 20, seed = 1) {
  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = n_boot, seed = seed))
  stats::sd(b$replicates$estimate[b$replicates$parameter == "early.log_mu"])
}

test_that("a vector-interface fit actually resamples", {
  d <- fixture()
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = as.integer(d$dead), data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  # The regression: this was exactly 0.
  expect_gt(sd_log_mu(fit), 0)
})

test_that("vector and formula interfaces bootstrap identically", {
  # The strong form, and it has to compare the REPLICATES, not a summary of
  # them: two different sets of replicate estimates can share an SD, so
  # comparing sd() alone would pass on outputs that differ.
  d <- fixture()
  fv <- suppressWarnings(hazard(
    time = d$int_dead, status = as.integer(d$dead), data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))
  ff <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  bv <- suppressWarnings(hzr_bootstrap(fv, n_boot = 20, seed = 1))
  bf <- suppressWarnings(hzr_bootstrap(ff, n_boot = 20, seed = 1))

  # Same model, same data, same seed: every replicate estimate must match.
  expect_equal(bv$replicates, bf$replicates, tolerance = 1e-8)
  expect_equal(bv$n_success, bf$n_success)
  # And they must not be trivially equal by both being constant.
  expect_gt(stats::sd(bv$replicates$estimate[
    bv$replicates$parameter == "early.log_mu"]), 0)
})

test_that("time_lower and time_upper are resampled alongside time", {
  # Interval-censored rows carry their bounds in separate vectors; resampling
  # `time` but not the bounds would silently pair row i's time with row j's
  # interval.
  d <- fixture()
  status <- as.integer(d$dead)
  idx <- which(status == 0L)[1:20]
  status[idx] <- 2L
  lower <- pmax(d$int_dead - 0.1, 0)

  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = status,
    time_lower = lower, time_upper = d$int_dead, data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  expect_gt(sd_log_mu(fit), 0)
})

test_that("a partially-stored vector-interface fit is refused, not half-rewired", {
  # Rewiring some vector arguments and not others is worse than rewiring none:
  # the rewired ones follow the resample while the rest evaluate against the
  # original data, pairing row i's time with row j's status. Silent corruption
  # producing plausible numbers.
  d <- fixture()
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = as.integer(d$dead), data = d,
    dist = "multiphase", phases = phases_fixed(),
    fit = TRUE, control = list(n_starts = 1, maxit = 200)))

  # Simulate an object fitted by an older version that did not store `status`.
  fit$data$status <- NULL
  expect_error(suppressWarnings(hzr_bootstrap(fit, n_boot = 5, seed = 1)),
               "status")
})

# A vector fit made WITHOUT `data =` stores its evaluated vectors but no data
# frame. hzr_bootstrap() sized the resample with nrow() of that NULL frame, so
# every such fit was refused with a message naming the vectors 'NA', 'NA'
# (#259, #312).

avc_fixture <- function() {
  e <- new.env()
  utils::data("avc", package = "TemporalHazard", envir = e)
  stats::na.omit(e$avc[, c("int_dead", "dead", "age")])
}

no_data_weibull <- function(d) {
  hazard(time = d$int_dead, status = d$dead, dist = "weibull",
         theta = c(0.1, 1), fit = TRUE)
}

test_that("a vector fit made without `data =` resamples like the formula fit (#259)", {
  d <- avc_fixture()
  vf <- no_data_weibull(d)
  ff <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
               theta = c(0.1, 1), fit = TRUE)
  # The case under test: the stored vectors are all there is to resample.
  expect_null(vf$data$frame)

  bv <- hzr_bootstrap(vf, n_boot = 10, seed = 1)
  bf <- hzr_bootstrap(ff, n_boot = 10, seed = 1)
  expect_equal(bv$n_success, 10L)
  # Nothing failed: the tally is present and empty, not NULL.
  expect_identical(bv$failure_reasons, stats::setNames(integer(0), character(0)))

  # The known positive the equality below needs: replicates vary on every
  # parameter, so matching replicates cannot be two constant sets.
  params <- unique(bv$replicates$parameter)
  expect_length(params, 2L)
  for (p in params) {
    expect_gt(stats::sd(bv$replicates$estimate[bv$replicates$parameter == p]),
              0, label = p)
  }
  # Same seed, same rows drawn: every replicate matches the formula path.
  expect_equal(bv$replicates, bf$replicates, tolerance = 1e-8)
  # And a different seed draws different rows.
  b2 <- hzr_bootstrap(vf, n_boot = 10, seed = 2)
  expect_false(isTRUE(all.equal(b2$replicates$estimate,
                                bv$replicates$estimate)))
})

test_that("a multiphase vector fit made without `data =` resamples every free parameter (#312)", {
  d <- avc_fixture()
  int_dead <- d$int_dead
  dead <- d$dead
  f <- suppressWarnings(hazard(
    time = int_dead, status = dead, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)))
  expect_null(f$data$frame)

  b <- suppressWarnings(hzr_bootstrap(f, n_boot = 5L, seed = 1L))
  expect_equal(b$n_success, 5L)
  free <- c("early.log_mu", "early.log_t_half", "early.nu", "constant.log_mu")
  expect_true(all(free %in% b$replicates$parameter))
  for (p in free) {
    expect_gt(stats::sd(b$replicates$estimate[b$replicates$parameter == p]),
              0, label = p)
  }
})

test_that("select mode on a vector fit made without `data =` is refused before seeding", {
  # The candidates would be read unresampled from the environment, so there
  # is nothing for the screen to resample them with.
  vf <- no_data_weibull(avc_fixture())
  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  msg <- tryCatch(hzr_bootstrap(vf, n_boot = 2, seed = 1, scope = ~ age),
                  error = conditionMessage)
  expect_match(msg, "^hzr_bootstrap\\(\\)")
  expect_match(msg, "`scope`", fixed = TRUE)
  expect_match(msg, "`data =`", fixed = TRUE)
  expect_no_match(msg, "NA", fixed = TRUE)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("a no-data vector fit missing a stored vector names it, not NA", {
  vf <- no_data_weibull(avc_fixture())
  one <- vf
  one$data$status <- NULL
  msg <- tryCatch(hzr_bootstrap(one, n_boot = 2, seed = 1),
                  error = conditionMessage)
  expect_match(msg, "`status`", fixed = TRUE)
  expect_no_match(msg, "`time`", fixed = TRUE)
  expect_no_match(msg, "NA", fixed = TRUE)
  expect_match(msg, "refit", fixed = TRUE)

  # With `time` gone there is no row count; the stored `status` must not be
  # reported missing along with it.
  no_time <- vf
  no_time$data$time <- NULL
  msg <- tryCatch(hzr_bootstrap(no_time, n_boot = 2, seed = 1),
                  error = conditionMessage)
  expect_match(msg, "`time`", fixed = TRUE)
  expect_no_match(msg, "`status`", fixed = TRUE)

  both <- vf
  both$data$time <- NULL
  both$data$status <- NULL
  msg <- tryCatch(hzr_bootstrap(both, n_boot = 2, seed = 1),
                  error = conditionMessage)
  expect_match(msg, "`time`, `status`", fixed = TRUE)
  expect_no_match(msg, "NA", fixed = TRUE)
})

test_that("a weighted vector fit made without `data =` resamples its weights with the rows", {
  # The oracle is a hand-drawn replicate, not another bootstrap: comparing
  # two bootstraps would share the resampling code under test.
  d <- avc_fixture()
  w <- seq(0.5, 1.5, length.out = nrow(d))
  fit <- hazard(time = d$int_dead, status = d$dead, weights = w,
                dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  expect_null(fit$data$frame)

  b <- hzr_bootstrap(fit, n_boot = 6, seed = 1)
  expect_equal(b$n_success, 6L)
  for (p in unique(b$replicates$parameter)) {
    expect_gt(stats::sd(b$replicates$estimate[b$replicates$parameter == p]),
              0, label = p)
  }

  # Replicate 1 is the first draw after set.seed(1).
  set.seed(1)
  idx <- sample.int(nrow(d), size = nrow(d), replace = TRUE)
  by_hand <- hazard(time = d$int_dead[idx], status = d$dead[idx],
                    weights = w[idx], dist = "weibull", theta = c(0.1, 1),
                    fit = TRUE)
  expect_equal(b$replicates$estimate[b$replicates$replicate == 1L],
               unname(by_hand$fit$theta), tolerance = 1e-8)
})

test_that("a `data =` that is not a data frame is refused, naming its class", {
  d <- avc_fixture()
  fit <- hazard(time = d$int_dead, status = d$dead, data = as.list(d),
                dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  msg <- tryCatch(hzr_bootstrap(fit, n_boot = 2, seed = 1),
                  error = conditionMessage)
  expect_match(msg, "must be a data frame", fixed = TRUE)
  expect_match(msg, "this fit's `data` is a list", fixed = TRUE)
  expect_no_match(msg, "NA", fixed = TRUE)
})

test_that("a data.frame subclass is not refused: is.data.frame, not an exact class test", {
  # A tibble is a data.frame subclass, and it bootstraps correctly, so the
  # refusal above must not catch one. `class(x) == "data.frame"`, or an
  # exact-class inherits() test, would. The class is stamped on here rather
  # than taken from tibble, which is not a dependency of this package.
  d <- avc_fixture()
  sub <- structure(d, class = c("tbl_df", "tbl", "data.frame"))
  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = sub,
                dist = "weibull", theta = c(0.1, 1, 0), fit = TRUE)
  expect_s3_class(fit$data$frame, "tbl_df")

  b <- hzr_bootstrap(fit, n_boot = 5, seed = 1)
  expect_equal(b$n_success, 5L)
  expect_gt(stats::sd(b$replicates$estimate[b$replicates$parameter == "age"]),
            0)
})

test_that("vectors that do not match the rows of `data =` are refused as such", {
  # The vectors are stored, so "not stored ... refit" would be false, and
  # refitting the same call would be refused again.
  d <- avc_fixture()
  fit <- hazard(time = d$int_dead, status = d$dead, data = d[1:100, ],
                dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  msg <- tryCatch(hzr_bootstrap(fit, n_boot = 2, seed = 1),
                  error = conditionMessage)
  n <- nrow(d)
  expect_match(msg, paste0("`time` (", n, "), `status` (", n, ")"),
               fixed = TRUE)
  expect_match(msg, "(100 rows)", fixed = TRUE)
  expect_no_match(msg, "not stored", fixed = TRUE)
  expect_no_match(msg, "NA", fixed = TRUE)
})

# A separate fixture, so the rows the earlier tests keep are not changed by
# na.omit() over an extra column.
avc_fixture_mal <- function() {
  e <- new.env()
  utils::data("avc", package = "TemporalHazard", envir = e)
  stats::na.omit(e$avc[, c("int_dead", "dead", "age", "mal")])
}

mp_phases <- function(early_formula = NULL) {
  list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                      formula = early_formula),
    constant = hzr_phase("constant")
  )
}

test_that("a saved fit whose phase formula was ignored is refused before seeding", {
  skip_on_cran() # multiphase fit
  # Saved before #299: built without `data`, so `~ mal` was never used. Its
  # stored call cannot be refit, and every replicate used to fail, returning
  # no replicates and no error.
  d <- avc_fixture_mal()
  hollow <- suppressWarnings(hzr_saved_before_299(
    mp_phases(~ mal), time = d$int_dead, status = d$dead,
    dist = "multiphase", fit = TRUE,
    control = list(n_starts = 1L, conserve = FALSE)))
  expect_false(any(grepl("mal", names(hollow$fit$theta))))

  set.seed(7)
  before <- get(".Random.seed", envir = globalenv())
  msg <- tryCatch(hzr_bootstrap(hollow, n_boot = 3L, seed = 1L),
                  error = conditionMessage)
  expect_match(msg, "^hzr_bootstrap\\(\\): cannot resample this fit: ")
  expect_match(msg, "phase 'early' has a formula, `~mal`, that the fit ignored",
               fixed = TRUE)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("stepping-only refusals do not reach the bootstrap: a 3-level factor still resamples", {
  skip_on_cran() # multiphase fit and replicates
  # hzr_stepwise() refuses to step a phase that inherits a factor with more
  # than two levels. A bootstrap refits the exact call, which is fine, so it
  # must take only the ignored-formula check, not the whole blocker.
  d <- avc_fixture_mal()
  d$ageg <- cut(d$age, 3, labels = c("lo", "mid", "hi"))
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ ageg, data = d, dist = "multiphase",
    phases = mp_phases(), fit = TRUE,
    control = list(n_starts = 1L, conserve = FALSE)))
  # The known positive: the whole blocker would refuse this fit.
  expect_match(.hzr_inherit_blocker(fit),
               "a term that expands to more than one column", fixed = TRUE)

  b <- suppressWarnings(hzr_bootstrap(fit, n_boot = 3L, seed = 1L))
  expect_equal(b$n_success, 3L)
  ageg <- grep("ageg", unique(b$replicates$parameter), value = TRUE)
  expect_gt(length(ageg), 0L)
  for (p in ageg) {
    expect_gt(stats::sd(b$replicates$estimate[b$replicates$parameter == p]),
              0, label = p)
  }
})

# Replace a fit's stored refit with one that cycles, per replicate, through a
# non-finite objective, a success and an error. Deterministic, and it reaches
# both failure branches without depending on which resamples happen to fail.
with_flaky_refit <- function(fit) {
  state <- new.env()
  state$calls <- 0L
  flaky <- function(...) {
    state$calls <- state$calls + 1L
    switch(state$calls %% 3L + 1L,
           stop("synthetic refit failure"),
           list(fit = list(objective = NaN, theta = fit$fit$theta)),
           fit)
  }
  env <- new.env(parent = fit$call_env %||% globalenv())
  assign("flaky_refit", flaky, envir = env)
  fit$call[[1L]] <- as.name("flaky_refit")
  fit$call_env <- env
  fit
}

test_that("failure_reasons tallies every failed replicate, by reason, without warning on partial failure", {
  vf <- with_flaky_refit(no_data_weibull(avc_fixture()))
  expect_no_warning(b <- hzr_bootstrap(vf, n_boot = 6L, seed = 1L))
  expect_equal(b$n_success, 2L)
  expect_equal(b$n_failed, 4L)
  expect_identical(
    b$failure_reasons[sort(names(b$failure_reasons))],
    c("non-finite objective (did not converge)" = 2L,
      "synthetic refit failure" = 2L)
  )
  expect_equal(sum(b$failure_reasons), b$n_failed)
})

test_that("a bootstrap whose every replicate fails warns, and still returns its reasons", {
  # A real refit failure: the captured `theta` no longer holds finite values,
  # so hazard() refuses every replicate's refit.
  d <- avc_fixture()
  th <- c(0.1, 1)
  vf <- hazard(time = d$int_dead, status = d$dead, dist = "weibull",
               theta = th, fit = TRUE)
  expect_true(exists("th", envir = vf$call_env, inherits = FALSE))
  assign("th", c(NA_real_, 1), envir = vf$call_env)

  expect_warning(
    b <- hzr_bootstrap(vf, n_boot = 3L, seed = 1L),
    "no replicate succeeded out of n_boot = 3. The most common failure (3 of 3)",
    fixed = TRUE
  )
  expect_equal(b$n_success, 0L)
  expect_equal(b$n_failed, 3L)
  expect_equal(sum(b$failure_reasons), b$n_failed)
  expect_length(b$failure_reasons, 1L)
  expect_match(names(b$failure_reasons), "theta", fixed = TRUE)
  expect_equal(nrow(b$replicates), 0L)
})

test_that("a refit error with an empty message is still counted", {
  # R cannot index a tally by the name "", so a bare stop() used to vanish
  # from failure_reasons and the tally no longer summed to n_failed.
  vf <- no_data_weibull(avc_fixture())
  env <- new.env(parent = vf$call_env %||% globalenv())
  assign("silent_refit", function(...) stop(""), envir = env)
  vf$call[[1L]] <- as.name("silent_refit")
  vf$call_env <- env

  b <- suppressWarnings(hzr_bootstrap(vf, n_boot = 3L, seed = 1L))
  expect_equal(b$n_failed, 3L)
  expect_equal(sum(b$failure_reasons), b$n_failed)
  expect_identical(b$failure_reasons,
                   c("error with an empty message" = 3L))
})

# Select mode (`scope =`) catches its replicate errors in a second handler. If
# that handler dropped the condition, a failed replicate would still count,
# but its reason would read "non-finite objective" for every failure (#333).
# The per-replicate base refit evaluates the stored call; hzr_stepwise()'s own
# refits and the up-front validation screen call hazard() directly, so
# replacing the call's function fails the replicates and nothing else.
test_that("select mode tallies each failed replicate under its own reason (#333)", {
  d <- avc_fixture()
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  state <- new.env()
  state$calls <- 0L
  flaky <- function(...) {
    state$calls <- state$calls + 1L
    switch(state$calls %% 3L + 1L,
           stop("synthetic base refit failure"),
           list(fit = list(objective = NaN, theta = base$fit$theta)),
           {
             mc <- match.call()
             mc[[1L]] <- quote(hazard)
             eval(mc, parent.frame())
           })
  }
  env <- new.env(parent = base$call_env %||% globalenv())
  assign("flaky_base", flaky, envir = env)
  base$call[[1L]] <- as.name("flaky_base")
  base$call_env <- env

  b <- suppressWarnings(hzr_bootstrap(base, n_boot = 6L, seed = 1L,
                                      scope = ~ age, criterion = "wald"))
  expect_equal(b$mode, "select")
  expect_equal(b$n_success, 2L)
  expect_equal(b$n_failed, 4L)
  expect_identical(
    b$failure_reasons[sort(names(b$failure_reasons))],
    c("base refit did not converge" = 2L,
      "synthetic base refit failure" = 2L)
  )
  expect_equal(sum(b$failure_reasons), b$n_failed)
})

test_that("a refit that returns a bare vector is a failed replicate, not a crash (#333)", {
  vf <- no_data_weibull(avc_fixture())
  env <- new.env(parent = vf$call_env %||% globalenv())
  assign("atomic_refit", function(...) c(1, 2), envir = env)
  vf$call[[1L]] <- as.name("atomic_refit")
  vf$call_env <- env

  b <- suppressWarnings(hzr_bootstrap(vf, n_boot = 3L, seed = 1L))
  expect_equal(b$n_success, 0L)
  expect_equal(b$n_failed, 3L)
  expect_identical(b$failure_reasons,
                   c("refit returned a numeric, not a fit object" = 3L))
})
