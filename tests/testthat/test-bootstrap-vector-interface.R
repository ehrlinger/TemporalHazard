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
  expect_match(msg, "is a list", fixed = TRUE)
  expect_no_match(msg, "NA", fixed = TRUE)
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
