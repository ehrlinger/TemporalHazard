# hzr_bootstrap(seed = ) and the caller's random number stream.
#
# A seeded bootstrap used to call set.seed() and leave the global generator
# there, so a script that seeded once at the top lost its stream at the first
# hzr_bootstrap(seed = ) call: every later draw followed from the bootstrap's
# seed instead. The seed now applies for the duration of the call only.

seed_fit <- function() {
  withr::local_seed(11)
  hazard(time = stats::rexp(60, 0.5), status = rep(1L, 60),
         theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
}

test_that("a seeded bootstrap restores the caller's stream", {
  fit <- seed_fit()

  set.seed(1)
  expected <- stats::runif(3)

  set.seed(1)
  hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_identical(stats::runif(3), expected)
})

test_that("a seeded bootstrap leaves an unseeded session unseeded", {
  fit <- seed_fit()
  withr::local_preserve_seed()

  suppressWarnings(rm(".Random.seed", envir = globalenv()))
  hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the same seed still gives the same replicates", {
  fit <- seed_fit()
  a <- hzr_bootstrap(fit, n_boot = 5, seed = 42)
  b <- hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_identical(a$replicates, b$replicates)
  # The replicates must differ from one another, or identical runs prove
  # nothing (AGENTS.md: assert a computation varied). Per parameter: pooled,
  # the sd of mu near 0.3 and nu near 1 is positive even for five identical
  # replicates (#611).
  sds <- tapply(a$replicates$estimate, a$replicates$parameter, stats::sd)
  expect_length(sds, 2L)
  expect_true(all(sds > 0))
})

test_that("seed = k gives the replicates set.seed(k) then seed = NULL gives (#611)", {
  # Pins "the replicates for a given seed are unchanged" to the plain
  # set.seed() semantics, and that seed = NULL draws from the caller's stream.
  fit <- seed_fit()
  withr::local_preserve_seed()
  set.seed(42)
  ref <- hzr_bootstrap(fit, n_boot = 5)
  got <- hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_identical(got$replicates, ref$replicates)
})

test_that("the caller's stream is restored when the bootstrap fails (#611)", {
  fit <- seed_fit()
  withr::local_preserve_seed()
  # Called after seeding; failing there leaves the seeded state behind unless
  # the restore runs on the error exit.
  local_mocked_bindings(.hzr_bootstrap_param_names = function(...) {
    stop("planted failure")
  })
  set.seed(1)
  expected <- stats::runif(3)
  set.seed(1)
  expect_error(hzr_bootstrap(fit, n_boot = 5, seed = 42), "planted failure")
  expect_identical(stats::runif(3), expected)
})

test_that("a seed that is not one whole number is refused, stream untouched (#611)", {
  fit <- seed_fit()
  withr::local_preserve_seed()
  for (bad in list(1.7, c(1, 2), TRUE, NA_real_, Inf, 2^31, "a")) {
    set.seed(1)
    expected <- stats::runif(3)
    set.seed(1)
    expect_error(hzr_bootstrap(fit, n_boot = 5, seed = bad),
                 "`seed` must be NULL or a single whole number",
                 info = deparse(bad))
    expect_identical(stats::runif(3), expected, info = deparse(bad))
  }
  # Negative and boundary whole numbers are valid seeds.
  expect_no_error(hzr_bootstrap(fit, n_boot = 2, seed = -1))
  expect_no_error(hzr_bootstrap(fit, n_boot = 2, seed = .Machine$integer.max))
})
