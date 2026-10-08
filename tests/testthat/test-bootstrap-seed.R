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
  # nothing (AGENTS.md: assert a computation varied).
  expect_gt(stats::sd(a$replicates$estimate), 0)
})
