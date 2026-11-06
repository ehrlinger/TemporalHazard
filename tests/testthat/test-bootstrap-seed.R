# hzr_bootstrap(seed = ) and the caller's random number stream.
#
# `seed = k` is set.seed(k) followed by the bootstrap's own draws, and it
# leaves the session's stream there (#653). #611 restored the caller's stream
# instead, which meant writing .Random.seed into .GlobalEnv, and CRAN policy
# forbids that. The documentation says so; these tests pin it.

seed_fit <- function() {
  withr::local_seed(11)
  hazard(time = stats::rexp(60, 0.5), status = rep(1L, 60),
         theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
}

test_that("a seeded bootstrap leaves the stream where set.seed(seed) would (#653)", {
  fit <- seed_fit()
  withr::local_preserve_seed()

  set.seed(1)
  untouched <- stats::runif(3)

  # The same call twice, from different caller streams, ends in the same
  # place: the stream follows from `seed`, not from the caller's.
  set.seed(1)
  hzr_bootstrap(fit, n_boot = 5, seed = 42)
  after_1 <- stats::runif(3)
  set.seed(2)
  hzr_bootstrap(fit, n_boot = 5, seed = 42)
  after_2 <- stats::runif(3)
  expect_identical(after_2, after_1)

  # And that place is set.seed(42) followed by the unseeded bootstrap's draws.
  set.seed(42)
  hzr_bootstrap(fit, n_boot = 5)
  expect_identical(stats::runif(3), after_1)

  # Known negative: the caller's own stream was not restored.
  expect_false(identical(after_1, untouched))
})

test_that("a seeded bootstrap does not write the stream back on exit (#653)", {
  fit <- seed_fit()
  withr::local_preserve_seed()

  set.seed(42)
  hzr_bootstrap(fit, n_boot = 5)
  expected <- get(".Random.seed", envir = globalenv())

  # An unseeded session is left seeded, at the bootstrap's own state, rather
  # than put back to having no stream.
  suppressWarnings(rm(".Random.seed", envir = globalenv()))
  hzr_bootstrap(fit, n_boot = 5, seed = 42)
  expect_identical(get(".Random.seed", envir = globalenv()), expected)
})

test_that("a call that stops before resampling leaves the stream alone (#653)", {
  fit <- seed_fit()
  withr::local_preserve_seed()
  # Fails after every refusal, just before the seed and the replicate loop.
  local_mocked_bindings(.hzr_bootstrap_param_names = function(...) {
    stop("planted failure")
  })
  set.seed(1)
  before <- get(".Random.seed", envir = globalenv())
  expect_error(hzr_bootstrap(fit, n_boot = 5, seed = 42), "planted failure")
  expect_identical(get(".Random.seed", envir = globalenv()), before)
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
