# hzr_theta_names(): the theta order, before any fit runs (#186).
#
# theta is positional and its entries are not on a common scale -- the late
# phase logs mu and tau but carries gamma, alpha and eta naturally. Wrapping
# the wrong element in log() produces a fit, not an error, so a shipped
# template needs to compare a hand-written theta0 against the real thing
# rather than against a comment describing it.

test_that("the returned order IS the order a fit uses", {
  skip_on_cran()
  # THE CONTRACT, and the only assertion here that can catch drift. A
  # hardcoded expected vector would encode the same assumption the function
  # encodes, and would pass with both of them wrong. Compare to what the
  # optimizer actually names its theta.
  specs <- list(
    list(early = hzr_phase("cdf"), late = hzr_phase("g3")),
    list(early = hzr_phase("cdf"), const = hzr_phase("constant")),
    list(a = hzr_phase("hazard"), b = hzr_phase("constant"),
         c = hzr_phase("cdf"))
  )
  set.seed(1)
  n <- 120
  tt <- stats::rexp(n, 0.3) + 0.05
  st <- rep(1L, n)

  for (ph in specs) {
    fit <- suppressWarnings(hazard(
      time = tt, status = st, dist = "multiphase", phases = ph,
      fit = TRUE, control = list(n_starts = 1)))
    # Guard the guard: an unfitted object has no theta names, and every
    # comparison below would then be NULL == NULL.
    expect_false(is.null(names(fit$fit$theta)))
    expect_identical(hzr_theta_names(ph), names(fit$fit$theta))
    expect_length(hzr_theta_names(ph), length(fit$fit$theta))
  }
})

test_that("the documented layout is what comes back", {
  # The order the help page prints. Stated explicitly so a change to it is a
  # visible test edit rather than a silent doc drift.
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("g3"))
  expect_identical(
    hzr_theta_names(ph),
    c("early.log_mu", "early.log_t_half", "early.nu", "early.m",
      "late.log_mu", "late.log_tau", "late.gamma", "late.alpha", "late.eta")
  )
})

test_that("covariates land in their phase's block, not at the end", {
  ph <- list(early = hzr_phase("cdf"), const = hzr_phase("constant"))
  got <- hzr_theta_names(ph, covariates = list(early = c("age", "sex")))
  expect_identical(
    got,
    c("early.log_mu", "early.log_t_half", "early.nu", "early.m",
      "early.age", "early.sex", "const.log_mu")
  )
  # A phase absent from `covariates` has none.
  expect_identical(hzr_theta_names(ph), c(
    "early.log_mu", "early.log_t_half", "early.nu", "early.m", "const.log_mu"))
})

test_that("unnamed phases get the same labels hazard() gives them", {
  skip_on_cran()
  # .hzr_validate_phases() auto-names phase_1, phase_2. Re-deriving those here
  # instead of reusing that validator would let the two spellings drift.
  ph <- list(hzr_phase("cdf"), hzr_phase("constant"))
  set.seed(2)
  n <- 100
  fit <- suppressWarnings(hazard(
    time = stats::rexp(n, 0.3) + 0.05, status = rep(1L, n),
    dist = "multiphase", phases = ph, fit = TRUE,
    control = list(n_starts = 1)))
  expect_identical(hzr_theta_names(ph), names(fit$fit$theta))
  expect_true(all(grepl("^phase_[12]\\.", hzr_theta_names(ph))))
})

test_that("a bad `phases` is refused the way hazard() refuses it", {
  expect_error(hzr_theta_names(hzr_phase("cdf")), "not a single hzr_phase")
  expect_error(hzr_theta_names(list()), "non-empty list")
  expect_error(hzr_theta_names(list(hzr_phase("cdf"), "nope")),
               "not an hzr_phase")
})

test_that("a covariates entry naming no such phase is an error, not ignored", {
  # Silently dropping it would return a vector short by exactly the covariates
  # the caller believed they had asked for -- and the caller is using this to
  # length-check a theta0, so a short answer is worse than no answer.
  ph <- list(early = hzr_phase("cdf"))
  expect_error(hzr_theta_names(ph, covariates = list(late = "age")),
               "no such phase")
  expect_error(hzr_theta_names(ph, covariates = list(early = 1L)),
               "character vector")
  expect_error(hzr_theta_names(ph, covariates = list("age")), "named list")
})

test_that("a covariate matrix with no colnames still yields conformable names", {
  # Latent defect the shared helper fixes. .hzr_phase_cov_names() falls back to
  # positional labels; passing character(0) instead produced FEWER names than
  # parameters, and `names(x) <- <short>` pads with NA rather than erroring, so
  # the misalignment was silent.
  ph <- list(early = hzr_phase("cdf"))
  x <- matrix(0, nrow = 3, ncol = 2)          # deliberately no dimnames
  expect_null(colnames(x))
  cov_names <- .hzr_phase_cov_names(ph, list(early = 2L), list(early = x))
  expect_identical(cov_names$early, c("x1", "x2"))
  nms <- .hzr_theta_names_list(ph, cov_names)
  # 1 log_mu + 3 shape + 2 covariates
  expect_length(nms, 6L)
  expect_false(anyNA(nms))
})
