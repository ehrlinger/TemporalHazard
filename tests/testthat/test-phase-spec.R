# test-phase-spec.R — Tests for hzr_phase() and internal helpers
#
# Validates construction, validation, printing, parameter name generation,
# starting value extraction, and phase list validation.

# ============================================================================
# hzr_phase() construction
# ============================================================================

test_that("hzr_phase('cdf') creates valid object with defaults", {
  p <- hzr_phase("cdf")
  expect_s3_class(p, "hzr_phase")
  expect_equal(p$type, "cdf")
  expect_equal(p$t_half, 1)
  expect_equal(p$nu, 1)
  expect_equal(p$m, 0)
  expect_null(p$formula)
})

test_that("hzr_phase('hazard') stores custom starting values", {
  p <- hzr_phase("hazard", t_half = 5, nu = 0.8, m = -0.5)
  expect_equal(p$type, "hazard")
  expect_equal(p$t_half, 5)
  expect_equal(p$nu, 0.8)
  expect_equal(p$m, -0.5)
})

test_that("hzr_phase('constant') stores NA for shape parameters", {
  p <- hzr_phase("constant")
  expect_equal(p$type, "constant")
  expect_true(is.na(p$t_half))
  expect_true(is.na(p$nu))
  expect_true(is.na(p$m))
})

test_that("hzr_phase('constant') ignores supplied shape values", {
  p <- hzr_phase("constant", t_half = 99, nu = 99, m = 99)
  expect_true(is.na(p$t_half))
  expect_true(is.na(p$nu))
  expect_true(is.na(p$m))
})

test_that("hzr_phase stores one-sided formula", {
  p <- hzr_phase("cdf", formula = ~ age + nyha)
  expect_true(inherits(p$formula, "formula"))
})

test_that("type is matched by match.arg", {
  # "c" is ambiguous (matches both "cdf" and "constant"), so use longer prefixes
  expect_equal(hzr_phase("cd")$type, "cdf")
  expect_equal(hzr_phase("h")$type, "hazard")
  expect_equal(hzr_phase("con")$type, "constant")
})


# ============================================================================
# Input validation
# ============================================================================

test_that("hzr_phase rejects invalid type", {
  expect_error(hzr_phase("unknown"))
})

test_that("hzr_phase rejects non-positive t_half", {
  expect_error(hzr_phase("cdf", t_half = 0), "t_half must be a positive")
  expect_error(hzr_phase("cdf", t_half = -1), "t_half must be a positive")
})

test_that("hzr_phase rejects non-finite shape parameters", {
  expect_error(hzr_phase("cdf", nu = Inf), "nu must be a finite")
  expect_error(hzr_phase("cdf", m = NA), "m must be a finite")
  expect_error(hzr_phase("cdf", t_half = NaN), "t_half must be a positive")
})

test_that("hzr_phase rejects m < 0 AND nu < 0", {
  expect_error(hzr_phase("cdf", nu = -1, m = -1),
               "both m < 0 and nu < 0")
})

test_that("hzr_phase rejects two-sided formula", {
  expect_error(hzr_phase("cdf", formula = y ~ x),
               "one-sided")
})

test_that("hzr_phase rejects non-formula formula argument", {
  expect_error(hzr_phase("cdf", formula = "not a formula"),
               "formula must be a one-sided formula")
})


# ============================================================================
# is_hzr_phase()
# ============================================================================

test_that("is_hzr_phase returns TRUE for phases", {
  expect_true(is_hzr_phase(hzr_phase("cdf")))
  expect_true(is_hzr_phase(hzr_phase("constant")))
})

test_that("is_hzr_phase returns FALSE for non-phases", {
  expect_false(is_hzr_phase("cdf"))
  expect_false(is_hzr_phase(list(type = "cdf")))
  expect_false(is_hzr_phase(42))
  expect_false(is_hzr_phase(NULL))
})


# ============================================================================
# print.hzr_phase()
# ============================================================================

test_that("print.hzr_phase produces expected output for cdf", {
  p <- hzr_phase("cdf", t_half = 2.5, nu = 1.5, m = 0.3)
  out <- capture.output(print(p))
  expect_true(any(grepl("hzr_phase.*cdf.*early", out)))
  expect_true(any(grepl("t_half", out)))
})

test_that("print.hzr_phase produces expected output for constant", {
  p <- hzr_phase("constant")
  out <- capture.output(print(p))
  expect_true(any(grepl("constant.*flat", out)))
  # No shape line for constant
  expect_false(any(grepl("t_half", out)))
})

test_that("print.hzr_phase shows formula when present", {
  p <- hzr_phase("hazard", formula = ~ age + shock)
  out <- capture.output(print(p))
  expect_true(any(grepl("covariates.*age.*shock", out)))
})


# ============================================================================
# .hzr_phase_n_shape()
# ============================================================================

test_that(".hzr_phase_n_shape returns 3 for cdf and hazard", {
  expect_equal(.hzr_phase_n_shape(hzr_phase("cdf")), 3L)
  expect_equal(.hzr_phase_n_shape(hzr_phase("hazard")), 3L)
})

test_that(".hzr_phase_n_shape returns 0 for constant", {
  expect_equal(.hzr_phase_n_shape(hzr_phase("constant")), 0L)
})


# ============================================================================
# .hzr_phase_theta_names()
# ============================================================================

test_that("theta names for cdf phase without covariates", {
  p <- hzr_phase("cdf")
  nms <- .hzr_phase_theta_names(p, "early")
  expect_equal(nms, c("early.log_mu", "early.log_t_half",
                       "early.nu", "early.m"))
})

test_that("theta names for hazard phase with covariates", {
  p <- hzr_phase("hazard")
  nms <- .hzr_phase_theta_names(p, "late", c("age", "nyha"))
  expect_equal(nms, c("late.log_mu", "late.log_t_half",
                       "late.nu", "late.m",
                       "late.age", "late.nyha"))
})

test_that("theta names for constant phase", {
  p <- hzr_phase("constant")
  nms <- .hzr_phase_theta_names(p, "const")
  expect_equal(nms, "const.log_mu")
})

test_that("theta names for constant phase with covariates", {
  p <- hzr_phase("constant")
  nms <- .hzr_phase_theta_names(p, "const", c("age"))
  expect_equal(nms, c("const.log_mu", "const.age"))
})


# ============================================================================
# .hzr_phase_n_params()
# ============================================================================

test_that("n_params correct for cdf without covariates", {
  expect_equal(.hzr_phase_n_params(hzr_phase("cdf")), 4L)  # mu + 3 shape
})

test_that("n_params correct for cdf with 2 covariates", {
  expect_equal(.hzr_phase_n_params(hzr_phase("cdf"), n_covariates = 2), 6L)
})

test_that("n_params correct for constant without covariates", {
  expect_equal(.hzr_phase_n_params(hzr_phase("constant")), 1L)  # mu only
})

test_that("n_params correct for constant with 3 covariates", {
  expect_equal(.hzr_phase_n_params(hzr_phase("constant"), n_covariates = 3), 4L)
})


# ============================================================================
# .hzr_phase_start()
# ============================================================================

test_that("start values for cdf phase", {
  p <- hzr_phase("cdf", t_half = 2, nu = 1.5, m = 0.3)
  sv <- .hzr_phase_start(p)
  expect_equal(length(sv), 4L)
  expect_equal(sv[1], log(0.1))         # log_mu default
  expect_equal(sv[2], log(2))           # log_t_half
  expect_equal(sv[3], 1.5)             # nu
  expect_equal(sv[4], 0.3)             # m
})

test_that("start values for constant phase", {
  p <- hzr_phase("constant")
  sv <- .hzr_phase_start(p)
  expect_equal(length(sv), 1L)
  expect_equal(sv[1], log(0.1))
})

test_that("start values include zero betas for covariates", {
  p <- hzr_phase("hazard", t_half = 5, nu = 1, m = 0)
  sv <- .hzr_phase_start(p, n_covariates = 3)
  expect_equal(length(sv), 7L)  # log_mu + 3 shape + 3 betas
  expect_equal(sv[5:7], c(0, 0, 0))
})

test_that("start values respect mu_start", {
  p <- hzr_phase("cdf")
  sv <- .hzr_phase_start(p, mu_start = 0.5)
  expect_equal(sv[1], log(0.5))
})

test_that(".hzr_phase_start rejects non-positive mu_start", {
  expect_error(.hzr_phase_start(hzr_phase("cdf"), mu_start = 0),
               "mu_start must be a positive")
  expect_error(.hzr_phase_start(hzr_phase("cdf"), mu_start = -1),
               "mu_start must be a positive")
})


# ============================================================================
# .hzr_validate_phases()
# ============================================================================

test_that("validate_phases accepts named list of phases", {
  phases <- list(
    early    = hzr_phase("cdf"),
    constant = hzr_phase("constant"),
    late     = hzr_phase("hazard")
  )
  result <- .hzr_validate_phases(phases)
  expect_equal(names(result), c("early", "constant", "late"))
})

test_that("validate_phases auto-names unnamed phases", {
  phases <- list(hzr_phase("cdf"), hzr_phase("constant"))
  result <- .hzr_validate_phases(phases)
  expect_equal(names(result), c("phase_1", "phase_2"))
})

test_that("validate_phases auto-names partially named phases", {
  phases <- list(early = hzr_phase("cdf"), hzr_phase("constant"))
  result <- .hzr_validate_phases(phases)
  expect_equal(names(result), c("early", "phase_2"))
})

test_that("validate_phases rejects empty list", {
  expect_error(.hzr_validate_phases(list()),
               "non-empty list")
})

test_that("validate_phases rejects non-list", {
  expect_error(.hzr_validate_phases("not a list"),
               "non-empty list")
})

test_that("validate_phases rejects bare hzr_phase (not wrapped in list)", {
  expect_error(.hzr_validate_phases(hzr_phase("cdf")),
               "not a single hzr_phase")
})

test_that("validate_phases rejects non-phase elements", {
  expect_error(.hzr_validate_phases(list(hzr_phase("cdf"), "bad")),
               "not an hzr_phase")
})

test_that("validate_phases rejects duplicate names", {
  phases <- list(early = hzr_phase("cdf"), early = hzr_phase("hazard"))
  expect_error(.hzr_validate_phases(phases),
               "unique")
})

test_that("validate_phases rejects a phase named 'total'", {
  # 'total' is the key .hzr_multiphase_cumhaz() writes the accumulator under,
  # so a phase of that name loses its own contribution vector entirely.
  phases <- list(total = hzr_phase("cdf"), late = hzr_phase("constant"))
  expect_error(.hzr_validate_phases(phases), "reserved")
})

test_that("hazard() rejects a phase named 'total' before it can collide", {
  # The guard runs during validation, so the call fails before the optimizer is
  # reached even with fit = TRUE. Without it this returned a populated fit whose
  # 'total' phase had lost its contribution vector to the accumulator.
  dat <- data.frame(time = c(1, 2, 3, 4, 5), status = c(1, 1, 0, 1, 0))
  expect_error(
    hazard(survival::Surv(time, status) ~ 1, data = dat, dist = "multiphase",
           fit = TRUE,
           phases = list(total = hzr_phase("cdf"), late = hzr_phase("constant"))),
    "reserved")
})

test_that("hzr_theta_names() rejects a phase named 'total'", {
  # Exported, and validates through the same path, so the reservation reaches it.
  expect_error(hzr_theta_names(list(total = hzr_phase("cdf"))), "reserved")
})

test_that("'total' is still the name the cumhaz accumulator is stored under", {
  # Cross-check for the guard above. Asserting only that "reserved" is thrown
  # restates the guard's own string: rename the accumulator and that test stays
  # green while the guard protects a name that no longer collides. This reads
  # the key off the returned object instead, so it fails when the two drift.
  phases <- .hzr_validate_phases(
    list(early = hzr_phase("cdf"), late = hzr_phase("constant")))
  tt <- c(0.5, 1, 2, 3, 5)
  theta <- unlist(lapply(phases, .hzr_phase_start))
  cc <- stats::setNames(as.list(rep(0L, length(phases))), names(phases))
  xl <- stats::setNames(vector("list", length(phases)), names(phases))
  contrib <- .hzr_multiphase_cumhaz(tt, theta, phases, cc, xl, per_phase = TRUE)

  expect_equal(setdiff(names(contrib), names(phases)), "total")
  expect_length(contrib, length(phases) + 1L)
})


# ============================================================================
# Round-trip: start values match theta names in length
# ============================================================================

test_that("start values and theta names have matching length", {
  cases <- list(
    list(phase = hzr_phase("cdf"),      name = "early", ncov = 0L),
    list(phase = hzr_phase("cdf"),      name = "early", ncov = 3L),
    list(phase = hzr_phase("hazard"),   name = "late",  ncov = 2L),
    list(phase = hzr_phase("constant"), name = "const", ncov = 0L),
    list(phase = hzr_phase("constant"), name = "const", ncov = 1L)
  )

  for (cc in cases) {
    sv  <- .hzr_phase_start(cc$phase, n_covariates = cc$ncov)
    nms <- .hzr_phase_theta_names(cc$phase, cc$name,
                                   if (cc$ncov > 0) paste0("x", seq_len(cc$ncov))
                                   else character(0))
    np  <- .hzr_phase_n_params(cc$phase, n_covariates = cc$ncov)

    expect_equal(length(sv), length(nms),
                 label = paste("start/names length for", cc$name,
                               "ncov=", cc$ncov))
    expect_equal(length(sv), np,
                 label = paste("start/n_params for", cc$name,
                               "ncov=", cc$ncov))
  }
})


# ============================================================================
# Fixed parameter support
# ============================================================================

test_that("hzr_phase stores fixed parameter names", {
  ph <- hzr_phase("cdf", t_half = 1, nu = 2, m = 0, fixed = "shapes")
  expect_equal(sort(ph$fixed), c("m", "nu", "t_half"))

  ph2 <- hzr_phase("hazard", t_half = 5, nu = 1, m = 0, fixed = c("nu", "m"))
  expect_equal(sort(ph2$fixed), c("m", "nu"))

  ph3 <- hzr_phase("cdf", t_half = 1, nu = 1, m = 0)
  expect_equal(ph3$fixed, character(0))
})

test_that("hzr_phase rejects invalid fixed names", {
  expect_error(
    hzr_phase("cdf", t_half = 1, nu = 1, m = 0, fixed = "mu"),
    "Invalid fixed"
  )
  expect_error(
    hzr_phase("cdf", t_half = 1, nu = 1, m = 0, fixed = "log_mu"),
    "Invalid fixed"
  )
})

test_that("fixed is ignored for constant phases with warning", {
  expect_warning(
    ph <- hzr_phase("constant", fixed = "shapes"),
    "ignored for constant"
  )
  expect_equal(ph$fixed, character(0))
})

test_that(".hzr_phase_free_mask works correctly", {
  phases <- list(
    early    = hzr_phase("cdf", t_half = 1, nu = 1, m = 0,
                          fixed = "shapes"),
    constant = hzr_phase("constant"),
    late     = hzr_phase("hazard", t_half = 5, nu = 1, m = 0,
                          fixed = c("nu", "m"))
  )
  cov_counts <- c(early = 0L, constant = 0L, late = 2L)

  mask <- .hzr_phase_free_mask(phases, cov_counts)

  # early: log_mu=T, log_t_half=F, nu=F, m=F
  # constant: log_mu=T
  # late: log_mu=T, log_t_half=T, nu=F, m=F, beta1=T, beta2=T
  expected <- c(TRUE, FALSE, FALSE, FALSE,   # early
                TRUE,                          # constant
                TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)  # late
  expect_equal(mask, expected)
  expect_equal(sum(mask), 6)  # 6 free params
  expect_equal(sum(!mask), 5) # 5 fixed params
})

test_that(".hzr_phase_free_mask is all TRUE when nothing fixed", {
  phases <- list(
    early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0),
    const = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, const = 0L)
  mask <- .hzr_phase_free_mask(phases, cov_counts)
  expect_true(all(mask))
  expect_equal(length(mask), 5)  # 4 early + 1 const
})
