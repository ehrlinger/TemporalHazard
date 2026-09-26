# The score criterion reads each candidate from `data` by position and scores
# it against the fit's stored rows. Only the row COUNT was checked, so the
# same rows in another order scored every candidate against the wrong
# observations and entered a different variable, with no warning (#487).
#
# Each test asserts what the user reads: the variable entered, or the error.

ro_avc <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)
}

ro_weibull <- function(d) {
  suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age, data = d,
                          dist = "weibull", theta = c(0.5, 1, 0), fit = TRUE))
}

ro_multiphase <- function(d) {
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
}

ro_screen <- function(fit, data, scope, max_steps = 1L, slentry = 0.30) {
  suppressWarnings(hzr_stepwise(fit, scope = scope, data = data,
                                direction = "forward", criterion = "score",
                                slentry = slentry, max_steps = max_steps,
                                trace = FALSE))
}

ro_scope <- ~ com_iv + opmos + mal + nyha + inc_surg + orifice
ro_refusal <- "does not hold the rows the model was fitted on"

test_that("a frame sorted by time is refused, not screened (#487)", {
  d <- ro_avc()
  fit <- ro_weibull(d)
  srt <- d[order(d$int_dead), ]
  # Known positive: the same rows, in a different order.
  expect_identical(nrow(srt), nrow(d))
  expect_false(identical(srt$int_dead, d$int_dead))
  expect_equal(srt[order(as.integer(rownames(srt))), ], d,
               ignore_attr = TRUE)

  # The original order enters com_iv; before #487 the sorted frame entered
  # mal (Q = 1.92, p = 0.17) instead.
  expect_identical(ro_screen(fit, d, ro_scope)$steps$variable, "com_iv")
  expect_error(ro_screen(fit, srt, ro_scope), ro_refusal)
  # The message names the columns that moved, so the reader can see why.
  expect_error(ro_screen(fit, srt, ro_scope), "columns `study`, `status`")
})

test_that("a shuffled frame is refused on a multiphase constant phase (#487)", {
  skip_on_cran() # a multiphase fit plus a screen
  d <- ro_avc()
  base <- ro_multiphase(d)
  set.seed(487)
  perm <- sample(nrow(d))
  sh <- d[perm, ]
  expect_equal(sh[order(perm), ], d, ignore_attr = TRUE)
  expect_true(all(perm != seq_along(perm)))

  sc <- list(early = NULL,
             constant = ~ age + com_iv + opmos + mal + nyha + inc_surg)
  # Original order enters opmos; before #487 the shuffled frame entered
  # com_iv (Q = 4.04) with no warning.
  expect_identical(ro_screen(base, d, sc)$steps$variable, "opmos")
  expect_error(ro_screen(base, sh, sc), ro_refusal)
})

test_that("a column only `data` has is a candidate, not a mismatch (#487)", {
  # A candidate derived after the fit is absent from the fit's frame. It
  # cannot be compared, and must not be read as a reordering.
  d <- ro_avc()
  fit <- ro_weibull(d)
  d2 <- d
  d2$com_iv_x2 <- 2 * d$com_iv
  expect_false("com_iv_x2" %in% names(fit$data$frame))
  expect_identical(
    ro_screen(fit, d2, ~ com_iv_x2 + opmos + mal)$steps$variable,
    "com_iv_x2"
  )
})

test_that("an aligned frame still screens past the first step (#487)", {
  skip_on_cran() # a multiphase fit plus a three-step screen
  # Every accepted step refits; the next step then checks `data` against the
  # REFIT's stored frame, which must agree with it too.
  d <- ro_avc()
  base <- ro_multiphase(d)
  sc <- list(early = NULL,
             constant = ~ age + com_iv + opmos + mal + nyha + inc_surg)
  sw <- ro_screen(base, d, sc, max_steps = 3L, slentry = 0.99)
  expect_identical(sw$steps$variable, c("opmos", "inc_surg", "mal"))
})

test_that("rows hazard() dropped at time 0 do not read as a reordering (#487)", {
  # The fit's frame is the frame after the time-0 drop (#374, #476), so it is
  # shorter than the frame the caller holds; hzr_stepwise() trims that frame
  # first, and the check must then pass.
  d <- ro_avc()
  d$int_dead[c(3L, 40L)] <- 0
  fit <- ro_weibull(d)
  expect_identical(fit$data$dropped_time_zero, 2L)
  trimmed <- d[-c(3L, 40L), ]
  expect_identical(
    ro_screen(fit, d, ro_scope)$steps$variable,
    ro_screen(fit, trimmed, ro_scope)$steps$variable
  )
  expect_identical(ro_screen(fit, d, ro_scope)$steps$variable, "com_iv")
})

test_that("a lookup-only `data` of another length is not compared (#487)", {
  skip_on_cran() # a multiphase fit plus a screen
  # On the vector interface `data` may serve only to look names up, and is
  # then stored at its own length. It is not the fit's rows, so it says
  # nothing about the order of the frame the screen is given.
  set.seed(3)
  n <- 300
  D <- data.frame(tt = stats::rexp(n, 0.3) + 0.01,
                  ev = stats::rbinom(n, 1, 0.75),
                  x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  D$tt <- D$tt * exp(-0.8 * D$x1)
  fit <- suppressWarnings(hazard(
    time = D$tt, status = D$ev, data = D[1:5, ], dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                  const = hzr_phase("constant")),
    fit = TRUE
  ))
  expect_identical(nrow(fit$data$frame), 5L)
  sw <- ro_screen(fit, D, list(early = NULL, const = ~ x1 + x2))
  expect_identical(sw$steps$variable, "x1")
})

test_that("the bootstrap select mode screens each replicate unrefused (#487)", {
  skip_on_cran() # a bootstrap of stepwise screens
  # hzr_bootstrap() refits the base on each resample and screens that same
  # resample, so its frames agree by construction and must not be refused.
  d <- ro_avc()
  fit <- ro_weibull(d)
  bs <- suppressWarnings(hzr_bootstrap(fit, n_boot = 3, seed = 1,
                                       scope = ~ com_iv + opmos + mal,
                                       criterion = "score"))
  expect_identical(bs$n_failed, 0L)
  expect_identical(bs$n_success, 3L)
  expect_true("com_iv" %in% bs$replicates$parameter)
})
