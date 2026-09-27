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

# A screen that counts the unverified-order warnings it raises and muffles
# the rest. ro_screen() suppresses every warning, so it cannot see them.
ro_screen_w <- function(fit, data, scope) {
  n <- 0L
  sw <- withCallingHandlers(
    hzr_stepwise(fit, scope = scope, data = data, direction = "forward",
                 criterion = "score", max_steps = 1L, trace = FALSE),
    hzr_score_rows_unverified = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning")
  )
  list(sw = sw, n_unverified = n)
}

test_that("a lookup-only `data` of another length is checked by time (#487)", {
  skip_on_cran() # a multiphase fit plus four screens
  # On the vector interface `data` may serve only to look names up, and is
  # then stored at its own length, so it is not the fit's rows and cannot be
  # compared with the screen's `data`. The fit's event times can.
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
  sc <- list(early = NULL, const = ~ x1 + x2)
  set.seed(4)
  sh <- D[sample(nrow(D)), ]
  expect_false(identical(sh$tt, D$tt))

  # `tt` holds the fit's times in order: checked, and screened silently.
  r <- ro_screen_w(fit, D, sc)
  expect_identical(r$sw$steps$variable, "x1")
  expect_identical(r$n_unverified, 0L)
  # The same rows shuffled: `tt` holds the times out of order, so refused.
  expect_error(ro_screen(fit, sh, sc), "column `tt` holds the fit's event")

  # Without `tt` nothing can be checked, and the reader is told, once.
  r <- ro_screen_w(fit, D[, -1L], sc)
  expect_identical(r$sw$steps$variable, "x1")
  expect_identical(r$n_unverified, 1L)
  expect_identical(ro_screen_w(fit, sh[, -1L], sc)$n_unverified, 1L)
  # The warning says why; the screen's other warnings are muffled outside it.
  suppressWarnings(expect_warning(
    hzr_stepwise(fit, scope = sc, data = sh[, -1L], direction = "forward",
                 criterion = "score", max_steps = 1L, trace = FALSE),
    "has 5 rows, not the fit's 300",
    class = "hzr_score_rows_unverified"
  ))
})

test_that("a vector fit without `data =` is not screened silently (#487)", {
  skip_on_cran() # a multiphase fit plus three screens
  # The reviewer's repro: no stored frame at all, so a shuffled avc entered
  # com_iv (Q = 4.04, p = 0.044) where the original enters opmos (Q = 8.13),
  # with no error or warning.
  d <- ro_avc()
  fit <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  expect_null(fit$data$frame)
  sc <- list(early = NULL,
             constant = ~ age + com_iv + opmos + mal + nyha + inc_surg)
  set.seed(487)
  sh <- d[sample(nrow(d)), ]
  expect_false(identical(sh$int_dead, d$int_dead))

  r <- ro_screen_w(fit, d, sc)
  expect_identical(r$sw$steps$variable, "opmos")
  expect_identical(r$n_unverified, 0L)
  expect_error(ro_screen(fit, sh, sc), "column `int_dead` holds the fit's")
  # The Wald criterion refits on the vector interface, whose response is the
  # fit's stored vectors, so it misread the shuffle the same way (com_iv for
  # opmos); it is checked too.
  expect_error(
    suppressWarnings(hzr_stepwise(fit, scope = sc, data = sh,
                                  direction = "forward", criterion = "wald",
                                  max_steps = 1L, trace = FALSE)),
    "column `int_dead` holds the fit's"
  )

  no_time <- setdiff(names(d), "int_dead")
  expect_identical(ro_screen_w(fit, sh[, no_time], sc)$n_unverified, 1L)
})

test_that("a vector fit WITH `data =` is checked under Wald and AIC (#487)", {
  skip_on_cran() # two multiphase fits plus eight screens
  # The fit stores `data` as its frame, at the fit's length, so the time
  # check defers to the frame comparison. The score test ran that; the Wald
  # and AIC screens never reach it, and their vector-interface refits pair
  # the stored response with `data` read by position: a shuffle entered
  # com_iv (Wald p = 0.098) where the original enters opmos (p = 0.021).
  d <- ro_avc()
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                         fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  fv <- suppressWarnings(hazard(time = d$int_dead, status = d$dead, data = d,
                                dist = "multiphase", phases = phases,
                                fit = TRUE))
  expect_identical(nrow(fv$data$frame), nrow(d))
  expect_null(fv$call$formula)
  sc <- list(early = NULL, constant = ~ age + com_iv + opmos + mal + inc_surg)
  set.seed(487)
  sh <- d[sample(nrow(d)), ]
  expect_false(identical(sh$int_dead, d$int_dead))
  screen <- function(fit, data, cr) {
    suppressWarnings(hzr_stepwise(fit, scope = sc, data = data,
                                  direction = "forward", criterion = cr,
                                  max_steps = 1L, trace = FALSE))
  }
  for (cr in c("wald", "aic")) {
    expect_identical(screen(fv, d, cr)$steps$variable, "opmos")
    expect_error(screen(fv, sh, cr), ro_refusal)
  }

  # An unweighted formula fit refits from `data` alone, response and all, so
  # a shuffle is self-consistent there and still enters opmos, unrefused.
  ff <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                                dist = "multiphase", phases = phases,
                                fit = TRUE))
  expect_null(ff$data$weights)
  for (cr in c("wald", "aic")) {
    expect_identical(screen(ff, sh, cr)$steps$variable, "opmos")
  }
})

test_that("a weighted formula fit is checked under Wald and AIC (#487)", {
  skip_on_cran() # a multiphase fit plus five screens
  # A formula refit rebuilds its response from `data`, but takes the base
  # fit's stored `weights`, in the fit's row order. A shuffle then weighted
  # the wrong rows: Wald entered opmos (p = 0.0019) where the original
  # enters inc_surg (p = 0.0025), and AIC entered nothing, unrefused.
  d <- ro_avc()
  w <- ifelse(d$opmos > stats::median(d$opmos), 3, 1)
  fw <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, weights = w,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  expect_false(is.null(fw$call$formula))
  expect_identical(length(fw$data$weights), nrow(d))
  sc <- list(early = NULL, constant = ~ age + com_iv + opmos + mal + inc_surg)
  set.seed(487)
  sh <- d[sample(nrow(d)), ]
  expect_false(identical(sh$int_dead, d$int_dead))
  screen <- function(data, cr) {
    suppressWarnings(hzr_stepwise(fw, scope = sc, data = data,
                                  direction = "forward", criterion = cr,
                                  max_steps = 1L, trace = FALSE))
  }
  expect_identical(screen(d, "wald")$steps$variable, "inc_surg")
  for (cr in c("wald", "aic")) {
    expect_error(screen(sh, cr), ro_refusal)
  }
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
  # com_iv enters every replicate, each at its own estimate: three screens
  # on three resamples, not one result repeated.
  hit <- bs$replicates$parameter == "com_iv"
  expect_equal(bs$replicates$replicate[hit], 1:3)
  expect_identical(length(unique(bs$replicates$estimate[hit])), 3L)
})
