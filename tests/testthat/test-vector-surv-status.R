# ---------------------------------------------------------------------------
# A Surv() object passed as `status` on the vector interface (#226)
#
# survival::Surv() and this package code censoring with different integers:
#
#   TemporalHazard   -1 left   0 right   1 event   2 interval
#   Surv "left"                0 left    1 event
#   Surv "interval"            0 right   1 event   2 left   3 interval
#
# The formula path always translated. The vector path unclassed the Surv and
# took its second column unchanged, so a left-censored row was fitted as
# right-censored, and for "interval" and "counting" the second column is not
# even the status (it is time2 / stop). No error, no warning. Both paths now
# read the Surv through one helper, so these tests assert the two agree.
# ---------------------------------------------------------------------------

surv_left_data <- function() {
  set.seed(226)
  n  <- 150
  # An early-hazard component, so the multiphase early phase has something
  # to estimate; without it that fit's Hessian is ill-conditioned.
  tt <- ifelse(runif(n) < 0.3, rexp(n, rate = 4), rexp(n, rate = 0.3))
  ev <- rep(1, n)
  lc <- sample.int(n, 35)
  ev[lc] <- 0                         # Surv "left": 0 = left-censored
  tt[lc] <- tt[lc] + 0.5              # observed only as "before this time"
  data.frame(tt = tt, ev = ev)
}

surv_interval_data <- function() {
  set.seed(2261)
  n  <- 150
  tt <- rexp(n, rate = 0.3)
  ev <- ifelse(tt >= 6, 0, 1)
  tt <- pmin(tt, 6)
  lo <- tt
  hi <- rep(NA_real_, n)
  idx_event <- which(ev == 1)
  iv <- idx_event[1:25]               # event known only within a window
  lf <- idx_event[26:40]              # event known only to precede a time
  lo[iv] <- pmax(tt[iv] - 0.5, 1e-3)
  hi[iv] <- tt[iv] + 0.5
  ev[iv] <- 3                         # Surv "interval" code
  lo[lf] <- tt[lf] + 0.5
  ev[lf] <- 2                         # Surv "left" code under "interval"
  data.frame(lo = lo, hi = hi, ev = ev)
}

surv_mp_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}

test_that("the #226 reproducer stores left-censored rows as -1", {
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  s  <- survival::Surv(c(1, 2, 3, 4, 5, 6), c(1, 0, 1, 1, 0, 1),
                       type = "left")
  f <- hazard(time = c(1, 2, 3, 4, 5, 6), status = s, dist = "multiphase",
              phases = ph, fit = FALSE)

  expect_equal(f$data$status, c(1, -1, 1, 1, -1, 1))
  expect_equal(f$data$time_upper, c(1, 2, 3, 4, 5, 6))
  expect_null(f$data$time_lower)
})

test_that("objective = 'sas' refuses a Surv left-censored row on both paths", {
  ph <- list(early = hzr_phase("cdf"), late = hzr_phase("constant"))
  d  <- data.frame(tt = c(1, 2, 3, 4, 5, 6), ev = c(1, 0, 1, 1, 0, 1))

  err_vector <- tryCatch(
    hazard(time = d$tt, status = survival::Surv(d$tt, d$ev, type = "left"),
           dist = "multiphase", phases = ph, fit = FALSE, objective = "sas"),
    error = conditionMessage
  )
  err_formula <- tryCatch(
    hazard(survival::Surv(tt, ev, type = "left") ~ 1, data = d,
           dist = "multiphase", phases = ph, fit = FALSE, objective = "sas"),
    error = conditionMessage
  )

  expect_match(err_vector, "does not support left-censored rows",
               fixed = TRUE)
  expect_identical(err_vector, err_formula)
})

test_that("Surv 'interval' codes 0/1/2/3 map to 0/1/-1/2 on the vector path", {
  lo <- c(1, 2, 3, 4)
  hi <- c(NA, NA, NA, 6)
  s  <- survival::Surv(lo, hi, c(0, 1, 2, 3), type = "interval")
  f  <- hazard(time = lo, status = s, theta = c(mu = 0.5, nu = 1),
               dist = "weibull", fit = FALSE)

  expect_equal(f$data$status,     c(0, 1, -1, 2))
  expect_equal(f$data$time_lower, c(0, 0, 0, 4))
  expect_equal(f$data$time_upper, c(1, 2, 3, 6))
})

test_that("Surv 'interval' rows contribute to the vector-path likelihood", {
  d  <- surv_interval_data()
  iv <- d$ev == 3
  fit_all <- hazard(time = d$lo,
                    status = survival::Surv(d$lo, d$hi, d$ev,
                                            type = "interval"),
                    theta = c(mu = 0.5, nu = 1), dist = "weibull", fit = TRUE)
  dk <- d[!iv, ]
  fit_drop <- hazard(time = dk$lo,
                     status = survival::Surv(dk$lo, dk$hi, dk$ev,
                                             type = "interval"),
                     theta = c(mu = 0.5, nu = 1), dist = "weibull", fit = TRUE)

  expect_equal(sum(fit_all$data$status == 2), 25L)
  expect_true(is.finite(fit_all$fit$objective))
  expect_true(is.finite(fit_drop$fit$objective))
  # Each interval row adds log(S(l) - S(u)) < 0, so dropping 25 of them has
  # to move the log-likelihood by a visible amount.
  expect_gt(abs(fit_all$fit$objective - fit_drop$fit$objective), 1)
})

# ---------------------------------------------------------------------------
# Interface parity on real fits
# ---------------------------------------------------------------------------

surv_parity_fits <- function(type, dist) {
  if (type == "left") {
    d <- surv_left_data()
    frm <- survival::Surv(tt, ev, type = "left") ~ 1
    s <- survival::Surv(d$tt, d$ev, type = "left")
    tm <- d$tt
  } else {
    d <- surv_interval_data()
    frm <- survival::Surv(lo, hi, ev, type = "interval") ~ 1
    s <- survival::Surv(d$lo, d$hi, d$ev, type = "interval")
    tm <- d$lo
  }
  args <- if (dist == "weibull") {
    list(dist = "weibull", theta = c(mu = 0.3, nu = 1))
  } else {
    list(dist = "multiphase", phases = surv_mp_phases(),
         control = list(n_starts = 1, maxit = 2000))
  }
  list(
    formula = do.call(hazard, c(list(frm, data = d, fit = TRUE), args)),
    vector  = do.call(hazard, c(list(time = tm, status = s, fit = TRUE), args))
  )
}

for (type in c("left", "interval")) {
  for (dist in c("weibull", "multiphase")) {
    test_that(sprintf("Surv '%s' gives the same %s fit on both interfaces",
                      type, dist), {
      fits <- surv_parity_fits(type, dist)
      ff <- fits$formula
      fv <- fits$vector

      # The data must carry the censoring it claims before a fit is compared:
      # a vector path that recoded nothing would still fit.
      expected_left <- if (type == "left") 35L else 15L
      expect_equal(sum(ff$data$status == -1), expected_left)

      expect_identical(fv$data$status, ff$data$status)
      expect_identical(fv$data$time, ff$data$time)
      expect_identical(fv$data$time_lower, ff$data$time_lower)
      expect_identical(fv$data$time_upper, ff$data$time_upper)

      expect_true(is.finite(ff$fit$objective))
      expect_true(is.finite(fv$fit$objective))
      expect_equal(fv$fit$objective, ff$fit$objective, tolerance = 1e-8)
      expect_equal(unname(coef(fv)), unname(coef(ff)), tolerance = 1e-6)
    })
  }
}

test_that("right and counting Surv give identical data on both interfaces", {
  d <- data.frame(start = c(0, 1, 0, 2), stop = c(1, 3, 2, 5),
                  ev = c(1, 0, 0, 1))

  ff <- hazard(survival::Surv(stop, ev) ~ 1, data = d,
               theta = c(mu = 0.5, nu = 1), dist = "weibull")
  fv <- hazard(time = d$stop, status = survival::Surv(d$stop, d$ev),
               theta = c(mu = 0.5, nu = 1), dist = "weibull")
  expect_equal(fv$data$status, c(1, 0, 0, 1))
  expect_identical(fv$data[c("time", "status", "time_lower", "time_upper")],
                   ff$data[c("time", "status", "time_lower", "time_upper")])

  ff <- hazard(survival::Surv(start, stop, ev) ~ 1, data = d,
               theta = c(mu = 0.5, nu = 1), dist = "weibull")
  fv <- hazard(time = d$stop,
               status = survival::Surv(d$start, d$stop, d$ev),
               theta = c(mu = 0.5, nu = 1), dist = "weibull")
  # Column 2 of a counting Surv is `stop`, not the status.
  expect_equal(fv$data$status, c(1, 0, 0, 1))
  expect_equal(fv$data$time_lower, c(0, 1, 0, 2))
  expect_identical(fv$data[c("time", "status", "time_lower", "time_upper")],
                   ff$data[c("time", "status", "time_lower", "time_upper")])
})

test_that("an interval2 Surv, left-censored row included, matches the formula", {
  # Surv() converts "interval2" to "interval". A left-censored row (lo = NA)
  # stores its upper bound in time1, so `time` is 5 there, not NA.
  d  <- data.frame(lo = c(NA, 2, 3, 1), hi = c(5, 2, NA, 4))
  s  <- survival::Surv(d$lo, d$hi, type = "interval2")
  ff <- hazard(survival::Surv(lo, hi, type = "interval2") ~ 1, data = d,
               theta = c(mu = 0.5, nu = 1), dist = "weibull")
  fv <- hazard(time = unclass(s)[, 1], status = s,
               theta = c(mu = 0.5, nu = 1), dist = "weibull")

  expect_equal(fv$data$status,     c(-1, 1, 0, 2))
  expect_equal(fv$data$time,       c(5, 2, 3, 1))
  expect_equal(fv$data$time_lower, c(0, 0, 0, 1))
  expect_equal(fv$data$time_upper, c(5, 2, 3, 4))
  expect_identical(fv$data[c("time", "status", "time_lower", "time_upper")],
                   ff$data[c("time", "status", "time_lower", "time_upper")])
})

test_that("a vector Surv fit bootstraps with its Surv-supplied bounds", {
  # The bounds come from the Surv, so they never appear in the stored call.
  # hzr_bootstrap() used to rewire only arguments named there, and refitted
  # every replicate with no bounds: 100% success over a different model, or
  # over nothing (sd exactly 0). Compare replicates, not a summary of them.
  di <- surv_interval_data()
  set.seed(2262)
  stop_t <- rexp(150, rate = 0.3) + 0.2
  dc <- data.frame(st = stop_t * runif(150, 0, 0.5), sp = stop_t,
                   e = rbinom(150, 1, 0.7))
  cases <- list(
    interval = list(
      d = di, frm = survival::Surv(lo, hi, ev, type = "interval") ~ 1,
      vec = function(d) {
        hazard(time = lo,
               status = survival::Surv(lo, hi, ev, type = "interval"),
               data = d, theta = c(mu = 0.3, nu = 1),
               dist = "weibull", fit = TRUE)
      }),
    counting = list(
      d = dc, frm = survival::Surv(st, sp, e) ~ 1,
      vec = function(d) {
        hazard(time = sp, status = survival::Surv(st, sp, e),
               data = d, theta = c(mu = 0.3, nu = 1),
               dist = "weibull", fit = TRUE)
      })
  )
  for (nm in names(cases)) {
    d  <- cases[[nm]]$d
    fv <- cases[[nm]]$vec(d)
    ff <- hazard(cases[[nm]]$frm, data = d, theta = c(mu = 0.3, nu = 1),
                 dist = "weibull", fit = TRUE)
    bv <- suppressWarnings(hzr_bootstrap(fv, n_boot = 10, seed = 1))
    bf <- suppressWarnings(hzr_bootstrap(ff, n_boot = 10, seed = 1))

    expect_equal(bv$n_success, 10L, label = nm)
    expect_equal(bv$replicates, bf$replicates, tolerance = 1e-8, label = nm)
    sds <- tapply(bv$replicates$estimate, bv$replicates$parameter, stats::sd)
    expect_true(all(sds > 0), label = nm)
  }
})

test_that("an unsupported Surv type is refused identically on both paths", {
  st <- factor(c("a", "b", "a"), levels = c("censor", "a", "b"))
  d  <- data.frame(tt = c(1, 2, 3))
  d$st <- st

  err_vector <- tryCatch(
    hazard(time = d$tt, status = survival::Surv(d$tt, st),
           theta = c(mu = 0.5, nu = 1), dist = "weibull"),
    error = conditionMessage
  )
  err_formula <- tryCatch(
    hazard(survival::Surv(tt, st) ~ 1, data = d,
           theta = c(mu = 0.5, nu = 1), dist = "weibull"),
    error = conditionMessage
  )
  expect_match(err_vector, "Unsupported Surv() type: mright", fixed = TRUE)
  expect_identical(err_vector, err_formula)
})

test_that("a Surv status whose times disagree with 'time' or a bound stops", {
  s <- survival::Surv(c(1, 2, 3), c(1, 0, 1), type = "left")
  expect_error(
    hazard(time = c(1, 2, 4), status = s, theta = c(mu = 0.5, nu = 1),
           dist = "weibull"),
    "'time' does not match the times in the Surv object", fixed = TRUE
  )
  expect_error(
    hazard(time = c(1, 2, 3), status = s, time_upper = c(1, 2, 9),
           theta = c(mu = 0.5, nu = 1), dist = "weibull"),
    "'time_upper' does not match the Surv object", fixed = TRUE
  )
  sc <- survival::Surv(c(0, 1), c(2, 3), c(1, 0))
  expect_error(
    hazard(time = c(2, 3), status = sc, time_lower = c(0, 0.5),
           theta = c(mu = 0.5, nu = 1), dist = "weibull"),
    "'time_lower' does not match the Surv object", fixed = TRUE
  )
  # A bound the Surv does not define is the caller's to supply: under
  # "left" that is `time_lower`, the counting-process entry time.
  f <- hazard(time = c(1, 2, 3), status = s, time_lower = c(0, 0.5, 0),
              theta = c(mu = 0.5, nu = 1), dist = "weibull")
  expect_equal(f$data$time_lower, c(0, 0.5, 0))
  expect_equal(f$data$status, c(1, -1, 1))
})
