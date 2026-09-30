# With no stored frame to compare, hzr_stepwise() checks the row order of
# `data` by finding a column equal to the fit's event times. With ties in those
# times, rows reordered WITHIN a tie leave that column unchanged, so #487 could
# only warn, and warned identically whether or not the rows had moved. The
# order is now checked against every per-row input the fit stores: status,
# interval bounds, weights and covariate design columns (#515).
#
# The check does not need the exact fitted order. If every per-row input is
# found in `data` in order, any remaining reordering swaps rows identical in
# everything the likelihood reads, and changes nothing.

rt_avc <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc)
}

rt_phases <- list(
  early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                       fixed = "shapes"),
  constant = hzr_phase("constant")
)

rt_scope <- list(early = NULL,
                 constant = ~ age + com_iv + opmos + mal + nyha + inc_surg)

rt_refusal <- "does not hold the rows the model was fitted on"

# Counts the unverified-order warnings a one-step screen raises, keeps their
# text, and muffles the rest.
rt_screen <- function(fit, data, criterion) {
  n <- 0L
  msg <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(fit, scope = rt_scope, data = data, direction = "forward",
                 criterion = criterion, max_steps = 1L, trace = FALSE),
    hzr_score_rows_unverified = function(w) {
      n <<- n + 1L
      msg <<- c(msg, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning")
  )
  list(sw = sw, n_unverified = n, msg = msg)
}

# Reorders rows at random within each group of identical `key` rows.
rt_shuffle_within <- function(key, seed) {
  set.seed(seed)
  g <- interaction(key, drop = TRUE)
  perm <- seq_len(nrow(key))
  for (lev in levels(g)) {
    p <- which(g == lev)
    if (length(p) > 1L) perm[p] <- p[sample.int(length(p))]
  }
  perm
}

test_that("a stored covariate settles ties, and a moved one is refused (#515)", {
  # A vector fit with a covariate matrix and no stored frame. The check is
  # called directly: it runs once per screen, before any refit.
  # `z` is age made distinct on every row, so with the times it tells every
  # row apart.
  d <- rt_avc()
  d$tt <- ceiling(d$int_dead / 12)
  d$z <- d$age + seq_len(nrow(d)) / 1000
  fit <- suppressWarnings(hazard(time = d$tt, status = d$dead,
                                 x = as.matrix(d["z"]), dist = "weibull",
                                 theta = c(0.5, 1, 0), fit = TRUE))
  expect_null(fit$data$frame)
  d <- d[, setdiff(names(d), "int_dead")]
  chk <- function(data) .hzr_check_data_row_order(fit, data, score = TRUE)

  # Known positive: in order, with the time, status and covariate.
  expect_no_condition(chk(d))
  # Without the status, time and z together still tell every row apart.
  expect_gt(anyDuplicated(d[c("tt", "dead")]), 0L)
  expect_identical(anyDuplicated(d[c("tt", "z")]), 0L)
  expect_no_condition(chk(d[names(d) != "dead"]))
  # Without the covariate, rows tied in time and status still differ in z.
  expect_warning(chk(d[names(d) != "z"]), "covariate `z`",
                 class = "hzr_score_rows_unverified")

  # Reordered within (time, status): the covariate moved, and is named.
  perm <- rt_shuffle_within(d[c("tt", "dead")], seed = 515)
  expect_identical(d$dead[perm], d$dead)
  expect_false(identical(d$z[perm], d$z))
  expect_error(chk(d[perm, ]), "column `z` holds the fit's covariate `z`")

  # A multiphase phase design is checked the same way.
  mp <- fit
  mp$data$x <- NULL
  mp$fit$x_list <- list(constant = as.matrix(d["z"]))
  expect_error(.hzr_check_data_row_order(mp, d[perm, ], score = TRUE),
               "phase `constant` covariate `z`")
  expect_no_condition(.hzr_check_data_row_order(mp, d, score = TRUE))
})

test_that("tied times with status in order are verified, not warned (#515)", {
  skip_on_cran() # a multiphase fit plus six screens
  # The issue's repro: fitted on rows sorted by time, then screened on rows
  # sorted by time and age. 27 rows move, all within ties, and 2 of them pair
  # a time with the wrong status.
  d <- rt_avc()
  ds <- d[order(d$int_dead), ]
  d2 <- d[order(d$int_dead, d$age), ]
  fit <- suppressWarnings(hazard(time = ds$int_dead, status = ds$dead,
                                 dist = "multiphase", phases = rt_phases,
                                 fit = TRUE))
  expect_null(fit$data$frame)
  expect_lt(length(unique(fit$data$time)), nrow(ds))
  expect_identical(d2$int_dead, ds$int_dead)
  expect_identical(sum(d2$dead != ds$dead), 2L)

  for (cr in c("score", "wald", "aic")) {
    # Known positive: the rows the fit was made on, in its order.
    r <- rt_screen(fit, ds, cr)
    expect_identical(r$n_unverified, 0L)
    expect_identical(r$sw$steps$variable, "opmos")
    # The same rows reordered within ties: refused, naming the status column.
    expect_error(rt_screen(fit, d2, cr), rt_refusal)
    expect_error(rt_screen(fit, d2, cr), "column `dead` holds the fit's status")
  }
})

test_that("an input `data` does not hold leaves the order unverified (#515)", {
  skip_on_cran() # a multiphase fit plus five screens
  # Discrete times (14 distinct in 305 rows) and weights passed as a vector.
  d <- rt_avc()
  d$tt <- ceiling(d$int_dead / 12)
  w <- ifelse(d$opmos > stats::median(d$opmos), 3, 1)
  fit <- suppressWarnings(hazard(time = d$tt, status = d$dead, weights = w,
                                 dist = "multiphase", phases = rt_phases,
                                 fit = TRUE))
  expect_null(fit$data$frame)
  expect_identical(length(unique(fit$data$time)), 14L)
  d <- d[, setdiff(names(d), "int_dead")]

  # `data` holds the times and the status in order, but not the weights, and
  # rows tied in time and status still differ in weight: unverified, and the
  # warning says which input was missing.
  r <- rt_screen(fit, d, "score")
  expect_identical(r$n_unverified, 1L)
  expect_match(r$msg, "weights", fixed = TRUE)

  # With the weights as a column, the order is verified.
  dw <- d
  dw$w <- w
  r_w <- rt_screen(fit, dw, "score")
  expect_identical(r_w$n_unverified, 0L)
  expect_identical(r_w$sw$steps$variable, r$sw$steps$variable)

  # A shuffle within tied times moves status with it: refused, with or without
  # the weights column.
  perm <- rt_shuffle_within(d["tt"], seed = 487)
  expect_gt(sum(perm != seq_along(perm)), 250L)
  expect_identical(d$tt[perm], d$tt)
  expect_false(identical(d$dead[perm], d$dead))
  expect_error(rt_screen(fit, d[perm, ], "score"), rt_refusal)
  expect_error(rt_screen(fit, dw[perm, ], "aic"), rt_refusal)
})

test_that("rows identical in every fit input may swap, and change nothing (#515)", {
  skip_on_cran() # a multiphase fit plus four screens
  # The check accepts a `data` whose rows are the fit's rows up to swaps
  # among rows identical in time, status and weight. The property that makes
  # that safe: the screen's answer is then exactly the aligned answer.
  d <- rt_avc()
  d$tt <- ceiling(d$int_dead / 12)
  d$w <- ifelse(d$opmos > stats::median(d$opmos), 3, 1)
  fit <- suppressWarnings(hazard(time = d$tt, status = d$dead, weights = d$w,
                                 dist = "multiphase", phases = rt_phases,
                                 fit = TRUE))
  d <- d[, setdiff(names(d), "int_dead")]
  perm <- rt_shuffle_within(d[c("tt", "dead", "w")], seed = 515)
  sw_d <- d[perm, ]
  # Known positive: many rows moved, and candidates moved with them, while
  # every fit input stayed in place.
  expect_gt(sum(perm != seq_along(perm)), 200L)
  expect_false(identical(sw_d$opmos, d$opmos))
  expect_identical(sw_d[c("tt", "dead", "w")], d[c("tt", "dead", "w")],
                   ignore_attr = TRUE)

  for (cr in c("score", "wald")) {
    a <- rt_screen(fit, d, cr)
    b <- rt_screen(fit, sw_d, cr)
    expect_identical(b$n_unverified, 0L)
    expect_identical(b$sw$steps$variable, a$sw$steps$variable)
    expect_equal(b$sw$steps$p_value, a$sw$steps$p_value, tolerance = 1e-6)
  }
})
