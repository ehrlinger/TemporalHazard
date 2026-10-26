# A multiphase fit saved before attr(x_list, "from_formula") is classified by
# its stored design. Under time_windows a phase formula can use columns NAMED
# like the expanded global ones (x1_w1, x1_w2) that hold other covariates;
# names alone read it as inherited and predict(newdata = ) and hzr_gof()
# swapped its covariates for the global ones (#658). The fit with its record
# is the known positive for each comparison.

build_legacy_windows_fixture <- function() {
  set.seed(658)
  n <- 400
  # Not centred, so a swapped design moves the mean curve as well.
  d <- data.frame(x1 = stats::rnorm(n), x1_w1 = stats::rnorm(n, mean = 1),
                  x1_w2 = stats::rnorm(n, mean = -1))
  tt <- stats::rweibull(n, shape = 1.2,
                        scale = 6 * exp(-0.5 * d$x1 - 0.4 * d$x1_w1))
  cc <- stats::runif(n, 2, 12)
  d$time <- pmin(tt, cc)
  d$status <- as.integer(tt <= cc)
  own <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1, data = d,
    dist = "multiphase", time_windows = 2,
    phases = list(
      early    = hzr_phase("cdf", t_half = 2, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ x1_w1 + x1_w2),
      constant = hzr_phase("constant", formula = ~ 1)
    ),
    fit = TRUE
  ))
  inh <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ x1, data = d,
    dist = "multiphase", time_windows = 2,
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  strip <- function(f) {
    attr(f$fit$x_list, "from_formula") <- NULL
    f
  }
  list(d = d, own = own, own_legacy = strip(own), inh = inh,
       inh_legacy = strip(inh))
}

# Fitted once per file, not once per block: the blocks only read it.
legacy_windows_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- build_legacy_windows_fixture()
    cache
  }
})

test_that("legacy time_windows fit: inheritance is read from the design's values", {
  fx <- legacy_windows_fixture()
  expect_true(isTRUE(fx$own$fit$converged))
  expect_true(isTRUE(fx$inh$fit$converged))
  expect_identical(colnames(fx$own$fit$x_list$early), c("x1_w1", "x1_w2"))
  expect_false(.hzr_phase_inherits_global(fx$own, "early"))
  expect_false(.hzr_phase_inherits_global(fx$own_legacy, "early"))
  # A legacy phase that did inherit is still read as inheriting.
  expect_true(.hzr_phase_inherits_global(fx$inh_legacy, "early"))
  expect_true(.hzr_phase_inherits_global(fx$inh_legacy, "constant"))
})

test_that("legacy time_windows fit: predict(newdata = ) keeps a phase formula's own covariates", {
  fx <- legacy_windows_fixture()
  nd <- fx$d[, c("x1", "x1_w1", "x1_w2")]
  for (t in c(1, 6)) {
    ndt <- cbind(nd, time = t)
    rec <- predict(fx$own, newdata = ndt, type = "cumulative_hazard")
    # The own-formula design is not windowed, so the stored design is the
    # design at any time: an oracle independent of predict(newdata = ).
    oracle <- .hzr_multiphase_cumhaz(
      rep(t, nrow(nd)), fx$own$fit$theta, fx$own$fit$phases,
      fx$own$fit$covariate_counts, fx$own$fit$x_list
    )
    expect_equal(rec, oracle, tolerance = 1e-12)
    expect_equal(predict(fx$own_legacy, newdata = ndt,
                         type = "cumulative_hazard"),
                 oracle, tolerance = 1e-12)
    # A legacy fit that inherited still gets the windowed global design.
    expect_equal(predict(fx$inh_legacy, newdata = ndt,
                         type = "cumulative_hazard"),
                 predict(fx$inh, newdata = ndt, type = "cumulative_hazard"),
                 tolerance = 1e-12)
  }
})

test_that("legacy time_windows fit: hzr_gof()'s curve keeps a phase formula's own covariates", {
  fx <- legacy_windows_fixture()
  rec <- suppressWarnings(hzr_gof(fx$own))
  leg <- suppressWarnings(hzr_gof(fx$own_legacy))
  for (v in c("par_cumhaz", "par_surv", "par_cumhaz_early")) {
    expect_equal(leg[[v]], rec[[v]], tolerance = 1e-12, info = v)
  }
  inh_rec <- suppressWarnings(hzr_gof(fx$inh))
  inh_leg <- suppressWarnings(hzr_gof(fx$inh_legacy))
  expect_equal(inh_leg$par_cumhaz, inh_rec$par_cumhaz, tolerance = 1e-12)
})

test_that("legacy time_windows fit: an ignored phase formula still reads as inherited", {
  # A pre-#299 vector fit wrote a phase formula it could not evaluate and took
  # the window-expanded global x instead. Without the record, its design is
  # the inherited one in names AND values, so it must still inherit.
  fx <- legacy_windows_fixture()
  d <- fx$d
  ignored <- suppressWarnings(hzr_saved_before_299(
    time = d$time, status = d$status, x = cbind(x1 = d$x1),
    time_windows = 2, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ x1),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  expect_true(isTRUE(ignored$fit$converged))
  expect_identical(colnames(ignored$fit$x_list$early), c("x1_w1", "x1_w2"))
  expect_false(is.null(ignored$fit$phases$early$formula))
  legacy <- ignored
  attr(legacy$fit$x_list, "from_formula") <- NULL
  expect_true(.hzr_phase_inherits_global(legacy, "early"))

  ndt <- data.frame(x1 = d$x1, time = 6)
  expect_equal(predict(legacy, newdata = ndt, type = "cumulative_hazard"),
               predict(ignored, newdata = ndt, type = "cumulative_hazard"),
               tolerance = 1e-12)
  expect_equal(suppressWarnings(hzr_gof(legacy))$par_cumhaz,
               suppressWarnings(hzr_gof(ignored))$par_cumhaz,
               tolerance = 1e-12)
})
