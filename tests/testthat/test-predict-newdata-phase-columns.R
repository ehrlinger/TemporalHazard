# predict(newdata = ) on a multiphase fit that has BOTH a global covariate
# and phase-formula covariates.  A phase without its own formula inherits the
# global design, so at newdata it must be rebuilt from the global formula's
# columns -- not from every non-time column of newdata, which made the two
# kinds of phase impossible to satisfy with one newdata.

.pn_data <- function(n = 500, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rbinom(n, 1, 0.5),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  # A genuine two-phase process, so the early phase is identified and the
  # Hessian inverts (se.fit needs it).  Early: H1(t) = mu1 t / (t + 0.3).
  mu1 <- 0.6 * exp(0.8 * d$z)
  u <- stats::rexp(n)
  t1 <- ifelse(u < mu1, 0.3 * u / (mu1 - u), Inf)
  t2 <- stats::rexp(n, 0.15 * exp(0.4 * d$x + 0.5 * (d$g == "b")))
  tt <- pmin(t1, t2)
  cens <- stats::runif(n, 1, 8)
  d$time <- pmin(tt, cens) + 1e-3
  d$status <- as.integer(tt <= cens)
  d
}

.pn_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(1)
      cache <<- hazard(
        survival::Surv(time, status) ~ x + g, data = .pn_data(),
        dist = "multiphase",
        phases = list(
          early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                               fixed = "shapes", formula = ~ z),
          constant = hzr_phase("constant")
        ),
        fit = TRUE, control = list(n_starts = 1, maxit = 500)
      )
    }
    cache
  }
})

# Structural reference: H_j(t | x) = exp(x beta_j) H0_j(t).  H0_j comes from a
# time-only newdata (every covariate at 0); beta_j is read from theta by name.
.pn_reference <- function(fit, time, x, z, g) {
  th <- fit$fit$theta
  base <- predict(fit, newdata = data.frame(time = time),
                  type = "cumulative_hazard", decompose = TRUE)
  early <- base$early * exp(th[["early.z"]] * z)
  constant <- base$constant *
    exp(th[["constant.x"]] * x +
          th[["constant.gb"]] * (g == "b") + th[["constant.gc"]] * (g == "c"))
  list(early = early, constant = constant, total = early + constant)
}

.pn_time <- c(0.1, 0.5, 1, 2, 4)

test_that("the fixture fits, with covariates that move each phase", {
  fit <- .pn_fit()
  expect_true(isTRUE(fit$fit$converged))
  th <- fit$fit$theta
  # A covariate effect of 0 would make the reference below unable to tell
  # the right design from a dropped column.
  expect_gt(abs(th[["early.z"]]), 0.3)
  expect_gt(abs(th[["constant.x"]]), 0.1)
  expect_gt(abs(th[["constant.gb"]]), 0.05)
  expect_true(all(is.finite(diag(vcov(fit))[c("early.z", "constant.x")])))
})

test_that("predict(newdata) matches exp(x beta_j) H0_j(t) per phase", {
  fit <- .pn_fit()
  ref <- .pn_reference(fit, .pn_time, x = 0.7, z = 1, g = "b")
  nd <- data.frame(time = .pn_time, x = 0.7, z = 1, g = "b")

  total <- predict(fit, newdata = nd, type = "cumulative_hazard")
  expect_equal(unname(total), ref$total, tolerance = 1e-10)

  dec <- predict(fit, newdata = nd, type = "cumulative_hazard",
                 decompose = TRUE)
  expect_equal(dec$early, ref$early, tolerance = 1e-10)
  expect_equal(dec$constant, ref$constant, tolerance = 1e-10)

  surv <- predict(fit, newdata = nd, type = "survival")
  expect_equal(unname(surv), exp(-ref$total), tolerance = 1e-10)
})

test_that("each row of newdata uses its own covariates", {
  fit <- .pn_fit()
  nd <- data.frame(time = 1, x = c(-1, 0, 1), z = c(0, 1, 0),
                   g = c("a", "b", "c"))
  got <- predict(fit, newdata = nd, type = "cumulative_hazard")
  want <- vapply(seq_len(nrow(nd)), function(i) {
    .pn_reference(fit, 1, nd$x[i], nd$z[i], nd$g[i])$total
  }, numeric(1))
  expect_equal(unname(got), want, tolerance = 1e-10)
  # The rows differ, so a design that ignored a covariate would fail here.
  expect_gt(max(want) / min(want), 1.2)
})

test_that("extra and reordered newdata columns give the same answer", {
  fit <- .pn_fit()
  nd <- data.frame(time = .pn_time, x = 0.7, z = 1, g = "b")
  base <- predict(fit, newdata = nd, type = "cumulative_hazard")

  reordered <- nd[, c("g", "z", "time", "x")]
  expect_equal(predict(fit, newdata = reordered, type = "cumulative_hazard"),
               base, tolerance = 1e-12)

  extra <- cbind(nd, unused = 99, also_unused = "q")
  expect_equal(predict(fit, newdata = extra, type = "cumulative_hazard"),
               base, tolerance = 1e-12)

  # A factor carrying the fit's levels is the same design as its label.
  as_factor <- nd
  as_factor$g <- factor("b", levels = c("a", "b", "c"))
  expect_equal(predict(fit, newdata = as_factor, type = "cumulative_hazard"),
               base, tolerance = 1e-12)
})

test_that("a newdata missing a global covariate is an error naming it", {
  fit <- .pn_fit()
  expect_error(
    predict(fit, newdata = data.frame(time = 1, z = 1, g = "a"),
            type = "cumulative_hazard"),
    "missing the global covariate\\(s\\) 'x'"
  )
})

# The design helper needs only the stored design, so these use unfitted
# objects: the fit is irrelevant to how the columns are rebuilt.
.pn_unfitted <- function(formula, data = .pn_data()) {
  hazard(formula, data = data, dist = "multiphase",
         phases = list(
           early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                fixed = "shapes", formula = ~ z),
           constant = hzr_phase("constant")
         ))
}

test_that("data-dependent terms keep the fit's parameters at newdata", {
  d <- .pn_data()
  rows <- c(1, 7, 42)

  # scale() and poly() are computed from the data they see.  At newdata they
  # must reuse the fit's centre, scale and basis, or three rows would be
  # rescaled among themselves -- a wrong design with no error.
  for (f in list(survival::Surv(time, status) ~ scale(x),
                 survival::Surv(time, status) ~ poly(x, 2))) {
    obj <- .pn_unfitted(f, d)
    nd <- data.frame(time = 1, x = d$x[rows], z = 0)
    got <- TemporalHazard:::.hzr_global_design(obj, nd)
    expect_equal(unname(got), unname(obj$data$x[rows, , drop = FALSE]),
                 tolerance = 1e-12)
    # A single row, where recomputing gives NaN or an error.
    one <- TemporalHazard:::.hzr_global_design(obj, nd[1, ])
    expect_equal(unname(one), unname(obj$data$x[rows[1], , drop = FALSE]),
                 tolerance = 1e-12)
  }
})

test_that("a variable from the formula's environment is not a newdata column", {
  cutoff <- 0.5
  obj <- .pn_unfitted(survival::Surv(time, status) ~ I(x > cutoff))
  got <- TemporalHazard:::.hzr_global_design(
    obj, data.frame(time = 1, x = c(0, 1), z = 0))
  expect_equal(unname(got[, 1]), c(0, 1))
})

test_that("se.fit = TRUE is centred on the structural reference", {
  fit <- .pn_fit()
  ref <- .pn_reference(fit, .pn_time, x = 0.7, z = 1, g = "b")
  nd <- data.frame(time = .pn_time, x = 0.7, z = 1, g = "b")

  se <- predict(fit, newdata = nd, type = "cumulative_hazard", se.fit = TRUE)
  expect_equal(se$fit, ref$total, tolerance = 1e-10)
  expect_true(all(is.finite(se$se.fit) & se$se.fit > 0))
  expect_true(all(se$lower < se$fit & se$fit < se$upper))

  sd <- predict(fit, newdata = nd, type = "cumulative_hazard", se.fit = TRUE,
                decompose = TRUE)
  expect_equal(sd$fit[sd$component == "early"], ref$early, tolerance = 1e-10)
  expect_equal(sd$fit[sd$component == "constant"], ref$constant,
               tolerance = 1e-10)
  expect_true(all(is.finite(sd$se.fit) & sd$se.fit > 0))
})
