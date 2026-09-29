# hzr_deciles() took each subject's expected events as its cumulative hazard
# at exit, and counted events unweighted (#491). A left-truncated subject is
# at risk only from its entry time, so its expected count is H(exit) -
# H(entry), and a weighted fit conserves weighted events, sum(w * H) =
# sum(w * d). hzr_gof() already did both; a correctly specified fit read as
# miscalibrated in hzr_deciles(). Both fixtures conserve events, so the
# decile totals must match the observed events.

.dec_avc <- local({
  data(avc, package = "TemporalHazard")
  na.omit(avc)
})

.dec_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}

test_that("hzr_deciles() subtracts the entry-time cumulative hazard (#491)", {
  avc <- .dec_avc
  set.seed(2)
  late <- stats::runif(nrow(avc)) < 0.4
  entry <- ifelse(late, avc$int_dead * stats::runif(nrow(avc), 0.1, 0.9), 0)
  fit <- hazard(
    time       = avc$int_dead,
    status     = avc$dead,
    time_lower = entry,
    x          = as.matrix(avc[, c("age", "mal")]),
    dist       = "multiphase",
    phases     = .dec_phases(),
    fit        = TRUE,
    control    = list(n_starts = 1, maxit = 1000, conserve = TRUE)
  )
  expect_true(isTRUE(fit$spec$control$conserve_applied))

  nd_exit  <- data.frame(age = avc$age, mal = avc$mal, time = avc$int_dead)
  nd_entry <- data.frame(age = avc$age, mal = avc$mal, time = entry)
  h_exit  <- predict(fit, newdata = nd_exit,  type = "cumulative_hazard")
  h_entry <- predict(fit, newdata = nd_entry, type = "cumulative_hazard")
  # The entry term is not negligible, so omitting it would be caught.
  expect_gt(sum(h_entry), 10)

  dec <- hzr_deciles(fit, time = 12)
  ov <- attr(dec, "overall")
  expect_equal(ov$total_events, sum(avc$dead))
  expect_equal(ov$total_expected, sum(h_exit - h_entry), tolerance = 1e-8)
  expect_equal(ov$total_expected / ov$total_events, 1, tolerance = 1e-6)
  # The groups carry the same totals as the summary.
  expect_equal(sum(dec$expected), ov$total_expected, tolerance = 1e-10)
  expect_equal(sum(dec$events), ov$total_events)
})

test_that("hzr_deciles() weights events and expected events (#491)", {
  avc <- .dec_avc
  set.seed(1)
  w <- stats::runif(nrow(avc), 0.5, 2)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data    = avc,
    weights = w,
    dist    = "multiphase",
    phases  = .dec_phases(),
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 1000, conserve = TRUE)
  )
  expect_true(isTRUE(fit$spec$control$conserve_applied))
  # The weights matter: the unweighted and weighted event counts differ.
  expect_gt(abs(sum(w * avc$dead) - sum(avc$dead)), 5)

  dec <- hzr_deciles(fit, time = 12)
  ov <- attr(dec, "overall")
  expect_equal(ov$total_events, sum(w * avc$dead), tolerance = 1e-10)
  expect_equal(ov$total_expected / ov$total_events, 1, tolerance = 1e-6)
  expect_equal(sum(dec$events), ov$total_events, tolerance = 1e-10)
  expect_equal(sum(dec$expected), ov$total_expected, tolerance = 1e-10)
  # n stays a head count.
  expect_equal(sum(dec$n), nrow(avc))
})

.dec_weibull <- function(avc, w) {
  hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data = avc, weights = w, dist = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit = TRUE
  )
}

test_that("hzr_deciles() chi-square does not move with the weights' scale", {
  avc <- .dec_avc
  set.seed(1)
  w <- stats::runif(nrow(avc), 0.5, 2)
  d1  <- hzr_deciles(.dec_weibull(avc, w), time = 12)
  d10 <- hzr_deciles(.dec_weibull(avc, 10 * w), time = 12)
  # Known positive: the weighted counts do scale, by exactly 10.
  expect_equal(attr(d10, "overall")$total_events /
                 attr(d1, "overall")$total_events, 10, tolerance = 1e-10)
  # The test statistic, its p-value and the rates do not.
  expect_gt(attr(d1, "overall")$chi_sq, 0.5)
  expect_equal(attr(d10, "overall")$chi_sq / attr(d1, "overall")$chi_sq, 1,
               tolerance = 1e-5)
  # The two fits agree to optimizer precision, not bit for bit, which moves
  # the p-values in the fifth digit. Scale-dependence moved them by orders of
  # magnitude.
  expect_equal(d10$p_value, d1$p_value, tolerance = 1e-3)
  expect_equal(d10$observed_rate, d1$observed_rate, tolerance = 1e-8)
  expect_equal(d10$expected_rate, d1$expected_rate, tolerance = 1e-5)
})

test_that("hzr_deciles() with unit weights is the unweighted result", {
  avc <- .dec_avc
  d0 <- hzr_deciles(.dec_weibull(avc, NULL), time = 12)
  d1 <- hzr_deciles(.dec_weibull(avc, rep(1, nrow(avc))), time = 12)
  expect_gt(attr(d0, "overall")$chi_sq, 0.5)
  expect_equal(as.data.frame(unclass(d1)), as.data.frame(unclass(d0)),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(attr(d1, "overall")$chi_sq, attr(d0, "overall")$chi_sq,
               tolerance = 1e-8)
})

test_that("an interval's lower bound is not taken as an entry time", {
  avc <- .dec_avc
  set.seed(3)
  ic <- sample(which(avc$dead == 1), 30)
  lo <- avc$int_dead
  lo[ic] <- avc$int_dead[ic] * 0.7
  st <- avc$dead
  st[ic] <- 2L
  fit <- hazard(
    time = avc$int_dead, status = st, time_lower = lo,
    x = as.matrix(avc[, c("age", "mal")]), dist = "weibull",
    theta = c(0.01, 0.5, 0, 0), fit = TRUE
  )
  nd_hi <- data.frame(age = avc$age, mal = avc$mal, time = avc$int_dead)
  nd_lo <- data.frame(age = avc$age, mal = avc$mal, time = lo)
  h_hi <- predict(fit, newdata = nd_hi, type = "cumulative_hazard")
  h_lo <- predict(fit, newdata = nd_lo, type = "cumulative_hazard")
  # The lower bounds carry real cumulative hazard, so subtracting them shows.
  expect_gt(sum(h_lo[ic]), 1)

  dec <- hzr_deciles(fit, time = 12, status = avc$dead)
  expect_equal(attr(dec, "overall")$total_expected, sum(h_hi),
               tolerance = 1e-8)
})
