# hzr_bootstrap(id =) resamples patients, not rows (#667).
#
# Time-varying covariate data carry one row per patient-segment: the
# counting-process layout of SAS `lcensor <start>; time <stop>;`. Resampling
# rows splits a patient across replicates, so a replicate can hold two of a
# patient's three segments, or one of them twice. SAS %HAZBOOT draws patients
# under `_TVC=1, _IDVAR=` and keeps every segment of each one drawn. These
# tests check that per patient and per replicate, never on a summary.

# Several segments per patient, a covariate that changes between them.
tvc_fixture <- function(n_pat = 60, seed = 11) {
  set.seed(seed)
  rows <- lapply(seq_len(n_pat), function(i) {
    k <- sample(1:3, 1)
    cuts <- sort(stats::runif(k, 0, 5))
    start <- c(0, utils::head(cuts, -1))
    data.frame(pid = i, start = start, stop = cuts,
               x = stats::rnorm(k), event = c(rep(0L, k - 1L),
                                              stats::rbinom(1, 1, 0.6)))
  })
  do.call(rbind, rows)
}

tvc_fit <- function(d) {
  hazard(survival::Surv(start, stop, event) ~ x, data = d,
         dist = "weibull", theta = c(mu = 0.3, nu = 1, beta = 0),
         fit = TRUE)
}

test_that(".hzr_bootstrap_index() draws whole patients", {
  id <- rep(c("a", "b", "c", "d", "e"), times = c(3, 1, 2, 4, 2))
  rows_of <- split(seq_along(id), id)
  set.seed(3)
  for (r in 1:50) {
    idx <- TemporalHazard:::.hzr_bootstrap_index(length(id), 1, id)
    # Every patient drawn k times contributes exactly k copies of each of
    # its rows, so each of its rows appears the same number of times.
    for (p in names(rows_of)) {
      counts <- vapply(rows_of[[p]], function(j) sum(idx == j), integer(1))
      expect_length(unique(counts), 1L)
    }
    # Five patients drawn, whatever their row counts.
    n_drawn <- sum(vapply(rows_of, function(rr) sum(idx == rr[1]),
                          integer(1)))
    expect_identical(n_drawn, 5L)
  }
})

test_that(".hzr_bootstrap_index() row-level default is unchanged", {
  set.seed(9)
  a <- TemporalHazard:::.hzr_bootstrap_index(20L, 1, NULL)
  set.seed(9)
  b <- sample.int(20L, size = 20L, replace = TRUE)
  expect_identical(a, b)
})

test_that("fraction scales patients, not rows, under id =", {
  id <- rep(1:10, times = rep(c(1, 5), 5))   # 30 rows, 10 patients
  set.seed(4)
  idx <- TemporalHazard:::.hzr_bootstrap_index(length(id), 0.5, id)
  drawn <- vapply(split(seq_along(id), id), function(rr) sum(idx == rr[1]),
                  integer(1))
  expect_identical(sum(drawn), 5L)
})

test_that("id = accepts a column name or a row-length vector, and refuses bad ids", {
  d <- tvc_fixture()
  fit <- tvc_fit(d)
  b_col <- hzr_bootstrap(fit, n_boot = 3, seed = 1, id = "pid")
  b_vec <- hzr_bootstrap(fit, n_boot = 3, seed = 1, id = d$pid)
  expect_identical(b_col$replicates, b_vec$replicates)

  expect_error(hzr_bootstrap(fit, n_boot = 2, id = "no_such_col"),
               "no_such_col")
  expect_error(hzr_bootstrap(fit, n_boot = 2, id = d$pid[-1]),
               "one value per row")
  bad <- d$pid
  bad[2] <- NA
  expect_error(hzr_bootstrap(fit, n_boot = 2, id = bad), "missing")
})

test_that("patient resampling changes the replicates", {
  d <- tvc_fixture()
  fit <- tvc_fit(d)
  b_row <- suppressWarnings(hzr_bootstrap(fit, n_boot = 5, seed = 2))
  b_id  <- hzr_bootstrap(fit, n_boot = 5, seed = 2, id = "pid")
  expect_identical(b_id$n_success, 5L)
  # Same seed, different resampling unit: the replicates must differ.
  expect_false(isTRUE(all.equal(b_row$replicates$estimate,
                                b_id$replicates$estimate)))
})

test_that("entry times without id = warn once; with id = they do not", {
  d <- tvc_fixture()
  fit <- tvc_fit(d)
  expect_warning(hzr_bootstrap(fit, n_boot = 2, seed = 1),
                 "id =", fixed = TRUE)
  expect_no_warning(hzr_bootstrap(fit, n_boot = 2, seed = 1, id = "pid"))

  # No entry times: no warning.
  d1 <- d[!duplicated(d$pid), ]
  d1$start <- 0
  fit1 <- tvc_fit(d1)
  expect_no_warning(hzr_bootstrap(fit1, n_boot = 2, seed = 1))
})

test_that("the vector interface takes id = as a vector", {
  d <- tvc_fixture()
  # Intercept only: a vector fit given `x` directly cannot be bootstrapped
  # at all (its design is not resampled with the rows), with or without id.
  fit <- hazard(time = d$stop, status = d$event, time_lower = d$start,
                dist = "weibull", theta = c(mu = 0.3, nu = 1), fit = TRUE)
  b <- hzr_bootstrap(fit, n_boot = 3, seed = 1, id = d$pid)
  expect_identical(b$n_success, 3L)
  expect_error(hzr_bootstrap(fit, n_boot = 2, id = "pid"), "data")
})

test_that("the entry-time warning uses hazard()'s genuine-entry rule", {
  has_entry <- TemporalHazard:::.hzr_bootstrap_has_entry
  expect_false(has_entry(list(time = 1:3, status = c(0, 1, 0),
                              time_lower = NULL)))
  # Interval bounds on status-2 rows are not entries.
  expect_false(has_entry(list(time = 1:3, status = c(2, 2, 2),
                              time_lower = c(0.5, 1, 2))))
  # time_lower equal to time, or 0, means no entry.
  expect_false(has_entry(list(time = 1:3, status = c(0, 1, 0),
                              time_lower = c(0, 2, 3))))
  expect_true(has_entry(list(time = 1:3, status = c(0, 1, 0),
                             time_lower = c(0, 1.5, 3))))
})
