# Under objective = "sas" an interval-censored row contributes PROC HAZARD's
# interval-mean-hazard term, so the maximised value is not a log-likelihood.
# print(), summary() and the stepwise trace called it one (#544, #556): a
# reader comparing it with another model's log-likelihood, or taking an LR
# or AIC from it, got a finite, plausible and wrong number.

sas_label_data <- function(n = 120, intervals = TRUE) {
  withr::local_seed(544)
  st <- rep(if (intervals) c(1, 0, 2) else c(1, 0, 0), length.out = n)
  lo <- seq(0.10, 2.40, length.out = n)
  up <- lo + seq(0.20, 4.00, length.out = n)
  tt <- ifelse(st == 2, up, lo + 0.4)
  data.frame(tt = tt, st = st, lo = ifelse(st == 2, lo, 0),
             x = stats::rbinom(n, 1, 0.5))
}

sas_label_fit <- function(d, objective) {
  suppressWarnings(hazard(
    time = d$tt, status = d$st, time_lower = d$lo, time_upper = d$tt,
    dist = "multiphase", objective = objective, fit = TRUE, data = d,
    phases = list(early = hzr_phase("cdf", t_half = 0.8, nu = 0.4, m = -0.6),
                  constant = hzr_phase("constant")),
    control = list(n_starts = 1)))
}

printed <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")

test_that("a SAS objective is not printed or summarised as a log-likelihood (#544)", {
  skip_on_cran()
  d <- sas_label_data()
  fs <- sas_label_fit(d, "sas")
  fl <- sas_label_fit(d, "likelihood")
  # The premise: the two objectives differ on these data, and the SAS value
  # is not the log-likelihood at its own estimates.
  ll_at_sas <- as.numeric(hzr_evaluate(fl, coef(fs))$logLik)
  expect_gt(abs(fs$fit$objective - ll_at_sas), 1)

  out <- printed(fs)
  expect_match(out, "SAS objective:", fixed = TRUE)
  expect_false(grepl("log-lik:", out, fixed = TRUE))

  s <- summary(fs)
  expect_identical(s$log_lik, NA_real_)
  expect_identical(s$objective, "sas")
  expect_identical(s$objective_value, fs$fit$objective)
  out_s <- printed(s)
  expect_match(out_s, "SAS objective:", fixed = TRUE)
  expect_false(grepl("log-lik:", out_s, fixed = TRUE))

  # Known positive: the likelihood fit keeps its label and its value.
  expect_match(printed(fl), "log-lik:", fixed = TRUE)
  sl <- summary(fl)
  expect_identical(sl$log_lik, fl$fit$objective)
  expect_identical(sl$objective, "likelihood")
  expect_match(printed(sl), "log-lik:", fixed = TRUE)
})

test_that("without interval-censored rows the SAS objective is the log-likelihood (#544)", {
  skip_on_cran()
  d <- sas_label_data(intervals = FALSE)
  fs <- sas_label_fit(d, "sas")
  fl <- sas_label_fit(d, "likelihood")
  expect_equal(fs$fit$objective, fl$fit$objective, tolerance = 1e-8)
  expect_match(printed(fs), "log-lik:", fixed = TRUE)
  expect_identical(summary(fs)$log_lik, fs$fit$objective)
})

test_that("interval rows the fit did not read do not change the label (#544)", {
  skip_on_cran()
  # Every interval row given weight 0: the fitted SAS objective reads none of
  # them, so it is the log-likelihood and keeps that label.
  d <- sas_label_data()
  w <- ifelse(d$st == 2, 0, 1)
  fit <- suppressWarnings(hazard(
    time = d$tt, status = d$st, time_lower = d$lo, time_upper = d$tt,
    weights = w, dist = "multiphase", objective = "sas", fit = TRUE,
    phases = list(early = hzr_phase("cdf", t_half = 0.8, nu = 0.4, m = -0.6),
                  constant = hzr_phase("constant")),
    control = list(n_starts = 1)))
  expect_true(any(fit$data$status == 2))
  expect_match(printed(fit), "log-lik:", fixed = TRUE)
  expect_identical(summary(fit)$log_lik, fit$fit$objective)
})

test_that("interval rows a phase design dropped do not change the label (#544)", {
  skip_on_cran()
  # A phase covariate missing on every interval row: the multiphase design
  # drops those rows (fit$fit$rows_used), so the SAS objective reads none.
  d <- sas_label_data()
  d$z <- ifelse(d$st == 2, NA, stats::rnorm(nrow(d)))
  fit <- suppressWarnings(hazard(
    time = d$tt, status = d$st, time_lower = d$lo, time_upper = d$tt,
    data = d, dist = "multiphase", objective = "sas", fit = TRUE,
    phases = list(early = hzr_phase("cdf", t_half = 0.8, nu = 0.4, m = -0.6,
                                    formula = ~ z),
                  constant = hzr_phase("constant")),
    control = list(n_starts = 1)))
  # The premise: the rows are stored, and the fit did not read them.
  expect_true(any(fit$data$status == 2))
  expect_false(any(fit$fit$rows_used[fit$data$status == 2]))
  expect_match(printed(fit), "log-lik:", fixed = TRUE)
  expect_identical(summary(fit)$log_lik, fit$fit$objective)
})

test_that("summary()'s new fields come last (#544)", {
  skip_on_cran()
  s <- summary(sas_label_fit(sas_label_data(), "sas"))
  n <- length(names(s))
  expect_identical(names(s)[c(n - 1L, n)], c("objective", "objective_value"))
})

test_that("the formula interface's interval rows are seen too (#544)", {
  skip_on_cran()
  d <- sas_label_data()
  d$lo2 <- ifelse(d$st == 2, d$lo, d$tt)
  d$hi2 <- ifelse(d$st == 0, Inf, d$tt)
  fit <- suppressWarnings(hazard(
    survival::Surv(lo2, hi2, type = "interval2") ~ 1, data = d,
    dist = "multiphase", objective = "sas", fit = TRUE,
    phases = list(early = hzr_phase("cdf", t_half = 0.8, nu = 0.4, m = -0.6),
                  constant = hzr_phase("constant")),
    control = list(n_starts = 1)))
  # The premise: Surv's interval code (3) was stored as this package's 2.
  expect_equal(as.numeric(fit$data$status), d$st)
  expect_match(printed(fit), "SAS objective:", fixed = TRUE)
  expect_identical(summary(fit)$log_lik, NA_real_)
})

test_that("an AIC screen on a SAS objective says it is not an AIC (#544)", {
  skip_on_cran()
  d <- sas_label_data()
  fs <- sas_label_fit(d, "sas")
  classes <- character()
  withCallingHandlers(
    suppressMessages(hzr_stepwise(fs, scope = list(constant = ~ x), data = d,
                                  criterion = "aic", trace = FALSE)),
    warning = function(w) {
      classes <<- c(classes, class(w)[1L])
      invokeRestart("muffleWarning")
    })
  expect_identical(sum(classes == "hzr_stepwise_sas_objective"), 1L)
  # Known positive: a Wald screen on the same fit does not warn so.
  expect_no_warning(
    withCallingHandlers(
      suppressMessages(hzr_stepwise(fs, scope = list(constant = ~ x),
                                    data = d, criterion = "wald",
                                    trace = FALSE)),
      warning = function(w) {
        if (!inherits(w, "hzr_stepwise_sas_objective")) {
          invokeRestart("muffleWarning")
        }
      }),
    class = "hzr_stepwise_sas_objective")
})

test_that("the stepwise trace labels a SAS objective as one (#556)", {
  skip_on_cran()
  d <- sas_label_data()
  fs <- sas_label_fit(d, "sas")
  sw <- suppressWarnings(suppressMessages(hzr_stepwise(
    fs, scope = list(constant = ~ x), data = d, criterion = "wald",
    trace = FALSE)))
  trace <- paste(sw$trace_msg, collapse = "\n")
  expect_match(trace, "SAS objective = ", fixed = TRUE)
  expect_false(grepl("logLik = ", trace, fixed = TRUE))
})
