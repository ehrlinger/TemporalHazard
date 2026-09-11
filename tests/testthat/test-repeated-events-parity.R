# Parity with SAS for hz.ce_cardioversion_repeated.ehb (maze, permanent AF).
#
# The job builds bd_card, runs %repeat on it, adjusts one line, and fits the
# result. bd_card was never saved, so .hzr_derive_bd_card() rebuilds it from
# the permanent datasets. Every expected value below comes from the job's own
# .log (each DATA step's output shape) or .lst (the counts HAZARD printed) --
# none from the R code.
#
# The data are PHI and live only on the study volume. This test skips when the
# volume is not mounted, which is every CI runner, so a green CI run is NOT
# evidence of parity. Only a local run with the volume mounted is.

skip_if_no_maze_datasets <- function() {
  d <- .hzr_maze_datasets_dir() # nolint: object_usage_linter.
  testthat::skip_if(is.na(d), "maze study datasets not mounted")
  d
}

test_that("hz.ce_cardioversion_repeated.ehb: %repeat reproduces the SAS log and HAZARD's counts", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  dir <- skip_if_no_maze_datasets()

  bd_card <- .hzr_derive_bd_card(dir)
  # The rebuild, against the .log: WORK.BD_CARD has 3246 observations and 357 variables.
  expect_equal(dim(bd_card), c(3246L, 357L))

  id <- "ccfid"
  tm <- "iv_event"
  fu <- "iv_end"
  ind <- "ce_card"
  s1 <- .hzr_re_stage1(bd_card)
  s2 <- .hzr_re_stage2(s1, id, tm, ind)
  s3 <- .hzr_re_stage3(s2, id, tm, fu, ind)
  s4 <- .hzr_re_stage4(s3, id, tm, fu, ind)
  s5 <- .hzr_re_stage5(s4, ind)
  s6 <- .hzr_re_stage6(s5, id, tm, fu)
  s7 <- .hzr_re_stage7(s6, id, tm, ind)

  # Row coverage first: each DATA step's output shape, from .log lines 196, 211,
  # 229, 244, 252, 267 and 290. Stages 6 and 7 carry two fewer columns than SAS
  # because lag_iv and number, the macro's loop counters, are not returned.
  shapes <- rbind(dim(s1), dim(s2), dim(s3), dim(s4), dim(s5), dim(s6), dim(s7))
  expect_equal(shapes[, 1], c(3246L, 709L, 709L, 963L, 963L, 963L, 962L))
  expect_equal(shapes[, 2], c(358L, 358L, 358L, 360L, 361L, 367L - 2L, 367L - 2L))

  # The exported function is the same seven stages, and this data trips none of
  # its input warnings. Those warnings name subjects, so they are counted here,
  # never printed; and the comparison with the stages reports a bare FALSE
  # rather than a diff of patient rows.
  warned <- testthat::capture_warnings(events <- hzr_repeated_events(bd_card, id, tm, fu, ind))
  expect_length(warned, 0L)
  expect_true(isTRUE(all.equal(events, `row.names<-`(s7, NULL))))

  # The job's line 65, after the macro. It moves an event time that does not
  # exceed its start; it removes no rows (.log line 315: 962 in, 962 out).
  nudge <- which(events$iv_start >= events$iv_event)
  events$iv_event[nudge] <- events$iv_event[nudge] + 0.0001141553

  # HAZARD's own tallies, from the .lst.
  expect_false(anyNA(events$ce_card))
  expect_equal(nrow(events), 962L)
  expect_equal(sum(events$ce_card == 1), 388L)
  expect_equal(sum(events$ce_card != 1), 574L)
  expect_equal(sum(events$iv_start > 0), 387L)

  # The ranges, compared as ratios. The minima are smaller than any sensible
  # tolerance, which would make expect_equal() compare absolutely -- an
  # assertion that could not fail.
  left <- events$iv_start[events$iv_start > 0]
  expect_equal(min(events$iv_event) / 0.0001140795, 1, tolerance = 1e-6)
  expect_equal(max(events$iv_event) / 12.99137, 1, tolerance = 1e-6)
  expect_equal(min(left) / 0.0001140795, 1, tolerance = 1e-6)
  expect_equal(max(left) / 9.716832, 1, tolerance = 1e-6)
})

test_that("hz.ce_cardioversion_repeated.ehb: the result does not depend on input row order", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  dir <- skip_if_no_maze_datasets()

  bd_card <- .hzr_derive_bd_card(dir)
  # Ties are what make order matter. Without them the shuffle below proves
  # nothing, so require them. The job feeds the macro from a PROC SQL join,
  # whose output order is not guaranteed.
  expect_gt(sum(duplicated(bd_card[c("ccfid", "iv_event")])), 0)

  # The columns the job's two fits read: the segments, the indicator, and the
  # stratified model's covariates. Six other columns do move with order --
  # other event types' (ce_cva, ce_embol, dt_cva, dt_embol, iv_cva, iv_embol),
  # which neither fit reads.
  used <- c("ccfid", "iv_event", "iv_start", "iv_seg", "ce_card", "event", "event_no",
            "rcensor", "renewal", "first", "last", "maze_prc", "iso_pvi")
  canon <- function(x) {
    x <- x[used]
    x <- x[do.call(order, c(unname(as.list(x)), method = "radix")), , drop = FALSE]
    row.names(x) <- NULL
    x
  }
  # Its warnings name subjects, and would print without failing the test.
  run <- function(d) {
    warned <- testthat::capture_warnings(out <- hzr_repeated_events(d, "ccfid", "iv_event", "iv_end", "ce_card"))
    expect_length(warned, 0L)
    out
  }
  base <- canon(run(bd_card))
  # `canon()` keeps ccfid, so a failure must not print a diff.
  for (k in 1:5) {
    shuffled <- withr::with_seed(k, bd_card[sample(nrow(bd_card)), , drop = FALSE])
    expect_true(isTRUE(all.equal(canon(run(shuffled)), base)), label = paste("shuffle", k))
  }
})

# The fits. The job's PROC HAZARD statements go through hzr_translate_sas(),
# and the calls it emits are evaluated over the rebuilt input, so the model
# specification comes from the translator rather than from this file.

cardioversion_job <- function(dir) {
  path <- file.path(dirname(dir), "distributions", "hz.ce_cardioversion_repeated.ehb.sas")
  testthat::skip_if_not(file.exists(path), "cardioversion job not found beside the datasets")
  # STEEPEST, in both blocks, is SAS's steepest-descent warm-up (issue #145),
  # which the translator leaves out. The other two rows come from its
  # fail-closed scan of steps that may change the %repeat output EVENTS
  # before a fit reads it. The first is the job's line 65 nudge, a genuine
  # rewrite, which this file applies itself in cardioversion_events(). The
  # second is %HAZPLOT(IN=EVENTS, ...) before the stratified fit. It only
  # reads EVENTS, but a macro call naming the output stops by design, since
  # a false stop costs a deleted chunk and a miss fits the wrong data.
  # Nothing else may go missing.
  expect_warning(job <- hzr_translate_sas(path), "STEEPEST")
  expect_equal(
    job$untranslated$construct,
    c("EVENTS changed after %repeat", "STEEPEST", "EVENTS changed after %repeat", "STEEPEST")
  )
  job
}

# The fit input: %repeat and then the job's line 65, as in the first test,
# which checks it against the listing's tallies.
cardioversion_events <- function(dir) {
  events <- hzr_repeated_events(.hzr_derive_bd_card(dir), "ccfid", "iv_event", "iv_end", "ce_card")
  nudge <- which(events$iv_start >= events$iv_event)
  events$iv_event[nudge] <- events$iv_event[nudge] + 0.0001141553
  events
}

# Evaluates one emitted status chunk and fit chunk, changing only the control
# entries named in `control`. SAS names are case-insensitive, and the
# translator writes them in upper case.
fit_translated <- function(job, events, fit_slot, status_slot, control) {
  env <- new.env()
  env$EVENTS <- stats::setNames(events, toupper(names(events)))
  eval(job$calls[[status_slot]], env)
  cl <- job$calls[[fit_slot]][[3L]]
  cl$control[names(control)] <- control
  eval(cl, env)
}

# Row coverage: the fit saw what the listing's Initial Summary counted. The
# two status counts sum to 962, so no row carries any other code.
expect_listing_rows <- function(fit) {
  d <- fit$data
  expect_equal(length(d$time), 962L)
  expect_equal(sum(d$status == 1), 388L)
  expect_equal(sum(d$status == 0), 574L)
  expect_equal(sum(d$time_lower > 0), 387L)
}

# `sas` is the listing's Parameter Estimate Summary in R's theta order, named
# by R's parameter names. SAS estimates every shape on the log scale, so for
# the rows marked `logged` R's estimate is logged and its standard error
# converted by the delta method, se(log|x|) = se(x) / |x|.
#
# Tolerances, all measured 2026-09-10 against a fit that stopped short (model
# 2 under default control, 0.013 below the listing's log likelihood), which
# fails every one of them:
# - log likelihood: half a unit in the listing's third decimal;
# - estimates: 1e-3 of the listing's standard error. Both optimizers stop on a
#   tolerance on a flat surface, so the printed digits are not all fixed by
#   the data; a thousandth of a standard error is agreement for any purpose.
# - standard errors: 2e-3 relative. SAS's Hessian is numeric.
expect_matches_listing <- function(fit, sas, ll) {
  f <- fit$fit
  # Not a hollow fit: .hzr_optim_generic() clamps a non-finite objective to
  # 1e10 and still reports convergence, so the log likelihood is checked
  # against the listing, never taken from `converged`.
  expect_true(f$converged)
  expect_lt(abs(f$objective - ll), 5e-4)
  expect_equal(names(f$theta), row.names(sas))
  expect_true(is.matrix(f$vcov))
  se <- sqrt(diag(f$vcov))
  expect_true(all(is.finite(se) & se > 0))

  # The log likelihood at the listing's own estimates, with no optimizer
  # involved, so a failure here is the likelihood and not the search. log|x|
  # drops the sign; every logged parameter in this job is positive.
  expect_true(all(f$theta[sas$logged] > 0))
  d <- fit$data
  ll_at_sas <- .hzr_logl_multiphase( # nolint: object_usage_linter.
    theta = ifelse(sas$logged, exp(sas$est), sas$est),
    time = d$time, status = d$status, time_lower = d$time_lower,
    time_upper = d$time_upper, weights = d$weights,
    phases = f$phases, covariate_counts = f$covariate_counts, x_list = f$x_list
  )
  expect_lt(abs(ll_at_sas - ll), 5e-4)

  est <- ifelse(sas$logged, log(abs(f$theta)), f$theta)
  se <- ifelse(sas$logged, se / abs(f$theta), se)
  expect_lt(max(abs(est - sas$est) / sas$se), 1e-3)
  expect_lt(max(abs(se / sas$se - 1)), 2e-3)
}

sas_table <- function(names, logged, est, se) {
  data.frame(logged = logged, est = est, se = se, row.names = names)
}

test_that("hz.ce_cardioversion_repeated.ehb: the first model reproduces the listing", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  dir <- skip_if_no_maze_datasets()
  job <- cardioversion_job(dir)

  # One start, at the job's PARMS. SAS's iteration 0 already sits at the
  # optimum, and the listing's LL is from that point. The default five starts
  # return this same fit from start 1, but the perturbed starts end at worse
  # local optima (LL -267.914, -268.175, -274.368) and their Hessians warn --
  # noise that says nothing about this fit.
  expect_no_warning(
    fit <- fit_translated(job, cardioversion_events(dir), "fit", "status", list(n_starts = 1L))
  )
  expect_listing_rows(fit)
  # CONSERVE: "Conservation of events: Invoked at each iteration".
  expect_true(fit$spec$control$conserve_applied)

  # .lst lines 100 and 116-125.
  sas <- sas_table(
    c("phase_1.log_mu", "phase_1.log_t_half", "phase_1.nu", "phase_1.m",
      "phase_2.log_mu", "phase_2.log_tau", "phase_2.gamma", "phase_2.alpha", "phase_2.eta"),
    logged = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
    #      E0        E2        E3        E4        L0         L1        L2        L3        L4
    est = c(-1.19976, -4.56749, -1.76693, 1.761105, -1.40628, -0.609944, 1.863922, 0.916978,
            -1.99124),
    se = c(0.09400547, 0.06371024, 0.2689927, 0.3089678, 0.2012175, 0.3764006, 0.6892761,
           0.3171719, 0.8600162)
  )
  expect_matches_listing(fit, sas, ll = -267.885)
})

test_that("hz.ce_cardioversion_repeated.ehb: the stratified model reproduces the listing", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  dir <- skip_if_no_maze_datasets()
  job <- cardioversion_job(dir)
  events <- cardioversion_events(dir)
  expect_false(anyNA(events[c("maze_prc", "iso_pvi")]))

  # reltol is tightened on purpose. Under the default (1e-5) BFGS stops at LL
  # -242.2675 and reports convergence; R's log likelihood at the listing's own
  # estimates is -242.2541, so the likelihoods agree and the default optimizer
  # tolerance stopped 0.013 short on a flat ridge in the late shapes.
  ctl <- list(n_starts = 1L, reltol = 1e-12, maxit = 5000L)
  expect_no_warning(fit <- fit_translated(job, events, "fit_2", "status_2", ctl))
  expect_listing_rows(fit)
  # NOCONSERVE: "Conservation of events: Not invoked".
  expect_false(fit$spec$control$conserve_applied)

  # .lst lines 2760 and 2772-2793.
  sas <- sas_table(
    c("phase_1.log_mu", "phase_1.log_t_half", "phase_1.nu", "phase_1.m",
      "phase_1.MAZE_PRC", "phase_1.ISO_PVI",
      "phase_2.log_mu", "phase_2.log_tau", "phase_2.gamma", "phase_2.alpha", "phase_2.eta",
      "phase_2.MAZE_PRC", "phase_2.ISO_PVI"),
    logged = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
               FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    #      E0        E2        E3        E4        MAZE_PRC  ISO_PVI
    est = c(-1.02384, -4.56614, -1.77136, 1.770403, -0.54922, 0.1504561,
            #  L0        L1         L2        L3         L4        MAZE_PRC   ISO_PVI
            -1.16547, -0.692292, 2.014376, 0.7332091, -2.11007, -0.938336, -0.461191),
    se = c(0.1162237, 0.06331743, 0.2646996, 0.3068885, 0.1807944, 0.2205966,
           0.2346658, 0.3799748, 0.8113678, 0.2687002, 0.9540579, 0.1622002, 0.2525617)
  )
  expect_matches_listing(fit, sas, ll = -242.254)
})

test_that("hzr_translate_sas() on the cardioversion job builds EVENTS, then stops at line 65", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("haven")
  dir <- skip_if_no_maze_datasets()
  sas <- file.path(dirname(dir), "distributions", "hz.ce_cardioversion_repeated.ehb.sas")
  testthat::skip_if_not(file.exists(sas), "hz.ce_cardioversion_repeated.ehb.sas not on the volume")

  job <- suppressWarnings(hzr_translate_sas(sas))
  expect_equal(names(job$calls)[1:5], c("data", "repeated", "rewrite", "status", "fit"))

  bd_card <- .hzr_derive_bd_card(dir)
  names(bd_card) <- toupper(names(bd_card))
  env <- new.env(parent = globalenv())
  env$BD_CARD <- bd_card
  eval(job$calls$data, env)
  expect_no_warning(eval(job$calls[["repeated"]], env))
  # Shape first, then values: the .log's 962 rows, and 357 input columns plus
  # the function's eight.
  expect_equal(dim(env$EVENTS), c(962L, 365L))
  expect_equal(sum(env$EVENTS$CE_CARD == 1), 388L)
  # The job's line 65 is job code, not macro code; the document stops on it.
  expect_error(eval(job$calls$rewrite, env), "IV_EVENT=IV_EVENT+0.0001141553", fixed = TRUE)
})
