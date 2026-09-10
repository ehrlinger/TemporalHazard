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
  # its input warnings.
  expect_no_warning(events <- hzr_repeated_events(bd_card, id, tm, fu, ind))
  expect_equal(events, `row.names<-`(s7, NULL))

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
  run <- function(d) hzr_repeated_events(d, "ccfid", "iv_event", "iv_end", "ce_card")
  base <- canon(run(bd_card))
  for (k in 1:5) {
    shuffled <- withr::with_seed(k, bd_card[sample(nrow(bd_card)), , drop = FALSE])
    expect_equal(canon(run(shuffled)), base, label = paste("shuffle", k))
  }
})
