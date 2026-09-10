# The record of what a fit did not do (#242, following #197): builder,
# validator, formatter and output check. Every expected cause string below is
# copied from inst/dev/DEGRADED-RECORD-DESIGN.md, never from the package's own
# constants -- a test that reads the table it is checking cannot fail.

vc <- diag(2)

fake_fit <- function(degraded = character(0), causes = character(0),
                     vcov = vc, weak = NULL, dist = "weibull",
                     control = list()) {
  structure(list(spec = list(dist = dist, control = control),
                 fit = list(vcov = vcov, weak = weak),
                 degraded = degraded, degraded_causes = causes),
            class = "hazard")
}

record <- function(...) {
  .hzr_degraded_record(control = list(), dist = "weibull", ...)
}

test_that("a clean single-distribution fit records nothing", {
  rec <- record(vcov = vc, weak = NULL)
  expect_identical(rec$degraded, character(0))
  expect_length(rec$degraded_causes, 0L)
})

test_that("each capability is recorded with the cause it was given", {
  cases <- list(
    list(args = list(vcov = NULL, weak = NA, fitted = FALSE),
         want = c(fitting = "not requested (fit = FALSE)",
                  standard_errors = "model not fitted",
                  weak_direction_check = "model not fitted")),
    list(args = list(vcov = NULL, weak = NA, fitted = FALSE,
                     not_fitted_cause = paste0("no starting values (theta = NULL); ",
                                               "the optimizer did not run")),
         want = c(fitting = "no starting values (theta = NULL); the optimizer did not run",
                  standard_errors = "model not fitted",
                  weak_direction_check = "model not fitted")),
    list(args = list(vcov = vc, weak = NA, imported = TRUE),
         want = c(fitting = "imported from SAS output; not fitted in R",
                  weak_direction_check = "imported from SAS output; no R Hessian")),
    list(args = list(vcov = NULL, weak = NA, imported = TRUE),
         want = c(fitting = "imported from SAS output; not fitted in R",
                  standard_errors = "covariance not imported",
                  weak_direction_check = "imported from SAS output; no R Hessian")),
    list(args = list(vcov = NA, weak = NA,
                     reasons = list(se = "Hessian not invertible",
                                    weak = "standard errors unavailable")),
         want = c(standard_errors = "Hessian not invertible",
                  weak_direction_check = "standard errors unavailable")),
    # State forces the entries; a missing reason must not remove them.
    list(args = list(vcov = NA, weak = NA),
         want = c(standard_errors = "cause not recorded",
                  weak_direction_check = "cause not recorded"))
  )
  for (i in seq_along(cases)) {
    rec <- do.call(record, cases[[i]]$args)
    expect_identical(rec$degraded, names(cases[[i]]$want), info = i)
    expect_identical(rec$degraded_causes, cases[[i]]$want, info = i)
  }
})

test_that("every conserve_disabled_reason becomes the spec's cause", {
  want <- c(
    not_requested         = "not requested (conserve = FALSE)",
    unsupported_censoring = "status outside {0, 1} (left or interval censoring)",
    single_phase          = "fewer than two phases",
    no_events             = "no events",
    setup_failed          = "setup did not complete"
  )
  for (r in names(want)) {
    rec <- .hzr_degraded_record(
      vcov = vc, weak = NULL, dist = "multiphase",
      control = list(conserve_applied = FALSE, conserve_disabled_reason = r))
    expect_identical(rec$degraded_causes,
                     c(conservation_of_events = want[[r]]), info = r)
  }
  # A reason the optimizer adds later still produces an entry.
  rec <- .hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = FALSE, conserve_disabled_reason = "new"))
  expect_identical(rec$degraded_causes[["conservation_of_events"]],
                   "cause not recorded")
  # Applied CoE, and CoE on a non-multiphase fit, record nothing.
  expect_length(.hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = TRUE))$degraded, 0L)
  expect_length(.hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "weibull",
    control = list(conserve_applied = FALSE))$degraded, 0L)
})

test_that("the conserved phase's variance is recorded only when CoE ran", {
  rec <- .hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = TRUE),
    reasons = list(conserved_variance = "numDeriv not installed"))
  expect_identical(rec$degraded_causes,
                   c(conserved_phase_variance = "numDeriv not installed"))
  rec <- .hzr_degraded_record(
    vcov = vc, weak = NULL, dist = "multiphase",
    control = list(conserve_applied = TRUE),
    reasons = list(conserved_variance = NA_character_))
  expect_length(rec$degraded, 0L)
})

test_that("entries come out in canonical order", {
  rec <- .hzr_degraded_record(
    vcov = NA, weak = NA, dist = "multiphase",
    control = list(conserve_applied = FALSE,
                   conserve_disabled_reason = "unsupported_censoring"),
    reasons = list(weak = "standard errors unavailable",
                   se = "Hessian not invertible"))
  expect_identical(rec$degraded, c("standard_errors", "weak_direction_check",
                                   "conservation_of_events"))
})

test_that("the validator accepts a record that matches the object", {
  expect_true(.hzr_validate_degraded(fake_fit(), fitted = TRUE))
  ok <- fake_fit(c("standard_errors", "weak_direction_check"),
                 c(standard_errors = "x", weak_direction_check = "y"),
                 vcov = NA, weak = NA)
  expect_true(.hzr_validate_degraded(ok, fitted = TRUE))
})

test_that("the validator rejects each way the record can go wrong", {
  bad <- list(
    list(fake_fit(c("standard_errors", "standard_errors"),
                  c(standard_errors = "x", standard_errors = "x"), vcov = NA),
         "duplicate"),
    list(fake_fit("coffee", c(coffee = "x")), "unknown capability"),
    list(fake_fit(c("weak_direction_check", "standard_errors"),
                  c(weak_direction_check = "y", standard_errors = "x"),
                  vcov = NA, weak = NA),
         "canonical order"),
    list(fake_fit("standard_errors", c(weak_direction_check = "x"), vcov = NA),
         "names\\(degraded_causes\\)"),
    list(fake_fit("standard_errors", c(standard_errors = ""), vcov = NA),
         "empty cause"),
    list(fake_fit(vcov = NA), "'standard_errors' is absent"),
    list(fake_fit("standard_errors", c(standard_errors = "x")),
         "'standard_errors' is listed"),
    list(fake_fit(weak = NA), "'weak_direction_check' is absent"),
    list(fake_fit(dist = "multiphase", control = list(conserve_applied = FALSE)),
         "'conservation_of_events' is absent"),
    list(fake_fit("conserved_phase_variance", c(conserved_phase_variance = "x")),
         "was not applied")
  )
  for (b in bad) {
    expect_error(.hzr_validate_degraded(b[[1]], fitted = TRUE), b[[2]],
                 info = b[[2]])
  }
  expect_error(.hzr_validate_degraded(fake_fit(), fitted = FALSE),
               "'fitting' is absent")
  expect_error(.hzr_validate_degraded(fake_fit(), fitted = TRUE, imported = TRUE),
               "'fitting' is absent")
})

test_that("the formatter prints none, the entries, or not recorded", {
  expect_identical(.hzr_format_not_done(character(0), character(0)),
                   "  Not done in this run: none")
  expect_identical(
    .hzr_format_not_done(c("standard_errors", "weak_direction_check"),
                         c(standard_errors = "a", weak_direction_check = "b")),
    c("  Not done in this run:", "    standard_errors: a",
      "    weak_direction_check: b"))
  expect_identical(
    .hzr_format_not_done(NULL, NULL),
    "  Not done in this run: not recorded (object built before this record existed)")
})

test_that("the output check passes a faithful block and names each fault", {
  obj <- fake_fit(c("standard_errors", "weak_direction_check"),
                  c(standard_errors = "a", weak_direction_check = "b"),
                  vcov = NA, weak = NA)
  good <- c("hazard object", "  Not done in this run:",
            "    standard_errors: a", "    weak_direction_check: b",
            "  other: 1")
  chk <- function(lines, o = obj) .hzr_check_not_done_output(lines, o)

  expect_identical(chk(good), character(0))
  expect_identical(chk(good[-(2:4)]), "heading missing")
  expect_identical(chk(c(good, "  Not done in this run: none")),
                   "heading printed more than once")
  expect_true("entry missing: weak_direction_check" %in% chk(good[-4]))
  expect_true("'none' printed over a non-empty record" %in%
                chk("  Not done in this run: none"))
  expect_true("cause differs for standard_errors" %in%
                chk(sub(": a$", ": z", good)))
  expect_true("entry not in the record: fitting" %in%
                chk(append(good, "    fitting: q", after = 2)))
  expect_identical(chk(good[c(1, 2, 4, 3, 5)]), "entries out of order")

  expect_identical(chk("  Not done in this run: none", fake_fit()), character(0))

  legacy <- fake_fit()
  legacy$degraded <- NULL
  legacy$degraded_causes <- NULL
  expect_identical(
    chk("  Not done in this run: not recorded (object built before this record existed)",
        legacy),
    character(0))
  expect_identical(chk("  Not done in this run: none", legacy),
                   "an object with no record was not reported as unrecorded")
})
