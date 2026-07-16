# Unit tests for the bh.dead .lst parsers.
#
# Inputs below reproduce the SAS listing FORMAT with entirely INVENTED numbers.
# Never paste real listing values here: this repository is public, and a
# renamed variable with verbatim values is not anonymised. Round, obviously
# fake numbers make that impossible to get wrong by accident.

# They run from a source checkout (devtools::load_all() / devtools::test()) and
# SKIP against an installed package, because the parsers they exercise live in
# inst/dev/, which is .Rbuildignore'd and therefore absent from the tarball.
# They never touch the secure volume and contain no study values.

.bhdead_helper <- system.file("dev", "bhdead-parity", "parse-bhdead-lst.R",
                              package = "TemporalHazard")
if (nzchar(.bhdead_helper)) {
  source(.bhdead_helper, local = TRUE)
}

.fake_bh_lines <- c(
  "                        Study Title Redacted                       3",
  "                                                  Early Phase",
  "",
  "                           Obs    _NAME_         N     PCT          MIN        MAX        MEAN        STD",
  "",
  "                             1    E0           500    100.0    -100.000    100.000     -1.0000     10.000",
  "                             2    VAR_A        250     50.0      -2.000      2.000      0.2000      2.000",
  "                             3    VAR_B        100     20.0      -3.000      3.000      0.3000      3.000",
  "                        Study Title Redacted                       4",
  "                                                  Constant Phase",
  "",
  "                           Obs    _NAME_         N     PCT          MIN        MAX        MEAN        STD",
  "",
  "                             1    C0           500    100.0    -200.000    200.000     -2.0000     20.000",
  "                             2    VAR_A        150     30.0      -4.000      4.000      0.4000      4.000"
)

.fake_hz_lines <- c(
  "          Phase        Parameter       Specified        Used       |  Theta",
  "                       THALF             1.100000      1.100000    |  E2",
  "                       NU                 -0.5000       -0.5000    |  E3",
  "                       M                        0             0    |  E4",
  "                       MUE                 0.8000        0.8000    |  E0",
  "          Constant:    MUC               0.050000      0.050000    |  C0",
  "                                          Estimates for Model Parameters",
  "                                          Phase       Parameter     Fixed?     Estimate",
  "                                                      THALF         No            1.200000",
  "                                                      NU            No           -0.600000",
  "                                                      M             Yes                  0",
  "                                                      MUE                        0.850000",
  "                                          Constant:   MUC                       0.055000"
)

test_that(".hzr_parse_bhdead_phase extracts the Early table", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  d <- .hzr_parse_bhdead_phase(.fake_bh_lines, "Early Phase")
  expect_identical(d$name, c("E0", "VAR_A", "VAR_B"))
  expect_identical(d$n, c(500L, 250L, 100L))
  expect_equal(d$pct, c(100.0, 50.0, 20.0))
  expect_equal(d$mean[2], 0.2)
  expect_equal(d$sd[3], 3.0)
})

test_that(".hzr_parse_bhdead_phase does not bleed across phase boundaries", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  d <- .hzr_parse_bhdead_phase(.fake_bh_lines, "Constant Phase")
  expect_identical(d$name, c("C0", "VAR_A"))
  expect_equal(d$pct, c(100.0, 30.0))
})

test_that(".hzr_parse_bhdead_shape reads specified and converged values", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  s <- .hzr_parse_bhdead_shape(.fake_hz_lines)
  expect_equal(s$specified[["thalf"]], 1.1)
  expect_equal(s$specified[["mue"]], 0.8)
  expect_equal(s$converged[["thalf"]], 1.2)
  expect_equal(s$converged[["nu"]], -0.6)
  expect_equal(s$converged[["m"]], 0)
  expect_equal(s$converged[["muc"]], 0.055)
})

test_that(".hzr_parse_bhdead_phase does not absorb a trailing unrelated table", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  # The terminal phase table has no following label to bound it. A trailing
  # table that happens to match the data-row shape must not bleed in -- its
  # Obs sequence restarts at 1 rather than continuing the run.
  lines <- c(
    "Constant Phase", "", "Obs _NAME_ N PCT MIN MAX MEAN STD", "",
    "  1    C0    500  100.0   -1.000   1.000   0.1000   1.000",
    "", "Some Other Unrelated Table",
    "Obs _NAME_ N PCT MIN MAX MEAN STD", "",
    "  1    XYZ_9   500  100.0   -9.000   9.000   9.0000   9.000"
  )
  d <- .hzr_parse_bhdead_phase(lines, "Constant Phase")
  expect_identical(d$name, "C0")
})

test_that(".hzr_parse_bhdead_phase spans a paginated table via a contiguous Obs run", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  # The phase label and column header repeat on each page, but the Obs
  # sequence continues uninterrupted across the page break.
  lines <- c(
    "Early Phase", "", "Obs _NAME_ N PCT MIN MAX MEAN STD", "",
    "  1    E0      500  100.0   -1.000   1.000   0.1000   1.000",
    "  2    VAR_A   250   50.0   -2.000   2.000   0.2000   2.000",
    "Study Title Redacted    2",
    "Early Phase", "", "Obs _NAME_ N PCT MIN MAX MEAN STD", "",
    "  3    VAR_B   100   20.0   -3.000   3.000   0.3000   3.000",
    "  4    VAR_C    50   10.0   -4.000   4.000   0.4000   4.000"
  )
  d <- .hzr_parse_bhdead_phase(lines, "Early Phase")
  expect_identical(d$name, c("E0", "VAR_A", "VAR_B", "VAR_C"))
})

test_that(".hzr_build_bhdead_fixture returns a schema-valid fixture", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  hz <- tempfile(fileext = ".lst"); bh <- tempfile(fileext = ".lst")
  on.exit(unlink(c(hz, bh)), add = TRUE)
  writeLines(.fake_hz_lines, hz); writeLines(.fake_bh_lines, bh)
  fix <- .hzr_build_bhdead_fixture(
    hz, bh,
    meta = list(resampl = 1000L, sle = 0.12, sls = 0.10, n_obs = 100L,
                captured_on = "2026-07-16", source = "synthetic")
  )
  expect_identical(.hzr_validate_bhdead_fixture(fix), fix)
  expect_identical(fix$early$name[1], "E0")
})
