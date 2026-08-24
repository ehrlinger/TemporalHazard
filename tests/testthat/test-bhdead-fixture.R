# Unit tests for the bh.dead parity fixture schema/loader.
# These construct synthetic fixtures in-memory: they never touch the secure
# volume and contain no study values.
#
# They run from a source checkout (devtools::load_all() / devtools::test()) and
# SKIP against an installed package, because the helper they exercise lives in
# inst/dev/, which is .Rbuildignore'd and therefore absent from the tarball.

.bhdead_helper <- system.file("dev", "bhdead-parity", "bhdead-fixture.R",
                              package = "TemporalHazard")
if (nzchar(.bhdead_helper)) {
  source(.bhdead_helper, local = TRUE)
}

.fake_phase <- function(names_vec) {
  data.frame(
    name = names_vec,
    n    = rep(10L, length(names_vec)),
    pct  = rep(50, length(names_vec)),
    min  = rep(-1, length(names_vec)),
    max  = rep(1, length(names_vec)),
    mean = rep(0, length(names_vec)),
    sd   = rep(1, length(names_vec)),
    stringsAsFactors = FALSE
  )
}

.fake_fixture <- function() {
  list(
    shape = list(
      specified = c(thalf = 1, nu = -1, m = 0, mue = 1, muc = 1),
      converged = c(thalf = 1, nu = -1, m = 0, mue = 1, muc = 1)
    ),
    early    = .fake_phase(c("E0", "VAR_A", "VAR_B")),
    constant = .fake_phase(c("C0", "VAR_A", "VAR_B")),
    meta = list(resampl = 1000L, sle = 0.12, sls = 0.10, n_obs = 100L,
                captured_on = "2026-07-16", source = "synthetic")
  )
}

test_that(".hzr_validate_bhdead_fixture accepts a well-formed fixture", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  fix <- .fake_fixture()
  expect_identical(.hzr_validate_bhdead_fixture(fix), fix)
})

test_that(".hzr_validate_bhdead_fixture reports every missing field at once", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  fix <- .fake_fixture()
  fix$meta$sle <- NULL
  fix$shape$converged <- NULL
  expect_error(.hzr_validate_bhdead_fixture(fix), "sle")
  expect_error(.hzr_validate_bhdead_fixture(fix), "converged")
})

test_that(".hzr_validate_bhdead_fixture rejects a phase missing a column", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  fix <- .fake_fixture()
  fix$early$pct <- NULL
  expect_error(.hzr_validate_bhdead_fixture(fix), "early")
})

test_that(".hzr_bhdead_candidates drops the intercept row and lowercases", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  expect_identical(.hzr_bhdead_candidates(.fake_fixture()), c("var_a", "var_b"))
})

test_that(".hzr_load_bhdead_fixture returns NULL when the file is absent", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  expect_null(.hzr_load_bhdead_fixture(file.path(tempdir(), "nope.rds")))
})

test_that(".hzr_load_bhdead_fixture round-trips a valid fixture", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  p <- file.path(tempdir(), "bhdead-test.rds")
  on.exit(unlink(p), add = TRUE)
  saveRDS(.fake_fixture(), p)
  expect_identical(.hzr_load_bhdead_fixture(p), .fake_fixture())
})

test_that("path helpers honour their environment variables", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  withr::with_envvar(c(TEMPORALHAZARD_BHBLT = "/tmp/x.sas7bdat"), {
    expect_identical(.hzr_bhblt_path(), "/tmp/x.sas7bdat")
  })
  withr::with_envvar(c(TEMPORALHAZARD_BHDEAD_FIXTURE = "/tmp/y.rds"), {
    expect_identical(.hzr_bhdead_fixture_path(), "/tmp/y.rds")
  })
})

test_that("path helpers return \"\" when their environment variable is unset", {
  testthat::skip_if_not(nzchar(.bhdead_helper),
                        "bh.dead parity helpers not installed (inst/dev is .Rbuildignore'd)")
  withr::with_envvar(c(TEMPORALHAZARD_BHBLT = NA), {
    expect_identical(.hzr_bhblt_path(), "")
  })
  withr::with_envvar(c(TEMPORALHAZARD_BHDEAD_FIXTURE = NA), {
    expect_identical(.hzr_bhdead_fixture_path(), "")
  })
})
