# The fixture is synthetic: real SAS `outhaz` layout, invented numbers. See
# data-raw/outhaz_fixture.R for the layout it reproduces and why the values
# are not from a real fit.

fixture <- function() {
  system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
}

test_that("hzr_read_outhaz returns the four documented elements", {
  out <- hzr_read_outhaz(fixture())
  expect_named(out, c("estimates", "status", "vcov", "flags"))
})

# THE LITERALS BELOW ARE 18 SIGNIFICANT DIGITS, NOT 17. IEEE754 guarantees a
# double round-trips through 17, but R's parser does not deliver that here: it
# reads 17-digit C0 one ULP low, and `tolerance = 0` sees the difference. Any
# literal added to these tests must be checked to round-trip -- write it with
# sprintf("%.18g") and confirm identical(as.numeric(s), x) -- or the test
# asserts a number the fixture does not contain.
test_that("estimates carry full double precision, not the .lst's 7 figures", {
  out <- hzr_read_outhaz(fixture())
  expect_equal(out$estimates[["THALF"]], 0.031415926535897934, tolerance = 0)
  expect_equal(out$estimates[["C0"]], -3.32192809488736263, tolerance = 0)
})

test_that("status marks free parameters 1 and fixed parameters 0", {
  out <- hzr_read_outhaz(fixture())
  expect_equal(unname(out$status[["THALF"]]), 1L)
  expect_equal(unname(out$status[["M"]]), 0L)
})

test_that("vcov covers only free parameters and is symmetric", {
  out <- hzr_read_outhaz(fixture())
  expect_setequal(rownames(out$vcov), c("THALF", "NU", "E0", "C0"))
  expect_equal(rownames(out$vcov), colnames(out$vcov))
  # SAS writes the two triangles independently, so the block is symmetric to
  # about 1e-18 and not bit-identical. Compare, do not identical().
  expect_equal(out$vcov, t(out$vcov))
  expect_equal(out$vcov[["THALF", "THALF"]], 0.101321183642337789,
               tolerance = 0)
})

test_that("flags are separated from parameters, not returned as estimates", {
  out <- hzr_read_outhaz(fixture())
  expect_equal(unname(out$flags[["G1FLAG"]]), 2)
  expect_false("G1FLAG" %in% names(out$estimates))
  expect_false("G1FLAG" %in% names(out$status))
})

test_that("a file that is not an outhaz dataset errors, naming the columns", {
  tmp <- withr::local_tempfile(fileext = ".rds")
  saveRDS(data.frame(a = 1), tmp)
  expect_error(hzr_read_outhaz(tmp), "_NAME_")
})
