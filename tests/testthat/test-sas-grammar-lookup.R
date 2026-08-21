test_that("context decides which token a colliding spelling resolves to", {
  expect_equal(.hzr_sas_token("M", "HAZARD", "PARM"), "M")
  expect_equal(.hzr_sas_token("M", "HAZARD", "PHOP"), "MOVE")
  expect_equal(.hzr_sas_token("NOS", "HAZPRED", "HZPP"), "NOSURV")
})

test_that("aliases resolve to their canonical token", {
  expect_equal(.hzr_sas_token("MI", "HAZARD", "HZRP"), "MAXITER")
  expect_equal(.hzr_sas_token("QUASI", "HAZARD", "HZRP"), "QUASINEWTON")
  expect_equal(.hzr_sas_token("PARMS", "HAZARD", "STMT"), "PARAMETERS")
})

test_that("an unknown keyword returns NA rather than guessing", {
  expect_true(is.na(.hzr_sas_token("NOTAKEYWORD", "HAZARD", "HZRP")))
  # Right spelling, wrong context is still unknown.
  expect_true(is.na(.hzr_sas_token("SLENTRY", "HAZARD", "HZRP")))
})
