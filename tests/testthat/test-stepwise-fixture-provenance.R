# Provenance of the AVC stepwise SAS capture (#540).
#
# The capture fixes all three early shapes (THALF, NU and M). The 64-bit (LP64)
# builds of the HAZARD C binary compute the wrong likelihood for exactly that
# case (ehrlinger/hazard#47): on this job today's affected binary reports
# -201.096 for the base model, where R and the capture both give -211.545. So
# a capture from an affected build carries log-likelihoods R cannot reproduce.
#
# This test pins R's log-likelihood at the capture's own printed final
# estimates to the capture's final logLik. Printed to three decimals, the two
# agree to 0.001; the affected build is off by about 10. The whole-model
# tolerance of 10 in test-stepwise-parity.R compares two fits and cannot
# catch the difference.

test_that("R reproduces the AVC stepwise capture's logLik at its estimates (#540)", {
  fix <- .hzr_load_stepwise_fixture("avc-forward-wald")
  if (is.null(fix)) {
    skip("stepwise-avc-forward-wald.rds fixture not found")
  }
  data(avc, envir = environment())
  avc <- na.omit(avc)

  # The capture's fixed shapes come from its PARMS statement
  # (inst/extdata/stepwise-fixtures/payload/avc-forward-wald.sas), which the
  # fixture does not record.
  shapes <- c(t_half = 0.1512095, nu = 1.438652, m = 1)
  est <- fix$final$coef
  early <- setdiff(est$variable[est$phase == "early"], "e0")
  constant <- setdiff(est$variable[est$phase == "constant"], "c0")
  phases <- list(
    early = hzr_phase("cdf", t_half = shapes[["t_half"]],
                      nu = shapes[["nu"]], m = shapes[["m"]],
                      fixed = c("t_half", "nu", "m"),
                      formula = reformulate(early)),
    constant = hzr_phase("constant", formula = reformulate(constant))
  )
  spec <- hazard(survival::Surv(int_dead, dead) ~ 1, data = avc,
                 dist = "multiphase", phases = phases, fit = FALSE)

  nm <- hzr_theta_names(phases, covariates = list(early = early,
                                                  constant = constant))
  key <- paste(est$phase, ifelse(est$variable %in% c("e0", "c0"), "log_mu",
                                 est$variable), sep = ".")
  theta <- stats::setNames(rep(NA_real_, length(nm)), nm)
  theta[c("early.log_t_half", "early.nu", "early.m")] <-
    c(log(shapes[["t_half"]]), shapes[["nu"]], shapes[["m"]])
  # Every printed estimate must land on a parameter, and every parameter
  # must be filled: a partial match would evaluate somewhere else. The set
  # comparison ignores repeats, so a duplicated row is refused first: the
  # assignment below would keep only its last value (Copilot, #561).
  expect_identical(anyDuplicated(key), 0L)
  expect_setequal(key, setdiff(nm, c("early.log_t_half", "early.nu",
                                     "early.m")))
  theta[key] <- est$estimate
  expect_false(anyNA(theta))

  ll <- as.numeric(hzr_evaluate(spec, theta)$logLik)
  expect_equal(fix$final$logLik, -182.659)
  expect_lt(abs(ll - fix$final$logLik), 0.01)
})
