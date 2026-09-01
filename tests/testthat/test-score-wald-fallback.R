# Wald fallback for candidates the score statistic cannot test (#130).
#
# The score criterion computes SAS HAZARD's Q exactly: Q = grad^2 * I22, the
# reciprocal Schur complement of the OBSERVED information at beta = 0
# (src/vars/q1.c).  When a candidate's true effect is far from zero the
# log-likelihood is convex there, the Schur complement turns negative, and the
# statistic is undefined -- so the criterion declines candidates in proportion
# to how predictive they are.
#
# SAS has the same defect.  q1.c documents it ("IT IS POSSIBLE THAT THE
# PROGRAM WILL RETURN A NEGATIVE Q VALUE ... THE USER SHOULD USE THE MORE
# EXPENSIVE Q2 AS AN ALTERNATIVE") and dqstat.c declines the candidate with
# p = 1.  Q2 is referenced once in the C tree and never implemented.
#
# So the statistic stays bit-faithful and only the *handling* diverges: a
# candidate the score cannot test is refit and tested by Wald, which is what
# the unbuilt Q2 was for.  These tests pin down both halves -- that the
# statistic is unchanged, and that the screen no longer misses strong
# candidates.

planted <- function(seed = 11, n = 400, beta = 0.9) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- stats::rnorm(n)
  data.frame(
    tt = stats::rexp(n, rate = 0.3 * exp(beta * x1)) + 0.01,
    ev = stats::rbinom(n, 1, 0.85),
    x1 = x1, x2 = x2, x3 = x3
  )
}

# Formula interface deliberately: hzr_stepwise() refuses a vector-interface
# base fit up front (.hzr_refit_blocker(), #159), because every candidate
# refit would fail.  Nothing about the Wald fallback depends on the
# interface -- the score criterion sees the same fitted MLE either way.
planted_fit <- function(D) {
  suppressWarnings(hazard(
    survival::Surv(tt, ev) ~ 1, data = D, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                  const = hzr_phase("constant")),
    fit = TRUE))
}

test_that("the score statistic itself is unchanged for an untestable candidate", {
  skip_on_cran()
  # Parity guard. The divergence from SAS is in the handling, not in Q, so
  # .hzr_score_q() must still decline x1 with the same reason as before.
  D <- planted()
  q <- .hzr_score_q(planted_fit(D), var = "x1", phase = "const", data = D)
  expect_true(is.na(q$stat))
  expect_identical(q$reason, "information_indefinite")
})

test_that("the planted strong effect is selected, not the noise variables", {
  skip_on_cran()
  # x1 has beta_true = 0.9 (LR ~ 178, p ~ 0); x2 and x3 are pure noise. Before
  # the fallback the screen entered x2 and x3 and never tested x1 -- a
  # confident, plausible, wrong selection.
  D <- planted()
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + x2 + x3),
    data = D, direction = "both", slentry = 0.05))

  steps <- as.data.frame(sw)
  entered <- steps$variable[toupper(steps$action) == "ENTER"]
  expect_true("x1" %in% entered)
  expect_false("x2" %in% entered)
  expect_false("x3" %in% entered)
  # The two assertions above do NOT test the fallback's narrowness, despite an
  # earlier comment here saying so: x2 and x3 have Wald p-values of 0.0997 and
  # 0.149 on this fixture, so a fallback that refit every candidate would
  # still see both declined by `slentry`. The count is what can fail --
  # refitting everything makes it 3.
  expect_identical(sw$criteria$n_wald_fallbacks, 1L)
})

test_that("the Wald fallback is reported rather than silent", {
  skip_on_cran()
  D <- planted()
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + x2 + x3),
    data = D, direction = "both", slentry = 0.05))
  # A criterion that quietly switched itself for some candidates would be the
  # same class of defect as the one being fixed.
  expect_true(sw$criteria$n_wald_fallbacks >= 1L)
})

test_that("only untestable-because-strong candidates fall back", {
  skip_on_cran()
  # A constant column is declined for a reason a refit cannot rescue, so it
  # must not cost one. Paying a refit per degenerate candidate would give back
  # the whole point of the score criterion.
  D <- planted()
  D$flat <- 1
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + flat),
    data = D, direction = "both", slentry = 0.05))
  entered <- as.data.frame(sw)$variable[toupper(as.data.frame(sw)$action) == "ENTER"]
  expect_false("flat" %in% entered)
  expect_true("x1" %in% entered)
  # The two assertions above stayed green when "constant" and "collinear" were
  # added to .hzr_score_fallback_reasons -- `flat` still does not enter, it
  # just costs a refit on the way out. Only the count sees the widening.
  expect_identical(sw$criteria$n_wald_fallbacks, 1L)
})

test_that("a Wald-decided row says so, so its p-value can be recomputed", {
  skip_on_cran()
  # THE REPORTING DEFECT. A rescued candidate carries a Wald z on a row still
  # labelled criterion = "score", where every neighbouring row carries a 1-df
  # chi-square Q. `df` cannot tell them apart -- both read 1 -- so a reader
  # recomputing pchisq(stat, df) on this fixture gets 1.1e-04 against a true
  # 2.5e-50, wrong by 46 orders of magnitude.
  D <- planted()
  sw <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1 + x2 + x3),
    data = D, direction = "both", slentry = 0.05))
  steps <- as.data.frame(sw)
  row <- steps[toupper(steps$action) == "ENTER" & steps$variable == "x1", ]
  expect_equal(nrow(row), 1L)
  expect_identical(row$criterion, "score")
  expect_identical(row$stat_type, "wald_z")
  # The point of the column: recomputing from `stat` under the distribution
  # `stat_type` names reproduces the reported p-value.
  #
  # Compared on the LOG scale, and guarded by p_value > 0 first. An
  # expect_equal(tolerance = 1e-8) on the natural scale cannot fail here:
  # all.equal switches from relative to absolute comparison once the target
  # falls below the tolerance, and at 2.5e-50 everything -- including a
  # p_value that underflowed to exactly 0 -- is within 1e-8 of everything
  # else. A hollow pass, in the test written to catch hollow passes.
  expect_gt(row$p_value, 0)
  expect_equal(log(row$p_value),
               stats::pnorm(-abs(row$stat), log.p = TRUE) + log(2),
               tolerance = 1e-8)
  # And the chi-square reading -- the one `df` alone invites -- does not.
  # Also on the log scale: the two differ by ~46 orders of magnitude, so the
  # gap is what should be asserted, not a difference that rounds to 1.1e-04
  # whatever the second term does.
  expect_gt(
    stats::pchisq(row$stat, row$df, lower.tail = FALSE, log.p = TRUE) -
      log(row$p_value),
    30
  )
})

test_that("fallback reasons are classified, not lumped together", {
  # Every reason eligible for a refit must be one .hzr_score_reason_text()
  # knows, or the warning would print a bare code.
  expect_false(any(
    .hzr_score_reason_text(.hzr_score_fallback_reasons) %in%
      .hzr_score_fallback_reasons
  ))
  # Degenerate causes must stay out of the fallback set: a refit cannot make a
  # constant or collinear column testable.
  expect_false("collinear" %in% .hzr_score_fallback_reasons)
  expect_false("constant" %in% .hzr_score_fallback_reasons)
  expect_false("non_numeric" %in% .hzr_score_fallback_reasons)
  expect_true("information_indefinite" %in% .hzr_score_fallback_reasons)
})

test_that("a fallback refit that fails is recorded and warned, not silent", {
  skip_on_cran()
  # Before this, a failed fallback refit hit `next` and recorded nothing: no
  # warning, no refit_failures token, no stopped_refit_failed. The resulting
  # row was byte-identical to one from a run where no refit was attempted at
  # all, so the only signal left -- an information_indefinite count -- meant
  # the opposite of what a reader would take it to mean.
  D <- planted()
  base <- planted_fit(D)

  testthat::local_mocked_bindings(
    .hzr_refit_with_scope = function(...) stop("forced refit failure")
  )

  w <- testthat::capture_warnings(
    sw <- hzr_stepwise(fit = base, scope = list(const = ~ x1 + x2 + x3),
                       data = D, direction = "both", slentry = 0.05)
  )

  # The candidate the score could not test is x1; its refit was attempted
  # and failed, so it must appear by name.
  expect_true(any(grepl("Wald-fallback refit failed", w)))
  expect_true(any(grepl("x1", sw$criteria$refit_failures, fixed = TRUE)))
  expect_gt(sw$criteria$n_refit_failures, 0L)
  # And nothing may be reported as a successful substitution.
  expect_identical(sw$criteria$n_wald_fallbacks, 0L)
})

test_that("an ordinary score row reads score_q, not wald_z", {
  skip_on_cran()
  # The score half of "stat_type varies". A weaker planted effect is scorable,
  # so it enters on Q rather than being rescued.
  #
  # The Wald half of this test used to live here, running the same fixture
  # under criterion = "wald". It failed in CI on Linux and Windows with zero
  # ENTER rows -- not because of the threshold, but because the candidate
  # REFIT failed there ("1 candidate refit FAILED and could not be tested --
  # x1@const"), which the Wald path needs and the score path does not. A
  # multiphase refit is not a stable thing to hang a reporting assertion on,
  # so the Wald side is pinned in test-candidate-score.R against the scorer
  # itself, where no refit is involved.
  D <- planted(beta = 0.35)
  sw_s <- suppressWarnings(hzr_stepwise(
    fit = planted_fit(D), scope = list(const = ~ x1), data = D,
    direction = "forward", criterion = "score", slentry = 0.2))
  st_s <- as.data.frame(sw_s)
  st_s <- st_s[toupper(st_s$action) == "ENTER", ]
  expect_gte(nrow(st_s), 1L)
  expect_identical(unique(st_s$stat_type), "score_q")
  expect_identical(unique(st_s$df), 1L)
  # And this arm really was scored, not rescued -- otherwise "score_q" could
  # quietly become "wald_z" and the assertion above would be vacuous.
  expect_identical(sw_s$criteria$n_wald_fallbacks, 0L)
  # Q is a 1-df chi-square, and reads as one.
  expect_equal(stats::pchisq(st_s$stat[1], 1, lower.tail = FALSE),
               st_s$p_value[1], tolerance = 1e-8)
})

test_that("the bootstrap aggregates Wald fallbacks instead of dropping them", {
  skip_on_cran()
  # Each replicate runs under suppressWarnings(), so a select-mode run in
  # which the fallback fired everywhere reported nothing at all -- and the
  # bootstrap is where a wholesale substitution matters most, since these
  # entries drive the pooled selection frequencies.
  D <- planted()
  # No `data =`: hzr_bootstrap() supplies the resampled frame itself, and
  # passing one collides with it through `...`.
  boot <- suppressWarnings(hzr_bootstrap(
    planted_fit(D), scope = list(const = ~ x1 + x2 + x3),
    n_boot = 3L, direction = "both", criterion = "score",
    slentry = 0.05, verbose = FALSE))
  expect_true(is.numeric(boot$n_wald_fallbacks))
  expect_gte(boot$n_wald_fallback_replicates, 1L)
  # The total counts entries, not replicates, so it can only be the larger.
  expect_gte(boot$n_wald_fallbacks, boot$n_wald_fallback_replicates)
})
