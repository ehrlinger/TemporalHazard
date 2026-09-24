# `$scope$force_in` and `$scope$force_out` record the caller's strings, and a
# name that pinned nothing is listed there as though it had been applied
# (#451). The names as given are worth keeping -- they are the record of what
# was asked for -- so the RESOLVED identities are recorded alongside them
# rather than replacing them.
#
# Four fields, four questions, none overlapping:
#   force_in             what did I ask for
#   force_in_resolved    what did those names resolve to
#   unresolved$force_in  which of them resolved to nothing
#   frozen               what did the loop hold in that I never named

srp_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
}

test_that("a resolved pin is recorded as its identity, not as written (#451)", {
  skip_on_cran() # a backward screen
  d <- srp_data()
  names(d)[names(d) == "age"] <- "_X1"
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ `_X1` + mal,
                                 data = d, dist = "weibull",
                                 theta = c(0.1, 1, 0, 0), fit = TRUE))
  sw <- suppressWarnings(hzr_stepwise(fit, data = d, direction = "backward",
                                      criterion = "wald", slstay = 1e-9,
                                      force_in = "_X1", trace = FALSE))
  # As given: exactly the caller's string.
  expect_identical(sw$scope$force_in, "_X1")
  # Resolved: the model's own label for that column.
  expect_identical(sw$scope$force_in_resolved, "`_X1`")
  # The two must DIFFER here, or the test cannot tell them apart.
  expect_false(identical(sw$scope$force_in, sw$scope$force_in_resolved))
  # Known positive: the pin actually held, so this is a live screen.
  expect_false("`_X1`" %in% sw$steps$variable[sw$steps$action == "drop"])
})

test_that("a pin that resolved to nothing is absent from the resolved field (#451)", {
  skip_on_cran() # a backward screen
  d <- srp_data()
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age + mal,
                                 data = d, dist = "weibull",
                                 theta = c(0.1, 1, 0, 0), fit = TRUE))
  sw <- suppressWarnings(hzr_stepwise(fit, data = d, direction = "backward",
                                      criterion = "wald", slstay = 1e-9,
                                      force_in = c("age", "nosuch"),
                                      trace = FALSE))
  expect_identical(sw$scope$force_in, c("age", "nosuch"))       # as given
  expect_identical(sw$scope$force_in_resolved, "age")           # what held
  expect_identical(sw$scope$unresolved$force_in, "nosuch")      # what did not
  # The three fields answer three different questions about one argument.
  expect_false("nosuch" %in% sw$scope$force_in_resolved)
})

test_that("the frozen set does not leak into the resolved pins (#451)", {
  skip_on_cran() # an oscillating two-way screen
  # THE test that distinguishes this design from the rejected one. The
  # internal `effective_force_in` is unique(c(force_in_id, frozen)), and
  # recording THAT would put variables into force_in_resolved that the
  # caller never named, duplicating $scope$frozen. Recording the resolved
  # identities keeps the two separate.
  d <- srp_data()
  base <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                  data = d, dist = "weibull",
                                  theta = c(0.1, 1), fit = TRUE))
  sw <- suppressWarnings(hzr_stepwise(base, data = d, scope = c("age", "mal"),
                                      direction = "both", criterion = "wald",
                                      slentry = 0.9, slstay = 0.5,
                                      max_move = 0L, trace = FALSE))
  # Known positive: something really did freeze, or this cannot fail.
  expect_gt(length(sw$scope$frozen), 0L)
  # The caller named no pins at all, so both fields must be empty even
  # though the loop froze variables.
  expect_identical(sw$scope$force_in, character())
  expect_identical(sw$scope$force_in_resolved, character())
  expect_false(any(sw$scope$frozen %in% sw$scope$force_in_resolved))
})

test_that("force_out is recorded both ways too (#451)", {
  skip_on_cran() # a forward screen
  d <- srp_data()
  names(d)[names(d) == "age"] <- "_X1"
  base <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                  data = d, dist = "weibull",
                                  theta = c(0.1, 1), fit = TRUE))
  sw <- suppressWarnings(hzr_stepwise(base, data = d, scope = c("_X1", "mal"),
                                      direction = "forward", criterion = "wald",
                                      slentry = 0.9, force_out = "_X1",
                                      trace = FALSE))
  expect_identical(sw$scope$force_out, "_X1")
  expect_identical(sw$scope$force_out_resolved, "`_X1`")
  # Known positive: the exclusion held.
  expect_false("`_X1`" %in% sw$steps$variable)
})
