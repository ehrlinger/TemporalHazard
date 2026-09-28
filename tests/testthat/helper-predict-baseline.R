# predict() at every covariate 0: the baseline a structural reference is
# built from, H(t | x) = exp(x beta) H0(t). A time-only newdata gives it, and
# since #522 predict() warns that it did so. Tests that want the baseline call
# this, which asserts that warning was given exactly once and passes any
# other warning through.
predict_baseline <- function(...) {
  n <- 0L
  value <- withCallingHandlers(
    predict(...),
    hzr_predict_covariates_zero = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  testthat::expect_identical(n, 1L)
  value
}
