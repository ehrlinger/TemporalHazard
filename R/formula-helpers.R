#' Parse Surv() formula for hazard modeling
#'
#' Extracts time, status, time_lower, time_upper, and predictors from a formula
#' of the form `Surv(time, status) ~ x1 + x2 + ...`.
#' Supports right-censored, left-censored, interval-censored, and
#' counting-process (start-stop) data.
#'
#' For counting-process (start-stop) data, use `Surv(start, stop, event)`.
#' The start times are returned as `time_lower` and stop times as `time`,
#' enabling the likelihood to compute `H(stop) - H(start)` per epoch.
#'
#' `Surv()` and this package code censoring status differently, so the
#' returned `status` is translated, not passed through:
#'
#' | Meaning   | TemporalHazard | `Surv` "left" | `Surv` "interval" |
#' | --------- | -------------- | ------------- | ----------------- |
#' | left      | `-1`           | `0`           | `2`               |
#' | right     | `0`            | --            | `0`               |
#' | event     | `1`            | `1`           | `1`               |
#' | interval  | `2`            | --            | `3`               |
#'
#' Under `type = "interval"`, `Surv()` reuses the `time2` column to hold the
#' status of any non-interval row, so an upper bound is read only where the
#' row is genuinely interval-censored.
#'
#' @param formula A formula object with Surv() on the LHS.
#' @param data A data frame containing variables referenced in the formula.
#' @return A list with elements: time, status, time_lower, time_upper, x
#' @keywords internal
.hzr_parse_formula <- function(formula, data) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }

  # Parse the LHS (should be Surv(...))
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  # Extract Surv() call
  if (!is.call(lhs) || !grepl("Surv$", deparse(lhs[[1L]]))) {
    stop("Formula LHS must be a Surv() call.", call. = FALSE)
  }

  # Evaluate the entire Surv() call in the context of data
  # Make sure Surv is available in the evaluation environment
  env_with_surv <- list2env(c(as.list(data), list(Surv = survival::Surv)))
  surv_obj <- eval(lhs, envir = env_with_surv)

  # Surv() returns a Surv object; extract the matrix and attributes
  if (!inherits(surv_obj, "Surv")) {
    stop("Formula LHS must return a Surv object.", call. = FALSE)
  }

  surv_type <- attr(surv_obj, "type")
  surv_mat <- unclass(surv_obj)

  if (surv_type == "right") {
    # Format: [time, status]
    time <- surv_mat[, 1L]
    status <- surv_mat[, 2L]
    time_lower <- NULL
    time_upper <- NULL
  } else if (surv_type == "left") {
    # Format: [time, status], where Surv codes 1 = event, 0 = left-censored.
    # TemporalHazard codes left-censoring as -1; passing 0 through would
    # silently read these rows as right-censored.
    time <- surv_mat[, 1L]
    status <- ifelse(surv_mat[, 2L] == 0, -1, 1)
    time_lower <- NULL
    time_upper <- surv_mat[, 1L]
  } else if (surv_type == "interval") {
    # Format: [time1, time2, status], where Surv codes
    #   0 = right-censored, 1 = event, 2 = left-censored, 3 = interval.
    # Map onto TemporalHazard's 0 / 1 / -1 / 2.
    surv_status <- surv_mat[, 3L]
    status <- c(0, 1, -1, 2)[surv_status + 1L]
    time <- surv_mat[, 1L]

    # Surv stores the status in `time2` for every non-interval row, so that
    # column is a sentinel except where surv_status == 3.
    is_interval <- surv_status == 3
    time_upper <- ifelse(is_interval, surv_mat[, 2L], time)

    # `time_lower` doubles as the counting-process entry time for status
    # 0 and 1, where the likelihood forms H(stop) - H(start).  Surv
    # "interval" carries no left truncation, so entry time is 0 outside the
    # interval rows; using `time` there would cancel those rows out.
    time_lower <- ifelse(is_interval, surv_mat[, 1L], 0)
  } else if (surv_type == "counting") {
    # Start-stop (counting process) format: Surv(start, stop, event)
    # Used for repeating events / epoch-decomposed longitudinal data.
    # Each epoch contributes H(stop) - H(start) to the likelihood.
    time_lower <- surv_mat[, 1L]  # entry (start) time
    time <- surv_mat[, 2L]        # exit (stop) time
    time_upper <- NULL
    status <- surv_mat[, 3L]
  } else {
    stop("Unsupported Surv() type: ", surv_type, call. = FALSE)
  }

  # Parse RHS (predictors)
  x <- NULL
  if (!is.null(rhs)) {
    # Reconstruct as a formula for model.matrix()
    rhs_formula <- formula(paste("~", deparse(rhs)))
    tryCatch({
      x <- stats::model.matrix(rhs_formula, data = data)
      # Remove intercept column if present
      if (ncol(x) > 0 && colnames(x)[1L] == "(Intercept)") {
        x <- x[, -1L, drop = FALSE]
      }
      if (ncol(x) == 0) {
        x <- NULL
      }
    }, error = function(e) {
      stop("Failed to parse formula RHS: ", e$message, call. = FALSE)
    })
  }

  list(
    time = time,
    status = status,
    time_lower = time_lower,
    time_upper = time_upper,
    x = x,
    surv_type = surv_type
  )
}
