#' Read a SAS `outhaz` estimate dataset
#'
#' `PROC HAZARD`'s `outhaz=` dataset stores the converged estimates and the
#' asymptotic variance-covariance matrix at full double precision, where the
#' printed `.lst` carries about seven significant figures. For any quantity the
#' dataset holds, it is the better parity reference: print precision stops
#' being the binding constraint and optimizer convergence takes over.
#'
#' The log-likelihood is *not* stored here; take it from the `.lst`.
#'
#' @param path Path to a `.sas7bdat` written by `outhaz=`, or to an `.rds`
#'   holding the same data frame.
#' @return A plain list -- **not** a `hazard` object -- with `estimates`
#'   (named numeric), `status` (named integer, 1 free / 0 fixed), `vcov`
#'   (matrix over free parameters, dimnames set) and `flags` (named numeric of
#'   model-structure flags). When no parameter is free, `vcov` is a 0x0 matrix
#'   rather than `NULL`, so check its dimensions rather than `is.null()`.
#'
#' @section Experimental:
#' This function is experimental and its return shape is expected to change.
#' Because the result is an unclassed list, it has no `predict()` method: the
#' `hzr_translate_sas(librefs = )` path emits `predict()` against it, which
#' does not work today (#151). The `_STATUS_` coding is asserted against a
#' synthetic fixture, so a real `OUTHAZ=` file using a different convention
#' would yield an empty `vcov` alongside a fully populated `estimates`.
#' @export
#' @examples
#' f <- system.file("extdata", "outhaz-fixture.rds", package = "TemporalHazard")
#' if (nzchar(f)) str(hzr_read_outhaz(f))
hzr_read_outhaz <- function(path) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    d <- readRDS(path)
  } else {
    if (!requireNamespace("haven", quietly = TRUE)) {
      stop("hzr_read_outhaz() needs the 'haven' package to read a .sas7bdat.",
           call. = FALSE)
    }
    d <- as.data.frame(haven::read_sas(path))
  }

  if (!all(c("_NAME_", "_EST_", "_STATUS_") %in% names(d))) {
    stop("hzr_read_outhaz(): '", path, "' is not an outhaz dataset -- it must ",
         "carry columns _NAME_, _EST_ and _STATUS_. Found: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  }

  nm <- as.character(d[["_NAME_"]])
  est <- stats::setNames(as.numeric(d[["_EST_"]]), nm)
  st_raw <- suppressWarnings(as.integer(d[["_STATUS_"]]))

  # Rows with NA status are model-structure flags (G1FLAG, FIXDEL0, ...),
  # not parameters. They carry no estimate row in the vcov block.
  is_param <- !is.na(st_raw)
  status <- stats::setNames(st_raw[is_param], nm[is_param])
  flags <- est[!is_param]

  free <- names(status)[status == 1L]
  vcov <- as.matrix(d[match(free, nm), free, drop = FALSE])
  dimnames(vcov) <- list(free, free)
  storage.mode(vcov) <- "double"

  list(estimates = est[is_param], status = status, vcov = vcov, flags = flags)
}
