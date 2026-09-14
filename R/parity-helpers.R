#' @importFrom utils write.csv read.table
#' @keywords internal
NULL

# parity-helpers.R -- Helpers for cross-validating R output against the legacy C binary
#
# STATUS: Stub / scaffold (M1 infrastructure, awaiting C binary integration)
#
# PURPOSE
# -------
# These functions wrap the compiled HAZARD C binary so the test suite can invoke
# it programmatically, capture its output, and compare parameter estimates with
# those produced by the R implementation.  This cross-validation is the primary
# correctness signal during the SAS->R migration.
#
# CURRENT STATE
# -------------
# .hzr_get_hazard_binary() -- locates binary via env var / option / package fallback
# .hzr_run_hazard_binary() -- writes working files, builds a .sas control file,
#                             and invokes the legacy executable as:
#                             hazard.exe < prefix.sas > prefix.lst 2>&1
# .hzr_parse_hazard_output() -- parses parameter tables from the .lst output.
#
# ENABLING C PARITY TESTS
# -----------------------
# 1. Prefer external path via TEMPORAL_HAZARD_BIN or options(temporal_hazard.binary=...)
#    for local acceptance testing.
# 2. Optionally place the compiled binary at inst/bin/hazard (chmod +x) for local use.
# 3. Complete .hzr_run_hazard_binary(): actual system() call + output parser.
# 4. Update test-parity-core.R to call .hzr_run_hazard_binary() instead of
#    loading .rds fixtures where true C-binary comparison is desired.
#
# Until the binary is wired up, parity tests use synthetic R-generated fixtures
# produced by golden_fixtures.R.

# Parity Testing Helpers
#
# Functions to generate golden fixtures from the C reference hazard binary
# and validate R implementation against them.

#' Get path to compiled hazard binary
#'
#' @noRd
.hzr_get_hazard_binary <- function(binary = c("hazard", "hazpred")) {
  binary <- match.arg(binary)
  env_var <- if (identical(binary, "hazpred")) {
    "TEMPORAL_HAZPRED_BIN"
  } else {
    "TEMPORAL_HAZARD_BIN"
  }
  option_name <- if (identical(binary, "hazpred")) {
    "temporal_hazard.hazpred_binary"
  } else {
    "temporal_hazard.binary"
  }

  external_path <- Sys.getenv(env_var, unset = "")
  if (!nzchar(external_path)) {
    external_path <- getOption(option_name, default = "")
  }

  if (nzchar(external_path)) {
    if (!file.exists(external_path)) {
      stop(
        "Binary path from ", env_var, "/options(", option_name, "=...) does not exist: ",
        external_path,
        call. = FALSE
      )
    }
    external_path
  } else {
    # Fallback to package-local binary in inst/bin.
    # Return NA if not available (allows graceful skipping in tests)
    local_candidates <- c(binary, paste0(binary, ".exe"))
    for (candidate in local_candidates) {
      bin_path <- system.file("bin", candidate, package = "TemporalHazard")
      if (nzchar(bin_path) && file.exists(bin_path)) {
        return(bin_path)
      }
    }
    NA
  }
}

#' Sanitise variable names for legacy SAS control files
#'
#' @param x Character vector of candidate names
#'
#' @return Uppercase names limited to SAS-friendly identifier characters
#'
#' @noRd
.hzr_as_legacy_names <- function(x) {
  if (length(x) == 0L) {
    character()
  } else {
    x <- toupper(gsub("[^A-Za-z0-9_]", "_", x))
    x <- gsub("^[^A-Za-z_]+", "X", x)
    make.unique(x, sep = "_")
  }
}

#' Infer how many leading theta elements are baseline shape parameters
#'
#' @param dist Distribution name
#' @param control Optional control list
#'
#' @return Integer scalar
#'
#' @noRd
.hzr_shape_parameter_count <- function(dist, control = list()) {
  if (!is.null(control$shape_param_count)) {
    as.integer(control$shape_param_count)
  } else {
    switch(
      dist,
      weibull = 2L,
      exponential = 1L,
      loglogistic = 2L,
      lognormal = 2L,
      0L
    )
  }
}

#' Split theta into baseline and covariate components for legacy output
#'
#' @param theta Numeric parameter vector
#' @param x_names Legacy covariate names used in the data extract
#' @param dist Distribution name
#' @param control Optional control list
#'
#' @return List with entries baseline and covariates
#'
#' @noRd
.hzr_split_theta_for_legacy <- function(theta, x_names, dist, control = list()) {
  if (is.null(theta) || !length(theta)) {
    list(baseline = numeric(), covariates = numeric())
  } else {
    theta_names <- names(theta)
    if (is.null(theta_names)) {
      theta_names <- rep("", length(theta))
    }
    theta_names_upper <- .hzr_as_legacy_names(
      ifelse(nzchar(theta_names), theta_names, paste0("TH", seq_along(theta)))
    )
    x_names_upper <- .hzr_as_legacy_names(x_names)

    matched_covariates <- theta_names_upper %in% x_names_upper
    if (any(matched_covariates)) {
      baseline <- theta[!matched_covariates]
      covariates <- theta[matched_covariates]
      names(baseline) <- theta_names_upper[!matched_covariates]
      names(covariates) <- theta_names_upper[matched_covariates]
      list(baseline = baseline, covariates = covariates)
    } else {
      shape_count <- .hzr_shape_parameter_count(dist = dist, control = control)
      shape_count <- min(shape_count, length(theta))
      baseline <- theta[seq_len(shape_count)]
      covariates <- if (length(theta) > shape_count) {
        theta[seq.int(shape_count + 1L, length(theta))]
      } else {
        numeric()
      }

      names(baseline) <- theta_names_upper[seq_len(shape_count)]
      if (length(covariates)) {
        covariate_names <- if (length(x_names_upper) >= length(covariates)) {
          x_names_upper[seq_along(covariates)]
        } else {
          theta_names_upper[seq.int(shape_count + 1L, length(theta))]
        }
        names(covariates) <- covariate_names
      }

      list(baseline = baseline, covariates = covariates)
    }
  }
}

#' Format a named numeric vector as SAS assignments
#'
#' @param values Named numeric vector
#'
#' @return Length-one character scalar
#'
#' @noRd
.hzr_format_legacy_assignments <- function(values) {
  if (!length(values)) {
    ""
  } else {
    parts <- paste0(
      names(values), "=",
      format(unname(values), digits = 10, scientific = TRUE, trim = TRUE)
    )
    paste(parts, collapse = ", ")
  }
}

#' Build PROC HAZARD option tokens from control settings
#'
#' @param control Optional control list
#'
#' @return Character vector of PROC option tokens
#'
#' @noRd
.hzr_format_proc_options <- function(control = list()) {
  options <- character()

  if (isTRUE(control$nocov)) {
    options <- c(options, "NOCOV")
  }
  if (isTRUE(control$nocor)) {
    options <- c(options, "NOCOR")
  }
  if (isTRUE(control$p)) {
    options <- c(options, "P")
  }
  if (isTRUE(control$conserve)) {
    options <- c(options, "CONSERVE")
  }
  if (isTRUE(control$quasi)) {
    options <- c(options, "QUASI")
  }
  if (!is.null(control$condition)) {
    options <- c(options, paste0("CONDITION=", control$condition))
  }
  if (!is.null(control$outhaz)) {
    options <- c(options, paste0("OUTHAZ=", .hzr_as_legacy_names(control$outhaz)))
  }

  options
}

#' Build a %HAZARD control file
#'
#' @param data_name Logical SAS dataset token to embed in PROC HAZARD
#' @param data_file Path to the temporary flat data file
#' @param time_var Legacy time variable name
#' @param status_var Legacy event variable name
#' @param x_names Legacy covariate names
#' @param dist Distribution name
#' @param theta Optional starting parameter vector
#' @param control Optional control list
#'
#' @return Character vector of SAS/control lines
#'
#' @noRd
.hzr_build_hazard_sas_input <- function(
    data_name,
    data_file,
    time_var,
    status_var,
    x_names = character(),
    dist,
    theta = NULL,
    control = list()) {

  sas_lines <- control$sas_lines %||% NULL
  if (!is.null(sas_lines)) {
    return(as.character(sas_lines))
  }

  proc_options <- .hzr_format_proc_options(control)
  proc_line <- paste(
    "PROC HAZARD",
    paste0("DATA=", data_name),
    paste(proc_options, collapse = " ")
  )
  proc_line <- trimws(gsub("[[:space:]]+", " ", proc_line))
  proc_line <- paste0(proc_line, ";")

  theta_parts <- .hzr_split_theta_for_legacy(theta = theta, x_names = x_names, dist = dist, control = control)
  parms_line <- if (length(theta_parts$baseline)) {
    paste0("     PARMS ", .hzr_format_legacy_assignments(theta_parts$baseline), ";")
  } else {
    NULL
  }

  fix_tokens <- control$fix %||% character()
  if (length(fix_tokens) && !is.null(parms_line)) {
    fix_tokens <- paste0("FIX", .hzr_as_legacy_names(fix_tokens))
    parms_line <- paste0(sub(";$", "", parms_line), " ", paste(fix_tokens, collapse = " "), ";")
  }

  phase_assignments <- control$phase_assignments %||% NULL
  covariate_lines <- character()
  if (length(theta_parts$covariates)) {
    if (is.null(phase_assignments)) {
      covariate_lines <- paste0(
        "     EARLY ",
        .hzr_format_legacy_assignments(theta_parts$covariates),
        ";"
      )
    } else if (is.list(phase_assignments)) {
      for (phase_name in names(phase_assignments)) {
        vars <- .hzr_as_legacy_names(unlist(phase_assignments[[phase_name]], use.names = FALSE))
        values <- theta_parts$covariates[names(theta_parts$covariates) %in% vars]
        if (length(values)) {
          covariate_lines <- c(
            covariate_lines,
            paste0("     ", toupper(phase_name), " ", .hzr_format_legacy_assignments(values), ";")
          )
        }
      }
    }
  }

  c(
    "%HAZARD(",
    paste("* TEMPORAL_HAZARD_DATA_FILE:", normalizePath(data_file, winslash = "/", mustWork = FALSE), ";"),
    paste("* DIST:", toupper(dist), ";"),
    proc_line,
    paste0("     TIME ", time_var, ";"),
    paste0("     EVENT ", status_var, ";"),
    parms_line,
    covariate_lines,
    ");"
  )
}

#' Build a %HAZPRED control file
#'
#' @param data_name Logical SAS dataset token to embed in PROC HAZARD
#' @param input_data_name Optional IN= dataset token for predictions
#' @param output_data_name Optional OUT= dataset token for predictions
#' @param time_var Legacy time variable name
#' @param control Optional control list
#'
#' @return Character vector of SAS/control lines
#'
#' @noRd
.hzr_build_hazpred_sas_input <- function(
    data_name,
    input_data_name,
    output_data_name,
    time_var,
    control = list()) {

  sas_lines <- control$sas_lines %||% NULL
  if (!is.null(sas_lines)) {
    return(as.character(sas_lines))
  }

  predict_types <- control$predict_types
  if (is.null(predict_types)) {
    predict_types <- c("SURVIVAL", "HAZARD")
  }
  predict_types <- toupper(predict_types)
  proc_parts <- c(
    "PROC HAZARD",
    paste0("DATA=", data_name)
  )
  if (!is.null(input_data_name)) {
    proc_parts <- c(proc_parts, paste0("IN=", input_data_name))
  }
  if (!is.null(output_data_name)) {
    proc_parts <- c(proc_parts, paste0("OUT=", output_data_name))
  }

  c(
    "%HAZPRED(",
    paste0(trimws(gsub("[[:space:]]+", " ", paste(proc_parts, collapse = " "))), ";"),
    paste0("     TIME ", time_var, ";"),
    paste0("     PREDICT ", paste(predict_types, collapse = " "), ";"),
    ");"
  )
}

#' Execute legacy binary using stdin/stdout redirection
#'
#' @param bin_path Path to executable
#' @param sas_file Path to control input file
#' @param lst_file Path to listing output file
#'
#' @return List with command string, status code, stdout text, and output lines
#'
#' @noRd
.hzr_run_legacy_binary <- function(bin_path, sas_file, lst_file) {
  command <- paste(
    shQuote(bin_path),
    "<",
    shQuote(sas_file),
    ">",
    shQuote(lst_file),
    "2>&1"
  )

  status_code <- system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
  output_text <- if (file.exists(lst_file)) readLines(lst_file, warn = FALSE) else character()

  list(
    call_string = command,
    status_code = unname(as.integer(status_code)),
    stdout = paste(output_text, collapse = "\n"),
    output_text = output_text
  )
}

#' Run hazard C binary on input data
#'
#' Invokes the compiled \code{hazard} executable with prepared input and parses output.
#' The legacy executable is called using shell redirection as
#' \code{hazard.exe < prefix.sas > prefix.lst 2>&1}.
#'
#' @param time Numeric vector of follow-up times
#' @param status Numeric vector of event indicators (0/1)
#' @param x Optional design matrix (univariable if NULL; multivar if provided)
#' @param theta Starting parameter vector (names must match C code expectations)
#' @param dist Distribution name ("weibull", etc.)
#' @param control Optional control list (nocov, nocor, condition, etc.)
#' @param data_name Character name for temporary data file prefix
#' @param binary Which legacy executable to run: \code{"hazard"} or \code{"hazpred"}
#'
#' @return List with:
#'   - \code{call_string}: Full command-line invocation for reproducibility
#'   - \code{stdout}: Raw standard output from hazard binary
#'   - \code{parsed}: Parsed results (data frame with estimates, SE, z, p-values)
#'   - \code{info}: Metadata (df, condition, convergence flags)
#'
#' @noRd
.hzr_run_hazard_binary <- function(
    time,
    status,
    x = NULL,
    theta = NULL,
    dist = "weibull",
    control = list(),
    data_name = "golden_test",
    binary = c("hazard", "hazpred")) {

  binary <- match.arg(binary)
  bin_path <- .hzr_get_hazard_binary(binary = binary)
  if (is.na(bin_path)) {
    return(list(
      call_string = NA_character_,
      stdout = "",
      parsed = .hzr_parse_hazard_output(character(), theta = theta),
      info = list(
        dist = dist,
        n = length(time),
        events = sum(status),
        binary = binary,
        binary_path = NA_character_,
        status_code = NA_integer_,
        sas_file = NA_character_,
        lst_file = NA_character_
      )
    ))
  }

  temp_dir <- tempdir()

  prefix <- file.path(temp_dir, data_name)
  data_file <- file.path(temp_dir, paste0(data_name, ".csv"))
  sas_file <- file.path(temp_dir, paste0(data_name, ".sas"))
  lst_file <- file.path(temp_dir, paste0(data_name, ".lst"))

  time_var <- .hzr_as_legacy_names(control$time_var %||% "TIME")
  status_var <- .hzr_as_legacy_names(control$status_var %||% "STATUS")
  dataset_name <- .hzr_as_legacy_names(control$data_name %||% toupper(data_name))

  df <- data.frame(
    time,
    as.integer(status)
  )
  names(df) <- c(time_var, status_var)

  if (!is.null(x)) {
    if (is.vector(x)) {
      x_names <- .hzr_as_legacy_names(colnames(as.matrix(x)) %||% "X1")
      df[[x_names[1L]]] <- x
    } else {
      x_names <- .hzr_as_legacy_names(colnames(x) %||% paste0("X", seq_len(ncol(x))))
      for (i in seq_len(ncol(x))) {
        df[[x_names[i]]] <- x[, i]
      }
    }
  } else {
    x_names <- character()
  }

  # Persist a simple flat data extract alongside the control program.
  write.csv(df, data_file, row.names = FALSE)

  if (!is.null(control$sas_file)) {
    sas_file <- control$sas_file
  } else {
    sas_lines <- if (identical(binary, "hazpred")) {
      .hzr_build_hazpred_sas_input(
        data_name = dataset_name,
        input_data_name = .hzr_as_legacy_names(control$input_data_name %||% dataset_name),
        output_data_name = .hzr_as_legacy_names(control$output_data_name %||% paste0(dataset_name, "_OUT")),
        time_var = time_var,
        control = control
      )
    } else {
      .hzr_build_hazard_sas_input(
        data_name = dataset_name,
        data_file = data_file,
        time_var = time_var,
        status_var = status_var,
        x_names = x_names,
        dist = dist,
        theta = theta,
        control = control
      )
    }
    writeLines(sas_lines, sas_file, useBytes = TRUE)
  }

  exec_result <- .hzr_run_legacy_binary(
    bin_path = bin_path,
    sas_file = sas_file,
    lst_file = lst_file
  )

  parsed <- if (!is.na(exec_result$status_code) && exec_result$status_code == 0L) {
    .hzr_parse_hazard_output(exec_result$output_text, theta = theta)
  } else {
    structure(
      data.frame(
        parameter = character(),
        estimate = numeric(),
        std_err = numeric(),
        z_stat = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
      ),
      output_text = paste(exec_result$output_text, collapse = "\n"),
      parse_status = "binary_failed"
    )
  }

  list(
    call_string = exec_result$call_string,
    stdout = exec_result$stdout,
    parsed = parsed,
    info = list(
      dist = dist,
      n = length(time),
      events = sum(status),
      binary = binary,
      binary_path = bin_path,
      status_code = exec_result$status_code,
      prefix = prefix,
      data_file = data_file,
      sas_file = sas_file,
      lst_file = lst_file
    )
  )
}

#' Parse hazard binary output
#'
#' Parses tabular output from the C hazard binary into a structured data frame.
#' Expects columns: parameter, estimate, StdErr (or SE), z-value (or z), p-value (or pval).
#' The parser is flexible and case-insensitive for column matching.
#'
#' @param output_text Character vector (lines) from binary output file
#' @param theta Optional parameter vector for context/validation
#'
#' @return Data frame with columns: parameter, estimate, std_err, z_stat, p_value
#' @keywords internal
.hzr_parse_hazard_output <- function(output_text, theta = NULL) {
  if (length(output_text) == 0 || all(nchar(trimws(output_text)) == 0)) {
    # Empty output: return empty data frame with proper structure
    return(structure(
      data.frame(
        parameter = character(),
        estimate = numeric(),
        std_err = numeric(),
        z_stat = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
      ),
      output_text = "",
      parse_status = "empty"
    ))
  }

  # Join text into single string and parse as table
  text <- paste(output_text, collapse = "\n")
  text_trimmed <- trimws(text)

  # Try to read as whitespace-separated table (skip comment lines)
  tryCatch(
    {
      lines <- strsplit(text_trimmed, "\n")[[1]]
      # Find header line (contains keywords like 'Parameter', 'Estimate', etc.)
      header_idx <- grep("(?i)(parameter|estimate|std|err|z[-_]?stat|p[-_]?val)", lines, perl = TRUE)[1]

      if (is.na(header_idx)) {
        warning("Could not identify header line in hazard output", call. = FALSE)
        return(structure(
          data.frame(
            parameter = character(),
            estimate = numeric(),
            std_err = numeric(),
            z_stat = numeric(),
            p_value = numeric(),
            stringsAsFactors = FALSE
          ),
          output_text = text,
          parse_status = "no_header"
        ))
      }

      # Read the table from R data (skip header line)
      data_start <- header_idx + 1
      data_lines <- lines[data_start:length(lines)]
      data_text <- paste(data_lines, collapse = "\n")

      # Parse as fixed-format or space-separated table
      df <- read.table(
        text = data_text,
        header = FALSE,
        colClasses = "character",
        fill = TRUE,
        na.strings = c("NA", ".", "NaN", "Inf"),
        stringsAsFactors = FALSE
      )

      if (nrow(df) == 0) {
        return(structure(
          data.frame(
            parameter = character(),
            estimate = numeric(),
            std_err = numeric(),
            z_stat = numeric(),
            p_value = numeric(),
            stringsAsFactors = FALSE
          ),
          output_text = text,
          parse_status = "empty_table"
        ))
      }

      # Map columns based on position and header keywords
      # Standard format: Parameter, Estimate, StdErr, Z, P
      # This is flexible and should handle variations

      # Assume first column is parameter names, next 4 are numeric
      param_col <- df[[1]]
      est_col <- suppressWarnings(as.numeric(df[[2]]))
      se_col <- suppressWarnings(as.numeric(df[[3]]))
      z_col <- suppressWarnings(as.numeric(df[[4]]))
      p_col <- suppressWarnings(as.numeric(df[[5]]))

      result <- data.frame(
        parameter = as.character(param_col),
        estimate = est_col,
        std_err = se_col,
        z_stat = z_col,
        p_value = p_col,
        stringsAsFactors = FALSE
      )

      structure(
        result,
        output_text = text,
        parse_status = "success"
      )
    },
    error = function(e) {
      warning("Failed to parse hazard binary output: ", e$message, call. = FALSE)
      structure(
        data.frame(
          parameter = character(),
          estimate = numeric(),
          std_err = numeric(),
          z_stat = numeric(),
          p_value = numeric(),
          stringsAsFactors = FALSE
        ),
        output_text = text,
        parse_status = paste("error:", e$message)
      )
    }
  )
}

#' Generate golden fixture for univariable Weibull model
#'
#' Runs the C hazard binary on a dataset and stores results as test fixture.
#'
#' @param fixture_name Character identifying the fixture (e.g., "avc_death_hz")
#' @param time Numeric vector of follow-up times
#' @param status Numeric vector of event indicators
#' @param x Optional design matrix for multivariable models
#' @param theta Optional starting parameters
#' @param control Optional control list
#' @param output_dir Directory to store the fixture. Required; the function
#'   never writes to a default path. Pass \code{tempdir()} for a throwaway
#'   run or an explicit project path (e.g. \code{"inst/fixtures"}) to persist.
#'
#' @return Invisibly returns the parsed output of the reference binary.
#'
#' @noRd
.hzr_generate_golden_fixture <- function(
    fixture_name,
    time,
    status,
    x = NULL,
    theta = NULL,
    control = list(),
    output_dir) {

  stopifnot(
    "`output_dir` must be a single non-empty path" =
      is.character(output_dir) && length(output_dir) == 1L && nzchar(output_dir)
  )
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Run binary
  result <- .hzr_run_hazard_binary(
    time = time,
    status = status,
    x = x,
    theta = theta,
    control = control,
    data_name = fixture_name
  )

  # Save fixture as RDS for reproducible testing
  fixture_file <- file.path(output_dir, paste0(fixture_name, ".rds"))
  saveRDS(result, fixture_file)

  message("Generated golden fixture: ", fixture_file)

  invisible(result$parsed)
}
