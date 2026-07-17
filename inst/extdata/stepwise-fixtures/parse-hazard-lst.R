# parse-hazard-lst.R -- Parse a HAZARD C-binary .lst listing into the three
# files consumed by .hzr_build_stepwise_fixture().
#
# Usage (from the project root, after devtools::load_all()):
#
#   source("inst/extdata/stepwise-fixtures/parse-hazard-lst.R")
#
# Inputs:
#   inst/extdata/stepwise-fixtures/payload/capture/avc-forward-wald.lst
#   inst/extdata/stepwise-fixtures/payload/capture/stepwise_meta.txt
#
# Outputs (written in-place to payload/capture/):
#   stepwise_trace.csv   -- one row per accepted step
#   stepwise_final.csv   -- final coefficient table
#   stepwise_meta.txt    -- updated with logLik= entry
#
# Then builds:
#   inst/fixtures/stepwise-avc-forward-wald.rds

lst_file    <- "inst/extdata/stepwise-fixtures/payload/capture/avc-forward-wald.lst"
meta_file   <- "inst/extdata/stepwise-fixtures/payload/capture/stepwise_meta.txt"
out_trace   <- "inst/extdata/stepwise-fixtures/payload/capture/stepwise_trace.csv"
out_final   <- "inst/extdata/stepwise-fixtures/payload/capture/stepwise_final.csv"
fixture_out <- "inst/fixtures/stepwise-avc-forward-wald.rds"

stopifnot(file.exists(lst_file), file.exists(meta_file))

lines <- readLines(lst_file, warn = FALSE)

# ---- 1. Step trace ----------------------------------------------------------
# Banner format (lines after the *---* separator):
#   |Step   N. VAR   added to PHASE PHASE                  |
#   |                   with Z =   Z_VAL and P = P_VAL     |
# or
#   |Step   N. VAR   deleted from PHASE PHASE              |
#   |                   with Z =   Z_VAL and P = P_VAL     |
# Step 0 is the intercept-only baseline; skip it.

step_idx <- grep("^\\|Step\\s+\\d+\\.", lines)
step_rows <- vector("list", length(step_idx))
n_steps <- 0L

for (i in step_idx) {
  line <- lines[i]
  step_num <- as.integer(sub("^\\|Step\\s+(\\d+)\\..*", "\\1", line))
  if (step_num == 0L) next

  if (grepl("added to", line, ignore.case = TRUE)) {
    action <- "enter"
    rest   <- sub("^\\|Step\\s+\\d+\\.\\s+", "", line)
    variable <- trimws(sub("\\s+added to.*", "", rest))
    phase_str <- trimws(sub(".*added to\\s+(\\S+)\\s+PHASE.*", "\\1", line,
                             ignore.case = TRUE))
  } else if (grepl("deleted from", line, ignore.case = TRUE)) {
    action <- "drop"
    rest   <- sub("^\\|Step\\s+\\d+\\.\\s+", "", line)
    variable <- trimws(sub("\\s+deleted from.*", "", rest))
    phase_str <- trimws(sub(".*deleted from\\s+(\\S+)\\s+PHASE.*", "\\1", line,
                             ignore.case = TRUE))
  } else {
    next
  }

  phase <- tolower(phase_str)

  # The next line carries "with Z = ... and P = ..."
  z_line <- if (i + 1L <= length(lines)) lines[i + 1L] else ""
  z_val  <- as.numeric(sub(".*with Z\\s*=\\s*([+-]?[0-9.]+).*", "\\1", z_line))

  p_raw  <- trimws(sub(".*and P\\s*=\\s*(.+?)\\s*\\|.*", "\\1", z_line))
  p_val  <- if (startsWith(p_raw, "<")) {
    as.numeric(trimws(sub("^<\\s*\\.?", "0.", p_raw)))
  } else {
    as.numeric(p_raw)
  }

  n_steps <- n_steps + 1L
  step_rows[[n_steps]] <- data.frame(
    step_num = step_num,
    action   = action,
    variable = variable,
    phase    = phase,
    stat     = z_val^2,   # Wald chi-square = Z^2 (recorded for reference; the
                          # current parity test compares final-model agreement,
                          # not per-step statistics)
    df       = 1L,
    p_value  = p_val,
    stringsAsFactors = FALSE
  )
}

trace_df <- do.call(rbind, step_rows[seq_len(n_steps)])
# SAS uppercases all variable names; R's avc uses lowercase.
trace_df$variable <- tolower(trace_df$variable)
utils::write.csv(trace_df, out_trace, row.names = FALSE)
message("Wrote ", out_trace, " (", nrow(trace_df), " steps)")

# ---- 2. Final coefficient table ---------------------------------------------
# Lives in the "Parameter Estimate Summary" block under "Final Results:".
# Format (fixed-width columns):
#
#   Phase      Parameter    Estimate       Std error             Z    Prob>|Z|
#   ----------------------...
#   Early:     E0             -2.20788      0.4338312         -5.089    <.0001
#              AGE          -0.00544498    0.002572203         -2.117    0.0343
#              ...
#              ------...
#   Constant:  C0             -10.6589       1.726197         -6.175    <.0001
#              ...
#   ------...
#
# The inner --------... separates phases.
# The outer (longer) --------... ends the table.

param_summary_idx <- grep("Parameter Estimate Summary", lines, fixed = TRUE)[1]
stopifnot(!is.na(param_summary_idx))

# Header line: "Phase  Parameter  Estimate  Std error  Z  Prob>|Z|"
col_header_idx <- grep(
  "Phase\\s+Parameter\\s+Estimate\\s+Std error",
  lines[param_summary_idx:length(lines)]
)[1] + param_summary_idx - 1L

# Data starts after the dashes line immediately below the header
data_start <- col_header_idx + 2L   # +1 dashes, +1 first data line

# Stop when we hit "Estimates for Model Parameters" (next section)
stop_idx <- grep("Estimates for Model Parameters", lines, fixed = TRUE)
stop_idx <- stop_idx[stop_idx > data_start][1L]

final_rows    <- list()
current_phase <- NA_character_

for (i in seq(data_start, stop_idx - 1L)) {
  line    <- lines[i]
  trimmed <- trimws(line)

  if (nchar(trimmed) == 0L) next
  if (grepl("^-{5,}", trimmed))  next   # separator line

  # Phase label on this line?
  phase_lbl <- regmatches(
    line, regexpr("(Early|Constant):", line, perl = TRUE, ignore.case = TRUE)
  )
  if (length(phase_lbl) > 0L && nchar(phase_lbl) > 0L) {
    current_phase <- tolower(sub(":", "", phase_lbl))
    # Strip the "Early:  " or "Constant:  " prefix before parsing numeric fields
    line <- sub(".*(?:Early|Constant):\\s+", "", line, ignore.case = TRUE,
                perl = TRUE)
  }

  # Expect: PARAM  estimate  std_error  z_stat  prob
  # Prob may be "<.0001"
  m <- regexpr(
    paste0("^\\s*(\\S+)\\s+([+-]?[0-9.Ee+-]+)\\s+([0-9.Ee+-]+)\\s+",
           "([+-]?[0-9.Ee+-]+)\\s+([<]?[.0-9]+)"),
    line, perl = TRUE
  )
  if (m < 0L) next

  tok <- strsplit(trimws(regmatches(line, m)), "\\s+")[[1]]
  if (length(tok) < 5L) next

  p_raw <- tok[5L]
  p_val <- if (startsWith(p_raw, "<")) {
    as.numeric(sub("^<\\.?", "0.", p_raw))
  } else {
    as.numeric(p_raw)
  }

  final_rows[[length(final_rows) + 1L]] <- data.frame(
    variable  = tok[1L],
    phase     = current_phase,
    estimate  = as.numeric(tok[2L]),
    std_error = as.numeric(tok[3L]),
    z_stat    = as.numeric(tok[4L]),
    p_value   = p_val,
    stringsAsFactors = FALSE
  )
}

final_df <- do.call(rbind, final_rows)
final_df$variable <- tolower(final_df$variable)
utils::write.csv(final_df, out_final, row.names = FALSE)
message("Wrote ", out_final, " (", nrow(final_df), " parameters)")

# ---- 3. Augment meta.txt with logLik ----------------------------------------
# Grab the LAST "Log likelihood = ..." line in the file (= final model's).
ll_val <- as.numeric(trimws(sub(
  ".*Log likelihood\\s*=\\s*", "",
  tail(grep("Log likelihood\\s*=", lines, value = TRUE), 1L)
)))
message("Final log-likelihood: ", ll_val)

meta_lines <- readLines(meta_file)
# Correct the misleading metadata: this is a two-phase MULTIPHASE model, and
# SAS's SELECTION always scores with Q (the score test), never a Wald. The old
# capture recorded dist=weibull / criterion=wald, but "weibull" is only the
# early phase's SHAPING family and "wald" described the R side it was once
# compared against -- neither describes what SAS computed. Record the truth:
# dist=multiphase, shaping=weibull, criterion=score.
force_kv <- function(lines, key, value) {
  lines <- lines[!grepl(paste0("^", key, "\\s*="), lines)]
  c(lines, paste0(key, "=", value))
}
meta_lines <- force_kv(meta_lines, "dist", "multiphase")
meta_lines <- force_kv(meta_lines, "shaping", "weibull")
meta_lines <- force_kv(meta_lines, "criterion", "score")
meta_lines <- force_kv(meta_lines, "logLik", ll_val)
writeLines(meta_lines, meta_file)
message("Updated ", meta_file)

# ---- 4. Build the .rds fixture ----------------------------------------------
devtools::load_all(quiet = TRUE)

fix <- TemporalHazard:::.hzr_build_stepwise_fixture(
  trace_csv = out_trace,
  final_csv = out_final,
  meta_txt  = meta_file,
  out_path  = fixture_out
)
message("Built fixture: ", fixture_out)
message("\nStep trace:")
print(fix$steps)
message("\nFinal coefficients:")
print(fix$final$coef)
message("\nlogLik: ", fix$final$logLik)
