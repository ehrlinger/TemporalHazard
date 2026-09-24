# Regenerate tests/testthat/fixtures/setg1-oracle.csv (#424, #421 item 3).
#
# Runs the HAZARD C binary on one early-phase PARMS statement per grid row and
# records what it PRODUCES, on two datasets. The binary is the oracle for the
# class the translator must give each job. Run from the package root:
#
#   Rscript data-raw/setg1-oracle.R [path/to/hazard]
#
# The binary reads its control file on stdin and its data from
# $TMPDIR/hzr.<JOBID>.<JOBIX>.dta, an XPORT file, with OBSCOUNT set to the row
# count. Without the data file every job dies at "open input file", and the
# classes cannot be told apart.
#
# The verdict is read from the fatal-exit reason plus the FULL set of fitted
# markers, never from one message string (the method of stream D's
# ~/.cache/th-gate/D/setg1_oracle_synthetic.R, ported here):
#   "runs"      every fitted marker is present;
#   "refused"   a fatal SYNTAX or SEMANTIC exit;
#   "no_result" a fatal ERROR exit, or a raised error code with the fitted
#               markers missing;
#   "unclear"   anything else.
#
# The two datasets (#468 review):
#   "ref"   the package's `avc` data (complete cases); `stable_ref` says whether
#           the class held with its times multiplied by 10 and by 0.1;
#   "synth" stream D's dataset "synthetic-n200-seed5", ported unchanged from
#           D's generator: seeded, no study data.
# Each dataset has its own controls, which must land on their SPECIFIC
# verdicts: a positive control that reads "runs" and a negative control that
# does not run. The script stops if either fails on either dataset.
#
# `class` is the verdict the translator must give. The translator speaks only
# about the shapes SETG1 singles out, so:
#   ref "runs"      -> "runs": a job PROC HAZARD fits on the reference data is
#                      not warned about, even where it fails to converge on
#                      other data -- that is fitting, not a verdict on the job;
#   ref "refused"   -> "refused", and it must be refused on synth too;
#   ref "no_result" -> "no_result" when synth also gives none, and
#                      "may_not_fit" when synth fits: the failure depends on
#                      the data, so the translator may only say SAS MAY stop.

args <- commandArgs(TRUE)
bin <- if (length(args)) args[[1L]] else
  path.expand("~/Documents/GitHub/hazard/dist/bin/hazard")
stopifnot(file.exists(bin))
work <- tempfile("setg1-oracle-")
dir.create(work)

grid <- c(
  known_positive = "MUE=0.2 THALF=1 NU=1 M=1",
  SETG1910_neg = "MUE=0.2 THALF=-1 FIXTHALF NU=1 M=1",
  SETG1910_zero = "MUE=0.2 THALF=0 FIXTHALF NU=1 M=1",
  thalf_neg_free = "MUE=0.2 THALF=-1 NU=1 M=1",
  thalf_zero_free = "MUE=0.2 THALF=0 NU=1 M=1",
  thalf_half_free = "MUE=0.2 THALF=-.5",
  SETG1940 = "MUE=0.2 THALF=1 M=-1 NU=-1 FIXM FIXNU",
  SETG1950 = "MUE=0.2 THALF=1 M=0 NU=0 FIXM FIXNU",
  SETG1960 = "MUE=0.2 THALF=1 M=1 NU=0 FIXM FIXNU",
  SETG1920 = "MUE=0.2 THALF=1 M=-1 NU=-1 FIXM FIXNU FIXMNU1",
  SETG1930 = "MUE=0.2 THALF=1 M=0 NU=0 FIXM FIXNU FIXMNU1",
  SETG1900 = "MUE=0.2 THALF=1 NU=1 M=1 DELTA=-2 FIXDELTA",
  SETG1901 = "MUE=0.2 THALF=1 NU=1 M=1 DELTA=2 FIXDELTA",
  mnu_zero_free = "MUE=0.2 THALF=1 M=0 NU=0",
  g1flag4_mpos = "MUE=0.2 THALF=1 M=1 NU=0",
  g1flag4_mpos_fixnu = "MUE=0.2 THALF=1 M=1 NU=0 FIXNU",
  g1flag4_mneg = "MUE=0.2 THALF=1 M=-1 NU=0",
  g1flag4_mneg_fixnu = "MUE=0.2 THALF=1 M=-1 NU=0 FIXNU",
  # g1flag 4 with M fixed: SETG1 reaches the same branch, but M does not move.
  g1flag4_mneg_both = "MUE=0.2 THALF=1 M=-1 NU=0 FIXM FIXNU",
  g1flag4_mzero_fixnu = "MUE=0.2 THALF=1 M=0 NU=0 FIXNU",
  # Known negatives: jobs PROC HAZARD runs.
  nu_zero_fixm = "MUE=0.2 THALF=1 M=1 NU=0 FIXM",
  thalf_fixed_pos = "MUE=0.2 THALF=0.5 FIXTHALF NU=1 M=1",
  mnu_both_fixed_pos = "MUE=0.2 THALF=1 NU=2 M=1 FIXM FIXNU",
  shape_defaults = "MUE=0.2"
)

# --- ref: the package's avc data ---
load("data/avc.rda")
a <- avc[stats::complete.cases(avc), ]
D0 <- data.frame(TT = a$int_dead, DEAD = a$dead,
                 AGE = as.numeric(scale(a$age)))
# --- synth: stream D's dataset, unchanged. Change NOTHING without the label. ---
SYNTH_LABEL <- "synthetic-n200-seed5"
set.seed(5)
n <- 200L
D1 <- data.frame(TT = stats::rexp(n, 0.25) + 0.01,
                 DEAD = stats::rbinom(n, 1, 0.8),
                 AGE = as.numeric(scale(stats::rnorm(n))))

# Controls per dataset. `M=1 NU=1` does not fit cleanly on synth (D measured
# "unclear"), so each dataset carries a positive control that does.
controls <- list(
  ref = c(positive = "MUE=0.2 THALF=1 NU=1 M=1",
          negative = "MUE=0.2 THALF=1 M=-1 NU=-1 FIXM FIXNU"),
  synth = c(positive = "MUE=0.2 THALF=0.5 M=-1 NU=1 FIXM FIXNU FIXMNU1",
            negative = "MUE=0.2 THALF=0.5 M=-1 NU=-1 FIXM FIXNU FIXMNU1")
)

fitted_markers <- c("Estimates for Model Parameters", "Final Results",
                    "Parameter Estimate Summary", "Log likelihood",
                    "Optimization terminated")

run_one <- function(parms, ix, D) {
  haven::write_xpt(D, file.path(work, sprintf("hzr.J424.%s.dta", ix)),
                   version = 5, name = "HZRCALL")
  ctl <- file.path(work, sprintf("hzr.J424.%s.sas", ix))
  lst <- file.path(work, sprintf("hzr.J424.%s.lst", ix))
  writeLines(c(paste0("(   PROC HAZARD DATA=D; EVENT DEAD; TIME TT; PARMS ",
                      parms, ";   )"),
               sprintf("; OBSCOUNT %d ;", nrow(D)), "JOBID J424;",
               sprintf("JOBIX %s;", ix)), ctl)
  system(paste0("TMPDIR=", shQuote(work), " ", shQuote(bin), " < ",
                shQuote(ctl), " > ", shQuote(lst), " 2>&1"))
  out <- readLines(lst, warn = FALSE)
  fatal <- sub("^fatal exit: ", "",
               unlist(regmatches(out, gregexpr("fatal exit: [A-Z]+", out))))
  raised <- sub("^\\[hazard ERROR [^]]*\\] ", "",
                out[grepl("^\\[hazard ERROR", out)])
  nmark <- sum(vapply(fitted_markers,
                      function(m) any(grepl(m, out, fixed = TRUE)), TRUE))
  verdict <- if (nmark == length(fitted_markers)) "runs" else
    if (length(fatal) && fatal[[1L]] %in% c("SYNTAX", "SEMANTIC")) "refused" else
      if (length(fatal) || length(raised)) "no_result" else "unclear"
  list(
    class = verdict,
    reason = if (length(raised)) raised[[1L]] else
      if (verdict == "runs") "fitted (all markers)" else "",
    version = grep("C-Version", out, value = TRUE)
  )
}

for (ds in names(controls)) {
  D <- if (ds == "ref") D0 else D1
  pos <- run_one(controls[[ds]][["positive"]], paste0("P", toupper(ds)), D)
  neg <- run_one(controls[[ds]][["negative"]], paste0("N", toupper(ds)), D)
  cat(ds, "controls: positive", pos$class, "| negative", neg$class, "\n")
  stopifnot(pos$class == "runs", neg$class != "runs")
}

rows <- lapply(seq_along(grid), function(k) {
  base <- run_one(grid[[k]], sprintf("X%d", k), D0)
  synth <- run_one(grid[[k]], sprintf("S%d", k), D1)
  scaled <- vapply(c(10, 0.1), function(s) {
    D <- D0
    D$TT <- D$TT * s
    run_one(grid[[k]], sprintf("%s%d", if (s > 1) "A" else "B", k), D)$class
  }, character(1))
  class <- switch(
    base$class,
    runs = "runs",
    refused = if (synth$class == "refused") "refused" else
      stop("refused on the reference data only: ", grid[[k]]),
    no_result = if (synth$class == "runs") "may_not_fit" else "no_result",
    stop("unclear verdict on the reference data: ", grid[[k]])
  )
  data.frame(id = names(grid)[[k]], parms = grid[[k]], class = class,
             binary_ref = base$class, reason_ref = base$reason,
             stable_ref = all(scaled == base$class),
             binary_synth = synth$class, reason_synth = synth$reason,
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
version <- run_one(grid[[1L]], "V1", D0)$version[[1L]]
dest <- "tests/testthat/fixtures/setg1-oracle.csv"
writeLines(c(
  "# Generated by data-raw/setg1-oracle.R; do not edit by hand.",
  paste0("# Binary: ", sub(path.expand("~"), "~", bin, fixed = TRUE), " (",
         trimws(sub("^\\$Note:", "", version)), ")"),
  paste0("# Binary file date: ", format(file.mtime(bin), "%Y-%m-%d")),
  paste0("# Source checkout when generated: ",
         system2("git", c("-C", dirname(dirname(dirname(bin))), "rev-parse",
                          "HEAD"), stdout = TRUE)),
  paste0("# Datasets: ref = avc (complete cases); synth = ", SYNTH_LABEL,
         " (stream D's generator, ported)"),
  "# Controls passed on each dataset: positive reads runs, negative does not.",
  paste0("# Generated: ", format(Sys.Date()))
), dest)
suppressWarnings(utils::write.table(tab, dest, append = TRUE, sep = ",",
                                    row.names = FALSE, qmethod = "double"))
print(tab[, c("id", "class", "binary_ref", "binary_synth", "stable_ref")],
      row.names = FALSE)
