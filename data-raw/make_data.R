## ── Convert inst/extdata CSVs to lazy-loaded data/ objects ──────────────────
##
## Run this script from the package root:
##   source("data-raw/make_data.R")
##
## It reads each CSV from inst/extdata/, creates the corresponding data.frame,
## and saves it as a compressed .rda file in data/.
##
## All CSVs originate from SAS exports where "." denotes missing values.

dir.create("data", showWarnings = FALSE)

.read_sas_csv <- function(path, ...) {
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", "."), ...)
}

# ── AVC: Atrioventricular Canal Repair (Cleveland Clinic, 1977–1993) ─────────
avc <- .read_sas_csv("inst/extdata/avc.csv")
stopifnot(nrow(avc) == 310)
usethis::use_data(avc, overwrite = TRUE, compress = "xz")

# ── CABGKUL: Primary Isolated CABG (KU Leuven, 1971–2008) ──────────────────
cabgkul <- .read_sas_csv("inst/extdata/cabgkul.csv")
stopifnot(nrow(cabgkul) == 5880, sum(cabgkul$dead) == 545)
usethis::use_data(cabgkul, overwrite = TRUE, compress = "xz")

# ── OMC: Open Mitral Commissurotomy ─────────────────────────────────────────
omc <- .read_sas_csv("inst/extdata/omc.csv")
stopifnot(nrow(omc) == 339)
usethis::use_data(omc, overwrite = TRUE, compress = "xz")

# ── TGA: Transposition of the Great Arteries (BCH/CHOP) ────────────────────
tga <- .read_sas_csv("inst/extdata/tga.csv")
stopifnot(nrow(tga) == 470)
usethis::use_data(tga, overwrite = TRUE, compress = "xz")

# ── Valves: Heart Valve Replacement ─────────────────────────────────────────
valves <- .read_sas_csv("inst/extdata/valves.csv")
stopifnot(nrow(valves) == 1533)
usethis::use_data(valves, overwrite = TRUE, compress = "xz")

message("All 5 datasets saved to data/")

# ── uslife2023: US Life Table 2023, all-interval-censored SAS parity anchor ──
#
# Derived from /studies/general/uslife/table2023/datasets/built.sas7bdat by
# `IF D_ALL=0 THEN DELETE`, keeping the three columns the parity fit uses.
# The CSV in inst/extdata/ is the durable source: it must remain reproducible
# after the SAS license lapses and the .sas7bdat becomes unreadable.
#
# The gates below are the figures printed by the SAS job's own .lst
# (distributions/hz.icall.lst), not values recomputed from this CSV -- they
# fail loudly if the CSV is ever regenerated from a different dataset.
uslife2023 <- .read_sas_csv("inst/extdata/uslife2023.csv")
stopifnot(
  nrow(uslife2023) == 124L,
  identical(names(uslife2023), c("age_l", "age_u", "d_all")),
  # Every interval is exactly one year wide: this is what makes the fixture
  # the clean anchor, since log(u - l) = 0 switches the width term off.
  all(uslife2023$age_u - uslife2023$age_l == 1),
  abs(sum(uslife2023$d_all) - 100000.0125) < 5e-5,
  abs(min(uslife2023$d_all) - 0.2352) < 5e-5,
  abs(max(uslife2023$d_all) - 3620.335) < 5e-4
)
usethis::use_data(uslife2023, overwrite = TRUE, compress = "xz")
