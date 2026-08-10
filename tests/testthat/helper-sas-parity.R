# Shim. The SAS .lst parsers themselves live at
# inst/sas-parity/helper-sas-parity.R in the source tree -- installed as
# <library>/TemporalHazard/sas-parity/helper-sas-parity.R, since R strips the
# inst/ prefix -- so that they ship with an installed package. testthat
# auto-sources this file, which loads them.
#
# Why they moved out of tests/: R CMD INSTALL skips tests/ unless
# --install-tests is passed, so a plain install.packages() or
# remotes::install_github() left the parsers unreachable. Downstream study
# pipelines that check their own SAS .lst output against these parsers had no
# way to get at them short of cloning the repository. Everything under inst/
# is installed unconditionally, so system.file() now resolves them.
#
# Keeping this shim -- rather than sourcing from inst/ inside each test file
# -- preserves testthat's helper auto-sourcing: tests see the parser functions
# exactly as they did before, and no test file needed to change.

.hzr_sas_parity_helper <- system.file(
  "sas-parity", "helper-sas-parity.R",
  package = "TemporalHazard"
)

if (!nzchar(.hzr_sas_parity_helper) || !file.exists(.hzr_sas_parity_helper)) {
  stop("helper-sas-parity.R: SAS parsers not found. Expected them at ",
       "sas-parity/helper-sas-parity.R inside the installed package ",
       "(inst/sas-parity/ in the source tree), reachable via ",
       "system.file(\"sas-parity\", \"helper-sas-parity.R\", ",
       "package = \"TemporalHazard\"). Reinstall the package if this ",
       "persists.", call. = FALSE)
}

# local = TRUE evaluates into this file's environment, which testthat has set
# to the test environment -- so the parsers land where the tests can see them.
source(.hzr_sas_parity_helper, local = TRUE)
