* The SAS-parity data resolvers no longer fall back to a hard-coded study
  path; set `HAZARD_MAZE_DATASETS` and `HAZARD_ACHALASIA_DIR` (once, in
  `~/.Renviron`) to run the volume-gated parity tests and
  `tools/release-suite.R`. The repository is public, and a study's folder name
  says which study it is.
