* The examples for `hzr_translate_sas()`, `hzr_gof()` and `hzr_deciles()` are no longer
  wrapped in `\donttest{}`. Each runs in well under a second, so `R CMD check` now runs them.
* `hzr_bootstrap()` gains a quick example (ten resamples) that runs during `R CMD check`.
  The longer examples stay in `\donttest{}`.
* The second reference in the `DESCRIPTION` now reads Rajeswaran et al. (2018), the
  authors (year) form CRAN asks for. The year, 2018, is the print issue
  (Statistical Methods in Medical Research 27(1)); the article appeared online in 2016.
