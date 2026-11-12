* `hzr_stepwise()`: when one of a multiphase candidate refit's two starts
  refuses because a shape parameter ran to a limit and the other start's fit
  is kept, the refusal is no longer dropped. It is kept on the fit as
  `$fit$refit_refusal`, and a warning of class `hzr_limit_refit_refused`
  names the limit and that start's objective (#736).
* `hzr_bootstrap()` with a `scope` now counts a replicate as restarted when
  the refit of its base model was restarted, not only when its final fit
  was; every candidate in that replicate was measured against the base
  (#736).
* `control$limit_restart` is now validated on an unfitted call too
  (`fit = FALSE`), whose control a later stepwise or bootstrap refit uses
  (#736).
* A hard limit refusal now reports the objective from the user's own start.
  When a restart that also ran to a limit went higher, that value is named
  separately: a fit at a limit can reach a degenerate likelihood, such as a
  spike on tied times, which the message used to call "the best objective
  reached". The refusal also carries `$start_objective` (#736).
* The limit line and messages now say "no start clear of a limit did better,
  original kept" where they said "none better", which was wrong when a start
  at a limit did better, and describe a flat direction as
  "early.nu: a flat direction (not identified)" rather than "ran to" one
  (#736).
* Documentation: the README no longer says the shapes are always estimated
  jointly beside an example that fixes them, and the architecture vignette no
  longer says the multi-start offsets are random draws; they do not use R's
  random number stream (#736).
* **The limit layer now catches the `m < 0` joint limit (#734).** There
  `nu` and `t_half` run to infinity together and a cdf or hazard phase
  becomes a power law in t over the data. No boundary record flags it, and
  when the Hessian was still usable, such fits came back with
  `converged = TRUE` and a clean `$fit$limit`. On a fit nothing else flagged, a free `t_half` outside the
  observed times, or in a phase with `m < 0`, is now moved across the data's
  whole log-span with the phase's other parameters refitted. A loss in the
  log-likelihood of less than about 0.23 means `t_half` is not identified; a
  gain means the fit was not a maximum, and fires the check too. As
  with a flat direction, this tries the restart grid and keeps the original
  unless a start clear of the limit does better, with a usable Hessian of its
  own. It never refuses. A move that cannot be computed does not count as
  a flat span: the other direction is still tried, and when neither is
  measured flat the record and messages say "an unchecked span".
