* **`hzr_deciles(event_time =)` on a multiphase `time_windows` fit used the
  exit window's coefficient.** The #629 fix re-expanded the design at the
  horizon only. A caller's `event_time` kept the design expanded at each
  subject's stored exit time, so a row whose `event_time` fell in another
  window was charged that window's coefficient in the expected events, with
  no warning. The design is now expanded at `event_time` whenever it differs
  from the stored times. The default `event_time` is unaffected.
* **`hzr_bootstrap()`'s warning on dropped replicates now says what the
  summary is conditional on.** Dropping the replicates whose likelihood is
  flat or still rising leaves a summary over only the resamples in which the
  model is identified, which understates the uncertainty and moves the other
  parameters' intervals too. On `avc` with a constant phase carrying `mal`,
  17 of 30 replicates were dropped and the `constant.mal` interval narrowed
  from about [-19, 12] to [-0.1, 1.8]. The warning says so, and when 20% or
  more of the replicates are dropped it says to treat the summary as
  unreliable.
