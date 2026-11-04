- **The score criterion no longer declines every candidate when a fixed
  shape sits on its boundary and rows are interval-censored (#712).** With a
  left- or interval-censored row the score test falls back to a numeric
  Hessian, and it differentiated every parameter, fixed shapes included. A
  shape fixed at its boundary, such as an early phase's `nu = 0`, made that
  matrix non-finite, so every candidate was declined as `nuisance_singular`
  and `hzr_stepwise()` ran the screen on Wald refits, a different screen from
  SAS's. It now differentiates only the parameters the test reads, which
  changes no value it uses. This changes the result of a score screen on
  such a fit.
