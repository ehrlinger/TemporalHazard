* A test now pins the speed-up from 1.2.14's CoE gradient cache (#631). It
  counts the per-phase work and checks that a call matching the cache key
  does none, and that changing any one of `theta`, the times or the entry
  times alone forces it. The existing test compared gradients, which a
  recomputed value matches too, so it passed with caching switched off; it
  also passed with `theta`, the times or the entry times dropped from the
  key. The new test fails on each of those four.
