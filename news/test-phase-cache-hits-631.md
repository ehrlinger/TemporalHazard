* Two tests now pin the speed-up from 1.2.14's CoE gradient cache (#631).
  One counts the per-phase work and checks that a call matching the cache
  key does none, and that changing any one of `theta`, the times or the
  entry times alone forces it. The other runs a real CoE fit and checks that
  the second gradient evaluation at each point reuses the first's cache and
  does no work. The existing test compared gradients, which a recomputed
  value matches too, so it passed with caching switched off, with `theta`,
  the times or the entry times dropped from the key, and with the CoE
  wrapper handing its second evaluation a fresh cache. The new tests fail on
  each of those.
