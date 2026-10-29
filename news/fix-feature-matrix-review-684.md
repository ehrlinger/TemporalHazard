* The feature test matrix (#682) now reads its Reference column from
  comparisons cited by name, so tests that only check a reference file is
  well formed are no longer counted. Its tests column counts only tests that
  ran. The generator refuses a run saved without provenance, a run of another
  commit or dependency set, and a run whose generated test names differ from
  the tree (#684 review).
