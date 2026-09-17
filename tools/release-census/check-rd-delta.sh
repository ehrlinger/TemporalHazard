#!/bin/bash
# CRAN Manuals-chapter checks on the Rd files that CHANGED since a ref.
#
# Usage: bash tools/release-census/check-rd-delta.sh <old-ref>
#
# Raw Unicode in Rd is the defect only the PDF-manual step catches, and
# `R CMD check --no-manual` skips it silently -- so a text-level sweep of the
# changed Rd is worth having before the manual build is run.
set -eu

# `grep` on this machine is ugrep, which does NOT read [^\x00-\x7F] as a
# byte class -- it matched pure ASCII, and the known positive below caught
# that. perl is used instead, and the known positive stays as the guard.
# Prints the offending lines AND exits non-zero when there are none, because
# `perl -ne` exits 0 whether or not it printed anything -- testing its exit
# status alone reported a hit on every file, which the known positive caught.
nonascii() {
  out="$(perl -ne 'print "$.:$_" if /[^\x00-\x7F]/' "$1")"
  [ -n "$out" ] || return 1
  printf '%s\n' "$out"
}

OLD="${1:?old ref required}"

CHANGED="$(git diff --name-only "$OLD..HEAD" -- man/ | grep '\.Rd$' || true)"
ALLRD="$(ls man/*.Rd | wc -l | tr -d ' ')"
NCH="$(echo "$CHANGED" | grep -c . || true)"
echo "Rd files in man/: $ALLRD; changed since $OLD: $NCH"
echo

# KNOWN POSITIVE: the non-ASCII detector must fire on a known non-ASCII byte.
KP="$(mktemp)"
printf 'beta-hat: \xce\xb2\n' > "$KP"
if ! nonascii "$KP" > /dev/null; then
  echo "KNOWN POSITIVE FAILED: non-ASCII detector missed a Greek beta" >&2
  rm -f "$KP"
  exit 1
fi
printf 'plain ascii only\n' > "$KP"
if nonascii "$KP" > /dev/null; then
  echo "KNOWN POSITIVE FAILED: non-ASCII detector fired on pure ASCII" >&2
  rm -f "$KP"
  exit 1
fi
rm -f "$KP"
echo "known positive: non-ASCII detector distinguishes both cases -- PASS"
echo

if [ -z "$CHANGED" ]; then
  echo "no Rd changed; nothing to check"
  exit 0
fi

echo "=== raw non-ASCII in CHANGED Rd (the PDF-manual failure mode) ==="
HITS=0
for f in $CHANGED; do
  [ -f "$f" ] || continue
  if nonascii "$f" > /dev/null; then
    echo "  $f:"
    nonascii "$f" | head -5 | sed 's/^/    /'
    HITS=$((HITS + 1))
  fi
done
[ "$HITS" -eq 0 ] && echo "  (none)"
echo

echo "=== \\dontrun in CHANGED Rd (CRAN prefers \\donttest) ==="
D=0
for f in $CHANGED; do
  [ -f "$f" ] || continue
  if grep -q '\\dontrun' "$f"; then echo "  $f"; D=$((D + 1)); fi
done
[ "$D" -eq 0 ] && echo "  (none)"
echo

echo "=== CHANGED Rd for an exported object with no \\value ==="
V=0
for f in $CHANGED; do
  [ -f "$f" ] || continue
  # data and package docs legitimately differ; report them rather than skip.
  if grep -q '^\\docType{data}' "$f"; then
    grep -q '^\\format' "$f" || { echo "  $f (data doc, no \\format)"; V=$((V + 1)); }
    continue
  fi
  if grep -q '^\\usage' "$f" && ! grep -q '^\\value' "$f"; then
    echo "  $f"
    V=$((V + 1))
  fi
done
[ "$V" -eq 0 ] && echo "  (none)"
echo

# Non-ASCII inside \enc{...}{...} is the CRAN-sanctioned form: it carries an
# ASCII fallback, so it is NOT a defect. man/hzr_log1pexp.Rd and
# man/hzr_log1mexp.Rd both spell Maechler that way and are correct. Only RAW
# non-ASCII is what the PDF-manual step fails on.
echo "=== whole-package sanity: raw non-ASCII in any Rd (\\enc{} is fine) ==="
A=0
for f in man/*.Rd; do
  if nonascii "$f" > /dev/null; then
    if nonascii "$f" | grep -q -F '\enc{'; then
      echo "  $f (inside \\enc{} -- correct, has an ASCII fallback)"
    else
      echo "  $f  <-- RAW non-ASCII: this is the PDF-manual failure mode"
      A=$((A + 1))
    fi
  fi
done
if [ "$A" -eq 0 ]; then echo "  (no raw non-ASCII)"; fi
exit "$A"
