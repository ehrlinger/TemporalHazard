#!/bin/bash
# API surface inventory diff between two refs. Deliverable 2 support.
# Usage: bash tools/release-census/api-inventory.sh <old-ref> <new-ref>
# Run from the repo/worktree root. Prints to stdout; no side effects.
set -eu

OLD="${1:?old ref required}"
NEW="${2:?new ref required}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

get_ns() { # ref -> file
  if [ "$2" = "WORKTREE" ]; then cp NAMESPACE "$1"; else git show "$2:NAMESPACE" > "$1"; fi
}

get_ns "$WORK/old.ns" "$OLD"
get_ns "$WORK/new.ns" "$NEW"

for side in old new; do
  grep '^export(' "$WORK/$side.ns" | sed 's/^export(//; s/)$//' | sort > "$WORK/$side.exp"
  grep '^S3method(' "$WORK/$side.ns" | sed 's/^S3method(//; s/)$//' | sort > "$WORK/$side.s3"
  # Normalise imports: roxygen writes importFrom() either one-per-line or as a
  # multi-line list, so compare the (package, symbol) pairs, not the layout.
  tr -d ' \t' < "$WORK/$side.ns" | tr '\n' ' ' | tr -s ' ' \
    | sed 's/) */)\n/g' | grep -E '^(import|useDynLib|exportPattern|exportClass|exportMethod)' \
    | sed 's/,$//' | tr -d ' ' \
    | awk -F'[(,)]' '{ n = 0; for (i = 3; i <= NF; i++) if ($i != "") { print $1 "(" $2 "," $i ")"; n++ } if (n == 0) print $1 "(" $2 ")" }' \
    | sort -u > "$WORK/$side.other"
done

echo "== counts =="
printf "%s: exports=%s S3methods=%s\n" "$OLD" "$(wc -l < "$WORK/old.exp" | tr -d ' ')" "$(wc -l < "$WORK/old.s3" | tr -d ' ')"
printf "%s: exports=%s S3methods=%s\n" "$NEW" "$(wc -l < "$WORK/new.exp" | tr -d ' ')" "$(wc -l < "$WORK/new.s3" | tr -d ' ')"

# KNOWN POSITIVE for the comm logic: a name present only on one side must show up.
# We prove the comparison can report a removal by feeding it a synthetic sentinel.
echo "== known positive (comm can detect a removal) =="
cp "$WORK/old.exp" "$WORK/kp.exp"
echo "zzz_sentinel_export_that_never_existed" >> "$WORK/kp.exp"
sort -o "$WORK/kp.exp" "$WORK/kp.exp"
KP="$(comm -23 "$WORK/kp.exp" "$WORK/new.exp")"
if printf '%s\n' "$KP" | grep -qx 'zzz_sentinel_export_that_never_existed'; then
  echo "PASS: removal detector fired on the sentinel"
else
  echo "FAIL: removal detector did not fire; got: [$KP]" >&2
  exit 1
fi

echo "== exports REMOVED ($OLD -> $NEW) =="
comm -23 "$WORK/old.exp" "$WORK/new.exp" | sed 's/^/  -/'
echo "(end)"
echo "== exports ADDED ($OLD -> $NEW) =="
comm -13 "$WORK/old.exp" "$WORK/new.exp" | sed 's/^/  +/'
echo "(end)"
echo "== S3 methods REMOVED =="
comm -23 "$WORK/old.s3" "$WORK/new.s3" | sed 's/^/  -/'
echo "(end)"
echo "== S3 methods ADDED =="
comm -13 "$WORK/old.s3" "$WORK/new.s3" | sed 's/^/  +/'
echo "(end)"
echo "== other NAMESPACE directives (import/useDynLib/etc) =="
if diff "$WORK/old.other" "$WORK/new.other"; then echo "  (identical)"; fi
echo "(end)"
