#!/bin/bash
# CRAN Cookbook spot-checks on the exports ADDED since a reference ref.
#
# Usage: bash tools/release-census/check-new-exports.sh <old-ref>
#
# Scope is deliberately narrow: only what CHANGED. A full-package Cookbook
# audit is the release session's job, not this one's.
#
# Checks, per added export:
#   \value present            -- CRAN rejects an exported object without one
#   \examples present         -- and whether they are wrapped in \dontrun{}
#   \dontrun vs \donttest     -- \dontrun is the one CRAN objects to
#   \title in Title Case
#   software names quoted     -- 'SAS', 'C', 'R' etc in \title/\description
set -eu

OLD="${1:?old ref required}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

git show "$OLD:NAMESPACE" | grep '^export(' | sed 's/^export(//; s/)$//' | sort > "$WORK/old.exp"
grep '^export(' NAMESPACE | sed 's/^export(//; s/)$//' | sort > "$WORK/new.exp"
ADDED="$(comm -13 "$WORK/old.exp" "$WORK/new.exp")"

if [ -z "$ADDED" ]; then
  echo "no exports added since $OLD -- nothing to check"
  exit 0
fi

echo "exports added since $OLD:"
echo "$ADDED" | sed 's/^/  +/'
echo

# Known positive: the \value detector must fire on a file that has none.
printf '\\name{zzz}\n\\title{No Value Here}\n' > "$WORK/novalue.Rd"
if grep -q '\\value' "$WORK/novalue.Rd"; then
  echo "KNOWN POSITIVE FAILED: \\value detector fires on a file without one" >&2
  exit 1
fi
printf '\\name{zzz}\n\\value{something}\n' > "$WORK/hasvalue.Rd"
if ! grep -q '\\value' "$WORK/hasvalue.Rd"; then
  echo "KNOWN POSITIVE FAILED: \\value detector misses a file that has one" >&2
  exit 1
fi
echo "known positive: \\value detector distinguishes both cases -- PASS"
echo

FAILED=0
for fn in $ADDED; do
  echo "=== $fn ==="
  # Find the Rd that documents it: \name or an \alias.
  RD="$(grep -l -e "^\\\\alias{$fn}" -e "^\\\\name{$fn}" man/*.Rd 2>/dev/null | head -1 || true)"
  if [ -z "$RD" ]; then
    echo "  FAIL: no man/*.Rd documents $fn"
    FAILED=$((FAILED + 1))
    continue
  fi
  echo "  Rd: $RD"

  if grep -q '^\\value' "$RD"; then
    echo "  \\value        : present"
  else
    echo "  \\value        : MISSING  <-- CRAN rejects this"
    FAILED=$((FAILED + 1))
  fi

  if grep -q '^\\examples' "$RD"; then
    if grep -q '\\dontrun' "$RD"; then
      echo "  \\examples     : present, but uses \\dontrun  <-- prefer \\donttest"
      FAILED=$((FAILED + 1))
    elif grep -q '\\donttest' "$RD"; then
      echo "  \\examples     : present, \\donttest (runnable, not run by CRAN)"
    else
      echo "  \\examples     : present and run by CRAN"
    fi
    # A Suggests used in an example must be guarded.
    if grep -qE 'ggplot2|haven|withr|numDeriv|scales|knitr' "$RD"; then
      if grep -q 'requireNamespace' "$RD"; then
        echo "  Suggests use  : guarded by requireNamespace()"
      else
        echo "  Suggests use  : a Suggests package appears UNGUARDED in the Rd"
        FAILED=$((FAILED + 1))
      fi
    fi
  else
    echo "  \\examples     : MISSING"
    FAILED=$((FAILED + 1))
  fi

  TITLE="$(sed -n 's/^\\title{\(.*\)}$/\1/p' "$RD" | head -1)"
  echo "  \\title        : $TITLE"
  # Title Case check: flag a lower-case word that is not a common minor word.
  BAD="$(echo "$TITLE" | tr ' ' '\n' \
        | grep -E '^[a-z]' \
        | grep -v -E '^(a|an|the|and|or|but|for|nor|of|to|in|on|at|by|from|with|as|per|via|into|over)$' \
        | tr '\n' ' ' || true)"
  if [ -n "$BAD" ]; then
    echo "  Title Case    : check these lower-case words: $BAD"
  else
    echo "  Title Case    : ok"
  fi

  # Software names must be quoted in DESCRIPTION; in Rd it is a house-style
  # consistency point rather than a CRAN rule, so this reports, not fails.
  UNQ="$(grep -o -E '(^|[^'"'"'])\b(SAS|PROC HAZARD|HAZARD)\b' "$RD" | head -3 | tr '\n' ' ' || true)"
  if [ -n "$UNQ" ]; then
    echo "  software names: occurrences to eyeball for quoting: $UNQ"
  fi
  echo
done

echo "checks failing: $FAILED"
# Non-zero when any check failed, so the caller can gate on it. Capped at 1:
# an exit status is taken modulo 256, so a count could wrap to "success".
[ "$FAILED" -eq 0 ]
