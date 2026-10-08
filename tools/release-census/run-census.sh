#!/bin/bash
# Behavioural census, end to end: install two refs into two private libraries,
# run the fixed battery under each, compare.
#
# Usage:
#   bash tools/release-census/run-census.sh <old-ref> <old-ver> <new-ref> <new-ver> <workdir>
#
# Example (1.3.0 release gate):
#   bash tools/release-census/run-census.sh v1.2.9 1.2.9 origin/main 1.2.11 "$WORK"
#
# Holds a heavy-gate slot for the whole run and releases it on EVERY exit
# path, including failure. No timing is recorded: six sessions share this
# CPU, so every wall-clock number from this programme is contended.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

OLD_REF="${1:?old ref}"
OLD_VER="${2:?old version}"
NEW_REF="${3:?new ref}"
NEW_VER="${4:?new version}"
WORK="${5:?work dir}"

mkdir -p "$WORK"
LIB_OLD="$WORK/lib-$OLD_VER"
LIB_NEW="$WORK/lib-$NEW_VER"
LOG="$WORK/census.log"

OLD_SHA="$(git rev-parse "$OLD_REF")"
NEW_SHA="$(git rev-parse "$NEW_REF")"

{
  echo "=== behavioural census ==="
  echo "old: $OLD_REF = $OLD_SHA (expect version $OLD_VER)"
  echo "new: $NEW_REF = $NEW_SHA (expect version $NEW_VER)"
  echo "work: $WORK"
} | tee "$LOG"

# --- heavy gate ------------------------------------------------------------
SLOT=""
release_slot() {
  if [ -n "$SLOT" ] && [ -d "$SLOT" ]; then
    bash "$HERE/gate-slot.sh" release "$SLOT" | tee -a "$LOG"
  fi
}
# Release on exit; a signal exits, so the EXIT trap releases once and the
# census does not carry on without its slot.
trap release_slot EXIT
# While a slot is being waited for, the acquiring child is a separate
# process: a signal must stop it too, or it could later take a slot that no
# one releases. Its own EXIT trap drops its ticket.
ACQ_PID=""
on_signal() {
  if [ -n "$ACQ_PID" ]; then kill "$ACQ_PID" 2> /dev/null || true; fi
  exit 130
}
trap on_signal INT TERM

bash "$HERE/gate-slot.sh" acquire release-census "census $OLD_VER vs $NEW_VER" \
  > "$WORK/slot.path" &
ACQ_PID=$!
wait "$ACQ_PID"
ACQ_PID=""
SLOT="$(cat "$WORK/slot.path")"
echo "holding gate slot: $SLOT" | tee -a "$LOG"

# --- installs --------------------------------------------------------------
# Each ref goes into its OWN library so the two versions cannot see each
# other. census-run.R then refuses to proceed unless the package it loaded
# came from the library it was handed, at the version it was promised.
bash "$HERE/install-ref.sh" "$OLD_SHA" "$LIB_OLD" "$OLD_VER" 2>&1 | tee -a "$LOG"
bash "$HERE/install-ref.sh" "$NEW_SHA" "$LIB_NEW" "$NEW_VER" 2>&1 | tee -a "$LOG"

# --- runs ------------------------------------------------------------------
# The SAME cases file under both, so a case cannot drift between the two.
CASES="$HERE/census-cases.R"

Rscript "$HERE/census-run.R" "$LIB_OLD" "$OLD_VER" "$CASES" \
  "$WORK/census-$OLD_VER.rds" 2>&1 | tee "$WORK/run-$OLD_VER.log" | tee -a "$LOG"

Rscript "$HERE/census-run.R" "$LIB_NEW" "$NEW_VER" "$CASES" \
  "$WORK/census-$NEW_VER.rds" 2>&1 | tee "$WORK/run-$NEW_VER.log" | tee -a "$LOG"

# --- compare ---------------------------------------------------------------
# Exits non-zero when a planted change went undetected. No `|| true`: an
# empty or failed comparison is a broken harness, not a pass.
set +e
Rscript "$HERE/census-compare.R" \
  "$WORK/census-$OLD_VER.rds" "$WORK/census-$NEW_VER.rds" \
  "$WORK/census-report.txt" 2>&1 | tee -a "$LOG"
CMP=${PIPESTATUS[0]}
set -e

{
  echo ""
  echo "=== census-compare exit status: $CMP ==="
  echo "old sha: $OLD_SHA"
  echo "new sha: $NEW_SHA"
  echo "report : $WORK/census-report.txt"
} | tee -a "$LOG"

exit "$CMP"
