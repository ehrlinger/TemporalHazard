#!/bin/bash
# Behavioural census, end to end: install two refs into two private libraries,
# run the fixed battery under each, compare.
#
# Usage:
#   bash tools/release-census/run-census.sh <old-ref> <old-ver> <new-ref> <new-ver> <workdir>
#
# Example (1.3.1 cut, #683):
#   bash tools/release-census/run-census.sh v1.3.0 1.3.0 \
#     design/limit-restart-683 1.3.0 "$WORK"
#
# Every path is keyed by ROLE (old/new), never by version. Two refs can carry
# the same DESCRIPTION Version (v1.3.0 and a 1.3.1 branch before its bump
# both read 1.3.0); version-keyed paths made them share one library, so the
# "old" run executed the new code and the census compared new code with
# itself. Each install is stamped with its commit sha and the runs and the
# comparison refuse a pair whose shas are equal.
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
LIB_OLD="$WORK/lib-old"
LIB_NEW="$WORK/lib-new"
LOG="$WORK/census.log"

OLD_SHA="$(git rev-parse --verify "$OLD_REF^{commit}")"
NEW_SHA="$(git rev-parse --verify "$NEW_REF^{commit}")"

# Same commit on both sides: the census would compare code with itself and
# report "no change" over nothing.
if [ "$OLD_SHA" = "$NEW_SHA" ]; then
  echo "REFUSING: old ($OLD_REF) and new ($NEW_REF) are the same commit $OLD_SHA" >&2
  exit 1
fi
# A reused work dir would let a stale library or RDS from another pair stand
# in for this one; insist on a fresh one.
for p in "$LIB_OLD" "$LIB_NEW" "$WORK/census-old.rds" "$WORK/census-new.rds"; do
  if [ -e "$p" ]; then
    echo "REFUSING: $p already exists; use a fresh work dir" >&2
    exit 1
  fi
done

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
  if [ -n "$ACQ_PID" ]; then
    kill "$ACQ_PID" 2> /dev/null || true
    wait "$ACQ_PID" 2> /dev/null || true
    # The acquire may have finished just before the kill: if it printed a
    # slot, that slot is ours and the EXIT trap must release it.
    if [ -z "$SLOT" ] && [ -s "$WORK/slot.path" ]; then
      SLOT="$(cat "$WORK/slot.path")"
    fi
  fi
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

# The installed packages must carry DIFFERENT build shas. This reads what
# was installed, not what was asked for, so a library that one install wrote
# twice is caught here and not after the comparison has printed "identical".
stamp() { sed -n 's/^CensusSha: *//p' "$1/TemporalHazard/DESCRIPTION" | head -1; }
GOT_OLD="$(stamp "$LIB_OLD")"
GOT_NEW="$(stamp "$LIB_NEW")"
if [ -z "$GOT_OLD" ] || [ -z "$GOT_NEW" ] || [ "$GOT_OLD" = "$GOT_NEW" ]; then
  echo "REFUSING: installed build shas old='$GOT_OLD' new='$GOT_NEW'" | tee -a "$LOG" >&2
  exit 1
fi
echo "installed build shas: old $GOT_OLD, new $GOT_NEW" | tee -a "$LOG"

# --- runs ------------------------------------------------------------------
# The SAME cases file under both, so a case cannot drift between the two.
CASES="$HERE/census-cases.R"

Rscript "$HERE/census-run.R" "$LIB_OLD" "$OLD_VER" "$CASES" \
  "$WORK/census-old.rds" "$OLD_SHA" 2>&1 | tee "$WORK/run-old.log" | tee -a "$LOG"

Rscript "$HERE/census-run.R" "$LIB_NEW" "$NEW_VER" "$CASES" \
  "$WORK/census-new.rds" "$NEW_SHA" 2>&1 | tee "$WORK/run-new.log" | tee -a "$LOG"

# --- compare ---------------------------------------------------------------
# Exits non-zero when a planted change went undetected. No `|| true`: an
# empty or failed comparison is a broken harness, not a pass.
set +e
Rscript "$HERE/census-compare.R" \
  "$WORK/census-old.rds" "$WORK/census-new.rds" \
  "$WORK/census-report.txt" "$CASES" 2>&1 | tee -a "$LOG"
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
