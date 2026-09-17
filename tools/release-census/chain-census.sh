#!/bin/bash
# Run a SECOND census only after a first one has finished and released its
# heavy-gate slot.
#
# Usage:
#   bash tools/release-census/chain-census.sh <first-workdir> <first-pid> \
#        <old-ref> <old-ver> <new-ref> <new-ver> <second-workdir>
#
# WHY A CHAIN AND NOT A SECOND TICKET NOW: the gate README is explicit --
# "One ticket per JOB, not per session. Take the next job's ticket only when
# the previous job has finished and its slot is released, so a chained runner
# holds exactly one ticket." Taking a second ticket while the first job is
# still queued would put one session twice in a FIFO built to stop exactly
# that. So this waits, then run-census.sh takes a fresh ticket at the real
# epoch of that moment -- which is also strictly behind every ticket taken
# while we were waiting, so it cannot jump anyone.
#
# It also guarantees the two runs are never concurrent, which matters beyond
# fairness: two of our own heavy R jobs on a loaded machine only slow each
# other, and neither result would be worth a timing anyway (we record none).
set -eu

HERE="$(cd "$(dirname "$0")" && pwd)"

FIRST_WORK="${1:?first workdir}"
FIRST_PID="${2:?first pid}"
OLD_REF="${3:?old ref}"
OLD_VER="${4:?old version}"
NEW_REF="${5:?new ref}"
NEW_VER="${6:?new version}"
SECOND_WORK="${7:?second workdir}"

mkdir -p "$SECOND_WORK"
CHAINLOG="$SECOND_WORK/chain.log"

say() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $*" | tee -a "$CHAINLOG"; }

say "waiting for the first census (pid $FIRST_PID) to finish"

# Wait on BOTH conditions, not either: the sentinel line proves the compare
# ran, and the dead pid proves the trap that releases the slot has fired.
# Keying on the log alone could start the second run while the first still
# held its slot.
while :; do
  DONE_LOG=0
  DONE_PID=0
  if grep -q "census-compare exit status" "$FIRST_WORK/orchestrator.log" 2>/dev/null; then
    DONE_LOG=1
  fi
  if ! kill -0 "$FIRST_PID" 2>/dev/null; then
    DONE_PID=1
  fi
  if [ "$DONE_LOG" = 1 ] && [ "$DONE_PID" = 1 ]; then
    say "first census finished and its process is gone"
    break
  fi
  # A dead process with no sentinel means the first run died without
  # comparing. Say so loudly rather than chaining on top of a broken run.
  if [ "$DONE_PID" = 1 ] && [ "$DONE_LOG" = 0 ]; then
    say "ABORT: first census process $FIRST_PID is gone but never wrote a"
    say "ABORT: 'census-compare exit status' line. It failed before comparing."
    say "ABORT: not chaining a second run on top of that."
    exit 1
  fi
  sleep 30
done

# Belt and braces: do not start while a slot still names this session.
for s in /private/tmp/claude-504/th-heavy-gate-slots/slot1 \
         /private/tmp/claude-504/th-heavy-gate-slots/slot2; do
  if [ -d "$s" ] && grep -q "release-census" "$s/owner.txt" 2>/dev/null; then
    say "a gate slot still names release-census; waiting for its release"
    while [ -d "$s" ] && grep -q "release-census" "$s/owner.txt" 2>/dev/null; do
      sleep 30
    done
    say "slot released"
  fi
done

say "starting the $OLD_VER-baseline census (it takes its own ticket now)"
bash "$HERE/run-census.sh" "$OLD_REF" "$OLD_VER" "$NEW_REF" "$NEW_VER" \
  "$SECOND_WORK" 2>&1 | tee -a "$CHAINLOG"
RC=$?
say "second census orchestrator exit: $RC"
exit "$RC"
