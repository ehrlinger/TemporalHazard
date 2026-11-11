#!/bin/bash
# Heavy-gate slot acquisition, per ~/.cache/th-gate/SLOTS.md.
#
#   bash gate-slot.sh acquire <slug> <description>   -> prints the slot dir held
#   bash gate-slot.sh release <slotdir>
#
# The gate root is TH_GATE_ROOT, defaulting to this machine's shared root,
# ~/.cache/th-gate/slots, which every session here uses. The old default,
# /private/tmp/claude-504/th-heavy-gate-slots, was wiped by the 2026-09-22
# reboot; a job that still used it took a slot in an empty queue nobody else
# read, and so ran outside the two-job limit. Set TH_GATE_ROOT on any other
# machine.
#
# FIFO by ticket. Polls every 120 s and never shortens the interval: the
# ticket order IS the fairness mechanism, and racing it is what the queue
# was added to stop. Never removes another session's ticket or slot.
set -eu

# The root is configurable (TH_GATE_ROOT); the default is this programme's.
GATE="${TH_GATE_ROOT:-$HOME/.cache/th-gate/slots}"
QUEUE="$GATE/queue"
SLOTS="$GATE"

cmd="${1:?acquire|release}"

if [ "$cmd" = "release" ]; then
  SLOT="${2:?slot dir required}"
  # Only this gate's slots, and only one this tool acquired: releasing by
  # path alone removed a slot another job held (2026-09-17).
  case "$SLOT" in
    "$SLOTS/slot1"|"$SLOTS/slot2") ;;
    *) echo "refusing to release $SLOT: not a slot of $SLOTS" >&2; exit 1 ;;
  esac
  if ! grep -qx 'owner: release-census' "$SLOT/owner.txt" 2>/dev/null; then
    echo "refusing to release $SLOT: not held by release-census" >&2
    exit 1
  fi
  rm -rf "$SLOT"
  echo "released $SLOT"
  exit 0
fi

SLUG="${2:?slug required}"
DESC="${3:?description required}"
HEADSHA="$(git rev-parse HEAD)"
WT="$(pwd)"

mkdir -p "$QUEUE"
# Microsecond arrival time, so two tickets in the same second still sort in
# arrival order (the queue sorts on this first field, numerically); the pid
# keeps the name unique.
T="$(perl -MTime::HiRes=time -e 'printf "%.6f", time')-$SLUG-$$"
echo "release-census | $WT | $HEADSHA | $DESC" > "$QUEUE/$T"
echo "took ticket $T" >&2

# Drop the ticket on any exit path that is not a successful acquire.
ACQUIRED=""
cleanup() {
  if [ -z "$ACQUIRED" ]; then rm -f "$QUEUE/$T"; fi
}
trap cleanup EXIT

i=0
while :; do
  i=$((i + 1))
  HEAD_TICKET="$(ls "$QUEUE" | sort -t- -k1,1n | head -1)"
  if [ "$HEAD_TICKET" = "$T" ]; then
    for s in slot1 slot2; do
      # From mkdir to ACQUIRED the slot exists but nothing would release it:
      # the EXIT trap removes only the ticket, and `release` refuses a slot
      # without a matching owner.txt. An interrupt in that window left the
      # slot held until someone removed it by hand. So a signal there is
      # RECORDED, not acted on, and handled once the slot is fully ours:
      # the slot is removed and the acquire exits non-zero.
      SIGNALLED=""
      trap 'SIGNALLED=1' INT TERM HUP
      if mkdir "$SLOTS/$s" 2>/dev/null; then
        {
          echo "owner: release-census"
          echo "worktree: $WT"
          echo "head: $HEADSHA"
          echo "job: $DESC"
          echo "start: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
          echo "start_epoch: $(date +%s)"
        } > "$SLOTS/$s/owner.txt"
        rm -f "$QUEUE/$T"
        ACQUIRED="$SLOTS/$s"
        trap - INT TERM HUP
        if [ -n "$SIGNALLED" ]; then
          rm -rf "$ACQUIRED"
          echo "signalled while acquiring; released $ACQUIRED" >&2
          exit 130
        fi
        echo "acquired $SLOTS/$s after $i poll(s)" >&2
        echo "$ACQUIRED"
        exit 0
      fi
      trap - INT TERM HUP
      if [ -n "$SIGNALLED" ]; then exit 130; fi
    done
    echo "poll $i: head of queue, both slots held; waiting" >&2
  else
    echo "poll $i: waiting behind $HEAD_TICKET" >&2
  fi
  # No GNU timeout on this machine; sleep is the only timer, and 120 s is
  # the interval the README fixes.
  sleep 120
done
