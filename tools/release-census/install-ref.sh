#!/bin/bash
# Install one git ref of TemporalHazard into a private library tree.
#
# Usage: bash tools/release-census/install-ref.sh <ref> <libdir> <expected-version>
#
# The tree is a clean `git archive` export, not the working tree: AGENTS.md
# requires that, and a worktree's `.git` is a FILE, so a working-tree build
# drags it into the tarball. Nothing is committed or checked out here.
set -eu

REF="${1:?ref required}"
LIB="${2:?library dir required}"
EXPECT="${3:?expected version required}"

mkdir -p "$LIB"
SRC="$(mktemp -d)"
trap 'rm -rf "$SRC"' EXIT

mkdir -p "$SRC/tree"
git archive "$REF" | tar -x -C "$SRC/tree"

GOT="$(sed -n 's/^Version: *//p' "$SRC/tree/DESCRIPTION" | head -1)"
if [ "$GOT" != "$EXPECT" ]; then
  echo "REFUSING: $REF has DESCRIPTION Version '$GOT', expected '$EXPECT'" >&2
  exit 1
fi
echo "exported $REF at version $GOT"

# --no-docs skips the help database; the census calls functions, not ?fn.
R CMD INSTALL --no-docs --library="$LIB" "$SRC/tree"

# Prove the install landed where we asked and at the version we asked for.
# An install that silently went elsewhere would leave the census reading the
# stale 1.1.0 copy in the system library.
if [ ! -f "$LIB/TemporalHazard/DESCRIPTION" ]; then
  echo "REFUSING: no TemporalHazard/DESCRIPTION under $LIB after install" >&2
  exit 1
fi
INSTALLED="$(sed -n 's/^Version: *//p' "$LIB/TemporalHazard/DESCRIPTION" | head -1)"
if [ "$INSTALLED" != "$EXPECT" ]; then
  echo "REFUSING: installed version '$INSTALLED' != expected '$EXPECT'" >&2
  exit 1
fi
echo "INSTALL OK: $REF -> $LIB (version $INSTALLED)"
