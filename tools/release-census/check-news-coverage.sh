#!/bin/bash
# Does NEWS.md explain the upgrade, read as a user of <old-ref> would read it?
#
# Usage: bash tools/release-census/check-news-coverage.sh <old-ref> <old-version>
#
# Method: every issue number a commit between <old-ref> and HEAD claims to
# close is looked for in the NEWS text ABOVE the <old-version> heading -- that
# is, in the sections a user upgrading from <old-version> would read. An issue
# that was closed by a commit and is named nowhere in that text is a
# CANDIDATE finding: either the change was not user-visible (fine, and say
# so), or the upgrade is under-described (a finding).
#
# This is a candidate list, not a verdict. Each hit must be read against the
# commit before it is called a finding: a refactor or a test-only change
# legitimately has no NEWS entry.
set -eu

OLD="${1:?old ref required}"
OLDVER="${2:?old version required}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# The NEWS text a user of $OLDVER reads: everything above that heading.
awk -v v="# TemporalHazard $OLDVER" '$0 == v { exit } { print }' NEWS.md > "$WORK/news-since.txt"
NEWS_LINES="$(wc -l < "$WORK/news-since.txt" | tr -d ' ')"
echo "NEWS text above the '$OLDVER' heading: $NEWS_LINES lines"

# Issue numbers NEWS cites in that text.
grep -o -E '#[0-9]+' "$WORK/news-since.txt" | sort -u > "$WORK/news-issues.txt"
echo "distinct issue numbers cited in that NEWS text: $(wc -l < "$WORK/news-issues.txt" | tr -d ' ')"

# Issue numbers the commits claim to close. Merge-commit PR numbers are
# EXCLUDED on purpose: "Merge pull request #350" is a PR number, and NEWS
# cites issues, so counting them would fabricate dozens of false misses.
git log --no-merges --format='%H%x09%s%x09%b' "$OLD..HEAD" > "$WORK/commits.tsv"
echo "non-merge commits in $OLD..HEAD: $(wc -l < "$WORK/commits.tsv" | tr -d ' ')"

grep -o -i -E '(closes|close|closed|fixes|fix|fixed|resolves|resolve) +#[0-9]+' "$WORK/commits.tsv" \
  | grep -o -E '#[0-9]+' | sort -u > "$WORK/commit-issues.txt"
echo "distinct issues claimed closed by commits: $(wc -l < "$WORK/commit-issues.txt" | tr -d ' ')"
echo

# KNOWN POSITIVE: the comparison must report an issue that NEWS does not
# mention. A sentinel number that cannot be in NEWS is injected and must
# come out the other side.
cp "$WORK/commit-issues.txt" "$WORK/kp.txt"
echo "#9999999" >> "$WORK/kp.txt"
sort -o "$WORK/kp.txt" "$WORK/kp.txt"
if ! comm -23 "$WORK/kp.txt" "$WORK/news-issues.txt" | grep -q '^#9999999$'; then
  echo "KNOWN POSITIVE FAILED: the missing-from-NEWS detector did not fire" >&2
  exit 1
fi
echo "known positive: missing-from-NEWS detector fires on a sentinel -- PASS"
echo

echo "=== issues closed since $OLD but NOT cited in the NEWS a $OLDVER user reads ==="
MISSING="$(comm -23 "$WORK/commit-issues.txt" "$WORK/news-issues.txt")"
if [ -z "$MISSING" ]; then
  echo "  (none)"
else
  for iss in $MISSING; do
    n="${iss#\#}"
    echo ""
    echo "  $iss -- commits claiming to close it:"
    grep -i -E "(closes|close|closed|fixes|fix|fixed|resolves|resolve) +#$n\b" "$WORK/commits.tsv" \
      | cut -f1,2 | sed 's/^/    /' | cut -c1-140
  done
fi
echo ""
echo "=== issues cited in NEWS with no commit claiming to close them ==="
echo "(normal: NEWS often cites an issue for context rather than as the fix)"
comm -13 "$WORK/commit-issues.txt" "$WORK/news-issues.txt" | tr '\n' ' '
echo ""
