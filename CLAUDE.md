@AGENTS.md

# Claude Code specifics

`AGENTS.md`, imported above, is the operational contract and applies in full. It is written to
be tool neutral so that Codex and other agents read the same rules. Only the Claude Code
affordances live here.

## Before you touch code

`AGENTS.md` says to orient on the public API surface before editing. In Claude Code the way to
do that is the codemap: it lives in the Obsidian vault under `Claude/repomaps/` and is read via
the `read-codemap` skill (`/codemap temporal_hazard`). If the codemap looks stale, say so and
offer to refresh it (`/regenerate-codemap`) rather than working from a guess.

## Prose

`AGENTS.md` names the voice and the reader personas. In Claude Code, apply the
`ehrlinger-writing` skill: it carries that voice, the persona menu and the project context,
read from the vault sources. This repo has no `.claude/house-style.md` — the skill is the
only route.
