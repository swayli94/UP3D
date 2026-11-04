---
name: openspec-archive-change
description: Archive an OpenSpec change only after the project's devflow verification, review, and explicit user approval gates pass. Use for "openspec archive" or "opsx archive".
allowed-tools: Bash
license: MIT
compatibility: Requires OpenSpec CLI and the project's devflow adapter.
metadata:
  author: openspec
  version: "1.1"
---

# Archive through devflow

This is an adapter, not a second archive workflow. Read `AGENTS.md`, the protocol printed by `devflow protocol execution`, and the active Change artifacts.

1. Resolve the named active Change; do not infer an archived Change or ask to archive an ambiguous one.
2. Inspect `devflow resume --change <id>`. If `openspec/project.json` declares `devflowCommand` (for example `["conda", "run", "-n", "env", "devflow"]`), use it in place of `devflow`.
3. If tasks, current-platform and required-platform verification, independent review, or specification validation are incomplete, stop and report the exact blocker. Do not narrow the required profile or treat an old Run or review as current.
4. Present the current commit, verification scope and evidence, independent review, specification sync impact, and archive destination to the user.
5. Proceed only after the user explicitly approves this archive. Invoke `devflow archive --change <id> --approval "<user's explicit approval>"`.
6. Report the actual archive directory and the regenerated roadmap view separately from any requested merge. A successful archive does not imply user approval to merge or push main.

Never run native `openspec archive`, move a Change directory manually, edit task markers, or use `--no-validate` or automatic confirmation to bypass this adapter. The devflow adapter alone delegates the final validated OpenSpec archive after all project gates pass.
