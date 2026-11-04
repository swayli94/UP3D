---
name: openspec-apply-change
description: Implement an OpenSpec change through the project's devflow execution gates. Use for "openspec apply", "opsx apply", or continuing implementation.
allowed-tools: Bash
license: MIT
compatibility: Requires OpenSpec CLI and the project's devflow adapter.
metadata:
  author: openspec
  version: "1.1"
---

# Apply through devflow

This is an adapter, not a second implementation workflow. Read `AGENTS.md`, the protocol printed by `devflow protocol execution`, the active Change artifacts, and `openspec/project.json` before acting.

1. Select the named Change. If it is missing or ambiguous, list active Changes and ask the user.
2. Read the current specifications and all Change artifacts. Use OpenSpec status, instructions, show, and validate only for their read-only specification duties.
3. Enter the declared branch and worktree with the project's devflow CLI: `devflow start --change <id>`, or `devflow resume --change <id>` for existing work. If `openspec/project.json` declares `devflowCommand` (for example `["conda", "run", "-n", "env", "devflow"]`), use it in place of `devflow`.
4. Select each task with `devflow run --change <id> --task <task-id>`. Respect dependencies, resource ownership, and the declared retry and failure targets.
5. Implement only that action. Execute its checks immediately with `devflow verify --change <id> --task <task-id>`, which reuses Runs that still apply to the current content; use `--checks` only to force a rerun. Preserve the resulting Runs.
6. Mark completion only through `devflow complete --change <id> --task <task-id>`, supplying a Run, operator evidence, or checkpoint approval when required. Never edit a task checkbox directly.
7. On failure, preserve output and use `devflow fail` to reopen the declared recovery task. Stop at the retry limit or whenever project guidance requires a user decision.
8. Use `devflow review` for independent review and `devflow handoff` to persist complete cross-machine state.

Do not use native apply instructions to edit task markers, and do not invoke native OpenSpec archive. These entries do not bypass Run, independent review, user approval, or archive gates. Report implementation, verification, review, approval, merge, and archive as separate states.
