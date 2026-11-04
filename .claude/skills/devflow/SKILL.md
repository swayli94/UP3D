---
name: devflow
description: Execute one change with evidence gates, recovery, review, handoff and archive controls.
---

# Devflow Claude Code Adapter

Read the canonical execution protocol printed by `devflow protocol execution` (prefixed by `devflowCommand` when `openspec/project.json` declares one), which matches the installed engine, and the project's governance in `AGENTS.md`. Those documents are authoritative; this adapter does not duplicate them.

Use the deterministic `devflow` CLI for state transitions and the project configuration in `openspec/` for actual commands, cases, checks, evidence and resources. Do not claim the formal skill is unavailable when the adapter and CLI are present.
