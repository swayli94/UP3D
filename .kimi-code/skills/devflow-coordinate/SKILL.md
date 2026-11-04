---
name: devflow-coordinate
description: Coordinate dependency-ready OpenSpec Changes across multiple agents using one Change per worktree, static roadmap views, and the canonical devflow execution protocol. Use for graph-guided project decomposition, parallel worker assignment, progress refresh, and multi-Change integration coordination.
---

# Devflow Coordination Kimi Code Adapter

Read the canonical guide printed by `devflow protocol coordination` and the protocol printed by `devflow protocol execution` (prefixed by `devflowCommand` when `openspec/project.json` declares one), plus the consuming project's `AGENTS.md`. Those documents are authoritative; this adapter does not duplicate their workflow.

Coordinate at Change granularity: one Change, one worktree, one responsible worker. Agents author the decomposition through OpenSpec; the engine does not generate it or schedule agent processes. Return assignment packets when worker launch is unavailable or unauthorized. Never create a second task store, bypass evidence and approval gates, or merge, push, or archive without explicit authorization.
