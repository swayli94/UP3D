---
name: devflow-refactor
description: Assess, adopt, plan, or execute behavior-preserving refactoring in an existing brownfield or legacy repository using characterization baselines, bounded OpenSpec Changes, and devflow evidence gates. Use for restructuring, module extraction, legacy modernization, dependency replacement, or requests to refactor an already developed repository.
---

# Devflow Refactoring Claude Code Adapter

Read the canonical guide printed by `devflow protocol refactoring` and the protocol printed by `devflow protocol execution` (prefixed by `devflowCommand` when `openspec/project.json` declares one), plus the consuming project's `AGENTS.md`. Those documents are authoritative; this adapter does not duplicate their workflow.

Default to read-only assessment for a broad or unspecified request. Use `openspec-propose` for one bounded planning slice and `openspec-apply-change` only after a separate explicit apply request. Use the deterministic `devflow` CLI for execution and state transitions. Never turn this skill into a whole-repository rewrite, a second checkbox store, or a bypass of specification, evidence, review, or approval gates.
