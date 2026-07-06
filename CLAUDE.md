# pyFP3D — Project Instructions for Claude Code

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in supersonic
zones, wake cut + Kutta condition for lift. Target: wings at M∞ 0.3–0.87,
workstation-scale (minutes for 1–3 M nodes).

## Document map (read the relevant one before coding)

- [docs/roadmap.md](docs/roadmap.md) — **active tracker**: phase order (Track P:
  P0–P7 solver, Track M: M0–M2 meshing), gate checklists, progress ledger.
  "What phase are we in" and "what gate is open" live here, nowhere else.
- [docs/design.md](docs/design.md) — theory & numerics reference: equations (§2–§3),
  wake/Kutta (§4), BCs (§5), discretization (§6), Numba kernel rules (§7), solver
  strategy (§8), V0–V6 validation ladder (§10), risks/mitigations (§12).
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) — layout, per-module status, and
  **"Known gaps"**: read it before touching the G1.6 sphere-Cp problem (formerly
  G1.2; P1 gates renumbered 2026-07-06, mapping in roadmap.md) — it is already
  root-caused (curved-wall/flat-facet variational crime; tiered fix routes in
  design.md §5.1); do not re-propose h-refinement, recovery tweaks, or Nitsche.
- `cases/reference_data/` — ground truth, never edit.

## Hard rules and current phase

@docs/agent-rules.md

## Workflow

1. Before coding: find the open gate in roadmap.md's current phase and plan against
   its acceptance criterion. Every visual gate needs a headless artifact
   (`artifacts/<gate_id>/*.png` + `summary.csv`; matplotlib `Agg`, PyVista
   off-screen — never GUI-only checks).
2. After any kernel or assembly change, run the primary regression first:
   `pytest tests/test_v0_freestream.py`
3. Full suite is fast (~10 s): `pytest tests/`
4. Numba debugging: `PYFP3D_NOJIT=1` swaps `@njit` for identity — print/pdb work.
5. When a gate closes: tick it in roadmap.md, update the progress ledger and the
   "Current phase" line in docs/agent-rules.md, keep the commit phase-scoped.

Gate IDs are `G<phase>.<n>` per roadmap.md Track P numbering (P2 = wake/Kutta,
P3 = subsonic compressible, P4 = transonic, P5 = ONERA M6, P6 = Newton/performance);
design.md §11 mirrors the same order.
