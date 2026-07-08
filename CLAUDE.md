# pyFP3D — Project Instructions for Claude Code

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in supersonic
zones, wake cut + Kutta condition for lift. Target: wings at M∞ 0.3–0.87,
workstation-scale (minutes for 1–3 M nodes).

## Document map (read the relevant one before coding)

- [docs/roadmap.md](docs/roadmap.md) — **active tracker**: phase order (Track P:
  P0–P8 solver, Track M: M0–M2 meshing), gate checklists, progress ledger.
  "What phase are we in" and "what gate is open" live here, nowhere else.
- [docs/design.md](docs/design.md) — theory & numerics reference: equations (§2–§3),
  wake/Kutta (§4), BCs (§5), discretization (§6), Numba kernel rules (§7), solver
  strategy (§8), V0–V6 validation ladder (§10), risks/mitigations (§12).
- [docs/demo_report.md](docs/demo_report.md) — **evidence dossier** for completed
  phases (P0, P1-partial, P2, P3, P4, P5, M0, M1): one self-checking demo per phase under
  `cases/demo/<phase>/` with committed figures + measured gate numbers. When a
  phase closes, add its demo + report section here.
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) — layout, per-module status, and
  **"Known gaps"**: read it before touching the G1.6 sphere-Cp problem (formerly
  G1.2; P1 gates renumbered 2026-07-06, mapping in roadmap.md) — it is already
  root-caused, and the G1.3/G1.4 oracles ruled out boundary-data corrections
  (design.md §5.1.2); do not re-propose h-refinement, recovery tweaks, Nitsche,
  or flux-data corrections. Open route: Option C gate re-spec + curved elements.
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
3. Full suite: `pytest tests/` (140 passed + 4 skipped + 2 xfailed since P5, ~5 min;
   the always-on coarse transonic smoke is ~170 s of it, the G3.2 medium-mesh
   nested Picard solve ~45 s, the rule-7 sweep's M6 coarse+medium cut_wake
   ingest ~15 s. The M6 .msh files are gitignored — the 13 M1 tests skip until
   you run `cases/meshes/onera_m6/generate_onera_m6.py` (~30 s). The heavy
   transonic gates only run under `PYFP3D_TRANSONIC_GATES=1`, shown as 4 skipped.)
4. Numba debugging: `PYFP3D_NOJIT=1` swaps `@njit` for identity — print/pdb work.
5. When a gate closes: tick it in roadmap.md, update the progress ledger and the
   "Current phase" line in docs/agent-rules.md, keep the commit phase-scoped.
6. **Cost caution — do not recompute expensive artifacts casually.** Some
   evidence is committed precisely because regenerating it is slow: the P4 heavy
   demo figures (`cases/demo/p4_transonic/run_demo.py` under
   `PYFP3D_TRANSONIC_GATES=1`) are ~40 min of Picard; the medium G4.1 gate is
   ~17 min; the P5 M6 medium from-scratch continuation+polish is ~45–75 min
   (its solution npz `cases/demo/p5_onera_m6/results/medium_solution.npz` is a
   LOCAL gitignored cache like the .msh — the demo re-solves it when absent, or
   under `PYFP3D_P5_RESOLVE=1`; the committed PNG/CSV are the evidence); the
   ONERA M6 medium/fine `.msh` are minutes to regenerate. Treat the
   committed baseline as authoritative and only rerun the heavy part when a real
   solver/mesh/reference change would move those numbers AND you will commit the
   refresh. For routine edits, verify on the cheap coarse path. Prefer reading a
   committed CSV/PNG over recomputing it.

Gate IDs are `G<phase>.<n>` per roadmap.md Track P numbering (P2 = wake/Kutta,
P3 = subsonic compressible, P4 = transonic, P5 = ONERA M6, P6 = consistent/
differentiable artificial-density flux, P7 = Newton/performance);
design.md §11 mirrors the same order.
