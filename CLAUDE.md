# pyFP3D — Project Instructions for Claude Code

3D unstructured-mesh **full-potential** transonic flow solver (Python + Numba):
one scalar φ per node, Galerkin P1 tets, artificial-density upwinding in supersonic
zones, wake cut + Kutta condition for lift. Target: wings at M∞ 0.3–0.87,
workstation-scale (minutes for 1–3 M nodes).

## Document map (read the relevant one before coding)

- [docs/roadmap.md](docs/roadmap.md) — **active tracker**: phase order (Track P:
  P0–P12 solver, Track M: M0–M4 meshing, Track B: level-set wake — **IN PROGRESS**,
  B1–B5 + B7 closed / B6 open — and the designed-not-started Track V viscous
  coupling), gate checklists, progress ledger. "What phase are we in" and "what
  gate is open" live here, nowhere else. Track B numerics live in a separate
  spec, [docs/design_track_b.md](docs/design_track_b.md) (it supersedes DN1).
- [docs/design.md](docs/design.md) — theory & numerics reference: equations (§2–§3),
  wake/Kutta (§4), BCs (§5), discretization (§6), Numba kernel rules (§7), solver
  strategy (§8), V0–V6 validation ladder (§10), risks/mitigations (§12).
- [docs/demo_report.md](docs/demo_report.md) — **evidence dossier** for completed
  phases (P0, P1-partial, P2, P3, P4, P5, P6, P7, P8 + its capability assessment,
  P9, P10-partial G10.2/G10.3, M0, M1, and Track B B3/B5/B6/B7): one self-checking
  demo per phase under
  `cases/demo/<phase>/` with committed figures + measured gate numbers.
  When a phase closes, add its demo + report section here.
- [docs/discussion_notes/](docs/discussion_notes/) — **discussion & reference
  material only, NEVER a coding spec.** Design notes for future tracks (DN1
  level-set wake, DN2/DN6 VII coupling, DN4/DN5 Newton) +
  [PLAN.md](docs/discussion_notes/PLAN.md), the cross-track integration view
  (Chinese). When writing code, plan against **roadmap.md gates + design.md
  numerics only**; a discussion-note idea becomes actionable only after it is
  merged into roadmap.md/design.md (as Track B/V were on 2026-07-10). If a note
  contradicts roadmap/design, roadmap/design win. Sync PLAN.md when a gate
  closes.
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
3. Full suite: `pytest tests/` (**276 passed + 17 skipped + 2 xfailed since B7
   2026-07-12**, measured 719.29 s @16 threads: +92 Track B tests (B1 dual-mesh,
   B2 multivalued, B3 lifting, B4 TE-Kutta, B5 far-field, B6 transonic + LS
   Newton, B7 ONERA M6 3D) over the 184+8+2 P10/G10.2 baseline; B6 added
   `test_b6_transonic.py` (9 fast + 2 gated) and `test_b6_newton.py` (2 fast +
   2 gated); B7 added `test_b7_onera_m6.py` (6 fast + 5 gated — the gated ones
   are the M6 3D solves, ~20 min each, hence 12 → 17 skipped by default);
   some skip when the gitignored wake-free meshes aren't generated locally —
   M3 medium (~40 s) and the M4 ONERA M6 family (~12 s);
   ~5 min — G8.3 measured 301.66 s; the always-on coarse transonic smoke is ~170 s
   of it, the G3.2 medium-mesh nested Picard solve ~45 s, the rule-7 sweep's M6
   coarse+medium cut_wake ingest ~15 s. The M6 .msh files are gitignored — the 13
   M1 tests skip until you run `cases/meshes/onera_m6/generate_onera_m6.py`
   (~30 s). The heavy transonic/Newton gates (P4 medium + G4.3 sweep, P5, gated
   P8 G8.1/G8.2 + FD pocket, the gated B6 M0.80 dual-mesh + LS-Newton runs, and
   the gated B7 M6 3D dual-mesh solves) only run under
   `PYFP3D_TRANSONIC_GATES=1`, and make up most of the 17 skipped.)
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
P3 = subsonic compressible, P4 = transonic, P5 = ONERA M6, P6 = surface-Cp
recovery, P7 = differentiable flux at frozen selection, P8 = fully-coupled
Newton/performance, P9 = grid-convergence/accuracy-gap discrimination (evidence
phase), P10 = Newton generality & continuation efficiency, P11 = curved wall
elements, P12 = backlog; Track-P renumbers 2026-07-08 and 2026-07-11 (two
same-day insertions) — docs before those dates use the then-current IDs, so
pre-2026-07-11 "P9 curved walls"/"P10 backlog" read as P11/P12). Track V gates are
`GV<phase>.<n>` (phases V1–V4, distinct from validation-case IDs V0–V6);
design.md §11 mirrors the same order.
