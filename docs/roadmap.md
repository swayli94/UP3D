# pyFP3D Development Roadmap
## Active implementation roadmap — track index

**This file (plus the per-track files below) is authoritative for phase status
and gates.** Since 2026-07-15 the phase entries, gate checklists and progress
ledgers live in ONE FILE PER TRACK under [roadmap/](roadmap/) — split **verbatim**
from this file (nothing was reworded in the split). Any reference elsewhere of
the form "roadmap.md Track X / phase entry / ledger" resolves through the table:

| Track | File | Status (one line — details + gates + ledger in the track file) |
|-------|------|----------------------------------------------------------------|
| **P — solver** | [roadmap/track_p.md](roadmap/track_p.md) | P0–P9 ✓ (P1: only G1.6 open, strict xfail) · P10 ◐ (G10.2/G10.3 ✓, G10.1 open) · P11 conditional, not opened (only the G1.6 rationale survives) · P12 backlog · P13 ◐ (G13.1 ✓, G13.2 conforming ✓, G13.3 subsonic ✓ / transonic NEGATIVE-open) · **P14 ◐ 2026-07-17** (pressure-equality Kutta estimator, from A2: G14.1–G14.6 ✓ — S1+S2 both gone, M0.84 roughness 0.0970→0.0043/0.0365→0.0024, all-station TE gap 0.2206→0.0040/0.1585→0.0024; ★ **cross-model V14.6**: conforming-pressure lands ON the level-set answer (0.2776/0.2823 vs 0.2772/0.2813 = 0.17%/0.36%; the probe path was 4.5%/4.3% below) ⇒ the conforming-vs-LS lift disagreement *was* the Kutta form; ★ **G14.7 XFAIL, user arbitration open** — the swap moves cl_KJ +4.85% off the probe-path G8.2 locks and closes 69% of P9's 0.019 gap, pre-registered) |
| **M — meshing** | [roadmap/track_m.md](roadmap/track_m.md) | M0, M1(+M1b), M3, M4, M5 ✓ · M2 ◐ (wing-body mesh ✓ 2026-07-13, **body re-spec'd + regenerated 2026-07-16**: 5 root chords, wing centered, 2-diameter ellipsoid nose, graded skin; solver leg = B9) |
| **B — level-set wake** | [roadmap/track_b.md](roadmap/track_b.md) | B1–B5, B7, B8, B11–B13, B15 ✓ · B6 ◐ (coarse gate ✓; medium quantitative closed by GB15.4) · B14 designed-not-scheduled · B10 shelved · **B9 (wing-body LS solve, M∞0.5) = NEXT** |
| **V — viscous coupling** | [roadmap/track_v.md](roadmap/track_v.md) | designed 2026-07-09/10, zero implementation |
| **A — verification & analysis** | [roadmap/track_a.md](roadmap/track_a.md) | created 2026-07-15 · **A1 ✓ CLOSED 2026-07-16** (GA1.1–GA1.5: 4-driver timing instrumentation + conforming-vs-level-set × Picard-vs-Newton cost benchmark; 3-D Newton is precond-bound, the 2.5-D seed headline does not transfer) · **A2 ✓ CLOSED 2026-07-17** (TE/Kutta fidelity attribution, GA2.1–GA2.5: **S1** the conforming Γ(z) jitter is a measurement artifact of the per-station probe-difference Kutta target, not flow content (fixed-Γ discriminator D=7.33/25.70 coarse/medium; closure |F|/|Γ| ≤ 0.6%); **S2** the TE Cp jump = potential-jump Kutta form error (34×/133× vs level-set) + a P1 recovery artifact (both paths); fix routed to **P14**, which built it and confirmed the attribution — see the Track-P row) |

[design.md](design.md) remains the design reference for equations, numerics and
architecture; Track B numerics live in [design_track_b.md](design_track_b.md)
(supersedes DN1). Evidence for closed phases: [demo_report.md](demo_report.md)
(index) + [demo_report/](demo_report/) (per track). Human-readable status
snapshot + document map: [overview.md](overview.md). Analysis/review reports
(non-normative): [analysis/](analysis/).
(`docs/discussion_notes/` — the DN1–DN6 design-note sources and the PLAN.md
integration view — was deleted 2026-07-14, commit 0e4895a; historical copies via
`git show 8aa4aee:docs/discussion_notes/<file>`.)

**Gate-ID / renumbering conventions.** Gate IDs are `G<phase>.<n>` (Track V:
`GV<phase>.<n>`; Track A: `GA<phase>.<n>`). Phases were renumbered twice on Track P (2026-07-08 and
2026-07-11) and twice on Track B (2026-07-12 and 2026-07-13); documents dated
before those days use the then-current IDs. The authoritative mapping notes are
kept verbatim inside the affected phase entries in the track files.

## 0. Working rules

- A phase is complete only when its medium-mesh gate passes and the full coarse regression suite stays green.
- After any kernel change, run `pytest tests/test_v0_freestream.py` before any broader validation.
- Use SciPy/PyAMG for linear algebra, keep Numba kernels SoA-only, and validate against full-potential references before Euler.
- Keep state in git plus this file; update the progress ledger when a gate closes.
- Every visual gate must have a headless path: generate PNG/CSV artifacts by script (no GUI-only checks).

## 0.1 Headless Linux feasibility and tooling policy

The current plan is feasible on Linux without a graphical desktop, with one constraint:
all visualization checks must be script-driven and artifact-based.

Required execution modes:
- Mesh generation/inspection: use `gmsh` CLI (`-2/-3`, `-format msh4`) and scriptable quality reports.
- Solver/testing: run via `pytest` and Python CLIs only.
- Visualization artifacts: prefer `matplotlib` for plots and `pyvista` off-screen rendering; use `pvpython`/`pvbatch` if ParaView is installed headless.
- CI compatibility: every gate that references a figure should accept generated PNG + numeric CSV summary as evidence.

Non-goals for v1.0:
- No dependency on interactive ParaView GUI for gate closure.

---


## Agent Workflow Bindings

- **agent-rules.md** (kept < 30 lines): points to design.md + this file, states
  the current phase, the hard rules from §0, and "do not edit files under
  `cases/reference_data/`".
- **Session protocol:** per work item — plan mode first (Claude reads the
      relevant design.md sections + this plan's phase entry, proposes an
  implementation plan, human review), then implement, then gates, then commit.
- **Skill sedimentation checkpoints:** end of P1 and end of P2, distill
  recurring numba pitfalls and the test-first procedure observed so far into
  `.claude/skills/fp3d-dev/SKILL.md`. Do not write the skill before P1 —
  capture real failure modes, not anticipated ones.
- **Subagents (from P4 onward):**
  - `validation-runner` (read-only + bash): runs the full gate suite,
    returns a residual/convergence summary table — keeps hundreds of lines of
    pytest output out of the main context.
  - `derivation-checker` (read-only): cross-checks kernel implementations
      against design.md §2–§6 formulas; invoked before closing P4 and P7.
- **Reference data discipline:** digitized FP reference Cp/shock data lives
  in `cases/reference_data/` with provenance notes (source figure, digitizing
  method); gates compare against these files, never against numbers embedded
  in test code.

---

## Progress ledger

Moved per-track (2026-07-15): each `roadmap/track_*.md` ends with its own
"Progress ledger" section (status legend: ✓ closed · ◐ partially closed /
in progress · ☐ open or not started · ⊘ shelved).
