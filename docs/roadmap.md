# pyFP3D Development Roadmap
## Active implementation roadmap — track index

**This file (plus the per-track files below) is authoritative for phase status
and gates.** Since 2026-07-15 the phase entries, gate checklists and progress
ledgers live in ONE FILE PER TRACK under [roadmap/](roadmap/) — split **verbatim**
from this file (nothing was reworded in the split). Any reference elsewhere of
the form "roadmap.md Track X / phase entry / ledger" resolves through the per-track bullets:

- **P — solver** — [roadmap/track_p.md](roadmap/track_p.md) — P0–P9 ✓ (P1: only G1.6 open, strict xfail) · P10 ◐ (G10.2/G10.3 ✓, G10.1
  open) · **P11 ✓ CLOSED 2026-07-19** (user-directed, opened+closed same day; sphere leg): curved wall elements measured NEGATIVE (medium
  11.56%→11.33% = the G1.4 oracle ceiling; superparametric O(h) risk fired) and **G1.6 RE-ATTRIBUTED** —
  the order collapse was the fixed-bulk sweep floor (far-only refinement: 3.17×, order 1.89) and 11.6% ≈ intrinsic P1 capability at
  h=0.08 (structured control ~2nd order with flat facets); G1.6 xfail stays, route fork (Option C re-spec / P2 wall layer / accept) =
  user's call · P12 backlog · P13 ◐ (G13.1 ✓, G13.2 conforming ✓, G13.3 subsonic ✓ / transonic NEGATIVE-open) · **P14 ✓ CLOSED
  2026-07-17** (pressure-equality Kutta estimator, from A2; opened + closed same day; G14.1–G14.7 ✓; demo 28 PASS): S1+S2 both gone —
  M0.84 roughness 0.0970→0.0043/0.0365→0.0024, all-station TE gap 0.2206→0.0040/0.1585→0.0024, TE spike 0.1143→0.0533;
  ★ **the conforming path now matches the level-set path** (cross-model V14.6: cl_p/cl_KJ 0.2776/0.2823 vs 0.2772/0.2813 = 0.15%/0.34%;
  the probe path was 4.5%/4.3% below). G14.7 re-specced at close from the probe G8.2 locks to the level-set oracle (the +4.85% cl_KJ move
  is the finding — 69% of P9's 0.019 gap was Kutta-estimator bias)
- **M — meshing** — [roadmap/track_m.md](roadmap/track_m.md) — M0, M1(+M1b), M2, M3, M4, M5 ✓ — **M2 ✓ (solver leg closed by B9
  2026-07-17; ledger erratum 2026-07-19: A3 claimed this line was fixed while it still read ◐)** (wing-body mesh ✓ 2026-07-13, **body
  re-spec'd + regenerated 2026-07-16**: 5 root chords, wing centered, 2-diameter ellipsoid nose, graded skin)
- **B — level-set wake** — [roadmap/track_b.md](roadmap/track_b.md) — B1–B5, B7, B8, B11–B32 ✓ (**erratum 2026-07-19:
  this row had silently dropped B16/B17/B18 — see the track file for those entries**) · B6 ◐ (coarse gate ✓;
  medium quantitative closed by GB15.4, regressed under the first B20 re-baseline, **restored by B21**) · **B21 ✓ CLOSED 2026-07-19**
  (executes the Kimi-inspection N1 finding): `freeze_side_state` captured the frozen selection on the UNPATCHED side field —
  the one B20 consumer the patch missed; aligning it restores the committed M6-medium M0.84 ramp (γ 0.088343, res 9e-14, 515 s) ⇒
  **GB20.7's "real capability loss" verdict is OVERTURNED** (the loss was B20's own patch gap, not the fix's intrinsic cost) · **B17 ✓
  CLOSED 2026-07-18** (resolves GB16.4: the far-field aux pin must carry jump=γ, not 0 — a BC-modelling error, not a non-convergence;
  `farfield_aux="pin_gamma"` new default both solvers; medium 0.2117 Picard / 0.2114 Newton post-B20) · **B16 ✓ CLOSED 2026-07-18**
  (wing-body LS-Newton churn = near-singular far-field aux block, cond1 9.1e18 → 8.7e6 pinned) · **B25 ✓ CLOSED 2026-07-19 (the cure)**
  (`inboard_clip` moves the wake sheet's inboard boundary to the fuselage surface / symmetry plane = conforming fragment topology —
  medium α=3.06 junction pocket corrM **14.66→0.63**, n_sup 88→0, cl_p +0.38% within [A, oracle];
  default None bit-identical) · **B24 ✓ CLOSED 2026-07-19 (negative)** (the pocket FOLLOWS the free edge — hypothesis re-confirmed —
  but the waterline-extension variants trade the singularity for equal-or-worse forms; the (b)-1 route is CLOSED) · **B23 ✓ CLOSED
  2026-07-19** (junction discriminator: the pocket is lift/wake-coupled — α=0 clean at both levels —
  attributed to the wake inboard FREE-EDGE singularity, NOT G1.6 faceting; P11 close-out input) · **B27 ✓ CLOSED 2026-07-20** (B18 demo
  refresh: the pocket-healed LS reaches the SAME ceiling site as conforming — LS+clip coarse M0.84 reached / medium M0.7625 vs conforming
  M0.84/M0.79; conforming + LS A/C legs bit-reproduce the committed B18/B26 anchors (336/336,
  `b27_b18_demo_refresh/g27_consistency.csv`); cross-model upgraded M0.5 (2.6%) + M0.65 medium (2.4% PASS) + M0.75 (2.5%)) · **B29 ✓
  CLOSED 2026-07-20** (flat-fragment = the wing-body LS production config, user-adjudicated B28 §6: B18 C side = clip+flat sheet,
  M0.5 anchors 0.2115/0.2184; medium ceiling 0.7625→**0.775**; cross-model **0.5/1.1/1.1 %** M0.5/0.65/0.75; GB18.5 live flat
  cl_fus 0.0382 vs conf 0.0423 — the B26 tilted ×2 out-band reading retired; demo 8/8) · **B28 ✓ CLOSED 2026-07-20**
  (cl_fus decoupling + GB9.4 RE-SPEC — the "fuselage spurious lift" label retired: cl_fus = physical carryover + wake-sheet POSITION
  sensitivity, NOT the G1.6 error; gate re-spec'd to out-band cross-model ≤15%, medium 7.0% PASS, b9 demo 8/8) · **B30 ✓ CLOSED
  2026-07-21** ((b)-class ceiling attribution: the conforming M0.80+ stall and the LS+clip-medium 0.775 death are the SAME mechanism —
  the wing-tip P13 free-edge singularity + high-M Newton, NOT a wake-model pocket; named the C-class tip cure) · **B31 ✓ CLOSED
  2026-07-22** (C-class wing-tip cure: production pressure+taper CURES the conforming 0.83 dying level via the Gamma-pin row blend
  in `newton.py` (frozen weld-sign, FD-verified); LS-side C-class CLOSED negative — C1 inboard backflow / C3 coarse divergence,
  `outboard_fringe` retained default-inert) · **B32 ✓ CLOSED 2026-07-22** (② weld-sign per-step refresh rolled back — ill-posed
  switching system, B31 frozen semantics restored bit-identical; ① conforming tip_taper adopted — wing-body medium ceiling M0.79 →
  **M0.84 reached** (cl_p 0.2738, 0 clamps), cl_p cost ≈ −1.3%, demo 8/8) · **B26 ✓
  CLOSED 2026-07-20 (B26-A)** (the junction pocket WAS the LS wing-body transonic ceiling limiter — healed by B25's `inboard_clip`:
  medium 0.50→0.7625, coarse 0.82→0.84 reached; death class flips (a)→(b), wing-tip P13; the A-side anchor divergence = the B21/B22
  freeze-capture effect) · **B18 ✓ CLOSED 2026-07-18** (wing-body transonic: conforming reaches M0.84/M0.79; **erratum 2026-07-20:
  the "LS junction-limited — G1.6 class" reading is RETIRED — the pocket is the B23 free-edge singularity, healed;
  see B26/B27**) · **B14 ✓ CLOSED 2026-07-17** (`precond="schur"` Schur+AMG structural preconditioner:
  the A1 precond bottleneck is gone, M6 medium M0.84 43.6%→2.6% (same-session A/B; A1's earlier independent read was 42.6%), ramp
  1.43×/subsonic 2.08×; small-scale slower, fine route still unbuilt) · B10 shelved · **B20 ✓ CLOSED 2026-07-18**:
  mixed-side plain elements read their density from the MAIN field (measured behind a temporary `plain_density` knob, REMOVED on
  adoption); the assembly now agrees with the diagnostic (`element_mach2` has defaulted to main since 2026-07-14). 2.5-D **+0.0000 %**;
  M6 coarse ramp side m0.7875-not-converged → main **M0.84 converged**; ★★ **GB20.5 the B19 hypothesis SPLITS** —
  B18 wing-body @M0.5 side res 6.8e-5/82 clamped/Mmax 3.920 (clamped, non-converged) vs main **res 1.1e-13/6 clamped**/Mmax 5.220
  (genuine) ⇒ the CONVERGENCE churn was largely the contamination, but the junction POCKET is REAL and the literal hypothesis is
  **REFUTED** (unclamping reveals a genuine M≈5.2 = G1.6 faceted geometry). ★★ **ADOPTED PERMANENTLY, knob REMOVED** (user-arbitrated):
  the side reading is an internal inconsistency (one equation, two velocity fields, in an uncut element), so it is not left behind a
  switch; 3-D committed LS numbers re-based, 2.5-D untouched · **B19 ✓ CLOSED 2026-07-18**: the LS-Newton Jacobian is now EXACT in 3-D —
  TWO defects (row/column DOF maps, then the gradient FACTORS: the residual is `ρ̃(read-field grad)·V·(scatter-field grad·B_a)` and both
  factors had used the side gradient); targeted FD probe 1.146e-01 → 1.334e-08 with the ε discriminator FLIPPING from independent to
  ~1/ε; R bit-identical (0.000e+00) ⇒ no converged result moves; ★ RECORDED NEGATIVE: **no convergence gain** (γ/M_max identical, same
  steps, +3.6 % wall — the plateau is the B15 churn limit cycle, which an exact derivative cannot fix).
  Leg B measured the residual's own asymmetry: the side field manufactures a **spurious supersonic state** (q² 3.22 vs 1.34) on 252
  elements ⇒ routed to a new phase, NOT adopted · **B9 ✓ CLOSED 2026-07-17 (RE-SPEC'D, user-approved): wing-body cross-model validation —
  LS (Picard) + conforming (NEW capability, Newton) agree to 0.4%/0.6% at medium M0.5; GB9.4 fuselage-lift XFAIL ⇒ G1.6 error
  (**corrected by B28 2026-07-20**: wake-sheet POSITION sensitivity, not an error — gate re-spec'd to out-band cross-model ≤15%,
  medium gap 7.0% PASS, demo 8/8);
  LS Newton diverges = neumann blockage not the solver**
- **V — viscous coupling** — [roadmap/track_v.md](roadmap/track_v.md) — designed 2026-07-09/10 · **V1 ◐ OPENED 2026-07-22**
  (gates re-spec'd at opening against the B32/A4 state; same day re-phased: V1 IBL3 core / V2 transpiration channel /
  V3 loose coupling incl. the fuselage body-of-revolution smoke; V4 optional; V5 tight coupling with the GV5.0 M6
  subsonic bridge and GV5.3 anchored on committed Cp; V6 wake sheet; wing-body VII deferred)
- **A — verification & analysis** — [roadmap/track_a.md](roadmap/track_a.md) — created 2026-07-15 · **A1 ✓ CLOSED 2026-07-16**
  (GA1.1–GA1.5: 4-driver timing instrumentation + conforming-vs-level-set × Picard-vs-Newton cost benchmark;
  3-D Newton is precond-bound, the 2.5-D seed headline does not transfer) · **A2 ✓ CLOSED 2026-07-17** (TE/Kutta fidelity attribution,
  GA2.1–GA2.5: **S1** the conforming Γ(z) jitter is a measurement artifact of the per-station probe-difference Kutta target, not flow
  content (fixed-Γ discriminator D=7.33/25.70 coarse/medium; closure |F|/|Γ| ≤ 0.6%); **S2** the TE Cp jump = potential-jump Kutta form
  error (34×/133× vs level-set) + a P1 recovery artifact (both paths); fix routed to **P14**, which built it and confirmed the
  attribution — see the Track-P row) · **A3 ✓ CLOSED 2026-07-18** (GA3.1–GA3.6: response to the 2026-07-17 independent Kimi inspection —
  ★★ **C1 VERIFIED**: the LS Newton Jacobian is not dR/dφ on mixed-side plain elements (M6 coarse M0.70:
  targeted probe 1.146e-01 vs control 6.33e-10, and **eps-independent** — 1.532e-01 at eps 1e-6/1e-7/1e-8, max/min 1.00 ⇒ a missing term,
  not FD noise) ⇒ the LS Newton is a **quasi**-Newton in 3-D; **R is untouched so every converged state and gate number stands**;
  RECORDED not fixed, the kernel fix is its own phase. Plus C2/C3 backported to the conforming Newton, reader C4/C5 (silently dropped
  surface groups ⇒ Γ(root) silently pinned to 0), C6/C7/P1/T1/T2/F0, 17/17 docs findings dispositioned, and the close-out ritual extended
  to five surfaces + a backport check. No `pyfp3d/` numerics change)

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
`GV<phase>.<n>`; Track A: `GA<phase>.<n>`). Phases were renumbered three times on Track P (2026-07-08, then
twice on 2026-07-11) and three times on Track B (twice on 2026-07-12, once on
2026-07-13); documents dated
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
