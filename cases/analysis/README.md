# Analysis & verification cases (Track A)

Cross-cutting *analysis* studies that measure the solver/mesh machinery rather
than add to it — profiling, method A/B comparisons, cost accounting. The
counterpart to `cases/demo/` (per-phase capability demos), but for Track A work.

Each `run_*.py` regenerates every figure/CSV in its `results/` folder and
self-checks the study's gate numbers (exit code 0 = all checks pass). Heavy 3-D
parts run only under `PYFP3D_TRANSONIC_GATES=1` and cache to gitignored `.npz`;
the committed PNG/CSV are the evidence. Roadmap gates:
[docs/roadmap/track_a.md](../../docs/roadmap/track_a.md); conclusions:
[docs/demo_report/track_a.md](../../docs/demo_report/track_a.md).

- `a1_solver_bottleneck/` — A1 conforming-vs-level-set × Picard-vs-Newton cost benchmark (GA1.1–GA1.5) —
  `python cases/analysis/a1_solver_bottleneck/run_a1.py` — ~5 min (+ gated 3-D `run_a1_m6.py` ~1 h)
- `a2_te_kutta_fidelity/` — A2 TE/Kutta fidelity attribution (GA2.1–GA2.5; S1 jitter + S2 TE Cp gap) —
  `python cases/analysis/a2_te_kutta_fidelity/run_a2.py` — ~90 s (+ gated `run_a2_interventions.py`)
- `c1_ls_jacobian_fd/` — **A3 / GA3.5 + B19 + B20 + B21** — is the LS Newton Jacobian exact on a 3D mesh?
  (kimi code-review C1; `run_check.py` measures it before/after the B19 fix, `run_legb_probe.py` measures B19 Leg B's residual asymmetry;
  `run_legb_apply.py` / `run_legb_beforeafter.py` / `run_legb_b18.py` are B20's before/after — HISTORICAL post-adoption,
  the knob is gone; `run_gb207_recipe.py` = the GB20.7 freeze_tol sweep, verdict overturned;
  `run_n1_freeze_capture.py` = **B21's discriminator** — post-N1-fix the committed M6-medium ramp reaches M0.84 again,
  `results/n1_freeze_fix_sweep.csv`) — `python cases/analysis/c1_ls_jacobian_fd/run_check.py` ·
  `.../run_legb_probe.py` · `.../run_n1_freeze_capture.py` — ~4 min each; the two ramp sweeps ~18 min / ~90 min
- `p14_te_pressure_diag/` — P14 pressure-Kutta diagnostics (20 checks; feeds the P14 demo) —
  `python cases/analysis/p14_te_pressure_diag/run_diag.py` — heavy
- `b9_fuselage_guardrail/` — B9 fuselage-Cp guardrail (GB9.6, RECORDED — no pass/fail) —
  `python cases/analysis/b9_fuselage_guardrail/run_guardrail.py` — heavy
- `b23_junction_discriminator/` — B23 wing-body junction discriminator (P11 close-out mandate:
  GB9.4/GB20.5 get their own discriminator — D1 alpha sweep / D2 junction-only refinement / D3 fairing A-B
  (recorded low-cost infeasible) / D4 canonical crease / D5 pocket localization / W2 fuselage-lift decomposition;
  verdict + LS-workflow route in `VERDICT.md`, pre-registration in `PRE_REGISTRATION.md`) —
  `python cases/analysis/b23_junction_discriminator/run_d1.py` · `.../run_d2.py` · `.../run_d4.py` ·
  `.../run_d5_localize.py` · `.../run_w2_decomp.py` — heavy (LS solves ~15 min each; D4 h01 leg ~10 min)
- `b24_wake_inboard_end/` — B24 wake inboard free-end fix A/B (executes the B23 verdict:
  TE-polyline inboard end extended along the fuselage waterline so the near field has no free edge;
  VERDICT: free-edge hypothesis re-confirmed but both extension variants replace the pocket with an equal-or-worse singularity —
  route closed, falls back to B23 §(b)-2 P13-inboard; side product = the R8 levelset panel-selection fix for multi-segment TEs) —
  `python cases/analysis/b24_wake_inboard_end/run_e1.py` — heavy (reuses the b23 meshes + control data)
- `b25_inboard_fragment_clip/` — B25 wake inboard fragment clip (executes B23 §(b)-2 as the conforming-topology fix:
  the LS sheet's inboard clip moves from the junction station (q>=0) to the wake-surface/fuselage intersection trace —
  waterline -> tail -> symmetry plane, so the sheet's boundary never sits in the fluid interior;
  weld/taper v1 REJECTED as unphysical + B8 dead-oracle; pre-registered in `PRE_REGISTRATION.md`) —
  `python cases/analysis/b25_inboard_fragment_clip/run_f1.py` — medium (reuses the b23 meshes + control data)
- `b26_ls_transonic_ceiling/` — B26 pocket-healed LS transonic ceiling re-measure (same-code A/C with the only
  variable = `inboard_clip`, B18 recipe verbatim: **medium ceiling 0.50 → 0.7625, coarse 0.82 → 0.84 reached**;
  C-side death class flips (a)→(b), dying peak at the wing TIP = P13 class ⇒ the LS ceiling is co-located with conforming;
  T1: the A-side anchor divergence = the B21/B22 freeze-capture effect; verdict + pre-registration in the dir) —
  `python cases/analysis/b26_ls_transonic_ceiling/run_g1.py` — heavy (4 ramps, ~69 min solve)
- `b27_b18_demo_refresh/` — B27 B18 demo refresh (LS legs revived) + consistency proof: conforming legs
  bit-reproduce the B18 anchors, LS A/C legs bit-reproduce the B26 anchors (**336/336 rows bit-identical** in
  `results/g27_consistency.csv`); cross-model upgraded to M0.65 (2.4% PASS) + M0.75 (2.5%); refreshed demo evidence lands in
  `cases/demo/b18_wingbody_transonic/results/` (checks.csv **8/8 PASS**) — (regeneration =
  `cases/demo/b18_wingbody_transonic/run_demo.py`, gated) — heavy (~1 h 39 min full re-solve)

*(Rows for a2/b9/p14 added in A3 2026-07-18: they existed on disk but the
table still listed only a1. Note two rows are NOT Track A — `b9_*` and
`p14_*` are analysis studies belonging to Track B / Track P phases that live
here because they are studies, not capability demos.)*
