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
- `a4_ue_error_band/` — A4 wall edge-velocity (u_e) error-band study = the Track-V IBL input-quality
  prerequisite (analytic ground truth: cylinder u_e=2Usinθ + sphere u_e=1.5Usinθ, coarse+medium ×
  linear/quadratic recovery). Medium smooth-wall band ≈ 2.5% peak-relative / 0.04·U∞ max-norm / O(h);
  LE/stagnation band worst-relative (the IBL du_e/ds zone); NACA0012 16° TE clears the sub-6° quadratic
  guard (corrects the audit's blanket TE claim). VERDICT + CSV/PNG in the dir —
  `python cases/analysis/a4_ue_error_band/run.py` — ~30 s (cheap M0 solves, no cache)
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
- `b28_cl_fus_flat_sheet/` — B28 cl_fus out-band flat-vs-tilted decoupling + GB9.4 re-spec (decisive leg: flat-fragment
  out-band cl_fus 0.0326 vs conforming oracle 0.0351 = 7.25% <= 15% TOL ⇒ **F1: sheet-POSITION sensitivity, not an error**;
  the "fuselage spurious lift" label retired; verdict + pre-registration in the dir) —
  `python cases/analysis/b28_cl_fus_flat_sheet/run_f2.py` — medium (~40 min)
- `b32_tip_taper_adoption/` — B32 ② weld-sign freeze fix + ① conforming taper production adoption (user-adjudicated from
  the B31 verdict) — **GB32.1 ✗ ROLLED BACK**: per-step weld-sign refresh turns the fixed system into a state-dependent
  switching system (ill-posed — healthy 0.82 seed diverged; B31 frozen semantics restored bit-identical; the 0.84
  fresh-seed hazard is handled by the F2 healthy-seed pattern) / **GB32.2 ✓ ADOPTED**: b18 CONF legs on tip_taper
  (vanish_smooth 0.05·b_semi), conforming medium 0.50–0.79 strict re-solved 0 clamps, **ceiling climb 0.79 → 0.84
  REACHED** (cl_p 0.2738, 0 clamps), demo checks 8/8 PASS, cl_p cost ≈ −1.3 % in the F3 band
  (`results/g2_adoption_cost.csv`), cross-model gaps improved to 0.3 % (M0.65) / 0.2 % (M0.75) / GB32.3 VERDICT in dir —
  `python cases/analysis/b32_tip_taper_adoption/run_g1.py` · `.../run_g2.py` — solves cached, CSV/PNG regenerate in minutes
- `b31_tip_termination/` — B31 C-class wing-tip cure (sheet-termination re-spec) + LS step-semantics companion
  (GB31.1 tip atlas PASS: all 8 dying-level clamps = cap_wall at the sheet tip edge / GB31.2a taper factorial
  ⇒ **taper is the cause** (0.83 die→cure pair, 0 clamps) / GB31.2b production pressure+taper port ⇒ **✓
  dying level 0.83 cured, 0.84 converged from a healthy seed** (chain-seed failure diagnosed as the fixable
  weld-sign freeze hazard, F2 probe; cl_p cost −3.00% flagged) / GB31.3 LS ladder **C1 ✗ (inboard backflow
  −19.5%) / C3 ✗ (coarse divergence) ⇒ LS-side C-class closed (negative)**, remaining route B10 roll-up /
  GB31.4 step-semantics closed; verdict + pre-registration in the dir, library changes default-off
  bit-identical) — `python cases/analysis/b31_tip_termination/run_g1.py` · `.../run_g2.py` · `.../run_g2b.py` ·
  `.../run_g2c.py` · `.../run_g3.py` — all solves cached, CSV/PNG regenerate in minutes
- `b30_transonic_ceiling/` — B30 (b)-class wing-body transonic ceiling attribution + López dissipation lever
  (GB30.1 anchor PASS 4/4 / GB30.2 dying-level clamp census: both paths 100% tip-localized O(10) ⇒ **B30-SAME** /
  GB30.3 `upwind_c=2.0` climb + `[1.8,1.6]` staging ⇒ **LS L1-◐ (one rung, elevated dissipation only) / CONF L1-✗
  (true Newton stall)** ⇒ dissipation NOT the constraint; **C-class tip cure named next candidate**, user's call;
  verdict + pre-registration in the dir, zero `pyfp3d/` change) —
  `python cases/analysis/b30_transonic_ceiling/run_g1.py` · `.../run_g2.py` · `.../run_g3.py` — all solves cached,
  census/CSV/PNG regenerate in minutes

*(Rows for a2/b9/p14 added in A3 2026-07-18: they existed on disk but the
table still listed only a1. Note two rows are NOT Track A — `b9_*` and
`p14_*` are analysis studies belonging to Track B / Track P phases that live
here because they are studies, not capability demos.)*
