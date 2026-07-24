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

- `v1_ibl3_standalone/` — **Track V / V1** GV1.1 standalone IBL3 verification (prescribed u_e, no FP coupling:
  laminar/turbulent flat plates + Falkner–Skan decelerating branch + refinement family, FE vs the closure's own
  2-D ODE marches; pre-registered, **9 PASS / 2 FAIL** — (a) ×2 = closure-family fixed point H*≈2.7083 ≠ Blasius,
  (e) first-run FAIL = under-damped streamwise 2h grid mode at outflow ⇒ fixed by the D-HB streamwise-tensor
  stabilization (ε_s=0.02, order ≈1.0 restored); SUPG/upwind remains the V3+ upgrade route;
  VERDICT + PRE_REGISTRATION + CSV/PNG in the dir, design record `docs/design_track_v.md` §9) —
  `python cases/analysis/v1_ibl3_standalone/run.py` — ~7 min (exit 1 = honest FAIL present)
- `v2_transpiration_channel/` — **Track V / V2** GV2.1 transpiration channel verification (δ*→ṁ = ∇_Γ·(ρ_e u_e δ*)
  wall-RHS channels in all five FP drivers, `None` ⇒ legacy path bit-identical; pre-registered, **23 PASS / 0 FAIL /
  16 RECORDED** — (a) MMS cylinder-blowing convergence strict-decreasing order 1.65/1.64 ≥ 1.0, (b) five-driver ṁ=0
  bit-identity, (c) FD Jacobian 6.6e-09–7.2e-08 < 1e-5; VERDICT + PRE_REGISTRATION + CSV/PNG in the dir) —
  `python cases/analysis/v2_transpiration_channel/run.py` — heavy (exit 1 = honest FAIL present)
- `v3_loose_coupling/` — **Track V / V3** GV3.1/GV3.2 loose viscous–inviscid coupling on NACA0012 2.5-D strip
  (M0.5/α2°/Re3e6 vs committed XFOIL reference `cases/reference_data/naca0012_viscous_xfoil/`; pre-registered,
  **2 PASS / 4 FAIL / 23 RECORDED** — PASS: Δcl ratio 0.542 ∈ [0.5, 2.0], GV3.2 loop convergence 5 outer iters ω=1.0
  (transonic M0.72 record point Newton 4 iters, no retuning); FAIL (binding band): cf first-post-transition station
  +44% (XFOIL e^N ramp vs our instantaneous switch, all other stations ≤15%), δ* H-family offset +13–27% on the lower
  side = closure-family difference; VERDICT + PRE_REGISTRATION + CSV/PNG in the dir, design record
  `docs/design_track_v.md` §10) — `python cases/analysis/v3_loose_coupling/run.py` — heavy (exit 1 = honest FAIL
  present)
- `v3_fuselage_smoke/` — **Track V / V3** GV3.3 fuselage body-of-revolution smoke (Track V's only fuselage-alone
  item: full-2π BoR, M0.3/α0/Re3e6 per body length, non-lifting Picard + loose coupling; pre-registered, **0 PASS /
  2 FAIL / 7 RECORDED** — closed-body scheme stabilized through three debug rounds: tail-band Dirichlet pin +
  transpiration masking on the pinned band + FP non-convergence guard; 10/10 outer iterations no numerical event,
  mid-body x/L ∈ [0.34, 0.82] axisymmetry excellent (σ/μ(δ*) 0.018–0.068, crossflow ~1e-6); FAILs localized:
  (a) σ/μ worst 0.5533 at x/L = 0.940 (post-trip ring + tail cone), (b) crossflow 0.2631/0.2295 tail-cone; loop
  NOT converged — tail-cone ṁ_max ×5.7 over k = 5→10 = measured loose-coupling stern instability, V4 decision
  input; VERDICT + PRE_REGISTRATION + CSV/PNG in the dir, mesh generator `cases/meshes/fuselage_bor/`, design
  record `docs/design_track_v.md` §10) — `python cases/analysis/v3_fuselage_smoke/run.py --levels coarse` —
  ~17 min (exit 1 = honest FAIL present)
- `v5_1b_scaled_newton/` — **Track V / V5** GV5.1b scaled + damped augmented Newton (the GV5.1
  follow-up designed on the IBL-floor diagnosis; pre-registered 8b7793f before the first execution;
  the GV5.1 amended protocol verbatim — HEAD-regen loose-converged seeds, wiring guard PASS both
  legs), **1 PASS / 1 FAIL / 7 RECORDED** — the machinery is delivered and exact (solver-internal
  row/column equilibration + Levenberg damping + floor-reached stop, flags default OFF = legacy path
  bit-identical; band (a) suite + tight fleet 28 green); the FAIL is the medium live-seed e2
  identity 1.96e-10 vs a ≤ 1e-10 threshold chosen at implementation time (NOT pre-registered) =
  SuperLU pivot-order roundoff through cond(J) ~ 1e10, user adjudication requested (VERDICT §3);
  band (b): the amended seeds sit INSIDE the 10× floor band from iter 0 (F_BL = 1.00× the floor) —
  no above-band contraction segment exists by construction → pre-registered fallback: medium
  termination floor_reached at iter 5 (replacing GV5.1's 10-step λ-collapse crawl) at the same
  merit (9.074e-11 ≈ 9.025e-11), coarse ends below GV5.1 and still descending
  (2.044e-10 < 2.068e-10), the k=1 standalone descends markedly deeper (F_BL 3.268e-6, −31 % vs
  the k1seed; merit 2.3× below); band (c) coarse 10 vs 8 NOT met, medium 5 vs 10 met (degenerate
  band-entry iter 0); μ rejection-retries 0 on all three runs — scaling the active ingredient; the
  window question REFRAMED to an above-band-seed protocol (candidate GV5.1c); VERDICT +
  PRE_REGISTRATION + CSVs in the dir, design record `docs/design_track_v.md` §14) —
  `python cases/analysis/v5_1b_scaled_newton/run.py` — coarse leg ~5 min / medium ~27 min /
  k=1 ~2 min under external load (2026-07-24 measurement; the default full run exits 1 = the
  medium live-check FAIL above, see VERDICT §3; `--levels coarse` alone exits 0; a loose-regen
  wiring-guard failure raises RuntimeError = recipe error; `--levels`/`--no-k1` for partial
  re-runs)
- `v5_ibl_floor/` — **Track V / V5** IBL-floor diagnosis (the GV5.1 follow-up; RECORDED diagnostic study, no
  pass/fail bands; pre-registered 53bf904 before the first execution): dense SVD of J_BL,BL at the coarse +
  medium loose-converged states and the coarse k=1 fixture, **14 RECORDED** — the near-null cluster PERSISTS
  at the converged states (S1 500/1236 <1e-6·σmax cond 1.3e11; S2 1082/2460 cond 4.0e13; the s1/s3 spectra
  overlap curve-for-curve), carried by the turbulent (A, Ψ) variables mid-chord → TE; the raw cond 4e10–4e13
  is MOSTLY a scaling artifact (row+col equilibration → 2e4/7e5/1e7, sub-1e-6 count 501/500/1082 → 0/0/2,
  no exact null directions) with a genuine scaled (A, Ψ) stiffness 1e5–1e7 remaining = the GV5.1b/GV5.4
  target; the F_BL floor lives in the TE band (B, δ) equations essentially entirely INSIDE J's range
  (left-null alignment ≤ 7.7e-3); closure-floor active set EMPTY (DELTA_MIN sensitivity identically zero —
  the Q6(b) substitute); eps_diff ×4 moves the floor ≤ 6 % (not an artificial-viscosity truncation); the
  pseudo-time controller bottoms out with the residual frozen at 3.154e-6 from iter 0 = a formulation floor
  globalization alone cannot pass; findings `results/findings.md` (+ summary.csv + CSV/PNG, design record
  `docs/design_track_v.md` §13) —
  `python cases/analysis/v5_ibl_floor/run.py` — ~50 min under external load (2026-07-24 measurement;
  minutes-scale unloaded; exit 0 always — RECORDED study with no bands; a loose-regen wiring-guard failure
  raises RuntimeError = recipe error; `--states`/`--phases` for partial re-runs)
- `v5_m6_bridge/` — **Track V / V5** GV5.0 M6 subsonic loose-coupling bridge (RECORDED entry check: ONERA M6
  coarse+medium conforming, M0.5/α3.06, forced x_tr/c 0.05, Re_MAC 11.72e6; V3 loose driver + new
  `viscous/coupling.py::build_wing_case` 3-D wing IBL case — LE-band laminar pin per local x/c, both TE natural
  outflow, root symmetry natural, tip band z > 0.95·b_semi pinned + ṁ-masked; pre-registered,
  **16 RECORDED / 0 FAIL** — bridge answer: the loose loop does NOT converge ≤10 at either level; coarse =
  root-upper-TE separation patch (H 4–5.5) feedback runaway ṁ_max ×12.4 (GV3.3-stern class), medium = patch
  resolved away but bounded δ* limit cycle 2–12 %/k; ΔCL DOWN both estimators (coarse −5.2 %/−4.8 %, medium
  −2.4 %/−2.1 % input-limited); crossflow first live 3-D exercise max|B|/|A| ≤ 0.072; δ*(z) CSVs feed GV5.3's
  bands; VERDICT + PRE_REGISTRATION + CSV/PNG in the dir) —
  `python cases/analysis/v5_m6_bridge/run.py --levels coarse` — coarse ~25 min / medium ~6.5 h (exit 1 = honest
  FAIL present)
- `v5_tight_coupling/` — **Track V / V5** GV5.1 tight-coupling exactness + convergence (augmented (φ, Γ, BL)
  Newton on the NACA0012 2.5-D strip, coarse + medium; pre-registered incl. the Addendum-2 amended seed
  (loose-converged state, user-adjudicated), **9 PASS / 1 FAIL / 36 RECORDED** — (a) FD exactness PASS both
  levels (worst sweet-spot coarse 2.246e-8 seed / 2.244e-8 endpoint, medium 5.074e-9 seed+endpoint; masked
  0/1236 + 0/2460 rows; veps omission ≤ 3.0e-8 scaled, decision 5); (b) quadratic tail HONEST FAIL (medium
  binding): the polish runs 10 iterations un-converged, F_BL pinned from iter 0 at the loose-final IBL floor
  (medium 1.708e-6 / coarse 3.11e-6) = the intrinsic floor of the steady IBL residual on the
  cond(J_BL,BL) ~ 4e10 near-null manifold, NOT a coupling defect; (c) N_aug ≤ 2 not met standalone nor as
  polish (N_total 14/13 vs loose 4/5); finding: the committed GV3.1 medium fixed point is NOT reproducible
  (IBL-floor trajectory scatter — three code/env → three fixed points; diagnosis
  `results/gv5_1_medium_seed_diagnosis.md`; HEAD-regen seed user-accepted); VERDICT + PRE_REGISTRATION +
  CSVs in the dir, design record `docs/design_track_v.md` §12) —
  `python cases/analysis/v5_tight_coupling/run.py --levels coarse` — coarse ~1 min / medium ~3 min
  (exit 1 = honest FAIL present)

*(Rows for a2/b9/p14 added in A3 2026-07-18: they existed on disk but the
table still listed only a1. Note two rows are NOT Track A — `b9_*` and
`p14_*` are analysis studies belonging to Track B / Track P phases that live
here because they are studies, not capability demos. `v1_ibl3_standalone/`
is Track V's GV1.1 gate case, likewise a study; `v2_transpiration_channel/`,
`v3_loose_coupling/`, `v3_fuselage_smoke/`, `v5_m6_bridge/`,
`v5_1b_scaled_newton/`, `v5_ibl_floor/` and `v5_tight_coupling/` are
Track V's GV2.1 / GV3.x / GV5.x gate cases (the `v5_ibl_floor/` one a
RECORDED diagnosis, no bands), same status.)*
