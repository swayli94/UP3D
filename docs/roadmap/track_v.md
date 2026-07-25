# pyFP3D Roadmap — Track V (viscous–inviscid interaction V1–V6)

> Split from `docs/roadmap.md` on 2026-07-15. **No longer verbatim**: re-specified
> at Track V opening 2026-07-22, and re-phased the same day (the oversized V1
> split into V1 core / V2 channel / V3 coupling; old V2–V4 → V4–V6). Global
> working rules, gate-ID conventions and the track index live in
> [roadmap.md](../roadmap.md); the human-readable status snapshot is
> [overview.md](../overview.md).

## Track V — Viscous–inviscid interaction (designed 2026-07-09/10; **V1 ✓ CLOSED 2026-07-22 · GV1.1 9P/2F** · **V2 ✓ CLOSED 2026-07-22 · GV2.1 23P/0F** · **V3 ✓ CLOSED 2026-07-22 · GV3.1/3.2 2P/4F/23R · GV3.3 0P/2F/7R** · **V4 ⊘ SKIPPED 2026-07-22** · **V5 ✓ CLOSED 2026-07-25 · GV5.0 ✓ 16R/0F · GV5.1 ✓ 9P/1F/36R · IBL-floor diag ✓ 2026-07-24 14R · GV5.1b ✓ 2026-07-24 2P/0F/7R (1P/1F/7R as executed; (a)-medium cond-aware PASS adjudicated 2026-07-24) · GV5.1c ✓ 2026-07-24 2P/1F/7R (the above-band window read: NO slope-2 above the floor; mid-range stall) · GV5.1d ✓ 2026-07-24 2P/1F/7R (the near-band window read: NO basin adjacent to the floor either; the stall extends down to 24× floor) · GV5.5 ✓ 2026-07-24 2P/1F/9R (the TE-outflow row replacement does NOT break the floor — m2 5554×/245998× the floor, the pre-registered "worse" clause; flag stays default-OFF) · GV5.2 ✓ 2026-07-25 band-(b) FAIL + the loose-recipe transonic-limit anatomy (1/4 legs converged; every computed shock 0.06–0.10 c aft of the experimental bracket) · GV5.3 ✓ 2026-07-25 0P/1F/17R (band (b) honest FAIL — the loose viscous Cp does NOT move toward the committed M6 experiment: 1/5 unmasked stations + a pooled increase; band (a) RECORDED input-limited, Δcl_KJ −2.20 % DOWN under the A4 floor) · GV5.4 ✓ 2026-07-25 0P/1F/17R (augmented step 7.53× the inviscid step RECORDED, above the ≤ ~2× reference band; block preconditioner NOT-WORKING honest FAIL — block-Jacobi diverges, exact-BL Schur stalls 1/4; Schur-aware reduced-space preconditioning = the registered follow-up)**)

Deliverable: `pyfp3d/viscous/` — Drela IBL3 6-equation integral boundary layer
(δ, A, B, Ψ, C_τ1, C_τ2; surface Galerkin P1 FE on wall + wake sheet — **no
streamline integration**) coupled to the FP solver through a **transpiration
BC** (no mesh motion, RHS-only; δ* = 0 bit-identical to the inviscid path),
progressing loose → tight coupling.

> **★ 2026-07-22 re-spec (Track V opening, user-directed).** The phase entries
> and gates below were re-specified against the shipped code state BEFORE any
> implementation. Evidence base: the 2026-07-20 review
> [20260720-2015-wingbody-trackv-review.md](../inspection/20260720-2015-wingbody-trackv-review.md)
> §3 (hooks audit: ~70 % of the solver-side hooks exist), the 2026-07-22
> pre-Track-V audit
> [20260722-0335-b28-b32-audit-pre-trackv.md](../inspection/20260722-0335-b28-b32-audit-pre-trackv.md)
> §5, and **A4** (u_e input error band,
> `cases/analysis/a4_ue_error_band/VERDICT.md`). The pre-2026-07-22 GV1.1–GV1.3
> sketch is **superseded** by the gates below (Track V had zero implementation,
> so no historical gate result is affected). Binding design = THIS file. The
> historical notes DN2/DN6 (recover via `git show 8aa4aee:docs/discussion_notes/...`)
> carry two known traps: DN6 §10.2's streamline-integration module layout
> (`streamline.py` etc.) is obsolete — the binding route is IBL3 surface FE —
> and DN6 §8.3's "wall-term placeholder in `kernels/residual.py`" has **no code
> counterpart** (that kernel is pure volume assembly). The transpiration RHS
> structural template is `solve/wall_correction.py::assemble_wall_flux_correction_rhs`
> (P11-era infra, default-inert, annotated negative for its own G1.3 purpose but
> kept precisely as this template).
>
> **★ 2026-07-22 re-phase + reference pin (same day, user-directed).** The
> oversized V1 was split three ways — **V1** IBL3 core (standalone), **V2**
> transpiration channel through all drivers, **V3** loose coupling (2.5-D
> ladder + fuselage smoke) — and the original V2–V4 renumbered **V4–V6** (all
> gate IDs re-mapped accordingly; still zero implementation, so no gate result
> is affected). GV5.3 is anchored on the committed M6 experiment **Cp** data
> (no experimental CL value is committed). Binding reference on hand: Drela
> 2013 = AIAA 2013-2437, DOI 10.2514/6.2013-2437, local copy
> `docs/references/Drela_2013_IBL3_general_configurations.pdf` (symlink to the
> uploaded PDF; `docs/references/` is gitignored ⇒ local-only, re-fetch on a
> fresh clone). The BLWF58 method document (TsAGI; optional V4 only) is NOT on
> hand — non-blocking given V4's concrete skip criterion; fetch it before
> opening V4, or skip.

**Prerequisite state at opening (all measured, none pending):**

- **P6 ✓** smoothed wall tangential gradients
  (`post/surface.py::smooth_wall_tangential_gradients`) = the u_e / du_e/ds input path.
- **A4 ✓ (2026-07-22)** u_e inviscid input error band, analytic ground truth:
  medium smooth-wall ≈ **2.5 % peak-relative / 0.04·U∞ max-norm / 0.012·U∞ rms, O(h)**;
  **LE/stagnation band 4–7 % @ medium** (= the IBL seeding / du_e/ds zone, the least
  trustworthy input zone); linear-vs-quadratic recovery has no universal winner
  (~1 % region-dependent) → per-zone choice, LE band linear+smoothed; the sub-6°
  quadratic-recovery guard does **not** fire on NACA0012 (16° TE) — thin-TE
  (RAE-class) airfoils must re-check (GV5.2 clause).
- **P8/P14 ✓** conforming coupled (φ,Γ) Newton + pressure-Kutta (V5's augmentation
  base); **B14 ✓** `solve/schur_ls.py` (block-preconditioning structural prototype).
- **Solver hooks in place**: Picard `body_source_rhs` (`solve/picard.py`), LS Newton
  `b_base` slot (`solve/newton_ls.py`), wall-RHS assembly template
  (`solve/wall_correction.py`), wall triangle adjacency + area/normal infra
  (`post/surface.py`), colored parallel assembly (`mesh/coloring.py`).
- **Missing, to be built** (work-size order): the `viscous/` package itself (the
  Track-P-sized bulk); the conforming-Newton external-RHS channel (small, touches
  the core solver); the committed XFOIL viscous reference (external generation,
  `cases/reference_data/naca0012_incompressible/generate_panel_reference.py`
  precedent — XFOIL is NOT a repo dependency); the LS-path sheet-surface
  integration (V6 only, design adjudication GV6.0).
- **Reference data committed**: `cases/reference_data/rae2822_experiment/`
  (M0.725/α2.55 + M0.73/α3.19, Re 6.5e6), `naca0012_experiment/`,
  `onera_m6_experiment/` (Cp only — no CL values; see GV5.3).

> **Naming disambiguation.** Track V phase IDs V1–V6 are development *phases*;
> they are unrelated to the §10 validation-ladder case IDs V0–V6 in design.md
> (e.g. "V6 < 1%" remains the Γ-consistency metric — with phases now also
> running to V6, always read "V\<n\>" against context). Gates here are
> GV<phase>.<n>.

### V1 — IBL3 solver core (standalone verification) ✓ CLOSED 2026-07-22 · GV1.1 (9 PASS / 2 FAIL)

**Deliverable:**

- `viscous/surface_mesh.py` — compact wall-surface DOF numbering on the existing
  wall triangulation (`post/surface.py::wall_triangle_adjacency`); per-node local
  Cartesian basis (Drela 2013 §III.B — in-plane rotation invariance absorbs the
  TE kink, no special TE equations). **Data-layout design points recorded at
  build time**: (1) reserve the wake-sheet unknowns (V6 continuation — same 6
  equations, wake closures); (2) a master-map hook so an IBL surface mesh built
  on the *uncut* wall can be fed from cut-mesh (LS) solutions; (3) single-group
  (`wall`) scope — the `wall`+`fuselage` seam is wing-body, out of V1 scope.
- `viscous/closures.py` — Drela 2013 wall + wake closure fits; laminar +
  turbulent + **forced transition** (free-transition e^N is a recorded follow-up,
  NOT gated in V1 — the XFOIL comparisons of V3 must force transition at matched
  x_tr/c on both sides, or the gate compares transition models, not BLs).
- `viscous/ibl3.py` — 6-equation nonlinear surface Galerkin P1 FE: residual +
  analytic Jacobian, Numba kernels per design.md §7 (njit, no object mode,
  `PYFP3D_NOJIT=1` debuggable).

**Gates:**

- [x] **GV1.1 standalone IBL3 verification** (prescribed u_e, no FP coupling):
  (a) laminar flat plate → Blasius: H within ±2 % of 2.59, δ*(x) ∝ √x;
  (b) turbulent flat plate: C_f(Re_θ) within ±5 % of the closure's own reference
  correlation; (c) prescribed decelerating u_e: separation indicator (H rise)
  at the self-similar reference location, band pre-registered; (d) quasi-2D
  invariant: crossflow unknowns (B, Ψ, C_τ2) ≈ 0 (structural lock); (e) surface
  refinement ×2: error drops, measured order recorded.
  **EXECUTED 2026-07-22 → 9 PASS / 2 FAIL** (pre-registered, no re-spec;
  VERDICT + evidence `cases/analysis/v1_ibl3_standalone/`):
  (a) FAIL ×2 — H +3.76 % at outflow and δ* exponent 0.5288, both the
  closure family's own fixed point H*≈2.7083 ≠ Blasius (Stage-2 finding,
  pre-registered as known risk); (b) PASS 0.07 %; (c) PASS P1/P2/P3;
  (d) PASS machine-zero lock, both regimes; (e) first execution FAIL —
  strict decrease broken 100×16→200×32 by an under-damped streamwise 2h
  grid mode at the outflow strip (growth ∝1/h; isotropic D-HB loses inside
  ε∈[0.001,0.01]) → **fixed same-day by the D-HB streamwise-tensor
  follow-up** (anisotropic ν_s = ε_s·max(q)·h_e along s1, ε_s=0.02
  calibrated): errH 4.31e-4→2.21e-4→1.12e-4, order [0.96, 0.99] — PASS;
  SUPG/upwind of the defect convection remains the V3+ upgrade route
  (design doc §9 item 4).
  **V1 CLOSED 2026-07-22 (user-directed)**: (a) ×2 accepted as recorded
  FAIL — closure-family physics, FE matches the same-closure 2-D march to
  <1e-4; no open technical items in V1 scope.

**Prereq:** P6 ✓ + A4 ✓ (both done). V1 touches no wing-body wound and is
independent of the LS-side (b)-class work — parallelizable.

### V2 — Transpiration channel through all three drivers ✓ CLOSED 2026-07-22 · GV2.1 (23 PASS / 0 FAIL / 16 RECORDED)

**Deliverable** (solver plumbing; IBL-independent — parallelizable with V1):

- `viscous/transpiration.py` — δ* → ṁ = ∇_Γ·(ρ_e u_e δ* ê) and the wall-RHS
  assembly (template: `solve/wall_correction.py::assemble_wall_flux_correction_rhs`);
  u_e extraction per-zone per A4 (LE band linear+smoothed, elsewhere quadratic).
  GV2.1 exercises the channel with manufactured ṁ; the δ* → ṁ operator gets its
  first live exercise in V3's coupled gates.
- **Conforming-Newton external-RHS channel** in `solve/newton.py`
  (`R_free -= (Tᵀ b_ext)[free]` class of change): small but touches the core
  solver ⇒ bit-identity when absent + FD Jacobian check (project discipline).
- **Compressible-Picard RHS threading**: `solve_laplace` already has
  `body_source_rhs`, but the compressible Picard (`solve_subsonic` /
  `solve_subsonic_lifting`) does **not** — threading the wall RHS through it is
  needed by GV3.1/GV3.3; the LS path uses the existing `b_base` slot.

**Gates:**

- [x] **GV2.1 transpiration channel exactness**: (a) manufactured blowing on the
  M0 cylinder (Fourier-mode ṁ has an analytic exterior Laplace solution): φ
  error O(h) vs analytic; (b) ṁ = 0 **bit-identical** on ALL drivers (Picard,
  conforming Newton with channel absent, LS `b_base`) — the pre-2026-07-22
  sketch's "δ* = 0 bit-identical" clause lives here now; (c) conforming Newton
  Jacobian stays EXACT under lagged ṁ (FD check).
  **EXECUTED 2026-07-22 → 23 PASS / 0 FAIL / 16 RECORDED** (pre-registered,
  no re-spec; VERDICT + evidence `cases/analysis/v2_transpiration_channel/`):
  (a) relmax 2.1572e-02 > 6.8738e-03 > 2.2062e-03 strict decrease, measured
  orders 1.650/1.640 ≥ 1.0 — the transpiration sign convention pinned by the
  analytic match (a flipped sign lands at O(2)); (b) bit-identical on all
  five legs (`solve_laplace` / `solve_subsonic` / `solve_subsonic_lifting` /
  `solve_newton_lifting` / `solve_multivalued_newton`); (c) J_ff/B
  bit-invariant under lagged b_ext, residual identity 0.0 exact, FD
  6.6e-09–7.2e-08 < 1e-5. No implementation fixes needed during execution.
  **V2 CLOSED 2026-07-22.**

**Prereq:** P6 ✓ + A4 ✓. No BL solve is involved (manufactured ṁ), so V2 is
logically independent of V1; V3 needs both.

### V3 — Loose coupling (2.5-D ladder + fuselage smoke) ✓ CLOSED 2026-07-22 · GV3.1/GV3.2 (2 PASS / 4 FAIL / 23 RECORDED) · GV3.3 (0 PASS / 2 FAIL / 7 RECORDED)

**Deliverable:**

- `viscous/coupling.py` — loose driver: FP solve → u_e → IBL3 solve → ṁ → RHS →
  FP re-solve, with under-relaxation on δ*.
- **Committed XFOIL reference**: generation script + CSV under
  `cases/reference_data/` (δ*/c, C_f, cl at matched M/Re/α/x_tr) — listed here
  because the 2026-07-20 review flagged it as a gate-blocking external artifact.
- **Body-of-revolution smoke mesh** (GV3.3): a standalone fuselage-alone family
  reusing `meshgen/fuselage.py` (`FuselageParams` + `add_fuselage_solid`, full-2π
  revolve) + an M0-style far field; surface tagged `wall` so the single-group
  surface machinery applies unchanged; wake-free, non-lifting (α = 0) — no
  Kutta/wake dependency anywhere in the case.

**Gates:**

- [x] **GV3.1 coupled 2.5-D NACA0012 subsonic vs committed XFOIL reference**
  (matched M/Re/α, forced transition, attached): δ*/c and C_f band OUTSIDE the
  LE/stagnation zone (band pre-registered at execution **and quoted alongside
  the A4 input band ≈2.5 % medium** — viscous-model error and inviscid-input
  error reported separately, never summed silently); viscous Δcl < 0
  (direction) with magnitude vs XFOIL's own viscous decrement (band
  pre-registered). LE-band pointwise comparison is RECORDED, not gated
  (input-limited 4–7 % per A4).
  **EXECUTED 2026-07-22 → 1 PASS / 4 FAIL (+ 18 RECORDED; honest FAIL)** —
  pre-registered + two same-day addenda (cf compared in the freestream frame —
  XFOIL DUMP cf is freestream-normalized, our OUT_CF1 is local; Dirichlet
  inflow-band stations labeled `pinned`, excluded from all statistics);
  VERDICT + evidence `cases/analysis/v3_loose_coupling/`, XFOIL reference +
  generation script `cases/reference_data/naca0012_viscous_xfoil/`.
  **PASS Δcl**: cl 0.2844 → 0.2719, XFOIL's own decrement 0.0230, ratio
  0.542 ∈ [0.5, 2.0]. **FAIL cf** upper/lower worst +43.7 %/+44.8 % at the
  FIRST post-trip station x/c = 0.055 only — XFOIL's e^N intermittency ramps
  cf over a finite run while our flags switch instantaneously; every other
  banded station ≤ 15 % (coarse passes, worst 13.5 %). **FAIL δ*** upper/lower
  worst 27.9 %/27.6 % at x/c = 0.074 — cf matches there (±7 %) so θ is
  consistent; the bias is H in APG on the lower side (+13…+27 %), a
  closure-family difference the 25 % band mostly but not fully absorbs.
- [x] **GV3.2 loose-loop convergence**: ‖Δδ*‖/‖δ*‖ < 1e-3 in ≤ 10 outer
  iterations on the GV3.1 case, under-relaxation factor recorded honestly; one
  transonic-attached 2.5-D point (M ~0.70–0.75) run and RECORDED (iteration
  count + relaxation), not pass/fail — the DN6-predicted near-separation
  divergence risk is measured here, and feeds the V4 skip decision.
  **EXECUTED 2026-07-22 → PASS** (same VERDICT): medium 5 outer iterations at
  **ω = 1.0** (coarse 4); RECORDED transonic point M 0.72 (Newton driver): 4
  iterations, no per-case tuning, cl 0.3764, IBL residual floor 3.2e-6.
- [x] **GV3.3 fuselage smoke** (added 2026-07-22, user-directed) — the minimal
  genuinely-3-D closed-surface transpiration exercise, and Track V's **only
  fuselage-alone item** (no junction, no wake): body of revolution at α = 0,
  subsonic Picard (non-lifting) + loose coupling, forced transition, Re
  recorded. Asserted (bands pre-registered): (a) azimuthal δ* scatter at fixed
  x-stations within band — the surface-FE scatter measure on an unstructured
  triangulation, the genuinely-3-D analogue of GV1.1(d); (b) crossflow unknowns
  (B, Ψ, C_τ2) ≈ 0 (axisymmetric flow on a genuinely 3-D surface). RECORDED,
  not pass/fail: (c) nose/tail stagnation bands seeded per the A4 LE-band
  discipline (linear+smoothed u_e); (d) transpiration on/off fuselage Cp delta
  + tail-cone adverse-gradient H rise — an indicated tail separation is
  recorded and masked, not chased. Headless artifacts per CLAUDE.md §Workflow
  rule 1.
  **EXECUTED 2026-07-22 → 0 PASS / 2 FAIL / 7 RECORDED (honest FAIL)** —
  pre-registered; VERDICT + evidence `cases/analysis/v3_fuselage_smoke/`,
  smoke-mesh generator `cases/meshes/fuselage_bor/`. Three debug rounds to a
  stable closed-body scheme (exactly-singular Newton at the aft pole → tail
  transpiration-sink runaway → Goldstein separation crash ⇒ FINAL: tail-band
  Dirichlet pin + transpiration masking on the pinned band + FP
  non-convergence guard; airfoil path untouched): 10/10 outer iterations, no
  numerical event; mid-body x/L ∈ [0.34, 0.82] axisymmetry excellent
  (σ/μ(δ*) 0.018–0.068, crossflow ratio ~1e-6). **FAIL (a)** σ/μ worst
  0.5533 at x/L = 0.940 (12/63 window stations over band, localized in the
  post-trip ring 0.20–0.33 and the tail cone ≥ 0.82); **FAIL (b)**
  max|B|/max|A| = 0.2631, max|Cτ2|/max|Cτ1| = 0.2295 (maxima tail-cone).
  Loop NOT converged — tail-cone ṁ_max ×5.7 over k = 5→10 = the measured
  loose-coupling stern instability. Medium not executed (pre-registered
  optional; coarse verdict decisive).

**Prereq:** V1 + V2.

### V4 — Quasi-simultaneous coupling ⊘ SKIPPED 2026-07-22 (user-directed; criterion met on GV3.2)

**Deliverable:** Hilbert-integral surface surrogate (`viscous/hilbert.py`);
the BLWF58 method document is the reference description of the approach (NOT
on hand — see the header reference pin; fetch it before opening V4, or skip).

**Gates:**

- [ ] GV4.1 ≥ 30 % fewer coupling iterations than the V3 loose loop on the same
  ladder, OR converges a case the loose loop cannot (near-separation robustness).
- **Skip criterion (concrete):** if GV3.2 passes at ≤ 10 iterations including
  the recorded transonic point without per-case tuning, V4 is SKIPPED — record
  the decision in the ledger and move to V5.
  **Decision inputs measured 2026-07-22:** GV3.2 PASSED at 4–5 iterations
  including the transonic point (M 0.72, 4 iterations, no tuning) — the skip
  criterion is MET by its letter. Counter-evidence from GV3.3 (added
  2026-07-22, after the criterion was written): the loose loop does NOT
  converge at the closed-body stern (tail-cone ṁ growth ×5.7 over k = 5→10)
  — a live case for GV4.1's "converges a case the loose loop cannot" clause,
  relevant if V5/V6 want closed-body geometries. **DECIDED 2026-07-22
  (user-directed): V4 SKIPPED** — the criterion is met by its letter
  (GV3.2: 4–5 iterations incl. transonic M 0.72, no per-case tuning).
  The GV3.3 stern instability is logged as the **reopen trigger**: if
  V5's augmented Newton stalls, or closed-body viscous cases enter scope
  before V5 lands, V4 reopens as the fallback.

**Prereq:** V3.

### V5 — Tight coupling: augmented Newton ✓ CLOSED 2026-07-25 · GV5.0 (16 RECORDED / 0 FAIL) · GV5.1 (9 PASS / 1 FAIL / 36 RECORDED) · GV5.1b/1c/1d · GV5.5 · GV5.2 · GV5.3 · GV5.4 (0 PASS / 1 FAIL / 17 RECORDED) — all five gates executed, close-out pending user adjudication

**Deliverable:** augmented (φ, Γ, BL) Newton on the P8/P14 machinery; coupling
blocks J_φ,BL (∂ṁ/∂BL through the transpiration assembly) and J_BL,φ (∂u_e/∂φ
through the recovery operator chain); GMRES + block preconditioning (AMG on the
φ block / ILU on the BL block; `solve/schur_ls.py` is the structural prototype —
**note** the BL block is O(6 × wall nodes), far bigger than the LS aux thin
band, so exact Schur elimination may not pay: measure, don't assume).

**Gates:**

- [x] **GV5.0 M6 subsonic loose-coupling bridge** (RECORDED, entry check; added
  2026-07-22, user-directed; **EXECUTED 2026-07-23, 16 RECORDED / 0 FAIL** —
  evidence `cases/analysis/v5_m6_bridge/`: the loose loop does NOT converge
  ≤10 at either level — coarse: root-upper-TE separation patch (H 4–5.5)
  feedback runaway, ṁ_max ×12.4 (GV3.3-stern class); medium: refinement
  removes the patch (0 TE nodes H>3.5), runaway gone, but a bounded
  δ* limit cycle (2–12 %/k) never meets tol_ds 1e-3; ΔCL DOWN both
  estimators (coarse −5.2 %/−4.8 %, medium −2.4 %/−2.1 % = input-limited
  under the A4 2.5 % floor); crossflow small (max|B|/|A| ≤ 0.072); tip
  mask validated; δ*(z) CSVs = GV5.3's band feed) — runs on the **V3
  loose driver** (no augmented
  Newton), scheduled here so the 2.5-D → transonic-3-D jump is bridged and the
  crossflow content (Ψ, B equations) gets its **first live 3-D exercise**
  before GV5.3: ONERA M6 (existing `cases/meshes/onera_m6/` family, coarse +
  medium), conforming path, M0.5 / α 3.06° (the committed subsonic convention),
  forced transition, Re recorded (e.g. the M6 experiment 11.72e6). Outputs
  RECORDED, not pass/fail (no δ*(z) truth data exists for M6): δ*(z) spanwise
  distribution at fixed x/c stations (committed CSV + PNG), crossflow-magnitude
  field, ΔCL viscous−inviscid (expected DOWN, direction recorded), 3-D
  loose-loop iteration count vs the GV3.2 2.5-D count; wing-tip band masked
  (tip_taper r_c = 0.05·b_semi). Feeds GV5.3's band pre-registration.
- [x] **GV5.1 exactness + convergence**: both coupling blocks FD-verified
  (project Jacobian discipline; the B19/B31 FD-gate pattern) + quadratic tail on
  the GV3.1 case; outer iterations ≤ half the V3 loose loop. **EXECUTED
  2026-07-23, 9 PASS / 1 FAIL / 36 RECORDED** (pre-registered incl. Addenda
  1–2, amended seed = the loose-converged state, user-adjudicated; evidence
  `cases/analysis/v5_tight_coupling/`, VERDICT + diagnosis committed) — (a) FD
  exactness **PASS** both levels (worst sweet-spot coarse 2.246e-8 seed /
  2.244e-8 endpoint, medium 5.074e-9 seed+endpoint; masked 0/1236 + 0/2460;
  veps omission ≤ 3.0e-8 scaled, decision 5); (b) quadratic tail **HONEST
  FAIL** (medium binding): the polish runs 10 iterations un-converged — F_BL
  pinned from iter 0 at the loose-final IBL floor (medium 1.708e-6, coarse
  3.11e-6), lam → 0, no slope-2 regime (medium p = 0.02/0.50/16.07; F_φ
  resolved at iter 1, 1.16e-7) — mechanism: the IBL steady residual has an
  intrinsic floor on the cond(J_BL,BL) ~ 4e10 near-null manifold (the
  standalone pseudo-time solve stalls there too), NOT a tight-coupling
  defect; (c) N_aug ≤ 2 not met standalone AND as polish (N_polish = 10,
  N_total 14/13 vs loose 4/5). Finding: the committed GV3.1 medium fixed
  point is NOT reproducible — IBL-floor trajectory scatter, three code/env →
  three fixed points cl 0.2217/0.2719/0.2814 (diagnosis
  `results/gv5_1_medium_seed_diagnosis.md`; HEAD-regen seed user-accepted,
  wiring guard |dcl_k0| ≤ 1e-8 PASS 1.309e-9). **Pre-registered
  FD note (added 2026-07-22, user-directed):** J_BL,φ runs through the u_e
  recovery chain, whose per-zone choice (LE band linear+smoothed, elsewhere
  quadratic — per A4) is a discrete, non-differentiable switch: at zone-boundary
  nodes use one-sided differences with a pre-registered tolerance, or make the
  zone switch a smooth weighting; the choice is recorded in the gate.
- [x] **GV5.2 2-D transonic VII vs experiment — EXECUTED 2026-07-25,
  band (b) FAIL + the loose-recipe transonic-limit anatomy
  (`cases/analysis/v5_2_rae2822/`)**: RAE2822
  (`cases/reference_data/rae2822_experiment/`, M0.725/α2.55 + M0.73/α3.19,
  Re 6.5e6): shock location within a pre-registered band of experiment + Cp RMS
  recorded. The NEW 2.5-D RAE2822 mesh family
  (`cases/meshes/rae2822_2.5d/`, coarse 5560 nodes / medium 20790, Cook
  Table 6.1 ordinates) + the A4 TE-wedge pre-check: wedge 9.46°/9.92°
  mesh-crease vs 12.91° ordinate fit ⇒ the ~6° guard clears, no fallback
  (RECORDED). Band (b): every computed shock 0.06–0.10 c DOWNSTREAM of the
  bracket (medium P1 terminal 0.6288 vs [0.495, 0.580]; medium P2 loop
  runaway k = 4, §6 RECORDED) ⇒ **FAIL**; the loose recipe converges only
  1/4 legs (coarse P1) at these points — the displacement-thickness
  feedback is too weak at M ≥ 0.725 (motivates the tight/augmented path).
- [x] **GV5.3 M6 wing direction+magnitude check — EXECUTED 2026-07-25,
  band (b) honest FAIL + band (a) RECORDED input-limited
  (`cases/analysis/v5_3_m6_cp/`)** (re-anchored 2026-07-22, user-directed:
  `cases/reference_data/onera_m6_experiment/` holds Cp only — no experimental
  CL value is committed, so this gate does NOT use the external "experiment ≈
  0.26–0.27" figure). (a) CL moves **down** from the same-mesh k = 0 inviscid
  baseline (anchored to **cl_KJ 0.2823 medium, P14 pressure-Kutta** via the
  W1 wiring guard): Δcl_KJ −2.20 % medium / −1.03 % coarse, direction DOWN
  both estimators but under the A4 input floor (2.5 % medium) ⇒ RECORDED,
  flagged input-limited (most of the medium move arrives with the late
  k = 9–10 separation-patch event). (b) Viscous Cp at the committed 7 span
  stations (`experiment-Cp.dat`, TEST 2308 M0.8395/α3.06): Cp RMS-to-experiment
  does **NOT** decrease vs the same-mesh inviscid baseline — **1/5** unmasked
  stations improved and the pooled RMS INCREASED 0.1288 → 0.1299 (coarse 0/5,
  pooled +0.0024) ⇒ honest **FAIL** (direction verdict; every |ΔRMS| < 0.05
  flagged input-limited per the pre-registered A4 Cp-scale annotation; the
  tip-masked η = 0.96/0.99 stations recorded-only, no anomaly). Neither level
  converges ≤ 10 outer (the GV5.0/GV5.2 loose-loop signature); the FP rescue
  chain load-bearing (10/21, 10/22 stall-accepts); the first-execution medium
  k = 0 wiring-guard fire root-caused to a driver short-circuit and fixed
  under addendum #1. Reading: the 3-D counterpart of GV5.2 — further reads
  belong to the tight/augmented path. (The pre-P14 "0.245
  vs 0.288" framing is superseded — see scope guards.)
- [x] **GV5.4 cost — EXECUTED 2026-07-25 (0 PASS / 1 FAIL / 17
  RECORDED, `cases/analysis/v5_4_cost/`)**: augmented step 22.93 s vs the
  in-session inviscid anchor 3.05 s/step = **7.53×, above the ≤ ~2×
  reference band** (recorded either way per the registration; 4/5
  augmented steps carry capped-GMRES work, annotated; 8-thread session
  walls, not comparable to 16-thread ledger entries). Block
  preconditioner NOT-WORKING at medium (D5): rung-1 block-Jacobi
  diverges (rel_res 5.75e4 — the φ–BL off-diagonal coupling too
  strong); rung-2 exact-BL Schur converges 1/4 steps (2.66e-8 @277 it,
  then stalls vs rtol 1e-8 @300 cap — pure AMG-φ cannot see the
  J_hB·J_BB⁻¹·J_Bh Schur correction). "Measure before Schur" answered:
  splu(J_BL,BL) setup 1.8 s (elimination cheap — Krylov convergence is
  the bottleneck), AMG-φ 0.2 s / ILU-BL 2.5 s. System 124,216 DOFs;
  W1/W2/W3 PASS (cl_p 0.26429 vs the addendum-#4 P14-locked 0.2646 =
  0.116 %; FD median φ 8.7e-12 / Γ 7.3e-12 / BL 0). Reading:
  wing-scale augmented Newton needs a stronger (Schur-aware, (A,Ψ)
  structured) reduced-space preconditioner before cost reads into the
  ≤ ~2× band — the honest negative the registration asked to record;
  the EW-forcing variant registered-not-opened (user adjudication).
  Design record `docs/design_track_v.md` §20.
- [x] **GV5.5 TE-band (B, δ) formulation — breaking the IBL floor
  (STANDALONE ITEM, registered 2026-07-24, user-directed; OPENED and
  EXECUTED 2026-07-24 — 2 PASS / 1 FAIL / 9 RECORDED)**. Target: the
  steady-IBL residual floor (max-norm coarse 3.154e-6 / medium 1.710e-6,
  the committed loose-final floors) localized by the committed diagnosis
  (`cases/analysis/v5_ibl_floor/` findings Q5 = design doc §13 item 3) in
  the **TE-band (B, δ) equations**, lying essentially entirely inside J's
  range — a formulation floor (Q7: the pseudo-time controller bottoms out
  with the residual frozen), not a solver limitation, and not crossable by
  globalization alone (GV5.1b: scaling + damping delivered and exact, μ
  inert, floor intact). Route chosen at opening per the registration:
  (a) TE natural-outflow discretization work FIRST, as the row-level
  variant **V1 = TE-outflow row replacement** (first-order extrapolation:
  δ-carrier row 6i+0 `R = δ_i − δ_up`, H-carrier row 6i+2
  `R = H_i − H_up`, exact Jacobian rows, CSR in-pattern guard, default-OFF
  flag `te_extrapolate`; `te_outflow_pairs` supplied by the case layer).
  **Outcome (cases/analysis/v5_5_te_floor/VERDICT.md): V1 does NOT break
  the floor** — the variant system sees the amended seed at residual
  9.8/4.8 (the replaced rows measure the natural TE jump), the pseudo-time
  stalls (all steps rejected, cfl → 1e-3 floor), and the BINDING m2
  (original-system residual at the V1 terminal) lands 5554× (coarse) /
  245998× (medium, vs the seed's own flag-OFF floor per the scatter
  clause) ABOVE the floor — the pre-registered "worse" clause; the damage
  peaks at the LE suction zone (x_c ≈ 0.027, F_B/F_Psi), not the TE.
  Band (a) FD PASS both levels; V0 control coarse bit-close, medium on
  the 4th fixed point (scatter clause handled as pre-registered); guards:
  plate H bands PASS flag-ON, loose smoke flag-ON coarse RED (cl_rel
  2.62% > 2.5%, cap-hit) / medium marginal PASS (2.49%); tight polish
  secondary read no floor break either (7.32e-5 / diverged 3.98 vs
  committed finals 3.07e-6 / 1.708e-6). The flag stays default-OFF
  (legacy paths bit-identical); the escalation ladder (upwind
  boundary-flux (a)-variant / closure regularization (b)) stays
  registered-not-opened — opening = user's adjudication. Executed under
  the temporary 8-thread session constraint. Prereq: none
  beyond the committed diagnosis; the GV5.1c window read is informative
  but NOT binding for the opening.

**Prereq:** P8 ✓ + P14 ✓ + V3. **Wing-body VII is explicitly OUT of V5 scope**
(scope guards below).

### V6 — Wake-sheet IBL correction ☐ (continuation of V1's data layout, not an independent solver)

Same 6 equations with wake closure relations; the wake unknowns were reserved in
V1's layout. δ*_wake enters as the wake-sheet RHS mass source; TE thickness
continuity δ_wake(TE) = δ*_upper + δ*_lower.

**Gates:**

- [x] **GV6.0 design adjudication (BEFORE code, user-adjudicated — OPENED
  2026-07-25, RULED 2026-07-25 (user); document
  `cases/analysis/v6_0_design_adjudication/DESIGN_ADJUDICATION.md`)**: the
  LS-path sheet-source mechanism. **Ruling: Option A — V6 closes
  conforming-only; the LS sheet-source leg = recorded follow-up (its own
  pre-registration if/when opened). δ*_wake producer = (i) prescribed
  (TE-continuity construction + straight-wake relaxation, pinned in the
  GV6.1 pre-registration); the solved wake IBL (the "same 6 equations"
  vision) = recorded follow-up producer (ii), opened only if the GV6.2
  measured effect is significant vs the A4 input band.** Branch push + PR
  deferred to bundle the adjudication + the GV6.1 pre-registration (user).
  (Survey basis: `pyfp3d/wake/` has no sheet-surface integration machinery;
  the conforming path needs NO new mechanism — explicit `wake_minus/plus`
  faces + slave→master folding IS the weak-form flux channel, currently
  unconsumed; the path-independent physics gaps = wake closures, TE
  confluence, the wake IBL state.)
- [ ] **GV6.1 conforming sheet source**: δ*_wake enters via the
  `constraints/wake.py` reduce RHS (Tᵀ b_wake); TE thickness continuity
  asserted; δ*_wake = 0 **bit-identical**.
- [ ] **GV6.2 measured effect**: wake-IBL on/off cl (and TE-region Cp) delta on
  the GV3.1 case, direction-checked against XFOIL's wake modelling; RECORDED
  with the A4 input band quoted.

**Design constraints (unchanged from DN2 §4.5):** TE kink absorbed by Drela
local-basis adaptation; **straight wake + mass-transpiration relaxation, no
geometric relaxation**.

**Prereq:** V3 (+ GV6.0 adjudication for any LS leg).

### Scope guards (re-based 2026-07-22; the DN2 §9 / DN6 §13–14 envelope stands)

- **Validity envelope**: attached / mildly-shocked flow (M_shock ≲ 1.3); not
  massive or shock-induced separation. The M6 M0.84 shock sits at the envelope
  edge; the wing-TIP singularity zone (local M_max 2.5+) is OUTSIDE it — IBL
  seeding/comparison bands must mask the tip band (on the conforming path the
  production `tip_taper` r_c = 0.05·b_semi band is the natural mask).
- **VII does not close the inviscid-discretization CL gap — updated numbers**:
  after P13/P14 the remaining inviscid gap to the FP reference is ≈ 0.5 %
  (cl_KJ 0.2823 medium / 0.2866 fine vs 0.288), so the remaining delta to
  experiment is now genuinely viscous-dominated; viscosity moves CL **down**.
  **No experimental CL value is committed** (`onera_m6_experiment/` holds Cp
  only) ⇒ the M6 viscous gate GV5.3 anchors on the committed Cp data, not on
  any external CL figure. (History, superseded: the pre-P14 "0.245 vs 0.288 →
  sharp-TE/LE P1 floor → P9/P11" attribution — P11 measured NEGATIVE 2026-07-19,
  and P14/G14.7 showed 69 % of the old 0.019 gap was Kutta-estimator bias. See
  the 2026-07-20 review §3.)
- **Input-error discipline (A4, standing rule for every V-gate)**: every
  viscous-vs-reference comparison quotes the A4 inviscid u_e input band (medium
  ≈ 2.5 % peak-relative; LE/stagnation 4–7 %) alongside the viscous
  discrepancy. Tight LE-band comparisons are input-limited by construction.
  Only the **unchosen** G1.6 route (b) — isoparametric P2 wall layer — could
  raise the input band to O(h²): if a V-gate fails for input reasons alone,
  that is the recorded escalation route, NOT a viscous-model fix (G1.6 Option C
  close-out note, 2026-07-22).
- **Wing-body VII (applying V5/V6 to the M2 wing-body) is DEFERRED** until the
  LS-side wing-tip (b)-class is cured or explicitly accepted (B30 attribution;
  B31 cured the conforming side via tip_taper, LS C-class closed negative
  C1/C3) — otherwise wing-body viscous gates re-enter "viscous model vs known
  inviscid wound" attribution confusion. The M6 **wing** gate GV5.3 is
  unblocked (conforming, tip_taper production since B32).
- **Reynolds number and transition are new physical inputs** (the FP solver has
  neither): Re enters only the closures; V1 ships forced transition, free
  transition (e^N) is a recorded follow-up. Gate comparisons must match Re and
  x_tr explicitly.
- The V1–V3 core sequence is parallelizable with the remaining Track-B/LS work
  (it depends only on P6 + A4 and touches no wing-body wound), but the full
  sequence is a large, self-contained solver effort (6-equation nonlinear
  surface FE + closures + coupling): **budget it like a Track-P phase, not a
  side task** — the 2026-07-22 three-way split exists precisely to give the
  effort Track-P-sized checkpoints.

---


## Progress ledger

### Track V — viscous–inviscid interaction

Track status: **◐ IN PROGRESS — V1 OPENED 2026-07-22** (gates re-spec'd at
opening against the shipped B32/A4 state; same day, user-directed: V1 split
three ways → V1 core / V2 channel / V3 coupling, old V2–V4 → V4–V6, GV5.3
re-anchored on the committed M6 experiment Cp, Drela 2013 reference pinned —
still zero implementation, so no gate result is affected). Design 2026-07-09/10
(DN2 + DN6, historical; two known DN6 traps annotated in the header re-spec
block). Validity envelope: attached / mildly-shocked flow. **VII does not close
the inviscid-discretization CL gap** — the inviscid baseline is now clean to ≈
0.5 % (P13/P14), so the remaining delta to experiment is the viscous target
(assessed on the committed Cp data); direction is DOWN.

- V1 — **✓ CLOSED 2026-07-22 · GV1.1 9 PASS / 2 FAIL** — IBL3
  solver core shipped (`viscous/surface_mesh.py`, `closures.py`, `ibl3.py`;
  wake unknowns reserved in the data layout). GV1.1 verdict + evidence:
  `cases/analysis/v1_ibl3_standalone/VERDICT.md`; implementation record:
  `docs/design_track_v.md` §9. (a) ×2 accepted as recorded FAIL at closing
  (user-directed) = closure-family fixed point, FE matches the same-closure
  2-D march to <1e-4; (e) first-run
  FAIL = streamwise 2h grid mode → fixed by the D-HB streamwise-tensor
  stabilization (ε_s=0.02), PASS; (b)(c)(d) PASS. Prereqs P6 ✓ + A4 ✓;
  no wing-body contact.
- V2 — **✓ CLOSED 2026-07-22 · GV2.1 23 PASS / 0 FAIL / 16 RECORDED** —
  transpiration channel through all three drivers shipped
  (`viscous/transpiration.py`; conforming-Newton external-RHS channel +
  compressible-Picard RHS threading; LS rides the existing `b_base`). Gate
  GV2.1 (manufactured-blowing exactness + ṁ=0 bit-identity on all drivers +
  FD) verdict + evidence: `cases/analysis/v2_transpiration_channel/VERDICT.md`
  — (a) orders 1.650/1.640 ≥ 1.0, sign pinned analytically; (b) five legs
  bit-identical; (c) Jacobian bit-invariant + FD exact under lagged ṁ.
- V3 — **✓ CLOSED 2026-07-22 · GV3.1/GV3.2 2 PASS / 4 FAIL / 23 RECORDED ·
  GV3.3 0 PASS / 2 FAIL / 7 RECORDED (honest FAILs)** — loose coupling
  shipped (`viscous/coupling.py`: CouplingCase builders + run_loose_coupling
  outer loop; committed XFOIL reference
  `cases/reference_data/naca0012_viscous_xfoil/`; BoR smoke-mesh generator
  `cases/meshes/fuselage_bor/`). Gate verdicts + evidence:
  `cases/analysis/v3_loose_coupling/VERDICT.md` (GV3.1/3.2 — PASS Δcl ratio
  0.542 ∈ [0.5, 2.0] vs XFOIL's own decrement, PASS loop convergence 5 iters
  ω = 1.0 (transonic M 0.72 record: 4 iters, no tuning); FAILs localized:
  cf +44 % at the first post-trip station only (XFOIL e^N ramp vs
  instantaneous switch), δ* H-family offset ≤ 27.9 % at x/c = 0.074) and
  `cases/analysis/v3_fuselage_smoke/VERDICT.md` (GV3.3 — closed-body scheme
  stabilized through three debug rounds: tail-band pin + transpiration
  masking + FP guard; mid-body axisymmetry excellent, FAILs at the post-trip
  ring and the tail cone; loop NOT converged — measured stern instability,
  V4 decision input). Also fixed en route: IBL3 local-basis crossflow
  leakage (25.9/0.15 → 1.8e-4/1.6e-3, `viscous/ibl3.py`). V4 skip criterion
  MET by its letter (GV3.2), counter-evidence logged (GV3.3 stern).
- V4 — **⊘ SKIPPED 2026-07-22 (user-directed)** — quasi-simultaneous
  coupling (`viscous/hilbert.py`, BLWF58 reference — NOT on hand). Skip
  criterion met by its letter on GV3.2 (4–5 iterations incl. transonic
  M 0.72, no tuning). **Reopen trigger** (logged from GV3.3): V5's
  augmented Newton stalls, or closed-body viscous cases enter scope
  before V5 lands.
- V5 — ✓ CLOSED 2026-07-25 (all five gates executed; close-out pending user
  adjudication) — tight coupling: augmented (φ, Γ, BL) Newton on
  P8/P14. **Entry check GV5.0 ✓ EXECUTED 2026-07-23** (16 RECORDED / 0 FAIL;
  `cases/analysis/v5_m6_bridge/`): the bridge answer is that the loose loop is
  NOT sufficient on the 3-D lifting wing — coarse runs away on a root-upper-TE
  separation patch (ṁ_max ×12.4, the GV3.3-stern/Veldman class), medium
  resolves the patch away but sits in a bounded unconverged δ* cycle; ΔCL DOWN
  both estimators at both levels (medium −2.4 % input-limited); crossflow
  small; tip mask validated; δ*(z) CSVs feed GV5.3's bands. New machinery:
  `viscous/coupling.py::build_wing_case` (LE-band laminar pin per local x/c,
  both TE natural outflow, root symmetry natural, tip band z > 0.95·b_semi
  pinned + ṁ-masked via the GV3.3 machinery) + `tests/test_v5_wing_case.py`
  (5). **GV5.1 ✓ EXECUTED 2026-07-23** (9 PASS / 1 FAIL / 36 RECORDED;
  `cases/analysis/v5_tight_coupling/`, VERDICT + PRE_REGISTRATION Addenda
  1–2): the exact augmented Newton is delivered and FD-verified at both
  levels (worst sweet-spot 2.2e-8 coarse / 5.1e-9 medium; new machinery
  `viscous/tight.py` + `viscous/tight_driver.py`, `tests/v5_state.py` + 3
  tight test files); the quadratic tail is HONEST FAIL — F_BL pins at the
  IBL floor from iteration 0 (medium 1.708e-6 / coarse 3.11e-6), an
  intrinsic floor of the steady IBL residual on the cond(J_BL,BL) ~ 4e10
  near-null manifold (the standalone pseudo-time solve stalls there too),
  NOT a coupling defect; N_aug ≤ 2 not met standalone nor as polish
  (N_total 14/13 vs loose 4/5). Finding: the committed GV3.1 medium fixed
  point is not reproducible (IBL-floor trajectory scatter; diagnosis
  committed; HEAD-regen seed user-accepted). **IBL-floor follow-up
  diagnosis ✓ EXECUTED 2026-07-24** (14 RECORDED, no bands;
  `cases/analysis/v5_ibl_floor/`): the near-null cluster PERSISTS at the
  loose-converged states, carried by the turbulent (A, Ψ) variables
  mid-chord → TE; the raw cond 4e10–4e13 is MOSTLY a scaling artifact
  (row+col equilibration → 2e4/7e5/1e7, sub-1e-6 count 501/500/1082 →
  0/0/2 — no exact null directions) but a genuine scaled (A, Ψ)
  stiffness of 1e5–1e7 remains = the real GV5.1b/GV5.4 target; the F_BL
  floor lives in the TE band (B, δ) equations essentially entirely
  INSIDE J's range; the closure-floor active set is EMPTY (the
  floor-active-null hypothesis dead), eps_diff ×4 moves the floor ≤ 6 %
  (not an artificial-viscosity truncation), and the pseudo-time
  controller bottoms out with the residual frozen = a formulation floor
  that globalization alone cannot pass. **GV5.1b ✓ EXECUTED 2026-07-24**
  (2 PASS / 0 FAIL / 7 RECORDED adjudicated; 1 PASS / 1 FAIL / 7
  RECORDED as executed, preserved in commit 1c55906;
  `cases/analysis/v5_1b_scaled_newton/`,
  VERDICT + PRE_REGISTRATION committed 8b7793f): the scaled + damped
  Newton machinery is delivered and exact — solver-internal row/column
  equilibration + Levenberg diagonal damping + a floor-reached stop,
  flags default OFF (legacy path bit-reproduces the committed histories;
  new `tests/test_v5_tight_scaled.py` (8), tight fleet 28 passed twice).
  Band (a) suite PASS both levels; the medium live-seed e2 identity
  reads 1.96e-10 vs a ≤ 1e-10 threshold chosen at implementation time
  (NOT pre-registered) = SuperLU pivot-order roundoff through
  cond(J) ~ 1e10, the backward-error floor — **adjudicated PASS
  2026-07-24 (user) under the cond-aware read** tol = max(1e-10,
  10·κ₁(J)·eps), a ~1e-5-class bound at κ₁ ~ 1e10, a ~4-decade margin
  (VERDICT §3; run.py now computes the tolerance live from a κ₁
  one-norm estimate). Band (b): the amended seeds sit INSIDE the
  10× floor band from iter 0 (F_BL = 1.00× the floor), no above-band
  contraction segment exists by construction → the pre-registered
  fallback: medium terminates floor_reached at iter 5 (replacing
  GV5.1's 10-step λ-collapse crawl) at the same merit
  (9.074e-11 ≈ 9.025e-11); coarse ends below GV5.1 and still
  descending (merit 2.044e-10 < 2.068e-10); the k=1 standalone
  descends markedly deeper (F_BL 3.268e-6, −31 % vs the k1seed; merit
  2.3× below). Band (c): coarse 10 vs 8 NOT met, medium 5 vs 10 met
  (degenerate band-entry iter 0). μ rejection-retries = 0 across all
  three runs — the scaling is the active ingredient, the damping arm
  inert at these states. The window question is REFRAMED, not
  answered: it needs an above-band seed (early loose iterate /
  perturbed δ*) = candidate GV5.1c; breaking the floor itself = the
  TE-band (B, δ) formulation work, queued. VERDICT
  `cases/analysis/v5_1b_scaled_newton/VERDICT.md`, design record
  `docs/design_track_v.md` §14. **GV5.1c ✓ EXECUTED 2026-07-24**
  (2 PASS / 1 FAIL / 7 RECORDED; `cases/analysis/v5_1c_above_band_window/`,
  VERDICT + PRE_REGISTRATION committed 1e90d59 pre-execution; the
  above-band-seed window read, user-directed): the above-band seeds
  delivered as pre-registered (the amended seed + δ×(1+ε) at the free
  BL nodes, ε = 1e4 by the deterministic calibration bisection → seed
  F_BL 3.219e-1 coarse / 1.819e-1 medium ≈ 1e4× the floor band) and
  the pre-floor slope-2 window is MEASURED — **no quadratic regime
  anywhere above the floor**: the clean-descent steps are
  line-search-capped halvings (λ = 0.5 → p = 1.00 by construction,
  the backtracking cap, not Newton asymptotics) and the trajectory
  STALLS mid-range (F_BL ~ 3e-2 → 1.3e-2 / 2.2e-2 over 10
  iterations), never reaching the band (4262× / 12867× the floor at
  the cap); binding medium median p = 0.56 → honest FAIL (coarse
  1.00 recorded); regression slopes 0.75/0.62; μ rejection-retries
  0 again (the line search carries all the globalization). Band (a)
  PASS both levels with the cond-aware e2 tolerance pre-registered
  (e2 2.06e-9 / 2.40e-9 vs 3.9e-2 / 5.2e-2). New finding: the
  tight-Newton obstacle is not only the formulation floor — a
  mid-range descent barrier sits 3–4 decades above it; whether a
  quadratic basin exists ADJACENT to the floor = the near-band-seed
  follow-up question (candidate GV5.1d, user adjudication).
  Executed under the temporary 8-thread session constraint (runner
  default 16; wall times flagged non-comparable); the medium fixed
  point scattered AGAIN at 8 threads (a 4th fixed point cl
  0.28245999, unperturbed F_BL 1.824e-6 = 1.07× floor; coarse
  bit-identical). Design record `docs/design_track_v.md` §15.
  **GV5.1d ✓ EXECUTED 2026-07-24** (2 PASS / 1 FAIL / 7 RECORDED;
  `cases/analysis/v5_1d_near_band_window/`, VERDICT + PRE_REGISTRATION
  committed pre-execution; the near-band seed, user-directed): the
  seeds calibrated INTO the near-band windows as pre-registered (T1 =
  [1e-4, 1e-3]; coarse ε = 10 → F_BL 1.711e-4 = 5.42× the band,
  medium ε = 56 → 6.02e-4 = 35×; the T2 escalation never fired — ≥ 3
  above-band triples on both T1 legs) and the near-band window is
  MEASURED — **no quadratic basin adjacent to the floor either**:
  coarse halves once (the λ = 0.5 cap) then crawls (λ → 6e-5,
  ≤ 0.03 dex/step) to 7.59e-5 = 24× the floor, never entering the
  band; medium's FIRST accepted step moves F_BL AWAY from the band
  (6.0e-4 → 9.8e-4 — the merit bought by block rebalance, not BL
  descent) then crawls to 8.43e-4 = 493× the floor; binding medium
  median p = 1.17 → honest FAIL (coarse 0.35 recorded); regression
  slopes 0.15/0.88; μ rejection-retries 0 for the third time. Band
  (a) PASS both levels (cond-aware e2, ~12-decade margin). The stall
  is NOT a mid-range barrier with a basin below it: the flat/ragged
  merit neighborhood extends DOWN to within ~1.5 decades of the
  floor — basin hunting is exhausted (GV5.1b/1c/1d), and **GV5.5 is
  now the only registered open route for the floor itself**.
  Executed under the temporary 8-thread session constraint; medium on
  the same 4th fixed point as GV5.1c (cl 0.28245999; coarse
  bit-identical). Design record `docs/design_track_v.md` §16. V5
  stays **OPEN**; the V4-reopen trigger stays parked.
- **GV5.5 TE-band (B, δ) formulation — EXECUTED 2026-07-24 (2 PASS /
  1 FAIL / 9 RECORDED, `cases/analysis/v5_5_te_floor/`)**: route (a)
  variant V1 = TE-outflow row replacement (first-order extrapolation on
  the δ-carrier row 6i+0 and the H-carrier row 6i+2, exact Jacobian rows,
  default-OFF flag `te_extrapolate`) does **NOT** break the floor — the
  amended seed sits at variant residual 9.8 (coarse) / 4.8 (medium), the
  pseudo-time stalls with all steps rejected, and the BINDING m2
  (original-system residual at the V1 terminal) = 1.752e-2 = **5554×**
  the floor (coarse) / 4.487e-1 = **245998×** the seed's own flag-OFF
  floor (medium; the 8-thread scatter clause fired — the 4th fixed point
  cl 0.28245999 again, handled as pre-registered) — the pre-registered
  "worse" clause. Band (a) FD PASS both levels (≤1.8e-7); the damage
  peaks at the **LE suction zone** (x_c ≈ 0.027, F_B/F_Psi rows), not
  the TE; tight-polish secondary read no floor break either (coarse
  7.32e-5 vs committed 3.07e-6; medium diverged 3.98). Guards: plate H
  bands PASS flag-ON; loose smoke flag-ON coarse RED (cl_rel 2.62% >
  2.5%, cap-hit) / medium marginal PASS (2.49%). The flag stays
  default-OFF (legacy bit-identical); the escalation ladder (upwind
  boundary-flux (a)-variant / closure regularization (b)) stays
  registered-not-opened — opening = user's adjudication. Design record
  `docs/design_track_v.md` §17. Next = **GV5.2/GV5.3/GV5.4** (user
  sequencing 2026-07-24: GV5.1d → GV5.5 → GV5.2–5.4). Remaining: RAE2822
  transonic VII vs
  committed experiment (GV5.2; needs the 2.5-D RAE2822 mesh family + A4
  TE-wedge pre-check), M6 CL-down + Cp-RMS-down vs committed experiment Cp
  (GV5.3 — anchored on committed data only, no external CL figure), cost
  recorded (GV5.4). Block precond: AMG-φ / ILU-BL, `schur_ls.py` prototype;
  BL block is NOT thin — measure before Schur. Wing-body VII out of scope
  (deferred, see scope guards).
- **GV5.2 RAE2822 transonic VII vs committed experiment — EXECUTED
  2026-07-25 (band (b) FAIL + 3 recipe-limit RECORDED + 2
  outside-envelope RECORDED, `cases/analysis/v5_2_rae2822/`)**: the
  loose GV3.1 recipe (ω = 1.0, ≤ 10 outer, tol_ds = 1e-3) with the GV3.2
  Newton-driver protocol at the two dataset-labeled points (P1
  M 0.725/α 2.55, P2 M 0.73/α 3.19, Re 6.5e6, x_tr/c = 0.03). Band (a):
  TE wedge 9.46° coarse / 9.92° medium mesh-crease (A4 method) vs 12.91°
  ordinate fit; quadratic recovery available ⇒ no fallback. **Band (b)
  FAIL** (medium binding): P1 terminal x_shock 0.6288 outside
  [0.495, 0.580] (leg non-converged at the k = 10 cap — ds_change_rel
  oscillates, mdot grows, IBL capped every outer), P2 loop runaway at
  k = 4 (mdot_max = 1.59, the GV3.3 class; §6 recipe-limit RECORDED);
  coarse recorded: P1 converged (7 outer) but 0.6122 out of band, P2
  capped + IBL exactly-singular warning, 0.6520 out. Every computed
  shock sits 0.06–0.10 c DOWNSTREAM of the experimental bracket,
  worsening with Mach/α and NOT improving with mesh refinement (coarse →
  medium P1 moves further aft) — the miss is not a coarse-mesh artifact.
  (c) Cp RMS RECORDED 0.185/0.146, 0.176/0.118, 0.265/0.129
  (upper/lower) — dominated by the shock displacement + the over-deep LE
  suction; lower side markedly better. Outside-envelope RECORDED:
  coarse P2 M_peak 1.365, medium P1 1.306 (> 1.3). Execution mechanics
  (all pre-registered in addenda BEFORE each (re-)execution): the
  RAE2822 reflex camber broke two global-y assumptions — `cut_wake`'s
  Kutta probe gained a TE-wedge bisector-normal fallback (fires only
  where the old code raised) and the runner's Cp side split switched to
  the outward-normal idiom; the FP driver gained a cheap→deep rescue
  chain (strict 1e-10 → `solve_newton_transonic` Mach continuation →
  honesty-guarded stall acceptance) against the M ≥ 0.725 shock-cell
  plateaus — load-bearing everywhere (2 continuation cold starts per
  level, 2/10 → 10/22 stall-accepts per leg). Reading: the
  displacement-thickness feedback through the loose update is too
  weak/slow at these transonic points — the next transonic-VII reads
  should come from the tight/augmented path, not further loose-loop
  tuning. Design record `docs/design_track_v.md` §18. Next =
  **GV5.3/GV5.4** (user sequencing 2026-07-24). Executed under the
  temporary 8-thread session constraint (~43 min for 4 legs).
- **GV5.3 M6 wing direction+magnitude check vs committed Cp — EXECUTED
  2026-07-25 (band (b) honest FAIL + band (a) RECORDED input-limited,
  0P/1F/17R, `cases/analysis/v5_3_m6_cp/`)**: the loose GV3.1 recipe on
  the GV5.0 wing case at TEST 2308 (M 0.8395/α 3.06, Re_MAC 11.72e6,
  x_tr/c 0.05), the P14 transonic FP recipe verbatim (NEWTON_M6_RECIPE
  imported; the k = 0 inviscid baseline anchored to the committed P14
  cl_KJ 0.2823 medium / 0.2688 coarse via the W1 wiring guard, PASS both
  levels). **Band (a) RECORDED input-limited**: Δcl_KJ −2.20 % medium /
  −1.03 % coarse — direction DOWN both estimators (cl_p −2.40 %/−1.35 %,
  the expected viscous decambering) but under the A4 2.5 % floor; most
  of the medium move arrives with the late k = 9–10 separation-patch
  event (ds_max 0.0021 → 0.0039, mdot 0.007 → 0.106, bounded).
  **Band (b) honest FAIL** (medium binding): the viscous Cp does NOT
  move toward the committed 7-station experiment — 1/5 unmasked stations
  improved (only the root η = 0.20), pooled RMS 0.1288 → 0.1299
  INCREASED (coarse 0/5, pooled +0.0024); every |ΔRMS| < 0.05 flagged
  input-limited per the pre-registered A4 Cp-scale annotation — a
  DIRECTION verdict; the tip-masked η = 0.96/0.99 stations recorded-only
  (no anomaly). (c) RECORDED: neither level converges ≤ 10 outer (the
  GV5.0/GV5.2 signature — medium sits in a tight limit cycle k = 1–8
  before the k = 9 event); the IBL residual floor 1.9e-6 medium; the
  pre-registered FP rescue chain load-bearing (10/21 coarse, 10/22
  medium stall-accepts; the unconditional P14 ramp after addendum #1).
  Execution narrative: the first execution's medium k = 0 FAILED W1
  (cl 0.226 — the GV5.0 bridge's M0.5 short-circuit returned the
  half-converged M0.70 probe seed and the ramp never ran); root-caused
  by a measured diagnostic (the P14 ramp from the same failed seed lands
  on the anchored branch to 9 digits — NOT an 8-thread branch scatter)
  and fixed under addendum #1; the loop itself converged 9 outer even
  from the poisoned seed. Reading: the 3-D counterpart of GV5.2 — at
  transonic conditions the loose displacement-thickness feedback is too
  weak to repair the inviscid-family mismatch (shallow LE suction, aft
  shock); further reads belong to the tight/augmented path, not more
  loose-loop tuning. No follow-up opened (user adjudication). Design
  record `docs/design_track_v.md` §19. Next = **GV5.4** (user sequencing
  2026-07-24). Executed under the temporary 8-thread session constraint
  (coarse 1721 s + medium 12479 s).
- **GV5.4 augmented-step cost on M6 medium — EXECUTED 2026-07-25 (0
  PASS / 1 FAIL / 17 RECORDED, `cases/analysis/v5_4_cost/`)**: the
  tight/augmented Newton path (the GV5.1b scaled+damped driver with an
  injectable `step_solve` solve callback — a library change, default
  `None` = splu bit-identical) measured on the 124,216-DOF W2 system
  (62,820 φ + 166 Γ + 61,230 BL). **Band (a) RECORDED**: augmented
  step 22.93 s vs the in-session inviscid anchor 3.05 s/step = **7.53×
  — above the ≤ ~2× reference band** (the registration records the
  number either way; 4/5 augmented steps carry capped-GMRES work,
  annotated in the VERDICT; walls under the temporary 8-thread session
  constraint, not comparable to 16-thread ledger entries). **Band (b)
  FAIL (D5) — the block preconditioner does NOT work at medium**:
  rung-1 block-Jacobi (AMG-φ + ILU-BL) diverges (rel_res 5.75e4 — the
  φ–BL off-diagonal coupling too strong for a block-diagonal
  approximation); rung-2 exact-BL Schur (AMG-φ on the Schur operator)
  converges 1/4 steps (2.66e-8 @277 iterations, then stalls
  2.07e-7/6.13e-5/2.52e-6 vs rtol 1e-8 at the 300 cap — a pure AMG-φ
  cycle cannot represent the J_hB·J_BB⁻¹·J_Bh Schur correction). The
  registered "measure before Schur" question answered: splu(J_BL,BL)
  setup is only 1.8 s (the elimination is cheap — Krylov convergence
  is the bottleneck, not the factorization); AMG-φ setup 0.2 s, ILU-BL
  2.5 s. Guards: W1 cl_p 0.26429 vs the addendum-#4 P14-locked anchor
  0.2646 (0.116 %), W2/W3 PASS (FD median φ 8.7e-12 / Γ 7.3e-12 /
  BL 0); the IBL floor trajectory intact (merit pinned 9.35e-9, |F_BL|
  2.888e-6, steps 2–5 λ ~ 1e-4 strictly decreasing). Execution
  narrative: four addenda — #1 W1 cl_p tolerance scoped medium-only,
  #2 the M0.70 seed chain not raising on non-convergence, #3 the seed
  chain = the A1 conf_newton verbatim, #4 the W1 anchor re-pinned from
  the stale A1 0.26918 to the P14 probe G8.2 lock 0.2646. Reading: the
  pre-registered honest negative — wing-scale augmented Newton needs a
  stronger (Schur-aware, (A,Ψ)-structured) reduced-space preconditioner
  before the augmented-step cost can read into the ≤ ~2× band; the
  EW-forcing variant stays registered-not-opened (user adjudication).
  Design record `docs/design_track_v.md` §20. Next = **V5 close-out —
  all five gates executed; user adjudication**. Executed under the
  temporary 8-thread session constraint.
- V6 — ☐ — wake-sheet IBL correction, a continuation of V1's data layout (wake
  unknowns reserved). GV6.0 LS sheet-source design adjudication BEFORE code
  (conforming needs no new mechanism; may close conforming-only), GV6.1
  conforming sheet source + δ*_wake=0 bit-identity, GV6.2 measured on/off
  effect. Straight wake + mass-transpiration relaxation, no geometric
  relaxation.
