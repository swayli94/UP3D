# pyFP3D Roadmap — Track V (viscous–inviscid interaction V1–V4)

> Split verbatim from `docs/roadmap.md` on 2026-07-15 (content unchanged; only
> this header and the ledger heading were added). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md).

## Track V — Viscous–inviscid interaction (designed 2026-07-09/10; **V1 ◐ OPENED 2026-07-22**)

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

**Prerequisite state at opening (all measured, none pending):**

- **P6 ✓** smoothed wall tangential gradients
  (`post/surface.py::smooth_wall_tangential_gradients`) = the u_e / du_e/ds input path.
- **A4 ✓ (2026-07-22)** u_e inviscid input error band, analytic ground truth:
  medium smooth-wall ≈ **2.5 % peak-relative / 0.04·U∞ max-norm / 0.012·U∞ rms, O(h)**;
  **LE/stagnation band 4–7 % @ medium** (= the IBL seeding / du_e/ds zone, the least
  trustworthy input zone); linear-vs-quadratic recovery has no universal winner
  (~1 % region-dependent) → per-zone choice, LE band linear+smoothed; the sub-6°
  quadratic-recovery guard does **not** fire on NACA0012 (16° TE) — thin-TE
  (RAE-class) airfoils must re-check (GV3.2 clause).
- **P8/P14 ✓** conforming coupled (φ,Γ) Newton + pressure-Kutta (V3's augmentation
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
  integration (V4 only, design adjudication GV4.0).
- **Reference data committed**: `cases/reference_data/rae2822_experiment/`
  (M0.725/α2.55 + M0.73/α3.19, Re 6.5e6), `naca0012_experiment/`,
  `onera_m6_experiment/`.

> **Naming disambiguation.** Track V phase IDs V1–V4 are development *phases*;
> they are unrelated to the §10 validation-ladder case IDs V0–V6 in design.md
> (e.g. "V6 < 1%" remains the Γ-consistency metric). Gates here are GV<phase>.<n>.

### V1 — IBL3 solver core + loose coupling (2.5-D ladder) ◐ OPENED 2026-07-22

**Deliverable** (Track-P-sized; budget it like one):

- `viscous/surface_mesh.py` — compact wall-surface DOF numbering on the existing
  wall triangulation (`post/surface.py::wall_triangle_adjacency`); per-node local
  Cartesian basis (Drela 2013 §III.B — in-plane rotation invariance absorbs the
  TE kink, no special TE equations). **Data-layout design points recorded at
  build time**: (1) reserve the wake-sheet unknowns (V4 continuation — same 6
  equations, wake closures); (2) a master-map hook so an IBL surface mesh built
  on the *uncut* wall can be fed from cut-mesh (LS) solutions; (3) single-group
  (`wall`) scope — the `wall`+`fuselage` seam is wing-body, out of V1 scope.
- `viscous/closures.py` — Drela 2013 wall + wake closure fits; laminar +
  turbulent + **forced transition** (free-transition e^N is a recorded follow-up,
  NOT gated in V1 — XFOIL comparisons must force transition at matched x_tr/c on
  both sides, or the gate compares transition models, not BLs).
- `viscous/ibl3.py` — 6-equation nonlinear surface Galerkin P1 FE: residual +
  analytic Jacobian, Numba kernels per design.md §7 (njit, no object mode,
  `PYFP3D_NOJIT=1` debuggable).
- `viscous/transpiration.py` — δ* → ṁ = ∇_Γ·(ρ_e u_e δ* ê) and the wall-RHS
  assembly (template: `solve/wall_correction.py::assemble_wall_flux_correction_rhs`);
  u_e extraction per-zone per A4 (LE band linear+smoothed, elsewhere quadratic).
- `viscous/coupling.py` — loose driver: FP solve → u_e → IBL3 solve → ṁ → RHS →
  FP re-solve, with under-relaxation on δ*.
- **Conforming-Newton external-RHS channel** in `solve/newton.py`
  (`R_free -= (Tᵀ b_ext)[free]` class of change): small but touches the core
  solver ⇒ bit-identity when absent + FD Jacobian check (project discipline).
  Picard: `solve_laplace` already has `body_source_rhs`, but the compressible
  Picard (`solve_subsonic` / `solve_subsonic_lifting`) does **not** — threading
  the wall RHS through it is a small V1 deliverable (needed by GV1.3/GV1.5);
  LS path uses the existing `b_base` slot.
- **Committed XFOIL reference**: generation script + CSV under
  `cases/reference_data/` (δ*/c, C_f, cl at matched M/Re/α/x_tr) — listed here
  because the 2026-07-20 review flagged it as a gate-blocking external artifact.
- **Body-of-revolution smoke mesh** (GV1.5): a standalone fuselage-alone family
  reusing `meshgen/fuselage.py` (`FuselageParams` + `add_fuselage_solid`, full-2π
  revolve) + an M0-style far field; surface tagged `wall` so the single-group
  surface machinery applies unchanged; wake-free, non-lifting (α = 0) — no
  Kutta/wake dependency anywhere in the case.

**Gates:**

- [ ] **GV1.1 standalone IBL3 verification** (prescribed u_e, no FP coupling):
  (a) laminar flat plate → Blasius: H within ±2 % of 2.59, δ*(x) ∝ √x;
  (b) turbulent flat plate: C_f(Re_θ) within ±5 % of the closure's own reference
  correlation; (c) prescribed decelerating u_e: separation indicator (H rise)
  at the self-similar reference location, band pre-registered; (d) quasi-2D
  invariant: crossflow unknowns (B, Ψ, C_τ2) ≈ 0 (structural lock); (e) surface
  refinement ×2: error drops, measured order recorded.
- [ ] **GV1.2 transpiration channel exactness**: (a) manufactured blowing on the
  M0 cylinder (Fourier-mode ṁ has an analytic exterior Laplace solution): φ
  error O(h) vs analytic; (b) ṁ = 0 **bit-identical** on ALL drivers (Picard,
  conforming Newton with channel absent, LS `b_base`) — the GV1.3-sketch
  "δ* = 0 bit-identical" clause lives here now; (c) conforming Newton Jacobian
  stays EXACT under lagged ṁ (FD check).
- [ ] **GV1.3 coupled 2.5-D NACA0012 subsonic vs committed XFOIL reference**
  (matched M/Re/α, forced transition, attached): δ*/c and C_f band OUTSIDE the
  LE/stagnation zone (band pre-registered at execution **and quoted alongside
  the A4 input band ≈2.5 % medium** — viscous-model error and inviscid-input
  error reported separately, never summed silently); viscous Δcl < 0
  (direction) with magnitude vs XFOIL's own viscous decrement (band
  pre-registered). LE-band pointwise comparison is RECORDED, not gated
  (input-limited 4–7 % per A4).
- [ ] **GV1.4 loose-loop convergence**: ‖Δδ*‖/‖δ*‖ < 1e-3 in ≤ 10 outer
  iterations on the GV1.3 case, under-relaxation factor recorded honestly; one
  transonic-attached 2.5-D point (M ~0.70–0.75) run and RECORDED (iteration
  count + relaxation), not pass/fail — the DN6-predicted near-separation
  divergence risk is measured here, and feeds the V2 skip decision.
- [ ] **GV1.5 fuselage smoke** (added 2026-07-22, user-directed) — the minimal
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
  recorded and masked, not chased. Headless artifacts per workflow rule 1.

**Prereq:** P6 ✓ + A4 ✓ (both done). V1 touches no wing-body wound and is
independent of the LS-side (b)-class work — parallelizable.

### V2 — Quasi-simultaneous coupling ☐ (optional; decision fed by GV1.4)

**Deliverable:** Hilbert-integral surface surrogate (`viscous/hilbert.py`);
the BLWF58 method document is the reference description of the approach.

**Gates:**

- [ ] GV2.1 ≥ 30 % fewer coupling iterations than the V1 loose loop on the same
  ladder, OR converges a case the loose loop cannot (near-separation robustness).
- **Skip criterion (concrete):** if GV1.4 passes at ≤ 10 iterations including
  the recorded transonic point without per-case tuning, V2 is SKIPPED — record
  the decision in the ledger and move to V3.

**Prereq:** V1.

### V3 — Tight coupling: augmented Newton ☐

**Deliverable:** augmented (φ, Γ, BL) Newton on the P8/P14 machinery; coupling
blocks J_φ,BL (∂ṁ/∂BL through the transpiration assembly) and J_BL,φ (∂u_e/∂φ
through the recovery operator chain); GMRES + block preconditioning (AMG on the
φ block / ILU on the BL block; `solve/schur_ls.py` is the structural prototype —
**note** the BL block is O(6 × wall nodes), far bigger than the LS aux thin
band, so exact Schur elimination may not pay: measure, don't assume).

**Gates:**

- [ ] **GV3.0 M6 subsonic loose-coupling bridge** (RECORDED, entry check; added
  2026-07-22, user-directed) — runs on the **V1 loose driver** (no augmented
  Newton), scheduled here so the 2.5-D → transonic-3-D jump is bridged and the
  crossflow content (Ψ, B equations) gets its **first live 3-D exercise**
  before GV3.3: ONERA M6 (existing `cases/meshes/onera_m6/` family, coarse +
  medium), conforming path, M0.5 / α 3.06° (the committed subsonic convention),
  forced transition, Re recorded (e.g. the M6 experiment 11.72e6). Outputs
  RECORDED, not pass/fail (no δ*(z) truth data exists for M6): δ*(z) spanwise
  distribution at fixed x/c stations (committed CSV + PNG), crossflow-magnitude
  field, ΔCL viscous−inviscid (expected DOWN, direction recorded), 3-D
  loose-loop iteration count vs the GV1.4 2.5-D count; wing-tip band masked
  (tip_taper r_c = 0.05·b_semi). Feeds GV3.3's band pre-registration.
- [ ] **GV3.1 exactness + convergence**: both coupling blocks FD-verified
  (project Jacobian discipline; the B19/B31 FD-gate pattern) + quadratic tail on
  the GV1.3 case; outer iterations ≤ half the V1 loose loop.
- [ ] **GV3.2 2-D transonic VII vs experiment**: RAE2822
  (`cases/reference_data/rae2822_experiment/`, M0.725/α2.55 + M0.73/α3.19,
  Re 6.5e6): shock location within a pre-registered band of experiment + Cp RMS
  recorded. **Needs a NEW 2.5-D RAE2822 mesh family** (small meshgen task, the
  2.5-D extrusion machinery exists) **+ the A4 TE-wedge pre-check**: if the RAE
  TE wedge < ~6° the quadratic-recovery guard fires ⇒ TE-band u_e falls back to
  linear+smoothed, recorded in the gate.
- [ ] **GV3.3 M6 wing direction+magnitude check**: CL moves **down** from the
  converged inviscid baseline — **cl_KJ 0.2823 (medium, P14 pressure-Kutta) /
  0.2866 (fine, P13 tapered)** — toward experiment ≈ 0.26–0.27. The expected
  move ≈ 0.02 (≈ 7 %) is safely above the A4 input floor (2.5 %), so this gate
  is attributable. (The pre-P14 "0.245 vs 0.288" framing is superseded — see
  scope guards.)
- [ ] **GV3.4 cost (RECORDED)**: augmented step wall-time ≤ ~2× the inviscid
  Newton step on M6 medium with the block preconditioner working; measured
  number recorded either way.

**Prereq:** P8 ✓ + P14 ✓ + V1. **Wing-body VII is explicitly OUT of V3 scope**
(scope guards below).

### V4 — Wake-sheet IBL correction ☐ (continuation of V1's data layout, not an independent solver)

Same 6 equations with wake closure relations; the wake unknowns were reserved in
V1's layout. δ*_wake enters as the wake-sheet RHS mass source; TE thickness
continuity δ_wake(TE) = δ*_upper + δ*_lower.

**Gates:**

- [ ] **GV4.0 design adjudication (BEFORE code, user-adjudicated)**: the LS-path
  sheet-source mechanism. `pyfp3d/wake/` has **no sheet-surface integration
  machinery** — the options are zero-isosurface polygon integration (new
  geometry code) or a volume-band approximation (deviates from the "sheet RHS"
  formulation). The conforming path needs NO new mechanism (explicit
  `wake_minus/plus` faces + slave→master folding IS the weak-form flux channel).
  V4 may close conforming-only; the LS leg then becomes a recorded follow-up.
- [ ] **GV4.1 conforming sheet source**: δ*_wake enters via the
  `constraints/wake.py` reduce RHS (Tᵀ b_wake); TE thickness continuity
  asserted; δ*_wake = 0 **bit-identical**.
- [ ] **GV4.2 measured effect**: wake-IBL on/off cl (and TE-region Cp) delta on
  the GV1.3 case, direction-checked against XFOIL's wake modelling; RECORDED
  with the A4 input band quoted.

**Design constraints (unchanged from DN2 §4.5):** TE kink absorbed by Drela
local-basis adaptation; **straight wake + mass-transpiration relaxation, no
geometric relaxation**.

**Prereq:** V1 (+ GV4.0 adjudication for any LS leg).

### Scope guards (re-based 2026-07-22; the DN2 §9 / DN6 §13–14 envelope stands)

- **Validity envelope**: attached / mildly-shocked flow (M_shock ≲ 1.3); not
  massive or shock-induced separation. The M6 M0.84 shock sits at the envelope
  edge; the wing-TIP singularity zone (local M_max 2.5+) is OUTSIDE it — IBL
  seeding/comparison bands must mask the tip band (on the conforming path the
  production `tip_taper` r_c = 0.05·b_semi band is the natural mask).
- **VII does not close the inviscid-discretization CL gap — updated numbers**:
  after P13/P14 the remaining inviscid gap to the FP reference is ≈ 0.5 %
  (cl_KJ 0.2823 medium / 0.2866 fine vs 0.288), so the ~0.02 delta to
  experiment (≈ 0.26–0.27) is now genuinely viscous-dominated; viscosity moves
  CL **down**. (History, superseded: the pre-P14 "0.245 vs 0.288 → sharp-TE/LE
  P1 floor → P9/P11" attribution — P11 measured NEGATIVE 2026-07-19, and
  P14/G14.7 showed 69 % of the old 0.019 gap was Kutta-estimator bias. See the
  2026-07-20 review §3.)
- **Input-error discipline (A4, standing rule for every V-gate)**: every
  viscous-vs-reference comparison quotes the A4 inviscid u_e input band (medium
  ≈ 2.5 % peak-relative; LE/stagnation 4–7 %) alongside the viscous
  discrepancy. Tight LE-band comparisons are input-limited by construction.
  Only the **unchosen** G1.6 route (b) — isoparametric P2 wall layer — could
  raise the input band to O(h²): if a V-gate fails for input reasons alone,
  that is the recorded escalation route, NOT a viscous-model fix (G1.6 Option C
  close-out note, 2026-07-22).
- **Wing-body VII (applying V3/V4 to the M2 wing-body) is DEFERRED** until the
  LS-side wing-tip (b)-class is cured or explicitly accepted (B30 attribution;
  B31 cured the conforming side via tip_taper, LS C-class closed negative
  C1/C3) — otherwise wing-body viscous gates re-enter "viscous model vs known
  inviscid wound" attribution confusion. The M6 **wing** gate GV3.3 is
  unblocked (conforming, tip_taper production since B32).
- **Reynolds number and transition are new physical inputs** (the FP solver has
  neither): Re enters only the closures; V1 ships forced transition, free
  transition (e^N) is a recorded follow-up. Gate comparisons must match Re and
  x_tr explicitly.
- V1 is parallelizable with the remaining Track-B/LS work (it depends only on
  P6 + A4 and touches no wing-body wound), but it is a large, self-contained
  solver effort (6-equation nonlinear surface FE + closures): **budget it like
  a Track-P phase, not a side task**.

---


## Progress ledger

### Track V — viscous–inviscid interaction

Track status: **◐ IN PROGRESS — V1 OPENED 2026-07-22** (gates re-spec'd at
opening against the shipped B32/A4 state; the pre-2026-07-22 GV1.1–GV1.3
sketch is superseded — Track V had zero implementation, so no historical gate
result is affected). Design 2026-07-09/10 (DN2 + DN6, historical; two known
DN6 traps annotated in the header re-spec block). Validity envelope: attached /
mildly-shocked flow. **VII does not close the inviscid-discretization CL gap**
— the inviscid baseline is now clean to ≈ 0.5 % (P13/P14), so the ~0.02 delta
to experiment is the viscous target; direction is DOWN.

- V1 — **◐ OPENED 2026-07-22** — IBL3 solver core + loose coupling, 2.5-D ladder + fuselage smoke
  (`viscous/surface_mesh.py`, `closures.py`, `ibl3.py`, `transpiration.py`, `coupling.py`;
  conforming-Newton external-RHS channel + compressible-Picard RHS threading; committed XFOIL reference;
  body-of-revolution smoke mesh). Gates GV1.1 (standalone IBL3 vs analytic/self-similar), GV1.2
  (transpiration exactness + ṁ=0 bit-identity on all three drivers + FD), GV1.3 (coupled NACA0012 vs
  XFOIL, forced transition, A4 band quoted), GV1.4 (loose loop ≤ 10 iters; transonic point RECORDED →
  V2 skip decision), GV1.5 (fuselage body-of-revolution smoke — axisymmetry + crossflow≈0 asserted,
  rest RECORDED; Track V's only fuselage-alone item, no junction/wake). Prereqs P6 ✓ + A4 ✓;
  no wing-body contact.
- V2 — ☐ (optional) — quasi-simultaneous coupling (`viscous/hilbert.py`, BLWF58 reference). Gate GV2.1;
  concrete skip criterion wired to GV1.4.
- V3 — ☐ — tight coupling: augmented (φ, Γ, BL) Newton on P8/P14. Entry check GV3.0 = M6 subsonic M0.5
  loose-coupling bridge (RECORDED; V1 driver, first live 3-D crossflow exercise, δ*(z) + ΔCL direction +
  3-D iteration count; bridges 2.5-D → transonic 3-D, feeds GV3.3's bands); FD-verified coupling blocks
  (GV3.1), RAE2822 transonic VII vs committed experiment (GV3.2; needs the 2.5-D RAE2822 mesh family +
  A4 TE-wedge pre-check), M6 CL-down direction+magnitude vs inviscid 0.2823/0.2866 (GV3.3), cost
  recorded (GV3.4).
  Block precond: AMG-φ / ILU-BL, `schur_ls.py` prototype; BL block is NOT thin — measure before Schur.
  Wing-body VII out of scope (deferred, see scope guards).
- V4 — ☐ — wake-sheet IBL correction, a continuation of V1's data layout (wake unknowns reserved).
  GV4.0 LS sheet-source design adjudication BEFORE code (conforming needs no new mechanism; may close
  conforming-only), GV4.1 conforming sheet source + δ*_wake=0 bit-identity, GV4.2 measured on/off effect.
  Straight wake + mass-transpiration relaxation, no geometric relaxation.
