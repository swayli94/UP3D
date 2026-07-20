# pyFP3D Roadmap — Track B (level-set embedded wake B1–B27)

> Split verbatim from `docs/roadmap.md` on 2026-07-15 (content unchanged; only
> this header and the ledger heading were added). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md).

## Track B — Level-set embedded wake

- (designed 2026-07-07; IN PROGRESS — B1 ✓ B2 ✓ B3 ✓ B4 ✓ B5 ✓ B7 ✓; B6 ◐ since 2026-07-12:
  coarse gate met + LS Newton delivered, medium closure open; **B8 ✓ CLOSED 2026-07-14 characterized-not-cured (user-arbitrated;
  both constraint-side cures measured negative; B9 unblocked)**; **B11 ✓ CLOSED 2026-07-14 — LS-path infrastructure:
  unified post-processing + GMRES/AMG scaling (the deferred §5.3 escape from the splu wall), NEW appended after B10**;
  **B12 ✓ CLOSED 2026-07-14 — lagged-LU direct-reuse for LS Newton (M6 medium Newton 2.18× via 1 factorization vs 7), NEW appended after
  B11**; **B13 ✓ CLOSED 2026-07-14 — lagged-LU on the Picard outer loop (M6 medium lifting 6.55× 447.6→68.3 s, end-to-end ~3× ~330→112 s)**;
  **B14 ✓ CLOSED 2026-07-17 — Schur-eliminated-aux + AMG structural preconditioner (`precond="schur"`);
  the A1 precond bottleneck GONE (M6 medium M0.84 43.6% → 2.6%, ramp 1.43× / subsonic 2.08×), γ = the committed GB15.4;
  fine-scale route remains the unbuilt designed use-case**; **B15 ✓ CLOSED 2026-07-15 — LS Newton transonic ramp + N5 freeze-selection:
  the Picard shock plateau is GONE (M6 medium M0.84 2304.7 s bounded-stall → 657 s all-levels-converged, 3.5×;
  NACA coarse M0.80 5.6× and strict), plus FOUR errata proving the conforming N5 recipe is not mechanically portable**;
  **B9 ✓ CLOSED 2026-07-17 (RE-SPEC'D, user-approved) — wing-body cross-model validation: LS (Picard) + conforming (NEW capability, Newton)
  AGREE to 0.4%/0.6% at medium M0.5; GB9.4 fuselage-lift XFAIL ⇒ G1.6 fuselage-Cp error; LS Newton diverges on the wing-body = the `neumann`
  far-field blockage, not the solver**; **B16 ✓ CLOSED 2026-07-18 (churn fix; lift-convergence OPEN) —
  LS Newton far-field BC generalisation: the wing-body churn is a near-singular far-field aux block (cond1 O(1e19), the 8 rows |R|≈84
  reproduced), `farfield_aux="pin"` (default) drops it to 8.7e6 ⇒ freestream Newton reaches res 5.88e-14 / 0 limited at COARSE (lift matches
  conforming 0.1%) where legacy churns at 7.95; neumann byte-identical. ★★ BUT GB16.4 XFAIL —
  the MEDIUM Newton-pin STALLS at res 7e-6, lift 22% below Picard/conforming (which agree, B9's 0.4%);
  the {Newton,Picard,conforming} triangle does not close ⇒ UNRESOLVED non-convergence, open follow-up (**RESOLVED by B17**)**;
  **B17 ✓ CLOSED 2026-07-18 (NEW, user-directed; resolves GB16.4) — the far-field pin must carry jump=γ, not 0:
  B16's freestream jump=0 REMOVED the outflow wake circulation (a BC-modelling error, NOT a non-convergence —
  an independent Picard-pin converges to the same medium 0.169 the Newton-pin "stalls" at);
  `farfield_aux="pin_gamma"` (new default, both solvers) closes the triangle MONOTONE to conforming (coarse 0.2087, medium 0.2117
  Picard/0.2114 Newton — 0.2115 pre-B20; the medium Newton trajectory now converges to |R|~1e-13 post-B20); vortex brackets from +2.5%;
  B9 coarse-12.8% erratum'd as far-field contamination**; **B18 ✓ CLOSED 2026-07-18 (NEW, user-directed;
  executes the GB16.6 debt) — wing-body transonic (M0.84): CONFORMING reaches it (coarse M0.84 0.2617, medium M0.79 0.2579, clean cl(M)
  rise), LEVEL-SET is junction-limited (post-B20 re-baseline: coarse ceiling ~M0.55/Mmax 1.31, medium dies ~M0.5 with a GENUINE unclamped
  Mmax 5.22 — the junction pocket WORSENS with refinement, closed-negative; ★ attribution CORRECTED by GB20.5:
  the pocket is the G1.6/GB9.4 faceted-geometry error, NOT mixed-plain contamination); no common transonic Mach at medium so the trustworthy
  cross-model stays M0.5 (2.6%), plus a post-B20 coarse M0.6 transonic point (0.2178 vs 0.2174, 0.2%); GB18.1 PASS + GB18.2–5 RECORDED;
  no pyfp3d/ change**; **B19 ✓ CLOSED 2026-07-18 (NEW, user-directed; executes A3/C1) — the LS-Newton Jacobian is now EXACT in 3-D (two
  defects: DOF column maps + gradient factors; probe 1.146e-01 → 1.33e-08, ε-discriminator flipped);
  R bit-identical, NO convergence gain (GB19.4 negative — the plateau is B15 selection churn);
  Leg B measured the mixed-plain side-density contamination (spurious supersonic q² 3.22 vs 1.34) and routed it to B20**;
  **B20 ✓ CLOSED 2026-07-18, ADOPTED PERMANENTLY (user-arbitrated; knob removed), re-baselined 2026-07-19 —
  mixed-side plain elements read their density from the MAIN field; every moved 3-D number went B20's way except ONE:
  the M6 medium ramp ceiling fell M0.84 → ≈M0.675 (GB20.7 called it "a REAL capability loss" — **overturned by B21 the next day**:
  the loss was B20's own patch gap in `freeze_side_state`)**; **B21 ✓ CLOSED 2026-07-19 (NEW;
  executes the Kimi-inspection N1 finding) — `freeze_side_state` missed the B20 patch (froze selections the live system would not make, 3-D
  only; 83+9 probe-confirmed); the one-line fix RESTORES the M6-medium M0.84 ramp (γ 0.088343, res 9e-14, 515 s) ⇒ GB20.7's "real capability
  loss" verdict OVERTURNED — the loss was B20's own patch gap; 3-D capture-consistency test lock added (verified failing pre-fix)**;
  **B22 ✓ CLOSED 2026-07-19 (NEW; executes B21's follow-up + Kimi N3/§2/§5) — B15 demo refreshed 20/20 + B14 7/7 on the B21 state (coarse
  ramp γ 0.084931 disclosed); N3 CLOSED via gated absolute anchor locks (`test_b22_ls_3d_anchors.py`);
  re-baseline erratum checklist now process (CLAUDE.md step 5 / discipline #11); next-phase priorities analysis recommends P11 (user's
  call)**; **B23 ✓ CLOSED 2026-07-19 (NEW) — wing-body junction discriminator: the spurious supersonic pocket is lift/wake-coupled not
  geometric (α=0 clean at both levels, pocket grows superlinearly with α), attributed to the wake inboard free-edge singularity;
  P11 close-out input delivered**; **B24 ✓ CLOSED 2026-07-19 (NEW; negative) — waterline-extension route CLOSED:
  the pocket follows the free edge (hypothesis re-confirmed) but both extension variants (B1 flush / B3 offset-cone) trade the singularity
  for equal-or-worse forms; decision-tree exit 3 ⇒ back to B23 (b)-2 (P13 free-edge treatment rescoped to the inboard end)**;
  **B25 ✓ CLOSED 2026-07-19 (NEW; cure C-A) — inboard fragment clip gives the free edge a legal home (= conforming fragment topology):
  medium α=3.06 junction pocket corrM 14.66 → 0.63 with all physical guardrails clean; `inboard_clip` in `wake/cut_elements.py` +
  `meshgen/fuselage.py:make_inboard_clip`, default None bit-identical**; **B26 ✓ CLOSED 2026-07-20 (NEW; B26-A ceiling-lift) —
  post-cure LS transonic ceiling re-measured: medium 0.50 → 0.7625, coarse 0.82 → 0.84 reached;
  both death causes flip (a)-pocket-rejection → (b)-high-M Newton stall with the peak at the WING TIP (P13-class, not the pocket)**;
  **B27 ✓ CLOSED 2026-07-20 (NEW) — B18 demo refreshed on the B25/B26 state: GB27.1/27.2/27.3(0.65) PASS, 27.3(0.75)/27.4/27.5 RECORDED as
  pre-registered; checks.csv 8/8 PASS, 336/336 bit-identical; B18 façade re-worded from "LS junction-limited (closed-negative)" to
  "post-cure LS ceiling co-located with conforming"** — B10 = curved wake shelved)

> **★ Track-B renumber 2026-07-12 (user-directed).** TWO renumbers landed the
> same day. **(1)** A new **B4 — TE control-volume / implicit-Kutta
> re-derivation** was INSERTED (B3's emergent circulation converges to the wrong
> value; design_track_b.md §9). **(2)** The half-integer IDs were then
> regularized away — the far-field A/B and the M6 3D gate become full phases and
> everything after shifts up by one. **Net mapping from the pre-2026-07-12
> scheme: old B4 (transonic) → B6; old B4.5 (M6 3D) → B7; old B5 (multi-wake) →
> B8; old B6 (curved wake, shelved) → B9; and the far-field A/B (once B3.5, then
> B4.5) → B5.** Docs written before 2026-07-12 use the old IDs; docs written
> between the two same-day renumbers use the interim B4.5/B5/B5.5/B6/B7 IDs.
>
> **★ Track-B renumber 2026-07-13 (user-directed).** A new **B8 — level-set
> tip-edge desingularization (row-blend tip taper)** was INSERTED (the LS
> analogue of P13/G13.2's conforming taper; G13.2 finding (8)). Everything after
> shifts up by one: old **B8 (multi-wake) → B9**; old **B9 (curved wake,
> shelved) → B10**. Docs written before 2026-07-13 use the pre-insertion IDs
> (B8 = multi-wake, B9 = curved wake).

Deliverable: `wake/` — a level-set wake representation + multivalued (CutFEM-style)
elements + implicit Kutta (TE duplication + wake least-squares condition; penalty
Kutta demoted to an optional diagnostic — design_track_b.md D2), replacing the
conforming embedded wake sheet + master–slave Γ elimination.
**Purpose (user-arbitrated 2026-07-11): mesh/geometry workflow capability** —
no pre-embedded wake surface, α sweeps without remeshing, blunt-TE anchoring,
multi-wake / wake–fuselage intersections for M2, structural elimination of the
st133-class Kutta-probe failures. NOT solver speed: the original
kill-the-Γ-secant efficiency motivation is obsolete post-P8 Newton
(design_track_b.md §1); efficiency criteria below are non-regression guards only.
Design record: DN1 (historical; `discussion_notes/` deleted 2026-07-14 —
`git show 8aa4aee:docs/discussion_notes/20260707_1505_levelset_wake_design.md`),
superseded by design_track_b.md. **Status (2026-07-12): B1 + B2 + B3 + B4 CLOSED.**
The level-set path now produces LIFT with an implicit Kutta (no Γ secant, no
master–slave Γ): Γ matches the conforming solver **within 1% on the same mesh**
at M=0 and M=0.5, cl lands inside the committed [PG, KT] bracket, and the
**wake-free** M3 mesh (no `wake` tag at all) reproduces the embedded-mesh
circulation to 0.3% — the workflow payoff. **★ The B4 finding (design_track_b.md
§9):** the wake LS is STRUCTURALLY blind to a constant jump (Σ_c ∇N_c = ∇(1) = 0
⇒ residual identically zero, measured 1.9e-16), so "g₂ IS the discrete Kutta" is
FALSE — Γ needs its own condition. B4 supplies it: the **nonlinear TE
pressure-equality (Bernoulli) Kutta** |q_u|² = |q_l|², factorized exactly as
(q_u+q_l)·(q_u−q_l) = 0 and linearized by freezing the mean, with q recovered on
the TE's **WALL-ADJACENT** upper/lower control volumes. **B5 closed
2026-07-12 — far-field verdict: option a (Dirichlet+vortex) stays the default**
(domain-robust to <1%; the López Neumann outlet truncates O(Γ/R) and needs a
2–4× larger domain; M6 leg folded into B7). Next = B6 (transonic) / B7
(M6 3D).
**Numerics reference (2026-07-11):** [design_track_b.md](../design_track_b.md) —
theory/implementation analysis cross-checked against the López dissertation;
supersedes DN1 as the Track B numerics spec (key deltas: 3D wake BC uses the
López g₁+g₂ two-component LS form, NOT the Núñez full-vector form; Kutta is
implicit via TE duplication + the wake LS condition, penalty demoted to
optional; on-wake nodes need the ε side-shift DN1 missed). Gate re-specs were
**user-arbitrated 2026-07-11 and merged into the phase entries below** (B3 re-spec;
B5 far-field A/B NEW; B4 re-anchored post-P4-erratum, medium at M0.7875;
B7 M6 3D gate NEW); design_track_b.md §7 is the arbitration record.


### B1 — Level-set wake + cut-element identification ✓ (closed 2026-07-11)
**Deliverable:** Level-set wake + cut-element identification (`wake/levelset.py`, `wake/cut_elements.py`)
**Gate (dual-mesh: M0/M1 wake-embedded + M3/M4 wake-free):**
- [x] **CLOSED 2026-07-11** (`tests/test_b1_cut_elements.py`, 34 passed, 2.5D coarse+medium of both families + 3D M6 of both families):
      (a) M0 embedded — every conforming sheet node ε-shifted "+" (D4 stress test at scale), census **exactly** == `cut_wake`'s minus-side
      element star (`cut_elems ∪ te_lower_elems`, element-by-element), TE nodes == `wc.te_nodes`; (b) M3 wake-free —
      generic cuts, gap-free corridor TE→far field at α=0 AND after `update_direction` to α=4° **on the same mesh**;
      (c) **M1/M4 ONERA M6 (3D)** — swept TE polyline: census is a strict **superset** of the conforming minus-star (0 missing, +2.9%
      extras, all tip-edge straddlers — the sheet's tip edge conforms to element edges on the M1 mesh but passes THROUGH elements for the
      level set; expected and explained), spanwise clip verified (nothing cut wholly outboard of the tip); M4 wake-free —
      same census structure with no `wake` tag, spanwise-gap-free sheet, α re-aim 0°→3.06° on the same mesh. Delivered:
      TE-**polyline** ruled level set (D9) with an **oblique (v, d̂, n̂) frame** — ★ a swept TE is not perpendicular to the wake direction,
      so an orthogonal span projection leaks the downstream distance into the spanwise coordinate and wrongly clips ~60% of the true M6 cut
      set (measured, then fixed; regression-pinned); **spanwise clip** (crossings must satisfy 0 ≤ q ≤ span_length —
      outboard of the tip the sheet has ended, Γ(tip)=0; without it the level set re-creates P5's far-field branch-ray artifact);
      downstream-crossing test (excludes the ahead-of-LE sign-change region); `te_lower_elems` recorded for B2's López-fig-3.6c aux
      assignment

### B2 — Multivalued FE assembly ✓ (closed 2026-07-11)
**Deliverable:** Multivalued FE assembly (`wake/multivalued.py::MultivaluedOperator`, `kernels/cut_assembly.py`; `solve/picard_ls.py::solve_multivalued_laplace` non-lifting driver — parallel to the conforming path, which stays byte-untouched)
**Gates (dual-mesh):**
- [x] **CLOSED 2026-07-11** (`tests/test_b2_multivalued.py`, 17 passed on coarse+medium of both 2.5D families and 3D M6 coarse of both
      families; some medium/M6 parametrizations skip in CI where the meshes are gitignored). Key design:
      a cut element is the SAME P1 element matrix assembled twice with `dofs_upper`/`dofs_lower`, expressed as a sparse **redirection** of
      the single-valued matrix — on a cut element the entries whose two nodes are on OPPOSITE sides move their column main(b)→aux(b)
      (`multivalued_redirection_coo`); everything else is byte-identical to `PicardOperator.assemble_matrix()`.
      Aux rows carry the B2 **continuity ("weld") closure** aux_k = main_j (`continuity_closure_coo`), which makes the extended (n_total =
      n_main + n_ext) system reduce EXACTLY to the single-valued one — proven directly (`test_extended_matrix_folds_to_stiffness`:
      folding aux→main recovers the stiffness matrix to 1e-13). The extended matrix is structurally nonsymmetric (weld rows), solved by
      sparse-direct LU (`spsolve`); GMRES+AMG is the B3+ scaling path (design_track_b.md §5.3).
      B3 replaces the weld block with the g₁+g₂ wake LS (implicit Kutta), at which point [φ] becomes nonzero.
- [x] V0 freestream (φ = U·x, full Dirichlet) < 1e−12 on the cut mesh: 2.5D M0/M3 α=0 and α=4° = **0.0** (exact linear field); 3D M6 M1/M4 = **1.1e−14 / 3.4e−14**
- [x] V1 MMS slope ≥ 1.9: cube cut in generic position (8° tilted half-plane), 3-level slope **1.94**
- [x] Laplace α = 0 gives cl ≈ 0: TE jump = 0 (the weld forbids a jump) ⇒ cl_KJ = 0, and the main potential matches the single-valued `solve_laplace` oracle to **~3e−11** — on both mesh types (dual-mesh rule)

### B3 — Lifting solve with implicit Kutta ✓ (closed 2026-07-12)
**Deliverable:** Lifting solve with **implicit Kutta** — no Γ secant, no master–slave Γ constraint:
the TE jump is carried by the multivalued aux DOFs, the g₁+g₂ wake LS convects it downstream, and its VALUE is set by B4's nonlinear TE
pressure-equality Kutta. Γ is a RESULT. `kernels/cut_assembly.py` (`mass_conservation_coo` with per-side ρ per D10, `wake_ls_coo`,
`te_kutta_coo`), `wake/multivalued.py` (`closure="wake_ls"`, side potentials/densities, TE control volumes),
`solve/picard_ls.py::solve_multivalued_lifting`. Far field = Dirichlet freestream + vortex on the **MAIN** DOFs, aux **FREE**.
**Gates (dual-mesh) — `tests/test_b3_lifting.py`, 6 passed:**
- [x] **V3 M0.5 α=2°: cl inside the committed [PG, KT] bracket** — cl_KJ **0.2828** (medium) inside [PG 0.2788, KT 0.2919], read from `cases/reference_data/naca0012_m05/cl_reference.csv` (the same file the conforming G3.2 gate reads). Holds on BOTH the wake-embedded M0 family and the **wake-free M3** family.
- [x] **Same-mesh A/B vs conforming** (only possible on the embedded family): Γ within **0.1–0.7%** at M=0 and M=0.5, coarse/medium/fine (0.1177/0.1191/0.1197 vs 0.1175/0.1200/0.1202). *The old "within 1% same-mesh" clause — retired on 2026-07-12 as unmeetable — is in fact MET, now that B4 supplies a real Kutta.*
- [x] **Wake-free M3 mesh** (no `wake` tag, generic cuts — the actual workflow form): Γ within **0.3%** of the embedded-mesh result; medium cl inside [PG, KT]. *(M3 coarse cl 0.2773 sits 0.5% under the PG edge — a coarse-mesh accuracy artifact; the gate lives on medium per hard rule 2.)*
- [x] Γ emerges (no outer Γ loop) and the wake jump is CONVECTED, not decaying (pinned: far-field aux DOFs must stay free)

### B4 — TE control-volume / implicit-Kutta re-derivation ✓ (NEW + closed 2026-07-12, user-directed)
**The B3 blocker and its fix (design_track_b.md §9).** Two structural facts:
**(1) The wake LS CANNOT pin Γ.** Its residual is identically zero for any spatially-constant jump, because Σ_c ∇N_c = ∇(Σ_c N_c) = ∇(1) = 0
(partition of unity) — measured **1.9e-16**. ⇒ design_track_b.md §2.3/D2's "the g₂ on the TE-adjacent wake element IS the discrete Kutta
condition" is **FALSE and retired**; the López dissertation has no explicit Kutta anywhere ("Kutta" never appears in its method chapter).
**(2) Γ was therefore pinned by a single, WRONG equation** — the TE aux row (lower-side mass conservation), whose control volume is up/down asymmetric on a symmetric airfoil (TE fan 9 upper / 6 lower / 3 cut, because the ε shift sends every on-sheet node "+"). It over-circulated **+42%** (Γ → 0.168 vs 0.120), mesh-convergently.
**★ Deliverable — the nonlinear TE pressure-equality (Bernoulli) Kutta.** Symmetrizing the control volume is NOT available (the mesh is
naturally asymmetric at the TE — user-arbitrated), so the condition is instead a POINTWISE PHYSICAL statement needing no symmetry:
|q_u|² = |q_l|², factorized **exactly** as (q_u+q_l)·(q_u−q_l) = 0 and linearized by freezing the mean s̄ = q_u+q_l at the previous iterate
— a row that is LINEAR in φ, re-linearized each Picard outer (same cadence as the density lag, no new outer loop), converging to the exact
nonlinear condition. It sits on the TE **aux** row; the displaced lower-side mass-conservation entries are re-routed onto the TE **main**
row, which then carries the TOTAL (upper+lower) balance — so mass stays conserved and no side is arbitrarily robbed of its equation.
**Why it is non-degenerate where g₂ is not:** q_u and q_l are recovered on DIFFERENT element sets, so q_u−q_l is not a jump gradient and
does not vanish for a constant jump.
**★ The control volumes must be WALL-ADJACENT** (elements carrying a wall face — the upper/lower body surface at the TE), not the whole
element fan: the Kutta condition is about the SURFACE velocities. Measured (Γ, coarse/medium/fine vs conforming 0.1175/0.1200/0.1202):
full-fan recovery **0.1407/0.1355/0.1329** (+11–15%, interior and wake elements pollute the average) → wall-adjacent
**0.1177/0.1191/0.1197** (**<1%**).
**Interfaces:** `solve_multivalued_lifting(..., te_kutta="pressure")` (default); `te_kutta="mass"` keeps the old B3 row for the before/after
contrast. `kernels/cut_assembly.py::te_kutta_coo`, `wake/multivalued.py::{_build_te_control_volumes, te_velocities}`.
**The D2 penalty-Kutta fallback is no longer needed** — this route has no penalty weight and no tuning parameter (s̄ is solved for, not
calibrated).
**Gates — `tests/test_b4_te_control_volume.py`, 8 passed (~29 s):**
- [x] the wake-LS constant-jump null space is pinned numerically (1.9e-16) — the reason a separate TE condition is structurally required
- [x] the TE control volumes are wall-adjacent (every element carries a wall face)
- [x] the below-TE fan is never cut (the López p.57 ε-shift trap; regression pin)
- [x] **NACA0012 α=2° incompressible: emergent Γ within 5% of conforming** — measured **0.1–0.7%** on coarse and medium; and the old `te_kutta="mass"` row is still >30% out (before/after contrast held honest)
- [x] visual artifact `artifacts/EXPORT_TE_DIAGNOSIS/b4_te_kutta.png` + `summary.csv` (wall-adjacent control volumes + the three Γ-vs-h curves against conforming)

### B5 — Far-field A/B: Dirichlet+vortex vs Neumann outlet ✓ (was B4.5, orig B3.5; NEW 2026-07-11; closed 2026-07-12, user-arbitrated)
**Deliverable:** option a (spherical Dirichlet + vortex on the main DOFs) vs option b (López-style Neumann outlet, no vortex; domain re-calibrated per the dissertation §4.1.4 — note López uses **10²–10⁷ chord** domains vs pyFP3D's 15c) — design_track_b.md §5.4/D7. **Interface:** `solve_multivalued_lifting(farfield="vortex"|"neumann"|"freestream")`; default `"vortex"`.
**Note (2026-07-12):** far-field truncation is **NOT** the B3 over-circulation cause — imposing the conforming true solution's far-field trace as Dirichlet still yields Γ = 0.1721. Confirmed after B4 closed: with B4's correct Kutta, the far-field question decouples cleanly and is purely an O(Γ/R) truncation study.
**Gate:**
- [x] **CLOSED 2026-07-12 — verdict: option a (Dirichlet+vortex) STAYS the default.** López-style domain-size re-calibration on the NACA
      dual-mesh families (M0 embedded + M3 wake-free), coarse, M0.5 α2°, far-field radius R ∈ {15,30,60,120}c (demo
      `cases/demo/b4p5_farfield/`, `tests/test_b45_farfield.py` 10 passed). **Measured:** option a is **domain-robust** —
      Γ within **0.45%** (M0) / **1.09%** (M3) of the truth across 15→120c, and **0.25%** of the conforming solver at 15c;
      option b truncates the O(Γ/R) point-vortex tail — Γ **−4.07%** below a at 15c, halving each doubling of R (−2.0% at 30c, −0.99% at
      60c, −0.50% at 120c), so it meets the B3 ±2% band only at **R ≥ ~30c** and <1% at **R ≥ 60c** (a 2–4× larger domain, 4× the tets at
      equal near-body h — consistent with why López needs 10²–10⁷c domains); freestream-Dirichlet is the crudest at every R (and DIVERGES on
      the compact 15c M0). Both mesh families give bit-for-bit the same story. Since the O(Γ/R) truncation is geometry-universal (a 3D wing
      truncates the same horseshoe tail), this decides the far-field default for the M6 B-path too.
      **M6 leg folded into B7** (user-arbitrated 2026-07-12): the level-set B-path *solve* on M6 needs the 3D wake-BC machinery that is B7's
      deliverable — and, separately measured here, the span-uniform option-a vortex without the P5 Γ(z) taper recreates the branch-ray
      artifact on M6, itself B7 machinery.

### B6 — Transonic + Mach continuation on the level-set path ◐ IN PROGRESS 2026-07-12 (coarse gate ✓; medium = bounded Picard state, isolated fold solution deferred to LS Newton) (was B5, orig B4)
**Deliverable:** Transonic + Mach continuation on the level-set path (~~inherits `damping_theta`~~ — **measured 2026-07-12:
the inherited stabilizer does NOT transplant as-is**, see findings below). Delivered so far (design_track_b.md §10):
per-side artificial density on the cut elements with a same-side-restricted upstream walk (`MultivaluedOperator.element_rho_tilde`, D10 —
subcritically an exact no-op; the M0.80 blow-up cells sit in the pocket ABOVE the airfoil, zero on the wake strip, so the shock machinery is
isomorphic to conforming); **supersonic-zone-LOCALIZED damping** (`damping_scope="supersonic"`) —
the P4 whole-field θ·diag form is a Jacobi smoother that throttles the smooth global circulation mode, which on the B path is a SOLUTION
mode (conforming keeps Γ outside the damped matrix as a secant unknown; measured: Γ crawls 0.0005→0.017 in 160 outers vs undamped
convergence in 35); `solve_multivalued_transonic` Mach ramp with **no Γ secant** (a level = one warm-started Picard solve —
the st133-class per-station closure failure is structurally impossible); `post/surface_ls.py` (D11 wall Cp + shock extraction on the B
path); +9 suite tests (`tests/test_b6_transonic.py`) incl. the two recorded negative results.
**★ Fold findings (2026-07-12, coarse M0.80 α1.25 vs the G8.1 anchors — Newton shock 0.658/cl_p 0.459/Γ≈0.2295;
conforming Picard's own committed state is a STALL at Γ 0.1819/shock 0.604):** (1) the **live option-a Γ→far-field-vortex feedback has loop
gain > 1 near the fold** — Γ climbs monotonically THROUGH the conforming-Picard value AND the Newton value at flat residual ~5e-5, then
blows up (M_max 37); under-relaxation cannot fix monotone gain > 1 (1+ω(λ−1) > 1 ∀ω>0; ω_γ=0.1 measured = slower divergence).
(2) the per-level lagged vortex + polish (P5 pattern) measures the outer map g(Γ_ff) = 0.1366→0.2189→0.2884→0.4514:
**no fixed point below the isentropic-validity ceiling** ⇒ with a live/lagged vortex this discretization at coarse M0.80 sits PAST the fold
— the P8 conforming-MEDIUM phenomenon one mesh earlier (the LS path lifts a few % higher at equal h, see (4)).
(3) **the López Neumann outlet (B5 option b, no Γ feedback) removes the loop and CONVERGES to near the Newton solution** on BOTH mesh
families (physical, 0 lim/flr, M_max 1.39): M0 embedded Γ 0.2114 (−7.9% of Newton), shock 0.644, cl_p 0.4154;
**M3 wake-free Γ 0.2315 (+0.9% of Newton!), shock 0.678, cl_p 0.4556 (−0.7% of Newton 0.459)** —
both far closer to the truth than the conforming Picard's own stall (Γ 0.1819 = −21%, shock 0.604, cl_p 0.357);
structurally why the dissertation runs all transonic cases on the outlet form. ⇒ **B6 transonic recipe = `farfield="neumann"`** (B5's
subsonic option-a default verdict unaffected). (4) ★ **INVERSION — the LS Picard tracks the conforming NEWTON truth to ≤1%;
the deviator is the conforming Picard itself.** The raw Picard-vs-Picard gap grows with pocket strength (M0.5 +0.2% / M0.65 +0.5% / M0.70
+4.9% / M0.75 +10.5% coarse; medium M0.70 +7.4%) — but same-mesh conforming NEWTON arbitration shows the conforming PICARD under-circulates
by −4.1% (coarse M0.70, Newton Γ 0.1151) / −8.4% (coarse M0.75, Newton 0.1377) / −6.6% (medium M0.70, Newton 0.1190) —
the P4-erratum bias (frozen-Γ inner solves + budgeted secant early-stop) quantified at weak shocks —
while the LS Picard sits at **+0.6% / +1.0% / +0.25%** of the Newton truth respectively, converging TOWARD it under refinement (no
early-stoppable Γ outer exists on the LS path: Γ is a solution mode, converged with the field to residual ~1e-7).
**Gates (dual-mesh; re-anchored 2026-07-11 P4-erratum aware; ★ BASELINE CHANGED 2026-07-12, USER-ARBITRATED: the reference is the same-mesh conforming NEWTON truth, NOT the conforming Picard — the Picard stall under-circulates 4–8% at these shock strengths (finding 4), so it was never a valid A/B target; this aligns the B6 reference with the G8.1 anchor):**
- [x] **coarse M0.80 α1.25° inside the G8.1 Newton-lock bands** with the B6 neumann recipe: **MET** —
      M0 embedded shock 0.644 / cl_p 0.4154 / Γ 0.2114 (−7.9% of Newton 0.2295); M3 wake-free shock 0.678 / cl_p 0.4556 (−0.7% of Newton
      0.459) / Γ 0.2315 (+0.9%); both inside shock 0.658±0.06 and Γ ±10%, physical (0 lim/flr, M_max 1.39).
      (The conforming Picard itself misses these locks by −21% cl / −0.054 shock — the deprecated baseline.)
- [~] **medium M0.7875** (the G8.1 re-specced case): **PARTIAL — bounded Picard-quality state reached, isolated fold solution deferred to
      the LS Newton.** With the fold recipe (C=2.0, θ=0.5 localized, dm=0.025 into the fold) the LS neumann solve on M0 embedded stays
      BOUNDED and physical at every level (M_max 1.44 vs Newton 1.404, **0 limited / 0 floored**, shock captured) —
      but it STALLS at Γ 0.2146, **−18.8% of the same-mesh conforming Newton truth 0.2643**, residual parked at ~5e-6.
      This is the FP non-uniqueness fold: a Picard method — conforming OR level-set — does not reach the isolated solution here (exactly why
      G8.1 re-specced the *conforming* path to Newton locks at medium; P4-erratum / P8 fold record), and the heavier dissipation needed to
      keep the finer-mesh shock bounded depresses the lift further. Earlier attempt with the coarse recipe (C=1.5, θ=0.2) **diverged**
      (M_max 60, 1103 limited) — the medium shock needs more dissipation than coarse. ⇒ the quantitative medium gate needs the **LS Newton**
      (post-B6 re-derivation, design_track_b.md §5.5, explicitly deferred). The dm=0.025 M3 wake-free medium leg is not yet measured (timed
      out behind M0).
- [x] **same-mesh A/B vs the conforming NEWTON solution within 2%** (re-baselined): **MET** where both converge — LS Picard is +0.5% (coarse M0.65) / +0.6% (M0.70) / +1.0% (M0.75) / +0.25% (medium M0.70) of the Newton Γ, and −7.9%/+0.9% (M0/M3 coarse M0.80, the fold, Picard-quality band). The old "Picard-vs-Picard ±2%" reading is retired (measured the conforming stall bias, not an LS error).
- [x] fold discipline applies (per-mesh locks, no cross-mesh convergence claims) — enforced: coarse M0.80 and medium M0.7875 are separate anchors, the dual-mesh spread widens at the fold (M0/M3 straddle Newton) and is reported per mesh, never as a convergence claim.

**★ B6-Newton (post-B6 re-derivation, design_track_b.md §5.5/§10.6; 2026-07-12) — the LS Newton that the medium fold needs.**
`solve/newton_ls.py::solve_multivalued_newton`: exact Jacobian = Picard matrix + per-side Terms 2/3 (P7 sensitivities through the DOF
indirection) + the EXACT quadratic TE-Kutta derivative; wake-LS rows linear (no correction); no Γ DOF (no Woodbury);
nonsymmetric → splu. **FD-verified 1.3e-9** (`tests/test_b6_newton.py`; the Terms-2/3 row-map —
drop non-TE aux, reroute TE aux → TE main — was FD-caught at 1e-4 before the fix). **★ Reaches machine-converged, terminal-QUADRATIC
discrete fold solutions (0 lim/flr) where the Picard only stalled:** coarse M0.80 M0 |R| 9.4e-13 Γ 0.2124 (−7.4% of conforming Newton) / M3
3.2e-11 Γ 0.2322 (+1.2%); **medium M0.7875 M3 wake-free (the workflow mesh) |R| 1.5e-12 Γ 0.2292** —
the fold is a genuine discrete solution on the B path, closing the "is it a solution?" question the Picard stall left open.
**Two honest gaps remain (recorded):** (1) M0-embedded medium live-Newton limit-cycles at |R|~3e-6 (bounded/physical but not
machine-converged) — the P8/N5 near-tie churn in LS form; fix = wire in the N5 frozen-selection/refresh (§5.5 says it transplants, interface
`freeze` reserved); (2) the converged LS fold lift sits ~13% below the conforming-Newton truth at medium (a real discretization difference,
both machine-converged) — apportionment (B5 neumann O(Γ/R) −4% vs cut-integration O(h) vs artificial-density mesh-dependence) is the
recorded next investigation (candidate: a vortex-far-field LS Newton + C-sweep). ⇒ B6 stays IN PROGRESS:
coarse gate met, LS Newton delivered + fold reachable, medium quantitative closure needs those two items.

### B7 — ONERA M6 3D gate ✓ CLOSED 2026-07-12 (was B5.5, orig B4.5; NEW 2026-07-11, user-arbitrated)
**Deliverable:** the 3D-only machinery — TE-polyline ruled level set (D9), g₂ spanwise-free wake BC (D1), tip Γ→0 — is untestable on the 2.5D meshes of B1–B6
**Gate:**
- [x] **M6 coarse vs the P5/P8 baseline: Γ(z) distribution, cl_KJ, and shock positions within A/B bands — MET on BOTH families** (dual-mesh rule), M∞ 0.84 / α 3.06°, `farfield="neumann"`, Mach ramp 0.60→0.84 @ dm 0.04. Demo `cases/demo/b7_onera_m6/` **35/35 PASS**; `tests/test_b7_onera_m6.py` (6 fast + 5 gated).

**Measured (M1 wake-embedded / M4 wake-free; solve 22.7 / 18.4 min coarse):**

| | M1 embedded | M4 wake-free | P5 conforming Picard | P8 conforming **Newton** |
|---|---|---|---|---|
| cl_KJ | **0.2765** (+2.7% of Newton) | **0.2710** (**+0.7%**) | 0.24788 (−8.6% of Newton) | **0.2692** |
| cl_p (3D) | 0.2716 | 0.2656 | 0.24194 | 0.2560 |
| V6 consistency | 1.77% | 1.97% | 2.40% | — |
| shocks η .44/.65/.90 | 0.635/0.588/0.449 | 0.634/0.584/0.454 | 0.596/0.570/0.425 | 0.596/0.541/0.362 |
| Γ root → tip | 0.1076 → −0.0003 | 0.1055 → +0.0003 | 0.097 → 0.0206 | — |
| M_max, limited/floored | 1.453†, **0/0** | 1.368, **0/0** | 1.398, 0/0 | 2.13 |

† M_max re-read 2026-07-14 with the honest `element_mach2(mixed_plain="main")`
(default since the B8-backlog flip): **M1 1.453 (side) → 1.392 (main)** — the
committed 1.453 was itself a beyond-tip mixed-plain artifact cell (honest value
closer to P5's conforming 1.398); M4 and both 2.5D B6 states are **bit-identical**
under either reading. Gate bands unchanged. Artifact:
`cases/demo/b8_tip_taper_ls/results/mmax_reread.csv` (+ `run_b8_mmax_reread.py`).

**★ Finding 1 — the B6 lift INVERSION reproduces in 3D, on the first try.** Gated against the conforming **Newton** truth (the B6
user-arbitrated baseline, not the conforming Picard), the level-set Picard lands **+2.7% (M1) / +0.7% (M4)** of cl_KJ 0.2692 —
while the conforming Picard (P5) sits **−8.6%** below it. This is the same structure B6 measured in 2D and for the same reason:
the LS path has **no early-stoppable Γ outer** (the implicit Kutta makes Γ a solution mode converged with the field), whereas the conforming
Picard's frozen-Γ inner solves + budgeted per-station secant under-circulate (the P4-erratum bias;
P8 independently measured it at +7.9% for M6 medium). Gating B7's lift on P5 would have *penalised the B path for being closer to the truth*
— hence the Newton anchor. Note the **wake-free workflow mesh (M4) is the more accurate of the two**, which is the outcome Track B exists to
deliver.

**★ Finding 2 — the 3D far field: `farfield="neumann"`, and the P5 Γ(z) taper is NOT needed on the B path.**
The B-path vortex (`picard_ls._farfield_main`) is a **span-uniform** 2D point vortex whose branch cut is the ray y=0, x>0 *at every z*.
On M6 that is wrong in two independent, separately-measured ways (demo `farfield_decision.png`;
gated `test_farfield_vortex_is_contraindicated_in_3d`), both showing up as a spurious near-sonic spot at the **outlet, where the sheet
leaves the domain** (max local Mach there, M∞ 0.5):
  - **(a) non-coplanarity** — the α-aimed sheet has climbed to y ≈ x·tan α ≈ 0.5 by the outlet, far off the vortex's y=0 cut, so the outlet carries a prescribed Γ jump **no cut supports**. This is B3's recorded coplanarity rule, now in 3D. Outlet M **0.958** vs neumann **0.513**.
  - **(b) span-uniformity** — re-aiming the sheet coplanar (direction (1,0,0)) *shrinks but does not remove* it (outlet M **0.825**): one scalar Γ cannot match Γ(z)→0, and outboard of the tip there is no cut at all. This is exactly P5's branch-ray artifact, whose conforming fix was the Γ(z) taper.
  ⇒ **neumann carries no vortex, so neither defect can exist** — the taper is unnecessary on the B path rather than merely unimplemented. Cost: B5's O(Γ/R) outlet truncation (a few % of lift on a compact domain), which is why the bands are A/B bands, not <1% bands.

**★ Finding 3 — Γ(z) comes out spanwise-SMOOTH with no smoothing applied** (unplanned; it became visible the moment the real P5 curve was
overlaid — `gamma_of_z.png`). Normalised RMS second difference of Γ(z): **0.0079 (M1) / 0.0091 (M4) vs 0.0970 for the conforming P5 —
an 11–12× reduction.** The conforming path runs a **separate secant per TE station**, so its Γ(z) carries station-to-station jitter (this is
the very defect P5's `INVESTIGATION_gamma_smoothing.md` chased, concluding that spanwise-Γ *smoothing* moves Γ **away** from the
self-consistent value, and it is the same machinery whose single-station failure — st133, 32% under-circulated —
cost P5 an entire investigation). The implicit Kutta has **no per-station loop to be noisy in**:
Γ is one solution mode of the coupled system. So Track B does not merely *fix* the P5 spanwise-Γ problem —
it makes the problem **structurally impossible**.

**★ Finding 4 — the 3D-only machinery works, and is cheap.** Γ(z) decays monotonically root→tip and reaches **~3e-4 at the tip** on both
families — the spanwise clip delivers Γ(tip)=0 *discretely*, the level-set analogue of the conforming free-edge rule.
The swept TE-polyline oblique frame and the g₂ spanwise-free jump gradient (structurally untestable on quasi-2D meshes) needed **no new
code** — B1's fixes held. Cost was far below the plan's risk estimate: the per-outer `spsolve` on ~12k 3D DOFs is ~0.6 s, so a full 7-level
continuation is **~20 min**, not hours.

**Honest caveats (recorded, not chased):**
1. **Convergence semantics = the recorded transonic Picard tail, not `tol_residual`.** The top Mach levels exhaust the 600-outer budget and
   park at |R| ~ 4–6e-6 (M1: levels 0.72–0.84; M4: 0.68/0.76/0.84). The field is **bounded and physical at every level** (0 limited / 0
   floored throughout) and every gate metric is in band, so the gate is asserted on *bounded + in-band*, not on `converged` —
   the same P4/B6 engineering-converged regime. The cure is the LS Newton.
2. **LS Newton on M6 = DEFERRED.** `solve/newton_ls.py` uses a plain `splu`, and P8/N6 measured true-3D LU fill at ~100× the 2.5D cost (it needed lagged-LU). Porting `direct_refactor_every` to `newton_ls` is the follow-up; B6 already demonstrated the fold capability in 2D.
3. **Shock positions sit ~0.02–0.04 c aft of P5** (in band vs P5's ±0.06, and consistent with the higher circulation). Against the P8 *Newton* shocks the η=0.90 station is 0.087 aft — recorded; a shock-position A/B against Newton needs the deferred LS Newton to be a like-for-like comparison.
4. **Coarse only** (medium/fine M6 deferred: cost + fold risk).

### B8 — Level-set tip-edge desingularization ✓ CLOSED 2026-07-14 as CHARACTERIZED-NOT-CURED (user-arbitrated; row-blend 2026-07-13 + re-spec round 2026-07-14 both NEGATIVE with mechanisms pinned; B9 unblocked)
**Motivation:** P13/G13.2 fixed the tip/wake-edge singularity on the **conforming**
path with a spanwise loading taper `Γ_eff(z) = F(z)·Γ_Kutta(z)`
(`constraints/wake.py::tip_taper_factors`, applied on the per-station Kutta
target). That mechanism **cannot be ported literally** to the level-set path:
the LS path has **no Γ DOF** (Γ is a solution mode of the implicit Kutta) and its
TE Kutta row is the **homogeneous** pressure-equality condition
`s·(q_u−q_l) = 0` (`kernels/cut_assembly.py::te_kutta_coo`, RHS ≡ 0), so scaling
that row by F is an algebraic **no-op** (G13.2 finding (8)). The clean analogue —
also finding (8) — is a **convex BLEND** of the pressure-equality row with B2's
continuity weld `continuity_closure_coo` (Ŵ: `aux_k − main_j = 0`), per TE node:

    F_i · [ s·(q_u − q_l) ]  +  (1 − F_i) · [ φ_aux(i) − φ_main(i) ]  = 0

F=1 inboard ⇒ full pressure Kutta (bit-identical to today); F=0 at the tip ⇒
weld ⇒ jump `[φ]=0` at that node ⇒ the tip is unloaded — the structural analogue
of `Γ_eff→0`. The blend is **not** a no-op because the weld row is not
proportional to the pressure row (q_u, q_l come from different element sets, so
`q_u−q_l` is not a jump gradient). **This is a different MODEL from
`Γ=F·Γ_Kutta`** (finding (8)): r_c is calibrated independently and the two-path
comparison is a **physics A/B**, not a model-identity check. F(z) reuses the same
`tip_taper_factors`, fed the per-TE-node spanwise arclength `cm.q[cm.te_nodes]`
with `z_tip = cm.span_length`.
**Deliverable:** `tip_taper` (per-TE-node F array, default None ⇒ ones ⇒
bit-identical) threaded through `MultivaluedOperator.assemble_matrix` (blend in
the `closure="wake_ls"`, `te_kutta="pressure"` branch), the Picard drivers
(`solve_multivalued_lifting`/`_transonic`), and the LS Newton
(`solve_multivalued_newton`: blend residual + Jacobian on the TE aux rows).
Subsonic only — transonic M0.84 convergence stays with G13.3's round-ladder.
**Gates:**
- [x] `tip_taper=None` ⇒ B3/B4/B5/B6/B7 bit-identical; `F≡1` == current
  pressure-Kutta path (bitwise); `F≡0` ⇒ single-valued reduction (Γ≈0 / cl≈0).
  ✓ `tests/test_b8_tip_taper_ls.py` (13 passed); B-suite 59 passed / 0 failed.
- [x] blend-is-not-a-no-op unit test (contrast with the scaling no-op that
  motivates the model). ✓
- [ ] ❌ **mechanism probe — NOT MET, and the reason is the finding (below).**
  (★ user re-affirmed 2026-07-14: this unchecked line is KEPT BY DESIGN as the
  honest record of the characterized-not-cured closure — do not "fix" it.)
- [x] **two-path physics A/B** — ✓ RUN, and it is the decisive measurement:
  the taper bounds the CONFORMING edge and does NOT bound the level-set one.

**★★ RESULT 2026-07-13 — THE ROW BLEND DOES NOT CLOSE B8, BECAUSE ITS PREMISE IS
WRONG.** Demo `cases/demo/b8_tip_taper_ls/` (**12/12**, M6 coarse+medium, M∞0.5,
α3.06, `upwind_c=0` — no limiter, no shock; artifacts `b8_taper_ls.csv/.png`,
`checks_b8.csv`). The blend is **correctly implemented and behaves exactly as its
model predicts**: it converges cleanly (0 lim / 0 flr at every r_c), it
**UNLOADS the tip circulation far past the conforming criterion** (Γ_last ~ h^q
with **q = 4.73**, criterion q ≥ 1), and it is **perfectly LOCAL** (inboard Γ
**+0.01%**, cl_KJ **+0.03%**). **And yet the tip-edge peak STILL DIVERGES under
every taper**: p = **+1.341** untapered → **+1.37 / +1.41 / +1.58 / +1.37**
tapered (larger r_c is *worse*). Three findings:

1. **G13.2's DISCRETE mechanism does NOT transfer.** There, `p ≈ 1 − q` (the
   outermost TE station sheds Γ_last as a concentrated vortex over the last cell
   ⇒ edge ~ Γ_last/h ⇒ q ≥ 1 kills it). Here **q = 4.73 yet p = +1.37**, nowhere
   near 1 − q = −3.73. **Killing Γ_last does not kill the peak.**
2. **The lift cost is ~0 (+0.03%, vs the conforming taper's −1.74%) because there
   is NOTHING TO UNLOAD:** the level-set **implicit Kutta already drives
   Γ(tip) → 0 emergently** (B7 measured ±3e-4). The conforming path *needs* the
   taper because its free-edge rule leaves Γ_last ~ √h (q = 0.44). **The
   level-set path never had that disease.**
3. **★ MECHANISM — where the peak actually lives.** The peak cell is **OUTBOARD
   of the geometric tip** (z/b = **1.0118**, dx = +0.061), it is a
   **`beyond_tip` element** — one the **SPANWISE CLIP refuses to cut**
   (`cut_elements.py`: a crossing needs `q ≤ span_length`) — it is the **SAME
   element tapered or not** (elem 93977), and it is **NOT a small-cut sliver**
   (volume **0.71×** the median, and not even a cut element, so the CutFEM
   small-cut instability is ruled out). ⇒ **the level-set tip singularity lives
   in how the embedded sheet TERMINATES, not in the circulation it sheds.**
   (The untapered p = +1.341 reproduces G13.1's level-set exponent 1.34 exactly,
   so the metric is measuring the right object.)

⇒ **The two paths' tip singularities are DIFFERENT OBJECTS.** Finding (8)'s
"clean analogue" is a faithful analogue of the conforming *model*, but the
level-set path does not have the conforming path's disease — the analogue treats
a patient that is not ill. **The shipped machinery (`tip_taper` on the LS path)
is correct, tested and bit-identical by default; it is simply not the cure.**
**B8 needs a RE-SPEC aimed at the sheet TERMINATION** (the spanwise clip /
beyond-tip zone) — candidate directions: a graded/faded sheet termination, a
ghost-penalty-style stabilization of the clip boundary, or extending the sheet
past the tip. **User arbitration required before re-speccing.**

**Caveat recorded (cost boundary):** the LS path has **no `precond` option** —
`solve_multivalued_lifting` is hardcoded to sparse-direct `spsolve` (a deliberate
B2 decision; GMRES+AMG is the deferred B3+ scaling path). M6 **medium** costs
**~484 s / solve** at 67,426 extended dofs (~1.2 GB RSS). **M6 fine (~450k dofs)
on the LS path would hit the same splu wall P9 hit on the conforming Newton, with
no AMG escape hatch** ⇒ this demo is coarse+medium **by necessity**, and any
fine-mesh LS work needs the deferred GMRES+AMG path first.

**★★ RE-SPEC ROUND 2026-07-14 (user re-spec doc → diagnosis-first → span blend
RAN → NEGATIVE, mechanism pinned). Two findings supersede parts of the
2026-07-13 verdict above:**

1. **★ THE COMMITTED LS EXPONENT p = +1.34 WAS A ×5 METRIC ARTIFACT** (demo
   `run_b8_termination_diagnosis.py`, artifact `b8_termination_diagnosis.csv`;
   solves cached WITH `phi_ext` this time). `element_mach2` reads **mixed-side
   PLAIN elements** (exactly the beyond-tip class) from the aux-substituted
   SIDE field, but a plain element's **assembly is single-valued on MAIN dofs**
   (`mass_conservation_coo` scatters `el[plain]`) — the B6 `own_side_field`
   disease in the one element class own_side_field cannot fix (neither side
   field is the assembled one). At the verdict element (medium, elem 93977):
   **side 1.532 vs main 0.309**. The **HONEST untapered exponent is +0.620**
   (+0.367 with a V/median ≥ 0.1 sliver filter) — **the SAME OBJECT as the
   conforming +0.52, NOT a stronger one** ⇒ the 2026-07-13 "different objects"
   verdict holds for the *mechanism* (TE taper can't reach it) but its
   "LS 1.34 ≥ conforming 0.52" magnitude comparison is RETIRED, and **G13.1's
   LS-vs-conforming exponent comparison needs the same erratum** (its
   conforming numbers are metric-clean and stand). Fix shipped **opt-in**:
   `element_mach2(mixed_plain="main")` — default `"side"` keeps every
   committed diagnostic (B6/B7 M_max locks) bit-identical; **flipping the
   default + re-reading the B6/B7 M_max gate numbers is a recorded user
   arbitration item.** Related recorded (not fixed): the same aux-mixed side
   field feeds `element_densities`, so junk density weights DO enter the
   matrix on mixed-side plain elements (measured rho_up min 0.43 vs physical
   ~0.87 at M0.5 medium; no NaN) — a second arbitration item, since fixing it
   moves every committed LS gate number.
2. **The honest residual singularity is the sheet TERMINATION carrying a
   FINITE jump:** the last cut ring holds |δ| ≈ **0.026** (M6, h- AND
   TE-taper-INDEPENDENT — Γ_last→0 under the TE blend while the ring jump
   does not move; the two are **decoupled**, which is *why* the TE row blend
   measured no effect), dropped to zero across one single-valued element.
   Also recorded en route: the **untapered** emergent Γ(tip)→0 property
   degrades under refinement (Γ_last 0.00011 coarse → 0.00218 medium — B7's
   ±3e-4 was coarse-only).

**★★ THE SPAN BLEND (wake-LS-row termination softening) RAN → NEGATIVE ON
LOCALITY, mechanism pinned** (demo `run_b8_span_blend.py` **8/8**, artifacts
`b8_span_blend.csv/.png`, `checks_b8_span_blend.csv`; machinery
`MultivaluedOperator(span_blend=(form, r_blend))` — per non-TE cut node j,
`w_j·[wake-LS row] + (1−w_j)·s_j·[φ_aux−φ_main] = 0`, w_j from the same
`tip_taper_factors`, s_j = the row's own LS magnitude for h-invariance,
beyond-tip straddler nodes get w=0 at any r_blend; default None bit-identical,
`tests/test_b8_span_blend.py` 11 passed, B-suite 116 passed / 9 skipped).
The blend **hits its target** — the termination-ring jump is WELDED (0.026 →
0.001/0.0003 at r_blend 0.05/0.08 b) — **and the price disqualifies the
model**: (a) **~20% GLOBAL lift loss** (−19.8/−20.2/−21.8% at r_blend
0.03/0.05/0.08 medium), **UNIFORM in z** — Γ(z) scales down root-to-tip
including where the blended rows are bitwise identical to baseline — and
essentially **r_blend-INSENSITIVE** (2-point spread); (b) **component
isolation** (coarse): the straddler weld ALONE costs −13.3%, the inboard
smooth blend ALONE −10.8% ⇒ not a defect of either piece: **the implicit
Kutta has no per-station target — Γ is ONE global solution mode, and ANY
δ-pin on the sheet near the tip re-levels it** (G13.2 finding-(2)'s
fixed-point amplification acting GLOBALLY; the conforming secant keeps it
per-station, which is why the same F(z) costs −1.6% there and −20% here);
(c) the loss **GROWS under refinement** (rb0.08: −16.9% → −21.8%), so it
**CONTAMINATES the exponent**: p_unload ≈ −0.10 of the apparent
+0.37 → +0.05 reduction (corrected ~+0.15 — suggestive, NOT certifiable
under a 20% flow distortion, and moot at this cost); (d) a **narrow blend
(~2 elements, rb0.03) is ACTIVELY harmful** (ring jump ×2.26, honest
p +1.31). ⇒ **every constraint-side route is now measured: TE rows (no
effect), wake-LS rows (globally amplified damage). If a cure is still
wanted, it must change the FUNCTION SPACE at the termination** (how the
spanwise clip ends the multivalued region — e.g. sub-element termination of
the aux-DOF set), not add sheet-side constraints.

**★ ARBITRATED 2026-07-14 (user): B8 CLOSED as CHARACTERIZED-NOT-CURED.**
The honest LS tip exponent (+0.62 / +0.37 no-sliver) is the same object at
the same magnitude as the conforming +0.52 that every closed conforming gate
lives with ⇒ **B9 (wing-body LS solve, M∞0.5) is UNBLOCKED**. Recorded
BACKLOG (not scheduled): (a) flip the `element_mach2` default to
`mixed_plain="main"` + re-read the B6/B7 M_max locks + the G13.1 LS-exponent
erratum; (b) fix the `element_densities` mixed-plain junk weights (moves
every committed LS gate number — needs its own A/B). The function-space
termination re-spec stays a candidate only if a future gate actually needs
the LS tip edge bounded. All shipped machinery (`span_blend`, `mixed_plain`)
is default-inert and stays. Evidence: demo_report "Track B / B8 re-spec"
section.

### B9 — Wing-body cross-model validation (LS + conforming), M∞ 0.5 ✓ CLOSED 2026-07-17 (was B8 2026-07-13; orig B6→B5) ★★ RE-SPEC 2026-07-17 (user-approved — multi-element leg moved out, conforming wing-body leg added)

**★★ RESULT 2026-07-17 — the two wake models AGREE on the wing-body to 0.4% at
medium.** Demo `cases/demo/b9_wingbody/` (7 PASS + 1 XFAIL). Medium M0.5:
conforming-pressure Newton cl_p 0.2173 / cl_kj 0.2188 vs level-set Picard
0.2165 / 0.2175 = **0.4% / 0.6%** (GB9.5 PASS, pre-registered < 1%), the
wing-body analogue of the wing-alone P14 cross-model 0.17%/0.36%; section Cp
overlays pointwise at all 4 stations; ~~the coarse 12.8%/8.7% gap was pure
resolution~~ (**B17 erratum 2026-07-18: the coarse gap is largely far-field
CONTAMINATION, not resolution — the coarse legacy far-field aux carry |jump|=53
garbage; pinning them jump→γ lifts coarse LS to 0.2087 ≡ conforming. The medium
headline stands**). GB9.1 (conforming Newton converges coarse+medium, probe 2/5 +
pressure 3/4 steps, 0 lim/flr) ✓, GB9.2 (LS Picard converges both) ✓, GB9.3
(junction TE-CV wing-side, `tests/test_b9_wingbody_{conforming,ls}.py`) ✓.
- ★ **GB9.4 XFAIL (documented-open, like G1.6): the "fuselage carries no lift"
  premise is FALSE as measured.** Fuselage pressure-lift is 16% (conf) / 20%
  (LS) of the wing's at medium, and the LS value GROWS with refinement
  (0.164 → 0.205) while conforming stays flat — the signature of the G1.6
  smooth-wall flat-facet natural-BC error on the fuselage (GB9.6's subject),
  NOT clean physics. Band NOT moved after the fact (house rule); recorded as
  the open state a wing-body body-surface claim carries until P11/Option C.
  **★★ B28 erratum (2026-07-20): the "spurious lift ⇒ G1.6" reading is
  measured FALSE.** The pre-registered flat-fragment decoupling
  (`cases/analysis/b28_cl_fus_flat_sheet/`) drags the LS sheet GEOMETRY flat
  at y=0 (`WakeLevelSet(sheet_direction=...)`, default None bit-identical)
  while physics convects with the flow: out-band cl_fus 0.0326 vs the
  conforming oracle 0.0351 (**7.25% ≤ 15% TOL**), where the flow-aligned
  sheet gives 0.0504 (35.3%) — wake-sheet POSITION sensitivity, not a
  discretization error. Decomposition: physical carryover baseline (~10% of
  cl_p, both paths, refinement-flat) + pocket imprint (B25, cured) + sheet
  position (B28). GB9.4 re-spec'd per B23 §(c) — see the ledger entry; B9
  demo re-run 8/8 PASS (medium gap 7.0%).
- ★ **LS uses PICARD, not Newton — measured, not assumed (user-questioned).**
  The committed LS Newton recipes (lagged-LU, `precond="schur"` B14, N5 freeze,
  the B15 Mach ramp) ALL diverge/churn on the subsonic wing-body. Diagnostic at
  the converged Picard state: the wake/Kutta rows are PERFECT
  (max|R[te_aux]| = 1.8e-8) but 8 far-field FLUID rows carry |R| ≈ 84 — a
  localized inconsistency in the **`farfield="freestream"` Newton path**, which
  was never exercised (every committed LS Newton run uses `neumann`, and the
  fuselage BLOCKAGE makes `neumann` unbounded, res → 1e43). So the failure is
  the far-field BC layer, orthogonal to the linear solver (Schur/spsolve) and
  the wake model. LS Picard with `freestream` converges cleanly (res 3e-7).
  A `newton_ls.py` freestream bug (Dirichlet values left None) was fixed in
  passing; the residual inconsistency at the outer boundary is a recorded
  follow-up (subsonic LS Newton for a blunt body), not a B9 blocker.
- ★ **Conforming wing-body — the NEW capability.**
  `onera_m6_wingbody_mesh(embed_wake=True)`: the fuselage is built as TWO
  π-revolves (`add_fuselage_solid_split`) so the y=0 meridians are genuine seam
  edges — a single periodic revolve surface with the waterline imprinted to the
  degenerate tail pole is UNMESHABLE (every 2D algorithm fails "1D mesh not
  forming a closed loop"). The sheet passes through the body, fragment-trims to
  exposed TE + waterline + aft symmetry + tip/downstream, and embeds; Netgen is
  OFF for embed (segfaults). `cut_wake` / constraints / P14 pressure-Kutta ALL
  unchanged; the waterline duplicates via the existing boundary-edge rule.
  `embed_wake=False` bit-identical (n_tets 65621 exact). Generation-time
  `cut_wake` ingest gate = the crack detector (all free nodes at the tip).
  Coarse 90099 tets, medium 679391.
- **GB9.6 RECORDED** (`cases/analysis/b9_fuselage_guardrail/`): isolated
  body-of-revolution h-sweep α=0, azimuthal Cp scatter (reference-free) median
  DECAYS 0.0036/0.0022/0.0010 (coarse/medium/fine h_body) but max GROWS
  0.042/0.096/0.117 (nose/tail poles) — the G1.6 error class, quantified.

**Deliverable (re-spec'd):** run BOTH wake models with Newton on the M2
wing-body geometry (coarse + medium, M∞ 0.5, α 3.06° — the committed M6
subsonic convention), each with its best-known recipe, and compare spanwise
lift Γ(z)/cl(z), section Cp, convergence histories, and the A1 timing
breakdown. The conforming leg is a NEW capability (no conforming wing-body
mesh/generator existed; `cut_wake` ValueErrors on the wake-free family) —
the design routes it through a wake-embedded variant of the M2 family
(`onera_m6_wingbody_conforming/`), with the sheet fragment-trimmed to the
fuselage waterline and extended to the symmetry plane aft of the body
(simple-connectivity) and to the far-field sphere (branch-cut Dirichlet
data), reusing the existing `cut_wake` boundary-edge duplication rule
unchanged (the waterline duplicates by the same mechanism as the wing-alone
symmetry root edge).

**Original gates (pre-2026-07-17), kept for the audit trail:**
- ~~two-element cl's plausible~~ — SUPERSEDED (multi-element moved out of B9,
  not scheduled anywhere)
- fuselage carries no lift → now GB9.4
- ★ fuselage surface-Cp guardrail (user-arbitrated 2026-07-14: isolated
  body-of-revolution subsonic h-sweep, error magnitude RECORDED, no pass/fail
  line; caveat carried by every wing-body surface-pressure claim until
  P11/Option C) → now GB9.6, kept verbatim in intent (user re-confirmed
  2026-07-17)

**Gates (pre-registered 2026-07-17, BEFORE any wing-body solve was run):**
- [x] **GB9.1** ✓ conforming wing-body capability: wake-embedded mesh variant
      passes the family quality gates + a generation-time `cut_wake` ingest
      gate (all free nodes at the tip z≈B_SEMI — the crack detector; innermost
      station z≈0.15; waterline wake∩fuselage nodes all duplicated); M0.5
      probe Newton → probe-seeded pressure Newton (P14 recipe) both converge
      on coarse AND medium with 0 limited / 0 floored. **✓ probe 2/5 +
      pressure 4/3 steps, 0/0; meshes 90099 / 679391 tets.**
- [x] **GB9.2** ✓ LS wing-body capability: converges on coarse AND medium;
      TE-node census = the M2 locks (76 coarse / 150 medium) at α 3.06°
      (TE identification is aim-independent; the cut-element census is NOT —
      recorded, not locked). **★ RE-SPEC of the SOLVER within this gate: the
      committed LS Newton recipes (lagged-LU, Schur, N5 freeze, Mach ramp) all
      diverge/churn on the subsonic wing-body — the failure is the
      `neumann` far field (fuselage blockage → unbounded), not the solver; LS
      uses its proven subsonic solver PICARD with `farfield="freestream"` (res
      3e-7). See the RESULT block above.**
- [x] **GB9.3** ✓ junction TE-CV verification (the M2 open item, BOTH paths):
      the junction station's upper/lower Kutta control volumes contain only
      wing-side elements — LS: audit `multivalued.py::_build_te_control_volumes`
      fans against the wing/fuselage face sets (wing-only `wall_nodes` is the
      wiring); conforming: `te_pressure.py::TEControlVolumes` builds adjacency
      from `boundary_faces["wall"]` only — assert construction + nonempty
      junction fans + finite implied targets. **✓ both paths, coarse+medium
      (`tests/test_b9_wingbody_{ls,conforming}.py`).**
- [~] **GB9.4** ✗ XFAIL (documented-open): fuselage carries no lift
      |cl_fus(pressure)| ≤ 0.05 · cl_p(wing) at medium. **MEASURED 16% (conf) /
      20% (LS), LS GROWS with refinement 0.164→0.205 ⇒ the G1.6 fuselage-Cp
      discretization error (GB9.6), not clean physics. Band NOT moved (house
      rule); the caveat a wing-body body-surface claim carries until
      P11/Option C.** Structural half census-locked (zero TE stations inboard
      of the junction). **★★ B28 RE-SPEC (2026-07-20, executes B23 §(c)):**
      the ≤5% band is retired by erratum (physical carryover is not an
      error); the gate is now **out-band cross-model consistency** —
      |cl_fus_out(conf) − cl_fus_out(LS)| ≤ 0.15·|cl_fus_out(conf)| at
      medium (hard), coarse RECORDED, plus a RECORDED band/out/poles
      decomposition (`cross_model_m05.csv` extra columns). TOL=15% = 1.5×
      the out-band refinement noise (8-10%), pre-registered. Re-run: medium
      gap **7.0% PASS**, coarse 10.0% RECORDED — demo 8/8.
- [x] **GB9.5** ✓ cross-model agreement at medium M0.5: conforming-pressure vs
      LS agree to **< 1%** on cl_p(wing) AND on cl_KJ under the SAME
      exposed-span reducer (trapezoid over actual stations + tip closure, NO
      root flat-extension; the conforming flat carry-through value recorded
      as a supplementary column). Band justification: the wing-alone
      cross-model precedent is 0.17%/0.36% at M0.84 with shocks and different
      mesh families; wing-body is same-geometry subsonic but the junction
      discretization is new — 1% = the G14.7 house band. Coarse RECORDED.
      Same-extractor discipline throughout (A2/V14.6): section Cp via
      `section_cp_curve`/`section_cp_curve_levelset` (both walk `wall` only),
      LS TE gap re-measured through the same all-station sweep, never the
      mvop's own CVs. **✓ MEDIUM cl_p 0.4% / cl_kj 0.6% (conf 0.2173/0.2188 vs
      LS 0.2165/0.2175); coarse 12.8%/8.7% RECORDED (~~pure resolution~~ — B17
      erratum: largely far-field CONTAMINATION; pin_gamma lifts coarse LS to
      0.2087 ≡ conforming).**
- [x] **GB9.6** ✓ RECORDED fuselage surface-Cp guardrail (carried from the
      2026-07-14 arbitration): isolated body-of-revolution h-sweep
      (h_body 0.060/0.030/0.015 = the wing-body coarse/medium/fine skin
      resolutions), α=0 non-lifting; primary metric = azimuthal Cp scatter
      per x-bin (reference-free absolute discretization error at α=0),
      secondary = self-convergence deltas (stated as SELF-convergence, with
      the G1.6 sphere calibration quoted — this error class partially hides
      from Richardson). No pass/fail line; the caveat stands until
      P11/Option C. **✓ azimuthal scatter median DECAYS 0.0036/0.0022/0.0010,
      max GROWS 0.042/0.096/0.117 (nose/tail poles) — the G1.6 class.**
- Convergence histories + timing breakdowns: RECORDED-not-gated (no committed
  wing-body anchors exist; these runs become the first).

**Scope guards (carried + updated):** subsonic M∞ 0.5 ONLY (M0.84 excluded —
G13.3 transonic NEGATIVE stands); meshes = the M2 family and its new
conforming variant, coarse + medium only, both gitignored (regenerate before
running); `wall_tag` stays `"wall"` on both paths (widening it to include
`fuselage` would mint spurious Kutta stations along the waterline);
conforming recipe = `M6_NEWTON_KW` (farfield_spanwise_gamma=True,
precond="direct", direct_refactor_every=1000) with the P14 probe-seeded
pressure Newton; LS recipe = lagged-LU primary (A1 anchor lineage) +
`precond="schur"` as a secondary timing row at medium (B14: 2.08× at M6
medium M0.5), γ agreement |Δγ| ≲ 1e-8 per the GB14.4 precedent.

### B10 — Curved wake / free wake ⊘ (was B9 2026-07-13; orig B7→B6; SHELVED 2026-07-10)
**Deliverable:** Curved wake / free wake — **SHELVED 2026-07-10** (DN2 §4.5.6: loading error of a straight wake is O(θ²) ≈ 0.1%; per-update CutElementMap/DOF rebuild cost; discrete cut-set jumps conflict with Newton; López precedent). `update_direction()` interface capability retained.
**Gate:** — (shelved; no gate)

### B11 — LS-path infrastructure: unified post-processing + GMRES/AMG scaling ✓ (NEW 2026-07-14, user-directed; appended after B10, no renumber; CLOSED 2026-07-14)
**Deliverable:** two long-standing LS-path infrastructure gaps closed together
(a B9 enabler — the wing-body medium LS solve would otherwise hit the same
splu wall).

**(1) Unified post-processing.** `post/surface.py` (conforming) and
`post/surface_ls.py` (level-set) now share private cores — `surface._cp_from_q2`
(the per-triangle isentropic/Bernoulli Cp branch), `surface._pressure_force`
(the `-(cp·area)·n_out/s_ref` integral + lift/drag projection),
`section_cut._wall_plane_crossings` + `_resolve_station` + `_section_curve_dict`
(the triangle plane-cut loop + station resolve + chord/x_le normalize),
`surface_ls._d11_wall_state` (the D11 two-sided q² selection, formerly
duplicated twice inside surface_ls). The three near-duplicate blocks
(Cp+D11, the copy-pasted section-cut loop, the force integral) collapse to
one implementation each. New upper layer `post/unified.py` dispatches by
keyword — `wall_cp` / `wall_forces` / `section_cp`, taking `phi=` (conforming)
or `mvop=,phi_ext=` (level-set) — so one call site serves both paths; outputs
are `np.array_equal` to the legacy functions by construction. **Every legacy
public function keeps its name/signature** (14+ demos, 10+ test files), and the
extraction preserved float op order, so the bit-identity locks
(`test_b7_onera_m6.py::test_d11_upper_surface_equals_the_main_based_section`,
the shock `x_shock` asserts, `test_post_surface.py`) pass unchanged. Bonus:
`section_cp_curve_levelset` / `wall_cp_levelset` gain the opt-in `smooth_passes`
(the conforming G6.1 gradient smoothing), default 0 = bit-identical.

**(2) GMRES+AMG on the LS solvers (the deferred design_track_b.md §5.3
landing).** `solve_multivalued_laplace` / `_lifting` / `_newton` grow
`precond=None|"ilu"|"amg"` (None = the pre-B11 sparse-direct `spsolve`,
bit-identical default) plus `linear_rtol=1e-10`, `gmres_restart`,
`gmres_maxiter`, `amg_rebuild_every`; `solve_multivalued_transonic` inherits
them through `**kwargs` (zero code); `newton_ls` adds `seed_precond`. "ilu"/"amg"
run `solve/linear.solve_gmres` on the fused nonsymmetric matrix (the escape from
the M6-fine splu wall). **★ ILU is the effective escape** (`precond="ilu"`,
spilu on the real fused A_free): converges at 434 iters coarse, |Δγ| < 1e-8,
warm-started per outer. **★ AMG does NOT converge on the lifting operator
(measured, honest result).** `_amg_surrogate_preconditioner` builds SA-AMG on an
SPD surrogate (the single-valued Picard block + SPD springs tying each aux dof
to its coincident host so SA aggregates them, §5.3 "把 N_ext 个辅助 DOF 当普通节点
处理") — this works for the `continuity`-closure (Laplace/B2) system, but on the
`wake_ls`-closure lifting/transonic/Newton operator the aux rows are the g₁+g₂
wake-LS + nonlinear TE-Kutta rows (convection-like, not SPD springs), the
surrogate cannot model them, and **GMRES stalls at the restart cap** (coarse
M0.5 lifting: γ 0.0033 vs 0.139, all 80 outers stalled, 455 s vs ILU's 2.7 s).
So AMG stays wired for the SPD Laplace case + as the recorded §5.3 knob, and
**ILU is the shipped lifting escape**. The **Núñez symmetric row assignment
(§5.3 fallback) stays not-prebuilt** — the route that would restore genuine AMG
applicability if ILU's fill ever becomes the bottleneck at fine scale.
**lagged-LU (`direct_refactor_every`) is OUT of B11 scope** (recorded roadmap
follow-up below): B11 ships the iterative escape, not the direct-reuse one.

**Gates:**
- [x] **G11.1 bit-identity:** full suite green with zero committed-number changes
      (D11 array-equal lock + every shock lock pass untouched); the refactor
      adds only the 9 `test_b11_post_unified.py` tests; `precond` default None
      pinned by `test_b11_linear_ls.py::test_precond_default_is_none`.
- [x] **G11.2 unified == legacy:** `np.array_equal` on cp/q2/section outputs,
      exact-float on cl, both paths (`test_b11_post_unified.py` + demo
      `run_b11_unified_post.py` self-check column max|Δcp| = 0.0).
- [x] **G11.3 GMRES correctness:** ILU reproduces spsolve on the coarse
      (+ demo medium) 2.5D meshes — Laplace/lifting/Newton |Δγ| < 1e-8 subsonic,
      converged, 0 stalls (`test_b11_linear_ls.py`; gated transonic-forwarding
      smoke |Δγ| < 1e-6). AMG reproduces spsolve on the SPD Laplace only;
      measured to STALL on the wake_ls lifting operator (455 s, non-converged)
      ⇒ ILU is the shipped lifting escape (an honest §5.3 finding).
- [x] **G11.4 scaling headline — the splu wall quantified:** M6 medium LS A/B
      CSV committed (`cases/demo/b11_ls_infra/results/m6_medium_ab.csv`, gated
      one-shot). **spsolve = 454.8 s at 67,426 dofs** (the splu wall; the P9
      catastrophe at 450 k fine). ★ **Honest finding:** the M6-medium 3D fused
      matrix resists cheap incomplete factorization — ILU-GMRES advances ~17 of
      26 outers, then even shifted-MILU at fill 20 goes singular at a hard
      outer; near-full fill (≈spsolve cost at this size) would be needed ⇒ **at
      67 k dofs spsolve is still the right tool, ILU is not advantageous there.**
      The escape's payoff is the FINE-scale regime where spsolve is impossible on
      memory (feasibility, extrapolated); the escape is *demonstrated to
      converge* at 2.5D medium (|Δγ| 7.5e-10, 0 stalls; `solver_ab.csv`).
      **Follow-ups:** the Núñez symmetric row assignment (§5.3 — would restore
      AMG and a cheaper factorization; still not scheduled) and the
      `direct_refactor_every` (lagged-LU) port into `newton_ls` (roadmap "LS
      Newton on M6 = DEFERRED") — the latter **executed as B12 (2026-07-14)**.

**Evidence:** tests `tests/test_b11_post_unified.py` (9) +
`tests/test_b11_linear_ls.py` (10 + 1 gated); demos
`cases/demo/b11_ls_infra/` (`run_b11_unified_post.py`, `run_b11_gmres_ls.py`,
`run_b11_m6_headline.py` [gated]). Conforming solver numerics byte-untouched;
no Numba kernel or COO-assembly path touched (pure SciPy/PyAMG + numpy).

### B12 — Lagged-LU direct-reuse for LS Newton (medium/M6-scale enabler) ✓ CLOSED 2026-07-14 (NEW, user-directed; appended after B11, no renumber; executes the B11/G11.4 recorded follow-up "LS Newton on M6 = DEFERRED")
**Deliverable:** make the level-set Newton solve affordable at medium/M6 sizes by
porting the conforming N6 **lagged-LU direct-reuse** mechanism into
`solve/newton_ls.py`.

**Why this and not the B11 iterative escapes.** B11 measured (G11.4 A/B,
`cases/demo/b11_ls_infra/results/`) that the iterative escapes fail on the fused
level-set matrix beyond coarse — so at medium/M6 sparse-direct is the only
converging tool, and the cost driver is then the **number of factorizations**:

| case | dofs | spsolve | ILU | AMG |
|---|---|---|---|---|
| 2.5D coarse lifting | 6,614 | ✓ 1.9 s | ✓ (\|Δγ\| 9.5e-9) | STALL |
| 2.5D **medium** lifting | 22,880 | ✓ 8.6 s | ✗ **diverges** γ=−137, 77 stalls | STALL |
| **M6 medium** lifting | 67,426 | ✓ 454.8 s / 26 outer / **17.5 s per factor** | ✗ **factor_failed** | STALL |

LS Newton with `precond=None` factorizes on **every** Newton step, so on M6
medium it pays that 17.5 s once per iteration.

**Implementation (`solve/newton_ls.py`, the ONLY production change).** Two new
kwargs on `solve_multivalued_newton`: `direct_refactor_every: int = 1` and
`direct_reuse_rtol: float = 1e-8`, active only on the `precond is None` branch.
`k=1` (default) = the byte-identical per-step `spsolve`. `k>1` = refactor the LU
(`spla.splu`) every k-th step and drive the intermediate steps with GMRES on the
FRESH Jacobian preconditioned by the stale (exact) LU, converged to
`direct_reuse_rtol`; a reuse step whose GMRES fails falls back to a refactor +
exact solve in the same iteration (robustness never below every-step-direct).
This is the N6 mechanism (`solve/newton.py::direct_refactor_every`) **minus the
Woodbury** — the level-set system has NO Γ DOF, so the step is a plain
`J_free d = −R_free` with no low-rank coupling. New monitor `n_refactor`.

**Gates:**
- [x] **G12.1 bit-identity — CLOSED 2026-07-14.** `direct_refactor_every=1`
      (default) is byte-identical to the pre-B12 `spsolve` path (same `phi_ext`,
      `n_refactor==0`, `n_gmres_total==0`); the two params default to 1 / 1e-8.
      `tests/test_b12_lagged_lu_ls.py::test_lagged_lu_param_defaults` +
      `test_default_bit_identical_to_spsolve`. B6/B11 Newton locks unchanged.
- [x] **G12.2 numerical equivalence (core gate) — CLOSED 2026-07-14.** On the
      coarse 2.5D mesh at M0.70 (upwind active), `direct_refactor_every` 2 and
      1000 both reach the spsolve default's converged γ (0.1778053693) to
      **bit-identity** (|Δγ| < 1e-8), 0 lim/flr, 0 GMRES stalls, and actually
      reuse the stale LU: **k=1000 refactors ONCE** over 6 Newton iters (vs 5
      spsolve factorizations at k=1), k=2 refactors 3× — `n_refactor < n_newton`
      and a reuse GMRES step ran. `tests/test_b12_lagged_lu_ls.py::
      test_lagged_lu_matches_spsolve` (k∈{2,1000}).
- [x] **G12.3 scaling headline (gated) — CLOSED 2026-07-14.** M6-medium subsonic
      (M0.5, α3.06, `farfield="neumann"`, 67,426 dofs) LS Newton A/B from a shared
      Picard seed. Both take **7 Newton steps to a genuine converged solution**
      (0 lim/flr); spsolve refactors all 7 (**145.6 s**), lagged-LU (k=1000)
      refactors **ONCE** + 30 reuse-GMRES iters (**66.7 s = 2.18× faster on the
      Newton phase**), 0 stalls, **γ bit-identical** (0.06685284, |Δγ| = 6.74e-13).
      Demo `cases/demo/b12_lagged_lu/run_b12_m6_newton.py` (**6/6 PASS**),
      `results/m6_newton_ab.csv` + `checks_m6_newton.csv`.
      **★ Honest boundary:** at 67k dofs one splu fits in memory ⇒ this is a REAL
      runnable medium-scale win, not an extrapolation. But lagged-LU still needs
      ≥1 in-memory splu, so it does **not** break the FINE-mesh memory wall (P9's
      26 GB / 4h39m is per-factorization, not per-count) — that regime remains
      the Núñez symmetric-row-assignment → AMG route (design_track_b §5.3, not
      prebuilt).

**Out of scope (recorded):** the LS-Newton Mach-continuation *ramp* wrapper
(`newton_ls` is single-Mach-level; the conforming `solve_newton_transonic`
analogue is a separate follow-up), and genuine AMG applicability (Núñez rows).

**Evidence:** `tests/test_b12_lagged_lu_ls.py` (4); demo
`cases/demo/b12_lagged_lu/` (6/6, G12.3 headline CSV). `solve/newton_ls.py` is
the only production change; default path byte-identical. Suite +4.

### B13 — Lagged-LU on the Picard outer loop (the post-B12 cost driver) ✓ CLOSED 2026-07-14 (NEW, user-directed; appended after B12, no renumber)
**Deliverable:** the B12 lagged-LU mechanism applied to the Picard OUTER loop
(`solve_multivalued_lifting`; transonic inherits via `**kwargs`) — after B12
the M6-medium cost driver is one 17.5 s spsolve per Picard outer (B11 lifting
headline 454.8 s / 26 outers; the B12 demo's Newton seed 263 s / 15 outers).
User goal arbitrated 2026-07-14: **compute speed at medium scale is the
objective; fine-mesh extension is optional** — which ranks lagged-LU (this)
above the structural preconditioner (B14, designed-not-scheduled below).
`solve_multivalued_laplace` is excluded — it is a single-shot solve, nothing to
amortize.

**External-analysis corrections (recorded; GLM analysis + comparison doc,
baseline f9d400a, both predating B12):** (1) "lagged-LU port = not done" was
true at their baseline; B12 landed it for Newton the same day. (2) The Schur
direction in both docs is inverted — the efficient elimination removes the
SMALL aux block (`K = J_mm − J_ma·J_aa⁻¹·J_am`, J_aa an n_ext×n_ext thin-strip
matrix ~8k at M6 medium), not the main block (which would need A_mm⁻¹ = an AMG
inner solve per application). (3) J_aa is not fully constant — wake-LS rows are
(§5.5) but the TE-Kutta rows (76–150) re-linearize each outer; refactoring the
thin strip is milliseconds, a non-issue. (4) 454.8 s is 26 outers, not one
splu (17.5 s each).

**Gates (GB13.x — deliberately NOT G13.x, which is P13's namespace; Track V's
GV prefix is the precedent):**
- [x] **GB13.1 bit-identity — CLOSED 2026-07-14.** `direct_refactor_every=1`
      (default) is byte-identical to the per-outer `spsolve` (same `phi_ext`,
      `n_refactor==0`, `n_gmres_total==0`); defaults pinned.
      `tests/test_b13_lagged_picard.py`.
- [x] **GB13.2 equivalence (core) — CLOSED 2026-07-14.** Coarse 2.5D lifting
      M0.5: k∈{4,1000} reach the spsolve γ to <1e-8, converged, 0 stalls,
      `n_refactor < n_outer` (k=1000 refactors ONCE); also under
      `farfield="neumann"` (the B6/B7 recipe) at k=8.
      ★ **Measured finding: `direct_reuse_rtol` must default 1e-10, NOT B12's
      1e-8** — a Picard fixed point is pinned only by its lag tolerances
      (1e-6), so an inexact reuse step SHIFTS the stopping point (|Δγ| 8e-8 at
      rtol 1e-8), whereas Newton's terminus is pinned by `tol_residual`
      regardless; 1e-10 restores <1e-8 agreement for ~1–2 extra Krylov iters
      on a near-exact preconditioner.
- [x] **GB13.3 M6-medium headline (gated) — CLOSED 2026-07-14.** B11-headline
      lifting (M0.5, α3.06, neumann, tol 1e-7, 67,426 dofs, 26 outers both):
      spsolve **447.6 s** (17.2 s/outer) vs lagged-LU (k=1000) **68.3 s**
      (2.63 s/outer) = **6.55× faster**, 2 refactors vs 26, γ bit-identical
      (0.06685270, |Δγ| 6.9e-13). The 1 GMRES "stall" is the designed safety
      net — an early outer's large density move exhausts the stale LU and
      triggers an extra refactor (hence 2, not 1), never a divergence.
      **End-to-end seed+Newton** (the B12 pipeline, both mechanisms on): seed
      **42 s** (1 refactor / 15 outers, was 263 s spsolve) + Newton **69.9 s**
      = **111.9 s total vs ~330 s post-B12 baseline (~3×)**, Newton γ 0.06685284
      in the B12 lock band. Demo `cases/demo/b13_lagged_picard/` (**6/6 PASS**),
      `m6_lifting_ab.csv` + `m6_end_to_end.csv`.
      **★ Honest boundary:** amortizes the factorization COUNT; still needs
      ≥1 in-memory splu ⇒ does NOT break the fine-mesh memory wall (that is
      B14's unique value).

**Result headline:** M6-medium lifting **447.6 s → 68.3 s (6.55×)**; end-to-end
seed+Newton **~330 s → 111.9 s (~3×)** — the LS medium workflow is now the same
order as conforming M6 medium (solve 140–240 s). **Evidence:**
`tests/test_b13_lagged_picard.py` (5); demo `cases/demo/b13_lagged_picard/`
(6/6, GB13.3 CSVs). `solve/picard_ls.py` is the only production change; default
path byte-identical. Suite +5.

**★ Workflow evidence riding on B12/B13 (2026-07-15, demo — not a gate):
`cases/demo/m6_medium_ls_workflow/` (10/10 self-checks; full record in
demo_report "M6 medium level-set WORKFLOW").** M6 medium wake-free LS solves
BOTH subsonic M0.5 (cl 0.212, strict-converged) and transonic M0.84 (cl 0.276,
M_max 2.455, bounded/engineering-converged, B7 semantics) at
conforming-comparable cost. Mesh A/B (wake-free vs embedded) cl within
0.62%/0.85%; method A/B (LS vs conforming) 0.47% subsonic / +10.65% transonic
(the B6/B7 conforming-Picard under-circulation, not an LS error). Load-bearing
recipe finding: the M6-medium transonic **Picard residual plateaus at
1e-5..1e-4** (P4/B6/N5 shock-position soft mode) — use `tol_residual=1e-5` or
every level burns its full budget (~1 h); a strict-converged transonic wants
the LS Newton ramp instead — which **B15 then delivered** (GB15.4 removes this
plateau, 3.5×).

### B14 — Schur-eliminated aux block + AMG(SPD Picard main block) ✓ CLOSED 2026-07-17 (opened + built same day, user-directed)

**Built.** `precond="schur"` on both LS drivers, in `pyfp3d/solve/schur_ls.py`
(`SchurReducedSystem` + `main_block_preconditioner` + `jaa_diagnostic`), shared
by `solve_multivalued_newton` and `solve_multivalued_lifting` (transonic
wrappers inherit via `**kwargs`). Per Newton step / Picard outer: factor the
aux thin-strip `lu_aa = splu(J_aa)` (n_ext-sized — 1004/3701 dofs at M6
coarse/medium, split+factor ≤ 19 ms), run GMRES on the matrix-free reduced
operator `K = J_mm − J_ma·J_aa⁻¹·J_am` preconditioned by AMG on the SPD
single-valued Picard block restricted to main-free (the conforming
`solve/newton.py` analogue), back-substitute the aux part exactly. The AMG
hierarchy is invalidated at the three freeze selection-epoch sites alongside
`lu_direct`. A stalled reduced GMRES falls back to a full fused spsolve in the
same step (`n_schur_fallback`, the lagged-LU safety-net pattern) — never
triggered in the whole campaign.

**Gates.**
- [x] **GB14.1** diagnostic-first — J_aa factors on every measured case; 1-norm
  cond estimate finite (2.5D coarse 5.1e8 / 2.5D medium 8.2e9 / M6 coarse 6.5e6
  / M6 medium 7.4e7). Measured, not assumed — the constant-jump null vector
  mixes main+aux columns and the TE-Kutta rows pin the level, so J_aa is
  generically nonsingular, and it is.
- [x] **GB14.2** correctness — 2.5D coarse lifting + Newton M0.7 land on the
  spsolve γ (\|Δγ\| 4.2e-11 / 2.0e-12), 0 stalls / 0 fallbacks. This is the
  exact operator where the B11 spring surrogate stalled to γ 0.0033 (vs 0.139):
  no springs, no bias, the circulation mode survives.
- [x] **GB14.3** the discriminating tier — 2.5D **MEDIUM** lifting, where ILU
  DIVERGED (γ = −136.99, 77 stalls, `b11_ls_infra/solver_ab.csv`): schur
  converges to γ 0.14137632, \|Δγ\| 9.3e-10, 0 stalls/fallbacks. "Passing there
  is what a real escape means."
- [x] **GB14.4** 3-D capability — ONERA M6 wake-free COARSE + MEDIUM, subsonic
  M0.5 lifting AND transonic M0.84 Newton ramp, all converge / target-reached;
  γ matches the same-session lagged-LU arm to \|Δγ\| ≤ 1.5e-8 and the committed
  GB15.4 state exactly (γ **0.088338**, M_max **2.4938**).
- [x] **GB14.5** inertness — default `precond=None` byte-identical
  (`np.array_equal` on `phi_ext`), `n_schur_fallback`/`n_gmres_total` inert at 0.

**★ Timing (RECORDED, not a gate — the design pre-declared the medium-scale
gain "uncertain"; it landed on the winning side).** M6 medium fresh
same-session A/B (identical lagged-LU seeds):

| M6 medium (63,100 dofs) | lagged-LU | schur | speedup | precond share |
|---|---|---|---|---|
| M0.5 subsonic lifting | 73.2 s | **35.2 s** | **2.08×** | 51.7% → **5.2%** |
| M0.84 transonic ramp | 671.2 s | **469.3 s** | **1.43×** | 43.6% → **2.6%** |

The A1 bottleneck (LS Newton M6 medium M0.84 = 42.6% precond, A1's own
measurement; the same-session B14 A/B reads 43.6%) is structurally
gone: the full-size `splu` factorizations lagged-LU still needed (2 per solve,
17.5 s each) are replaced by a thin-strip LU + AMG V-cycles. That beats the
user's stated <10% target on both regimes. γ is bit-close to the committed
GB15.4 (\|Δγ\| 2e-13 on the ramp).

**★ Honest limit — where schur is SLOWER.** At small scale it LOSES: 2.5D
coarse/medium and M6 coarse are 3–6× slower with schur (the direct solve there
is already trivially cheap — the reduced-GMRES iteration count costs more than
the tiny factorization). The win appears only at M6-medium size and grows with
the mesh, exactly as the design predicted ("marginal at medium, the unique
value is the FINE memory-bounded path"). The **fine-scale route (AMG O(n) +
thin-strip LU, no full-size splu that cannot fit in memory) remains the
designed, unbuilt use-case** — out of scope here by user direction
(coarse+medium only). Fallbacks (block-triangular; Núñez additive row
assignment) were NOT needed — the aux block factored and GMRES converged on
every case.

Evidence: `tests/test_b14_schur_ls.py` (9, incl. the gated GB14.3 medium
escape); demo `cases/demo/b14_schur_precond/` (**7/7 self-checks, incl. gated
M6 coarse+medium A/B**; the M6-medium states cache to gitignored
`results/*.npz`, committed evidence = the 5 PNGs + `schur_ab.csv` +
`jaa_diag.csv` + `checks.csv`).

<details><summary>Original design snapshot (2026-07-14, now built)</summary>

**Why not now (user-arbitrated 2026-07-14):** at medium scale its marginal gain
over B13 is uncertain (lagged-LU already amortizes to ~1–3 factorizations per
solve ≈ 35 s; Schur+AMG trades that for per-outer GMRES+AMG cost of the same
order), and its unique value — the FINE-scale memory-bounded path (AMG O(n) +
thin-strip LU, no full-size splu) — addresses a regime the user has declared
optional. **Trigger:** GB13.3 lands and medium is still too slow, or a real
M6/wing-body FINE campaign is scheduled.

**Design snapshot (ready to build):** new `precond="schur"` on the LS drivers.
Free dofs split main-free/aux (aux are never Dirichlet — the B3 load-bearing
fact). Per outer/Newton step: `lu_aa = splu(J_aa)` (n_ext×n_ext thin strip,
~8k at M6 medium, milliseconds; TE-Kutta rows re-linearize per step, so
refactor per step). Reduced operator matrix-free:
`K x = J_mm x − J_ma·lu_aa.solve(J_am x)`; reduced RHS
`r = b_m − J_ma·lu_aa.solve(b_a)`; back-substitution
`φ_a = lu_aa.solve(b_a − J_am φ_m)`. Preconditioner =
`build_amg_preconditioner(op.assemble_matrix(rho_own))` restricted to
main-free — the exact conforming analogue (AMG on the SPD Picard block,
constraints eliminated exactly), **with NO springs**: the B11 surrogate's
mismatch (springs bias the solution toward jump≈0, killing the global
circulation mode — γ 0.0033 vs 0.139) disappears structurally because no aux
dof survives into the preconditioned system. GMRES then faces "elliptic +
cut-strip-localized correction" — the operator shape the conforming path
already proved AMG-preconditionable. **Diagnostic-first gate:** J_aa
invertibility/conditioning (the constant-jump null vector mixes main+aux
columns ⇒ J_aa generically nonsingular, TE-Kutta pins the level — measure,
don't assume); the discriminating tier is **2.5D medium lifting, where ILU
DIVERGED (γ=−137)** — passing there is what "a real escape" means. Fallbacks:
block-triangular preconditioner; last resort = the Núñez additive symmetric
row assignment (§5.3 — a discretization change with penalty-weight
calibration, demoted to third line).

</details>

### B15 — LS Newton transonic ramp + N5 freeze-selection ✓ CLOSED 2026-07-15 (NEW, user-directed; appended after B14, no renumber)

**Why (the cost driver, measured):** the LS transonic **Picard**
(`solve_multivalued_transonic`) is a Mach ramp whose top levels park on the
**shock-position residual plateau** (the P4/B6/N5 soft mode) and burn their whole
outer budget there. On the 2026-07-15 M6-medium workflow solve the ramp is 7
levels (0.60→0.84, dm 0.04) and the embedded per-level cache shows levels **0.80
and 0.84 do NOT converge** — each runs its full 200-outer budget — which is the
bulk of the **24.5 min (embedded) / 38.4 min (wake-free)** wall clock. `tol_residual`
is already set to 1e-5 *above* the plateau (1e-7 would burn ~1 h). A Picard method
cannot do better: the plateau is intrinsic. **Newton has no shock-position soft
mode** — the demo's own note says "a strict-converged transonic wants the LS Newton
ramp". But `newton_ls` could not run a ramp: `freeze=` was a reserved no-op, the
convergence gate hard-requires 0 limited/floored (which shock limiter cells block),
and there was no Mach-ramp wrapper. B15 supplies all three.

**Gates**
- [x] **GB15.1 — frozen per-side selection + FD.** `MultivaluedOperator.newton_side_data(frozen=…)`
      + new `freeze_side_state` (per-side `(upstream, branch)` capture). The
      `kernels/upwind.py` frozen apparatus is reused **unmodified** — the per-side
      ops are plain walk-mode `UpwindOperator`s with a same-side-masked face graph
      — so this is wiring, not new numerics. Residual/Jacobian extracted into
      `LSNewtonSystem` so the solver and the FD gate share ONE assembly path.
      **FD: rel 6.7e-9** (eps 1e-5; clean round-off scaling 5.8e-8 / 6.0e-7 at
      1e-6 / 1e-7 ⇒ a true derivative), 96.9% of free rows kept by the ε-guard,
      on a real pocket (nu_max 0.785, 1118 elements on branches 1/2).
      Frozen sweep reproduces the live density **bitwise** at the freeze point.
- [x] **GB15.2 — the freeze cures the limit cycle, and does not move the answer.**
      ★ On **NACA coarse M0.75** the LIVE LS Newton **does not converge**: it parks
      in a genuine **period-6 limit cycle** (3.2e-7, 2.8e-7, 2.7e-7, 1.3e-6, 8.6e-7,
      4.3e-7, repeating) at |R|≈2.7e-7 — three orders above tol — with **0
      limited/floored** (a CLEAN stall = the assignment churn). Arming the freeze
      converts it to a converged solve: **22 steps → |R| 8.5e-13**, 0 reverts, and
      **γ 0.218809 vs the live cycle's 0.218804** ⇒ the freeze removes the churn, it
      does not select a different state.
      ★★ **TRIGGER ERRATUM (measured; the conforming recipe does NOT transfer):**
      `solve/newton.py` also freezes on `live_stalled`. Porting that verbatim makes
      the LS solver freeze a **still-MOVING** assignment — the LS live residual
      bounces ±2× for tens of steps *while still descending* (γ travels 0.183→0.243
      over that stretch: slow progress in a stiff direction, **not** a stall). The
      frozen step then diverges → revert → re-arm: **3 reverts, no convergence**, on
      a case (medium M0.75) the untouched live path converges (54 steps, 7.5e-12).
      With the stall trigger removed it freezes late (|R|<1e-6, assignment settled)
      and converges: **53 steps, 2.1e-12, 0 reverts, exactly the live γ 0.243305**,
      with `residual_unfrozen` 2.1e-12 confirming the LIVE selection agrees there.
      ⇒ **the freeze arms on `freeze_tol` ALONE.** Fail-safes added: disarm after
      `freeze_max_reverts` (3) so the freeze can only ever HELP, never cost
      convergence; the reported `n_limited`/`n_floored` are always re-read LIVE
      (a frozen finish shows 0 floored BY DESIGN and can never be its own evidence).
- [x] **GB15.3 — the Mach-ramp wrapper.** `solve_multivalued_newton_transonic`:
      upward `mach_schedule`, warm start from the last CONVERGED level only,
      dm-halving retry inserted BELOW a failed level and run STRICT, optional
      `upwind_c_post` staging, honest `target_reached`/`m_final` (the P13/G13.3
      erratum: never census a state whose ramp did not reach the target).
      `intermediate_tol` = loose stopping tol on intermediate levels, strict final.
      ★ **LS-specific: the freeze stays ARMED on loose intermediate levels**, unlike
      the conforming mask (`newton.py:888` sets `freeze_tol=None` there): the accept
      gate requires 0 limited/floored on EVERY route, and on a 0.60→0.84 ramp the
      shock forms MID-ramp, so those levels carry limiter cells and can only reach a
      0-clamped accept THROUGH a freeze. Loosen the tolerance, keep the mechanism.
      (The conforming fold contraindication — a loose level leaving an untracked Γ
      seed, G10.2 — has **no analogue here**: the LS path has no Γ DOF, Γ is a
      solution mode carried inside `phi_ext`.) A/B on NACA coarse M0.80 + fold.
- [x] **GB15.4 — M6 medium M0.84: the plateau is GONE, 3.5× faster.** ONERA M6
      medium wake-free (63,100 ext dofs), M0.84/α3.06. **Picard (committed):
      2304.7 s = 38.4 min, residual parked on the 1e-5..1e-4 plateau, top two
      levels burning their full 200-outer budget. Newton ramp (B15): 657.4 s =
      11.0 min (3.51×; committed `summary.csv` — an earlier draft's 672 s was a
      pre-CSV trial run), EVERY level converged to ~1e-11**, the freeze armed at
      every level with **0 reverts**. Per level: M0.60 ✓5 / 0.65 ✓19 / 0.70 ✓23 /
      0.75 ✓16 / 0.80 ✓19 / **0.84 ✓16 steps, |R| 6.9e-11**. Physics cross-checks
      against the committed Picard: **M_max 2.4938 vs 2.4549** (1.6%), **3 clamped
      cells of 330k vs Picard's ≤3**.
      ★ **HONEST LIMIT:** most levels accept via `assignment_cycle` — the FROZEN
      system converges to ~1e-11 and is accepted at the **assignment-discontinuity
      floor** (the live residual stops improving across refreshes). That is the
      N5 semantics the conforming path also uses; it is **NOT** a claim that the
      LIVE residual is below 1e-10. It beats the Picard plateau by 6–7 orders,
      but "live-strict solution" would be an over-claim.
      ⇒ Closes the deferred **B6-medium quantitative** and **B7-quantitative**
      items (the LS Newton on M6 was DEFERRED at B7).
      ★★ **ERRATUM TRAIL (2026-07-19, two reversals — read both).**
      (1) Under the first B20 re-baseline the SAME committed recipe reached
      only **M0.6625 (2/5 levels)**, γ 0.088338 → 0.071909, demo **17/20**;
      GB20.7 recorded it as "a REAL capability loss on the freeze_tol axis".
      (2) **B21 then found the actual mechanism and RESTORED the capability:**
      the stall was B20's own patch gap — `freeze_side_state` captured the
      frozen selection on the unpatched side field (Kimi N1) — and with the
      capture aligned the recipe reaches **M0.84 again**: γ **0.088343**
      (vs 0.088338 here), M_max **2.4818** (vs 2.4938), res 9.0e-14,
      **0 lim/1 flr** (vs 3 clamped), **515 s** (vs 657 s). So this gate's
      capability clause STANDS; its numbers carry the small B20+B21
      re-baseline above (demo refresh = recorded follow-up). What remains
      superseded: the M_max cross-check against the LS Picard was COMMON-MODE
      (both solvers shared the contaminated density; the conforming record is
      1.995 on a different mesh family — an open cross-model question).

★★ **FOUR ERRATA — porting the conforming N5 recipe is NOT mechanical** (the same
lesson B8 taught). Every one was forced out by measurement, none was foreseen:
1. **The TE polyline must come from the AUTHORITATIVE geometry.** A hand-rolled
   `x_te(0)=0.8059` vs `wing3d.x_te(0)=0.80611` — off by **2e-4**. `CutElementMap`
   finds TE nodes by matching the polyline onto WALL NODES (M2: the M6 TE endpoints
   are exact wall nodes), so 2e-4 matches **nothing** ⇒ **0 TE nodes ⇒ no Kutta ⇒ Γ
   unpinned ⇒ 340k limited cells + NaN** — and the solver **passed silently**.
   ⇒ Both LS solvers now **raise** on `te_nodes == 0`, pointing at `meshgen.wing3d`.
2. **★ `freeze_tol` must sit ABOVE the CHURN FLOOR, and that floor RISES with Mach**
   (measured: **<1e-6 at M0.60 → 8.6e-6 at M0.65 → 2.7e-4 at M0.70**). Below it, a
   discrete upwind-selection flip throws the residual back before the freeze can arm
   (measured: clean descent to 8.6e-6, then ×300 to **the same value 2.6e-3, twice** —
   the signature of a discrete flip) and the ramp dies. **Same law as "tol_residual
   must sit above the Picard plateau".**
3. **★ Residuals are NOT comparable across a SELECTION EPOCH.** The frozen phase
   drives `r_best` to 1.5e-11; after a refresh the residual legitimately returns to
   the live scale (2.6e-3), and the fail-fast (`res > 100*r_best`) reads a 1e8×
   "blow-up" and kills a perfectly healthy freeze-refresh cycle. ⇒ `r_best` is reset
   on every freeze / refresh / revert.
4. **★ The frozen phase's clamp count is STALE BY CONSTRUCTION and must not gate
   acceptance.** Under a freeze `n_floored` counts `branch==3` = the cells clamped
   **at the freeze point**; it never falls. Gating on `n_flr == 0` therefore
   **refuses a 7.8e-14 machine-precision solution forever** (measured at M0.70: the
   freeze cured the period-7 limit cycle, the floored cell **cleared itself** in the
   live field — final live `lim/flr = 0/0` — and the gate still would not fire).
   ⇒ The frozen phase need only be no worse than at the freeze; the **LIVE
   re-evaluation** in the honesty branch is the arbiter (and it is strict).

★ **New knob `freeze_max_clamped`** (default **0** = the conforming N5 rule, bit-identical).
At M6 medium M0.70 a **single** persistently-floored cell (of 330k) blocks the freeze at **any** `freeze_tol`.
The frozen sweep **represents a clamped cell exactly** (branch 3: `nu=0`, `rho=rho_floor`, `s_e=s_u=0` —
a flat clamp with zero derivative), so the 0-clamped precondition is stricter than the machinery needs;
relaxing it lets the freeze arm and the ramp completes.
⚠ **TWO CORRECTIONS to an earlier draft of this entry (2026-07-15, self-caught):**
  (a) **The clamped cells do NOT "clear themselves".** That was over-generalised from ONE isolated 80-step run at M0.70 (driven to 7.8e-14,
  ending 0/0). In the SHIPPED ramp — which accepts at `assignment_cycle` after ~23 steps — the cells **PERSIST**:
  M0.70 `0/1`, M0.75 `0/1`, M0.80 `1/1`, **M0.84 `1/2` = 3 clamped cells** (which is exactly the Picard's ≤3, so it is consistent, not
  alarming). The freeze proceeds **WITH** them present.
  (b) **The convergence semantics ARE relaxed** — the earlier "the convergence gate is untouched" was FALSE.
  With `freeze_max_clamped > 0` the `assignment_cycle` / `refresh_budget` accept routes do NOT re-check the clamp count, so the returned
  `converged=True` M0.84 state **carries 3 clamped cells of 330k**. Only the strict `tol` route still demands live 0-limited/0-floored.
  State this whenever the M6 number is quoted.
⚠ **P9/G9.1 is CITED, NOT RE-TESTED.** P9/G9.1 records that permanently-**limited** cells block the N5 freeze machinery on the CONFORMING
path; our blocker at M6 medium is mostly **floored** cells — the same *precondition*, a different clamp.
`freeze_max_clamped` exists **only on the LS path** (`newton.py` still has the hard 0-clamped rule), and whether relaxing it would unblock
G9.1's conforming fine mesh is an **UNTESTED HYPOTHESIS**, not a result. Do not cite B15 as having revived G9.1.

**Bit-identity:** `freeze_tol=None` (default) + `tol_residual_loose/rel=None` +
`accept_on_stall=False` ⇒ the pre-B15 live solver, byte-identical (locked).
Tests `tests/test_b15_ls_newton_freeze.py` (12 — 11 at closure + the 2026-07-15
errata lock `test_freeze_max_clamped_relaxes_the_convergence_semantics`).

Working rules (DN1 §9–§10):

- **No big-bang rewrite.** `solve/picard_ls.py` lives alongside
  `solve_subsonic_lifting`; the suite runs both paths parameterized; the default
  flips per-phase only after that phase's gate.
- **Dual-mesh testing (NEW 2026-07-11, user-directed;
  design_track_b.md §5.7).** Every B1–B6 gate runs on BOTH mesh types:
  (a) the existing wake-embedded meshes (M0/M1 — the "C-grid" analogue: nodes
  lie exactly on the wake plane, exercising the ε side-shift at scale, and
  enabling strict same-mesh A/B against the conforming path), and (b) the
  wake-free meshes (the "O-grid" analogue: generic cuts through generic
  elements — the actual workflow target): **M3** quasi-2D and **M4** ONERA M6.
  Where no conforming counterpart exists, acceptance compares against the (a)
  results and external references. The 3D pair (M1/M4) is not optional
  cosmetics: the swept TE and the wing tip carry machinery — oblique span
  frame, spanwise clip (Γ(tip)=0) — that the quasi-2D meshes cannot exercise
  at all, and B1 found a real defect in exactly that machinery.
- **Sequencing guard vs P8 (recorded 2026-07-10; wording updated 2026-07-11).**
  The P8 fully-coupled Newton is designed on the *conforming* wake (the
  Γ-Jacobian blocks come from `wake.py::self._h`; design.md §8.1), while B3's
  implicit Kutta removes the Γ DOF entirely. Land P8 on the conforming path
  first (done — P8 closed 2026-07-11); a level-set Newton is a post-B6
  re-derivation, not a parallel design (design_track_b.md §5.5: the wake-LS
  Jacobian blocks are constant in φ, no Γ elimination/Woodbury needed) —
  Track B blocks nothing in P7–P12.

---


### B16 — LS Newton far-field BC generalisation (far-field aux-DOF pin) ✓ CLOSED 2026-07-17 (NEW, user-directed; appended after B15, no renumber; executes the B9 recorded follow-up)

**Why (the B9 recorded follow-up, now measured — not prose).** B9 shipped its
LS leg on Picard because the LS Newton churns on the wing-body (medium Picard
**1458.9 s** vs the conforming Newton's 52.4 s), and recorded the diagnostic —
"`neumann` res→1e43, freestream Newton 8 far-field fluid rows |R|≈84" — as
prose only ("Findings recorded in memory", commit 555cfd8), which violates
session-discipline #3. B16 reproduces it as a committed artifact AND fixes it.

**Root cause (GB16.1, measured; corrects the pre-registered proposal).** A wake
level set has **no outflow clip** (`cut_elements.py`: a crossing need only be
downstream of the TE and within the span), so the sheet reaches the far-field
boundary and the outer nodes it crosses each carry an aux DOF governed only by a
**near-singular wake-LS row on a giant outer tet**. At the converged freestream
Picard state those aux hold garbage (coarse wing-body: |jump| **53.4** at x≥10
vs the physical Γ̄ **0.0586**); the Picard fixed point tolerates it (it solves
those rows to zero garbage-and-all), but the Newton residual reads it as an O(1)
inconsistency — **exactly the 8 far-field MAIN rows, max|R| = 84.457** (reproduced
to the digit). The aux-block conditioning is the single-number tell:
`jaa_diagnostic` cond1 = **9.1e18** (legacy free-aux — ABOVE the GB14.1 1e14
ceiling, i.e. genuinely singular; ⚠ doc erratum 2026-07-19: the prose's
6.36e18 was a pre-CSV trial value — the committed CSV reads 9.1e18 in both
epochs; onenormest estimates run-vary, the 12-order drop is the point) → **8.70e6**
(pinned). ⚠ **The proposal's
mechanism was wrong**: it attributed Picard's success to `closure="continuity"`
(weld) vs Newton's `wake_ls`, but the lifting/transonic Picard uses `wake_ls`
too (weld is only the Laplace seed) — both paths carry the same aux rows; the
difference is that Picard's fixed point absorbs the near-singular rows and
Newton's residual does not.

**Fix = `farfield_aux="pin"` (default) on `solve_multivalued_newton`, mode-adaptive
(user-arbitrated).** On a Dirichlet far field (freestream/vortex, where the
far-field MAIN DOFs are already Dirichlet) the far-field-BOUNDARY aux enter the
Dirichlet set at the branch value their host carries: freestream → the same
single-valued φ∞ (jump→0, consistent with the freestream BC already suppressing
circulation at a 25-MAC boundary); vortex → `main − side·γ` refreshed per step
(jump→γ, the conforming `lower_branch_mask` analogue, side from
`cm.node_side ∈ {±1}` — distinct from the B3 negative result, which pinned BOTH
sides to one branch and drained the circulation). **`neumann` is byte-identical
either way** (its outer aux ARE constrained by the wake-LS rows there, and every
committed LS-Newton anchor uses neumann), so the default flip is vacuous on all
committed evidence. New helper `farfield_aux_dofs(mesh, cm)` (`solve/picard_ls.py`).

**Gates**

- [x] **GB16.1 — diagnostic (evidence-debt repayment).** Committed script
  reproduces the churn on the coarse wing-body Picard state: max|R[free]| =
  **84.457**, 8 outer MAIN rows |R|>1 in the x∈[7,13] wake corridor, junk aux
  |jump| **53.4** at x≥10 vs Γ̄ 0.0586, cond1 **9.1e18 → 8.70e6** under the pin
  (doc erratum 2026-07-19: 6.36e18 was a pre-CSV trial value; the committed CSV reads 9.1e18 in both epochs).
  D8: the pin drives the outer jumps to **5.3e-15** — and ★ cures the 4 INTERIOR
  junk aux too (their wake-LS rows now anchor to clean Dirichlet data), so the
  R2 "interior aux unpinned" risk did not materialize.
- [x] **GB16.2 — wing-alone / neumann not perturbed.** naca coarse M0.7 neumann
  is `np.array_equal(phi_ext)` pin vs legacy (test); the B12 anchor
  (`onera_m6/medium` M0.5 neumann γ **0.06685284**) and the B15 ramp anchor
  (`onera_m6_wakefree/medium` M0.84 γ **0.088338**, M_max 2.4938, the
  freeze_max_clamped=8 / 3-clamped caveat unchanged) are the gated tier.
  (★ The B15 anchor VALUES are pre-B20; the re-baseline superseded them —
  see the GB15.4 erratum and GB20.7.)
- [x] **GB16.3 — wing-body freestream Newton converges (coarse ✓; medium
  gated).** Coarse M0.5: legacy churns (res **7.95**, 3690 limited) → pin reaches
  **res 5.88e-14, 0 limited** with the outer jumps at 5.3e-15. The BC layer is
  the fix. ⚠ **HONEST LIMIT (pre-registered branch, band NOT moved):** the pin
  state carries `n_flr=3`, so the strict `converged` flag (which needs 0 floored)
  does not fire — but D5 places those 3 cells at the wing-fuselage junction as
  the **B8 mixed-plain junk / G1.6 fuselage-Cp** class (M²_side 7.32 vs M²_main
  0.29; max honest M²_main 1.273 at the junction), the SAME root as **GB9.4**'s
  fuselage-lift xfail — a pre-existing issue orthogonal to the far-field BC fix,
  not created by B16. Recorded, not chased (G1.6 fix routes are closed negatives).
  ★ **And the coarse machine-converged lift is RIGHT:** cl_p(wing) **0.2086**
  vs the conforming **0.2089** (0.1%) — a properly converged Newton-pin agrees
  with the independent conforming path. (The medium is a different story — see
  GB16.4 and the ★★ block: the medium Newton STALLS and its lift is 22% low.)
- [ ] **GB16.4 — Newton-pin vs Picard vs conforming lift (XFAIL — the OPEN
  non-convergence).** The lift triangle does NOT close consistently, and the
  alignment FLIPS with resolution (all cl_p, wing): **coarse** the
  machine-converged Newton-pin **0.2086** matches conforming **0.2089** (0.1%)
  while LS Picard 0.1853 is the low outlier; **medium** the Newton-pin STALLS at
  res 7e-6 with **0.1690** — 22% BELOW both LS Picard **0.2165** and conforming
  **0.2173**, which agree (B9's 0.4% headline). ⇒ at least one path is NOT
  converged (see the ★★ block below). XFAIL: this is the open state, not a pass.
- [x] **GB16.5 — Schur/B14 compatible.** `SchurReducedSystem(...,
  n_aux_expected=mvop.n_ext − ff_aux_dofs.size)`: the pinned aux leave a
  contiguous free-aux tail; the split constructs and still fails loudly on the
  un-adjusted (legacy) count. `pytest tests/test_b14_schur_ls.py` stays green.
- [~] **GB16.6 — transonic wing-body stretch (RECORDED, gated).** wingbody
  freestream + a B15-style ramp to M0.84; records target_reached / death level /
  whether the death is the BC layer or the shock — no pass/fail. ★ **B18 erratum
  (2026-07-18): GB16.6 was spec'd "RECORDED" but NEVER implemented** (no demo
  code, no checks row — an evidence debt, discipline #3). **B18 executes it** and
  the answer is a NEGATIVE: the LS ramp dies at the wing-fuselage junction
  (coarse ~M0.575, medium at the first level ~M0.5, Mmax artifact growing with
  refinement — the G1.6/GB9.4 class), NOT at the BC layer or a wing shock. See B18.

★★ **RESOLVED BY B17 (2026-07-18) — it was NOT a non-convergence, it was a
BC-modelling error in the freestream pin.** See the **B17** section below. The
short version: `farfield_aux="pin"` on a freestream far field forces the outflow
wake jump to **0**, which REMOVES the circulation the wake physically carries out
— a resolution-dependent lift error, not a solver stall. Proof it is not a stall:
an independent **Picard** pin converges cleanly (res 7.5e-8) to the SAME medium
0.1691 the Newton-pin "stalls" at (0.1690) — both solvers agree per-BC. The coarse
"match to conforming" was a coincidence (jump=0 there cancelled the coarse legacy's
outer-tet garbage). The fix is `farfield_aux="pin_gamma"` (jump→γ): coarse
0.2087 / medium 0.2117 (Picard) and 0.2114 (Newton; 0.2115 pre-B20),
**monotone** to conforming, both solvers agreeing to 0.1%. The junction churn
that limited the medium residual floor at close-out was the B20 mixed-plain
contamination — post-B20 the same trajectory converges to |R| ~1e-13 (see the
★ note below the gates). The original (now-superseded) diagnosis is kept below
for the record.

★★ ~~UNRESOLVED NON-CONVERGENCE~~ (user-flagged 2026-07-18 — superseded by B17).
The far-field aux pin definitively fixes the **churn** — a self-contained result
on its own evidence (cond1 O(1e19)→8.7e6, res 84→**machine at coarse**, 0
limited vs legacy's 3690, and the coarse converged lift matches conforming to
0.1%). But comparing the LIFT across the three paths exposes a discrepancy the
pin does NOT resolve, and it FLIPS with resolution:
- **coarse:** the machine-converged Newton-pin (0.2086) ≈ conforming (0.2089);
  LS Picard (0.1853) is the low outlier.
- **medium:** LS Picard (0.2165) ≈ conforming (0.2173) — B9's headline — and the
  Newton-pin (0.1690, STALLED at res 7e-6, not machine) is the low outlier, 22%
  below.
The {Newton-pin, Picard, conforming} triangle therefore does **not** close
consistently ⇒ **at least one path is not converged.** Two live possibilities,
neither ruled out: **(a)** the medium Newton-pin is simply non-converged (it
stalls at 7e-6; a warm start from the *converged* B9 Picard state also failed to
converge within ~10 min, so it is NOT merely a shallow cold seed); **(b)** the
B9 medium LS-Picard≈conforming 0.4% agreement could itself be a non-converged
coincidence (both stopping near, but not at, a truly converged state). The
coarse evidence (converged Newton-pin ≈ conforming) favours (a), but does not
settle it. **UNRESOLVED — the medium-convergence / lift-consistency problem is
the open B16 follow-up; analysis deferred (user).** B16's churn fix stands; the
"the LS Newton now matches the other paths" claim does **not** — do not make it.

★ **New knob `farfield_aux` (default `"pin"`).** Defensible-as-default:
freestream Newton was NEVER exercised before B9 (the `ff_vals=None` bug made it
non-finite until 695baa0), vortex Newton has zero committed recipes, and neumann
is byte-identical — so "default leaves every committed anchor bit-identical"
holds vacuously. `"legacy"` is kept as the GB16.1 pathology reproduction switch.

**Bit-identity:** `farfield="neumann"` and `farfield_aux="legacy"` are
byte-identical to pre-B16 (locked: `test_neumann_bit_identical_pin_vs_legacy`,
and the b12/b13/b14/b15 neumann suites all pass unchanged). `schur_ls.py` changed
docstring/error text only (the assert body is preserved).
Tests `tests/test_b16_farfield_aux.py` (9, incl. gated GB16.3).

Working rules (DN1 §9–§10): the no-big-bang / dual-mesh / P8-sequencing block
above (B15) applies unchanged — B16 touches only the LS Newton far-field wiring
and one pure-query helper; the conforming path and the Picard drivers are
byte-untouched.

★ **B20 re-baseline number shifts (2026-07-19; the regenerated artifacts are
authoritative, the prose above keeps the pre-B20 values for the record).**
Coarse: legacy churn res 7.95 → **5.25** / 3690 limited → **11**; pin res
5.88e-14 → **1.49e-13**, floored 3 → **0** (those 3 were recorded above as
"B8/G1.6 class, not chased" — they were the mixed-plain contamination); outer
jumps 5.3e-15 → **1.8e-15**; max|R| 84.457 → **84.8** (so "reproduced to the
digit" now reads "to 3 figures"); neumann-unbounded diagnostic 2.6e8 → 5.8e5.
Medium: pin res 7.03e-6 → **1.79e-13** (the pin state now CONVERGES — the
churn that limited it was the contamination, GB20.5). Every gate verdict is
unchanged; GB16.4's XFAIL row in the artifact is superseded by B17.

---

### B17 — Far-field aux pin carries jump=γ, not 0 (resolves GB16.4) ✓ CLOSED 2026-07-18 (NEW, user-directed; appended after B16, no renumber; executes the B16 GB16.4 open follow-up)

**GB16.4 was NOT a non-convergence — it was a boundary-condition modelling error
in the B16 freestream pin.** The B16 pin, on `farfield="freestream"`, forced the
outflow wake potential-jump to **0**. Physically the wake carries its jump
[φ]=Γ out to the boundary; zeroing it **removes the outflow circulation**, a
resolution-dependent lift error. This was invisible at coarse (the jump=0 error
happened to cancel the coarse legacy's near-singular outer-tet garbage) and
un-masked at medium (where the legacy wake-LS already carries the jump correctly,
so the pin's error dominates: −22%).

**Decisive discriminator (E2): both solvers agree per-BC.** Giving the *Picard*
driver the same freestream pin (new `farfield_aux` knob on
`solve_multivalued_lifting`) makes medium Picard-pin converge cleanly (res 7.5e-8,
34 outers) to **cl_p 0.1691** — matching the "stalled" Newton-pin **0.1690** to
0.1%. Two independent solvers landing on the same value ⇒ it is a genuine
(BC-determined) state, **not** a Newton stall. So B16's possibility (a) was wrong;
the pin fixed point itself is wrong.

**The fix — `farfield_aux="pin_gamma"` (jump→γ), the new default.** aux = host
φ∞ − side·γ, refreshed with the live γ (Picard: per outer; Newton: per step) — the
same near-singular-aux Dirichlet cure B16 correctly identified, but with the
physical ring value. The triangle then closes and is **monotone** to conforming:

| far-field aux | coarse cl_p | medium cl_p | trend |
| --- | --- | --- | --- |
| conforming (P14 Newton, ref) | 0.2089 | 0.2173 | ↑ |
| legacy (free aux) | 0.1853 | 0.2165 | coarse polluted by \|jump\|=53 outer tet |
| pin jump=0 (B16) | 0.2086 | 0.1690 | ✗ non-monotone, kills outflow circ |
| **pin_gamma (B17)** | **0.2087** | **0.2117** (Picard) / **0.2114** (Newton; 0.2115 pre-B20) | ✓ monotone |

Newton-pin_gamma and Picard-pin_gamma agree to 0.1% at both resolutions; both
undershoot conforming by 0.1%/2.6% (far-field truncation). At B17 close-out the
medium Newton-pin_gamma still carried the wing-fuselage-junction churn (nlim 42
/ nflr 40, res 5.5e-5) — but the lift was correct regardless (γ stable 0.06420,
cl_p 0.2115). ★ **Post-B20 update (2026-07-19): that churn WAS the mixed-plain
contamination** — the re-baselined trajectory converges to **|R| ~1e-13, γ
0.064201, cl_p 0.2114** with the clamps gone (GB20.5); the "only limits the
residual floor" framing is superseded. ★ **B16 conflated two orthogonal issues:** the
far-field near-singular *conditioning* (which the pin cures, jump value
irrelevant) and the outflow *circulation* (which needs jump=γ). The junction churn
is a third, pre-existing issue that only limits the residual floor, not the lift.

- [x] **GB17.1 — coarse 4-way triangle + ring-jump collapse.** pin_gamma 0.2087
  ≈ conforming 0.2089 (0.1%); the ring jump collapses legacy **53.4** → pin **0** →
  pin_gamma **0.063** = γ; and the coarse **legacy** garbage (\|jump\|=53) is a 12%
  lift deficit ⇒ B9's "coarse 12.8% = resolution" was largely far-field
  **contamination** (see the B9 erratum).
- [x] **GB17.2 — post-processing NOT the cause.** cl_p (surface-pressure integral)
  and cl_KJ (circulation integral) move together (~22% both at medium) ⇒ the pin
  gap is a genuine flow-state change, not a post artifact. The user's suspicion
  (Cp "looks aligned" yet cl_p differs) is a Cp-axis scale illusion: the plotted
  sectional cl(z) is the **Γ-based** `2Γ/(u·c)`, and per-station ∫Cp differs
  24–44% while the Cp curves differ only ~0.03–0.05 on a ±1 axis.
- [x] **GB17.3 — pin jump=0 is a BC error, both solvers.** medium Picard-pin-jump0
  0.169 ≡ B16 Newton-pin-jump0 0.1690 ⇒ GB16.4 is a BC-modelling error, not a stall.
- [x] **GB17.4 — pin_gamma closes the triangle.** Newton 0.2114 (0.2115
  pre-B20) ≈ Picard 0.2117 (<1%), monotone coarse→medium toward conforming.
- [x] **GB17.5 — spanwise Γ(z) uniform offset removed (RECORDED).** the jump=0
  spanwise Γ is a uniform multiplicative deficit (~22% at every station); pin_gamma
  restores it.
- [x] **GB17.6 — vortex evaluation (RECORDED, user-requested).** `farfield="vortex"`
  does NOT close the residual gap — it **brackets** conforming from the other side
  (medium **+2.5%** vs pin_gamma's −2.6%) and its free far-field aux **churn at
  coarse** (res 3.2, \|jump\|=71, needs its own pin). The 2–3% is far-field
  truncation; **freestream pin_gamma stays recommended** (clean at both resolutions).

**Defaults (user-arbitrated 2026-07-18):** `farfield_aux="pin_gamma"` is the new
default on BOTH `solve_multivalued_newton` (was `"pin"`) and
`solve_multivalued_lifting` (was free/legacy). Safety: it acts **only** on
`farfield="freestream"`, and is inert (bit-identical to legacy) on vortex/neumann
— so every committed 2.5D NACA vortex/neumann Picard run and every neumann Newton
anchor is byte-untouched. The B9/B16 **freestream** Picard demos were pinned to
explicit `farfield_aux="legacy"` to keep their committed numbers reproducible; a
**B9 erratum** records that coarse 12.8% was contamination (its medium
legacy≈conforming headline stands — legacy happens to carry the jump correctly at
medium). B16's committed jump=0 numbers reproduce with explicit `farfield_aux="pin"`.
`"pin"` (jump=0) is kept as the diagnostic value.

Tests `tests/test_b17_farfield_pin_gamma.py` (6, ungated); demo
`cases/demo/b17_farfield_pin_gamma/` (Part 1–2 coarse ungated, Part 3 medium gated).
The B16 tests are unchanged except `test_farfield_aux_knob` (default is now
`pin_gamma`).

---

### B18 — Wing-body transonic (M0.84): conforming reaches it, level-set is junction-limited ✓ CLOSED 2026-07-18 (NEW, user-directed; appended after B17, no renumber; executes the GB16.6 debt)

**The wing-body transonic capability is asymmetric, and that is the finding.**
Subsonic M0.5 wing-body is done (B9/B17). Pushing the Mach up:

- **Conforming (Newton + pressure Kutta, Mach continuation) IS the wing-body
  transonic path.** coarse reaches **M0.84 (cl_p 0.2617**, Mmax 2.15, strict
  res 2.8e-12); medium reaches **M0.79 STRICT (cl_p 0.2579**, res 2.2e-14), with
  a clean transonic rise cl_p(M) = **0.2173 / 0.2321 / 0.2579** at M0.50/0.65/0.79.
  Medium **M0.80+ stalls** (res ~2–7e-6, 0 clamp) — NOT slivers (the medium mesh
  is clean, min_dihedral 9.75°, 0 tets < 5°; the coarse mesh has 27 slivers yet
  reaches 0.84), a sharper shock/junction interaction; recorded, not chased.
- **Level-set (B15 freeze-ramp + B17 pin_gamma) does NOT reach transonic on the
  wing-body.** The wing-fuselage junction carries a spurious supersonic pocket
  (M²≈1.27 already at M0.5) that **WORSENS with refinement**: at close-out,
  coarse ceiling **M0.575** (Mmax 1.44), medium dies at the FIRST transonic
  level ~**M0.5** (Mmax artifact 3.96, nlim 43/nflr 40). This is a
  **closed-negative discretization error** (discipline #8), characterized
  not chased — and it is the direct analogue of GB9.4's "LS fuselage lift GROWS
  0.164→0.205 with refinement."
  ★★ **ATTRIBUTION CORRECTED + numbers re-based by B20/GB20.5 (erratum
  2026-07-19).** B18's original guess "**B8 mixed-plain** class" is **measured
  FALSE**: removing the mixed-plain contamination made the medium case CONVERGE
  (res 6.8e-5 → 1.1e-13, clamps 82 → 6) and **UNCLAMPED the pocket, revealing a
  GENUINE Mmax 5.22** at subsonic freestream — i.e. the pocket is the
  **G1.6/GB9.4 faceted-geometry** error, not the mixed-plain density. Post-B20
  re-baselined numbers: coarse ceiling **~M0.55 (Mmax 1.31, m_last_conv=nan)**;
  medium dies ~M0.5 with **Mmax 5.22 genuine (nlim 3/nflr 3)**. The old
  Mmax 3.96 was a CLAMPED non-converged artifact. Conforming legs untouched.
  ★★★ **POCKET CURED + CEILING RE-MEASURED by B23–B27 (erratum 2026-07-20).**
  The "junction-limited closed-negative" headline is **retired**: B23 attributed
  the pocket to the wake inboard FREE-EDGE singularity (lift/wake-coupled, not
  faceted geometry — α=0 clean at both levels), B24 closed the waterline-extension
  route, B25's `inboard_clip` healed the pocket (medium α=3.06 corrM 14.66 →
  0.63, guardrails clean), and B26 re-measured the post-cure ceiling: **medium
  0.50 → 0.7625, coarse 0.82 → 0.84 reached**; the death cause flipped to
  (b)-class wing-tip P13 + high-M Newton stall (the same class as the conforming
  medium 0.80+ stall), no longer the junction. B27 refreshed the demo to 8/8
  checks PASS and extended the cross-model table (below). What survives of the
  original B18: the conforming legs (untouched, bit-reproduced by B27) and the
  GB9.4/G1.6 fuselage-lift characterization (cl_fus persists — now the primary
  G1.6 suspect).

Consequence (2026-07-18; **superseded 2026-07-20 by B27**): there WAS **no
common transonic Mach at medium** (LS could not leave 0.5), so the trustworthy
cross-model check STAYED **M0.5 (B9/B17: 2.6%; post-B20 LS value 0.2117 vs conf
0.2173)**. ★ B27 (2026-07-20): post-cure the cross-model gap is **flat across
Mach — M0.5 2.6 % / M0.65 2.4 % (conf 0.2321 vs LS+clip 0.2266, PASS ≤ 5 %) /
M0.75 2.5 % (conf 0.2483 vs LS+clip 0.2421, RECORDED no-threshold)** — one
~2.5 % band = the known B17 LS cl_p↔cl_kj convention gap, not physics; the
M0.5-only restriction is retired. ★ Erratum 2026-07-19: the coarse M0.60
transonic cross-model, originally SKIPPED, **exists in the re-baselined
artifact** (`cross_model.csv` row 3: conf 0.2178 vs LS 0.2174 = **0.2 %**,
increment comparison 2 %, under-resolved; the LS state is `ls_coarse_06.npz`
— quote it with the coarse strict-ceiling caveat above, not as a converged
capability claim). ★ **Retired 2026-07-20 (B27 §3.1)**: the cured C-side value
is 0.2133 (gap **2.1 %**) — the old A-side 0.2 % was a pocket-contaminated
coincidence landing near conforming; the row is superseded.

- [x] **GB18.1 — conforming transonic (PASS).** coarse M0.84 0.2617
  (proof-of-concept, under-resolved); medium M0.79 0.2579 strict, monotone cl(M)
  rise. The wing-body transonic deliverable.
- [~] **GB18.2 — LS transonic ceiling (RECORDED; superseded by B26).** coarse ~M0.575, medium ~M0.5
  at close-out (post-B20 re-baseline: coarse **~M0.55**, medium ~M0.5); the
  GB16.6 debt repaid as a negative (junction, not BC layer). **B26 (2026-07-20)
  re-answers POSITIVE: post-cure ceiling medium 0.7625 / coarse 0.84 reached.**
- [~] **GB18.3 — cross-model (RECORDED; superseded by B27).** M0.5 medium (2.6%) is the only
  trustworthy cross-model; the medium transonic cross-model is BLOCKED by the LS
  junction, and that is the finding. **B27 (2026-07-20): the transonic
  cross-model now exists — M0.65 2.4 % PASS, M0.75 2.5 % RECORDED.**
- [~] **GB18.4 — junction transonic characterization (RECORDED; re-answered by B27/GB27.4).** the spurious
  pocket grows with refinement (close-out: coarse Mmax 1.4 → medium 4.0;
  post-B20: coarse **1.31** → medium **5.22 genuine, unclamped**), GB9.4 sign;
  attribution corrected by GB20.5 (G1.6 geometry, NOT mixed-plain — see the
  erratum block above). **Final attribution (B23, 2026-07-19): the wake inboard
  free-edge singularity — cured C-side by B25; the residual ceiling limiter is
  the wing-tip P13 class + high-M Newton, per B26 committed peaks.**
- [~] **GB18.5 — fuselage lift at the medium transonic top (RECORDED).** cl_fus
  16% of wing cl_p at M0.79 — the G1.6 flat-facet natural-BC error persists into
  transonic (GB9.4 class).

★ **No `pyfp3d/` numerics change** — B18 is a pure demo/tests/docs phase using the
existing conforming `solve_newton_transonic` and LS `solve_multivalued_newton_transonic`
(the B15 ramp with B17 pin_gamma). ★ **Recipe note:** the conforming wing-body
medium ramp needs `freeze_tol` raised to the wing-body churn floor (1e-6 → 1e-5,
the B17 lesson) or it stalls at M0.80 with the wing-alone recipe. ★ **fine
excluded** (G13.3 negative + LS has no fine escape). Tests
`tests/test_b18_wingbody_transonic.py` (4, ungated); demo
`cases/demo/b18_wingbody_transonic/` (7 gated gates at the 2026-07-18 close-out:
1 PASS + 6 RECORDED — **superseded 2026-07-20 by the B27 refresh: checks.csv
8/8 PASS**, see the ★★★ erratum block above).

---


### B19 — LS Newton Jacobian exactness on mixed-side plain elements ✓ CLOSED 2026-07-18 (NEW, user-directed; appended after B18, no renumber; executes the A3/GA3.6 C1 finding)

**Trigger.** A3 verified the 2026-07-17 inspection's C1 (evidence:
`cases/analysis/c1_ls_jacobian_fd/`, demo_report/track_a.md §A3): on 3-D meshes
the LS Newton Jacobian is **not** the derivative of its residual on
**mixed-side plain elements** — targeted probe ‖Jv−FD‖/‖FD‖ **1.146e-01** vs
control **6.33e-10**, and **ε-independent** (1.532e-01 at ε 1e-6/1e-7/1e-8,
max/min 1.00) ⇒ a missing term, not FD noise. 3378 such elements on M6 coarse.
R is untouched by a Jacobian error, so every converged state and gate number
stands; what degrades is the convergence RATE (a quasi-Newton in 3-D).

**★ The phase is TWO problems, deliberately separated (user-arbitrated
2026-07-18). Do not merge the legs — merging them makes any result
unattributable.**

**Root cause (re-derived from the code, not transcribed).** One element has
TWO different DOF maps, and Terms 2/3 use one map for both roles:

- **Row map** — where the element's residual LANDS. `mass_conservation_coo`
  scatters plain elements onto **main** DOFs (`cut_assembly.py:173`), cut
  elements onto `dofs_upper`/`dofs_lower`.
- **Column map** — which DOFs the element's side FIELD READS.
  `side_potentials` (`multivalued.py:381-382`) is a **per-NODE** rule: a cut
  node takes the **aux** value on the opposite side. `newton_side_data`
  computes the gradient from that side field (`multivalued.py:562`
  `grad, q2 = self.op.velocities(phi_s)`), so ρ — and hence the residual —
  really does depend on those aux DOFs.

`_side_dofvecs` (`multivalued.py:509,513`) applies its override only to
`cm.cut_elems`; plain elements keep main connectivity on both sides. So
`newton_terms23_side_coo` (`cut_assembly.py:441-456`) scatters the aux
sensitivity onto a **main** column. **For cut elements the two maps coincide**
— which is exactly why the 2.5-D FD gate passes and why this survived to B18.

#### Leg A — make J exact for the R we have (low risk, do FIRST)

Split the two roles: rows keep `dofvec` (the mass-conservation scatter),
columns become a new per-side **read map** mirroring `side_potentials`
node-by-node, for ALL elements:

    read_u = where(node_side[el] == +1, el, ext_dof_of_node[el])
    read_l = where(node_side[el] == -1, el, ext_dof_of_node[el])
    # ext_dof_of_node < 0 (non-cut node) falls back to main

Term 2: rows `dofvec[e]`, cols `read[e]`. Term 3: rows `dofvec[e]`, cols
`read[u(e)]`. ★ **This is not a new rule — on cut elements `read_u`/`read_l`
reproduce `dofs_upper`/`dofs_lower` exactly** (assert it), so the fix
generalizes the existing per-element special case into the per-node rule it
should always have been. That equivalence is the strongest evidence the fix is
right, and it is gated below.

★★ **The index split is only HALF of it — a SECOND missing term, found by
measurement (2026-07-18).** After the row/column index split the targeted probe
still read 1.4697e-02 (from 1.146e-01) and was still ε-independent. A
column-wise FD localization, then a block isolation that rebuilt Terms 2 and 3
separately, showed `|FD23 − J23| ≡ |FD23 − J2|` **exactly** on most touched aux
columns ⇒ **Term 3 contributes nothing there** (those elements are subsonic,
s_u = 0) ⇒ the残 error is in **Term 2**, not the upstream coupling the first
row-classification had suggested. ★ *This is why the localization was run
instead of reasoning from the row/element classification: the first reading of
that classification pointed at Term 3 and was wrong.*

The mechanism, and it is the same duality one level down — **not just the DOF
indices but the GRADIENT FACTORS come from two different fields**. An element's
residual contribution is

    R_a(e) = rho_tilde( grad of the READ field ) * V_e * ( grad of the
             SCATTER field . B_a )

— the density is built from the side field, but the `K @ phi` that it weights
is contracted with the field the SCATTER map reads (main, for a plain element).
Terms 2/3 used the side gradient for **both** factors. So the fix is: the ROW
factor uses `grad_row` = the gradient over `phi_ext[dofvec]`, the COLUMN factor
keeps `grad` = the side field's. On cut elements `phi_ext[dofs_upper]` IS
`phi_up` over that element, so the two gradients coincide and nothing moves —
the same reason the index split was inert there.

- [x] **GB19.1 — the C1 measurement inverts ✓ PASS** (`results/c1_fd_probes.csv`
  post-fix, `c1_fd_probes_prefix.csv` pre-fix, `b19_three_states.csv` = the
  three-state ledger). Targeted probe **1.145684e-01 → 1.333699e-08**; control
  6.327479e-10 unchanged; global-free **2.47e-03 → 8.49e-10**.
  ★★ **The discriminator FLIPPED, and that is the real proof.** The same ε
  sweep that convicted the code now acquits it: pre-fix **1.532e-01 at every ε**
  (spread 1.00 = a missing term); post-fix **1.6e-09 / 2.1e-08 / 2.2e-07**
  (spread 131.5, scaling like 1/ε = pure FD roundoff, the numerical floor). An
  ε-independent error became an ε-sensitive one — exactly the transition a real
  fix must produce, and one no amount of tolerance-loosening could fake.
  ★ **The intermediate state is recorded, not erased:** after the index split
  ALONE the probe read **1.4697e-02** — 8× better, still ε-independent
  (2.913917e-02 at all three ε) ⇒ held open as PARTIAL rather than rounded into
  a pass, which is what forced the second defect to be found.
  ★ **Ruled out along the way:** `_side_element_sets` uses the *same*
  `node_side.max() == 1` rule as `mass_conservation_coo`, so the side-set
  assignment was never a second inconsistency.
- [x] **GB19.2 — R is BIT-IDENTICAL (the load-bearing claim).** On a fixed
  φ_ext, the assembled residual before and after must agree bit for bit on a
  3-D mesh. This is what licenses "no converged answer moves"; assert it, do
  not infer it.
- [x] **GB19.3 — `read_*` ≡ `dofs_upper`/`dofs_lower` on cut elements**, and
  the 2.5-D Jacobian is bit-identical (quasi-2D has no mixed-side plain
  elements ⇒ J itself must not move). Full ungated suite unchanged.
- [x] **GB19.4 — the 3-D convergence delta, recorded honestly (A/B COMPLETE —
  recorded NEGATIVE).** M6 coarse M0.5→0.70 freestream LS Newton,
  n_seed 30 / n_newton_max 40 / tol 1e-8: **does NOT reach tol** — it settles
  into a clean PERIODIC plateau at **~3.5e-07** (residual cycling 3.31e-07 …
  4.12e-07 with an obvious period), **0 limited / 0 floored**, γ 0.07212068,
  M_max 1.134, 62.8 s.
  ★ **That plateau is the B15 upwind selection/branch churn limit cycle, not a
  Jacobian defect** — with 0 limited/0 floored and a periodic residual, the
  iteration is bouncing between assignments, which is exactly the object B15's
  freeze machinery was built to cure (this run deliberately used no freeze, to
  observe the raw Newton). An exact Jacobian cannot fix a discontinuous
  selection: **the fix removes the wrong-derivative error, it does not remove
  the churn floor, and those are different problems.** Recording this so the
  phase does not get credited with a convergence improvement it did not make.
  ★★ **A/B COMPLETE — the exact Jacobian changes NOTHING measurable here, and
  that is the finding.** Same case, same driver, stash A/B:

  | | steps | converged | residual | γ | M_max | wall |
  |---|---|---|---|---|---|---|
  | pre-fix (J wrong) | 40 | False | 3.707261e-07 | 0.07212068 | 1.134235 | 60.6 s |
  | post-fix (J exact) | 40 | False | 3.973387e-07 | **0.07212068** | **1.134235** | 62.8 s |

  **γ identical to 8 decimals, M_max to 6** (as GB19.2's bit-identical R
  requires), the plateau sits at the same level — the two residuals differ only
  by which phase of the limit cycle step 40 lands on — and the step count is
  unchanged. Cost of exactness: **+3.6 % wall** (the extra `grad_row` einsum
  per side per step).
  ⇒ **The phase must NOT be credited with a convergence improvement.** The
  binding constraint on this case is the selection churn, which an exact
  derivative cannot touch. What B19 Leg A buys is **correctness** — the driver
  is now a Newton method rather than a quasi-Newton, so its behaviour is
  analysable and its future improvement (freeze, better preconditioning) rests
  on a true Jacobian. Whether that pays off anywhere is unmeasured, and this
  gate deliberately does not claim it.
- [x] **GB19.5 — a 3-D FD gate enters the suite** (`tests/test_b19_jacobian_3d.py`;
  the heavy FD gate itself gated, two structural checks ungated), so this
  element class can never regress unseen again. ★ This gate exists because the
  ROOT process failure was a structural blind spot — fixing the bug without
  fixing the blind spot would leave the trap armed.
  ★★ **ERRATUM, measured 2026-07-18 — the blind spot was MIS-STATED, by the C1
  write-up and by this author's first version of the test.** The claim was
  "quasi-2-D has no mixed-side plain elements". **It has 129** (coarse). What it
  has **zero** of is the subset that can reach an aux DOF: **0 of those 129
  touch a cut node**. A mixed-side plain element only mis-scatters when its
  side field actually READS an aux value, i.e. when one of its nodes is a cut
  node — that is the real invariant, and it is what the test now asserts. The
  coarser claim was wrong in fact while accidentally right in conclusion; it is
  corrected everywhere rather than quietly left standing (the P14 lesson: do
  not carry an attribution forward as if it were measured).

#### Leg B — is the R we have the R we want? (model question, do SECOND)

Leg A makes Newton exact for the CURRENT residual. It does not ask whether
that residual is right. For a mixed-side plain element the current R takes its
**stiffness from the main-field φ** but its **density from the side field** —
an asymmetry that was never a deliberate modelling decision on record.

- [x] **GB19.6 — characterized ✓, and the asymmetry is NOT benign** (zero-change
  probe `run_legb_probe.py`, `results/legb_density_gap.csv`; M6 coarse M0.70,
  Picard-400 state, |R| 1.95e-07).

  Scope first: of 3378 mixed-side plain elements only **252 actually read an
  aux DOF** — **0.1888 % of domain volume**. So this is a thin strip, and a
  small gap there would have argued for leaving R alone.

  **It is not small:**

  | metric | value |
  |---|---|
  | max \|ρ_side − ρ_main\| | **0.4474** (**45.3 %** relative) |
  | mean relative gap | 4.51 % (volume-weighted 0.0260 absolute) |
  | elements differing > 1 % / > 10 % | **154 / 39** of 252 |
  | **max q² read by the SIDE field** | **3.2229** (M ≈ 1.80, at the M_cap=3 limiter) |
  | **max q² read by the MAIN field** | **1.3379** (M ≈ 1.16) |

  ★★ **The last two rows are the finding.** On these elements the side field
  produces a **spurious supersonic state** where the main field is barely
  transonic — because the side field substitutes AUX values, which hold the
  flow on the OTHER side of the wake, into an element that is not cut and whose
  stiffness is assembled on main DOFs. The artificial-density switch then fires
  on that fictitious q², and the solver's own density — not merely a diagnostic
  — carries it.

  ★ **Same class, third bite.** B8 measured exactly this contamination as a ×5
  metric artifact in `element_mach2` (`element_mach2` reads mixed-side plain
  elements' side field; the honest exponent was +0.62, the artifact +1.34); A3
  found it in the Jacobian; this finds it in the residual's density.

  ★ **HYPOTHESIS, explicitly not a conclusion — a possible link to the
  junction/tip pockets.** B18 recorded that the level-set wing-body carries a
  spurious supersonic pocket that **worsens with refinement** (M² ≈ 1.27 at
  M0.5 coarse, M_max artifact 3.96 at medium), closed-negative as a
  discretization error. Spurious side-field supersonic states in the
  mixed-side plain strip are a *candidate* contributing mechanism with the
  right qualitative signature (more elements straddle as h falls). **This was
  measured on the wing-alone M6, NOT on the wing-body, and no causal link is
  demonstrated.** The named test would be: rerun the B18 medium level-set case
  with the main-field density on this element class and see whether M_max 3.96
  moves. **Do not repeat this hypothesis as a result** (the P14 lesson).

  **Routing: NOT adopted here, by design.** Changing the density source changes
  R ⇒ it moves every converged level-set answer, so it is a discretization
  change requiring its own phase with its own before/after on the committed
  cases (γ, cl_p, M_max, shock position, and the B18 pocket). Leg B's
  deliverable is the measured case for opening it — which the 45 % / spurious-
  supersonic numbers make, and which the 0.19 % volume share does not weaken,
  because the affected strip sits exactly where the tip and junction pathologies
  live.
- ★ **Prior art, same element class:** B8 measured a ×5 metric artifact from
  `element_mach2` reading mixed-side plain elements' side field
  (`b8-respec-negative-and-metric-artifact`). This class has now bitten twice —
  once in a metric, once in the Jacobian — which is itself the argument for
  asking the model question rather than only patching J.

### B20 — Mixed-side plain elements: main-field density (the B19 Leg B fix) ✓ CLOSED 2026-07-18, ADOPTED PERMANENTLY (user-arbitrated; knob removed); re-baselined 2026-07-19 (NEW, user-directed; executes the B19/GB19.6 finding)

**Trigger.** B19 Leg B (GB19.6) measured that on **mixed-side plain elements**
the residual takes its stiffness from the main field but its **density from the
side field**, and that this is not benign: the side field manufactures a
**spurious supersonic state** (q² 3.2229, M≈1.80, at the M_cap limiter) where
the main field reads 1.3379 (M≈1.16), with max density error **45.3 %** on the
252 aux-reading elements. The contamination reaches the SOLVER's density, not
just a diagnostic.

**★ The reporting layer already made this call.** `element_mach2` has carried
`mixed_plain="main"` as its DEFAULT since 2026-07-14 (user-arbitrated B8
backlog flip): the DIAGNOSTIC reads the main field for exactly this class,
calling the side reading a "×5 inflation … manufactured the LS tip exponent
p=+1.34 (honest +0.62)". **B20 makes the ASSEMBLY consistent with the
diagnostic** — a plain (single-valued, uncut) element is single-valued, so its
density belongs to the single-valued main field. The side-field substitution
only ever made sense for cut elements (which are genuinely two-valued and
assembled twice).

**Architecture — ★★ PERMANENT AND NON-OPTIONAL (user-arbitrated 2026-07-18).**
The fix was first built as a default-off `plain_density` knob so it could be
measured against the old behaviour (that A/B is the evidence below). **On the
strength of that evidence the user directed that it be hard-coded: no switch,
no default, no fallback.** The knob is REMOVED; mixed-side plain elements read
the main field unconditionally, at all three density sites
(`element_rho_tilde`, `newton_side_data`, `element_densities`) and in Leg A's
`_side_readvecs`.

**Why it is not a preference.** The decisive argument needs no physics: the
element's stiffness is contracted with the MAIN field (it scatters onto main
DOFs) while its density came from a DIFFERENT field — **one equation built from
two velocity fields, internally inconsistent**, regardless of which branch one
argues is "physically right". It is also uncut, so no jump passes through it.
And the reporting layer had already ruled the same way: `element_mach2` has
defaulted to `mixed_plain="main"` since 2026-07-14. B20 makes the ASSEMBLY
agree with the DIAGNOSTIC.

**Consequence, accepted with the decision:** the pre-B20 3-D level-set numbers
were solutions to a discretization carrying a known internal inconsistency.
They are re-based here; the old values stay traceable in git and in the
before/after CSVs.

- [x] **GB20.1 ✓ localized** (`results/legb_apply.csv`). Quasi-2-D NACA:
  `"main"` ≡ `"side"` **bit-identical (max|ΔR| = 0.000e+00)** — its 129
  mixed-side plain elements have **0** cut nodes, so main == side there.
  M6 coarse: R moves on **164 rows** (of 12k), 4 of them just outside the
  mixed-plain node set (legitimate Term-3 upstream coupling). ★ **A workspace
  aliasing bug was caught and fixed here first:** `PicardOperator.velocities`
  returns VIEWS into a shared buffer, so recomputing the main-field gradient
  clobbered the caller's side values in place — the first run moved 2940
  elements (not 129) and tripled subsonic Γ. Detaching with `.copy()` before
  the second `velocities` call fixed it; the 2940→164 drop is the tell.
  ★ *Caught by measuring an unexpected result, not by rationalizing it.*
- [x] **GB20.2 ✓ the Jacobian stays exact under `"main"`** — targeted probe
  **8.07e-09**, control **6.29e-10** (`results/legb_apply.csv`). Leg A composes
  with Leg B: under main density the mixed-plain columns return to main
  (`_side_readvecs` mode-aware) and the row/column gradient factors coincide,
  so J = dR/dφ still holds. The naive first cut (before the aliasing fix) read
  1.97e-02 — that was the same buffer bug, not a missing term.
- [x] **GB20.3 ✓ subsonic 2.5-D is UNAFFECTED** (`results/legb_subsonic_ab.csv`).
  NACA medium M0.5 embedded Γ **0.088144 → 0.088144 (+0.0000 %)**, wake-free
  0.088348 unchanged, B3 wake-free-vs-embedded 0.2313 % unchanged. Expected and
  honest: the quasi-2-D mesh has no aux-reading mixed-side plain elements (the
  class is a 3-D tip/junction phenomenon), so Leg B is a genuine no-op on every
  committed 2.5-D case. The change lives only in 3-D.
- [x] **GB20.4 ✓ transonic before/after — main density REACHES the target where
  side stalls** (`results/legb_transonic_ab.csv`, M6 coarse wing-alone ramp to
  M0.84): **side m_final 0.7875 NOT converged** (γ 0.077959, M_max 1.3157) vs
  **main m_final 0.84 CONVERGED** (γ 0.084812, M_max 1.4491). A genuine
  capability move on the wing-alone M6 coarse: removing the spurious-supersonic
  side-density contamination let the transonic ramp climb the last two levels
  and converge at the full target. (M_max is higher in main only because it is
  at M0.84 vs M0.7875 — a higher freestream, not a worse pocket.) Recorded as a
  positive signal, wing-alone/coarse; the wing-body medium is GB20.5.
- [x] **GB20.5 ✓ THE HYPOTHESIS TEST — a SPLIT result, and it is the most
  informative one** (`results/legb_b18_hypothesis.csv`, B18 medium LS wing-body,
  the side leg reproduces the committed baseline: side m_final 0.5, Mmax 3.920,
  nlim 42 / nflr 40, res 6.8e-5 ≈ committed 3.964 / 43 / 40).

  | | m_final | converged (res) | Mmax | nlim / nflr |
  |---|---|---|---|---|
  | side (committed) | 0.5 | **NO (6.8e-5)** | 3.920 | **42 / 40** |
  | main (Leg B) | 0.5 | **YES (1.1e-13)** | 5.220 | **3 / 3** |

  ★★ **The naive read — "Mmax 3.92 → 5.22, +33 %, worse" — is WRONG.** The side
  state is NOT converged: it churns at res 6.8e-5 with **82 cells on the
  limiter/floor**, so its 3.92 is a CLAMPED, propped-up value. The main state
  **converges to machine precision (1.1e-13) with only 6 clamped cells** — a
  genuine discrete solution.

  **So the hypothesis SPLITS in two, and both halves are findings:**
  - ★ **The CONVERGENCE pathology WAS substantially the mixed-plain side-density
    artifact.** B18's "dies at M0.5, 42/40 clamped, churns" is largely CURED by
    main density: res 6.8e-5 → **1.1e-13**, clamps 82 → **6**. The junction
    churn/limit-cycle that made the level-set wing-body untractable is, to a
    large degree, the spurious side-supersonic feedback GB19.6 measured.
  - ★ **The junction POCKET itself is REAL, and B19's literal hypothesis is
    REFUTED.** Removing the contamination did not remove the pocket — it
    UNCLAMPED it, revealing a genuine converged M≈5.2 spike at the wing-body
    junction (a subsonic-freestream local pocket = the G1.6/GB9.4 faceted-
    geometry discretization error, NOT the mixed-plain density). Main still
    cannot advance past M0.5.

  **Bottom line:** Leg B is not the cure for the junction pocket (that is
  geometry → curved elements, G1.6), but it IS a large part of the cure for the
  level-set wing-body's non-convergence. RECORDED, not pass/fail, as
  pre-registered. **Do not repeat "the pocket is mixed-plain contamination" —
  measured false.**
- [x] **GB20.6 ✓ full suite unchanged + adoption dossier delivered.** Suite
  **465 passed + 22 skipped + 2 xfailed** (1112.37 s @16) — identical to the
  B19 baseline. (Wording updated post-permanence: at measurement time the fix
  sat behind a default-off `plain_density` knob, which is what "suite
  unchanged" verified; the knob was then REMOVED on adoption — see the section
  header — and the re-baseline below re-verified 465+22+2 with the fix
  hard-coded.)

  **Adoption dossier (as presented for the user's decision; the user then
  adopted `main` PERMANENTLY, 2026-07-18).**

  | evidence | side (today) | main (Leg B) |
  |---|---|---|
  | 2.5-D subsonic Γ (B3/B4) | 0.088144 | **0.088144 — no-op** |
  | M6 coarse ramp → M0.84 | m 0.7875, **not converged** | **M0.84, converged** |
  | B18 wing-body medium @M0.5 | res 6.8e-5, **82 clamped** | **res 1.1e-13, 6 clamped** |
  | B18 junction Mmax | 3.920 (clamped) | 5.220 (genuine, unclamped) |
  | Jacobian exactness | exact (B19) | **exact** (8.07e-09) |

  **The case FOR adopting `main` as default:**
  1. ★ **It is the principled model, not a tuning choice.** A plain element is
     uncut and single-valued; reading its density from a side field that
     imports the OTHER side of the wake is simply wrong. The reporting layer
     already made this exact call in 2026-07-14 (`element_mach2`'s
     `mixed_plain="main"` default) — adoption makes assembly and diagnostic
     agree instead of contradicting each other.
  2. Every committed **2.5-D** result is bit-identical (the class is 3-D only),
     so B3/B4/B6/B11–B17's quasi-2-D locks do not move at all.
  3. It **improves 3-D convergence** materially (M6 coarse reaches the target;
     the wing-body converges to machine precision instead of churning).
  4. The Jacobian stays exact, so B19's work composes.

  **The cost / what it does NOT buy:**
  - It **re-bases the 3-D committed level-set numbers** (M6 transonic γ/M_max
    in B7/B15, and the B9-LS / B16 / B17 / B18 wing-body values). Regenerating
    that committed evidence is hours of heavy compute — the real price.
  - It does **NOT** fix the junction pocket (GB20.5: real, G1.6-class geometry).
  - The M6 coarse γ moves 0.0780 → 0.0848 (+8.8 %) with no independent
    reference to say which is closer to truth; "converges to the target" is an
    objective improvement, "more accurate" is NOT claimed.

#### B20 re-baseline of the committed 3-D level-set evidence (2026-07-19)

Making the fix permanent re-bases the 3-D level-set numbers. **Scope, measured
rather than assumed:** the default suite is **465+22+2 unchanged** and the
gated 3-D LS tests are **67/67 green** — *no test lock breaks at all*.

★ **Why the suite could not see it, and that is itself a gap:** the 3-D numbers
are not locked by tests, they live in demo evidence. B15's tests, for instance,
run entirely on the 2.5-D NACA mesh (which B20 provably cannot touch), while
its M6-medium γ = 0.088338 is a *demo* number no test asserts. Same shape as
B19's blind spot, one layer out: **a change that moves 3-D level-set physics
raises no alarm in the suite.** Recorded as an open process gap.

★ **Demo-cache trap (cost me one false "regenerated" result):** the heavy demos
cache solves to gitignored `results/*.npz` and a plain re-run silently REUSES
them. B7's first "re-run" reported 35 passed with one punctuation diff — it was
`[M1] cached / [M4] cached`. Correct procedure: **delete the level-set `.npz`
first, then verify the log contains zero `cached` lines.** Conforming caches
may legitimately stay (B20 cannot touch that path).

**What actually moved.**

*Unchanged:* all 2.5-D; the whole conforming path; **B9's cross-model headline
(LS cl_p 0.2165 / cl_kj 0.2175 — not one digit)**; the gate verdicts of B9,
B16, B17, B18. *(Moved one digit, disclosed not hidden: B17's medium Newton
0.2115 → 0.2114.)*

*Improved — every one in the direction B20 predicts:*

| | before | after |
|---|---|---|
| B7 M6 M0.84 M_max | 1.453 | **1.392** |
| B7 tip Γ(z)→0 | −0.0003 | **−0.0000** |
| B16 coarse legacy limited cells | **3690** | **11** |
| B16 coarse pin floored cells | 3 | **0** |
| B17 pin_gamma medium clamps (= the B18 side leg; ⚠ erratum 2026-07-19: this row was mislabeled "B16 medium pin" — B16's medium pin was already 0/0 in its own pre-B20 CSV) | 42/40 churn | **0/0** |
| B18 wing-body medium residual | 6.8e-5 | **1.1e-13** |
| M6 coarse ramp → M0.84 | 0.7875, not converged | **0.84, converged** |

★ B16's `flr 3 → 0` is worth naming: those three floored cells were recorded in
B16 as "the B8/G1.6 junction class, orthogonal to the BC fix, not chased". They
were this contamination.

*Regressed — ONE case:* **M6 medium M0.84** (B15's committed recipe; B14 fails
with it because it asserts the same GB15.4 physics). The ramp reaches
**M0.6625 (2/5 levels)** where it used to reach M0.84; γ 0.088338 → 0.071909.
B15 now scores 17/20, B14 5/7, and **every failing check in both traces to this
one case**.

★★ **A finding that changes how the old number reads.** B15's gate compared
M_max against `PICARD_M6 = dict(..., m_max=2.4549)` — the LS **Picard**. Both
LS solvers read the same contaminated density, so **that check was common-mode:
it verified the two LS solvers agreed with each other, never that M_max was
right.** The conforming record on M6 medium is **1.995**; the old LS 2.45–2.49
sat well above it. Separately, the new M_max 1.5822 is measured at **M0.6625,
not M0.84** — different freestream, not a like-for-like number, so that FAIL is
downstream of the reach failure, not an independent finding. **The real
regression is one clause: the ramp no longer reaches M0.84 on this mesh.**

- [x] **GB20.7 — ANSWERED: a REAL capability loss, not a recipe mismatch**
  (`results/gb207_recipe_sweep.csv`, 2026-07-19). `freeze_tol` swept 1e-3 →
  1e-6 on the committed M6-medium call:

  | freeze_tol | m_final | levels | accept routes |
  |---|---|---|---|
  | 1e-3 (committed) | 0.6625 | 2/5 | assignment_cycle |
  | 1e-4 | 0.6625 | 2/5 | assignment_cycle |
  | **1e-5** | **0.6750** | **3/6** | + **tol** |
  | 1e-6 | 0.6750 | 3/6 | + tol |

  ★ **The hypothesis was directionally right and quantitatively insufficient.**
  Lowering `freeze_tol` does help — the ceiling moves 0.6625 → 0.6750 and a
  level starts converging on `tol` instead of escaping via `assignment_cycle`,
  confirming 1e-3 was indeed calibrated above a churn floor that B20 lowered.
  But the ceiling moves by 0.0125, not to 0.84. **The recipe was a contributor,
  not the cause.**
  **Honest bound:** only `freeze_tol` was varied (`dm`, `n_newton_max`,
  `m_start` were not), so this is the evidence-based conclusion, not a proof
  that no recipe can reach M0.84.

  ★★ **Synthesis — the contamination was acting as an unintended stabiliser.**
  The same pattern appears on the wing-body (GB20.5): with the contamination
  removed the solver converges *beautifully* (0/0 clamps, |R| ~1e-13) but
  cannot climb as far. Before B20 the M6-medium ramp reached M0.84 — into a
  state whose M_max 2.45–2.49 sits ~25 % above the conforming reference
  (**1.995**) and was only ever validated against the equally-contaminated LS
  Picard. **The trade is real and should be stated as a trade:** the old code
  went further into states we now have reason to distrust; the new code stops
  earlier and what it produces is clean.

  **Post-B20 level-set transonic envelope (M6 wake-free):** coarse **M0.84
  converged** (an improvement — it was 0.7875 not-converged), medium **≈M0.675**
  (was M0.84). ⇒ **GB15.4's "reaches M0.84" clause is now a NEGATIVE and
  B14's "== the committed GB15.4 physics" clause is superseded**; both need a
  re-spec against the new envelope (the G14.7 precedent: re-lock against what
  is now the honest oracle). Left OPEN for the user rather than re-specced
  unilaterally, since it redefines a committed capability claim.

  *(Stall symptom record, from the pre-answer entry — kept because it is what
  the sweep was designed against: **0/0 clamped at every level including the
  failed ones**, |R| 9.2e-14 where it converges, freeze armed every level with
  zero reverts, converged levels accepted via `assignment_cycle`. It is not
  diverging and not being clamped — it simply cannot climb on the now-clean
  field.)*

  ★★ **OVERTURNED BY B21 (2026-07-19, same day — see the B21 section).** The
  freeze_tol sweep varies WHEN the freeze arms; it structurally cannot see
  WHAT the freeze captures. The Kimi inspection's N1 found `freeze_side_state`
  computed its capture on the UNPATCHED side field (the one consumer the B20
  patch missed), so every armed freeze locked a selection the live system
  would not make (probe: 83+9 element differences, all aux-touching
  mixed-plain). B21's one-line fix restores the committed recipe to **M0.84**
  (γ 0.088343, res 9e-14, 515 s) ⇒ this gate's "real capability loss" verdict
  and the "contamination was a stabiliser" synthesis are both RETIRED; the
  loss was B20's own incomplete patch. What survives: freeze_tol was a real
  (minor) contributor under the broken capture, and the common-mode M_max
  critique stands. The GB15.4/GB14.4 re-spec question is moot in its original
  form — replaced by a small numeric re-baseline (GB21.2).

★ **Backport check (the A3 rule):** the conforming path has no side/aux DOFs
and no mixed-side plain class, so there is nothing to backport — recorded, not
skipped.

---

### B21 — Freeze-capture alignment with the B20 density (the Kimi-inspection N1 fix) ✓ CLOSED 2026-07-19 (NEW; appended after B20, no renumber; executes the 2026-07-19 second-round inspection's N1 finding)

**Trigger.** The 2026-07-19 Kimi inspection
(`docs/inspection/20260719-0555-full-inspection.md`, N1) found — and confirmed
with a committed probe (`docs/inspection/20260719-n1-freeze-probe.py`) — that
`freeze_side_state` was the ONE consumer of the side q²/ρ path the B20 patch
missed: it captured the frozen (upstream, branch) selection on the UNPATCHED
side field while `newton_side_data` applies `_apply_main_density` first. The
docstring invariant "the frozen sweep reproduces the live density bitwise at
the freeze point" was broken on 3-D meshes: at the GB19.6 workload (M6 coarse,
seeded M0.70) the shipped capture differed from the live system's selection on
**83 upstream + 9 branch** entries, ALL on aux-touching mixed-plain elements.
Pre-B20 the capture and the live sweep were consistent (both unpatched); B20's
partial patch CREATED the inconsistency — so every armed freeze since the
re-baseline locked a selection the live system would not make.

**Fix.** One line: the capture path now runs `_apply_main_density` before
`freeze_upwind_state` (`pyfp3d/wake/multivalued.py::freeze_side_state`),
restoring the invariant by construction. 2.5-D bit-identical (0 aux-touching
mixed-plain elements — the B19 erratum invariant; `test_v0_freestream` + the
full b15/b19/b6 files green untouched).

- [x] **GB21.1 — the M6-medium discriminator: N1 WAS THE MECHANISM; GB20.7 is
  OVERTURNED.** `cases/analysis/c1_ls_jacobian_fd/run_n1_freeze_capture.py`
  re-runs the SAME committed B15 M6-medium call as the GB20.7 sweep, post-fix:

  | freeze_tol | reached | levels | γ | M_max | res | clamps | wall |
  |---|---|---|---|---|---|---|---|
  | 1e-3 (committed) | **M0.84 ✓** | 6/6 | **0.088343** | 2.4818 | 9.0e-14 | 0/1 | **515 s** |
  | 1e-5 (GB20.7's best) | **M0.84 ✓** | 6/6 | 0.088343 | 2.4818 | 7.8e-14 | 0/1 | 540 s |

  (`results/n1_freeze_fix_sweep.csv`.) The committed recipe reaches M0.84
  again — faster than pre-B20 (657 s) and cleaner (0 lim/1 flr vs 3 clamped;
  res ~9e-14; some levels now accept via `tol`, not only `assignment_cycle`),
  and the two freeze_tol decades agree on γ to 5e-7 (the recipe sensitivity
  GB20.7 measured is gone). ⇒ **GB20.7's "real capability loss" verdict is
  overturned: the loss was B20's own patch gap, not the fix's intrinsic
  cost.** GB20.7's sweep could not see this axis by construction — freeze_tol
  varies WHEN the freeze arms, N1 is WHAT it captures. The "contamination was
  an unintended stabiliser" synthesis is likewise RETIRED — the stabiliser
  story explained a symptom whose real cause was capture/live inconsistency.
  ★ **What survives of GB20.7:** the honest observation that freeze_tol was a
  contributor (1e-5 helped 0.6625→0.675 under the broken capture), and the
  common-mode critique of the old M_max validation (2.48 is still only
  LS-vs-LS checked; the conforming record 1.995 is a different mesh family —
  an open cross-model question, unchanged by B21).
- [x] **GB21.2 — capability restored ⇒ the GB15.4/GB14.4 re-spec question is
  MOOT in its original form.** GB15.4's "reaches M0.84" clause holds again;
  what remains is a small numeric re-baseline (γ 0.088338 → **0.088343**,
  M_max 2.4938 → **2.4818**, 3 clamped → 1 floored) — the B15/B14 demos and
  gated anchors should be refreshed against these values (the demo re-run is
  the recorded follow-up, not blocking: the solver call itself is the evidence
  here, committed as CSV). *(Executed by B22, 2026-07-19: B15 demo 20/20.)* **Post-B21 LS transonic envelope (M6 wake-free):
  coarse M0.84 converged AND medium M0.84 converged.**
- [x] **GB21.3 — the invariant is test-locked on 3-D** —
  `tests/test_b15_ls_newton_freeze.py::test_freeze_capture_matches_live_density_3d`
  (gated; M6 coarse at the probe's exact workload), with premise asserts
  (aux-touching mixed-plain > 0, nu_max > 0) so it cannot pass vacuously;
  **verified FAILING with the fix stashed** (the A3 discipline) and passing
  with it: 13 passed. This is the N2-adjacent lock for the freeze path; the
  broader N3 gap (3-D LS numbers unlocked) remains open.
- ★ **Backport check (the A3 rule):** the conforming `freeze_upwind_state`
  has ONE field — no side/aux, nothing to align. N/A, recorded.
- ★ **Fourth bite from one element class** (B8 diagnostic → B19 Jacobian →
  GB19.6 residual → B21 freeze capture). The standing rule is now written
  into design_track_b.md §21: **any consumer of the side-field q²/ρ must
  explicitly decide and record its density source** — grep for
  `velocities(phi_up)` / `velocities(phi_lo)` when touching this layer.

---

### B22 — Evidence refresh + 3-D LS anchor locks + re-baseline process rule ✓ CLOSED 2026-07-19 (NEW; appended after B21, no renumber; executes B21's recorded follow-up + the Kimi-inspection N3/§2/§5 items, user-directed "依次解决4条遗留问题")

**Scope.** Four leftovers from B21's close-out, in the recorded order: (1) the
B15/B14 demo refresh against the B21 state; (2) the N3 test-lock gap; (3) the
re-baseline erratum process rule; (4) the next-phase priority re-assessment.
**No `pyfp3d/` change** — demos re-run, tests added, docs/process only.

- [x] **GB22.1 — B15 demo refreshed: 20/20 PASS** (was 19/19 at close-out,
  17/20 under the first B20 re-baseline). LS npz caches deleted first, zero
  `cached` lines (the demo-cache trap). The refreshed committed record: M6
  medium M0.84 γ **0.088343**, M_max **2.4818**, |R| 9.048e-14, **0 lim/1
  flr**, 6/6 levels (M0.60 accepts via `tol`, the rest `assignment_cycle`),
  freeze armed everywhere with 0 reverts, **511 s = 4.51×** vs the committed
  Picard 2304.7 s. Cross-checked bit-consistent with B21's
  `n1_freeze_fix_sweep.csv` (γ and |R| to the printed digit — same-machine
  determinism measured, which is what licenses the tight anchor bands below).
- [x] **GB22.2 — B14 demo refreshed: 7/7 PASS** (was 5/7 under B20). Medium
  ramp lagged **505 s** vs schur **345 s = 1.47×**, precond share 45.8 % →
  **1.8 %**, γ **0.088343** on BOTH arms (|Δγ| 8.6e-13), 0 fallbacks; the
  GB14.4 anchor constants re-pinned to the B21 state (γ 0.088343 / M_max
  2.4818 / 515.3 s, provenance noted in the script). ★ Disclosed: the M6
  COARSE ramp also moved under B21 — γ 0.0848 (GB20.4-era) → **0.084931**,
  M_max **1.3684** (34.5 s) — the freeze capture patch touches every 3-D
  freeze-armed ramp, not just medium. ⚠ Found en route: B14's
  `demo_report.md` index row had been MISSING since its close-out (the D9
  close-out-debt class) — added, with the omission recorded in the row.
- [x] **GB22.3 — the N3 gap is CLOSED: gated absolute anchor locks**
  (`tests/test_b22_ls_3d_anchors.py`, +2 gated). Each test RE-SOLVES the
  committed B15 recipe (so it cannot pass vacuously) and asserts absolute
  anchors: `target_reached`, all levels converged, |R| < 1e-9, n_lim == 0,
  n_flr ≤ 1, γ within rtol 1e-4, M_max within rtol 1e-3 — coarse (γ
  0.08493098 / 1.3684, ~35 s) and medium (γ 0.088343 / 2.4818, ~9 min).
  Bands: the two B21 sweep variants agreed on γ to 5.5e-6 relative ⇒ 1e-4
  carries ~20× margin yet catches a B20-sized move (−18 %) by four orders.
  **This is the alarm that did not exist when the suite stayed green through
  two re-baselines in two days.**
- [x] **GB22.4 — re-baseline erratum checklist is now process** — CLAUDE.md
  workflow step 5 gained the ★ clause (grep the moved values across docs/,
  correct or annotate every old-section quote in the same commit) and
  agent-rules discipline **#11** points at it. Executes the 2026-07-19
  inspection's §2 recommendation; its D1/D2/D7/D8 were all products of this
  rule not existing.
- [~] **GB22.5 — next-phase priority re-assessment (RECORDED, decision =
  user).** `docs/analysis/next_phase_priorities_2026-07-19.md`: recommends
  **P11 (curved wall elements) first** — post-B20/B21 the wing-body line has
  no numerical suspects left, and the G1.6 faceted-geometry root now owns
  three refinement-worsening wounds (sphere Cp 11.6 %, GB9.4 fuselage lift
  16–20 %, LS junction pocket Mmax 5.22); LS fine route second (technically
  ready via B14, but no blocked question needs it yet); Track V after P11
  (viscous corrections on top of a refinement-worsening geometry error would
  inherit it), V1's 2.5-D ladder parallelizable; plus a cheap A-track
  M_max same-family cross-model check (the 2.4818-vs-1.995 gap is currently
  a cross-FAMILY comparison and must not be quoted as a defect).
- ★ **Backport check (the A3 rule):** no solver change in this phase — N/A,
  recorded.

### B23 — Wing-body junction discriminator: the pocket is the wake inboard FREE-EDGE singularity ✓ CLOSED 2026-07-19 (NEW; appended after B22, no renumber; executes the P11 close-out junction-discrimination requirement)

**Trigger.** B18/B20 left the junction pocket attributed to "G1.6 faceted
geometry", but P11 close-out forbids reusing that label for the junction
(crease geometry ≠ smooth-wall faceting) without a dedicated discrimination
experiment. Pre-registered campaign
(`cases/analysis/b23_junction_discriminator/`, PRE_REGISTRATION + VERDICT +
committed CSV/PNG).

- **D1 α-sweep (LS Picard M0.5, freestream+pin_gamma): the pocket is
  lift/wake-coupled, NOT geometric.** α=0 is clean at BOTH levels (medium
  Mmax 0.66 / coarse 0.64, cl_fus ≈ 0 self-check ✓, no pocket); the pocket
  appears with α and grows SUPERLINEARLY (medium α=3.06 corrM **14.66**, 104
  supersonic elements, peak BEHIND the fuselage at x=2.13, z≈z_junc; coarse
  same-type 1.14). Pre-registered branch 2 ⇒ wake/lift-coupling driven.
- **Attribution: the wake sheet's inboard FREE EDGE** — the sheet ends at the
  junction station (q≥0) instead of reaching the fuselage surface, trans-
  planting the P13-class free-edge singularity from the (handled) wing tip to
  the (unhandled) inboard end. W2 (fuselage spurious lift) characterized in
  the same campaign.
- **Routing:** (b)-1 waterline extension → B24 (closed negative); (b)-2 P13
  free-edge treatment rescoped to the inboard end → realized as B25.
- ★ **No `pyfp3d/` numerics change** — pure analysis phase.

### B24 — Wake inboard-end waterline extension: route CLOSED (negative) ✓ CLOSED 2026-07-19 (NEW; appended after B23, no renumber; pre-registered E1)

**Trigger.** B23's (b)-1 candidate: extend the TE polyline to the fuselage
waterline so the free end lands on the body. Pre-registered
(`cases/analysis/b24_wake_inboard_end/`).

- **The free-edge hypothesis is re-confirmed** — the pocket FOLLOWS the free
  edge: B1 (flush waterline) moves the peak past x_tail at all probed legs
  (medium α=2.0: 1.71@1.57 → 5.58@2.48). Mechanism certain.
- **But both extension variants trade the singularity for equal-or-worse
  forms:** B1 medium α=3.06 corrM **78.56 non-converged** (solver killed);
  B3 (offset-cone δ) migrates the pocket BACK into the near field on every
  leg (9.97–16.64 @ x≈2.15–2.32).
- **Decision-tree exit 3: "extension-class insufficient" — the (b)-1 route
  is CLOSED**; fall back to B23 (b)-2 (P13 free-edge treatment rescoped to
  the inboard end), realized by B25 as the fragment clip.
- ★ **No `pyfp3d/` numerics change** — pure analysis phase.

### B25 — Inboard fragment clip: the junction pocket CURED (C-A) ✓ CLOSED 2026-07-19 (NEW; appended after B24, no renumber; pre-registered F1, v2.1)

**The cure.** Give the free edge a legal home = the conforming fragment
topology: move the sheet's inboard boundary from "junction station q≥0" to
the fuselage surface / symmetry plane (`inboard_clip`, built by
`meshgen/fuselage.py:make_inboard_clip`, consumed in
`wake/cut_elements.py:122,190-204`). **Default None ⇒ bit-identical** (the
A side of every A/C comparison is the default path). Pre-registered
(`cases/analysis/b25_inboard_fragment_clip/`).

- [x] **Primary criteria all decisive (medium α=3.06):** corridor corrM
  **14.66 → 0.63** (≤1.3), corridor n_sup **88 → 0**, cl_p **+0.38 %**
  (within [A, oracle 0.2173] — the A1 direction clause), |Δγ| +0.37 %,
  root te_jump profile 0.28 %, α=0 inert (|Δcl_p| ≈ 0, peak 0.62@x=1.79
  side-identical), convergence 56 outer ≤ 1.5×A, no singular/clamped rows,
  strip-jump metric 1.16×γ (no anchoring pathology), sliver min dihedral
  11.0° / p05 37.3°.
- [~] **One secondary guardrail recorded non-blocking:** out-of-band
  fuselage-lift carryover +135 % (0.0214 → 0.0504) exceeds the literal 20 %
  band; same-code oracle comparison attributes it to the flat-vs-tilted
  sheet MODEL difference, not topology failure or pre-registration rot
  (note 2 of the VERDICT). Becomes a P11/curved-wall-element watch item
  (B27 §6.3).
- **Evidence chain closed:** B23 (attribution) → B24 (the pocket follows
  the free edge) → B25 (legal home cures it). α=2 legs same-type.
- Tests: `tests/test_b1_cut_elements.py::TestInboardFragmentClip` (4) +
  regression locks in `tests/test_m2_wingbody.py` / `test_b1_cut_elements.py`.

### B26 — Post-cure LS transonic ceiling re-measured: the pocket WAS the limiter (B26-A) ✓ CLOSED 2026-07-20 (NEW; appended after B25, no renumber; pre-registered G1)

**Trigger.** With the pocket cured (B25), does the B18 "junction-limited"
LS ceiling move? Same-code A/C (default vs `inboard_clip`), the B18 recipe
frozen verbatim, honest-stop ramp. Pre-registered
(`cases/analysis/b26_ls_transonic_ceiling/`).

- [x] **GB26 ceiling-lift (B26-A):** C medium m_last **0.50 → 0.7625**
  (climbs five loose rungs 0.55–0.75 all 0–1 clamp + 0.7625 strict res
  2.6e-11; dies 0.775); C coarse 0.82 → **0.84 REACHED** (strict res
  6.9e-11, freeze catches 6/1). A sides unchanged (medium 0.50 dies
  0.5125 (a)-class; coarse 0.82 dies 0.84).
- [x] **Death-cause flip (a) → (b):** the pre-cure medium died of
  (a)-class pocket rejection (8/3 > freeze_max_clamped=8 at 0.55 loose);
  post-cure BOTH levels die of **(b)-class high-M Newton stall with the
  peak at the WING TIP** (medium M4.18 @ z=1.20, corridor corrM 1.07
  clean; P13-class — the same class as the conforming medium 0.80+ stall),
  not the junction.
- [x] **cl_p same-trend vs conforming anchors:** C coarse 0.2542@0.84 vs
  conf 0.2617@0.84; C medium near-death 0.2475@0.775 vs conf 0.2579@0.79
  — 2–4 % low = the known B17 LS cl_p↔cl_kj convention gap, monotone
  transonic rise 0.1289@0.50 → 0.1491@0.7625 (+16 %).
- ~69 min total solve (T2 budget). The B18 LS ceiling story is hereby
  superseded — see the B18-entry ★★★ erratum.

### B27 — B18 demo refresh: LS legs resurrected, 8/8 checks PASS ✓ CLOSED 2026-07-20 (NEW; appended after B26, no renumber; pre-registered GB27.1–27.5)

**Trigger.** The B18 demo still told the retired "junction-limited
closed-negative" story. Full re-solve from no cache (~1 h 39 min, T2
budget) on the B25/B26 state; `cases/analysis/b27_b18_demo_refresh/` +
refreshed `cases/demo/b18_wingbody_transonic/` (checks.csv **8/8 PASS**).

- [x] **GB27.1 (PASS):** conforming legs bit-reproduce the committed B18
  anchors (0.2173/0.2321/0.2579/0.2617 + cross 0.2178 + Mmax 2.15 +
  reached flags) — B21/B22 inert on the conforming path, zero T1 drift.
- [x] **GB27.2 (PASS):** LS A/C ceiling legs bit-reproduce B26 committed
  (36 summary + 292 level items; `g27_consistency.csv` **336/336
  bit-identical**).
- [x] **GB27.3 (0.65 PASS / 0.75 RECORDED):** the transonic cross-model
  NOW EXISTS — M0.65: conf 0.2321 vs LS+clip 0.2266 = **2.4 % ≤ 5 %
  PASS**; M0.75: conf 0.2483 vs LS+clip 0.2421 = **2.5 %** (RECORDED,
  no-threshold per pre-registration). The gap is **flat across Mach**
  (M0.5 2.6 / M0.65 2.4 / M0.75 2.5 %) = one ~2.5 % B17 convention band,
  not physics.
- [~] **GB27.4 (RECORDED):** GB18.4 re-answered in-demo (pocket = B23
  inboard free-edge singularity, cured C-side; residual limiter = wing-tip
  P13 + high-M Newton, peaks quoted from B26 committed); GB18.5 refreshed
  (conf cl_fus 0.0423 = 16 % @0.79 live; C-side new-ceiling cl_fus 0.0781
  / out-band 0.0565 ≈ ×2 → P11 watch item — **CLOSED by B28 (2026-07-20):
  wake-sheet position sensitivity, not a lesion**).
- [~] **GB27.5 (RECORDED):** the T1 erratum is in the demo docstring —
  the A-side re-test climbing past the B18 committed anchors (dies
  0.50/0.55) is the **B21/B22 freeze-capture repair effect**, not physics
  drift; the pocket's true kill line A medium = 0.55 (Mmax 13.1 >
  freeze_max_clamped=8).
- **Independent observations:** coarse 0.60 cross — the cured C-side value
  0.2133 (gap 2.1 %) retires the old A-side 0.2 % row (a pocket-
  contaminated coincidence); `b18_sections_conf_medium.png` had been
  silently empty since an earlier `section_cp_curve` API drift — fixed
  here (NOT a B21/B22 regression); conf medium **0.75 NEW point 0.2483**
  strict (res 8.3e-11), cl(M) monotone 0.2173/0.2321/0.2483/0.2579, the
  0.80+ stall narrative unchanged.
- **Façade:** the B18 demo headline is re-worded from "LS junction-limited
  (closed-negative)" to "**post-cure LS ceiling co-located with
  conforming** (coarse 0.84 = 0.84; medium 0.7625 ≈ 0.79)"; conforming
  drops to cross-validator, and the Track V sheet-topology prerequisites
  are all in place (B25 §7).
- ★ **No `pyfp3d/` numerics change** — demo/tests/docs phase.

## Progress ledger

### Track B — level-set embedded wake

Track status: **◐ IN PROGRESS** — **B27 ✓ CLOSED 2026-07-20** (B18 demo refreshed on the B25/B26 state:
8/8 checks PASS, 336/336 bit-identical; transonic cross-model now exists — M0.65 2.4 % PASS / M0.75 2.5 % RECORDED;
façade re-worded "junction-limited" → "post-cure LS ceiling co-located with conforming") · **B26 ✓ CLOSED 2026-07-20** (post-cure LS ceiling
re-measured: medium 0.50 → **0.7625**, coarse 0.82 → **0.84 reached**; death cause flips (a)-pocket → (b)-wing-tip P13 + high-M Newton) ·
**B25 ✓ CLOSED 2026-07-19** (the CURE: `inboard_clip` gives the wake inboard free edge its conforming-fragment home —
medium α=3.06 pocket corrM **14.66 → 0.63**, guardrails clean; default None bit-identical) · **B24 ✓ CLOSED 2026-07-19** (negative:
the pocket follows the free edge but waterline-extension variants trade the singularity for worse —
(b)-1 route CLOSED) · **B23 ✓ CLOSED 2026-07-19** (junction discriminator: pocket is lift/wake-coupled, α=0 clean;
attribution = wake inboard FREE-EDGE singularity) · **B22 ✓ CLOSED 2026-07-19** (evidence refresh:
B15 demo 20/20 + B14 7/7 on the B21 state; **N3 closed** — gated absolute anchor locks on the M6 coarse/medium ramps;
re-baseline erratum checklist = CLAUDE.md step 5 / discipline #11; next-phase analysis recommends P11, user's call) · **B21 ✓ CLOSED
2026-07-19** (N1 freeze-capture alignment: `freeze_side_state` missed the B20 patch; the one-line fix RESTORES the M6-medium M0.84 ramp —
γ 0.088343, res 9e-14, 515 s — ⇒ **GB20.7 overturned**; 3-D capture lock added, verified failing pre-fix) · **B20 ✓ CLOSED 2026-07-18,
ADOPTED PERMANENTLY (user-arbitrated 2026-07-18; the A/B knob REMOVED), re-baselined 2026-07-19** (mixed-plain main-field density;
the apparent M6-medium regression was B20's own patch gap, resolved by B21) · **B19 ✓ CLOSED 2026-07-18** (LS-Newton Jacobian made exact in
3-D; Leg B measured the spurious-supersonic side-field contamination in the same element class and routed it to B20, which then adopted the
fix). — design 2026-07-07; B10 shelved 2026-07-10;
numerics spec [design_track_b.md](../design_track_b.md) (supersedes DN1) + gate
re-arbitration 2026-07-11; **B1 CLOSED 2026-07-11**, with M3/M4 delivered the
same day; next = B2 *(that opening timeline is HISTORICAL — the live status is
the ledger table below and the track line in agent-rules.md; as of 2026-07-20
B1–B9 and B11–B27 are closed, B6 ◐, B10 shelved)*. Purpose is user-arbitrated as **mesh/geometry workflow
capability, not solver speed** (the kill-the-Γ-secant efficiency motivation is
obsolete post-P8 Newton), so the efficiency criteria in the B-gates are
non-regression guards only. Coexistence strategy: a parallel `solve/picard_ls.py`
path with a per-phase default flip — the conforming-path solver numerics stay
byte-untouched. Sequencing guard: P8's Newton landed on the conforming wake
(closed), and a level-set Newton is a post-B6 re-derivation (simpler — the
wake-LS Jacobian blocks are constant in φ, no Γ elimination/Woodbury); Track B
blocks nothing in P7–P12, and M2 (wing-body) wants it.

- B28 — ✓ — 2026-07-20 — **cl_fus decoupling + GB9.4 RE-SPEC — the "fuselage spurious lift" label retired** (executes B23 §(c);
  pre-registered `cases/analysis/b28_cl_fus_flat_sheet/`). New `WakeLevelSet(sheet_direction=...)` knob decouples sheet GEOMETRY
  drag-out from physics convection (default None bit-identical; `TestSheetDirection` 4 + b1/m2/v0 73 green). Decisive leg
  (medium α=3.06, strict conv res 9.0e-8): flat-fragment out-band cl_fus **0.0326 vs conforming oracle 0.0351 = 7.25% ≤ 15% TOL**,
  flow-aligned sheet 0.0504 = 35.3% ⇒ **F1: sheet-POSITION sensitivity, not an error**. cl_fus = physical carryover (~10% cl_p,
  refinement-flat) + pocket imprint (B25 cured) + sheet position. GB9.4 re-spec'd: hard gate (medium) out-band cross-model ≤15%,
  coarse RECORDED + full band/out/poles decomposition; old ≤5% band retired by erratum. B9 demo LS leg switched to flat-fragment
  (legacy farfield kept), re-run **8/8 PASS** — GB9.4 medium gap **7.0%**, GB9.5 un-flipped (cl_p 0.5%/cl_kj 0.3%),
  cl_fus/wing LS 0.205→0.176; LS Picard 59/73 outer converged (220 budget). B26's "C-side out-band ×2" watch item CLOSED;
  M2 wiring suspicion excluded. Flat-fragment as PRODUCTION config for B18/GB18.5/GB27.4 left to the user.
- B27 — ✓ — 2026-07-20 — **B18 demo refresh — LS legs resurrected, checks.csv 8/8 PASS** (no `pyfp3d/` change;
  full re-solve ~1 h 39 min). **GB27.1 ✓** conforming legs bit-reproduce the committed B18 anchors (B21/B22 inert on conforming).
  **GB27.2 ✓** LS A/C legs bit-reproduce B26 committed (`g27_consistency.csv` **336/336 bit-identical**).
  **GB27.3 ✓/RECORDED** the transonic cross-model NOW EXISTS — M0.65 conf 0.2321 vs LS+clip 0.2266 = **2.4 % PASS (≤5 %)**;
  M0.75 conf 0.2483 vs LS+clip 0.2421 = **2.5 %** RECORDED; gap flat across Mach (2.6/2.4/2.5 %) = the B17 convention band.
  **GB27.4/27.5 RECORDED** (junction re-characterized in-demo; A-side climb past old anchors = the B21/B22 freeze-capture repair, not
  drift). Independent: coarse 0.60 cross old 0.2 % row retired (C-side 0.2133, gap 2.1 %);
  the silently-empty `b18_sections_conf_medium.png` (older `section_cp_curve` API drift) fixed;
  conf medium **0.75 NEW point 0.2483** strict, cl(M) monotone. Façade: "junction-limited (closed-negative)" → "**post-cure LS ceiling
  co-located with conforming** (coarse 0.84 = 0.84; medium 0.7625 ≈ 0.79)"; Track V sheet-topology prerequisites in place.
- B26 — ✓ — 2026-07-20 — **Post-cure LS transonic ceiling re-measured — the pocket WAS the limiter (B26-A).**
  Same-code A/C (default vs `inboard_clip`), B18 recipe frozen verbatim, ~69 min. C medium m_last **0.50 → 0.7625** (five loose rungs
  0.55–0.75 all 0–1 clamp + 0.7625 strict res 2.6e-11; dies 0.775); C coarse 0.82 → **0.84 REACHED** (strict res 6.9e-11).
  A sides unchanged. **Death-cause flip (a)→(b):** pre-cure medium died of (a)-class pocket rejection (8/3 > freeze_max_clamped=8);
  post-cure both levels die of **(b)-class high-M Newton stall, peak at the WING TIP** (medium M4.18 @ z=1.20;
  corridor corrM 1.07 clean) = the conforming 0.80+ stall class, no longer the junction. cl_p same-trend vs conforming (0.2542@0.84 vs
  0.2617; 2–4 % low = B17 convention gap). Supersedes the B18 LS ceiling story (★★★ erratum in the B18 entry).
- B25 — ✓ — 2026-07-19 — **Inboard fragment clip — the junction pocket CURED (C-A).** The wake sheet's inboard boundary moves from "junction
  station q≥0" to the fuselage surface / symmetry plane (= conforming fragment topology): `inboard_clip` via
  `meshgen/fuselage.py:make_inboard_clip` + `wake/cut_elements.py:122,190-204`; **default None ⇒ bit-identical**.
  Medium α=3.06 decisive leg: corridor corrM **14.66 → 0.63**, n_sup **88 → 0**, cl_p **+0.38 %** (within [A, oracle 0.2173]), — Δγ —
  0.37 %, root te_jump 0.28 %, α=0 inert, 56 outer ≤ 1.5×A, strip-jump 1.16×γ, sliver min dihedral 11.0°.
  One secondary guardrail (out-of-band cl_fus carryover +135 %) oracle-attributed to the flat-vs-tilted sheet model difference —
  recorded non-blocking, P11 watch item. Chain closed: B23 attribution → B24 follows-the-edge → B25 legal home.
  Tests `TestInboardFragmentClip` (4) + M2/B1 regression locks.
- B24 — ✓ — 2026-07-19 — **Wake inboard-end waterline extension — route CLOSED (negative).** E1:
  the pocket FOLLOWS the free edge (B1 flush moves the peak past x_tail on all legs — mechanism re-confirmed), but both extension variants
  trade the singularity for equal-or-worse forms (B1 medium α=3.06 corrM 78.56 NON-converged;
  B3 offset-cone migrates back near-field 9.97–16.64 @ x≈2.15–2.32 every leg). Decision-tree exit 3:
  "extension-class insufficient", (b)-1 CLOSED → back to B23 (b)-2, realized as B25. No `pyfp3d/` change.
- B23 — ✓ — 2026-07-19 — **Wing-body junction discriminator — the pocket is the wake inboard FREE-EDGE singularity** (P11 close-out
  requirement; pre-registered D-campaign). D1 α-sweep (LS Picard M0.5, freestream+pin_gamma):
  α=0 clean at BOTH levels (Mmax 0.66/0.64, cl_fus ≈ 0 self-check), the pocket appears with α and grows SUPERLINEARLY (medium α=3.06 corrM
  14.66, 104 supersonic elements, peak BEHIND the fuselage x=2.13 @ z≈z_junc) ⇒ **lift/wake-coupled, NOT faceted geometry**.
  Attribution: the sheet ends at the junction station (q≥0) instead of reaching the fuselage —
  the P13-class free-edge singularity transplanted inboard. Routes: (b)-1 → B24 (closed negative); (b)-2 → B25 (the cure).
  No `pyfp3d/` change.
- B22 — ✓ — 2026-07-19 — **Evidence refresh + 3-D LS anchor locks + re-baseline process rule** (executes B21's recorded follow-up + Kimi
  N3/§2/§5; no `pyfp3d/` change). **GB22.1 ✓** B15 demo **20/20** (caches deleted, zero `cached` lines):
  medium γ 0.088343 / M_max 2.4818 / — R — 9.048e-14 / 0 lim/1 flr / 6/6 levels / **511 s = 4.51×**, bit-consistent with
  `n1_freeze_fix_sweep.csv`. **GB22.2 ✓** B14 demo **7/7**: medium lagged 505 s vs schur **345 s = 1.47×**, precond → **1.8 %**, γ 0.088343
  both arms; anchors re-pinned; ★ coarse ramp also moved under B21 (γ 0.0848 → **0.084931**, M_max 1.3684);
  ⚠ B14's missing demo_report index row found + added. **GB22.3 ✓ N3 CLOSED**: `tests/test_b22_ls_3d_anchors.py` (+2 gated) re-solves the
  committed coarse (~35 s) + medium (~9 min) ramps and asserts m_final/γ/M_max/clamps absolutely (γ rtol 1e-4 = 20× the measured run-to-run
  spread, 4 orders under a B20-sized move) — the alarm that did not exist when the suite stayed green through two re-baselines.
  **GB22.4 ✓** CLAUDE.md step-5 re-baseline erratum checklist + agent-rules discipline #11.
  **GB22.5 RECORDED** `docs/analysis/next_phase_priorities_2026-07-19.md` — recommends **P11 first** (G1.6 owns three refinement-worsening
  wounds; wing-body line has no numerical suspects left), LS fine second, Track V after P11 (V1 ladder parallelizable);
  decision = user. Backport N/A.
- B21 — ✓ — 2026-07-19 — **Freeze-capture alignment with the B20 density (the Kimi-inspection N1 fix) —
  restores the M6-medium M0.84 ramp and OVERTURNS GB20.7.** `freeze_side_state` was the one consumer of the side q²/ρ path the B20 patch
  missed: it captured (upstream, branch) on the UNPATCHED side field while `newton_side_data` patches first —
  probe-confirmed **83 upstream + 9 branch** selection differences vs the live system (M6 coarse seeded M0.70, all aux-touching mixed-plain;
  `docs/inspection/20260719-n1-freeze-probe.py`). Pre-B20 both were consistently unpatched;
  **B20's partial patch CREATED the inconsistency** ⇒ every armed freeze locked a selection the live system would not make = the actual
  mechanism of the "GB20.7 capability loss". **Fix = one line** (capture runs `_apply_main_density`); 2.5-D bit-identical.
  **GB21.1 ✓** committed recipe reaches **M0.84 again**: γ **0.088343** (pre-B20 0.088338), M_max 2.4818, res 9.0e-14, 0 lim/1 flr, **515
  s** (pre-B20 657 s/3 clamped); freeze_tol 1e-3 and 1e-5 agree to 5e-7 (`n1_freeze_fix_sweep.csv`).
  The "contamination was a stabiliser" synthesis is RETIRED. **GB21.2 ✓** GB15.4's capability clause holds again —
  remaining is a small numeric re-baseline (B15/B14 demo refresh = recorded follow-up); post-B21 envelope:
  M6 coarse AND medium reach M0.84. **GB21.3 ✓** 3-D freeze-capture lock `test_freeze_capture_matches_live_density_3d` (gated,
  premise-asserted, verified FAILING pre-fix). Backport N/A (conforming has one field). ★ Fourth bite from the mixed-plain class ⇒
  design_track_b.md §21 rule: every side-field q²/ρ consumer must decide+record its density source.
- B20 — ✓ — 2026-07-18 — **Mixed-side plain elements now read their density from the MAIN field —
  PERMANENTLY, with no switch (user-arbitrated).** Executes GB19.6. ★★ **Why it is not optional:**
  the element's stiffness is contracted with the MAIN field (it scatters onto main DOFs) while its density came from a DIFFERENT field —
  **one equation built from two velocity fields, internally inconsistent** — and the element is UNCUT, so no wake jump passes through it at
  all. That argument needs no physics and admits no "preference" reading. ★ **The reporting layer had already ruled the same way** —
  `element_mach2` has defaulted to `mixed_plain="main"` since 2026-07-14 — so B20 makes the ASSEMBLY agree with the DIAGNOSTIC.
  The fix was first built as a default-off `plain_density` knob purely to measure the A/B below;
  on that evidence the user directed it be hard-coded and the knob REMOVED. ★★ **A workspace-aliasing bug was caught by measuring an
  unexpected result:** `PicardOperator.velocities` returns VIEWS into a shared buffer, so recomputing the main gradient inside the density
  path overwrote the caller's side values in place — quasi-2D "moved" 0.77, subsonic Γ tripled, the Jacobian degraded:
  all ONE bug (2940 elements clobbered vs the 129 in the mask), fixed with `.copy()`. *Accepting the plausible "Γ tripled ⇒ big effect"
  story would have recorded an aliasing bug as a physics finding.* **GB20.1 ✓** quasi-2D bit-identical (0.000e+00), M6 R moves on 164 of
  ~12k rows. **GB20.2 ✓** the Jacobian stays EXACT (targeted 8.07e-09 / control 6.29e-10; the gated 3-D FD gate passes 3/3 post-permanence)
  — B19 Leg A ∘ B20 compose. **GB20.3 ✓** 2.5-D subsonic Γ 0.088144 → 0.088144 (**+0.0000 %**):
  the class is 3-D only, so every committed quasi-2D lock is untouched. **GB20.4 ✓** M6 coarse ramp→M0.84:
  old **m 0.7875 NOT converged** → new **M0.84 CONVERGED** (γ 0.0780→0.0848; "converges to target" is the objective gain, "more accurate" is
  NOT claimed — no independent reference). **GB20.5 ✓ RECORDED — the hypothesis SPLITS.** B18 medium wing-body @M0.5:
  old res 6.8e-5 / **82 clamped** / Mmax 3.920 (a CLAMPED, NON-converged number) → new **res 1.1e-13 / 6 clamped** / Mmax 5.220 (a genuine
  converged solution). ⇒ ★ the **CONVERGENCE** pathology was largely this contamination (B18's churn/clamping substantially cured), but ★
  the junction **POCKET is REAL** and B19's literal hypothesis is **REFUTED** — removing the contamination UNCLAMPED it, revealing a genuine
  M≈5.2 spike at subsonic freestream = the **G1.6/GB9.4 faceted-geometry** error, not the mixed-plain density.
  It still cannot pass M0.5. **Do not repeat "the pocket is mixed-plain contamination" — measured false.** **Accepted cost:**
  the pre-B20 3-D level-set numbers were solutions to a discretization carrying a known internal inconsistency and are re-based;
  old values stay traceable in git and in the before/after CSVs. The three A/B scripts are marked HISTORICAL (the knob they toggle is gone;
  reproduce at commit 5369a84); the standing Jacobian check is `tests/test_b19_jacobian_3d.py`.
  Evidence `cases/analysis/c1_ls_jacobian_fd/results/legb_*.csv`.
- B19 — ✓ — 2026-07-18 — **The LS-Newton Jacobian is now exact in 3-D (Leg A), and the residual's own asymmetry is measured and routed (Leg
  B).** Executes the A3/GA3.6 C1 finding; opened + closed same day, user-directed, two deliberately separated legs.
  **★ Leg A was TWO defects, not one.** (1) **DOF maps**: Terms 2/3 used the mass-conservation SCATTER map for both rows and columns, but
  the columns must follow `side_potentials`' per-node READ map — they coincide on cut elements (asserted:
  `readvec` reproduces `dofs_upper`/`dofs_lower`) and diverge on mixed-side plain elements, a 3-D-only class.
  (2) **Gradient factors, the same duality one level down**: the residual is
  `rho_tilde(grad of the READ field) * V * (grad of the SCATTER field . B_a)`, so the ROW factor must use `grad_row` while the COLUMN factor
  keeps the side gradient; the code used the side gradient for both. ★ Fixing (1) alone left the probe at **1.4697e-02** —
  8× better but **still ε-independent** ⇒ recorded PARTIAL instead of rounded into a pass, which is what forced (2) out.
  ★ A column-wise FD localization plus a block isolation (`\|FD23−J23\| ≡ \|FD23−J2\|` exactly ⇒ Term 3 contributes nothing there) found (2)
  — **the first reading of the row classification had pointed at Term 3 and was wrong**; acting on it would have put a new bug into correct
  code. **GB19.1 ✓** targeted probe **1.145684e-01 → 1.333699e-08**, control 6.327479e-10 unchanged, global-free 2.47e-03 → 8.49e-10;
  ★★ **the ε discriminator FLIPPED** — pre-fix 1.532e-01 at every ε (spread 1.00 = a missing term), post-fix 1.6e-09/2.1e-08/2.2e-07 (spread
  131.5, ~1/ε = pure FD roundoff): an ε-independent error became ε-sensitive, the transition a real fix must produce and one no
  tolerance-loosening can fake. **GB19.2 ✓** `max\|ΔR\| = **0.000e+00**` bit-identical on a 3-D mesh, verified by `git stash` A/B after EACH
  fix ⇒ no converged level-set result can move; asserted, not inferred. **GB19.4 ✓ RECORDED NEGATIVE —
  the fix buys NO convergence**: pre/post 40 steps, not converged, γ **0.07212068** identical to 8 dp, M_max **1.134235** to 6 dp, same
  plateau (residuals differ only by limit-cycle phase), **+3.6 % wall** for the extra einsum.
  The plateau is the **B15 selection-churn limit cycle** (0 limited/0 floored, clean period) and an exact derivative cannot fix a
  discontinuous selection — **the phase must not be credited with a convergence improvement**;
  what it buys is correctness (a Newton method rather than a quasi-Newton). **GB19.5 ✓** `tests/test_b19_jacobian_3d.py` (+3:
  gated 3-D FD gate + 2 structural locks) closes the blind spot that let this live through B6–B18; ★★ **ERRATUM:
  the blind spot was mis-stated** by the C1 write-up AND by this author's first test — quasi-2-D has **129** mixed-side plain elements, not
  zero; what it has zero of is the subset that can READ an aux (**0 of 129 touch a cut node**), which is the real invariant.
  **GB19.6 ✓ (Leg B) — the residual's asymmetry is NOT benign.** Zero-change probe: only **252** elements (**0.1888 %** of volume) actually
  read an aux, but there max\|ρ_side − ρ_main\| = **0.4474 (45.3 %)**, 154/39 of 252 differ >1 %/>10 %, and decisively **the SIDE field
  reads q² up to 3.2229 (M≈1.80, at the M_cap limiter) where the MAIN field reads 1.3379** —
  a **spurious supersonic state** that the artificial-density switch then acts on, i.e. the contamination reaches the solver's density, not
  just a diagnostic. ★ Third bite from one element class (B8's ×5 `element_mach2` metric artifact → A3/B19 Jacobian → now the residual).
  ★ **HYPOTHESIS, not a result:** a candidate contributor to B18's refinement-worsening wing-body pocket (M_max artifact 3.96);
  measured on wing-alone M6, no causal link shown; named test = rerun B18 medium LS with main-field density and see if 3.96 moves.
  **NOT adopted** — changing the density source changes R and moves every converged answer ⇒ its own phase.
  Evidence `cases/analysis/c1_ls_jacobian_fd/` (`run_check.py`, `run_legb_probe.py`, `b19_three_states.csv`, `b19_convergence_ab.csv`,
  `legb_density_gap.csv`, `c1_fd_probes{,_prefix}.csv`).
- B1 — ✓ — 2026-07-11 — **B1 delivery (2026-07-11):** `pyfp3d/wake/levelset.py` (TE-**polyline** ruled straight wake per design_track_b.md
  D9, per-segment frames, `update_direction()` re-aims the wake without touching the mesh) + `pyfp3d/wake/cut_elements.py` (ε side-shift
  relative to local edge length (D4), **downstream-crossing test** excluding the ahead-of-LE sign-change region, TE-node flagging, below-TE
  fan recorded as `te_lower_elems` for B2's López-fig-3.6c aux assignment, per-node ext DOFs, López eq.
  3.33–3.34 `dofs_upper`/`dofs_lower` tables); imported by nothing in the shipped solver paths.
  Gate evidence (`tests/test_b1_cut_elements.py`, **34 passed**, the FULL dual-mesh matrix —
  2.5D M0/M3 coarse+medium AND 3D M1/M4 ONERA M6): M0 embedded — every conforming sheet node ε-shifted "+" (the D4 stress test at scale),
  census cross-validated EXACTLY against `cut_wake` (`cut_elems ∪ te_lower_elems` == the minus-side element star, element-by-element), TE
  nodes == `wc.te_nodes`; M3 wake-free — generic cuts, gap-free corridor TE→far field at α=0 AND re-aimed to α=4° **on the same mesh**;
  M1/M4 ONERA M6 — census a strict **superset** of the conforming minus-star (0 missing, +2.9% tip-edge straddlers:
  expected, since in an embedded method the sheet's tip EDGE need not conform), spanwise clip verified.
  ★ **Two 3D-only mechanisms found and fixed here** (both invisible on quasi-2D meshes): (1) the swept TE span axis is NOT perpendicular to
  the wake direction ⇒ q must come from the **oblique (v, d̂, n̂) frame** — an orthogonal projection leaks the downstream distance into the
  spanwise coordinate and wrongly clipped ~60% of the true M6 cut set (measured, fixed, regression-pinned);
  (2) the **spanwise clip** (crossings must satisfy 0 ≤ q ≤ span_length) is mandatory — without it the level set cuts the wake-plane
  extension beyond the tip, i.e. P5's far-field branch-ray artifact re-created (the conforming path gets the same semantics from its
  free-edge rule, Γ(tip)=0). Suite **218+8+2** (was 184+8+2; +34, some of which skip when the gitignored wake-free meshes aren't generated
  locally); conforming solver paths byte-untouched; all runs at the 8-thread cap alongside the in-flight P9 fine demo.
- B2 — ✓ — 2026-07-11 — **B2 delivery (2026-07-11):** multivalued (CutFEM-style) FE assembly.
  `pyfp3d/kernels/cut_assembly.py` (`multivalued_redirection_coo` + `continuity_closure_coo`) +
  `pyfp3d/wake/multivalued.py::MultivaluedOperator` (extended n_total = n_main + n_ext DOF assembly, TE-jump/Γ extraction) +
  `pyfp3d/solve/picard_ls.py::solve_multivalued_laplace` (non-lifting direct-LU driver, parallel to the conforming path).
  **Key simplification (design_track_b.md §2.5/D6):** a cut element is the same P1 element matrix assembled twice with B1's `dofs_upper`
  /`dofs_lower`; expressed as a sparse redirection of the single-valued matrix — only the entries whose two nodes are on OPPOSITE sides move
  their column main(b)→aux(b), everything else byte-identical to `PicardOperator.assemble_matrix()`.
  Aux rows carry the B2 continuity ("weld") closure aux_k = main_j, so the extended system reduces EXACTLY to the single-valued one
  (`test_extended_matrix_folds_to_stiffness`: fold recovers the stiffness matrix to 1e-13).
  Extended matrix is nonsymmetric ⇒ `spsolve`; GMRES+AMG deferred to B3+ scaling (design_track_b.md §5.3).
  Gate (`tests/test_b2_multivalued.py`, **17 passed**, coarse+medium both 2.5D families + 3D M6 coarse both families;
  medium/M6 skip in CI where gitignored): V0 freestream **0.0** (2.5D, α=0/4°) / **1e-14** (3D M6) < 1e−12;
  V1 MMS slope **1.94** ≥ 1.9 (generic-position cube cut); Laplace α=0 ⇒ TE jump = 0, cl_KJ = 0, main φ == single-valued oracle to
  **3e-11**. Suite **235+8+2** (was 218+8+2; +17, some medium/M6 skip in CI — the B2 commit message's "229" was measured before the medium
  parametrization was added, corrected 2026-07-12); conforming solver paths byte-untouched; 8-thread cap. Next = B3 (implicit Kutta:
  g₁+g₂ wake LS replaces the weld).
- B3 — ✓ — 2026-07-12 — **B3 CLOSED (with B4).** Lifting solve with implicit Kutta on the level-set path:
  no Γ secant, no master–slave Γ — the TE jump is carried by the multivalued aux DOFs, the g₁+g₂ wake LS convects it, and its VALUE comes
  from B4's nonlinear TE pressure-equality Kutta. Γ is a RESULT. Delivered: `kernels/cut_assembly.py` (`mass_conservation_coo` per-side ρ
  per D10, `wake_ls_coo`, `te_kutta_coo`), `wake/multivalued.py` (`closure="wake_ls"`, side potentials/densities, TE control volumes),
  `solve/picard_ls.py::solve_multivalued_lifting`. Far field = freestream + vortex on the **MAIN** DOFs, aux **FREE**.
  Gate (`tests/test_b3_lifting.py`, **6 passed**): V3 M0.5 α=2° cl_KJ **0.2828** (medium) INSIDE the committed [PG 0.2788, KT 0.2919]
  bracket (read from `cases/reference_data/naca0012_m05/cl_reference.csv`), on BOTH the wake-embedded M0 and the **wake-free M3** families;
  same-mesh A/B vs conforming Γ within **0.1–0.7%** at M=0 and M=0.5 (0.1177/0.1191/0.1197 vs 0.1175/0.1200/0.1202 on coarse/medium/fine);
  the **wake-free** mesh (no `wake` tag, generic cuts — the workflow form) reproduces the embedded-mesh Γ to **0.3%**;
  the wake jump is CONVECTED, not decaying. **Five correctness fixes landed here (all load-bearing):**
  (1) far-field **aux DOFs must stay FREE** (Neumann) — pinning them to the vortex lower branch drains the circulation (jump decays
  0.0147→0.001); the vortex goes on the **main** DOFs only. (2) The wake must be **coplanar with the vortex branch cut** (chord plane y=0,
  design.md §4) — aiming the level set along the freestream while the branch cut stays horizontal leaves an unsupported Dirichlet jump at
  the outlet ⇒ spurious velocity ⇒ density blow-up (all high-M cells at x≈15, NaNs). (3) The per-side cut-strip density **limit-cycles and
  must be under-relaxed** (`omega_rho`, default 0.5); full adoption diverges after ~80 outers (Γ 0.126→0.010, M_max→6.7).
  (4) **D11 is mandatory**: wall Cp from `phi_main` makes lower-surface TE triangles reference the TE's UPPER value ⇒ cl_pressure =
  **−3.35** (junk); the per-side `phi_up`/`phi_lo` mapping brings cl_p within 0.4% of cl_KJ.
  (5) **Compressibility is carried by the BULK density, NOT the far-field vortex** — PG-scaling the vortex (β<1) leaves Γ unchanged, while
  the bulk density raises it 0.1086→0.1256 (the correct 1/β direction).
- B4 — ✓ — 2026-07-12 — **B4 NEW + CLOSED same day (user-directed) — TE control-volume / implicit-Kutta re-derivation.**
  The B3 blocker was that the emergent Γ converged to the WRONG value (0.2074/0.1760/0.1704 vs conforming 0.1175/0.1200/0.1202 —
  mesh-convergent ⇒ a METHOD defect, +42%). **Root cause, two structural facts:** (1) **the wake LS CANNOT pin Γ** —
  its residual is identically zero for any spatially-constant jump because Σ_c ∇N_c = ∇(1) = 0 (partition of unity), measured **1.9e-16** ⇒
  design_track_b.md §2.3/D2's "g₂ IS the discrete Kutta condition" is **FALSE and retired** (the López dissertation has no explicit Kutta
  anywhere — the word never appears in its method chapter); (2) Γ was therefore pinned by a single, WRONG equation —
  the TE aux row (lower-side mass conservation), whose control volume is up/down **asymmetric** on a symmetric airfoil (TE fan 9 upper / 6
  lower / 3 cut, because the ε shift sends every on-sheet node "+"). **★ Fix = the NONLINEAR TE pressure-equality (Bernoulli) Kutta.**
  Symmetrizing the control volume is NOT available (the mesh is naturally asymmetric at the TE —
  user-arbitrated 2026-07-12), so the condition is a POINTWISE PHYSICAL statement needing no symmetry: — q_u — ² = — q_l —
  ², factorized **exactly** as (q_u+q_l)·(q_u−q_l)=0 and linearized by freezing the mean s̄ = q_u+q_l at the previous iterate ⇒ a row LINEAR
  in φ, re-linearized each Picard outer (same cadence as the density lag — no new outer loop) and converging to the exact nonlinear
  condition. It replaces the TE **aux** row; the displaced lower-side mass-conservation entries are re-routed onto the TE **main** row,
  which then carries the TOTAL (upper+lower) balance, so mass stays conserved and no side is arbitrarily robbed of its equation.
  **Why it is non-degenerate where g₂ is not:** q_u and q_l are recovered on DIFFERENT element sets, so q_u−q_l is NOT a jump gradient and
  does not vanish for a constant jump. **★ The control volumes must be WALL-ADJACENT** (elements carrying a wall face = the upper/lower body
  surface at the TE), not the whole fan — the Kutta condition is about SURFACE velocities:
  full-fan recovery gives Γ 0.1407/0.1355/0.1329 (+11–15%, interior and wake elements pollute the average), wall-adjacent gives
  **0.1177/0.1191/0.1197** (**<1%** of conforming). **The D2 penalty-Kutta fallback is no longer needed** —
  no penalty weight, no tuning parameter (s̄ is solved for, not calibrated). Also fixed en route:
  the ε shift was manufacturing **spurious cuts in the below-TE fan** (3 of 6 elements got a bogus UPPER copy BELOW the wake) —
  exactly the López p.57 warning; that fix alone restored mesh convergence (Γ went from a mesh-independent wrong 0.186 to the convergent
  0.207→0.176→0.170). Gate (`tests/test_b4_te_control_volume.py`, **8 passed**, ~29 s): LS null space pinned (1.9e-16);
  TE control volumes wall-adjacent; below-TE fan never cut; **emergent Γ within 5% of conforming (measured 0.1–0.7%)** while the old
  `te_kutta="mass"` row is still >30% out; visual artifact `artifacts/EXPORT_TE_DIAGNOSIS/b4_te_kutta.png`.
  Interfaces: `solve_multivalued_lifting(..., te_kutta="pressure")` (default), `te_kutta="mass"` retained for the before/after contrast.
- B5 — ✓ — 2026-07-12 — (was B4.5, orig B3.5) **NEW 2026-07-11 + CLOSED 2026-07-12 (user-arbitrated) — far-field A/B:
  option a (Dirichlet+vortex) STAYS the default.** `solve_multivalued_lifting` grew `farfield="vortex"` (default, option a:
  spherical Dirichlet freestream + PG vortex on the MAIN DOFs with the emergent Γ refreshed in each outer iter, aux FREE) / `"neumann"`
  (option b, López: inflow Dirichlet freestream + outflow Neumann outlet carrying the freestream flux ρ∞(u·n̂), NO vortex, NO Γ feedback) /
  `"freestream"` (Dirichlet freestream everywhere, crudest). Helpers `_farfield_split`/`_neumann_outlet_rhs` in `solve/picard_ls.py`.
  **López-style domain-size re-calibration** (the dissertation §4.1.4 method) on BOTH NACA families (M0 embedded + M3 wake-free), coarse,
  M0.5 α2°, R ∈ {15,30,60,120}c: option a is **domain-robust** (Γ within 0.45%/1.09% of the truth over 15→120c;
  0.25% of conforming at 15c), option b truncates the **O(Γ/R)** point-vortex tail (−4.07% at 15c → −0.50% at 120c, halving each doubling of
  R ⇒ meets the B3 ±2% band only at **R≥~30c**, <1% at **R≥60c** = 2–4× larger domain), freestream crudest at every R (DIVERGES on compact
  15c M0). Both families bit-for-bit agree. ⇒ option a stays default (compact 15c workflow); option b validated but domain-hungry.
  **M6 leg folded into B7** (the 3D B-path solve is B7 machinery; the span-uniform option-a vortex also recreates the P5 branch-ray artifact
  on M6 without the Γ(z) taper — B7). Evidence: demo `cases/demo/b4p5_farfield/` (`farfield_domain_study.png` + summary/checks CSV,
  self-checking), `tests/test_b45_farfield.py` **10 passed** (15c coarse locks + `_farfield_split`/RHS unit checks).
  Conforming path byte-untouched.
- B6 — ◐ — 2026-07-12 — (was B5, orig B4) **Transonic + Mach continuation on the level-set path —
  coarse gate MET, medium fold = LS-Newton (delivered).** Full detail in the B6 gate section above (§"B6 —
  Transonic…") + design_track_b.md §10/§10.6. Delivered: per-side artificial density with a same-side-restricted upstream walk (D10;
  subcritical exact no-op), **supersonic-zone-localized damping** (the P4 whole-field θ·diag throttles the implicit-Kutta circulation —
  a Jacobi-smoother-vs-solution-mode effect), `solve_multivalued_transonic` (Mach ramp, **no Γ secant**), `post/surface_ls.py` (D11
  wall-Cp/shock). **★ Gate baseline changed (user-arbitrated): same-mesh conforming NEWTON truth, not the conforming Picard** (which
  under-circulates 4–8% at these shocks). **coarse M0.80 MET** dual-mesh (M0 Γ 0.2124/−7.9%, M3 0.2322/+0.9%, shock 0.644/0.678, 0 lim/flr;
  demo `cases/demo/b6_transonic/` 14/14). **★ Fold findings:** live option-a Γ→vortex loop-gain>1 near the fold ⇒ transonic recipe =
  `farfield="neumann"`; and the raw Picard-vs-Picard A/B gap is the conforming Picard's own stall bias (Newton-arbitrated).
  **★ LS Newton (`solve/newton_ls.py`, design §5.5/§10.6): DELIVERED + FD-verified 1.3e-9**, reaches machine-converged terminal-quadratic
  discrete **fold** solutions (0 lim/flr): coarse M0.80 M0 — R — 9.4e-13 / M3 3.2e-11; **medium M0.7875 M3 wake-free (workflow mesh) —
  R — 1.5e-12** — closing the "is it a solution?" question the Picard stall left open. **Two honest gaps (open):**
  M0-embedded medium live-Newton limit-cycles at 3e-6 (P8/N5 near-tie churn → wire in frozen selection);
  converged LS fold lift ~13% below conforming-Newton (discretization difference to apportion — mesh + B5 neumann −4% + cut-O(h);
  user decided NOT to chase it now). Tests `tests/test_b6_transonic.py` (9 + 2 gated) + `tests/test_b6_newton.py` (2 + 2 gated).
- B7 — ✓ — 2026-07-12 — (was B5.5, orig B4.5) **ONERA M6 3D gate — CLOSED, dual-mesh, first try.**
  Full detail in the B7 gate section above + design_track_b.md §11. M∞0.84/α3.06 coarse, `farfield="neumann"`, ramp 0.60→0.84 @ dm 0.04;
  **M1 embedded** cl_KJ 0.2765 / shocks 0.635/0.588/0.449 / Γ 0.1076→−0.0003 / M_max 1.453 / **0 lim,flr** (22.7 min) and **M4 wake-free**
  cl_KJ 0.2710 / 0.634/0.584/0.454 / 0.1055→+0.0003 / 1.368 / **0 lim,flr** (18.4 min); V6 1.77%/1.97% (P5 coarse 2.40%);
  dual-mesh A/B 2.0%. **★ The B6 lift INVERSION reproduces in 3D:** against the conforming **NEWTON** truth (cl_KJ 0.2692, the B6-arbitrated
  baseline) the LS Picard is **+2.7% (M1) / +0.7% (M4)** while the conforming **Picard** (P5, 0.24788) is **−8.6%** —
  the LS path has no early-stoppable Γ outer (implicit Kutta ⇒ Γ is a solution mode), so gating on P5 would penalise the B path for being
  closer to the truth; the **wake-free workflow mesh is the more accurate of the two**. **★ 3D far field = neumann, and the P5 Γ(z) taper is
  structurally UNNECESSARY on the B path** (not merely unimplemented): the B-path vortex is span-uniform with a y=0,x>0 branch cut at every
  z, which misfires two independent ways — (a) the α-aimed sheet is not coplanar with that cut (B3's rule in 3D;
  outlet M 0.958 vs neumann 0.513) and (b) even re-aimed coplanar, one scalar Γ cannot match Γ(z)→0 (P5's branch-ray artifact;
  outlet M 0.825) — and neumann carries no vortex, so neither can exist. **★ Γ(z) comes out spanwise-SMOOTH with NO smoothing** (unplanned
  finding): normalised RMS 2nd difference 0.0079/0.0091 vs the conforming P5's 0.0970 — **11–12× smoother**.
  The conforming path runs a separate secant PER TE STATION (the machinery whose single-station failure, st133, cost P5 a whole
  investigation, and whose jitter `INVESTIGATION_gamma_smoothing.md` failed to smooth away); the implicit Kutta has no per-station loop —
  Γ is ONE solution mode ⇒ the P5 spanwise-Γ problem is not fixed but made **structurally impossible**.
  **★ The 3D-only machinery needed NO new solver code** (B1's oblique-frame + spanwise-clip fixes held): Γ(tip) → ~3e-4 discretely;
  the only gap was post-processing (`post/surface_ls.py`: `section_cp_curve_levelset` + `cl_pressure_3d_levelset`).
  Cost far under the risk estimate (~0.6 s/outer at ~12k 3D DOFs ⇒ ~20 min/solve, not hours). **Caveats (recorded, not chased):**
  top Mach levels park on the P4/B6 Picard residual tail ( — R — ~4–6e-6, 600-outer cap) —
  bounded + physical + in band at every level, so the gate asserts *bounded*, not `converged`;
  **LS Newton on M6 deferred** (plain splu; P8/N6's true-3D LU fill ⇒ needs lagged-LU); shocks sit 0.02–0.04c aft of P5 (in band) and η=0.90
  is 0.087 aft of the P8 Newton shock. Evidence: demo `cases/demo/b7_onera_m6/` (**35/35 PASS**, 4 figures + summary/farfield/checks CSV),
  `tests/test_b7_onera_m6.py` (6 fast + 5 gated). Conforming path byte-untouched.
- B8 — ✓ — 2026-07-14 — (NEW 2026-07-13, user-approved; **CLOSED 2026-07-14 as CHARACTERIZED-NOT-CURED, user-arbitrated —
  B9 unblocked**) **Level-set tip-edge desingularization (row-blend tip taper)** — the LS analogue of P13/G13.2's conforming taper.
  The conforming `Γ_eff(z)=F(z)·Γ_Kutta(z)` cannot be ported: the LS path has no Γ DOF and its TE Kutta row `s·(q_u−q_l)=0` is homogeneous ⇒
  scaling by F is a no-op (G13.2 finding (8)). Fix = a convex BLEND per TE node of the pressure-equality row with B2's continuity weld:
  `F·[s·(q_u−q_l)] + (1−F)·[φ_aux−φ_main] = 0` (F=1 inboard ⇒ pressure Kutta bit-identical; F=0 at tip ⇒ weld ⇒ jump=0 ⇒ tip unloaded).
  A DIFFERENT model from `Γ=F·Γ_Kutta` ⇒ r_c independently calibrated, two-path comparison is a physics A/B. **★★ RESULT 2026-07-13 —
  the blend does NOT close the gate, because its PREMISE is wrong** (demo `cases/demo/b8_tip_taper_ls/` **12/12**;
  M6 coarse+medium, M0.5, no limiter). The blend works exactly as its model says — converges 0 lim/flr, unloads the tip circulation far past
  the criterion (**Γ_last ~ h^4.73**, criterion q≥1), perfectly LOCAL (inboard Γ +0.01%, cl +0.03%) —
  **and the tip edge still DIVERGES** (p **+1.341** untapered → +1.37…+1.58 tapered; bigger r_c is worse).
  **(1)** G13.2's `p ≈ 1−q` does NOT transfer (q=4.73 yet p=+1.37) ⇒ killing Γ_last does not kill the peak.
  **(2)** Lift cost ~0 because **there is nothing to unload** — the LS implicit Kutta already drives Γ(tip)→0 emergently (B7: ±3e-4);
  the conforming path needs the taper only because its free-edge rule leaves Γ_last ~ √h (q=0.44). **(3) ★ MECHANISM:**
  the peak cell is **OUTBOARD of the geometric tip** (z/b=1.0118), a **`beyond_tip` element the SPANWISE CLIP refuses to cut**, the SAME
  element tapered or not, and **NOT a small-cut sliver** (V 0.71× median, not even cut) ⇒ **the LS tip singularity lives in how the embedded
  sheet TERMINATES, not in the circulation it sheds.** The two paths' tip singularities are DIFFERENT OBJECTS.
  Machinery shipped (correct, tested, bit-identical by default) but **B8 needs a RE-SPEC aimed at the sheet termination (spanwise clip /
  beyond-tip zone) — user arbitration required.** Cost caveat: the LS path has **no AMG option** (hardcoded `spsolve`, B2 decision) —
  M6 medium is 484 s/solve at 67k dofs; **fine would hit the splu wall with no escape hatch**.
  **★★ RE-SPEC ROUND 2026-07-14 (diagnosis-first, then span blend — see the B8 section):**
  (1) the committed p=+1.34 was a **×5 METRIC ARTIFACT** (`element_mach2` reads mixed-side plain/beyond-tip elements from the
  aux-substituted side field; assembly uses MAIN dofs — elem 93977: side 1.532 vs main 0.309);
  the **HONEST exponent is +0.62 (+0.37 no-sliver) = the SAME object as the conforming +0.52**;
  fix opt-in `element_mach2(mixed_plain="main")`, default bit-identical. (2) The honest residual object is the **termination ring's FINITE
  jump** ( — δ — ≈0.026, h- and TE-taper-independent — decoupled from Γ_last, which is why the TE blend measured nothing).
  (3) The **span blend of the wake-LS rows** (`MultivaluedOperator(span_blend=…)`, default None bit-identical, 11 tests, B-suite 116/9)
  **WELDS the ring (0.026→0.0003) but is NEGATIVE on locality: ~20% GLOBAL lift loss, uniform in z, r_blend-insensitive, h-GROWING** (⇒ its
  flat p at rb0.08 is confounded, corrected ~+0.15); component isolation: straddler weld alone −13.3%, inboard blend alone −10.8% ⇒ **the
  implicit Kutta has no per-station target — ANY sheet-side δ-pin re-levels the global Γ mode ~10×** (the conforming secant keeps the same
  F(z) at −1.6%). Demos `run_b8_termination_diagnosis.py` + `run_b8_span_blend.py` (8/8). **Both constraint-side routes now measured dead;
  any further cure must change the FUNCTION SPACE at the termination. ★ ARBITRATED 2026-07-14:
  CLOSED as characterized-not-cured (honest exponent = the conforming object every closed gate lives with) ⇒ B9 UNBLOCKED;
  backlog **EXECUTED 2026-07-14 (user-directed)**: `element_mach2` default flipped to `mixed_plain="main"` ("side" stays opt-in for
  reproducing committed diagnostics; demo repro scripts pin it explicitly), B6/B7 M_max re-read from the cached states WITHOUT re-solving
  (`run_b8_mmax_reread.py` + `mmax_reread.csv`: side values reproduce committed to 6 digits ⇒ reconstruction verified;
  **M1 1.453 side → 1.392 main** — the committed M_max was itself a beyond-tip artifact cell; M4 + both 2.5D states bit-identical;
  all gate bands unchanged), G13.1/P9 LS-exponent errata placed at the original claims in roadmap/demo_report, and the M2 LS-ingestion
  census landed as tests+CSV (`test_m2_wingbody.py` +7, `ls_ingest_census.csv` — the f3c7989 prose numbers 1,415/76 confirmed exactly at
  α=0, medium 29,108/150 added). Still recorded, NOT scheduled: the `element_densities` mixed-plain junk-weight fix.**
- B9 — ✓ — 2026-07-17 — (was B8 2026-07-13; orig B6→B5) **★★ RE-SPEC + CLOSED 2026-07-17 (user-approved):
  wing-body cross-model validation (LS + conforming), M∞ 0.5, α 3.06°, M2 coarse+medium.** Multi-element leg SUPERSEDED. ★ HEADLINE:
  the two wake models AGREE to **cl_p 0.4% / cl_kj 0.6%** at medium (conf-pressure 0.2173/0.2188 vs LS-Picard 0.2165/0.2175, GB9.5 PASS <
  1%; coarse 12.8% = resolution) — the wing-body analogue of the wing-alone P14 cross-model.
  Conforming leg = NEW capability (`onera_m6_wingbody_mesh(embed_wake=True)`: split fuselage two π-revolves + through-body sheet + Netgen
  off; `cut_wake`/constraints/P14-Kutta unchanged; embed_wake=False bit-identical). GB9.1 ✓ GB9.2 ✓ GB9.3 ✓ (junction TE-CV, tests) GB9.5 ✓
  GB9.6 ✓ RECORDED; **GB9.4 XFAIL** (fuselage lift 16-20%, resolution/model-sensitive ⇒ G1.6 fuselage-Cp error, band NOT moved).
  ★ LS uses PICARD not Newton — the committed LS Newton recipes (Schur/freeze/lagged-LU/ramp) all diverge on the wing-body because `neumann`
  is unbounded under the fuselage blockage; diagnostic: te_aux perfect (1.8e-8), 8 far-field fluid rows — R —
  ≈84 in the freestream-Newton path (never exercised). Demo `cases/demo/b9_wingbody/` 7 PASS + 1 XFAIL;
  guardrail `cases/analysis/b9_fuselage_guardrail/`. Closes the solver leg of Track M's M2.
- B10 — ⊘ SHELVED — 2026-07-10 — (was B9 2026-07-13; orig B7→B6) Curved wake / free wake. Recorded reasons (DN1 §8 / DN2 §4.5.6):
  the loading error of a straight wake is O(θ²) ≈ 0.1%; per-update CutElementMap/DOF rebuild cost;
  discrete cut-set jumps conflict with Newton; López precedent. The `update_direction()` interface capability is retained —
  it is what B1's α re-aim tests exercise.
- B11 — ✓ — 2026-07-14 — (NEW 2026-07-14, user-directed; appended after B10, no renumber) **LS-path infrastructure:
  unified post-processing + GMRES/AMG scaling.** Two gaps closed (a B9 enabler). **(1)** `post/surface.py` + `post/surface_ls.py` now share
  private cores (`_cp_from_q2`, `_pressure_force`, `_wall_plane_crossings`/`_resolve_station`/`_section_curve_dict`, `_d11_wall_state`)
  under a keyword-dispatched upper layer `post/unified.py` (`wall_cp`/`wall_forces`/`section_cp`, `phi=` conforming vs `mvop=,phi_ext=`
  level-set); every legacy function keeps its name/signature and outputs are `np.array_equal` (D11 lock + shock locks pass unchanged).
  **(2)** the deferred design_track_b.md §5.3 GMRES+AMG landing: `solve_multivalued_laplace`/`_lifting`/`_newton` grow
  `precond=None — "ilu" — "amg"` (None = the bit-identical `spsolve` default; transonic inherits via `**kwargs`), the escape from the
  M6-fine splu wall (roadmap "no precond option" caveat). **★ ILU is the effective escape** (spilu on the real fused matrix, 434 iters
  coarse, exact); **AMG (SA on an SPD Picard-block surrogate + aux↔host springs) converges only on the SPD Laplace/continuity system —
  on the `wake_ls` lifting operator its convection-like aux rows defeat the SPD surrogate and GMRES STALLS (measured:
  γ 0.0033 vs 0.139, 455 s, all outers stalled)**, so AMG stays a Laplace/§5.3 knob and ILU is shipped.
  Núñez symmetric row assignment stays not-prebuilt (§5.3). lagged-LU (`direct_refactor_every`) port to `newton_ls` = recorded out-of-scope
  follow-up → **executed by B12 (2026-07-14)**. Evidence: `tests/test_b11_post_unified.py` (9) + `tests/test_b11_linear_ls.py` (10 + 1
  gated); demos `cases/demo/b11_ls_infra/` (unified-post + GMRES A/B + gated M6-medium headline CSV). Conforming numerics byte-untouched;
  no Numba/COO path touched.
- B12 — ✓ — 2026-07-14 — (NEW 2026-07-14, user-directed; appended after B11, no renumber; executes the B11/G11.4 follow-up) **Lagged-LU
  direct-reuse for LS Newton (medium/M6-scale enabler).** B11 measured that the iterative escapes fail beyond coarse (ILU diverges at 2.5D
  medium lifting, `factor_failed`s at M6 medium; AMG stalls), so at medium/M6 sparse-direct is the only converging tool and the cost driver
  is the NUMBER of factorizations (17.5 s each at 67k dofs). `solve_multivalued_newton` gains `direct_refactor_every` (default 1 =
  bit-identical per-step `spsolve`) + `direct_reuse_rtol`: with `k>1` it refactors the LU every k-th Newton step and drives the intermediate
  steps with GMRES preconditioned by the stale (exact) LU — the N6 mechanism (`solve/newton.py`) ported **minus the Woodbury** (the LS
  system has no Γ DOF ⇒ plain `J_free d = −R_free`). **G12.1 (bit-identity) ✓ + G12.2 (equivalence/reuse) ✓** —
  coarse M0.70 k∈{2,1000} reach the spsolve γ (0.1778053693) to bit-identity, 0 stalls, k=1000 refactors ONCE over 6 Newton iters.
  **G12.3 (M6-medium subsonic A/B) ✓** — M6 medium M0.5 (67,426 dofs), 7 Newton steps both, spsolve refactors 7× (**145.6 s**) vs lagged-LU
  1× + 30 reuse-GMRES iters (**66.7 s = 2.18×**), γ bit-identical ( — Δγ — 6.7e-13), 0 stalls, 0 lim/flr. Honest boundary:
  a real medium-scale win (one splu fits at 67k), but does NOT break the FINE memory wall (still needs ≥1 in-memory splu;
  that's the Núñez→AMG route). `solve/newton_ls.py` is the only production change, default byte-identical.
  Evidence: `tests/test_b12_lagged_lu_ls.py` (4); demo `cases/demo/b12_lagged_lu/` (6/6).
- B13 — ✓ — 2026-07-14 — (NEW 2026-07-14, user-directed; appended after B12) **Lagged-LU on the Picard OUTER loop** —
  the post-B12 cost driver (one 17.5 s spsolve per outer; B11 lifting headline 447.6 s / 26 outers, Newton seed 263 s / 15).
  `solve_multivalued_lifting` gains `direct_refactor_every` (default 1 = bit-identical) + `direct_reuse_rtol` (**1e-10, NOT B12's 1e-8** —
  a Picard fixed point is pinned only by its 1e-6 lag tolerances, so an inexact reuse step SHIFTS the stopping point, measured — Δγ —
  8e-8 at 1e-8; Newton's terminus is pinned by tol_residual regardless); transonic inherits via `**kwargs`;
  laplace excluded (single-shot). User goal arbitrated: medium-scale speed is the objective, fine optional ⇒ this outranks the structural
  preconditioner (B14). **GB13.1 ✓ + GB13.2 ✓ + GB13.3 ✓** — M6-medium lifting **447.6 s → 68.3 s (6.55×)**, 2 refactors vs 26 outers, γ
  bit-identical ( — Δγ — 6.9e-13); end-to-end seed+Newton **~330 s → 111.9 s (~3×)**, seed 263→42 s.
  (1 GMRES stall = the designed safety-net refactor on an early large-density outer, not a divergence.)
  External-doc corrections recorded (Schur direction inverted in both external docs; 454.8 s = 26 outers not one splu).
  Evidence: `tests/test_b13_lagged_picard.py` (5); demo `cases/demo/b13_lagged_picard/` (6/6).
- B14 — ✓ — 2026-07-17 — (OPENED + CLOSED 2026-07-17, user-directed) **Schur-eliminated aux block + AMG(SPD Picard main block) —
  the structural preconditioner, BUILT.** `precond="schur"` on both LS drivers (`pyfp3d/solve/schur_ls.py::SchurReducedSystem`):
  per step/outer eliminate the SMALL aux thin-strip block exactly (`K = J_mm − J_ma·J_aa⁻¹·J_am`, `lu_aa = splu(J_aa)`, n_ext-sized —
  1004/3701 dofs at M6 coarse/medium, split+splu **≤19 ms**), GMRES on the reduced main-free operator preconditioned by AMG on the SPD
  single-valued Picard block, **NO springs** — the B11 surrogate's jump≈0 bias (γ 0.0033 vs 0.139) is structurally absent, and the
  circulation mode survives. Shared by `solve_multivalued_newton` + `solve_multivalued_lifting` (transonic wrappers inherit via kwargs;
  `solve_multivalued_laplace` out of scope). A stalled reduced GMRES falls back to a full fused spsolve in the same step
  (`n_schur_fallback`) — **0 fallbacks anywhere in the campaign**. **GB14.1 ✓** J_aa factors on all 4 measured cases, cond1
  5.1e8/8.2e9/6.5e6/**7.4e7** (finite; measured, not assumed). **GB14.2 ✓** 2.5D coarse lifting+Newton land on the spsolve γ to \|Δγ\|
  4.2e-11/2.0e-12, 0 fallbacks — on the exact operator where the B11 surrogate stalled to γ 0.0033.
  **GB14.3 ✓** the pre-registered discriminating tier — 2.5D MEDIUM lifting, where ILU DIVERGED to γ=−136.99 (77 stalls) —
  schur converges to γ 0.14137632, \|Δγ\| 9.3e-10, "a real escape". **GB14.4 ✓** 3-D capability:
  M6 wake-free COARSE + MEDIUM × M0.5 lifting + M0.84 ramp all converge/target-reached, γ matches the lagged-LU arm to \|Δγ\| ≤ 1.5e-8 and
  the committed GB15.4 state exactly (γ **0.088338**, M_max **2.4938**). **GB14.5 ✓** default `precond=None` byte-identical,
  `n_schur_fallback` inert. ★ **TIMING (RECORDED, not gated — the design said medium-scale gain was uncertain;
  it landed on the winning side):** M6 medium M0.5 lifting **73.2 → 35.2 s = 2.08×** (precond 51.7% → 5.2%);
  M6 medium M0.84 ramp **671.2 → 469.3 s = 1.43×** (precond **42.6%/43.6% → 2.6%**, i.e. the A1 bottleneck is GONE —
  beats the user's <10% target). ★ **Honest limit:** at SMALL scale schur is SLOWER (2.5D coarse/medium and M6 coarse:
  the direct solve is already trivially cheap, the extra Krylov iters cost more than the tiny factorization) —
  the speedup appears only at M6-medium size and grows with it, so the **remaining designed value is the fine memory-bounded path (AMG O(n)
  + thin-strip LU, no full-size splu), out of scope here (user: coarse+medium)**. Fallbacks (block-triangular;
  Núñez additive) NOT needed. Evidence: `tests/test_b14_schur_ls.py` (9, incl. gated GB14.3);
  demo `cases/demo/b14_schur_precond/` (**7/7 incl. gated M6 coarse+medium**).
- B15 — ✓ — 2026-07-15 — (NEW 2026-07-15, user-directed; appended after B14, no renumber) **LS Newton transonic ramp + N5 freeze-selection —
  the Picard shock-position PLATEAU is removed.** Root cause measured: the LS transonic Picard's top Mach levels never converge and burn
  their full outer budget on the plateau (M6 medium M0.84: levels 0.80/0.84 each run all 200 outers) — that IS the 24.5/38.4 min.
  Newton has no such soft mode, but `newton_ls` could not ramp (`freeze=` was a reserved no-op;
  the 0-clamped gate blocks shock limiter cells; no Mach-ramp wrapper). Delivered: `MultivaluedOperator.newton_side_data(frozen=…)` +
  `freeze_side_state` (the `kernels/upwind.py` frozen apparatus reused UNMODIFIED — the per-side ops are already walk-mode with a
  side-masked graph ⇒ wiring, not new numerics); `LSNewtonSystem` (residual+Jacobian in ONE code path shared with the FD gate);
  `solve_multivalued_newton_transonic`. **GB15.1 ✓** FD rel **6.7e-9**, frozen sweep reproduces live density BITWISE at the freeze point.
  **GB15.2 ✓** the freeze cures a genuine **period-6 limit cycle** (NACA coarse M0.75: live stuck at 2.7e-7 with 0 lim/flr = clean
  assignment churn → frozen **22 steps to 8.5e-13**, 0 reverts, γ 0.218809 vs the live cycle's 0.218804 ⇒ it removes churn, it does NOT move
  the solution). **GB15.3 ✓** NACA coarse M0.80/α1.25 (B6 gate; conforming-Newton truth M_max 1.408): Picard 41.9 s → — R —
  1.55e-5 with only **3/5 levels converged (not a solution)** vs Newton **7.5 s → 3.1e-12 strict (5.6×)**, M_max 1.3924 (−1.1% of truth);
  **+`intermediate_tol` 6.5 s with γ 0.212445 IDENTICAL to strict** (48→38 steps) ⇒ the loose-intermediate knob is FREE.
  **GB15.4 ✓** M6 medium M0.84 wake-free: committed Picard **2304.7 s (38.4 min)** on the 1e-5..1e-4 plateau → Newton ramp **657 s (11.0
  min) = 3.51×, ALL 6 levels converged to ~1e-11**, freeze armed everywhere, 0 reverts; M_max 2.4938 vs Picard 2.4549 (1.6%), 3 clamped
  cells of 330k vs ≤3. ⇒ closes the deferred **B6-medium quantitative** + **B7-quantitative** items. ★ **HONEST LIMIT:**
  5 of 6 levels accept via `assignment_cycle` — the FROZEN system converges to ~1e-11 and is accepted at the **assignment-discontinuity
  floor** (the N5 semantics the conforming path also uses); this is NOT a claim that the LIVE residual is <1e-10.
  6–7 orders better than the Picard plateau, but 'live-strict' would be an over-claim. Also open (unchanged by B15):
  the LS-vs-conforming **discretization gap** (γ −7.4% of the same-mesh conforming-Newton truth at NACA coarse M0.80;
  B6 recorded ~13%) — B15 makes it measurable strict-to-strict for the first time, it does not close it.
  ★★ **FOUR ERRATA (the conforming N5 recipe is NOT mechanically portable — the B8 lesson again; all four forced out by measurement):**
  (1) the TE polyline must come from the AUTHORITATIVE geometry — hand-rolled x_te off by **2e-4** matches ZERO wall nodes ⇒ **0 TE nodes ⇒
  no Kutta ⇒ 340k limited cells + NaN, passed SILENTLY** ⇒ both LS solvers now RAISE on `te_nodes == 0`;
  (2) **`freeze_tol` must sit ABOVE the CHURN FLOOR, which RISES with Mach** (<1e-6 @M0.60 → 8.6e-6 @M0.65 → **2.7e-4** @M0.70) —
  below it a discrete selection flip throws the residual back before the freeze can arm and the ramp DIES at M≈0.66 (**same law as
  'tol_residual above the Picard plateau'**); (3) **residuals are NOT comparable across a SELECTION EPOCH** —
  the frozen phase drives `r_best` to 1.5e-11, a refresh legitimately returns it to the live scale (2.6e-3), and the fail-fast reads a 1e8×
  blow-up and kills a healthy freeze-refresh cycle ⇒ `r_best` reset on freeze/refresh/revert;
  (4) **the frozen clamp count is STALE BY CONSTRUCTION** (`n_floored` = `branch==3` at the freeze point, never falls) ⇒ gating on
  `n_flr==0` **refuses a 7.8e-14 machine-precision solution forever** ⇒ the **LIVE re-evaluation** is the arbiter.
  New knob **`freeze_max_clamped`** (default 0 = the conforming rule): a SINGLE floored cell of 330k otherwise blocks the freeze at ANY
  `freeze_tol` — **the P9/G9.1 wall** — yet the frozen sweep represents a clamped cell exactly (branch 3), so the precondition was stricter
  than the machinery needs; relaxed, the ramp completes **WITH** the clamped cells (⚠ errata 2026-07-15: they PERSIST —
  the converged M0.84 state carries 3 of 330k, and the `assignment_cycle`/`refresh_budget` accept routes no longer re-check the clamp count;
  see the corrections block in the B15 entry). Defaults byte-identical (`freeze_tol=None`). Evidence:
  `tests/test_b15_ls_newton_freeze.py` (12, incl. the errata lock); demo `cases/demo/b15_ls_newton_ramp/` (**19/19 incl. gated M6**).
  ★★ **Erratum trail (2026-07-19): under the first B20 re-baseline the "reaches M0.84" clause went NEGATIVE** (M0.6625, γ 0.071909, demo
  17/20; GB20.7 called it a real capability loss) — **then B21 found the true mechanism (the B20 patch gap in `freeze_side_state`, Kimi N1)
  and RESTORED it: M0.84 again, γ 0.088343, res 9e-14, 0 lim/1 flr, 515 s**. The capability clause stands;
  numbers re-based slightly (demo refresh = recorded follow-up). The M_max 2.4938 cross-check remains common-mode-only (conforming record
  1.995, different mesh family).
- B16 — ✓ — 2026-07-17 — (NEW 2026-07-17, user-directed; appended after B15, no renumber; executes the B9 recorded follow-up) **LS Newton
  far-field BC generalisation — far-field aux-DOF pin.** Root cause (GB16.1, MEASURED — B9 had it as prose only):
  a wake level set has no outflow clip ⇒ the sheet reaches the far field and the outer nodes it crosses carry aux DOFs governed only by
  **near-singular wake-LS rows on giant outer tets**; at the freestream Picard state they hold garbage (coarse — jump —
  **53.4** at x≥10 vs Γ̄ 0.0586), which Picard's fixed point absorbs but the Newton residual reads as an O(1) inconsistency —
  **the 8 far-field MAIN rows, max\|R\| = 84.457** reproduced to the digit; aux-block cond1 **9.1e18** (legacy, ABOVE the 1e14 GB14.1
  ceiling; doc erratum 2026-07-19: 6.36e18 was a pre-CSV trial value) **→ 8.70e6** (pin). ⚠ **The proposal's mechanism was wrong** (Picard
  lifting uses `wake_ls` too, not the weld — the difference is fixed-point absorption vs Newton residual, not the closure).
  **Fix = `farfield_aux="pin"` (default), mode-adaptive (user-arbitrated):** on a Dirichlet far field the far-field-boundary aux enter the
  Dirichlet set at the host's branch value (freestream → φ∞, jump→0; vortex → main−side·γ, jump→γ, the conforming `lower_branch_mask`
  analogue, NOT the B3 both-sides pin); `neumann` byte-identical. Helper `farfield_aux_dofs` (`solve/picard_ls.py`);
  Schur adapted (`n_aux_expected = n_ext − n_pinned`, assert kept). **GB16.1 ✓** (diagnostic + D8:
  pin drives outer jumps 53.4→**5.3e-15** and ★ cures the 4 INTERIOR junk aux too ⇒ R2 risk void).
  **GB16.2 ✓** neumann `array_equal` pin vs legacy; B12 γ 0.06685284 / B15 γ 0.088338 anchors gated.
  **GB16.3 coarse ✓** legacy churns (res **7.95**, 3690 limited) → pin **res 5.88e-14, 0 limited**, and the coarse converged lift matches
  conforming (cl_p 0.2086 vs 0.2089, 0.1%); ⚠ pin carries `n_flr=3` at the wing-fuselage junction = the **B8 mixed-plain / G1.6
  fuselage-Cp** class (same root as **GB9.4**), orthogonal to the BC fix. **★★ GB16.4 XFAIL —
  ~~the OPEN non-convergence~~ RESOLVED BY B17 2026-07-18 (it was a freestream-pin BC-modelling error, jump=0 kills the outflow circulation;
  see the B17 row):** the {Newton-pin, LS-Picard, conforming} lift triangle does NOT close and FLIPS with resolution — coarse:
  Newton-pin 0.2086 ≈ conforming 0.2089, Picard 0.1853 low; medium: Picard 0.2165 ≈ conforming 0.2173 (B9's headline), Newton-pin 0.1690 low
  (22%, STALLED at res 7e-6) ⇒ at least one path is not converged. Two live possibilities (neither ruled out):
  the medium Newton-pin is non-converged (a warm start from the converged Picard also failed to converge in ~10 min ⇒ not merely a shallow
  seed), OR the B9 LS-Picard≈conforming 0.4% was itself a non-converged coincidence. UNRESOLVED, analysis deferred;
  the churn fix stands on the coarse machine-converged evidence, the "Newton now matches the other paths" claim does NOT.
  **GB16.5 ✓** Schur split with pinned aux constructs + fails loudly on the legacy count; `test_b14` green.
  **GB16.6** transonic stretch RECORDED (gated). ★ New knob `farfield_aux` default `"pin"` —
  defensible-as-default (freestream Newton never exercised pre-B9, vortex Newton zero committed recipes, neumann byte-identical ⇒ "default
  leaves every committed anchor bit-identical" holds vacuously); `"legacy"` is the pathology reproduction switch.
  Evidence: `tests/test_b16_farfield_aux.py` (9, incl. gated GB16.3); demo `cases/demo/b16_farfield_aux/` (5/5 coarse PASS + gated medium).
- B17 — ✓ — 2026-07-18 — (NEW 2026-07-18, user-directed; appended after B16, no renumber; executes the B16 GB16.4 open follow-up)
  **Far-field aux pin carries jump=γ, not 0 — resolves GB16.4.** GB16.4 was NOT a non-convergence:
  the B16 freestream pin forced the outflow wake jump to **0**, which REMOVES the circulation the wake physically carries out ⇒ a
  resolution-dependent lift error (invisible at coarse, where jump=0 cancelled the legacy outer-tet garbage;
  −22% at medium, where legacy already carries the jump). **Decisive discriminator (GB17.3):**
  giving the *Picard* driver the same freestream pin (new `farfield_aux` knob on `solve_multivalued_lifting`) makes medium Picard-pin
  converge cleanly (res 7.5e-8) to cl_p **0.1691** — matching the "stalled" Newton-pin 0.1690 to 0.1% ⇒ two independent solvers on the same
  value ⇒ a genuine BC-determined state, not a Newton stall. **Fix = `farfield_aux="pin_gamma"` (jump→γ, the new default):**
  aux = host φ∞ − side·γ, refreshed with the live γ — same near-singular-aux Dirichlet cure, physical ring value.
  Triangle then closes MONOTONE to conforming (cl_p wing): coarse conf 0.2089 / legacy 0.1853 / pin0 0.2086 / **pin_gamma 0.2087**;
  medium conf 0.2173 / legacy 0.2165 / pin0 0.1690 / **pin_gamma 0.2117 (Picard) = 0.2115 (Newton)** —
  both solvers agree 0.1%, undershoot conforming 0.1%/2.6% (far-field truncation). **GB17.1 ✓** ring jump collapses legacy 53.4 → pin 0 →
  pin_gamma 0.063=γ; coarse legacy garbage ( — jump — =53) IS a 12% deficit ⇒ **B9 "coarse 12.8% = resolution" was contamination** (erratum
  added). **GB17.2 ✓** cl_p (surface integral) and cl_KJ (circulation) move together ⇒ the gap is a flow-state change NOT a post artifact;
  the user's "Cp aligns yet cl_p differs" is a Cp-axis scale illusion, and the plotted sectional cl(z) is Γ-based `2Γ/(u·c)`.
  **GB17.3 ✓** pin jump=0 = BC error, both solvers 0.169. **GB17.4 ✓** pin_gamma closes triangle (Newton≈Picard <1%, monotone).
  **GB17.5** spanwise Γ(z) uniform offset removed (RECORDED). **GB17.6** vortex does NOT close the gap —
  BRACKETS conforming from above (medium +2.5% vs pin_gamma −2.6%) and churns at coarse (free aux);
  freestream pin_gamma stays recommended (RECORDED, user-requested). ★ **B16 conflated two orthogonal issues** —
  far-field near-singular *conditioning* (pin cures, jump value irrelevant) vs outflow *circulation* (needs jump=γ);
  the wing-fuselage-junction churn (nlim 42/nflr 40 at medium, G1.6/GB9.4 class) survives but limits only the residual floor, not the lift.
  **Defaults (user-arbitrated):** `pin_gamma` new default on BOTH solvers, acts only on freestream, inert (bit-identical to legacy) on
  vortex/neumann ⇒ every committed 2.5D vortex/neumann Picard + neumann Newton anchor byte-untouched;
  B9/B16 freestream Picard demos pinned to explicit `legacy`; B16 jump=0 reproduces with explicit `"pin"`.
  Evidence: `tests/test_b17_farfield_pin_gamma.py` (6, ungated); demo `cases/demo/b17_farfield_pin_gamma/` (3 coarse PASS + gated medium).
  ★ **B20 re-baseline (2026-07-19):** medium Newton 0.2115 → **0.2114**; the junction churn WAS the mixed-plain contamination —
  the re-baselined medium Newton trajectory converges to ** — R — ~1e-13** (γ 0.064201) with the clamps gone, so "limits only the residual
  floor" is superseded.
- B18 — ✓ — 2026-07-18 — (NEW 2026-07-18, user-directed; appended after B17, no renumber; executes the GB16.6 debt) **Wing-body transonic
  (M0.84) — conforming reaches it, level-set is junction-limited.** ★★ **[superseded 2026-07-20 by B26/B27 —
  the "junction-limited" story is RETIRED: pocket = the B23 free-edge singularity, healed by B25;
  the post-cure LS ceiling is co-located with conforming; demo refreshed 8/8 PASS. Conforming legs stand.]**
  The capability is asymmetric and that IS the finding. **Conforming** (Newton + pressure Kutta, Mach continuation) is the wing-body
  transonic path: **coarse M0.84 cl_p 0.2617** (Mmax 2.15, strict), **medium M0.79 cl_p 0.2579** strict (res 2.2e-14), clean cl_p(M) rise
  **0.2173/0.2321/0.2579** at M0.50/0.65/0.79; medium **M0.80+ stalls** (res ~2e-6, 0 clamp — NOT slivers:
  medium mesh clean 0-tets<5°, coarse has 27 yet reaches 0.84; a sharper shock/junction interaction, recorded not chased). ★ recipe:
  the conforming wing-body medium ramp needs `freeze_tol` raised to the wing-body churn floor (1e-6→1e-5, the B17 lesson) or it stalls at
  M0.80. **Level-set** (B15 freeze-ramp + B17 pin_gamma) does NOT reach transonic on the wing-body:
  the wing-fuselage junction spurious supersonic pocket (**G1.6/GB9.4/B8 mixed-plain**, M²≈1.27 already at M0.5) **WORSENS with refinement**
  — coarse ceiling **M0.575** (Mmax 1.44), medium dies at the FIRST transonic level ~**M0.5** (Mmax artifact 3.96, nlim 43/nflr 40);
  the direct analogue of GB9.4's fuselage-lift-grows-with-refinement. Closed-negative discretization error (discipline #8), characterized
  not chased. ⇒ **no common transonic Mach at medium** (LS can't leave 0.5), so the trustworthy cross-model stays M0.5 (2.6%, B9/B17);
  a coarse M0.60 transonic cross-model was skipped (LS coarse ceiling 0.575<0.60). **GB18.1 ✓ PASS** (conforming transonic).
  **GB18.2/3/4/5 RECORDED** (LS ceiling; cross-model M0.5-only; junction worsens with refinement coarse Mmax 1.4→medium 4.0;
  fuselage lift 16% of wing @M0.79, GB9.4 class). ★ **repays the GB16.6 evidence debt** (spec'd RECORDED but never implemented;
  B18 executes it as a negative). ★ **NO pyfp3d/ numerics change** — pure demo/tests/docs on existing `solve_newton_transonic` +
  `solve_multivalued_newton_transonic`. fine excluded (G13.3). Evidence: `tests/test_b18_wingbody_transonic.py` (4, ungated);
  demo `cases/demo/b18_wingbody_transonic/` (7 gates: 1 PASS + 6 RECORDED). ★★ **B20/GB20.5 erratum (2026-07-19):**
  the junction-pocket attribution "B8 mixed-plain" is measured **FALSE** — it is the **G1.6 faceted-geometry** error (removing the
  contamination CONVERGED the medium case, res 6.8e-5 → 1.1e-13 / clamps 82 → 6, and UNCLAMPED a genuine **Mmax 5.22**;
  the old 3.96 / nlim 43 was a clamped artifact). Post-B20: coarse ceiling **~M0.55 (Mmax 1.31)**;
  and the coarse M0.6 cross-model point, originally "skipped", now EXISTS in the re-baselined artifact (conf 0.2178 vs LS 0.2174 = 0.2 %).

