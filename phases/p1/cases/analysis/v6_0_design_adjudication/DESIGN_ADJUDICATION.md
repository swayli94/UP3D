# GV6.0 — Design adjudication: the V6 wake-sheet IBL correction (BEFORE code, user-adjudicated)

Registered in `docs/roadmap/track_v.md` (the V6 section): *"GV6.0 design
adjudication (BEFORE code, user-adjudicated): the LS-path sheet-source
mechanism. `pyfp3d/wake/` has no sheet-surface integration machinery — the
options are zero-isosurface polygon integration (new geometry code) or a
volume-band approximation (deviates from the 'sheet RHS' formulation). The
conforming path needs NO new mechanism (explicit `wake_minus/plus` faces +
slave→master folding IS the weak-form flux channel). V6 may close
conforming-only; the LS leg then becomes a recorded follow-up."*

This document is the adjudication request. It states the code-survey facts,
prices the options, and asks for the user's ruling. **No solver/library code
is written before the ruling** (the GV6.0 ordering constraint).

Date: 2026-07-25. Author: Kimi (Track V V6 opening).

---

## 1. Code-survey facts (committed-state, with file:line)

### 1.1 The conforming channel exists but is currently unconsumed

- `mesh/wake_cut.py:408-414,434-435` — `cut_wake()` splits the original
  `wake` face group into `boundary_faces["wake_minus"]` (original node ids)
  and `boundary_faces["wake_plus"]` (slave ids), also stored as
  `wc.wake_faces_minus/plus`; topology assert gives each face exactly one
  owner (`wake_cut.py:586-591`).
- **No consumer exists.** Grep across `pyfp3d/`, `tests/`, `cases/*.py`:
  outside `mesh/wake_cut.py` itself nothing reads `wake_minus/plus`. Kutta/TE
  conditions go through probe nodes (`constraints/wake.py:196-212`) and wall
  faces (`constraints/te_pressure.py:108`); the weak mass-flux continuity is
  carried implicitly by the master–slave fold T. The registered phrase
  "slave→master folding IS the weak-form flux channel" is true of the fold
  itself, but the explicit faces are **data without a consumer** — GV6.1
  connects the two halves; it does not switch on a live channel.
- The RHS channel is threaded through all three drivers and gated:
  - Picard: `solve/picard.py:609-617,659` (`b_zero` → `con.reduced_rhs(b, Γ)`)
  - Newton: `solve/newton.py:372-378` (`external_rhs` → R_red = TᵀR_vol)
  - Tight: `viscous/tight.py:10-12,40-43` (R_red = Tᵀ(R_inv − b(m(U))))
  - Gate precedent: GV2.1 (`tests/test_v2_transpiration.py:247-267` exercises
    a nonzero RHS through exactly this Tᵀ route on a wake-cut mesh, with the
    GV2.1(a) analytic sign pin and the GV2.1(b) zeros bit-identity).
- `viscous/transpiration.py:66-114` `assemble_transpiration_rhs(nodes, faces,
  m_dot)` is face-group-agnostic (docstring: faces may be any cut-mesh face
  group, TE-duplicated slaves included) — feeding it
  `wc.wake_faces_minus/plus` yields b_wake with no new kernel; slave-node
  loads fold into master rows automatically under Tᵀ (that IS the sheet
  flux). `transpiration_from_delta_star` (`transpiration.py:165-174`, δ*→ṁ =
  div_Γ(ρ_e u_e δ*)) is likewise SurfaceMesh-generic. One caveat: the
  function's sign convention is pinned to "blowing out of the body" (module
  docstring `transpiration.py:19-27`, the negated load at `:114`); the wake
  sheet is an internal slit loaded from both sides, so the sign for the
  sheet source needs its own analytic pin (GV2.1(a)-style MMS, new test).

### 1.2 The LS path genuinely has no sheet-integration machinery

- `pyfp3d/wake/` = `levelset.py` (implicit ruled-surface (s,d,q) evaluation +
  normals), `cut_elements.py` (classification only — no polygon extraction),
  `multivalued.py` (extended-DOF assembly). Grep for
  `integr|area|quadrature|marching|isosurf` inside `pyfp3d/wake/`: **zero
  hits**. The only "marching" in the repo is `post/section_cut.py:81-87`
  (z=const post-processing sections, unrelated).
- The existing implicit-Kutta wake-LS mechanism is **volume-weighted**
  (Π_e = ½ V_e[...] at element centroids, `multivalued.py:89-93`,
  `kernels/cut_assembly.py:185-265`) — the whole LS path integrates over
  volumes, never on the s=0 sheet. The registration's two LS options are
  accurately stated: (i) zero-isosurface polygon integration = entirely new
  geometry code (marching triangles on tets, per-sheet quadrature rules);
  (ii) a volume-band approximation = a deviation from the "sheet RHS"
  formulation (and it would need its own band-thickness convergence study to
  be interpretable at all).

### 1.3 The physics-layer gaps (path-independent — owed by ANY option)

These are the real V6 work items; they exist regardless of the GV6.0 ruling:

1. **No wake closure relations.** `viscous/closures.py` implements the wall
   profile family only (laminar Bernstein corrections + turbulent Spalding
   wall law + the Coles *wall-profile* wake component at `closures.py:216` —
   that "wake" is the turbulent wall-profile term, not a free-shear wake
   closure). Symmetric-wake shape-factor evolution, wake dissipation/lag
   relations: absent.
2. **No TE confluence logic.** The two TE copies (upper/lower) are separate
   IBL outflow points with no coupling (`coupling.py:17-21`); the V6
   thickness-continuity condition δ_wake(TE) = δ*_upper + δ*_lower exists
   only as roadmap text. (GV5.5's `te_extrapolate` is a one-sided
   first-order outflow BC, default OFF — not a confluence.)
3. **No wake IBL state.** The "reserved wake unknowns" of V1 are a
   layout-isomorphism declaration (`surface_mesh.py:28-32,261-264`: one
   SurfaceMesh instance per boundary group, same 6-equation block) — there
   is no second SurfaceMesh, no wake station structure, and the IBL3 solver
   takes exactly one SurfaceMesh (`ibl3.py:1182-1185`).
4. **No δ*_wake/u_e,wake producer anywhere** (grep `delta_wake|ds_wake`
   = 0 hits).

### 1.4 The GV6.2 target case is conforming-native

GV6.2 measures the on/off effect on **the GV3.1 case** — NACA0012 quasi-2.5-D,
M 0.5 / α 2° / Re 3.0e6, forced x_tr/c 0.05, run on
`cases/meshes/naca0012_2.5d/*.msh` **cut with `cut_wake()`** (conforming
wake-cut), against the committed XFOIL reference
(`cases/reference_data/naca0012_viscous_xfoil/`, whose DUMP tables carry
wake rows — usable for the direction check). The science question of V6 —
*how large is the wake-sheet displacement effect on cl / TE-region Cp* — is
therefore answerable on the conforming path alone. The LS leg is a
solver-completeness item, not a science item.

---

## 2. The registered options, priced

### Option A — V6 closes conforming-only; the LS leg becomes a recorded follow-up

- Conforming work (GV6.1): wake SurfaceMesh (V1 layout, per the registration)
  + wake station structure (each TE station seeds a straight-wake line per
  the design constraint: **straight wake, mass-transpiration relaxation, no
  geometric relaxation**) + wake closures + δ*_wake→ṁ→b_wake→Tᵀ wiring +
  TE-continuity assert + δ*_wake = 0 bit-identity gate + sign-pin MMS test.
  All solver-side pieces exist; the new code is physics (closures, wake
  case builder) plus wiring.
- LS work: none in V6; recorded follow-up (its own pre-registration if/when
  opened).
- Risk: the wing gates (GV5.x) ran on the conforming path, so closing V6
  conforming-only leaves the LS path without the wake correction — but the
  LS path also lacks the entire tight-coupling machinery of V5, so this is
  consistent with the project state, not a new asymmetry.

### Option B — LS leg via zero-isosurface polygon integration (new geometry code)

- Marching-triangles on cut tets + sheet quadrature + a new sheet-load
  assembler + validation gates for the geometry itself (area normals,
  orientation, span blending at the tip clip) — a geometry subsystem
  comparable in class to `cut_elements.py` itself, before any physics.
- Estimated cost: an order of magnitude more code than Option A's wiring;
  new failure modes (slivers, on-sheet degeneracies) that have nothing to do
  with the wake physics being tested.

### Option C — LS leg via a volume-band approximation

- Deviates from the registered "sheet RHS" formulation (the registration
  itself flags this). Interpretability requires a band-thickness convergence
  study; the result would be a *different* numerical object from the
  conforming sheet source, complicating any conforming-vs-LS comparison.

### Pricing summary

| Option | Science (GV6.2 answer) | Code cost | Risk |
|---|---|---|---|
| A | full | wiring + physics closures | low |
| B | full + LS completeness | new geometry subsystem | high, orthogonal to physics |
| C | partial (approximate object) | medium + convergence study | medium; deviates from formulation |

## 3. Secondary design point — how δ*_wake is produced (ruling requested too)

The registration fixes the *channel* (GV6.1: δ*_wake enters via the
`constraints/wake.py` reduce RHS (Tᵀ b_wake); TE thickness continuity
asserted; δ*_wake = 0 bit-identical) but not the *producer*. Two readings:

- **(i) Prescribed wake-δ\* (recommended for GV6.1/GV6.2)**: δ*_wake(TE) =
  δ*_upper + δ*_lower by construction (the continuity assert is then a
  construction identity, verified at runtime), relaxed downstream along the
  straight wake by a transpiration-relaxation model (no geometric
  relaxation, per the design constraint; the relaxation law pinned in the
  GV6.1 pre-registration, XFOIL-referenced). This closes GV6.1 + GV6.2 and
  answers the science question at minimal cost.
- **(ii) Solved wake IBL**: the full V6 vision ("same 6 equations with wake
  closure relations") — a second SurfaceMesh with a wake station structure
  and new free-shear wake closures (Drela 2013 / XFOIL practice; local copy
  `docs/references/Drela_2013_IBL3_general_configurations.pdf`). Substantially
  larger: new closures, a wake solver instance, TE-confluence coupling into
  the wall IBL, and its own verification ladder.

Recommendation: take **(i) for GV6.1/GV6.2** and register **(ii)** as the
follow-up "solved wake IBL" item (to be opened only if the GV6.2 measured
effect is large enough to matter against the A4 input band).

## 4. Recommendation

**Option A + producer (i)**: close V6 conforming-only with a prescribed
(TE-continuity + straight-wake relaxation) δ*_wake entering through the
existing Tᵀ RHS channel; register the LS sheet-source leg and the solved
wake-IBL producer as separate recorded follow-ups. This answers the V6
science question (the GV6.2 on/off delta on the GV3.1 case, direction-checked
against XFOIL's wake modelling, with the A4 input band quoted) with no new
geometry code and no deviation from the registered formulation.

## 5. Ruling requested (user)

1. GV6.0 main ruling: **A** (recommended) / B / C.
2. δ*_wake producer: **(i)** prescribed (recommended) / (ii) solved wake IBL.
3. If A+(i): GV6.1 opens with its own pre-registration (closures/relaxation
   law, the sign-pin MMS, the bit-identity gate, the TE-continuity assert)
   committed before the first code change, per Track-V discipline.

## 6. RULING (user, 2026-07-25)

- GV6.0 main ruling: **Option A** — V6 closes conforming-only; the LS
  sheet-source leg becomes a recorded follow-up (its own pre-registration
  if/when opened).
- δ*_wake producer: **(i) prescribed** — TE-continuity δ_wake(TE) =
  δ*_upper + δ*_lower by construction + a straight-wake relaxation law
  (pinned in the GV6.1 pre-registration, XFOIL-referenced); the solved wake
  IBL ("same 6 equations with wake closure relations") is registered as the
  follow-up producer (ii), to be opened only if the GV6.2 measured effect
  is significant against the A4 input band.
- Process: the GV6.1 pre-registration is committed on this branch before
  the first code change; branch push + PR deferred to bundle the
  adjudication document and the GV6.1 pre-registration together (user,
  2026-07-25).

---

*Survey basis: two read-only codebase surveys dated 2026-07-25 (wake /
constraints / transpiration machinery and viscous / case-layout / reference
data), cited inline above. No files were modified for this document.*

---

## Addendum 2026-07-25 — V6 close-out adjudication (user)

The GV6.2 measured effect came in NOT significant vs the A4 input band
(0 PASS / 0 FAIL / 24 RECORDED,
`bench/studies/v6_2_measured_effect/VERDICT.md`: on/off Δ-cl +0.00015
(+0.0547 %) = 0.022× the 2.5 % band; TE-region max |ΔCp| 0.00250 =
0.051× the propagated δCp_A4 0.0493; L-robust over L_rel ∈ {0.5, 1.0,
2.0} c; XFOIL wake direction agrees, rate/TE-anchor model-form
differences recorded). Per §6's producer-(ii) clause ("opened only if
the GV6.2 measured effect is significant vs the A4 input band"), the
condition is NOT met — adjudicated 2026-07-25 (user): **producer (ii)
NOT opened; V6 ✓ CLOSED**. The LS sheet-source leg remains a recorded
follow-up (its own pre-registration if/when opened).
