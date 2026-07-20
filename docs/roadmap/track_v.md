# pyFP3D Roadmap — Track V (viscous–inviscid interaction V1–V4)

> Split verbatim from `docs/roadmap.md` on 2026-07-15 (content unchanged; only
> this header and the ledger heading were added). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md).

## Track V — Viscous–inviscid interaction (IBL coupling; designed 2026-07-09/10; NOT started)

Deliverable: `viscous/` — Drela IBL3 6-equation integral boundary layer
(δ, A, B, Ψ, C_τ1, C_τ2; surface Galerkin FE on wall + wake sheet — no streamline
integration) coupled to the FP solver through a **transpiration BC** (no mesh
motion), progressing loose → tight coupling. Design records: DN2 + DN6
(historical; `discussion_notes/` deleted 2026-07-14 — recover via
`git show 8aa4aee:docs/discussion_notes/20260707_2118_ibl_viscous_coupling_design.md`
and `.../20260709_0145_3d_vii_implementation_analysis.md`; the load-bearing
content is summarized here and in design.md §5/§8 notes; design.md §11 is a
pointer to the roadmap/ files since 2026-07-15). Design
complete; **no implementation exists yet**.

> **Naming disambiguation.** Track V phase IDs V1–V4 are development *phases*;
> they are unrelated to the §10 validation-ladder case IDs V0–V6 in design.md
> (e.g. "V6 < 1%" remains the Γ-consistency metric). Gates here are GV<phase>.<n>.
> (The historical DN6 note internally also used "V2 = 3D crossflow" in one
> diagram — the phase entries below are the binding ones.)

### V1 — IBL3 solver + loose coupling ☐
**Deliverable:** IBL3 solver + loose coupling: `viscous/ibl3.py`, `viscous/transpiration.py` (δ*→ṁ), `viscous/coupling.py`; wall-Neumann blowing source in `kernels/residual.py`; wake-sheet mass source in `constraints/wake.py`
**Gates (sketch):**
- [ ] GV1.1 2D NACA0012 subsonic δ* vs XFOIL
- [ ] GV1.2 loose loop converges in 5–10 iterations
- [ ] GV1.3 δ* = 0 bit-identical to the inviscid path
**Prereq:** P6 ✓ (smoothed wall gradient is the IBL edge-velocity input)

### V2 — Quasi-simultaneous coupling ☐ (optional)
**Deliverable:** Quasi-simultaneous coupling (Hilbert-integral potential surrogate, `viscous/hilbert.py`) — **optional**, skip if V1 converges fast enough
**Gates (sketch):**
- [ ] GV2.1 30–50% fewer coupling iterations
- [ ] improved robustness near separation
**Prereq:** V1

### V3 — Tight coupling: augmented Newton ☐
**Deliverable:** Tight coupling: augmented (φ, Γ, δ, A, B, Ψ, C_τ1, C_τ2) Newton; J_φ,δ* = ∂ṁ/∂δ* and J_δ*,φ blocks; GMRES + block preconditioning (AMG potential / ILU BL)
**Gates (sketch):**
- [ ] GV3.1 quadratic convergence
- [ ] GV3.2 2D transonic VII shock shift vs experiment
- [ ] GV3.3 M6: CL moves **down** from the converged inviscid value toward experiment ≈ 0.26–0.27 (direction + magnitude, DN6 §13.3)
**Prereq:** P8 + V1

### V4 — Wake-sheet IBL correction ☐ (a continuation of V1, not an independent phase)
**Deliverable:** Wake-sheet IBL correction — a continuation of V1, **not** an independent phase (same 6 equations, wake closure relations; reserve the wake unknowns in V1's data layout)
**Gate sketch (design constraints):**
- [ ] δ*_wake enters as the wake-sheet RHS mass source
- [ ] TE kink absorbed by Drela local-basis adaptation (DN2 §4.5.1–4.5.3)
- [ ] straight wake + mass-transpiration relaxation, no geometric relaxation (DN2 §4.5.6)
**Prereq:** V1

Scope guards (DN2 §9, DN6 §13–14):

- Validity envelope: attached / mildly-shocked flow; not for massive separation
  or shock-induced separation.
- **VII does not close the inviscid-discretization CL gap.** The M6 CL 0.245 vs
  the FP-reference 0.288 is attributed to the sharp-TE/LE P1 wall-gradient floor
  → **P11** (curved elements), with **P9** first discriminating how much of the
  gap is plain resolution (post-2026-07-11 numbering; this bullet said "P9 =
  curved elements" under the old IDs). Viscosity moves CL *down* from the
  (accurate) inviscid value toward experiment — do not book the 0.245→0.288
  recovery to Track V.
- V1 is parallelizable with P7/P8 (depends only on P6), but it is a large,
  self-contained solver effort (6-equation nonlinear surface FE + closures);
  budget it like a Track-P phase, not a side task.

---


## Progress ledger

### Track V — viscous–inviscid interaction

Track status: **☐ NOT STARTED** — design complete 2026-07-09/10 (DN2 + DN6,
historical notes, see the Track V header above); no implementation exists yet. Validity envelope: attached /
mildly-shocked flow, not massive or shock-induced separation. **VII does not
close the inviscid-discretization CL gap** — viscosity moves CL *down* toward
experiment (≈ 0.26–0.27 on M6), so the 0.245-vs-0.288 gap belongs to P9
(discriminating how much of it is resolution) and P11 (curved walls), never to
Track V.

- V1 — ☐ — IBL3 solver + loose coupling (`viscous/ibl3.py`, `transpiration.py`, `coupling.py`;
  wall-Neumann blowing source in `kernels/residual.py`; wake-sheet mass source in `constraints/wake.py`).
  Depends only on **P6 ✓** (the smoothed wall gradient is the IBL edge-velocity input) and is parallelizable with P7/P8 —
  but it is a Track-P-sized effort (6-equation nonlinear surface FE + closures), not a side task. Gates GV1.1–GV1.3.
- V2 — ☐ (optional) — Quasi-simultaneous coupling (Hilbert-integral potential surrogate, `viscous/hilbert.py`).
  Skip if V1's loose loop converges fast enough. Gate GV2.1.
- V3 — ☐ — Tight coupling: augmented (φ, Γ, δ, A, B, Ψ, C_τ1, C_τ2) Newton on top of the **P8** Jacobian machinery;
  GMRES + block preconditioning (AMG potential / ILU BL). Gates GV3.1–GV3.3 — including the direction check:
  M6 CL moves **down** from the converged inviscid value toward experiment ≈ 0.26–0.27.
- V4 — ☐ — Wake-sheet IBL correction — a continuation of V1, NOT an independent phase (same 6 equations + wake closure
  relations; reserve the wake unknowns in V1's data layout). Straight wake + mass-transpiration relaxation, no geometric
  relaxation.
