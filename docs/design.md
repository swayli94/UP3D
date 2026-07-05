# pyFP3D — Design Reference
## A 3D Unstructured-Grid Full-Potential Transonic External-Flow Solver in Python + Numba

Status: pre-implementation design. Intended as the repository-level reference for
incremental ("vibe-coding") development. Every phase in §11 has an explicit
acceptance test; do not advance phases until the gate passes.

---

## 1. Scope and positioning

**Target problem.** Steady transonic external flow over wings and wing-body
configurations, freestream Mach 0.3–0.87, attached or mildly shocked flow
(local normal-shock Mach ≲ 1.3). Isentropic, irrotational, inviscid.

**What FP buys over Euler.** One scalar unknown φ per node instead of five
conservative variables; elliptic-dominated operator amenable to strong implicit
solvers (Newton + AMG-preconditioned Krylov); 1–2 orders of magnitude cheaper
per solution. Ideal as the aerodynamic engine inside optimization / surrogate
loops, and as the inviscid core for later viscous–inviscid coupling (the same
role TSFOIL plays in 2D — this solver is the 3D unstructured analogue of the
pyTSFoil TSD→FP upgrade path, now free of small-disturbance and Cartesian-grid
restrictions).

**What FP cannot do.** Rotational post-shock flow (shock is replaced by an
isentropic compression satisfying a mass-conserving jump, not
Rankine–Hugoniot); shock position error grows with shock strength; known
non-uniqueness of FP solutions for airfoils in a narrow band around
M∞ ≈ 0.82–0.85 at low lift (Steinhoff–Jameson). These are accepted model
limitations, documented in §12, not implementation targets.

**Mesh class.** Unstructured tetrahedral (optionally mixed prism/tet later),
vertex-centered. The wake must be a conforming internal surface in the mesh
(see §4, §12 — this is the single strongest constraint imposed on the mesh
generator).

---

## 2. Governing equations and normalization

Velocity potential Φ with **V** = ∇Φ. Nondimensionalize lengths by reference
chord c, velocity by U∞, density by ρ∞. With φ the nondimensional potential and
q = |∇φ|:

**Conservative full-potential equation**

    ∇·(ρ ∇φ) = 0                                   (2.1)

**Isentropic density law**

    ρ(q²) = [ 1 + (γ−1)/2 · M∞² (1 − q²) ]^{1/(γ−1)}          (2.2)

**Local Mach number**

    M²(q²) = q² M∞² / [ 1 + (γ−1)/2 · M∞² (1 − q²) ]           (2.3)

**Speed of sound (nondim, squared)**

    a² = 1/M∞² + (γ−1)/2 · (1 − q²)                            (2.4)

Useful derived quantities (implement as pure numba functions, unit-tested
against hand values):

- Critical speed q\*² where M = 1:  q\*² = [2 + (γ−1)M∞²] / [(γ+1)M∞²].
- dρ/d(q²) = −(M∞²/2) ρ^{2−γ} = −ρ/(2a²) · … ; in practice implement
  ρ′(q²) = −(M∞²/2) · ρ^{(2−γ)} directly from (2.2).
- Mass flux ρq along a streamline is maximized exactly at M = 1. This
  non-monotonicity is *the* mathematical reason central discretizations fail
  in supersonic zones and upwinding of ρ is required (§3).

**Pressure coefficient (exact isentropic, not linearized):**

    Cp = 2/(γ M∞²) · [ ρ^γ − 1 ]                               (2.5)

with ρ from (2.2). Never use the linearized Cp = −2u; this solver's whole point
is to be valid where linearization fails.

**Character of (2.1).** Quasi-linear second order; in non-conservative form
(a²δij − uiuj) ∂²φ/∂xi∂xj = 0 — elliptic where M < 1, hyperbolic where M > 1.
Type-dependent discretization is mandatory; we do it through density upwinding
rather than operator switching (no Murman–Cole style switch logic on
unstructured grids).

---

## 3. Transonic treatment: artificial density (flux biasing)

Follow the Hafez–South–Murman / Holst artificial-compressibility family: in
supersonic regions replace ρ in (2.1) by an upwinded density

    ρ̃ = ρ − ν · Δℓ · (∂ρ/∂ℓ)_upwind                            (3.1)

where ℓ is the local streamline direction, Δℓ a local mesh length, and the
switching function

    ν = C · max( 0 , 1 − M_c² / M² )                            (3.2)

with cutoff Mach M_c ≈ 0.95–1.0 and constant C ∈ [1, 2] (C, M_c are solver
parameters — same calibration philosophy as the pyTSFoil parameter set; expose
them in the config from day one).

**Unstructured implementation (element-based, first order in the upwind
term):**

For each element e with velocity **V**_e = ∇φ|_e (constant on linear tets):

1. Identify the *upstream neighbor* element u(e): among face-neighbors of e,
   the one whose face normal (outward from e) is most anti-aligned with
   **V**_e. Precompute face adjacency once; recompute u(e) each nonlinear
   iteration (cheap: 4 dot products per tet).
2. Set ρ̃_e = ρ_e − ν_e (ρ_e − ρ_{u(e)}). This is (3.1) with Δℓ·∂ρ/∂ℓ ≈
   ρ_e − ρ_upstream, robust and grid-independent in form.
3. Use ρ̃_e in the element flux/residual (§6). Subsonic elements (ν = 0)
   reduce identically to the central Galerkin scheme — second order accurate.

Properties to preserve (write asserts/tests for these):
- ν = 0 ⇒ scheme is symmetric and the Picard matrix is SPD (§8).
- Monotone shock capture over 2–3 cells; no expansion shocks (upwinding of ρ
  provides the requisite dissipation only in compression through sonic).
- Freestream preservation: uniform φ = x must give machine-zero residual on
  any mesh, including across the wake cut. This is the first regression test
  after every kernel change.

---

## 4. Circulation, wake model, Kutta condition in 3D

Lifting flow ⇒ φ is multivalued; the domain must be cut by a **wake surface**
W emanating from the trailing edge and extending to the far field.

**Mesh requirement.** W is a conforming internal surface: nodes on W are
duplicated (upper copy i⁺, lower copy i⁻), elements above attach to i⁺, below
to i⁻. Gmsh: model the wake as an embedded surface and duplicate with a
plugin/crack step, or generate the duplication in the solver's mesh
preprocessor from a tagged internal surface (recommended — keep mesh
generation vanilla, do node splitting in `mesh/wake_cut.py`).

**Jump conditions on W:**

    [φ] ≡ φ⁺ − φ⁻ = Γ(s)      (constant along streamwise wake lines)   (4.1)
    [ρ ∂φ/∂n] = 0             (mass flux continuous)                    (4.2)

where s parameterizes the spanwise direction. (4.2) is enforced naturally by
assembling wake elements as interior (fluxes from both sides sum into the
constraint system); (4.1) by a master–slave linear constraint:

    φ_{i⁺} = φ_{i⁻} + Γ(s_i)                                            (4.3)

Implementation: eliminate the i⁺ DOFs (fold their rows/columns into i⁻ plus a
right-hand-side contribution proportional to Γ). Keep the elimination map
precomputed; Γ enters only the RHS, so Γ updates do not require matrix
re-assembly in Picard iterations.

**Kutta condition.** Pressure equality at the TE upper/lower ⇒ in potential
form, per spanwise station j:

    Γ_j^{new} = φ_{TE,j}^{upper} − φ_{TE,j}^{lower}                     (4.4)

evaluated one node off the TE (or extrapolated to the TE), updated in the
outer nonlinear loop with under-relaxation ω_Γ ≈ 0.7–1.0. Γ(s) between TE
stations: piecewise linear in the spanwise parameter; wake lines inherit the Γ
of the TE station they emanate from — the wake mesh preprocessor must build
this station→wake-line map. Convergence monitor: ‖ΔΓ‖∞ alongside the residual
norm.

Wake geometry: planar, aligned with freestream (or the chord plane) — standard
FP practice; force-free wake relaxation is out of scope (error is second order
in loading for attached flow).

---

## 5. Boundary conditions

| Boundary | Condition | Implementation |
|---|---|---|
| Solid wall (wing, body) | ρ ∂φ/∂n = 0 | Natural (do-nothing) in the Galerkin weak form — zero surface integral. No penalty, no ghost cells. |
| Far field | φ → φ∞ + φ_vortex | Dirichlet on outer nodes: φ∞ = x cosα cosβ + y sinβ + z sinα cosα-corrected uniform flow, plus a compressible horseshoe-vortex correction with total circulation ΣΓ_j Δs_j (Prandtl–Glauert-scaled). At R ≳ 25–50 chords the correction is small but including it lets the domain shrink to ~15 chords. |
| Wake | (4.1)–(4.2) | Constraint elimination, §4. |
| Symmetry plane (half-model) | ∂φ/∂n = 0 | Natural, free. |

One subtlety: with all-Neumann walls and Dirichlet far field the system is
well-posed; if a pure-Neumann variant is ever used (e.g. channel flows), pin
one node.

---

## 6. Spatial discretization

**Galerkin FEM on linear tetrahedra** (P1). On P1 tets this is algebraically
equivalent to a vertex-centered median-dual finite-volume scheme, but the FEM
view gives the cleanest implementation: no dual-mesh metric construction, no
face loops — only an element loop.

Weak form: find φ (satisfying Dirichlet + wake constraints) s.t. for all test
functions N_i,

    R_i = Σ_e  ρ̃_e ( ∇φ_e · ∇N_i ) V_e  = 0                            (6.1)

with ∇φ_e = Σ_k φ_k ∇N_k|_e. Precompute per element: volume V_e and the 4×3
shape-gradient matrix B_e (constant). Then per nonlinear evaluation, per
element: one 4×3 gemv for ∇φ_e, density law, upwind lookup, one 3×4 gemv for
the scatter — trivially numba-friendly.

**Picard (frozen-density) matrix:**

    A_ij = Σ_e ρ̃_e (∇N_i · ∇N_j) V_e                                    (6.2)

— a weighted stiffness matrix. SPD when ν ≡ 0. With upwinding, ρ̃_e couples to
the upstream element's DOFs; in **Picard** we deliberately ignore that coupling
in A (treat ρ̃_e as a frozen scalar) — A stays symmetric, and the upwind
physics enters through the residual/RHS lag. This is the classical, very robust
scheme; its convergence rate degrades with shock strength, motivating the
Newton option in §8.

**Newton Jacobian (Phase 5+):** differentiate (6.1) w.r.t. φ_k:

    ∂R_i/∂φ_k = Σ_e V_e [ ρ̃_e ∇N_i·∇N_k
                + (∂ρ̃_e/∂q²_e) · 2(∇φ_e·∇N_k)(∇φ_e·∇N_i) ]              (6.3)

plus the upwind chain ∂ρ̃_e/∂ρ_{u(e)} coupling e's test functions to u(e)'s
DOFs (widens the stencil by one layer of elements; the sparsity map must be
built from the element + upwind-neighbor graph). The exact Jacobian is
nonsymmetric and indefinite in supersonic zones — GMRES + ILU, not CG.

Order of accuracy: 2nd in subsonic regions, 1st locally at captured shocks —
standard and acceptable.

---

## 7. Numba kernel architecture

Design rules (hard-won numba constraints — encode these in the project skill):

1. **Struct-of-arrays only.** All mesh/solution state as flat, contiguous,
   explicitly-typed `np.float64` / `np.int32` arrays. No Python objects, no
   dicts, no dataclass instances inside `@njit`. A thin Python-level `Mesh`
   dataclass *holds* the arrays; kernels receive arrays as arguments.
2. **Two-pass or colored assembly for `prange`.** Scatter-add races are the
   central parallelization hazard:
   - Residual: greedy **element coloring** (no two same-color elements share a
     node); outer serial loop over colors, inner `prange` over that color's
     elements. Coloring computed once at startup (~8–12 colors typical for
     tets).
   - Matrix: precompute `elem_to_csr[e, 4, 4] → nnz index` map once from the
     symbolic sparsity pattern; assembly writes are then to disjoint or
     colored locations identically.
3. **`@njit(cache=True, fastmath=True)`** on all kernels; `parallel=True` only
   where profiled to matter (element loops; not the small BC loops).
4. **No allocation inside hot kernels.** Preallocate residual, ∇φ workspace
   (per-color or per-thread), density arrays; pass in.
5. **Pure functions for physics.** `density(q2, minf2, gamma)`,
   `mach2(q2, ...)`, `cp(...)` as scalar njit'd functions — unit-testable in
   isolation, reused by post-processing.
6. Linear algebra stays in **SciPy/PyAMG land** (compiled, no need to rewrite):
   `scipy.sparse.csr_matrix` + `pyamg.smoothed_aggregation_solver` as a
   preconditioner for `scipy.sparse.linalg.cg` (Picard/subsonic) or `gmres`
   (Newton/transonic). Optional matrix-free numba matvec later; not Phase 1.

**Module layout**

```
pyfp3d/
├── config.py            # dataclass: M∞, α, γ, ν-params (C, M_c), tolerances, relaxations
├── mesh/
│   ├── reader.py        # meshio front-end: gmsh/.su2/.cgns → raw arrays
│   ├── wake_cut.py      # node duplication along tagged wake surface; TE-station→wake-line map
│   ├── metrics.py       # numba: B_e, V_e, face adjacency, quality checks
│   └── coloring.py      # greedy element coloring for parallel assembly
├── physics/
│   └── isentropic.py    # numba scalars: ρ(q²), M²(q²), a², Cp, q*², ρ′
├── kernels/
│   ├── gradient.py      # ∇φ_e per element
│   ├── upwind.py        # upstream-element id, ν_e, ρ̃_e
│   ├── residual.py      # colored assembly of R (6.1)
│   └── jacobian.py      # Picard matrix (6.2); later Newton (6.3)
├── constraints/
│   ├── dirichlet.py     # far-field values incl. vortex correction
│   └── wake.py          # master–slave elimination, Γ handling, Kutta update (4.4)
├── solve/
│   ├── linear.py        # CG/GMRES + AMG/ILU wrappers, tolerances
│   ├── picard.py        # outer loop: ρ̃ update → Γ update → linear solve → relax
│   ├── newton.py        # Phase 5: exact Jacobian, line search / pseudo-transient
│   └── continuation.py  # Mach/α ramping for hard cases
├── post/
│   ├── surface.py       # nodal Cp, sectional cl, forces & moments
│   ├── trefftz.py       # induced drag from Γ(s) (far-field), cross-check vs near-field
│   └── vtk_out.py       # meshio VTK export: φ, V, M, ρ, Cp
└── cases/               # regression decks: sphere, naca0012_extruded, oneraM6
```

---

## 8. Nonlinear solution strategy

**Baseline: relaxed Picard (Phases 1–4).**

```
assemble sparsity, coloring, wake elimination map      # once
φ ← farfield initial guess (uniform flow)
repeat:
    ∇φ_e, q²_e, ρ_e            # element sweep
    u(e), ν_e, ρ̃_e             # upwind sweep
    A(ρ̃), b(Dirichlet, Γ)      # colored assembly (or RHS-only if ρ̃ lagged further)
    solve A δφ = −R             # CG+AMG (sym) / GMRES+ILU (if unsym pieces added)
    φ ← φ + ω δφ                # ω ≈ 1 subsonic; 0.7–0.9 transonic
    Γ_j ← Γ_j + ω_Γ (Kutta target − Γ_j)
until ‖R‖₂/‖R₀‖₂ < 1e−10 and ‖ΔΓ‖∞ < tol
```

Expected behavior: subcritical cases converge in O(10–30) Picard iterations;
transonic O(100–300). Each linear solve with AMG is O(N) — total cost minutes
for ~1–3 M nodes on a workstation, which is the design target.

**Accelerations, in order of implementation value:**
1. AMG preconditioner reuse across Picard iterations (re-setup every k≈5–10).
2. Inexact solves: loose linear tolerance (1e−2 relative) early, tighten near
   convergence (Eisenstat–Walker style schedule).
3. Mach continuation: converge M∞ − 0.05 first, restart. Essential near the
   FP non-uniqueness band and for shocked ONERA M6.
4. Full Newton with (6.3) + pseudo-transient continuation (add V_e/Δτ mass
   term) — quadratic terminal convergence; only after Picard is bulletproof.

---

## 9. Post-processing and forces

- Nodal q², M, Cp via volume-weighted element averages (or superconvergent
  patch recovery later).
- Surface forces: integrate Cp·n over wall triangles → C_L, C_D_pressure,
  C_m. **Warning:** near-field pressure drag in FP is contaminated by the
  isentropic shock model and discretization; report it, but treat Trefftz-plane
  induced drag from Γ(s) (+ a wave-drag estimate from shock-swept entropy proxy
  if ever needed) as the trustworthy decomposition.
- Sectional cl(η) from spanwise Γ: cl = 2Γ/(c(η)) — direct cross-check against
  surface integration; the two must agree to ~1% on a converged solution
  (excellent bug detector for wake/Kutta code).
- VTK output every k iterations for shock-position monitoring.

---

## 10. Verification & validation ladder

Each item is a pytest in `cases/`; keep them runnable on coarse meshes in CI.

| # | Case | Checks | Gate |
|---|---|---|---|
| V0 | Freestream preservation, arbitrary tet mesh + wake cut | ‖R(φ=x)‖∞ | < 1e−12 |
| V1 | Laplace MMS (ν=0, ρ=1), manufactured φ = sin·cos products | L2 error vs h | slope ≥ 1.9 |
| V2 | Incompressible sphere | Cp vs 1 − (9/4)sin²θ | max err < 2% (medium mesh) |
| V3 | Subcritical extruded NACA0012, M 0.5, α 2° (span-periodic or large AR + symmetry) | Cp vs 2D FP/panel reference; cl | Δcl < 2% |
| V4 | Transonic NACA0012 M 0.80 α 1.25° (extruded) | shock x/c vs published FP (Holst AGARD FP workshop results — note FP ≠ Euler shock position) | shock Δx/c < 0.03 |
| V5 | ONERA M6, M 0.84, α 3.06° | λ-shock topology, section Cp at η = 0.44/0.65/0.90 vs FP references (FLO-class / TRANAIR); CL | qualitative λ + CL within FP-literature scatter |
| V6 | Γ-consistency | sectional cl from Γ vs surface integration | < 1% |

Reference data caveat: validate against *full-potential* published results
first (Holst's Progress in Aerospace Sciences 2000 review collects them), Euler
second — otherwise model error and code error are confounded.

---

## 11. Development roadmap (vibe-coding phases)

> **Superseded numbering.** The phase order below (P0–P6) is the original design-time
> sketch. [docs/roadmap.md](roadmap.md) Track P is the *active* phase order and takes
> precedence where the two disagree — notably, roadmap.md's P2 is wake cut/circulation/
> Kutta (critical phase) and P3 is subsonic compressible, the reverse of the ordering
> below. Treat this section as theory/gate reference only (§10 gates, phase content
> descriptions); for "what phase are we in" and gate sequencing, follow roadmap.md.

Each phase is a self-contained PR-sized unit with its gate from §10.

- **P0 — Mesh infrastructure.** meshio reader, metrics (B_e, V_e), adjacency,
  coloring, VTK writer. Gate: metrics unit tests (ΣV_e = volume of unit cube;
  ∇(linear field) exact), coloring validity assert.
- **P1 — Laplace solver.** ρ ≡ 1, Dirichlet far field, natural walls, CG+AMG.
  Gates: V0 (without wake), V1, V2.
- **P2 — Subsonic compressible, non-lifting.** Density law + Picard loop, no
  upwinding, no wake. Gate: sphere at M∞ = 0.3 vs Prandtl–Glauert-corrected V2;
  convergence in <30 Picard its.
- **P3 — Lift: wake cut + Kutta.** Node duplication, constraint elimination,
  Γ update, vortex far-field correction. Gates: V0 *with* wake cut, V3, V6.
- **P4 — Transonic: artificial density.** Upwind element search, ν switch,
  ρ̃; relaxation + Mach continuation. Gates: V4, then V5.
- **P5 — Performance & robustness.** Newton (6.3), pseudo-transient, AMG
  reuse, profiling (target: ONERA M6 medium mesh < 5 min single node).
- **P6 — Extensions (backlog).** Mixed prism/tet; embedded-boundary wake
  alternative; VII coupling hook (transpiration BC ∂φ/∂n = d(u_e δ*)/ds —
  reuses the IBL work from pyTSFoil); adjoint via the Newton Jacobian
  transpose (nearly free once (6.3) exists — high value for the MDO thread).

---

## 12. Known risks and mitigations

1. **Wake-conforming mesh is the #1 practical risk.** Everything else is
   textbook; getting Gmsh to embed a wake surface and the preprocessor to
   split nodes robustly (TE seam, tip edge, wake–far-field intersection) will
   consume real effort. Mitigate: build `wake_cut.py` against a trivially
   simple extruded-airfoil mesh first; write topological asserts (each wake
   face has exactly one ⁺ and one ⁻ element; TE nodes not duplicated).
2. **FP non-uniqueness** (M∞ ≈ 0.82–0.85, low α, conventional airfoils):
   continuation direction can select different branches. Document; always ramp
   Mach upward from a subcritical converged state.
3. **Strong shocks** (M_local > 1.3): isentropic jump under-predicts shock
   pressure rise, shock sits aft of Euler. Not a bug; state validity envelope
   in the README and warn at runtime when max(M) exceeds 1.35.
4. **Picard stall on shocked cases**: symptoms — residual plateau with shock
   oscillating one cell. Mitigations in order: raise C in (3.2), lower ω,
   continuation, then Newton+PTC.
5. **Numba compile-time friction** during vibe coding: keep kernels small and
   argument lists explicit; `cache=True`; a `PYFP3D_NOJIT=1` env switch that
   swaps `njit` for identity decorator makes pdb-debugging possible.

---

## 13. Core references

- Holst, T.L., "Transonic flow computations using nonlinear potential
  methods," *Progress in Aerospace Sciences* 36 (2000) — the canonical review;
  contains the artificial-density family, workshop validation data.
- Hafez, M., South, J., Murman, E., "Artificial compressibility methods for
  numerical solutions of transonic full potential equation," *AIAA J.* 17(8),
  1979.
- Jameson, A., "Iterative solution of transonic flows over airfoils and
  wings," *CPAM* 27 (1974); FLO-series reports for wing FP baselines.
- Steinhoff, J., Jameson, A., "Multiple solutions of the transonic potential
  flow equation," *AIAA J.* 20(11), 1982 — non-uniqueness.
- Neel, R.E., *Advances in the Computation of Transonic Full Potential Flows
  on Unstructured Grids*, PhD thesis, Virginia Tech, 1997 — closest prior art
  to this exact solver concept (vertex-based unstructured FP).
- Caughey, D.A., Jameson, A., "Basic advances in the finite-volume method for
  transonic potential-flow calculations" — FV/FEM equivalence viewpoint.
- Drela, M., *Flight Vehicle Aerodynamics*, MIT Press 2014 — wake/Kutta and
  Trefftz-plane treatment in potential methods.
