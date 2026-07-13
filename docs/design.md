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

### 3.1 The shipped P4 operator, and its two distinct defects (P6 design pass)

The P4 implementation (`kernels/upwind.py`) is the *element-upwind* realization
of López's FE artificial compressibility (dissertation Eq. 3.19–3.27, §13):

    ρ̃_e = ρ_e − μ_s (ρ_e − ρ_{u(e)}),   μ_s = C · max(0, μ_e, μ_{u(e)})     (3.3)
    μ = 1 − M_c²/M²          (López 3.21; UP3D guards the denominator as
                              max(M², M_c²) so μ ≡ 0 subcritically, bit-for-bit)

with `C = upwind_c` (López's μ_c), `M_c = m_crit`. The upstream element `u(e)`
is chosen by a **discrete multi-hop directional walk** through face neighbours
(`upstream_elements`, up to 4 hops, most-anti-streamwise neighbour each hop,
target reach 0.8·element streamwise extent). ρ̃_e is a **frozen per-element
scalar** that weights both the residual (6.1) and the Picard matrix (6.2). Two
hard clamps guard positivity: a q² speed limiter at M_cap = 3 (`limit_q2_field`,
López's M²_max ~ O(3), §3.4) and a ρ̃ floor at 0.05 (`rho_tilde_sweep`).

This operator sits at the centre of **two separate issues** the roadmap
originally ran together — and an N1 investigation (2026-07-08) established, with
measurement, that they have **different root causes and different fixes**:

- **The surface-Cp sawtooth (accuracy — gate G6.1) is a wall-Cp gradient
  RECOVERY artifact, NOT the flux.** The ≈2-cell serration in the supersonic-run
  wall Cp is produced by the *per-triangle constant P1 gradient recovery*
  (`triangle_tangential_gradients`, no averaging) on the sliver prism-split wall
  triangulation: adjacent sliver triangles give oscillating gradients even for a
  *smooth* φ. Decisive evidence (both on the coarse NACA0012 M0.80 solution):
  (i) recomputing the wall Cp of the **same** walk solution with a normal-gated
  edge-neighbour-smoothed gradient drops the G6.1 sawtooth metric **~330×**
  (0.0758 → 0.00023, 39 → 1 slope reversals) with the shock at 0.604→0.607 and
  cl_p within 0.3 % — the φ field was already smooth; (ii) a **smoother
  artificial-density flux** (the streamline kernel of §3.2) does **not** reduce
  the metric at all (equal-or-worse at every calibration). So G6.1 is fixed in
  post-processing (§9 `smooth_passes`), not by changing the flux. *This
  overturns the earlier attribution (a ρ̃ selection-flip) recorded in
  demo_report §P4 — that attribution is refuted.*
- **The flux's non-differentiability in φ (blocks the exact Newton Jacobian) is
  a real, separate defect — but it belongs to P7, not G6.1.** The integer `u(e)`,
  the `max(0, μ_e, μ_{u(e)})` switch, and the two hard clamps are non-smooth in
  φ; the exact Newton Jacobian (§6.3) needs ∂ρ̃/∂φ well-defined **at frozen
  selection** — which the shipped **walk** already provides (next paragraph).
  Deriving + FD-verifying that ∂ρ̃/∂φ, not a sawtooth fix, is the **P7**
  deliverable and the **P8 Newton prerequisite**; the §3.2 kernel is an *optional*
  smoother alternative, not the requirement.

The PDF-verified insight (López Appendix B, §6.3) sharpens what the P7 flux must
provide: **it does not require a differentiable *selection*.** López freezes
`u(e)` within each Newton step and differentiates only *through* ρ_{u(e)}, μ, ρ
(Eqs. B.3–B.8); his own `μ_c·max(0,μ,μ_up)` is only C⁰ yet converges
quadratically (measure-zero non-smooth locus). And the two clamps are inactive
at a converged physical solution (P5 medium: **0 floored / 0 limited**), so their
non-smoothness never enters the converged Jacobian — smoothing them (§3.2) is a
transient-robustness option, not a correctness requirement.

### 3.2 A differentiable upwind density (P7 — the P8 Newton prerequisite)

> Scope note (re-scoped 2026-07-08): this kernel is an **optional** smoother /
> ~10×-faster Picard flux, **not** the mandated P7 Newton prerequisite. The P7
> prerequisite is the *frozen-selection* differentiability of the shipped **walk**
> (§3.1, §6.3) — its ∂ρ̃/∂φ at fixed u(e) is well-defined as-is and sparse
> (~+1 element/row). This kernel is a P8 Newton-flux **candidate** only if its
> denser Jacobian (§6.3) measures net-favourable. It is not the G6.1 sawtooth fix
> (§9) and stays opt-in (`UpwindOperator(weighted=True, mode="kernel")`); shipped
> default = the P4 walk.

The naive form — replacing the single integer neighbour with a **smooth
face-neighbour weighting**

    ρ_up(e) = Σ_f w_f(e) · ρ_{nb(f)} / Σ_f w_f(e),
    w_f(e) = max(0, −V_e · n̂_f)^p         (n̂_f = outward unit face normal)     (3.4)

is C¹ in V_e = ∇φ but **transiently unstable** on the sliver tets: an N1 probe
(a fast frozen-Γ Picard from the walk's converged coarse field) blows it up
(M_max 20–40 in ~60 iterations) from a *good* field for every p and reach-gain.
The cause is **not** dissipation magnitude (at the converged field it gives 1.3×
the walk's) — the *near* ring reaches only ≈0.3·extent, and averaging over it
samples a too-local, too-central region that destroys the upwind character
stabilising the hyperbolic zone. (Two dead ends ruled out en route:
face-normal inflow weighting `max(0,−V·n̂_f)` produces *negative*-reach
anti-dissipative cells on slivers — use centroid-displacement alignment
`max(0,−V̂·d̂_{e,nb})`; and a reach multiplier `G·(ρ_e−ρ_up)` cannot rescue it,
since `G·0 = 0` when the near blend can't reach upstream.)

The **shipped opt-in operator** is therefore a **multi-ring streamline-Gaussian
kernel** with genuine reach: over a precomputed depth-3 BFS neighbourhood,
weight each cell c by a Gaussian centred one streamwise extent upstream,

    up_c = −V̂·(c_c − c_e),  perp²_c = |c_c − c_e|² − up_c²
    w_c  = exp(−½[(up_c − reach·ext)²/σ_s² + perp²_c/σ_p²]),   up_c > 0
    ρ_up = Σ_c w_c ρ[c] / Σ_c w_c                                       (3.4′)

This keeps the walk's stabilising **reach** (samples genuinely-upstream cells,
not the near ring), is C^∞ in V̂ = ∇φ/|∇φ| at fixed neighbourhood (Newton-ready),
and — no per-iteration walk to recompute — runs ~10× faster per density
iteration than the walk. Probed stable: reach ≈ 0.8–1.2·extent, σ_s,σ_p ≈
0.25–0.4·extent converge to the walk's *own* physical state (M_max 1.37,
0 floored/limited) where the near blend diverged; reach ≈ 0.6 is too short.
Preserved: the exact ν ≡ 0 subcritical no-op (G4.2), the positivity guards, and
the (M²−1)/M² dissipation floor. The shock-point operator stays a **hard**
max(ν_e, ν_up) on the smoothly-blended upstream Mach — exactly 0 subcritically,
so no smooth `max_ε` is used (max_ε(0,0)=ε≠0 would break the no-op).

**It does not reduce the surface-Cp sawtooth** (that is a recovery artifact,
§3.1/§9): calibrated to reproduce the walk's coarse shock/cl (reach 1.0, C≈2.5:
shock 0.619 vs 0.604, cl_KJ 0.356 vs 0.364) it gives an *equal-or-worse* G6.1
metric. Its value is as an **optional** Picard-speed / smoother flux; the default
P7/P8 Newton flux is the sparse frozen **walk** (§6.3). If ever promoted to the
Newton flux it needs `upwind_c` recalibration (≈2.0–2.5 vs the walk's 1.5)
because the kernel's effective dissipation differs from the walk's. Shipped
default = the P4 walk.

**Optional — smooth density clamp (only if Newton stalls on clamp-touching
transients).** Replace the hard M_cap clamp by
σ = ½(1 + tanh((M² − M²_max)/ε)), ρ_clamped = (1−σ)ρ(M²) + σ ρ(M²_max), which is
C¹ across the cap. López used a **hard** clamp and still converged quadratically
(the clamp is active only at isolated wingtip-TE geometric singularities —
dissertation Fig. 3.5, §3.4 "only local parts of the boundary are affected"), and
UP3D's converged states do not touch it, so this is scoped as a fallback, not a
prerequisite.

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

**The duplication must include the trailing-edge nodes.** In the continuum
[φ](TE) = Γ — that is the Kutta condition itself — so the cut (and the
node doubling) extends onto the wall. Keeping the TE node single-valued
tapers the jump to zero across the first wake cell, which is equivalent to
a point vortex of strength Γ parked at the TE: wall suction ~ (Γ/2πr)²
integrates to a spurious force ~ Γ²/h that *diverges* under refinement
(measured in P2: −0.27 out of an expected cl ≈ 0.6 on the coarse NACA0012
mesh, from six TE-adjacent triangles; with the TE doubled the same case
gives cl = 0.6012 vs Kutta–Joukowski 0.6).

**Kutta condition.** Pressure equality at the TE upper/lower ⇒ in potential
form, per spanwise station j:

    Γ_j^{new} = φ_{TE,j}^{upper} − φ_{TE,j}^{lower}                     (4.4)

evaluated one node off the TE (or extrapolated to the TE), updated in the
outer nonlinear loop with under-relaxation ω_Γ ≈ 0.7–1.0. On the linear
(Laplace) driver the map Γ → Γ_target is affine with measured slope
b ≈ 0.93 (the smooth-flow jump one node off the TE is nearly Γ itself), so
plain relaxation converges slowly; the P2 driver applies a per-station
secant (Aitken) step from the second update on and converges in 2–3
updates. The under-relaxed form remains the fallback and the baseline for
the nonlinear (P3+) outer loop. Γ(s) between TE
stations: piecewise linear in the spanwise parameter; wake lines inherit the Γ
of the TE station they emanate from — the wake mesh preprocessor must build
this station→wake-line map. Convergence monitor: ‖ΔΓ‖∞ alongside the residual
norm.

Wake geometry: planar, aligned with freestream (or the chord plane) — standard
FP practice; force-free wake relaxation is out of scope (error is second order
in loading for attached flow).

**Designed alternative — Track B level-set wake (2026-07-07; not started).** A
level-set + multivalued-element (CutFEM-style) wake representation with an
implicitly enforced Kutta condition (Γ emerges from the solution — no Γ DOF,
no secant) would replace the conforming sheet + master–slave elimination above,
so meshes need no pre-embedded wake surface (α sweeps without remeshing;
multi-wake for wing–body). Numerics reference:
[design_track_b.md](design_track_b.md) (2026-07-11 López-dissertation
cross-check; supersedes DN1
`discussion_notes/20260707_1505_levelset_wake_design.md` as the Track B design
spec — notably the original penalty Kutta is demoted to an optional diagnostic
in favor of TE duplication + the wake least-squares condition); phase plan
B1–B5 in roadmap.md Track B. Curved/free wake (B6) is
**shelved** (2026-07-10, DN2 §4.5.6): the straight-wake loading error is
O(θ²) ≈ 0.1%, geometry updates force per-step cut-set/DOF rebuilds that conflict
with Newton, and under VII the wake instead carries a mass-transpiration
relaxation (§5 note). The conforming path in this section stays the shipped
default throughout P7–P10.

### 4.1 Tip / wake-edge singularity (P13)

The rigid **planar** wake of the previous section terminates at a free edge —
the swept tip TE corner in 3D, from which the sheet extends downstream but has
no continuation *beyond* the tip (M1: tip nodes stay single-valued, Γ(tip) = 0
discretely). That free edge is **singular**, and the singularity is a property
of the wake **model**, not of any particular discretization or wall geometry.

**Bound Γ vs trailing vorticity.** The regularity condition at the tip is
Γ(tip) = 0 — the upper/lower speeds equalize, so there is no *bound* vortex
there. This is **necessary but not sufficient**. The physical object shed into
the wake is the **trailing vortex sheet** of strength γ(s) = −dΓ/ds (s the
spanwise parameter): a wing unloads towards the tip, and it is the *unloading
rate* |dΓ/ds| — not Γ itself — that peaks at the tip (for the near-elliptic
loading of a finite wing, dΓ/ds → −∞ at the tip in the continuum; on the B7
level-set solution's smooth Γ(z) it is measured ~10× larger at the tip than
mid-span while Γ itself falls 0.10 → 0). A **flat** sheet of trailing vorticity
that simply ends induces a crossflow that turns the free edge — the same
mechanism as the leading edge of a flat plate at incidence — so the edge
velocity behaves as **1/√r** (r the distance from the edge), a **flat-plate-edge
(square-root) singularity**. It is *not* the 1/r of a concentrated line vortex
(which would need the sheet's vorticity collected onto a line — precisely the
roll-up the rigid planar model omits). A refinement study distinguishes them by
the exponent p in peak-velocity ~ h^(−p): 1/√r ⇒ p ≈ 0.5, 1/r ⇒ p = 1.

**Evidence (P13/G13.1, 2026-07-13).** Probed subsonically (M∞ 0.5, no shock /
artificial density / limiter) on the ONERA M6 tip box, the conforming
three-point (coarse/medium/fine) log-log exponent of peak local Mach vs 1/h is
**p ≈ 0.59** — consistent with 1/√r, refuting the earlier "1/r-type" phrasing.
(**Read §4.2 before using this continuum picture to reason about the fix.** The
exponent is right, but the *mechanism the solver actually sees* is discrete —
the outermost TE station shedding its retained Γ over one cell — and it is that
discrete picture, not this continuum one, that predicts which tapers work.)
The peak diverges on **both** the conforming and the Track-B level-set paths
(the representation change does not blunt it), while the surrounding field
stays mesh-converged (a *localized* edge singularity). Under enough refinement
the edge overspeed trips the speed limiter / density floor even subsonically, so
the fine mesh stops being a discrete solution (the M∞ 0.5 analogue of the P9/G9.1
transonic finding). Demo: `cases/demo/p13_tip_edge_singularity/`.

### 4.2 The fix: a spanwise loading taper (P13/G13.2)

Curved *wall* elements (P11) cannot remove a wake-sheet edge, and Track B's
level-set changes the wake *representation*, not the model. But the cure turned
out to be **neither roll-up nor a vortex core** — both of which earlier drafts of
this section proposed. The measured mechanism is **discrete**:

**The solver never sees the continuum edge.** It sees the **outermost TE
station**, which retains circulation Γ_last and sheds it as a *concentrated
vortex over the last cell* — the sheet's free-edge nodes are single-valued, so
the jump falls from Γ_last to 0 across one element. That shed vortex induces a
velocity ~ Γ_last / h. Hence, with Γ_last ~ h^q,

    edge peak ~ h^(q−1)        ⇒        p ≈ 1 − q                       (4.5)

and the regularization criterion is **q ≥ 1**. This predicts the untapered
baseline exactly: Γ ~ √u with u_last ~ h gives Γ_last ~ √h, i.e. q = 0.44
(measured), so p ≈ 0.56 — against a measured p = 0.52. It is *not* the continuum
edge exponent that matters, and *not* a viscous core: it is how fast the last
station's circulation vanishes under refinement.

**The fix follows immediately.** Taper the accepted circulation toward the tip:

    Γ_eff(z) = F(z) · Γ_Kutta(z),   F = smoothstep((z_tip − z)/r_c)      (4.6)

applied to the per-station Kutta target (`constraints/wake.py::tip_taper_factors`;
`solve_newton_lifting(tip_taper=…)`, default `None` = bit-identical). F is
geometry-only (independent of φ), so the Newton Jacobian keeps its structure —
the Kutta row is simply scaled by the same F_j. Shipped model: `vanish_smooth`
with **r_c = 0.05 · b_semi**.

**F must have COMPACT support.** A tanh taper (F(tip)=½) also reaches q≈1 and does
regularize, but tanh never reaches 1: its tails unload the wing to η≈0.77,
costing −7.4% lift, −4.9% inboard circulation, and breaking TE pressure closure
at mid-outboard span where there is no singularity to fix. A tip model must be
*local*.

**The taper is amplified, not applied.** Γ is a fixed point of Γ = F·Γ_Kutta(Γ),
and the Kutta map has slope b ≈ 0.93 (§4, P2 measurement), so

    Γ/Γ*  =  F (1 − b) / (1 − F b)                                       (4.7)

— F = 0.8 yields 0.21×, not 0.8×. This is why the measured q (≈3.3) far exceeds
the naive exponent of F, and why r_c must be kept small.

**Result.** The edge peak goes flat (p: 0.592 → 0.009 over coarse/medium/fine),
the M6 fine mesh becomes a **genuine discrete solution** (0 limited/floored) —
what P9/G9.1 could not achieve — and at M∞0.84 the medium M_max drops 2.13 →
1.725. **Price:** a deliberate tip unloading worth ~1–1.6% of cl (scales with
r_c), localized inboard of η≈0.95. Evidence:
`cases/demo/p13_tip_edge_singularity/`.

**Caveat (P13/G13.3): a THIRD singularity, on the WALL.** Removing the tip
singularity did **not** deliver 3D grid convergence — the M6 lift sequence is
still not in the asymptotic range (its increments *grow* under uniform
refinement, the signature of a singularity still being resolved). A three-region
box study on the self-similar ladder localizes what is left, and it is **not**
mid-span/LE resolution (an earlier reading of the spanwise "slosh", since
retracted — that slosh was an artifact of comparing two different mesh families,
see below):

| region | coarse_ss → medium → fine | p | verdict |
|---|---|---|---|
| **tip-cap edge (WALL)** | 0.662 → 0.824 → **1.015** | **+0.321** | **DIVERGES** |
| wake free edge (the §4.2 fix) | 0.536 → 0.565 → 0.570 | +0.045 | bounded |
| wing p99 (control) | 0.642 → 0.628 → 0.629 | −0.014 | bounded |

The §4.2 fix holds and the ordinary wing field is converged; what diverges is the
sharp convex edge where `meshgen/wing3d.py`'s **flat tip cap** meets the
upper/lower surfaces. The real ONERA M6 has a **rounded** cap — the flat one is a
documented deliberate simplification — so this is an edge singularity of exactly
the kind §4.1–§4.2 treat, only on the **body** instead of the wake. **It is not a
P11 (curved wall *element*) problem:** isoparametric elements cannot regularize a
genuinely sharp geometric *edge*; the geometry itself is what is wrong, so the fix
is Track M (round the cap). Evidence:
`cases/demo/p13_tip_edge_singularity/run_g133_ladder.py` (5/5).

**Mesh-ladder prerequisite (fixed 2026-07-13).** Any M6 three-point study must run
on `generate_onera_m6.RICHARDSON_LADDER = (coarse_ss, medium, fine)`. The shipped
`coarse` is **not on the refinement ray**: `_level_params` derived every length
from `h_wall` except the far field, which was clamped (`h_far = min(2.5,
120·h_wall)`), and the clamp bit only at `coarse` — so coarse→medium refined the
far field 1.39× against the wall's 2×. `coarse_ss` is that level re-cut without
the clamp; the ladder now refines by exactly 2.000 in every length scale. This
invalidates every earlier M6 three-point Richardson, **including P9/G9.1's**.

---

## 5. Boundary conditions

| Boundary | Condition | Implementation |
|---|---|---|
| Solid wall (wing, body) | ρ ∂φ/∂n = 0 | Natural (do-nothing) in the Galerkin weak form — zero surface integral. No penalty, no ghost cells. |
| Far field | φ → φ∞ + φ_vortex | Dirichlet on outer nodes: φ∞ = x cosα cosβ + y sinβ + z sinα cosβ (uniform flow at incidence α, sideslip β), plus a compressible horseshoe-vortex correction with total circulation ΣΓ_j Δs_j (Prandtl–Glauert-scaled). At R ≳ 25–50 chords the correction is small but including it lets the domain shrink to ~15 chords. |
| Wake | (4.1)–(4.2) | Constraint elimination, §4. |
| Symmetry plane (half-model) | ∂φ/∂n = 0 | Natural, free. |

One subtlety: with all-Neumann walls and Dirichlet far field the system is
well-posed; if a pure-Neumann variant is ever used (e.g. channel flows), pin
one node.

**Designed extension — Track V transpiration BC (not started).** The VII
coupling (Drela IBL3, roadmap.md Track V) replaces the homogeneous wall
condition by a blowing source ρ ∂φ/∂n = d(ρ_e u_e δ*)/ds on wall facets and adds
the analogous δ*_wake mass source on the wake sheet — RHS-only at the Picard
level, no mesh motion, δ* = 0 bit-identical to the inviscid path (design: DN2
§3.2/§7–8, DN6 §8–9).

### 5.1 Boundary-flux correction for the wall geometric error (G1.6 candidate fix routes; gates G1.3–G1.5 + DP1)

Context (gate G1.6, formerly G1.2, incompressible sphere): the medium-mesh Cp error is ~11.6%
against a <2% target, root-caused (see PROJECT_STRUCTURE.md "Known gaps") to a
geometric/variational inconsistency, not to the surface gradient recovery. After
integration by parts, the P1 Galerkin wall term ⟨v, ∇φ·ñ⟩ is dropped as a
natural BC — which enforces "zero flux through the **flat facet** (normal ñ)",
whereas the physical condition is "zero flux through the **true curved surface**
(normal n)". The facet normal deviates from the true normal by O(h), producing a
first-order geometric error in wall velocity/Cp; the raw nodal potential shows
sub-first-order wall convergence, consistent with this mechanism. The previously
recorded conclusion was "a true fix needs curved/isoparametric wall elements, a
separately-scoped effort". The literature survey below found intermediate routes
with a far smaller footprint; they should be verified first, before deciding
whether the curved-element effort is still needed. Options are ordered by
implementation footprint. (Standing prohibitions remain: do **not** re-propose
further h-refinement or recovery-scheme tweaks — both are ruled out with
evidence.)

> **Verification outcome (2026-07-06, see §5.1.2):** the G1.3/G1.4 oracle
> experiments ruled Option A out — on body-fitted meshes the boundary-data
> defect it corrects is (near-)zero and the measured ceiling is ≈ 11.3% vs
> the < 2% target. DP1 took the "> 5%" branch: Option C + separately-scoped
> curved elements. §§ Option A/B below are kept as the design record.

**Option A (recommended, implement first — RULED OUT by the G1.3/G1.4
oracles, see §5.1.2): true-normal weak-flux correction
(lagged flux correction).**

Provenance (two independent lines of work, same practical recipe):

- Krivodonova & Berger (JCP 2006): on straight-sided meshes, impose the solid
  wall condition using the normal of the *physical* geometry rather than the
  *computational* geometry; for the Euler equations this dramatically improves
  solution quality without curved meshes. Ciallella, Gaburro, Lorini &
  Ricchiuto (Appl. Math. Comput. 2023) extend the same polynomial correction to
  general 2D/3D boundary conditions.
- The Shifted Boundary Method (Main & Scovazzi 2018; the Atallah, Canuto &
  Scovazzi analysis series) and the earlier Bramble–Dupont–Thomée (1972)
  boundary-value corrections: correct the boundary-condition **data** onto the
  approximate boundary so as to cancel the geometric consistency error.

Mathematical form. At a point x on a wall facet, let p(x) be the closest-point
projection of x onto the true surface (the sphere), n = n(p(x)) the true unit
outward normal, and ñ the facet unit outward normal. Decompose

    ñ = (ñ·n) n + t,    t := ñ − (ñ·n) n    (t lies in the true tangent plane, |t| = O(h))

so that ∇φ·ñ = (ñ·n)(∇φ·n) + ∇φ·t. The physical boundary condition
∇φ·n(p(x)) = 0 eliminates the first term; the wall term is no longer zero but
⟨v, ∇φ·t⟩. Move it to the right-hand side in a **lagged (Picard)** fashion:

    ∫_Ω ∇φ^{k+1}·∇v dV = ⟨ v,  ∇φ^k · ( ñ − (ñ·n) n ) ⟩_{Γ_h,wall}
    (all other boundary conditions unchanged)

Implementation notes (binding for the eventual implementation):

- The stiffness matrix is **completely unchanged**: it stays SPD, and the AMG
  hierarchy and preconditioner are reused as-is. The only change is one
  RHS-assembly loop over wall facets.
- New geometric data required: closest-point projections and true normals at
  facet quadrature points. For the sphere these are analytic,
  n = (x − c)/|x − c|; for the wing phase later they come from CAD/analytic
  geometry or a fine reference surface — design the interface as a replaceable
  `closest_point_normal(x)` callback, precomputed into SoA arrays (per the
  agent-rules Numba hard constraints: no Python-object operations in hot loops).
- Quadrature: 3-point Gauss on each facet (or whatever rule the existing surface
  integrals use); ∇φ^k is the piecewise-constant gradient of the adjacent
  element.
- Fixed-point loop: the correction is O(h), so the iteration is contractive. For
  the pure-Laplace validation, run a fixed 3–5 outer iterations (repeated solves
  with the same matrix — cheap); once in the full-potential regime, fold it into
  the existing ρ-Picard loop at zero marginal cost.
- Freestream-preservation check: uniform flow φ = U∞·x does not satisfy the
  sphere wall condition, so a nonzero correction term on wall-bearing meshes is
  expected; but it must be confirmed that the V0 freestream gate
  (`tests/test_v0_freestream.py`, which uses wall-free / all-far-field
  configurations) never triggers this correction path.
- Prerequisite: consistent wall-facet winding — the existing winding assert in
  `_wall_vertex_normals` is a precondition for this correction and must not be
  removed.

#### 5.1.1 Cylinder pre-study (gate G1.3, formerly G1.2-a0)

`cases/meshes/cylinder_2.5d/` is the **designated rapid testbed** for the
Option A/B route, to be exercised before the sphere (gate G1.3 precedes
G1.4). Rationale: it exhibits the **same** curved-wall variational crime,
already quantified — max |Cp err| 0.091 (coarse) → 0.045 (medium), ~O(h)
(`tests/test_m0_cylinder.py`); every geometric ingredient Option A needs is
available in closed form; the meshes are cheap (6.9k / 17.3k tets); and the
diagnostic is a one-dimensional curve Cp(θ).

Closed forms for the cylinder (radius a, axis z, freestream U along x, with
r_xy² = x² + y²):

    φ_exact = U x (1 + a²/r²),    p(x) = a (x, y, 0)/r_xy,    n = (x, y, 0)/r_xy

The `closest_point_normal(x)` callback of §5.1 gets one analytic implementation
each for the cylinder and the sphere (n = (x − c)/|x − c|); the interface design
itself is unchanged.

Caveats (all three are binding; none may be dropped when citing this pre-study):

1. **Necessary, not sufficient.** The cylinder has single curvature, the sphere
   double curvature. **The G1.6 gate closes only on the sphere.** Cylinder
   results serve solely as prerequisite evidence for entering the sphere
   experiment G1.4.
2. **Spanwise-noise floor.** The quasi-2D mesh carries the O(h) spanwise noise
   inherent to the 3-tet prism split (max |w|/U∞ ≈ 2.9e-2 on coarse), which
   pollutes in-plane gradient recovery at the same O(h). The cylinder
   acceptance criterion is therefore: corrected error significantly below the
   uncorrected one **and** the Cp-error convergence order recovering from
   sub-first to ≈ first order — **no absolute threshold**. The oracle run
   measures this floor as a by-product; its magnitude feeds back into the
   G2.5(b) re-spec.
3. The cylinder suction peak Cp = −3 is stronger than the sphere's −1.25 — a
   harsher stress test; percentage figures must not be compared across the two
   geometries.

Visualization: on the single-layer mesh, slicing degenerates to a trivial
operation — the symmetry plane carries its own 2D triangulation, so fields plot
directly via tripcolor. The helper is to be defined with the final interface
prototype of P2's `post/section_cut.py` (signature reserves the z = const
parameter); P2 then adds the general 3D interpolation path.

Distinction from the rejected Nitsche/penalty prototype (record this explicitly
to prevent misclassification as a repeat experiment): the earlier attempt
changed the **enforcement mechanism** for the same condition on the same wrong
(flat-facet) geometry; Option A corrects the boundary condition's **data
itself** — the closest-point projection brings in the true normal. That is the
substantive contribution of the "shift" in SBM, independent of the enforcement
mechanism.

Theoretical expectations and ceiling (recorded to manage expectations):

- The SBM analyses show that for P1 elements — piecewise-constant gradient, zero
  Hessian, hence no second-order Taylor expansion available — the naive shifted
  Neumann condition loses one order in L², but the H¹ seminorm retains its
  optimal first order. Cp is a gradient quantity controlled by H¹, so this
  ceiling is not an obstacle for this project.
- This project uses **body-fitted** meshes (vertices lie on the true surface):
  the geometric gap has thickness O(h²), far better than the O(h) gap of
  unfitted SBM; the first-order normal-rotation correction is therefore the
  dominant error term.
- Expectation: Cp recovers close to first-order convergence, with a good chance
  the medium-mesh error drops below 2%; the exact ceiling is measured by the
  oracle experiment (gate G1.4).

#### 5.1.2 Pre-study outcome (2026-07-06): Option A ruled out on body-fitted meshes

The G1.3 cylinder pre-study and the G1.4 sphere oracle run the same day
falsified the expectation above (full evidence: roadmap G1.3/G1.4 entries;
`tests/test_wall_correction_cylinder.py`; `artifacts/G1.3/`; oracle results in
`cases/demo/p1_laplace/results/`, demo absorbed there 2026-07-07). The
mechanism, in hindsight elementary:

- For a **harmonic** potential and **body-fitted** wall vertices, the region
  between a flat facet and the true curved surface is closed by exactly those
  two surfaces; the true surface carries zero flux, so by the divergence
  theorem the exact solution's **net flux through every facet is exactly
  zero**. The natural-BC consistency defect ⟨v, ∇φ_exact·ñ⟩ is only a
  first-moment residue — measured ~2e-5 max on the coarse cylinder, shrinking
  ~O(h⁴). The SBM/BDT corrections earn their keep on *unfitted* boundaries
  with O(h) gaps; with an O(h²) body-fitted gap there is nothing to correct.
- Option A's t-form ⟨v, ∇φ_exact·t⟩ additionally assembles to machine zero on
  the cylinder (exact adjacent-facet cancellation on the uniformly spaced
  circle); the assembly was verified against a hand-computed single-facet
  case, so the zero is real.
- Measured ceilings: cylinder Cp error unchanged
  (9.10e-2/4.49e-2/2.22e-2 max over coarse/medium/fine, slope 1.02, corrected
  ≡ uncorrected); sphere medium-mesh max |Cp err| 0.1156 → 0.1133 with the
  full consistency defect restored — an **11.3% ceiling** vs the < 2% target.
- Bonus finding that re-frames the cylinder case: ~76% of its Cp error at
  every level is reproduced by feeding the *exact* potential through the
  surface recovery — on the quasi-2D single-layer strip the sliver-shaped
  wall triangulation dominates the recovery error, and the wall nodal φ
  converges at a healthy ~1.2 order. The cylinder therefore does **not**
  exhibit the sphere's sub-first-order pathology and is no longer the
  designated testbed for it (§5.1.1's "same variational crime" rationale is
  superseded on this point; the meshes and the analytic-Cp validation remain
  in use as the M0 end-to-end check).

Consequences (DP1 "> 5%" branch, recorded in the roadmap): the sphere's
remaining error is attributed to the domain perturbation (the missing
facet/surface slivers in ∫_Ωh vs ∫_Ω) plus P1 gradient approximation — i.e.
exactly the terms Gap-SBM's gap integrals model (Option B below), which is
why B survives as optional pre-study material for the curved-element design
pass; but with O(h²) gap thickness its payoff is expected second-order-small.
The accuracy route is curved/isoparametric wall elements as a
separately-scoped effort — **now the P9 phase** (G9.1 sphere-Cp < 2%; it also
owns the residual V6 < 1% floor, G9.2), and gate G1.6 is to be redefined per
Option C.

**Option B (escalation if Option A falls short — per §5.1.2 now optional
pre-study material only): Gap-SBM gap correction.**

Collins, Li, Lozinski & Scovazzi, "Gap-SBM" (arXiv:2508.09613, 2025) targets
exactly the P1 Neumann suboptimality above. It builds an approximate gap
geometry from the distance map between the true and computational boundaries,
extends the solution and test functions into the gap, and applies approximate
quadrature to the corrected variational form; no extra degrees of freedom, no
cut cells, no ghost-penalty terms, with proven optimal L² and H¹ convergence for
the Neumann problem. The concrete footprint is: wall-facet surface integrals
multiplied by a gap-thickness coefficient, plus a few distance-vector correction
terms at wall nodes — implementable as SoA arrays + a facet loop, Numba-friendly.
Caveats: the paper's analysis is 2D (the authors state 3D is a direct
extension), and in the body-fitted case the gap coefficient is O(h²), so these
extra terms are small to begin with — which is precisely the rationale for
"Option A first".

**Option C (pragmatic fallback): redefine the gate rather than the scheme.**

If curved elements ultimately still require their own effort, redefine the G1.6
acceptance criterion as "Cp error / convergence order relative to a
**geometry-consistent reference solution** (a high-accuracy reference on the
same polyhedral domain, e.g. BEM or an ultra-fine mesh)", stripping the
geometric model error out of the code-correctness verification. This conforms to
the §10 validation-ladder principle of not confounding model error with code
error. This option does not improve physical accuracy; it only adjusts the
verification yardstick.

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

**Newton Jacobian (needs the P7 differentiable flux; the Newton solve itself is
P8).** Differentiating (6.1) w.r.t. φ_k gives López's two-term form
(dissertation Eq. 3.24, §13):

    ∂R_i/∂φ_k = Σ_e V_e [ ρ̃_e ∇N_i·∇N_k                     ← Term 1 (Picard, SPD)
                + (∂ρ̃_e/∂φ_k) (∇φ_e·∇N_i) ]                  ← Term 2 (density)   (6.3)

Term 1 is exactly the frozen-density Picard matrix (6.2) — the existing colored
assembler produces it verbatim. Term 2 carries the density sensitivity; expanded
through the isentropic law and the switching function (López Appendix B, verified
against the PDF) it splits by flow regime, and in the supersonic accelerating
case reads

    ∂ρ̃_e/∂φ_k     = (∂ρ_e/∂φ_k)(1 − μ) − (ρ_e − ρ_up)(∂μ/∂φ_k)          (B.3)
    ∂ρ̃_e/∂φ_k^up  = +μ (∂ρ_up/∂φ_k^up)                                   (B.4)

with ∂ρ/∂φ = (∂ρ/∂u²)(∂u²/∂φ), ∂ρ/∂u² = −(ρ∞/2a∞²)[…]^{(2−γ)/(γ−1)} (B.8),
∂μ/∂φ = (**+**M_c²/M⁴)(∂M²/∂u²)(∂u²/∂φ) (B.9–B.12; **sign corrected 2026-07-10
by the P7 G7.3 FD gate** — for μ = 1 − M_c²/M², dμ/dM² = +M_c²/M⁴, and with
dM²/du² > 0 the chain is positive (μ grows with speed, as it must); the earlier
"−" here was a transcription typo. The shipped derivative
(`kernels/upwind.py::rho_tilde_sensitivities_sweep`, P7) uses the + sign and
matches a central difference of the shipped flux at frozen selection to
~4e-10 in every regime, `tests/test_p7_diff_flux.py`). **(B.4) is nonzero**: the
current element's residual depends on the *upstream* element's DOFs, so the
Newton stencil is **wider than Picard** — the sparsity map and element coloring
must be rebuilt from the element + upstream graph.

**Stencil width is NOT "one element-layer" in UP3D (correction 2026-07-08).**
López's *single-hop* upstream gives a one-layer-wider stencil; UP3D reaches
farther because the sliver prism-split tets need multi-hop reach (§3: a single
hop reaches only ~0.37 of the streamwise extent). The **default Newton flux is
the frozen walk** — differentiate the shipped upwind density at fixed selection,
NOT a replacement flux:
- **Frozen walk (default)** (`mode="walk"`, u(e) a single element ≤4 hops away):
  Term 3 couples e to exactly *one* upstream element — nnz increase is modest
  (~+1 element/row, closest to López) — but at graph-distance ≤4, so the
  coloring/CSR graph gains long-range e→u(e) edges (higher chromatic number).
  Differentiable at frozen selection as-is (isentropic ρ + branch-wise ∂ν/∂φ,
  B.3–B.6); `max_ε` optional only to damp active-set churn far from the solution
  (López's C⁰ `max` already reaches quadratic near it, where the selection and
  active set stop changing). This is the P7 deliverable.
- **Streamline-Gaussian kernel (optional)** (`mode="kernel"`): ρ_up is a weighted
  sum over the *whole depth-3 BFS neighbourhood*, so Term 3 couples e to *every*
  cell in that neighbourhood — a materially **denser** Jacobian (row nnz
  ~O(neighbourhood size), not +1); the weight dependence w_f(∇φ_e) adds only
  self-coupling. Smoother and frozen-selection-free, but a denser system — a
  Picard-speed path, promoted to the Newton flux only if it measures favourable.

Retire the old "~30 % memory / one layer" estimate; take a **measured**
nnz/GMRES-iteration/AMG-setup figure in N2/N3 before committing the wider
coloring/CSR rebuild.

Two facts from the PDF that shape the implementation:

- López keeps the switching-function derivative ∂μ/∂φ (B.3, B.6) — this is a
  *full* Jacobian ("strategy A"); dropping it (Term-1-plus-local-density-only)
  degrades quadratic to superlinear.
- The derivation treats the **upstream selection u(e) as frozen** while
  differentiating ρ_up; the selection is re-evaluated each Newton step but not
  differentiated. This is exactly why P6 need only make the flux differentiable
  at fixed selection (§3.1), not differentiate the selection itself.

The exact Jacobian (Term 1 + Term 2 + upstream coupling) is **nonsymmetric** and
indefinite in supersonic zones — GMRES + AMG (aggregation on the symmetric part)
or GMRES + ILU, not CG (§8).

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
│   ├── newton.py        # Phase 8: exact Jacobian, line search / pseudo-transient
│   └── continuation.py  # Mach/α ramping for hard cases
├── post/
│   ├── surface.py       # nodal Cp, sectional cl, forces & moments
│   ├── trefftz.py       # induced drag from Γ(s) (far-field), cross-check vs near-field
│   └── vtk_out.py       # meshio VTK export: φ, V, M, ρ, Cp
└── cases/               # regression decks: sphere, naca0012_extruded, oneraM6
```

---

## 8. Nonlinear solution strategy

**Baseline: relaxed Picard (Phases 1–5).**

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

### 8.1 Fully-coupled Newton with Γ as an unknown (P8)

The Picard/secant architecture nests three loops (Mach continuation → per-station
Γ secant → density Picard). The P5 medium failure root-caused the fragility to
**secant–density coupling**: the secant reads its Kutta mismatch off a not-yet-
converged density field, so a bad slope estimate over-shoots Γ, which forces more
density work, which pollutes the next mismatch — a positive feedback that
diverges when pushed (the P5 fix regularized it by early-stopping + a fixed-Γ
polish, not by curing it). López avoids this class of failure entirely by solving
one Newton system in which **Γ is an unknown alongside φ**, so φ and Γ converge
together and there is no secant reading a lagged residual.

UP3D's master–slave reduction makes this cheap because the Γ-Jacobian blocks are
almost already in the code. With φ_full = T·φ_red + g(Γ) (the wake jump) and the
far-field Dirichlet data also carrying Γ through the vortex correction, the
coupled step is

    ⎡ TᵀJ T        ∂R_red/∂Γ ⎤ ⎡ δφ_red ⎤     ⎡ −Tᵀ R      ⎤
    ⎢                        ⎥ ⎢        ⎥  =  ⎢            ⎥
    ⎣ ∂F/∂φ_red    −I        ⎦ ⎣ δΓ     ⎦     ⎣ −F         ⎦

where
- `TᵀJ T` is the reduced full Jacobian (6.3), Term 1 reusing `assemble_matrix`;
- `F_j(φ_red) = kutta_targets_j(φ_red) − Γ_j` is the Kutta closure residual
  (`wake.py::kutta_targets`, reusable verbatim), with ∂F_j/∂φ_red the sparse
  averaging matrix (±1/n_j on this station's upper/lower TE probe nodes);
- **∂R_red/∂Γ_j = TᵀJ g_j** (g_j = station-j slave indicator). At the Picard
  level TᵀA g_j = h_j, which `wake.py` **already materializes** as `self._h`; the
  exact Newton column uses the full J and can be lagged to A as a documented
  approximation.
- **The far-field column.** Γ enters φ_full a *second* way: the far-field vortex
  Dirichlet values are linear in Γ (`dirichlet.py`, `farfield_spanwise_gamma`
  taper or the span-uniform mean). So ∂R_red/∂Γ_j gets an extra
  −A_coupling·(∂vals_red/∂Γ_j) term. **This term is easy to miss** — it is not in
  the Picard code path (where it is folded silently into the RHS) — and must be
  added to the Newton column.

Globalization (López-verified, dissertation §4.4/§4.8): plain Newton from the
freestream diverges transonically, so wrap it in **Mach continuation + load
stepping**. Two corrections to the earlier internal schedule proposal, verified
against the PDF (Tables 4.7/4.8/4.13):

- Within a single case López holds M_crit and μ_c **fixed** and ramps only M∞;
  the harder the case, the *lower* M_crit and *higher* μ_c chosen for the whole
  ramp (Case 1: M_crit 0.99, μ_c 1.0; Case 4 at α=2°/M∞0.75: M_crit 0.90,
  μ_c 1.1). There is **no** per-step M_crit=0.99→0.90 sweep.
- ONERA M6 (Table 4.13): 12 steps, M_crit held **constant at 0.95**, μ_c held at
  2.0 while M∞ ramps 0.50→0.84 (steps 1–8), then M∞ fixed and μ_c **decreased**
  2.0→1.6 (steps 9–12) to sharpen the shock after the field is established. So
  μ_c scheduling is a *post-target dissipation reduction*, not a during-ramp
  increase.

Add an Eisenstat–Walker forcing schedule for the inner GMRES (loose early, tight
near convergence — the existing `forcing` param is a one-parameter stand-in), and
reuse the AMG hierarchy across a few Newton steps. Expected cost from López's
data: subsonic 5–10 Newton iterations total; transonic 4–9 per load step; ONERA
M6 ≈ 12 steps × 5–9 ≈ 60–110 total Newton iterations (Tables 4.5/4.6/4.9,
§4.7) — versus the ~10⁴ Picard iterations today.

An N-Γ split (keep the Γ secant outer loop, Newton only the density inner solve)
is retained **only as a fallback**: it preserves the secant and therefore the P5
coupling-instability risk.

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

### 9.1 Wall-Cp gradient recovery — the G6.1 sawtooth fix (P6)

The wall velocity is the in-plane tangential gradient of φ (natural BC, §5). The
default recovery is **per-triangle constant** (`triangle_tangential_gradients` /
`wall_cp_curve` / `wall_force_coefficients`, no averaging) — chosen so the sharp
TE crease is not polluted by mixing the two sides' tangent planes. Its cost on
the prism-split sliver wall triangulation is the ≈2-cell **surface-Cp sawtooth**:
adjacent sliver triangles give oscillating constant gradients even for a smooth
φ. This was misattributed to the artificial-density flux; it is a **recovery
artifact** (N1, 2026-07-08 — §3.1: nodal averaging of the *same* solution drops
the G6.1 metric ~330×; a smoother flux does not reduce it).

**Fix (gate G6.1): `smooth_wall_tangential_gradients` (`post/surface.py`),
threaded as `smooth_passes` through the Cp curve and the force integral.** A few
Jacobi passes average each wall triangle's gradient with its edge neighbours,
**gated by outward-normal alignment** (`n̂_out·n̂_out(nb) > cos_thresh`) so it never
averages across the sharp TE (the two TE triangles' outward normals are nearly
anti-parallel). Measured on the coarse NACA0012 M0.80 walk solution, one pass:
G6.1 metric 0.0758 → 0.00023 (330×, reversals 39 → 1), shock 0.604 → 0.607,
cl_p within 0.3 %, |Cp|_max unchanged (TE preserved) — and it is a fixed linear
operator in the gradient (differentiable). `smooth_passes = 0` is the
bit-identical default; the near-field cd_pressure shifts (~15 % on coarse) but is
the explicitly-untrusted FP quantity (use Trefftz for drag).

**Use it for the Cp curve, not the loads (G6.3, measured 2026-07-08).** On
ONERA M6 coarse the V6 consistency |CL_p − CL_KJ|/CL_KJ *worsens* with smoothing
(2.40 % → 3.35 % → 3.88 % for `smooth_passes` 0/1/2): the ±sawtooth largely
cancels in the surface integral, and the averaging instead smears the true LE
suction peak, moving CL_p ~1 % further below the trustworthy CL_KJ. So the whole
V6 floor is the sharp-TE/LE P1 wall gradient (→ P9 curved elements), not the
sawtooth. Recommendation: `smooth_passes > 0` for the reported Cp curve;
`smooth_passes = 0` for `wall_force_coefficients` (the param stays, opt-in, but
raw loads are more accurate).

---

## 10. Verification & validation ladder

Each item is a pytest in `cases/`; keep them runnable on coarse meshes in CI.

| # | Case | Checks | Gate |
|---|---|---|---|
| V0 | Freestream preservation, arbitrary tet mesh + wake cut | ‖R(φ=x)‖∞ | < 1e−12 |
| V1 | Laplace MMS (ν=0, ρ=1), manufactured trigonometric φ (implemented: sin πx · sin πy · sin πz) | L2 error vs h | slope ≥ 1.9 |
| V2 | Incompressible sphere | Cp vs 1 − (9/4)sin²θ | max err < 2% (medium mesh) |
| V3 | Subcritical extruded NACA0012, M 0.5, α 2° (span-periodic or large AR + symmetry) | Cp vs 2D FP/panel reference; cl | Δcl < 2% |
| V4 | Transonic NACA0012 M 0.80 α 1.25° (extruded) | shock x/c vs published FP (Holst AGARD FP workshop results — note FP ≠ Euler shock position) | shock Δx/c < 0.03 |
| V5 | ONERA M6, M 0.84, α 3.06° | λ-shock topology, section Cp at η = 0.44/0.65/0.90 vs FP references (FLO-class / TRANAIR); CL | qualitative λ + CL within FP-literature scatter |
| V6 | Γ-consistency | **2.5D (G2.4):** sectional cl = 2Γ/(Uc) vs surface-pressure integration. **3D (G5.2):** total CL_p vs the *spanwise-integrated* Kutta–Joukowski CL_KJ = 2∫Γ(z)dz/(U·S). | < 1% |

V6 caveat (2D vs 3D): the *sectional* identity cl(z) = 2Γ(z)/(U c) holds only
in 2D / quasi-2D flow. On a finite 3D wing the local flow at a section is not
2D (spanwise flow + induced downwash), so the pressure-integrated sectional cl
does **not** equal 2Γ(z)/(U c) station-by-station — only the *span-integrated*
total lift is exact: L = ρ∞U∞∫Γ(z)dz (far-field/Trefftz-plane Kutta–Joukowski).
So V6 is checked per-section only for the extruded quasi-2D case (G2.4,
`sectional_cl_from_gamma`); on ONERA M6 it is the integrated `cl_kj_3d`
(`post/surface.py`) vs total CL_p. The M6 spanwise `2Γ/c` plot is a *loading
distribution*, not a per-section consistency gate.

Reference data caveat: validate against *full-potential* published results
first (Holst's Progress in Aerospace Sciences 2000 review collects them), Euler
second — otherwise model error and code error are confounded.

---

## 11. Development roadmap (vibe-coding phases)

> Phase numbering follows [docs/roadmap.md](roadmap.md) Track P, the *active*
> tracker (gate checklists, progress ledger, and the parallel mesh Track M live
> there). This section summarizes phase content and maps it onto the §10
> verification ladder; for detailed gates and current status, follow roadmap.md.

Each phase is a self-contained PR-sized unit with its gate from §10.

- **P0 — Repo scaffolding + mesh infrastructure.** meshio reader, metrics
  (B_e, V_e), adjacency, coloring, VTK writer. Gate: metrics unit tests
  (ΣV_e = volume of unit cube; ∇(linear field) exact), coloring validity assert.
- **P1 — Laplace solver.** ρ ≡ 1, Dirichlet far field, natural walls, CG+AMG.
  Gates: V0 (without wake), V1, V2.
- **P2 — Lift: wake cut + Kutta, on Laplace (★ critical phase).** Node
  duplication, constraint elimination, Γ update, vortex far-field correction
  (incompressible form) — all the hard topology/constraint machinery lands
  against the linear operator. Gates: V0 *with* wake cut, V3 (incompressible
  variant), V6.
- **P3 — Subsonic compressible.** Density law + Picard loop with
  under-relaxation, no upwinding; Prandtl–Glauert-scaled vortex far field.
  Gates: sphere at M∞ = 0.3 vs Prandtl–Glauert-corrected V2; V3; convergence
  in <30 Picard its; P1/P2 gates stay green (ν ≡ 0 path identical to Laplace).
- **P4 — Transonic: artificial density.** Upwind element search, ν switch,
  ρ̃; relaxation + Mach continuation. Gate: V4.
- **P5 — 3D validation: ONERA M6.** Requires the swept-wing mesh (roadmap.md
  Track M1). Gates: V5; V6 consistency in 3D (the *span-integrated* CL_KJ =
  2∫Γdz/(U·S) vs total CL_p — the per-section 2Γ/c identity is 2D-only, see §10).
> Track-P renumber (2026-07-08): the old P6 "differentiable flux (remove the
> sawtooth)" split, once the sawtooth was root-caused to *recovery* not *flux*,
> into P6 (recovery) + P7 (flux) + P8 (Newton, was P7) + P9 (curved walls, new)
> + P10 (backlog, was P8). Old gate IDs map G6.x → G6.x/G7.x, G7.x → G8.x.

- **P6 — Surface-pressure recovery (remove the Cp sawtooth)** (§3.1, §9.1).
  The sawtooth is a per-triangle wall-gradient **recovery** artifact on the
  sliver wall triangulation, NOT the artificial-density flux (proven: nodal
  smoothing of the same solution drops the metric ~330×; a smoother flux does
  not). Fix = normal-gated recovery smoothing in post-processing
  (`smooth_passes`, §9.1), preserving the sharp TE, applied to the Cp curve and
  the force integral. **Done** (2026-07-08). Gates G6.1–G6.4.
- **P7 — Differentiable artificial-density flux at frozen selection** (§3.1/§6.3;
  re-scoped 2026-07-08). The **P8 Newton prerequisite** = the shipped **walk**
  flux's ∂ρ̃/∂φ made differentiable at frozen u(e) (sparse, ~+1 element/row,
  closest to López); deliverable = derive + FD-verify it (`max_ε` optional). The
  streamline-Gaussian kernel (Eq. 3.4′, §3.2) is an **optional** ~10×-faster
  Picard flux with a denser Newton Jacobian, a P8 flux candidate only if measured
  net-favourable. Neither fixes the sawtooth (that is P6). Gates G7.1–G7.3.
- **P8 — Performance & robustness: fully-coupled Newton** (§3.2, §6.3, §8.1).
  Consumes the P7 flux. Full Jacobian with the nonzero upstream coupling (López
  Eq. B.4), fully-coupled (φ, Γ) solve replacing the P5-fragile Γ-secant, Mach
  continuation + load stepping, GMRES + AMG, Eisenstat–Walker, AMG reuse,
  profiling (target: ONERA M6 medium mesh < 5 min single node). Gates G8.1–G8.3.
- **P9 — Curved / isoparametric wall elements** (§5.1). The shared accuracy route
  for G1.6 sphere-Cp < 2% (Option C / DP1 "> 5%" branch) and the residual
  V6 < 1% floor after P6 (sharp-TE/LE P1 wall gradient). Gates G9.1–G9.2.
- **P10 — Extensions (backlog).** Mixed prism/tet; embedded-boundary wake
  alternative; VII coupling hook (transpiration BC ∂φ/∂n = d(u_e δ*)/ds —
  reuses the IBL work from pyTSFoil); adjoint via the Newton Jacobian
  transpose (nearly free once (6.3) exists — high value for the MDO thread).

Two further tracks are **designed but not started** (2026-07-07..10; roadmap.md
holds the authoritative phase tables and sequencing guards):

- **Track B — level-set embedded wake** (B1–B5; B6 curved/free wake shelved):
  replaces the §4 conforming wake with a level-set + multivalued-FE + penalty
  Kutta path; coexists with the shipped path, blocks nothing in P7–P10 (§4 note;
  DN1).
- **Track V — viscous–inviscid coupling** (V1 loose → V3 tight Newton; V4 wake
  IBL rides on V1): Drela IBL3 6-equation surface FE + the §5 transpiration BC;
  V1 needs only P6, V3 consumes P8. Phase IDs V1–V4 are *phases*, distinct from
  the §10 validation-case IDs V0–V6 (DN2/DN6).

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

Boundary-flux correction for curved walls on flat-facet meshes (§5.1):

- Krivodonova, L., Berger, M., "High-order accurate implementation of solid
  wall boundary conditions in curved geometries," *J. Comput. Phys.* 211
  (2006) 492–512.
- Ciallella, M., Gaburro, E., Lorini, M., Ricchiuto, M., "Shifted boundary
  polynomial corrections for compressible flows: high order on curved domains
  using linear meshes," *Appl. Math. Comput.* (2023).
- Main, A., Scovazzi, G., "The shifted boundary method for embedded domain
  computations. Part I: Poisson and Stokes problems," *J. Comput. Phys.* 372
  (2018) 972–995.
- Atallah, N.M., Canuto, C., Scovazzi, G., "The high-order Shifted Boundary
  Method and its analysis," *CMAME* 394 (2022) 114885.
- Collins, J.H., Li, C., Lozinski, A., Scovazzi, G., "Gap-SBM: A New
  Conceptualization of the Shifted Boundary Method with Optimal Convergence
  for the Neumann and Dirichlet Problems," arXiv:2508.09613 (2025).
- Bramble, J.H., Dupont, T., Thomée, V., "Projection methods for Dirichlet's
  problem in approximating polygonal domains with boundary-value corrections,"
  *Math. Comp.* 26 (1972) 869–879.
- Burman, E., Hansbo, P., Larson, M.G., "A cut finite element method with
  boundary value correction," *Math. Comp.* 87 (2018) 633–657.

FE full-potential, embedded wake, and fully-coupled Newton (§3.1–3.2, §6.3, §8.1):

- López Canalejo, I., *A Finite-Element Transonic Potential Flow Solver with an
  Embedded Wake Approach for Aircraft Conceptual Design*, PhD dissertation,
  Technische Universität München, 2021 (`docs/references/Dissertation_Inigo_Lopez.pdf`)
  — the closest prior art for the FE artificial-density flux and a full-Jacobian
  Newton solver. §3.3 artificial compressibility (Eq. 3.19–3.27, switching
  function μ_s = μ_c·max(0,μ,μ_up)); §3.4 limit velocity / M²_max ~ O(3) clamp;
  §3.5 embedded wake (2D/3D conditions, TE and small-cut treatment); Ch.4
  Newton convergence (Tables 4.5/4.6/4.9, strict quadratic) and load-stepping
  schedules (Tables 4.7/4.8 NACA0012, 4.13 ONERA M6); ONERA M6 loads Table 4.15
  (CL 0.288 KRATOS = Tranair); Appendix B Eq. (B.1)–(B.17) full sensitivity
  derivation incl. the nonzero upstream coupling (B.4/B.6) and switching-function
  derivatives (B.3/B.6). NACA0012 validation used a **sharp** TE (Eq. 4.2 modifies
  the last thickness coefficient to close the otherwise-blunt profile) with full
  quadratic Newton convergence — i.e. a sharp 2D TE is not itself a Newton
  obstacle; only the 3D wingtip-TE geometric singularity needs the density clamp.
- López Canalejo, I., Núñez, M., Baiges, J., Rossi, R., "An embedded approach for
  the solution of the full potential equation with finite elements," *CMAME* 388
  (2022) 114244 — the journal condensation of the dissertation Ch.3–4.
- Cai, X.-C., Keyes, D.E., Young, D.P., "Parallel Newton–Krylov–Schwarz
  algorithms for the transonic full potential equation," *SIAM J. Sci. Comput.*
  19 (1998) 246–265 — inexact Newton + GMRES for the FPE.
- Eisenstat, S.C., Walker, H.F., "Choosing the forcing terms in an inexact Newton
  method," *SIAM J. Sci. Comput.* 17 (1996) 16–32.
