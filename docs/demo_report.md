# Phase demo report — evidence for completed phases

**Scope.** One self-contained demo case per completed roadmap phase (Track P:
P0, P1-partial, P2, P3, P4-partial; Track M: M0, M1), designed as *evidence that the phase's
functionality works, is numerically stable, and physically sensible* — not
merely that tests pass. Each demo is a standalone script with built-in
acceptance checks against the roadmap gate criteria; its figures and CSVs are
committed under `cases/demo/<phase>/results/` so this report's numbers and
images always correspond to a reproducible run.

**Reproduce.** `python cases/demo/<phase>/run_demo.py` (headless, matplotlib
Agg). Exit code 0 = all checks pass; `results/checks.csv` holds the measured
value, criterion, and PASS/FAIL/XFAIL status per check. Status of this report:
generated 2026-07-07 from fresh runs; full pytest suite 117 passed + 2 xfailed
(~96 s since P3's medium-mesh compressible gate).

**Honesty rule.** P1 is *not* fully closed: gate G1.6 (sphere Cp < 2%) is
open, held as a strict xfail, and shown here as an XFAIL with its root cause —
the demos document the negative results (G1.3/G1.4 oracles, DP1) as evidence,
not as gaps to hide. Likewise P4 is *not* closed: the G4.1 medium-mesh gate
diverged on its first actual run (2026-07-07) — the coarse evidence and the
failure are both documented in §P4.

| Phase | Demo | Checks | Verdict |
|---|---|---|---|
| P0 mesh infrastructure | `cases/demo/p0_infrastructure/` | 4 PASS | closed, reproduced |
| P1 Laplace solver | `cases/demo/p1_laplace/` | 9 PASS + 1 XFAIL (G1.6) | closed gates reproduced; G1.6 open by design |
| P2 wake cut + Kutta | `cases/demo/p2_kutta_lifting/` | 11 PASS | closed, reproduced |
| M0 quasi-2D meshing | `cases/demo/m0_meshgen/` | 6 PASS | closed, reproduced |
| P3 subsonic compressible | `cases/demo/p3_subsonic/` | 14 PASS | closed, reproduced |
| P4 transonic artificial density | `cases/demo/p4_transonic/` | 10 PASS | **OPEN** — coarse evidence reproduced; G4.1 medium gate DIVERGED when first run (see §P4) |
| M1 swept-wing meshing (ONERA M6) | `cases/demo/m1_wing_mesh/` | 13 PASS | closed, reproduced |

---

## P0 — mesh infrastructure (G0.1–G0.4, closed)

**Purpose.** Establish that the geometric kernels every later phase builds on
are *exact*, so any later error is attributable to the PDE discretization,
never to geometry, I/O, or parallel-assembly bookkeeping.

**Case setup.** Analytic meshes (5-tet unit cube, 8-tet octahedron, Kuhn
cubes, scaled cube) plus the committed gmsh families.

**Key figures.**

![G0.1 volume exactness](../cases/demo/p0_infrastructure/results/g01_volume_exactness.png)
![G0.2 gradient exactness](../cases/demo/p0_infrastructure/results/g02_gradient_exactness.png)
![G0.3 coloring](../cases/demo/p0_infrastructure/results/g03_coloring.png)
![G0.4 VTK roundtrip](../cases/demo/p0_infrastructure/results/g04_vtk_roundtrip.png)

**Measured results.**

| Gate | Check | Measured | Criterion |
|---|---|---|---|
| G0.1 | max relative volume error, 5 analytic cases | 1.1e-16 | < 1e-12 |
| G0.2 | max linear-field gradient error (incl. random sliver tets) | 1.5e-13 | < 1e-12 |
| G0.3 | greedy coloring valid on 4 mesh families (24–45 colors) | all valid | validate_coloring |
| G0.4 | max write→read difference, 17k-tet gmsh mesh + fields | 0 (bit-exact) | < 1e-15 |

**Conclusion & analysis.** The P1-tet metric pipeline (volumes, Jacobians,
basis gradients) is exact to machine precision even on deliberately bad
(random, slivered) tets, the element coloring that licenses race-free
`prange` assembly is verified valid on every committed mesh family with
balanced class sizes, and VTK I/O is lossless. This is the foundation claim
of P0: downstream accuracy discussions can exclude geometry/infrastructure
as an error source.

---

## P1 — Laplace solver (G1.1, G1.2 closed; G1.3/G1.4 negative oracles; G1.6 OPEN)

**Purpose.** Show the three closed claims (consistency, order of accuracy,
solver scalability), show the flow physics is right, and document the open
G1.6 gap *with its evidence trail* so nobody re-litigates dead ends.

**Case setup.** MMS on structured cubes (φ = sin πx sin πy sin πz);
incompressible flow past the unit sphere on the committed `sphere_shell`
gmsh meshes (analytic φ and Cp = 1 − (9/4)sin²θ available in closed form).

**Key figures.**

![V0 freestream](../cases/demo/p1_laplace/results/v0_freestream.png)
![G1.1 MMS convergence](../cases/demo/p1_laplace/results/g11_mms_convergence.png)
![G1.2 AMG-CG scaling](../cases/demo/p1_laplace/results/g12_amg_cg_scaling.png)
![sphere flow field](../cases/demo/p1_laplace/results/sphere_flowfield.png)
![G1.6 open gate](../cases/demo/p1_laplace/results/g16_sphere_cp_open_gate.png)
![G1.4 oracle negative result](../cases/demo/p1_laplace/results/g14_oracle_negative_result.png)

**Measured results.**

| Gate | Check | Measured | Criterion | Status |
|---|---|---|---|---|
| V0 | max interior residual, φ = x, 4 mesh types (largest: 62k-tet NACA) | 8.8e-14 | < 1e-12 | PASS |
| G1.1 | MMS L2 slope over 4 levels (n = 4…32) | 1.96 | ≥ 1.9 | PASS |
| G1.2 | CG+AMG iterations 8 → 11 → 14 over 64× nodes (growth) | 1.75× | < 2× | PASS |
| V2-sanity | sphere solve residual (medium, 95k tets) | 5.7e-10 | < 1e-8 | PASS |
| V2-sanity | wall speed at stagnation poles / equator (exact 0 / 1.5) | 0.20 / 1.45 | sane | PASS |
| **G1.6** | **max wall Cp error, medium sphere** | **11.6%** | **< 2%** | **XFAIL (open)** |
| G1.4 | Option A t-form correction RHS magnitude (medium) | 6.8e-5 | near-zero data defect | PASS (negative result) |
| G1.4 | best exact-gradient correction moves max Cp err by | 0.23 pp | ineffective | PASS (negative result) |
| DP1 | best oracle-corrected max Cp err | 11.3% | > 5% branch confirmed | PASS |

**Conclusion & analysis.** The solver is *consistent* (freestream to machine
zero on every mesh type, including quasi-2D symmetry rows), *second-order
accurate* (MMS slope 1.96), and *scalable* (AMG-CG iteration count nearly
flat over a 64× node increase — the linear algebra will not be the bottleneck
at the 1–3 M-node target). The sphere flow field is physically correct:
stagnation at the poles, suction band at the equator (measured wall speed
1.45 vs exact 1.5), fore-aft symmetry.

The G1.6 gap (11.6% max wall Cp vs 2% gate) is *root-caused, not mysterious*:
the natural BC is satisfied on flat polyhedral facets instead of the true
curved sphere (variational crime). The oracle experiment reproduces the
recorded ceiling exactly (uncorrected 0.1156 → t-form 0.1164 / full-flux
0.1133): even feeding the *exact analytic gradient* into a boundary-data
correction moves the error by only 0.23 percentage points, because on
body-fitted meshes there is almost no boundary-data defect to fix (t-form
RHS max ~7e-5 on the sphere; exactly machine zero, ~1e-17, on the G1.3
cylinder). Hence DP1's "> 5%" branch: boundary-data corrections (Option A),
h-refinement, recovery tweaks, and Nitsche are ruled out with evidence; the
sanctioned route is the Option C gate re-spec plus separately-scoped
curved/isoparametric wall elements.

---

## P2 — wake cut, circulation, Kutta (G2.1–G2.5, closed)

**Purpose.** Show that lift is produced by machinery that is exact where it
must be exact (the cut), convergent where it iterates (the Kutta loop), and
physically consistent where it can be cross-checked (three independent lift
routes).

**Case setup.** NACA0012 at α = 4°, incompressible, on the M0 quasi-2D
meshes (16k / 62k tets), wake cut + master–slave elimination + secant Kutta
loop; reference = Hess–Smith panel solution
(`cases/reference_data/naca0012_incompressible/`, cl(4°) = 0.482556).

**Key figures.**

![G2.1/G2.2 cut exactness](../cases/demo/p2_kutta_lifting/results/g21_g22_cut_exactness.png)
![G2.3 Kutta convergence](../cases/demo/p2_kutta_lifting/results/g23_kutta_convergence.png)
![G2.3 Cp vs panel](../cases/demo/p2_kutta_lifting/results/g23_cp_vs_panel.png)
![G2.4 lift cross-check](../cases/demo/p2_kutta_lifting/results/g24_cl_crosscheck.png)
![lifting flow field](../cases/demo/p2_kutta_lifting/results/lifting_flowfield.png)
![G2.5 spanwise decay](../cases/demo/p2_kutta_lifting/results/g25_spanwise_decay.png)

**Measured results.**

| Gate | Check | Measured | Criterion |
|---|---|---|---|
| G2.1 | max free-dof residual, φ = x on cut mesh (coarse/medium) | 8.4e-13 / 5.6e-13 | < 1e-12 |
| G2.1 | max wake-master row residual | 6.9e-16 / 5.3e-16 | < 1e-13 |
| G2.2 | max \|[φ] − Γ\| for prescribed Γ = 0.3, every wake pair | 7.2e-16 | < 1e-12 |
| G2.3 | Kutta updates to convergence (secant) | 2 | < 20 |
| G2.3 | cl = 0.47858 vs panel 0.48256 (medium; coarse −3.0%) | −0.82% | \|err\| < 2% |
| G2.4 | cl from Γ (0.478646) vs cl from pressure (0.478576) | 0.015% | < 1% |
| — | max speed one cell off the TE (no vortex-like suction spike) | 0.96 U∞ | Kutta physically enforced |
| G2.5a | spanwise gradient of interpolated φ = x on cut mesh | 8.9e-15 | < 1e-12 |
| G2.5b | p99 \|w\|/U∞ 4.82e-3 → 2.35e-3 (RMS 1.62e-3 → 8.2e-4) | ratio 2.05 | ≥ 1.8 (1st order) |

**Conclusion & analysis.** The wake-cut machinery adds nothing spurious: with
Γ = 0 the cut mesh preserves freestream to the same machine zero as an uncut
mesh, and a prescribed jump is reproduced to 1e-16 at every wake node pair.
The secant-accelerated Kutta loop converges in 2 updates — essential, since
the measured relaxation-map slope b ≈ 0.93 means plain under-relaxation would
need O(100) outer solves. Physically, lift is over-determined and consistent:
surface-pressure integration, Kutta–Joukowski from the circulation the loop
found, and the independent panel reference agree (Γ-vs-pressure 0.015%,
vs panel −0.82% on the medium mesh, converging from −3.0% on coarse) — the
circulation is not a tuning knob, it *is* the lift the pressure field
carries. The mid-span field shows smooth flow off the trailing edge and the
branch cut carrying the constant jump Γ = 0.239 in the perturbation
potential. Quasi-2D consistency holds at the level the discretization
permits: the 3-tet prism split makes exact spanwise invariance of a *solved*
field impossible (the roadmap's re-specced G2.5b), and the demo reproduces
the honest criterion — field-wide spanwise noise decays at clean first order
(ratio 2.05 at h-ratio 2), with the un-gated single-element max confined to
the leading-edge gradient peak, not the wake. This also documents why **TE
nodes are duplicated**: a single-valued TE would taper [φ] to zero across
the first wake cell — a point vortex parked at the TE whose spurious suction
~Γ²/h *diverges* under refinement (roadmap P2 assert block).

---

## P3 — subsonic compressible (G3.1–G3.3, closed)

**Purpose.** Show that the density machinery (isentropic ρ(q²) in a nested
Picard outer loop) adds *compressibility and nothing else*: it reproduces the
classical corrections where they apply, degenerates bit-exactly to the
validated Laplace solvers at M∞ → 0, and its assembly rewrite (the retired P1
tech debt) is a pure performance change, not a numerics change.

**Case setup.** (a) Sphere at M∞ = 0.3 (medium shell, 95k tets), freestream
Dirichlet, compressible vs incompressible on the *same mesh with the same
quadratic surface recovery* — the comparison G3.1 prescribes, which cancels
the known G1.6 flat-facet wall bias. (b) NACA0012 at M∞ = 0.5, α = 2° on the
M0 quasi-2D meshes, nested Picard (outer density update, inner secant Kutta
at frozen ρ), PG-scaled vortex far field; reference =
`cases/reference_data/naca0012_m05/` (the P2-verified panel solution under
Prandtl–Glauert and Kármán–Tsien corrections). (c) M∞ = 0 bit-identity
against `solve_laplace_lifting`.

**Key figures.**

![G3.1 sphere Cp + convergence](../cases/demo/p3_subsonic/results/g31_sphere_cp_and_convergence.png)
![G3.2 Cp + Picard convergence](../cases/demo/p3_subsonic/results/g32_naca_cp_and_convergence.png)
![G3.2 cl bracket](../cases/demo/p3_subsonic/results/g32_cl_bracket.png)

**Measured results.**

| Gate | Check | Measured | Criterion |
|---|---|---|---|
| P3-debt | colored fast path vs P1 serial reference kernels | 5.7e-16 rel | < 1e-13 |
| P3-debt | reassembly determinism | bitwise equal | bitwise |
| P3-debt | hot reassembly speedup (medium NACA, warm JIT) | ~160× | > 2× |
| G3.1 | sphere Cp peak −1.33591 vs PG-corrected −1.34020 | 0.32% | < 2% |
| G3.1 | non-lifting Picard (density lag < 1e-10) | 11 iters, monotone | converges |
| G3.2 | cl = 0.28437 vs [PG 0.27877, KT 0.29186] | inside bracket | PG ≤ cl ≤ KT |
| G3.2 | cl vs PG/KT midpoint 0.28532 | −0.33% | < 2% |
| G3.2 | density (outer) iterations | 15 | < 30 |
| G3.2 | residual history 5.9e-6 → 5.0e-11 | strictly monotone | monotone |
| G3.2 | max local Mach | 0.726 | < 1 (ν ≡ 0 regime) |
| G3.3 | A(ρ(M=0)), Kutta φ and Γ vs P2 driver | bitwise equal | bitwise |
| G3.3 | P1/P2 gates with ρ machinery in the loop | 117 passed + 2 xfailed | all green |

**Conclusion & analysis.** The compressibility machinery lands where the
classical corrections say it must: the sphere's suction-peak amplification
matches Prandtl–Glauert to 0.32% once the mesh's own (known, G1.6) wall bias
is cancelled by same-mesh comparison, and the airfoil cl falls inside the
[PG, KT] correction bracket, −0.33% from its midpoint — with the mesh-family
trend pointing toward the KT side in the continuum limit, as expected for a
12%-thick section. Convergence behaves like the theory: the density lag
contracts geometrically (rate ≈ 0.3) and the residual decays strictly
monotonically to the linear-solver floor. Two negative results are recorded
as design constraints: interleaving one Γ update per density iteration (the
literal design.md §8 pseudocode) injects 10× residual spikes down the whole
history — nesting the P2 secant Kutta loop at frozen ρ (where the Kutta map
is exactly affine) restores monotonicity; and a loose *relative* inner CG
tolerance combined with warm starts false-converges (CG exits at x0 without
computing the density correction), so the shipped inexact option is a
forcing term `atol = η‖b − A x₀‖` (η = 0.05 ≈ 2× faster; default η = 0).
The M∞ → 0 limit is bit-exact, not approximately-equal: ρ ≡ 1.0 bitwise from
the density law, β = 1 reduces the PG vortex bit-exactly, the P1/P2 drivers
share the rewritten assembly, and pyamg's unseeded spectral-radius RNG — the
one source of run-to-run scatter, measured at 2e-11 between *identical*
solves — is now pinned in `solve/linear.py::build_amg_preconditioner`.

---

## P4 — transonic artificial density (G4.2/G4.3 closed; G4.1 OPEN — medium gate diverged)

**Purpose.** Show that the artificial-density upwinding produces sharp,
monotone, correctly-placed shocks; that it is an *exact* no-op below
critical Mach; and record the scheme-hardening evidence trail that P4's
"main risk" phase actually consumed (four instability mechanisms found
and fixed with measurements).

**Case setup.** NACA0012 quasi-2D family, M∞ = 0.80, α = 1.25° — the
canonical transonic benchmark (strong upper shock, weak lower shock).
Pipeline: Mach continuation (0.70 → 0.75 → 0.80), each supercritical
level closing the Kutta condition by an outer per-station secant around
frozen-Γ pseudo-time density solves. Reference:
`cases/reference_data/naca0012_m080/` (Euler anchor ~0.60c/0.35c +
documented conservative-FP aft-shift band; provenance in its README —
an open digitized FP table for this case was not retrievable, and the
README says so rather than inventing one).

**Key figures.**

![upwind reach evidence](../cases/demo/p4_transonic/results/p4_upwind_reach.png)
![G4.1 Cp and shock](../cases/demo/p4_transonic/results/g41_cp_shock.png)

**Measured results (coarse evidence run; the G4.3 coarse sweep is in
artifacts/G4.3; the G4.1 medium run is in artifacts/G4.1 and FAILED —
see the addendum below).**

| Gate | Check | Measured | Criterion |
|---|---|---|---|
| G4.2 | ν at M∞ = 0.5 with machinery active | max ν = 0.0 exactly | ν ≡ 0 |
| G4.2 | φ/Γ vs the P3 code path (upwind_c = 0) | bitwise identical | bit-identical |
| scheme | single-hop upstream reach (M0.70 pocket) | median 0.37 extents | documented root cause |
| scheme | multi-hop walk reach | median 1.00 extents | ~1 streamwise cell |
| G4.1 | upper shock x/c | 0.599 | 0.62 ± 0.03 (ref band) |
| G4.1 | shock monotone / expansion shock | monotone, none | required |
| G4.1 | shock sharpness | 1 deduped station (2–3 raw cells) | 2–3 cells |
| G4.1 | lower weak shock x/c | 0.362 | ~0.35 (reported) |
| G4.1 | M_max | 1.363 | physical, no limited cells |
| G4.1 | Γ closure | secant \|F\| = 9.3e-5 in 8 evals | < 2e-4 |

**Conclusion & analysis.** The shipped scheme is the design.md §3
artificial density plus four evidence-forced hardenings, each of which
was found by a measured failure: (1) a **multi-hop upstream walk** —
on the prism-split meshes a single face-neighbor hop reaches only ~1/3
of an element's streamwise extent, putting the effective ν below the
(M²−1)/M² stability bound for M ≳ 1.2 (measured blow-up at M0.75 for
every ω ∈ 0.3–0.9, C ∈ 1.5–2.0, while marginally-stable M0.70 crawled);
(2) a **speed limiter** q² ≤ q²(M=3) — without it, transients run to the
vacuum limit and the positivity guards then stabilize *spurious* dead-
cell solutions (measured: off-body supersonic blobs acting as fake
blockage); (3) a **shock-point operator** ν = max(ν_e, ν_upstream) — the
first subsonic cell downstream of the shock is otherwise purely central
on the field's largest jump; and (4) a **pseudo-transient term**
diag(m/Δτ) (design.md §8 acceleration 4 pulled forward) — the exact-
solve Picard limit-cycles at M0.80 shock strength for every relaxation
tried (φ, ρ̃, or both), while Δτ ≈ 1e-3–3e-3 bounds the update and
yields the physical field (M_max 1.36–1.46); Δτ-ramping (SER) re-ignites
the instability, so Δτ stays fixed and the sharp-shock residual settles
into a slowly-decaying engineering tail (cl drift < 1e-3 per hundreds of
iterations) instead of 1e-10 — the documented P4 convergence semantics,
with Newton (P6) as the designed cure. Γ closure had to move OUTSIDE the
density loop entirely: nested exact-Kutta runs away (Γ 0.115 → 4.99) and
damped interleaving limit-cycles, because the transonic target map's
slope crosses 1 where relaxed fixed-point updates provably diverge — the
outer secant on density-converged evaluations converges in ~4–8 warm-
started evaluations. The coarse result lands where the references say it
should: upper shock 0.599 (Euler anchor 0.60–0.63 band), weak lower
shock 0.362 (~0.35), monotone with no expansion shock, and the
subcritical bit-identity G4.2 guarantees P3 behavior is untouched.

**Addendum (2026-07-07, audit): the G4.1 MEDIUM gate FAILED on its first
actual run — P4 is NOT closed.** The phase had been declared closed with
the medium result left as an unfilled placeholder; running the gate
(`PYFP3D_TRANSONIC_GATES=1`, 2h43m wall) produced a diverged solve:
M_max 30.1 (vs the physical 1.36 on coarse), 423 q²-limited + 271
density-floored cells, Kutta |F| 2.5e-3 vs tol 2e-4, a spurious "shock"
at x/c 0.802, and sign-inconsistent cl (pressure −0.171 vs KJ +0.212).
Both supercritical continuation levels exhausted their 12 Γ evaluations
without secant convergence (n_picard_total 19331; each frozen-Γ eval
burns its full 800-iteration budget by design — the tol_rho early-exit
is unreachable at transonic). Evidence:
`artifacts/G4.1/summary_medium.csv`, `medium_gate_pytest.log`,
`v4_1_cp_shock_medium.png`. Interpretation: the coarse-calibrated
`TRANSONIC_DEFAULTS` (fixed Δτ = 2e-3, dm = 0.05) do not transfer to the
medium mesh — the pseudo-transient damping diag(m_lumped/Δτ) scales as
h³ while the operator stiffness does not, so the stabilizer weakens
under refinement (hypothesis, not yet verified). Candidate routes:
mesh-scaled Δτ, finer Mach continuation steps, or P6 Newton. The gate is
re-opened in roadmap.md; the coarse evidence above and G4.2/G4.3 stand.

---

## M0 — quasi-2D meshing pipeline (closed; acceptance link = G2.5)

**Purpose.** Show the mesh-side evidence for M0: the pipeline (vanilla-Gmsh
planar mesh → single-layer extrusion → globally consistent min-global-index
prism→3-tet split) produces topologically sound, refinable meshes that the
solver actually converges on. M0's formal acceptance was G2.5 — the Track M
↔ Track P link demonstrated in the P2 demo above.

**Case setup.** All seven committed meshes in `cases/meshes/`
(naca0012_2.5d coarse/medium, cylinder_2.5d coarse/medium/fine,
sphere_shell coarse/medium); end-to-end solve on the cylinder family
against the analytic Cp = 1 − 4 sin²θ.

**Key figures.**

![mesh gallery](../cases/demo/m0_meshgen/results/mesh_gallery.png)
![topology asserts](../cases/demo/m0_meshgen/results/topology_asserts.png)
![mesh quality](../cases/demo/m0_meshgen/results/mesh_quality.png)
![cylinder Cp convergence](../cases/demo/m0_meshgen/results/cylinder_cp_convergence.png)

**Measured results.**

| Check | Measured | Criterion |
|---|---|---|
| topology asserts (tags, quad-split consistency, wake cut) on all 7 meshes | all pass | hard rule 7 |
| min tet volume across families | 1.3e-6 > 0 | no degenerate/inverted tets |
| isotropic sphere family max aspect ratio | 3.5 | < 5 (quasi-2D far-field anisotropy is by design) |
| cylinder solve residual, all 3 levels | ≤ 3.1e-11 | < 1e-8 |
| max wall Cp error coarse → medium → fine | 9.1% → 4.5% → 2.2% | monotone, ≥ 25%/level |
| Cp error slope vs h | 1.02 | ~O(h), documented curved-wall limit |

**Conclusion & analysis.** Every committed mesh passes the full topology
assert battery (agent-rules hard rule 7), including the wake-cut asserts on
the NACA family — the wake is a single conforming interior sheet from TE to
far field, and the prism split is globally consistent so no lateral quad is
cracked. Element quality is deliberate rather than accidental: the isotropic
sphere family stays under aspect 3.5, while the quasi-2D families carry
their large aspect ratios exactly where the single-layer design puts them
(far field: in-plane coarsening at fixed dz). The pipeline's meshes are not
just valid but *useful*: the cylinder case solves to 1e-11 residuals at
every level and its wall-Cp error falls monotonically at the expected O(h)
(slope 1.02) — first order because of the same flat-facet curved-wall
recovery limit root-caused at G1.6, i.e. a documented solver-side limit, not
a meshing defect. Combined with G2.3–G2.5 running on these meshes (P2 demo),
M0's deliverable is demonstrated end to end.

---

## M1 — swept/tapered wing meshing, ONERA M6 (closed; consumed by P5)

**Purpose.** Show the mesh-side evidence for M1: a scripted, refinable
ONERA M6 half-wing tet mesh whose chord-plane wake sheet — swept from the
sharp (foilmod zero-thickness) TE, ending exactly at the tip, reaching the
spherical far field at 15 MAC — is ingested by the P2 solver preprocessor
with the topology asserts green. The new mesh-side machinery is
`pyfp3d/meshgen/wing3d.py` (OCC ruled loft + `occ.fragment` +
`mesh.embed`); the new solver-side machinery is wake_cut.py's handling of
a swept TE (per-node stations, off-plane Kutta-probe fallback) and of the
sheet's interior FREE edge at the tip (single-valued nodes ⇒ Γ(tip) = 0
discretely).

**Case setup.** The `cases/meshes/onera_m6` family (coarse 55.5k /
medium 350.7k tets; fine 2.513M validated at generation time). The .msh
files are large and gitignored — regenerate coarse+medium with
`generate_onera_m6.py` (~30 s) before running this demo; the committed
per-level stats CSVs and inspection PNGs are the persistent evidence.
Solver axis convention: chord x, lift y, span z.

**Key figures.**

![wing + wake gallery](../cases/demo/m1_wing_mesh/results/wing_wake_gallery.png)
![tip cut planes](../cases/demo/m1_wing_mesh/results/tip_cut_planes.png)
![mesh quality](../cases/demo/m1_wing_mesh/results/mesh_quality.png)

**Measured results.**

| Check | Measured | Criterion |
|---|---|---|
| tags + P2 topology asserts through cut_wake, coarse & medium | pass | M1 gate "same asserts" |
| per-node TE stations on the swept TE | 83 (coarse) / 166 (medium) | == n_TE_nodes |
| tip free-edge nodes single-valued, at z = b | 106 / 208, none duplicated | wake-tip semantics |
| wake-tip closure: tip edge one open chain from the exact tip TE corner | pass (both levels) | no cracks / self-intersections |
| Kutta probe pairs found on the unstructured TE | 83 / 166, y>0 upper, y<0 lower | design.md (4.4) fallback |
| min dihedral coarse/medium/fine | 7.5° / 11.0° / 3.5° | ≥ 2° |
| max aspect ratio coarse/medium/fine | 9.3 / 6.9 / 6.5 | ≤ 60 |
| refinement ladder (one h_wall parameter, 2×) | 55.5k → 350.7k → 2513k tets | monotone ~2³/level |
| freestream residual on the CUT coarse mesh | 4.3e-14 | < 1e-10 (G2.1 analogue) |

**Conclusion & analysis.** The M1 gate items are all measured green: the
solver preprocessor ingests the family unchanged (same read_mesh/cut_wake
call as the 2.5D cases), the quality report is comfortably inside bounds
on all three levels, and the family is one script with one parameter. Two
findings worth recording: (1) for a sheet that ends *inside* the domain,
`occ.fragment` alone does not make the tet mesh conform — it stitches the
shared TE edge and the boundary trims, but `gmsh.model.mesh.embed` must
still be called on the trimmed sheet face; (2) the sheet's tip edge is an
interior free edge whose node stars are NOT split by the sheet, so the
duplication map must exclude them — which is also the physically correct
discrete statement Γ(tip) = 0 (the trailing jump vanishes at the tip).
Both are documented in the wake_cut.py module docstring; the free-edge
path is exactly inert on the quasi-2D meshes (their sheets have no free
edges), which the unchanged P2/M0 test battery confirms. Mesh sizes are
runtime-driven per the P4 lesson (solver wall time is the binding
constraint): coarse is the P5 development mesh, medium the gate mesh.

---

## Cross-phase summary

- **Functionality**: every closed gate's headline number is reproduced from
  scratch by the demos (MMS slope 1.96, CG 8→11→14, cl −0.82%, Γ-lift
  cross-check 0.015%, spanwise ratio 2.05, cylinder slope 1.02).
- **Numerical stability**: machine-zero consistency checks (V0, G2.1, G2.2,
  G2.5a) hold on the largest committed meshes, the linear solver is
  mesh-independent, and the Kutta outer loop converges in 2 updates.
- **Physical soundness**: stagnation/suction structure on the sphere,
  smooth TE flow, and three mutually independent lift routes agreeing —
  physics cross-checks, not code self-consistency.
- **Open item**: G1.6 remains open (11.6% vs 2%); its Option C acceptance
  re-spec is the open P1 task, with Option A conclusively ruled out by the
  oracle demos.
- **P3 additions**: compressibility validated against the classical
  correction band (sphere 0.32% vs PG; airfoil cl inside [PG, KT]);
  M∞ → 0 is bit-identical to the P1/P2 Laplace drivers; assembly is
  colored-`prange` with precomputed geometry (~160× hot reassembly) and
  every solve is bit-reproducible run-to-run (seeded AMG setup).
