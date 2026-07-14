# Phase demo report — evidence for completed phases

**Scope.** One self-contained demo case per completed roadmap phase (Track P:
P0, P1-partial, P2, P3, P4; Track M: M0, M1), designed as *evidence that the phase's
functionality works, is numerically stable, and physically sensible* — not
merely that tests pass. Each demo is a standalone script with built-in
acceptance checks against the roadmap gate criteria; its figures and CSVs are
committed under `cases/demo/<phase>/results/` so this report's numbers and
images always correspond to a reproducible run.

**Reproduce.** `python cases/demo/<phase>/run_demo.py` (headless, matplotlib
Agg). Exit code 0 = all checks pass; `results/checks.csv` holds the measured
value, criterion, and PASS/FAIL/XFAIL status per check. Status of this report:
updated 2026-07-07 (P4 re-closed same day, see §P4 addenda) from fresh runs;
full pytest suite 136 passed + 2 skipped + 2 xfailed (~5 min; the always-on
coarse transonic smoke is ~170 s of it).

**Honesty rule.** P1 is *not* fully closed: gate G1.6 (sphere Cp < 2%) is
open, held as a strict xfail, and shown here as an XFAIL with its root cause —
the demos document the negative results (G1.3/G1.4 oracles, DP1) as evidence,
not as gaps to hide. P4 was briefly re-opened the same day it was first
declared closed (its medium-mesh gate had never actually been run) — the
divergence, root cause, and fix are all documented in §P4's addenda rather
than silently corrected away.

| Phase | Demo | Checks | Verdict |
|---|---|---|---|
| P0 mesh infrastructure | `cases/demo/p0_infrastructure/` | 4 PASS | closed, reproduced |
| P1 Laplace solver | `cases/demo/p1_laplace/` | 9 PASS + 1 XFAIL (G1.6) | closed gates reproduced; G1.6 open by design |
| P2 wake cut + Kutta | `cases/demo/p2_kutta_lifting/` | 11 PASS | closed, reproduced |
| M0 quasi-2D meshing | `cases/demo/m0_meshgen/` | 6 PASS | closed, reproduced |
| P3 subsonic compressible | `cases/demo/p3_subsonic/` | 14 PASS | closed, reproduced |
| P4 transonic artificial density | `cases/demo/p4_transonic/` | 10 PASS | closed, reproduced (re-closed 2026-07-07 — see Addendum 3) |
| M1 swept-wing meshing (ONERA M6) | `cases/demo/m1_wing_mesh/` | 13 PASS | closed, reproduced |
| P5 3D validation (ONERA M6) | `cases/demo/p5_onera_m6/` | 16 PASS | closed 2026-07-08 (V6 < 1% deferred to P9) |
| P6 surface-pressure recovery | `cases/demo/p6_surface_recovery/` | 6 PASS (incl. gated M6) | closed 2026-07-08 (sawtooth = recovery artifact) |
| P7 differentiable walk flux | `cases/demo/p7_diff_flux/` | 7 PASS (incl. gated converged-field) | closed 2026-07-10 (FD 3–5e-10) |
| P8 fully-coupled Newton | `cases/demo/p8_newton/` | 15 PASS (parts 2–3 gated) | closed 2026-07-11 (G8.1 + G8.2 + G8.3) |
| P8 capability assessment | `cases/demo/p8_capability/` | 36 PASS (full matrix gated) | **evaluation demo, not a gate** (2026-07-11) |
| P10 (partial) G10.2 continuation tolerance | `cases/demo/p10_newton_usability/` | split A/B verdict | G10.2 + G10.3 closed 2026-07-11; phase stays open (G10.1) |
| P9 grid-convergence & accuracy-gap discrimination | `cases/demo/p9_grid_discrimination/` | 11 PASS + 3 XFAIL | closed 2026-07-11 (G9.3 verdict awaits arbitration) |
| **Track B** B1 cut-element identification | `tests/test_b1_cut_elements.py` (test-only) | 34 PASS | closed 2026-07-11 |
| **Track B** B2 multivalued assembly | `tests/test_b2_multivalued.py` (test-only) | 17 PASS | closed 2026-07-11 |
| **Track B** B3 + B4 lifting + TE Kutta | `cases/demo/b3_levelset_lifting/` | 13 demo PASS (+6, +8 tests) | closed 2026-07-12 |
| **Track B** B5 far-field A/B | `cases/demo/b4p5_farfield/` | 9 demo PASS (+10 tests) | closed 2026-07-12 |
| **Track B** B6 transonic (level-set) + LS Newton | `cases/demo/b6_transonic/` | 14 demo PASS (+9, +2 tests; +2, +2 Newton) | ◐ IN PROGRESS 2026-07-12 — coarse M0.80 gate met, LS Newton reaches quadratic fold solutions, medium closure open |

> Track-P renumber (2026-07-08, then 2026-07-11 ×2): P6 = surface recovery;
> P7 = differentiable flux (Newton prereq); P8 = fully-coupled Newton;
> P9 = grid-convergence discrimination; P10 = Newton generality/continuation;
> P11 = curved wall elements; P12 = backlog.
> **Track-B renumber (2026-07-12 ×2):** a new **B4** (TE control volume) was
> inserted, then the half-integer IDs were regularized away — the far-field gate
> is now **B5** (was B3.5, then B4.5; its demo dir keeps the old `b4p5_` name),
> transonic is **B6**, ONERA M6 3D is **B7**. See roadmap.md for the full mapping.

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

## P4 — transonic artificial density (closed 2026-07-07: G4.1/G4.2/G4.3 all green)

> **★ ERRATUM (2026-07-11, P8/N5 finding — user-approved: record, do not
> re-open).** The "converged" states this section reports are Picard STALL
> states, not solutions of the discrete equations: the P8 exact-Newton
> residual at the committed coarse M0.80 state is **2.2e-4**, and Newton
> walks from it to the true discrete solution (coarse: shock 0.658,
> cl 0.459, M_max 1.408; on medium the family steepens into the FP
> non-uniqueness fold with no reachable solution at M0.80). The gates below
> stand as **Picard-quality/robustness gates** (the machinery is the
> production warm-start engine); Newton-era physical acceptance lives in
> the G8.1 regression locks. Full evidence: roadmap P4 ledger erratum +
> demo_report §P8.

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

**Measured results (coarse evidence run, re-measured 2026-07-07 under the
`damping_theta` fix — see Addendum 3; the committed coarse/medium and
G4.3-sweep figures + CSVs live in `cases/demo/p4_transonic/results/`,
regenerated by `run_demo.py` heavy mode — see the supplementary analysis
subsection; the G4.1 medium run now PASSES).**

| Gate | Check | Measured | Criterion |
|---|---|---|---|
| G4.2 | ν at M∞ = 0.5 with machinery active | max ν = 0.0 exactly | ν ≡ 0 |
| G4.2 | φ/Γ vs the P3 code path (upwind_c = 0) | bitwise identical | bit-identical |
| scheme | single-hop upstream reach (M0.70 pocket) | median 0.37 extents | documented root cause |
| scheme | multi-hop walk reach | median 1.00 extents | ~1 streamwise cell |
| G4.1 coarse | upper shock x/c | 0.604 | 0.62 ± 0.03 (ref band) |
| G4.1 coarse | shock monotone / expansion shock | monotone, none | required |
| G4.1 coarse | shock sharpness | 1 deduped station (2–3 raw cells) | 2–3 cells |
| G4.1 coarse | lower weak shock x/c | 0.358 | ~0.35 (reported) |
| G4.1 coarse | M_max | 1.373 | physical, no limited cells |
| G4.1 coarse | Γ closure | secant \|F\| = 1.50e-4 in 9 evals | < 2e-4 |
| G4.1 medium | upper shock x/c | 0.633 | 0.62 ± 0.03 (ref band) |
| G4.1 medium | cl_pressure / cl_KJ | 0.349 / 0.354 (sign-consistent) | physical, consistent |
| G4.1 medium | M_max | 1.366 | physical, no limited cells |
| G4.1 medium | Γ closure | secant \|F\| = 1.23e-4, 16m39s wall | < 2e-4 (was: diverged, 2h43m) |

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
with Newton (P7) as the designed cure. Γ closure had to move OUTSIDE the
density loop entirely: nested exact-Kutta runs away (Γ 0.115 → 4.99) and
damped interleaving limit-cycles, because the transonic target map's
slope crosses 1 where relaxed fixed-point updates provably diverge — the
outer secant on density-converged evaluations converges in ~4–9 warm-
started evaluations. The coarse result lands where the references say it
should: upper shock 0.604 (Euler anchor 0.60–0.63 band), weak lower
shock 0.358 (~0.35), monotone with no expansion shock, and the
subcritical bit-identity G4.2 guarantees P3 behavior is untouched. (These
are the post-fix numbers, re-measured under the `damping_theta` default —
see Addendum 3; they are consistent with the original 0.599/0.362 to
within the run-to-run Γ-secant noise documented in Addendum 2.)

**Addendum (2026-07-07, audit): the G4.1 MEDIUM gate FAILED on its first
actual run — P4 is NOT closed.** *[Superseded the same day by Addendum 3
below, which lands the fix and re-closes P4; retained verbatim as the
honest failure trail, not as current status.]* The phase had been declared closed with
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
mesh-scaled Δτ, finer Mach continuation steps, or P7 Newton. The gate is
re-opened in roadmap.md; the coarse evidence above and G4.2/G4.3 stand.

**Addendum 2 (2026-07-07, diagnosis follow-up — informal ad-hoc scripts,
not a committed gate artifact): root cause verified, one route ruled
out, three working fixes measured.** The pseudo-transient damping ratio
(m_lumped/Δτ)/diag(K) at the shipped Δτ = 2e-3 was measured directly
rather than only inferred from scaling arguments: wall-node median 0.035
on coarse falls to 0.0092 on medium (~4× weaker) — confirming the
addendum-1 hypothesis. A standalone frozen-Γ density driver, built to
isolate the density solve from the Γ secant entirely, reproduces the
M0.75 blow-up in 50 iterations (M_max → 47), so the divergence lives in
the pseudo-time-stabilized density iteration itself, not in the outer
Kutta root-find. Rescaling Δτ by h² (5e-4 on medium) restores M0.75
stability (M_max 1.22, zero limited/floored cells) — **but the same
state still diverges stepping to M0.80**: the damping needed also grows
with shock strength, so no single Δτ-vs-mesh-size rescaling suffices on
its own. A finer Mach-continuation sub-step (intermediate M0.775, still
at Δτ = 5e-4) **also diverges at M0.80 — the "finer dm" candidate route
is ruled out**; this is not a per-step transient-overshoot problem.
Three candidates were measured to stabilize the M0.80 step from a
converged M0.75 state (500 iterations each, zero limited/floored cells
throughout): global Δτ = 2e-4 (M_max 1.37); **local damping
θ·diag(A_free)** at θ = 0.2 (M_max 1.37, the recommended fix — mesh- and
shock-strength-independent by construction, unlike the mass-lumped
global form that caused the original failure); upwind_c raised 1.5→2.0
at Δτ = 5e-4 (M_max 1.34, more dissipative, would need a shock-position
re-validation against G4.3). None of the three alone closes the gate as
specced: Kutta mismatch |F| still drifts to 4e-4–9e-4 under a fixed
500–800-iteration budget, mirroring a drift already visible on
coarse — |F| rises monotonically from 2.9e-4 at iteration 100 to 8.9e-4
at iteration 800 of a *single* eval, i.e. the eval's own Kutta target
moves inside the fixed budget rather than converging to it and
stalling. This also explains why coarse itself is slower than the
design.md §8 O(100–300)-iteration expectation: the eval path solves
every inner CG to rtol = 1e-10 (`forcing = 0`, bypassing the
P3-validated η ≈ 0.05 forcing-term acceleration already shipped
elsewhere) and rebuilds the AMG hierarchy every 4 iterations; medium
profiling attributes 64% of eval wall-time to CG (22 CG-iterations per
outer step vs 3 on coarse — itself a symptom of the weak damping's poor
conditioning) and 27% of coarse eval wall-time to AMG rebuilds. No
equation-level bug was found in `picard.py`/`upwind.py`/`wake.py` during
this session; one piece of negative evidence supports that: an
independent from-scratch reimplementation of the frozen-Γ inner loop
diverged immediately until it included the same per-iteration
`h_j = TᵀA(ρ̃)g_j` recompute that `WakeConstraint.update_matrix` already
performs, confirming that recompute is load-bearing rather than
removable overhead. Recommended before the next medium-gate attempt:
switch to local θ·diag(A_free) damping (validate θ on coarse first —
G4.2 bit-identity and the G4.1 coarse shock position must survive),
enable eval-path forcing plus a wider AMG rebuild interval, and add an
adaptive per-eval exit on |F| drift (or re-match tol_gamma to the
measured drift floor) instead of the fixed-iteration budget. Separately
identified P5 blocker: `WakeConstraint.update_matrix`'s per-station
`h_j` loop is one sparse matvec per wake station — inert on the
single-station 2.5D meshes but adds ~166 extra matvecs per density
iteration on the ONERA M6 medium mesh (166 stations, M1 delivery); it
must be batched into one `A @ G` sparse product before P5 solves start.

**Addendum 3 (2026-07-07, same day as Addendum 2: fix landed, medium gate
closed, P4 re-closed).** The recommended fix from Addendum 2 was
implemented as-is and closed the gate on the first attempt, with none of
the secondary mitigations (eval-path forcing, wider AMG rebuild interval,
adaptive |F|-drift exit) needed. `solve/picard.py::solve_subsonic_lifting`
gained `damping_theta`: D = θ·diag(A_free), rebuilt every outer iteration
from that iteration's own (upwinded) operator, mutually exclusive with the
retired global `pseudo_dt` (both given raises `ValueError`).
`solve/continuation.py::TRANSONIC_DEFAULTS` now defaults to
`damping_theta = 0.2` in place of `pseudo_dt = 2e-3`. Calibration order
followed the Addendum 2 recommendation: G4.2 bitwise no-op re-verified
first (both new params default `None`, so any caller that doesn't pass
them — including the G4.2 test — is untouched), then the G4.1 coarse
gate re-measured (shock 0.599 → 0.604, within the Γ-secant noise; the
converged Γ/cl shifted a touch more, 0.170 → 0.182 / cl_p 0.334 → 0.357,
since the two stabilizers reach slightly different points on the same
bounded engineering-convergence tail — the shock position and all gated
criteria are unaffected, and cl is reported-not-gated for G4.1), then a
reduced-budget stability probe on the medium mesh (`max_gamma_evals=4,
n_picard_eval=150`, 77 s wall) to catch a repeat divergence cheaply before
committing to the full budget: it came back stable and, unexpectedly,
needed only 1 Γ eval per supercritical level (vs 12 exhausted before) —
suggesting the Kutta-mismatch drift Addendum 2 flagged as a likely second
blocker was largely a symptom of the poorly-conditioned global damping,
not an independent problem. The full medium gate then passed outright
(committed `cases/demo/p4_transonic/results/g41_summary_medium.csv`),
**16m39s wall** (vs the divergent
attempt's 2h43m) — upper shock x/c 0.633 (band 0.62 ± 0.03), Kutta
|F| = 1.23e-4 < tol 2e-4, M_max 1.366, zero limited/floored cells,
cl_pressure/cl_KJ sign-consistent at 0.349/0.354 (the divergent run had
−0.171/+0.212), n_picard_total 12931 (vs 19331 — fewer Γ evals needed,
not merely faster ones: per-iteration wall time also fell from ~0.51 s to
~0.18 s, consistent with the Addendum 2 profiling that blamed the old
global damping's poor conditioning for elevated CG iteration counts on
medium). The G4.3 ten-case sweep was re-run under the new default and
stays green — all converged, zero limited cells, smooth trends — with one
difference recorded as evidence rather than a regression: the
M0.82/α = 1.25° corner's cl moved 0.389 → 0.458 and its Kutta |F| (1.92e-4)
sits closest to the 2e-4 tolerance of the ten cases; G4.3 gates on
convergence and physicality, not an exact cl, so this does not fail it,
but a future θ-sensitivity study should know this corner is the closest
to the boundary. The separately-identified P5 blocker was fixed in the
same session: `constraints/wake.py::WakeConstraint.update_matrix` now
builds one sparse indicator matrix `G` (one column per station) and
computes all stations' `h_j` vectors as a single `T^T @ (A @ G)` product;
verified bit-identical to the old per-station loop on the real ONERA M6
coarse mesh (83 stations). The full default suite (136 passed + 2 skipped
+ 2 xfailed, ~5 min) and the coarse transonic smoke were both re-run
green after both changes landed. Net effect: P4 is closed, and P5 starts
without either the transonic stability gap or the wake per-station cost
that this diagnosis session identified as its two open risks.

### P4 supplementary analysis — mesh-refinement study (G4.1 coarse ↔ medium), the surface-Cp sawtooth, and the G4.3 sweep

This subsection collects the two heavier `PYFP3D_TRANSONIC_GATES=1`
evidence figures. They are committed under
`cases/demo/p4_transonic/results/` and regenerated by that phase's
`run_demo.py` in heavy mode (`PYFP3D_TRANSONIC_GATES=1 python
cases/demo/p4_transonic/run_demo.py`, ~40 min for the medium solve + the
10-case sweep); the always-on default demo run produces only the coarse
figure. They answer a question the coarse Cp figure above raises on
sight: **why is the transonic surface Cp visibly serrated (a ~2-cell
sawtooth along the supersonic run) when the P1/P2/P3 curves on the *same*
mesh, drawn by the *same* extractor, are smooth?**

> **⚠️ ATTRIBUTION CORRECTED (P6 N1, 2026-07-08).** The analysis below
> attributes the sawtooth to the artificial-density flux (the integer
> upstream-walk selection). **That attribution is wrong.** The sawtooth is a
> per-triangle wall-gradient **recovery** artifact on the sliver wall
> triangulation, not the flux. Decisive evidence (design.md §3.1/§9.1;
> `cases/demo/p6_surface_recovery/`): nodal/edge-neighbour smoothing of the *same*
> walk solution's wall gradient drops the sawtooth metric ~330× (0.0758 →
> 0.00023), while a *smoother artificial-density flux* (the P6 streamline
> kernel) does not reduce it at all. G6.1 fixes it in post-processing
> (`smooth_wall_tangential_gradients`, `smooth_passes`), not by changing the
> flux. The paragraphs below are kept as the original (superseded) reasoning;
> read them with this correction in mind — in particular the "requires
> changing the spatial operator" conclusion at the end is refuted.

**[SUPERSEDED — see the correction banner above.]** The sawtooth was thought
to be a real feature of the artificial-density flux. The original mechanism
argument, kept for the record:

- `wall_cp_curve` (`post/section_cut.py`) is deliberately triangle-wise:
  each wall triangle crossed by the section plane contributes *one* point
  carrying *that triangle's own* piecewise-constant tangential velocity,
  with **no nodal averaging** (so the sharp-TE crease needs no
  special-casing). P1, P2, and P3 all use this identical function
  (`p3_subsonic/run_demo.py:205`, `p2_kutta_lifting/run_demo.py:233`), and
  their curves are smooth — so the extractor is not the source.
- The difference is entirely in the field being sampled. In the subsonic
  phases the artificial density is off (ν ≡ 0 — the exact G4.2 no-op), so
  φ solves a *symmetric elliptic* density-weighted Laplace problem and its
  per-triangle tangential velocity varies smoothly element to element. In
  the transonic pocket the shipped scheme adds
  ρ̃_e = ρ_e − ν_e (ρ_e − ρ_u(e)) with **ν = max(ν_e, ν_u)**, where the
  upstream element **u(e) is a discrete, integer-valued directional walk**
  through face neighbours (`upwind.py::upstream_elements`). Both u(e) and
  the `max` are *discontinuous* functions of the local velocity direction
  and the local mesh geometry: two adjacent supersonic wall elements
  routinely select geometrically different upstream cells (different ρ),
  so ρ̃ carries a mesh-scale, near-odd-even (≈2h wavelength) perturbation
  riding on top of the smooth density field. That perturbation propagates
  through the density-weighted stiffness into φ and hence into the raw
  per-triangle wall velocity — displayed undamped by the no-averaging
  extractor. It is the unstructured-mesh signature of a scalar
  artificial-viscosity/artificial-density scheme whose dissipation
  operator is only C⁰-rough in its stencil selection.
- **It is confined to the supersonic pocket, exactly where ν > 0.** On the
  coarse upper surface the serration runs from just aft of the LE peak to
  the shock at x/c ≈ 0.60 and then *stops* — the post-shock recompression
  (subsonic, ν ≡ 0) is smooth, and so is the entire lower surface aft of
  its small front pocket. Cp\* (sonic) = −0.435 marks the boundary: the
  jitter lives precisely on the portion of the curve below (more negative
  than) Cp\*. This is the tell that it is the upwind flux, not noise from
  the Γ-secant or the residual tail.

**Mesh-refinement evidence (the sawtooth is O(h) and decays under
refinement).** The medium mesh (61.8k tets, mean edge 0.085) resolves the
same case at ~3.8× the coarse cell count (16.4k tets, mean edge 0.178):

![G4.1 coarse Cp/shock](../cases/demo/p4_transonic/results/g41_cp_shock_coarse.png)
![G4.1 medium Cp/shock](../cases/demo/p4_transonic/results/g41_cp_shock_medium.png)

| Quantity | Coarse (16.4k tets) | Medium (61.8k tets) | Reading |
|---|---|---|---|
| Upper shock x/c | 0.604 | 0.633 | both in the 0.62 ± 0.03 band; drifts aft with resolution (correct direction — coarse dissipation smears the shock forward) |
| Shock sharpness | 1 station, sawtooth pocket | 1 station, near-smooth pocket, steeper jump | serration amplitude drops sharply; shock stays crisp |
| M_max | 1.373 | 1.366 | physical, essentially converged |
| Upper supersonic cells | 58 | 121 | ~2× — finer sampling of the same pocket |
| cp_min (upper) | −1.122 | −1.111 | converging |
| Lower weak shock x/c | 0.358 | 0.364 | stable |
| cl_pressure / cl_KJ | 0.357 / 0.364 | 0.349 / 0.354 | sign-consistent, moves < 3% under refinement |
| Kutta \|F\| | 1.50e-4 | 1.23e-4 | both < 2e-4 tol |

The visual and quantitative story agree: refining h roughly halves the
element scale, the serration amplitude falls with it, the shock sharpens
and moves slightly aft to its better-resolved position, and every
integrated quantity (cl, M_max, shock location) moves < 3%. That is the
signature of a **bounded, mesh-convergent discretization artifact**, not
a physical oscillation, an instability, or a convergence failure — the
monotone-shock and no-expansion-shock gates pass on both meshes.

**This is a *spatial-discretization* defect and is essentially orthogonal
to the Newton work (now P7).** It must not be conflated with the
sharp-shock residual tail, which is a *nonlinear/temporal convergence*
problem that P7 Newton is designed to cure (consistent linearization →
residual to ~1e-10 instead of the engineering tail). P7 Newton changes
*how the nonlinear system is solved*, not the discrete flux itself: run
over the *same* rough ρ̃ operator (discrete integer-valued u(e) +
`max(ν_e, ν_u)`), Newton would converge the serration *more* tightly, not
remove it — indeed it could sharpen the sawtooth once the
θ-damping/under-relaxation that partly masks it is gone. (There is a real
coupling in the other direction, though: the same integer-walk + `max`
non-differentiability that causes the sawtooth also blocks an *exact*
Newton Jacobian, which is why the fix is scheduled as P6, ahead of P7.)
Removing the sawtooth requires changing the *spatial* operator — this is
exactly the new **P6 phase** (roadmap.md P6): (a) a
directionally-consistent / smoother upwind-density flux — a continuous
streamwise-projected density bias in place of the discrete face-neighbour
hop selection and the `max` switch; (b) curved / isoparametric wall
elements (separately scoped, design.md §5.1.2); or (c) mesh refinement,
already shown here to be O(h) convergent but explicitly *rejected* as the
accepted fix (it only reduces amplitude, never removes the oscillation,
and is prohibitive in 3D). It is also explicitly *not* a Δτ/θ-damping
problem (Addendum 3) — the damping change left the shock position and the
serration character unchanged.

**G4.3 parameter sweep (10 cases, re-run under the `damping_theta`
default).** Two α-lines (0.0° and 1.25°) across M∞ ∈ {0.74, 0.76, 0.78,
0.80, 0.82}, gating on convergence + physicality (not an exact cl):

![G4.3 sweep dashboard](../cases/demo/p4_transonic/results/g43_sweep_dashboard.png)

| Trend | Evidence | Reading |
|---|---|---|
| Shock migrates aft with M∞ | upper x/c 0.244 → 0.527 (α = 0°), 0.346 → 0.723 (α = 1.25°) | monotone, physically correct — stronger freestream pushes the shock back |
| Lift builds with M∞ | cl 0.246 → 0.458 (α = 1.25°); ≈ 0 for α = 0° (−6e-4 → −8e-4) | correct: the symmetric α = 0° case carries no lift, the lifting line climbs monotonically |
| All 10 converged | `converged = True`, `n_limited = 0` every row | no limited/floored cells anywhere in the envelope |
| Kutta closure | \|F\| 5.5e-6 → 1.9e-4 | all < 2e-4; the M0.82/α = 1.25° corner (1.9e-4) sits closest to tol |
| M_max stays physical | 1.03 (M0.74/α0) → 1.42 (M0.82/α1.25) | rises smoothly with M∞ and α, never near the limiter cap |

The α = 0° line is the useful control: zero lift to 6e-4 across the whole
range confirms the wake/Kutta machinery injects no spurious circulation
into a symmetric case even as the shock strengthens. The one flagged
datum — the M0.82/α = 1.25° cl of 0.458 (shifted from 0.389 under the old
global damping, and its \|F\| the closest to the 2e-4 boundary) — is
recorded as evidence of a different-but-converged, still-physical state on
the same bounded engineering-convergence tail, since G4.3 gates on
convergence and physicality rather than an exact cl value. A future
θ-sensitivity study should treat that corner as the envelope edge.

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

## P5 — 3D validation: ONERA M6 (closed 2026-07-08; V6 < 1% deferred to P9 curved elements)

**Purpose.** First 3D transonic validation on the ONERA M6 half wing
(M∞ = 0.84, α = 3.06°): the swept-wake / symmetry / far-field pipeline
computes, the tip does not diverge, the λ-shock signature is present, and the
3D post-processing (sectional Cp at η = 0.44/0.65/0.90 plus an inboard
η = 0.20 panel, spanwise Γ(η), the 3D Kutta–Joukowski consistency check,
planform Cp map) is exercised against the user-committed **viscous** AGARD
experiment as a qualitative overlay. This section records the closed-out
evidence **and the two-investigation diagnosis history that closing the
medium gate required** — the first root-cause attribution was wrong, and the
record of how it was overturned is part of the evidence.

**Case setup.** `cases/meshes/onera_m6/{coarse,medium}.msh`. New code:
`post/surface.py::planform_area` + `cl_kj_3d`, `post/section_cut.py::
section_cp_curve` (`wall_cp_curve` refactor verified bit-identical);
`solve/continuation.py` forwards `rtol` (**rtol=1e-7 ~5.5× faster, M_max
identical to 5 digits**; default 1e-10 keeps the P4 path bit-identical).
Final calibrated recipe: `seed40 / eval300 / gamma10 / rtol1e-7` **+ the two
3D closure ingredients landed 2026-07-08** — `farfield_spanwise_gamma=True`
(Γ(z)-tapered vortex far field, 0 at/beyond the sheet tip) and
`n_kutta_polish=4` (fixed-Γ Kutta-closure polish, `omega_rho_polish=0.5`).
All runs cap `NUMBA_NUM_THREADS=16`. Reference:
`reference_data/onera_m6_experiment/` (AGARD AR-138 Test 2308, viscous —
qualitative overlay, **not** a point-wise gate for the inviscid FP solver).

**Key figures (post-fix).**

![sectional Cp vs AGARD experiment](../cases/demo/p5_onera_m6/results/g51_sections_coarse.png)
![spanwise circulation and loading](../cases/demo/p5_onera_m6/results/g52_spanwise_coarse.png)
![upper-surface Cp planform map](../cases/demo/p5_onera_m6/results/g51_surface_cp_coarse.png)
![medium sectional Cp](../cases/demo/p5_onera_m6/results/g51_sections_medium.png)

**Measured results (post-fix; demo 16/16 PASS).**

| Check | Coarse (55.5k) | Medium (350.7k) | Criterion |
|---|---|---|---|
| physical + tip stable (M_max, floored/limited) | **1.398, 0/0 PASS** | **1.995, 0/0 PASS** | M_max < cap, zero floored/limited |
| upper shock x/c η=0.44/0.65/0.90 | 0.596 / 0.570 / 0.425 | 0.594 / 0.526 / 0.345 (~1 cell) | present, monotone, forward-migrating |
| Γ root → tip | 0.097 → 0.0206 | 0.097 → 0.0151 | smooth band-mean decay, Γ_tip → 0 |
| V6 consistency \|CL_p−CL_KJ\|/CL_KJ | 2.40% | 1.82% | **< 3% floor bound (re-spec; < 1% post-P6)** |
| CL (pressure) | 0.2419 | 0.2453 | reported, coarse→medium convergence |
| Kutta \|F\|_max | 5.3e-4 | 5.8e-4 | reported (28× tighter than pre-fix 1.6e-2) |
| demo self-checks | **8/8 PASS** | **8/8 PASS** | — |

**Diagnosis history — why the gate stayed open for two investigations.**
The pre-fix medium solve failed `physical`: M_max 5.204, 8 floored / 4
limited; 26/350 718 cells M > 2, split **18 on the wing at the outboard TE**
(z/b 0.80–0.81) + **8 at the far-field sphere** beyond the tip
(`diagnose_medium.py` reproduces every audit below from the cache).

*Investigation 1 (A–E, `INVESTIGATION_gamma_smoothing.md`):* the cluster
co-locates with the steepest spanwise-Γ roll-off, so a spanwise-Γ-noise
hypothesis was tested via a Gaussian smoother — and **refuted** (flattening
\|dΓ/dz\| 7× at fixed Γ left the spike; in-continuation smoothing made it
worse). The investigation then attributed the cluster to a sharp-TE P1
discretization singularity (G1.6 family) on the strength of refinement-
sharpening (coarse 1.47 → medium 5.20) + well-shaped cells, and concluded
"NOT a wake/Kutta change". **That attribution was wrong.**

*Investigation 2 (T1–T4, `INVESTIGATION_kutta_closure.md`) — the actual root
cause:* the per-station Kutta audit showed the mismatch was a **single-station
anomaly**: st133 (z/b = 0.801) was left **32% under-circulated** (Γ 0.0431 vs
its own target 0.0592; \|F\|/Γ = 37% there vs ≤ 5% at all 165 other stations).
T1: no tet contains both master+slave wake nodes (0/350 718 — jump-in-cell
ruled out). T2: better density convergence at the same Γ leaves the defect
bit-identical (under-convergence ruled out). T3 (decisive): setting ONLY
st133's Γ to its own Kutta target collapses the cluster 18 → 0 cells (band
M_max 3.10 → 1.16) on the same mesh — **a 1/h singularity cannot depend on
Γ**; the refinement-sharpening was the finer mesh resolving a sharper
overspeed around the same Γ *error*. Root cause: the continuation's
per-station secant does not converge at the top Mach level on the 3D mesh
(pushing to 16 evals diverges, M_max ≈ 29 — the 10-eval budget was
early-stopping regularization); the ρ̃-floor + limiter then froze the
under-circulated station's TE overspeed into a spurious M ≈ 3 state. The
8 far-field cells were the independent span-uniform-2D-vortex branch-ray
artifact beyond the tip (no wake cut exists there to carry the prescribed
jump). A latent probe-assignment degeneracy (adjacent stations share Kutta
probe nodes on the swept unstructured TE; st133/134 share their upper probe)
is why st133 specifically stalled — recorded as a known-robustness item.

**Fix + verification.** `farfield_spanwise_gamma=True` removes the far-field
cluster at the boundary-data level; `n_kutta_polish=4` closes the stuck
station by a secant-free damped fixed-point iteration (\|F\| halves per step
to 5.8e-4; omega_rho < 1 *inside* the secant is NOT the fix — measured to
diverge). Both default off; the full default suite is bit-identical
(140 + 4 + 2). From-scratch medium: M_max 1.995 — the mild bounded
tip-TE-corner P1 overshoot, the only surviving trace of the singularity
family — with **0 floored / 0 limited**.

**G5.2 V6 re-spec (user-approved).** V6 is a systematic O(h)
CL_p-below-CL_KJ offset (coarse 2.40% → medium 1.82%, sign-constant) from the
sharp-TE/LE P1 wall gradient and the P4 surface-Cp sawtooth — both **P6**
targets; removing the M>2 clusters left V6 unchanged, proving it independent
of the wake/far-field defects. It is therefore reported against a 3% floor
bound here, and the original < 1% acceptance moves to a post-P6 re-measure
(roadmap P5/P6).

---

## P6 — surface-pressure recovery (G6.1–G6.4, closed 2026-07-08)

**Purpose.** Remove the non-physical ≈2-cell **sawtooth** in the supersonic-run
surface Cp. The N1 investigation root-caused it — decisively — to a *per-triangle
wall-gradient recovery artifact* on the sliver prism-split wall triangulation,
**not** the artificial-density flux: on the same solution, nodal/edge-neighbour
smoothing of the wall gradient drops the metric ~330×, whereas a smoother flux
(the P7 streamline kernel) does not reduce it at all. The fix is a normal-gated
recovery smoothing in post-processing (`smooth_wall_tangential_gradients`,
`smooth_passes`; design.md §9.1), applied to the Cp curve and the force integral,
gated so it never averages across the sharp TE.

**Case setup.** The P4 (NACA0012 M∞ 0.80, α 1.25°) and P5 (ONERA M6 M∞ 0.84,
α 3.06°) cases are re-run through the walk solver; the recovery is a
post-processing change, so the solve (φ, Γ, shock, M_max, cl_KJ) is unchanged —
only the reported Cp/forces move. `smooth_passes = 1`.

**Key figures.**

![G6.1 wall Cp raw per-triangle vs recovery-smoothed (coarse NACA0012 M0.80)](../cases/demo/p6_surface_recovery/results/g61_cp_raw_vs_smoothed_coarse.png)
![G6.1 the finding in one panel: a smoother FLUX (kernel) does not remove the sawtooth; smoothing the RECOVERY does](../cases/demo/p6_surface_recovery/results/g61_flux_vs_recovery_coarse.png)
![G6.1 ONERA M6 coarse section Cp at eta=0.65, raw vs smoothed](../cases/demo/p6_surface_recovery/results/g61_m6_section_cp_coarse.png)

**Measured results (coarse NACA0012 M0.80, walk solution).**

| quantity | raw per-triangle | recovery-smoothed |
|---|---|---|
| sawtooth metric (upper) | 0.0758 | **0.00023** (330×) |
| slope reversals (upper) | 39 | 1 |
| upper shock x/c | 0.604 | 0.607 |
| cl_pressure | 0.3572 | 0.3562 (−0.3%) |
| cd_pressure | 0.01758 | 0.02014 (~15%; near-field FP drag, untrusted — §9) |
| cl_KJ / M_max (solution) | unchanged | unchanged |

**Conclusion & analysis.** The sawtooth is eliminated on the coarse mesh with no
refinement and no solver change — a recovery, not a flux, problem. The shock and
cl_KJ are untouched; cl_p shifts −0.3% and the (explicitly-untrusted) near-field
cd_p shifts ~15%. The one-panel figure is the evidence that the flux is not the
cause: the walk flux and the smoother streamline-kernel flux serrate identically
(raw recovery), while smoothing the same walk solution's recovery is clean. This
**overturns** the P4-supplementary attribution (see the correction banner in the
P4 section). On ONERA M6 coarse (η=0.65) the same smoothing drops the section-Cp
metric 0.145 → 0.046 (3D confirmation, gated part). **G6.3 nuance (measured):**
smoothing helps the *Cp curve* but not the *loads* — the M6 V6 consistency
|CL_p−CL_KJ|/CL_KJ *worsens* 2.40% → 3.35% under smoothing (the ±sawtooth cancels
in the integral; averaging smears the LE peak, moving CL_p ~1% further below the
trustworthy CL_KJ). So the whole V6 floor is the sharp-TE/LE P1 wall gradient →
tracked to **P9** (curved wall elements), and `smooth_passes` is recommended for
the Cp curve only, `= 0` for `wall_force_coefficients`. Demo: `cases/demo/p6_surface_recovery/` (6/6 PASS incl. the gated M6
part); the differentiable flux itself is **P7** (Newton prerequisite), not a
sawtooth fix.

---

## P7 — frozen-selection differentiable walk flux (G7.3, closed 2026-07-10)

**Purpose.** Deliver the P8 fully-coupled-Newton prerequisite: the exact
sensitivity `∂ρ̃/∂φ` of the shipped P4 **walk** artificial density at **frozen
upstream selection** u(e) (design.md §3.1/§6.3, López Appendix B.3–B.6). Scope
decision (with the user): the forward flux is byte-untouched — no `max_ε`, no
flux replacement — so G7.1/G7.2 hold by construction (locked by V0/G4.2 re-runs
+ a forward-path regression guard) and the phase is exactly the derivative +
its finite-difference verification.

**What was built.** `physics/isentropic.py::mach_squared_derivative_wrt_q_sq`
(dM²/dq² = M∞²[1+(γ−1)/2·M∞²]/D², strictly positive) and
`kernels/upwind.py::rho_tilde_sensitivities_sweep` +
`UpwindOperator.rho_tilde_sensitivities` (walk mode only): branch-wise
`(s_e, s_u) = (∂ρ̃_e/∂q²_e, ∂ρ̃_e/∂q²_{u(e)})` — subsonic `(ρ'_e, 0)`;
accelerating ν=ν_e `(ρ'_e(1−ν) − (ρ_e−ρ_u)ν'_e, ν·ρ'_u)` (B.3+B.4);
shock-point ν=ν_u `(ρ'_e(1−ν), ν·ρ'_u − (ρ_e−ρ_u)ν'_u)`; floored /
self-upstream flat branches → 0, exactly mirroring `rho_tilde_sweep`'s clamp.
The DOF chain `∂q²/∂φ_k = 2∇φ·∇N_k` stays with the caller — P8's Term-2/Term-3
assembly consumes `(s_e, s_u)` as the physics factor.

**Verification method.** Directional (JVP) central difference against the
*shipped* `rho_tilde_sweep` with u(e) held frozen — the forward flux is reused
verbatim, so the check verifies exactly the derivative P8 relies on.

**Key figures.**

![V7.1 analytic vs FD scatter + per-regime rel-err histogram (constructed multi-regime field)](../cases/demo/p7_diff_flux/results/v71_fd_scatter_constructed.png)
![V7.3 frozen-selection regimes over the real supersonic pocket (converged G4.1 coarse field)](../cases/demo/p7_diff_flux/results/v73_regime_map_converged.png)
![V7.4 FD accuracy on the converged G4.1 field](../cases/demo/p7_diff_flux/results/v74_fd_scatter_converged.png)

**Measured results (gate < 1e-6).**

| field | max rel err | regimes exercised |
|---|---|---|
| structured cube, 4 fields (`tests/test_p7_diff_flux.py`) | **3–5e-10** | subsonic / accelerating / shock-point / self-upstream / floored |
| constructed multi-regime, NACA coarse 16.4k elements | **3.5e-9** | 4.1k subsonic / 6.0k accelerating / 6.3k shock-point |
| **converged G4.1 M0.80 field** (the P8 target state) | **5.7e-9** | pocket = 1189 accelerating + 977 shock-point, M_max 1.3729 |

**Conclusion & analysis.** The derivative is exact to FD-noise level in every
frozen-selection branch, including on the real converged transonic field — the
P8 Newton Term-2/Term-3 physics factor is in place, sparse (~+1 upstream
element/row), with the forward P4/P5/P6 paths bit-identical (suite 165 passed +
4 skipped + 2 xfailed). Two findings worth their record: (1) **sign
arbitration** — the FD gate settled the design.md §6.3 chain to
`dμ/dM² = +M_c²/M⁴` (the doc's "−" was a transcription typo, fixed); (2) the
**C⁰ kink locus is real but measure-zero in practice** — an FD probe straddling
the max(ν_e,ν_u) tie or the switch threshold reads a branch *average* (~1e-5
apparent error, not a bug); on generic fields only 0.04 % of elements (2/16.4k
on the converged field) sit inside the ε-neighbourhood, but symmetry-degenerate
(separable) fields on structured/prism-split meshes park whole element slabs
exactly on the tie — the measured trap is documented in the test docstrings,
and any future FD check must use generic/noise-broken fields. Demo:
`cases/demo/p7_diff_flux/` (7/7 PASS incl. the gated converged-field part).

---

## P8 — fully-coupled Newton (closed 2026-07-11: G8.1 + G8.2 + G8.3)

**Purpose.** Replace the Picard/secant iteration with a fully-coupled
(φ_red, Γ) Newton on the exact Jacobian (design.md §6.3 at frozen selection,
§8.1 coupled system): quadratic convergence to the actual discrete solution,
Kutta closed as an unknown (no secant–density coupling — the P5 instability
class), and the speed to retire the ~10⁴-iteration Picard budgets.

**What was built.** N2: `kernels/jacobian.py::assemble_newton_jacobian` —
Terms 1+2 fused on the shared Picard CSR pattern, Term 3 (upstream coupling,
graph-distance ≤ 4) as active-set COO rebuilt per step; JVP FD-verified to
~1e-10. N3/N4: `solve/linear.py::solve_gmres` + `solve/newton.py` — one shared
`eval_residual` path, exact δΓ elimination with the far-field vortex column
FD-guarded, GMRES+AMG and **direct (splu + Woodbury)** linear paths. N5: the
transonic robustness chain — **direct exact steps** (the shock-position soft
mode stiffens under refinement; η-accurate Krylov steps stall: measured
frozen-system residual flat at 3.7e-6 with GMRES converging to η),
**stall-adaptive freeze** of the upwind assignment with **active-set refresh**
(the 2.5D prism-split mesh parks ~10³ elements in the max(ν_e,ν_u) near-tie
band — the P7 kink trap in Newton form; live Newton limit-cycles there,
measured branch flips 300–800/step), two-cycle acceptance with the honest
`residual_unfrozen` floor, and freeze-revert / level-fail-fast /
best-of-tried-line-search safety nets.

**★ Baseline findings (user-arbitrated 2026-07-11; roadmap P4 erratum).**
(1) The P4 Picard "engineering-converged" states are **not discrete
solutions**: the coupled Newton residual at the committed coarse M0.80/α1.25
state is **2.2e-4**, and Newton started from it walks in 6 quadratic steps to
the true solution — **shock 0.658, cl 0.459, M_max 1.408** (dissipation-scan
robust, continuation-path independent to the last bit). (2) On the medium
mesh the solution family steepens into the **FP non-uniqueness fold**:
Newton-converged M0.775 → shock 0.570/cl 0.396 (residual 1.8e-13), M0.7875 →
shock 0.674/cl 0.523 (7.9e-11); **no reachable isolated solution at M0.80**
(M_max ≈ 1.45, beyond the isentropic validity envelope — conservative FP
over-lifts strong-shock cases vs Euler, Holst PAS 2000). G8.1 was therefore
re-specced to coarse M0.80 + medium M0.7875 with regression-lock physics
bands; `cases/reference_data/naca0012_m080/` untouched (hard rule 6).

**Key figures.**

![V8.1a coupled-Newton convergence, coarse (subsonic + M0.80 final level)](../cases/demo/p8_newton/results/v81a_convergence_coarse.png)
![V8.1b medium M0.7875 convergence + V8.2-lite runtime breakdown](../cases/demo/p8_newton/results/v81b_convergence_medium.png)

The sawtooth in both convergence plots is the freeze-refresh cycle working as
designed: each frozen phase collapses quadratically to ~3e-13, the live
re-evaluation jumps to the current assignment-staleness level, and the refresh
contracts it (measured stale counts 693 → 81 → 2 → 0 at M0.70 medium) until
the assignment is self-consistent or its intrinsic discontinuity floor
(~1.3e-7 on medium, reported honestly) is reached.

**Measured results (demo 15/15 PASS, `cases/demo/p8_newton/`).**

| check | value | criterion |
|---|---|---|
| subsonic cl Newton vs P3 Picard | 1.4e-7 | < 0.5 % |
| coarse M0.80 terminal residual / quadratic drops | 3.0e-13; 1.3e-3, 1.8e-3 | < 1e-9; both < 3e-2 |
| coarse shock / cl (regression lock) | 0.6581 / 0.4590 | 0.658 ± 0.012 / 0.459 ± 0.010 |
| coupled Kutta closure \|F\| | 8.3e-17 | machine (secant era: ~1e-4) |
| medium M0.7875 terminal residual / drops | 7.8e-11; 1.7e-3, 1.0e-3 | < 1e-9; both < 3e-2 |
| medium shock / cl (regression lock) | 0.6738 / 0.5234 | 0.674 ± 0.012 / 0.523 ± 0.010 |
| assignment-discontinuity floor (honesty) | 1.3e-7 | < 1e-5, reported |
| medium gate run end-to-end | ~100 s | Picard G4.1 medium: 16m39s, non-solution |

**N6 addendum (2026-07-11) — ONERA M6 + performance, G8.2/G8.3 closed.**
The M6 medium (63k nodes / 351k tets) Newton run at M0.84/α3.06 is
**249.2 s end to end** (mesh+cut 7.3 s, solve 239.8 s, forces + 3 section
shocks 2.1 s) against the 300 s gate — vs 4539 s for the P5 Picard recipe;
the coarse mesh takes 42 s. Both meshes have a reachable isolated Newton
solution at the full M0.84 (the FP-fold contingency planned for this run
never triggered): every continuation level converges with zero dm-halving,
the frozen phases end terminal-quadratic (medium final level
2.6e-7 → 2.1e-10 → 7.0e-15), 0 limited/floored, coupled Kutta |F| ~2e-16.

Two ingredients close the runtime gate: the **lagged-LU direct mode**
(`direct_refactor_every` — on a true-3D mesh the LU fill makes each splu
~18.6 s at 63k dofs, ~100× the thin 2.5D cost, and the every-step-direct
N5 recipe spent 1606 s, 97% in splu; the lagged mode refactors once per
level and drives the steps between with GMRES on the fresh coupled operator
preconditioned by the stale LU at rtol 1e-8, falling back to refactor +
exact Woodbury if GMRES fails — same solution, 6.4× faster) and the P5
**dm=0.05 Mach schedule** (the M6 family is far from the NACA-medium fold).

![V8.1c M6 medium convergence + V8.2 runtime breakdown](../cases/demo/p8_newton/results/v82_m6_medium.png)

**P5-caveat measurement** (the recorded follow-up to the P4 erratum): the
committed P5 Picard states are not discrete solutions either, but the
failure is milder in degree — Newton residual 8.6e-6 coarse / 7.6e-6 medium
(Kutta |F| ~5.5e-4) vs the P4 stall's 2.2e-4. The Newton true solutions:

| quantity | P5 Picard (committed) | Newton true solution |
|---|---|---|
| cl_p coarse | 0.2419 | **0.2560 (+5.8%)** |
| cl_p medium | 0.2453 | **0.2646 (+7.9%)** |
| cl_KJ medium | 0.2499 | 0.2692 |
| shocks η44/65/90 (medium) | 0.594/0.526/0.345 | 0.596/0.541/0.362 |
| M_max (medium) | 1.995 | 2.13 |
| Kutta \|F\| | 5.8e-4 (secant+polish) | ~2e-16 (coupled unknown) |

The under-convergence lives in the circulation/lift, not the shock
positions. cl_KJ 0.2692 narrows the inviscid-vs-Tranair/KRATOS (0.288) gap
assigned to P9 from 0.043 to 0.019. The P5 gates stand as Picard-quality
gates (roadmap P5 ledger note).

**G8.3**: the default regression suite is **301.66 s (5m02s)** at
NUMBA_NUM_THREADS=16 — 182 passed + 8 skipped + 2 xfailed; every heavy
transonic/M6 gate sits behind PYFP3D_TRANSONIC_GATES=1.

**Conclusion.** P8 closed: G8.1 terminal quadratic convergence on both gate
cases with the FD-verified Jacobian (rel ~1e-10 on the converged pocket),
G8.2 M6 medium end to end in 249.2 s < 5 min, G8.3 CI budget 5m02s < 10 min.
The production path for a 3D transonic case is now: Picard warm levels +
coupled Newton finish, ~18× faster than the Picard recipe and converging to
the actual discrete solution.

---

## P8 capability assessment — cross-case evaluation demo (2026-07-11, NOT a gate)

**Purpose.** Post-P8 stock-taking requested by the user: run the production P8
Newton solver over the geometry × mesh matrix and measure, in one reproducible
place, (a) convergence behaviour (residual, Kutta closure, circulation → KJ
lift, per Mach level — using the new `gamma_history`/`level_results` solver
instrumentation added for this demo, additive keys only, suite bit-unchanged
at 182+8+2), (b) section-Cp accuracy against the available references, and
(c) end-to-end cost — as the evidence base for choosing the next track (curved
walls vs Track V viscous vs Track B level-set wake; on 2026-07-11, after this
demo, the user inserted the P9 discrimination phase and the P10
Newton-usability phase — curved walls are now P11, backlog P12). This demo asserts
convergence quality and regression locks but does NOT close or claim any
roadmap gate. Demo: `cases/demo/p8_capability/` (part 1 NACA coarse always;
the full matrix under `PYFP3D_TRANSONIC_GATES=1`, ~23 min with the 16-thread
cap `NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16`).

**Case matrix and measured results (36/36 PASS, `results/checks.csv` +
`results/summary.csv`).**

| case | mesh (nodes/tets) | condition | levels / Newton steps | final ‖R‖∞ | Kutta ‖F‖ | cl_p (cl_KJ) | shock x/c | end-to-end |
|---|---|---|---|---|---|---|---|---|
| NACA sub | coarse 5.6k/16.4k | M0.50/α2.00 | 1 / 2 | 4.7e-13 | 0 | 0.2776 (0.2781) | — | 3.6 s |
| NACA sub | medium 20.9k/61.8k | M0.50/α2.00 | 1 / 2 | 2.1e-13 | 0 | 0.2844 (0.2844) | — | 13.6 s |
| NACA tr | coarse | M0.78/α1.00 | 5 / 29 | 8.3e-11 | 0 | 0.2626 (0.2658) | 0.486 | 4.1 s |
| NACA tr | medium | M0.78/α1.00 | 9 / 252 | 2.0e-13 | 0 | 0.3238 (0.3257) | 0.555 | 54.1 s |
| NACA tr (fold attempt) | coarse | M0.78/α1.25 | 5 / 35 | 2.2e-11 | 0 | 0.3399 (0.3445) | 0.522 | 4.4 s |
| NACA tr (fold attempt) | medium | M0.78/α1.25 | 7 / 194 | 2.1e-11 | 0 | 0.4339 (0.4372) | 0.602 | 44.2 s |
| ONERA M6 | coarse 11.0k/55.5k | M0.84/α3.06 | 4 / 35 | 6.9e-12 | 1.7e-16 | 0.2560 (0.2621) | 0.600/0.573/0.429 | 13.4 s |
| ONERA M6 | medium 63.2k/350.7k | M0.84/α3.06 | 4 / 47 | 7.0e-15 | 2.1e-16 | 0.2646 (0.2692) | 0.596/0.541/0.362 | 256.5 s |

All eight runs: 0 limited / 0 floored, Kutta closed to machine precision
(secant era: ~1e-4), terminal super-linear collapse (assessment band 5e-2 on
the best consecutive drop pair; the G8.1 gate cases keep their 3e-2 in the
gated tests — the one 3.7e-2 pair, NACA coarse α1.0, is a warm start already
at 6e-6 leaving only a 2-step tail: 1.64e-7 → 8.3e-11). Subsonic reference:
corrected 2D panel bracket [PG 0.2788, KT 0.2919] — medium 0.2844 inside the
bracket, −0.3% of the midpoint (P3 G3.2 semantics); coarse −2.7%. M6 medium
reproduces every G8.2 regression lock (cl_p 0.2646, shocks 0.596/0.541/0.362,
M_max 2.129, 257 s < 300 s).

**★ Fold-zone grid-sensitivity finding (the user's contingency ladder
exhausted).** The NACA transonic pair was specified SAME-condition
M0.78/α1.25 on both meshes with an α→1.0 fallback if the FP fold interferes.
Both meshes converge cleanly at BOTH alphas — but to points far apart on the
fold-steep solution family: α1.25 coarse shock 0.522/cl 0.3399 vs medium
0.602/0.4339 (Δcl 0.094); after the rule-mandated rerun, α1.0 STILL gives
0.486/0.2626 vs 0.555/0.3238 (Δcl 0.061 > the 0.05 comparability band). This
is the G8.1 fold finding in grid form: the measured family slope dcl/dM ≈
6–10 in the M0.775–0.80 zone means O(h) discretization differences act like
an O(0.01) M∞ shift — same-condition grid comparison is intrinsically
ill-conditioned near the fold, and NONE of the four states is a solver
failure (all are true discrete solutions, terminal-quadratic, 0 lim/flr).
The demo therefore regression-locks each mesh's own Newton solution
(±0.012 shock / ±0.010 cl, G8.1 semantics) instead of asserting a
grid-convergence band, and reports both attempts (`summary.csv`; the α1.25
pair is the dashed overlay in the Cp figure). Contrast: away from the fold
the grid behaviour is benign — M6 M0.84 coarse→medium moves cl_p by only
0.0086 and the η0.44 shock by 0.004. The medium α1.0 run also exercised the
N5 robustness chain for real: 9 levels / 252 Newton steps with dm-halving
retries visible in the convergence figure, still finishing at 2.0e-13.

**Key figures.**

![NACA convergence](../cases/demo/p8_capability/results/naca_convergence.png)

![M6 convergence](../cases/demo/p8_capability/results/m6_convergence.png)

![NACA Cp](../cases/demo/p8_capability/results/naca_cp_sections.png)

![M6 Cp](../cases/demo/p8_capability/results/m6_cp_sections.png)

![timing](../cases/demo/p8_capability/results/timing.png)

Presentation notes: Cp curves are plotted with the P6 normal-gated recovery
smoothing (1 pass); forces and shock/regression locks are measured on raw
curves (G6.3/G8.2 protocol). The M6 lift-convergence panels make the V6 gap
directly visible: cl_KJ (from Γ) settles ~2% above the pressure-integrated
cl_p on both meshes. The M6 coarse tip section (η0.90) smears the double
shock the medium mesh resolves — expected at 11k nodes.

**Timing.** Newton end-to-end (mesh+cut / solve / post): NACA medium 54.1 s
(1.0/52.6/0.5), M6 coarse 13.4 s (1.2/10.8/1.4), M6 medium 256.5 s
(7.1/239.5/9.9) — vs the RECORDED Picard ledger baselines (not rerun, cost
rule): NACA medium M0.80 G4.1 999 s to a state the P4 erratum showed is NOT
a discrete solution, and M6 medium P5 solve 4539 s with Kutta |F| 5.8e-4 —
~17.7× the Newton end-to-end that closes Kutta to 2e-16.

**Honest capability boundaries (assessment).**

1. **No non-lifting Newton path**: `solve_newton_lifting` structurally
   requires the wake cut + Kutta/Γ block — a sphere cannot even build the
   `(mesh_cut, wc)` pair (`cut_wake` raises without a wake group).
   Non-lifting bodies run Picard (`solve_laplace`/`solve_subsonic`, P1/P3
   demos) and carry the OPEN G1.6 flat-facet sphere-Cp gap (~11.6% vs the
   2% gate, refinement-saturating ~3.6% at 7M tets; root-caused geometry
   variational crime → P11 curved walls; the missing Newton entry itself is now
   P10/G10.1). The sphere is deliberately absent from this
   Newton demo (user arbitration 2026-07-11).
2. **V6 lift floor (P11 curved walls; attribution under P9 test)**: M6 medium cl_KJ 0.2692 vs the Tranair/KRATOS
   inviscid reference 0.288 — the remaining 0.019 is the sharp-TE/LE P1
   wall O(h) floor (P5/P8 evidence), visible in this demo as the ~2%
   cl_KJ-over-cl_p offset.
3. **Fold-zone conditions are certification-hostile** (finding above):
   single-mesh numbers near the FP fold are not mesh-converged engineering
   data; report them with the family-slope context or move off the fold.
4. **cd_p untrusted** (FP pressure drag on P1 walls, ~15% recovery
   sensitivity — P6 record); not asserted anywhere in this demo.
5. **Viscous effects absent**: the AGARD overlay is qualitative; VII
   (Track V, designed) moves CL down toward experiment (~0.26–0.27) and
   does NOT close the 0.288 inviscid gap (that belongs to P11 curved walls,
   pending the P9 discrimination).

**Development outlook (evidence from this demo; the user arbitrates the
order).**

- **Curved walls (now P11)** attacks the only two measured ACCURACY gaps — G1.6
  (11.6% sphere Cp) and V6 (cl_KJ 0.019 below the inviscid references) —
  both root-caused to the flat-facet P1 wall. Highest leverage on
  reference-matching; medium-high risk (own effort, roadmap).
- **Track V (VII/IBL)** adds the missing physics for absolute CL/CD against
  experiment (the M6 Cp overlays in this demo show exactly the inviscid
  offsets it would shrink); V1 is parallelizable with the accuracy phases, V3 consumes the
  P8 Newton machinery.
- **Track B (level-set wake)** buys geometry flexibility (M2 wing-body);
  no accuracy payoff on the current case set.
- **Discrete adjoint (backlog, now P12)** is cheap to open now: the exact P8 Jacobian +
  `reduce_operator` machinery is the transpose seed; note the fold finding
  also warns that fold-zone gradients will be ill-conditioned.
- The `gamma_history`/`level_results` instrumentation added here is the
  natural hook for all of these (convergence dashboards, VII coupling
  monitors, adjoint checkpointing).

---

## P10 (partial) — G10.2 level-adaptive intermediate continuation tolerance (closed 2026-07-11; G10.1 still open, so the phase stays open)

**Demo:** `cases/demo/p10_newton_usability/run_ab_g102.py` (one-shot A/B
evidence, not a suite test; 16-thread timing protocol). Result: **34 PASS +
6 XFAIL** — the XFAILs are the *documented negative result* for the
fold-zone case, not open defects.

**What shipped** (`pyfp3d/solve/newton.py`): `solve_newton_transonic`
gained the opt-in `intermediate_tol` (default None = bit-identical,
suite-locked); `solve_newton_lifting` gained the loose acceptance
criteria `tol_residual_loose` / `tol_residual_rel` / `accept_on_stall`
(all default off) and reports `accept_reason`
("tol"/"loose_tol"/"rel_drop"/"stall") per solve and per level. Loose
acceptance keeps the 0-limited/0-floored and ‖F‖ guards, requires ≥ 1
Newton step at the level, and never applies inside a frozen phase; only
ORIGINAL-SCHEDULE intermediate levels run loose — the final level and
every dm-halving retry level keep the full strict tol-1e-10 +
freeze/honesty machinery. Suite +2 (`tests/test_p10_continuation.py`);
baseline **184 passed + 8 skipped + 2 xfailed**.

**A/B on the two committed recipes** (`results/summary.csv`,
`results/levels.csv`, `results/*_ab.png`):

| case | default | adaptive (`intermediate_tol=1e-5`) | verdict |
|---|---|---|---|
| ONERA M6 medium M0.84/α3.06 (`NEWTON_M6_RECIPE`) | 239.5 s, 47 steps (intermediate levels 11+12+12) | **140.3 s (+41.4%)**, 18 steps (3+1+2 loose, ending ~1e-5) — final level IDENTICAL: 12 steps, ‖R‖ 7.8e-15, cl 0.2646 / M_max 2.129 / shocks 0.596-0.541-0.362 equal to 4 digits, all G8.2 locks PASS | **PROMOTED** into `NEWTON_M6_RECIPE` (≥ 15% criterion met; gated G8.2 now ~145 s) |
| NACA0012 medium M0.7875/α1.25 (`NEWTON_TRANSONIC_RECIPE`, fold zone) | 79.1 s, 403 steps, 6 dm-halvings, converges 7.8e-11 | ends UNCONVERGED at 4.3e-6: cl 0.369 vs lock 0.523, shock 0.535 vs 0.674 | **NEGATIVE result recorded** — recipe unchanged |

**The fold-zone failure mechanism (the pre-registered P8 trap, now
measured).** Round 1 exposed a degenerate path: warm-started levels ENTER
below any absolute threshold (~2-9e-6 here), so the naive ‖R‖ ≤ 1e-5
clause accepted intermediate levels with ZERO Newton steps — the ramp
becomes a level skip. Two hardenings were added and kept (≥ 1 step
required; dm-halving retries strict), after which the loose ramp does 1–4
steps/level — but near the fold (dcl/dM ≈ 6–10) that still never tracks
the circulation: the final level arrives with an untracked Γ seed and
stalls at the ~5e-6 live-churn floor for the full 60-step budget, and so
do the STRICT retry levels warm-started from the same loose state (round
2: 0.7812 and 0.7781 both 60 steps, no convergence — a strictly
re-converged level cannot repair the seed within budget either).
Conclusion: **loose intermediates are contraindicated in fold zones**;
away from folds (M6 class) they are pure profit. This is the P8 N2–N4
"warm-start only from CONVERGED levels" warning in its G10.2 form, and it
is why the promotion is per-recipe.

**G10.3 no-ramp direct-solve feasibility (closed 2026-07-11; verdict: KEEP
the Mach ramp).** `run_g103_noramp.py` (9 PASS; `results/g103_noramp.csv`,
`g103_summary.csv`, `g103_noramp_convergence.png`): single-level Newton at
the target M∞, seedings s1 = Picard-5 / s2 = Picard-40, four cases.
Measured answer to "what is the ramp actually buying":

| case | s1 (Picard-5) | s2 (Picard-40) |
|---|---|---|
| M6 coarse M0.84 | **A** — same solution, 5.9 s vs ramp 8.0 s, but clamped transient (peak 11) | **A**, clamp-free, 11.0 s (−37% vs ramp) |
| M6 medium M0.84 | **A** — same solution (cl 0.2646, ‖R‖ 6.6e-15), **79.0 s vs ramp 141.4 s (+44%)**, but clamped transient (peak 45; final 0/0) | **A**, clamp-free, 132.0 s (+6.6%) |
| NACA coarse M0.80 | A — 4.0 s (clamped transient + 1 freeze revert; single-case exception) | **C** — the un-continued Picard-40 seed itself diverges (1625 lim/138 flr, M_max at the 3.0 cap) |
| NACA medium M0.7875 (fold) | **C** — stalls at 4.6e-6, cl 0.449 ≠ lock 0.523 | **C** — seed diverges (9934 limited) |

Findings: (1) **far from the fold, branch selection is not the binding
constraint** — no-ramp converges to the identical solution under both
seedings; the ramp's measurable value there is a clamp-free transient.
(2) No seeding satisfies the pre-registered promotion rule (s1 fails the
clamp-free clause despite +44%; s2 is clamp-free but under 20%), so
**no recipe change** — the +44%-but-clamped observation is recorded for
user arbitration if the clamp-free clause is ever to be relaxed.
(3) **Deep Picard seeds are actively harmful without a ramp** (both s2
fold-zone divergences) — consistent with the P4/P5 record that Picard
itself needs continuation at supercritical M∞. (4) Fold-zone no-ramp
fails on medium regardless of seeding — the ramp stays, everywhere, as
shipped. Instrumentation added: `clamp_history` on `solve_newton_lifting`
(additive key).

## P9 — grid-convergence & accuracy-gap discrimination (closed 2026-07-11; G9.3 verdict awaits user arbitration)

**Demo:** `cases/demo/p9_grid_discrimination/run_demo.py` — **11 PASS + 3
XFAIL**, where the XFAILs ARE the result (the fine-mesh failure, with its
root cause), not open defects. The M6 fine solution is npz-cached locally
(gitignored, like P5's); the committed CSVs/PNGs are the evidence.

**The phase asked:** is the 0.019 M6 lift gap (cl_KJ 0.2692 vs
Tranair/KRATOS 0.288) resolution, or a sharp-TE/LE flat-facet floor?
**The answer:** the question was mis-posed — the gap is **unsplittable as
posed**, because the 3D sequence does not converge, and the reason it does
not converge is a *different* defect than either candidate.

### G9.2 (PASS, clean) — a sharp TE imposes NO lift floor

NACA0012 2.5D, M0.5/α2, coarse→medium→fine, all three converged
(|R| ≤ 4.4e-11). Error vs the corrected-panel PG–KT midpoint:

| level | tets | cl | \|error\| |
|---|---|---|---|
| coarse | 16 386 | 0.2776 | 2.71% |
| medium | 61 788 | 0.2844 | 0.33% |
| fine | 239 022 | 0.2853 | **0.03%** |

Monotone to well inside the ±1% clause. **The 2D leg of the
"sharp-TE/LE P1 wall gradient" attribution is gone.**

### G9.1 (INVALID) — ★ the fine mesh is not a discrete solution

| level | tets | converged | cl_KJ | M_max (unlimited) | cells over M_cap=3 |
|---|---|---|---|---|---|
| coarse | 55 531 | yes (6.9e-12, 0/0) | 0.2621 | 1.40 | 0 |
| medium | 350 718 | yes (7.1e-15, 0/0) | 0.2692 | 2.13 | 0 |
| fine | 2 513 255 | **NO** (1.1e-5, 1 limited) | *0.2393* | **7.93** | **9** |

The fine Newton limit-cycles at ‖R‖ ~ 1e-5 for its entire 60-step budget:
permanently speed-limited cells block the N5 freeze machinery (which
requires 0 limited/floored to engage), so the assignment churn never
freezes. Its cl_KJ (0.2393, italic above) is a **limit-cycle artifact, not
a lift** — with only two discrete solutions there is no three-point
Richardson, so the pre-registered bands (≥ 0.283 / ≤ 0.278) **cannot
fire**. Nothing was fabricated: `verdict.csv` records the extrapolant and
both attribution shares as `n/a`.

**★ Where the singularity is — this is the phase's real finding.** All 9
capped cells sit at:

- **z/b = 0.998–1.000** — the wing tip,
- **x − x_TE = +0.002 … +0.017** — *aft* of the trailing edge (so **not on
  the wing at all**),
- **|y| < 0.003** — in the chord plane, i.e. **on the wake sheet**,

that is, at the **free tip edge of the rigid planar wake sheet**, exactly
where `Γ(tip) = 0` is enforced (M1: tip free edges stay single-valued).
This is the classical **vortex-sheet-edge singularity**: a flat sheet of
trailing vorticity (strength −dΓ/ds, largest at the tip where the wing
unloads fastest — **not** the bound Γ, which correctly → 0 there) that
simply *ends* induces a **1/√r flat-plate-edge** velocity at its free edge
(**corrected from the earlier "1/r-type"** — P13/G13.1 measured the
conforming refinement exponent p ≈ 0.59, i.e. 1/√r; a 1/r line vortex would
give p = 1), whereas the real flow rolls the sheet up into a tip vortex.
P5's "bounded tip-TE-corner
P1 overshoot" (M_max 1.995 on medium, recorded then as *the only surviving
singularity trace*) is the **same object**, seen at a resolution too coarse
to reveal that it is unbounded: refinement makes it **worse** (1.40 → 2.13
→ 7.93), not better.

**⇒ It is a WAKE-MODEL defect, not a wall-element defect.** Curved or
isoparametric *wall* elements cannot remove the edge of a wake sheet.

### G9.3 — attribution verdict (user arbitrates)

1. The 0.019 gap **cannot be split** by grid convergence: the 3D sequence
   does not converge.
2. **2D sharp TE: exonerated** (G9.2). **3D blocker: the rigid planar
   wake's tip edge** (G9.1).
3. **Recommendation:** **P11 (curved walls) is not supported by P9 as the
   3D-lift fix.** Its remaining justification is **G1.6** (sphere-Cp on a
   *smooth curved* wall — a different, still-valid mechanism). The 3D
   accuracy route now points at the **tip/wake treatment**: wake **roll-up**
   or an explicit **tip-vortex** model — which must land before *any* 3D
   grid-convergence claim is possible.

   ★ **Correction 2026-07-12 (this was over-stated as "P9 corroborates the
   Track B route").** Track B's level-set wake changes the wake
   *representation*, not the rigid-planar-sheet *model*: it keeps the same
   sheet ending at the tip with Γ(tip)→0. The `cases/demo/p13_tip_edge_singularity/`
   probe (subsonic M0.5, no limiter — the clean geometric test) measures the
   tip-edge peak Mach **diverging under refinement on BOTH the conforming and
   level-set paths** (same mesh, coarse→medium: ×1.38 conforming, ×2.28
   level-set — ★ **ERRATUM 2026-07-14**: the LS ×2.28 is the `element_mach2`
   mixed-plain ×5 metric artifact retired at the B8 re-spec diagnosis; honest
   LS exponent +0.62 ≈ conforming +0.52, both still diverge, see §B8), while
   the wing control stays flat. So Track B does **not** fix
   this singularity; only a genuinely new wake model (roll-up / tip vortex)
   does, and **no current Track B phase does that** (B9 free-wake is shelved,
   and is about O(θ²) deflection rather than roll-up). What Track B *does*
   eliminate is the separate **secant/Kutta-closure** family of M6 defects (the
   P5 st133 stall, the swept-TE probe-sharing degeneracy, the Γ(z) spanwise
   jitter, the far-field branch-ray artifact) — see the B7 section.

### Solver-path findings recorded en route (feed P10; no default changed)

- **`precond="direct"` does not scale to the fine mesh.** A *single* `splu`
  at ~450k dofs ran **4 h 39 min without returning** (RSS 26 GB) and was
  killed — vs **18.6 s** per factorization on medium (63k). True-3D LU fill
  is the wall; this is the N6 finding one mesh level further.
- **The `precond="amg"` fallback is valid *and faster than direct at every
  size* — once the Eisenstat–Walker forcing is tightened to η = 1e-8.**
  Validated against the G8.2 locks on medium *before* spending the fine
  budget: **66 s vs 141 s** direct, same solution to 4 digits (cl 0.2646,
  M_max 2.129, shocks 0.596/0.541/0.362), terminal quadratic in the frozen
  phase, 0 GMRES stalls; coarse **8 s vs 42 s**. N5's "Krylov steps stall on
  the shock-position soft mode" is a property of the **loose default
  forcing**, not of AMG — a candidate P10 recipe change.
- **The fine mesh's cold Picard seed overshoots the LE into the density
  floor.** A cold M0.70 Picard-5 seed lands 4036 speed-limited + 1847
  density-floored cells; level-0 Newton then stalls at ‖R‖ ~ 6e-2, and since
  level 0 cannot dm-halve, the *whole* solve breaks (a M0.50 cold seed still
  floors 658). Fixed at **continuation-path level only** — `m_start = 0.30`
  (deep subcritical, where a crude seed cannot reach the floor) plus
  `n_picard_seed = 12` — giving a clean 0/0 start; the ramp then carries it
  up in 13 levels, **5294 s (88 min)**, inside the 2 h budget. Path changes
  are safe by G8.2's continuation-path independence.

## Track B — level-set embedded wake (B1 ✓ B2 ✓ B3 ✓ B4 ✓ B5 ✓ B7 ✓, closed 2026-07-11/12; B6 ◐)

**What the track replaces.** The conforming path represents the wake as a *mesh
surface*: the sheet is embedded in the geometry, its nodes are duplicated by the
preprocessor, and Γ is a global unknown eliminated by a master–slave constraint
and chased by a secant loop. Track B removes all of that. The wake becomes a
**level set** evaluated on an *unmodified* mesh; elements the sheet passes through
get a second set of DOFs (multivalued / CutFEM-style); the jump is convected by a
wake least-squares condition; and Γ is no longer an unknown at all — it is a
**result**, pinned by a Kutta condition at the TE.

**Purpose (user-arbitrated 2026-07-11): mesh/geometry workflow capability, NOT
solver speed.** No pre-embedded wake surface, α sweeps without remeshing, blunt-TE
anchoring, multi-wake/wake–fuselage intersections, and the structural elimination
of the P5 st133-class Kutta-probe failures. The original "kill the Γ-secant for
speed" motivation is obsolete post-P8 Newton.

**Dual-mesh rule (the acceptance discipline).** Every gate runs on **both** mesh
families:

| family | meshes | role |
|---|---|---|
| **wake-embedded** (the "C-grid" analogue) | M0 quasi-2D, M1 ONERA M6 | nodes lie *exactly* on the sheet ⇒ stresses the ε side-shift at scale, and enables a strict **same-mesh A/B against the conforming solver** |
| **wake-free** (the "O-grid" analogue) | M3 quasi-2D, M4 ONERA M6 | **no `wake` tag at all**; the level set makes generic cuts through generic elements — **the actual workflow target**, where no conforming counterpart exists |

**Evidence map.**

| gate | evidence | checks | verdict |
|---|---|---|---|
| B1 cut-element identification | `tests/test_b1_cut_elements.py` | 34 PASS | closed 2026-07-11 (test-only; no demo dir) |
| B2 multivalued FE assembly | `tests/test_b2_multivalued.py` | 17 PASS | closed 2026-07-11 (test-only; no demo dir) |
| B3 lifting solve + implicit Kutta | `cases/demo/b3_levelset_lifting/` + `tests/test_b3_lifting.py` | 13 demo PASS + 6 PASS | closed 2026-07-12 |
| B4 TE control volume / Kutta | same demo + `tests/test_b4_te_control_volume.py` | 8 PASS | closed 2026-07-12 |
| B5 far-field A/B | `cases/demo/b4p5_farfield/` + `tests/test_b45_farfield.py` | 9 demo PASS + 10 PASS | closed 2026-07-12 |
| B6 transonic (level-set) | `cases/demo/b6_transonic/` + `tests/test_b6_transonic.py` | 14 demo PASS + 9 PASS (+2 gated) | ◐ IN PROGRESS 2026-07-12 — coarse M0.80 gate met, medium M0.7875 fold deferred to LS Newton |

Numerics spec: [design_track_b.md](design_track_b.md) (supersedes DN1;
B6 findings in §10). **B6 in progress** (transonic on the level-set path);
**next = B7** (ONERA M6 3D).

### B6 — transonic on the level-set path (IN PROGRESS 2026-07-12)

B6 carries the multivalued implicit-Kutta solver into the transonic regime.
Three measured findings, each overturning a "transplant the conforming recipe"
default (design_track_b.md §10):

1. **Per-side artificial density** (D10): a cut element has two velocity
   states, so `rho_tilde` is evaluated once per side and the upstream walk
   runs on a same-side-restricted face graph (the wake is a slip line).
   Subcritically an exact no-op; the M0.80 blow-up cells sit in the pocket
   ABOVE the airfoil (zero on the wake strip), so the shock machinery is
   isomorphic to conforming.
2. **Damping must be LOCALIZED to the ν>0 zone.** The P4 whole-field
   θ·diag stabilizer is a Jacobi smoother — near-transparent-yet-throttling
   to smooth global modes. The implicit Kutta makes Γ a smooth global
   SOLUTION mode (conforming keeps it as an outer secant unknown, outside the
   damped matrix), so global damping throttles it: Γ crawls 0.0005→0.017 in
   160 outers vs undamped convergence in 35. `damping_scope="supersonic"`.
3. **Near the FP fold the option-a Γ→vortex feedback has loop gain > 1** —
   Γ climbs monotonically through both the conforming-Picard stall and the
   Newton truth, then blows up; the López **Neumann outlet** (B5 option b,
   no Γ feedback) removes the loop and converges. ⇒ B6 transonic recipe =
   `farfield="neumann"`.

**★ A/B inversion (the headline result):** the raw same-mesh Picard-vs-Picard
Γ gap grows with pocket strength (+10.5% at coarse M0.75), which looks like an
LS error — but same-mesh conforming **Newton** arbitration shows the
*conforming Picard* under-circulates 4–8% (its P4-erratum stall bias,
quantified at weak shocks), while the LS Picard sits within **+0.25–1.0%** of
the Newton truth and converges toward it under refinement. **User arbitrated
(2026-07-12): the B6 gate baseline is the same-mesh conforming Newton truth,
not the conforming Picard.**

Coarse M0.80 α1.25 (Neumann), vs the G8.1 Newton truth (Γ 0.2295 / shock 0.658
/ cl_p 0.459; the conforming Picard stall is Γ 0.1819 = −21%):

| mesh | Γ | vs Newton | shock | cl_p | M_max | lim/flr |
|---|---|---|---|---|---|---|
| M0 wake-embedded | 0.2114 | −7.9% | 0.644 | 0.4154 | 1.39 | 0/0 |
| M3 wake-free | 0.2315 | **+0.9%** | 0.678 | 0.4556 | 1.39 | 0/0 |

**Medium M0.7875 = the FP fold.** With more dissipation (C=2.0, θ=0.5
localized, dm=0.025 — the coarse recipe diverges) the LS solve stays bounded
and physical (M_max 1.44, 0 lim/flr) but stalls at Γ −18.8% of the same-mesh
Newton truth (0.2643): a Picard method — conforming or level-set — does not
reach the isolated fold solution (why G8.1 re-specced the conforming path to
Newton locks). The quantitative medium gate needs the **LS Newton** (post-B6
re-derivation, design_track_b.md §5.5, explicitly deferred).

![B6 stabilizer story: throttle / runaway / converge](../cases/demo/b6_transonic/results/stabilizer_story.png)
![B6 transonic Cp + shock, dual-mesh vs conforming](../cases/demo/b6_transonic/results/transonic_cp_shock.png)
![B6 A/B gap vs Mach](../cases/demo/b6_transonic/results/ab_gap_vs_mach.png)

> **Track-B renumber 2026-07-12.** Two renumbers landed the same day: a new **B4**
> (TE control volume) was inserted, then the half-integer IDs were regularized
> away. The far-field A/B gate documented below as **B5** was called *B3.5*, then
> *B4.5*, in earlier docs — including the demo directory name
> `cases/demo/b4p5_farfield/`, which is kept as-is so the committed paths stay
> stable. See roadmap.md for the full mapping.

---

### B1 — level-set wake + cut-element identification (closed 2026-07-11)

**Evidence:** `tests/test_b1_cut_elements.py` — **34 passed**, across the full
dual-mesh matrix (2.5D M0/M3 coarse+medium *and* 3D M1/M4 ONERA M6). No demo
directory: B1 delivers no solve, only a geometric predicate, so its evidence is
the test matrix rather than a figure.

**Deliverable.** `wake/levelset.py` + `wake/cut_elements.py`: a TE-**polyline**
ruled level set, and the census of elements it cuts. The mesh is never modified.

| check | measured | why it matters |
|---|---|---|
| M0 (embedded): cut census vs the conforming wake | **exactly** == `cut_wake`'s minus-side element star, element by element | the level set reproduces the conforming topology it replaces |
| M0: on-sheet nodes | every one ε-shifted **"+"** | D4 side-shift stress test at scale |
| M3 (wake-free): corridor TE → far field | gap-free at α=0 **and** after `update_direction` to α=4° **on the same mesh** | α sweeps without remeshing — the workflow payoff |
| M1/M4 (ONERA M6, 3D): census vs conforming | strict **superset**: **0 missing**, **+2.9%** extras, all tip-edge straddlers | expected — the sheet's tip *edge* need not conform in an embedded method |
| M1/M4: spanwise clip | verified (nothing cut wholly outboard of the tip) | encodes Γ(tip) = 0 |

**★ Two 3D-only mechanisms found and fixed here — both invisible on quasi-2D
meshes.** This is the concrete justification for the dual-mesh rule:

1. **The swept TE span axis is not perpendicular to the wake direction.** The
   spanwise coordinate must be measured in the **oblique (v, d̂, n̂) frame**. An
   orthogonal projection leaks the downstream distance into the spanwise
   coordinate and wrongly clipped **~60% of the true M6 cut set** (measured, then
   fixed, now regression-pinned).
2. **The spanwise clip is mandatory** (crossings must satisfy 0 ≤ q ≤
   span_length). Without it the level set cuts the wake-plane *extension beyond
   the tip* — i.e. it re-creates P5's far-field **branch-ray artifact**. The
   conforming path gets the same semantics for free from its free-edge rule.

**The meshes the rule needs** (Track M deliverables M3/M4, built the same day).
Left: the M3 wake-free quasi-2D layer — note there is no wake line in the
topology, only a size-field corridor. Right: the M4 wake-free ONERA M6 corridor.

![M3 wake-free quasi-2D layer](../cases/meshes/naca0012_wakefree_2.5d/coarse_layer.png)
![M4 wake-free ONERA M6 wake corridor](../cases/meshes/onera_m6_wakefree/coarse_wake_corridor.png)

M3 coarse: 29,250 tets, corridor median edge 0.0595 vs an h_wake target of 0.06.
M4 coarse/medium: 50,605 / 329,645 tets — **within 6–9% of the M1 counts at equal
h_wall**, which is precisely what makes the B7 A/B against the P5/P8 baseline a
*controlled* comparison.

---

### B2 — multivalued FE assembly (closed 2026-07-11)

**Evidence:** `tests/test_b2_multivalued.py` — **17 passed** (coarse+medium of both
2.5D families, M6 coarse of both 3D families; some parametrizations skip in CI
where the meshes are gitignored). Test-only, like B1: B2 delivers the assembly, not
yet a lifting solve.

**The key design insight.** A cut element is *the same P1 element matrix assembled
twice*, once with `dofs_upper` and once with `dofs_lower`. That is expressible as a
sparse **column redirection** of the ordinary single-valued matrix: on a cut
element, entries whose two nodes lie on **opposite** sides move their column
main(b) → aux(b) (`multivalued_redirection_coo`). Everything else is byte-identical
to `PicardOperator.assemble_matrix()`. There is **one mesh and one extra DOF per cut
node** — *not* two meshes (López fig. 3.6's "two meshes" is a visualization only).

At B2 the aux rows carry a **continuity ("weld") closure** aux_k = main_j, which
forces the jump to zero. That makes the extended system a strict generalization of
the single-valued one — and that is exactly what B2 proves:

| gate | measured | criterion |
|---|---|---|
| extended matrix folds back to the single-valued stiffness | **1e-13** | the weld closure degenerates *exactly* |
| V0 freestream (φ = U·x) on the cut mesh, 2.5D M0/M3, α=0 and α=4° | **0.0** (exact linear field) | < 1e−12 |
| V0 freestream, 3D M6 M1/M4 | **1.1e−14 / 3.4e−14** | < 1e−12 |
| V1 MMS convergence slope (cube cut by an 8°-tilted half-plane, generic position) | **1.94** | ≥ 1.9 |
| Laplace at α = 0 ⇒ cl ≈ 0, main potential vs the single-valued `solve_laplace` oracle | **~3e−11** | the weld forbids a jump ⇒ cl_KJ = 0 |

**Recorded consequence:** the extended matrix is structurally **nonsymmetric** (the
weld rows), so CG is inapplicable — B2 solves by sparse-direct LU (`spsolve`), and
GMRES+AMG is the B3+ scaling path.

---

### B3 — lifting solve with implicit Kutta (closed 2026-07-12)

**Demo:** `cases/demo/b3_levelset_lifting/` — `python run_demo.py`, ~1 min,
**13/13 PASS**. NACA0012 medium, incompressible, α = 0 and α = 4, **on the same
mesh with the same level set**. The mesh topology knows nothing about the wake.

B3 replaces B2's weld closure with the real thing: the TE jump is carried by the
multivalued aux DOFs, the g₁+g₂ wake LS convects it downstream, and its **value**
is set by B4's TE Kutta condition. **Γ is a RESULT** — no secant, no master–slave
constraint, no Γ unknown.

| check | measured | criterion |
| --- | --- | --- |
| Γ at α = 0 — M0 embedded / M3 wake-free | −3.89e−4 / −4.15e−4 | \|Γ\| < 1e−3 (symmetric ⇒ no circulation) |
| Γ at α = 4 — M0 embedded | **0.2384** (conforming 0.2393, **0.4%**) | > 0.2 |
| Γ at α = 4 — M3 wake-free | **0.2339** | > 0.2 |
| aux DOFs / main DOFs | **9.5%** | < 15% — the enrichment is a thin strip |
| cut tets | 2982 of 61788 (**4.8%**) | (context: the level set touches ONE element layer) |
| TE Kutta control volumes | 2 upper / 2 lower, wall-adjacent | both non-empty |
| jump drift TE → far field | **0.0%** | < 10% (no drain) |
| [φ] near the TE vs the reported Γ | 0.0% | < 10% |
| cl_p vs conforming (α = 4, M0) | 0.4770 vs 0.4786 | within 3% |
| cl_p vs cl_KJ = 2Γ (M0) | 0.4770 vs 0.4769 | within 5% (D11 mapping correct) |
| cl_p wake-free M3 vs conforming | 0.4674 vs 0.4786 | within 5% |
| M3 mesh has a `wake` tag | **False** | topology knows nothing about the wake |
| wake-free Γ vs embedded Γ (α = 4) | **0.2339 vs 0.2384** (1.9%) | within 5% (generic cuts reproduce it) |

The **gate** itself is the compressible one (`tests/test_b3_lifting.py`, 6 passed):
at M0.5, α = 2°, cl_KJ = **0.2828** (medium) sits **inside** the committed
[PG 0.2788, KT 0.2919] bracket read from `cases/reference_data/naca0012_m05/cl_reference.csv`
— the same file the conforming G3.2 gate reads — on **both** mesh families. Same-mesh
A/B vs conforming: Γ within **0.1–0.7%** (0.1177/0.1191/0.1197 vs 0.1175/0.1200/0.1202
on coarse/medium/fine).

**Figures.**

**1. The lift, on both mesh families.** Speed (own-side) and the perturbation
potential drawn **per element** — i.e. exactly as the multivalued DOFs store it. At
α = 4 a crisp branch cut carries [φ] = Γ; at α = 0 the field is flat with **no jump
at all**. The M3 panel exposes the coarser wake-free triangulation that the level
set cuts through generically.

![B3 flow field, lift vs no-lift, M0 embedded](../cases/demo/b3_levelset_lifting/results/flowfield_lift_vs_nolift_m0.png)
![B3 flow field, lift vs no-lift, M3 wake-free](../cases/demo/b3_levelset_lifting/results/flowfield_lift_vs_nolift_m3.png)

**2. How the jump survives to the far field.** LEFT: the nodal [φ] at every cut node
vs downstream distance is **flat at Γ from the TE (d = 0) out to the far field
(d ≈ 15 c)** — the g₁+g₂ wake LS convects it unchanged, and the far-field aux DOFs
are left **FREE** so it exits rather than being drained. *(Pinning them to the
vortex's lower branch was measured to decay the jump 0.0147 → 0.001 — i.e. to drain
the circulation. This is a load-bearing fix, not a detail.)* RIGHT: the **storage** —
the MAIN dof holds the node's own-side value, the AUX dof the other side, and the gap
between them is exactly Γ, all the way out.

![B3 wake jump convection and storage](../cases/demo/b3_levelset_lifting/results/wake_jump_m0.png)

**3. Surface Cp on both families**, using the **D11 per-side DOF mapping** (solid =
M0 embedded, dotted = M3 wake-free, grey dashed = conforming; Cp axis inverted,
suction up). Lower-surface TE triangles *must* read the TE's AUX value — reading
`phi_main` alone gives cl_p = **−3.35**, junk. At α = 0 upper and lower collapse; the
M3 cl_p lands within 2.3% of conforming despite being a coarser, wake-free mesh.

![B3 wall Cp, both mesh families](../cases/demo/b3_levelset_lifting/results/wall_cp.png)

**4. The dual-mesh rule made visible** — the same level-set path on the
wake-**embedded** M0 mesh (which *has* a `wake` tag, its wake nodes lying exactly on
the sheet) and on the wake-**free** M3 mesh (**no `wake` tag anywhere**, generic cuts
through generic elements). Γ agrees to **1.9%**. This is the payoff: **lift on a mesh
that never had a wake embedded**, where no conforming counterpart exists at all.

![B3 dual-mesh: embedded vs wake-free](../cases/demo/b3_levelset_lifting/results/dual_mesh_embedded_vs_free.png)

---

### B4 — TE control volume / implicit-Kutta re-derivation (NEW + closed 2026-07-12)

**Evidence:** `tests/test_b4_te_control_volume.py` — **8 passed** (~29 s); the
B3 demo above is shared. B4 was **inserted** into the track mid-flight because B3's
emergent Γ converged to the **wrong value** — and the reason turned out to be
structural, not a bug.

**★ The finding: the wake LS CANNOT pin Γ.** Its residual is **identically zero for
any spatially-constant jump**, because Σ_c ∇N_c = ∇(Σ_c N_c) = ∇(1) = 0 (partition of
unity) — measured **1.9e−16**. Therefore design_track_b.md §2.3/D2's claim that "the
g₂ on the TE-adjacent wake element **is** the discrete Kutta condition" is **FALSE and
retired**. (Cross-checked against the source: the López dissertation has no explicit
Kutta condition anywhere — the word never appears in its method chapter.) **Γ needs
its own equation.**

Without one, Γ was being pinned by a single *wrong* equation — the TE aux row
(lower-side mass conservation), whose control volume is up/down **asymmetric** on a
symmetric airfoil (the TE fan is 9 upper / 6 lower / 3 cut, because the ε shift sends
every on-sheet node "+"). It over-circulated by **+42%**, *mesh-convergently* — which
is the signature of a method defect rather than a discretization error.

**★ The fix: the nonlinear TE pressure-equality (Bernoulli) Kutta.** Symmetrizing the
control volume is **not available** — the mesh is naturally asymmetric at the TE
(user-arbitrated) — so the condition is instead a **pointwise physical statement that
needs no symmetry**:

> |q_u|² = |q_l|², factorized **exactly** as (q_u + q_l)·(q_u − q_l) = 0, and
> linearized by freezing the mean s̄ = q_u + q_l at the previous iterate.

That yields a row **linear in φ**, re-linearized once per Picard outer (the same
cadence as the density lag — **no new outer loop**), converging to the exact nonlinear
condition. It replaces the TE **aux** row; the displaced lower-side mass-conservation
entries are re-routed onto the TE **main** row, which then carries the **total**
(upper + lower) balance — so mass stays conserved and no side is arbitrarily robbed of
its equation. **Why it is non-degenerate where g₂ is not:** q_u and q_l are recovered on
**different element sets**, so q_u − q_l is *not* a jump gradient and does not vanish for
a constant jump.

**★ The control volumes must be WALL-ADJACENT.** The Kutta condition is about *surface*
velocities, so q must be recovered on the elements carrying a **wall face** (the
upper/lower body surface at the TE), **not** the whole element fan. This is the single
most consequential detail in B4, and it is measurable:

| Γ recovery (α = 2°, incompressible) | coarse | medium | fine | vs conforming |
|---|---:|---:|---:|---:|
| conforming reference | 0.1175 | 0.1200 | 0.1202 | — |
| old B3 `te_kutta="mass"` row | 0.2074 | 0.1760 | 0.1704 | **+42%** (wrong equation) |
| pressure Kutta, **full-fan** control volume | 0.1407 | 0.1355 | 0.1329 | **+11–15%** (interior + wake elements pollute the average) |
| pressure Kutta, **wall-adjacent** ✓ | **0.1177** | **0.1191** | **0.1197** | **< 1%** |

The wall-adjacent control volumes are the highlighted elements in the level-set region
figure — which also shows *where* the level set acts at all: **one** layer of elements
(4.8% of the tets) plus the below-TE fan. The mesh is never modified. Note the cut layer
sits just **below** the sheet: the ε side-shift sends on-sheet nodes "+", so the sheet
effectively lies at y = −ε. **That bias is exactly what B4's TE condition had to be made
immune to.**

![B4 level-set region and TE control volumes](../cases/demo/b3_levelset_lifting/results/levelset_region_m0.png)

**Gate checks.** LS constant-jump null space pinned numerically (1.9e−16 — the reason a
separate TE condition is *structurally* required); TE control volumes verified
wall-adjacent (every element carries a wall face); the below-TE fan is never cut (the
López p.57 ε-shift trap, regression-pinned — the ε shift had been manufacturing
**spurious cuts** there, giving 3 of 6 elements a bogus UPPER copy *below* the wake);
emergent Γ within 5% of conforming (**measured 0.1–0.7%**), while the old
`te_kutta="mass"` row is still >30% out, so the before/after contrast stays honest.

**Consequence: the D2 penalty-Kutta fallback is no longer needed** — this route has **no
penalty weight and no tuning parameter** (s̄ is solved for, not calibrated). Interfaces:
`solve_multivalued_lifting(..., te_kutta="pressure")` (default), with `te_kutta="mass"`
retained purely for the contrast. Derivation in [design_track_b.md §9](design_track_b.md).

---

### B5 — far-field A/B: Dirichlet+vortex vs Neumann outlet (closed 2026-07-12)

**Demo:** `cases/demo/b4p5_farfield/` (directory name predates the renumber) —
`python run_demo.py` redraws and self-checks from the committed `summary.csv`
(**9/9 PASS**); `PYFP3D_B45_RESOLVE=1 python run_demo.py` re-solves the whole study
from scratch (~15 min, threads capped). `tests/test_b45_farfield.py` (**10 passed**,
~20 s) holds the cheap 15c locks.

**Question.** The level-set lifting path needs a far-field BC, and two self-consistent
options exist (design_track_b.md §5.4):

- **option a (vortex)** — spherical Dirichlet freestream **+ a PG point vortex** on the
  far-field MAIN DOFs, with the emergent Γ refreshed into the vortex each outer
  iteration. pyFP3D's compact **15c** domain is calibrated *for* this correction.
- **option b (Neumann)** — the **López** form: inflow Dirichlet freestream (**no
  vortex**), outflow a Neumann outlet carrying the freestream flux ρ∞(u·n̂).

Option b is attractive for the workflow (no Γ-into-far-field feedback, simplest α
sweep), but **with no vortex it truncates the O(Γ/r) far-field tail**, so its domain
must grow — which is why the dissertation uses 10²–10⁷-chord domains.
New interface: `solve_multivalued_lifting(farfield="vortex"|"neumann"|"freestream")`,
default `"vortex"`; helpers `_farfield_split`/`_neumann_outlet_rhs` in
`solve/picard_ls.py`. Conforming path byte-untouched.

**Method — a López-style domain-size re-calibration** (the dissertation §4.1.4 method).
Coarse NACA0012, M0.5, α = 2°, on **both** Track B mesh families, with the far-field
radius swept over R ∈ {15, 30, 60, 120}c. Γ vs R, one panel per family, against the
conforming reference and its ±2% B3 band:

![B5 far-field domain-size study](../cases/demo/b4p5_farfield/results/farfield_domain_study.png)

**Result** (M0 embedded shown; the M3 wake-free family agrees to the third digit):

| R/c | conforming | option a (vortex) | option b (Neumann) | b − a | freestream − a |
|----:|-----:|-----:|-----:|-----:|-----:|
| 15  | 0.1391 | 0.1394 | 0.1337 | **−4.07%** | −7.52% (**diverges**) |
| 30  | 0.1389 | 0.1392 | 0.1364 | −2.01% | −3.87% |
| 60  | 0.1389 | 0.1397 | 0.1383 | −0.99% | −1.96% |
| 120 | 0.1391 | 0.1388 | 0.1381 | −0.50% | −0.99% |

- **Option a is domain-robust:** Γ stays within **0.45%** (M0) / **1.09%** (M3) of the
  truth across the whole 15→120c sweep, and within **0.25%** of the conforming solver
  at 15c.
- **Option b truncates O(Γ/R):** −4.07% at 15c, and the error **halves every time R
  doubles** — the textbook point-vortex far-field decay, visible as a clean straight
  line in the figure. It therefore meets the B3 ±2% band only at **R ≥ ~30c**, and gets
  below 1% only at **R ≥ 60c** — a 2–4× larger domain, ~4× the tets at equal near-body
  h. This is exactly why López needs 10²–10⁷-chord domains.
- **Freestream-Dirichlet** (no vortex, whole boundary) is crudest at every R and
  **diverges** on the compact 15c M0 mesh (M_max 5.9): a lifting body cannot sit in a
  tight box without either the far-field vortex or an outlet.

**★ Verdict — option a stays the DEFAULT.** For pyFP3D's compact 15c workflow the vortex
correction pays for itself. Option b is **validated** as an alternative but is
domain-hungry, so its workflow simplicity does not pay at pyFP3D's scale. Because the
O(Γ/R) truncation is **geometry-universal** (a 3D wing truncates the same
horseshoe-vortex tail), this also decides the far-field default for the M6 B-path — so
the gate did not need its own M6 run to be conclusive.

**The M6 leg is folded into B7** (user-arbitrated 2026-07-12): running the level-set
B-path *solve* on M6 needs the 3D wake-BC machinery that is B7's deliverable — and,
separately measured here, the span-uniform option-a vortex **without** the P5 Γ(z)
taper recreates the **branch-ray artifact** on M6, itself B7 machinery.

---

## B7 — ONERA M6 3D gate (closed 2026-07-12; `cases/demo/b7_onera_m6/`, 35/35 PASS)

**Why the phase exists.** B1–B6 all ran on quasi-2D meshes, where the wake sheet is a
flat strip: no sweep, no tip, no spanwise direction. Three pieces of the level-set path
were therefore *structurally* untested — the **TE-polyline ruled level set** (D9; its
per-segment (v, d̂, n̂) frame is oblique, and B1 already found a real defect there), the
**spanwise clip** `0 ≤ q ≤ span_length` (what makes Γ(tip)=0 *discretely* — the LS
analogue of the conforming free-edge rule), and the **g₂ spanwise-free jump gradient**
(the trailing-vortex DOF). B7 runs the full transonic B-path solve on ONERA M6 and A/Bs
it against the committed P5/P8 conforming baseline.

**Setup.** M∞ 0.84 / α 3.06°, coarse, `farfield="neumann"`, Mach ramp 0.60→0.84 @
dm 0.04, dual-mesh (M1 wake-embedded + M4 wake-free — M4 is within 6–9% of M1's tet
count at equal h_wall, which is what makes the comparison controlled).

| | **M1** embedded | **M4** wake-free | P5 conforming **Picard** | P8 conforming **Newton** |
|---|---|---|---|---|
| cl_KJ | **0.2765** (+2.7% of Newton) | **0.2710** (**+0.7%**) | 0.24788 (**−8.6%**) | **0.2692** (truth) |
| cl_p (3D) | 0.2716 | 0.2656 | 0.24194 | 0.2560 |
| V6 consistency | 1.77% | 1.97% | 2.40% | — |
| shocks η .44/.65/.90 | 0.635/0.588/0.449 | 0.634/0.584/0.454 | 0.596/0.570/0.425 | 0.596/0.541/0.362 |
| Γ root → tip | 0.1076 → **−0.0003** | 0.1055 → **+0.0003** | 0.097 → 0.0206 | — |
| M_max, limited/floored | 1.453†, **0/0** | 1.368, **0/0** | 1.398, 0/0 | 2.13 |
| solve wall time | 22.7 min | 18.4 min | — | — |

> † **Re-read 2026-07-14** (B8-backlog `element_mach2` default flip to
> `mixed_plain="main"`): M1 **1.453 (side) → 1.392 (main)** — the committed
> M_max was itself a beyond-tip mixed-plain artifact cell; honest value sits
> closer to P5's conforming 1.398. M4 and both 2.5D B6 states are bit-identical
> under either reading; all gate bands unchanged. Artifact:
> `cases/demo/b8_tip_taper_ls/results/mmax_reread.csv`.

**★ Finding 1 — the B6 lift inversion reproduces in 3D, first try.** Against the
conforming **Newton** truth (the B6 user-arbitrated baseline), the level-set Picard sits
**+2.7% (M1) / +0.7% (M4)**, while the conforming **Picard** (P5) sits **−8.6%** below
it. Same mechanism as B6's 2D finding: the LS path has **no early-stoppable Γ outer**
(the implicit Kutta makes Γ a *solution mode*, converged with the field), whereas the
conforming Picard's frozen-Γ inner solves plus budgeted per-station secant
under-circulate (the P4-erratum bias; P8 measured it independently at +7.9% for M6
medium). Gating B7's lift on P5 would have **penalised the B path for being closer to
the truth** — hence the Newton anchor. Note the **wake-free workflow mesh (M4) is the
more accurate of the two**, which is the outcome Track B exists to deliver.

**★ Finding 2 — the 3D far field, and why the P5 Γ(z) taper is *structurally
unnecessary* here** (`farfield_decision.png`). The B-path vortex
(`picard_ls._farfield_main`) is a **span-uniform** 2D point vortex whose branch cut is
the ray y=0, x>0 *at every z*. On a 3D wing it misfires two independent ways — both
measured, both appearing as a spurious near-sonic spot at the **outlet, where the sheet
leaves the domain** (max local Mach there, M∞ 0.5):

| far field | outlet M_max | mechanism |
|---|---|---|
| `neumann` (no vortex) | **0.513** | — |
| vortex, sheet re-aimed to y=0 (coplanar) | 0.825 | (b) only |
| vortex, sheet α-aimed (the default) | **0.958** | (a) + (b) |

- **(a) non-coplanarity** — the α-aimed sheet has climbed to y ≈ x·tan α ≈ 0.5 by the
  outlet (x≈10c), far off the vortex's y=0 cut, so the outlet carries a prescribed Γ
  jump **no cut supports**. This is B3's recorded coplanarity rule, in 3D.
- **(b) span-uniformity** — re-aiming coplanar *shrinks but does not remove* it: one
  scalar Γ cannot match Γ(z)→0, and outboard of the tip there is **no cut at all**.
  This is precisely P5's branch-ray artifact, whose conforming fix was the Γ(z) taper.

⇒ **neumann carries no vortex, so neither defect can exist**: the taper is *unnecessary*
on the B path, not merely unimplemented. Price: B5's O(Γ/R) outlet truncation (a few % of
lift on a compact domain) — which is why the gate uses A/B bands, not <1% bands.

**★ Finding 3 — Γ(z) comes out spanwise-smooth, with no smoothing applied** (unplanned — it
became visible the moment the *real* P5 curve was overlaid in `gamma_of_z.png`). Normalised RMS
second difference of Γ(z): **0.0079 (M1) / 0.0091 (M4) vs 0.0970 for the conforming P5 — an
11–12× reduction.** The conforming path solves a **separate secant per TE station**, so its Γ(z)
carries station-to-station jitter — the very defect P5's `INVESTIGATION_gamma_smoothing.md`
chased (concluding that spanwise-Γ *smoothing* moves Γ **away** from the self-consistent value),
and the same machinery whose single-station failure (st133, 32% under-circulated) cost P5 an
entire investigation. The implicit Kutta has **no per-station loop to be noisy in**: Γ is one
solution mode of the coupled system. Track B therefore does not *fix* the P5 spanwise-Γ problem
— it makes it **structurally impossible**.

**★ Finding 4 — the 3D-only machinery needed no new solver code.** B1's oblique-frame and
spanwise-clip fixes held: Γ(z) decays monotonically root→tip and reaches **~3e-4 at the
tip** on both families (`gamma_of_z.png`). The only code gap was post-processing (the TE
node is multivalued, and `section_cp_curve` takes a single nodal field):
`post/surface_ls.py` gained `section_cp_curve_levelset` (D11 per-side plane cut — the
**upper surface is bit-identical** to the `main`-based curve, so every gate shock metric
is unaffected; the lower surface is where D11 bites, and reading `main` there is the junk
that gave B3's cl_pressure = −3.35) and `cl_pressure_3d_levelset` (planform-area
normalisation, pairing with `cl_kj_3d` for V6). Cost came in far under the plan's risk
estimate: the per-outer `spsolve` on ~12k 3D DOFs is ~0.6 s, so a 7-level continuation is
**~20 min, not hours**.

**Honest caveats (recorded, not chased).** (1) **Convergence semantics = the recorded
transonic Picard tail, not `tol_residual`**: the top Mach levels exhaust the 600-outer
budget at |R| ~4–6e-6 (M1 levels 0.72–0.84; M4 0.68/0.76/0.84). The field is **bounded and
physical at every level** (0 limited / 0 floored throughout) and every gate metric is in
band, so the gate asserts *bounded + in-band*, not `converged` — the same P4/B6
engineering-converged regime. (2) **LS Newton on M6 deferred**: `newton_ls.py` uses a
plain `splu`, and P8/N6 measured true-3D LU fill at ~100× the 2.5D cost (it needed
lagged-LU); porting `direct_refactor_every` is the follow-up. (3) Shocks sit 0.02–0.04 c
aft of P5 (in band, and self-consistent with the higher circulation); against the P8
*Newton* shocks the η=0.90 station is 0.087 aft — a like-for-like shock A/B needs the
deferred LS Newton. (4) Coarse only.

**Figures:** `gamma_of_z.png` (Γ(z) both families vs the committed P5 curve; tip→0),
`section_cp.png` (upper/lower Cp at η = 0.44/0.65/0.90), `shock_planform.png` (swept
shock line, forward migration toward the tip), `farfield_decision.png` (the table above).
`p5_gamma_baseline.csv` is committed alongside so the A/B curve reproduces without the
gitignored P5 solution cache. Tests: `tests/test_b7_onera_m6.py` (6 fast + 5 gated).

---

## P13/G13.1 — Tip / wake-edge singularity: characterization (`cases/demo/p13_tip_edge_singularity/`, 10/10, 2026-07-13)

*(This is the evidence for P13/G13.1 — roadmap.md Track P P13, design.md §4.1.
It began 2026-07-12 as the "wake MODEL vs REPRESENTATION" probe below; the
2026-07-13 update added the conforming **fine** third point, the refinement-rate
fit, and the dΓ/dz mechanism, closing G13.1.)*

**★ Rate + mechanism (P13/G13.1, 2026-07-13).** With the conforming **fine** M6
mesh added as a third point (coarse/medium/fine = 55.5k / 350.7k / 2.51M tets),
the tip-box peak Mach goes **0.712 → 0.981 → 1.510**, a log-log exponent
**p = 0.59** (peak ~ h^−p) — squarely in the flat-plate-edge band [0.4, 0.65],
i.e. **1/√r, not 1/r** (a 1/r concentrated line vortex would give p = 1). The
**driver is the trailing vorticity dΓ/dz**, not the bound circulation Γ: Γ → 0
at the tip (a necessary-not-sufficient regularity condition), but the *unloading
rate* |dΓ/dz| is largest at the tip (~10× mid-span on B7's smooth Γ(z)), and a
terminating flat vortex sheet cannot regularize its own free edge — exactly a
flat-plate leading/trailing edge at incidence. Within the same tip box the
**p95/mean stay flat** (0.573 → 0.562 → 0.525 / ~0.49) while the peak diverges —
the localized-edge signature. **★ And the conforming fine M∞0.5 solve does NOT
converge** (limited/floored cells, ~1.4k NaN): the tip singularity trips the
speed limiter / density floor **even subsonically**, so the fine mesh is not a
discrete solution — the exact M∞0.5 analogue of what G9.1 found transonically.
This **corrects the committed "1/r-type" phrasing** (this file above and
roadmap.md) to 1/√r and supplies the dΓ/dz mechanism. Figure `tip_edge_growth.png`
now overlays the measured p, the 1/√r (p=0.5) and 1/r (p=1) guide slopes.

**Original probe (2026-07-12): wake MODEL vs REPRESENTATION.**

**Why.** P9/G9.1 found the M6 transonic solve does not converge under refinement —
the unlimited local Mach at the wake sheet's free tip edge climbs 1.40 → 2.13 → 7.93
(coarse/medium/fine) and 9 cells cross the M_cap=3 limiter on the fine mesh, blocking
the Newton freeze machinery so the fine "solution" is a limit-cycle artifact. P9 called
it a vortex-sheet-edge singularity of the rigid planar wake and pointed at "the
tip/wake treatment (Track B / …)". One doc over-stated that as "precisely what Track B
exists to fix". This demo settles whether the Track B level-set **representation**
removes it, or whether it is a wake-**model** property present on any discretization.

**Method (the point is that it is cheap and clean).** Probe SUBSONICALLY (M∞ 0.5): no
shock, no artificial density, no speed limiter — nothing to confound the pure
potential-flow 1/r edge signal. Measure the peak local Mach in the P9 tip-edge box
(z/b > 0.95, at/just aft of the swept tip TE corner, chord plane) coarse→medium, on the
conforming path (`solve/newton.py`) and the level-set path (`solve/picard_ls.py`, both
the wake-embedded M1 and wake-free M4 meshes). **conforming-M1 vs level-set-M1 use the
IDENTICAL onera_m6 mesh** — a true same-mesh A/B of the representation change. A wing
box (x < x_TE, z/b < 0.95) is the control: the real, bounded flow.

**Result.**

| path | tip-edge M_max (levels) | growth c→m | exponent p | tip-box p95 (c→m) |
|---|---|---|---|---|
| conforming (M1) | 0.712 → 0.981 → **1.510** (fine ✗conv) | ×1.38 | **0.59** (3-pt) | ×0.98 (flat) |
| level-set (M1, same mesh) | 0.672 → 1.532 | ×2.28 | 1.34 (2-pt) | ×0.98 (flat) |
| level-set (M4, wake-free) | 0.661 → 1.151 | ×1.74 | 0.89 (2-pt) | ×0.93 (flat) |

> ★ **ERRATUM 2026-07-14 (B8 re-spec diagnosis).** The two level-set rows were
> read through `element_mach2`'s then-default mixed-plain "side" handling,
> which inflates beyond-tip mixed-side plain cells ×5 (elem 93977: side 1.532
> vs main 0.309). The HONEST LS exponent is **+0.62** (+0.37 excluding the
> straddler sliver) — the same object and magnitude as the conforming +0.52.
> The "LS ≥ conforming" magnitude comparison is RETIRED; the qualitative claim
> (the edge diverges on both paths) stands. See §B8
> (`run_b8_termination_diagnosis.py`, `b8_termination_diagnosis.csv`).

The tip-edge **peak** diverges on **all three** paths while the same-box
**p95/mean stay flat** — the signature of a *localized* edge singularity (only
the few cells at the very corner grow), seen with **zero transonic machinery**, so
it is a genuine potential-flow feature, not a shock/limiter artifact. (The bulk
"wing" interior is flat coarse→medium too; only the fine conforming wing *max* is
polluted by a separate sharp-TE edge cell — another P1 edge feature — so the clean
same-box control plotted is the tip-box p95.) The conforming three-point exponent
**p = 0.59** puts the growth in the **1/√r flat-plate-edge** band, not 1/r. The
level-set representation does **not** remove it: the honest LS exponent (+0.62,
per the erratum above — the raw "1.34 ≥ conforming 0.52" comparison is retired)
is the same object and magnitude as the conforming +0.52, and the LS peak
sits in the **+2.9% straddler cells** at/just beyond the geometric tip (z/b ≈ 1.01),
where the jump terminates mid-element.

**⇒ It is a WAKE-MODEL defect.** Track B changes the wake *representation*, not the
model (B7 keeps the same rigid planar sheet ending at the tip with Γ(tip)→0), so it
does **not** fix G9.1's cause. The model-level fix is wake **roll-up / an explicit tip
vortex**, which **no current Track B phase does** (B9 free-wake is shelved, and is about
O(θ²) deflection rather than roll-up). This corrects the earlier "P9 corroborates the
Track B route" framing. What Track B *does* structurally eliminate is the separate
**secant/Kutta-closure** family (P5 st133, probe-sharing, Γ(z) jitter, branch-ray) — see
the B7 section. Figure `tip_edge_growth.png` (log-log peak Mach vs 1/h, tip edge vs wing
control); heavy solves cache to `results/*.npz` (gitignored, ~20 min to regenerate).

---

## P13/G13.2 — Tip-edge desingularization: the spanwise loading taper (`cases/demo/p13_tip_edge_singularity/run_taper_probe.py` + `run_taper_physics.py`, 2026-07-13)

**The fix.** Taper the accepted circulation toward the tip,
`Γ_eff(z) = F(z)·Γ_Kutta(z)`, on the per-station Kutta target
(`constraints/wake.py::tip_taper_factors`; `solve_newton_lifting(tip_taper=…)`,
default `None` = bit-identical; reaches the transonic driver via `newton_kw`).
**Shipped model (user-arbitrated): `vanish_smooth` (smoothstep, compact
support), r_c = 0.05·b_semi.**

**★ The mechanism is DISCRETE — and it is neither roll-up nor a vortex core**
(both of which this report and design.md §4.1 previously proposed). The solver
never sees the continuum edge: it sees the **outermost TE station**, which
retains `Γ_last` and sheds it as a *concentrated vortex over the last cell*
(free-edge nodes are single-valued, so the jump falls to 0 in one element),
inducing `~Γ_last/h`. With `Γ_last ~ h^q`:

> **edge peak ~ h^(q−1)  ⟹  p ≈ 1 − q,  criterion q ≥ 1**

This predicts the baseline *exactly*: Γ~√u with u_last~h gives Γ_last~√h, i.e.
q = 0.44 measured ⟹ p_pred 0.56 against **p_meas 0.52**.

**★ The taper is amplified, not applied.** Γ is a fixed point of
`Γ = F·Γ_Kutta(Γ)`, and the Kutta map has slope b≈0.93 (P2), so
`Γ/Γ* = F(1−b)/(1−F·b)` — a taper of 0.8 yields **0.21×**, not 0.8× (test-locked).
This is why the measured q ≈ 3.3 far exceeds the naive exponent, and why r_c must
stay small.

**Result (M∞0.5, strict OFF-BODY edge box).**

| | coarse | medium | fine | p |
|---|---|---|---|---|
| untapered | 0.712 | 0.981 | 1.510 | **+0.592** |
| tapered | 0.567 | 0.565 | 0.570 | **+0.009** |

**★ The M6 fine mesh is now a GENUINE DISCRETE SOLUTION** (converged, 0 limited /
0 floored / 0 NaN) — exactly what P9/G9.1 could not achieve. Transonic **M0.84**:
coarse/medium converge 0-limited and the medium `M_max` drops **2.13 → 1.725**
(P5's "bounded tip-TE-corner overshoot" was the same object).

**The price, and that it is LOCAL** (`run_taper_physics.py`). cl_KJ falls
**−1.1…−1.6 %** (scaling with r_c: −0.70/−1.58/−3.27 % at r_c = 0.03/0.05/0.08 b).
But Γ(η), cl(η) and the sectional Cp are **unchanged inboard of η≈0.95** (inboard
circulation −0.51 %), and TE pressure closure at η=0.90 stays at the baseline
value (0.232 vs 0.218). The taper even makes cl *more* mesh-convergent than the
untapered baseline (+0.2 % vs +0.7 % coarse→medium — the untapered case is still
gaining spurious tip lift).

**★ The tanh form is disqualified — but not for the reason first argued.** It
*does* regularize (q ≈ 1.00, exactly the marginal case its s=0 predicts). It is
rejected for **unbounded support**: tanh never reaches 1, so it depresses F over
**57 of 83 stations** (inboard to η=0.77), costing **−7.4 % lift**, **−4.9 %
inboard circulation**, and **breaking TE pressure closure at η=0.90** (gap 0.972
vs baseline 0.218) — where there is no singularity to fix. A tip model must be
*local*; the tanh silently re-rigs the whole wing.

**★ A metric trap, found the hard way.** G13.1's tip box (z/b>0.95, dx ≥ −0.05)
admits WING cells. Once the edge is regularized, the box's max **migrates onto the
ordinary wing suction peak** at z/b≈0.95 and stops measuring the singularity —
making a *working* fix look like it made p *worse*. The edge must be measured
off-body (dx>0, z/b>0.98).

**Two pre-existing artifacts, confirmed taper-independent** (both user-flagged on
first sight of the figures): the Cp **sawtooth** is the **P6/G6.1 wall-gradient
recovery** artifact (`smooth_passes`, default 0 — slope reversals 40→2 when
enabled; the *raw* count is identical for all three tip models); the **Γ(η)
jitter** is conforming **Kutta-probe placement** on the unstructured swept TE (P5
known item; RMS d² 0.042 vs P5's 0.097 — the level-set path is 11–12× smoother by
construction, B7); and the **baseline TE Cp gap of 0.14–0.22** is the conforming
**potential-jump** Kutta (design.md 4.4) being only an approximation of true
pressure equality, plus the sharp-TE P1 floor — which is precisely why Track B's
B4 had to introduce the explicit nonlinear |q_u|²=|q_l|² Kutta.

**Open: the level-set clause.** Not a mechanical port — the LS path has **no Γ
DOF** and its TE row is **homogeneous** (`s·(q_u−q_l)=0`), so scaling it by F is a
**no-op**. The clean analogue blends the pressure-equality row with B2's
continuity weld (`F·K̂ + (1−F)·Ŵ = 0`), which is a *different model* needing its
own r_c calibration. Designed, not implemented.

**★ Honest gap — G13.3 is still blocked, by further causes.** All three M6 meshes
are now discrete solutions, so a Richardson is *mechanically* possible — but the
lift sequence is **not in the asymptotic range**: cl_KJ 0.2001 → 0.2005 → 0.2121,
i.e. increments **+0.2 % then +5.8 %**, which *grow*. **No extrapolation was run**
(P9/G9.3 discipline: report `n/a`, never fabricate). Removing the tip singularity
**revealed** the remaining problems rather than curing them — see the G13.3
section below, which localizes them (and **retracts** this section's first
reading, that the residual growth was "mid-span Γ ⇒ wing/LE resolution": the wing
interior turns out to be *converged*).

Demos: `run_taper_probe.py` (**16 passed + 2 xfailed**) and `run_taper_physics.py`
(**10 passed + 1 xfailed**). The xfails are the *documented negatives*, and they
are the point: `tanh_half` sits exactly on the criterion's knife edge (q = 1.00)
**and** is disqualified for unbounded support (it unloads the wing to η = 0.77);
`vanish_smooth` at r_c = 0.03 is under-regularized (p = 0.207 > 0.20) — which is
*why* the shipped r_c is 0.05. Tests: `tests/test_p13_tip_taper.py` (15).

---

## P13/G13.3 — The third singularity: a flat tip cap (`cases/demo/p13_tip_edge_singularity/run_g133_ladder.py`, 5/5, 2026-07-13)

G13.2 fixed the wake's free tip edge, and Track M (M1b) fixed the mesh ladder —
and the M6 lift sequence *still* is not asymptotic. Growing increments under
**uniform** refinement are the signature of a singularity still being resolved,
so a third object had to be there. A three-region box study on the self-similar
ladder `RICHARDSON_LADDER = (coarse_ss, medium, fine)` finds it, and it is **on
the wall, not in the wake**:

| region | coarse_ss → medium → fine | p | verdict |
|---|---|---|---|
| **tip-cap edge (WALL)** | 0.662 → 0.824 → **1.015** | **+0.321** | **DIVERGES** |
| wake free edge (the G13.2 fix) | 0.536 → 0.565 → 0.570 | +0.045 | bounded |
| wing p99 (control) | 0.642 → 0.628 → 0.629 | −0.014 | bounded |

The control matters: the ordinary wing field is **mesh-converged**, so this is a
*localized* divergence (a singularity), not "the whole solution is under-resolved"
— and it retracts the earlier "mid-span Γ / wing-LE resolution" reading. The
G13.2 taper also **holds** across the full ladder.

**Root cause: a documented, deliberate geometry simplification.**
`meshgen/wing3d.py` builds a **flat** tip cap where the real ONERA M6 has a
**rounded** one. A flat cap meets the upper and lower surfaces at a sharp convex
edge — in potential flow, an edge singularity of exactly the kind P13 exists to
remove, only on the **body** instead of the wake. **This is not a P11 (curved wall
*element*) problem:** isoparametric elements cannot regularize a genuinely sharp
geometric *edge*. The geometry itself is wrong, so the fix belongs to **Track M**
(round the cap), and neither P11 nor the level-set port is required.

**⇒ Three distinct defects blocked 3D grid convergence, and they were different
objects:** (1) the wake free tip edge (p = 0.59) → **fixed** by the G13.2 taper;
(2) the `h_far` mesh-ladder clamp → **fixed** by Track M M1b; (3) the flat tip-cap
wall edge (p = +0.32) → **fixed** by Track M **M5** — the geometry is corrected
(seam crease q = −0.92 vs the flat cap's −0.00) AND the flow follows: on the round
ladder the cap-surface exponent drops to +0.09 (bounded) and the lift Richardson
becomes definable (p = 2.31, cl_KJ→0.2050). See the P13/G13.3 (subsonic) section
below. All three defects are addressed at M0.5. **At M0.84 the story is
different** — the rounded tip lets flow wrap around and accelerate, which
*amplifies* the (pre-existing, shared) sharp tip-TE singularity until the fine
mesh's Mach ramp dies at M = 0.75 and never reaches M0.84. A sub-problem at the
trailing edge, not a defect of the cap; the P9 band verdict is not earned on
either geometry. See the P13/G13.3 (transonic) section below.

**✓ Evidence status (audit 2026-07-13; RESTORED same day).** This report's G13.2
transonic claim (cl_KJ 0.2593 → 0.2652 → **0.2866** at M∞0.84, M_max 2.818, 0
limited/floored) was found **prose only** — no committed script, CSV or cached
solve, a repo-wide search finding those numbers in the `.md` files and nowhere
else, while a P11 ledger status had been changed on its strength. It was **re-run
from scratch and REPRODUCES to 4 digits**, and now has a committed artifact:
demo `cases/demo/p13_tip_edge_singularity/run_g132_transonic.py` (5/5 PASS),
`results/g132_transonic.csv`:

| level | n_tets | cl_KJ | cl_p | M_max | over M_cap | lim | flr | conv | wall_s |
|---|---|---|---|---|---|---|---|---|---|
| coarse | 55 531 | 0.2593 | 0.2534 | 1.394 | 0 | 0 | 0 | ✓ | 9 |
| medium | 350 718 | 0.2652 | 0.2608 | 1.725 | 0 | 0 | 0 | ✓ | 149 |
| fine | 2 513 255 | 0.2866 | 0.2835 | 2.818 | 0 | 0 | 0 | ✓ | 2679 |

All three are genuine discrete solutions; the census is G9.1's own (unlimited
Mach field + M_cap count), so it is a strict A/B — the only change is `tip_taper`.
⇒ "the flat-cap fine mesh is a discrete solution" and "cl_KJ reaches 0.2866" are
real, not prose. **Two honest limits remain, and the M5 work above resolves the
first:** (a) 0.2866 is on the **FLAT-cap** sequence, which is not asymptotic (the
flat tip cap diverges — G13.3), so it is a REPORTED single point, not a Richardson
value; the definitive P9-band verdict needs the M0.84 Richardson on the **round**
ladder. (b) The rerun exposed the recipe trap that lost the original evidence: the
fine mesh (~450k dofs) needs `precond="amg"` + tight EW forcing (η=1e-8) +
`m_start=0.30, n_picard_seed=12`; the medium recipe's `precond="direct"` is P9's
4h39m/26GB splu trap (killed at 1h16m/24GB on the first attempt). The docs' "38
min" was in the right ballpark — the amg rerun took 44.6 min at RSS 3.9 GB.

---

## Track M / M5 — the tip cap, rounded (`cases/demo/m5_round_tip/`, 9/9 PASS, 2026-07-13)

G13.3 named the defect; this is the fix. `meshgen/wing3d.py` grew
**`tip_cap="round"`** (default `"flat"` ⇒ every existing family bit-identical, so
the P5 / P8-G8.2 / B7 / M1 locks are untouched), which closes the wing with the
**half body of revolution swept by the tip section about its own chord line** —
`{√(y² + (z−B_SEMI)²) ≤ t(x)}`, with `t(x)` the tip section's local half
thickness. It is an OCC revolve of the tip section's upper half-face about an
edge *of* that face (the degenerate-at-the-axis kind, like a sphere from a half
disc), fused onto the loft.

★ **Why this construction and not a fillet, a loft, or curved elements.** The cap
radius **vanishes at the LE and the TE**, so the cap degenerates to a *point* at
each — which means the **TE line, the wake sheet, the tip TE corner, the Kutta
stations and B_SEMI are all unchanged**. Only the tip *wall* moves. That is what
makes the A/B against M1 controlled, and it is why the fix needed **no solver
code at all**. (A constant-radius fillet is geometrically impossible here — the
section thickness goes to zero at the TE — and P11's isoparametric *elements*
cannot regularize a sharp *edge* in the first place.) The revolved solid was
verified analytically and asserted at generation to protrude past the wing by at
most **4 × 10⁻⁶ m**, three orders below the finest edge size, and never to reach
aft of the local TE — so it can neither cut the wake sheet nor spawn slivers.

★ **The gate metric is the SEAM CREASE ANGLE, because the solver never sees the
CAD.** It only ever sees the triangulation, so "the geometry is round" is not the
claim that has to be true — "the *mesh* has no edge" is. Measure the turning
angle between the outward normals of adjacent wall triangles across the
tip-section seam at `z = B_SEMI` (the locus that *is* the sharp edge in M1), away
from the LE and TE — both sharp **by design** in either family, the TE because it
carries the Kutta condition (`post/surface.py::wall_crease_angles`):

| seam crease (p99) | coarse | medium | fine | exponent q |
|---|---|---|---|---|
| **flat cap (M1)** | **91.9°** | **92.1°** | **92.1°** | **−0.00** — does not move |
| **round cap (M5)** | 46.8° | 24.5° | **13.7°** | **−0.92** — O(h) |

(round cap, *max* rather than p99: 46.8 → 25.0 → 18.1°, q = −0.68 — see below.)

That contrast is the whole gate. A sharp convex edge creases by its own turning
angle: halving `h` **resolves** it better and removes nothing, which is exactly
why refinement made the flow singularity *worse* (p = +0.321). Facets
approximating a **smooth** surface crease by O(h · curvature) and halve when `h`
halves — the discrete statement of "there is no edge in the limit".

Two honesty notes on the metric. (i) The flat cap's **median** seam crease is
**0.00°**: the rest of the seam is already smooth, so this metric is reading the
edge and nothing else — it is not a diffuse mesh-quality number. (ii) The round
cap's **max** decays more slowly (q = −0.68) than its p99 (−0.92), because the
max is a single worst facet at the *thin* end of the seam window, where the cap
radius is smallest and the local facet size is set by `h_edge` (the TE/LE
refinement) rather than by `h_tip`. Both **decay**, which is the claim; the flat
cap's decays not at all.

**★ `h_tip = 0.25 · h_wall` is load-bearing, not a tuning knob.** The cap radius
is only `TIP_CAP_RADIUS` = 22 mm (1.9 % of the semi-span). At `h_wall` the coarse
cap would be about *one element* wide — i.e. the mesh would quietly discretize the
rounded cap back into a flat one and the geometry change would do nothing. `h_tip`
scales with the level, so the cap costs the same *fraction* at every level and the
refinement ray is preserved.

Other measured items (demo `run_demo.py`, 9/9; `tests/test_m5_round_tip.py`, 19):
cap resolved to its apex (wall `z_max` 1.218465 / 1.218466 vs the analytic apex
1.218467); quality inside the M1 bounds (min dihedral 4.05° / 5.18°, max aspect
16.2 / 6.4); `cut_wake` keeps the M1 semantics (85 TE stations; the tip TE corner
is still a free-edge node ⇒ Γ(tip) = 0 discretely); G2.1 freestream on the cut
mesh < 1e-10; **tip TE corner offset exactly 0.0** and the wake sheet still in the
chord plane, ending at B_SEMI. **Cost ×1.29 / ×1.28 tets** (59,359 / 448,197 vs
M1's `coarse_ss` 46,067 / `medium` 350,718) — level-independent, as a self-similar
ladder requires. Note the comparison is against M1's `coarse_ss`, not its shipped
`coarse`: the latter still carries the M1b `h_far` clamp and would report a
spurious ×1.07.

The **flow** consequence is P13/G13.3's, measured on the round ladder below.

---

## P13/G13.3 (subsonic) — rounding the tip cap restores 3D grid convergence (`cases/demo/p13_tip_edge_singularity/run_g133_roundtip.py`, 9/9, 2026-07-13)

The M5 geometry fix, tested where it matters: a strict A/B of the box study and
the three-point Richardson on the round ladder vs the flat one, at M∞0.5/α3.06
(subsonic, so the edge signal is geometric — no limiter/shock in the way), tip
taper on both, only `tip_cap` differing. All six levels converged, 0 limited /
0 floored.

**The cap edge singularity is gone.** Measuring the fluid just outboard of the
tip (`z > B_SEMI`), the peak-Mach exponent falls from **+0.327 (flat)** to
**+0.091 (round, bounded)**; the tip region with the design-sharp LE/TE excluded
(chord-frac 0.05–0.90) is **converged, p = −0.006**. The wake free edge stays
bounded (+0.071, the G13.2 fix holds) and the wing interior stays converged
(−0.013).

**★ The three-point Richardson G9.1 could never run is now earned.** Round-cap
cl_KJ **0.2159 → 0.2073 → 0.2055**, increments **−3.95% then −0.88% (shrinking)**
⇒ the sequence is asymptotic; observed order **p = 2.31**, extrapolated
**cl_KJ(h→0) = 0.2050**. The flat ladder's cl was **non-monotone**
(0.2015 → 0.2005 → 0.2121) — no Richardson is definable from it, which is exactly
G9.1's failure reproduced and now removed.

| region | flat (M1) | round (M5) | verdict |
|---|---|---|---|
| cap surface (`z > B_SEMI`) | p = +0.327 | **p = +0.091** | bounded |
| tip, LE/TE excluded | — | **p = −0.006** | converged |
| wake free edge (G13.2) | +0.045 | +0.071 | bounded |
| wing interior (control) | −0.014 | −0.013 | converged |
| cl_KJ Richardson | non-monotone → n/a | **p = 2.31, cl→0.2050** | earned |

**★ Honest caveat — the metric trap (G13.1 finding 6).** The *broad* G13.1 tip
box `(z/b>0.98) & (dx<0)` still shows a divergence (p = +0.38), but its maximum
has **migrated**: on the fine round mesh it sits at **chord-frac 0.999 — the
zero-thickness trailing edge**, which is sharp *by design* (it carries the Kutta
condition), present in *both* families, and is not something any tip-cap change
removes or should. It is a local, integrable feature and does not spoil the
integrated lift — which is why the lift Richardson is clean. Once the edge you
fixed is gone, a broad max-in-box metric latches onto the next-sharpest feature;
the cap-surface and TE-excluded boxes above are the honest measures of the fix.

**Scope.** This is the subsonic (M0.5) leg — it proves the *mechanism* (the
geometry fix restores grid convergence and the Richardson is now definable).
Firing P9's pre-registered decision bands (cl_KJ∞ ≥ 0.283 resolution / ≤ 0.278
floor) is a *transonic* question (M0.84), run next.

---

## P13/G13.3 (transonic) — the round cap AMPLIFIES the tip-TE singularity, and the fine mesh never reaches M0.84 (`cases/demo/p13_tip_edge_singularity/run_g133_roundtip_transonic.py`, 5/5 + `..._locate.py`, 2026-07-13/14)

The M0.84/α3.06 three-point Richardson on the round ladder, same recipe as the
committed flat-cap run (`run_g132_transonic.py`: tip taper, `precond="amg"` +
tight EW forcing + the `m_start`/Picard fine guards). It is a **negative result**,
and the demo's checks assert it truthfully.

| level | n_tets | ramp reached | last converged | M_max | over M_cap | cl_KJ (flat-norm) |
|---|---|---|---|---|---|---|
| coarse | 59,359 | **M0.84** ✓ | 0.8400 | 1.51 | 0 | 0.2769 |
| medium | 448,197 | **M0.84** ✓ | 0.8400 | 2.00 | 0 | 0.2763 |
| **fine** | 3,257,273 | **M0.75** ✗ | 0.7375 | n/a | n/a | **n/a** |

**★ The finding is sharper than "the fine did not converge": the fine mesh never
reaches M0.84 at all.** Its Mach continuation **breaks down at M = 0.75** — it
dm-halves until dm < dm_min and gives up (one cell in the density floor, residual
stalling ~8e-6) — so it never gets to M0.80, let alone M0.84. **There is no M0.84
fine state to census.** Only two of three points exist, so there is no three-point
Richardson and no P9 band verdict.

**★ The site is the SHARP TIP TE — and the round cap does not create it, it heats
it.** The committed **flat** M0.84 run *completes* the ramp at the same refinement
level and converges (M_max 2.818, 0 over cap, cl_KJ 0.2866), which at first seems
to exonerate the shared trailing edge — if the TE were the culprit, the flat cap
would die too. It does not. But the field-saving rerun
(`run_g133_roundtip_transonic_locate.py`) settles it: of the **20 fastest cells on
the failed fine field, 20 are on the sharp tip TE** (z/b 0.97–0.99, chord-frac
≈ 0.998) and **none is on the new cap surface** (z/b > 1). The rounded tip lets
flow **wrap around the tip and accelerate**, which raises the velocity at that
same, pre-existing sharp TE at *every* level (M_max 1.51 / 2.00 round vs
1.39 / 1.73 flat) — and at fine resolution that pushes it past the limiter, where
the flat cap's cooler tip flow stays under. **⇒ the cap did not add a
singularity; it amplified the tip-TE one that was always there.** This also
*confirms* the subsonic box study's site (chord-frac 0.999), rather than
contradicting it.

**⚠ Method erratum, self-caught — three numbers retracted.**
`solve_newton_transonic` returns the **failed level's state** (`converged=False`,
at *that* level's Mach, not at `m_inf`). An earlier version of this demo censused
that state at `M_INF = 0.84` anyway — applying the **wrong freestream Mach** to a
M ≈ 0.75 velocity field — and so reported a spurious **"M_max 3.97 / 5 cells over
M_cap / cl_KJ 0.2415 at M0.84"**. **Those three numbers are retracted.** `solve()`
now records `target_reached` / `m_final` / `m_last_converged` and **refuses to
census a state whose ramp did not reach the target**; every census column reads
`n/a` for the fine level in `g133rt_transonic.csv`, and `get()` rejects any cache
that predates the guard.

**Net G13.3 picture.** Rounding the cap is a clean **subsonic** fix (Richardson
p = 2.31). Its **transonic** cost is that the amplified tip-TE flow becomes
unsolvable at fine resolution — a sub-problem at the *trailing edge*, not a defect
of the cap. And the P9 transonic band verdict is **not earned on either
geometry**: the flat cap's solves converge but its sequence *rises*
(0.2593 → 0.2652 → 0.2866, non-asymptotic, on the clamped ladder and polluted by
the flat cap edge), so 0.2866 is a single reported point, not a Richardson
extrapolation; the round cap has no third point at all. The
"0.019 gap = resolution ⇒ P11 refuted" conclusion therefore has **no clean
asymptotic discrete-solution basis on either geometry**.

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

---

## Track B / B8 — level-set tip-edge desingularization (row-blend tip taper)

**Status: ◐ IMPLEMENTED, GATE NOT MET — a NEGATIVE result, and the reason is the
finding.** Demo `cases/demo/b8_tip_taper_ls/run_b8_taper_ls.py` (**12/12 PASS** —
every check asserts a *measured fact*, including the ones that record the
failure). Artifacts: `results/b8_taper_ls.csv`, `results/b8_taper_ls.png`,
`results/checks_b8.csv`. Tests `tests/test_b8_tip_taper_ls.py` (13).
Case: ONERA M6 coarse+medium, M∞ 0.5, α 3.06°, `upwind_c = 0` (no limiter, no
shock — the clean geometric probe, as G13.1 established).

### What was built

P13/G13.2 killed the tip/wake-edge singularity on the **conforming** path with a
spanwise loading taper `Γ_eff(z) = F(z)·Γ_Kutta(z)`. That cannot be ported to the
level-set path: it has **no Γ DOF** (the implicit Kutta makes Γ a solution mode)
and its TE Kutta row `s·(q_u − q_l) = 0` is **homogeneous** (RHS ≡ 0), so scaling
it by F is an algebraic **no-op** (G13.2 finding (8)). B8 implements finding (8)'s
prescribed analogue — a convex **blend** of the pressure-equality row with B2's
continuity weld, per TE node:

```
F_i · [ s·(q_u − q_l) ]  +  (1 − F_i) · [ φ_aux − φ_main ]  =  0
```

`F=1` inboard ⇒ pure pressure Kutta (**bit-identical**); `F=0` at the tip ⇒ weld
⇒ jump = 0 ⇒ tip unloaded. Shipped in `kernels/cut_assembly.py`
(`te_kutta_coo(weights=F)` + new `te_weld_coo`), `wake/multivalued.py`
(`assemble_matrix(tip_taper=…)`), and threaded through
`solve_multivalued_lifting` / `_transonic` / `solve_multivalued_newton` (blended
residual **and** Jacobian). Default `None` ⇒ every existing B3–B7 result is
bitwise unchanged (B-suite **59 passed / 0 failed**).

### The measured result

| variant | edge peak coarse→medium | **p** | **q** (Γ_last ~ h^q) | Δ inboard Γ | Δ cl_KJ |
|---|---|---|---|---|---|
| untapered | 0.672 → 1.532 | **+1.341** | 0.44 | — | — |
| `vanish_smooth` r_c=0.03 | 0.675 → 1.564 | +1.37 | | | |
| `vanish_smooth` r_c=0.05 | 0.681 → 1.619 | +1.41 | | | |
| `vanish_smooth` r_c=0.08 | 0.702 → 1.856 | +1.58 | | | |
| `vanish_linear` r_c=0.05 | 0.678 → 1.569 | **+1.367** | **4.73** | **+0.01%** | **+0.03%** |

The blend **works exactly as its model predicts**: it converges cleanly (0
limited / 0 floored at every r_c), it **unloads the tip circulation far past the
conforming criterion** (q = 4.73, criterion q ≥ 1), and it is **perfectly local**
(inboard Γ +0.01%, cl_KJ +0.03%). **And the tip-edge peak still diverges under
every taper** — larger r_c is *worse*.

The untapered **p = +1.341 reproduces G13.1's level-set exponent (1.34)** on the
same meshes, so the off-body metric (dx > 0, z/b > 0.98 — G13.2 finding (6)'s
trap-free box) is measuring the right object.

### Why — three findings

1. **G13.2's DISCRETE mechanism does not transfer.** There, `p ≈ 1 − q`: the
   outermost TE station sheds its retained `Γ_last` as a concentrated vortex over
   the last cell, so the edge velocity ~ `Γ_last/h`, and `q ≥ 1` kills it. Here
   **q = 4.73 yet p = +1.37** — nowhere near `1 − q = −3.73`. **Killing `Γ_last`
   does not kill the peak.**

2. **The lift cost is ~0 (+0.03%, vs the conforming taper's −1.74%) because there
   is nothing to unload.** The level-set **implicit Kutta already drives
   Γ(tip) → 0 emergently** (B7 measured ±3e-4). The conforming path *needs* the
   taper only because its free-edge rule leaves `Γ_last ~ √h` (q = 0.44). **The
   level-set path never had that disease.**

3. **★ MECHANISM — where the peak actually lives.** The peak cell is **outboard
   of the geometric tip** (z/b = **1.0118**, dx = +0.061); it is a **`beyond_tip`
   element** — one the **spanwise clip refuses to cut** (`cut_elements.py`: a
   crossing needs `q ≤ span_length`); it is the **same element tapered or not**
   (elem 93977); and it is **not a small-cut sliver** (volume **0.71×** the
   median, and not even a cut element — so the CutFEM small-cut instability is
   ruled out). ⇒ **the level-set tip singularity lives in how the embedded sheet
   TERMINATES, not in the circulation it sheds.**

### Two-path physics A/B (the decisive measurement)

Read from the committed conforming numbers
(`cases/demo/p13_tip_edge_singularity/results/taper_probe.csv` — **same meshes,
same M0.5/α3.06**; read, never recomputed, per the cost rule):

| | conforming `Γ = F·Γ_Kutta` | level-set row blend |
|---|---|---|
| untapered | p = **+0.521** (diverges) | p = **+1.341** (diverges) |
| `vanish_linear` r_c=0.05 | p = **+0.103** (**bounded**) | p = **+1.367** (**diverges**) |
| lift cost | **−1.74%** | **+0.03%** |

**The same F(z), the same meshes, the same condition: the conforming taper bounds
its edge; the level-set row blend does not.**

### Verdict

**The two paths' tip singularities are DIFFERENT OBJECTS.** Finding (8)'s "clean
analogue" is a faithful analogue of the conforming *model* — but the level-set
path does not have the conforming path's disease, so the analogue treats a
patient that is not ill. **The shipped machinery is correct, tested, and
bit-identical by default; it is simply not the cure.** **B8 needs a re-spec aimed
at the sheet TERMINATION** (the spanwise clip / beyond-tip zone) — candidate
directions: a graded/faded sheet termination, ghost-penalty-style stabilization
of the clip boundary, or extending the sheet beyond the tip. **User arbitration
required before re-speccing.**

### Cost boundary (recorded)

The level-set path has **no `precond` option** — `solve_multivalued_lifting` is
hardcoded to sparse-direct `spsolve` (a deliberate B2 decision to decouple
"is the assembly correct" from preconditioner convergence; GMRES + AMG is the
deferred B3+ scaling path). M6 **medium costs ~484 s / solve** at 67,426 extended
DOFs (~1.2 GB RSS). **M6 fine (~450k DOFs) on the level-set path would hit the
same splu wall P9 hit on the conforming Newton — with no AMG escape hatch.** So
this demo is coarse+medium **by necessity**, and any fine-mesh level-set work
needs the deferred GMRES+AMG path first.

## Track B / B8 re-spec — termination diagnosis + span blend (CLOSED as CHARACTERIZED-NOT-CURED, user-arbitrated 2026-07-14)

Demos `cases/demo/b8_tip_taper_ls/run_b8_termination_diagnosis.py` (diagnosis;
artifact `results/b8_termination_diagnosis.csv`) and `run_b8_span_blend.py`
(**8/8 PASS**; artifacts `results/b8_span_blend.csv/.png`,
`checks_b8_span_blend.csv`). Same condition throughout: M6 coarse+medium,
M∞0.5, α3.06, `farfield="neumann"`, `upwind_c=0`. The user's re-spec proposal
(span-blend of the wake-LS rows) was reviewed against the code first; the
review found the B8 verdict number standing on an unaudited metric chain, so a
**diagnosis ran before any implementation** — and it split the verdict in two.

### Diagnosis finding 1 — the committed LS exponent p = +1.341 was a ×5 METRIC ARTIFACT

`element_mach2` reports **mixed-side PLAIN elements** — exactly the
`beyond_tip` class where B8's verdict cell lives — from the aux-substituted
SIDE field (`side_potentials` puts the aux value at "−" cut nodes), but a
plain element's **assembly is single-valued on MAIN dofs**
(`mass_conservation_coo` scatters `el[plain]`). The diagnostic field is not
the assembled field: this is the **B6 `own_side_field` disease in the one
element class `own_side_field` cannot fix** (neither side field is the
assembled one there). Measured at the verdict element (medium, elem 93977,
2/4 nodes carry aux DOFs): **side 1.532 vs main 0.309**.

| untapered tip-edge box peak | coarse | medium | p |
|---|---|---|---|
| committee (`element_mach2`, the committed metric) | 0.672 | 1.532 | **+1.341** |
| **honest** (main field on mixed-side plain cells) | 0.672 | 0.984 | **+0.620** |
| honest, no-sliver (V/median ≥ 0.1) | 0.551 | 0.691 | **+0.367** |

⇒ **the honest LS exponent is the SAME OBJECT as the conforming +0.52** — the
2026-07-13 "LS 1.34 ≥ conforming 0.52" magnitude comparison is RETIRED (and
G13.1's LS-vs-conforming exponent comparison carries the same erratum; its
conforming numbers are metric-clean and stand). Fix shipped **opt-in**:
`element_mach2(mixed_plain="main")` — the default `"side"` stays bit-identical
because the **B6/B7 M_max gate locks read through this function** (flipping
the default + re-reading those numbers is a recorded backlog item, per the
2026-07-14 arbitration). Related recorded, NOT fixed: the same aux-mixed side
field feeds `element_densities`, so junk density weights DO enter the matrix
on mixed-side plain elements (measured rho_up min 0.43 vs physical ~0.87,
M0.5 medium; no NaN) — also backlog, since fixing it moves every committed LS
gate number.

### Diagnosis finding 2 — the honest residual object: a FINITE termination-ring jump, decoupled from the TE

The last cut ring carries |δ| ≈ **0.026** (max over termination-box aux
nodes), **h-independent** (0.0262 coarse / 0.0256 medium) and **TE-taper
independent** (0.0256 → 0.0258 under `vanish_linear` r_c=0.05 — while the
same taper drives Γ_last to exactly 0). **The ring jump and the TE jump are
decoupled — which is precisely why the 2026-07-13 TE row blend measured no
effect** (q = 4.73 yet p unchanged). Also recorded: the *untapered* emergent
Γ(tip)→0 property degrades under refinement (Γ_last 0.00011 coarse → 0.00218
medium; B7's ±3e-4 was coarse-only).

### The span blend — machinery + result: it WELDS its target, and the price disqualifies the model

`MultivaluedOperator(span_blend=(form, r_blend))`: per non-TE cut node j,

    w_j · [wake LS row]  +  (1 − w_j) · s_j · [φ_aux − φ_main]  =  0

with `w_j = tip_taper_factors(q_j, span_length, form, r_blend)` (the same
F(z) family, row-level per node), `s_j` = the row's own LS 1-norm (LS entries
are O(h), a bare weld O(1) — the normalization keeps the blend itself
h-invariant), and beyond-tip straddler nodes (q > span_length) welded at any
r_blend. Default `None` bit-identical (`tests/test_b8_span_blend.py`, 11
passed; B-suite 116 passed / 9 skipped; full suite 350 + 17 + 2). The blend
needs **no solver plumbing**: it lives in the cached `_ls_coo`, so Picard,
transonic and LS-Newton all inherit it through `assemble_matrix`.

| variant | ring |δ|max (medium) | honest no-sliver p | Δcl_KJ medium |
|---|---|---|---|
| none (baseline) | 0.0256 | +0.367 | — |
| vanish_smooth rb0.03 | **0.0580 (×2.26 — WORSE)** | **+1.311 (worse)** | −19.8% |
| vanish_smooth rb0.05 | 0.0010 (0.04×) | +0.388 | −20.2% |
| vanish_smooth rb0.08 | 0.0003 (0.01×) | +0.048 (confounded) | −21.8% |

Four measured facts (all asserted by `checks_b8_span_blend.csv`):

1. **The lift cost is GLOBAL**: Γ(z) scales down ~0.8× **uniformly
   root-to-tip**, including where the blended rows are bitwise identical to
   baseline (test-locked) — the global circulation mode re-levels, this is
   not local tip unloading. And it is **r_blend-insensitive** (2-point spread
   over a 1.7× dose range).
2. **Component isolation** (coarse): the straddler weld ALONE costs
   **−13.3%**, the inboard smooth blend ALONE **−10.8%** — both components of
   a sheet-side δ-pin are amplified, so this is not a defect of either piece.
   ⇒ **the implicit Kutta has no per-station target: Γ is ONE global solution
   mode, and ANY constraint intervention on the sheet near the tip re-levels
   it** — G13.2 finding (2)'s fixed-point amplification (slope b ≈ 0.93 ⇒
   ~10×) acting GLOBALLY, where the conforming secant keeps it per-station.
   That is why the same F(z) costs −1.6% on the conforming path and −20% here.
3. **The loss GROWS under refinement** (rb0.08: −16.9% coarse → −21.8%
   medium), so it **contaminates the exponent**: p_unload ≈ −0.10 of the
   apparent +0.37 → +0.05 reduction; the corrected ~+0.15 hints at a real
   partial regularization but is **not certifiable** from a 2-point ladder
   under a 20% global flow distortion — and is moot at this cost.
4. **A narrow blend (~2 elements, rb0.03) is ACTIVELY harmful** — it
   re-creates the Heaviside it was meant to remove, steeper (ring jump ×2.26,
   honest p +1.31).

### Verdict (user-arbitrated 2026-07-14)

**Both constraint-side routes are now measured dead** — TE rows (2026-07-13:
no effect on the peak) and wake-LS rows (this: target welded, ~20% global
10×-amplified lift damage). Any further cure must change the **FUNCTION
SPACE** at the termination (how the spanwise clip ends the multivalued
region), not add sheet-side constraints. **B8 is CLOSED as
CHARACTERIZED-NOT-CURED**: the honest LS tip exponent (+0.62 / +0.37
no-sliver) is the same object, at the same magnitude, as the conforming
+0.52 that every closed conforming gate lives with — so **B9 (wing-body LS
solve, M∞0.5) is UNBLOCKED**. Recorded backlog (arbitration items 2–3): the
`mixed_plain` default flip + B6/B7 M_max re-reads + G13.1-LS erratum; the
`element_densities` mixed-plain junk-weight fix. All shipped machinery
(`span_blend`, `mixed_plain`) is default-inert and stays.

## Track B / B11 — LS-path infrastructure: unified post-processing + GMRES/AMG scaling (`cases/demo/b11_ls_infra/`, 2026-07-14)

Two long-standing LS-path infrastructure gaps, closed together (a B9 enabler).

**(1) Unified post-processing** (`run_b11_unified_post.py`, 4/4 PASS).
`post/surface.py` (conforming) and `post/surface_ls.py` (level-set) now share
private cores — `_cp_from_q2` (the isentropic/Bernoulli Cp branch),
`_pressure_force` (the `-(cp·area)·n_out/s_ref` integral + lift/drag),
`_wall_plane_crossings` / `_resolve_station` / `_section_curve_dict` (the
triangle plane-cut + station-resolve + chord/x_le), `_d11_wall_state` (the D11
two-sided q² selection) — under a keyword-dispatched upper layer `post/unified.py`
(`wall_cp` / `wall_forces` / `section_cp`, `phi=` conforming vs `mvop=,phi_ext=`
level-set). The three near-duplicate blocks (Cp+D11, the copy-pasted section-cut
loop, the force integral) collapse to one implementation each. The demo solves one
NACA0012 level-set case and extracts its wall Cp through BOTH the unified entry
and the legacy functions on BOTH dispatch forms: **max|Δcp| = 0.0 exactly** on
all four (LS wall_cp, LS section_cp, conforming wall_forces, conforming
section_cp). Every legacy public function keeps its name/signature; the extraction
preserved float op order, so the bit-identity locks
(`test_b7_onera_m6.py::test_d11_upper_surface_equals_the_main_based_section`, the
shock `x_shock` asserts, `test_post_surface.py`) pass unchanged. Evidence:
`results/cp_unified_overlay.png`, `summary.csv`, `checks_unified_post.csv`;
`tests/test_b11_post_unified.py` (9 passed).

**(2) GMRES escape from the splu wall** (`run_b11_gmres_ls.py`, 4/4 PASS; the
deferred design_track_b.md §5.3 landing). `solve_multivalued_laplace` / `_lifting`
/ `_newton` grow `precond=None|"ilu"|"amg"` (None = the pre-B11 `spsolve`,
bit-identical default; `solve_multivalued_transonic` inherits via `**kwargs`).
★ **ILU is the effective escape** — spilu on the real fused nonsymmetric matrix,
warm-started per outer:

| mesh | precond | γ | \|Δγ\| vs spsolve | GMRES iters | stalls | wall |
|------|---------|-----|-----|-----|-----|-----|
| coarse | spsolve | 0.139418 | 0 (baseline) | — | — | 1.9 s |
| coarse | ilu | 0.139418 | 3.9e-10 | 434 | 0 | 2.6 s |
| medium | spsolve | 0.141376 | 0 (baseline) | — | — | 8.6 s |
| medium | ilu | 0.141376 | 7.5e-10 | 148 | 0 | 18.2 s |

(ILU is *slower* than `spsolve` at these small 2.5D sizes — expected: the win is
at M6-fine scale where the single splu factor's fill blows up, the P9 4h39m/26 GB
wall. The point here is that ILU-GMRES CONVERGES to the same solution, so the
escape route works.) The medium ILU needed a robustness ladder in
`build_ilu_preconditioner` (gentler drop + MODIFIED-ILU `SMILU_2`): the default
`drop_tol=1e-4` left the incomplete factor *exactly singular* on the fused matrix
(its aux weld / wake-LS / TE-Kutta rows carry zero/negative diagonals); MILU
compensates the dropped mass onto the diagonal. The default first attempt is
unchanged, so the conforming Newton caller's ILU success path is bit-identical.

★ **AMG does NOT converge on the lifting operator — an honest §5.3 finding.**
`_amg_surrogate_preconditioner` builds SA-AMG on an SPD surrogate (the
single-valued Picard block + SPD springs tying each aux dof to its coincident
host so SA aggregates them). This converges for the SPD `continuity`-closure
(Laplace) system, but on the `wake_ls`-closure lifting operator the aux rows are
the g₁+g₂ wake-LS + nonlinear TE-Kutta rows (convection-like, not SPD springs),
the surrogate cannot model them, and **GMRES stalls at the restart cap**:
measured coarse M0.5 lifting **γ 0.0033 vs 0.139, all 80 outers stalled, 455 s**
vs ILU's 2.7 s. So AMG stays wired for the Laplace case + as the recorded §5.3
knob, and **ILU is the shipped lifting escape**. The Núñez symmetric row
assignment (which would restore genuine AMG applicability) stays not-prebuilt.
Evidence: `results/solver_ab.csv`, `gmres_ls_ab.png`, `checks_gmres_ls.csv`;
`tests/test_b11_linear_ls.py` (8 passed + 1 gated). lagged-LU
(`direct_refactor_every`) port to `newton_ls` = recorded out-of-scope follow-up.

**(3) M6 medium headline — the splu wall quantified** (`run_b11_m6_headline.py`,
gated one-shot; G11.4; `results/m6_medium_ab.csv`). The M6 medium level-set
solve at M0.5/α3.06 (neumann far field), **67,426 extended dofs**:

| precond | γ | wall | per-outer | outers | note |
|---------|-----|------|-----------|--------|------|
| spsolve | 0.0668527 | **454.8 s** | 17.5 s | 26 | the splu wall, converged |
| ilu | n/a | 306.7 s (17 outers) | — | — | factor went singular at a hard outer |

The spsolve baseline is the splu wall, quantified: **454.8 s / 67 k dofs** — and
at the M6 FINE ~450 k dofs it is the P9 catastrophe (a single splu ran 4h39m /
26 GB without returning, killed). ★ **Honest finding: the M6-medium 3D fused
matrix resists cheap incomplete factorization.** ILU-GMRES factored and advanced
~17 of the 26 outers, but at a hard outer even the diagonal-shifted MODIFIED-ILU
(`A + 1e-3·mean|diag|·I`, `SMILU_2`) at fill 20 produced an exactly-singular
incomplete factor — near-full fill (the first-run `fill×4≈40`, which does
complete but is ~spsolve cost at this size) would be needed. So **at 67 k dofs
spsolve is still the right tool; ILU is not advantageous there.** The ILU
escape's value is the FINE-scale regime where spsolve is *impossible* on memory
(where even an expensive high-fill ILU is bounded and the only option) — that is
a feasibility argument, extrapolated, not run. The escape is *demonstrated to
converge* at 2.5D medium above (|Δγ| 7.5e-10, 148 iters, 0 stalls). The demo
records this outcome rather than crashing (the `build_ilu_preconditioner`
robustness ladder — gentler drop → MILU → diagonal shift — was added here after
the M6 medium exposed the singular-factor failure that the 2.5D meshes do not).

**Net (B11):** the unified post-processing is bit-identical and shipped; the
GMRES/ILU escape is wired into every LS driver (`precond=None|"ilu"|"amg"`,
default bit-identical) and *works* at 2.5D; AMG on the SPD surrogate is honestly
measured to stall on the `wake_ls` operator (Laplace-only); and the M6-medium
splu wall is quantified with the escape's advantageous regime placed at fine
scale. The Núñez symmetric row assignment (which would restore genuine AMG and a
cheaper factorization) and the `newton_ls` lagged-LU port are the recorded
follow-ups.

## Track B / B12 + B13 — lagged-LU direct-reuse (the medium/M6 scaling escape) (`cases/demo/b12_lagged_lu/`, `cases/demo/b13_lagged_picard/`, 2026-07-14)

B11's finding — the iterative escapes are coarse-only on the fused LS matrix
(ILU diverges at 2.5D medium lifting, `factor_failed`s at M6 medium; AMG stalls)
— means that at medium/M6 sizes sparse-direct is the only converging tool and
the cost driver is the **number of factorizations** (17.5 s per splu at 67 k
dofs). B12/B13 port the conforming N6 **lagged-LU** mechanism (refactor the LU
every k-th step, drive the steps in between with GMRES preconditioned by the
stale *exact* LU) — **minus the Woodbury**, since the LS system has no Γ DOF and
its step is a plain `J_free d = −R_free`. Default `direct_refactor_every=1` is
byte-identical per-step `spsolve`.

**B12 (`solve_multivalued_newton`, demo 6/6):** M6-medium Newton (M0.5, 67,426
dofs, 7 steps) — spsolve refactors 7× (**145.6 s**), lagged-LU (k=1000)
refactors **once** + 30 reuse-GMRES iters (**66.7 s, 2.18×**), γ bit-identical
(|Δγ| 6.7e-13), 0 stalls.

**B13 (`solve_multivalued_lifting`, demo 6/6):** the Picard OUTER loop is the
post-B12 driver. M6-medium lifting (26 outers) — spsolve **447.6 s** vs lagged-LU
**68.3 s (6.55×)**, 2 refactors vs 26, γ bit-identical (|Δγ| 6.9e-13); the
end-to-end seed+Newton pipeline **~330 s → 111.9 s (~3×)**. ★ Measured: the
Picard `direct_reuse_rtol` must be 1e-10 (not B12's 1e-8) — a Picard fixed point
is pinned only by its 1e-6 lag tolerances, so an inexact reuse step SHIFTS the
stopping point. **Honest boundary:** both amortize the factorization COUNT; they
still need ≥1 in-memory splu, so they do NOT break the FINE memory wall (that is
the B14 Schur design / Núñez route, designed-not-scheduled).

## M6 medium level-set WORKFLOW — methods × meshes × regimes (`cases/demo/m6_medium_ls_workflow/`, 2026-07-15)

The capability B11/B12/B13 unlock, made concrete: the ONERA M6 **medium**
level-set solve, both **subsonic (M∞0.5)** and **transonic (M∞0.84)**, on the
**wake-free workflow mesh** (the wake built analytically from the TE polyline,
nothing embedded), at a wall-clock now comparable to the conforming path. The
demo cross-checks three axes — **mesh** (level-set wake-free M4 vs wake-embedded
M1), **method** (level-set vs conforming), **regime** (M0.5 vs M0.84) — and is
self-checking (`checks.csv`, **10/10 PASS**). Solves cache to gitignored
`results/*.npz` (P5 policy); the committed evidence is 4 PNGs + `summary.csv` +
`checks.csv`.

### It solves — both regimes, both meshes

`summary.csv` (cl_kj / M_max / γ_root / γ_tip):

| regime | LS wake-free | LS embedded | conforming |
|---|---|---|---|
| subsonic M0.5 | 0.2116 / 1.15 | 0.2129 / 0.98 | 0.2126 / (P3-class) |
| transonic M0.84 | 0.2765 / 2.455 | 0.2789 / 2.195 | 0.2499 / 1.995 (P5 Picard) |

- **Mesh A/B (LS wake-free vs embedded):** cl within **0.62 % (sub) / 0.85 %
  (trans)** — the dual-mesh agreement B7 measured at coarse (~2 %) tightens at
  medium.
- **Method A/B (LS vs conforming):** **0.47 % subsonic**; **+10.65 % transonic**
  — the documented B6/B7 inversion (the conforming *Picard* under-circulates
  4–8 % at these shocks; the LS Picard sits closer to the conforming *Newton*
  truth, so the gap is the conforming baseline's, not the LS path's).
- **Γ(tip) → 0** on every LS solution (0.013–0.027 of Γ_max — the spanwise clip
  / free-edge rule, discretely); **M_max bounded** < M_cap = 3 (2.455 / 2.195).

### ★ The transonic solve is bounded/engineering-converged, not strict — and *why* was diagnosed live

The M6-medium transonic ramp was first run with `tol_residual=1e-7` and *looked*
stuck (~1 h, no completion). A per-level streaming diagnostic (the P4/P5/G13.3
method: ramp M0.60→0.84 manually, print each level's residual trajectory + M_max
+ lim/flr) showed it was **not a breakdown**: M0.64 drove the residual to
5.8e-7 by outer 59 then **plateaued flat at 3.5e-7** for 240 more outers, with
**M_max 1.515, 0 lim/flr** — the P4/B6/N5 transonic Picard residual floor (the
shock-position soft mode), stiffest at the shock-forming levels. `tol=1e-7` was
below the plateau, so every level burned its full budget on the floor. Re-run
above the plateau (`tol=1e-5`, the B7 engineering-converged level), the ramp
completes: M_max climbs **1.39 → 1.52 → 1.64 → 1.78 → 1.95 → 2.17 → 2.455**, all
< 3, γ rising physically 0.069 → 0.087, ≤ 3 clamped cells of 329 k. This is the
**B7 gate semantics — assert bounded, not `converged`** — the same status as the
committed P5 medium (M_max 1.995) and B7 coarse.

### The four visualizations (physical reasonableness)

- **`spanwise_loading.png`** — Γ(z) root→tip. Subsonic: the three curves overlap
  (elliptic-ish, → 0 at tip). Transonic: the two LS curves overlap and are
  **smooth**; the conforming curve carries the P5 spanwise Kutta-probe jitter and
  sits below (under-circulated) — a direct picture of B7's "LS Γ(z) is 11–12×
  smoother" and the loading inversion.
- **`section_cp.png`** — Cp(x/c) at η = 0.44/0.65/0.90. LS and conforming track
  each other; the transonic shock (Cp shoulder → drop, moving forward with η) and
  the −Cp* line are clear; the conforming TE Cp dip (its potential-jump Kutta)
  and the LE recovery sawtooth (P6/G6.1, `smooth_passes=0`) are visible and
  documented.
- **`wake_potential.png`** — perturbation potential φ′ on the η=0.65 slice: red
  above / blue below, and the **jump across the wake line persists downstream**
  (convected, not decaying) — the implicit-Kutta wake-LS signature.
- **`tip_mach.png`** — local Mach on a chord-plane slab over the outer span
  (Mach-1-centred diverging map): subsonic is near-uniformly subsonic with a thin
  LE line; transonic shows a **supersonic pocket along the whole outer-span LE**
  up to and around the tip.

**Net:** M6 medium on the level-set wake-free workflow mesh is solvable and
physically sound in **both** regimes, mesh-independent (0.6–0.9 %) and
method-consistent (0.5 % subsonic; the transonic gap is the conforming Picard
baseline's, per B6/B7), at conforming-comparable cost thanks to B12/B13. The
transonic state is the bounded/engineering-converged B7-class solution; a
strict-converged transonic would want the LS Newton ramp (newton_ls + B12
lagged-LU — deferred, the residual plateau is exactly what Newton removes).

---

## Track B / B15 — LS Newton transonic ramp + N5 freeze-selection (`cases/demo/b15_ls_newton_ramp/`, 2026-07-15)

The direct answer to the sentence that closes the section above ("a strict-converged
transonic would want the LS Newton ramp — the residual plateau is exactly what Newton
removes"). **13/13 self-checks PASS** ungated (`results/checks.csv`); part 3 (ONERA M6
medium) is gated.

### The cost is the Picard plateau, not the linear algebra

The 24.5 / 38.4 min of the M6-medium M0.84 level-set solve is **not** an inner-solver
cost: the ramp's top Mach levels (0.80, 0.84) simply **never converge** and burn their
full 200-outer budget on the **shock-position residual plateau** (the P4/B6/N5 soft
mode). `tol_residual` is already pinned *above* the plateau at 1e-5 — 1e-7 would make
every level burn its budget (~1 h). No Picard tuning escapes this; the plateau is
intrinsic to the method. Newton has no such soft mode.

Before B15 the LS Newton could not run a ramp at all: `freeze=` was a **reserved
no-op**, the convergence gate hard-requires 0 limited/floored (which shock limiter
cells block), and there was **no Mach-ramp wrapper**.

### GB15.1 — the frozen per-side selection is exact (FD gate)

The `kernels/upwind.py` frozen apparatus is reused **unmodified**: the level-set
per-side operators are already walk-mode `UpwindOperator`s with a side-masked face
graph, so B15 is *wiring*, not new numerics. New: `newton_side_data(frozen=…)` +
`MultivaluedOperator.freeze_side_state`. Residual and Jacobian were extracted into
**`LSNewtonSystem`** so the solver and the FD gate share **one** assembly path — a
Jacobian must be the derivative of the residual the solver actually evaluates, and an
FD test that re-implements the assembly only tests its own copy.

- **FD: rel 6.7e-9** (eps 1e-5), with clean round-off scaling (5.8e-8 / 6.0e-7 at
  1e-6 / 1e-7) — that scaling *is* the evidence it is a true derivative, not a
  coincidence. ε-guard excludes only **3.1%** of free rows, on a real pocket
  (`nu_max` 0.785; 1,118 elements on branches 1/2).
- The frozen sweep reproduces the live density **bitwise** at the freeze point.
- A clean freeze has **0 floored by design** (no floor is re-applied on branches 0–2)
  — which is exactly what unblocks the 0-clamped convergence gate at a shock.

### GB15.2 — the freeze cures a genuine limit cycle, and does not move the answer

`results/gb152_freeze_limit_cycle.png`. On **NACA coarse M∞0.75** the LIVE LS Newton
**does not converge**: it parks in a **period-6 limit cycle** (3.2e-7, 2.8e-7, 2.7e-7,
1.3e-6, 8.6e-7, 4.3e-7 — repeating) at |R|≈2.7e-7, three orders above tol, with **0
limited / 0 floored** — a *clean* stall, i.e. the upwind assignment churn.

| | converged | steps | \|R\| final | γ |
|---|---|---|---|---|
| live selection | ✗ | 60 (budget) | **2.84e-7** | 0.218804 |
| frozen selection | **✓** | **22** | **8.49e-13** | **0.218809** |

The freeze removes the churn; it does **not** select a different state (γ agrees to 5
digits), and it took **0 reverts**.

### ★ GB15.2 erratum — the conforming freeze recipe does NOT transfer

`solve/newton.py` arms the freeze on `(r < freeze_tol) **or** live_stalled`. Ported
verbatim, this **breaks** the level-set solver. Measured on medium M0-embedded M0.75:
the LS live residual **bounces ±2× for tens of steps while still descending** (γ
travels 0.183 → 0.243 across that stretch — slow progress in a stiff direction, **not**
a stall). `live_stalled` misfires at ~5e-6, the freeze locks a **still-moving**
assignment, the frozen step diverges → revert → re-arm: **3 reverts, no convergence** —
on a case the *untouched live path converges* (54 steps, 7.5e-12).

With the stall trigger removed, the same solve freezes late (|R| < 1e-6, assignment
settled) and converges: **53 steps, 2.1e-12, 0 reverts, landing on exactly the live
γ = 0.243305**, with `residual_unfrozen = 2.1e-12` confirming the LIVE selection agrees
there (an honest pass). ⇒ **the LS freeze arms on `freeze_tol` alone.**

Two fail-safes ship with it: `freeze_max_reverts` (default 3) **disarms** a
misbehaving freeze so it can only ever HELP, never cost convergence; and the reported
`n_limited`/`n_floored` are always re-read **LIVE** — a frozen finish shows 0 floored
*by design* and can never be its own evidence.

### GB15.3 — the ramp beats the Picard, and `intermediate_tol` is free

`results/gb153_ramp_vs_picard.png`. **NACA coarse M∞0.80 / α1.25** — the B6 gate
condition, whose same-mesh **conforming-Newton truth** is shock 0.658 / cl_p 0.459 /
**M_max 1.408**.

| method | wall | γ | \|R\| | levels converged | iterations |
|---|---|---|---|---|---|
| Picard | 44.0 s | 0.190374 | **1.55e-5** (stalled) | **3/5** | 962 outers |
| Newton ramp | **8.1 s (5.4×)** | 0.212445 | **3.10e-12** | 5/5 | 48 steps |
| Newton + `intermediate_tol` | **6.8 s (6.5×)** | **0.212445** | 3.10e-12 | 5/5 | **38 steps** |

- The **plateau is gone**: Picard's top two levels burn their budgets and it ends at
  1.55e-5 — *not a solution*. Newton reaches 3.1e-12, 0 limited / 0 floored.
- **`intermediate_tol` is free**: the final γ is **identical to 6 digits**, while
  Newton steps drop 48 → 38.
- **M_max 1.3924 vs the conforming-Newton truth 1.408 (−1.1%)** ⇒ physically
  consistent.
- **Honest boundary:** γ is **−7.4%** of that same truth. This is the **LS-vs-conforming
  discretization gap** (B6 recorded it as ~13% while the LS side was a Picard stall).
  B15 makes it measurable **strict-to-strict for the first time** — it does **not**
  close it. That remains open.

### The LS-specific ramp mask (differs from conforming on purpose)

The freeze stays **ARMED on loose intermediate levels**, where the conforming mask
(`newton.py:888`) sets `freeze_tol=None`. Reason: *every* accept route — loose ones
included — requires 0 limited/floored, and on a 0.60→0.84 ramp the shock forms
**mid-ramp**, so those levels carry limiter cells and can only reach a 0-clamped accept
**through** a freeze. Loosen the *tolerance*, keep the *mechanism*. The conforming fold
contraindication (a loose level leaving an untracked **Γ** seed, G10.2) has **no
analogue** here: the level-set path has no Γ DOF — Γ is a solution mode carried inside
`phi_ext`.

### Seed trap (found and fixed)

`_seed_from_picard` did **not** forward B13's `direct_refactor_every`, so the Newton
warm start paid a **full sparse factorization on every seed outer** — ≈17 s × `n_seed`
at M6 medium, i.e. **~11 min of pure seed** on a 40-outer warm start, dwarfing the
Newton solve it was seeding. The seed keeps the *lifting* default
`direct_reuse_rtol=1e-10` (**not** Newton's 1e-8): a Picard fixed point is pinned only
by its lag tolerances, so an inexact reuse step would *shift* where it stops (B13).

**Defaults unchanged:** `freeze_tol=None` ⇒ the pre-B15 live solver, byte-identical
(locked). Tests `tests/test_b15_ls_newton_freeze.py` (12 — 11 at closure + the
2026-07-15 errata lock `test_freeze_max_clamped_relaxes_the_convergence_semantics`).

### ★ GB15.4 — ONERA M6 medium M0.84: the plateau is gone, 3.5× faster

The target case: the one the committed Picard needs **38.4 min** to leave
*bounded-but-not-converged*. The Picard baseline is **not re-run** — it is committed
evidence (`cases/demo/m6_medium_ls_workflow/results/summary.csv`).

| | Picard (committed) | **Newton ramp (B15)** |
|---|---|---|
| wall clock | 2304.7 s (**38.4 min**) | **657.4 s (11.0 min) = 3.51×** (committed `summary.csv`; an earlier draft's 672 s was a pre-CSV trial run) |
| residual | **1e-5…1e-4 plateau** — top two levels burn their full 200-outer budget | **~1e-11, every level converged** |
| M_max | 2.4549 | 2.4938 (1.6% apart) |
| clamped cells | ≤3 / 329k | 3 / 330k |

Per level, freeze armed everywhere, **0 reverts**: M0.60 ✓5 steps · 0.65 ✓19 · 0.70
✓23 · 0.75 ✓16 · 0.80 ✓19 · **0.84 ✓16 steps, |R| 6.9e-11**.

**Honest limit (stated, not buried):** most levels accept via `assignment_cycle` —
the **frozen** system converges to ~1e-11 and is accepted at the
**assignment-discontinuity floor** (the live residual stops improving across
refreshes). That is the N5 semantics the conforming path also uses. It beats the
Picard plateau by 6–7 orders of magnitude, but calling it a "live-strict solution"
would be an over-claim.

### ★★ Four errata — porting the conforming N5 recipe is NOT mechanical

Every one was forced out by measurement; none was foreseen. This is the same lesson
B8 taught, and it is the most reusable output of B15.

1. **The TE polyline must come from the authoritative geometry.** A hand-rolled
   `x_te(0)=0.8059` vs `wing3d.x_te(0)=0.80611` — off by **2e-4**. `CutElementMap`
   matches the polyline onto **wall nodes** (M2: the M6 TE endpoints are exact wall
   nodes), so 2e-4 matches **nothing** ⇒ **0 TE nodes ⇒ no Kutta ⇒ Γ unpinned ⇒ 340k
   limited cells and NaN** — and the solver **passed silently**. Both LS solvers now
   **raise** on `te_nodes == 0`.
2. **`freeze_tol` must sit ABOVE the churn floor — and that floor RISES with Mach**
   (**<1e-6** at M0.60 → **8.6e-6** at M0.65 → **2.7e-4** at M0.70). Set below it, a
   discrete upwind-selection flip throws the residual back before the freeze can arm
   (clean descent, then ×300 to **the same value 2.6e-3, twice** — the signature of a
   discrete flip) and the ramp dies. **The same law as "`tol_residual` must sit above
   the Picard plateau".**
3. **Residuals are not comparable across a SELECTION EPOCH.** The frozen phase drives
   `r_best` to 1.5e-11; after a refresh the residual legitimately returns to the live
   scale (2.6e-3) and the fail-fast reads a 1e8× "blow-up", killing a healthy
   freeze-refresh cycle. `r_best` is now reset on every freeze / refresh / revert.
4. **The frozen phase's clamp count is stale by construction** — `n_floored` counts
   `branch==3`, the cells clamped *at the freeze point*, and never falls. Gating
   acceptance on `n_flr == 0` **refuses a 7.8e-14 machine-precision solution forever**
   (measured at M0.70: the freeze cured the period-7 limit cycle and the floored cell
   **cleared itself** — final live `lim/flr = 0/0` — yet the gate would not fire). The
   **live** re-evaluation is the arbiter.

**New knob `freeze_max_clamped`** (default **0** = the conforming N5 rule,
bit-identical): at M6 medium M0.70 a **single** persistently-floored cell of 330k
blocks the freeze at **any** `freeze_tol`. The frozen sweep represents a clamped cell
*exactly* (branch 3: `nu=0`, `rho=rho_floor`, `s_e=s_u=0`), so the 0-clamped
precondition is stricter than the machinery needs; relaxing it lets the freeze arm and
the ramp complete.

⚠ **CORRECTIONS (2026-07-15, self-caught after the first draft):**
1. **The clamped cells do NOT "clear themselves".** Over-generalised from ONE isolated
   80-step M0.70 run (driven to 7.8e-14, ending 0/0). In the SHIPPED ramp (which accepts at
   `assignment_cycle` after ~23 steps) they **PERSIST**: M0.70 `0/1`, M0.75 `0/1`, M0.80
   `1/1`, **M0.84 `1/2` = 3 clamped cells of 330k** — exactly the Picard's ≤3, so consistent,
   but the freeze proceeds **WITH** them, it does not dissolve them.
2. **The convergence semantics ARE relaxed** — "the convergence gate is untouched" was
   FALSE. With `freeze_max_clamped > 0` the `assignment_cycle` / `refresh_budget` accept
   routes do NOT re-check the clamp count, so the returned `converged=True` M0.84 state
   **carries 3 clamped cells**. Only the strict `tol` route still demands live
   0-limited/0-floored. Quote the M6 number WITH this caveat.
3. **P9/G9.1 is CITED, NOT RE-TESTED.** Its record is about permanently-**limited** cells on
   the CONFORMING path; ours are mostly **floored** — same precondition, different clamp.
   `freeze_max_clamped` exists **only on the LS path** (`newton.py` keeps the hard rule).
   Whether relaxing it would unblock G9.1's conforming fine mesh is an **UNTESTED
   HYPOTHESIS**. B15 has NOT revived G9.1.

