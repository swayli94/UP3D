# pyFP3D Roadmap — Track P (solver phases P0–P14)

> Split verbatim from `docs/roadmap.md` on 2026-07-15 (content unchanged; only
> this header and the ledger heading were added). Global working rules, gate-ID
> conventions and the track index live in [roadmap.md](../roadmap.md); the
> human-readable status snapshot is [overview.md](../overview.md).

## Track P — Solver phases

### P0 — Repo scaffolding + mesh infrastructure
**Deliverables:** package skeleton per design.md §7 layout; `mesh/reader.py`
(meshio → SoA arrays + boundary tags), `mesh/metrics.py` (B_e, V_e, face
adjacency), `mesh/coloring.py`, `post/vtk_out.py`; `agent-rules.md`; CI running
pytest on coarse meshes.
**Gates:**
- [x] G0.1 ΣV_e equals analytic volume (unit cube, sphere shell) to 1e−12 — `tests/test_mesh_volume.py`
- [x] G0.2 element-wise ∇(a·x+b) exact to machine precision on random tets — `tests/test_mesh_gradient.py`
- [x] G0.3 coloring validity: no two same-color elements share a node — `tests/test_mesh_coloring.py`
- [x] G0.4 round-trip: read → write VTK → fields visualize correctly — `tests/test_io_vtk.py`
**Visual test examples:**
- V0.1 export per-element volume error heatmap for unit cube/sphere shell; expect no spatial bias (random sign, machine-level magnitude).
- V0.2 visualize random linear field recovery: plot exact vs reconstructed gradient components on slices; expect pixel-level overlap.
- V0.3 plot element-color index in 3D; inspect adjacent elements and confirm no same-color node-sharing conflicts.
- V0.4 store all visual checks as reproducible artifacts (`artifacts/<gate_id>/*.png` + `summary.csv`) so they can run on headless CI.
**Effort:** 2–3 sessions. **Risk:** low (pure infrastructure).

### P1 — Laplace solver (ρ ≡ 1)
**Deliverables:** `kernels/residual.py`, Picard-degenerate driver (single
linear solve), Dirichlet far field, natural walls; `solve/linear.py`
(CG + pyamg smoothed aggregation).

**Renumbering note (2026-07-06):** the former G1.2 sub-items (G1.2-a0/-a/-b/-c)
were promoted to first-class gates and the whole P1 list renumbered in workflow
order. Mapping (old → new): G1.1 → G1.1 (unchanged); G1.3 → G1.2;
G1.2-a0 → G1.3; G1.2-a → G1.4; G1.2-c → DP1 (decision point, not a gate);
G1.2-b → G1.5; G1.2 → G1.6. Git history and documents dated before 2026-07-06
(including the shared P1 review report) use the old IDs.

**Gates (numbering = workflow order; G1.3 → G1.4 → DP1 → G1.5 → G1.6):**
- [x] G1.1 = V1: MMS on cube, L2 convergence slope ≥ 1.9 over 3 mesh levels — `tests/test_laplace_mms.py`, slope ≈ 1.94–1.96 (sin πx·sin πy·sin πz manufactured solution, consistent quadrature load vector; a harmonic-polynomial exact solution was tried first and rejected — it's reproduced to machine precision at every h by this structured mesh's stiffness stencil, giving no convergence signal).
- [x] G1.2 (formerly G1.3) AMG-preconditioned CG: iteration count roughly mesh-independent — `tests/test_laplace_cg_iterations.py`, iterations 8→11→14 over an 8×/level node-count increase (n=8,16,32 cube), well inside a 2× growth-ratio bound.
- [x] G1.3 (formerly G1.2-a0) — **Cylinder oracle pre-study** — **COMPLETED
      2026-07-06 with a NEGATIVE result: acceptance criterion NOT met**
      (`tests/test_wall_correction_cylinder.py`, 9 regression tests locking
      in the measured facts + the acceptance itself as a `strict=True` xfail;
      artifacts in `artifacts/G1.3/` per spec). Spec was: assemble Option A's
      correction RHS on the wall facets from the **exact gradient + exact
      normal** (closed-form: φ_exact = U x (1 + a²/r²), n = (x, y, 0)/r_xy),
      solve once, require corrected error significantly below uncorrected and
      ≈ first-order convergence. Measured outcome (coarse/medium/fine,
      6.9k/17.3k/50.2k tets):
      (i) Option A's t-form RHS assembles to **machine zero** (~1e-17) —
      adjacent-facet contributions cancel exactly on the uniformly spaced
      circle, even though single quadrature points carry O(h) flux defects
      (assembly verified against a hand-computed single-facet case);
      (ii) the **full** consistency defect ⟨N_i, ∇φ_exact·ñ⟩ (7-point
      quadrature, upper bound for any boundary-data correction) is only
      ~2e-5 max and shrinks ~O(h⁴): for a harmonic potential with
      **body-fitted wall vertices** the exact net flux through each facet is
      *exactly zero* (divergence theorem over the closed facet/true-arc
      sliver) — there is nothing for a boundary-data correction to correct;
      (iii) corrected Cp error is unchanged: 9.10e-2/4.49e-2/2.22e-2 max at
      slope 1.02, identical to uncorrected;
      (iv) **error anatomy overturns the cylinder's "same variational crime"
      designation**: exact-potential-through-recovery reproduces ~76% of the
      total Cp error at every level (the quasi-2D single-layer sliver-strip
      recovery dominates), and wall nodal φ error converges at a healthy
      ~1.2 order — the cylinder does NOT exhibit the sphere's decreasing
      sub-first-order pathology;
      (v) by-product for G2.5(b): spanwise-noise floor max |w|/U∞ =
      2.88e-2 / 1.50e-2 / 7.40e-3, unaffected by the correction.
- [x] G1.4 (formerly G1.2-a) — **Oracle ceiling experiment** (sphere) —
      **RUN 2026-07-06** (immediately, since the G1.3 mechanism made the
      outcome predictable and the run is cheap;
      demo now absorbed into `cases/demo/p1_laplace/run_demo.py`, results in
      `cases/demo/p1_laplace/results/`). Exact-gradient + exact-normal correction on the
      sphere_shell coarse/medium meshes, both the t-form and the full-flux
      variant: medium-mesh max |Cp err| = 0.1156 uncorrected → 0.1164
      (t-form) / **0.1133 (full-flux)**. The Option A ceiling is ≈ **11.3%**
      — same mechanism as G1.3: with body-fitted vertices the BC-data defect
      is near-zero; the sphere's error lives in the domain perturbation
      (missing facet/surface slivers) + P1 approximation, not in the BC data.
- **DP1 (formerly G1.2-c) — Fix-route decision point** — **DECIDED
      2026-07-06: the "> 5%" branch applies** (measured ceiling ≈ 11.3%):
      boundary-data corrections (Option A, and by the G1.3(ii) flux argument
      any correction of this family) are ruled out on body-fitted meshes;
      curved/isoparametric wall elements remain the only accuracy-improving
      route, as a separately-scoped effort, and G1.6 is to be redefined per
      Option C's geometry-consistent-reference yardstick. Note for Option B
      (Gap-SBM): its gap terms target exactly the missing-sliver domain
      perturbation identified above, so it is the one intermediate route the
      G1.3/G1.4 evidence does *not* kill — but with O(h²) gap thickness its
      expected payoff is second-order-small; treat it as optional pre-study
      material for the curved-element design pass, not as a gate.
      (Original branch spec, for the record: < 2% → G1.5 rollout; 2–5% →
      Option B pre-study; > 5% → curved elements + Option C redefinition.)
- ~~G1.5 (formerly G1.2-b) — Lagged implementation + h-sweep~~ — **VOID per
      DP1** (only applicable had the G1.4 ceiling been < 2%; it is 11.3%).
      The RHS-correction infrastructure built for G1.3/G1.4
      (`solve/wall_correction.py`) stays: assembly-verified, reusable for
      Gap-SBM-style facet integrals if Option B material is ever picked up.
- [ ] G1.6 = V2 (formerly G1.2): incompressible sphere, max |Cp − (1 − 9/4 sin²θ)| < 2% (medium) — **NOT MET.** `tests/test_laplace_sphere.py` has the real gate as a `strict=True` xfail. Best implementation so far (`post/surface.py::wall_tangential_gradient_quadratic`, quadratic surface-patch recovery) gives ~11.6% max error on medium (down from ~12.0% with the earlier linear recovery). **Root-caused**, not just hypothesized (see PROJECT_STRUCTURE.md "Known gaps" for the full evidence trail): a clean, single-variable h_min sweep plus an oracle test (exact analytic potential fed straight through the recovery step, bypassing the FEM solve) show the recovery scheme is a minor contributor at every mesh size tested — the dominant error is the volume PDE solve's own accuracy at the wall, caused by the natural (zero-flux) BC being satisfied on the flat polyhedral wall-facet approximation instead of the true curved sphere (a geometric/variational-crime inconsistency, evidenced by sub-first-order, decreasing convergence order of the raw nodal potential itself). A direct Nitsche/penalty fix was tried and rejected (made error and CG iterations both worse with increasing penalty strength — see PROJECT_STRUCTURE.md for why). Fix routes researched and tiered 2026-07-06 (design.md §5.1, Option A/B/C) and **resolved the same day by G1.3/G1.4 + DP1 above**: boundary-data corrections are ruled out by measurement (oracle ceiling ≈ 11.3%); the sanctioned route is now Option C — redefine this gate's acceptance against a geometry-consistent reference solution (drafting that criterion is the open next action) — with curved/isoparametric wall elements as the separately-scoped physical-accuracy effort. Do not re-propose h-refinement, post-processing tweaks, or further boundary-data corrections — all three are ruled out with evidence.
**Visual test examples:**
- V1.1 log-log plot of L2 error vs h for 3 meshes with fitted slope; visually confirm near-2nd-order straight-line trend.
- V1.2 CG convergence history overlay for coarse/medium/fine; curves should have similar shape and iteration count.
- V1.3 sphere Cp contour and meridian line plot against analytic curve; check front stagnation point Cp peak and rear symmetry.
**Effort:** 2–3 sessions. **Risk:** low.

### P2 — Wake cut, circulation, Kutta — on Laplace  ★ critical phase
Incompressible lifting flow over the extruded NACA0012 (mesh M0). All the
hard topology and constraint machinery lands here, against a linear operator.
**Deliverables:** `mesh/wake_cut.py` (node duplication from tagged wake
surface, TE-station → wake-line map), `constraints/wake.py` (master–slave
elimination, RHS-only Γ dependence), Kutta update loop, vortex far-field
correction (incompressible form); `post/section_cut.py`: given z = const,
extract (a) a field slice (φ, |V|, later M/ρ/Cp) and (b) wall-surface Cp(x/c)
split into upper/lower curves — headless per §0.1, writes
`artifacts/<gate_id>/*.png` (matplotlib) + numeric CSV, no GUI; needed for
V2.3 visuals and reused heavily in P3/P4 debugging. The degenerate single-layer
implementation (symmetry-plane tripcolor) is pre-warmed by the G1.3
(cylinder oracle pre-study) visualization work: the interface is defined with
its final signature (z = const parameter) already at the G1.3 stage, so P2 only
adds the general z = const interpolation path.
**Topology asserts (run at preprocess time, every mesh):**
- each wake face has exactly one ⁺-side and one ⁻-side element
- TE nodes **ARE duplicated** (**re-spec 2026-07-06** — the original "TE
  nodes are NOT duplicated" was implemented first and is quantitatively
  wrong: in the continuum [φ](TE) = Γ — that *is* the Kutta condition — so
  a single-valued TE node tapers the jump to zero over the first wake cell,
  which is equivalent to parking a point vortex of strength Γ at the TE.
  Measured on the coarse NACA0012 mesh with Γ = 0.3 prescribed: peak wall
  |V| = 4.6 U∞ at the TE, a spurious suction force of −0.27 out of an
  expected cl ≈ 0.6 from just the 6 TE-adjacent wall triangles, and the
  spike scales like Γ²/h — it *diverges* under refinement. With the TE
  duplicated: cl = 0.6012 vs the Kutta–Joukowski 0.6 (+0.2%). Evidence:
  `mesh/wake_cut.py` module docstring); wake–symmetry-plane edges handled;
  no orphan duplicated nodes
- duplicated-node count equals wake-sheet-node count (TE included, per the
  re-spec above)
**Gates:**
- [x] G2.1 = V0-with-cut: φ = x on the cut mesh, ‖R‖∞ < 1e−12 with Γ = 0
      *(closed 2026-07-06: ‖R‖∞ = 8.4e-13 over all free non-wall reduced
      dofs, 6.9e-16 on the folded wake-master rows — synthetic strip +
      NACA0012 coarse, `tests/test_p2_wake_cut.py`)*
- [x] G2.2 jump consistency: prescribe Γ = const analytically, verify [φ]
      reproduced exactly and residual clean at the cut
      *(closed 2026-07-06: [φ] = Γ to 1e-13 abs at every wake pair,
      free-dof residual < 1e-8 · scale; artifacts/G2.2)*
- [x] G2.3 = V3(incompressible): NACA0012 α = 4°, cl vs 2D panel/XFOIL-inviscid
      reference, Δcl < 2%; Kutta loop converges in < 20 updates
      *(closed 2026-07-06: medium mesh cl_p = 0.47858 vs Hess–Smith panel
      reference 0.482556 (`cases/reference_data/naca0012_incompressible/`,
      provenance in its README) → **−0.82%**; coarse −3.0% (recorded, not
      gated). Kutta loop: **2 updates** via per-station secant acceleration
      on the affine Γ-map — plain ω-relaxation alone would need O(100)
      updates because the measured map slope is b ≈ 0.93;
      `solve/picard.py::solve_laplace_lifting` docstring. artifacts/G2.3)*
- [x] G2.4 = V6: sectional cl from Γ vs surface-pressure integration, < 1%
      *(closed 2026-07-06: medium 0.478646 vs 0.478576 → **0.01%**; coarse
      0.32%. artifacts/G2.4)*
- [x] G2.5 quasi-2D consistency (criterion (b) **re-specified 2026-07-06**
      per the evidence note below): on the M0 single-layer mesh,
      (a) φ = x freestream: max over elements of |(∇φ_e)_z|/U∞ < 1e-12;
      (b) converged incompressible lifting solution of G2.3: the
      field-wide spanwise noise — 99th percentile of |(∇φ_e)_z|/U∞ over
      elements — decays at ≥ ~1st order across the mesh family, AND the
      V2.5 mid-plane |w| heatmap shows no coherent stripe along wake/TE.
      The single-element max is recorded but not gated (it sits in the
      leading-edge peak-gradient region and tracks local h·|∇²φ|, not the
      split asymmetry).
      *(closed 2026-07-06: (a) 8.5e-15 / 8.9e-15 coarse/medium;
      (b) p99 4.82e-3 → 2.35e-3 (ratio 2.05, h ratio 2 — clean 1st order),
      RMS 1.62e-3 → 8.2e-4 (1.97), max 6.4e-2 → 4.5e-2 at the LE
      ((0.014, −0.024) → (−0.0001, −0.0074)), heatmap stripe-free;
      artifacts/G2.5, `tests/test_p2_kutta_naca0012.py`)*
      **Evidence note (2026-07-06, from M0 delivery — criterion (b) needs
      re-spec before P2 closes it:** (a) holds at machine precision (verified,
      `tests/test_m0_extrude.py`, `tests/test_m0_cylinder.py`), but for any
      *solved* field (b) is unachievable as written: a 3-tet prism split is
      necessarily asymmetric under the z-mirror (on a lateral quad face,
      ∫N_i dS is S/3 or S/6 depending on the diagonal), so a z-invariant
      field is not a solution of the discrete system and the discrete
      minimizer carries O(h) spanwise noise. Measured on the delivered M0
      meshes (non-lifting solves): cylinder max|w|/U∞ = 2.9e-2 (coarse) →
      1.5e-2 (medium), NACA0012 α=0 5.4e-2 (coarse) — clean ~O(h) decay, 10
      orders above 1e-12. Suggested replacement for (b): |w| decreasing at
      ≥ 1st order over the mesh family AND no coherent stripe along wake/TE
      in the V2.5 heatmap (the bug-detector intent survives; the machine-zero
      wording does not). Achieving literal 1e-12 would require a
      z-mirror-symmetric subdivision, which needs Steiner points and violates
      the M0 3-tet spec. The G1.3 cylinder oracle pre-study (2026-07-06,
      same quasi-2D mesh family) measured the floor for re-setting criterion
      (b)'s threshold: max |w|/U∞ = 2.88e-2 / 1.50e-2 / 7.40e-3 over
      coarse/medium/fine (clean ~O(h)), identical with and without the wall
      correction (`artifacts/G1.3/summary.csv`).
**Visual test examples:**
- V2.1 on cut mesh with Γ=0, display residual magnitude on wake-adjacent cells; expect no hot stripe along wake.
- V2.2 with prescribed constant Γ, plot φ+ and φ- on paired wake nodes and their difference [φ]; expect a flat spanwise profile at target Γ.
- V2.3 NACA0012 Cp(x/c) on upper/lower surfaces versus panel/XFOIL reference; inspect TE pressure matching trend from Kutta updates.
- V2.4 spanwise sectional cl comparison chart (Γ-derived vs pressure-integrated); expect near-overlap within 1% band.
- V2.5 heatmap of |w| (spanwise velocity magnitude) on the mid-plane slice; expect uniform machine-level noise, no structure along the wake or TE. Free bug detector for asymmetric prism subdivision, wake-duplication holes, and assembly bugs.
**Effort:** 4–6 sessions — budget generously; this is the phase where the
project's main risk is retired. **Risk:** high. If G2.1 fails intermittently
across meshes, stop and harden `wake_cut.py` before touching anything else.

### P3 — Subsonic compressible
**Deliverables:** `physics/isentropic.py` (unit-tested scalar functions),
Picard outer loop with density update + under-relaxation, Prandtl–Glauert-
scaled vortex far field. Also retires the P1 assembly tech debt *before* the
outer loop makes assembly hot (Picard reassembles A every iteration): the P1
kernels in `kernels/residual.py` are correct but serial and recompute each
element's Jacobian/volume per call with in-loop allocations — precompute
per-element B_e/V_e once, remove hot-loop allocations, and wire
`mesh/coloring.py` into a colored-`prange` assembly per design.md §7 rules
2/4 (agent-rules hard rule #3), with the V0 freestream check and a
bit-identical-assembly regression as the safety net.
**Gates:**
- [x] G3.1 sphere M∞ = 0.3: Cp peak vs PG-corrected incompressible, < 2%
      — **closed 2026-07-07**, measured **0.32%** on medium (0.33% coarse):
      Cp_peak −1.33591 vs PG-corrected −1.34020, both solves on the same
      mesh + same quadratic recovery so the G1.6 flat-facet bias cancels;
      non-lifting Picard 11 iterations, monotone to the CG floor
      (`tests/test_p3_subsonic.py`, `artifacts/G3.1/`).
- [x] G3.2 NACA0012 M∞ = 0.5 α = 2°: cl vs 2D FP/panel reference, < 2%;
      Picard converges < 30 iterations, monotone residual — **closed
      2026-07-07**. Reference spec note: a corrected-panel reference has an
      inherent model band — PG 0.278774 vs Kármán–Tsien 0.291861 (±2.3%),
      wider than the gate tolerance — so the reference value is the
      PG/KT **midpoint 0.285318** with an additional inside-the-bracket
      assert (`cases/reference_data/naca0012_m05/README.md`). Measured:
      medium cl = 0.284372 → **−0.33%** from the midpoint and inside
      [PG, KT]; **15 density iterations** (nested Picard, see below),
      residual strictly monotone 5.9e-6 → 5.0e-11; max local M = 0.726
      (subcritical, ν ≡ 0 regime) (`tests/test_p3_naca0012_m05.py`,
      `artifacts/G3.2/`). Two solver-structure findings with evidence:
      (1) Γ interleaved one-update-per-density-iteration (the literal §8
      pseudocode) injects 10× residual spikes down the whole history —
      the shipped loop NESTS the P2 secant Kutta iteration at frozen ρ
      inside each density update, which restores a strictly monotone
      residual and re-validates the secant's exact-affine assumption;
      (2) a loose *relative* inner CG tolerance is not a usable inexact
      scheme once warm starts are on (CG exits at x0 without computing the
      density correction — measured false convergence); the opt-in
      acceleration is a forcing-term `atol = η‖b − A x0‖` (η = 0.05 ≈ 2×
      faster, bounded residual-tail noise; default η = 0 keeps the gate's
      monotone criterion as written).
- [x] G3.3 all P1/P2 gates still green with ρ machinery in the loop (ν ≡ 0
      path must be bit-identical to Laplace when M∞ → 0) — **closed
      2026-07-07**: full suite 117 passed + 2 xfailed (~96 s; the P1/P2
      drivers now run the SAME colored assembly as the Picard loop, so the
      guarantee is structural); at M∞ = 0 the assembled matrix bits AND
      the full secant-Kutta solve (φ and Γ) equal `solve_laplace_lifting`
      bitwise, and the non-lifting path equals `solve_laplace` bitwise
      (`tests/test_p3_subsonic.py::TestG33BitIdenticalLaplaceLimit`).
      Enabling fix: pyamg's spectral-radius estimate starts from an
      unseeded `np.random` vector, making even two identical P2 solves
      differ at 2e-11 — `solve/linear.py::build_amg_preconditioner` now
      pins (and restores) the seed, so every solve in the code base is
      bit-reproducible run-to-run.
**Visual test examples:**
- V3.1 sphere Cp contour at M∞=0.3 and line cut against PG-corrected baseline; check compressibility amplification pattern is symmetric.
- V3.2 NACA0012 surface Mach and Cp plots at M∞=0.5; verify suction peak shift and physically smooth recovery to TE.
- V3.3 residual-vs-iteration plot in semilog scale; verify monotone decrease and no sawtooth instability from density update.
**Effort:** 2–3 sessions. **Risk:** low-medium.

### P4 — Transonic: artificial density
**Deliverables:** `kernels/upwind.py` (upstream-element search via face
adjacency, ν switch per design.md (3.2), ρ̃), Mach continuation in
`solve/continuation.py`, shock monitors (max local M, shock x/c extraction).
**Gates:**
- [x] G4.1 = V4: NACA0012 M∞ = 0.80 α = 1.25° (extruded): monotone shock over
      2–3 cells, no expansion shock, shock Δx/c < 0.03 vs published FP results
      — **CLOSED 2026-07-07** (re-closed the same day it was found open by
      audit: divergence → root cause verified → local-damping fix landed →
      medium gate re-measured green; full trail in the "diagnosis
      follow-up" and "fix landed" notes below). Reference spec note
      (provenance in `cases/reference_data/naca0012_m080/README.md`): an
      open digitized *FP* table for this case was not retrievable, so the
      reference is the **Euler anchor 0.60–0.63c + documented
      conservative-FP aft-shift band, gated as x/c = 0.62 ± 0.03**; the weak
      lower shock (~0.35c) is reported, not gated. Coarse evidence green,
      re-measured under the fix (`cases/demo/p4_transonic/`, `damping_theta`
      now the default): upper/lower shock 0.604 / 0.358, cl_pressure 0.357,
      M_max 1.373, monotone jump, ≤ 3 stations, no expansion shock,
      Γ-secant |F| = 1.50e-4 < tol 2e-4, zero limited/floored cells.
      **Medium measured (artifacts/G4.1/summary_medium.csv, fix landed):
      PASS in 16m39s wall** (vs the divergent attempt's 2h43m) — upper
      shock x/c **0.633** (in the 0.62 ± 0.03 band), lower 0.364,
      cl_pressure 0.349 / cl_KJ 0.354 (sign-consistent, unlike the divergent
      run's −0.171 / +0.212), M_max 1.366, Kutta |F| = 1.23e-4 < tol 2e-4,
      zero limited/floored cells, n_picard_total 12931 (vs 19331 divergent —
      fewer Γ evals needed, not merely faster ones: better conditioning cut
      CG from ~22 to a few iterations/outer). (`tests/test_p4_transonic.py`,
      `artifacts/G4.1/`; the medium gate + G4.3 sweep still run only under
      `PYFP3D_TRANSONIC_GATES=1`, excluded from the default suite — now
      minutes rather than hours of Picard, but P7 Newton remains the target
      for the O(seconds) speed the roadmap ultimately wants.)
      **Scheme-hardening evidence trail** (each item forced by a measured
      failure; details in demo_report §P4): (1) multi-hop upstream walk —
      single-hop reach on prism-split meshes is only ~0.37 of the
      streamwise extent (median), putting effective ν below the (M²−1)/M²
      bound for M ≳ 1.2: blow-up at M0.75 for every ω ∈ 0.3–0.9 ×
      C ∈ 1.5–2.0, while the walk restores median reach 1.00 and a clean
      contraction; (2) q² limiter at M_cap = 3 — without it transients
      reach the vacuum limit and the positivity guards stabilize *spurious*
      dead-cell attractors (measured off-body supersonic blobs at |y|≈0.6
      acting as blockage, Cp garbage); (3) shock-point operator
      ν = max(ν_e, ν_up) — the first subsonic cell downstream of the shock
      is otherwise purely central on the field's largest jump;
      (4) pseudo-transient diag(m_lumped/Δτ), Δτ ≈ 1e-3–3e-3 (design.md §8
      acceleration 4 pulled forward) — exact-solve Picard limit-cycles at
      M0.80 for every φ/ρ̃ relaxation tried; SER Δτ-ramping re-ignites the
      instability, so Δτ stays fixed and convergence is **engineering
      semantics**: physical M_max, zero limited cells, Kutta mismatch below
      eval noise, cl drift < 1e-3 per hundreds of iterations, residual in a
      bounded slowly-decaying shock-cell tail (NOT 1e-10; Newton in P7 is
      the designed cure — documented in `solve/picard.py` and
      `solve/continuation.py`).
      **Diagnosis follow-up (2026-07-07, later session — hypothesis
      VERIFIED, one candidate route RULED OUT, three stabilizer candidates
      measured; ad-hoc diagnostic scripts, not yet a committed gate):**
      damping ratio (m_lumped/Δτ)/diag(K) at Δτ=2e-3 measured directly:
      wall-node median 0.035 (coarse) → 0.0092 (medium, ~4× weaker). A
      standalone frozen-Γ density driver (bypassing the Γ secant entirely)
      reproduces the M0.75 blow-up in 50 iterations (M_max → 47), proving
      the divergence is in the density solve, not the Kutta loop. Δτ=5e-4
      (h²-rescaled) stabilizes M0.75 (M_max 1.22, zero limited cells) but
      **M0.80 still diverges from that state** — the needed damping also
      grows with shock strength, so no single Δτ rescaling by mesh size
      alone suffices. A finer continuation step (intermediate M0.775, still
      at Δτ=5e-4) **also diverges at M0.80 — the "finer dm" candidate route
      is ruled out**; this is not a per-step transient-overshoot problem.
      Three candidates DO stabilize the M0.80 step from a converged M0.75
      state (500 iterations each, zero limited/floored cells throughout):
      (a) global Δτ=2e-4 (M_max 1.37); (b) **local pseudo-time damping
      θ·diag(A_free)**, θ=0.2, replacing the mass-lumped global form —
      mesh- and shock-strength-independent by construction, the recommended
      fix; (c) upwind_c 1.5→2.0 at Δτ=5e-4 (M_max 1.34, more dissipative —
      would move the shock and needs a G4.3 re-validation). Under all
      three, Kutta mismatch |F| still drifts to 4e-4–9e-4 (> tol_gamma =
      2e-4) — the same drift character measured on coarse (|F| 2.9e-4 at
      100 iterations inside one eval → 8.9e-4 at the full 800, i.e. the
      eval's OWN target is budget-dependent, not converging to a fixed
      point within the budget) — so re-closing the gate also needs an
      adaptive per-eval exit or a re-matched tol_gamma, not just a
      stability fix. Separately (why coarse itself is slow, not only
      medium): the eval path solves each inner CG exactly (`forcing=0`,
      unlike the P3-validated η≈0.05 forcing-term acceleration) and
      rebuilds AMG every 4 iterations; medium profiling shows CG at 64% of
      eval wall-time (22 CG-iterations/outer vs 3 on coarse — itself a
      symptom of the weak damping's poor conditioning), AMG rebuild at 27%
      of coarse eval wall-time. No equation-level bug was found in
      picard.py/upwind.py/wake.py during this session.

      **Fix landed and verified (2026-07-07, same day as the diagnosis
      above):** `solve/picard.py::solve_subsonic_lifting` gained a
      `damping_theta` parameter — D = θ·diag(A_free) recomputed fresh every
      outer iteration from that iteration's own (upwinded) operator,
      mutually exclusive with the retired `pseudo_dt` global form (passing
      both raises). `solve/continuation.py::TRANSONIC_DEFAULTS` now
      defaults to `damping_theta = 0.2` instead of `pseudo_dt = 2e-3`. G4.2
      bit-identity is untouched (both new params default `None`, no-op for
      any caller that doesn't pass them). Calibrated on coarse first per
      the plan above: G4.1 coarse shock position holds (0.604 vs the prior
      0.599) and G4.2 stays bitwise. The medium gate then **passed outright
      on the first attempt with the new default and NO other changes**
      (no forcing, no wider `amg_rebuild_every`, no adaptive |F|-drift exit
      needed) — see the G4.1 gate entry above for the numbers. A
      reduced-budget stability probe run first (`max_gamma_evals=4,
      n_picard_eval=150`, 77 s) had already shown the mismatch drift
      predicted above did not materialize at this θ: only 1 Γ eval was
      needed at both M0.75 and M0.80, vs 12 exhausted before. The G4.3
      10-case sweep was re-run under the new default and stays green (all
      converged, zero limited cells); one measured difference worth
      recording as evidence rather than a regression: the M0.82/α=1.25°
      corner's cl moved 0.389 → 0.458 and its Kutta |F| sits at 1.92e-4,
      close to the 2e-4 tolerance — the sweep gates on convergence and
      physicality, not an exact cl, so this does not fail it, but it is the
      one corner closest to the tolerance boundary under the new damping.
      Separately, the P5 blocker identified in this same diagnosis --
      `constraints/wake.py::WakeConstraint.update_matrix`'s per-station
      `h_j` loop -- is also closed: batched into a single `T^T @ (A @ G)`
      sparse product (G = one indicator column per station) instead of one
      sparse matvec per station, verified bit-identical to the old loop on
      the real ONERA M6 coarse mesh (83 stations). Full default suite
      unaffected by either change: 136 passed + 2 skipped + 2 xfailed,
      ~5 min.
- [x] G4.2 subcritical no-op: with M∞ = 0.5, ν ≡ 0 everywhere and results
      bit-identical to P3 — **closed 2026-07-07**: max ν = 0.0 exactly (the
      switch is `max(0, ·)` and subcritical states never enter it), and the
      full solve with the P4 machinery active (upstream walk + ρ̃ sweep in
      the loop, upwind_c = 1.5) is **bitwise identical** in φ and Γ to the
      literal P3 code path (upwind_c = 0 bypasses the sweep)
      (`tests/test_p4_upwind.py::test_g42_subcritical_noop_bitwise`,
      `artifacts/G4.2/summary.csv`).
- [x] G4.3 robustness sweep: M∞ ∈ {0.74…0.82} × α ∈ {0°, 1.25°} all converge
      with one parameter set (C, M_c, ω fixed) — **closed 2026-07-07**: all
      10 cases with `TRANSONIC_DEFAULTS` (C = 1.5, M_c = 0.95, Δτ = 2e-3,
      ω_seed = 0.9, budgets n_picard_eval = 800 × ≤ 12 Γ evals — budgets are
      load-bearing: the frozen-Γ evals never meet tol_rho early-exit at
      transonic, so each runs its full n_picard_eval; cost is exactly
      levels × evals × budget), zero limited/floored cells, smooth trends and no
      outlier: α = 0 gives cl ≈ −7e-4 (symmetry check) with the shock
      marching 0.244 → 0.528 over M 0.74 → 0.82; α = 1.25° gives shock
      0.345 → 0.698, cl 0.245 → 0.389, M_max 1.19 → 1.40
      (`artifacts/G4.3/summary.csv` + V4.3 dashboard; 22 min wall on
      coarse).
**Visual test examples:**
- V4.1 contour plot of local Mach on airfoil mid-span plane with shock marker; check one dominant shock, 2-3 cell thickness, and no downstream expansion shock. Reuse `post/section_cut.py` (P2 deliverable) for the slice/Cp extraction rather than ad-hoc plotting.
- V4.2 compare M∞=0.5 fields from P3 and P4 with a difference heatmap; expect machine-zero map confirming ν no-op.
- V4.3 parameter sweep dashboard (small multiples): residual convergence + shock x/c; verify smooth trend with Mach/alpha and no outlier divergence.
**Effort:** 3–5 sessions. **Risk:** medium — expect Picard stall tuning
(design.md §12.4 mitigation ladder: raise C → lower ω → continuation).

### P5 — 3D validation: ONERA M6 (CLOSED 2026-07-08; V6 < 1% deferred to P9 curved elements — G6.3 confirmed it is the sharp-edge P1 floor, not the sawtooth)
**Gates** (status 2026-07-08 — closed after the second re-diagnosis; demo
16/16 PASS, `cases/demo/p5_onera_m6/`):
- [x] G5.1 = V5: M∞ = 0.84 α = 3.06°: λ-shock topology; section Cp at
      η = 0.44/0.65/0.90 within FP-literature scatter; CL reported with mesh
      convergence (3 levels).
      *Coarse (55.5k): physical M_max 1.398, 0 floored/limited; upper shock
      present/monotone/forward-migrating x/c 0.596→0.570→0.425; CL 0.2419.
      Medium (350.7k): physical M_max 1.995, 0 floored/limited; shocks
      0.594→0.526→0.345 (~1 cell wide); CL 0.2453; Kutta |F| 5.8e-4 (28×
      tighter than the pre-fix 1.6e-2). Section Cp = qualitative viscous-AGARD
      overlay (inviscid FP shock sits aft — the documented tendency); 3-level
      CL convergence scoped coarse+medium (no fine transonic per the cost
      decision).*
- [x] G5.2 spanwise Γ(η) smooth to the tip; V6 consistency **reported at the
      discretization floor (< 3% bound); the original < 1% is RE-SPECCED
      (2026-07-08) as a post-P6 target** — measured V6 is a systematic O(h)
      CL_p-below-CL_KJ offset (coarse 2.40%, medium 1.82%) from the sharp-TE/LE
      P1 wall gradient + the P4 surface-Cp sawtooth (both P6 targets), NOT a
      wake/far-field defect: removing the M>2 clusters left V6 unchanged
      (T4/polish, INVESTIGATION_kutta_closure.md).
      *Coarse: Γ 0.097→0.0206 smooth; medium: Γ 0.097→0.0151 smooth (the
      pre-fix st133 dip was the closure failure, healed).*
**Visual test examples:**
- V5.1 ONERA M6 wing surface pressure/Mach map with iso-lines; visually confirm lambda shock legs and spanwise migration.
- V5.2 sectional Cp plots at eta=0.44/0.65/0.90 with literature bands; check each section sits inside expected scatter envelope.
- V5.3 spanwise Γ(eta) and CL mesh-convergence chart (coarse/medium/fine); expect smooth tip decay and monotone CL convergence.
**Effort:** 2–4 sessions (solver) + Track M1 meshing effort. **Risk:** medium;
failures here usually trace back to wake–tip topology (Track M) or far-field
distance, not the field operator.

**CLOSED — how the medium `physical` failure was actually solved (second
re-diagnosis, 2026-07-08, T1–T4; this OVERTURNS the earlier "TE
discretization singularity, NOT a wake/Kutta change" conclusion).** The
pre-fix medium solve failed `physical` (M_max 5.204, 8 floored / 4 limited;
26/350 718 cells M > 2 = 18 on the wing at the outboard TE, z/b 0.80–0.81,
plus 8 at the far-field sphere beyond the tip). Cell-level audits + four
targeted experiments localized and then *refuted* every candidate in turn
(full lab record: `cases/demo/p5_onera_m6/INVESTIGATION_kutta_closure.md`;
all audits reproducible via `python cases/demo/p5_onera_m6/diagnose_medium.py`):
- **T1 straddle census**: no tet anywhere contains both a master and a slave
  wake node (0/350 718) — the "Γ-jump compressed into one P1 cell" hypothesis
  is dead; the cut topology is clean.
- **T2 stabilization test**: at the SAME (cached) Γ, under-relaxing the
  density (omega_rho 0.5) cuts the eval drho limit cycle ~4× (0.18 → 0.047)
  yet the defect stays bit-identical → **not under-convergence**; the defect
  is a genuine fixed point of that Γ.
- **T3 single-station Γ scan (decisive)**: the per-station Kutta audit showed
  the mismatch was a SINGLE-station anomaly — st133 (z/b = 0.801, the
  steepest-|dΓ/dz| station) was left **32% under-circulated** (Γ 0.0431 vs
  its own Kutta target 0.0592; |F|/Γ = 37% there vs ≤ 5% at all 165 other
  stations). Setting ONLY that station's Γ to its target collapses the whole
  cluster (band M_max 3.10 → 1.16, 18 → 0 cells, monotone in Γ₁₃₃) — same
  mesh, same h, same TE elements. **A 1/h discretization singularity cannot
  depend on the circulation value → the TE-singularity attribution was
  wrong**; the refinement-sharpening signal was the finer mesh resolving a
  sharper overspeed around the same Γ *error*.
- **Root cause**: the continuation's per-station secant does NOT converge Γ
  at the top Mach level on the 3D medium mesh — pushing it harder (16 evals)
  diverges outright (measured M_max ≈ 29); the 10-eval budget was
  early-stopping regularization that left st133 under-circulated. The Γ
  deficit drives a real TE overspeed which the ρ̃-floor + speed limiter then
  freeze into a spurious M ≈ 3 state (the `q2_at_mach` "guards stabilize
  spurious states" mechanism). The swept-TE probe-assignment degeneracy
  (adjacent stations share Kutta probe nodes; st133/134 share their upper
  probe — `diagnose_medium.py` audit) is the latent reason st133 specifically
  stalls; recorded as a known-robustness item, not fixed for this gate.
- **Far-field 8 cells** (incl. the pre-fix M_max = 5.20 cell): confirmed
  independent far-field-BC artifact — the span-uniform 2D vortex prescribes a
  ~Γ_mean jump across its branch ray **beyond the tip**, where no wake cut
  exists to carry it (the cells hug y = 0⁺ at z/b 1.0–1.15 on the sphere).

**Fix (recipe-level, every non-3D path bit-identical; landed 2026-07-08):**
1. `farfield_spanwise_gamma=True` (`constraints/dirichlet.py::
   farfield_dirichlet`): Γ(z)-tapered vortex far field — per-station Γ
   interpolant, exactly 0 at/beyond the sheet tip (the first-order horseshoe
   form design.md §5 always intended). Kills the far-field cluster.
2. `n_kutta_polish=4` (`solve/continuation.py::solve_transonic_lifting`): a
   **fixed-Γ Kutta-closure polish** after the Mach continuation — apply the
   measured Kutta target, re-solve with under-relaxed density
   (`omega_rho_polish=0.5`), repeat. Secant-free (no overshoot), contractive:
   |F| halves each step down to 5.8e-4. NOTE omega_rho < 1 *inside* the
   active-secant continuation is NOT the fix — it destabilizes the top-Mach
   secant (measured M_max → 29); only the polish under-relaxes.
   Both parameters default off (0 / False) — P2/P3/P4 and all 2.5D paths are
   bit-identical (full default suite re-verified 140 + 4 + 2).
From-scratch medium verification: physical M_max 1.995 (the mild bounded
tip-TE-corner P1 overshoot — the only surviving trace of the singularity
family), 0 floored / 0 limited, Kutta |F| 5.8e-4; demo 16/16 PASS.

**Two dead routes stay dead** (do not re-attempt): spanwise-Γ *smoothing*
(A–E record, `INVESTIGATION_gamma_smoothing.md` — it moves Γ AWAY from the
self-consistent value; the polish drives the stuck station TOWARD its own
target, the opposite operation), and TE-element-level treatments *for this
gate* (the T3 scan showed the gate-failing amplitude was never the TE
discretization; the residual tip-corner M ≈ 2.0 overshoot remains a real but
bounded P6/curved-element accuracy item).

**Deferred**: V6 < 1% → **P9** (curved wall elements). Re-measured under the P6
recovery smoothing (G6.3, 2026-07-08): the sawtooth is *not* the cause — V6
*worsens* 2.40% → 3.35% with smoothing (the ±sawtooth cancels in the surface
integral; smoothing smears the LE peak), so the entire V6 floor is the
sharp-TE/LE P1 wall gradient, which only curved/isoparametric wall elements
(P9) can remove.

### P6 — Surface-pressure recovery (remove the Cp sawtooth) ★ accuracy phase

> **Track-P renumber (2026-07-08).** The old P6 "differentiable
> artificial-density flux (remove the sawtooth)" conflated three now-separate
> deliverables. Split by the N1 finding into: **P6** (this — the sawtooth,
> which is a *recovery* artifact), **P7** (the differentiable flux, a Newton
> prerequisite), **P8** (fully-coupled Newton + performance, was P7), **P9**
> (curved wall elements, new), **P10** (backlog, was P8). Gate-ID mapping:
> old G6.x → new **G6.x (recovery)** + **G7.x (flux)**; old G7.x (Newton) →
> **G8.x**; new **G9.x** (curved). Docs/history before 2026-07-08 use the old
> IDs.

**Motivation.** The shipped transonic solve shows a bounded ≈2-cell
**sawtooth** in the supersonic-run surface Cp (see demo_report §P4 and
`cases/demo/p4_transonic/results/g41_cp_shock_coarse.png`). **N1 root cause
(2026-07-08, ★ below):** it is a **per-triangle wall-gradient RECOVERY
artifact** on the sliver prism-split wall triangulation — the piecewise-constant
P1 gradient of a *smooth* φ oscillates cell-to-cell — **not** the
artificial-density flux (a smoother flux does not reduce it; smoothing the
recovery drops the metric ~330×). Removing a non-physical oscillation from the
reported surface pressure and loads is a first-class accuracy requirement (it
corrupts Cp, sectional loads, and any pressure-based objective), so it gets its
own phase; the fix is a post-processing recovery, cheap and refinement-free.

**Deliverable (landed 2026-07-08, `30ab85e`).** A normal-gated wall-gradient
recovery smoothing (`post/surface.py::smooth_wall_tangential_gradients`) applied
in post-processing to the Cp curve **and** the force integral (`smooth_passes`
on `wall_cp_curve` / `section_cp_curve` / `wall_force_coefficients`); the metric
is `post/section_cut.py::cp_oscillation_metric` (shock-robust, sign-alternating).
Design: design.md §9.1. The differentiable *flux* work moved to P7.

**★ ROOT-CAUSE CORRECTION (N1, 2026-07-08): the sawtooth is a wall-Cp
gradient-RECOVERY artifact, not the artificial-density flux.** Decisive
evidence (design.md §3.1/§9.1): nodal/edge-neighbour smoothing of the *same*
walk solution's wall gradient drops the G6.1 metric ~330×, while a smoother
flux (the streamline kernel) does not reduce it at all. So G6.1 is fixed in
**post-processing recovery** (`smooth_wall_tangential_gradients`, normal-gated,
preserves the sharp TE — `smooth_passes` on `wall_cp_curve` / `section_cp_curve`
/ `wall_force_coefficients`), and the differentiable **flux** (streamline
kernel) is moved to **P7** as an *optional* Picard-speed path (re-scoped
2026-07-08: the P7 Newton prerequisite is the frozen-walk ∂ρ̃/∂φ, **not** the
kernel; the kernel does not address the sawtooth either). Metric now shock-robust (sign-alternating: counts only
slope-reversal points, so the monotone shock is excluded).

**Gates (repointed 2026-07-08):**
- [x] G6.1 surface-Cp smoothness: the sign-alternating sawtooth metric
      (`post/section_cut.py::cp_oscillation_metric`) on the **G4.1 coarse**
      supersonic-run wall Cp drops far below the raw per-triangle baseline under
      the recovery smoothing. Measured (coarse NACA0012 M0.80 walk solution,
      1 pass): **0.0758 → 0.00023** (330×, reversals 39 → 1). *(The original
      "coarse ≤ medium baseline" spec is obsolete: the fix is a recovery, not a
      flux, so it works at every mesh level; the acceptance is smoothed ≪ raw
      with the TE preserved.)*
- [x] G6.2 physics preserved: the smoothing is post-processing — the solve
      (φ, Γ, shock, M_max, cl_KJ) is **unchanged**; the smoothed shock stays in
      the 0.62 ± 0.03 band (0.607 coarse), smoothed cl_p within a few % of the
      raw (−0.3 % coarse), |Cp|_max unchanged (TE not polluted). Near-field
      cd_pressure shifts (~15 % coarse) — the explicitly-untrusted FP quantity,
      not gated (design.md §9). P4/P5 gate numbers (shock/M_max/cl_KJ) are
      untouched; only the reported cl_p/cd_p change and are re-recorded.
- [x] G6.3 V6 re-measured under smoothing (M6 coarse, 2026-07-08): the sawtooth
      is **not** what inflates V6. `smooth_passes` 0→1→2 gives V6 = |CL_p−CL_KJ|/
      CL_KJ = **2.40% → 3.35% → 3.88%** (CL_p 0.2419→0.2396→0.2383, CL_KJ
      0.2479) — smoothing moves CL_p slightly *further below* CL_KJ because it
      smears the LE suction peak; the ±sawtooth largely cancels in the integral.
      So the whole V6 floor is the sharp-TE/LE P1 wall gradient → **P9** (curved
      elements); V6<1% stays deferred there. **Consequence:** `smooth_passes>0`
      is for the reported **Cp curve** (removes the sawtooth); for the **force
      integral** keep `smooth_passes=0` (raw CL_p is closer to the trustworthy
      CL_KJ). The `smooth_passes` param on `wall_force_coefficients` stays
      (opt-in, default 0) but is *not* recommended for loads.
- [x] G6.4 no regression: `smooth_passes = 0` bit-identical to the current
      per-triangle path; full default suite green (157 passed).
**Demo:** `cases/demo/p6_surface_recovery/` (6/6 PASS incl. the gated M6 part) — re-runs the P4
(NACA M0.80) + P5 (M6, gated) cases showing raw→smoothed Cp, physics
preservation, and the kernel negative result. See demo_report §P6.
**Visual test examples:**
- V6.1 coarse wall Cp raw per-triangle vs recovery-smoothed on one figure;
  the smoothed curve loses the sawtooth, shock crisp, TE preserved.
- V6.2 one-panel finding: walk-flux-raw + kernel-flux-raw (both serrated) vs
  walk-flux-smoothed (clean) — a smoother flux does not help; the recovery does.
- V6.3 (gated) ONERA M6 section Cp at η=0.65 raw vs smoothed.
**Effort:** done (recovery + demo + tests). **Risk:** low (post-processing,
`smooth_passes=0` default keeps everything bit-identical).

### P7 — Differentiable artificial-density flux at frozen selection (Newton prerequisite)

> **Re-scoped 2026-07-08 (was "ship the streamline-Gaussian kernel").** Following
> the analysis to its end: the Newton prerequisite is **not** a new flux — it is
> the *existing walk flux made differentiable at frozen selection*. The kernel is
> demoted to an optional Picard-speed path. See the motivation.

**Motivation.** The exact Newton Jacobian (P8, design.md §6.3 / López
Appendix B) needs ∂ρ̃/∂φ **at frozen upstream selection** — López re-selects
u(e) each Newton step but never differentiates the selection (design.md §3.1).
At frozen u(e) the shipped P4 **walk** flux is *already* differentiable: ρ_e and
ρ_up follow the C¹ isentropic law, and the switch derivative ∂ν/∂φ exists on each
branch (López B.3–B.6). Its only kinks are the two `max` —
`ν_e = C·max(0, 1−M_c²/M_e²)` at M_e = M_c and the shock-point `max(ν_e, ν_up)` —
and **López keeps exactly this C⁰ switch and still converges quadratically** (the
active set freezes near the solution; design.md §3.1). **So the Newton
prerequisite is small and sparse: derive + finite-difference-verify ∂ρ̃/∂φ for
the shipped walk flux at frozen selection — the walk is kept, not replaced.** Its
upstream coupling stays ~+1 element/row (the walk's single u(e)) — the sparse
structure closest to López's one-hop stencil. This phase does **not** touch the
sawtooth (a recovery artifact, P6).

**Optional `max_ε` (robustness only, not required for correctness).** A smooth
`max_ε` on the **inner** `max(ν_e, ν_up)` removes Jacobian chatter if the active
set churns *far* from the solution during load stepping; it is not needed for
terminal quadratic convergence (López's C⁰ `max` suffices). **Never smooth the
outer `max(0, ·)`** — that clamp is what makes ν ≡ 0 (ρ̃ = ρ, bit-identical to
P3) subcritically; smoothing it would leak ν into subsonic cells and break
G4.2/G7.1. The outer clamp's kink is inactive at the converged state (0
floored/limited on P5), so it never enters the converged Jacobian.

**Kernel demoted to an optional Picard-speed path (2026-07-08).** The multi-ring
streamline-Gaussian kernel built during P6/N1
(`UpwindOperator(weighted=True, mode="kernel")`, depth-3 BFS + Gaussian, C^∞ in
∇φ, ~10× faster per density iteration than the walk — `tests/test_p6_weighted_flux.py`,
params threaded through `solve/picard.py`+`solve/continuation.py`) is **not** the
Newton prerequisite: N1 proved it (a) does not fix the sawtooth and (b) is
unnecessary for differentiability (the frozen walk above suffices). Its cost for
Newton is a **materially denser Jacobian** — Term 3 couples each element to its
*whole* depth-3 neighbourhood, not +1 element (design.md §6.3). It stays shipped
opt-in (`mode="kernel"`, walk the default) as a Picard accelerator, and is a
**candidate** P8 Newton flux *only if* N2 measures its dense Jacobian
net-favourable vs the sparse walk. Its `upwind_c`≈2.0–2.5 recalibration +
shock-ladder re-validation travel with it — deferred, not a Newton gate.

**Gates (closed 2026-07-10; scope decision: `max_ε` NOT implemented — the
forward walk flux is byte-untouched, so G7.1/G7.2 hold by construction and the
whole phase is the derivative + its verification):**
- [x] G7.1 subcritical no-op — **held by construction (2026-07-10)**: no
      `max_ε` was applied (López's C⁰ `max` already reaches quadratic
      convergence, so smoothing is not needed for P8; it remains an option IF
      a Newton transient ever chatters); the walk flux is byte-identical to
      P4 and the no-op stays locked by the G4.2 bitwise test
      (`test_p4_upwind.py`, re-run green).
- [x] G7.2 solution preserved — **held by construction (2026-07-10)**: no flux
      change of any kind; V0 + G4.2 re-run green, plus a forward-path
      regression guard (`test_p7_diff_flux.py::
      test_forward_walk_flux_unchanged_by_sensitivity_call`).
- [x] G7.3 differentiability (core gate) — **CLOSED 2026-07-10**: exact
      branch-wise sensitivities `(s_e, s_u) = (∂ρ̃_e/∂q²_e, ∂ρ̃_e/∂q²_{u(e)})`
      of the shipped walk flux at frozen u(e) (López B.3–B.6, + the ρ̃-floor
      flat branch → 0, + self-upstream → (ρ'_e, 0)) implemented in
      `kernels/upwind.py::rho_tilde_sensitivities_sweep` /
      `UpwindOperator.rho_tilde_sensitivities` (walk mode only; new physics
      scalar `physics/isentropic.py::mach_squared_derivative_wrt_q_sq`), and
      verified as a directional (JVP) central difference against the *shipped*
      `rho_tilde_sweep` at frozen selection. Measured: **max rel err 3–5e-10**
      (gate < 1e-6; pure FD-noise level) across subsonic / accelerating /
      shock-point (ν_u > ν_e) / self-upstream / floored regimes on the
      structured cube (`tests/test_p7_diff_flux.py`, 8 tests, green under
      `PYFP3D_NOJIT=1` too); **3.5e-9** on a constructed multi-regime field on
      the NACA coarse mesh (16.4k elements, all regimes populated); **5.7e-9
      on the real converged G4.1 M0.80 field** (supersonic pocket 2166
      elements = 1189 accelerating + 977 shock-point, M_max 1.3729 — the P8
      Newton target state). **Sign finding:** the FD gate arbitrated the
      design.md §6.3 ∂μ/∂φ chain — `dμ/dM² = +M_c²/M⁴`; the doc's earlier "−"
      was a transcription typo, corrected 2026-07-10.
      **Kink caveat (recorded for P8):** ρ̃ is C⁰ (not C¹) exactly AT the
      max(ν_e, ν_u) tie and at the switch threshold M² = M_c² — the
      measure-zero locus of design.md §3.1. An FD probe straddling a kink
      returns a branch average by construction, so tests/demo exclude the
      kink's ε-neighbourhood (0.04 % of elements on generic fields, 2/16.4k on
      the converged field). Measured trap worth remembering:
      symmetry-degenerate fields (separable φ on structured/prism-split
      meshes) park whole element slabs exactly on the tie — the FD then reads
      ~1e-5 errors that are branch averages, not derivative bugs
      (`test_p7_diff_flux.py` docstrings).
**Demo:** `cases/demo/p7_diff_flux/` — 7/7 PASS incl. the gated converged-field
part (V7.1 FD scatter/histogram, V7.2 regime map, V7.3 regimes over the real
supersonic pocket, V7.4 converged-field FD accuracy; committed PNG/CSV).
**Effort:** 1 session, as planned. The kernel flux stays opt-in with its
calibration deferred (a P8 measurement item, not on the Newton path).

### P8 — Performance & robustness: fully-coupled Newton
**Deliverables:** exact Newton Jacobian (design.md §6.3 / López Eq. 3.24 +
Appendix B) with widened sparsity, **fully-coupled Newton with Γ as an
unknown** (design.md §8.1), Mach continuation + load stepping, GMRES + AMG
path, AMG setup reuse, Eisenstat–Walker inexact-solve schedule, profiling
report. Consumes the **P7** differentiable flux: the exact Jacobian needs
∂ρ̃/∂φ well-defined at frozen selection — the shipped **walk** flux already
provides it (P7; the walk is kept, its ∂ρ̃/∂φ FD-verified), so the default
Newton flux is the sparse frozen walk. The streamline-Gaussian kernel is an
optional denser-Jacobian alternative, promoted only if N2 measures it favourable.

**Design-pass decision (2026-07-08, design.md §8.1 — López dissertation
Ch.3–4 + Appendix B, verified against the PDF; supersedes the earlier
internal note where it conflicts).**
- **Full Jacobian ("strategy A"), not an approximate one.** López keeps the
  switching-function derivative ∂μ/∂φ (B.3/B.6) and the nonzero upstream
  coupling (B.4) → strict quadratic convergence (Tables 4.5/4.6/4.9). Dropping
  Term 3 → only superlinear. The upstream coupling widens the stencil; recolor
  and extend the CSR scatter map accordingly. **Correction (2026-07-08): the
  stencil is NOT "one element-layer wider" in UP3D** — that is López's
  single-hop result; UP3D's upstream is multi-hop (sliver tets, §3). The width
  depends on which P7 flux P8 differentiates (design.md §6.3):
  - **default — frozen walk** → Term 3 adds ~+1 upstream element/row (sparse,
    closest to López), but at graph-distance ≤4 (long-range coloring/CSR edges);
    ∂ρ̃/∂φ FD-verified in P7, `max_ε` optional for churn robustness;
  - *optional — streamline-Gaussian kernel* → Term 3 couples each element to its
    **whole depth-3 BFS neighbourhood** → a materially **denser** Jacobian
    (smoother, frozen-selection-free), only if N2 measures it net-favourable.
  N2 **measures** nnz/GMRES-iter/AMG-setup before committing the coloring/CSR
  rebuild; the old "~30 % memory" estimate is retired.
- **Fully-coupled (φ_red, Γ), not the Γ-secant.** The P5 medium instability was
  secant–density coupling; solving Γ inside the Newton system removes the
  secant and the instability. The Γ-Jacobian blocks are nearly in the code
  already: `∂R_red/∂Γ_j = TᵀJ g_j` (≈ `wake.py::self._h` at the Picard level),
  and `F_j = kutta_targets_j − Γ_j` with a sparse ±1/n_j `∂F/∂φ_red`. **Do not
  omit the far-field column:** Γ also enters the vortex Dirichlet data
  (`dirichlet.py`), adding −A_coupling·(∂vals_red/∂Γ_j) — this term is folded
  silently into the Picard RHS and is easy to miss under Newton. The N-Γ split
  (keep secant, Newton only the density inner) is a **fallback only** (retains
  the P5 risk).
- **Load stepping — two PDF corrections to the earlier note.** (i) Within a
  case, López holds M_crit and μ_c **fixed** and ramps only M∞ (Tables 4.7/4.8);
  harder cases get lower M_crit / higher μ_c for the whole ramp — there is **no**
  per-step M_crit 0.99→0.90 sweep. (ii) ONERA M6 (Table 4.13): 12 steps, M_crit
  **constant 0.95**, μ_c held at 2.0 while M∞ ramps 0.50→0.84, then μ_c
  **decreased** 2.0→1.6 at fixed M∞ to sharpen the shock — μ_c scheduling is a
  post-target dissipation reduction, not a during-ramp increase.
- **Expected cost (López §4.7):** subsonic 5–10 Newton iters total; transonic
  4–9 per load step; ONERA M6 ≈ 12 × 5–9 ≈ 60–110 total (vs ~10⁴ Picard today).
  CL reference for M6 is 0.288 (KRATOS = Tranair, Table 4.15) — UP3D's current
  0.245 is a separate sharp-TE/P1 accuracy gap tracked via V6 → P9 (curved
  elements), not a Newton target.

**Sub-phase order** (Newton needs the P7 flux first; N0 optional):
- **N0 (optional)** smooth density clamp — only if a Newton transient stalls on
  the clamp (design.md §3.2); skip by default (converged states are clamp-free).
- **N1** = the P7 frozen-selection differentiable **walk** flux (∂ρ̃/∂φ
  FD-verified; `max_ε` optional). The streamline-Gaussian kernel is the optional
  denser-Jacobian alternative, not the default N1 deliverable.
- **N2 ✓ (2026-07-10)** Newton Jacobian assembly (Term 1 reuse + Term 2 local + Term 3 upstream)
  with finite-difference verification and the G4.2 subcritical bit-identity.
  **Default = the sparse frozen walk** (~+1 upstream element/row); if the kernel
  is being considered, N2 first **measures** its Jacobian nnz, GMRES iterations,
  and AMG setup vs the walk before committing any wider coloring/CSR rebuild.
  *Landed: Terms 1+2 fused on the shared Picard pattern, Term 3 as active-set
  COO rebuilt per Newton step (no recolor/wider-CSR machinery needed — measured
  ms-scale at pocket sizes); JVP-vs-FD rel ~1e-10; forward paths byte-untouched.*
- **N3 ✓ (2026-07-10)** GMRES + AMG (aggregation on the symmetric part) in `solve/linear.py`.
- **N4 ✓ (2026-07-10)** fully-coupled Newton driver `solve/newton.py` (Γ Jacobian incl. the
  far-field column) + subsonic verification (cl matches P3 < 0.5 %).
  *Landed: exact δΓ elimination, far-field column FD-guarded, Γ/cl match P3 to
  1.9e-8 / ~1e-7, terminal orders [2.57, 1.79] on the cold start; details in the
  ledger.*
- **N5 ✓ (2026-07-11)** transonic + load stepping → **G8.1** quadratic
  convergence (case set re-specced to the solvable domain, see the gate).
  *Landed: the transonic robustness chain — `precond="direct"` exact steps
  (splu + Woodbury on the rank-n_st coupling; the shock-position soft mode
  leaves η-accurate Krylov steps stalled — measured frozen-system residual
  flat at 3.7e-6 with GMRES converging to η), stall-adaptive freeze of the
  upwind assignment (`UpwindOperator.freeze_upwind_state` + frozen sweeps;
  the 2.5D prism family parks ~1e3 elements in the max(ν_e,ν_u) near-tie
  band — the P7 kink trap in Newton form — and live Newton limit-cycles
  there), active-set refresh with two-cycle acceptance + honest
  `residual_unfrozen` floor reporting, freeze-revert/level-fail-fast/
  best-of-tried line-search safety nets. Recipe:
  `tests/test_p8_newton.py::NEWTON_TRANSONIC_RECIPE` (dm 0.025, dm_min
  0.003, freeze_tol 1e-6, refresh 8, direct, budget 60). Medium G8.1 run
  ~100 s vs Picard's 16m39s-to-a-non-solution.*
- **N6 ✓ (2026-07-11)** ONERA M6 + performance → **G8.2/G8.3** closed.
  *Landed: (1) `direct_refactor_every`/`direct_reuse_rtol` in
  `solve_newton_lifting` — the lagged-LU direct mode: refactor every k-th
  step, in between GMRES on the FRESH coupled operator preconditioned by
  the stale LU at tight rtol 1e-8 (orders below the soft-mode-stall η;
  affordable only because the stale LU is near-exact), GMRES failure falls
  back to refactor + exact Woodbury in the same iteration — robustness
  never below the every-step-direct path; default 1 = bit-identical N5
  behavior (suite-locked). Why: true-3D LU fill makes splu ~100× the thin
  2.5D cost at equal dofs (measured 18.6 s vs ~0.2 s at ~6e4 dofs; the
  every-step-direct M6 medium run was 1606 s with 97% in splu; profile in
  the N6 session). (2) `NEWTON_M6_RECIPE`
  (tests/test_p8_newton.py): N5 chain + refactor_every 1000 (once per
  level + on-demand) + the P5 dm=0.05 schedule (M6 is far from the
  NACA-medium fold — all levels converge with zero dm halving) +
  farfield_spanwise_gamma. (3) **M0.84/α3.06 has a reachable isolated
  Newton solution on BOTH M6 meshes** (the plan's fold contingency never
  triggered): coarse 42 s / medium 249 s end to end, terminal quadratic
  in the frozen phases (medium final level …2.6e-7→2.1e-10→7.0e-15),
  0 limited/floored, coupled Kutta |F| ~2e-16. (4) **P5-caveat measured**
  (the agent-rules recorded item): Newton residual at the committed P5
  Picard states = 8.6e-6 coarse / 7.6e-6 medium (|F| 5.3e-4/5.8e-4 — far
  better than the P4 stall's 2.2e-4, but not solutions either); the true
  solutions sit at **cl_p 0.2560 coarse (+5.8%) / 0.2646 medium (+7.9%)**
  with shock positions essentially unchanged (0.596/0.541/0.362 vs P5's
  0.594/0.526/0.345) and M_max 2.13 vs 1.995; cl_KJ 0.2692 — the
  inviscid-vs-Tranair(0.288) gap assigned to P9 narrows 0.043→0.019.
  Under-convergence lives in the circulation/lift, not the shocks.
  Solution-identity note: coarse is continuation-path independent to
  1e-15; medium's two paths (dm 0.025 vs 0.05) agree to |dφ| 2.6e-4,
  |dΓ| 3.7e-5, cl to 5 digits — the N5 assignment-discontinuity floor
  semantics (live floor ~5e-5, stale 84 vs 92), covered by the lock
  bands.*

**Gates:**
- [x] G8.1 **(re-specced 2026-07-11, user-approved — case set moved to the
      solvable domain)** Newton terminal quadratic convergence on **coarse
      M0.80/α1.25 AND medium M0.7875/α1.25** (consecutive ≥1.5-digit residual
      collapses ending < 1e-9; measured coarse …6.99e-6→1.47e-7→8.2e-11, medium
      frozen phases to 2.6e-13, cf. López Table 4.9); Term-2/Term-3 Jacobian
      finite-difference-verified (rel err < 1e-6) incl. the supersonic pocket
      (gated `test_jacobian_fd_converged_pocket` on the converged Newton coarse
      M0.80 field, rel ~1e-10); G4.2 subcritical no-op still bit-identical to
      the P3 path (suite-locked). **Why the original "G4.1 case (M0.80 medium)"
      was moved:** the Newton runs exposed that (i) the P4 Picard states are
      NOT discrete solutions (Newton residual 2.2e-4 at the committed coarse
      state — see the P4 ledger erratum) and (ii) on the medium mesh the true
      solution family steepens into the FP non-uniqueness fold (cl 0.396 at
      M0.775 → 0.523 at M0.7875, measured at 1e-13/8e-11 residuals; no
      reachable isolated solution at M0.80, M_max ≈ 1.45 beyond the isentropic
      validity envelope — design.md §12 risks 2/3). Physics acceptance is a
      REGRESSION LOCK around the measured Newton solutions (coarse shock
      0.658/cl 0.459; medium 0.674/0.523) — the Euler-anchored G4.1 band does
      not bind conservative FP at this shock strength (Holst PAS 2000; the
      reference CSV is untouched per hard rule 6).
      Closed 2026-07-11: `tests/test_p8_newton.py` gated ×2 +
      `tests/test_p8_jacobian.py` gated FD, demo `cases/demo/p8_newton/`.
- [x] G8.2 ONERA M6 medium mesh < 5 min single node, end to end — closed
      2026-07-11: **249.2 s** (mesh+cut 7.3 s / solve 239.8 s / forces+
      3-section shocks 2.1 s; P5 Picard was 4539 s to a state with Newton
      residual 7.6e-6). Timing protocol: the CLAUDE.md 16-thread cap must
      cover BLAS/OMP too — `NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
      OPENBLAS_NUM_THREADS=16` (A/B measured 2026-07-11: without the
      BLAS/OMP caps the identical run is ~333 s, oversubscription on the
      16C/32T box costs ~33%; gated test 252 s with caps). Physics regression locks
      around the measured Newton solution: cl_p 0.2646 ± 0.005, shocks
      η44/65/90 = 0.596/0.541/0.362 ± 0.02, M_max 2.134 ± 0.05. Gated
      `test_g82_m6_medium_newton_end_to_end` + demo part 3.
- [x] G8.3 full regression suite runtime < 10 min (CI budget) — closed
      2026-07-11: **301.66 s (5m02s)**, 182 passed + 8 skipped + 2 xfailed
      (all heavy transonic/M6 gates behind PYFP3D_TRANSONIC_GATES=1).
**Demo:** `cases/demo/p8_newton/` — Newton convergence order + runtime
breakdown (V8.1a/b/c, V8.2 M6 breakdown, G8.2 CSV).
**Visual test examples:**
- V8.1 Newton convergence plot of ||R|| and estimated convergence order p_k; confirm terminal quadratic regime (p_k near 2).
- V8.2 runtime breakdown chart (assembly/linear solve/postprocess) for ONERA M6 medium mesh; verify bottleneck location matches profiling report.
- V8.3 regression dashboard trend over commits (runtime + key errors); confirm speedup does not degrade physical metrics.
**Effort:** 3–5 sessions. **Risk:** medium.

### P9 — Grid-convergence & accuracy-gap discrimination (evidence phase; NEW 2026-07-11)

> **Track-P renumbers 2026-07-11 (user-directed, two same-day insertions):**
> this evidence phase is inserted as P9, and a second insertion the same day
> added P10 (Newton generality & continuation efficiency). NET mapping from
> the pre-2026-07-11 IDs: curved walls P9 → **P11**, backlog P10 → **P12**.
> Docs written before 2026-07-11 use the old IDs (as with the 2026-07-08
> renumber); the intermediate one-commit state (curved walls = P10) appears
> only in commit 467737e.

**Motivation.** Before committing to the large curved-walls effort (now P11),
discriminate what the remaining external lift gap actually is. The P8
capability demo (`cases/demo/p8_capability/`, 2026-07-11) showed the
attribution of the M6 gap (cl_KJ 0.2692 vs Tranair/KRATOS 0.288, i.e. 0.019)
to the sharp-TE/LE flat-facet P1 wall is **inference-grade, not oracle-grade**
— unlike G1.6, where the exact-potential recovery oracle + the nodal-φ
convergence-order degradation + the refinement-saturation sweep nail the
variational crime. Two demo data points argue for caution: (a) 2D sharp-TE
subsonic lift converges cleanly to the corrected-panel bracket (coarse −2.7%
→ medium −0.3% of the PG–KT midpoint), and (b) M6 cl_KJ is still RISING under
refinement (0.2621 → 0.2692) with no sign of the sphere-style saturation.
If the gap Richardson-extrapolates away, P11's lift case collapses to
"improves the constant"; the sphere-Cp (G1.6) case for P11 stands regardless.

**Gates (decision bands PRE-REGISTERED here to forbid post-hoc fitting):**
- [x] G9.1 **M6 three-point grid study** at M0.84/α3.06: coarse/medium/fine
      Newton (`NEWTON_M6_RECIPE`; fine = the M1 2513k-tet mesh, gitignored,
      regenerated on demand; fine solve robustness/perf at ~450k dofs is the
      phase's technical unknown — lagged-LU one factorization per level first,
      `precond="amg"` fallback; budget ≤ ~2 h, run ONCE, npz cached locally
      like P5). Richardson-extrapolate cl_KJ and cl_p with the measured
      per-level h-ratio from tet counts (~1.85/1.93); validity requires
      monotone three-point sequences. **Verdict bands on extrapolated
      cl_KJ∞:** ≥ 0.283 ⇒ resolution-dominated (gap closable by refinement);
      ≤ 0.278 ⇒ saturating floor confirmed (P11 confirmed for lift);
      in between ⇒ inconclusive, report as-is. Artifact:
      `cases/demo/p9_grid_discrimination/results/` (three-point values,
      observed order p, extrapolants, verdict CSV + convergence PNG).
      **CLOSED 2026-07-11 — outcome: the pre-registered bands CANNOT FIRE,
      and the reason is the phase's real finding.** (a) **Solver path
      (recorded, feeds P10):** `precond="direct"` does NOT scale — a SINGLE
      `splu` at ~450k dofs ran **4 h 39 min without returning** (RSS 26 GB,
      killed) vs 18.6 s on medium; the pre-registered `precond="amg"`
      fallback works and is FASTER than direct at every size **once the
      Eisenstat–Walker forcing is tightened to η = 1e-8** (validated against
      the G8.2 locks on medium before spending the fine budget: 66 s vs
      141 s direct, same solution to 4 digits, 0 GMRES stalls; coarse 8 s vs
      42 s). N5's "Krylov steps stall on the shock-position soft mode" is a
      property of the LOOSE default forcing, not of AMG. (b) **Seeding
      (recorded):** the fine mesh's cold Picard-5 seed OVERSHOOTS the LE into
      the density floor (M0.70 cold: 4036 limited / 1847 floored, level-0
      Newton stalls at |R|~6e-2 and cannot dm-halve ⇒ the whole solve breaks;
      M0.50 cold: still 658 floored). Fixed at continuation-path level only:
      m_start 0.30 + n_picard_seed 12 ⇒ clean 0/0 start, ramp carries it up
      (13 levels, 5294 s = 88 min, inside the 2 h budget). (c) **★ THE
      FINDING — the fine mesh is NOT a discrete solution, because a
      singularity DIVERGES under refinement:** max local Mach (unlimited)
      **1.40 → 2.13 → 7.93** across coarse/medium/fine, with **0 / 0 / 9**
      cells crossing the M_cap = 3 speed limiter. Permanently-limited cells
      block the N5 freeze machinery (which requires 0 limited), so the fine
      Newton limit-cycles at |R| ~ 1e-5 for the full budget; its cl_KJ
      (0.2393) is a limit-cycle artifact, NOT a lift. With only two discrete
      solutions (coarse 0.2621, medium 0.2692) there is no three-point
      Richardson ⇒ **verdict INVALID; the ≥0.283 / ≤0.278 bands cannot
      fire.** (d) **★ WHERE the singularity sits (this redirects P11):** all
      9 capped cells are at **z/b = 0.998–1.000** (the tip), **x − x_TE =
      +0.002…+0.017** (AFT of the trailing edge, i.e. NOT on the wing), and
      **|y| < 0.003** (in the chord plane) — i.e. **on the WAKE SHEET at its
      free tip edge**, exactly where `Γ(tip) = 0` is enforced (M1: tip free
      edges stay single-valued). This is the classic **vortex-sheet-edge
      singularity** of a rigid planar wake (the real flow rolls up into a tip
      vortex). The **driver is the trailing vorticity dΓ/dz, NOT the bound
      circulation Γ**: Γ→0 correctly at the tip (Γ(tip)=0 is a
      necessary-not-sufficient regularity condition), but the *unloading rate*
      |dΓ/dz| is largest at the tip (measured ~10× mid-span on B7's smooth Γ(z)),
      and a terminating flat vortex sheet has a free-edge crossflow singularity
      like a flat-plate edge — **1/√r (flat-plate-edge type), NOT 1/r** (a 1/r
      concentrated line vortex would give refinement exponent p=1; P13/G13.1
      measures the conforming peak-Mach exponent **p≈0.59**, `cases/demo/
      p13_tip_edge_singularity/`). P5's "bounded tip-TE-corner P1 overshoot" (M_max 1.995 medium, recorded
      as the only surviving singularity trace) is the SAME object seen at a
      resolution too coarse to reveal that it is unbounded. **It is a WAKE
      MODEL defect, not a wall-element defect — curved wall elements (P11) do
      not remove a wake-sheet edge.** (The 9 cells' `|y| < 0.003` chord-plane
      locus, unpersisted in the P9 census CSV, was independently re-measured
      2026-07-12 from the fine cache: |y| 0.0009–0.0027, all 9 — confirmed.)
      ★ **The "fix" is roll-up / a tip vortex, NOT Track B as such** (correction
      2026-07-12; the earlier "Track B / tip treatment" phrasing over-promised).
      Track B changes the wake **representation, not the model** — it keeps the
      same rigid planar sheet ending at the tip — so it does not remove this
      singularity. Measured decisively by `cases/demo/p13_tip_edge_singularity/`
      (subsonic M0.5, NO limiter — the clean geometric probe): the tip-edge peak
      Mach **diverges under refinement on BOTH the conforming and level-set
      paths** (same mesh coarse→medium: ×1.38 conforming, ×2.28 level-set —
      ★ ERRATUM 2026-07-14: the ×2.28 LS reading is the `element_mach2`
      mixed-plain ×5 metric artifact retired at the B8 re-spec (honest LS
      exponent +0.62 ≈ conforming +0.52; both still diverge) — the
      LS peak sitting in its +2.9% straddler cells, i.e. no blunting), while the
      wing control stays flat. The only fix is wake roll-up / an explicit tip
      vortex; **no current Track B phase does that** (B10 free-wake is SHELVED and
      is about O(θ²) deflection, not roll-up).
- [x] G9.2 **2.5D sharp-TE lift oracle**: NACA0012 fine mesh (generated on
      demand), subsonic M0.5/α2 Newton, coarse/medium/fine vs the
      corrected-panel reference — PASSES iff the |error| vs the PG–KT midpoint
      decreases monotonically coarse→medium→fine AND fine is within ±1%.
      A pass certifies a sharp TE alone imposes NO lift floor at these
      resolutions in 2D — the V6 mechanism would then have to be 3D-specific
      (swept-LE suction resolution, tip, planar-wake model), sharpening what
      P11 must actually fix. Artifact: same demo results/.
      **CLOSED 2026-07-11 — PASS, cleanly:** |error| vs the PG–KT midpoint
      **2.71% → 0.33% → 0.03%** (coarse/medium/fine, all three converged to
      |R| ≤ 4.4e-11), fine well inside the ±1% clause. **A sharp trailing
      edge imposes NO lift floor at these resolutions** — the 2D sharp TE is
      EXONERATED, and with it the "sharp-TE/LE P1 wall gradient" story for
      the 3D gap loses its 2D leg.
- [x] G9.3 **attribution verdict recorded**: split the 0.019 gap into
      resolution / geometry-floor / unattributed shares from G9.1 + G9.2,
      write the P11 (curved walls) go/no-go recommendation into the progress ledger +
      demo_report (user arbitrates the decision).
      **RECORDED 2026-07-11 — the gap is UNSPLITTABLE as posed, and the
      question turned out to be mis-posed.** (1) The 0.019 gap CANNOT be
      split by grid convergence: the 3D sequence does not converge, because
      the fine mesh is not a discrete solution (G9.1c). No resolution /
      floor shares can be reported honestly — `verdict.csv` records both as
      `n/a`, not as numbers. (2) What the evidence DOES say: the 2D sharp TE
      is exonerated (G9.2, error → 0.03%), and the 3D blocker is a
      **divergent vortex-sheet-edge singularity on the rigid planar wake at
      the tip** (G9.1d) — a WAKE-MODEL object, not a wall-geometry one.
      (3) **Recommendation to the user (arbitration pending):** P11 (curved
      wall elements) is **NOT supported by P9 as the fix for the 3D lift
      gap** — its remaining justification is G1.6 (sphere-Cp on a SMOOTH
      curved wall), which is a different, still-valid mechanism. The 3D
      accuracy path now points at the **tip/wake treatment: wake roll-up or
      an explicit tip-vortex model**, which must land before ANY 3D
      grid-convergence claim is possible. ★ **NB (correction 2026-07-12):
      Track B's level-set wake is NOT that fix** — it changes the wake
      representation, not the rigid-planar-sheet model, and
      `cases/demo/p13_tip_edge_singularity/` measured the edge singularity
      diverging on the level-set path too (B10 free-wake, the only roll-up-
      adjacent phase, is shelved and is about O(θ²) deflection). Suggested
      order: decide the tip/wake route (a genuinely new wake model) first;
      keep P11 scoped to G1.6 only.

**Non-goals:** no solver/numerics changes; no reference-data edits (hard rule
6); fine meshes stay gitignored; heavy runs are one-shot demo scripts with
committed CSV/PNG evidence, NOT suite tests (G8.3 CI budget untouched).
**Demo:** `cases/demo/p9_grid_discrimination/`. **Effort:** 1–2 sessions.
**Risk:** low-medium (the ~450k-dof fine Newton solve is the unknown).

### P10 — Newton generality & continuation efficiency (NEW 2026-07-11, user-proposed; NOT started)

**Motivation.** Two solver-level items raised by the user after reviewing the
P8 capability demo. (1) The P8 Newton driver has NO non-lifting entry —
`solve_newton_lifting` structurally requires the wake cut + Kutta/Γ block, so
sphere-class bodies run Picard only (the boundary recorded in demo_report
§capability; promoted here out of the backlog). (2) Intermediate
Mach-continuation levels are converged far deeper than their warm-start role
requires. **Framing correction recorded with the proposal:** in the coupled
P8 Newton there is NO Kutta outer loop — Γ is updated EVERY Newton step
(δΓ = K δφ + F) and ‖F‖ sits at machine zero throughout (visible in the
capability-demo convergence figures); the two-level structure in those plots
is the MACH RAMP. The user's efficiency intuition ("don't converge the
residual deeply between the big circulation moves; tighten only once the
circulation has settled") therefore lands on the continuation levels: an
intermediate level only seeds the next level, yet today it runs to tol 1e-10
plus the freeze/refresh polish. Estimated from the committed capability-demo
data: M6 medium spends ≈ 3–4 steps × ~5.1 s on each intermediate level's
terminal tail (~20–25% of the 240 s solve); the NACA fold-zone medium run
spends most of its 252 steps in intermediate-level churn.

**Gates:**
- [ ] G10.1 **non-lifting Newton entry**: `solve_newton` on the RAW (uncut)
      mesh — `NewtonWorkspace` generalized to `wc=None` (no WakeConstraint /
      Kutta row / Γ unknown; reduction = identity; far field without the
      vortex basis). Acceptance: sphere_shell medium at M0.3 matches
      `solve_subsonic` (same discretization) to |Δφ|∞ < 1e-8 and peak-Cp
      diff < 1e-6, terminal quadratic order observed, m_inf = 0 converges in
      ONE step ≡ `solve_laplace`; suite green with additive tests only.
      NOT an accuracy gate: G1.6 (sphere-Cp 11.6%) is untouched — same
      discretization, same open xfail.
- [x] G10.2 **level-adaptive intermediate tolerance** in
      `solve_newton_transonic`: new opt-in knob (default None = current
      behaviour BIT-IDENTICAL, suite-locked). Candidate rule (the user's
      proposal made precise): accept an intermediate level once
      (relative residual drop ≥ 1e3 from the level start OR ‖R‖ ≤ 1e-5 OR
      the live-stall detector fires) AND 0 limited/floored AND ‖F‖ within
      tol — i.e. treat intermediate-level stall as "advance to the next Mach
      level" instead of "freeze and polish"; the FINAL level (and any
      `upwind_c_post` stages) keeps tol 1e-10 + the full freeze/honesty
      machinery. Acceptance = committed A/B on (a) the G8.2 M6 medium recipe
      and (b) the G8.1 NACA medium recipe: ALL regression locks intact,
      robustness not degraded (level count and dm-halvings not worse),
      end-to-end deltas reported; promote into the recipes only if (a)
      improves ≥ 15% — otherwise record the negative result and keep the
      default. **Known risk to test explicitly (P8 trap):** the N2–N4 record
      warns "warm-start only from CONVERGED levels" — a loosely-converged
      seed is a new, unvalidated category and may cost more downstream steps
      than it saves; that is exactly what the A/B decides.
      **CLOSED 2026-07-11 with a SPLIT A/B verdict**
      (`cases/demo/p10_newton_usability/`, ledger entry): (a) M6 medium
      locks intact + identical final level, **239.5→140.3 s (+41.4%)** ⇒
      `intermediate_tol=1e-5` promoted into `NEWTON_M6_RECIPE`;
      (b) NACA medium M0.7875 = the pre-registered trap MEASURED — loose
      intermediates leave the fold-zone ramp with an untracked Γ seed and
      the run ends unconverged (cl 0.369 vs lock 0.523) even after two
      hardenings (loose acceptance requires ≥1 Newton step; dm-halving
      retry levels run STRICT — both kept in the shipped rule), so
      `NEWTON_TRANSONIC_RECIPE` keeps the strict default: **loose
      intermediates are contraindicated near the fold.**
- [x] G10.3 **no-ramp direct-solve feasibility** (NEW 2026-07-11,
      user-directed after reviewing the G10.2 A/B; scheduled BEFORE P9):
      can `solve_newton_lifting` converge AT the target M∞ from a cold
      start — no Mach continuation at all? The G10.2 data motivates the
      question: on M6 medium the loose intermediate levels contribute only
      1–3 Newton steps each yet still cost ≈ 89 s of the 140.3 s adaptive
      solve (one splu + seed per level), so the remaining value of the
      ramp is exactly (i) FP branch selection near the non-uniqueness
      fold (design.md §12 risk 2 — the reason the ramp exists) and
      (ii) the Γ/active-set seed quality. **Protocol (evidence study,
      one-shot demo script, no default-path changes):** run the
      single-level Newton at target M∞ with the current machinery
      (precond direct + lagged LU, freeze_tol 1e-6, budget 60) under two
      seedings — (s1) the standard short Picard seed AT the target M∞
      (n_picard_seed 5, today's subsonic-entry default) and (s2) a
      deeper Picard seed (n_picard_seed ~40, the P5-era warm start) — on
      four cases: M6 coarse + medium M0.84/α3.06 (far from fold; the
      cases where skipping the ramp could pay) and NACA coarse M0.80 +
      medium M0.7875 α1.25 (fold zone; the branch-selection risk cases).
      **Pre-registered outcome classes per case:** (A) SAME-solution
      convergence — terminal quadratic, 0 lim/flr, ALL of that case's
      ramp-solution regression locks met ⇒ no-ramp is valid there;
      (B) WRONG-branch convergence — converges cleanly but outside the
      locks ⇒ hard evidence that the ramp is doing branch selection,
      record cl/shock of the off-branch solution; (C) no convergence
      (stall/divergence/clamped) ⇒ negative result. **Promotion rule:**
      adopt no-ramp into `NEWTON_M6_RECIPE` only if BOTH M6 cases land
      in class (A) with end-to-end ≥ 20% below the current post-G10.2
      recipe; the Mach ramp remains the shipped default for the 2.5D
      fold-zone recipe REGARDLESS of outcome (G10.2's measured
      contraindication already covers it), and remains available as the
      fallback everywhere. **Traps to test explicitly:** the Picard seed
      itself at a supercritical M∞ may need continuation (P4/P5 record)
      — report the seed's limited/floored counts; the N5 freeze-cap
      guard must not lock a mid-transient garbage assignment (report
      freeze reverts); a class-(A)-looking result with nonzero
      lim/flr at any point of the path still counts as (C) for
      promotion. Artifact: `cases/demo/p10_newton_usability/`
      `run_g103_noramp.py` + committed CSV/PNG (per-case class, steps,
      wall, lock deltas).
      **CLOSED 2026-07-11 — verdict: KEEP the Mach ramp** (all 8 runs in
      `results/g103_noramp.csv`; `clamp_history` instrumentation added,
      additive key). Measured: **far from the fold BOTH seedings are
      class (A)** — no-ramp M6 converges to the SAME solution (cl
      0.2646, shocks/M_max equal to the ramp locks), so branch selection
      is NOT the binding constraint there; but s1 (Picard-5) transits
      CLAMPED states (peak 45 lim+flr on medium; final 0/0) at
      141.4→79.0 s (+44%), failing the clamp-free clause, while s2
      (Picard-40, clamp-free) gains only +6.6% medium / −37% coarse,
      failing ≥ 20% — no promotion path satisfies the pre-registered
      rule. **Fold zone: class (C) on medium under BOTH seedings**
      (s1 stalls at 4.6e-6 with cl 0.449 ≠ 0.523; s2's un-continued
      Picard seed itself diverges into clamp land, 9934 limited —
      matching the P4/P5 "Picard needs continuation too" record; the
      deep seed is actively HARMFUL without a ramp). One recorded
      curiosity: NACA coarse M0.80 s1 is class (A) in 4.0 s (clamped
      transient + 1 freeze revert) — a single-case exception that does
      not change the rule. The +44%-but-clamped s1 observation is left
      for the USER to arbitrate if the clamp-free clause is ever to be
      relaxed; until then the ramp stays everywhere.

**Sequencing note (updated 2026-07-11):** user-directed order is
G10.2 (✓) → **G10.3 → then P9** — both G10.2 and a positive G10.3 cut the
G9.1 fine-mesh (~450k-dof) run cost; G10.1 has no ordering constraint and
may land whenever convenient. **Demo:** `cases/demo/p10_newton_usability/`
(G10.2 A/B + G10.3 no-ramp study both committed). **Effort:** 2–3 sessions.
**Risk:** low (G10.1) / medium (G10.2 — robustness interplay with the freeze
machinery, measured; G10.3 — branch selection without continuation is the
known FP hazard, hence the pre-registered class-(B) detector).

### P11 — Curved / isoparametric wall elements ★ accuracy phase (renumbered from P9 via P10 on 2026-07-11; opening CONDITIONAL on the G9.3 verdict)
**Motivation.** The shared accuracy route for two open items that a flat-facet
P1 wall cannot close: (i) **G1.6** incompressible sphere-Cp < 2% (root-caused to
the natural BC enforced on the flat polyhedral wall instead of the true curved
surface — design.md §5.1/§5.1.2, DP1 "> 5%" branch; boundary-data corrections
ruled out with evidence), and (ii) the **residual V6 < 1%** floor left after P6
removes the sawtooth — the remaining CL_p-below-CL_KJ gap is attributed to the
sharp-TE/LE P1 wall-gradient error (design.md §9.1; M6 cl_KJ 0.2692 vs the
Tranair/KRATOS 0.288 after the P8 Newton re-measurement) — **P9 (G9.1/G9.2)
tests this attribution before this phase opens**.
**Deliverables:** curved/isoparametric wall boundary elements (or the
geometry-consistent Option-C re-spec of G1.6, design.md §5.1 Option C), scoped
so the wall integral/gradient sees the true surface. **Gates:** G11.1 G1.6
sphere-Cp < 2% (medium); G11.2 V6 < 1% on the M6 medium gate. **Demo:**
`cases/demo/p11_curved_walls/` (future). **Effort:** large (own effort per DP1).
**Risk:** medium-high — isoparametric assembly is a real extension.

### P12 — Backlog (renumbered from P10 via P11 on 2026-07-11; post-v1.0, ordered by value to the MDO thread)
1. Discrete adjoint = transpose of the P8 Newton Jacobian + Kutta/wake constraint
   adjoint terms → gradients for shape optimization.
2. VII hook: transpiration BC ∂φ/∂n = d(u_e δ*)/ds on wall faces (reuses the
   pyTSFoil IBL line of work).
3. Mixed prism/tet elements; parameter calibration campaign for (C, M_c, ω,
   ω_Γ) with the stochastic-Kriging BO framework.

(The non-lifting Newton entry, briefly listed here on 2026-07-11, is promoted
to P10/G10.1 the same day.)

### P13 — Tip / wake-edge singularity: characterization + wake-model fix ★ accuracy phase (NEW 2026-07-13, user-approved; appended, no renumber)
**Motivation (direct descendant of P9).** P9/G9.1 localized the 3D
grid-convergence blocker to a **divergent singularity on the rigid planar wake
at its free tip edge** (the fine M6 mesh is not a discrete solution). Two things
that P9 left loose and this phase closes/records: (1) the singularity is driven
by the **trailing vorticity dΓ/dz** (the *unloading rate*, largest at the tip:
~10× mid-span on B7's smooth Γ(z)), **not** the bound circulation Γ (which
correctly → 0 at the tip — Γ(tip)=0 is a necessary-not-sufficient regularity
condition); a terminating flat vortex sheet has a free-edge crossflow
singularity like a flat-plate edge, i.e. **1/√r**, not the "1/r-type" earlier
docs wrote (a 1/r line vortex would refine at exponent p=1). (2) The fix is a
**wake-MODEL change** (roll-up / explicit tip vortex), **not** a wall-element
change (P11 — curved walls do not remove a sheet edge) and **not** Track B as
such (which changes the wake *representation*, not the model — measured: the
edge singularity diverges on the level-set path too). Implementation of the fix
is handed to **Track B as a rescope of the shelved B10** (level-set naturally
represents a movable/curved sheet via `update_direction()` — from "O(θ²)
deflection" to "roll-up").
**Gates (PRE-REGISTERED, P9-style):**
- [x] **G13.1 — Characterization (CLOSED 2026-07-13).** Prove it is a genuine
      wake-MODEL singularity, not a format artifact. Demo
      `cases/demo/p13_tip_edge_singularity/` (10/10 PASS), probed SUBSONICALLY
      (M∞0.5 — no shock/artificial-density/limiter to confound the geometric
      edge signal). **Measured:** (a) conforming three-point (coarse/medium/fine
      = 55.5k/350.7k/2.51M tets) log-log exponent of tip-box peak Mach vs 1/h
      **p = 0.59 ∈ [0.4, 0.65] ⇒ 1/√r (flat-plate edge)**; (b) the tip-box
      **p95/mean stay flat** (0.573→0.562→0.525 / ~0.49) while the **peak
      diverges** 0.712→0.981→1.510 — a localized-edge signature; (c) the edge
      diverges on **both** paths (conforming ×1.38 and level-set ×2.28
      coarse→medium on the same mesh; LS exponent 1.34 ≥ conforming 0.52 —
      ★ **ERRATUM 2026-07-14 (B8 re-spec):** these LS numbers are a **×5
      metric artifact** of `element_mach2`'s mixed-plain "side" reading
      (beyond-tip mixed-side plain cells; elem 93977 side 1.532 vs main 0.309);
      the HONEST LS exponent is **+0.62** (+0.37 no-sliver) — the SAME object
      and magnitude as the conforming +0.52, so "LS ≥ conforming" is RETIRED
      while "diverges on both paths" STANDS; see the B8 re-spec block +
      `run_b8_termination_diagnosis.py`) ⇒
      MODEL, not representation; peak sits **aft of the TE** in the chord plane
      (z/b 0.999, x−x_TE +0.006, |y|<0.001) ⇒ not a wall feature; (d) ★ the
      conforming **fine M0.5 solve does NOT converge** (limited/floored cells,
      ~1.4k NaN) — the tip singularity trips the limiter even subsonically,
      the exact M0.5 analogue of G9.1's transonic finding. Per the pre-registered
      risk clause, a drifting/off-band p would still close G13.1 (the core
      claims — singularity exists, model-not-representation, aft-of-TE — do not
      depend on p being exactly 0.5); it landed cleanly at 0.59.
- [x] **G13.2 — tip-edge desingularization (CONFORMING path CLOSED 2026-07-13;
      level-set clause OPEN).** The fix is a **spanwise loading taper**:
      `Γ_eff(z) = F(z) · Γ_Kutta(z)` applied to the per-station Kutta target
      (`constraints/wake.py::tip_taper_factors`, `solve_newton_lifting(
      tip_taper=…)`; default `None` = **bit-identical**). Shipped model
      (**user-arbitrated 2026-07-13**): `form="vanish_smooth"` (smoothstep,
      **compact support**), **r_c = 0.05 · b_semi**.
      **★ THE MECHANISM IS DISCRETE, AND IT IS NOT ROLL-UP, NOT A VORTEX CORE,
      AND NOT THE CONTINUUM EDGE EXPONENT.** The solver never sees the continuum
      edge: it sees the **outermost TE station**, which retains circulation
      `Γ_last` and sheds it as a CONCENTRATED VORTEX over the last cell (free-edge
      nodes are single-valued, so the jump falls to 0 in one element), inducing
      `~ Γ_last/h`. With `Γ_last ~ h^q` the edge peak grows as `h^(q−1)`, i.e.
      **p ≈ 1 − q**, and the criterion is **q ≥ 1**. This predicts the baseline
      exactly (measured q = 0.44, since Γ~√u and u_last~h ⇒ Γ_last~√h; p_pred 0.56
      vs p_meas 0.52). ⇒ **design.md §4.1's "roll-up / cored tip vortex" framing is
      superseded** (2nd correction to that paragraph).
      **★ The taper is AMPLIFIED, not applied:** Γ is a fixed point of
      `Γ = F·Γ_Kutta(Γ)` and the Kutta map has slope **b ≈ 0.93** (P2), so
      `Γ/Γ* = F(1−b)/(1−F·b)` — a taper of 0.8 gives **0.21×**, not 0.8×
      (test-locked). This is why the measured `q ≈ 3.3` far exceeds the naive
      `s+½`, and why the taper MUST be kept compact.
      **Measured (M∞0.5, strict off-body edge box):** edge peak
      **0.712 → 0.981 → 1.510 (p = +0.592)** untapered vs
      **0.567 → 0.565 → 0.570 (p = +0.009)** tapered — the singularity is GONE
      across a 45× cell-count range. **★ The M6 FINE mesh is now a GENUINE
      DISCRETE SOLUTION** (converged, **0 limited / 0 floored / 0 NaN**), which is
      exactly what G9.1 could not achieve.
      **★★ TRANSONIC M0.84 — G9.1's EXACT SCENARIO, CLEARED.** Where G9.1's fine
      mesh had **unlimited M_max 7.93 with 9 cells over M_cap=3, permanently
      limited, Newton limit-cycling, and a cl_KJ (0.2393) it had to declare an
      ARTIFACT**, the tapered fine solve is a **genuine discrete solution**:
      **M_max 2.818, 0 cells over M_cap, 0 limited, 0 floored, converged**
      (44.6 min).
      > ✓ **EVIDENCE RESTORED 2026-07-13 (was PROSE ONLY; now REPRODUCED with a
      > committed artifact).** The audit found this whole M0.84 paragraph existed
      > only as text — no script, no CSV, no cached solve — while a P11 ledger
      > status had been changed on its strength. It was re-run from scratch and
      > **reproduces to 4 digits**: demo
      > `cases/demo/p13_tip_edge_singularity/run_g132_transonic.py` (5/5 PASS),
      > artifact `results/g132_transonic.csv`. Measured coarse/medium/fine:
      > cl_KJ **0.2593 / 0.2652 / 0.2866**, cl_p 0.2534 / 0.2608 / 0.2835,
      > M_max **1.394 / 1.725 / 2.818**, all **0 over M_cap / 0 limited /
      > 0 floored / converged**. The census is G9.1's own (unlimited Mach field +
      > M_cap count), so it is a strict A/B; the only change is `tip_taper`.
      > ★ The rerun exposed WHY the evidence was lost: the first attempt reused
      > the medium recipe's `precond="direct"`, which at the fine mesh's ~450k
      > dofs is P9's documented 4h39m/26GB splu trap (killed at 1h16m/24GB). The
      > demo now uses P9's validated `precond="amg"` + tight EW forcing
      > (η=1e-8) + the fine seed guard `m_start=0.30, n_picard_seed=12`; that ran
      > clean at RSS 3.9 GB. **⇒ the numbers below now rest on reproducible
      > evidence, not prose.** (The sequence is still NOT asymptotic — increments
      > +2.27%, +8.05% — so no Richardson is extrapolated; that is G13.3's
      > flat-tip-cap finding, addressed on the geometry side by Track M / M5.
      > ★ Wording arbitration 2026-07-14: the verdicts these numbers support
      > read "STRONGLY INDICATED, NOT EARNED" — see below.)
      Coarse/medium likewise 0-limited, and the medium `M_max` drops
      **2.13 → 1.725** (P5's "bounded tip-TE-corner overshoot" was the same
      object). **★ AND THE LIFT GAP CLOSES WITH RESOLUTION:** cl_KJ
      **0.2593 → 0.2652 → 0.2866** (coarse/medium/fine), against the
      Tranair/KRATOS reference **0.288** — and since the taper *removes* ~1.5 %,
      the untapered equivalent is ≈ **0.291**. This lands on the
      **resolution-dominated** side of P9's PRE-REGISTERED band
      (cl_KJ∞ ≥ 0.283 ⇒ resolution-dominated; ≤ 0.278 ⇒ geometry floor ⇒ P11's
      lift case). ⇒ ★ **WORDING DOWNGRADED 2026-07-14 (user-arbitrated): the
      resolution-dominated reading of the 0.019 gap is STRONGLY INDICATED — NOT
      EARNED.** The 0.2866 point is a genuine fine discrete solution at 99.5 %
      of the reference, but the pre-registered band fires on the **Richardson
      extrapolant**, not a single fine value, and no admissible extrapolant
      exists (this sequence is non-asymptotic; the round ladder has no M0.84
      fine state — G13.3 transonic NEGATIVE). **P11's lift justification is
      correspondingly strongly indicated moot, not refuted** (P9/G9.3
      discipline). **Caveats (do not over-read):** the
      sequence is **NOT asymptotic** (increments +2.3 %, +8.1 % — they GROW) so
      **no Richardson was run**, and the mesh family is **not self-similar**
      (below), so `coarse` cannot legally enter one. The *direction* is
      unambiguous; the *extrapolated number* is not yet earned.
      **★ MESH-FAMILY DEFECT FOUND EN ROUTE (invalidates every past 3-point M6
      Richardson attempt, independent of any physics):** `generate_onera_m6.py`
      derives all sizes from `h_wall` but CLAMPS the far field —
      `h_far = min(2.5, 120·h_wall)` — so at coarse the clamp bites
      (2.5 instead of the self-similar 3.6). coarse→medium therefore refines the
      far field by only **1.39×** while the wall refines by **2×**; medium→fine is
      a clean 2× throughout. **`coarse` is not on the same refinement ray**, which
      is why the coarse→medium spanwise loading "slosh" (−6.3 % root vs +5.4 %
      mid-span, with a near-cancelling total) never looked like convergence: it
      was comparing two different mesh families. **G9.1's bands could never have
      fired legitimately even had its fine mesh converged.** Fix (cheap, not yet
      done): drop the clamp or re-ladder the family so the ratio is uniform.
      **Cost (the model's price):** `cl_KJ` **−1.1 … −1.6 %** (scales with r_c:
      −0.70 / −1.58 / −3.27 % at r_c = 0.03 / 0.05 / 0.08 b). It is **LOCAL**:
      Γ(η), cl(η) and sectional Cp are UNCHANGED inboard of η≈0.95 (inboard
      circulation −0.51 %), and TE pressure closure at η = 0.90 stays at the
      baseline value (0.232 vs 0.218). The taper also makes cl MORE
      mesh-convergent than the untapered baseline (+0.2 % vs +0.7 % coarse→medium
      — the untapered case is still gaining spurious tip lift).
      **★ The originally-proposed `tanh` form (F(tip)=½) is DISQUALIFIED — but not
      for the reason first argued.** It DOES regularize (q ≈ 1.00 — exactly the
      marginal case its s=0 predicts). It is rejected because it has **unbounded
      support**: tanh never reaches 1, so it depresses F over **57 of 83 stations**
      (inboard to η = 0.77), costing **−7.4 % lift**, **−4.9 % inboard
      circulation**, and **breaking TE pressure closure at η = 0.90** (gap 0.972 vs
      baseline 0.218) — where there is no singularity to fix. A tip model must be
      *local*; the tanh silently re-rigs the whole wing.
      **★ METRIC TRAP (cost a wrong conclusion once):** G13.1's tip box
      (z/b>0.95, dx ≥ −0.05) admits WING cells. Once the edge is regularized its
      max MIGRATES onto the ordinary wing suction peak at z/b≈0.95 and stops
      measuring the singularity — making a *working* fix look like it made p
      *worse*. The edge must be measured OFF-BODY (dx > 0, z/b > 0.98).
      Evidence: `cases/demo/p13_tip_edge_singularity/run_taper_probe.py` (the
      falsification ladder + r_c sweep) and `run_taper_physics.py` (Γ/cl/Cp
      distributions); tests `tests/test_p13_tip_taper.py` (15).
      **OPEN clause — level-set port.** NOT a mechanical port: the LS path has
      **no Γ DOF** and its TE row (`kernels/cut_assembly.py::te_kutta_coo`) is
      **homogeneous** (`s·(q_u−q_l) = 0`), so scaling it by F is a NO-OP. The
      clean analogue is to BLEND the pressure-equality row with B2's continuity
      weld (`[φ]=0`), normalized: `F·K̂ + (1−F)·Ŵ = 0` (pure Kutta at F=1, zero
      jump at F=0). That is a *different model* from `Γ = F·Γ_Kutta`, so it needs
      its own r_c calibration and the two-path A/B is a PHYSICS comparison, not a
      model-identity check. Designed, not implemented.
- [ ] **G13.3 — 3D grid-convergence closure (PARTIALLY ANSWERED; formal
      Richardson still BLOCKED).** ★ **The substantive question P9 posed is now
      ANSWERED, on the resolution side** (see G13.2's transonic entry): with the
      tapered fine mesh a genuine discrete solution at M0.84, cl_KJ reaches
      **0.2866** (untapered-equivalent ≈ 0.291) against the Tranair/KRATOS
      **0.288** ⇒ the **resolution-dominated reading is STRONGLY INDICATED —
      NOT EARNED as a P9-band verdict** (★ wording downgraded 2026-07-14,
      user-arbitrated: the band fires on the Richardson extrapolant, and none
      exists on either geometry — see the transonic entry below). **P11's lift
      case: strongly indicated moot, not refuted.**
      ★ **But the FORMAL Richardson remains blocked, for TWO reasons, and the tip
      is neither of them:** (a) the sequence is **not asymptotic** — increments
      **GROW** (M0.84: +2.3 %, +8.1 %; M0.5: −0.5 %, +5.8 %) ⇒ **no extrapolation
      was run** (P9/G9.3 discipline: report `n/a`, never fabricate); (b) the M6
      **mesh family was not self-similar** — `h_far = min(2.5, 120·h_wall)` clamps
      at coarse, so coarse→medium refined the far field 1.39× against the wall's
      2× ⇒ **`coarse` was not on the refinement ray and could never legally enter
      a 3-point Richardson** (this also explains the coarse→medium spanwise
      "slosh": −6.3 % root vs +5.4 % mid-span, near-cancelling in total — it was
      comparing two different mesh families, not measuring convergence).
      **Removing the tip singularity REVEALED both of these rather than curing
      them.** **(b) is now FIXED** (mesh ladder, below); **(a) survives the fix,
      and the reason is a THIRD singularity.**
      ★ **The growth is REAL FIELD RESOLUTION — NOT a Kutta-probe artifact
      (hypothesis raised and REFUTED 2026-07-13).** Decisive test: compare the
      probe-based `cl_KJ` against the **probe-free, pressure-integrated `cl_p`**
      across the three meshes. **Both rise, both with growing increments** —
      cl_KJ 0.2001 → 0.2005 → 0.2121 (+0.19 %, +5.76 %) vs
      cl_p 0.1875 → 0.1939 → 0.2091 (+3.41 %, +7.87 %) — so the conforming
      Kutta-probe extraction is EXONERATED (its geometry is healthy too:
      `probe_distance/h` is constant at 0.48–0.49 across all three levels, i.e.
      the probes scale correctly with h). The *solution field itself* is simply
      under-resolved. ⇒ **The level-set port is NOT the route to G13.3** (it
      remains a parity/robustness item, not an accuracy one), and the earlier
      suspicion that B7's smoother Γ(z) might be load-bearing for 3D accuracy is
      **not supported**.
      ★ **BONUS — the P5-era V6 item is behaving exactly as re-specced.** The
      V6 gap (cl_p below cl_KJ) converges cleanly at ~O(h):
      **6.30 % → 3.29 % → 1.41 %** (ratios 1.9, 2.3) — a genuine O(h)
      discretization floor that REFINES AWAY, now at **1.41 % on fine**, closing
      on the < 1 % target deferred since P5/G5.2. The discretization is
      consistent; the meshes were simply too coarse.
      **★ MESH LADDER FIXED 2026-07-13** (`generate_onera_m6.py`): `_level_params`
      grew `clamp_h_far` (default True ⇒ `coarse`/`medium`/`fine` stay
      **BIT-IDENTICAL**, so the P5 / P8-G8.2 / B7 / M1 locks are untouched), plus
      a new **`coarse_ss`** level — the ONLY level the clamp ever touched — recut
      without it (h_far 2.5 → 3.6, 46,067 tets, quality within M1 bounds:
      min dihedral 4.15°, max aspect 15.5). `RICHARDSON_LADDER =
      ("coarse_ss", "medium", "fine")` refines by **exactly 2.000 in EVERY length
      scale** (h_wall/h_edge/h_wake/h_far), locked by 3 new tests in
      `tests/test_m1_onera_m6.py`.
      **★★ BUT THE SEQUENCE IS STILL NOT ASYMPTOTIC — AND THE REASON IS A THIRD
      SINGULARITY, ON THE WALL.** On the self-similar ladder the increments still
      GROW (M0.5: cl_KJ 0.2015 → 0.2005 → 0.2121, cl_p 0.1888 → 0.1939 → 0.2091 —
      the cl_p increments grow ~3× while h halves). Growing increments under
      *uniform* refinement are the signature of a singularity still being
      resolved. Localized by a box study on the tapered solutions, now a
      committed self-checking demo —
      **`cases/demo/p13_tip_edge_singularity/run_g133_ladder.py` (5/5 PASS)**,
      artifacts `results/g133_boxstudy.csv` + `g133_ladder.csv` +
      `g133_ladder.png`, run on the SELF-SIMILAR ladder:

      | region | coarse_ss → medium → fine | p | verdict |
      |---|---|---|---|
      | **tip-cap edge (WALL)** | 0.662 → 0.824 → **1.015** | **+0.321** | **DIVERGES** |
      | wing p99 (control) | 0.642 → 0.628 → 0.629 | −0.014 | bounded |
      | wake free edge (G13.2-fixed) | 0.536 → 0.565 → 0.570 | +0.045 | bounded |

      *(Numbers corrected by the 2026-07-13 audit. The first-written table mixed
      the `coarse` and `coarse_ss` ladders and mis-stated the wake-edge row as
      0.567 → 0.579 → 0.649 / p = +0.107; re-measured from the caches it is
      0.536 → 0.565 → 0.570 / p = +0.045. The tip-cap fine peak is 1.015, not
      0.967 — a box-definition difference. **All three VERDICTS are unchanged**,
      and the demo now pins them.)*

      ⇒ The wake edge we fixed stays bounded and the wing interior is converged,
      but the **flat tip cap's sharp wall edge diverges**. **Root cause: a
      DOCUMENTED DELIBERATE GEOMETRY SIMPLIFICATION** — `meshgen/wing3d.py`
      builds a **flat tip cap** where the real ONERA M6 has a **rounded** one
      ("the real rounded tip cap is a deliberate simplification — standard for FP
      validation meshes"). A flat cap meets the upper/lower surfaces at a sharp
      convex edge, which in potential flow is an edge singularity of exactly the
      kind P13 exists to remove — only on the BODY rather than the wake. **Note
      this is NOT a P11 (curved wall ELEMENT) problem:** isoparametric elements
      cannot regularize a genuinely sharp geometric edge; the GEOMETRY is what is
      wrong.
      **⇒ NEXT for G13.3: round the tip cap** (Track M, `wing3d.py`), regenerate
      the ladder, then re-run the three-point Richardson. Neither P11 nor the
      level-set port is required.
**Non-goals:** fine meshes stay gitignored; heavy runs are one-shot cached
artifacts (the conforming fine M0.5 AMG solve is ~10–34 min, the M0.84 fine
continuation ~1 h).

---

### P14 — Probe-free conforming Kutta target: wall-adjacent-CV pressure-equality estimator ✓ CLOSED 2026-07-17 (user-directed open + close same day; G14.1–G14.7 all ✓; the conforming path now matches the level-set path on lift, Γ(z) smoothness, and TE Cp closure)

> **Origin.** Track A / A2 (`cases/analysis/a2_te_kutta_fidelity/`, closed
> 2026-07-17; [roadmap/track_a.md](track_a.md) A2, dossier
> [demo_report/track_a.md](../demo_report/track_a.md#track-a--a2--tekutta-fidelity-z-jitter--te-cp-jump-casesanalysisa2_te_kutta_fidelity-2026-071617))
> settled the cause of two long-standing conforming symptoms and routed the
> fix here. **A2 measured; A2 implemented nothing.** This phase is the
> implementation, and it is a change to the conforming **Kutta enforcement
> mechanism**, which is Track P's remit (Track A adds no physics).

**The defect A2 pinned (two symptoms, one estimator).** The conforming Kutta
target is a per-station potential-jump read at probe nodes one edge off the TE
(`constraints/wake.py::kutta_targets`, `mesh/wake_cut.py::_kutta_probe_nodes`).
Two consequences, both measured (A2 GA2.2/GA2.4):
- **S1 — Γ(z) jitter.** The probe estimator is a single-node sampler of a
  field that varies sharply at the swept, unstructured TE (probes off-plane by
  O(h); 35+41 of 166 stations share a probe). The fixed-Γ discriminator
  D = 7.33 (coarse) / 25.70 (medium) proved the estimator *manufactures* the
  jitter from a smooth field — it is a measurement-operator artifact, not flow
  content (closure |F|/max|Γ| ≤ 0.6%, so not unclosed stations either).
- **S2 — TE Cp jump.** Enforcing equal *potential* jump is only an
  approximation of equal *pressure* (|q_u|²=|q_l|²), so the upper/lower Cp land
  apart at the TE: same-estimator gap 0.318/0.228 (conforming) vs 0.009/0.002
  (level-set), 34×/133×; P13's prose 0.14–0.22 confirmed as the sweep median.

**Design.** Replace the probe-difference potential-jump target with a
**wall-adjacent control-volume recovered-velocity / pressure-equality**
estimator — the B4 objects (`wake/multivalued.py::_build_te_control_volumes` +
`te_velocities`, which recover each surface's velocity over a *consistent*
TE-adjacent element set rather than a single off-TE node) ported to the
conforming cut mesh. This kills the single-node sampling degeneracy (S1) and,
being equal-pressure not equal-potential, closes the TE Cp gap (S2) — one
estimator swap, both symptoms. Because pressure equality is nonlinear in Γ, the
per-station secant becomes a per-station nonlinear closure, which the transonic
path already approximates (`kutta_per_outer=1`); the subsonic path's affine
secant reasoning is replaced by the recovered-velocity residual.

**Gates (firmed at open 2026-07-17 — user-arbitrated TWO-TIER structure).**
Mapping from the provisional list: provisional G14.1/G14.2/G14.3 become tier-2
G14.5/G14.6/G14.7 with their numbers unchanged; provisional G14.4 keeps its ID;
tier 1 is the NEW subsonic milestone — a deliverable stop point if the M0.84
nonlinear per-station closure proves to need its own damping study (the spec's
recorded risk). **Wiring scope (user-arbitrated same day): the pressure
estimator wires into the coupled Newton drivers (`solve_newton_lifting` /
`solve_newton_transonic`) + `solve_laplace_lifting` (cheap B4-parity
verification) ONLY; `solve_subsonic_lifting`'s inner Kutta secant and
`continuation.py`'s Γ-secant/polish stay probe-based (the P5-era warm-start
engines — excluded to keep the blast radius down; the estimator MEASUREMENT
function is solver-independent, which is all the discriminator rerun needs.**

Tier 1 — subsonic M0.5 milestone:
- **G14.1 ✓ (2026-07-17) — estimator validity.** CV construction asserts all
  pass (two-sided non-empty wall-adjacent fans via EXACT wall-face ownership,
  probe-membership side identity, zero far-field-Dirichlet contact) + dF/dΓ
  non-degenerate and FD-exact at committed states + implied Γ* within the B4
  band: NACA incompressible pressure-vs-probe Γ **0.074%** (demo) / 0.01%
  (Stage-D at the converged probe state). Evidence:
  `cases/analysis/p14_te_pressure_diag/results/` +
  `cases/demo/p14_pressure_kutta/results/checks.csv` +
  `tests/test_p14_te_pressure.py` (15). One firming correction: the
  "uniform-sign dF/dΓ diagonal" clause holds at CONVERGED states (Stage-D,
  0 flips at all five) but NOT at a rough Picard-5 seed (measured: flipped
  stations on M6 medium M0.5) — σ-freeze demotes it to a recorded count
  (`kutta_sigma_sign_flips`; σ is merit weighting only, the per-step exact D
  carries the true signs).
- **G14.2 ✓ (2026-07-17) — S1/S2 A/B at M0.5.** M6 M0.5 Newton same-run A/B,
  committed `m05_ab.csv` + `m05_ab.png`: Γ(z) roughness probe→pressure
  **0.1203 → 0.0052** (coarse) / **0.0504 → 0.0024** (medium — BELOW the LS
  band 0.003–0.009); all-station raw TE Cp gap median **0.2278 → 0.0045** /
  **0.1603 → 0.0026** (51×/62×). Both symptoms gone in one swap, at M0.5,
  before any transonic tuning.
- **G14.3 ✓ (2026-07-17) — subsonic lift: cross-read + recorded move.**
  **Re-specced at the first tier-1 run — the provisional "< 1%" band was
  wrong in principle, twice** (first firmed < 1%, then < 3%, both failed
  honestly): the two closures agree POINTWISE to the probe's own O(h) reading
  bias (cross-read at the pressure-converged state: **3.67% coarse → 1.05%
  medium**, O(h) decay confirmed), but the Kutta map's near-unity slope
  (b ≈ 0.93, P2 record) amplifies estimator bias by 1/(1−b) ≈ 14× into the
  converged circulation — so the converged lift MUST move by the amplified
  probe-bias correction. Measured move: cl_p/cl_KJ **+2.1% coarse, +4.5%
  medium** at M0.5 (grows under refinement because the pressure closure
  converges while the probe path's amplified O(h) bias shrinks more slowly).
  Final gate form: (a) cross-read < 5%/2% (coarse/medium) + O(h) decay —
  catches a self-consistent-but-garbage closure (the B8 failure mode);
  (b) |Δcl| < 8% with the spanwise dΓ committed (`dgamma_*_m05.csv`).
- **G14.4 ✓ (2026-07-17) — inert by default.** New estimator behind
  `kutta_estimator` (default "probe"). Inertness is **structural, not just
  measured**: `kutta_blocks` returns the literal pre-P14 objects `(ws.K, F)`
  on the probe path, so the driver's matvec/rhs/Woodbury/step-map expressions
  execute unchanged; the CV object is built only when "pressure" is requested
  (the probe `__init__` does no extra work); the Laplace target swap is an
  `if` around the same `kutta_targets` call. Full suite **421 passed + 18
  skipped + 2 xfailed** (16m55s @8 threads) = the 406 baseline + P14's 15 new
  tests, **zero regressions** — every committed P4/P5/P8/P13 conforming lock
  green.

**Tier-1 solver-recipe finding (recorded 2026-07-17): seed the pressure
Newton from the probe solution.** The quadratic pressure row has a smaller
Newton basin than the affine probe row: on M6 medium M0.5 a Picard-5 cold
seed wanders to cl +16% and fail-fasts at step 29, while the probe-seeded
solve converges in **3 quadratic steps** (|F| 6e-3 → 2.8e-4 → 2e-6 → 5e-11,
26 s). Demo recipe: A/B legs seed pressure from the probe solve; the M0.84
ramp seeds its level 0 (M0.70) from a probe Newton solve at the same Mach
(later levels warm-start from the previous pressure level as usual). This is
the spec's "nonlinear closure may need its own damping" risk landing in its
mild form — a seeding rule, not a damping scheme.

**Pre-registered interpretation note for G14.7 (written BEFORE the tier-2
runs; FIRED — see G14.7 above).** G14.7's band ("< 1–2% vs the G8.2 locks")
was written on the assumption the estimator swap should not move loading.
Tier 1 measured that assumption wrong at M0.5 (+4.5% medium, mechanism above).
The G8.2 locks are PROBE-path locks; if the M0.84 pressure lift moves upward
by a similar amplified correction, cl_KJ moves TOWARD the Tranair/KRATOS 0.288
reference — i.e. into the P9 "0.019 gap" that P9 closed as
unsplittable-as-posed. A G14.7 FAIL-as-written in that direction is therefore
a candidate ACCURACY finding, not a defect; the gate is reported against its
pre-registered band and the verdict is user-arbitrated either way.
*(Outcome: it fired in exactly that direction — +4.85% cl_KJ, 69% of the gap
closed. The note is kept verbatim as the audit trail that the band was not
retrofitted to the result.)*

Tier 2 — transonic M0.84 (the A2 regime; provisional numbers unchanged):
- **G14.5 ✓ (2026-07-17) — S1 removed.** Pressure-estimator M0.84 Newton ramp
  (probe-seeded level 0), both committed M6 levels converged with 0
  limited/floored (coarse 11 steps |R| 7.0e-12; medium 12 steps |R| 5.6e-15).
  Γ(z) roughness **0.0043** (coarse, from A2's probe 0.0970 = 23×) and
  **0.0024** (medium, from 0.0365 = 15×) — both AT/BELOW the level-set band
  0.003–0.009. A2's fixed-Γ discriminator rerun on the NEW estimator:
  **D = 1.80** vs the probe's 7.33 (pre-registered confirm > 3 / refute < 1.5).
  **Honest reading: D = 1.80 lands in A2's INCONCLUSIVE zone, not at 1.0** —
  the pressure estimator still regenerates ~1.8× the roughness of a smooth
  input, i.e. it is 4× cleaner than the probe but not a perfect measurement
  operator. The gate's "→ O(1)" is met in the sense that mattered (the
  jitter-manufacturing mechanism is gone: 0.0043 absolute is LS-grade), and
  the residual 1.8 is recorded, not claimed away.
- **G14.6 ✓ (2026-07-17) — S2 removed.** All-station raw TE Cp gap median
  **0.0040** (coarse) and **0.0024** (medium) — both < 0.02 on the PRIMARY
  clause, raw recovery; the pre-registered fallback (≤ 3× LS band,
  smooth_passes=1) was not needed. **★ Baseline erratum (caught 2026-07-17,
  same day, on a user question — corrected in the demo + every doc):** the
  first write-up quoted these against A2's **0.318/0.228** and claimed
  **80×/95×**. Those are A2's *section-last-point* numbers; this demo measures
  the *all-station sweep*, whose A2 probe baseline is **0.2206/0.1585**
  (`a2_te_gap.csv`). The correct factors are **55× / 67×** — still the effect,
  but the metric had to match. Same trap on the other side: `a2_te_gap.csv`'s
  LS rows read ≈ 0 because they evaluate the LS's OWN control volumes (its
  own constraint residual, the "cannot be used as A/B" caveat), so V14.6 now
  re-measures the LS wall through this demo's sweep instead: **LS = 0.0047**
  (medium), i.e. the conforming pressure path is ~2× BELOW the level-set path
  on its own metric.
- **★ Correction to P14's own claim, and a partial correction to A2's S2
  decomposition (measured 2026-07-17 on a user question; V14.7,
  `te_spike_medium_m084.csv`).** Every earlier P14 write-up asserted that the
  TE Cp **spike** (A2's `spike_metric`: the last section point's deviation from
  its own x/c ∈ [0.85,0.97] trend — a COMMON-MODE metric, unlike the
  differential gap) would be **untouched**, because A2 read it as "a P1
  recovery artifact, present on the level-set path too ⇒ wake-model
  independent". That was a *prediction from A2's reasoning, never measured on
  the pressure path*. Measured now, medium M0.84 raw (smooth_passes=0, mean
  over A2's 4 η × 2 sides): conforming probe **0.1143** → conforming pressure
  **0.0533** (2.1×), which is also **below the level-set path's 0.0743**.
  **What A2 got right:** a spike remains (~0.05), and it is shared — the LS
  path has one too. **What needs correcting:** the conforming path's EXCESS
  over LS was *also* Kutta-form error, not recovery. That is mechanically
  sensible: a wrong Kutta gives a genuinely wrong TE flow, and the last point
  then departs from the upstream trend for a physical reason, which a
  common-mode metric cannot separate from a recovery artifact. Corroborating
  detail: P6 normal-gated smoothing no longer helps on the pressure path
  (0.0533 → 0.0660 → 0.0626 over 0/1/2 passes; A2 measured 0.147 → 0.081 on
  the probe path) — consistent with the Kutta-form component being gone and
  only the genuine recovery floor left for smoothing to chew on.
  **Still open:** that residual ~0.05 floor (the honest remainder of A2's
  second S2 component) — P14 does not close it.
- **G14.7 ✓ CLOSED — RE-SPECCED at close (user-arbitrated 2026-07-17): lift
  matches the level-set oracle.** The gate opened against the G8.2 **probe**
  locks on the premise that a Kutta *estimator* swap should not move loading.
  Tiers 1–2 measured that premise wrong for a measured reason (below), the
  gate XFAILed as written (band NOT moved after the fact), and the user
  arbitrated the recorded verdict: **accept the move as the finding, re-lock
  against the level-set independent oracle.** Re-specced acceptance —
  conforming pressure vs level-set, medium M0.84: cl_p **0.2776 vs 0.2772
  (0.15%)**, cl_KJ **0.2823 vs 0.2813 (0.34%)**, both **< 1%** (demo
  `medium_m084_lift_vs_levelset` PASS; V14.6 confirms live). **Recorded
  context — the move off the OLD probe locks:** cl_p **+4.92%**, cl_KJ
  **+4.85%**. The mechanism was pre-registered before the run (note below): the
  two closures agree pointwise to the probe's own O(h) reading bias (cross-read
  at the M0.84 pressure state 2.52% coarse / **0.79% medium** — a shifted
  closure, not a wandered solution), which the Kutta map's b ≈ 0.93 amplifies
  1/(1−b) ≈ 14× into Γ, so the converged lift MUST move. Direction: |cl_KJ −
  0.288| goes **0.0188 → 0.0057, 69% of P9's "0.019 gap" closed** by an
  estimator swap. **What closing this gate does NOT claim:** not a
  grid-convergence result (P9: the M6 fine mesh is not a discrete solution, so
  no Richardson exists here), not a re-opening of "the 0.019 gap is
  resolution" (still *strongly indicated, NOT earned* — 2026-07-14 wording
  arbitration), and not proof the pressure lift is *absolutely* right — what is
  established is that the conforming and level-set paths now AGREE, i.e. a
  measurable share of the old gap was **Kutta-estimator bias** that P9 could
  not see because both its meshes used the same estimator.
  **★ Independent corroboration (V14.6, added 2026-07-17 on a user question;
  `cross_model_medium_m084.csv`).** The level-set path has ALWAYS used
  pressure-equality Kutta on wall-adjacent CVs (B4). If the conforming lift
  move is really the Kutta *form*, the pressure path must land on the LS
  answer — and it does, from a different wake model (multivalued aux DOFs, no
  Γ DOF), a different DOF space, and a *different mesh family*
  (`onera_m6_wakefree`):

  | medium M0.84 | cl_p | cl_KJ | roughness | all-station TE gap |
  |---|---|---|---|---|
  | conforming **probe** (G8.2 lock) | 0.2646 | 0.2692 | 0.0365 | 0.1585 |
  | conforming **pressure** (P14) | **0.2776** | **0.2823** | **0.0024** | **0.0024** |
  | **level-set** Newton (B15/A1 cache) | 0.2772 | 0.2813 | 0.0033 | 0.0047 |

  The two independently-implemented wake models now agree to **0.17% / 0.36%**
  on cl_p / cl_KJ, where the probe path sat **4.5% / 4.3% below** the LS one —
  i.e. the long-standing conforming-vs-LS lift disagreement WAS the Kutta
  form, and it is gone. Caveats kept: different mesh families (not a
  same-mesh A/B), and the LS state carries 1 limited / 2 floored cells (the
  B15 `freeze_max_clamped` caveat) while the pressure state has 0/0.
  **This cross-model agreement is the physical-acceptance basis on which the
  user closed G14.7** (2026-07-17): the modified conforming result matches an
  independent implementation that already used the correct Kutta form, so the
  computed result is judged reasonable.

**Dead routes (A2, do not re-propose):** spanwise-Γ smoothing (moves Γ off the
self-consistent value — P5 `INVESTIGATION_gamma_smoothing.md`); full-element-fan
recovery at the TE (B4: +11–15% wrong, wall-adjacent < 1%); better probe-picking
alone (a band-aid on the sampling, not a fix of the mechanism).

**Risks / unknowns (as designed):** porting B4's control volumes onto the
conforming cut mesh has to cope with the duplicated TE nodes and the explicit
per-station Γ DOF (the level-set path has neither); the tip region (P13) and
the wall-adjacent-CV construction there are unexamined; a nonlinear per-station
closure may need its own damping. Diagnostic-first at open: build the
conforming TE control volumes and verify the recovered two-sided velocity is
non-degenerate in Γ before wiring the residual.

**Stage-D diagnostic — GO (2026-07-17, `cases/analysis/p14_te_pressure_diag/`,
20/20 checks, committed `results/diag_checks.csv` + `diag_states.csv`).**
Retired at the diagnostic: (a) CV construction is clean on NACA coarse + M6
coarse/medium — two-sided wall-adjacent fans never empty (min fan 1 element on
the NACA family, 2 on M6; recorded), exact wall-face ownership, probe-membership
side identity (slave="+"=upper ties to the `upper_hint` convention
structurally), zero far-field-Dirichlet contact (⇒ dF/dΓ carries only the
slave-jump chain, no V_red term); (b) **non-degeneracy is measured, not
assumed**: dF/dΓ is tridiagonal (bandwidth 1 in station index), uniform-sign,
|D_jj| ∈ [29, 216] across states with cond(D) 4.2–6.4 (NOT always strictly
diagonally dominant — freestream/coarse margins go to −0.27, which is why the
implied-target solve is dense, dominance being only sufficient); central-FD of
D at fixed φ_red is exact to roundoff (7.1e-11 / 1.2e-10 — F is exactly
quadratic in Γ, validating the whole row-builder chain); (c) **estimator
plausibility on the committed converged probe states**: implied Γ* deviates
from the closed Γ by 0.01% (NACA incompressible — the B4 band) and
0.88–2.15% median (M6 M0.84 states — the expected S2-scale offset); (d) **S1
preview**: measuring Γ*(z) on the SAME cached fields gives roughness
0.0226/0.0081/0.0074 (P5-coarse/P5-medium/Newton-medium) vs the probe-target
curve's 0.0965/0.0389/0.0365 — the pressure estimator does not manufacture the
jitter, and the medium numbers already sit in the LS band. Remaining risks →
tier 2: the nonlinear closure's behaviour inside an actual transonic Newton
ramp (damping), and the tip station under a live pressure closure.

**Trigger (designed-not-started, like B14):** open when a conforming-path
accuracy campaign needs smooth spanwise loading or a closed TE (e.g. a Track V
viscous coupling that reads sectional Cp, or an MDO objective on Γ(z)), or on
user direction. Not on the critical path to B9 (level-set, which has neither
symptom by construction).

---


## Progress ledger

### Track P — solver

| Phase | Status | Closed on | Notes |
|-------|--------|-----------|-------|
| P0 | ✓ | 2026-07-06 | `mesh/reader.py`, `metrics.py`, `coloring.py`, `physics/isentropic.py`, `post/vtk_out.py` implemented; G0.1–G0.4 unit tests pass. Three latent bugs found by manual audit and fixed, each now locked in by a regression test (`test_mesh_adjacency.py`, `test_mesh_reader_roundtrip.py`, `test_laplace_picard.py`): `metrics.py::build_face_adjacency` crashed under `@njit` (reflected-list dict values), `reader.py::write_mesh` dropped all named boundary tags (`.msh` writer ambiguity + only handled a legacy `"all_triangles"` block), `solve/picard.py::solve_laplace` reported a `residual_norm` dominated by Dirichlet-row flux imbalance instead of the free-dof residual. None were caught by the pre-existing test suite because nothing exercised those paths. Second audit (2026-07-06), same fix-plus-regression-test pattern: `tests/conftest.py` now writes gate artifacts to the persistent `artifacts/` dir instead of a deleted tempdir (`test_conftest_artifacts.py`); `metrics.py::element_gradients` raises on degenerate tets with a scale-relative threshold instead of silently returning zero gradients (`test_metrics_degenerate.py`); `reader.py` keeps the default "bulk" volume tag for meshes without named 3D groups (`test_mesh_reader_roundtrip.py`) — details in PROJECT_STRUCTURE.md. Closed 2026-07-06: the two remaining blockers were resolved by the M0 delivery — the mesh family exists on real geometry (cylinder + NACA0012) and the full coarse regression suite runs against those case meshes (87 passed, 2 xfailed); G0.1–G0.4 all green (`test_mesh_volume.py`, `test_mesh_gradient.py`, `test_mesh_coloring.py`, `test_io_vtk.py`). |
| P1 | ☐ (in progress; G1.1, G1.2 closed; G1.3, G1.4 completed 2026-07-06 with negative results, DP1 decided "> 5%" branch, G1.5 void; open: G1.6 pending its Option C re-spec) | | `kernels/residual.py`, `solve/linear.py`, `solve/picard.py`, `post/surface.py`, `tests/mesh_utils.py`, `tests/test_post_surface.py`, `cases/meshes/sphere_shell/{coarse,medium}.msh` implemented/committed. G1.1 (`test_laplace_mms.py`) and G1.3 (now renumbered G1.2; `test_laplace_cg_iterations.py`) pass. G1.2 (now renumbered G1.6; `test_laplace_sphere.py`) is a `strict=True` xfail, now root-caused rather than just hypothesized: added `post/surface.py::wall_tangential_gradient_quadratic` (quadratic tangential patch recovery, well-conditioned via SVD-based `lstsq` + rank-deficiency fallback), which improved medium-mesh max error from ~12.0% to ~11.6% — a real but modest gain, because a controlled investigation (clean single-variable h_min sweep + an oracle exact-potential-in recovery test + a rejected Nitsche/penalty prototype) showed the recovery scheme was never the dominant error source. The dominant source is the volume PDE solve's own accuracy at the wall: the natural BC is satisfied on the flat polyhedral wall-facet approximation instead of the true curved sphere, a geometric/variational-crime inconsistency evidenced by sub-first-order, decreasing convergence order of the raw nodal potential itself (not just its gradient). Closing G1.2 needs curved/isoparametric wall boundary elements — a separately-scoped effort — not more h-refinement or post-processing. See PROJECT_STRUCTURE.md "Known gaps" for the full evidence trail. Second audit (2026-07-06): `post/surface.py::_wall_vertex_normals` now raises on inconsistent wall-triangle winding instead of silently averaging cancelling normals into garbage tangent planes (`test_post_surface.py::test_inconsistent_wall_winding_raises`). 52 tests total (51 passed, 1 xfailed), full suite runs in ~10 s. 2026-07-06: sphere-Cp fix-route research complete — three-tier Option A/B/C plan (lagged true-normal flux correction / Gap-SBM / gate redefinition) with oracle-first verification order defined; see design.md §5.1. Same day: P1 gates renumbered in workflow order (old → new: G1.3 → G1.2, G1.2-a0 → G1.3, G1.2-a → G1.4, G1.2-c → DP1, G1.2-b → G1.5, G1.2 → G1.6; see the renumbering note in the P1 section). Same day, after renumbering: G1.3 cylinder oracle pre-study and G1.4 sphere oracle ceiling both completed with **negative results** (delivered: `solve/wall_correction.py` assembly-verified correction infrastructure, `post/section_cut.py` degenerate single-layer interface, `tests/test_wall_correction_cylinder.py` 10 tests, the sphere-oracle demo — since 2026-07-07 absorbed into `cases/demo/p1_laplace/run_demo.py` — cylinder fine.msh 50.2k tets, `artifacts/G1.3/` + oracle results in `cases/demo/p1_laplace/results/`): boundary-data corrections have (near-)zero lever on body-fitted meshes — exact per-facet net flux is zero by the divergence theorem — so the Option A ceiling is ≈ 11.3% on the medium sphere vs the < 2% target; DP1 "> 5%" branch taken (Option C gate re-spec + separately-scoped curved elements); the cylinder case additionally shown to be recovery-dominated (~76%), not crime-dominated, and de-designated as the G1.6 pathology testbed. See the G1.3/G1.4/DP1 gate entries for the full evidence. |
| P2 | ✓ | 2026-07-06 | Delivered: `mesh/wake_cut.py` (per-node flood-fill side classification — no planarity assumption, works for the M1 swept wake; ⁺-side slaves appended after original nodes so the reduced dof space is exactly the original node set; stations group TE nodes by (x,y), collapsing a quasi-2D extrusion to the single scalar Γ of the M0 spec; preprocess-time topology asserts), `constraints/wake.py` (master–slave elimination via sparse T, A_red = TᵀAT assembled once, Γ enters RHS-only through precomputed per-station vectors; folding the slave rows into masters enforces weak flux continuity (4.2)), `constraints/dirichlet.py` (far-field freestream + incompressible 2D vortex correction with the branch cut ON the wake sheet, so eliminated ⁺-side far-field wake nodes are automatically consistent), `solve/picard.py::solve_laplace_lifting` (Kutta outer loop, matrix+AMG hierarchy built once, secant-accelerated), `post/surface.py` (triangle-wise wall force integration, owner-tet-oriented normals, KJ sectional cl), `post/section_cut.py` (general marching-tets z=const path + sectional wall Cp(x/c) curves), Hess–Smith panel reference `cases/reference_data/naca0012_incompressible/` (two independent lift routes agree to 0.09%, lift slope 6.91/rad vs thickness-corrected 6.90). **One spec deviation with evidence: TE nodes ARE duplicated** — the originally specified single-valued TE produces a spurious TE suction ~Γ²/h that diverges under refinement (measured −0.27 of cl 0.6 on coarse; see the re-specced topology-assert block). Gates: G2.1 8.4e-13 (wake-master rows 6.9e-16); G2.2 [φ]−Γ < 1e-13; G2.3 medium cl −0.82% vs panel, Kutta 2 updates (secant; measured map slope b≈0.93 would need O(100) plain relaxed updates); G2.4 0.01%; G2.5 re-specced criterion (b) closed (p99 ratio 2.05 = 1st order, stripe-free heatmap). Suite: 100 passed + 2 xfailed (G1.6 strict xfail unchanged), ~38 s. Artifacts: `artifacts/G2.{1,2,3,4,5}/` PNG+CSV. |
| P3 | ✓ | 2026-07-07 | Delivered: **assembly tech debt retired** — `mesh/metrics.py::precompute_element_geometry` (B_e/V_e once per mesh), `mesh/coloring.py` numba-jitted greedy coloring (same visit order ⇒ identical assignment to the old pure-Python loop, which was ~seconds per call on real meshes) + `color_partition_csr`, `kernels/gradient.py` (prange velocity sweep, zero-alloc), `kernels/jacobian.py` (symbolic CSR pattern + `elem_to_csr` scatter map + colored-prange matrix kernel + `PicardOperator` per-mesh workspace), `kernels/residual.py::assemble_residual_colored`; the public `assemble_stiffness_matrix` now delegates to the fast path (P1/P2 drivers run the same code as the Picard loop) with the old serial kernels retained as the regression reference — fast-vs-reference 5.7e-16 rel, bit-deterministic across calls/threads (within a color no two elements share a node, so accumulation order is fixed by the color sequence), hot reassembly ~160× faster on the medium NACA mesh (`tests/test_p3_assembly.py`, demo part 1). **Subsonic compressible solver**: `physics/isentropic.py::density_field/mach_squared_field` (array sweeps of the §2 scalars; ρ ≡ 1.0 *bitwise* at M∞ = 0 — the G3.3 anchor), `solve/picard.py::solve_subsonic` (non-lifting density Picard) and `solve_subsonic_lifting` (nested: outer density update, inner P2 secant Kutta at frozen ρ; AMG reuse every 4 outers; opt-in forcing-term inexact solves η‖b−Ax₀‖, default off — see the G3.2 gate entry for why interleaved Γ updates and relative loose tolerances were both rejected with measurements), PG-scaled vortex far field (`constraints/dirichlet.py`, β = √(1−M∞²) stretches only the atan2 argument so the wake-jump/branch-cut structure is untouched and β = 1 reduces bit-exactly), `constraints/wake.py::WakeConstraint.update_matrix` (T topological, rebuilt never; A_red + h_j per density iteration), compressible Cp in `post/surface.py::wall_force_coefficients` + `post/section_cut.py::wall_cp_curve` (`m_inf` param, isentropic (2.5)), `solve/linear.py::build_amg_preconditioner` (seeded AMG setup — repeatable solves; see G3.3). Reference data: `cases/reference_data/naca0012_m05/` (PG + Kármán–Tsien corrected panel cl/Cp with provenance + verification trail). Gates: G3.1 **0.32%** (< 2%); G3.2 **cl −0.33%** from the PG/KT midpoint and inside the bracket, **15 iterations** (< 30), strictly monotone residual; G3.3 matrix/φ/Γ **bitwise** at M∞ = 0 + full suite green. Suite: 117 passed + 2 xfailed, ~96 s (G3.2's medium-mesh nested solve is ~45 s of it). Demo: `cases/demo/p3_subsonic/` (14 checks PASS) + docs/demo_report.md §P3. Known non-P3 fix bundled: `tests/test_p2_wake_cut.py` topology sweep now skips surface-only mesh assets (the new `cessna_surface.msh` broke `read_mesh` in the hard-rule-7 sweep — pre-existing on main). |
| P4 | ✓ | 2026-07-07 | **★ ERRATUM (2026-07-11, P8/N5 finding — user-approved: record, do not re-open):** the P4 "engineering-converged" Picard states are **NOT solutions of the discrete equations**. Measured with the P8 exact Newton machinery: the coupled Newton residual at the committed coarse M0.80/α1.25 Picard state is **2.2e-4** (Kutta |F| 1.5e-4), and Newton started FROM that state leaves it in 6 quadratic steps for the true discrete solution — **shock x/c 0.658, cl_p 0.459, M_max 1.408** (coarse; residual 8e-11, 0 limited/floored, dissipation-scan robust upwind_c 2.0→1.5 and continuation-path independent to the last bit). The heavily damped pseudo-time Picard iteration stalls at a non-solution whose shock position (0.604) happens to sit in the Euler-anchored band; the documented "bounded, slowly decaying residual tail" IS this stall. This also retro-explains the G4.3 M0.82 corner's damping-dependent cl 0.389→0.458 ("different converged state" — the less-damped path got closer to the actual solution). On the MEDIUM mesh the true solution family steepens into the FP non-uniqueness fold (Newton-converged: M0.775 shock 0.570/cl 0.396 @1.8e-13; M0.7875 shock 0.674/cl 0.523 @7.9e-11; **no reachable isolated solution at M0.80**, M_max≈1.45 beyond the isentropic validity envelope — conservative FP over-lifts strong-shock cases vs Euler, Holst PAS 2000; design.md §12 risks 2/3). Consequences: `cases/reference_data/naca0012_m080/` is UNTOUCHED (hard rule 6; its README already flags the Euler-anchor provenance); the P4 gates stand as **Picard-quality/robustness gates** (the machinery they lock — upwinding, damping, continuation — is the production warm-start engine and remains bit-locked by G4.2/G4.3); physical acceptance for Newton-era results lives in the G8.1 regression locks. P5's M6 results are not invalidated by construction (different mesh/flow, weaker section shocks), but their Picard-convergence caveat is now known to be of the same kind — the P8/N6 M6 Newton run will measure it. Original record follows. **Closed same day it was found open by audit** (opened AM, re-closed PM): the medium G4.1 gate had diverged on its first actual run (M_max 30.1, Kutta \|F\| 2.5e-3, spurious shock 0.802, cl_p/cl_KJ sign-inconsistent, 19331 total Picard iterations across 12 exhausted Γ evals × 2 levels). Root-cause verified same day (diagnosis follow-up in the G4.1 gate entry): the shipped global mass-lumped pseudo-time damping diag(m_lumped/Δτ) weakens ~4× coarse→medium at fixed Δτ, and the damping needed also grows with shock strength (finer Mach-continuation steps alone ruled out). Fix landed same day: `solve/picard.py::solve_subsonic_lifting`'s new `damping_theta` param (D = θ·diag(A_free), θ=0.2, recomputed every outer iteration from that iteration's own operator — mesh/shock-independent by construction), wired as `solve/continuation.py::TRANSONIC_DEFAULTS`'s new default, mutually exclusive with the retired `pseudo_dt`. Medium gate now **PASSES in 16m39s** (vs 2h43m divergent): upper shock x/c 0.633 (band 0.62±0.03), Kutta \|F\| 1.23e-4 < 2e-4 tol, M_max 1.366, zero limited/floored cells, cl_pressure/cl_KJ sign-consistent (0.349/0.354), n_picard_total 12931. G4.2 bit-identity and the G4.1 coarse shock position (0.604 vs prior 0.599) both re-verified first, per plan; the G4.3 10-case sweep re-run and stays green (one recorded, non-gating difference: the M0.82/α=1.25° corner's cl moved 0.389→0.458 under the new damping). Full default suite unaffected: 136 passed + 2 skipped + 2 xfailed, ~5 min. Same session, the P5 blocker flagged during the diagnosis is also closed: `constraints/wake.py::WakeConstraint.update_matrix`'s per-station `h_j` loop is now one batched `T^T @ (A @ G)` sparse product instead of one matvec per station, verified bit-identical on the real ONERA M6 coarse mesh (83 stations) — removes ~166 extra matvecs/density-iteration on the M6 medium mesh. Transonic convergence semantics unchanged (engineering-converged, not 1e-10; P7 Newton remains the designed cure for the residual tail). Demo `cases/demo/p4_transonic/` re-run, 10 PASS with the new default; demo_report.md §P4 updated with two new addenda. |
| P5 | ✓ | 2026-07-08 | **★ Convergence-quality note (2026-07-11, P8/N6 measurement — record only, gates stand):** the P5 Picard states are NOT discrete solutions (P4-erratum in kind, milder in degree): coupled Newton residual at the committed states 8.6e-6 coarse / 7.6e-6 medium (Kutta \|F\| ~5.5e-4); the Newton true solutions sit at cl_p 0.2560/+5.8% coarse, 0.2646/+7.9% medium with shock positions essentially unchanged and M_max 2.13 — the under-convergence is in circulation/lift. cl_KJ 0.2692 (vs P5 0.2499) narrows the P9 inviscid-vs-Tranair(0.288) gap 0.043→0.019. P5 gates stand as Picard-quality gates (the warm-start engine); Newton-era physics acceptance lives in the G8.2 regression locks. Original record follows. **Closed the same day the medium `physical` failure was RE-diagnosed a second time — the earlier "TE discretization singularity, NOT a wake/Kutta change" conclusion was OVERTURNED by four targeted experiments (T1–T4, `cases/demo/p5_onera_m6/INVESTIGATION_kutta_closure.md`).** Infrastructure (unchanged from the in-progress entry): `post/surface.py::planform_area` + `cl_kj_3d` (CL_KJ=2∫Γdz/(U·S)); `post/section_cut.py::section_cp_curve` (`wall_cp_curve` refactor bit-identical); `solve/continuation.py` forwards `rtol`/`maxiter` (**rtol=1e-7 ~5.5× faster, M_max identical to 5 digits**; default 1e-10 keeps P4 bit-identical); viscous AGARD `reference_data/onera_m6_experiment/` as qualitative overlay; all runs cap `NUMBA_NUM_THREADS=16`. **Second re-diagnosis (T1–T4):** T1 straddle census — 0 tets contain both master+slave wake nodes (cut topology clean). T2 — under-relaxing the eval density at the SAME cached Γ improves drho 4× but the 18-cell outboard-TE M>2 cluster stays bit-identical (not under-convergence). T3 (decisive) — the per-station Kutta mismatch was a SINGLE-station anomaly (st133, z/b=0.801, left 32% under-circulated, |F|/Γ=37% vs ≤5% at all 165 other stations); setting only that station's Γ to its own Kutta target collapses the cluster 18→0 cells (band M_max 3.10→1.16), so the amplitude was a **Kutta-closure failure**, not the TE discretization (a 1/h singularity cannot depend on Γ). Root cause: the per-station secant does not converge at the top Mach level on the 3D mesh (pushing to 16 evals diverges, M_max≈29 — the 10-eval budget was early-stopping regularization); the Γ deficit drives a real TE overspeed that the ρ̃-floor/limiter freeze into a spurious M≈3 state. The far-field 8-cell cluster (incl. the pre-fix M_max=5.20 cell) is the independent span-uniform-2D-vortex branch-ray artifact beyond the tip. **Fix (recipe-level, defaults off, all non-3D paths bit-identical):** `farfield_spanwise_gamma=True` (Γ(z)-tapered vortex far field, 0 at/beyond the sheet tip — `constraints/dirichlet.py`) + `n_kutta_polish=4` (fixed-Γ Kutta-closure polish after the continuation: apply the measured target, re-solve with `omega_rho_polish=0.5`, repeat; secant-free and contractive, |F| halves per step). **Gate numbers (from-scratch, demo 16/16 PASS):** COARSE M_max 1.398, 0 floored/limited, shocks x/c 0.596/0.570/0.425 (η 0.44/0.65/0.90), Γ 0.097→0.0206, V6 2.40%, CL 0.2419, |F| 5.3e-4. MEDIUM M_max 1.995 (bounded tip-TE-corner P1 overshoot — the only surviving singularity trace), **0 floored/limited**, shocks 0.594/0.526/0.345 (~1 cell), Γ 0.097→0.0151 (st133 dip healed), V6 1.82%, CL 0.2453 (coarse→medium 0.2419→0.2453), |F| 5.8e-4 (28× tighter than pre-fix). **G5.2 V6 re-spec (user-approved 2026-07-08):** V6 is a systematic O(h) CL_p-below-CL_KJ discretization floor (sharp-TE/LE P1 wall gradient + P4 Cp sawtooth — both P6 targets; removing the M>2 clusters left V6 unchanged), so it is REPORTED against a 3% floor bound and **<1% is deferred to post-P6 (re-measure after G6.1)**. Known-robustness item (recorded, not fixed): swept-TE Kutta probe assignment shares probe nodes between adjacent stations (st133/134 share their upper probe; 35 upper/41 lower of 166 stations — `diagnose_medium.py` audit) — the latent reason st133 specifically stalled. Dead routes stay dead: spanwise-Γ smoothing (A–E) and TE-element treatments for this gate. Tests: `tests/test_p5_onera_m6.py` 4 fast + 2 gated (updated to the polish recipe + re-specced V6 bound); demo `cases/demo/p5_onera_m6/` (16 checks) + refreshed committed PNG/CSV evidence (the solution npz caches are LOCAL/gitignored like the .msh, now storing residual/drho histories; the demo re-solves absent caches). Full default suite green: **140 passed + 4 skipped + 2 xfailed**. |
| P6 | ✓ | 2026-07-08 | **Surface-pressure recovery (sawtooth removal). ★ N1 root-cause correction (2026-07-08): the surface-Cp sawtooth is a per-triangle wall-gradient RECOVERY artifact, NOT the artificial-density flux.** Decisive evidence: on the *same* coarse NACA0012 M0.80 walk solution, nodal/edge-neighbour smoothing of the wall gradient drops the G6.1 sawtooth metric **330×** (0.0758→0.00023, 39→1 reversals) with shock 0.604→0.607 and cl_p −0.3 % and |Cp|_max unchanged; whereas the smoother streamline-kernel **flux** does not reduce the metric at all (equal-or-worse at every reach/C, even with the shock-robust metric and the solution matched). This **overturns** the earlier ρ̃-selection-flip attribution (demo_report §P4). **G6.1 fix (landed):** `post/surface.py::smooth_wall_tangential_gradients` — a few Jacobi passes averaging each wall triangle's gradient with its edge neighbours, **gated by outward-normal alignment so it never averages across the sharp TE**; threaded as `smooth_passes` through `wall_cp_curve` / `section_cp_curve` / `wall_force_coefficients` (default 0 = bit-identical). Linear in the gradient (differentiable). Metric re-specced shock-robust (sign-alternating second difference — counts only slope-reversal points, so the monotone shock is not counted). Tests: `tests/test_p6_recovery.py` (4), `tests/test_p6_cp_metric.py` (7). Demo `cases/demo/p6_surface_recovery/` re-runs the P4 (NACA M0.80) + P5 (M6, gated) cases showing raw-vs-smoothed Cp + the kernel negative result. **The kernel (differentiable flux) is moved to P7 as an optional Picard-speed path** (opt-in `mode="kernel"`, `upwind_c`≈2.0–2.5, ~10× faster/iter, stable, C^∞ in ∇φ; **re-scoped 2026-07-08** — the P7 Newton prerequisite is the frozen-walk ∂ρ̃/∂φ, sparse Jacobian, NOT the kernel, whose Newton Jacobian is denser) — it does not fix the sawtooth. **P4/P5 impact:** solutions unchanged (shock/M_max/cl_KJ identical, post-processing only); only reported cl_p (−0.3 % coarse) / cd_p (~15 % coarse, untrusted FP quantity) change under smoothing, re-recorded in the P6 demo. **Earlier N1 progress (2026-07-08):** delivered the G6.1 metric (`post/section_cut.py::cp_oscillation_metric` + `tests/test_p6_cp_metric.py`, 6 tests) and a differentiable upwind-density operator (`kernels/upwind.py`, `tests/test_p6_weighted_flux.py`, 6 tests + G4.2 no-op green). **Key implementation finding (design.md §3.2):** the literal Eq. 3.4 near-neighbour blend is **transiently unstable** on the prism-split sliver tets (a fast frozen-Γ probe from the walk's converged coarse field blows up to M_max 20–40 for every p/gain; diagnostic: not dissipation magnitude — the near ring reaches only ~0.3·extent and averaging destroys the upwind character; face-normal weighting also gives anti-dissipative negative-reach cells → use centroid-displacement weighting). The **working operator is a multi-ring streamline-Gaussian kernel** (Eq. 3.4′: depth-3 BFS neighbourhood, Gaussian weight centred on the point ~1 streamwise extent upstream) — genuine reach + smooth, C^∞ in ∇φ, no per-iteration walk so ~10× faster/iteration; probed stable, converging to the walk's own M_max 1.37 / 0 floored-limited state. **Calibration still open:** at reach 1.0/σ 0.35 the kernel converges to a *different* solution (coarse shock 0.604→0.641, cl_KJ 0.364→0.414 +14%) and the raw G6.1 metric doesn't drop — the metric must exclude the shock foot (isolate the sawtooth) and the kernel dissipation must be matched (reach/C) to reproduce the walk's shock/cl before the smoothness gain is measured. **Shipped default stays the P4 walk** (`upwind_weighted=False`, all committed P4/P5 gates bit-identical); the kernel is opt-in (`mode="kernel"`); params threaded through `solve/picard.py` + `solve/continuation.py`; demo metric-instrumented. **Design-pass background (2026-07-08):** **New phase (added 2026-07-07; numbered ahead of Newton by dependency order).** Remove the non-physical ≈2-cell surface-Cp sawtooth introduced by the P4 artificial density: replace the discrete integer-walk upstream selection u(e) + `max(ν_e, ν_u)` switch with a directionally-consistent, C¹ upwind-density operator that keeps the G4.2 subcritical no-op and the (M²−1)/M² dissipation floor. Root-cause analysis + coarse/medium evidence in demo_report §P4 "supplementary analysis". **Design pass (2026-07-08, design.md §3.1–3.2; audit of `kernels/upwind.py` + López dissertation Appendix B verified against the PDF):** the phase has **two independent defects** — (A) the sawtooth is a *selection-flip* spatial-consistency artifact (adjacent supersonic cells pick different integer u(e); O(h), present in Picard); (B) non-differentiability blocks the P7 Jacobian. Key finding: (B) does NOT need a differentiable *selection* — López freezes u(e) per Newton step and differentiates only through ρ_up/μ/ρ (Appendix B), and his own `max(0,μ,μ_up)` is C⁰ yet converges quadratically; the two hard clamps are inactive at the converged P5 state (0 floored/limited) so their non-smoothness never enters the converged Jacobian. **Selected route:** the streamline-projected weighted upwind density (design.md Eq. 3.4, `ρ_up=Σ_f w_f ρ_{nb(f)}/Σ_f w_f`, `w_f=max(0,−V_e·n̂_f)^p`) — kills (A) via a smooth spatial blend and (B) via C¹-in-∇φ at fixed neighbour set, in one operator; smooth the inner `max` with `max_ε`; the smooth density clamp (note-4 "N0") is **optional** (deferred unless a Newton transient stalls on it). (This dense history predates the 2026-07-08 Track-P renumber; the flux role it describes is now **P7**, the Newton role **P8**.) Closed gates: G6.1 sawtooth metric 0.0758→0.00023 (330×) under recovery smoothing, G6.2 physics preserved (post-processing; shock 0.607, cl_p −0.3%, TE ok), G6.4 no-regression (`smooth_passes=0` bit-identical, suite 157 passed). Open G6.3: re-measure V6 under smoothing (residual floor → P9). Demo `cases/demo/p6_surface_recovery/` (renamed from p6_diff_flux). |
| P7 | ✓ | 2026-07-10 | **CLOSED — frozen-selection ∂ρ̃/∂φ of the shipped walk flux, FD-verified at noise level.** Scope as re-scoped 2026-07-08 (core G7.3 + demo; `max_ε` NOT implemented — forward flux byte-untouched, so G7.1/G7.2 hold by construction, locked by V0/G4.2 re-runs + a forward-path regression guard). Delivered: `physics/isentropic.py::mach_squared_derivative_wrt_q_sq` (dM²/dq², strictly positive); `kernels/upwind.py::rho_tilde_sensitivities_sweep` + `UpwindOperator.rho_tilde_sensitivities` (walk mode only, kernel raises NotImplementedError) — exact branch-wise `(s_e, s_u) = (∂ρ̃_e/∂q²_e, ∂ρ̃_e/∂q²_u)` at frozen u(e): subsonic (ρ'_e, 0); accelerating ν=ν_e (López B.3/B.4): (ρ'_e(1−ν) − (ρ_e−ρ_u)ν'_e, ν ρ'_u); shock-point ν=ν_u: (ρ'_e(1−ν), ν ρ'_u − (ρ_e−ρ_u)ν'_u); floored/self-upstream flat branches → 0 exactly mirroring `rho_tilde_sweep`'s clamp. The DOF chain `∂q²/∂φ_k = 2∇φ·∇N_k` stays with the caller (P8 Term-2/-3 assembly). **G7.3 measured:** JVP-vs-FD (against the SHIPPED `rho_tilde_sweep` at frozen u) max rel err **3–5e-10** on the cube in every regime (`tests/test_p7_diff_flux.py`, 8 tests, also green under PYFP3D_NOJIT=1), **3.5e-9** constructed multi-regime NACA-coarse field, **5.7e-9 on the real converged G4.1 M0.80 field** (pocket 1189 accelerating + 977 shock-point, M_max 1.3729). **Sign correction:** design.md §6.3 ∂μ/∂φ chain fixed to `dμ/dM² = +M_c²/M⁴` (FD-arbitrated; the "−" was a typo). **Kink findings recorded for P8:** (i) the max(ν_e,ν_u) tie + switch threshold are the C⁰ measure-zero locus — FD probes straddling them read branch averages (~1e-5), NOT derivative bugs; tests/demo exclude the ε-neighbourhood (0.04 % generic, 2/16.4k converged); (ii) trap: separable fields on structured/prism-split meshes park whole slabs exactly on the tie (same-cell tets share gradients) — use generic/noise-broken fields for any future FD check; (iii) self-upstream elements (u=e, inflow boundary) have zero upwind jump — no kink, ∂ρ̃/∂q² = ρ'_e. Demo `cases/demo/p7_diff_flux/` 7/7 PASS (V7.1–V7.4, committed). Suite: 165 passed + 4 skipped + 2 xfailed. Forward paths untouched: P4/P5/P6 numbers bit-identical. (Pre-close design record: **Differentiable artificial-density flux at frozen selection — the P8 Newton prerequisite** (split out of the old P6 on 2026-07-08; **re-scoped later that day**). **Re-scope:** the Newton prerequisite is the *existing walk flux made differentiable at frozen selection*, NOT a new flux — at frozen u(e) the walk is already C¹ (isentropic ρ_e/ρ_up + branch-wise ∂ν/∂φ, López B.3–B.6); its only kinks are the two `max`, and López keeps exactly that C⁰ switch and still converges quadratically (design.md §3.1). So the deliverable is **derive + FD-verify ∂ρ̃/∂φ for the shipped walk** (upstream coupling ~+1 element/row, sparse, closest to López), with `max_ε` on the **inner** max as optional churn-robustness (never the outer clamp — that would break the subcritical no-op). **The multi-ring streamline-Gaussian kernel** (`UpwindOperator(weighted=True, mode="kernel")`, depth-3 BFS + Gaussian, C^∞ in ∇φ, ~10× faster/iter, `tests/test_p6_weighted_flux.py`, params threaded through `solve/picard.py`+`solve/continuation.py`) is **demoted to an optional Picard-speed path** — N1 proved it neither fixes the sawtooth (P6 recovery artifact) nor is needed for differentiability, and it gives a materially **denser** Newton Jacobian (whole depth-3 neighbourhood); it is a P8 Newton-flux candidate only if N2 measures it net-favourable. Gates: G7.1 subcritical no-op (inner-max `max_ε` only), G7.2 solution preserved (smoothing negligible; no `upwind_c` recalibration), G7.3 (core) FD-verify ∂ρ̃/∂φ at frozen selection. Kernel `upwind_c`≈2.0–2.5 calibration deferred, off the Newton path.) |
| P8 | ✓ | 2026-07-11 | **N6 + G8.2/G8.3 landed 2026-07-11 — PHASE CLOSED.** ONERA M6 Newton: M0.84/α3.06 reachable on both meshes (no fold), coarse 42 s / **medium 249.2 s end to end (G8.2 < 300 s)** via the lagged-LU direct mode (`direct_refactor_every`: true-3D splu is 18.6 s/refactor at 63k dofs — 97% of the 1606 s every-step run; reuse = stale-LU-preconditioned GMRES on the fresh coupled operator at rtol 1e-8, refactor fallback on GMRES failure, default 1 bit-identical) + `NEWTON_M6_RECIPE` (dm 0.05, refactor_every 1000, spanwise-Γ far field); suite **5m02s** (G8.3 < 10 min, 182+8+2). **P5-caveat measured:** Newton residual at the committed P5 Picard states 8.6e-6/7.6e-6 (coarse/medium; Kutta |F| ~5.5e-4) — better than the P4 stall's 2.2e-4 but not solutions; true solutions cl_p 0.2560/+5.8% coarse, **0.2646/+7.9% medium** (cl_KJ 0.2692 — the P9 inviscid-vs-0.288 gap narrows 0.043→0.019), shocks essentially unchanged (0.596/0.541/0.362), M_max 2.13; under-convergence is in circulation/lift, not shock position. Medium solution identity across continuation paths: |dφ| 2.6e-4 / cl to 5 digits (assignment-floor semantics, stale 84 vs 92 — lock bands cover it). Evidence: gated `test_g82_m6_medium_newton_end_to_end`, demo part 3 (V8.1c/V8.2 + g82_m6_medium.csv). Prior: **N5 + G8.1 landed 2026-07-11.** Transonic Newton robustness chain (see the N5 sub-phase entry): direct exact steps (splu+Woodbury — η-accurate Krylov steps stall on the refining shock-position soft mode), stall-adaptive freeze of the upwind assignment + active-set refresh with two-cycle acceptance and honest `residual_unfrozen` floor reporting (the 2.5D prism family parks ~1e3 elements in the max(ν_e,ν_u) near-tie band — live Newton limit-cycles at 3e-9…1e-5 depending on level), freeze-revert/fail-fast/best-of-tried safety nets; recipe `NEWTON_TRANSONIC_RECIPE`. **Two baseline findings en route (user-arbitrated 2026-07-11):** (1) the P4 Picard states are not discrete solutions (P4 ledger erratum — Newton residual 2.2e-4 at the committed coarse state; true coarse M0.80 solution = shock 0.658/cl 0.459/M_max 1.408); (2) the medium mesh steepens into the FP non-uniqueness fold (M0.775 cl 0.396 → M0.7875 cl 0.523, no reachable solution at M0.80) ⇒ **G8.1 re-specced to coarse M0.80 + medium M0.7875** with regression-lock physics bands (gate entry). G8.1 evidence: gated `test_g81_terminal_quadratic_{coarse_m080,medium_m07875}` + gated Newton-field FD pocket + demo `cases/demo/p8_newton/` (V8.1a/V8.1b + V8.2-lite runtime breakdown; medium gate run ~100 s vs Picard's 16m39s non-solution). Kutta closed to machine precision by the coupled solve (|F| ~1e-16 vs secant 1e-4). Prior: **Subsonic milestone landed 2026-07-10 (N2 Jacobian + N3 GMRES + N4 coupled driver; session-scoped by design — transonic tuning deferred to N5 with fresh gated G4.1 runs).** Delivered: `kernels/jacobian.py::assemble_newton_jacobian` — exact (6.3) at frozen selection, Terms 1+2 fused on the SHARED Picard CSR pattern/coloring (Term-2 footprint = Term-1; Term 2 added only when s_e≠0 so masked/limited elements reduce to the Picard matrix **bitwise** — a fused expression FMA-contracts differently under fastmath, measured), Term 3 as **active-set COO** (16 entries per s_u≠0 element, rebuilt every Newton step so upstream-selection churn cannot corrupt a reused pattern; `newton_nnz`/`n_term3_active` recorded = the N2 measurement; **no recolor/wider-CSR machinery needed** — COO build is ms-scale at pocket sizes); speed-limiter consistency via lim-masking s_e/s_u (flat clamp ⇒ derivative 0; bitwise no-op at 0 limited). `wake.py::reduce_operator(A)→(TᵀAT,TᵀAG)` pure extraction (`update_matrix` delegates, bit-identical); on J the H column is the exact wake-jump ∂R_red/∂Γ. `solve/linear.py::solve_gmres` + `build_ilu_preconditioner` (Jacobian nonsymmetric supersonic — asserted in tests). `solve/newton.py`: `NewtonWorkspace` (ONE shared eval_residual state-reconstruction path; affine far-field basis `vals_red(Γ)=vals0+V_red·Γ` extracted by unit-Γ probing with a machine-precision linearity guard; sparse Kutta row K with shared-probe rows correct by construction) + `solve_newton_lifting` — **exact δΓ elimination** ((2,2)=−I): `(J_ff+B·K)δφ=−R−B·F`, `B = J_red[free,dir]@V_red + H_J[free,:]` (**the easy-to-miss far-field vortex column is in and FD-guarded** by `test_gamma_column_fd`), GMRES on the low-rank LinearOperator preconditioned by AMG on the SPD Term-1 block (rebuilt every 2 steps), Eisenstat–Walker choice-2 forcing, safety-only backtracking (no damping_theta anywhere — Picard stabilizer), optional consistent `ptc_dtau`; convergence refused while limited/floored>0. `solve_newton_transonic` = upward-only Mach-continuation skeleton (warm start from last CONVERGED level, dm halving, `upwind_c_post` for López's post-ramp μ_c 2.0→1.6; interfaces final, **untuned — N5**). **Measured (NACA coarse M0.5/α2):** seeded Newton 2 steps 4.2e-7→3.8e-9→**4.7e-13**, cold-start 3 steps 4.5e-4→…→4.2e-11 with observed orders **[2.57, 1.79]** and terminal step gaining >3 digits; Γ matches the P3 Picard/secant to **1.9e-8**, φ to 2.0e-8, cl to ~1e-7 (same discretization); m_inf=0 → ONE step ≡ P2 Laplace. **Jacobian JVP vs frozen-selection residual FD: rel ~1.3–1.5e-10** (tol 1e-6) in all regimes incl. long-range Term 3 (938 active on the constructed cube pocket); Γ-column and Kutta-row FD green; forward Picard paths byte-untouched (structural bit-guard test). Tests: `test_p8_jacobian.py` (7 + 1 gated converged-pocket FD for the G8.1 clause), `test_p8_newton.py` (8); all green under JIT and PYFP3D_NOJIT. Suite **180 passed + 5 skipped + 2 xfailed** (baseline 165+4+2 + 15 new + 1 new gated skip, zero regressions). FD-kink protocol inherited from P7 (generic noise-broken fields, element-level ε-guard lifted to residual rows). Prior design record: design pass done 2026-07-08. Performance & robustness → **fully-coupled Newton** (was P6; renumbered 2026-07-07 so the differentiable-flux phase precedes it — the exact Jacobian is only well-defined on a C¹ flux). **Design pass (2026-07-08, design.md §8.1; López dissertation Ch.3–4 + Appendix B verified against the PDF, supersedes the earlier internal note where it conflicts):** (1) **full Jacobian "strategy A"** — keep the switching-function derivative ∂μ/∂φ (Eq. B.3/B.6) and the nonzero upstream coupling (B.4) → strict quadratic (Tables 4.5/4.6/4.9); dropping Term 3 → only superlinear. **Stencil is NOT "one element-layer wider" in UP3D** (that is López single-hop): the width depends on which P7 flux P8 differentiates — frozen walk ~+1 upstream element/row at graph-distance ≤4, vs kernel coupling the whole depth-3 BFS neighbourhood (denser); N2 measures both (design.md §6.3), old "~30% memory" retired. (2) **fully-coupled (φ_red, Γ)** solve, not the Γ-secant (the secant–density coupling was the P5 medium instability); the Γ-Jacobian blocks are nearly in-code — `∂R_red/∂Γ_j=TᵀJ g_j` (≈`wake.py::self._h`), `F_j=kutta_targets_j−Γ_j` with sparse ±1/n_j `∂F/∂φ`, **plus the easy-to-miss far-field column** −A_coupling·∂vals_red/∂Γ_j (Γ enters the vortex Dirichlet data too). (3) **load-stepping — two PDF corrections to the earlier note:** within a case M_crit/μ_c are held **fixed** and only M∞ ramps (Tables 4.7/4.8, no per-step M_crit sweep); ONERA M6 (Table 4.13) holds M_crit constant 0.95, μ_c at 2.0 during the M∞ ramp, then decreases μ_c 2.0→1.6 at fixed M∞. (4) **also refuted:** the note's premise that "López used a blunt TE" — Eq. 4.2 *sharpens* the NACA0012 TE and he still converged quadratically, so a sharp 2D TE is not itself a Newton obstacle. Expected ≈60–110 Newton iters on M6 (vs ~10⁴ Picard); CL ref 0.288 (KRATOS=Tranair, Table 4.15). Sub-phases N0(optional)/N1 differentiable flux/N2 Jacobian/N3 GMRES+AMG/N4 coupled driver/N5 transonic+load-stepping/N6 M6+perf. **N1 = the P7 frozen-selection differentiable walk flux** (∂ρ̃/∂φ FD-verified, `max_ε` optional) — the **default** Newton flux is the sparse frozen walk (~+1 upstream element/row, closest to López). The **streamline-Gaussian kernel** (`UpwindOperator(weighted=True, mode="kernel")`, C^∞, ~10× faster/iter, built during P6) is the **optional** denser-Jacobian alternative, promoted only if N2 measures it net-favourable — it does not fix the sawtooth (P6 recovery artifact, §9.1) and is not required for differentiability. (Renumbered old-P7 → **P8** on 2026-07-08 when the flux became its own P7.) Gates G8.1 (terminal quadratic on G4.1 + FD-verified Jacobian), G8.2 (M6 medium < 5 min), G8.3 (suite < 10 min). |
| P9 | ✓ | 2026-07-11 | **Grid-convergence & accuracy-gap discrimination — evidence phase, NEW 2026-07-11 (user-directed renumbers, net same-day mapping: curved walls → P11, backlog → P12).** Purpose: test the 0.019 M6 lift-gap attribution (cl_KJ 0.2692 vs Tranair/KRATOS 0.288 → "sharp-TE/LE P1 wall floor") BEFORE the large P10 effort — the P8 capability demo showed it is inference-grade (2D sharp-TE lift converges cleanly to the panel bracket; M6 cl_KJ still rising 0.2621→0.2692 under refinement, no sphere-style saturation). Gates with PRE-REGISTERED bands: G9.1 M6 coarse/medium/fine Newton + Richardson (cl_KJ∞ ≥ 0.283 ⇒ resolution-dominated; ≤ 0.278 ⇒ floor confirmed; else inconclusive); G9.2 NACA 2.5D coarse/medium/fine subsonic lift oracle (monotone error decrease + fine within ±1% of the PG–KT midpoint); G9.3 attribution verdict + P11 (curved walls) go/no-go recorded here (user arbitrates). No solver changes; fine meshes gitignored; heavy runs one-shot (CI budget untouched). **CLOSED 2026-07-11 — demo `cases/demo/p9_grid_discrimination/` (11 PASS + 3 documented XFAIL). The pre-registered bands could NOT fire, and why is the phase's finding.** **G9.2 PASS (clean):** 2D sharp-TE lift error vs the PG–KT midpoint **2.71% → 0.33% → 0.03%** (all three converged) — **a sharp TE imposes NO lift floor**; the 2D leg of the "sharp-TE/LE P1 wall" story is gone. **G9.1 INVALID (no 3-point Richardson):** the fine mesh is **not a discrete solution** — max local Mach (unlimited) **1.40 → 2.13 → 7.93** with **0/0/9** cells crossing M_cap=3; permanently-limited cells block the N5 freeze machinery (needs 0 limited), so fine limit-cycles at |R|~1e-5 for its full budget and its cl_KJ 0.2393 is an ARTIFACT, not a lift (coarse 0.2621 / medium 0.2692 are the only discrete solutions). **★ The divergent singularity sits ON THE WAKE SHEET at its free tip edge** (all 9 capped cells: z/b 0.998–1.000, x−x_TE +0.002…+0.017 i.e. AFT of the TE, |y|<0.003 in the chord plane) — where Γ(tip)=0 is enforced: the classic **vortex-sheet-edge singularity of a rigid planar wake** (real flow rolls up into a tip vortex). P5's "bounded tip-TE-corner P1 overshoot" (M_max 1.995) is the SAME object, just under-resolved. **⇒ a WAKE-MODEL defect, NOT a wall-element defect — curved elements do not remove a sheet edge. ★ Correction 2026-07-12: the fix is wake roll-up / a tip vortex, NOT Track B (which changes the wake REPRESENTATION, not the model — `cases/demo/p13_tip_edge_singularity/` measured the tip-edge peak Mach diverging under refinement on BOTH the conforming and level-set paths at M0.5; no current Track B phase does roll-up, B9 is shelved).** **G9.3 verdict:** the 0.019 gap is **UNSPLITTABLE as posed** (resolution/floor shares recorded as `n/a`, not fabricated); **P11 is NOT supported by P9 as the 3D-lift fix** — its remaining case is G1.6 (smooth curved wall, different mechanism). 3D accuracy now points at the tip/wake route (Track B), which must land before any 3D grid-convergence claim. **User arbitrates.** **Solver-path findings en route (feed P10):** (a) `precond="direct"` does not scale — one splu at 450k dofs ran **4h39m without returning** (26 GB, killed) vs 18.6 s on medium; (b) the `precond="amg"` fallback is valid AND faster at every size **once EW forcing is tightened to η=1e-8** (medium 66 s vs 141 s direct, same solution to 4 digits, 0 stalls; coarse 8 s vs 42 s) — N5's "Krylov stalls" is a property of the LOOSE forcing, not of AMG; (c) the fine cold Picard-5 seed overshoots the LE into the density floor (M0.70: 4036 lim/1847 flr, level-0 breaks) — fixed path-only with m_start 0.30 + n_picard_seed 12 (13 levels, 5294 s, in budget). |
| P10 | ◐ (G10.2 ✓ 2026-07-11; G10.1 open) | | **Newton generality & continuation efficiency (NEW 2026-07-11, user-proposed).** **G10.2 CLOSED 2026-07-11 — SPLIT A/B verdict** (demo `cases/demo/p10_newton_usability/`, run under the 16-thread timing protocol): `solve_newton_transonic(intermediate_tol=…)` opt-in, default None bit-identical (suite-locked; `solve_newton_lifting` grew `tol_residual_loose`/`tol_residual_rel`/`accept_on_stall` + `accept_reason` reporting; level_results record per-level accept_reason). Shipped rule = the pre-registered candidate + two A/B-measured hardenings: (1) loose acceptance requires ≥ 1 Newton step at the level (round-1 finding: warm-started levels ENTER below any absolute threshold — zero-step acceptance degenerates the ramp into a level skip); (2) dm-halving retry levels run STRICT (the halving cascade is the robustness fallback). **(a) M6 medium: all G8.2 locks intact, final level converges identically (12 steps, |R| 7.8e-15, cl/M_max/shocks equal to 4 digits), solve 239.5→140.3 s (+41.4%, intermediate levels 35→6 Newton steps: 1–3 loose steps each ending ~1e-5) ⇒ `intermediate_tol=1e-5` PROMOTED into `NEWTON_M6_RECIPE`** (gated G8.2 test now runs the adaptive path, ~145 s). **(b) NACA medium M0.7875: NEGATIVE result recorded — the P8 "warm-start only from CONVERGED levels" trap measured in G10.2 form:** near the fold (dcl/dM ~6–10) the loose ramp's 1–4-step levels never track Γ/shock, the final level and even STRICT halving-retry levels stall at the ~5e-6 live-churn floor for 60 steps each (round 2: cl 0.369 vs lock 0.523, unconverged) ⇒ `NEWTON_TRANSONIC_RECIPE` unchanged (strict); **loose intermediates are contraindicated in fold zones.** Suite +2 (`tests/test_p10_continuation.py`): default-path accept_reason lock + subsonic-ramp adaptive path (Γ matches strict to 1e-6, steps not worse). Baseline 184+8+2. **G10.3 CLOSED 2026-07-11 — verdict KEEP the ramp** (no-ramp single-level Newton at target M∞, 2 seedings × 4 cases, `run_g103_noramp.py`): far from the fold both seedings reach the SAME solution (class A — branch selection not binding there), but Picard-5 transits clamped states (peak 45 lim+flr, final 0/0, locks pass) at +44% (141.4→79.0 s) failing the pre-registered clamp-free clause, and clamp-free Picard-40 gains only +6.6% medium / −37% coarse failing ≥20%; fold-zone medium is class C under BOTH seedings (s1 stalls cl 0.449≠0.523; s2's un-continued Picard seed diverges, 9934 limited — deep seeds are HARMFUL without a ramp, the P4/P5 record); NACA coarse s1 class-A in 4 s recorded as a single-case exception. The +44%-but-clamped observation is recorded for user arbitration of the clamp-free clause; `clamp_history` added (additive). Remaining: And G10.1 non-lifting Newton entry (`wc=None` workspace, raw mesh; sphere M0.3 matches solve_subsonic to \|Δφ\|∞ < 1e-8, m_inf=0 one-step ≡ Laplace; NOT an accuracy gate — G1.6 untouched), no ordering constraint. Framing correction recorded: the coupled Newton has NO Kutta outer loop (Γ updates every step, ‖F‖ machine-zero) — the capability-demo two-level structure is the Mach ramp, so the user's "loose-then-tight" idea applies to continuation levels. |
| P11 | ☐ (**lift case STRONGLY INDICATED MOOT — NOT EARNED** (wording downgraded 2026-07-14, user-arbitrated; was "REFUTED with evidence 2026-07-13") — P13/G13.2 + G13.3-transonic NEGATIVE; G1.6 case stands as the only remaining justification) | | ★ **2026-07-13 (P13/G13.2); wording downgraded 2026-07-14 (user-arbitrated) — the lift justification is STRONGLY INDICATED MOOT by a real fine solution, NOT refuted per P9 discipline.** Once the tip taper made the M6 **fine** mesh a genuine discrete solution at M0.84 (G9.1's own scenario, where fine had been a limit-cycle artifact), cl_KJ reached **0.2866** (untapered-equivalent ≈ 0.291) against the Tranair/KRATOS **0.288** — 99.5 % of the reference, on the resolution side of P9's ≥ 0.283 threshold. **But the pre-registered band fires on the Richardson EXTRAPOLANT**, and no admissible extrapolant exists: the flat-cap sequence is non-asymptotic (increments grow) and the round ladder has no M0.84 fine state (G13.3 transonic NEGATIVE) ⇒ the *direction* is established with evidence, the *verdict* is not earned on either geometry. **P11's remaining valid case is G1.6 alone** (sphere-Cp on a SMOOTH curved wall — a different mechanism, untouched by any of this). Prior entry follows. **Curved / isoparametric wall elements (renumbered from P9 via P10 on 2026-07-11; opening CONDITIONAL on G9.3):** **P9 outcome (2026-07-11): G11.2's premise is REFUTED as stated** — the 2D sharp TE imposes no lift floor (G9.2, error → 0.03%), and the 3D blocker is a divergent vortex-sheet-edge singularity on the rigid planar wake at the tip (G9.1), which curved WALL elements cannot remove. **G11.1 (G1.6 sphere-Cp on a smooth curved wall) is unaffected and remains the phase's valid justification.** The 3D lift/accuracy route now points at the tip/wake model (Track B). Original scope follows: the shared accuracy route for **G11.1** G1.6 sphere-Cp < 2% (design.md §5.1 Option C / DP1 "> 5%" branch) and **G11.2** the residual V6 < 1% floor left after P6 removes the sawtooth (attributed to the sharp-TE/LE P1 wall gradient — P9 tests this; M6 cl_KJ 0.2692 vs Tranair/KRATOS 0.288). Large own effort per DP1. |
| P12 | ☐ | | Backlog (renumbered from P10 via P11 on 2026-07-11; originally P8): discrete adjoint (transpose of the P8 Newton Jacobian), VII transpiration BC (now expanded into Track V), mixed prism/tet + (C, M_c, ω) BO calibration. (Non-lifting Newton promoted to P10/G10.1 on 2026-07-11.) |
| P13 | ◐ (G13.1 ✓ + **G13.2 conforming ✓** + **G13.3 subsonic ✓** 2026-07-13; G13.2 level-set clause open; G13.3 transonic RAN → NEGATIVE (the round FINE mesh never reaches M0.84 — its Mach ramp dies at M=0.75; site = the sharp tip TE, which the rounded cap AMPLIFIES rather than creates; P9-band verdict NOT earned on either geometry)) | | **Tip / wake-edge singularity (NEW 2026-07-13, user-approved; descendant of P9/G9.1; appended, no renumber).** **G13.2 CLOSED (conforming path):** the fix is a **spanwise loading taper** `Γ_eff(z)=F(z)·Γ_Kutta(z)` on the per-station Kutta target (`constraints/wake.py::tip_taper_factors`; `solve_newton_lifting(tip_taper=…)`, default `None` = bit-identical). Shipped model: `vanish_smooth` (smoothstep, COMPACT support), r_c = 0.05·b_semi. ★ **NOT roll-up, NOT a vortex core, and NOT via Track B** — the earlier plan recorded in this row ("hand G13.2 to Track B as a rescope of the shelved B9") was **superseded**: the mechanism is DISCRETE (the outermost TE station sheds its retained Γ_last over ONE cell ⇒ p ≈ 1 − q, criterion q ≥ 1), so a Kutta-target taper fixes it with no wake-model rewrite. Edge peak p: **+0.592 → +0.009**; the M6 fine mesh becomes a genuine discrete solution; cost ≈ −1.1…−1.6 % cl, LOCAL (inboard of η≈0.95 unchanged). Tests `tests/test_p13_tip_taper.py` (15); demos `cases/demo/p13_tip_edge_singularity/` (G13.1 10/10, G13.2 probe 16+2 xfail, G13.2 physics 10+1 xfail, **G13.3 5/5**). **G13.3 SUBSONIC ✓ — the flat tip cap was a THIRD singularity, on the WALL, and rounding it (Track M / M5) RESTORES 3D grid convergence.** With the tip taper AND the self-similar ladder, lift increments still GREW; a box study localized the remaining divergence to the **flat tip cap's sharp wall edge** (p = +0.321) while the wake edge (+0.045) and wing interior (−0.014) stayed bounded. Track M / M5 (`wing3d.py::tip_cap="round"`, below) replaces it with a half-body-of-revolution cap. Re-run on the round ladder (demo `run_g133_roundtip.py`, **9/9**, M∞0.5, A/B vs the flat ladder): **(1)** the **cap-surface** exponent falls from **+0.327 (flat) → +0.091 (round, bounded)** and the tip region excluding the design-sharp LE/TE is **converged (p = −0.006)** ⇒ the cap edge singularity is GONE; **(2) ★ the three-point Richardson G9.1 could NEVER run is now EARNED** — round cl_KJ 0.2159 → 0.2073 → 0.2055, increments SHRINK (−3.95%, −0.88%), observed order **p = 2.31**, **cl_KJ(h→0) = 0.2050**; the flat ladder's cl was non-monotone (0.2015 → 0.2005 → 0.2121, no Richardson). **(3) HONESTY / metric trap (G13.1 finding 6):** the *broad* G13.1 tip box still "diverges" (p = +0.38), but its max MIGRATED to the **design-sharp tip TE** (fine peak at chord-frac 0.999, present in BOTH families, carries the Kutta condition) — a different, integrable feature that no tip-cap change removes and that does not spoil the integrated lift (which converged). **G13.3 TRANSONIC — RAN, NEGATIVE RESULT (2026-07-13, demo `run_g133_roundtip_transonic.py`, 6/6):** the M0.84/α3.06 three-point Richardson on the ROUND ladder was run (amg + tight EW forcing + m_start/Picard guards). **coarse + medium are clean discrete solutions** (M_max 1.51 / 2.00, 0 over cap, converged; cl_KJ flat-normalised 0.2769 → 0.2763, nearly flat), **but the FINE solve is NOT** — M_max 3.97, **5 cells over M_cap**, 1 floored, not converged; its cl_KJ 0.2415 is a limit-cycle **artifact** (cf. G9.1's 0.2393 non-lift), not a convergence point. — **and it is worse than "the fine did not converge": the fine mesh never REACHES M0.84 at all.** Its Mach continuation **breaks down at M = 0.75** (last converged level M = 0.7375; dm-halves to below dm_min and gives up, 1 cell in the density floor), so **no M0.84 fine state exists** and no census of one is possible. ⇒ **only two of three points exist; no three-point Richardson at M0.84; the P9 band verdict is STILL NOT EARNED.** ★★ **THE SITE IS THE SHARP TIP TE — the round cap does NOT create it, it HEATS it.** The committed **flat** M0.84 run (`run_g132_transonic.py`, 39e1ded, same taper) *completes* the ramp at the same refinement level and converges (M_max 2.818, 0 over cap, cl_KJ 0.2866), which at first seems to exonerate the shared TE. It does not: the field-saving rerun `run_g133_roundtip_transonic_locate.py` puts **20 of the 20 fastest cells ON the sharp tip TE** (z/b 0.97–0.99, chord-frac ≈ 0.998) and **NONE on the new cap surface** (z/b > 1). The rounded tip lets flow **wrap around the tip and accelerate**, raising the velocity at that same, pre-existing sharp TE at **every** level (M_max 1.51 / 2.00 round vs 1.39 / 1.73 flat) — and at fine that pushes it past the limiter where the flat cap's cooler tip flow stays under. **⇒ the cap did not add a singularity; it AMPLIFIED the tip-TE one that was always there.** The round cap remains a clean **subsonic** fix (Richardson p = 2.31); its transonic cost is a *new sub-problem at the tip TE*, not a defect of the cap geometry. ⚠ **METHOD ERRATUM (self-caught):** an earlier version of the transonic demo censused the *failed* (M ≈ 0.75) state at M_INF = 0.84 — the WRONG freestream Mach — and reported a spurious "M_max 3.97 / 5 cells over M_cap / cl_KJ 0.2415 at M0.84". **Those numbers are RETRACTED.** `solve()` now records the level sequence and refuses to census a state whose ramp did not reach the target (`target_reached`, `m_final`, `m_last_converged` in `g133rt_transonic.csv`; every census column reads `n/a` for the fine level). ⚠ **P9 transonic verdict, both geometries:** flat gives converged solves but the sequence 0.2593→0.2652→**0.2866** RISES (not asymptotic; on the CLAMPED non-self-similar ladder AND polluted by the flat cap edge), so 0.2866 is a single-point value, not a Richardson extrapolation; round has no third point at all. **Neither earns the "0.019 gap = resolution ⇒ P11 refuted" conclusion** — it has no clean asymptotic discrete-solution basis on either geometry. **Historical G13.1 record follows.** **Tip / wake-edge singularity — characterization + wake-model fix (NEW 2026-07-13, user-approved; direct descendant of P9/G9.1; appended, no renumber).** **G13.1 CLOSED 2026-07-13 — demo `cases/demo/p13_tip_edge_singularity/` 10/10 PASS.** Probed SUBSONICALLY (M∞0.5, no transonic machinery) to isolate the geometric edge signal. **(1) It is a real 1/√r flat-plate-edge singularity, not "1/r":** the conforming three-point (coarse/medium/fine 55.5k/350.7k/2.51M tets) log-log exponent of tip-box peak Mach vs 1/h is **p = 0.59 ∈ [0.4,0.65]** (peak 0.712→0.981→1.510; a 1/r line vortex would give p=1). **(2) Driver = trailing vorticity dΓ/dz, not bound Γ:** Γ→0 at the tip (necessary-not-sufficient), but |dΓ/dz| is ~10× larger at the tip than mid-span (B7's smooth Γ(z)); the *unloading rate*, not the loading, is what a terminating flat sheet cannot regularize. **(3) MODEL not representation:** the edge diverges on BOTH the conforming (×1.38) and level-set (×2.28) paths on the SAME mesh (LS exponent 1.34 ≥ conforming 0.52) — Track B's representation change does not blunt it; tip-box p95/mean stay FLAT (0.573→0.562→0.525 / ~0.49) while max diverges (localized-edge signature); peak sits AFT of the TE in the chord plane ⇒ not a wall feature (curved walls/P11 cannot fix it). **(4) ★ The conforming FINE M0.5 solve does NOT converge** (limited/floored, ~1.4k NaN cells) — the tip singularity trips the limiter even subsonically, the exact M0.5 analogue of G9.1's transonic non-solution. **This corrects the committed docs' "1/r-type" (roadmap :1036, demo_report :1333) to 1/√r, and adds the dΓ/dz mechanism.** **G13.2 (future):** roll-up / explicit tip-vortex model ⇒ bounded tip peak (p→0) + M6 fine a real discrete solution; implementation handed to **Track B as a rescope of the shelved B9** (level-set naturally represents a movable sheet via `update_direction()`; from O(θ²) deflection → roll-up). **G13.3 (future):** with the model in place, redo the M6 3D three-point Richardson (G9.1's blocked task) under the P9 pre-registered bands. No solver/numerics changes; fine meshes gitignored; the ~34 min conforming fine M0.5 AMG solve is a one-shot cached artifact. |
| P14 | ✓ CLOSED 2026-07-17 (**opened + closed same day, user-directed**; G14.1–G14.7 all ✓; demo `cases/demo/p14_pressure_kutta/` 28 PASS, diagnostic `cases/analysis/p14_te_pressure_diag/` 20/20; wiring scope = coupled Newton + `solve_laplace_lifting` only). **Result: the conforming path now matches the level-set path** — M0.84 medium cl_p/cl_KJ agree to 0.15%/0.34% (V14.6 cross-model), Γ(z) roughness 0.0970→0.0043/0.0365→0.0024 (at/below the LS band), all-station TE Cp gap 0.2206→0.0040/0.1585→0.0024, TE spike 0.1143→0.0533 (below LS). G14.7 re-specced at close from the probe G8.2 locks to the level-set oracle (the +4.85% cl_KJ move is the finding — 69% of P9's 0.019 gap was Kutta-estimator bias); residual ~0.05 shared P1 recovery floor and D=1.80 discriminator recorded, not closed. | | **Probe-free conforming Kutta target: wall-adjacent-CV pressure-equality estimator (NEW 2026-07-17, from A2; appended, no renumber).** A2 (`cases/analysis/a2_te_kutta_fidelity/`, closed 2026-07-17) proved two conforming symptoms are one estimator's fault: **S1** the Γ(z) jitter is a measurement-operator artifact of the per-station probe-difference potential-jump Kutta target (fixed-Γ discriminator D=7.33/25.70 coarse/medium — the estimator regenerates the jitter from a smooth field; closure |F|/max|Γ| ≤ 0.6%, so not unclosed stations); **S2** the TE Cp jump is that same target being equal-*potential* not equal-*pressure* (same-estimator gap 34×/133× vs level-set; P13's 0.14–0.22 confirmed). **Fix (this phase):** replace it with a wall-adjacent control-volume recovered-velocity / pressure-equality estimator — the B4 objects (`wake/multivalued.py::_build_te_control_volumes` + `te_velocities`) ported to the conforming cut mesh; kills the single-node sampling (S1) and closes the gap (S2) in one swap; per-station closure becomes nonlinear (as the transonic `kutta_per_outer=1` path already approximates). Provisional gates: G14.1 Γ(z) roughness → LS band + D→O(1); G14.2 section TE gap → LS band; G14.3 lift/V6 move <1–2% vs the G8.2 locks; G14.4 inert by default (flag off = bit-identical). Dead routes (A2): spanwise-Γ smoothing, full-fan TE recovery (+11–15%), better-probe-picking alone. Risks: conforming duplicated TE nodes + explicit Γ DOF (LS has neither), tip region unexamined, nonlinear closure may need damping — diagnostic-first at open. **Not on the B9 critical path** (level-set has neither symptom by construction). **Trigger:** a conforming accuracy campaign needing smooth Γ(z) or a closed TE (e.g. Track V sectional-Cp coupling, an MDO Γ(z) objective), or user direction. |

