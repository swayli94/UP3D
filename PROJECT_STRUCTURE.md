# Project Structure

This document describes the directory layout and initialization status of pyFP3D.

## Directory Tree

```
pyfp3d/                    # Main package
├── __init__.py           # Package entry point, NOJIT mode flag
├── mesh/                 # Mesh I/O, topology, metrics, coloring
│   ├── __init__.py
│   ├── reader.py         # [P0] meshio → SoA arrays + boundary tags
│   ├── metrics.py        # [P0] volumes, gradients, face adjacency; ✓ [P3] adds
│   │                       #   precompute_element_geometry (B_e/V_e once per mesh);
│   │                       #   ✓ [P4] precompute_face_normals (adjacency-ordered outward)
│   ├── coloring.py       # [P0] element graph coloring; ✓ [P3] greedy loop numba-jitted
│   │                       #   (same visit order => identical assignment) + color_partition_csr
│   └── wake_cut.py       # ✓ [P2] wake-sheet node duplication (per-node flood-fill side
│                           #   classification, no planarity assumption; TE nodes ARE
│                           #   duplicated -- see docstring + roadmap P2 assert re-spec;
│                           #   (x,y)-grouped spanwise stations -> quasi-2D collapses to
│                           #   one scalar Γ; Kutta probe nodes; preprocess topology asserts)
├── meshgen/              # ✓ [M0] Quasi-2D mesh generation (roadmap Track M)
│   ├── __init__.py
│   ├── extrude.py        # ✓ single-layer extrusion; globally consistent
│   │                       #   prism→3-tet split (min-global-index diagonal rule);
│   │                       #   assert_quad_split_consistency (M0 preprocessor assert)
│   └── planar.py         # ✓ vanilla-Gmsh 2D builders: cylinder annulus, NACA0012
│                           #   with wake line embedded via gmsh.model.mesh.embed
│                           #   (gmsh imported lazily; solver tests don't need it)
├── constraints/          # ✓ [P2] Constraint machinery
│   ├── __init__.py
│   ├── wake.py           # ✓ master–slave elimination (A_red = TᵀAT once; Γ RHS-only via
│   │                       #   precomputed per-station vectors); kutta_targets() (per-station
│   │                       #   mean of probe jumps)
│   └── dirichlet.py      # ✓ far-field freestream + 2D vortex correction (branch cut ON
│                           #   the wake sheet; eliminated ⁺-side far-field wake nodes
│                           #   automatically consistent); ✓ [P3] Prandtl-Glauert scaling
│                           #   (beta stretches only the atan2 argument; beta=1 bit-exact)
├── physics/              # Physics constants and constitutive relations
│   ├── __init__.py
│   └── isentropic.py     # ✓ [P0] ρ(q²), M(q²), a(q²), Cp, etc.; [P2] adds
│                           #   pressure_coefficient_incompressible (Bernoulli limit);
│                           #   ✓ [P3] density_field/mach_squared_field array sweeps
│                           #   (ρ ≡ 1.0 BITWISE at M∞=0 -- the G3.3 anchor)
├── kernels/              # Element-wise assembly kernels (Numba-jitted)
│   ├── __init__.py
│   ├── gradient.py       # ✓ [P3] prange velocity sweep: grad(phi)_e + q²_e from B_e,
│   │                       #   zero-alloc (outputs preallocated by PicardOperator)
│   ├── jacobian.py       # ✓ [P3] Picard matrix (6.2): symbolic CSR pattern, elem_to_csr
│   │                       #   scatter map, colored-prange data kernel, PicardOperator
│   │                       #   per-mesh workspace (B_e/V_e/coloring/pattern/buffers once);
│   │                       #   accumulation order fixed by color sequence => bit-deterministic
│   │                       #   across calls/threads; [P6] exact Newton (6.3) lands here
│   ├── residual.py       # [P1] serial reference kernels (KEPT, regression-tested against)
│   │                       #   + ✓ [P3] assemble_residual_colored; assemble_stiffness_matrix
│   │                       #   now delegates to the fast path (P1/P2 drivers share it) → [P6] Newton
│   └── upwind.py         # ✓ [P4] artificial density (3.1)-(3.2): MULTI-HOP upstream walk
│                           #   (single-hop reaches only ~1/3 extent on prism-split meshes --
│                           #   measured dissipation starvation), shock-point operator
│                           #   nu = max(nu_e, nu_up), rho_tilde floor, UpwindOperator workspace;
│                           #   exact bitwise no-op below M_crit (G4.2)
├── solve/                # Linear and nonlinear solvers
│   ├── __init__.py
│   ├── linear.py         # [P1] Dirichlet elimination + CG/PyAMG preconditioner (done);
│   │                       #   ✓ [P3] build_amg_preconditioner pins pyamg's spectral-radius
│   │                       #   RNG seed -- ALL solves bit-reproducible run-to-run (G3.3)
│   ├── picard.py         # [P1] Laplace driver (done); ✓ [P2] solve_laplace_lifting():
│   │                       #   Kutta outer loop, matrix+AMG built once, Γ updates RHS-only,
│   │                       #   secant (Aitken) acceleration -> 2 updates on the linear driver;
│   │                       #   ✓ [P3] solve_subsonic (non-lifting density Picard) +
│   │                       #   solve_subsonic_lifting (NESTED: outer ρ update, inner P2
│   │                       #   secant Kutta at frozen ρ -- interleaved Γ steps rejected with
│   │                       #   measurements; AMG reuse every 4 outers; opt-in forcing-term
│   │                       #   inexact inner solves, default off; M∞=0 bitwise ≡ P2);
│   │                       #   ✓ [P4] upwind_c/m_crit (ρ̃ in matrix+residual), q² limiter
│   │                       #   (m_cap), pseudo-transient diag(m/Δτ), omega_rho, kutta_per_outer,
│   │                       #   phi_init/gamma_init continuation seeds
│   ├── wall_correction.py # ✓ [P1/G1.3] true-normal weak-flux correction RHS (Option A);
│   │                       #   assembly-verified; correction itself RULED OUT by the
│   │                       #   G1.3/G1.4 oracles (design.md §5.1.2) -- kept as reusable
│   │                       #   facet-integral infrastructure (e.g. Gap-SBM terms)
│   ├── continuation.py   # ✓ [P4] Mach continuation + transonic Γ closure: outer per-station
│   │                       #   SECANT on F(Γ) = kutta_target(density-converged φ at fixed Γ) − Γ
│   │                       #   around frozen-Γ pseudo-time solves (nested/interleaved Kutta both
│   │                       #   fail at transonic -- measured; see module docstring)
│   └── newton.py         # [P6] Newton method (PLANNED, file not created yet)
└── post/                 # Post-processing
    ├── __init__.py
    ├── vtk_out.py        # [P0] Write .vtu for ParaView; also the PNG/CSV gate-artifact helpers
    │                       #      (export_error_heatmap, export_matplotlib_plot) live here, not
    │                       #      in a separate artifacts.py
    ├── shock.py          # ✓ [P4] shock monitors: Cp* sonic-crossing detection, shock x/c,
    │                       #      jump-width cell count, monotonicity + expansion-shock detectors
    ├── surface.py        # [P1] nodal_gradient_recovery() (volume-weighted, for interior fields)
    │                       #      ✓ [P3] m_inf param: exact isentropic Cp (2.5) in the force integral
    │                       #      and wall_tangential_gradient() (surface-only, for wall Cp --
    │                       #      the accurate one; see "Known gaps" for why it still isn't
    │                       #      accurate *enough* to close G1.6); ✓ [P2] adds triangle-wise
    │                       #      wall force integration (owner-tet-oriented outward normals,
    │                       #      no nodal averaging across the sharp TE) and KJ sectional cl
    └── section_cut.py    # ✓ [G1.3→P2] z = const section extraction: degenerate single-layer
                            #      path + [P2] general marching-tets interpolation path and
                            #      wall_cp_curve() sectional Cp(x/c) upper/lower split

cases/                     # Test cases and reference data
├── meshes/               # Mesh families (coarse/medium/fine)
│   ├── sphere_shell/     # [P1] Gmsh sphere-shell case for gate G1.6 (coarse/medium generated);
│   │                       #   [P3] reused for G3.1 (compressible-vs-PG same-mesh comparison)
│   ├── cylinder_2.5d/    # ✓ [M0] Single-layer extruded cylinder-flow test case
│   │                       #   (generate_cylinder.py; coarse 6.9k / medium 17.3k tets
│   │                       #   committed; analytic Cp = 1 - 4 sin^2(theta) validation);
│   │                       #   G1.3 pre-study ran here (fine.msh added, 50.2k tets); found
│   │                       #   recovery-dominated (~76%), NOT sphere-crime-dominated --
│   │                       #   de-designated as G1.6 testbed (design.md §5.1.2)
│   ├── naca0012_2.5d/    # ✓ [M0] Single-layer extruded NACA0012 + embedded wake sheet
│   │                       #   (generate_naca0012.py, one parameter h_wall per level;
│   │                       #   coarse 16.4k / medium 61.8k tets committed, fine on demand)
│   └── onera_m6/         # [M1] Swept wing -- not started
├── reference_data/       # Ground truth (DO NOT EDIT)
│   ├── naca0012_incompressible/  # ✓ [P2] Hess–Smith panel reference (generator script +
│   │                             #   cl_reference.csv / cp_alpha4.csv / convergence.csv +
│   │                             #   README provenance; two independent lift routes agree
│   │                             #   to 0.09%, panel-count converged)
│   ├── naca0012_m05/     # ✓ [P3] the same panel solution under Prandtl-Glauert AND
│   │                             #   Karman-Tsien corrections (G3.2 reference = PG/KT
│   │                             #   midpoint + inside-bracket assert; README provenance)
│   └── naca0012_m080/    # ✓ [P4] transonic shock reference: Euler anchor (~0.60c upper,
│                                 #   ~0.35c lower) + documented conservative-FP aft-shift band;
│                                 #   README records that no open FP table was retrievable
├── demo/                 # ✓ Per-phase evidence demos (docs/demo_report.md; one
│   ├── README.md         #   self-checking run_demo.py + committed results/ per phase)
│   ├── _common.py        #   shared chart style + CheckList acceptance recorder
│   ├── p0_infrastructure/  # G0.1-G0.4: volume/gradient exactness, coloring, VTK I/O
│   ├── p1_laplace/         # V0/G1.1/G1.2 + G1.4 oracle (absorbed) + G1.6 open-gate XFAIL
│   ├── p2_kutta_lifting/   # G2.1-G2.5: cut exactness, Kutta, cl vs panel, spanwise decay
│   ├── p3_subsonic/        # G3.1-G3.3: assembly-debt evidence, sphere-vs-PG, cl bracket,
│   │                       #   monotone nested Picard, M=0 bit-identity (14 checks)
│   ├── p4_transonic/       # G4.1-G4.3: subcritical bitwise no-op, upwind-reach evidence,
│   │                       #   coarse M0.80 shock quality vs reference band (10 checks)
│   └── m0_meshgen/         # mesh gallery, hard-rule-7 topology matrix, cylinder convergence
└── test_*.py             # [Deprecated] Integration tests (use tests/ now)

tests/                     # Unit and gate tests
├── conftest.py           # ✓ Pytest fixtures: artifacts_dir (persistent, PYFP3D_ARTIFACTS_DIR
│                           #   overridable), mesh_dir, etc.
├── test_conftest_artifacts.py       # ✓ Regression test: gate artifacts persist in artifacts/
├── test_metrics_degenerate.py       # ✓ Regression test: degenerate-tet guard in metrics.py
├── mesh_utils.py         # ✓ [P1] Dependency-free structured-cube + sphere-shell mesh generators
├── __init__.py
├── test_v0_freestream.py # ✓ [P0/P1] Primary regression test (incl. cut-free residual check)
├── test_mesh_*.py        # [P0] Gates G0.1–G0.4
├── test_mesh_adjacency.py           # ✓ [P0] Regression test for build_face_adjacency fix
├── test_mesh_reader_roundtrip.py    # ✓ [P0] Regression test for write_mesh tag-loss fix
├── test_laplace_mms.py              # ✓ [P1] Gate G1.1 -- PASSES
├── test_laplace_cg_iterations.py    # ✓ [P1] Gate G1.2 (formerly G1.3) -- PASSES
├── test_laplace_sphere.py           # ✓ [P1] Gate G1.6 (formerly G1.2) -- strict xfail, see "Known gaps"
├── test_laplace_picard.py           # ✓ [P1] Regression test for solve_laplace residual_norm fix
├── test_m0_extrude.py               # ✓ [M0] Prism-split unit tests (pure numpy, no Gmsh)
├── test_m0_cylinder.py              # ✓ [M0] Cylinder-flow validation (analytic Cp, spanwise)
├── test_wall_correction_cylinder.py # ✓ [P1] Gate G1.3 -- completed, acceptance NOT met
│                                     #   (negative result locked in; acceptance = strict xfail)
├── test_m0_naca0012.py              # ✓ [M0] NACA0012 family topology/wake-sheet/ingestion
├── test_p2_wake_cut.py              # ✓ [P2] Cut topology unit tests (synthetic strip, no Gmsh),
│                                     #   G2.1 + G2.2, assert-fires-on-broken-cut, hard-rule-7
│                                     #   sweep over every wake-tagged mesh in cases/meshes/
├── test_p2_kutta_naca0012.py        # ✓ [P2] Gates G2.3/G2.4/G2.5 + V2.1–V2.5 artifacts
├── test_subsonic_*.py    # [P3] Gates G3.1–G3.3
└── test_transonic_*.py   # [P4] Gates G4.1–G4.3

artifacts/                 # Gate outputs (auto-generated, gitignored)
├── G0.1/                 # Volume conservation heatmap
├── G0.2/                 # Gradient recovery plots
├── G0.3/                 # Element coloring 3D render
└── ...

pyproject.toml            # ✓ Project metadata and dependencies
setup.py                  # ✓ Legacy setup (pyproject.toml preferred)
CLAUDE.md                 # ✓ Claude Code project instructions (doc map + workflow; imports docs/agent-rules.md)
```

## Implementation Status

### ✓ Complete (P0)
- **pyfp3d/physics/isentropic.py** — All physics scalars, numba-jitted
- **pyfp3d/mesh/reader.py** — Mesh I/O (meshio → SoA), mesh validation, tagged round-trip write
- **pyfp3d/mesh/metrics.py** — Geometry: volumes, gradients, face adjacency (Numba)
- **pyfp3d/mesh/coloring.py** — Element graph coloring for @prange
- **pyfp3d/post/vtk_out.py** — VTK writer (point fields) + PNG/CSV artifact helpers
- **tests/conftest.py** — Pytest fixtures
- **tests/test_v0_freestream.py** — Smoke tests + regression baseline ✓
- **tests/test_mesh_volume.py** — Gate G0.1 (volume conservation) ✓
- **tests/test_mesh_gradient.py** — Gate G0.2 (gradient recovery) ✓
- **tests/test_mesh_coloring.py** — Gate G0.3 (element coloring) ✓
- **tests/test_io_vtk.py** — Gate G0.4 (VTK round-trip) ✓
- **tests/test_mesh_adjacency.py**, **tests/test_mesh_reader_roundtrip.py** — regression tests for
  two bugs found by manual audit (see below)
- **pyproject.toml** — Build metadata and dependencies
- **CLAUDE.md** — Claude Code project instructions, auto-loaded each session (replaces the
  former `.copilot-instructions.md`, whose content largely duplicated design.md/roadmap.md
  and had drifted; details now live only in the authoritative docs)

Three latent bugs found by manual code audit (not caught by the existing suite, because nothing
exercised these code paths) have been fixed, each now with a regression test:
- `mesh/metrics.py::build_face_adjacency` crashed under `@njit` (reflected-list dict values are
  not valid in numba nopython mode) — rewritten around a `numba.typed.Dict` keyed by sorted face,
  storing only the packed first-owner index instead of a growing list. (`test_mesh_adjacency.py`)
- `mesh/reader.py::write_mesh` silently dropped every named boundary group (`wall`, `farfield`,
  ...), writing only a legacy `"all_triangles"` block, and `.msh` was ambiguous between meshio's
  `ansys`/`gmsh` writers (defaulted to `ansys`, discarding all tag data). Now writes each boundary
  group as its own tagged `triangle` block plus `gmsh:physical`/`field_data`, explicitly via the
  `gmsh22` writer. (`test_mesh_reader_roundtrip.py`)
- `solve/picard.py::solve_laplace` reported `residual_norm` over *all* nodes, including Dirichlet
  (far-field) rows whose natural-BC flux imbalance is O(1) and never shrinks — swamping the actual
  free-dof residual (which was already converging to ~1e-10). Now restricted to free dofs and
  correctly nets out `body_source_rhs` when present. (`test_laplace_picard.py`)

A second manual audit (2026-07-06) fixed four more latent issues, again each locked in by a
regression test:
- `tests/conftest.py::artifacts_dir` handed out a `tempfile.TemporaryDirectory`, so every gate
  PNG/CSV was deleted at teardown and the repo `artifacts/` directory stayed permanently empty —
  violating the CLAUDE.md/roadmap rule that every visual gate leaves inspectable headless
  artifacts. Now defaults to the persistent (gitignored) `artifacts/`, overridable via
  `PYFP3D_ARTIFACTS_DIR`. (`test_conftest_artifacts.py`)
- `mesh/metrics.py::element_gradients` silently returned zero gradients for |det J| < 1e-20 — an
  absolute threshold that both let coplanar tets corrupt assembly without any error and zeroed
  perfectly well-shaped tiny elements (edge ~1e-7 ⇒ det ~1e-21). Now raises `ValueError` with a
  scale-relative threshold (|det J| < 1e-12 × product of edge norms). (`test_metrics_degenerate.py`)
- `post/surface.py::_wall_vertex_normals` silently averaged cancelling normals when wall-triangle
  winding is inconsistent (or a crease is razor-sharp), producing garbage tangent planes
  downstream. Now raises with a diagnostic when |Σ area·n̂| / Σ area < 0.05 at any vertex.
  (`test_post_surface.py::test_inconsistent_wall_winding_raises`)
- `mesh/reader.py::read_mesh` rebuilt `tag_names` from `field_data` unconditionally (meshio meshes
  always *have* that attribute), clobbering the default `["bulk"]` with `[""]` for meshes carrying
  no named 3D groups. Now rebuilds only when a dim-3 entry exists.
  (`test_mesh_reader_roundtrip.py::test_untagged_volume_keeps_default_bulk_name`)

Numba pitfall hit while landing the degenerate-tet guard (captured here for the end-of-P1 skill
sedimentation checkpoint): `@njit(cache=True)` disk-cache invalidation tracks only the cached
function's *own* source, not its jitted callees — after editing `metrics.py::element_gradients`,
the stale cached `kernels/residual.py::assemble_residual` (with the old callee behavior baked in)
could still be loaded depending on compile order, making a test pass or fail depending on which
test file ran first. Symptom: behavior differs between `pytest tests/test_x.py` and the full
suite. Fix during development: delete `pyfp3d/**/__pycache__/*.nb[ic]` after editing any function
that other cached kernels call, then re-run.

Doc-only fixes in the same pass: design.md §5 far-field formula typo (`z sinα cosβ`, was garbled
"cosα-corrected"), design.md §10 V1 / roadmap G1.1 manufactured-solution wording (sin·sin·sin, not
sin·cos), `coloring.py` docstring (pure-Python preprocessing, not "@njit-compatible"; architecture
is design.md §7, not §3), `isentropic.py::validate_physics_bounds` docstring bounds aligned with
the implementation, a dead `nodal_gradient_recovery_spr` cross-reference in `post/surface.py`
replaced with the real wall-recovery functions, and stale `docs/PROJECT_STRUCTURE.md` paths (this
file lives at the repo root).

### ✓ P1 gates G1.1 and G1.2 closed; G1.3–G1.6 open (P1 renumbered 2026-07-06, see roadmap.md; below)
- **pyfp3d/kernels/residual.py** — Laplace residual (6.1) + SPD stiffness matrix (6.2) assembly
- **pyfp3d/solve/linear.py** — Dirichlet elimination (principal-submatrix reduction) + CG+PyAMG
- **pyfp3d/solve/picard.py** — `solve_laplace()` driver (P1's Picard loop degenerates to one solve)
- **pyfp3d/post/surface.py** — `nodal_gradient_recovery()` (volume-weighted, interior fields),
  `wall_tangential_gradient()` (surface-only linear recovery, for wall Cp), and
  `wall_tangential_gradient_quadratic()` (surface-only quadratic patch recovery; a real but modest
  improvement, see below)
- **tests/mesh_utils.py**, **cases/meshes/sphere_shell/** — structured-cube (MMS) and sphere-shell
  (G1.6) mesh generators; sphere-shell coarse/medium `.msh` + inspection PNGs are committed
- **tests/test_laplace_mms.py** — Gate G1.1 (MMS convergence) ✓ — L2 slope ≈ 1.94–1.96 with a
  sin·cos manufactured solution and a proper 4-point quadrature-consistent load vector. (A
  harmonic-polynomial exact solution was tried first and rejected: this codebase's structured
  Kuhn-triangulated cube reproduces harmonic quadratics to machine precision at *every* h, giving
  zero convergence-order signal — the same reason central finite differences are exact for
  quadratics.)
- **tests/test_laplace_cg_iterations.py** — Gate G1.2 (formerly G1.3; CG+AMG mesh-independence) ✓ — iterations
  8→11→14 across an 8×/level node-count increase (n=8,16,32 cube), comfortably under a 2× cap.

**G1.6 (formerly G1.2; incompressible sphere Cp) is still open** — `tests/test_laplace_sphere.py::test_sphere_cp_medium_mesh`
is a `strict=True` xfail against the real <2% criterion, not a loosened threshold:
- The original `nodal_gradient_recovery` (volume-weighted average of the one-sided tets touching
  each wall node) gave ~26% max / 9% mean Cp error on the medium mesh — systematically low,
  because tangential velocity physically decays moving away from the wall and every incident tet
  sits on the inward side.
- A local least-squares ("SPR-style") extrapolation fix was tried and **rejected**: on a real
  graded mesh, some nodes' 1-ring element patches are nearly coplanar, and extrapolating a linear
  model through an ill-conditioned patch blew one node's Cp error up to 429 (found on the medium
  mesh). Don't resurrect this approach without a robust conditioning safeguard.
- The fix that shipped, `wall_tangential_gradient()`, computes the surface-tangential gradient
  directly from the wall's own P1 triangulation (physically exact for a natural-BC wall, where the
  normal component is zero by construction) instead of extrapolating the volume gradient. This
  improved medium-mesh max error to ~12%, and is a genuine, well-conditioned fix (no
  extrapolation, only interpolation/averaging) — but a mesh-refinement sweep (h_min = 0.08 → 0.04
  → 0.02 → 0.015, up to ~1.2M nodes / 7.4M tets) only reaches ~3.6% max error and visibly
  saturates rather than continuing to converge, so the remaining gap is *not* simply "refine the
  mesh more."
- **Root-caused (this session)** by isolating each candidate error source on a *clean*,
  single-variable h_min refinement sweep (h_min = 0.08 → 0.05 → 0.03 → 0.02, everything else in
  `generate_sphere_shell.py` held fixed — the earlier sweep above changed h_max/r_out/dist_max
  simultaneously with h_min, confounding the picture):
  - **Recovery scheme ruled out as the dominant cause.** An oracle test feeds the *exact* analytic
    potential straight into the recovery step, bypassing the FEM solve entirely. The recovery
    operator's own bias measured this way is 9.3%→0.7% (linear) or 0.5%→0.005% (quadratic, see
    below) across the same h_min range where the *full* FEM pipeline gives 12%→4.2% — i.e.
    recovery is a small fraction of the total error at every mesh size tested, for both schemes.
  - **Under-refined bulk/far mesh is a minor contributor, not the dominant one.** Tightening
    h_max/dist_max at a fixed h_min=0.03 helps a bit (5.6% → 4.3% max error) but plateaus by
    ~4% even with a much finer far mesh (0.5 vs 3.0 h_max, 4.5M vs 1.4M tets) — refining the bulk
    mesh alone doesn't close the gap either.
  - **Confirmed dominant cause: the volume PDE solve's own accuracy next to the wall**, not
    anything in post-processing. The raw nodal potential φ itself (not just its recovered
    gradient) has the same sub-first-order error at the wall as the derived Cp, and the
    convergence *order* measured on the clean sweep is decreasing as h shrinks (0.88 → 0.56 → 0.42
    for max nodal φ error) rather than settling to a fixed rate — a signature of a genuine
    geometric/consistency error, not plain discretization error. Mechanism: the natural
    (zero-flux) BC is satisfied on the flat *polyhedral* wall-facet approximation (Γ_h), not the
    true curved sphere (Γ); solving `Δφ_h = 0` on the (slightly wrong) domain bounded by Γ_h
    pollutes the whole solution through ellipticity, not just the boundary nodes. This is a
    textbook "variational crime" for Neumann conditions on curved boundaries meshed with flat
    (non-isoparametric) elements.
  - **A direct fix was tried and rejected**: a Nitsche/penalty term added to the stiffness matrix,
    weakly forcing each wall-adjacent tet's own volumetric gradient toward zero along the *true*
    (here, analytically known) surface normal, swept over penalty strength β. Result: error and CG
    iteration count both got *worse* monotonically with increasing β (e.g. medium mesh max error
    12%→17%→40%→98%→211% for β=1→10→100→1000→1e4). Diagnosis: a P1 tet spanning from the wall
    inward necessarily has a nonzero radial gradient component representing the interior falloff
    of tangential velocity (the exact solution's normal derivative is zero only exactly *at* the
    wall, not throughout the adjacent tet's finite thickness) — that's correct FEM behavior, not a
    BC violation, so this penalty fights the physically-correct solution instead of correcting an
    inconsistency. Don't resurrect this approach; a correct fix needs to change how the boundary
    integral/geometry itself is represented (see below), not add a volumetric penalty.
  - **Implemented as a genuine (if modest) improvement**: `wall_tangential_gradient_quadratic()`
    fits a local quadratic model per wall node (in its own reconstructed tangent plane, over its
    1-ring — expanded to 2-ring, then falling back to a 2-parameter linear fit, if the patch is
    rank-deficient for the 6-parameter fit; uses `np.linalg.lstsq`'s SVD-based minimum-norm
    solution throughout, never a normal-equations solve or 3D extrapolation, so it can't hit the
    ill-conditioned blowup the earlier volume-based SPR attempt did). This is exact for a locally
    quadratic field (vs. linear recovery's exactness only for locally linear fields) and cuts the
    recovery-only oracle error by roughly 20x, but — consistent with recovery not being the
    dominant error source — only trims medium-mesh total error from ~12.0% to ~11.6%. Adopted as
    the default for the G1.6 test since it's a strict, low-risk improvement, but it does not (and
    was never going to) close the gate alone. See `tests/test_post_surface.py` for regression
    coverage locking in both facts (recovery-only accuracy, and the fact that it's still not
    enough).
- **Fix routes researched and tiered (2026-07-06, see design.md §5.1 for the full writeup)**:
  verify **Option A** first — the true-normal weak-flux correction (lagged/SBM-style: correct the
  natural-BC *data* via closest-point projection onto the true surface; RHS-only, stiffness matrix
  unchanged, AMG hierarchy fully reused). An oracle experiment (exact analytic potential + exact
  normals in the correction RHS, gate G1.4) measures the accuracy ceiling of this
  first-order normal correction, then the DP1 decision point picks the route:
  Option B (Gap-SBM gap terms) if the ceiling lands at 2–5%, or curved/isoparametric wall
  elements — now demoted to the option of last resort, taken up as its own scoped effort only if
  both A and B fall short — with G1.6 redefined per Option C's geometry-consistent-reference
  yardstick in that case. Still do **not** propose further h-refinement or recovery-scheme
  tweaks; both remain ruled out with evidence. `cases/meshes/cylinder_2.5d/` is confirmed to
  share the same root cause with quantified error (max |Cp err| 0.091 coarse → 0.045 medium,
  ~O(h)) and is the designated fix-route testbed — see design.md §5.1.1 / gate G1.3; the
  gate criterion itself stays on the sphere.

### ✓ G1.3 + G1.4 completed 2026-07-06 — negative results; DP1 "> 5%" branch taken

The Option A (true-normal weak-flux correction) verification chain ran to completion the same
day it was renumbered, and falsified the route (full evidence: roadmap G1.3/G1.4/DP1 entries,
design.md §5.1.2, `artifacts/G1.3/` + oracle results in `cases/demo/p1_laplace/results/`):

- **Delivered code**: `pyfp3d/solve/wall_correction.py` (closest_point_normal callback interface
  with analytic cylinder/sphere implementations, domain-outward facet orientation from the
  owning tet, 3-point edge-midpoint facet quadrature, RHS assembly verified against a
  hand-computed single-facet case); `pyfp3d/post/section_cut.py` (P2's final interface,
  degenerate single-layer path); `tests/test_wall_correction_cylinder.py` (10 tests);
  the sphere-oracle experiment (now absorbed into `cases/demo/p1_laplace/run_demo.py`,
  2026-07-07); cylinder `fine.msh` (50.2k tets); shared cylinder-case helpers moved into
  `tests/mesh_utils.py`.
- **Core finding**: for a harmonic potential with body-fitted wall vertices, the exact net flux
  through every flat facet is exactly zero (divergence theorem over the facet/true-surface
  sliver), so boundary-DATA corrections have (almost) nothing to correct. Measured: the full
  consistency defect is ~2e-5 (coarse cylinder, ~O(h⁴)); Option A's t-form assembles to machine
  zero on the cylinder; corrected Cp errors unchanged. Sphere ceiling: 0.1156 → 0.1133 (11.3%)
  vs the < 2% target → DP1 "> 5%" branch: Option C gate re-spec + separately-scoped curved
  elements. SBM/BDT-style corrections earn their keep on *unfitted* O(h)-gap boundaries, not
  here.
- **Cylinder re-framed**: ~76% of its Cp error at every level is the surface *recovery* on the
  quasi-2D sliver-strip wall triangulation (exact-potential-through-recovery oracle), and wall
  nodal φ converges at a healthy ~1.2 order — no sphere-style sub-first-order pathology. The
  "same variational crime" designation is retracted; the case remains the M0 end-to-end check.
- **G2.5(b) reference numbers** (by-product): solved-field spanwise floor max |w|/U∞ =
  2.88e-2 / 1.50e-2 / 7.40e-3 over coarse/medium/fine, identical with/without correction.

### ✓ M0 mesh-side items delivered (2026-07-06)

- **pyfp3d/meshgen/extrude.py** — single-layer quasi-2D extrusion: one cell layer in z, both
  planes tagged `symmetry`; every prism split into 3 tets with the globally consistent
  min-global-index diagonal rule (Dompierre et al. "indirect" subdivision — since top-layer
  indices are bottom + n_2d, the smallest index on any lateral quad is always the smaller
  *bottom* node, so both sides of every shared quad pick the same diagonal by construction);
  tagged interior edge sheets (the wake) split by the same rule so they coincide exactly with
  tet faces. `assert_quad_split_consistency()` is the M0 preprocessor assert: single-owner tet
  faces must equal the tagged boundary set, interior-sheet faces must have exactly two owner
  tets; unit tests prove it fires on a deliberately broken split (`tests/test_m0_extrude.py`).
- **pyfp3d/meshgen/planar.py** — vanilla-Gmsh 2D builders (Distance+Threshold grading like the
  sphere case): `cylinder_annulus_2d()` and `naca0012_wake_2d()` (closed-TE NACA0012, circular
  far field r=15c centered at mid-chord, wake line TE→farfield embedded with
  `gmsh.model.mesh.embed` so triangle edges conform to it). Gmsh imported lazily.
- **cases/meshes/naca0012_2.5d/** (the M0 deliverable) — `generate_naca0012.py`, one parameter
  (h_wall) per level, everything else derived: coarse 16.4k / medium 61.8k tets committed
  (on-target vs the ~15k/60k spec; fine ~240k on demand). Per-level stats CSV + layer-inspection
  PNG written at generation time (headless, roadmap §0.1).
- **cases/meshes/cylinder_2.5d/** (extra validation case) — quasi-2D circular-cylinder flow,
  no wake (non-lifting), analytic solution phi = x(1 + a²/r²), Cp = 1 − 4 sin²θ. The simplest
  end-to-end check of pipeline + P1 solver: measured max |Cp err| 9.1% (coarse, 6.9k tets) →
  4.5% (medium, 17.3k tets) with quadratic surface recovery — same curved-wall/flat-facet
  variational crime as the G1.6 sphere, converging ~O(h), fully expected on this geometry.
- **tests/test_m0_extrude.py / test_m0_cylinder.py / test_m0_naca0012.py** — 21 tests covering
  the M0 gate items: reader ingestion, tags, quad-split consistency, wake-sheet topology (one
  connected planar interior sheet TE→farfield, both z-planes, nodes not duplicated), symmetry
  planes planar/disjoint from wall, cylinder Cp vs analytic, spanwise-gradient behavior.

**G2.5(b) finding (evidence in roadmap.md G2.5 note — the gate criterion needs re-spec before
P2 closes it):** the *interpolated* freestream phi = x is spanwise-zero to machine precision
(G2.5(a) verified), but the *solved* field is not and cannot be: a 3-tet prism split is
necessarily asymmetric under the z-mirror (on a lateral quad face, ∫N_i dS is S/3 or S/6
depending on the diagonal direction), so a z-invariant field does not satisfy the discrete
equations and the discrete minimizer carries O(h) spanwise velocity noise. Measured (non-lifting
solves): cylinder max|w|/U∞ 2.9e-2 (coarse) → 1.5e-2 (medium); NACA0012 α=0 5.4e-2 (coarse) —
clean ~O(h) decay, 10 orders above the literal 1e-12 criterion. Literal machine-zero would need
a z-mirror-symmetric subdivision (requires Steiner points, violating the 3-tet M0 spec). The
tests assert the honest behavior instead: machine-zero for the interpolant, small-and-decreasing
under refinement for solved fields.

### ✓ P2 delivered, M0 closed (2026-07-06)

- **pyfp3d/mesh/wake_cut.py** — wake-sheet node duplication: ⁺-side (upper-hint) copies
  appended after the original nodes, so the reduced dof space of the wake constraint is
  exactly the original node set; per-node flood-fill side classification (adjacency through
  non-wake faces only — no planarity assumption, ready for the M1 swept wake); spanwise
  stations group TE nodes by (x, y), so a quasi-2D extrusion collapses to the single scalar
  Γ of the M0 spec while a swept TE gets per-node stations; Kutta probe nodes (wall, one
  edge off each TE node, per side); preprocess-time topology asserts (roadmap P2 list,
  proven to fire on a deliberately broken cut).
- **THE P2 spec deviation (evidence-backed): TE nodes ARE duplicated.** The roadmap
  originally asserted the opposite; implemented first and measured: a single-valued TE node
  tapers [φ] from Γ to 0 across the first wake cell ≡ a point vortex of strength Γ parked
  at the TE → wall suction ~ (Γ/2πr)², a spurious TE force ~ Γ²/h that *diverges* under
  refinement (coarse NACA0012, Γ=0.3 prescribed: peak wall |V| 4.6 U∞, −0.27 spurious out
  of cl ≈ 0.6, from 6 triangles). With the TE doubled: cl = 0.6012 vs Kutta–Joukowski 0.6.
  Roadmap assert block re-specced accordingly; design.md §4 records the theory (the TE jump
  IS the Kutta condition).
- **pyfp3d/constraints/wake.py** — master–slave elimination: φ_full = T φ_red + g(Γ),
  A_red = TᵀAT assembled once (SPD preserved), Γ RHS-only via precomputed per-station
  vectors h_j = TᵀA g_j; folding slave rows into masters is exactly the weak mass-flux
  continuity (4.2). `kutta_targets()` = per-station mean of probe jumps (also filters the
  O(h) spanwise noise on quasi-2D meshes).
- **pyfp3d/constraints/dirichlet.py** — far-field freestream at incidence + incompressible
  2D point-vortex correction with the branch cut ON the wake sheet: an eliminated ⁺-side
  far-field wake node (master + Γ) automatically equals the upper-branch vortex value, so
  no special Dirichlet casing is needed.
- **solve/picard.py::solve_laplace_lifting()** — Kutta outer loop with the matrix, Dirichlet
  split and AMG hierarchy built once; Γ updates are RHS-only. Secant (Aitken) acceleration
  from the second update: the linear map's measured slope is b ≈ 0.93, so plain ω-relaxation
  would need O(100) updates; the secant hits the affine fixed point in 2.
- **post/surface.py** — triangle-wise wall force integration (`wall_force_coefficients`):
  in-plane tangential gradient IS the wall velocity (natural BC), Cp per triangle, outward
  normals oriented by the owning tet (no winding/star-shape assumptions), no nodal averaging
  across the sharp-TE crease; `sectional_cl_from_gamma` (KJ cross-check).
- **post/section_cut.py** — general z = const marching-tets path (linear-exact, unit-tested
  at an off-node plane) + `wall_cp_curve()`: triangle-wise sectional Cp(x/c) split
  upper/lower at the intersection-segment midpoints.
- **cases/reference_data/naca0012_incompressible/** — Hess–Smith panel reference (constant
  sources + single vortex + Kutta), same closed-TE coordinate set as the mesh so G2.3 is
  method-vs-method; cl(4°) = 0.482556 at N=800, Cp-integration vs Kutta–Joukowski agree to
  0.09%, lift slope 6.91/rad vs thickness-corrected 6.90. Provenance in its README.
- **Gates** (all green, `tests/test_p2_wake_cut.py` + `tests/test_p2_kutta_naca0012.py`,
  artifacts/G2.{1..5}/): G2.1 ‖R‖∞ = 8.4e-13 (folded wake-master rows 6.9e-16); G2.2
  [φ] − Γ < 1e-13; G2.3 medium cl_p = 0.47858 → −0.82% vs panel (coarse −3.0%), Kutta
  converged in 2 updates; G2.4 Γ-cl vs pressure-cl 0.01% (coarse 0.32%); G2.5 closed under
  the re-specced criterion (b) — p99 |w|/U∞ 4.82e-3 → 2.35e-3 (ratio 2.05 at h ratio 2),
  max recorded (LE peak-gradient region, not wake), stripe-free mid-plane heatmap.
- **M0 closed** with it: wake-cut topology asserts sweep every wake-tagged mesh in
  cases/meshes/ (hard rule 7 test) and the G2.5 acceptance link is green.

### ⏳ Next
- **P3 (subsonic compressible)** — Picard density outer loop + PG-scaled vortex far field;
  retire the P1 assembly tech debt first (precompute B_e/V_e, colored `prange` assembly per
  design.md §7; the outer loop makes assembly hot).
- ~~G1.3/G1.4 oracle experiments~~ — DONE 2026-07-06 with negative results (see the G1.3+G1.4
  section above); DP1 decided the "> 5%" branch.
- **G1.6 re-spec per Option C** — draft the geometry-consistent-reference acceptance criterion
  (design.md §5.1 Option C): compare against a high-accuracy reference on the *same polyhedral
  domain* (BEM or ultra-fine reference), separating geometric model error from code error. This
  is now the open P1 item; the curved/isoparametric-element effort for physical accuracy is
  separately scoped (backlog).

## Quick Start

### 1. Install dependencies
```bash
cd /home/lrz/code/UP3D
pip install -e ".[dev]"
```

### 2. Run smoke tests
```bash
pytest tests/test_v0_freestream.py -xvs
```

Expected output:
```
tests/test_v0_freestream.py::test_import_pyfp3d PASSED
tests/test_v0_freestream.py::test_import_physics PASSED
tests/test_v0_freestream.py::test_isentropic_stagnation PASSED
tests/test_v0_freestream.py::test_isentropic_freestream PASSED
tests/test_v0_freestream.py::test_pressure_coefficient_bounds PASSED
```

### 3. Run physics module directly (self-test)
```bash
python -m pyfp3d.physics.isentropic
```

Expected output:
```
=== Isentropic Physics Self-Test ===

Freestream Mach M∞ = 0.5

At stagnation (q² = 0):
  ρ = 1.000000 (expected 1.0)
  Cp = 1.000000 (expected 1.0)

At freestream (q² = 1):
  ρ = 1.000000
  M = 0.500000 (expected 0.5)
  Cp = 0.000000 (expected 0.0)

Critical speed q*² = 0.923077 where M = 1.0:
  Computed M = 1.000000 (expected 1.0)

✓ All checks passed!
```

## Next Steps (P1 Completion)

P0 (mesh I/O, metrics, coloring, VTK writer, gates G0.1–G0.4) is done. G1.1 and G1.2 (formerly
G1.3) are done; G1.3 and G1.4 completed 2026-07-06 with negative results and DP1 took the
"> 5%" branch (see above). To close P1, G1.6 remains, pending its Option C re-spec:

1. ~~**Run the G1.3 cylinder oracle pre-study and the G1.4 sphere oracle experiment, then follow
   the DP1 decision point**~~ — DONE 2026-07-06, negative result: the Option A ceiling is
   ≈ 11.3% (medium sphere) because on body-fitted meshes the boundary-data defect it corrects is
   (near-)zero — see the "G1.3 + G1.4 completed" section above and design.md §5.1.2. DP1 took
   the "> 5%" branch. The new open item: **draft the G1.6 Option C acceptance re-spec**
   (geometry-consistent reference on the same polyhedral domain); the curved/isoparametric-
   element effort for physical accuracy is separately scoped.

2. ~~**Create test meshes** (M0, parallel track)~~ — DONE 2026-07-06 (single-layer re-spec
   targets ~15k/60k/240k tets, not the older 30k/150k/700k): `pyfp3d/meshgen/` +
   `cases/meshes/naca0012_2.5d/` + `cases/meshes/cylinder_2.5d/`, topology validated by
   `tests/test_m0_*.py`. See "M0 mesh-side items delivered" above.

## References

- **Design & Theory:** [docs/design.md](docs/design.md)
- **Roadmap & Gates:** [docs/roadmap.md](docs/roadmap.md)
- **Agent Rules:** [docs/agent-rules.md](docs/agent-rules.md)
- **Claude Code Instructions:** [CLAUDE.md](CLAUDE.md)

---

**Last updated:** 2026-07-06  
**Status:** P2 closed 2026-07-06 (wake cut / circulation / Kutta on the M0 mesh: gates
G2.1–G2.5 green — medium-mesh cl −0.82% vs the Hess–Smith panel reference, Γ-vs-pressure lift
cross-check 0.01%, Kutta converged in 2 secant updates; one evidence-backed spec deviation:
TE nodes ARE duplicated, see the P2 section above), and M0 closed with it (topology asserts
sweep the whole mesh family; G2.5 criterion (b) re-specced to p99-of-|w| ≥ 1st-order decay,
measured ratio 2.05). **Next: P3 (subsonic compressible).** P0 gates G0.1–G0.4 green; the full coarse regression suite now runs against real case
meshes (cylinder + NACA0012); P0 closed 2026-07-06 (G0.1–G0.4 ticked in the roadmap). P1 (gates renumbered
in workflow order 2026-07-06, old → new: G1.3 → G1.2, G1.2-a0 → G1.3, G1.2-a → G1.4,
G1.2-c → DP1, G1.2-b → G1.5, G1.2 → G1.6): G1.1 (MMS) and G1.2 (CG+AMG mesh-independence)
closed; G1.3 (cylinder oracle) and G1.4 (sphere oracle ceiling) completed 2026-07-06 with
negative results — boundary-data corrections have (near-)zero lever on body-fitted meshes
(exact per-facet flux is zero by the divergence theorem), Option A ceiling ≈ 11.3% — so DP1
took the "> 5%" branch and G1.5 is void; G1.6 (sphere Cp) stays open with a `strict=True`
xfail tracking the real 2% criterion, its Option C re-spec (geometry-consistent reference) is
the open P1 item, and the curved/isoparametric-element effort is separately scoped. The M0
cylinder case was re-framed by G1.3: its ~O(h) wall-Cp error (9.1% → 4.5% → 2.2% max) is
~76% surface-recovery on the quasi-2D sliver strip, not the sphere pathology.
89 tests total (87 passed, 2 xfailed), full suite ~18 s.
