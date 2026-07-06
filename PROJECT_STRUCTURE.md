# Project Structure

This document describes the directory layout and initialization status of pyFP3D.

## Directory Tree

```
pyfp3d/                    # Main package
├── __init__.py           # Package entry point, NOJIT mode flag
├── mesh/                 # Mesh I/O, topology, metrics, coloring
│   ├── __init__.py
│   ├── reader.py         # [P0] meshio → SoA arrays + boundary tags
│   ├── metrics.py        # [P0] volumes, gradients, face adjacency
│   ├── coloring.py       # [P0] element graph coloring
│   └── wake_cut.py       # [P2] node duplication from wake surface
├── physics/              # Physics constants and constitutive relations
│   ├── __init__.py
│   └── isentropic.py     # ✓ [P0] ρ(q²), M(q²), a(q²), Cp, etc. (complete)
├── kernels/              # Element-wise assembly kernels (Numba-jitted)
│   ├── __init__.py
│   ├── residual.py       # [P1] Laplace residual + stiffness assembly (done) → [P3] isentropic → [P6] Newton
│   └── upwind.py         # [P4] Artificial compressibility
├── solve/                # Linear and nonlinear solvers
│   ├── __init__.py
│   ├── linear.py         # [P1] Dirichlet elimination + CG/PyAMG preconditioner (done)
│   ├── picard.py         # [P1] Laplace driver (single linear solve, no outer loop) (done)
│   └── newton.py         # [P6] Newton method
└── post/                 # Post-processing
    ├── __init__.py
    ├── vtk_out.py        # [P0] Write .vtu for ParaView; also the PNG/CSV gate-artifact helpers
    │                       #      (export_error_heatmap, export_matplotlib_plot) live here, not
    │                       #      in a separate artifacts.py
    └── surface.py        # [P1] nodal_gradient_recovery() (volume-weighted, for interior fields)
                            #      and wall_tangential_gradient() (surface-only, for wall Cp --
                            #      the accurate one; see "Known gaps" for why it still isn't
                            #      accurate *enough* to close G1.2)

cases/                     # Test cases and reference data
├── meshes/               # Mesh families (coarse/medium/fine)
│   ├── sphere_shell/     # [P1] Gmsh sphere-shell case for gate G1.2 (coarse/medium generated)
│   ├── naca0012_2.5d/    # [M0] Extruded NACA0012 + wake surface -- not started
│   └── onera_m6/         # [M1] Swept wing -- not started
├── reference_data/       # Ground truth (DO NOT EDIT)
│   ├── sphere_incomp_cp.csv
│   ├── naca0012_cl.csv
│   └── ...
└── test_*.py             # [Deprecated] Integration tests (use tests/ now)

tests/                     # Unit and gate tests
├── conftest.py           # ✓ Pytest fixtures: artifacts_dir, mesh_dir, etc.
├── mesh_utils.py         # ✓ [P1] Dependency-free structured-cube + sphere-shell mesh generators
├── __init__.py
├── test_v0_freestream.py # ✓ [P0/P1] Primary regression test (incl. cut-free residual check)
├── test_mesh_*.py        # [P0] Gates G0.1–G0.4
├── test_mesh_adjacency.py           # ✓ [P0] Regression test for build_face_adjacency fix
├── test_mesh_reader_roundtrip.py    # ✓ [P0] Regression test for write_mesh tag-loss fix
├── test_laplace_mms.py              # ✓ [P1] Gate G1.1 -- PASSES
├── test_laplace_cg_iterations.py    # ✓ [P1] Gate G1.3 -- PASSES
├── test_laplace_sphere.py           # ✓ [P1] Gate G1.2 -- strict xfail, see "Known gaps"
├── test_laplace_picard.py           # ✓ [P1] Regression test for solve_laplace residual_norm fix
├── test_wake_*.py        # [P2] Gates G2.1–G2.4
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

### ✓ P1 gates G1.1 and G1.3 closed; G1.2 open (see below)
- **pyfp3d/kernels/residual.py** — Laplace residual (6.1) + SPD stiffness matrix (6.2) assembly
- **pyfp3d/solve/linear.py** — Dirichlet elimination (principal-submatrix reduction) + CG+PyAMG
- **pyfp3d/solve/picard.py** — `solve_laplace()` driver (P1's Picard loop degenerates to one solve)
- **pyfp3d/post/surface.py** — `nodal_gradient_recovery()` (volume-weighted, interior fields),
  `wall_tangential_gradient()` (surface-only linear recovery, for wall Cp), and
  `wall_tangential_gradient_quadratic()` (surface-only quadratic patch recovery; a real but modest
  improvement, see below)
- **tests/mesh_utils.py**, **cases/meshes/sphere_shell/** — structured-cube (MMS) and sphere-shell
  (G1.2) mesh generators; sphere-shell coarse/medium `.msh` + inspection PNGs are committed
- **tests/test_laplace_mms.py** — Gate G1.1 (MMS convergence) ✓ — L2 slope ≈ 1.94–1.96 with a
  sin·cos manufactured solution and a proper 4-point quadrature-consistent load vector. (A
  harmonic-polynomial exact solution was tried first and rejected: this codebase's structured
  Kuhn-triangulated cube reproduces harmonic quadratics to machine precision at *every* h, giving
  zero convergence-order signal — the same reason central finite differences are exact for
  quadratics.)
- **tests/test_laplace_cg_iterations.py** — Gate G1.3 (CG+AMG mesh-independence) ✓ — iterations
  8→11→14 across an 8×/level node-count increase (n=8,16,32 cube), comfortably under a 2× cap.

**G1.2 (incompressible sphere Cp) is still open** — `tests/test_laplace_sphere.py::test_sphere_cp_medium_mesh`
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
    the default for the G1.2 test since it's a strict, low-risk improvement, but it does not (and
    was never going to) close the gate alone. See `tests/test_post_surface.py` for regression
    coverage locking in both facts (recovery-only accuracy, and the fact that it's still not
    enough).
- **What would actually close this gate**: genuine curved/isoparametric boundary elements (curve
  the geometric mapping of wall-adjacent tets to match the true surface, so the natural BC's
  implicit "do nothing" trick is applied on Γ instead of Γ_h) — a properly-derived shape/geometry
  correction, not a bolt-on penalty. This is a substantially larger, separately-scoped effort
  (touches `mesh/metrics.py` element-Jacobian machinery and `kernels/residual.py` assembly, not
  just `post/surface.py`) and should get its own design pass before implementation, rather than
  being rushed in alongside other P1/P2 work.

### ⏳ Next
- Scope a curved/isoparametric wall-boundary treatment as its own design item (see root-cause
  writeup above) — this is the concrete next step to actually close G1.2, not further h-refinement
  or more post-processing tweaks (both now ruled out as sufficient, with evidence).
- Generate NACA0012 mesh family (M0 gate) — not started.

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

P0 (mesh I/O, metrics, coloring, VTK writer, gates G0.1–G0.4) is done. G1.1 and G1.3 are done. To
close P1, only G1.2 remains:

1. **Implement a curved/isoparametric wall-boundary treatment** (see "Known gaps" above for the
   full root-cause investigation) — mesh refinement and surface-recovery improvements (both tried
   this session; see below) are now ruled out as sufficient, with evidence. The confirmed cause is
   a geometric/variational-crime inconsistency from enforcing the natural BC on the flat
   polyhedral wall approximation instead of the true curved sphere; closing G1.2 needs the
   boundary geometry itself represented correctly (isoparametric wall elements or an equivalent
   shape correction in `mesh/metrics.py`/`kernels/residual.py`), which deserves its own design
   pass before implementation.

2. **Create test meshes** (M0, parallel track)
   - Gmsh script: extruded NACA0012 + embedded wake surface
   - Generate coarse (30k), medium (150k), fine (700k) families
   - Validate topology with asserts

## References

- **Design & Theory:** [docs/design.md](../docs/design.md)
- **Roadmap & Gates:** [docs/roadmap.md](../docs/roadmap.md)
- **Agent Rules:** [docs/agent-rules.md](../docs/agent-rules.md)
- **These Instructions:** [.copilot-instructions.md](../.copilot-instructions.md)

---

**Last updated:** 2026-07-06  
**Status:** P0 closed (G0.1–G0.4 green). P1: G1.1 (MMS) and G1.3 (CG+AMG mesh-independence) closed;
G1.2 (sphere Cp) open with a `strict=True` xfail tracking the real 2% criterion. The saturation
cause is now root-caused (not just hypothesized) — a geometric/variational-crime inconsistency
from the natural BC being enforced on the flat polyhedral wall approximation rather than the true
curved sphere; recovery-scheme quality and bulk-mesh refinement were both tested and ruled out as
the dominant cause. A quadratic surface-patch recovery (`wall_tangential_gradient_quadratic`) was
implemented as a genuine but modest improvement (~12.0%→~11.6% max error on medium); closing the
gate for real needs curved/isoparametric wall elements — see "Known gaps" and "Next Steps" above.
45 tests total (44 passed, 1 xfailed), full suite ~10s.
